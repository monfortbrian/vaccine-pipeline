import os
import time
import logging
import requests as http_requests
from typing import Optional
from dataclasses import dataclass
from functools import lru_cache

from fastapi import Depends, HTTPException, Header, status
from jose import jwt, JWTError, ExpiredSignatureError

logger = logging.getLogger("kozi.auth")

# ─── Config ───────────────────────────────────────────────────────────────────
SUPABASE_JWT_SECRET = os.getenv("SUPABASE_JWT_SECRET", "")
SUPABASE_URL        = os.getenv("SUPABASE_URL", "")

if not SUPABASE_URL:
    logger.warning(
        "SUPABASE_URL is not set. ES256 token verification will fail. "
        "Add it in Railway: Settings → Variables → SUPABASE_URL"
    )

if not SUPABASE_JWT_SECRET:
    logger.warning(
        "SUPABASE_JWT_SECRET is not set. HS256 token verification will fail. "
        "Add it in Railway: Settings → Variables → SUPABASE_JWT_SECRET"
    )

# Custom inactivity session cap, must match INACTIVITY_LIMIT_MS in auth-provider.tsx
MAX_SESSION_AGE = 3600  # 1 hour (seconds)


# ─── Models ───────────────────────────────────────────────────────────────────
@dataclass
class UserClaims:
    sub: str                     # Supabase user UUID
    email: Optional[str] = None
    role: Optional[str] = None   # "authenticated"
    aal:  Optional[str] = None


# ─── JWKS fetch (cached — Supabase keys rotate rarely) ───────────────────────
@lru_cache(maxsize=1)
def _fetch_jwks() -> dict:
    """Fetch Supabase's public JWKS. Cached for process lifetime."""
    url = f"{SUPABASE_URL}/auth/v1/.well-known/jwks.json"
    try:
        resp = http_requests.get(url, timeout=10)
        resp.raise_for_status()
        logger.info(f"JWKS fetched from {url}")
        return resp.json()
    except Exception as e:
        logger.error(f"Failed to fetch JWKS from Supabase: {e}")
        raise


def _get_jwk_for_token(token: str) -> Optional[dict]:
    """Extract kid from token header, find matching key in JWKS."""
    try:
        header = jwt.get_unverified_header(token)
        kid = header.get("kid")
        alg = header.get("alg")
        logger.debug(f"Token header: alg={alg}, kid={kid}")
    except Exception as e:
        logger.warning(f"Could not read token header: {e}")
        return None

    if not kid:
        return None

    jwks = _fetch_jwks()
    for key in jwks.get("keys", []):
        if key.get("kid") == kid:
            return key
    return None


# ─── Core verification ────────────────────────────────────────────────────────
def _verify_token(token: str) -> UserClaims:
    try:
        header = jwt.get_unverified_header(token)
        alg = header.get("alg", "HS256")
    except JWTError as e:
        logger.warning(f"Could not parse token header: {e}")
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or malformed token",
            headers={"WWW-Authenticate": "Bearer"},
        )

    try:
        if alg == "ES256":
            # New Supabase ECC signing — verify against JWKS public key
            if not SUPABASE_URL:
                raise HTTPException(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail="Server misconfiguration: SUPABASE_URL not set",
                )
            jwk = _get_jwk_for_token(token)
            if not jwk:
                logger.warning("No matching JWK found for token kid")
                raise HTTPException(
                    status_code=status.HTTP_401_UNAUTHORIZED,
                    detail="Invalid or malformed token",
                    headers={"WWW-Authenticate": "Bearer"},
                )
            payload = jwt.decode(
                token,
                jwk,
                algorithms=["ES256"],
                options={
                    "require": ["sub", "exp", "iat"],
                    "verify_aud": False,
                },
            )

        else:
            # Legacy HS256 — verify against shared secret
            if not SUPABASE_JWT_SECRET:
                raise HTTPException(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail="Server misconfiguration: SUPABASE_JWT_SECRET not set",
                )
            payload = jwt.decode(
                token,
                SUPABASE_JWT_SECRET,
                algorithms=["HS256"],
                options={
                    "require": ["sub", "exp", "iat"],
                    "verify_aud": False,
                },
            )

    except ExpiredSignatureError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Session expired, please log in again",
            headers={"WWW-Authenticate": "Bearer"},
        )
    except HTTPException:
        raise
    except JWTError as e:
        logger.warning(f"JWT verification failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or malformed token",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Custom session age check via iat
    iat = payload.get("iat")
    if iat:
        age = time.time() - iat
        if age > MAX_SESSION_AGE:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Session expired, please log in again",
                headers={"WWW-Authenticate": "Bearer"},
            )

    return UserClaims(
        sub=payload["sub"],
        email=payload.get("email"),
        role=payload.get("role"),
        aal=payload.get("aal"),
    )


# ─── Helpers ──────────────────────────────────────────────────────────────────
def _extract_bearer(authorization: Optional[str]) -> Optional[str]:
    if not authorization:
        return None
    parts = authorization.split(" ", 1)
    if len(parts) != 2 or parts[0].lower() != "bearer":
        return None
    return parts[1].strip() or None


# ─── FastAPI dependencies ─────────────────────────────────────────────────────
async def require_user(
    authorization: Optional[str] = Header(None),
) -> UserClaims:
    """Protected endpoint, returns 401 if token missing or invalid."""
    token = _extract_bearer(authorization)
    if not token:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing Authorization header",
            headers={"WWW-Authenticate": "Bearer"},
        )
    return _verify_token(token)


async def optional_user(
    authorization: Optional[str] = Header(None),
) -> Optional[UserClaims]:
    """Optional auth, returns None if no/bad token, never raises."""
    token = _extract_bearer(authorization)
    if not token:
        return None
    try:
        return _verify_token(token)
    except HTTPException:
        return None