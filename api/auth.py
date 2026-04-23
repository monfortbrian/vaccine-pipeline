import os
import time
import logging
import requests
from typing import Optional
from dataclasses import dataclass

from fastapi import Depends, HTTPException, Header, status
from jose import jwt
from jose.exceptions import JWTError, ExpiredSignatureError

logger = logging.getLogger("kozi.auth")

# ─── Config ───────────────────────────────────────────────────────────────────
SUPABASE_URL = os.getenv("SUPABASE_URL", "").rstrip("/")
SUPABASE_KEY = os.getenv("SUPABASE_KEY", "")  # anon or service_role — used as apikey header

if not SUPABASE_URL:
    logger.warning("SUPABASE_URL is not set, JWT verification will fail")

JWKS_URL = f"{SUPABASE_URL}/auth/v1/keys"

# Custom session timeout must match INACTIVITY_LIMIT_MS in auth-provider.tsx
MAX_SESSION_AGE = 3600  # 1 hour (seconds)

# Cache JWKS in memory, re-fetch if key ID not found (key rotation)
_jwks_cache: Optional[dict] = None


# ─── Models ───────────────────────────────────────────────────────────────────
@dataclass
class UserClaims:
    sub: str                    # Supabase user UUID
    email: Optional[str] = None
    role: Optional[str] = None  # "authenticated"
    aal:  Optional[str] = None


# ─── JWKS fetch ───────────────────────────────────────────────────────────────
def _get_jwks(force_refresh: bool = False) -> dict:
    """
    Fetch Supabase JWKS with the apikey header.
    Supabase requires the anon/service key even for the public JWKS endpoint.
    """
    global _jwks_cache
    if _jwks_cache and not force_refresh:
        return _jwks_cache
    try:
        res = requests.get(
            JWKS_URL,
            headers={
                "apikey": SUPABASE_KEY,          # ← required by Supabase
                "Authorization": f"Bearer {SUPABASE_KEY}",
            },
            timeout=10,
        )
        res.raise_for_status()
        _jwks_cache = res.json()
        logger.info("JWKS fetched successfully")
        return _jwks_cache
    except Exception as e:
        logger.error(f"Failed to fetch JWKS from {JWKS_URL}: {e}")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Auth service unavailable — cannot fetch JWT keys",
        )


# ─── Token verification ───────────────────────────────────────────────────────
def _verify_token(token: str) -> UserClaims:
    try:
        # Get the key ID from the token header (unverified, just reading kid)
        headers = jwt.get_unverified_header(token)
        kid = headers.get("kid")

        jwks = _get_jwks()

        # Find the matching key
        key = next(
            (k for k in jwks.get("keys", []) if k.get("kid") == kid),
            None,
        )

        # If key not found, the key may have rotated - force refresh once
        if key is None:
            logger.info(f"Key {kid} not in cache, refreshing JWKS...")
            jwks = _get_jwks(force_refresh=True)
            key = next(
                (k for k in jwks.get("keys", []) if k.get("kid") == kid),
                None,
            )

        if key is None:
            logger.warning(f"JWT key ID {kid} not found in JWKS")
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Invalid token: signing key not found",
                headers={"WWW-Authenticate": "Bearer"},
            )

        # Decode and verify, Supabase ECC tokens don't always include audience
        payload = jwt.decode(
            token,
            key,
            algorithms=["ES256"],
            options={
                "require": ["sub", "exp", "iat"],
                "verify_aud": False,    # Supabase ECC tokens omit audience claim
            },
        )

        # Custom session timeout check via iat (issued-at)
        iat = payload.get("iat")
        if iat:
            age = time.time() - iat
            if age > MAX_SESSION_AGE:
                raise HTTPException(
                    status_code=status.HTTP_401_UNAUTHORIZED,
                    detail="Session expired — please log in again",
                    headers={"WWW-Authenticate": "Bearer"},
                )

        return UserClaims(
            sub=payload["sub"],
            email=payload.get("email"),
            role=payload.get("role"),
            aal=payload.get("aal"),
        )

    except HTTPException:
        raise  # re-raise our own HTTPExceptions unchanged

    except ExpiredSignatureError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Session expired — please log in again",
            headers={"WWW-Authenticate": "Bearer"},
        )

    except JWTError as e:
        logger.warning(f"JWT verification failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or malformed token",
            headers={"WWW-Authenticate": "Bearer"},
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
    """Protected endpoint — returns 401 if token missing or invalid."""
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
    """Optional auth returns None if no/bad token, never raises."""
    token = _extract_bearer(authorization)
    if not token:
        return None
    try:
        return _verify_token(token)
    except HTTPException:
        return None