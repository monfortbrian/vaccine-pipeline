import os
import time
import logging
from typing import Optional
from dataclasses import dataclass

from fastapi import Depends, HTTPException, Header, status
from jose import jwt, JWTError, ExpiredSignatureError

logger = logging.getLogger("kozi.auth")

# ─── Config ───────────────────────────────────────────────────────────────────
SUPABASE_JWT_SECRET = os.getenv("SUPABASE_JWT_SECRET", "")

if not SUPABASE_JWT_SECRET:
    logger.warning(
        "SUPABASE_JWT_SECRET is not set, all protected endpoints will return 401. "
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


# ─── Core verification ────────────────────────────────────────────────────────
def _verify_token(token: str) -> UserClaims:
    if not SUPABASE_JWT_SECRET:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Server misconfiguration: SUPABASE_JWT_SECRET not set",
        )

    try:
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
    except JWTError as e:
        logger.warning(f"JWT verification failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or malformed token",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Custom session age check via iat (issued-at timestamp)
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