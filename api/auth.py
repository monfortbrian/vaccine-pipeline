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


# CONFIG
SUPABASE_URL = os.getenv("SUPABASE_URL")
if not SUPABASE_URL:
    raise RuntimeError("SUPABASE_URL is not set")

JWKS_URL = f"{SUPABASE_URL}/auth/v1/keys"

# YOUR CUSTOM SESSION TIMEOUT (seconds)
MAX_SESSION_AGE = 3600  # 1 hour

# Cache JWKS to avoid calling Supabase every request
_jwks_cache = None


# MODELS
@dataclass
class UserClaims:
    sub: str
    email: Optional[str] = None
    role: Optional[str] = None
    aal: Optional[str] = None



# HELPERS
def _get_jwks():
    global _jwks_cache
    if _jwks_cache is None:
        try:
            res = requests.get(JWKS_URL, timeout=5)
            res.raise_for_status()
            _jwks_cache = res.json()
        except Exception as e:
            logger.error(f"Failed to fetch JWKS: {e}")
            raise HTTPException(
                status_code=500,
                detail="Auth service unavailable",
            )
    return _jwks_cache


def _extract_bearer(authorization: Optional[str]) -> Optional[str]:
    if not authorization:
        return None
    parts = authorization.split(" ", 1)
    if len(parts) != 2 or parts[0].lower() != "bearer":
        return None
    return parts[1].strip()



# CORE VERIFY
def _verify_token(token: str) -> UserClaims:
    try:
        jwks = _get_jwks()

        headers = jwt.get_unverified_header(token)

        key = next(
            (k for k in jwks["keys"] if k["kid"] == headers["kid"]),
            None
        )

        if not key:
            raise HTTPException(status_code=401, detail="Invalid token key")

        payload = jwt.decode(
            token,
            key,
            algorithms=["ES256"],
            audience="authenticated",
            options={"require": ["sub", "exp", "iat"]},
        )

        # CUSTOM SESSION CONTROL
        iat = payload.get("iat")
        if iat:
            age = time.time() - iat
            if age > MAX_SESSION_AGE:
                raise HTTPException(
                    status_code=401,
                    detail="Session expired (backend policy)",
                )

        return UserClaims(
            sub=payload["sub"],
            email=payload.get("email"),
            role=payload.get("role"),
            aal=payload.get("aal"),
        )

    except ExpiredSignatureError:
        raise HTTPException(
            status_code=401,
            detail="Session expired",
            headers={"WWW-Authenticate": "Bearer"},
        )

    except JWTError as e:
        logger.warning(f"JWT verification failed: {e}")
        raise HTTPException(
            status_code=401,
            detail="Invalid or malformed token",
            headers={"WWW-Authenticate": "Bearer"},
        )


# FASTAPI DEPENDENCIES

async def require_user(
    authorization: Optional[str] = Header(None),
) -> UserClaims:
    token = _extract_bearer(authorization)
    if not token:
        raise HTTPException(
            status_code=401,
            detail="Missing Authorization header",
        )
    return _verify_token(token)


async def optional_user(
    authorization: Optional[str] = Header(None),
) -> Optional[UserClaims]:
    token = _extract_bearer(authorization)
    if not token:
        return None
    try:
        return _verify_token(token)
    except HTTPException:
        return None