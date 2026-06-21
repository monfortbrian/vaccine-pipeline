"""
Structured console logger

Color-coded, level-aware, production-ready.
Import and use instead of logging.getLogger() directly.

Usage:
    from src.utils.logger import get_logger
    logger = get_logger("tope_deep.agents.Agent 3")
    logger.info("TCell Predictor started")
    logger.success("20 CTL epitopes predicted")   # custom level
    logger.warning("IEDB unavailable, MHCflurry fallback active")
    logger.error("NetMHCpan call failed: timeout")
"""

import logging
import sys
from typing import Optional

# ANSI color codes
_RESET  = "\033[0m"
_BOLD   = "\033[1m"
_DIM    = "\033[2m"

_COLORS = {
    "DEBUG":    "\033[36m",      # cyan
    "INFO":     "\033[37m",      # white
    "SUCCESS":  "\033[32m",      # green
    "WARNING":  "\033[33m",      # yellow
    "ERROR":    "\033[31m",      # red
    "CRITICAL": "\033[35m",      # magenta
}

_ICONS = {
    "DEBUG":    "·",
    "INFO":     "›",
    "SUCCESS":  "✓",
    "WARNING":  "⚠",
    "ERROR":    "✗",
    "CRITICAL": "!",
}

# Custom SUCCESS level (between INFO and WARNING)
SUCCESS_LEVEL = 25
logging.addLevelName(SUCCESS_LEVEL, "SUCCESS")


class _ColorFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        level    = record.levelname
        color    = _COLORS.get(level, _COLORS["INFO"])
        icon     = _ICONS.get(level, "›")
        name     = record.name.replace("tope_deep.", "")

        # Timestamp
        ts = self.formatTime(record, "%H:%M:%S")

        # Level badge
        badge = f"{color}{_BOLD}{icon} {level:<8}{_RESET}"

        # Logger name (dimmed)
        source = f"{_DIM}{name}{_RESET}"

        # Message
        msg = record.getMessage()
        if level in ("ERROR", "CRITICAL"):
            msg = f"{color}{_BOLD}{msg}{_RESET}"
        elif level == "SUCCESS":
            msg = f"{color}{msg}{_RESET}"
        elif level == "WARNING":
            msg = f"{_COLORS['WARNING']}{msg}{_RESET}"

        # Exception
        exc = ""
        if record.exc_info:
            exc = f"\n{_DIM}{self.formatException(record.exc_info)}{_RESET}"

        return f"  {_DIM}{ts}{_RESET}  {badge}  {source}  {msg}{exc}"


class _TOPELogger(logging.Logger):
    """Extended logger with success() method."""

    def success(self, msg: str, *args, **kwargs):
        if self.isEnabledFor(SUCCESS_LEVEL):
            self._log(SUCCESS_LEVEL, msg, args, **kwargs)


def get_logger(name: str, level: int = logging.INFO) -> _TOPELogger:
    """
    Return a color-coded TOPE_DEEP logger.

    Convention:
        tope_deep.api
        tope_deep.orchestrator
        tope_deep.agents.Agent 1   (Data Curator)
        tope_deep.agents.Agent 2   (Antigen Screener)
        tope_deep.agents.Agent 3   (TCell Predictor)
        ... etc
        tope_deep.validation
        tope_deep.storage
    """
    logging.setLoggerClass(_TOPELogger)

    logger = logging.getLogger(name)
    if logger.handlers:
        return logger  # already configured

    logger.setLevel(level)

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    handler.setFormatter(_ColorFormatter())

    logger.addHandler(handler)
    logger.propagate = False

    return logger  # type: ignore


def configure_root(level: int = logging.INFO):
    """
    Call once at startup (in main.py) to configure the root logger.
    Suppresses noisy third-party loggers.
    """
    # Root handler with color formatter
    root = logging.getLogger()
    root.setLevel(level)
    if not root.handlers:
        h = logging.StreamHandler(sys.stdout)
        h.setFormatter(_ColorFormatter())
        root.addHandler(h)

    # Suppress noisy libraries
    for noisy in ["httpx", "httpcore", "uvicorn.access", "hpack", "urllib3"]:
        logging.getLogger(noisy).setLevel(logging.WARNING)

    # Keep uvicorn.error visible but not uvicorn.access
    logging.getLogger("uvicorn").setLevel(logging.INFO)
    logging.getLogger("uvicorn.access").setLevel(logging.WARNING)