"""Load capture area presets from config."""

from __future__ import annotations

import json
from pathlib import Path

# Package root for bundled presets.json
_PACKAGE_DIR = Path(__file__).resolve().parent

# User overrides path (optional)
_USER_PRESETS_PATH = Path.home() / ".config" / "capture_visualiser" / "presets.json"


def _load_presets() -> list[dict]:
    """Load presets from bundled config or user override."""
    for path in (_USER_PRESETS_PATH, _PACKAGE_DIR / "presets.json"):
        if path.exists():
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
                return data.get("presets", [])
    return []


def get_presets() -> list[dict]:
    """Return list of preset dicts with id, name, description, width_mm, height_mm."""
    return _load_presets()
