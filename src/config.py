"""Load the YAML analysis configuration."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "config.yaml"


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Read the configuration file and resolve relative paths against the repository root."""
    path = Path(path) if path else DEFAULT_CONFIG
    with open(path, encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    for section in ("data", "output"):
        for key, value in config.get(section, {}).items():
            if key.endswith("_dir"):
                resolved = Path(value)
                config[section][key] = resolved if resolved.is_absolute() else ROOT / resolved
    return config
