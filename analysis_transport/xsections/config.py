"""Load and fingerprint the authoritative thin-target scan configuration."""

import hashlib
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 and earlier
    import tomli as tomllib


def load(path: Path):
    raw = path.read_bytes()
    values = tomllib.loads(raw.decode("utf-8"))
    values["config_sha256"] = hashlib.sha256(raw).hexdigest()
    values["config_path"] = str(path)
    return values
