#!/bin/bash
set -e

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
VENV="$SCRIPT_DIR/.build-venv"

if [ ! -d "$VENV" ]; then
    echo "Creating Python build environment..." >&2
    python3 -m venv "$VENV"
fi

. "$VENV/bin/activate"

python3 -m pip install --upgrade pip >/dev/null
python3 -m pip install -q "setuptools-scm>=8,<9"

python3 - <<'EOF'
from setuptools_scm import get_version

print(get_version(
    root="..",
    version_scheme="guess-next-dev",
    local_scheme="node-and-date",
    fallback_version="1.1.0.dev0",
))
EOF
