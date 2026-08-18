#!/bin/bash
set -e

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
VENV="$SCRIPT_DIR/.build-venv"

if [ ! -d "$VENV" ]; then
    echo "Creating Python build environment..." >&2
    python3 -m venv "$VENV"
fi

. "$VENV/bin/activate"

python3 -m pip install --upgrade pip >/dev/null
python3 -m pip install -q "setuptools-scm>=8,<9"
python -m setuptools_scm  --root "$REPO_ROOT" --config "$SCRIPT_DIR/python/pyproject.toml"
