#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WEB_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
QUADLET_DIR="${HOME}/.config/containers/systemd"
IMAGE="localhost/plaac-web:latest"

if ! command -v podman >/dev/null 2>&1; then
    echo "error: podman is not installed" >&2
    exit 1
fi

echo "Building $IMAGE..."
podman build --format docker -t "$IMAGE" "$WEB_DIR"

mkdir -p "$QUADLET_DIR"

install -m 0644 \
    "$SCRIPT_DIR/quadlet/plaac.container" \
    "$QUADLET_DIR/plaac.container"

loginctl enable-linger "$USER"

systemctl --user daemon-reload
systemctl --user restart plaac.service

systemctl --user --no-pager --full status plaac.service
