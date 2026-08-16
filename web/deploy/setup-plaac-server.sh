#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
QUADLET_DIR="${HOME}/.config/containers/systemd"
IMAGE="localhost/plaac-web:latest"

if ! command -v podman >/dev/null 2>&1; then
    if [ -f /etc/debian_version ] && command -v apt-get >/dev/null 2>&1; then
        echo "Podman is not installed; installing it..."
        sudo apt-get update
        sudo apt-get install -y podman
    else
        echo "error: podman is not installed and automatic installation is only supported on Debian" >&2
        exit 1
    fi
fi

echo "Building $IMAGE..."
podman build --cgroup-manager=cgroupfs --format docker -t "$IMAGE" -f "$REPO_DIR/web/Dockerfile" "$REPO_DIR"

mkdir -p "$QUADLET_DIR"

install -m 0644 \
    "$SCRIPT_DIR/quadlet/plaac.container" \
    "$QUADLET_DIR/plaac.container"

loginctl enable-linger "$USER"

systemctl --user daemon-reload
systemctl --user restart plaac.service

systemctl --user --no-pager --full status plaac.service
