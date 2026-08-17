#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
QUADLET_DIR="${HOME}/.config/containers/systemd"
IMAGE="localhost/plaac-web:latest"
DOMAIN="${DOMAIN:-}"

if [ -f /etc/debian_version ] && command -v apt-get >/dev/null 2>&1; then
    packages=()

    if ! command -v podman >/dev/null 2>&1; then
        packages+=(podman)
    fi

    if ! command -v caddy >/dev/null 2>&1; then
        packages+=(caddy)
    fi

    if [ "${#packages[@]}" -gt 0 ]; then
        echo "Installing host dependencies: ${packages[*]}..."
        sudo apt-get update
        sudo apt-get install -y "${packages[@]}"
    fi
else
    if ! command -v podman >/dev/null 2>&1; then
        echo "error: podman is not installed and automatic installation is only supported on Debian" >&2
        exit 1
    fi

    if ! command -v caddy >/dev/null 2>&1; then
        echo "error: caddy is not installed and automatic installation is only supported on Debian" >&2
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

echo "Configuring Caddy..."

sudo tee /etc/caddy/Caddyfile >/dev/null <<'EOF'
:80 {
    reverse_proxy 127.0.0.1:4567 {
        header_up Host localhost
        header_up -X-Forwarded-Host
    }
}
EOF

if [ -n "$DOMAIN" ]; then
    sudo tee -a /etc/caddy/Caddyfile >/dev/null <<EOF

$DOMAIN {
    reverse_proxy 127.0.0.1:4567 {
        header_up Host localhost
        header_up -X-Forwarded-Host
    }
}
EOF
fi

sudo systemctl enable caddy
sudo systemctl restart caddy

sudo systemctl --no-pager --full status caddy
