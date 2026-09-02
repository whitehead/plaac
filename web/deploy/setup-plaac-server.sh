#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
QUADLET_DIR="${HOME}/.config/containers/systemd"
CONTAINER_IMAGE="${1:-${IMAGE:-}}"
IMAGE="${CONTAINER_IMAGE:-localhost/plaac-web:latest}"
DOMAIN="${DOMAIN:-}"

if [ -f /etc/debian_version ] && command -v apt-get >/dev/null 2>&1; then
    echo "Installing host dependencies..."
    sudo apt-get update
    sudo apt-get install -y podman caddy python3-venv
else
    for command in podman caddy python3; do
        if ! command -v "$command" >/dev/null 2>&1; then
            echo "error: $command is not installed and automatic installation is only supported on Debian" >&2
            exit 1
        fi
    done

    if ! python3 -c 'import venv' >/dev/null 2>&1; then
        echo "error: Python venv support is not installed" >&2
        exit 1
    fi
fi

if [ -n "$CONTAINER_IMAGE" ]; then
    echo "Using container image $IMAGE..."
    podman pull "$IMAGE"
else
    PLAAC_VERSION=$("$REPO_DIR/cli/get_plaac_version.sh")
    PLAAC_GIT_SHA=$(git -C "$REPO_DIR" rev-parse HEAD)
    echo "Building $IMAGE..."
    podman build --build-arg "PLAAC_GIT_SHA=$PLAAC_GIT_SHA"  --build-arg "PLAAC_VERSION=$PLAAC_VERSION" --cgroup-manager=cgroupfs --format docker -t "$IMAGE" -f "$REPO_DIR/web/Dockerfile" "$REPO_DIR"
fi

mkdir -p "$QUADLET_DIR"

sed "s|@IMAGE@|$IMAGE|g" \
    "$SCRIPT_DIR/quadlet/plaac.container" |
    install -m 0644 /dev/stdin "$QUADLET_DIR/plaac.container"

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
        header_down Location http://localhost/ http://{http.request.host}/
    }
}
EOF

if [ -n "$DOMAIN" ]; then
    sudo tee -a /etc/caddy/Caddyfile >/dev/null <<EOF

$DOMAIN {
    reverse_proxy 127.0.0.1:4567 {
        header_up Host localhost
        header_up -X-Forwarded-Host
        header_down Location https://localhost/ https://$DOMAIN/
    }
}
EOF
fi

sudo systemctl enable caddy
sudo systemctl restart caddy

sudo systemctl --no-pager --full status caddy
