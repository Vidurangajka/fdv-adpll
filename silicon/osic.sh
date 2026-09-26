#!/usr/bin/env bash
# Run the IIC-OSIC-TOOLS container against silicon/ with sky130A selected.
# POSIX twin of osic.ps1 -- see there for why PDK is forced.
#
#   ./osic.sh bash layout/verify.sh     batch command in /work
#   ./osic.sh --shell                   interactive bash
#   ./osic.sh --gui | --stop            desktop at http://localhost:8080
#   ./osic.sh --check                   tool versions + PDK sanity
set -eu
IMAGE="${OSIC_IMAGE:-hpretl/iic-osic-tools:latest}"
WORK="$(cd "$(dirname "$0")" && pwd)"
GUI_NAME=fdv-osic
COMMON=(-e PDK=sky130A -v "$WORK:/work" -w /work)

docker info >/dev/null 2>&1 || { echo "Docker engine not reachable" >&2; exit 1; }

case "${1:-}" in
    --gui)
        docker rm -f "$GUI_NAME" >/dev/null 2>&1 || true
        docker run -d --name "$GUI_NAME" "${COMMON[@]}" -p 8080:80 "$IMAGE" >/dev/null
        echo "Desktop: http://localhost:8080  (default VNC password: abc123)" ;;
    --stop)
        docker rm -f "$GUI_NAME" >/dev/null ;;
    --shell)
        exec docker run --rm -it "${COMMON[@]}" "$IMAGE" --skip bash ;;
    --check)
        exec docker run --rm "${COMMON[@]}" "$IMAGE" --skip bash osic_check.sh ;;
    "")
        sed -n 2,8p "$0"; exit 1 ;;
    *)
        exec docker run --rm "${COMMON[@]}" "$IMAGE" --skip "$@" ;;
esac
