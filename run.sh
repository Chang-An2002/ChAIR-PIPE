#!/usr/bin/env bash
set -eo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec bash "$SCRIPT_DIR/submit_cpu10x_for_scChiatac.mm10.nobl.sh" "$@"
