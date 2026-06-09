#!/usr/bin/env sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
sh "$script_dir/run_tests.sh" "${1:-./tsoclip}"
