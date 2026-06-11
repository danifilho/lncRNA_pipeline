#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEF_DIR="${ROOT_DIR}/containers"
IMAGE_DIR="${DEF_DIR}/images"

if [ -n "${CONTAINER_ENGINE:-}" ]; then
  ENGINE="${CONTAINER_ENGINE}"
elif command -v apptainer >/dev/null 2>&1; then
  ENGINE="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  ENGINE="singularity"
else
  echo "Neither apptainer nor singularity is available in PATH." >&2
  exit 127
fi
mkdir -p "${IMAGE_DIR}"

BUILD_ARGS=()
if [ "${1:-}" = "--remote" ]; then
  BUILD_ARGS+=(--remote)
  shift
fi

if [ -n "${SINGULARITY_BUILD_OPTS:-}" ]; then
  # shellcheck disable=SC2206
  BUILD_ARGS+=(${SINGULARITY_BUILD_OPTS})
fi

cd "${ROOT_DIR}"

if [ "$#" -gt 0 ]; then
  DEFS=()
  for name in "$@"; do
    DEFS+=("${DEF_DIR}/${name%.def}.def")
  done
else
  DEFS=()
  while IFS= read -r def; do
    DEFS+=("${def}")
  done < <(find "${DEF_DIR}" -maxdepth 1 -name "*.def" | sort)
fi

for def in "${DEFS[@]}"; do
  if [ ! -f "${def}" ]; then
    echo "Definition file not found: ${def}" >&2
    exit 1
  fi
  image="${IMAGE_DIR}/$(basename "${def}" .def).sif"
  echo "Building ${image}"
  "${ENGINE}" build "${BUILD_ARGS[@]}" --force "${image}" "${def}"
done
