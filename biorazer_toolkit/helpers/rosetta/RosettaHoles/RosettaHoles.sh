#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
XML_PATH="${ROSETTA_HOLES_XML:-${SCRIPT_DIR}/RosettaHoles.xml}"

usage() {
    cat <<EOF
Usage: $(basename "$0") <input.pdb> [output_dir]

Required environment variables:
  ROSETTA_SCRIPTS_BIN   Path to the rosetta_scripts executable
  DALPHABALL_BIN        Path to the DAlphaBall executable required by HolesFilter

Optional environment variables:
  ROSETTA_HOLES_XML     Path to the RosettaScripts XML (default: ${XML_PATH})

Example:
  ROSETTA_SCRIPTS_BIN=/path/to/rosetta_scripts.linuxgccrelease \
  DALPHABALL_BIN=/path/to/DAlphaBall.gcc \
  $(basename "$0") model.pdb results
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage >&2
    exit 1
fi

if [[ -z "${ROSETTA_SCRIPTS_BIN:-}" ]]; then
    echo "ROSETTA_SCRIPTS_BIN is not set." >&2
    exit 1
fi

if [[ -z "${DALPHABALL_BIN:-}" ]]; then
    echo "DALPHABALL_BIN is not set." >&2
    exit 1
fi

INPUT_PDB="$(cd -- "$(dirname -- "$1")" && pwd)/$(basename -- "$1")"
OUTPUT_DIR="${2:-$(pwd)}"
SCOREFILE="${OUTPUT_DIR}/scores.sc"

if [[ ! -f "${INPUT_PDB}" ]]; then
    echo "Input PDB does not exist: ${INPUT_PDB}" >&2
    exit 1
fi

if [[ ! -f "${XML_PATH}" ]]; then
    echo "RosettaScripts XML does not exist: ${XML_PATH}" >&2
    exit 1
fi

if [[ ! -x "${ROSETTA_SCRIPTS_BIN}" ]]; then
    echo "rosetta_scripts executable is not executable: ${ROSETTA_SCRIPTS_BIN}" >&2
    exit 1
fi

if [[ ! -x "${DALPHABALL_BIN}" ]]; then
    echo "DAlphaBall executable is not executable: ${DALPHABALL_BIN}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

"${ROSETTA_SCRIPTS_BIN}" \
    -s "${INPUT_PDB}" \
    -parser:protocol "${XML_PATH}" \
    -out:path:all "${OUTPUT_DIR}" \
    -out:file:scorefile "${SCOREFILE}" \
    -holes:dalphaball "${DALPHABALL_BIN}" \
    -nstruct 1

echo "Wrote Holes score to ${SCOREFILE}"