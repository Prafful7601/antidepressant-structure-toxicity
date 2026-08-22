#!/bin/bash
# Download the MIMIC-IV Clinical Database Demo v2.2 and build it into DuckDB,
# reusing MIT-LCP/mimic-code's own buildmimic/duckdb pipeline rather than a
# hand-rolled loader (it creates the schema, loads hosp/ and icu/, derives
# the mimic-code "concepts" views used for severity scoring in task 6, and
# runs its own row-count check). This script then adds an independent, strict
# validation pass that exits non-zero on any mismatch (build_mimic.sh's own
# validate() step only warns).
#
# Usage:
#   ./scripts/build_database.sh
#
# Re-running is safe: mimic-code's build script records completed steps
# inside the database and skips them; the demo download is skipped if the
# files already exist (wget -N -c).
#
# To build against the FULL MIMIC-IV database instead of the demo, once
# PhysioNet access is granted: download the full hosp/ and icu/ data
# yourself (this script only fetches the open demo), point MIMIC_DATA_DIR at
# that directory, and change config.py's _DEFAULT_DB_PATH (see that file).
# This script itself does not need to change.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

# Pinned to the commit this pipeline was built and tested against, so a
# future mimic-code change can't silently alter the build.
readonly MIMIC_CODE_COMMIT="8bcbd190ca75670cd5281f9ead3611ae1cefb73e"
readonly MIMIC_CODE_DIR="${PROJECT_ROOT}/external/mimic-code"
readonly DEMO_VERSION="2.2"
readonly DATA_DIR="${MIMIC_DATA_DIR:-${PROJECT_ROOT}/data/raw}"
readonly DB_FILE="${PROJECT_ROOT}/data/mimic4_demo.duckdb"

echo "== 1/4: DuckDB CLI"
if ! command -v duckdb >/dev/null 2>&1; then
    if [ -x "${HOME}/.local/bin/duckdb" ]; then
        export PATH="${HOME}/.local/bin:${PATH}"
    else
        echo "ERROR: 'duckdb' CLI not found on PATH." >&2
        echo "  Install it (e.g. from https://github.com/duckdb/duckdb/releases," >&2
        echo "  pick duckdb_cli-osx-universal.zip for macOS) and re-run." >&2
        exit 1
    fi
fi
duckdb --version

echo "== 2/4: mimic-code (pinned commit ${MIMIC_CODE_COMMIT})"
if [ ! -d "${MIMIC_CODE_DIR}" ]; then
    git clone https://github.com/MIT-LCP/mimic-code.git "${MIMIC_CODE_DIR}"
fi
git -C "${MIMIC_CODE_DIR}" fetch --depth 1 origin "${MIMIC_CODE_COMMIT}"
git -C "${MIMIC_CODE_DIR}" checkout --quiet "${MIMIC_CODE_COMMIT}"

echo "== 3/4: MIMIC-IV demo v${DEMO_VERSION} data (hosp/, icu/)"
mkdir -p "${DATA_DIR}"
if [ -z "$(find "${DATA_DIR}/hosp" -maxdepth 1 -name '*.csv.gz' -print -quit 2>/dev/null)" ]; then
    # PhysioNet serves the whole open demo project as one zip; the archive's
    # top-level dir contains hosp/ and icu/ (plus ed/, which we don't need
    # here) straight from the published project, so no per-file listing or
    # wget recursion is required.
    zip_path="${PROJECT_ROOT}/data/mimic-iv-demo-${DEMO_VERSION}.zip"
    curl -fL --retry 5 --retry-connrefused -o "${zip_path}" \
        "https://physionet.org/content/mimic-iv-demo/get-zip/${DEMO_VERSION}/"
    tmp_extract=$(mktemp -d)
    unzip -q "${zip_path}" -d "${tmp_extract}"
    extracted_root="${tmp_extract}/mimic-iv-clinical-database-demo-${DEMO_VERSION}"
    cp -R "${extracted_root}/hosp" "${DATA_DIR}/hosp"
    cp -R "${extracted_root}/icu" "${DATA_DIR}/icu"
    rm -rf "${tmp_extract}" "${zip_path}"
else
    echo "   already downloaded, skipping"
fi

echo "== 4/4: build DuckDB database (schema, load, concepts, mimic-code's own row-count check)"
(
    cd "${MIMIC_CODE_DIR}/mimic-iv/buildmimic/duckdb"
    ./build_mimic.sh "${DATA_DIR}" "${DB_FILE}"
)

echo "== independent strict row-count validation (fails loudly on any mismatch)"
# Uses whatever Python is on PATH; activate the project's venv first (see
# README) so duckdb/pandas are importable.
python3 "${SCRIPT_DIR}/validate_row_counts.py"

echo "Database ready: ${DB_FILE}"
