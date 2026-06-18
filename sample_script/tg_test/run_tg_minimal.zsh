#!/usr/bin/env bash
set -euo pipefail

# Minimal local Tg workflow for RadonPy.
#
# Existing EQ result:
#   ./run_tg_minimal.zsh <DBID> [MPI] [OMP]
#
# From a monomer SMILES:
#   ./run_tg_minimal.zsh <DBID> '<SMILES>' [MPI] [OMP]
#
# Example:
#   ./run_tg_minimal.zsh P1 '*CC(C)*' 1 0

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
RADONPY_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"

DBID="${1:-}"
ARG2="${2:-}"

if [[ -z "${DBID}" ]]; then
  printf "%s\n" "Usage:" >&2
  printf "%s\n" "  ./run_tg_minimal.zsh <DBID> [MPI] [OMP]" >&2
  printf "%s\n" "  ./run_tg_minimal.zsh <DBID> '<SMILES>' [MPI] [OMP]" >&2
  printf "%s\n" "Example:" >&2
  printf "%s\n" "  ./run_tg_minimal.zsh P1 '*CC(C)*' 1 0" >&2
  exit 2
fi

SMILES=""
if [[ -n "${ARG2}" && ! "${ARG2}" =~ ^[0-9]+$ ]]; then
  SMILES="${ARG2}"
  MPI="${3:-1}"
  OMP="${4:-0}"
else
  MPI="${2:-1}"
  OMP="${3:-0}"
fi

mkdir -p "${SCRIPT_DIR}/.mplconfig" "${SCRIPT_DIR}/.cache"
export MPLCONFIGDIR="${SCRIPT_DIR}/.mplconfig"
export XDG_CACHE_HOME="${SCRIPT_DIR}/.cache"

set +u
source "/Users/shiqian/miniconda3/etc/profile.d/conda.sh"
conda activate radonpy_v1
set -u

export PYTHONPATH="${RADONPY_ROOT}:${PYTHONPATH:-}"
export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE

if [[ -x "/Users/shiqian/miniconda3/envs/radonpy_v1/bin/lmp" ]]; then
  export LAMMPS_EXEC="/Users/shiqian/miniconda3/envs/radonpy_v1/bin/lmp"
else
  export LAMMPS_EXEC="/Users/shiqian/Documents/app/lammps-22Jul2025/build/lmp"
fi

export RadonPy_DBID="${DBID}"
export RadonPy_Monomer_ID="${DBID}"
export RadonPy_TER_ID="${RadonPy_TER_ID:-CH3}"
export RadonPy_SMILES_TER="${RadonPy_SMILES_TER:-*C}"
export RadonPy_Charge="${RadonPy_Charge:-gasteiger}"
export RadonPy_OMP="${OMP}"
export RadonPy_MPI="${MPI}"
export RadonPy_GPU=0

# Small settings keep this as a runnable smoke example, not a production Tg job.
export RadonPy_NAtom="${RadonPy_NAtom:-200}"
export RadonPy_NChain="${RadonPy_NChain:-5}"
export RadonPy_RetryEQ="${RadonPy_RetryEQ:-0}"
export RadonPy_OMP_Psi4="${RadonPy_OMP_Psi4:-1}"
export RadonPy_MEM_Psi4="${RadonPy_MEM_Psi4:-2000}"
export RadonPy_Conf_MM_MPI="${RadonPy_Conf_MM_MPI:-0}"
export RadonPy_Conf_MM_OMP="${RadonPy_Conf_MM_OMP:-1}"
export RadonPy_Conf_MM_MP="${RadonPy_Conf_MM_MP:-0}"
export RadonPy_Conf_Psi4_OMP="${RadonPy_Conf_Psi4_OMP:-1}"
export RadonPy_Conf_Psi4_MP="${RadonPy_Conf_Psi4_MP:-0}"

cd "${SCRIPT_DIR}"

validate_smiles_arg() {
  [[ -z "${SMILES}" ]] && return

  python -c '
import os
import sys
from radonpy.core import poly, utils

smi = os.environ["RadonPy_SMILES"]
mol = utils.mol_from_smiles(smi)
if not poly.set_linker_flag(mol) or mol.GetIntProp("head_idx") == mol.GetIntProp("tail_idx"):
    sys.stderr.write(
        "Invalid RadonPy polymer SMILES: %s\n"
        "Use a repeat-unit SMILES with two linker atoms, for example: *CC(C)*\n" % smi
    )
    sys.exit(2)
'
}

check_existing_smiles() {
  [[ -z "${SMILES}" ]] && return
  [[ ! -f "${DBID}/analyze/qm_data.csv" ]] && return

  python -c '
import csv
import os
import sys

dbid = os.environ["RadonPy_DBID"]
requested = os.environ["RadonPy_SMILES"]
path = os.path.join(dbid, "analyze", "qm_data.csv")
with open(path, newline="") as fh:
    row = next(csv.DictReader(fh))
existing = row.get("smiles_1") or row.get("smiles_list")
if existing and existing != requested:
    sys.stderr.write(
        "%s already contains QM data for SMILES %r, not %r.\n"
        "Use a new DBID, or remove the old DBID directory before regenerating QM data.\n"
        % (dbid, existing, requested)
    )
    sys.exit(2)
'
}

check_existing_monomer_linkers() {
  [[ ! -d "${DBID}/analyze" ]] && return

  python -c '
import os
import sys
from radonpy.core import poly, utils

dbid = os.environ["RadonPy_DBID"]
save_dir = os.path.join(dbid, "analyze")
paths = [
    os.path.join(save_dir, "monomer_%s.json" % dbid),
    os.path.join(save_dir, "monomer_%s.pickle" % dbid),
]
mol = None
if os.path.isfile(paths[0]):
    mol = utils.JSONToMol(paths[0])
elif os.path.isfile(paths[1]):
    mol = utils.pickle_load(paths[1])

if mol is None:
    sys.exit(0)

if not poly.set_linker_flag(mol) or mol.GetIntProp("head_idx") == mol.GetIntProp("tail_idx"):
    sys.stderr.write(
        "%s/analyze contains a monomer without two RadonPy linker atoms.\n"
        "For propylene/polypropylene use a repeat-unit SMILES such as *CC(C)* with a fresh DBID.\n"
        % dbid
    )
    sys.exit(2)
'
}

ensure_terminal_obj() {
  local ter_json="${DBID}/analyze/ter_${RadonPy_TER_ID}.json"
  local ter_pickle="${DBID}/analyze/ter_${RadonPy_TER_ID}.pickle"

  python -c '
import os
from radonpy.core import calc, utils

dbid = os.environ["RadonPy_DBID"]
ter_id = os.environ.get("RadonPy_TER_ID", "CH3")
ter_smiles = os.environ.get("RadonPy_SMILES_TER", "*C")
charge = os.environ.get("RadonPy_Charge", "gasteiger")
save_dir = os.path.join(dbid, "analyze")
json_path = os.path.join(save_dir, f"ter_{ter_id}.json")
pickle_path = os.path.join(save_dir, f"ter_{ter_id}.pickle")
os.makedirs(save_dir, exist_ok=True)

if os.path.isfile(json_path):
    mol = utils.JSONToMol(json_path)
elif os.path.isfile(pickle_path):
    mol = utils.pickle_load(pickle_path)
else:
    mol = utils.mol_from_smiles(ter_smiles)

if any(not atom.HasProp("AtomicCharge") for atom in mol.GetAtoms()):
    calc.assign_charges(mol, charge=charge)

utils.pickle_dump(mol, pickle_path)
utils.MolToJSON(mol, json_path)
'
}

if [[ ! -f "${DBID}/analyze/results.csv" ]]; then
  if [[ -n "${SMILES}" ]]; then
    export RadonPy_SMILES="${SMILES}"
  fi
  validate_smiles_arg
  check_existing_smiles
  check_existing_monomer_linkers

  if [[ ! -f "${DBID}/analyze/monomer_${DBID}.json" && ! -f "${DBID}/analyze/monomer_${DBID}.pickle" ]]; then
    if [[ -z "${SMILES}" ]]; then
      printf "%s\n" "Missing ${SCRIPT_DIR}/${DBID}/analyze/results.csv and monomer object files" >&2
      printf "%s\n" "For a fresh run, provide a monomer SMILES, for example:" >&2
      printf "%s\n" "  ./run_tg_minimal.zsh ${DBID} '*CC(C)*' ${MPI} ${OMP}" >&2
      exit 2
    fi

    python ../qm.py
  fi

  ensure_terminal_obj
  export RadonPy_Monomer_Dir="${DBID}/analyze"
  python ../eq.py
fi

python ./tg.py
