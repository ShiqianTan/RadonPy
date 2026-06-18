#!/usr/bin/env bash
# Random copolymer Tg workflow: methyl acrylate (MA) + ethyl acrylate (MMA), 1:1
#
# Usage:
#   ./run_copoly_tg.zsh [MPI] [OMP]
#
# Examples:
#   ./run_copoly_tg.zsh            # 1 MPI process
#   ./run_copoly_tg.zsh 4 0        # 4 MPI, no OMP
#
# Steps:
#   1. Prepare MA monomer object  (Gasteiger charges, no Psi4)
#   2. Prepare MMA monomer object  (Gasteiger charges, no Psi4)
#   3. Equilibration MD           (builds 1:1 random copolymer amorphous cell)
#   4. Tg cooling MD
#
# Each step is skipped if its output files already exist (safe to resume).
#
# To use RESP charges from full QM instead of Gasteiger, run qm.py for each
# monomer manually (see tg_test/run_tg_minimal.zsh for the pattern) and
# place the results in MA/analyze/ and MMA/analyze/ before running this script.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
RADONPY_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
SAMPLE_DIR="${RADONPY_ROOT}/sample_script"

MPI="${1:-1}"
OMP="${2:-0}"

# ── Copolymer definition ────────────────────────────────
MA_ID="MA"
# MMA_ID="MMA"
MMA_ID="MMA"
# MA_SMILES="*C(C(=O)OC)*"      # methyl acrylate repeat unit
MA_SMILES="O=C(OC)C(*)(C(*))"      # methyl acrylate repeat unit
# MMA_SMILES="*CC(C(=O)OCC)*"     # ethyl acrylate repeat unit
MMA_SMILES="O=C(OC)C(*)(C(*))C"     # polymethyl acrylate (for RESP charge reference)
COPOLY_DBID="MA_MMA"
COPOLY_RATIO="0.498,1"
TER_ID="CH3"
TER_SMILES="*C"

# Small settings for a quick smoke test; increase for production:
#   RadonPy_NAtom=1000  RadonPy_NChain=10
export RadonPy_NAtom="${RadonPy_NAtom:-200}"
export RadonPy_NChain="${RadonPy_NChain:-5}"
# ────────────────────────────────────────────────────────

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

export RadonPy_OMP="${OMP}"
export RadonPy_MPI="${MPI}"
export RadonPy_GPU=0
export RadonPy_TER_ID="${TER_ID}"
export RadonPy_SMILES_TER="${TER_SMILES}"
export RadonPy_Charge="${RadonPy_Charge:-gasteiger}"
export RadonPy_RetryEQ="${RadonPy_RetryEQ:-0}"
export RadonPy_OMP_Psi4="${RadonPy_OMP_Psi4:-1}"
export RadonPy_MEM_Psi4="${RadonPy_MEM_Psi4:-2000}"
export RadonPy_Conf_MM_MPI="${RadonPy_Conf_MM_MPI:-0}"
export RadonPy_Conf_MM_OMP="${RadonPy_Conf_MM_OMP:-1}"
export RadonPy_Conf_MM_MP="${RadonPy_Conf_MM_MP:-0}"
export RadonPy_Conf_Psi4_OMP="${RadonPy_Conf_Psi4_OMP:-1}"
export RadonPy_Conf_Psi4_MP="${RadonPy_Conf_Psi4_MP:-0}"

cd "${SCRIPT_DIR}"

# ── Helper: generate terminal end-cap object ─────────────
ensure_terminal() {
  local dbid="$1"
  python - <<PYEOF
import os, sys
sys.path.insert(0, '${RADONPY_ROOT}')
from radonpy.core import calc, utils
dbid       = '${dbid}'
ter_id     = os.environ.get('RadonPy_TER_ID', 'CH3')
ter_smiles = os.environ.get('RadonPy_SMILES_TER', '*C')
charge     = os.environ.get('RadonPy_Charge', 'gasteiger')
save_dir   = os.path.join(dbid, 'analyze')
os.makedirs(save_dir, exist_ok=True)
json_path   = os.path.join(save_dir, 'ter_%s.json'   % ter_id)
pickle_path = os.path.join(save_dir, 'ter_%s.pickle' % ter_id)
if os.path.isfile(json_path):
    mol = utils.JSONToMol(json_path)
elif os.path.isfile(pickle_path):
    mol = utils.pickle_load(pickle_path)
else:
    mol = utils.mol_from_smiles(ter_smiles)
if any(not a.HasProp('AtomicCharge') for a in mol.GetAtoms()):
    calc.assign_charges(mol, charge=charge)
utils.pickle_dump(mol, pickle_path)
utils.MolToJSON(mol, json_path)
print('[terminal] %s ready: %s' % (ter_id, pickle_path))
PYEOF
}

# ── Step 1a: Prepare methyl acrylate (MA) monomer ────────
printf "\n=== Step 1a: Prepare MA monomer ===\n"
export RadonPy_DBID="${MA_ID}"
export RadonPy_Monomer_ID="${MA_ID}"
export RadonPy_SMILES="${MA_SMILES}"
python prepare_monomer.py

# ── Step 1b: Prepare ethyl acrylate (MMA) monomer ─────────
printf "\n=== Step 1b: Prepare MMA monomer ===\n"
export RadonPy_DBID="${MMA_ID}"
export RadonPy_Monomer_ID="${MMA_ID}"
export RadonPy_SMILES="${MMA_SMILES}"
python prepare_monomer.py

# ── Step 2: Equilibration MD for 1:1 MA/MMA copolymer ─────
printf "\n=== Step 2: Equilibration MD (random MA/MMA copolymer) ===\n"
if [[ ! -f "${COPOLY_DBID}/analyze/results.csv" ]]; then
  ensure_terminal "${COPOLY_DBID}"
  export RadonPy_DBID="${COPOLY_DBID}"
  export RadonPy_Monomer_ID="${MA_ID},${MMA_ID}"
  export RadonPy_Monomer_Dir="${MA_ID}/analyze,${MMA_ID}/analyze"
  export RadonPy_SMILES="${MA_SMILES},${MMA_SMILES}"
  export RadonPy_Copoly_Ratio="${COPOLY_RATIO}"
  export RadonPy_Copoly_Type="random"
  export RadonPy_Tacticity="atactic"
  export RadonPy_Temp=300.0
  export RadonPy_Press=1.0
  python "${SAMPLE_DIR}/eq.py"
else
  printf "[eq] Already done, skipping (found %s/analyze/results.csv)\n" "${COPOLY_DBID}"
fi

# ── Step 3: Tg cooling MD ─────────────────────────────────
printf "\n=== Step 3: Tg cooling MD ===\n"
export RadonPy_DBID="${COPOLY_DBID}"
python "${SAMPLE_DIR}/tg.py"

printf "\n=== Done! ===\n"
printf "Results: %s/%s/analyze/results.csv\n" "${SCRIPT_DIR}" "${COPOLY_DBID}"
