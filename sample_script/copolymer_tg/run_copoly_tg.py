#!/usr/bin/env python3
"""
Random copolymer Tg workflow — configurable two-monomer version.

Steps:
  1a. Prepare monomer 1 object  (Gasteiger charges by default)
  1b. Prepare monomer 2 object
  2.  Equilibration MD           (builds amorphous copolymer cell)
  3.  Tg cooling MD

Each step is skipped when its output already exists — safe to resume.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR   = Path(__file__).resolve().parent
RADONPY_ROOT = (SCRIPT_DIR / "../..").resolve()
SAMPLE_DIR   = RADONPY_ROOT / "sample_script"

if str(RADONPY_ROOT) not in sys.path:
    sys.path.insert(0, str(RADONPY_ROOT))


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="run_copoly_tg.py",
        description=(
            "Random copolymer Tg workflow.\n"
            "Prepares two monomers, runs equilibration MD, then Tg cooling MD."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
# MA + MMA, 1:1 ratio, 1 MPI process (quick smoke test)
  python run_copoly_tg.py MA "O=C(OC)C(*)(C(*))" MMA "O=C(OC)C(*)(C(*))C"

# Custom ratio (30 % MON1 / 70 % MON2), 4 MPI processes
  python run_copoly_tg.py MA "O=C(OC)C(*)(C(*))" MMA "O=C(OC)C(*)(C(*))C" \\
      --ratio 0.3,0.7 --mpi 4

# Production run with more atoms/chains
  python run_copoly_tg.py MA "O=C(OC)C(*)(C(*))" MMA "O=C(OC)C(*)(C(*))C" \\
      --natom 1000 --nchain 10 --mpi 8

# Alternating copolymer with RESP charges
  python run_copoly_tg.py MA "O=C(OC)C(*)(C(*))" MMA "O=C(OC)C(*)(C(*))C" \\
      --copoly-type alternating --charge resp

Notes
-----
  * SMILES must use * for polymer attachment points.
  * --ratio accepts two comma-separated weights that do NOT need to sum to 1;
    RadonPy normalises them internally (e.g. "1,3" means 25 %% / 75 %%).
  * For RESP charges, pre-run qm.py for each monomer and place results in
    <MON_ID>/analyze/ before calling this script, then pass --charge resp.
""",
    )

    # Positional — the four values everyone always needs
    p.add_argument("mon1_id",     metavar="MON1_ID",
                   help="Identifier for monomer 1 (e.g. MA)")
    p.add_argument("mon1_smiles", metavar="MON1_SMILES",
                   help="SMILES repeat unit of monomer 1 (use * for attachment points)")
    p.add_argument("mon2_id",     metavar="MON2_ID",
                   help="Identifier for monomer 2 (e.g. MMA)")
    p.add_argument("mon2_smiles", metavar="MON2_SMILES",
                   help="SMILES repeat unit of monomer 2 (use * for attachment points)")

    cop = p.add_argument_group("Copolymer settings")
    cop.add_argument("--ratio",       default="0.498,1", metavar="F1,F2",
                     help="Mole-fraction weights mon1:mon2 (default: %(default)s)")
    cop.add_argument("--copoly-type", default="random",
                     choices=["random", "alternating", "block"],
                     help="Sequence type (default: %(default)s)")
    cop.add_argument("--tacticity",   default="atactic",
                     choices=["atactic", "isotactic", "syndiotactic"],
                     help="Chain tacticity (default: %(default)s)")
    cop.add_argument("--temp",  default=300.0, type=float, metavar="K",
                     help="MD temperature in K (default: %(default)s)")
    cop.add_argument("--press", default=1.0,   type=float, metavar="atm",
                     help="MD pressure in atm (default: %(default)s)")

    ter = p.add_argument_group("Terminal group")
    ter.add_argument("--ter-id",     default="CH3", metavar="ID",
                     help="Terminal group identifier (default: %(default)s)")
    ter.add_argument("--ter-smiles", default="*C",  metavar="SMILES",
                     help="Terminal group SMILES (default: %(default)s)")

    sz = p.add_argument_group("System size")
    sz.add_argument("--natom",  default=int(os.environ.get("RadonPy_NAtom",  200)), type=int,
                    help="Target atoms per chain (default: %(default)s)")
    sz.add_argument("--nchain", default=int(os.environ.get("RadonPy_NChain",  5)),  type=int,
                    help="Number of chains (default: %(default)s)")

    par = p.add_argument_group("Parallelism")
    par.add_argument("--mpi", default=1, type=int, metavar="N",
                     help="MPI processes (default: %(default)s)")
    par.add_argument("--omp", default=0, type=int, metavar="N",
                     help="OpenMP threads per MPI rank, 0 = auto (default: %(default)s)")

    misc = p.add_argument_group("Miscellaneous")
    misc.add_argument("--charge",   default=os.environ.get("RadonPy_Charge", "gasteiger"),
                      choices=["gasteiger", "resp"],
                      help="Charge assignment method (default: %(default)s)")
    misc.add_argument("--retry-eq", default=0, type=int, metavar="N",
                      help="Max equilibration retries (default: %(default)s)")

    return p.parse_args()


# ── Environment ───────────────────────────────────────────────────────────────

def build_env(args: argparse.Namespace) -> dict[str, str]:
    lmp_conda = Path("/Users/shiqian/miniconda3/envs/radonpy_v1/bin/lmp")
    lmp_build = Path("/Users/shiqian/Documents/app/lammps-22Jul2025/build/lmp")
    lammps_exec = str(
        lmp_conda if lmp_conda.is_file() and os.access(lmp_conda, os.X_OK)
        else lmp_build
    )

    return {
        **os.environ,
        "PYTHONPATH":            f"{RADONPY_ROOT}:{os.environ.get('PYTHONPATH', '')}",
        "MPLCONFIGDIR":          str(SCRIPT_DIR / ".mplconfig"),
        "XDG_CACHE_HOME":        str(SCRIPT_DIR / ".cache"),
        "FLIB_FASTOMP":          "FALSE",
        "FLIB_CNTL_BARRIER_ERR": "FALSE",
        "LAMMPS_EXEC":           lammps_exec,
        "RadonPy_MPI":           str(args.mpi),
        "RadonPy_OMP":           str(args.omp),
        "RadonPy_GPU":           "0",
        "RadonPy_NAtom":         str(args.natom),
        "RadonPy_NChain":        str(args.nchain),
        "RadonPy_TER_ID":        args.ter_id,
        "RadonPy_SMILES_TER":    args.ter_smiles,
        "RadonPy_Charge":        args.charge,
        "RadonPy_RetryEQ":       str(args.retry_eq),
        "RadonPy_OMP_Psi4":      os.environ.get("RadonPy_OMP_Psi4",      "1"),
        "RadonPy_MEM_Psi4":      os.environ.get("RadonPy_MEM_Psi4",      "2000"),
        "RadonPy_Conf_MM_MPI":   os.environ.get("RadonPy_Conf_MM_MPI",   "0"),
        "RadonPy_Conf_MM_OMP":   os.environ.get("RadonPy_Conf_MM_OMP",   "1"),
        "RadonPy_Conf_MM_MP":    os.environ.get("RadonPy_Conf_MM_MP",    "0"),
        "RadonPy_Conf_Psi4_OMP": os.environ.get("RadonPy_Conf_Psi4_OMP", "1"),
        "RadonPy_Conf_Psi4_MP":  os.environ.get("RadonPy_Conf_Psi4_MP",  "0"),
    }


# ── Steps ─────────────────────────────────────────────────────────────────────

def run_step(script: Path, env: dict, label: str) -> None:
    print(f"\n=== {label} ===")
    result = subprocess.run(
        [sys.executable, str(script)],
        env=env,
        cwd=str(SCRIPT_DIR),
    )
    if result.returncode != 0:
        sys.exit(result.returncode)


def ensure_terminal(dbid: str, env: dict) -> None:
    """Build terminal end-cap object if not already present."""
    from radonpy.core import calc, utils

    ter_id     = env["RadonPy_TER_ID"]
    ter_smiles = env["RadonPy_SMILES_TER"]
    charge     = env["RadonPy_Charge"]
    save_dir   = SCRIPT_DIR / dbid / "analyze"
    save_dir.mkdir(parents=True, exist_ok=True)

    json_path   = save_dir / f"ter_{ter_id}.json"
    pickle_path = save_dir / f"ter_{ter_id}.pickle"

    if json_path.is_file():
        mol = utils.JSONToMol(str(json_path))
    elif pickle_path.is_file():
        mol = utils.pickle_load(str(pickle_path))
    else:
        mol = utils.mol_from_smiles(ter_smiles)

    if any(not a.HasProp("AtomicCharge") for a in mol.GetAtoms()):
        calc.assign_charges(mol, charge=charge)

    utils.pickle_dump(mol, str(pickle_path))
    utils.MolToJSON(mol, str(json_path))
    print(f"[terminal] {ter_id} ready: {pickle_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args        = parse_args()
    copoly_dbid = f"{args.mon1_id}_{args.mon2_id}"
    env         = build_env(args)

    (SCRIPT_DIR / ".mplconfig").mkdir(exist_ok=True)
    (SCRIPT_DIR / ".cache").mkdir(exist_ok=True)

    # Step 1a — monomer 1
    run_step(
        SCRIPT_DIR / "prepare_monomer.py",
        {**env,
         "RadonPy_DBID":       args.mon1_id,
         "RadonPy_Monomer_ID": args.mon1_id,
         "RadonPy_SMILES":     args.mon1_smiles},
        f"Step 1a: Prepare {args.mon1_id} monomer",
    )

    # Step 1b — monomer 2
    run_step(
        SCRIPT_DIR / "prepare_monomer.py",
        {**env,
         "RadonPy_DBID":       args.mon2_id,
         "RadonPy_Monomer_ID": args.mon2_id,
         "RadonPy_SMILES":     args.mon2_smiles},
        f"Step 1b: Prepare {args.mon2_id} monomer",
    )

    # Step 2 — equilibration (skip if already done)
    results_csv = SCRIPT_DIR / copoly_dbid / "analyze" / "results.csv"
    if results_csv.is_file():
        print(f"\n[eq] Already done, skipping (found {results_csv})")
    else:
        ensure_terminal(copoly_dbid, env)
        run_step(
            SAMPLE_DIR / "eq.py",
            {**env,
             "RadonPy_DBID":         copoly_dbid,
             "RadonPy_Monomer_ID":   f"{args.mon1_id},{args.mon2_id}",
             "RadonPy_Monomer_Dir":  f"{args.mon1_id}/analyze,{args.mon2_id}/analyze",
             "RadonPy_SMILES":       f"{args.mon1_smiles},{args.mon2_smiles}",
             "RadonPy_Copoly_Ratio": args.ratio,
             "RadonPy_Copoly_Type":  args.copoly_type,
             "RadonPy_Tacticity":    args.tacticity,
             "RadonPy_Temp":         str(args.temp),
             "RadonPy_Press":        str(args.press)},
            f"Step 2: Equilibration MD ({args.mon1_id}/{args.mon2_id} copolymer)",
        )

    # Step 3 — Tg cooling
    run_step(
        SAMPLE_DIR / "tg.py",
        {**env, "RadonPy_DBID": copoly_dbid},
        "Step 3: Tg cooling MD",
    )

    print(f"\n=== Done! ===")
    print(f"Results: {SCRIPT_DIR}/{copoly_dbid}/analyze/results.csv")


if __name__ == "__main__":
    main()
