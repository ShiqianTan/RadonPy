#!/usr/bin/env python3
"""
Fast monomer preparation without Psi4 QM.

Generates a monomer mol object using ETKDG geometry + Gasteiger charges
and saves the required files for eq.py (pickle, JSON, monomer_data CSV).

Usage (called from run_copoly_tg.zsh via env vars):
  RadonPy_DBID=MA RadonPy_Monomer_ID=MA RadonPy_SMILES='*CC(C(=O)OC)*' python prepare_monomer.py

Env vars:
  RadonPy_DBID        - Working directory name (required)
  RadonPy_Monomer_ID  - Monomer identifier used in filenames (defaults to DBID)
  RadonPy_SMILES      - Repeat-unit SMILES with * linker atoms (required)
  RadonPy_Charge      - Charge type: gasteiger (default) or zero
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RADONPY_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '../..'))
if RADONPY_ROOT not in sys.path:
    sys.path.insert(0, RADONPY_ROOT)

import matplotlib
matplotlib.use('Agg')

import pandas as pd
from radonpy.core import utils, calc, poly
from radonpy.sim import helper

dbid = os.environ.get('RadonPy_DBID')
monomer_id = os.environ.get('RadonPy_Monomer_ID', dbid)
smiles = os.environ.get('RadonPy_SMILES')
charge = os.environ.get('RadonPy_Charge', 'gasteiger')

if not dbid or not smiles:
    print('ERROR: RadonPy_DBID and RadonPy_SMILES must be set.', file=sys.stderr)
    sys.exit(1)

work_dir = os.path.join(SCRIPT_DIR, dbid)
save_dir = os.path.join(work_dir, 'analyze')
os.makedirs(save_dir, exist_ok=True)

data_csv = os.path.join(save_dir, 'monomer_%s_data.csv' % monomer_id)
if os.path.isfile(data_csv):
    print('[prepare_monomer] %s: already done, skipping.' % monomer_id)
    sys.exit(0)

print('[prepare_monomer] %s  SMILES=%s  charge=%s' % (monomer_id, smiles, charge))

mol = utils.mol_from_smiles(smiles)
if mol is None:
    print('ERROR: Cannot parse SMILES: %s' % smiles, file=sys.stderr)
    sys.exit(1)

if not poly.set_linker_flag(mol):
    print('ERROR: No linker atoms (*) found in SMILES: %s' % smiles, file=sys.stderr)
    sys.exit(1)

if mol.GetIntProp('head_idx') == mol.GetIntProp('tail_idx'):
    print('ERROR: Only one linker atom found. Use a repeat-unit SMILES such as *CC(C)*.', file=sys.stderr)
    sys.exit(1)

calc.assign_charges(mol, charge=charge)

utils.pickle_dump(mol, os.path.join(save_dir, 'monomer_%s.pickle' % monomer_id))
utils.MolToJSON(mol, os.path.join(save_dir, 'monomer_%s.json' % monomer_id))

monomer_data = {
    'monomer_ID': monomer_id,
    'smiles': smiles,
    'qm_method': 'MM+%s' % charge,
    'charge': charge,
    **helper.get_version(),
    'mol_weight': calc.molecular_weight(mol),
    'vdw_volume': calc.vdw_volume(mol),
}
pd.DataFrame([monomer_data]).set_index('monomer_ID').to_csv(data_csv)

print('[prepare_monomer] %s: done. Saved to %s' % (monomer_id, save_dir))
