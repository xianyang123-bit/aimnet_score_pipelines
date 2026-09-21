"""Strict prepared-receptor input contract for AIMNet scoring.

PDB charge columns are often blank; absence is never evidence of neutrality.
The sidecar binds an explicitly derived integer charge to the exact PDB bytes.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem


def read_pocket(path):
    path = Path(path)
    meta_path = path.with_suffix('.prep.json')
    if not meta_path.is_file():
        raise ValueError(f'Unprepared pocket: {path}. Missing {meta_path.name}. '
                         'Rebuild from a residue-labelled source structure with '
                         'prepare_aimnet_pocket.py; blank PDB charges cannot imply zero.')
    meta = json.loads(meta_path.read_text())
    if meta.get('schema') != 'aimnet_prepared_pocket_v1' or meta.get('status') != 'validated':
        raise ValueError('Pocket preparation did not pass validation')
    if hashlib.sha256(path.read_bytes()).hexdigest() != meta.get('pdb_sha256'):
        raise ValueError('Pocket PDB changed after preparation; regenerate its preparation record')
    charge = meta.get('net_charge')
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise ValueError('Preparation record must provide an explicit integer net_charge')
    xyz, numbers = [], []
    table = Chem.GetPeriodicTable()
    for line in path.read_text().splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        if line[17:20].strip() in ('POC', 'UNK'):
            raise ValueError('Dummy/unknown protein residue in prepared pocket')
        element = line[76:78].strip()
        if not element:
            raise ValueError('Explicit element symbols are required')
        numbers.append(table.GetAtomicNumber(element.title()))
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    xyz, numbers = np.asarray(xyz, np.float32), np.asarray(numbers, np.int64)
    if not len(numbers) or not np.isfinite(xyz).all() or not (numbers == 1).any():
        raise ValueError('Prepared pocket must have finite coordinates and explicit hydrogens')
    if len(numbers) != meta.get('n_atoms') or int((numbers == 1).sum()) != meta.get('n_hydrogens'):
        raise ValueError('Pocket atom counts disagree with preparation record')
    return xyz, numbers, charge
