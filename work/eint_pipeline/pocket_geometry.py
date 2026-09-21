"""Reject impossible overlaps for noncovalent protein-ligand scoring."""
import numpy as np
from scipy.spatial.distance import cdist


def check_contacts(ligand_xyz, ligand_z, pocket_xyz, pocket_z, threshold=1.6):
    ligand_xyz = np.asarray(ligand_xyz, dtype=float).reshape(-1, 3)
    pocket_xyz = np.asarray(pocket_xyz, dtype=float).reshape(-1, 3)
    ligand_z, pocket_z = np.asarray(ligand_z), np.asarray(pocket_z)
    if not np.isfinite(ligand_xyz).all() or not np.isfinite(pocket_xyz).all():
        raise ValueError('Nonfinite coordinates')
    distance = float(cdist(ligand_xyz[ligand_z != 1], pocket_xyz[pocket_z != 1]).min())
    if distance < threshold:
        raise ValueError(f'Severe protein-ligand heavy-atom overlap: {distance:.3f} A < {threshold} A. '
                         'Redock into the prepared receptor; covalent complexes require a separate protocol.')
    return distance
