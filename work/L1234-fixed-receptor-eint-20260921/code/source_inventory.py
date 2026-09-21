"""Read nearby nonprotein chemistry from the matched source structure."""
import numpy as np
from openmm import app,unit
from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
from scipy.spatial import cKDTree
from prepare_aimnet_pocket import AA,legacy_atoms,residue_key

def inventory(path, ref):
    source = app.PDBxFile(str(path))
    xyz = np.asarray(source.positions.value_in_unit(unit.angstrom))
    tree = cKDTree([r[2] for r in legacy_atoms(ref)])
    tables=[]
    with path.open() as f: PdbxReader(f).read(tables)
    comp=tables[0].getObj('chem_comp')
    names={comp.getValue('id',i):comp.getValue('name',i) for i in range(comp.getRowCount())}
    rows=[]
    for r in source.topology.residues():
        if r.name in AA or r.name in ('HOH','WAT'): continue
        atoms=list(r.atoms())
        distance=float(tree.query(xyz[[a.index for a in atoms]])[0].min())
        if distance<12:
            rows.append(dict(residue=residue_key(r), name=names.get(r.name), distance_to_reference_A=distance,
                             elements=sorted(set(a.element.symbol for a in atoms))))
    return rows
