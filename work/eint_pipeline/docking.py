"""Fresh graph-derived conformers; archived pose coordinates are discarded."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
SEED=20260921

def graph(m):
    params=Chem.RemoveHsParameters();params.removeDefiningBondStereo=True
    return Chem.MolToSmiles(Chem.RemoveHs(m,params),isomericSmiles=True)

def write_sdf(path,mols):
    w=Chem.SDWriter(str(path))
    for m in mols:w.write(m)
    w.close()

def fresh_conformer(m):
    # Serialize the molecular graph to ensure none of the old coordinates survive.
    smiles=graph(m)
    new=Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert new.GetNumConformers()==0
    assert new.GetNumAtoms()==m.GetNumAtoms() and Chem.GetFormalCharge(new)==Chem.GetFormalCharge(m)
    params=AllChem.ETKDGv3();params.randomSeed=SEED;params.numThreads=1;params.enforceChirality=True
    rc=AllChem.EmbedMolecule(new,params)
    fallback=False
    if rc!=0:
        params.useRandomCoords=True;new.RemoveAllConformers();rc=AllChem.EmbedMolecule(new,params);fallback=True
    assert rc==0,'Fresh ETKDG embedding failed'
    if AllChem.MMFFHasAllMoleculeParams(new):
        optimization=AllChem.MMFFOptimizeMolecule(new,mmffVariant='MMFF94s',maxIters=500);method='MMFF94s'
    else:optimization=AllChem.UFFOptimizeMolecule(new,maxIters=500);method='UFF'
    assert np.isfinite(new.GetConformer().GetPositions()).all()
    assert graph(new)==smiles
    for key in ('_Name','mol_id','source_index','layer','uniprot','paff','smiles'):
        if m.HasProp(key):new.SetProp(key,m.GetProp(key))
    return new,dict(graph_smiles=smiles,embedding_seed=SEED,random_coordinates_fallback=fallback,
                    initial_forcefield=method,initial_minimization_returncode=optimization,
                    n_atoms=new.GetNumAtoms(),n_hydrogens=sum(a.GetAtomicNum()==1 for a in new.GetAtoms()),formal_charge=Chem.GetFormalCharge(new))
