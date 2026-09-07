import json,numpy as np,torch
from pathlib import Path
from rdkit import Chem
from aimnet_safe_calculator import AIMNet2Calculator
from aimnet_t3_interaction import read_pocket,energies
K=23.060547830619
root=Path("/home/xianyang/aimnet2_score_pipelines/work/correction-20260906")
model="/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt"
calc=AIMNet2Calculator(model,device="cuda",compile_model=False)
direct=torch.jit.load(model,map_location="cuda").eval()
results=[]
for p in sorted((root/"t3").glob("L*/*/pocket.pdb"))[:3]:
    xyz,z,q=read_pocket(p)
    mol=next(m for m in Chem.SDMolSupplier(str(p.parent/"poses.sdf"),removeHs=False) if m)
    lx=np.asarray(mol.GetConformer().GetPositions(),np.float32)
    lz=np.array([a.GetAtomicNum() for a in mol.GetAtoms()]);lq=Chem.GetFormalCharge(mol)
    systems=[(np.concatenate([xyz,lx]),np.r_[z,lz],q+lq),(xyz,z,q),(lx,lz,lq)]
    safe=energies(calc,systems)
    dense=[]
    for x,n,c in systems:
        size=len(n)
        distances=np.linalg.norm(x[:,None,:].astype(np.float64)-x[None,:,:].astype(np.float64),axis=-1)
        np.fill_diagonal(distances,np.inf)
        def neighbors(cutoff):
            indices=[np.flatnonzero(row<cutoff) for row in distances]
            arr=np.full((size+1,max(1,max(map(len,indices)))),size,dtype=np.int32)
            for i,ids in enumerate(indices):arr[i,:len(ids)]=ids
            return torch.as_tensor(arr,device="cuda")
        inputs={"coord":torch.as_tensor(np.vstack([x,np.zeros((1,3),np.float32)]),device="cuda"),
                "numbers":torch.as_tensor(np.r_[n,0],device="cuda").int(),
                "charge":torch.tensor([c],device="cuda",dtype=torch.float32),
                "mol_idx":torch.zeros(size+1,device="cuda",dtype=torch.int32),
                "mol_sizes":torch.tensor([size],device="cuda",dtype=torch.int32),
                "nbmat":neighbors(5.0),"nbmat_lr":neighbors(float("inf"))}
        with torch.inference_mode():
            e=direct(inputs)["energy"]
        dense.append(float(e.reshape(-1)[0].cpu()))
    dense=np.array(dense)
    delta=float(((safe[0]-safe[1]-safe[2])-(dense[0]-dense[1]-dense[2]))*K)
    r={"target":str(p.parent),"safe_eV":safe.tolist(),"direct_dense_eV":dense.tolist(),"interaction_delta_kcal_mol":delta,"pass":abs(delta)<0.1}
    print(json.dumps(r),flush=True);results.append(r)
(root/"direct_dense_validation.json").write_text(json.dumps(results,indent=2))
assert all(r["pass"] for r in results)

