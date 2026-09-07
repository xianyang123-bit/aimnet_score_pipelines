import json,re,time,hashlib,sys
from pathlib import Path
from collections import Counter
import numpy as np,pandas as pd,torch
from rdkit import Chem
from scipy.stats import spearmanr,pearsonr
from scipy.optimize import linear_sum_assignment
from aimnet_safe_calculator import AIMNet2Calculator
root=Path(__file__).resolve().parent
data=root/"data"
def xyz(path):
 lines=path.read_text().splitlines();n=int(lines[0]);a=[l.split() for l in lines[2:] if l.strip()]
 assert len(a)==n
 z=np.array([Chem.GetPeriodicTable().GetAtomicNumber(v[0]) for v in a],np.int64)
 x=np.array([[float(c) for c in v[1:4]] for v in a],np.float32)
 assert np.isfinite(x).all()
 return x,z
def atoms(system):
 x,z=system
 return Counter((int(zz),*np.round(xx,4).tolist()) for xx,zz in zip(x,z))
jobs=[]
for group in ["cdk2","p38","tyk2"]:
 f=pd.read_csv(data/f"{group}.csv",dtype={"molecule_id":str})
 for row in f.to_dict("records"):
  key=row["molecule_id"]
  jobs.append(dict(dataset="KIN66",group=group,id=key,paths=[str(data/f"{group}_{kind}_{key}.xyz") for kind in ["complex","protein","ligands"]],charges=[0,0,0],reference=row))
f=pd.read_csv(root/"PLA15.csv",dtype={"molecule_id":str})
for row in f.to_dict("records"):
 key=row["molecule_id"];pdb=(data/f"{key}.pdb").read_text()
 q={k:int(v) for k,v in re.findall(r"REMARK (charge(?:_a|_b)?)\s+(-?\d+)",pdb)}
 assert q["charge"]==q["charge_a"]+q["charge_b"]
 jobs.append(dict(dataset="PLA15",group="PLA15",id=key,paths=[str(data/(key+s+".xyz")) for s in ["","_a","_b"]],charges=[q["charge"],q["charge_a"],q["charge_b"]],reference=row))
expanded=[];validation=[]
for j in jobs:
 systems=[xyz(Path(p)) for p in j["paths"]]
 x,z=systems[0];xf=np.concatenate([systems[1][0],systems[2][0]]);zf=np.concatenate([systems[1][1],systems[2][1]])
 assert Counter(z)==Counter(zf)
 dist=np.linalg.norm(x[:,None,:].astype(float)-xf[None,:,:].astype(float),axis=-1)
 dist[z[:,None]!=zf[None,:]]=1e6
 ia,ib=linear_sum_assignment(dist)
 deviation=float(dist[ia,ib].max())
 assert deviation<0.002,(j["id"],deviation)
 matched=np.empty_like(xf);matched[ib]=x[ia];n=len(systems[1][1])
 validation.append(dict(dataset=j["dataset"],group=j["group"],id=j["id"],max_fragment_coordinate_difference_A=deviation))
 expanded.append(dict(j,systems=systems,geometry="released_xyz"))
 expanded.append(dict(j,systems=[systems[0],(matched[:n],zf[:n]),(matched[n:],zf[n:])],geometry="matched_complex_coordinates"))
jobs=expanded
pd.DataFrame(validation).to_csv(root/"geometry_validation.csv",index=False)
print("Validated 81 systems and matched fragment atoms; testing both geometry conventions",flush=True)
rows=pd.read_csv(root/"predictions.csv").to_dict("records") if (root/"predictions.csv").exists() else []
info=json.loads((root/"model_metadata.json").read_text()) if (root/"model_metadata.json").exists() else {}
specs=[("current_wb97m","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt",None),
       ("aimnet2_2025_default","aimnet2-2025",None),
       ("aimnet2_2025_full_coulomb","aimnet2-2025","simple")]
specs += [(f"aimnet2_2025_member{i}",f"aimnet2-b973c-2025-d3_{i}",None) for i in range(1,4)]
for label,model,coulomb in specs:
 if sum(row["model"]==label for row in rows)==len(jobs):continue
 calc=AIMNet2Calculator(model,device="cuda",compile_model=False)
 if coulomb:calc.set_lrcoulomb_method(coulomb)
 info[label]=dict(model=model,metadata=calc.metadata,coulomb_method=calc._coulomb_method,coulomb_cutoff=calc._coulomb_cutoff)
 print(label,info[label],flush=True)
 for j in jobs:
  energies=[]
  for (x,z),q in zip(j["systems"],j["charges"]):
   with torch.inference_mode():
    out=calc({"coord":x.copy(),"numbers":z.copy(),"charge":np.array([q],np.float32)})
   energies.append(float(out["energy"].reshape(-1)[0].cpu()))
  e=(energies[0]-energies[1]-energies[2])*23.060547830619
  assert np.isfinite(e)
  rows.append(dict(model=label,geometry=j["geometry"],dataset=j["dataset"],group=j["group"],molecule_id=j["id"],n_atoms=len(j["systems"][0][1]),complex_charge=j["charges"][0],protein_charge=j["charges"][1],ligand_charge=j["charges"][2],complex_ev=energies[0],protein_ev=energies[1],ligand_ev=energies[2],interaction_kcal_mol=e,**{k:v for k,v in j["reference"].items() if k!="molecule_id"}))
  print(label,j["group"],j["id"],round(e,4),flush=True)
  pd.DataFrame(rows).to_csv(root/"predictions.csv",index=False)
 del calc
 torch.cuda.empty_cache()
(root/"model_metadata.json").write_text(json.dumps(info,indent=2,default=str))
frame=pd.DataFrame(rows)
frame=frame[frame.model!="aimnet2_2025_ensemble"]
members=["aimnet2_2025_default"]+[f"aimnet2_2025_member{i}" for i in range(1,4)]
ensemble=[]
for key,g in frame[frame.model.isin(members)].groupby(["dataset","group","molecule_id","geometry"]):
 assert len(g)==4
 row=g.iloc[0].to_dict();row["model"]="aimnet2_2025_ensemble"
 for col in ["complex_ev","protein_ev","ligand_ev","interaction_kcal_mol"]:row[col]=float(g[col].mean())
 ensemble.append(row)
frame=pd.concat([frame,pd.DataFrame(ensemble)],ignore_index=True)
frame.to_csv(root/"predictions.csv",index=False)
metrics=[]
for (model,dataset,geometry),f in frame.groupby(["model","dataset","geometry"]):
 for ref in ["int_b973c","int_wb97m","int_wb97m_dzvpd","int_wb97m_dzvpd_estimated","AIMNet2_2025_interaction","PLA15_reference_interaction"]:
  if ref not in f:continue
  v=f[np.isfinite(f[ref])].copy()
  if not len(v):continue
  delta=v.interaction_kcal_mol-v[ref]
  metrics.append(dict(model=model,dataset=dataset,geometry=geometry,reference=ref,n=len(v),mae=float(abs(delta).mean()),rmse=float(np.sqrt((delta**2).mean())),bias=float(delta.mean()),max_abs_error=float(abs(delta).max()),pearson=float(pearsonr(v.interaction_kcal_mol,v[ref]).statistic),spearman=float(spearmanr(v.interaction_kcal_mol,v[ref]).statistic)))
pd.DataFrame(metrics).to_csv(root/"metrics.csv",index=False)
print(pd.DataFrame(metrics).to_string(index=False))
