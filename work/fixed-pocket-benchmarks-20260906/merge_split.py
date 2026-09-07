from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from rdkit import Chem
r=Path(__file__).resolve().parent;d=r/"smoke/part3"
parts=[pd.read_csv(d/"completed_before_split.csv")]
parts[0]["execution_partition"]="original_worker9"
molpaths=[d/"completed_before_split.sdf"]
summaries=[]
for i in range(4):
 p=d/f"split{i}"
 summaries.append(json.loads((p/"aimnet2_score.summary.json").read_text()))
 assert summaries[-1]["n_errors"]==0
 f=pd.read_csv(p/"aimnet2_score.csv");f["execution_partition"]=f"split{i}";parts.append(f)
 molpaths.append(p/"aimnet2_score.minimized.sdf")
f=pd.concat(parts,ignore_index=True).sort_values("source_index")
j=next(j for j in json.loads((r/"jobs.json").read_text()) if j["key"]=="smoke/part3")
expected=pd.read_csv(j["interaction"])
assert len(f)==len(expected)==62 and set(f.mol_id)==set(expected.mol_id) and not f.mol_id.duplicated().any()
assert np.isfinite(f.aimnet2_composite_kcal_mol).all()
for c,rank in [("aimnet_interaction_kcal_mol","aimnet_interaction_rank"),("aimnet2_composite_kcal_mol","aimnet2_score_rank")]:
 f[rank]=rankdata(f[c])
mols={}
for p in molpaths:
 for m in Chem.SDMolSupplier(str(p),removeHs=False):
  assert m is not None
  key=m.GetIntProp("mol_id");assert key not in mols;mols[key]=m
assert set(mols)==set(f.mol_id)
f.to_csv(d/"aimnet2_score.csv",index=False)
with Chem.SDWriter(str(d/"aimnet2_score.minimized.sdf")) as w:
 for key in f.mol_id:w.write(mols[key])
s=summaries[0]
s.update(n_requested=62,n_scored=62,n_errors=0,
 n_optimization_converged=int(f.optimization_converged.sum()),
 n_complex_optimization_converged=int(f.complex_optimization_converged.sum()),
 elapsed_seconds_this_run=None,execution="16 preserved completed ligands plus 46 ligands split across four GPUs",
 output=str(d/"aimnet2_score.csv"),minimized_poses=str(d/"aimnet2_score.minimized.sdf"))
for key in ["n_actives","n_decoys","ef1_percent","ef1_top_n","ef1_actives_in_top_n","auroc","composite_vs_smina_spearman","composite_vs_interaction_spearman"]:
 s.pop(key,None)
(d/"aimnet2_score.summary.json").write_text(json.dumps(s,indent=2))
(r/"worker_9.json").write_text(json.dumps({"failures":[],"completed_via_split_job":"10340995","original_job_stopped_after_16_complete_ligands":True},indent=2))
print("Merged 62 unique smoke ligands and 62 minimized poses")
