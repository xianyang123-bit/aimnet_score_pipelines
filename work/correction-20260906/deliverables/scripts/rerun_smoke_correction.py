import sys,json,shutil,gzip
from pathlib import Path
import pandas as pd
from rdkit import Chem
import aimnet_chunk_smoke as chunk
import aimnet2_composite_smoke as comp
worker=int(sys.argv[1])
src=Path("/data/user_data/xianyang/t3-aimnet-full/full-dock-L4/B4RNM9")
root=Path("/home/xianyang/aimnet2_score_pipelines/work/correction-20260906/smoke")
d=root/f"part{worker}";d.mkdir(parents=True,exist_ok=True)
shutil.copy2(src/"pocket.pdb",d/"pocket.pdb")
shutil.copy2(src/"library.csv.gz",d/"library.csv.gz")
pose_src=src/"chunks/poses_000000_000250.sdf.gz"
with gzip.open(pose_src,"rb") as h: mols=[m for m in Chem.ForwardSDMolSupplier(h,removeHs=False) if m]
part=mols[worker::4]
with gzip.open(d/"poses.sdf.gz","wt") as h:
    with Chem.SDWriter(h) as w:
        for m in part:w.write(m)
sys.argv=["chunk","--target-dir",str(d),"--poses",str(d/"poses.sdf.gz"),"--out",str(d/"aimnet_interaction.csv"),"--batch-size","4"]
chunk.main()
sys.argv=["composite","--poses",str(d/"poses.sdf.gz"),"--interaction-csv",str(d/"aimnet_interaction.csv"),"--cpcm-model","/home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt","--gas-model","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt","--output",str(d/"aimnet2_score.csv"),"--max-poses","250","--refine-steps","300","--max-steps","1000","--torsion-force-constant","1","--device","cuda"]
comp.main()
if worker==0:
    shutil.copy2(pose_src,root/"poses.sdf.gz")
    shutil.copy2(src/"pocket.pdb",root/"pocket.pdb")
    shutil.copy2(src/"library.csv.gz",root/"library.csv.gz")
    inventory=[]
    for base in [Path("/data/user_data/xianyang/t3-aimnet-full"),Path("/data/user_data/xianyang/casf-2016/ligunity")]:
        for p in base.rglob("*"):
            if p.is_file() and "aimnet" in p.name.lower() and p.suffix in [".csv",".json"]:
                inventory.append({"path":str(p),"bytes":p.stat().st_size})
    (root.parent/"source_result_inventory.json").write_text(json.dumps(inventory,indent=2))
(root/f"worker_{worker}.json").write_text(json.dumps({"n_input":len(part),"completed":True}))

