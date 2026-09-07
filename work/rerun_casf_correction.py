import sys,json,gc,shutil,hashlib
from pathlib import Path
import pandas as pd
import torch
import aimnet_casf_interaction as interaction
import aimnet2_composite_smoke as composite
import analyze_casf_target as analyze
root=Path("/home/xianyang/aimnet2_score_pipelines/work/correction-20260906")
source=Path("/data/user_data/xianyang/casf-2016/ligunity/multi24")
worker=int(sys.argv[1]);workers=int(sys.argv[2])
tasks=pd.read_csv(source/"selected_targets.csv")
dest=root/"casf"; dest.mkdir(exist_ok=True)
if worker==0:shutil.copy2(source/"selected_targets.csv",dest/"selected_targets.csv")
jobs=[(r.target,source/r.target,dest/r.target,f"{r.target}_top20_allposes.sdf") for r in tasks.itertuples()]
jobs.append(("3ebp",Path("/data/user_data/xianyang/casf-2016/ligunity/rerank-3ebp-top20"),root/"casf-single","3ebp_retrieved_top20_allposes.sdf"))
fail=[]
for index in range(worker,len(jobs),workers):
    target,src,d,pose=jobs[index];d.mkdir(exist_ok=True)
    print("START",index,target,flush=True)
    try:
        for name in ["pocket.pdb",pose]:
            shutil.copy2(src/name,d/name)
        retrieval=src/"ligunity_all285.csv"
        if not retrieval.exists():retrieval=Path("/data/user_data/xianyang/casf-2016/ligunity/results")/f"{target}_ligunity_all285.csv"
        shutil.copy2(retrieval,d/"ligunity_all285.csv")
        for name in ["aimnet_interaction_bestpose.csv","aimnet2_composite_rerank.csv","comparison_summary.json"]:
            if (src/name).exists():shutil.copy2(src/name,d/(name+".before_20260906"))
        sys.argv=["interaction","--pocket",str(d/"pocket.pdb"),"--poses",str(d/pose),"--model","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt","--output-all",str(d/"aimnet_interaction_allposes.csv"),"--output-best",str(d/"aimnet_interaction_bestpose.csv"),"--output-best-sdf",str(d/"aimnet_interaction_bestpose.sdf"),"--batch-size","4"]
        interaction.main()
        sys.argv=["composite","--poses",str(d/"aimnet_interaction_bestpose.sdf"),"--interaction-csv",str(d/"aimnet_interaction_bestpose.csv"),"--cpcm-model","/home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt","--gas-model","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt","--output",str(d/"aimnet2_composite_rerank.csv"),"--summary",str(d/"aimnet2_composite_rerank.summary.json"),"--max-poses","20","--refine-steps","300","--max-steps","1000","--torsion-force-constant","1","--device","cuda"]
        composite.main()
        sys.argv=["analyze","--target-dir",str(d),"--target",target]
        analyze.main()
        print("DONE",index,target,flush=True)
    except Exception as e:
        fail.append({"task":index,"target":target,"error":repr(e)})
        print("FAILED",index,target,repr(e),flush=True)
    gc.collect();torch.cuda.empty_cache()
(root/f"casf_worker_{worker}.json").write_text(json.dumps({"failures":fail},indent=2))
if fail:sys.exit(1)

