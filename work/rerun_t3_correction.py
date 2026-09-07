import sys,json,gc
from pathlib import Path
import pandas as pd
import torch
import aimnet_t3_interaction as interaction
import aimnet2_composite_smoke as composite
root=Path("/home/xianyang/aimnet2_score_pipelines/work/correction-20260906/t3")
worker=int(sys.argv[1]);workers=int(sys.argv[2])
tasks=pd.read_csv(root/"tasks.tsv",sep="\t")
fail=[]
for index in range(worker,len(tasks),workers):
    row=tasks.iloc[index];d=Path(row.target_dir)
    print("START",index,row.layer,row.uniprot,flush=True)
    try:
        sys.argv=["interaction","--pocket",str(d/"pocket.pdb"),"--poses",str(d/"poses.sdf"),"--model","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt","--output",str(d/"aimnet_interaction.csv"),"--batch-size","4"]
        interaction.main()
        sys.argv=["composite","--poses",str(d/"poses.sdf"),"--interaction-csv",str(d/"aimnet_interaction.csv"),"--cpcm-model","/home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt","--gas-model","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt","--output",str(d/"aimnet2_score.csv"),"--summary",str(d/"aimnet2_score.summary.json"),"--max-poses","10","--refine-steps","300","--max-steps","1000","--torsion-force-constant","1","--device","cuda"]
        composite.main()
        print("DONE",index,flush=True)
    except Exception as e:
        fail.append({"task":index,"error":repr(e)})
        print("FAILED",index,repr(e),flush=True)
    gc.collect();torch.cuda.empty_cache()
(root/f"worker_{worker}.json").write_text(json.dumps({"failures":fail},indent=2))
if fail:sys.exit(1)

