import os,sys,json,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
worker=int(sys.argv[1]);jobs=json.loads((root/"jobs.json").read_text());fail=[]
(root/"claims").mkdir(exist_ok=True)
env={**os.environ,"OMP_NUM_THREADS":"2","MKL_NUM_THREADS":"2"}
for index,j in enumerate(jobs):
 try:(root/"claims"/str(index)).mkdir()
 except FileExistsError:continue
 d=root/j["key"];d.mkdir(parents=True,exist_ok=True)
 print("START",index,j["key"],flush=True)
 try:
  if "allposes" in j:
   cmd=[sys.executable,str(root/"aimnet_casf_interaction.py"),"--pocket",j["pocket"],"--poses",j["allposes"],
        "--model","aimnet2-2025","--output-all",str(d/"aimnet_interaction_allposes.csv"),"--output-best",j["interaction"],
        "--output-best-sdf",j["poses"],"--batch-size","4"]
   with (d/"allposes.log").open("w") as h:subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT,env=env,check=True)
  cmd=[sys.executable,str(root/"aimnet2_composite_smoke.py"),"--poses",j["poses"],"--pocket",j["pocket"],
       "--interaction-csv",j["interaction"],"--interaction-model","aimnet2-2025",
       "--gas-model","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt",
       "--cpcm-model","/home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt",
       "--output",str(d/j["output"]),"--max-poses",str(j["n"]),"--complex-max-steps","1000","--max-steps","1000"]
  with (d/"run.log").open("w") as h:subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT,env=env,check=True)
  summary=json.loads((d/j["output"]).with_suffix(".summary.json").read_text())
  assert summary["n_errors"]==0,summary
  print("DONE",index,j["key"],flush=True)
 except Exception as e:
  fail.append(dict(index=index,key=j["key"],error=repr(e)));print("FAILED",index,j["key"],repr(e),flush=True)
(root/f"worker_{worker}.json").write_text(json.dumps({"failures":fail},indent=2))
sys.exit(bool(fail))
