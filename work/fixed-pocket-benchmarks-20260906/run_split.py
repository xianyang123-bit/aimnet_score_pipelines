import os,sys,json,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parent
jobs=json.loads((root/"split_jobs.json").read_text())
worker=int(sys.argv[1]);workers=int(sys.argv[2]);fail=[]
for job in jobs[worker::workers]:
    d=root/job["key"];d.mkdir(parents=True,exist_ok=True)
    output=d/job["output"]
    if output.with_suffix(".summary.json").exists():
        print("ALREADY COMPLETE",job["key"],flush=True)
        continue
    cmd=[sys.executable,str(root/"aimnet2_composite_smoke.py"),"--poses",job["poses"],"--pocket",job["pocket"],"--interaction-csv",job["interaction"],
         "--cpcm-model","/home/xianyang/t3-aimnet-models/wb97m_cpcms_v2_0.jpt","--gas-model","/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt",
         "--output",str(output),"--max-poses",str(job["n"]),"--complex-max-steps","1000","--complex-fmax","0.002","--max-steps","1000","--fmax","0.002"]
    print("START",job["key"],flush=True)
    with (d/"run.log").open("w") as log:
        result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,"OMP_NUM_THREADS":"2","MKL_NUM_THREADS":"2"})
    if result.returncode:fail.append({"key":job["key"],"returncode":result.returncode})
    print("FINISHED",job["key"],result.returncode,flush=True)
(root/f"split_worker_{worker}.json").write_text(json.dumps({"failures":fail},indent=2))
sys.exit(bool(fail))
