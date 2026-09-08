import os,json,subprocess,time
from pathlib import Path
home=Path("/home/xianyang/aimnet2_score_pipelines/work")
data=Path("/data/user_data/xianyang")
full=data/"t3-full-screen-L123-20260908"
smoke=data/"t3-smoke-L123-20260908"
new=home/"t3-smoke-30-L123-20260908"
canonical=data/"t3-aimnet-full/vs-benchmark/dataset"
assert all((canonical/n).is_file() for n in ["actives.csv.gz","decoys.csv.gz","molecules.csv.gz"])
assert json.loads((new/"preparation_summary.json").read_text())["ready_targets"]==818
names=subprocess.check_output(["squeue","-u","xianyang","-h","-o","%j"],text=True).splitlines()
assert not any(n.startswith(("t3-full-","t3-smoke-","t3-fast-")) for n in names), "Cancelled-run jobs are still active"
candidates={}
def add(p,category,root):
 if p.is_symlink() or not p.is_file():return
 resolved=p.resolve()
 assert resolved.is_relative_to(root.resolve())
 s=p.stat()
 candidates[str(resolved)]={"path":str(resolved),"category":category,"bytes":s.st_size,"mtime_ns":s.st_mtime_ns,"root":str(root.resolve())}
for p in full.glob("L*/*/library.csv.gz"):add(p,"cancelled_full_library_preparation",full)
for root in [full,smoke,home/"t3-full-screen-L123-20260908",home/"t3-smoke-L123-20260908"]:
 if not root.exists():continue
 for base,dirs,files in os.walk(root,followlinks=False):
  dirs[:]=[d for d in dirs if not (Path(base)/d).is_symlink()]
  for name in files:
   p=Path(base)/name
   if "__pycache__" in p.parts and p.suffix==".pyc":add(p,"python_cache",root)
   elif name.endswith(".tmp"):add(p,"temporary_file",root)
   elif name=="input.sdf":add(p,"cancelled_docking_input",root)
totals={}
for r in candidates.values():
 t=totals.setdefault(r["category"],{"files":0,"bytes":0});t["files"]+=1;t["bytes"]+=r["bytes"]
plan={"candidates":list(candidates.values()),"totals":totals,
"preserved":"Original datasets, completed scores and poses, scripts, manifests, and the active 30-ligand run."}
(new/"cleanup_plan.json").write_text(json.dumps(plan,indent=2))
deleted=[]
for r in candidates.values():
 p=Path(r["path"])
 assert not p.is_symlink() and p.resolve().is_relative_to(Path(r["root"]))
 s=p.stat()
 if s.st_size!=r["bytes"] or s.st_mtime_ns!=r["mtime_ns"]:
  raise RuntimeError("File changed during cleanup; stopped: "+str(p))
 p.unlink()
 deleted.append(r)
report={"status":"complete","files_removed":len(deleted),"bytes_removed":sum(r["bytes"] for r in deleted),
        "totals":totals,"preserved":plan["preserved"],"completed":time.strftime("%Y-%m-%d %H:%M:%S")}
(new/"cleanup_report.json").write_text(json.dumps(report,indent=2))
print(json.dumps(report),flush=True)
