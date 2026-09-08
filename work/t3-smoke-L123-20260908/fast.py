#!/usr/bin/env python3
"""CPU docking queue plus concurrent GPU scoring; unchanged scientific protocol."""
from run import *
import functools, shlex, traceback
STATE = HOME / "fast"
CPU_WORKERS = 16
GPU_WORKERS = 8

def submit(mode, lane=None, dependency=None):
    args = ["sbatch", "--parsable", "--job-name=t3-fast-" + mode,
            "--output=" + str(STATE / "logs" / (mode + "-%j.out"))]
    if mode in ("gpu", "benchmark"):
        args += ["--partition=general", "--gres=gpu:1", "--cpus-per-task=4",
                 "--mem=32G", "--exclude=babel-x5-20",
                 "--time=" + ("01:00:00" if mode == "benchmark" else "2-00:00:00")]
    else:
        args += ["--partition=preempt", "--qos=preempt_cpu_qos", "--cpus-per-task=4",
                 "--mem=8G", "--time=2-00:00:00"]
    if dependency:
        args += ["--dependency=afterok:" + str(dependency)]
    cmd = [str(BIN / "python"), "-u", str(HOME / "fast.py"), mode]
    if lane is not None:
        cmd += [str(lane)]
    args += ["--wrap", shlex.join(cmd)]
    jid = subprocess.check_output(args, text=True).strip().split(";")[0]
    with (STATE / "jobs.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        jp = STATE / "jobs.json"
        jobs = json.loads(jp.read_text()) if jp.exists() else []
        jobs.append({"mode": mode, "lane": lane, "job_id": jid})
        atomic_json(jp, jobs)
    log("SUBMITTED", mode, lane, jid)
    return jid

@functools.lru_cache(maxsize=32)
def library(path):
    return pd.read_csv(Path(path) / "library.csv.gz")

def chunk_dir(task):
    return Path(task["target_dir"]) / "chunks" / f"{task['start']:08d}_{task['stop']:08d}"

def publish(task, metadata):
    atomic_json(chunk_dir(task) / "docked.json", metadata)
    atomic_json(STATE / "ready" / f"{task['id']:06d}.json", task)

def recover_docking(task):
    from rdkit import Chem
    cd = chunk_dir(task)
    frame = library(task["target_dir"]).iloc[task["start"]:task["stop"]]
    # score.log is opened only after successful SMINA exit in the original worker.
    for attempt in sorted(cd.glob("attempt-*"), reverse=True):
        if not (attempt / "score.log").exists() or not (attempt / "poses.sdf").exists():
            continue
        try:
            mols = [m for m in Chem.SDMolSupplier(str(attempt / "poses.sdf"), removeHs=False) if m]
            ids = [int(m.GetProp("mol_id")) if m.HasProp("mol_id") else int(m.GetProp("_Name").rsplit(":", 1)[-1]) for m in mols]
            if not ids or len(set(ids)) != len(ids) or not set(ids).issubset(set(frame.mol_id)):
                continue
            failed = [{"mol_id": int(mid), "error": "No pose in recovered completed docking"}
                      for mid in frame.mol_id if mid not in ids]
            publish(task, {"poses": str(attempt / "poses.sdf"), "failed": failed,
                           "reused": True, "source_attempt": str(attempt)})
            return True
        except Exception:
            continue
    return False

def dock_task_fast(task):
    import t3_full_dock as dock
    from rdkit import Chem
    cd = chunk_dir(task)
    cd.mkdir(parents=True, exist_ok=True)
    if (cd / "complete.json").exists():
        return
    if (cd / "docked.json").exists():
        atomic_json(STATE / "ready" / f"{task['id']:06d}.json", task)
        return
    if recover_docking(task):
        return
    td = Path(task["target_dir"])
    frame = library(str(td)).iloc[task["start"]:task["stop"]].copy()
    frame["uniprot"] = td.name
    ad = cd / ("dock-" + str(time.time_ns()))
    ad.mkdir()
    failed = []
    writer = Chem.SDWriter(str(ad / "input.sdf"))
    for row in frame.itertuples(index=False):
        try:
            mol = dock.embed_ligand(pd.Series(row._asdict()))
            if mol is None:
                raise ValueError("Embedding or SMILES parsing failed")
            writer.write(mol)
        except Exception as exc:
            failed.append({"mol_id": int(row.mol_id), "error": str(exc)})
    writer.close()
    poses = ad / "poses.sdf"
    if len(failed) < len(frame):
        box = json.loads((td / "box.json").read_text())
        cmd = [str(BIN / "smina"), "-r", str(td / "pocket.pdbqt"),
               "-l", str(ad / "input.sdf"), "--exhaustiveness", "8",
               "--num_modes", "1", "--cpu", "4", "--seed", "1", "-o", str(poses)]
        for i, axis in enumerate("xyz"):
            cmd += ["--center_" + axis, str(box["center"][i]), "--size_" + axis, str(box["size"][i])]
        with (ad / "dock.log").open("w") as h:
            proc = subprocess.run(cmd, stdout=h, stderr=subprocess.STDOUT, timeout=14400)
        if proc.returncode:
            raise RuntimeError(f"SMINA exit {proc.returncode}: {ad}/dock.log")
    publish(task, {"poses": str(poses), "failed": failed, "reused": False})
    log("DOCK_READY", task["id"], task["layer"], task["uniprot"], task["start"])

def score_task_fast(task):
    import t3_full_dock as dock
    from rdkit import Chem
    td = Path(task["target_dir"])
    start, stop = task["start"], task["stop"]
    cd = chunk_dir(task)
    attempt = cd / ("attempt-" + str(time.time_ns()))
    attempt.mkdir()
    frame = library(str(td)).iloc[start:stop].copy()
    frame["uniprot"] = td.name
    meta = json.loads((cd / "docked.json").read_text())
    pose_path = Path(meta["poses"])
    failed = meta["failed"]
    poses = [m for m in Chem.SDMolSupplier(str(pose_path), removeHs=False) if m] if pose_path.exists() else []
    pose_ids = [int(m.GetProp("mol_id")) if m.HasProp("mol_id") else int(m.GetProp("_Name").rsplit(":", 1)[-1]) for m in poses]
    if len(pose_ids) != len(set(pose_ids)) or not set(pose_ids).issubset(set(frame.mol_id)):
        raise RuntimeError("Docking returned duplicate or unexpected ligand IDs")
    energies = {}
    for mid, mol in zip(pose_ids, poses):
        energy = dock.affinity_from_mol(mol)
        if energy is not None and np.isfinite(energy):
            energies[mid] = energy
    frame["smina_affinity_kcal_mol"] = frame.mol_id.map(energies)
    if poses:
        out = attempt / "aimnet2_score.csv"
        command = [str(BIN / "python"), "-u", str(CODE / "aimnet2_composite_smoke.py"),
                   "--poses", str(pose_path), "--pocket", str(td / "pocket.pdb"),
                   "--interaction-model", "aimnet2-2025", "--gas-model", str(MODEL / "aimnet2_wb97m_0.jpt"),
                   "--cpcm-model", str(MODEL / "wb97m_cpcms_v2_0.jpt"),
                   "--output", str(out), "--max-poses", str(len(poses)),
                   "--complex-max-steps", "1000", "--max-steps", "1000", "--device", "cuda"]
        env = {**os.environ, "OMP_NUM_THREADS": "4", "MKL_NUM_THREADS": "4"}
        with (attempt / "score.log").open("w") as h:
            proc = subprocess.run(command, stdout=h, stderr=subprocess.STDOUT, env=env, timeout=14400)
        if proc.returncode:
            raise RuntimeError(f"AIMNet process failed ({proc.returncode}); see {attempt}/score.log")
        scored = pd.read_csv(out)
        if scored.mol_id.duplicated().any() or set(scored.mol_id) != set(pose_ids):
            raise RuntimeError("AIMNet ligand IDs do not match docking output")
        good = np.isfinite(scored.aimnet2_composite_kcal_mol)
        f = scored.loc[good]
        if not np.allclose(f.aimnet2_composite_kcal_mol,
                           f.aimnet_interaction_kcal_mol + f.desolvation_kcal_mol + f.local_strain_kcal_mol,
                           atol=1e-6, rtol=0):
            raise RuntimeError("Composite energy arithmetic failed")
        if not np.allclose(f.aimnet_interaction_kcal_mol,
                           (f.complex_bound_ev - f.pocket_ev - f.interaction_ligand_ev) * 23.060547830619,
                           atol=1e-6, rtol=0):
            raise RuntimeError("Interaction energy arithmetic failed")
        frame = frame.merge(scored, on="mol_id", how="left", validate="one_to_one")
    for column in METHODS.values():
        if column not in frame:
            frame[column] = np.nan
    if "error" not in frame:
        frame["error"] = ""
    frame["error"] = frame.error.fillna("")
    failed_map = {x["mol_id"]: x["error"] for x in failed}
    for i, row in frame.iterrows():
        if int(row.mol_id) in failed_map:
            frame.at[i, "error"] = "embedding: " + failed_map[int(row.mol_id)]
        elif int(row.mol_id) not in pose_ids:
            frame.at[i, "error"] = "docking: no valid pose returned"
    frame["scoring_ok"] = np.isfinite(frame.aimnet2_composite_kcal_mol)
    atomic_csv(cd / "results.csv.gz", frame)
    atomic_json(cd / "complete.json", {"n_sampled": len(frame), "n_docked": len(poses),
                "n_scored": int(frame.scoring_ok.sum()), "attempt": str(attempt),
                "n_failures": int((~frame.scoring_ok).sum())})
    log("CHUNK", td.parent.name, td.name, start, stop, "scored", int(frame.scoring_ok.sum()))

def maybe_summarize(task):
    td = Path(task["target_dir"])
    marker = HOME / "target_metrics" / (task["layer"] + "_" + task["uniprot"] + ".json")
    if marker.exists():
        return
    with (td / "aggregate.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        if marker.exists():
            return
        frame = library(str(td))
        if not all((td / "chunks" / f"{s:08d}_{min(s+CHUNK,len(frame)):08d}" / "complete.json").exists()
                   for s in range(0, len(frame), CHUNK)):
            return
        record = next(r for r in json.loads((HOME / "manifest.json").read_text())["targets"]
                      if r["layer"] == task["layer"] and r["uniprot"] == task["uniprot"])
        summarize_target(record, frame)

def prepare_fast():
    import itertools
    tasks = []
    records = json.loads((HOME / "manifest.json").read_text())["targets"]
    groups = [[r for r in records if r["layer"] == layer and r["status"] == "ready"]
              for layer in ("L1", "L2", "L3")]
    for group in itertools.zip_longest(*groups):
        for r in group:
            if r is None:
                continue
            for start in range(0, r["n"], CHUNK):
                tasks.append({"id": len(tasks), "layer": r["layer"], "uniprot": r["uniprot"],
                              "target_dir": r["target_dir"], "start": start, "stop": min(start+CHUNK,r["n"])})
    atomic_json(STATE / "tasks.json", tasks)
    completed, recovered = 0, 0
    for task in tasks:
        cd = chunk_dir(task)
        if (cd / "complete.json").exists():
            completed += 1
            maybe_summarize(task)
        elif cd.exists() and recover_docking(task):
            recovered += 1
    for f in HOME.glob("worker_*.json"):
        shutil.move(str(f), str(STATE / "legacy_status" / f.name))
    old_launch = HOME / "launch.json"
    shutil.copy2(old_launch, STATE / "launch_original.json")
    atomic_json(STATE / "migration.json", {"completed_chunks_preserved": completed,
                "docked_chunks_recovered": recovered, "n_chunks": len(tasks),
                "cpu_workers": CPU_WORKERS, "gpu_workers": GPU_WORKERS,
                "excluded_gpu_node": "babel-x5-20"})
    pilot_id = submit("benchmark")
    cpu = [{"lane": i, "job_id": submit("cpu", i)} for i in range(CPU_WORKERS)]
    gpu = [{"lane": i, "job_id": submit("gpu", i, dependency=pilot_id)} for i in range(GPU_WORKERS)]
    atomic_json(HOME / "launch.json", {"pilot": pilot_id, "workers": gpu, "cpu_workers": cpu,
                "implementation": "fast.py", "state": str(STATE)})
    log("MIGRATED", completed, recovered, len(tasks))

def cpu_worker(lane):
    tasks = json.loads((STATE / "tasks.json").read_text())[lane::CPU_WORKERS]
    deadline = time.monotonic() + 39 * 3600
    failures = []
    with (STATE / f"cpu_{lane}.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for task in tasks:
            if time.monotonic() > deadline:
                jid = submit("cpu", lane, dependency=os.environ.get("SLURM_JOB_ID"))
                atomic_json(STATE / f"cpu_{lane}.json", {"status": "continuation", "next_job_id": jid})
                return
            try:
                dock_task_fast(task)
            except Exception as exc:
                # Keep runtime failures separate and visible; they are not ligand labels.
                failure = {"task": task, "error": repr(exc)}
                failures.append(failure)
                atomic_json(STATE / "errors" / f"dock_{task['id']}.json", failure)
                atomic_json(STATE / f"cpu_{lane}.json", {"status": "failed", "failure": failure})
                raise
            atomic_json(STATE / f"cpu_{lane}.json", {"status": "running", "last_task": task["id"],
                        "job_id": os.environ.get("SLURM_JOB_ID"), "updated": time.time()})
        atomic_json(STATE / f"cpu_{lane}.json", {"status": "complete"})
        atomic_json(STATE / f"cpu_{lane}.done", {"complete": True})

def score_worker(worker_id):
    deadline = time.monotonic() + 36 * 3600
    completed = 0
    while time.monotonic() < deadline:
        ready = sorted((STATE / "ready").glob("*.json"))
        if not ready:
            if all((STATE / f"cpu_{i}.done").exists() for i in range(CPU_WORKERS)):
                return 0
            time.sleep(10)
            continue
        worked = False
        for path in ready:
            try:
                task = json.loads(path.read_text())
            except FileNotFoundError:
                continue
            cd = chunk_dir(task)
            with (cd / "score.lock").open("a") as lock:
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    continue
                if not path.exists():
                    continue
                try:
                    if not (cd / "complete.json").exists():
                        score_task_fast(task)
                    maybe_summarize(task)
                    path.unlink(missing_ok=True)
                    completed += 1
                    atomic_json(STATE / f"score_{worker_id}.json", {"status": "running",
                        "last_task": task["id"], "chunks_completed": completed,
                        "job_id": os.environ.get("SLURM_JOB_ID"), "updated": time.time()})
                except Exception as exc:
                    failure = {"task": task, "error": repr(exc), "job_id": os.environ.get("SLURM_JOB_ID")}
                    atomic_json(STATE / "errors" / f"score_{task['id']}.json", failure)
                    atomic_json(STATE / f"score_{worker_id}.json", {"status": "failed", "failure": failure})
                    raise
                worked = True
                break
        if not worked:
            time.sleep(5)
    return 75

def gpu_worker(lane):
    import torch
    if not torch.cuda.is_available():
        raise RuntimeError("GPU preflight failed before starting any scoring")
    test = torch.ones(1, device="cuda")
    torch.cuda.synchronize()
    log("GPU_READY", torch.cuda.get_device_name(0), os.environ.get("CUDA_VISIBLE_DEVICES"))
    config = json.loads((STATE / "benchmark.json").read_text())
    streams = config["streams_per_gpu"]
    children = [subprocess.Popen([str(BIN / "python"), "-u", str(HOME / "fast.py"), "score",
                                  str(lane * streams + i)],
                                 env={**os.environ, "OMP_NUM_THREADS": "4", "MKL_NUM_THREADS": "4"})
                for i in range(streams)]
    statuses = []
    try:
        while any(p.poll() is None for p in children):
            bad = [p.returncode for p in children if p.poll() is not None and p.returncode not in (0, 75)]
            if bad:
                raise RuntimeError("Scoring subprocess failed; inspect fast/errors and GPU logs")
            time.sleep(5)
        statuses = [p.returncode for p in children]
    finally:
        for p in children:
            if p.poll() is None:
                p.terminate()
        for p in children:
            p.wait()
    if any(code not in (0, 75) for code in statuses):
        raise RuntimeError("Scoring subprocess failed")
    if 75 in statuses:
        jid = submit("gpu", lane, dependency=os.environ.get("SLURM_JOB_ID"))
        atomic_json(STATE / f"gpu_{lane}.json", {"status": "continuation", "next_job_id": jid})
    else:
        atomic_json(STATE / f"gpu_{lane}.json", {"status": "complete"})

def benchmark():
    from rdkit import Chem
    import torch
    selftest()
    if not torch.cuda.is_available():
        raise RuntimeError("GPU preflight failed")
    # Reuse real completed poses; compare serial and concurrent scores on the same GPU.
    tasks = json.loads((STATE / "tasks.json").read_text())
    for task in tasks:
        cd = chunk_dir(task)
        if (cd / "complete.json").exists():
            meta = json.loads((cd / "complete.json").read_text())
            poses = Path(meta["attempt"]) / "poses.sdf"
            if poses.exists():
                mols = [m for m in Chem.SDMolSupplier(str(poses), removeHs=False) if m][:4]
                if len(mols) == 4:
                    break
    else:
        raise RuntimeError("No completed four-ligand docking chunk for benchmark")
    bd = STATE / ("benchmark-" + str(time.time_ns()))
    bd.mkdir()
    for name, group in [("serial", mols), ("parallel0", mols[::2]), ("parallel1", mols[1::2])]:
        writer = Chem.SDWriter(str(bd / (name + ".sdf")))
        for mol in group:
            writer.write(mol)
        writer.close()
    def command(name, count):
        return [str(BIN / "python"), "-u", str(CODE / "aimnet2_composite_smoke.py"),
                "--poses", str(bd / (name + ".sdf")), "--pocket", str(Path(task["target_dir"]) / "pocket.pdb"),
                "--interaction-model", "aimnet2-2025", "--gas-model", str(MODEL / "aimnet2_wb97m_0.jpt"),
                "--cpcm-model", str(MODEL / "wb97m_cpcms_v2_0.jpt"),
                "--output", str(bd / (name + ".csv")), "--max-poses", str(count),
                "--complex-max-steps", "1000", "--max-steps", "1000", "--device", "cuda"]
    env = {**os.environ, "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    begin = time.monotonic()
    with (bd / "serial.log").open("w") as h:
        subprocess.run(command("serial", 4), stdout=h, stderr=subprocess.STDOUT, env=env, check=True, timeout=1800)
    serial_time = time.monotonic() - begin
    begin = time.monotonic()
    handles = [(bd / f"parallel{i}.log").open("w") for i in range(2)]
    procs = [subprocess.Popen(command(f"parallel{i}", 2), stdout=handles[i],
                             stderr=subprocess.STDOUT, env=env) for i in range(2)]
    for p in procs:
        if p.wait(timeout=1800):
            raise RuntimeError("Concurrent scoring benchmark failed")
    for h in handles:
        h.close()
    parallel_time = time.monotonic() - begin
    a = pd.read_csv(bd / "serial.csv").set_index("mol_id").sort_index()
    b = pd.concat([pd.read_csv(bd / f"parallel{i}.csv") for i in range(2)]).set_index("mol_id").sort_index()
    columns = list(METHODS.values())[1:] + ["desolvation_kcal_mol", "local_strain_kcal_mol"]
    assert a.index.equals(b.index)
    assert np.isfinite(a[columns]).all().all() and np.isfinite(b[columns]).all().all()
    agreement = bool(np.allclose(a[columns], b[columns], atol=1e-5, rtol=1e-7))
    flags_agree = all(a[col].equals(b[col]) for col in ("optimization_converged", "complex_optimization_converged"))
    streams = 2 if agreement and flags_agree and parallel_time < .95 * serial_time else 1
    result = {"status": "passed", "n_ligands": 4, "serial_seconds": serial_time,
              "parallel_seconds": parallel_time, "speedup": serial_time / parallel_time,
              "max_abs_energy_difference": float(np.abs(a[columns]-b[columns]).to_numpy().max()),
              "streams_per_gpu": streams, "concurrent_scores_agree": agreement, "convergence_flags_agree": flags_agree,
              "gpu": torch.cuda.get_device_name(0), "directory": str(bd),
              "note": "Small controlled comparison; not an end-to-end completion-time estimate."}
    atomic_json(STATE / "benchmark.json", result)
    log("BENCHMARK_PASSED", result)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prepare", "cpu", "gpu", "score", "benchmark"])
    ap.add_argument("lane", nargs="?", type=int)
    args = ap.parse_args()
    if args.mode in ("cpu", "gpu", "score") and args.lane is None:
        ap.error("lane is required")
    if args.mode == "prepare":
        prepare_fast()
    elif args.mode == "benchmark":
        benchmark()
    elif args.mode == "cpu":
        cpu_worker(args.lane)
    elif args.mode == "gpu":
        gpu_worker(args.lane)
    else:
        sys.exit(score_worker(args.lane))
