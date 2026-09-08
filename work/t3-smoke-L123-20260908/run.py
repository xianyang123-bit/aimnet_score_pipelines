#!/usr/bin/env python3
"""Sample T3 L1-L3 docking and reconstructed AIMNet2 screening."""
import argparse, fcntl, gzip, hashlib, json, math, os, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

HOME = Path(__file__).resolve().parent
ROOT = Path("/data/user_data/xianyang/t3-smoke-L123-20260908")
DATA = Path("/data/user_data/xianyang/t3-aimnet-full/vs-benchmark/dataset")
BIN = Path("/home/xianyang/miniconda3/envs/t3-aimnet/bin")
MODEL = Path("/home/xianyang/t3-aimnet-models")
CODE = HOME / "protocol"
NWORKERS = 8
CHUNK = 25
METHODS = {"smina": "smina_affinity_kcal_mol",
           "pre_min_interaction": "pre_minimization_interaction_kcal_mol",
           "interaction": "aimnet_interaction_kcal_mol",
           "composite": "aimnet2_composite_kcal_mol"}
sys.path.insert(0, str(CODE))

def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + "." + str(os.getpid()) + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(tmp, path)

def atomic_csv(path, frame):
    path = Path(path)
    tmp = path.with_name(path.name + "." + str(os.getpid()) + ".tmp")
    frame.to_csv(tmp, index=False, compression="gzip" if path.suffix == ".gz" else None)
    os.replace(tmp, path)

def log(*args):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), *args, flush=True)

def metric(energy, labels):
    energy = np.asarray(energy, dtype=float)
    labels = np.asarray(labels, dtype=int)
    n, a = len(labels), int(labels.sum())
    if not n or not 0 < a < n:
        return {"auroc": None, "ef1": None, "top_n": None, "top_active_expected": None}
    k = max(1, math.ceil(n * .01))
    cutoff = np.partition(energy, k - 1)[k - 1]
    better, tied = energy < cutoff, energy == cutoff
    hits = float(labels[better].sum() + (k - better.sum()) * labels[tied].mean())
    ranks = rankdata(-energy, method="average")
    auc = float((ranks[labels == 1].sum() - a * (a + 1) / 2) / (a * (n - a)))
    return {"auroc": auc, "ef1": hits / k / (a / n),
            "top_n": k, "top_active_expected": hits}

def sbatch(mode, *extra, dependency=None):
    args = ["sbatch", "--parsable", "--job-name=t3-smoke-" + mode,
            "--output=" + str(HOME / "logs" / (mode + "-%j.out"))]
    if mode in ("pilot", "worker"):
        args += ["--partition=general", "--gres=gpu:1", "--cpus-per-task=4",
                 "--mem=32G", "--time=" + ("01:00:00" if mode == "pilot" else "2-00:00:00")]
    else:
        args += ["--partition=preempt", "--qos=preempt_cpu_qos", "--cpus-per-task=4",
                 "--mem=32G", "--time=02:00:00"]
    if dependency:
        args += ["--dependency=afterok:" + dependency]
    import shlex
    args += ["--wrap", shlex.join([str(BIN / "python"), "-u", str(HOME / "run.py"), mode, *map(str, extra)])]
    result = subprocess.check_output(args, text=True).strip().split(";")[0]
    log("SUBMITTED", mode, extra, result)
    return result

def prepare():
    ROOT.mkdir(parents=True, exist_ok=True)
    if (HOME / "manifest.json").exists():
        raise RuntimeError("Manifest already exists; use worker to resume.")
    original = Path("/home/xianyang/aimnet2_score_pipelines/work/t3-full-screen-L123-20260908/manifest.json")
    source_manifest = json.loads(original.read_text())
    records = []
    seed_base = 20260908
    for src in source_manifest["targets"]:
        source = Path(src["target_dir"])
        full = pd.read_csv(source / "library.csv.gz").sort_values("mol_id").reset_index(drop=True)
        if full.mol_id.duplicated().any() or not set(full.label) == {0, 1}:
            raise ValueError("Sampling requires unique IDs and both classes: " + str(source))
        n = min(250, len(full))
        active = full[full.label == 1]
        decoy = full[full.label == 0]
        na = max(1, min(len(active), n - 1, int(math.floor(n * len(active) / len(full) + .5))))
        nd = min(len(decoy), n - na)
        na = n - nd
        if not 0 < na <= len(active) or not 0 < nd <= len(decoy):
            raise ValueError("Cannot construct a two-class sample: " + str(source))
        seed_text = f"{seed_base}:{src['layer']}:{src['uniprot']}"
        seed = int.from_bytes(hashlib.sha256(seed_text.encode()).digest()[:8], "little")
        rng = np.random.default_rng(seed)
        selected = pd.concat([active.iloc[rng.choice(len(active), na, replace=False)],
                              decoy.iloc[rng.choice(len(decoy), nd, replace=False)]], ignore_index=True)
        selected = selected.iloc[rng.permutation(len(selected))].reset_index(drop=True)
        selected = selected.rename(columns={"source_index": "original_source_index"})
        selected.insert(0, "source_index", np.arange(len(selected)))
        assert len(selected) == n and not selected.mol_id.duplicated().any()
        assert int(selected.label.sum()) == na
        # Independent repeatability and source-label checks.
        again = np.random.default_rng(seed)
        expected_ids = set(active.iloc[again.choice(len(active), na, replace=False)].mol_id)
        expected_ids.update(decoy.iloc[again.choice(len(decoy), nd, replace=False)].mol_id)
        assert expected_ids == set(selected.mol_id)
        assert np.array_equal(selected.label.to_numpy(), full.set_index("mol_id").loc[selected.mol_id, "label"].to_numpy())
        td = ROOT / src["layer"] / src["uniprot"]
        td.mkdir(parents=True, exist_ok=True)
        atomic_csv(td / "library.csv.gz", selected)
        record = dict(layer=src["layer"], uniprot=src["uniprot"], n=n, n_actives=na, n_decoys=nd,
                      source_n=len(full), source_n_actives=len(active), source_n_decoys=len(decoy),
                      target_dir=str(td), status=src["status"], reason=src["reason"],
                      pocket_source=src["pocket_source"], n_chunks=math.ceil(n / CHUNK),
                      sampling_seed=str(seed), seed_base=seed_base)
        if record["status"] == "ready":
            for name in ("pocket.pdb", "pocket.pdbqt", "box.json"):
                shutil.copy2(source / name, td / name)
        atomic_json(td / "target.json", record)
        records.append(record)
        if len(records) % 100 == 0:
            log("SAMPLED", len(records), "/", len(source_manifest["targets"]))
    loads = [0] * NWORKERS
    for r in sorted((r for r in records if r["status"] == "ready"), key=lambda r: -r["n"]):
        lane = min(range(NWORKERS), key=lambda i: loads[i])
        r["lane"] = lane
        loads[lane] += r["n"]
    summary = pd.DataFrame(records).groupby("layer").agg(targets=("uniprot", "size"),
        sampled_ligands=("n", "sum"), actives=("n_actives", "sum"), decoys=("n_decoys", "sum"))
    manifest = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S"), data=str(DATA), root=str(ROOT),
                    scope="250-ligand stratified smoke test per T3 L1-L3 target",
                    source_manifest=str(original), sampling="Without replacement; proportional class allocation rounded to nearest integer, at least one of each class.",
                    seed_base=seed_base, seed_derivation="First 8 SHA256 bytes of seed:layer:uniprot, little endian; NumPy default_rng.",
                    chunk_size=CHUNK, workers=NWORKERS, lane_ligands=loads, targets=records,
                    layer_counts=json.loads(summary.to_json(orient="index")),
                    protocol="fixed_pocket_aimnet2025_minimized_v2",
                    failures="All sampled IDs retained. Missing energies rank tied last; scored-only metrics also reported.",
                    ties="Average-rank AUROC; expected active hits across the EF1 cutoff tie.",
                    finite_nonconverged="Retained and flagged, matching the smoke protocol.")
    atomic_json(HOME / "manifest.json", manifest)
    atomic_json(ROOT / "manifest.json", manifest)
    atomic_csv(HOME / "target_inventory.csv", pd.DataFrame(records))
    atomic_json(HOME / "progress.json", {"status": "prepared", "n_targets": len(records),
                "n_ready": sum(r["status"] == "ready" for r in records),
                "n_unavailable": sum(r["status"] != "ready" for r in records),
                "n_to_score": sum(r["n"] for r in records if r["status"] == "ready")})
    refresh_summary()
    log("PREPARATION_COMPLETE", summary.to_dict(orient="index"))
    launch()

def process_chunk(td, library, start, stop):
    import t3_full_dock as dock
    from rdkit import Chem
    cd = td / "chunks" / f"{start:08d}_{stop:08d}"
    cd.mkdir(parents=True, exist_ok=True)
    if (cd / "complete.json").exists():
        return
    attempt = cd / ("attempt-" + str(time.time_ns()))
    attempt.mkdir()
    frame = library.iloc[start:stop].copy()
    frame["uniprot"] = td.name
    writer = Chem.SDWriter(str(attempt / "input.sdf"))
    failed = []
    for row in frame.itertuples(index=False):
        try:
            mol = dock.embed_ligand(pd.Series(row._asdict()))
            if mol is None:
                raise ValueError("Embedding or SMILES parsing failed")
            writer.write(mol)
        except Exception as exc:
            failed.append({"mol_id": int(row.mol_id), "stage": "embedding", "error": str(exc)})
    writer.close()
    embedded = len(frame) - len(failed)
    pose_path = attempt / "poses.sdf"
    if embedded:
        box = json.loads((td / "box.json").read_text())
        command = [str(BIN / "smina"), "-r", str(td / "pocket.pdbqt"),
                   "-l", str(attempt / "input.sdf"), "--exhaustiveness", "8",
                   "--num_modes", "1", "--cpu", "4", "--seed", "1", "-o", str(pose_path)]
        for i, axis in enumerate("xyz"):
            command += ["--center_" + axis, str(box["center"][i]), "--size_" + axis, str(box["size"][i])]
        with (attempt / "dock.log").open("w") as h:
            proc = subprocess.run(command, stdout=h, stderr=subprocess.STDOUT, timeout=14400)
        if proc.returncode:
            raise RuntimeError(f"SMINA failed ({proc.returncode}); see {attempt}/dock.log")
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
        if not np.isfinite(scored.aimnet2_composite_kcal_mol).any():
            raise RuntimeError(f"All AIMNet scores failed; inspect {attempt}/score.log before resuming")
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

def summarize_target(record, library):
    td = Path(record["target_dir"])
    pieces = []
    for start in range(0, len(library), CHUNK):
        cd = td / "chunks" / f"{start:08d}_{min(start+CHUNK,len(library)):08d}"
        if not (cd / "complete.json").exists():
            raise RuntimeError("Cannot summarize an unfinished target")
        pieces.append(pd.read_csv(cd / "results.csv.gz"))
    frame = pd.concat(pieces, ignore_index=True)
    if len(frame) != len(library) or frame.mol_id.duplicated().any() or set(frame.mol_id) != set(library.mol_id):
        raise RuntimeError("Sample accounting failed")
    expected = library.set_index("mol_id").label
    if not np.array_equal(frame.label.to_numpy(), expected.loc[frame.mol_id].to_numpy()):
        raise RuntimeError("Label alignment failed")
    atomic_csv(td / "all_ligand_scores.csv.gz", frame)
    result = {k: record[k] for k in ("layer", "uniprot", "n", "n_actives", "n_decoys", "source_n", "source_n_actives", "source_n_decoys", "sampling_seed")}
    result.update(status="complete", methods={})
    for method, col in METHODS.items():
        energy = frame[col].to_numpy(float)
        labels = frame.label.to_numpy(int)
        valid = np.isfinite(energy)
        full = metric(np.where(valid, energy, np.inf), labels) if valid.any() else metric([], [])
        scored_only = metric(energy[valid], labels[valid])
        result["methods"][method] = {"n_scored": int(valid.sum()), "n_failed": int((~valid).sum()),
            "n_active_scored": int(labels[valid].sum()), "coverage": float(valid.mean()),
            "sample_failure_last": full, "scored_only": scored_only}
        ranked = frame[["mol_id", "source_index", "label", "paff", col, "error"]].copy()
        ranked["ranking_energy"] = np.where(valid, energy, np.inf)
        ranked["rank_average"] = rankdata(ranked.ranking_energy)
        ranked = ranked.sort_values(["ranking_energy", "mol_id"])
        atomic_csv(td / (method + "_ranking.csv.gz"), ranked)
    result["n_complex_converged"] = int(frame.get("complex_optimization_converged", pd.Series(dtype=bool)).fillna(False).sum())
    result["n_cpcm_converged"] = int(frame.get("optimization_converged", pd.Series(dtype=bool)).fillna(False).sum())
    atomic_json(td / "metrics.json", result)
    atomic_json(HOME / "target_metrics" / (record["layer"] + "_" + record["uniprot"] + ".json"), result)
    refresh_summary()
    return result

def refresh_summary():
    manifest = json.loads((HOME / "manifest.json").read_text())
    with (HOME / "summary.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        rows = []
        for r in manifest["targets"]:
            path = HOME / "target_metrics" / (r["layer"] + "_" + r["uniprot"] + ".json")
            row = {k: r[k] for k in ("layer", "uniprot", "n", "n_actives", "n_decoys", "source_n", "source_n_actives", "source_n_decoys", "sampling_seed")}
            row.update(status=r["status"] if r["status"] != "ready" else "pending", reason=r["reason"])
            failure_path = HOME / "target_failures" / (r["layer"] + "_" + r["uniprot"] + ".json")
            if failure_path.exists():
                row.update(status="failed", reason=json.loads(failure_path.read_text())["error"])
            if path.exists():
                row["reason"] = ""
                res = json.loads(path.read_text())
                row["status"] = "complete"
                for name, v in res["methods"].items():
                    row[name + "_n_scored"] = v["n_scored"]
                    row[name + "_coverage"] = v["coverage"]
                    for policy in ("sample_failure_last", "scored_only"):
                        for key, value in v[policy].items():
                            row[name + "_" + policy + "_" + key] = value
            rows.append(row)
        table = pd.DataFrame(rows)
        atomic_csv(HOME / "per_target_metrics.csv", table)
        layers = []
        for layer, group in table.groupby("layer"):
            completed = group[group.status == "complete"]
            lr = dict(layer=layer, n_targets=len(group), n_complete=len(completed),
                      n_unavailable=int((group.status == "unavailable").sum()),
                      n_pending=int((group.status == "pending").sum()),
                      n_failed=int((group.status == "failed").sum()),
                      n_sampled=int(group.n.sum()), n_sampled_completed=int(completed.n.sum()))
            for column in table.columns:
                if column.endswith(("_auroc", "_ef1")):
                    values = pd.to_numeric(completed[column], errors="coerce").dropna()
                    lr[column + "_macro_mean"] = float(values.mean()) if len(values) else None
                    lr[column + "_n_evaluable"] = len(values)
            layers.append(lr)
        atomic_csv(HOME / "per_layer_metrics.csv", pd.DataFrame(layers))
        atomic_json(HOME / "summary.json", {"layers": layers,
                     "all_available_targets_complete": all(x["n_pending"] == 0 and x["n_failed"] == 0 for x in layers),
                     "scope": "Macro means across completed, evaluable targets only; unavailable, failed and pending counts shown.",
                     "root": str(ROOT)})
        log("SUMMARY", [(x["layer"], x["n_complete"], x["n_pending"], x["n_unavailable"]) for x in layers])

def selftest():
    assert metric([0, 1], [1, 0])["auroc"] == 1.
    assert metric([1, 0], [1, 0])["auroc"] == 0.
    tie = metric([0, 0], [1, 0])
    assert tie["auroc"] == .5 and tie["ef1"] == 1.
    labels = np.array([1] + [0] * 2 + [1] * 57 + [0] * 190)
    assert len(labels) == 250
    assert np.isclose(metric(np.arange(250), labels)["ef1"], (1/3)/(58/250))
    assert metric([0, np.inf], [1, 0])["auroc"] == 1.
    assert metric([0, 1], [1, 1])["auroc"] is None
    log("METRIC_TESTS_PASSED")

def pilot():
    selftest()
    import torch
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA not available")
    manifest = json.loads((HOME / "manifest.json").read_text())
    r = min((r for r in manifest["targets"] if r["status"] == "ready"), key=lambda r: r["n"])
    source = Path(r["target_dir"])
    library = pd.read_csv(source / "library.csv.gz")
    trial = pd.concat([library[library.label == 1].head(1), library[library.label == 0].head(1)], ignore_index=True)
    td = ROOT / "pilot" / r["uniprot"]
    td.mkdir(parents=True, exist_ok=True)
    for name in ("pocket.pdb", "pocket.pdbqt", "box.json"):
        shutil.copy2(source / name, td / name)
    begun = time.time()
    process_chunk(td, trial, 0, 2)
    f = pd.read_csv(td / "chunks/00000000_00000002/results.csv.gz")
    if len(f) != 2 or not f.scoring_ok.all():
        raise RuntimeError("Pilot failed; full workers remain blocked")
    atomic_json(HOME / "pilot.json", {"status": "passed", "layer": r["layer"],
                "target": r["uniprot"], "n": 2, "elapsed_seconds": time.time() - begun,
                "gpu": torch.cuda.get_device_name(0),
                "note": "Two-ligand functionality check, not a performance estimate."})
    log("PILOT_PASSED")

def launch():
    path = HOME / "launch.json"
    if path.exists():
        raise RuntimeError("Already launched; inspect launch.json instead of submitting duplicates")
    pilot_id = sbatch("pilot")
    jobs = []
    for lane in range(NWORKERS):
        jobs.append({"lane": lane, "job_id": sbatch("worker", lane, dependency=pilot_id)})
    atomic_json(path, {"pilot": pilot_id, "workers": jobs})

def worker(lane):
    manifest = json.loads((HOME / "manifest.json").read_text())
    records = sorted((r for r in manifest["targets"] if r.get("lane") == lane), key=lambda r: r["n"])
    path = HOME / ("worker_" + str(lane) + ".json")
    deadline = time.monotonic() + 38 * 3600
    done, failures = 0, []
    with (HOME / ("worker_" + str(lane) + ".lock")).open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for r in records:
            td = Path(r["target_dir"])
            marker = HOME / "target_metrics" / (r["layer"] + "_" + r["uniprot"] + ".json")
            if marker.exists():
                done += 1
                continue
            try:
                library = pd.read_csv(td / "library.csv.gz")
                for start in range(0, len(library), CHUNK):
                    if time.monotonic() > deadline:
                        successor = sbatch("worker", lane, dependency=os.environ.get("SLURM_JOB_ID"))
                        atomic_json(path, {"status": "continuation_submitted", "job_id": os.environ.get("SLURM_JOB_ID"),
                                          "next_job_id": successor, "targets_complete": done, "failures": failures})
                        return
                    process_chunk(td, library, start, min(start + CHUNK, len(library)))
                    atomic_json(path, {"status": "running", "job_id": os.environ.get("SLURM_JOB_ID"),
                        "layer": r["layer"], "target": r["uniprot"], "next_ligand": min(start+CHUNK,len(library)),
                        "target_n": len(library), "targets_complete": done, "failures": failures, "updated": time.time()})
                summarize_target(r, library)
                done += 1
            except Exception as exc:
                failure = {"layer": r["layer"], "uniprot": r["uniprot"], "error": repr(exc),
                           "job_id": os.environ.get("SLURM_JOB_ID")}
                failures.append(failure)
                atomic_json(HOME / "target_failures" / (r["layer"] + "_" + r["uniprot"] + ".json"), failure)
                log("TARGET_FAILED", failure)
                refresh_summary()
                # Stop repeated infrastructure failures instead of failing the whole dataset silently.
                if len(failures) >= 3 and done == 0:
                    atomic_json(path, {"status": "failed", "targets_complete": done, "failures": failures})
                    raise RuntimeError("Three targets failed before any completed; inspect logs before resuming") from exc
        atomic_json(path, {"status": "complete" if not failures else "complete_with_failures",
                           "targets_complete": done, "failures": failures})

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["prepare", "pilot", "launch", "worker", "summary", "selftest"])
    parser.add_argument("lane", type=int, nargs="?")
    args = parser.parse_args()
    if args.mode == "worker":
        if args.lane is None or not 0 <= args.lane < NWORKERS:
            parser.error("worker requires lane 0..7")
        worker(args.lane)
    else:
        {"prepare": prepare, "pilot": pilot, "launch": launch,
         "summary": refresh_summary, "selftest": selftest}[args.mode]()
