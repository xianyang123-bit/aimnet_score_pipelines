"""Thirty-ligand T3 screening: five actives and twenty-five decoys per target."""
import importlib.util as _util
import sys as _sys
from pathlib import Path as _Path
_home = _Path(__file__).resolve().parent
_spec = _util.spec_from_file_location("_t3_small_base", _home / "legacy_run.py")
_base = _util.module_from_spec(_spec)
_spec.loader.exec_module(_base)
_base.HOME = _home
_base.ROOT = _Path("/data/user_data/xianyang/t3-smoke-30-L123-20260908")
_base.CODE = _home / "protocol"
_base.CHUNK = 10
_sys.path.insert(0, str(_base.CODE))
globals().update({k: v for k, v in vars(_base).items() if not k.startswith("_")})
SOURCE = Path("/home/xianyang/aimnet2_score_pipelines/work/t3-smoke-L123-20260908")

def prepare_small():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    if (HOME / "manifest.json").exists():
        raise RuntimeError("Already prepared; resume fast workers instead")
    original = json.loads((SOURCE / "manifest.json").read_text())
    ROOT.mkdir(parents=True, exist_ok=True)
    records = []
    cached_score_count, cached_pose_count = 0, 0
    for src in original["targets"]:
        origin = Path(src["target_dir"])
        frame = pd.read_csv(origin / "library.csv.gz").sort_values("mol_id").reset_index(drop=True)
        active, decoy = frame[frame.label == 1], frame[frame.label == 0]
        if len(active) < 5 or len(decoy) < 25:
            raise RuntimeError("Insufficient classes for 5+25 at " + str(origin))
        seed = int.from_bytes(hashlib.sha256(f"20260908:30:{src['layer']}:{src['uniprot']}".encode()).digest()[:8], "little")
        rng = np.random.default_rng(seed)
        selected = pd.concat([active.iloc[rng.choice(len(active), 5, replace=False)],
                              decoy.iloc[rng.choice(len(decoy), 25, replace=False)]], ignore_index=True)
        selected = selected.iloc[rng.permutation(30)].reset_index(drop=True)
        assert len(selected) == 30 and selected.mol_id.is_unique and int(selected.label.sum()) == 5
        repeat = np.random.default_rng(seed)
        ids = set(active.iloc[repeat.choice(len(active), 5, replace=False)].mol_id)
        ids.update(decoy.iloc[repeat.choice(len(decoy), 25, replace=False)].mol_id)
        assert ids == set(selected.mol_id)
        score_parts, pose_map = [], {}
        for cd in (origin / "chunks").glob("*"):
            if (cd / "complete.json").exists():
                scored = pd.read_csv(cd / "results.csv.gz")
                scored = scored[scored.mol_id.isin(ids) & np.isfinite(scored.aimnet2_composite_kcal_mol)]
                if len(scored):
                    score_parts.append(scored)
            pose_file = None
            if (cd / "docked.json").exists():
                pose_file = Path(json.loads((cd / "docked.json").read_text())["poses"])
            else:
                for attempt in sorted(cd.glob("attempt-*"), reverse=True):
                    if (attempt / "score.log").exists() and (attempt / "poses.sdf").exists():
                        pose_file = attempt / "poses.sdf"
                        break
            if pose_file is not None and pose_file.exists() and pose_file.stat().st_size:
                for mol in Chem.SDMolSupplier(str(pose_file), removeHs=False):
                    if mol is None:
                        continue
                    mid = int(mol.GetProp("mol_id")) if mol.HasProp("mol_id") else int(mol.GetProp("_Name").rsplit(":", 1)[-1])
                    if mid in ids:
                        pose_map[mid] = mol
        cache = pd.concat(score_parts, ignore_index=True) if score_parts else pd.DataFrame(columns=["mol_id"])
        if cache.mol_id.duplicated().any():
            raise RuntimeError("Duplicate cached score IDs at " + str(origin))
        cached_ids = set(cache.mol_id)
        # Storage order groups reusable scores; ligand selection never uses scores or completion.
        selected["_cached"] = selected.mol_id.isin(cached_ids)
        selected = selected.sort_values("_cached", ascending=False, kind="stable").drop(columns="_cached").reset_index(drop=True)
        selected = selected.rename(columns={"source_index": "previous_sample_index"})
        selected.insert(0, "source_index", np.arange(30))
        td = ROOT / src["layer"] / src["uniprot"]
        td.mkdir(parents=True, exist_ok=True)
        atomic_csv(td / "library.csv.gz", selected)
        if pose_map:
            writer = Chem.SDWriter(str(td / "cached_poses.sdf"))
            for mid in selected.mol_id:
                if mid in pose_map:
                    writer.write(pose_map[mid])
            writer.close()
        if len(cache):
            atomic_csv(td / "cached_scores.csv.gz", cache)
        record = dict(src)
        record.update(n=30, n_actives=5, n_decoys=25, n_chunks=3, target_dir=str(td),
                      sampling_seed=str(seed), seed_base=20260908, sampled_from_n=len(frame),
                      previous_target_dir=str(origin), n_cached_scores=len(cache), n_cached_poses=len(pose_map))
        if src["status"] == "ready":
            for name in ("pocket.pdb", "pocket.pdbqt", "box.json"):
                shutil.copy2(origin / name, td / name)
        for start in range(0, 30, CHUNK):
            subset = selected.iloc[start:start+CHUNK].copy()
            if not set(subset.mol_id).issubset(cached_ids):
                continue
            overlap = [c for c in cache.columns if c in subset.columns and c != "mol_id"]
            joined = subset.merge(cache.drop(columns=overlap), on="mol_id", validate="one_to_one")
            assert len(joined) == CHUNK and np.isfinite(joined.aimnet2_composite_kcal_mol).all()
            cd = td / "chunks" / f"{start:08d}_{start+CHUNK:08d}"
            cd.mkdir(parents=True, exist_ok=True)
            atomic_csv(cd / "results.csv.gz", joined)
            atomic_json(cd / "complete.json", {"n_sampled": CHUNK, "n_scored": CHUNK,
                        "n_failures": 0, "reused_scores": True, "source_target": str(origin), "attempt": str(td)})
        cached_score_count += len(cache)
        cached_pose_count += len(pose_map)
        atomic_json(td / "target.json", record)
        records.append(record)
    manifest = dict(original)
    manifest.update(root=str(ROOT), scope="Exactly 5 actives and 25 decoys per T3 L1-L3 target",
        source_manifest=str(SOURCE / "manifest.json"), targets=records, chunk_size=CHUNK,
        sampling="Without replacement from the preceding reproducible 250-ligand sample; fixed 5:25 class counts.",
        seed_derivation="First 8 SHA256 bytes of 20260908:30:layer:uniprot, little endian; NumPy default_rng.",
        n_cached_scores=cached_score_count, n_cached_poses=cached_pose_count)
    summary = pd.DataFrame(records).groupby("layer").agg(targets=("uniprot","size"), sampled_ligands=("n","sum"))
    manifest["layer_counts"] = json.loads(summary.to_json(orient="index"))
    atomic_json(HOME / "manifest.json", manifest)
    atomic_json(ROOT / "manifest.json", manifest)
    atomic_csv(HOME / "target_inventory.csv", pd.DataFrame(records))
    refresh_summary()
    atomic_json(HOME / "preparation_summary.json", {"targets":len(records),
        "ready_targets":sum(r["status"]=="ready" for r in records),
        "n_to_score":sum(r["n"] for r in records if r["status"]=="ready"),
        "cached_scores":cached_score_count,"cached_poses":cached_pose_count,
        "selection_checks":"Unique IDs, exact 5+25 counts and repeated seeded selection passed"})
    import fast
    fast.prepare_fast()

if __name__ == "__main__":
    prepare_small()
