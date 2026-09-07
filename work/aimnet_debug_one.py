#!/usr/bin/env python3
"""Multi-protein numerical interaction-energy diagnostics on prepared HiQBind inputs."""
import argparse, hashlib, json, time, os
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from scipy.spatial.distance import cdist, pdist
from aimnet.calculators import AIMNet2Calculator
from aimnet_debug_io import read_pocket, eval_systems
K = 23.060547830619

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", default="/home/xianyang/aimnet-score-reference/geometry")
    ap.add_argument("--reference", default="/home/xianyang/aimnet-score-reference/hiqbind.csv")
    ap.add_argument("--per-target", type=int, default=2)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reference-neighbors", action="store_true", help="Use explicit PyTorch distances to isolate neighbor backend errors")
    ap.add_argument("--tolerance-kcal", type=float, default=0.1)
    args = ap.parse_args()
    if args.reference_neighbors:
        from aimnet.calculators.calculator import AdaptiveNeighborList
        def reference_neighbors(self, positions, cell=None, pbc=None, batch_idx=None, fill_value=None):
            if cell is not None: raise ValueError("Reference neighbors only support nonperiodic inputs")
            n = len(positions)
            fill = n if fill_value is None else fill_value
            d = torch.cdist(positions,positions,compute_mode="donot_use_mm_for_euclid_dist")
            valid = (d < self.cutoff)
            valid.fill_diagonal_(False)
            if batch_idx is not None: valid &= batch_idx[:,None] == batch_idx[None,:]
            count = valid.sum(1).to(torch.int32)
            indices = torch.arange(n,device=positions.device,dtype=torch.int32).expand(n,n)
            indices = indices.masked_fill(~valid,fill).sort(dim=1).values
            return indices[:,:max(1,int(count.max()))].contiguous(), count, None
        AdaptiveNeighborList.__call__ = reference_neighbors
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    refs = pd.read_csv(args.reference)
    selected = []
    for directory in sorted(Path(args.geometry).iterdir()):
        if not directory.is_dir(): continue
        files = sorted(directory.glob("*_protein.pdb"))[:args.per_target]
        for p in files:
            ligand = p.with_name(p.name.replace("_protein.pdb", "_ligands.sdf"))
            selected.append((directory.name, p, ligand))
    models = ["aimnet2-2025", "/home/xianyang/t3-aimnet-models/aimnet2_wb97m_0.jpt"]
    rows, meta = [], {"torch": torch.__version__, "cuda": torch.cuda.is_available(),
        "device": torch.cuda.get_device_name(0), "tolerance_kcal_mol": args.tolerance_kcal,
        "selection": "first two lexicographic prepared complexes per protein; first SDF molecule",
        "charge_policy": "sum explicit PDB charge fields; ligand RDKit formal charge; complex sum",
        "caveat": "Blank PDB charge fields do not establish protonation correctness. Numerical checks, not affinity validation.",
        "models": {}, "reference_neighbors": args.reference_neighbors, "cuda_launch_blocking": os.environ.get("CUDA_LAUNCH_BLOCKING", "unset")}
    for model in models:
        label = Path(model).name
        try:
            calc = AIMNet2Calculator(model, device="cuda", compile_model=False)
            meta["models"][label] = {"metadata": dict(calc.metadata or {})}
            h2 = float(eval_systems(calc, [(np.array([[0,0,0],[.74,0,0]],np.float32), np.array([1,1]),0)])[0])
            meta["models"][label]["h2_eV"] = h2 if np.isfinite(h2) else None
        except Exception as e:
            meta["models"][label] = {"error": repr(e)}
            continue
        for target, protein, ligand in selected:
            row = {"model": label, "target": target, "complex": protein.name.replace("_protein.pdb",""),
                   "protein_path": str(protein), "ligand_path": str(ligand)}
            started = time.time()
            try:
                pxyz, pz, pq = read_pocket(protein)
                mols = [m for m in Chem.SDMolSupplier(str(ligand), removeHs=False) if m is not None]
                if not mols: raise ValueError("No valid ligand")
                m = mols[0]
                lxyz = np.asarray(m.GetConformer().GetPositions(), np.float32)
                lz = np.array([a.GetAtomicNum() for a in m.GetAtoms()])
                lq = int(Chem.GetFormalCharge(m))
                lines = [s for s in protein.read_text().splitlines() if s.startswith(("ATOM  ","HETATM"))]
                row.update(protein_atoms=len(pz), ligand_atoms=len(lz), ligand_records=len(mols),
                    protein_hydrogens=int((pz==1).sum()), ligand_hydrogens=int((lz==1).sum()),
                    protein_charge=int(pq), ligand_charge=lq, complex_charge=int(pq+lq),
                    explicit_pdb_charge_fields=sum(bool(s[78:80].strip()) for s in lines),
                    protein_sha256=hashlib.sha256(protein.read_bytes()).hexdigest(),
                    ligand_sha256=hashlib.sha256(ligand.read_bytes()).hexdigest(),
                    min_protein_distance_A=float(pdist(pxyz).min()),
                    min_ligand_distance_A=float(pdist(lxyz).min()),
                    min_contact_A=float(cdist(pxyz,lxyz).min()))
                if not np.isfinite(pxyz).all() or not np.isfinite(lxyz).all(): raise ValueError("Nonfinite coordinates")
                supported = set((calc.metadata or {}).get("implemented_species") or [])
                unsupported = sorted(set(map(int,np.r_[pz,lz])) - supported) if supported else []
                if unsupported: raise ValueError("unsupported elements: "+str(unsupported))
                def systems(px, lx):
                    return [(np.concatenate([px,lx]),np.r_[pz,lz],pq+lq),(px,pz,pq),(lx,lz,lq)]
                def energy(ss, batch=True):
                    return np.asarray(eval_systems(calc,ss) if batch else [eval_systems(calc,[s])[0] for s in ss],dtype=np.float64)
                def interaction(e): return float((e[0]-e[1]-e[2])*K)
                base = systems(pxyz,lxyz)
                separate = energy(base,False); batched = energy(base)
                row.update(complex_eV=float(separate[0]),protein_eV=float(separate[1]),ligand_eV=float(separate[2]),
                    interaction_kcal_mol=interaction(separate),
                    batch_delta_kcal_mol=interaction(batched)-interaction(separate),
                    max_batch_component_delta_kcal_mol=float(np.max(np.abs(batched-separate))*K))
                repeat = energy(base)
                row["repeat_delta_kcal_mol"] = interaction(repeat)-interaction(batched)
                angle = .731
                rot = np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]],np.float32)
                shifted = energy(systems(pxyz@rot + [3.25,-2.5,1.75],lxyz@rot + [3.25,-2.5,1.75]))
                row["rigid_transform_delta_kcal_mol"] = interaction(shifted)-interaction(batched)
                permuted = energy([(x[::-1].copy(),z[::-1].copy(),q) for x,z,q in base])
                row["permutation_delta_kcal_mol"] = interaction(permuted)-interaction(batched)
                swapped = energy(base[::-1])[::-1]
                row["batch_order_delta_kcal_mol"] = interaction(swapped)-interaction(batched)
                # Move ligand completely beyond protein bounding box. Diagnostic only:
                # global charge equilibration need not conserve separate fragment charges.
                separation = float(np.ptp(pxyz[:,0])+np.ptp(lxyz[:,0])+100)
                far = energy(systems(pxyz,lxyz-lxyz.mean(0)+pxyz.mean(0)+[separation,0,0]))
                row["far_interaction_kcal_mol"] = interaction(far)
                row["far_min_contact_A"] = float(cdist(pxyz,lxyz-lxyz.mean(0)+pxyz.mean(0)+[separation,0,0]).min())
                match = refs[(refs.target==target)&(refs.PDBID==row["complex"].split("_")[0])]
                if len(match)==1:
                    row["saved_reference_interaction_kcal_mol"] = float(match.iloc[0]["AIMNet2(Eint)"])
                    row["difference_from_saved_reference_kcal_mol"] = row["interaction_kcal_mol"]-row["saved_reference_interaction_kcal_mol"]
                keys=["batch_delta_kcal_mol","repeat_delta_kcal_mol","rigid_transform_delta_kcal_mol","permutation_delta_kcal_mol","batch_order_delta_kcal_mol"]
                row["finite"] = bool(np.isfinite(np.r_[separate,batched,repeat,shifted,permuted,swapped,far]).all())
                row["numerical_pass"] = row["finite"] and all(abs(row[k]) <= args.tolerance_kcal for k in keys)
                row["status"] = "pass" if row["numerical_pass"] else "fail"
            except Exception as e:
                row.update(status="error",error=repr(e))
            row["seconds"] = time.time()-started
            rows.append(row)
            pd.DataFrame(rows).to_csv(out/"results.csv",index=False)
            print(json.dumps(row,default=str),flush=True)
        del calc
        torch.cuda.empty_cache()
    meta.update(n_rows=len(rows),n_pass=sum(r["status"]=="pass" for r in rows),
                n_fail=sum(r["status"]=="fail" for r in rows),n_error=sum(r["status"]=="error" for r in rows))
    (out/"summary.json").write_text(json.dumps(meta,indent=2,default=str)+"\n")
    print(json.dumps(meta,default=str),flush=True)
if __name__=="__main__": main()

