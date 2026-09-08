import json,hashlib,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata
root=Path(__file__).resolve().parent
jobs=json.loads((root/"jobs.json").read_text())
for worker in range(8):
    assert not json.loads((root/f"worker_{worker}.json").read_text())["failures"]
def corr(x,y):
    ok=np.isfinite(x)&np.isfinite(y)
    return float(spearmanr(x[ok],y[ok]).statistic) if ok.sum()>2 and np.ptp(x[ok]) and np.ptp(y[ok]) else None
def metrics(f,col):
    f=f.sort_values(col,kind="stable");y=f.label.astype(int).to_numpy();n=len(y);a=int(y.sum());top=max(1,int(np.ceil(.01*n)))
    ranks=rankdata(-f[col].to_numpy())
    return dict(n=n,n_actives=a,ef1=float(y[:top].sum()/top/(a/n)) if a else None,
        auroc=float((ranks[y==1].sum()-a*(a+1)/2)/(a*(n-a))) if 0<a<n else None,
        active_top5=int(y[:5].sum()),active_top10=int(y[:10].sum()))
frames=[];records=[];changes=[]
for j in jobs:
    d=root/j["key"];p=d/j["output"];f=pd.read_csv(p);expected=pd.read_csv(j["interaction"])
    summary=json.loads(p.with_suffix(".summary.json").read_text())
    assert summary["n_errors"]==0,(j["key"],summary)
    assert len(f)==len(expected) and set(f.mol_id)==set(expected.mol_id) and not f.mol_id.duplicated().any(),j["key"]
    assert f.interaction_geometry.eq("fixed_pocket_aimnet2025_minimized_v2").all()
    assert np.isfinite(f[["aimnet_interaction_kcal_mol","aimnet2_composite_kcal_mol"]]).all().all()
    assert np.allclose((f.complex_bound_ev-f.pocket_ev-f.interaction_ligand_ev)*23.060547830619,f.aimnet_interaction_kcal_mol,atol=1e-6,rtol=0)
    assert np.allclose(f.aimnet_interaction_kcal_mol+f.desolvation_kcal_mol+f.local_strain_kcal_mol,f.aimnet2_composite_kcal_mol,atol=1e-6,rtol=0)
    assert f.interaction_model.eq("aimnet2-2025").all()
    assert np.allclose((f.gas_bound_ev-f.cpcm_bound_ev)*23.060547830619,f.desolvation_kcal_mol,atol=1e-6,rtol=0)
    f["benchmark"]=j["key"];frames.append(f)
    r=dict(benchmark=j["key"],n=len(f),complex_converged=int(f.complex_optimization_converged.sum()),cpcm_converged=int(f.optimization_converged.sum()))
    if j["key"].startswith("t3/"):
        r.update(layer=j["key"].split("/")[1],interaction_spearman=corr(-f.aimnet_interaction_kcal_mol,f.paff),
                 composite_spearman=corr(-f.aimnet2_composite_kcal_mol,f.paff),smina_spearman=corr(f.smina_score,f.paff),pre_min_interaction_spearman=corr(-f.pre_minimization_interaction_kcal_mol,f.paff))
    if "label" in f:
        for prefix,col in [("interaction","aimnet_interaction_kcal_mol"),("composite","aimnet2_composite_kcal_mol")]+([("retrieval","ligunity_rank")] if "ligunity_rank" in f else []):
            r.update({prefix+"_"+k:v for k,v in metrics(f,col).items()})
    if "allposes" in j:
        raw=pd.read_csv(d/"aimnet_interaction_allposes.csv")
        assert raw.error.fillna("").eq("").all()
        assert len(raw)==json.loads((d/"aimnet_interaction_bestpose.summary.json").read_text())["n_input_poses"]
        assert len(expected)==20
        best_values=raw.groupby("ligand_pdb").aimnet_interaction_kcal_mol.min()
        assert np.allclose(expected.aimnet_interaction_kcal_mol,best_values.reindex(expected.ligand_pdb).to_numpy(),atol=1e-6,rtol=0)
        r["n_allposes_scored"]=len(raw)
        oldbest=pd.read_csv(j["previous_best"])
        choices=expected[["ligand_pdb","pose_index"]].merge(oldbest[["ligand_pdb","pose_index"]],on="ligand_pdb",suffixes=("_new","_previous"))
        r["n_selected_poses_changed"]=int((choices.pose_index_new!=choices.pose_index_previous).sum())
        choices.to_csv(d/"pose_selection_changes.csv",index=False)
        r["pre_min_interaction_ef1"]=metrics(expected,"aimnet_interaction_kcal_mol")["ef1"]
    records.append(r)
    for col,name in [("aimnet_interaction_kcal_mol","interaction_ranking.csv"),("aimnet2_composite_kcal_mol","composite_ranking.csv")]:
        ranked=f.sort_values(col,kind="stable").copy();ranked["rank"]=range(1,len(ranked)+1);ranked.to_csv(d/name,index=False)
    oldpath=Path(j["previous_score"])
    if oldpath.exists():
        identifier="ligand_pdb" if "allposes" in j else "mol_id"
        cols=[identifier,"aimnet_interaction_kcal_mol","aimnet2_composite_kcal_mol"]
        aligned=f[cols].merge(pd.read_csv(oldpath)[cols],on=identifier,suffixes=("_new","_previous"))
        aligned["benchmark"]=j["key"];changes.append(aligned)
systems=pd.DataFrame(records);systems.to_csv(root/"per_system_metrics.csv",index=False)
full=pd.concat(frames,ignore_index=True);full.to_csv(root/"all_ligand_scores.csv.gz",index=False)
pd.concat(changes,ignore_index=True).to_csv(root/"protocol_comparison.csv.gz",index=False)
t3={}
for layer,g in systems[systems.benchmark.str.startswith("t3/")].groupby("layer"):
    t3[layer]=dict(n_systems=len(g),n_ligands=int(g.n.sum()),n_evaluable=int(g.composite_spearman.notna().sum()),
        mean_interaction_spearman=float(g.interaction_spearman.mean()),mean_composite_spearman=float(g.composite_spearman.mean()),
        mean_smina_spearman=float(g.smina_spearman.mean()),mean_pre_min_interaction_spearman=float(g.pre_min_interaction_spearman.mean()))
c=systems[systems.benchmark.str.startswith("casf/")]
casf=dict(n_pockets=len(c),n_ligands=int(c.n.sum()),pose_selection="All available poses rescored with AIMNet2(2025) member 0; lowest interaction-energy pose selected per ligand, then minimized",
    total_allposes_scored=int(c.n_allposes_scored.sum()),selected_poses_changed=int(c.n_selected_poses_changed.sum()),
    pre_min_interaction_ef1=float(c.pre_min_interaction_ef1.mean()),
    **{col:float(c[col].mean()) for col in ["interaction_ef1","composite_ef1","retrieval_ef1","interaction_active_top5","composite_active_top5","retrieval_active_top5","interaction_active_top10","composite_active_top10","retrieval_active_top10"]})
assert casf["n_pockets"]==24 and casf["n_ligands"]==480 and casf["total_allposes_scored"]==47896
assert systems.loc[systems.benchmark.eq("casf-single"),"n_allposes_scored"].iloc[0]==2000
assert sum(v["n_systems"] for v in t3.values())==93
smoke=full[full.benchmark.str.startswith("smoke/")].sort_values("source_index").copy()
assert len(smoke)==250 and not smoke.mol_id.duplicated().any()
assert len(full)==1667, len(full)
for col in ["aimnet_interaction_kcal_mol","aimnet2_composite_kcal_mol"]:smoke[col.replace("_kcal_mol","_rank")]=rankdata(smoke[col])
smoke["aimnet2_score_rank"]=rankdata(smoke.aimnet2_composite_kcal_mol)
smoke["aimnet_interaction_rank"]=rankdata(smoke.aimnet_interaction_kcal_mol)
smoke.to_csv(root/"smoke/aimnet2_score.csv",index=False)
for col,name in [("aimnet_interaction_kcal_mol","interaction_ranking.csv"),("aimnet2_composite_kcal_mol","composite_ranking.csv")]:
    ranked=smoke.sort_values(col,kind="stable").copy()
    ranked["rank"]=range(1,len(ranked)+1)
    ranked.to_csv(root/"smoke"/name,index=False)
for col in ["aimnet2_score_rank","aimnet_interaction_rank"]:
    full.loc[smoke.index,col]=smoke[col]
full.to_csv(root/"all_ligand_scores.csv.gz",index=False)
systems[systems.benchmark.str.startswith("t3/")].to_csv(root/"t3/per_target_metrics.csv",index=False)
systems[systems.benchmark.str.startswith("casf/")].to_csv(root/"casf/per_pocket_metrics.csv",index=False)
summary=dict(protocol="fixed_pocket_aimnet2025_minimized_v2",n_scored=len(full),n_system_jobs=len(jobs),
    n_complex_converged=int(full.complex_optimization_converged.sum()),n_cpcm_converged=int(full.optimization_converged.sum()),
    n_both_converged=int((full.complex_optimization_converged & full.optimization_converged).sum()),
    n_complex_energy_increased=int((full.complex_bound_ev>full.complex_initial_ev+1e-6).sum()),
    t3=t3,casf=casf,smoke={name:metrics(smoke,col) for name,col in [("pre_min_interaction","pre_minimization_interaction_kcal_mol"),("interaction","aimnet_interaction_kcal_mol"),("composite","aimnet2_composite_kcal_mol")]},
    excluded_inputs=["L3/P18564: zero-byte poses","L4/Q15726: zero-byte poses","L4/Q8TDU9: zero-byte poses"],
    validation="All input ligand IDs matched; no scoring errors; energy arithmetic passed at 1e-6 kcal/mol",
    convergence_policy="All finite scores included; nonconverged rows flagged",
    comparison_note="Interaction and fixed-pocket minimization use AIMNet2(2025) member 0. Desolvation uses the separate legacy wB97M gas and wB97M-CPCM pair. CASF comparisons also include refreshed pose selection.")
(root/"aggregate_summary.json").write_text(json.dumps(summary,indent=2))
(root/"smoke/aggregate_summary.json").write_text(json.dumps(summary["smoke"],indent=2))
lines=["# AIMNet2(2025) member-0 benchmark results","",f"Scored {len(full)} ligand cases. Complex minimization converged for {summary['n_complex_converged']}; isolated CPCM minimization converged for {summary['n_cpcm_converged']}.",
"","Composite = minimized complex interaction + bound desolvation + local ligand strain. Protein atoms remain fixed; all bound terms use the same minimized ligand geometry.",
"","## T3 affinity ranking","","| Layer | Systems | Ligands | Interaction Spearman | Composite Spearman |","|---|---:|---:|---:|---:|"]
for layer,v in t3.items():lines.append(f"| {layer} | {v['n_systems']} | {v['n_ligands']} | {v['mean_interaction_spearman']:.4f} | {v['mean_composite_spearman']:.4f} |")
lines+=["","Correlations compare negative energy with pAffinity. Undefined correlations are excluded from means.",
"","## CASF","",f"24 pockets, 480 freshly selected and minimized poses. Mean shortlist EF1: interaction {casf['interaction_ef1']:.4f}; composite {casf['composite_ef1']:.4f}; retrieval {casf['retrieval_ef1']:.4f}.",
f"Pre-minimization interaction EF1 was {casf['pre_min_interaction_ef1']:.4f}. Pose selection changed for {casf['selected_poses_changed']} of 480 ligands.",
"All 47,896 available poses were rescored with AIMNet2(2025), then one pose per ligand was selected for fixed-pocket minimization. The separate 20-ligand example is included in per_system_metrics.csv outside the 24-pocket average.",
"","## 250-ligand smoke test",""]
for name,v in summary["smoke"].items():lines.append(f"{name}: AUROC {v['auroc']:.4f}; EF1 {v['ef1']:.4f}.")
lines+=["","## Validation and coverage","",summary["validation"]+".",
"Three zero-byte T3 inputs remain excluded. Finite nonconverged scores are retained and flagged. The smoke set is nonrandom; T3 is active-only. These are public-checkpoint reconstructed scores.",
"","Each system includes component energies, minimized SDF poses, and interaction/composite ranking CSVs. protocol_comparison.csv.gz aligns old/new scores.", "", summary["comparison_note"], "", f"Both minimizations converged for {summary['n_both_converged']} of {len(full)} cases. {summary['n_complex_energy_increased']} complexes ended above their initial energy."]
previous=json.loads((root/"previous_checkpoint_summary.json").read_text())
lines += ["", "## Change from previous checkpoint", "",
"| Metric | Previous | AIMNet2(2025) member 0 |", "|---|---:|---:|"]
for layer, value in t3.items():
    lines.append(f"| T3 {layer} composite mean Spearman | {previous['t3'][layer]['mean_composite_spearman']:.4f} | {value['mean_composite_spearman']:.4f} |")
for name in ["interaction", "composite"]:
    lines.append(f"| CASF {name} mean shortlist EF1 | {previous['casf'][name+'_ef1']:.4f} | {casf[name+'_ef1']:.4f} |")
    lines.append(f"| Smoke {name} AUROC | {previous['smoke'][name]['auroc']:.4f} | {summary['smoke'][name]['auroc']:.4f} |")
(root/"report.md").write_text("\n".join(lines)+"\n")
files=[p for p in root.rglob("*") if p.is_file() and p.suffix not in (".log",".pyc") and "__pycache__" not in p.parts and p.name!="SHA256SUMS" and "test" not in p.relative_to(root).parts and not any(part in ("split0","split1","split2","split3") for part in p.parts) and not p.name.startswith("completed_before_split")]
(root/"SHA256SUMS").write_text("".join(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+str(p.relative_to(root))+"\n" for p in sorted(files)))
with zipfile.ZipFile(root.with_suffix(".zip"),"w",zipfile.ZIP_DEFLATED) as z:
    for p in files+[root/"SHA256SUMS"]:z.write(p,p.relative_to(root))
print(json.dumps(summary,indent=2))
