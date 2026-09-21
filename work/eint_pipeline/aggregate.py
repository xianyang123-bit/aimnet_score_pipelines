"""Aggregate all four layers from the unified fixed-receptor Eint run."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from run import case_path,load_tasks,save

def rho(frame,column):
    if len(frame)<3 or frame.paff.nunique()<2:return None
    x=float(spearmanr(-frame[column],frame.paff).statistic);return x if np.isfinite(x) else None
def aggregate(root):
    coverage=[];frames=[]
    for task in load_tasks(root):
        case=case_path(root,task);status={}
        details={}
        for stage in ('preparation','docking','scoring'):
            p=case/(stage+'_status.json')
            if p.exists():details[stage]=json.loads(p.read_text());status[stage]=details[stage].get('status')
        row=dict(index=task['index'],layer=task['layer'],uniprot=task['uniprot'],n_requested=task['n_ligands'],**status)
        prep=details.get('preparation',{})
        if prep.get('status')=='excluded_chemistry':
            residues=', '.join(x.get('residue','') for x in prep.get('blockers',[]))
            row['exclusion_reason']=prep.get('blocked_reason','excluded chemistry')+((': '+residues) if residues else '')
        elif prep.get('status')=='preparation_failed':row['exclusion_reason']=prep.get('error','preparation failed')
        elif prep.get('status')=='source_not_recovered':row['exclusion_reason']='Exact residue-labelled source structure was not recovered'
        elif prep.get('status')=='empty_original_cohort':row['exclusion_reason']='Original benchmark pose input was empty'
        else:row['exclusion_reason']=''
        score=case/'eint.csv'
        if score.exists():
            f=pd.read_csv(score);f=f[f.error.fillna('').eq('')].copy();frames.append(f);row.update(n_scored=len(f),n_converged=int(f.converged.sum()),eint_rho=rho(f,'aimnet_interaction_kcal_mol'),smina_rho=rho(f,'smina_affinity_kcal_mol'))
        else:row.update(n_scored=0,n_converged=0,eint_rho=None,smina_rho=None)
        coverage.append(row)
    cov=pd.DataFrame(coverage);cov.to_csv(root/'coverage.csv',index=False)
    scores=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
    if len(scores):
        scores['eint_rank']=scores.groupby(['layer','uniprot']).aimnet_interaction_kcal_mol.rank(method='average')
        scores['smina_rank']=scores.groupby(['layer','uniprot']).smina_affinity_kcal_mol.rank(method='average')
    scores.to_csv(root/'all_scores.csv',index=False)
    metrics=[]
    for layer,g in cov.groupby('layer'):
        usable=g[g.eint_rho.notna()]
        metrics.append(dict(layer=layer,n_systems=int((g.n_scored>0).sum()),n_evaluable=len(usable),n_ligands=int(g.n_scored.sum()),
                            mean_eint_rho=float(usable.eint_rho.mean()) if len(usable) else None,mean_smina_rho=float(usable.smina_rho.mean()) if len(usable) else None))
    pd.DataFrame(metrics).to_csv(root/'layer_metrics.csv',index=False)
    decomposition=float(abs((scores.complex_ev-scores.pocket_ev-scores.ligand_ev)*23.060547830619-scores.aimnet_interaction_kcal_mol).max())
    validation=dict(n_scoring_errors=int(scores.error.fillna('').ne('').sum()),n_converged=int(scores.converged.sum()),
                    n_nonconverged=int((~scores.converged).sum()),max_decomposition_error_kcal_mol=decomposition,
                    n_final_complex_energy_increases=int((scores.complex_ev>scores.initial_complex_ev+1e-5).sum()),
                    minimum_initial_heavy_contact_A=float(scores.initial_min_heavy_contact_A.min()),minimum_final_heavy_contact_A=float(scores.final_min_heavy_contact_A.min()))
    assert validation['n_scoring_errors']==0 and validation['n_final_complex_energy_increases']==0 and decomposition<1e-6
    summary=dict(n_requested_systems=len(cov),n_requested_ligands=int(cov.n_requested.sum()),n_scored_systems=int((cov.n_scored>0).sum()),
                 n_scored_ligands=int(cov.n_scored.sum()),layers=metrics,preparation_status=cov.preparation.value_counts(dropna=False).to_dict(),validation=validation)
    save(root/'summary.json',summary)
    lines=['# Unified L1–L4 fixed-receptor AIMNet2 Eint run','',
           'This run retains only the previous ligand all-atom relaxation in a fixed, fully protonated receptor followed by interaction-energy scoring. Fresh graph-derived ligand conformers were docked with SMINA before scoring. Composite, CPCM, ligand-strain, and hydrogen-only scoring are not part of this run.','',
           f"Scored **{summary['n_scored_ligands']} ligands across {summary['n_scored_systems']} chemically supported systems** from the requested 96-system, 917-ligand L1–L4 manifest. Unsupported chemistry, unavailable exact source structures, unsafe preparation, and the three original empty inputs are retained in `coverage.csv` with reasons.",'',
           '| Layer | Systems | Ligands | Mean Eint Spearman ρ | Mean SMINA Spearman ρ |','|---|---:|---:|---:|---:|']
    for x in metrics:lines.append(f"| {x['layer']} | {x['n_systems']} | {x['n_ligands']} | {x['mean_eint_rho']:+.4f} | {x['mean_smina_rho']:+.4f} |")
    lines += ['','Eint is `Ecomplex − Ereceptor − Eligand`, with receptor atoms fixed and the ligand minimized as a whole in the receptor field. Fragment energies use the final complex geometry. Lower Eint receives the better rank. Spearman values are unweighted means of within-target correlations between negative energy and experimental pAffinity.','',
              f"Validation: {validation['n_scoring_errors']} scoring errors; {validation['n_converged']}/{summary['n_scored_ligands']} ligand minimizations converged below 0.002 eV/Å; {validation['n_nonconverged']} finite nonconverged results are retained and flagged. No final complex energy increased. Maximum Eint decomposition error was {decomposition:.3g} kcal/mol. Minimum protein–ligand heavy-atom contact was {min(validation['minimum_initial_heavy_contact_A'],validation['minimum_final_heavy_contact_A']):.3f} Å.",'',
              '`all_scores.csv` contains ligand-level energies, ranks, convergence, contacts, and labels. `coverage.csv` contains every requested target and exclusion reason. `layer_metrics.csv`, `summary.json`, `protocol.json`, `tasks.json`, `code_sha256.json`, and `removed_L1_code.json` record aggregate results and provenance.','']
    (root/'README.md').write_text('\n'.join(lines))
    print(json.dumps(json.loads((root/'summary.json').read_text()),indent=2))
