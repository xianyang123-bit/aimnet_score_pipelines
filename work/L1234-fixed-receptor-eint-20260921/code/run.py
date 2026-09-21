#!/usr/bin/env python3
"""One manifest-driven, corrected-receptor Eint pipeline for L1, L2, L3 and L4."""
import argparse
import concurrent.futures
import hashlib
import json
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

BASE=Path('/home/xianyang/aimnet2_score_pipelines/work')
BIN=Path('/home/xianyang/miniconda3/envs/t3-aimnet/bin')
MODEL=Path('/home/xianyang/.cache/aimnet/aimnet2_2025_b973c_d3_0.pt')
MODEL_SHA='043ed5418a104e31f79462f8e5ebeca64a2d24422174f5d29f894d32271981b5'

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(path,data):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def get(url):
    req=urllib.request.Request(url,headers={'User-Agent':'AIMNet-corrected-pocket-benchmark/1.0'})
    with urllib.request.urlopen(req,timeout=45) as f:return f.read()
def load_tasks(root):return json.loads((root/'tasks.json').read_text())
def case_path(root,task):return root/task['layer']/task['uniprot']
def run_command(cmd,log,timeout=1800):
    with Path(log).open('w') as f:subprocess.run([str(x) for x in cmd],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)

def initialize(root):
    import pandas as pd
    from rdkit import Chem
    root.mkdir(exist_ok=False,parents=True);(root/'logs').mkdir()
    assert digest(MODEL)==MODEL_SHA,'Interaction checkpoint mismatch'
    tasks=[]
    archive=BASE/'t3-aimnet-affinity-balanced-poses-and-scores.tar.gz'
    with tarfile.open(archive) as tar:
        systems=sorted({tuple(Path(m.name).parts[:2]) for m in tar.getmembers() if m.isfile() and len(Path(m.name).parts)==3 and Path(m.name).name=='pocket.pdb'})
        assert len(systems)==96 and {x[0] for x in systems}=={'L1','L2','L3','L4'}
        for layer,up in systems:
            task=dict(index=len(tasks),layer=layer,uniprot=up);case=case_path(root,task);case.mkdir(parents=True)
            for name,dest in [('pocket.pdb','legacy_pocket.pdb'),('poses.sdf','original_poses.sdf'),('box.json','box.json')]:
                (case/dest).write_bytes(tar.extractfile(f'{layer}/{up}/{name}').read())
            old=BASE/'aimnet2025-benchmarks-20260906/t3'/layer/up/'aimnet2_score.csv'
            if old.exists():
                labels=pd.read_csv(old);shutil.copy2(old,case/'previous_score.csv')
                mols=list(Chem.SDMolSupplier(str(case/'original_poses.sdf'),removeHs=False))
                assert all(m is not None for m in mols)
                byid={int(m.GetProp('mol_id')):m for m in mols}
                assert labels.mol_id.is_unique and len(byid)==len(mols)
                for row in labels.itertuples():
                    m=byid[int(row.mol_id)]
                    assert abs(float(m.GetProp('paff'))-row.paff)<1e-7
                    assert Chem.GetFormalCharge(m)==row.formal_charge and m.GetNumAtoms()==row.n_atoms
                keep=['mol_id','paff','smiles','formal_charge','n_atoms']
                labels[keep].to_csv(case/'labels.csv',index=False)
                task.update(n_ligands=len(labels),mol_ids=labels.mol_id.astype(int).tolist())
            else:task.update(n_ligands=0,mol_ids=[])
            task['input_sha256']={p.name:digest(p) for p in case.iterdir() if p.is_file()}
            # Cache source structures and explicit chemistry decisions; never reuse scores.
            candidates=[BASE/'L1-corrected-rerank-20260921'/layer/up,BASE/'pocket-correction-20260921/recovered'/layer/up]
            for previous in candidates:
                if not previous.exists():continue
                for p in previous.glob('*.cif'):
                    if not (case/p.name).exists():shutil.copy2(p,case/p.name)
                if (previous/'uniprot.json').exists():shutil.copy2(previous/'uniprot.json',case/'uniprot.json')
            tasks.append(task)
    save(root/'tasks.json',tasks)
    review=json.loads((BASE/'L1-corrected-rerank-20260921/chemistry_review.json').read_text())
    save(root/'chemistry_review.json',{f'L1/{up}':value for up,value in review.items()})
    save(root/'protocol.json',dict(layers=['L1','L2','L3','L4'],n_systems=len(tasks),n_ligands=sum(t['n_ligands'] for t in tasks),
        input_archive=str(archive),archive_sha256=digest(archive),model=str(MODEL),model_sha256=MODEL_SHA,
        docking='Fresh ETKDGv3/MMFF conformers; SMINA exhaustiveness 8, seed 1, one pose; restored explicit ligand H',
        receptor='Validated complete capped protein residues, explicit H, pH 7.4, integer topology charge; all receptor atoms frozen',
        optimization='Previous ligand all-atom FIRE unchanged; maximum 1000 steps, max force 0.002 eV/A',
        score='Eint = complex - fixed receptor - ligand at its complex geometry; no separately relaxed fragments',
        unsupported='Explicit exclusions; no neutral-charge default or hydrogen-free receptor fallback',
        empty_inputs='The three original empty-pose systems remain documented without inventing replacement cohorts'))
    print(json.dumps(dict(root=str(root),systems=len(tasks),ligands=sum(t['n_ligands'] for t in tasks))),flush=True)

def recover(root,task):
    import numpy as np
    from openmm import app,unit
    from scipy.spatial import cKDTree
    from prepare_aimnet_pocket import AA,legacy_atoms
    from source_inventory import inventory
    case=case_path(root,task);report_path=case/'source_recovery.json'
    if not task['n_ligands']:
        save(report_path,dict(status='empty_original_cohort'));return
    if report_path.exists() and json.loads(report_path.read_text()).get('selected_source'):return
    ref=legacy_atoms(case/'legacy_pocket.pdb')
    cache=case/'uniprot.json'
    if not cache.exists():cache.write_bytes(get('https://rest.uniprot.org/uniprotkb/'+task['uniprot']+'.json'))
    data=json.loads(cache.read_text())
    ids={x['id'] for x in data.get('uniProtKBCrossReferences',[]) if x['database']=='PDB'}
    ids|={p.stem.upper() for p in case.glob('*.cif')}
    ids=sorted(ids,key=lambda p:(not (case/(p+'.cif')).exists(),p))
    report=json.loads(report_path.read_text()) if report_path.exists() else dict(status='searching',reference_sha256=digest(case/'legacy_pocket.pdb'),candidates=len(ids),results=[])
    checked={x['pdb_id'] for x in report['results'] if 'error' not in x}
    ids=[x for x in ids if x not in checked]
    def check(pdb_id):
        path=case/(pdb_id+'.cif')
        try:
            if not path.exists():path.write_bytes(get('https://files.rcsb.org/download/'+pdb_id+'.cif'))
            source=app.PDBxFile(str(path));atoms=[a for a in source.topology.atoms() if a.residue.name in AA and a.element.atomic_number!=1]
            xyz=np.asarray(source.positions.value_in_unit(unit.angstrom));tree=cKDTree(xyz[[a.index for a in atoms]])
            count=sum(len([atoms[i] for i in tree.query_ball_point(c,.02) if atoms[i].name==n and atoms[i].element.symbol.upper()==e.upper()])==1 for n,e,c in ref)
            return dict(pdb_id=pdb_id,matched_atoms=count,reference_atoms=len(ref),exact_match=count==len(ref),source=str(path),sha256=digest(path))
        except Exception as e:return dict(pdb_id=pdb_id,error=repr(e))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        for start in range(0,len(ids),2):
            report['results'].extend(ex.map(check,ids[start:start+2]));matches=[x for x in report['results'] if x.get('exact_match')]
            if matches:
                report.update(status='matched',selected_source=matches[0]['source'])
                report['heterogen_inventory']=inventory(Path(matches[0]['source']),case/'legacy_pocket.pdb')
            save(report_path,report)
            if matches:break
            print(task['layer'],task['uniprot'],'checked',len(report['results']),'/',report['candidates'],flush=True)
    if not report.get('selected_source'):report['status']='source_not_recovered';save(report_path,report)
    print(task['layer'],task['uniprot'],report['status'],flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage',choices=['init','pipeline','recover','prepare','dock','score','aggregate'])
    ap.add_argument('--run',type=Path,required=True);ap.add_argument('--index',type=int)
    ap.add_argument('--worker',type=int);ap.add_argument('--workers',type=int,default=8)
    args=ap.parse_args();root=args.run.resolve()
    if args.stage=='init':initialize(root);return
    if args.stage=='aggregate':
        from aggregate import aggregate
        aggregate(root);return
    tasks=load_tasks(root)
    if args.index is not None:selected=[tasks[args.index]]
    elif args.worker is not None and 0<=args.worker<args.workers:selected=tasks[args.worker::args.workers]
    else:ap.error('Supply --index or valid --worker/--workers')
    failures=[]
    for task in selected:
        case=case_path(root,task)
        from stages import prepare_target,dock_target,score_target
        stages=['recover','prepare','dock'] if args.stage=='pipeline' else [args.stage]
        for stage in stages:
            try:{'recover':recover,'prepare':prepare_target,'dock':dock_target,'score':score_target}[stage](root,task)
            except Exception as e:
                error=dict(stage=stage,layer=task['layer'],uniprot=task['uniprot'],error=repr(e))
                save(case/(stage+'_error.json'),error);failures.append(error);print(json.dumps(error),flush=True);break
    if failures:raise RuntimeError(f'{len(failures)} systems failed; see per-system error records')

if __name__=='__main__':main()
