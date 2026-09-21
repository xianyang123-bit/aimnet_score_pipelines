#!/usr/bin/env python3
"""Rank the 16 FEP series from their deposited aligned poses with fixed-receptor AIMNet2 Eint."""
import argparse,glob,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
FEP=Path('/home/xianyang/LigUnity-official/test_datasets/FEP')
MODEL=Path('/home/xianyang/.cache/aimnet/aimnet2_2025_b973c_d3_0.pt')
MODEL_SHA='043ed5418a104e31f79462f8e5ebeca64a2d24422174f5d29f894d32271981b5'
AA=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
METALS={'LI','NA','MG','AL','K','CA','MN','FE','CO','NI','CU','ZN','CD','HG'}
MODIFIED={'TPO','SEP','PTR','MSE','CSO','CSD','KCX','HYP','MLY','TYS'}
COFACTORS={'NAD','NAP','NDP','NAI','FAD','FMN','SAM','SAH','HEM','HEC','COH','PLP','ATP','ADP','GDP','GTP','F6P','FLC','POP'}

def save(path,data):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def taskdir(root,t):return root/t['system']
def protein_path(name):
 exact=FEP/'structures'/f'{name}_protein.pdb'
 if exact.exists():return exact
 x=list((FEP/'structures').glob(f'{name}*protein.pdb'));assert len(x)==1,(name,x);return x[0]
def ligand_path(name):
 x=list((FEP/'structures').glob(f'{name}*ligands.sdf'));assert len(x)==1,(name,x);return x[0]

def topology_key(m):
 from rdkit import Chem
 m=Chem.RemoveHs(Chem.Mol(m));rw=Chem.RWMol(m)
 for a in rw.GetAtoms():
  a.SetFormalCharge(0);a.SetIsAromatic(False);a.SetNoImplicit(True);a.SetNumExplicitHs(0);a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
 for b in rw.GetBonds():b.SetBondType(Chem.BondType.SINGLE);b.SetIsAromatic(False);b.SetStereo(Chem.BondStereo.STEREONONE)
 return Chem.MolToSmiles(rw,canonical=True,isomericSmiles=False)

def match_ligands(label,mols):
 from rdkit import Chem
 used=set();out=[];audit=[]
 keys=[topology_key(m) for m in mols]
 for j,lig in enumerate(label['ligands']):
  expected=-1.364*float(lig['act']);key=topology_key(Chem.MolFromSmiles(lig['smi']))
  dg=[i for i,m in enumerate(mols) if i not in used and m.HasProp('r_exp_dg') and abs(float(m.GetProp('r_exp_dg'))-expected)<0.005]
  exact=[i for i in dg if keys[i]==key];candidates=exact or dg
  if not candidates:raise ValueError(f"No deposited pose matches {label['pockets'][0]} ligand {j}")
  idx=min(candidates);used.add(idx);m=Chem.Mol(mols[idx])
  if not any(a.GetAtomicNum()==1 for a in m.GetAtoms()):raise ValueError(f'Pose {idx} has no explicit H')
  m.SetProp('_Name',f"{label['pockets'][0]}_{j:03d}");m.SetIntProp('mol_id',j);m.SetIntProp('source_index',idx)
  m.SetDoubleProp('paff',float(lig['act']));m.SetDoubleProp('source_r_exp_dg',float(m.GetProp('r_exp_dg')));m.SetProp('smiles',lig['smi'])
  out.append(m);audit.append(dict(mol_id=j,source_index=idx,n_candidates=len(candidates),topology_match=bool(exact),paff=float(lig['act']),source_r_exp_dg=float(m.GetProp('r_exp_dg')),formal_charge=Chem.GetFormalCharge(m),n_atoms=m.GetNumAtoms()))
 return out,audit

def make_reference(source,poses,out,cutoff=6.0):
 from openmm import app,unit
 from scipy.spatial import cKDTree
 src=app.PDBFile(str(source));atoms=list(src.topology.atoms());xyz=np.asarray(src.positions.value_in_unit(unit.angstrom))
 lig=np.concatenate([m.GetConformer().GetPositions()[[a.GetAtomicNum()>1 for a in m.GetAtoms()]] for m in poses])
 tree=cKDTree(lig);selected=set()
 for a in atoms:
  if a.residue.name in AA and a.element.atomic_number!=1 and tree.query(xyz[a.index])[0]<=cutoff:selected.add(a.residue)
 if not selected:raise ValueError('No protein residues within pocket cutoff')
 lines=[];serial=1
 for a in atoms:
  if a.residue not in selected or a.element.atomic_number==1:continue
  x,y,z=xyz[a.index];name=a.name[:4];el=a.element.symbol.upper()[:2]
  lines.append(f"ATOM  {serial:5d} {name:^4s} {a.residue.name:>3s} A{a.residue.index+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
  serial+=1
 out.write_text('\n'.join(lines)+'\nEND\n')
 return dict(cutoff_A=cutoff,n_residues=len(selected),n_heavy_atoms=serial-1)

def initial_contact(source,mol):
 from openmm import app,unit
 from scipy.spatial import cKDTree
 src=app.PDBFile(str(source));atoms=list(src.topology.atoms());xyz=np.asarray(src.positions.value_in_unit(unit.angstrom));p=[xyz[a.index] for a in atoms if a.residue.name in AA and a.element.atomic_number!=1]
 z=np.array([a.GetAtomicNum() for a in mol.GetAtoms()]);return float(cKDTree(p).query(mol.GetConformer().GetPositions()[z>1])[0].min())

def initialize(root):
 from rdkit import Chem
 if root.exists():raise FileExistsError(root)
 root.mkdir(parents=True);(root/'logs').mkdir();assert sha(MODEL)==MODEL_SHA
 labels=json.load(open(FEP/'fep_labels.json'));tasks=[]
 for label in labels:
  name=label['pockets'][0];case=root/name;case.mkdir();src=protein_path(name);sdf=ligand_path(name)
  mols=[m for m in Chem.SDMolSupplier(str(sdf),removeHs=False,sanitize=True) if m is not None];poses,audit=match_ligands(label,mols)
  writer=Chem.SDWriter(str(case/'poses.sdf'))
  for m in poses:writer.write(m)
  writer.close();shutil.copy2(src,case/'source_protein.pdb');meta=make_reference(src,poses,case/'legacy_pocket.pdb')
  contacts=[initial_contact(src,m) for m in poses];save(case/'mapping.json',dict(source_sdf=str(sdf),source_sha256=sha(sdf),records_in_source=len(mols),selected=len(poses),mapping=audit,reference=meta,initial_min_heavy_contact_A=dict(min=min(contacts),median=float(np.median(contacts)),max=max(contacts),n_below_1A=sum(x<1 for x in contacts))))
  tasks.append(dict(index=len(tasks),system=name,uniprot=label['uniprot'],n_ligands=len(poses),source_sdf=str(sdf),source_protein=str(src)))
 save(root/'tasks.json',tasks);save(root/'protocol.json',dict(dataset=str(FEP),n_systems=len(tasks),n_ligands=sum(t['n_ligands'] for t in tasks),pose='Deposited receptor-aligned FEP SDF coordinates; no redocking',mapping='Benchmark SMILES topology plus r_exp_dg-derived pAffinity; first source record retained for duplicate deposited protomer/mapping variants',receptor='Complete residues within 6 A of any series ligand, one-residue padding, ACE/NME caps, Amber14 protonation at pH 7.4, explicit H and validated integer charge',optimization='Existing ligand all-atom FIRE; receptor fixed; fmax 0.002 eV/A; 1000 steps',score='Eint = E(complex) - E(fixed receptor) - E(ligand at bound geometry)',model=str(MODEL),model_sha256=MODEL_SHA,unsupported='Nearby metals, modified residues and cofactors excluded rather than assigned neutral charge'))
 print(json.dumps(dict(root=str(root),systems=len(tasks),ligands=sum(t['n_ligands'] for t in tasks))))

def heterogen_review(source,reference):
 from openmm import app,unit
 from scipy.spatial import cKDTree
 from prepare_aimnet_pocket import legacy_atoms,residue_key
 src=app.PDBFile(str(source));atoms=list(src.topology.atoms());xyz=np.asarray(src.positions.value_in_unit(unit.angstrom));tree=cKDTree([x[2] for x in legacy_atoms(reference)])
 nearby=[]
 for r in src.topology.residues():
  if r.name in AA or r.name in ('HOH','WAT'):continue
  aa=list(r.atoms());dist=float(tree.query(xyz[[a.index for a in aa]])[0].min());elements=sorted({a.element.symbol.upper() for a in aa})
  if dist<4:nearby.append(dict(residue=residue_key(r),resname=r.name,distance_A=dist,elements=elements))
 blockers=[];ignored=[]
 for x in nearby:
  reason=None
  if x['resname'] in MODIFIED:reason='modified residue'
  elif set(x['elements'])&METALS:reason='metal or metal-containing species'
  elif x['resname'] in COFACTORS:reason='cofactor/substrate requires dedicated treatment'
  elif x['resname'] in ('ACE','NME','NMA'):reason='source terminal cap intersects selected pocket'
  if reason:blockers.append(dict(**x,reason=reason))
  else:ignored.append(x['residue'])
 return nearby,blockers,ignored

def prepare_one(root,t):
 from prepare_aimnet_pocket import prepare
 from prepared_pocket import read_pocket
 case=taskdir(root,t);status=case/'preparation_status.json'
 if status.exists():return
 nearby,blockers,ignored=heterogen_review(case/'source_protein.pdb',case/'legacy_pocket.pdb')
 if t['system']=='bace':
  from openmm import app
  from prepare_aimnet_pocket import residue_key
  ignored += [residue_key(r) for r in app.PDBFile(str(case/'source_protein.pdb')).topology.residues() if r.name=='TLA']
 save(case/'chemistry_review.json',dict(nearby=nearby,blockers=blockers,ignored_as_nonessential=sorted(set(ignored)),manual_note='BACE TLA is deposited tartaric acid and is intentionally omitted' if t['system']=='bace' else None))
 if blockers:save(status,dict(status='excluded_chemistry',blockers=blockers));print(t['system'],'excluded');return
 try:
  meta=prepare(case/'source_protein.pdb',case/'legacy_pocket.pdb',case/'prepared.pdb',padding_residues=1,ignored_heterogens=ignored,repair_missing_sidechains=True)
  _,_,q=read_pocket(case/'prepared.pdb');save(status,dict(status='prepared',n_atoms=meta['n_atoms'],n_hydrogens=meta['n_hydrogens'],net_charge=q,modeled_sidechain_atoms=len(meta['modeled_sidechain_atoms'])))
 except Exception as e:save(status,dict(status='preparation_failed',error=repr(e)))
 print(t['system'],json.loads(status.read_text())['status'])

def score_one(root,t):
 import pandas as pd,torch
 from rdkit import Chem
 from aimnet_interaction_calculator import InteractionCalculator
 from energy import EV_TO_KCAL_MOL,gas_energy,ligand_tensors,minimize_in_pocket
 from pocket_geometry import check_contacts
 from prepared_pocket import read_pocket
 case=taskdir(root,t);status=case/'scoring_status.json';prep=json.loads((case/'preparation_status.json').read_text()) if (case/'preparation_status.json').exists() else {}
 if prep.get('status')!='prepared':save(status,dict(status='not_prepared',preparation=prep));return
 poses=[m for m in Chem.SDMolSupplier(str(case/'poses.sdf'),removeHs=False) if m];pocket_xyz,pocket_z,pocket_charge=read_pocket(case/'prepared.pdb');calc=InteractionCalculator('aimnet2-2025',device='cuda',compile_model=False)
 pcoord=torch.as_tensor(pocket_xyz,device='cuda');pnumbers=torch.as_tensor(pocket_z,device='cuda');pe=gas_energy(calc,pcoord.unsqueeze(0),pnumbers.unsqueeze(0),torch.tensor([pocket_charge],device='cuda'))
 rows=[];writer=Chem.SDWriter(str(case/'minimized.sdf'))
 for i,mol in enumerate(poses):
  row=dict(system=t['system'],mol_id=mol.GetIntProp('mol_id'),source_index=mol.GetIntProp('source_index'),paff=mol.GetDoubleProp('paff'),formal_charge=Chem.GetFormalCharge(mol),pocket_charge=pocket_charge,error='')
  try:
   coord,numbers,charge=ligand_tensors(mol,'cuda');before=check_contacts(coord.squeeze(0).cpu().numpy(),numbers.squeeze(0).cpu().numpy(),pocket_xyz,pocket_z)
   refined,complex_ev,steps,maxforce,converged,initial_complex=minimize_in_pocket(calc,coord,numbers,charge,pcoord,pnumbers,pocket_charge,fmax=.002,max_steps=1000)
   ligand_ev=gas_energy(calc,refined,numbers,charge);initial_ligand_ev=gas_energy(calc,coord,numbers,charge);after=check_contacts(refined.squeeze(0).cpu().numpy(),numbers.squeeze(0).cpu().numpy(),pocket_xyz,pocket_z)
   row.update(aimnet_interaction_kcal_mol=(complex_ev-pe-ligand_ev)*EV_TO_KCAL_MOL,pre_minimization_interaction_kcal_mol=(initial_complex-pe-initial_ligand_ev)*EV_TO_KCAL_MOL,complex_ev=complex_ev,pocket_ev=pe,ligand_ev=ligand_ev,initial_complex_ev=initial_complex,minimization_steps=steps,final_max_force_ev_A=maxforce,converged=converged,initial_min_heavy_contact_A=before,final_min_heavy_contact_A=after)
   saved=Chem.Mol(mol)
   for j,pos in enumerate(refined.squeeze(0).cpu().numpy()):saved.GetConformer().SetAtomPosition(j,tuple(map(float,pos)))
   saved.SetDoubleProp('aimnet_interaction_kcal_mol',row['aimnet_interaction_kcal_mol']);writer.write(saved);writer.flush()
  except Exception as e:row['error']=repr(e)
  rows.append(row);pd.DataFrame(rows).to_csv(case/'eint.csv',index=False);print(t['system'],i+1,'/',len(poses),flush=True)
 writer.close();errors=sum(bool(x['error']) for x in rows);save(status,dict(status='ready' if not errors else 'completed_with_errors',n_scored=len(rows)-errors,n_errors=errors))

def aggregate(root):
 import pandas as pd
 from scipy.stats import spearmanr,pearsonr,kendalltau
 tasks=json.load(open(root/'tasks.json'));summary=[];all_frames=[];ranked=[]
 for t in tasks:
  case=taskdir(root,t);prep=json.loads((case/'preparation_status.json').read_text()) if (case/'preparation_status.json').exists() else {};score=json.loads((case/'scoring_status.json').read_text()) if (case/'scoring_status.json').exists() else {}
  row=dict(system=t['system'],n_requested=t['n_ligands'],preparation=prep.get('status'),scoring=score.get('status'),n_scored=score.get('n_scored',0),n_errors=score.get('n_errors',0),exclusion_or_error=json.dumps(prep.get('blockers',prep.get('error',''))))
  p=case/'eint.csv'
  if p.exists():
   raw=pd.read_csv(p);all_frames.append(raw);f=raw[raw.error.fillna('').eq('')].copy()
   if len(f):
    f['eint_rank']=f.aimnet_interaction_kcal_mol.rank(method='average',ascending=True)
    f['pre_eint_rank']=f.pre_minimization_interaction_kcal_mol.rank(method='average',ascending=True)
    f['experimental_rank']=f.paff.rank(method='average',ascending=False)
    f.sort_values('eint_rank').to_csv(case/'ranking.csv',index=False);ranked.append(f)
   for col,prefix in [('aimnet_interaction_kcal_mol','post'),('pre_minimization_interaction_kcal_mol','pre')]:
    if len(f)>=3 and f.paff.nunique()>1:
     row[prefix+'_spearman']=float(spearmanr(-f[col],f.paff).statistic);row[prefix+'_pearson']=float(pearsonr(-f[col],f.paff).statistic);row[prefix+'_kendall']=float(kendalltau(-f[col],f.paff).statistic)
  summary.append(row)
 out=pd.DataFrame(summary);out.to_csv(root/'system_metrics.csv',index=False)
 scores=pd.concat(all_frames,ignore_index=True) if all_frames else pd.DataFrame();scores.to_csv(root/'all_scores.csv',index=False)
 ranks=pd.concat(ranked,ignore_index=True).sort_values(['system','eint_rank']) if ranked else pd.DataFrame();ranks.to_csv(root/'ranked_ligands.csv',index=False)
 usable=out[out.post_spearman.notna()] if 'post_spearman' in out else out.iloc[0:0]
 ok=scores[scores.error.fillna('').eq('')].copy() if len(scores) else scores
 validation=dict(n_systems=len(tasks),n_prepared=int((out.preparation=='prepared').sum()),n_evaluable=len(usable),n_ligands_attempted=len(scores),n_ligands_scored=len(ok),n_errors=int(scores.error.fillna('').ne('').sum()) if len(scores) else 0,mean_post_spearman=float(usable.post_spearman.mean()) if len(usable) else None,median_post_spearman=float(usable.post_spearman.median()) if len(usable) else None,mean_pre_spearman=float(usable.pre_spearman.mean()) if len(usable) else None,n_converged=int(ok.converged.sum()) if len(ok) else 0,n_nonconverged=int((~ok.converged.astype(bool)).sum()) if len(ok) else 0)
 if len(ok):validation['max_decomposition_error_kcal_mol']=float(abs((ok.complex_ev-ok.pocket_ev-ok.ligand_ev)*23.060547830619-ok.aimnet_interaction_kcal_mol).max())
 save(root/'summary.json',validation);print(json.dumps(validation,indent=2))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['init','prepare','score','aggregate']);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--index',type=int);args=ap.parse_args();root=args.run.resolve()
 if args.stage=='init':initialize(root);return
 if args.stage=='aggregate':aggregate(root);return
 tasks=json.load(open(root/'tasks.json'));selected=tasks if args.index is None else [tasks[args.index]]
 for t in selected:(prepare_one if args.stage=='prepare' else score_one)(root,t)
if __name__=='__main__':main()
