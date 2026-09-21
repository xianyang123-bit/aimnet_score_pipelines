"""Preparation, fresh docking, and fixed-receptor Eint stages."""
import json
import shutil
from pathlib import Path

import numpy as np

from run import BASE,BIN,MODEL_SHA,case_path,digest,run_command,save

METALS={'LI','NA','MG','AL','K','CA','MN','FE','CO','NI','CU','ZN','CD','HG'}
MODIFIED={'TPO','SEP','PTR','MSE','CSO','CSD','KCX','HYP','MLY'}
COFACTOR_IDS={'NAD','NAP','NDP','NAI','FAD','FMN','SAM','SAH','HEM','HEC','COH','PLP','ATP','ADP','GDP','GTP','0WD'}
COFACTOR_WORDS=('COFACTOR','NICOTINAMIDE','DINUCLEOTIDE','FLAVIN','PORPHYRIN','HEME','COENZYME')

def chemistry_decision(root,task,report):
    key=f"{task['layer']}/{task['uniprot']}";manual=json.loads((root/'chemistry_review.json').read_text()).get(key)
    if manual:
        names={x.upper() for x in manual.get('exclude_names',[])}
        keys=[x['residue'] for x in report.get('heterogen_inventory',[]) if x['residue'].rsplit(':',1)[-1].upper() in names]
        if names and not keys:
            return dict(source='reviewed_L1_decision',blocked_reason='Reviewed heterogen names were not found in the matched source inventory',**manual)
        return dict(source='reviewed_L1_decision',excluded_residue_keys=keys,**manual)
    blockers=[];excluded=[]
    for item in report.get('heterogen_inventory',[]):
        residue=item['residue'];name=(item.get('name') or '').upper();resname=residue.rsplit(':',1)[-1].upper();elements={x.upper() for x in item.get('elements',[])}
        reason=None
        if resname in MODIFIED:reason='modified amino acid requires a dedicated force-field/protonation treatment'
        elif elements & METALS:reason='nearby metal or metal-containing species requires a dedicated treatment'
        elif resname in COFACTOR_IDS or any(x in name for x in COFACTOR_WORDS):reason='nearby cofactor requires validated retention and protonation'
        if reason:blockers.append(dict(**item,reason=reason))
        else:excluded.append(residue)
    if blockers:return dict(source='conservative_automatic_review',blocked_reason='Unsupported nearby chemistry',blockers=blockers,excluded_residue_keys=[])
    return dict(source='conservative_automatic_review',excluded_residue_keys=excluded,
                rationale='Nearby nonprotein species classified as replaceable deposited ligand/additive; no metal, modified amino acid, or recognized cofactor detected.')

def prepare_target(root,task):
    from prepare_aimnet_pocket import prepare
    from prepared_pocket import read_pocket
    case=case_path(root,task);status=case/'preparation_status.json'
    if status.exists() and json.loads(status.read_text()).get('status') in ('prepared','excluded_chemistry','empty_original_cohort','source_not_recovered'):return
    if not task['n_ligands']:save(status,dict(status='empty_original_cohort'));return
    recovery=json.loads((case/'source_recovery.json').read_text()) if (case/'source_recovery.json').exists() else {}
    if not recovery.get('selected_source'):save(status,dict(status='source_not_recovered'));return
    decision=chemistry_decision(root,task,recovery);save(case/'chemistry_decision.json',decision)
    if decision.get('blocked_reason'):save(status,dict(status='excluded_chemistry',**decision));return
    output=case/'prepared.pdb'
    if output.exists():output.unlink()
    if output.with_suffix('.prep.json').exists():output.with_suffix('.prep.json').unlink()
    try:
        meta=prepare(recovery['selected_source'],case/'legacy_pocket.pdb',output,padding_residues=decision.get('padding_residues',1),
                     ignored_heterogens=decision.get('excluded_residue_keys',[]),repair_missing_sidechains=True)
        _,_,charge=read_pocket(output)
        save(status,dict(status='prepared',n_atoms=meta['n_atoms'],n_hydrogens=meta['n_hydrogens'],net_charge=charge,
                         modeled_sidechain_atoms=len(meta['modeled_sidechain_atoms']),decision_source=decision['source']))
    except Exception as e:save(status,dict(status='preparation_failed',error=repr(e),decision=decision))
    print(task['layer'],task['uniprot'],json.loads(status.read_text())['status'],flush=True)

def dock_target(root,task):
    import pandas as pd
    from rdkit import Chem
    from docking import fresh_conformer,graph,write_sdf
    from pocket_geometry import check_contacts
    from prepared_pocket import read_pocket
    case=case_path(root,task);status=case/'docking_status.json';prep=case/'preparation_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='ready':return
    if not prep.exists() or json.loads(prep.read_text()).get('status')!='prepared':save(status,dict(status='not_prepared'));return
    labels=pd.read_csv(case/'labels.csv');source={int(m.GetProp('mol_id')):m for m in Chem.SDMolSupplier(str(case/'original_poses.sdf'),removeHs=False)}
    assert set(labels.mol_id)==set(source);pocket_xyz,pocket_z,pocket_charge=read_pocket(case/'prepared.pdb')
    box=json.loads((case/'box.json').read_text());dock=case/'docking';dock.mkdir(exist_ok=True);receptor=dock/'receptor.pdbqt'
    run_command([BIN/'obabel',case/'prepared.pdb','-O',receptor,'-xr'],dock/'receptor_conversion.log')
    final=[];records=[]
    for row in labels.itertuples():
        mid=int(row.mol_id);original=source[mid];fresh,record=fresh_conformer(original)
        inp=dock/f'{mid}.input.sdf';out=dock/f'{mid}.docked.sdf';hsdf=dock/f'{mid}.with_h.sdf';write_sdf(inp,[fresh])
        cmd=[BIN/'smina','-r',receptor,'-l',inp,'-o',out,'--exhaustiveness','8','--num_modes','1','--cpu','4','--seed','1']
        for key,values in [('center',box['center']),('size',box['size'])]:
            for axis,value in zip('xyz',values):cmd.extend([f'--{key}_{axis}',str(value)])
        run_command(cmd,dock/f'{mid}.smina.log');run_command([BIN/'obabel',out,'-O',hsdf,'-h'],dock/f'{mid}.hydrogenation.log')
        mols=list(Chem.SDMolSupplier(str(hsdf),removeHs=False));assert len(mols)==1 and mols[0] is not None;mol=mols[0]
        assert graph(mol)==graph(original) and mol.GetNumAtoms()==original.GetNumAtoms() and Chem.GetFormalCharge(mol)==Chem.GetFormalCharge(original)
        for key in ('_Name','mol_id','source_index','layer','uniprot','paff','smiles'):
            if original.HasProp(key):mol.SetProp(key,original.GetProp(key))
        contact=check_contacts(mol.GetConformer().GetPositions(),[a.GetAtomicNum() for a in mol.GetAtoms()],pocket_xyz,pocket_z)
        records.append(dict(mol_id=mid,smina_affinity_kcal_mol=float(mol.GetProp('minimizedAffinity')),minimum_heavy_contact_A=contact));final.append(mol)
    write_sdf(case/'poses.sdf',final);pd.DataFrame(records).to_csv(case/'docking.csv',index=False)
    save(status,dict(status='ready',n_docked=len(final),pocket_charge=pocket_charge,all_conformers_fresh=True));print(task['layer'],task['uniprot'],'docked',len(final),flush=True)

def score_target(root,task):
    import pandas as pd
    import torch
    from rdkit import Chem
    from aimnet_interaction_calculator import InteractionCalculator
    from energy import EV_TO_KCAL_MOL,gas_energy,ligand_tensors,minimize_in_pocket
    from pocket_geometry import check_contacts
    from prepared_pocket import read_pocket
    case=case_path(root,task);status=case/'scoring_status.json';dockstatus=case/'docking_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='ready':return
    if not dockstatus.exists() or json.loads(dockstatus.read_text()).get('status')!='ready':save(status,dict(status='not_docked'));return
    assert digest(Path('/home/xianyang/.cache/aimnet/aimnet2_2025_b973c_d3_0.pt'))==MODEL_SHA
    labels=pd.read_csv(case/'labels.csv').set_index('mol_id');dock=pd.read_csv(case/'docking.csv').set_index('mol_id')
    poses=list(Chem.SDMolSupplier(str(case/'poses.sdf'),removeHs=False));pocket_xyz,pocket_z,pocket_charge=read_pocket(case/'prepared.pdb')
    device='cuda';calc=InteractionCalculator('aimnet2-2025',device=device,compile_model=False)
    pcoord=torch.as_tensor(pocket_xyz,device=device);pnumbers=torch.as_tensor(pocket_z,device=device)
    pe=gas_energy(calc,pcoord.unsqueeze(0),pnumbers.unsqueeze(0),torch.tensor([pocket_charge],device=device));rows=[];writer=Chem.SDWriter(str(case/'minimized.sdf'))
    for i,mol in enumerate(poses):
        mid=int(mol.GetProp('mol_id'));row=dict(layer=task['layer'],uniprot=task['uniprot'],mol_id=mid,paff=float(labels.loc[mid].paff),pocket_charge=pocket_charge,
                                               smina_affinity_kcal_mol=float(dock.loc[mid].smina_affinity_kcal_mol),error='')
        try:
            coord,numbers,charge=ligand_tensors(mol,device);initial_contact=check_contacts(coord.squeeze(0).cpu().numpy(),numbers.squeeze(0).cpu().numpy(),pocket_xyz,pocket_z)
            refined,complex_ev,steps,maxforce,converged,initial_complex=minimize_in_pocket(calc,coord,numbers,charge,pcoord,pnumbers,pocket_charge,fmax=.002,max_steps=1000)
            ligand_ev=gas_energy(calc,refined,numbers,charge);initial_ligand_ev=gas_energy(calc,coord,numbers,charge)
            final_contact=check_contacts(refined.squeeze(0).cpu().numpy(),numbers.squeeze(0).cpu().numpy(),pocket_xyz,pocket_z)
            row.update(aimnet_interaction_kcal_mol=(complex_ev-pe-ligand_ev)*EV_TO_KCAL_MOL,
                       pre_minimization_interaction_kcal_mol=(initial_complex-pe-initial_ligand_ev)*EV_TO_KCAL_MOL,
                       complex_ev=complex_ev,pocket_ev=pe,ligand_ev=ligand_ev,initial_complex_ev=initial_complex,
                       minimization_steps=steps,final_max_force_ev_A=maxforce,converged=converged,
                       initial_min_heavy_contact_A=initial_contact,final_min_heavy_contact_A=final_contact)
            saved=Chem.Mol(mol)
            for j,pos in enumerate(refined.squeeze(0).cpu().numpy()):saved.GetConformer().SetAtomPosition(j,tuple(map(float,pos)))
            saved.SetDoubleProp('aimnet_interaction_kcal_mol',row['aimnet_interaction_kcal_mol']);writer.write(saved);writer.flush()
        except Exception as e:row['error']=repr(e)
        rows.append(row);pd.DataFrame(rows).to_csv(case/'eint.csv',index=False);print(task['layer'],task['uniprot'],i+1,'/',len(poses),flush=True)
    writer.close();errors=sum(bool(x['error']) for x in rows);save(status,dict(status='ready' if not errors else 'completed_with_errors',n_scored=len(rows)-errors,n_errors=errors))
