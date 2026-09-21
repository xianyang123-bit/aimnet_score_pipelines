#!/usr/bin/env python3
"""Recover a pocket from its original coordinate frame and prepare a capped cluster.

Requires a residue-labelled PDB/mmCIF and the original heavy-atom pocket PDB.
Uses complete standard residues, neutral ACE/NME peptide caps, explicit H, and
an integer net charge derived from a validated Amber topology. Refuses missing
heavy atoms, unknown chemistry, ambiguous mappings, and unavailable cap geometry.
This is a conservative baseline, not a protein pKa prediction method.
"""
import argparse
import hashlib
import io
import json
import random
import re
from pathlib import Path

import numpy as np
import openmm
from openmm import app, unit
from scipy.spatial import cKDTree

AA = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def residue_key(r):
    return f'{r.chain.id}:{r.id}:{r.insertionCode}:{r.name}'


def legacy_atoms(path):
    rows = []
    for l in Path(path).read_text().splitlines():
        if l.startswith(('ATOM  ', 'HETATM')) and l[76:78].strip() != 'H':
            rows.append((l[12:16].strip(), l[76:78].strip(),
                         [float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    if not rows:
        raise ValueError('Reference pocket has no heavy atoms')
    return rows


def prepare(source_path, reference_path, output_path, ph=7.4, padding_residues=1,
            ignored_heterogens=(), variants=None, seed=20260921, repair_missing_sidechains=False):
    out = Path(output_path)
    meta_path = out.with_suffix('.prep.json')
    if out.exists() or meta_path.exists():
        raise FileExistsError('Use new output paths; preparation never overwrites')
    source = (app.PDBxFile if str(source_path).lower().endswith(('.cif', '.mmcif')) else app.PDBFile)(str(source_path))
    atoms = list(source.topology.atoms())
    residues = list(source.topology.residues())
    coords = np.asarray(source.positions.value_in_unit(unit.angstrom))
    ref = legacy_atoms(reference_path)
    protein = [a for a in atoms if a.residue.name in AA and a.element.atomic_number != 1]
    if not protein:
        raise ValueError('Source has no recognized, residue-labelled protein atoms')
    tree = cKDTree(coords[[a.index for a in protein]])
    matches, max_distance = [], 0.0
    for name, element, xyz in ref:
        hits = [protein[i] for i in tree.query_ball_point(xyz, .02)
                if protein[i].name == name and protein[i].element.symbol.upper() == element.upper()]
        if len(hits) != 1:
            raise ValueError(f'Original coordinate mapping failed for {name} at {xyz}: {len(hits)} matches. '
                             'Use the exact source structure/assembly and coordinate frame.')
        a = hits[0]
        matches.append(a)
        max_distance = max(max_distance, float(np.linalg.norm(coords[a.index] - xyz)))
    if len({a.index for a in matches}) != len(matches):
        raise ValueError('Reference maps multiple atoms to one source atom')
    selected = {a.residue for a in matches}
    initial = set(selected)
    by_residue = {r: {a.name: a for a in r.atoms() if a.element.atomic_number != 1} for r in residues}
    # Only physically present peptide bonds qualify as cap geometry.
    previous, following = {}, {}
    for chain in source.topology.chains():
        rr = list(chain.residues())
        for r, s in zip(rr, rr[1:]):
            if r.name not in AA or s.name not in AA:
                continue
            c, n = by_residue[r].get('C'), by_residue[s].get('N')
            if c and n and 1.1 < np.linalg.norm(coords[c.index] - coords[n.index]) < 1.8:
                following[r], previous[s] = s, r
    for _ in range(padding_residues):
        selected |= {d[r] for r in list(selected) for d in (previous, following) if r in d}
    # Preserve disulfide partners rather than silently converting disulfides to thiols.
    source.topology.createDisulfideBonds(source.positions)
    disulfides = [(a, b) for a, b in source.topology.bonds() if a.name == b.name == 'SG']
    changed = True
    while changed:
        before = len(selected)
        for a, b in disulfides:
            if a.residue in selected or b.residue in selected:
                selected.update((a.residue, b.residue))
        # Fill short internal gaps so two fragments cannot borrow overlapping caps.
        for r in list(selected):
            gap, nxt = [], following.get(r)
            while nxt is not None and nxt not in selected and len(gap) < 3:
                gap.append(nxt)
                nxt = following.get(nxt)
            if nxt in selected and len(gap) <= 2:
                selected.update(gap)
        changed = len(selected) != before
    # Require a deliberate choice about nearby ligands/cofactors/ions.
    excluded = []
    selection_xyz = coords[[a.index for r in selected for a in r.atoms() if a.element.atomic_number != 1]]
    selection_tree = cKDTree(selection_xyz)
    for r in residues:
        if r.name in AA or r.name in ('HOH', 'WAT'):
            continue
        aa = list(r.atoms())
        distance = float(selection_tree.query(coords[[a.index for a in aa]])[0].min())
        if distance < 4.0:
            key = residue_key(r)
            excluded.append({'residue': key, 'distance_A': distance})
            if key not in ignored_heterogens:
                raise ValueError(f'Nearby nonprotein residue {key} ({distance:.2f} A). '
                                 'Review it: essential cofactors/metals need another protocol; '
                                 'explicitly --exclude-heterogen only a ligand intentionally replaced or nonessential species.')
    topology, positions, provenance = app.Topology(), [], []
    atom_map, caps = {}, []
    def add_residue(chain, src, name=None, cap_atoms=None):
        name = name or src.name
        r = topology.addResidue(name, chain, id=src.id, insertionCode=src.insertionCode)
        specs = cap_atoms or [(a.name, a.name) for a in src.atoms() if a.element.atomic_number != 1]
        for new_name, old_name in specs:
            if old_name not in by_residue[src]:
                raise ValueError(f'Missing backbone atom for cap: {residue_key(src)} {old_name}')
            old = by_residue[src][old_name]
            new = topology.addAtom(new_name, old.element, r)
            positions.append(coords[old.index] / 10)
            if not cap_atoms:
                atom_map[old] = new
            provenance.append({'source_atom_index': old.index, 'source_residue': residue_key(src),
                               'output_residue_index': r.index, 'output_atom_name': new_name, 'cap': bool(cap_atoms)})
        return r
    starts = sorted((r for r in selected if previous.get(r) not in selected), key=lambda r: r.index)
    visited = set()
    for start in starts:
        chain = topology.addChain(str(len(list(topology.chains())) + 1))
        if start in previous:
            p = previous[start]
            add_residue(chain, p, 'ACE', [('CH3', 'CA'), ('C', 'C'), ('O', 'O')])
            caps.append({'type': 'ACE', 'source': residue_key(p)})
        else:
            # Native N terminus only. An internal missing segment must not become a charged terminus.
            rr = list(start.chain.residues())
            if start is not rr[0]:
                raise ValueError(f'Unresolved chain break before {residue_key(start)}; no safe cap geometry')
        r = start
        while True:
            add_residue(chain, r)
            visited.add(r)
            nxt = following.get(r)
            if nxt in selected:
                r = nxt
                continue
            if nxt is not None:
                add_residue(chain, nxt, 'NME', [('N', 'N'), ('C', 'CA')])
                caps.append({'type': 'NME', 'source': residue_key(nxt)})
            elif 'OXT' not in by_residue[r]:
                raise ValueError(f'No cap geometry or explicit terminal OXT after {residue_key(r)}')
            break
    if visited != selected:
        raise ValueError('Cyclic/unsupported peptide connectivity requires explicit preparation')
    topology.createStandardBonds()
    for a, b in disulfides:
        if a in atom_map and b in atom_map:
            topology.addBond(atom_map[a], atom_map[b])
    positions = np.asarray(positions) * unit.nanometer
    modeled_atoms = []
    if repair_missing_sidechains:
        from pdbfixer import PDBFixer
        # Repair only side chains in the already selected, capped cluster.
        # Never invent unresolved loops, backbone atoms or terminal chemistry.
        buffer = io.StringIO()
        app.PDBFile.writeFile(topology, positions, buffer)
        fixer = PDBFixer(pdbfile=io.StringIO(buffer.getvalue()),
                         platform=openmm.Platform.getPlatformByName('Reference'))
        # Avoid the PDB interchange format rounding high-precision mmCIF coordinates.
        fixer.topology, fixer.positions = topology, positions
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        if fixer.missingTerminals:
            raise ValueError('Missing terminal atoms require explicit source preparation')
        for r, missing in fixer.missingAtoms.items():
            names = [a.name for a in missing]
            if r.name not in AA or set(names) & {'N', 'CA', 'C', 'O', 'OXT'}:
                raise ValueError(f'Missing backbone/cap atoms cannot be repaired automatically: {r} {names}')
            origin = next(p['source_residue'] for p in provenance if p['output_residue_index'] == r.index)
            modeled_atoms.extend(dict(source_residue=origin, output_residue_index=r.index,
                                      output_atom_name=n, cap=False, modeled=True, source_atom_index=None) for n in names)
        if modeled_atoms:
            before = {(a.residue.index, a.name): np.asarray(positions[a.index].value_in_unit(unit.angstrom))
                      for a in topology.atoms()}
            fixer.addMissingAtoms(seed=seed)
            after = {(a.residue.index, a.name): np.asarray(fixer.positions[a.index].value_in_unit(unit.angstrom))
                     for a in fixer.topology.atoms()}
            if set(before) - set(after) or any(np.linalg.norm(v-after[k]) > 1e-5 for k,v in before.items()):
                missing = sorted(set(before)-set(after))
                moved = [(k, float(np.linalg.norm(v-after[k]))) for k,v in before.items()
                         if k in after and np.linalg.norm(v-after[k]) > 1e-5]
                raise ValueError(f'Side-chain repair changed or removed existing source heavy atoms: missing={missing}, moved={moved[:10]}')
            topology, positions = fixer.topology, fixer.positions
            provenance.extend(modeled_atoms)
    ff = app.ForceField('amber14/protein.ff14SB.xml')
    modeller = app.Modeller(topology, positions)
    override = variants or {}
    unknown = set(override) - {residue_key(r) for r in selected}
    if unknown:
        raise ValueError(f'Unknown protonation override residues: {sorted(unknown)}')
    requested = []
    for r in topology.residues():
        origin = next(p['source_residue'] for p in provenance if p['output_residue_index'] == r.index)
        requested.append(override.get(origin) if r.name in AA else None)
    random.seed(seed)
    np.random.seed(seed)
    try:
        actual_variants = modeller.addHydrogens(ff, pH=ph, variants=requested,
                                              platform=openmm.Platform.getPlatformByName('Reference'))
    except ValueError as exc:
        match = re.search(r'residue (\d+)', str(exc))
        origin = 'unknown'
        if match:
            idx = int(match.group(1))
            origin = next((p['source_residue'] for p in provenance if p['output_residue_index'] == idx), 'unknown')
        raise ValueError(f'Receptor topology validation failed near source residue {origin}: {exc}') from exc
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    charge = sum(nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                 for i in range(system.getNumParticles()))
    if abs(charge - round(charge)) > 1e-4:
        raise ValueError(f'Nonintegral receptor charge: {charge}')
    residue_charges = []
    for r in modeller.topology.residues():
        q = sum(nonbonded.getParticleParameters(a.index)[0].value_in_unit(unit.elementary_charge) for a in r.atoms())
        residue_charges.append({'residue_index':r.index, 'name':r.name, 'charge_sum':float(q)})
    final_atoms = list(modeller.topology.atoms())
    final_xyz = np.asarray(modeller.positions.value_in_unit(unit.angstrom))
    heavy = final_xyz[[a.index for a in final_atoms if a.element.atomic_number != 1]]
    original_heavy = positions.value_in_unit(unit.angstrom)
    max_shift = float(np.linalg.norm(heavy - original_heavy, axis=1).max())
    if max_shift > 1e-5:
        raise ValueError(f'Heavy atoms moved during protonation: {max_shift} A')
    for a, b in modeller.topology.bonds():
        d = float(np.linalg.norm(final_xyz[a.index] - final_xyz[b.index]))
        lower, upper = (.7, 1.6) if 1 in (a.element.atomic_number, b.element.atomic_number) else (1.0, 2.4)
        if not lower < d < upper:
            raise ValueError(f'Invalid bond geometry {a.name}-{b.name}: {d:.3f} A')
    if cKDTree(final_xyz).query_pairs(.6):
        raise ValueError('Severe atom overlap in prepared pocket (<0.6 A)')
    text = io.StringIO()
    app.PDBFile.writeFile(modeller.topology, modeller.positions, text, keepIds=False)
    pdb_text = text.getvalue()
    meta = dict(schema='aimnet_prepared_pocket_v1', status='validated',
                source_structure=str(Path(source_path).resolve()), source_sha256=digest(source_path),
                reference_pocket=str(Path(reference_path).resolve()), reference_sha256=digest(reference_path),
                source_match_max_distance_A=max_distance, reference_atoms=len(ref),
                ph=ph, seed=seed, forcefield='amber14/protein.ff14SB.xml', openmm_version=openmm.__version__,
                charge_method='integer sum of Amber topology charges; AIMNet still predicts its own atomic charges',
                net_charge=int(round(charge)), topology_charge_sum=charge,
                residue_charge_sums=residue_charges,
                n_atoms=len(final_atoms), n_hydrogens=sum(a.element.atomic_number == 1 for a in final_atoms),
                selected_source_residues=[residue_key(r) for r in sorted(selected, key=lambda r:r.index)],
                original_pocket_residues=len(initial), padding_residues=padding_residues, caps=caps,
                protonation_variants=actual_variants, protonation_overrides=override,
                excluded_nearby_heterogens=excluded, water_policy='waters omitted; evaluate conserved waters separately',
                heavy_atom_max_displacement_A=max_shift, atom_provenance=provenance,
                modeled_sidechain_atoms=modeled_atoms,
                sidechain_repair='PDBFixer; existing heavy atoms fixed; no loops/backbone modeling' if modeled_atoms else None,
                limitations=['Standard-residue pH heuristic, not a pKa calculation',
                             'Fixed capped cluster; omitted protein and solvent response',
                             'Inspect caps and protonation near the binding site before production'])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(pdb_text)
    meta['pdb_sha256'] = digest(out)
    meta_path.write_text(json.dumps(meta, indent=2) + '\n')
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-structure', required=True)
    ap.add_argument('--reference-pocket', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--ph', type=float, default=7.4)
    ap.add_argument('--padding-residues', type=int, default=1)
    ap.add_argument('--exclude-heterogen', action='append', default=[])
    ap.add_argument('--variants-json', help='JSON map from chain:resid:insertion:resname to OpenMM protonation variant')
    ap.add_argument('--repair-missing-sidechains', action='store_true')
    args = ap.parse_args()
    if not 0 <= args.ph <= 14 or not 0 <= args.padding_residues <= 10:
        ap.error('Require pH 0..14 and padding-residues 0..10')
    meta = prepare(args.source_structure, args.reference_pocket, args.output, args.ph,
                   args.padding_residues, args.exclude_heterogen,
                   json.loads(Path(args.variants_json).read_text()) if args.variants_json else None,
                   repair_missing_sidechains=args.repair_missing_sidechains)
    print(json.dumps({k:meta[k] for k in ('status','n_atoms','n_hydrogens','net_charge','caps')}, indent=2))


if __name__ == '__main__':
    main()
