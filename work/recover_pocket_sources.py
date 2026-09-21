"""Look for exact source-coordinate matches in public structures for one target.

This does not align or replace the receptor. A candidate must contain every
reference heavy atom with the same name, element and coordinates (0.02 A).
"""
import argparse
import concurrent.futures
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
from openmm import app, unit
from scipy.spatial import cKDTree

from prepare_aimnet_pocket import AA, legacy_atoms, digest


def get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'AIMNet-pocket-source-audit/1.0'})
    with urllib.request.urlopen(req, timeout=40) as f:
        return f.read()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--uniprot', required=True)
    ap.add_argument('--layer', required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--archive', required=True)
    ap.add_argument('--limit', type=int, default=120)
    args = ap.parse_args()
    root = Path(args.root) / args.layer / args.uniprot
    root.mkdir(parents=True, exist_ok=True)
    ref = root / 'legacy_pocket.pdb'
    with tarfile.open(args.archive) as tar:
        matches = [m for m in tar.getmembers() if m.isfile() and
                   '/'.join(Path(m.name).parts[-3:]) == f'{args.layer}/{args.uniprot}/pocket.pdb']
        if len(matches) != 1:
            raise ValueError(f'Expected exactly one archived pocket: {len(matches)}')
        ref.write_bytes(tar.extractfile(matches[0]).read())
    data = json.loads(get(f'https://rest.uniprot.org/uniprotkb/{args.uniprot}.json'))
    ids = sorted({x['id'] for x in data.get('uniProtKBCrossReferences', []) if x['database'] == 'PDB'})
    (root/'uniprot.json').write_text(json.dumps(data))
    reference = legacy_atoms(ref)
    def check(pdb_id):
        path = root / f'{pdb_id}.cif'
        try:
            if not path.exists():
                path.write_bytes(get(f'https://files.rcsb.org/download/{pdb_id}.cif'))
            structure = app.PDBxFile(str(path))
            atoms = [a for a in structure.topology.atoms() if a.residue.name in AA and a.element.atomic_number != 1]
            xyz = np.asarray(structure.positions.value_in_unit(unit.angstrom))
            tree = cKDTree(xyz[[a.index for a in atoms]])
            count = 0
            for name, symbol, coordinate in reference:
                hits = [atoms[i] for i in tree.query_ball_point(coordinate, .02)
                        if atoms[i].name == name and atoms[i].element.symbol.upper() == symbol.upper()]
                count += len(hits) == 1
            result = dict(pdb_id=pdb_id, matched_atoms=count, reference_atoms=len(reference),
                          exact_match=count == len(reference), source=str(path), source_sha256=digest(path))
            if result['exact_match']:
                print('EXACT_MATCH', pdb_id, flush=True)
            return result
        except Exception as e:
            return {'pdb_id':pdb_id, 'error':repr(e)}
    print('Candidates', len(ids), 'checking', min(len(ids), args.limit), flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(check, ids[:args.limit]))
    report = dict(uniprot=args.uniprot, layer=args.layer, candidates=len(ids), checked=len(results),
                  reference_sha256=digest(ref), exact_matches=[r for r in results if r.get('exact_match')],
                  results=results)
    (root/'source_recovery.json').write_text(json.dumps(report, indent=2))
    print(json.dumps({k:v for k,v in report.items() if k != 'results'}, indent=2))


if __name__ == '__main__':
    main()
