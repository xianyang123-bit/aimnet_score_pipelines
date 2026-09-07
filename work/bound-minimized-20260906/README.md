# Fixed-pocket minimized composite
The updated aimnet2_composite_smoke.py minimizes ligand coordinates with the gas AIMNet2 model while keeping all pocket atoms fixed.
Eint = Ecomplex(bound) - Epocket(fixed) - Eligand(bound), converted from eV to kcal/mol.
Desolvation = Eligand_gas(bound) - Eligand_CPCM(bound).
Local strain = Eligand_CPCM(bound) - Eligand_CPCM(local minimum).
Composite = Eint + desolvation + local strain. The input CSV interaction is comparison metadata only.
Use --pocket pocket.pdb, --complex-max-steps 1000, --complex-fmax 0.002 and a new --output path.
The default pocket path is pocket.pdb beside the poses. --pocket-charge overrides the PDB net charge.
Legacy torsion/refinement CLI options are accepted but do not affect fixed-pocket minimization.
A minimized SDF and separate convergence flags for complex and free-ligand minimization are saved.
Existing output paths are rejected to prevent mixing scoring protocols.
Test results cover three smoke-test ligands, all converged. The previous 250-ligand, T3, and CASF ranking results have not been recomputed with this new protocol.
