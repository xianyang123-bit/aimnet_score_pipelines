# KIN66 and PLA15 interaction-energy validation

All 81 systems were evaluated as single points with the current safe-neighbor calculator. No additional minimization was performed. Model names and hashes are recorded; the main scoring pipeline was not modified.

Source: [authors' data deposit, version 2](https://doi.org/10.1184/R1/33043574.v2).

## B97-3c reference comparison

Errors are in kcal/mol. Positive bias means underbinding (interaction energy is too high).

| Model | Dataset | N | MAE | RMSE | Mean bias | Spearman |
|---|---|---:|---:|---:|---:|---:|
| Current wB97M checkpoint | KIN66 | 66 | 62.128 | 62.647 | +62.128 | 0.8056 |
| Current wB97M checkpoint | PLA15 | 15 | 29.429 | 31.037 | +29.429 | 0.9571 |
| AIMNet2(2025), member 0 | KIN66 | 66 | 2.764 | 3.494 | +2.610 | 0.9563 |
| AIMNet2(2025), member 0 | PLA15 | 15 | 4.168 | 5.397 | +3.003 | 0.9821 |
| AIMNet2(2025), 4-member mean | KIN66 | 66 | 2.754 | 3.225 | +2.676 | 0.9755 |
| AIMNet2(2025), 4-member mean | PLA15 | 15 | 4.025 | 4.770 | +2.949 | 0.9929 |

The current checkpoint systematically underbinds. Its RMSE against wB97M-D3(BJ)/def2-TZVPP is 70.835 kcal/mol on KIN66 and 45.736 on PLA15, so the discrepancy is not explained merely by comparing different DFT functionals.

## Reproduction and numerical checks

The authors' supplied AIMNet2(2025) predictions have B97-3c RMSE 3.126 on KIN66 and 4.768 on PLA15. Our four-member ensemble gives 3.225 and 4.770. Direct RMSE between our ensemble predictions and their model predictions is 0.135 and 0.041 kcal/mol, respectively. This is close reproduction, not exact identity.

The current checkpoint and AIMNet2(2025) default both use full nonperiodic Coulomb in this environment. Explicitly selecting full Coulomb repeats member 0 within 0.003 kcal/mol. All model evaluations used the verified explicit torch.cdist neighbors.

The released KIN66 fragment files differ from the complex coordinates by at most 0.000866 Angstrom after atom matching. We also evaluated fragments extracted at the exact complex coordinates. B97-3c RMSE shifts by less than 0.01 kcal/mol for current and member-0 models; PLA15 fragment coordinates match exactly. These file-precision differences do not explain the old-model errors.

KIN66 charges were zero for complex/protein/ligand as specified by the benchmark. PLA15 charges came from its PDB REMARK records; charge addition and fragment atom counts were checked. Atom matching used an element-constrained assignment with a 0.002 Angstrom maximum allowed difference.

The wB97M-V reference column contains only 20 directly converged KIN66 results. Its 50 estimated values are evaluated separately and never substituted for direct DFT values. Some systems have both entries.

## Implication for the scoring pipeline

Use a separately selected AIMNet2(2025) calculator for the interaction term. Retain a gas-phase reference compatible with the existing wB97M-CPCM model for the desolvation term; replacing that gas reference with B97-3c would mix incompatible energy baselines. This benchmark checks interaction-energy accuracy, not experimental affinity ranking.

## Files

- predictions.csv: every complex, protein, ligand, and interaction energy; charges and DFT references.
- metrics.csv: all models, geometry conventions, and reference levels, with MAE/RMSE/bias/correlations.
- geometry_validation.csv: coordinate matching diagnostics.
- model_metadata.json and provenance.json: runtime, models, hashes, sources.
- data/: released benchmark geometries and KIN66 reference tables; PLA15.csv is at the root.
- parity.png: predicted versus B97-3c interaction energies.
