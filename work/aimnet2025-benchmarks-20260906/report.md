# AIMNet2(2025) member-0 benchmark results

Scored 1667 ligand cases. Complex minimization converged for 1440; isolated CPCM minimization converged for 1490.

Composite = minimized complex interaction + bound desolvation + local ligand strain. Protein atoms remain fixed; all bound terms use the same minimized ligand geometry.

## T3 affinity ranking

| Layer | Systems | Ligands | Interaction Spearman | Composite Spearman |
|---|---:|---:|---:|---:|
| L1 | 24 | 231 | -0.0407 | -0.0877 |
| L2 | 24 | 236 | 0.0478 | 0.1018 |
| L3 | 23 | 230 | 0.2350 | 0.2174 |
| L4 | 22 | 220 | 0.1330 | 0.1169 |

Correlations compare negative energy with pAffinity. Undefined correlations are excluded from means.

## CASF

24 pockets, 480 freshly selected and minimized poses. Mean shortlist EF1: interaction 0.9722; composite 0.9722; retrieval 3.2024.
Pre-minimization interaction EF1 was 1.6389. Pose selection changed for 320 of 480 ligands.
All 47,896 available poses were rescored with AIMNet2(2025), then one pose per ligand was selected for fixed-pocket minimization. The separate 20-ligand example is included in per_system_metrics.csv outside the 24-pocket average.

## 250-ligand smoke test

pre_min_interaction: AUROC 0.5796; EF1 0.0000.
interaction: AUROC 0.5422; EF1 1.4368.
composite: AUROC 0.5216; EF1 0.0000.

## Validation and coverage

All input ligand IDs matched; no scoring errors; energy arithmetic passed at 1e-6 kcal/mol.
Three zero-byte T3 inputs remain excluded. Finite nonconverged scores are retained and flagged. The smoke set is nonrandom; T3 is active-only. These are public-checkpoint reconstructed scores.

Each system includes component energies, minimized SDF poses, and interaction/composite ranking CSVs. protocol_comparison.csv.gz aligns old/new scores.

Interaction and fixed-pocket minimization use AIMNet2(2025) member 0. Desolvation uses the separate legacy wB97M gas and wB97M-CPCM pair. CASF comparisons also include refreshed pose selection.

Both minimizations converged for 1313 of 1667 cases. 0 complexes ended above their initial energy.

## Change from previous checkpoint

| Metric | Previous | AIMNet2(2025) member 0 |
|---|---:|---:|
| T3 L1 composite mean Spearman | -0.0165 | -0.0877 |
| T3 L2 composite mean Spearman | 0.0441 | 0.1018 |
| T3 L3 composite mean Spearman | 0.0698 | 0.2174 |
| T3 L4 composite mean Spearman | 0.0177 | 0.1169 |
| CASF interaction mean shortlist EF1 | 0.7579 | 0.9722 |
| Smoke interaction AUROC | 0.4543 | 0.5422 |
| CASF composite mean shortlist EF1 | 0.9722 | 0.9722 |
| Smoke composite AUROC | 0.4472 | 0.5216 |
