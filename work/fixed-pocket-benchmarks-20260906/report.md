# Fixed-pocket minimized benchmark results

Scored 1667 ligand cases. Complex minimization converged for 1500; isolated CPCM minimization converged for 1506.

Composite = minimized complex interaction + bound desolvation + local ligand strain. Protein atoms remain fixed; all bound terms use the same minimized ligand geometry.

## T3 affinity ranking

| Layer | Systems | Ligands | Interaction Spearman | Composite Spearman |
|---|---:|---:|---:|---:|
| L1 | 24 | 231 | -0.0540 | -0.0165 |
| L2 | 24 | 236 | 0.0242 | 0.0441 |
| L3 | 23 | 230 | 0.0488 | 0.0698 |
| L4 | 22 | 220 | 0.0525 | 0.0177 |

Correlations compare negative energy with pAffinity. Undefined correlations are excluded from means.

## CASF

24 pockets, 480 previously selected poses. Mean shortlist EF1: interaction 0.7579; composite 0.9722; retrieval 3.2024.
Same selected poses as the previous scoring protocol; no all-pose minimization or reselection. The separate 20-ligand example is included in per_system_metrics.csv outside the 24-pocket average.

## 250-ligand smoke test

interaction: AUROC 0.4543; EF1 0.0000.
composite: AUROC 0.4472; EF1 1.4368.

## Validation and coverage

All input ligand IDs matched; no scoring errors; energy arithmetic passed at 1e-6 kcal/mol.
Three zero-byte T3 inputs remain excluded. Finite nonconverged scores are retained and flagged. The smoke set is nonrandom; T3 is active-only. These are public-checkpoint reconstructed scores.

Each system includes component energies, minimized SDF poses, and interaction/composite ranking CSVs. protocol_comparison.csv.gz aligns old/new scores.

Smoke input interaction previously used aimnet2-2025; current interaction uses the same legacy gas checkpoint as the other composite terms. Smoke comparisons therefore include a checkpoint change as well as geometry minimization.

Both minimizations converged for 1383 of 1667 cases. No complex energy increased during minimization.
