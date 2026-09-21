# Unified L1–L4 fixed-receptor AIMNet2 Eint run

This run retains only the previous ligand all-atom relaxation in a fixed, fully protonated receptor followed by interaction-energy scoring. Fresh graph-derived ligand conformers were docked with SMINA before scoring. Composite, CPCM, ligand-strain, and hydrogen-only scoring are not part of this run.

Scored **403 ligands across 41 chemically supported systems** from the requested 96-system, 917-ligand L1–L4 manifest. Unsupported chemistry, unavailable exact source structures, unsafe preparation, and the three original empty inputs are retained in `coverage.csv` with reasons.

| Layer | Systems | Ligands | Mean Eint Spearman ρ | Mean SMINA Spearman ρ |
|---|---:|---:|---:|---:|
| L1 | 14 | 133 | -0.1161 | +0.2322 |
| L2 | 13 | 130 | +0.0759 | -0.1100 |
| L3 | 9 | 90 | +0.0633 | +0.0728 |
| L4 | 5 | 50 | +0.0647 | +0.0854 |

Eint is `Ecomplex − Ereceptor − Eligand`, with receptor atoms fixed and the ligand minimized as a whole in the receptor field. Fragment energies use the final complex geometry. Lower Eint receives the better rank. Spearman values are unweighted means of within-target correlations between negative energy and experimental pAffinity.

Validation: 0 scoring errors; 383/403 ligand minimizations converged below 0.002 eV/Å; 20 finite nonconverged results are retained and flagged. No final complex energy increased. Maximum Eint decomposition error was 2.68e-09 kcal/mol. Minimum protein–ligand heavy-atom contact was 2.388 Å.

`all_scores.csv` contains ligand-level energies, ranks, convergence, contacts, and labels. `coverage.csv` contains every requested target and exclusion reason. `layer_metrics.csv`, `summary.json`, `protocol.json`, `tasks.json`, `code_sha256.json`, and `removed_L1_code.json` record aggregate results and provenance.
