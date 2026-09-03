# EDA 2.2 — Phase-Space Coverage Findings

## Simulated grid depth axis
Three discrete depth levels: 500, 2,000, 10,000 (nominal UMI target),
achieving a log10(mean UMIs/cell) range of [1.53, 4.14] (~34 to ~13,800)
across all non-null-control simulated runs, all simulators combined.

## Real dataset coverage
| Dataset | Technology | Mean UMIs/cell | Within simulated depth range? |
|---|---|---|---|
| pbmc68k | 10x droplet | 1,400 | Yes |
| baron | inDrop (droplet) | 5,828 | Yes |
| muraro | CEL-seq2 (full-length) | 22,272 | No — above ceiling |
| ts_lung | Smart-seq2 + 10x (mixed) | 142,233 | No — above ceiling |
| tasic2018 | Smart-seq2 (full-length) | 1,739,371 | No — above ceiling |
| segerstolpe | Smart-seq2 (full-length) | 374,906 | No — above ceiling |

Sparsity axis: all six real datasets fall within the simulated grid's
achieved sparsity range [43.81%, 99.86%] — no gap on this axis.

## Interpretation
The simulated parameter grid's depth axis was calibrated around
droplet-based (UMI-counting) sequencing depths. Full-length read-counting
protocols (CEL-seq2, Smart-seq2) produce mean per-cell counts 1-2 orders of
magnitude higher than UMI-based protocols for a fundamentally different
reason (counting all reads per transcript vs. one UMI per transcript), not
because those real datasets are unusually deep for their own technology.

## Scope implication (for manuscript limitations)
This project's simulated stress-testing grid, as currently designed,
directly covers the technical regime of droplet-based scRNA-seq (2 of 6
real validation datasets) but does not extend into the depth regime
occupied by full-length protocols (4 of 6 real validation datasets). Any
claims about simulated-to-real generalization at the highest real-dataset
depths (tasic2018, segerstolpe in particular) should be scoped as
extrapolation beyond the engineered grid, not interpolation within it.
Extending the depth axis to cover this range was not undertaken here and
would require an explicit scope decision (additional simulation runs at
higher depth_label values) rather than being addressed retroactively.
