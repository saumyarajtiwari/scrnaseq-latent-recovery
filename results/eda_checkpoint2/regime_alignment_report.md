# EDA 2.4 — Real vs. Simulated Regime Alignment Report

Nearest-neighbor matching performed in normalized (log10 depth, achieved
sparsity %) space against all 32,814 non-null-control simulated runs (see
`regime_alignment_analysis.R`). Distances are on a 0-1 normalized scale;
values under ~0.02 indicate a close, interpolated match, while larger
values indicate the real dataset falls outside the simulated grid's range
and the "nearest" match is necessarily an extrapolated edge-of-grid point,
not a true regime mirror (consistent with EDA 2.2's finding that 4 of 6
real datasets exceed the simulated depth ceiling).

## PBMC 68k
PBMC 68k is the most closely and confidently matched real dataset in this
collection (normalized distance = 0.0097), landing squarely within the
simulated grid at high sparsity (97.28%, vs. the matched regime's 98%
sparsity label) and moderate-high depth (log10 ≈ 3.15, matched to the
`depth_label = 10000` runs). Its nearest simulated regime is consistently
`scdesign3`, sparsity 0.98, depth 10000, high dropout, high separability
— indicating that PBMC 68k mirrors the project's high-sparsity,
high-dropout stress regime closely and can be treated as direct,
interpolated evidence that this specific failure regime occurs in real
10x droplet data.

## Baron (human pancreas)
Baron is the second most closely matched dataset (distance = 0.0167),
also falling inside the simulated grid's range at moderate sparsity
(89.22%) and high depth (log10 ≈ 3.77). Its nearest regime is
consistently `scdesign3`, sparsity 0.8, depth 10000, low dropout, medium
separability — a comparatively benign, well-resolved regime relative to
PBMC 68k's, consistent with Baron's fuller cell-type labeling (14 named
types, 0 QC-flagged cells) and its role in this project as a
"clean" droplet-based validation case rather than a stress case.

## Muraro (human pancreas)
Muraro falls outside the simulated grid's depth range (log10 ≈ 4.35,
above the grid's ceiling of 4.14) and its nearest match (distance =
0.096) is an extrapolated edge point at the grid's maximum depth label
(10000), moderate sparsity (0.7 label vs. Muraro's actual 72.01%),
`splatter` simulator, no dropout, low separability. This should be read
as evidence that Muraro's technical regime (CEL-seq2, higher per-cell
depth than any droplet-based simulation condition) sits beyond what the
current grid was designed to stress-test, rather than as a confirmed
regime match — a real coverage gap, not a validated alignment.

## Segerstolpe (human pancreas)
Segerstolpe shows the second-largest mismatch in this collection
(distance = 0.33), reflecting its extreme mean depth (log10 ≈ 5.57, over
an order of magnitude above the grid's ceiling) driven by Smart-seq2's
full-length read counting. Its nearest simulated points cluster at the
grid's maximum depth label with `splatter`, no dropout, varying
separability — none of which meaningfully represent Segerstolpe's actual
technical regime. This dataset should be explicitly framed in the
manuscript as falling outside the current simulated stress-testing
envelope rather than as validating any specific simulated failure mode.

## Tabula Sapiens — Lung
Tabula Sapiens Lung falls outside the simulated depth range (log10 ≈
5.15) with a moderately large mismatch (distance = 0.24), nearest-matched
to `scdesign3`, sparsity 0.7, depth 10000, no dropout, medium
separability — again an edge-of-grid extrapolation rather than a genuine
interpolated match. Given this dataset's real internal complexity (34
cell types, mixed Smart-seq2/10x sequencing technology, 4 donors), its
poor phase-space alignment likely reflects both its high effective depth
(driven substantially by its Smart-seq2 fraction) and its broader,
sparser cell-type distribution — a combination the current grid does not
directly represent.

## Tasic et al. 2018 (mouse cortex)
Tasic 2018 shows the largest mismatch of all six real datasets (distance
= 0.47), driven by its extreme mean depth (log10 ≈ 6.24, roughly two
orders of magnitude above the grid's ceiling) — a direct consequence of
Smart-seq2 full-length sequencing combined with this study's typically
long reads per cell. Its nearest simulated points are `splatter`, depth
10000, no dropout, high separability, but this pairing should not be
interpreted as a real alignment; Tasic 2018 sits further outside the
engineered stress-testing envelope than any other real dataset in this
project.

## Summary for manuscript framing
Two of six real datasets (PBMC 68k, Baron) provide genuine, interpolated
evidence that the project's simulated stress regimes — specifically
high-sparsity/high-dropout and moderate-sparsity/low-dropout droplet-based
conditions — occur in real biological data. The remaining four (Muraro,
Segerstolpe, Tabula Sapiens Lung, Tasic 2018) fall outside the simulated
grid's depth coverage entirely and cannot serve as regime-alignment
evidence in their current form; their inclusion validates the *sparsity*
axis and provides real-world biological/technical diversity (batch
complexity, tissue variety, rare cell types), but not the *depth* axis of
the engineered phase space. This asymmetry should be stated explicitly
in the manuscript, not implied as uniform validation across all six
datasets.
