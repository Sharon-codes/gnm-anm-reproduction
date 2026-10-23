# Published-result coverage and equivalence ledger

The articles contain six GNM figures and ten ANM figures. Every computational
analysis underlying those figures now has executable code and saved outputs.
This is analysis coverage, not a declaration of exact historical equivalence.
Paths below are relative to the repository root; files without repeated parent
paths are in `results/extended/`.

| Paper item | Executed evidence | Comparison and limit |
|---|---|---|
| GNM Table 1 | `results/gnm_benchmark.csv`, `results/data_audit.csv` | All 12 structures; ten counts agree. 1THO has 109 CA atoms versus 108 printed; 1ATN A+D has 629 versus 633. Printed g* constants do not match literal equations 3–5 calibration. |
| GNM Figures 1–3 | `results/gnm_profiles.png`, twelve profile CSVs | All prescribed profiles calculated. Pearson r/RMSE against crystallography recorded. Original full numerical curves unavailable; absolute normalization unresolved. |
| GNM Figure 4 | `gnm_figure4_correlations.png`, `5pti_correlation_pairs.csv`, `gnm_figure4_NMA_region_comparison.csv` | Original contour levels. Both marked positive centers and 3/4 marked negative centers lie inside matching outer contours. Missing terminal anticorrelation agrees with the reported limitation. Earlier all-atom NMA covariance unavailable. |
| GNM Figures 5–6 | `gnm_figures5_6_connectivity.png`, `gnm_chain_components.csv` | Both 1LTS C and 1OMF decomposed. Free path reproduces the U-shaped chain trend; literal Eq.8 gives the opposite trend; ring wording gives a flat trend. Independent area scaling and shared scale both exposed. Historical choice unconfirmed. |
| GNM historical timing | `results/gnm_benchmark.csv` solve timings | Current solve times measured. Original SGI Power Indigo hardware/code unavailable for a hardware-matched speed comparison. |
| ANM Figure 1 | `anm_structural_panels.png` | Actual 1AQB CA backbone and every 6 A contact drawn. Same construction; different rendering. |
| ANM Figure 2 | `anm_figure2_displacements.png`, `anm_figure2_contact_displacement.csv` | Equilibrium/displaced endpoints of an actual 13 A spring. Source is a geometric definition, not a numerical benchmark. |
| ANM Figure 3 | `anm_figure3_force_constants.png`, `anm_figure3_digitized_comparison.csv`, `force_constant_cutoff_sweep.csv` | All 12 visible spring-constant points digitized. Five GNM points within 1.6% at assumed 300 K. Seven ANM points within 2.7% using a separately labeled threefold scale alternative and 1e-6 zero threshold. Standard mechanical ANM scale fails. |
| ANM Figure 4 | `anm_figure4_spectra.png`, `anm_all_frequencies.csv`, `anm_spectra_summary.csv` | Full 519-mode spectra at body-text cutoffs 10/13/20 A; literal caption's 70 A also computed. Original bin boundaries, mass/time conversion and full histogram arrays unavailable. Shape comparison only; no Hz claim. |
| ANM Figure 5 | `results/anm_published_comparison.png`, `results/anm_published_figure5_landmarks.csv` | 13 A: experimental r=0.6578. Six vetted theory peak positions within one residue; five amplitudes differ by >1 A². One incorrectly selected experimental peak removed. Published 75 A² axis limit is clipping, not a measured terminal maximum. |
| ANM Figure 6 | `anm_structural_panels.png`, `principal_frame.json` | Computed B-factor colors, centered SVD axes, actual retinol atom sticks. CA representation instead of ribbon; author frame/origin/sign convention unreported. |
| ANM Figure 7 | `anm_figure7_anisotropy.png`, `anm_figure7_directional_variances.csv`, `1aqb_anm_covariance.npz` | All directional variances, isotropy departures, RMS displacements and covariance eigenvalue ratios. Trace/sum identities tested. Author numerical curves unavailable; historical curve equivalence unverified. |
| ANM Figure 8 | `anm_figure8_global_mode.png`, `anm_first12_mode_vectors.csv` | Second internal mode in documented principal axes. Maximum norm at residue 125 versus emphasized 94–99 loop. Axis/sign choices cannot fix this norm disagreement. |
| ANM Figure 9 | `anm_structural_panels.png`, bundled loop/retinol atom tables | Both views, four reported loop regions, two loops' heavy-atom sticks and retinol. Structural locations recreated; proposed ligand-gating mechanism not independently validated. |
| ANM Figure 10a | `anm_figure10a_conformations.png`, `anm_figure10a_conformations.csv`, `anm_mode_amplitude_comparison.csv` | Equilibrium and ± one calibrated modal RMS. Residue-95/residue-65 ratio 0.8207 versus published 8.7/6.8≈1.2794. One-RMS and two-RMS spans exposed; no amplitude fit. |
| ANM Figure 10b | `anm_figure10b_experimental_complex.png`, `1qab_RBP_TTR_heavy_atom_interface.csv` | Cited human RBP/TTR complex 1QAB, all six chains. Explicit 5 A heavy-atom interface summary is an added descriptive measurement. Pig/human numbering not asserted homologous. |
| ANM discussion of external MD/NMA/essential-dynamics studies | Mode outputs; methodological review | These are separate studies with different trajectories/force fields. Their original arrays are not available here; no independent MD validation claimed. |

## Comparison protocols

Primary networks, cutoffs and input policy were retained from the first
benchmark. Extra cutoff/normalization runs are diagnostics. Predictions were
never fitted to digitized peaks or displacement ratios.

Spring-constant points use the 846×620 Figure 3 raster. Logarithmic y ticks:
gamma=1 at pixel 374, gamma=0.1 at 526; gamma=10^((374-y)/152). Each pixel is
in the CSV. Three-pixel uncertainty corresponds to about 4.65% relative error,
not a statistical confidence interval. Temperature is not specified in the
figure's text: 300 K is a fixed assumption, not an estimated author temperature.

Physical prefactors are 3*kBT/gamma for GNM total fluctuations and kBT/gamma
for ANM Cartesian covariance. Multiplying inferred ANM gamma by three closely
matches the source plot but does not establish the original convention. It does
not repair profile/mode-shape differences. The source's 7 A ANM value also
requires discarding an extra near-zero mode (1e-6 instead of 1e-8).

GNM Figure 4 X centers use its 1003×936 raster: (31,20) is residue (0,0),
(895,825) is (58,58). One-residue uncertainty is assumed. Compare the nearest
resolved pair with the GNM outer contour; do not equate its magnitude with the
earlier NMA value. The helix marker is near the diagonal and is only location
evidence; a separate off-diagonal helix mean is reported.

Article URLs and PDF hashes/byte counts are in `config.json` and
`data/paper_provenance.json`. Original figures are not redistributed. Sparse
manual comparisons do not establish full-curve equivalence.

## Remaining historical evidence

The literal 1997 g* constants, original coordinate/selection choices, ANM
Figure 5 amplitudes and Figure 10a mode ratio remain unreconciled. Boundary,
scale and spectrum conventions are incompletely specified. ProDy agreement
does not resolve historical choices. Author numerical arrays, exact original
inputs or executable original code would permit direct tests; none accompanies
the retrieved articles. Historical ANM source-link lookups on 23 October 2026
did not yield downloadable code. No author has been contacted or asked to
endorse these results.
