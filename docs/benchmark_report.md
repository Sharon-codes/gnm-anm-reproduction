# Executed benchmark: interpretation and unresolved differences

This document distinguishes three questions: whether the formulas were
implemented correctly, whether current data support similar behavior, and
whether the exact historical published results were reproduced. Only the first
is fully supported by the present checks.

## Implementation verification

Eight tests passed: GNM graph/nullspace, ANM six rigid-body zeros on a rigid
test network, covariance versus an independent pseudoinverse, Hessian energy
curvature from finite differences, translation invariance, rotational
equivariance, explicit disconnected-network zero counts, and invalid-input
rejection. See `results/test_log.txt`.

Twelve GNM and six ANM runs were compared against ProDy 2.6.1 using identical
coordinates, spring constants, cutoffs, and matched eigenvalue-zero thresholds.
All 18 comparisons passed. Maximum relative covariance error was 7.74e-12.
The matrix maximum absolute errors and fluctuation errors are in
`results/prody_validation.csv`. This is a comparison to a community
implementation rather than proof that the original authors used identical
preprocessing.

**Numerical rank convention.** Our native absolute eigenvalue cutoff is 1e-8.
ProDy uses 1e-6. At 7 A on 1AQB these give 8 and 9 discarded modes,
respectively; the resulting covariance can differ dramatically. The initial
unmatched comparison is retained in
`results/prody_validation_initial_threshold_audit.csv`. Current validation
computes the independent covariance at ProDy's threshold for a fair numerical
comparison and also reports the native discrepancy. The primary ANM 13 A run
has six zero modes under both conventions and is not affected by this issue.

## Data audit

Original PDB downloads were made on 23 October 2026. URL and raw-file SHA256
are recorded in `data/provenance.json`. Fixed parsed tables have additional
SHA256 entries in `results/run_manifest.json`. Residues are not synthesized or
filled in. Only first-model ATOM C-alpha records are used, with one conformer
per residue and insertion codes retained. HETATM calcium ions named CA must
not be interpreted as protein C-alpha atoms.

| Entry | Paper description | Resolved C-alpha atoms used | Consequence |
|---|---:|---:|---|
| 1THO A | 108 | 109 | An insertion code is present; preserve it and disclose the count difference |
| 1ATN, A+D network | 633 total, 373 in evaluated A chain | 629 total, 371 in A | Original published counts do not equal current resolved protein C-alpha counts; do not pad the network |
| 1AQB A | 183 sequence residues | 175 | Residues 176-183 are explicitly missing in PDB REMARK 465; the published structure's sequence length is not its resolved coordinate count |

The other ten GNM cases match the counts printed in the 1997 table. The 1ATN
network includes both A and D chains, while its B-factor evaluation uses A
only, following the paper's multichain setup. Crystallographic symmetry mates
are not added. For 1OMF the actual PDB chain identifier is a blank space.

## GNM: all 12 original entries at 7 A

| Case | Pearson r with observed B |
|---|---:|
| 1LTS C | 0.6054 |
| 1CRN A | 0.6927 |
| 5PTI A | 0.8853 |
| 1ABA A | 0.6128 |
| 1VAA B | 0.5632 |
| 1THO A | 0.4202 |
| 5RSA A | 0.7235 |
| 2CCY A | 0.7350 |
| 3LZM A | 0.6168 |
| 1ULA A | 0.5721 |
| 1OMF blank | 0.4006 |
| 1ATN A, in A+D network | 0.7370 |
| Unweighted mean of per-structure r | 0.6304 |

These are new measurements of an existing model on the fixed inputs, not new
method results or published author scores. The primary paper does not tabulate
these r values. Each curve is rescaled to match the discrete sum of the
experimental B factors; RMSE after this rescaling is descriptive, not an
independent absolute-scale prediction. Pearson r does not depend on the fit.

**Published parameter check remains unresolved.** Table 1's g* values are
included in the CSV. Under the explicit convention Gamma=g*L and
B=(8*pi^2/3)*diag(L+)/g* from equations 3-5, our experimental-area fitted g*
does not reproduce those reported constants. Both the fitted multiplier and
its inverse are exposed, as are an alternate trapezoidal-area calculation and
RMSE with the published constant held fixed. No reciprocal substitution or
protein-specific cutoff tuning is used to manufacture agreement. The scale
notation/normalization and historical input choices need to be checked before
claiming the original absolute curves were reproduced.

## ANM: original 1AQB at 13 A

175 observed residues give 519 internal modes and six rigid-body zero modes.
Pearson r against crystallographic B factors is **0.6578**, and the
same-structure, mean-calibrated RMSE is recorded in the raw metrics table.
The predicted N-terminal B factor is approximately 306 A^2 versus the
experimental 76.66 A^2. The plot retains this large discrepancy in a full-range
panel and uses a separately labeled zoom for internal-peak comparisons.

Seven conspicuous theoretical local maxima were manually transcribed from the
original Figure 5 raster, without using them to fit the model. The embedded
image is 1494 by 850 pixels. Calibration used (x,y)=(180.5,741) for the origin,
x=470.5 for residue 40, and y=22 for B=75 A^2. The image pixel is treated in its
native coordinate system (downward y). One-residue and 1 A^2 error bars are
conservative transcription estimates, not statistical confidence intervals.
All picked pixels and calibrations are visible in the code/CSV; this small set
is not a digitization of the entire published curve and cannot establish
full-profile equivalence.

| Peak residue, approx. in paper | Approx. published B (A^2) | Reproduced peak residue | Reproduced B (A^2) |
|---|---:|---:|---:|
| 34 | 30.88 | 34 | 37.11 |
| 49 | 55.81 | 49 | 75.56 |
| 66 | 60.61 | 65 | 69.27 |
| 81 | 36.20 | 81 | 56.36 |
| 95 | 38.49 | 95 | 38.36 |
| 125 | 57.48 | 125 | 77.74 |
| 146 | 40.89 | 146 | 48.06 |

Locations agree to within one residue, but six of seven amplitudes differ
beyond the transcription estimates. Thus full quantitative Figure 5
reproduction has **not** been established. This comparison uses the prescribed
13 A run; no fit to the picked published amplitudes is made.

**Slow modes also need reconciliation.** The first internal mode is dominated
by residue 1, consistent with the described tail motion. In the current second
mode, the largest amplitude is at residue 125; residues 65,66,126,95 are also
among its top five. The paper emphasizes the 94-99 loop as the strongest
global-mode displacement. Coordinate rotations and eigenvector sign changes
cannot change which residue has the largest vector norm. The report therefore
does not assert that the biological global-mode interpretation is reproduced.

| Diagnostic cutoff (A) | Zero count, 1e-8 convention | Pearson r |
|---|---:|---:|
| 7 | 8 | 0.0510 |
| 10 | 6 | 0.6669 |
| 12 | 6 | 0.7252 |
| 13, primary | 6 | 0.6578 |
| 15 | 6 | 0.7404 |
| 20 | 6 | 0.7368 |

The 15 A diagnostic is not substituted for the paper's primary 13 A result.
The sensitivity itself warrants checking implementation/preprocessing history.
No claim of superiority or novelty follows from these correlations.

## What remains before an extension

1. Review the literal GNM scale definition against its table/plots and obtain
   clarification or historical inputs if possible. Do not infer an author error
   from this reproduction alone.
2. Reconcile the ANM Figure 5 amplitudes and second-mode ordering using a fully
   documented protocol; compare the original author implementation if available.
3. Complete the unattempted cross-correlation, anisotropy, and full slow-mode
   figure comparisons if the project is approved for that scope.
4. Freeze the accepted original baseline and define independent evaluation
   structures before testing any changes. Retain negative and discrepant results.

The present evidence supports an independently checked implementation and an
honest partial published-result reproduction. It does not support a new
contribution claim or a declaration of complete historical replication.
