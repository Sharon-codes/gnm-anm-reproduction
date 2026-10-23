# Executed benchmark: interpretation and unresolved differences

This document distinguishes three questions: whether the formulas were
implemented correctly, whether current data support similar behavior, and
whether the exact historical published results were reproduced. Only the first
is fully supported by the present checks.

## Implementation verification

Thirteen tests passed, including five added analysis checks (analytic chain covariance, boundary conditions, equilibrium-operator Hessian, principal-axis trace conservation, and correlation/thermal calibration). The original eight cover: GNM graph/nullspace, ANM six rigid-body zeros on a rigid
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

175 resolved residues give 519 internal modes and six rigid-body zeros.
Experimental B-factor Pearson r is **0.6578**. Six manually selected theoretical
peak landmarks remain after correcting an earlier selection that belonged to
the experimental curve. Five amplitudes differ by more than the approximate
1 A² transcription uncertainty. No parameter was fitted to these landmarks.
Their source pixels and calibration are in the CSV and code.

The predicted N-terminal B factor is approximately 306 A² versus observed
76.66 A². This is a model-versus-experiment discrepancy. The published graph
clips at 75 A², so its axis limit cannot be interpreted as a published model
peak or used to assert a model-versus-paper discrepancy. The revised plot
shows both the full computed range and the original axis range.

## Additional analyses now executed

GNM Figure 4 uses all original contour levels. Both approximate marked positive
regions and three of four marked negative regions are reproduced at the outer
contour level. Terminal anticorrelation is absent, as described in the paper.
Separate off-diagonal helix and inter-strand means are positive (0.2943 and
0.08319); an entire rectangular region's mean is not a contour maximum.

GNM Figures 5–6 now include full, nonbonded and chain-only networks for 1LTS C
and 1OMF. The free Rouse path gives the stated U-shaped chain profile. The
literal tridiagonal Eq.8 with endpoint diagonals 2 gives a central maximum;
the stated circle-graph interpretation gives a flat profile. Shared scale and
independent area calibration are both saved. This exposes the boundary and
scaling choices rather than silently changing the formula to match a picture.

All twelve visible ANM Figure 3 spring-constant points were digitized. Assuming
300 K, five GNM points agree within 1.6%. The standard ANM mechanical scaling
is approximately three times smaller than seven source points. A separately
labeled threefold alternative agrees within 2.7%, within the approximate 4.65%
digitization uncertainty. The 7 A point also requires the 1e-6 zero threshold.
This is evidence of a scale-convention difference, not proof of an author error
or proof of the original temperature/implementation choice.

Full spectra at 10,13,20 A and the caption's literal 70 A are saved. Units follow
sqrt(gamma*lambda) with an explicitly assumed 300 K; mass/time conversion is
unspecified, so they are not frequencies in Hz. Binning is documented as 25
bins rather than claimed identical to unrecovered historical bin boundaries.

Principal-axis variances, isotropy departures and RMS displacements reproduce
the Figure 7 computations. Absolute anisotropy departure and RMS displacement
correlate at 0.8593; contact count and RMS displacement correlate at -0.7649.
All trace/isotropy identities pass. The full covariance and frame are saved.

The first internal mode is tail-dominated. In the second mode, residue 125 has
the greatest norm. Its residue-95/residue-65 amplitude ratio is **0.8207**,
versus approximately **1.2794** (8.7/6.8) in Figure 10a. This ratio is invariant
to uniform scaling and coordinate rotations; fitting a global spring constant
cannot fix it. One-modal-RMS values are 0.5128 A at residue 65 and 0.4208 A at
95 under the same-structure thermal calibration. Both one- and two-RMS spans
are reported, without fitting to the source's absolute motion labels.

Actual loop/retinol atom coordinates support the two structural views. Figure
10b context uses the cited human RBP–transthyretin structure 1QAB. The additional
5 A heavy-atom interface summary is our explicit descriptive criterion, not an
original paper score or a binding-affinity prediction. Pig/human numbering is
not assumed equivalent.

The 12–15 A cutoff profiles are similar but not indistinguishable here:
pairwise profile correlations range from 0.9208 to 0.9804, and relative L2
differences from 0.1456 to 0.3730. These diagnostics are not replacements for
the prescribed 13 A primary result.

## Interpretation

All figure-related computational analyses now execute, with supporting tables,
arrays and images. The exact historical 1997 constants, ANM profile amplitudes
and second-mode ratio remain unreconciled. Historical original inputs/code or
author output arrays are needed to decide between preprocessing and convention
differences. The [coverage ledger](figure_coverage.md) records each limitation.
The evidence supports the implementation and several published qualitative or
convention-dependent numerical results; it does not support a declaration of
complete exact historical equivalence or a new contribution.
