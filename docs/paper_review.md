# Model understanding and original evaluation

## Why this pair

GNM and ANM form a coherent progression from a scalar contact-graph model to a
directional mechanical model of folded proteins. Their assumptions and matrix
operations can be inspected directly, without training a large foundation model.
They are proposals for a foundational reproduction exercise; no mentor approval
is implied. They are older foundational methods, not state-of-the-art AI models.

## Paper 1: GNM

Bahar I, Atilgan AR, Erman B (1997). *Direct evaluation of thermal fluctuations
in proteins using a single-parameter harmonic potential*. Folding & Design
2:173-181. DOI: 10.1016/S1359-0278(97)00024-2.
[Author-hosted full text](https://bahargroup.org/Faculty/bahar/publications/99.pdf).

**Input and output.** Input is the already folded structure, represented by one
C-alpha coordinate per resolved amino-acid residue. The contact graph determines
the relative mean-square displacement of each residue and the signed covariance
between residues. It does not determine the folded structure from a sequence.

**Construction.** Set A_ij=1 when two distinct C-alpha atoms are within 7 A and
zero otherwise. The Kirchhoff matrix is L=diag(A1)-A. A uniform spring strength
multiplies this matrix. Remove the translational zero eigenvalue and construct
the pseudoinverse C=sum_k u_k u_k^T/lambda_k over nonzero modes. Then C_ii gives
relative mean-square fluctuations; C_ij/sqrt(C_ii*C_jj) gives normalized
cross-correlations. A global thermal/spring factor supplies amplitude units.
The corresponding crystallographic relation is B_i=(8*pi^2/3)<|delta R_i|^2>.

In the notation of the 1997 paper's equations 3-5, Gamma=g*L and the prediction
would be (8*pi^2/3)*diag(L+)/g*. The table in our benchmark reports the fitted
g* under that explicit convention alongside the values printed in the paper.
The discrepancy remains open; it is not resolved by silently replacing a
reported parameter with its reciprocal.

**Assumptions.** The native coordinates define the equilibrium contact network.
Fluctuations are harmonic/Gaussian and isotropic, pair interactions are uniform,
and residue chemistry is largely discarded. No explicit solvent, ligand force
field, side-chain degrees of freedom, kinetic barriers, or time evolution is
modeled. Backbone neighbors enter through the same spatial contact rule.

**Original evaluation to inspect.** Table 1 specifies 12 PDB examples and
normalization parameters. Figures 1-3 compare predicted and experimental
C-alpha B-factor profiles; the authors normalize their theoretical curves to
the experimental areas. Figure 4 examines 5PTI fluctuation cross-correlations
against earlier detailed normal-mode analysis, and Figures 5-7 examine chain
connectivity. The final discussion gives timing examples on the original
workstation. Our local time measurements are not a hardware-matched replication
of those historical speed claims. The paper does not supply a per-protein
Pearson-r table to which our new r values can be equated.

**Implemented scope.** Original 12 entries, prescribed chain selections,
7 A contacts, full nonzero eigenspectrum, and per-residue B factors; all used
coordinates and B factors are bundled. Full published cross-correlation and
chain-connectivity figure reproduction has not been completed.

## Paper 2: ANM

Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin O, Bahar I (2001).
*Anisotropy of fluctuation dynamics of proteins with an elastic network model*.
Biophysical Journal 80:505-515. DOI: 10.1016/S0006-3495(01)76033-X.
[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1301252/).
[Author-hosted full text](https://bahargroup.org/Faculty/bahar/publications/143.pdf).

**Input and output.** The same kind of folded C-alpha structure is input, but
each residue has three displacement components. The output includes 3D mode
vectors, directional covariance blocks, and residue fluctuation magnitudes.

**Construction.** For contacts, V=(gamma/2)*sum_(i<j)(|r_i-r_j|-d_ij^0)^2.
At equilibrium, let n_ij=(r_i-r_j)/d_ij^0. The off-diagonal 3x3 Hessian block
is H_ij=-gamma*n_ij*n_ij^T and H_ii=-sum_(j!=i)H_ij. For a connected,
mechanically rigid 3D network there are six rigid-body zero modes. Additional
zero modes signal unconstrained internal motions. The inverse over retained
positive modes supplies covariance. Sum the three diagonal entries for each
residue to obtain total mean-square fluctuation; individual entries and block
eigenvectors describe directionality. Eigenvector signs are arbitrary.

**Assumptions and differences from GNM.** ANM retains equilibrium bond
directions and uses central springs, still with nonspecific interactions and a
local harmonic approximation. It cannot describe bond formation, folding,
unfolding, arbitrary large transitions, binding free energies, or mutation
effects as though it were a calibrated atomistic force field. A mode direction
is not an experimentally validated conformational pathway. Missing residues
and contacts can change the softest modes substantially.

**Original evaluation to inspect.** The application is pig plasma retinol
binding protein, PDB 1AQB. Figure 3 studies the spring scale versus cutoff;
Figure 4 compares vibrational spectra at different cutoffs; Figure 5 compares
the 13 A ANM B-factor profile with crystallography. Subsequent figures examine
directional anisotropy and slow collective motions, including loops 63-68 and
94-99 and regions 48-52 and 120-127. The authors describe the first mode as
tail-dominated and use the second mode to discuss the global motion.

**Implemented scope.** Primary 13 A network and full covariance; diagnostics
at 7,10,12,15,20 A; six-zero-mode check; experimental B-factor comparison;
seven approximate theoretical peak landmarks from Figure 5; first five mode
amplitude summaries. Full anisotropic Figure 7 and signed/principal-axis
Figure 8-10 reproductions and their supporting MD comparisons remain undone.

## Questions to answer without relying on generated summaries

1. Derive why graph row sums are zero, and explain why a direct matrix inverse
   is undefined. Demonstrate the rigid-body motions annihilated by each model.
2. Derive an ANM Hessian block by differentiating the central-spring energy.
   Explain why three scalar copies of a GNM matrix are not the ANM Hessian.
3. Explain why predicted B-factor Pearson correlation is unchanged when all
   spring constants are multiplied by one positive constant.
4. Explain why calibrating scale on measured B factors cannot validate the
   absolute fluctuation amplitude out of sample.
5. Distinguish crystal disorder from thermal motion, a static covariance from
   a time correlation, and a low-frequency mode from a simulated trajectory.
6. Explain the observed extra near-zero mode sensitivity at 7 A and why a
   consistent eigenvalue convention is required for a library comparison.
7. Account for each missing residue and insertion code before attributing a
   disagreement to an algorithm or proposing a modification.
8. Identify exactly which published figure/table each claimed reproduction
   uses, and which results remain only implementation validation.

The next decision should be whether this scope is suitable for mentorship and
how to reconcile the open baseline differences. Extensions should wait for a
fixed, reviewed reproduction protocol.
