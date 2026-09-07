# Extension: finite-rank surgery on actual weighted profile minimizers

Date: 2026-09-07. Status: derived from the independently audited signed
finite-rank construction; pending separate root/constructive audit.

Superseded implementation: the weighted section of
`flatify_adversary_2026_09_07_target_contraction_surgery.md` removes masking
and Schur bias entirely and has now been independently checked by the
constructive agent. The support restrictions here remain important.

Let C_ij=a_ij A_ij be a hollow weighted signing, with each nonzero amplitude
in [b,B], where 0<b<=B are fixed. Suppose Q(C)<=C0 n^(3/2).
Let G=UTU^T, U^TU=I_r, ||T||op<=1, and require G_ij=0 whenever a_ij=0,
for i!=j. This support condition is essential. Put t=theta sqrt(n).

Apply the same leverage mask Z and define G0=ZGZ and the PSD entrywise
majorant K=sum_l |lambda_l| |Zv_l||Zv_l|^T. Thus |G0_ij|<=K_ij<=mu.
For nonzero edges set

    p_ij=t/(2b) [K_ij+b A_ij G0_ij/a_ij].

At zero edges do nothing. Since b/a_ij<=1, these are nonnegative and at
most tmu/b, hence valid when tmu<=b. Their exact coefficient mean is

    Cbar=C-t offdiag(G0)-(t/b)(C circ K).

The zero-edge support condition ensures the displayed identity there too.
The audited diagonal-majorant argument applies to ANY real hollow C, not
only uniform signs, and gives Q(C circ K)<=2K_G mu Q(C). Thus the Schur
bias costs at most 2K_G theta C0 mu sqrt(n)/b after normalization.

Expected edit count is at most tnr/(2b). Centered coefficient magnitudes
are at most 2B and their total variance is at most 2tB^2 nr/b. Bernstein
and the full-cube union bound therefore give normalized rounding error

    2B sqrt(theta r a/b) n^(-3/4)+4B a/[3n^(3/2)],
    a=(n+2)log2.

Markov simultaneously ensures at most 2tnr/b edits for some realization.
The masking error is unchanged. Consequently an actual signing C' with
the SAME edge amplitudes and zero pattern satisfies

    Q(C')<=Q(C-theta sqrt(n) offdiag(G))+n^(3/2)e_n,

where

    e_n=theta sqrt(r/(mu n))+theta r/(2n)
       +2K_G theta C0 mu sqrt(n)/b
       +2B sqrt(theta r a/b)n^(-3/4)+4B a/[3n^(3/2)].

For r=o(sqrt(n)), choosing mu=(r/sqrt(n))^(1/3)/sqrt(n) makes e_n=o(1)
for fixed theta,b,B. Probabilities are eventually valid because
theta(r/sqrt(n))^(1/3)<=b. Exact or near optimality within this FIXED
weighted profile then yields the same minimax stationarity law for every
feature space whose entire signed operator ball respects the zero pattern.
The additional full-versus-offdiagonal payment remains theta r/(2n).

## Application to the reveal profile: precise support scope

In the central band R/u∈[delta/2,1-delta/2], u>=4, nonzero variances
are bounded below by min(delta/2,1/3) and above by 2/delta for large n.
The only structural missing block is between known positive and known
negative vertices. Hence arbitrary feature spaces supported on

* known-positive vertices together with unknown vertices, or
* known-negative vertices together with unknown vertices

obey the required zero-pattern condition. Their weighted cap minima are
therefore stationary under these sub-square-root-rank perturbations.
One may also use block-diagonal perturbations across the missing block.

This is an actual-profile operation, not a claim that a low-rank variance
innovation remains low rank after entrywise multiplication by its signing.
It does not control the signed reveal drift, and it does not assert
stationarity for finite-temperature pressure minimizers merely from cap
optimality. Those are separate obligations.
