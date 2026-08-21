# Active research state

Evidence cutoff: blank-slate direct campaign, ledger Section 10.145
(2026-08-21), started from commit `b5ec773`.

This is compact working context.  Use `ledger.md` and Git history only for
archive comparison or proof reconstruction.

## Exact problem

For a symmetric hollow sign matrix `A`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
\qquad
M_n=\min_A Q(A).
```

Determine whether `M_n/n^(3/2)` converges.  Convergence to any constant is a
solution.  Genuine nonconvergence requires fixed positive separation between
two infinite subsequences.

Let `N=binom(n,2)` and let

```math
\mathcal C_n^+
=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}
```

be the augmented cut code.  Under sign-to-bit identification,

```math
Q(a)=N-2d(a,\mathcal C_n^+),
\qquad M_n=N-2\rho(\mathcal C_n^+).
```

This is an antipodal covering-radius deficit, not one-sided Max-Cut.

## Rigorous frontier

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.
```

The exact recorded values for orders 3 through 14 are

```math
(3,4,4,5,9,10,12,13,17,18,20,21).
```

They are finite falsifiers, not asymptotic evidence.

## New reusable theorems from the blank-slate campaign

### Central radial moment proxy

Let `mu_A` be the uniform antipodal energy law and let `L_d(A)` be the least
central parity-lattice cap supporting a law with the same even moments through
degree `2d`.  Then

```math
Q(A)\operatorname{sech}\!\left(
 {\operatorname{arcosh}(w_A^{-1/2})\over d}
\right)
\le L_d(A)\le Q(A),
\qquad w_A=\mu_A\{\pm Q(A)\}\ge2^{1-n}.
```

For `d=alpha n`, the uniform factor tends to
`sech(log(2)/(2 alpha))`.  The LP dual consists of polynomials in `q^2`
nonnegative on a proposed central parity lattice.  This is an exact useful
certificate, but the required asymptotic limit of its minima is not a strict
reduction: its large-`alpha` oscillation vanishes exactly when the oscillation
of `M_n/n^(3/2)` does.

### Direct covering multiplicity obstruction

For `Z_r(b)=|C_n^+ intersect B(b,r)|`, the order-`K` Bonferroni sum is exactly

```math
|\{b:Z_r(b)>0\}|-(-1)^K
\sum_{Z_r(b)>0}\binom{Z_r(b)-1}K.
```

At `r=N/2-(c/2)n^(3/2)+O(1)` with `c<sqrt(log 2)`, an odd truncation capable
of proving noncoverage needs

```math
K\ge\exp((\log2-c^2)n+o(n)),
```

and an even truncation capable of proving coverage needs

```math
K\ge2^{n-(c/\log2)\sqrt n+o(\sqrt n)}.
```

Subexponential finite-replica inclusion--exclusion cannot improve the sphere
threshold.  At the necessary order the center array has linear affine rank;
full factorial moments invert to the full coset-multiplicity histogram.

### Rectangular projection and covariance correction

For sign rectangles,

```math
m\mu_k\le\min_B\|B\|_{\infty\to1}\le m\mu_k+k2^k,
\qquad \mu_k=\mathbb E|\varepsilon_1+\cdots+\varepsilon_k|.
```

If `T=B^T B` and
`eta_k=2^{-(k-2)} binom(k-2,floor((k-2)/2))`, then

```math
\|B\|_{\infty\to1}
\ge m\mu_k+{\eta_k^2\over m\mu_k}\sum_{p<q}T_{pq}^2.
```

The scalar rectangular constant is `0.3071059...`, below the current lower
frontier.  The correction also vanishes on every competitive signing: product
rounding gives `Q(A)>=|lambda|^3/(2(n-1))`, hence
`Q(A)=O(n^(3/2))` forces the normalized off-diagonal fourth-moment defect to
be `O(n^(-1/3))`.

### Bounded local stationarity is insufficient

For `n=2m`, the signing negative inside one `m`-set and positive elsewhere has
`Q=n^2/4` but no improving single-edge flip.  An explicit order-eight cap-12
signing has no improving one- or two-edge flip although a triple flip reaches
the exact optimum 10.  Exact-minimizer light/heavy cut witnesses therefore do
not yield a bounded-radius characterization.

## Strongest reusable older obstructions

1. Bounded moments, cycles, local/restriction profiles, fixed-level SOS, and
   separately paid scalar channels miss leading zero-entropy Boolean
   resonances.  Rich generic states reconstruct the full signed-coset response.
2. Action compactness preserves subsequential objectives but lacks lossless
   every-order exact-sign recovery.  Projective exchangeability and uniform
   mesoscopic sampling become iid and are quantitatively nonextremal.
3. Sign-near weighted recovery has a complete rounding theorem but its
   existence is equivalent to exact recovery and costs almost all edge bits.
4. Fixed-temperature pressure, bridge transport, posterior-state, sparse
   repair, and finite-feature composition branches are frozen under their
   archived no-go results.
5. Spectra alone fail: archived cospectral masks have different Boolean caps.
6. Arithmetic conference/Hadamard/Paley examples give selected constructions,
   not universal order behavior or nonconvergence.

## Current research status

The explicitly authorized blank-slate campaign is complete and is a
**STRIKE**.  It did not improve the rigorous interval, strictly reduce
convergence, or leave a class-A architecture.  Do not start another automatic
wave.

The best conditional direction is a nonperturbative coverage theorem for

```math
F_n(1)=\mathbb P_b\{Z_r(b)=0\},
```

which does not expand through subexponentially many raw factorial moments.
The exact convergence statement suggested by the pair transition is
`F_{n,1/2-epsilon}(1)=0` for every fixed positive `epsilon` and all large `n`.
Together with the known upper bound this is equivalent to convergence to
`1/2`, not a strict reduction.  The smaller first milestone is uniform
coverage at one fixed `c>0.336493364431...`, which would improve the lower
frontier.  This direction is independent of the frozen composition language,
but no simpler sufficient lemma has yet been found.

Resume only after explicit user authorization or a concrete imported theorem,
constructor, or counterexample that escapes the quantified obstruction.  Any
new candidate must state its exact implication to convergence/nonconvergence,
why its missing lemma contains less information than full Boolean/coset
optimization, and an immediate falsifier.
