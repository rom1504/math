# Solution-hidden benchmark: Viterbi / finite-state best-path response

Status: independently derived before literature lookup; proof audited
algebraically and accompanied by an exact finite verifier.

## Operational derivation

For a prefix of a finite hidden-state max-sum model, let

```math
a_i=\max\{\text{prefix score of a path ending at }i\}.        \tag{VI.1}
```

Every future induces a continuation value `g_i`, so the complete response is

```math
R_a(g)=\max_i(a_i+g_i).                                       \tag{VI.2}
```

Arbitrary one-step future scores realize every finite `g`.  Consequently the
coarsest exact absolute-score state is the full endpoint vector `a`.  If one
accumulated scalar is carried separately, the coarsest control state is its
projective class modulo `R 1`.

A block with score matrix `K` updates and composes by

```math
(U_Ka)_j=\max_i(a_i+K_{ij}),
\qquad
(K\star L)_{ik}=\max_j(K_{ij}+L_{jk}).                        \tag{VI.3}
```

Writing `m=max_i a_i` and `x=a-m1`, one stores the baseline `m` plus a shape
with `max x=0`.  This state and algebra were frozen before comparison with
the classical survivor-metric and max-plus transfer representations.

### Theorem VI.1 (exact response metrics and rate)

Under arbitrary future score probes,

```math
d_{abs}(a,b)=\|a-b\|_\infty,                                  \tag{VI.4}
```

while after optimal scalar calibration,

```math
d_{rel}([a],[b])={1\over2}\operatorname{osc}(a-b).             \tag{VI.5}
```

Both metrics are nonexpanded by every exact future block.  For absolute
scores in `[-B,B]^n`, the worst-case `epsilon`-response code has

```math
\log_2M_{abs}(\epsilon)=\Theta(n\log(B/\epsilon)).             \tag{VI.6}
```

For projective shapes of spread at most `D`, `0<epsilon<D/10`,

```math
\lfloor D/(5\epsilon)\rfloor^{n-1}
\le M_{rel}(\epsilon)
\le n(\lceil D/\epsilon\rceil+1)^{n-1}.                       \tag{VI.7}
```

#### Proof

The maximum is sup-norm Lipschitz, and a sufficiently negative future pins
any chosen coordinate, proving (VI.4).  Best approximation by a constant is
half the oscillation, and the same coordinate probes attain both endpoints,
proving (VI.5).  If `m<=a_i-b_i<=M`, monotonicity and additive homogeneity of
`U_K` keep every output difference in `[m,M]`, proving nonexpansiveness.
Grid covering and separated grid packing prove (VI.6).  For (VI.7), grid the
`n` faces `max x=0`; for the lower bound fix `x_n=0` and use a
`5epsilon`-grid in `[-D,0]^(n-1)`. `square`

One-time approximation therefore survives an arbitrary exact future without
growth.  Re-quantizing after each block can still accumulate fresh error.

### Theorem VI.2 (zero-temperature contraction dichotomy)

For an all-finite score matrix, the global projective Lipschitz coefficient
of `U_K` is

```math
\tau(K)=\begin{cases}
0,&K_{ij}=u_i+v_j\text{ for some }u,v,\\
1,&\text{otherwise}.
\end{cases}                                                   \tag{VI.8}
```

#### Proof

In the separable case the output is

```math
(U_Ka)_j=v_j+\max_i(a_i+u_i),
```

so its shape is independent of `a`.  Otherwise two rows have a nonconstant
difference across columns.  Set their input-score difference strictly
between two of those thresholds and suppress all other rows.  Some columns
then uniquely select the first row and others the second.  A sufficiently
small perturbation of one input row is copied to exactly the first group, so
the output oscillation change equals the input oscillation change.
Nonexpansiveness supplies the opposite inequality. `square`

Thus ordinary positivity or irreducibility gives no graded global forgetting
at zero temperature.  It is possible to have a positive, aperiodic Markov
transition whose log-score matrix has `tau=1`; Birkhoff contraction at
positive temperature tends to one as temperature tends to zero.

### Corollary VI.3 (near-rank-one memory reset)

If

```math
K_{ij}=u_i+v_j+E_{ij},\qquad |E_{ij}|\le\eta,                  \tag{VI.9}
```

then, for `gamma(a)=max_i(a_i+u_i)`,

```math
\|U_Ka-(\gamma(a)1+v)\|_\infty\le\eta.                       \tag{VI.10}
```

After this block all old relative-score coordinates may be discarded in
favor of the fixed shape `v` with future-uniform error `eta`.  The shape-code
rate drops from order `(n-1)log(D/epsilon)` to
`O((n-1)log(1+eta/epsilon))`; if `eta<=epsilon`, one shape codeword suffices.

This is a quantitative approximate-state conclusion beyond exact dynamic
programming.  Its hypothesis is a uniform log-likelihood row-mixing
condition, not ordinary stochastic mixing.

## Literature comparison

The independently predicted vector is the classical Viterbi survivor metric,
and (VI.3) is max-plus matrix multiplication.  Residual normalization agrees
with weighted determinization.  The strongest exact external match is
Merlet's rank-one max-plus memory-loss mechanism; stochastic Viterbi barriers
give a different, probabilistic path-stabilization mechanism.

Primary references:

- [Viterbi, *Error bounds for convolutional codes and an asymptotically
  optimum decoding algorithm*](https://doi.org/10.1109/TIT.1967.1054010)
- [Mohri, *Semiring frameworks and algorithms for shortest-distance
  problems*](https://cs.nyu.edu/~mohri/pub/jalc.pdf)
- [Mohri, weighted determinization and minimization](https://aclanthology.org/J97-2003/)
- [Merlet, *Memory loss property for products of random matrices in the
  max-plus algebra*](https://doi.org/10.1287/moor.1090.0434)
- [Lember--Koloydenko, *A constructive proof of the existence of Viterbi
  processes*](https://doi.org/10.1109/TIT.2010.2040897)

## Benchmark verdict

**Pass, independently predicted.**  Contextual response recovered the exact
classical state and composition before names were supplied.  More
importantly, it predicted the sharp approximate response rate and the
zero-temperature obstruction to ordinary mixing.  The benchmark reinforces
the revised law: static response geometry determines one-shot rate, while a
rank-one reset—not generic positivity—is what makes that rate reusable with
less memory.
