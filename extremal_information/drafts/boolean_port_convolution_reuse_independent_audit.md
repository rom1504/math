# Independent audit: Boolean-port convolution and reuse

**Verdict: PASS, with the occurrence semantics kept explicit.**  CR.1--CR.3
have the stated normalizations and constants.  The positive compiler applies
to an expanded tree of independently sampled occurrences, not to arbitrary
reuse of one DAG node state.

## 1. Convolution and the sharp reuse factor

With normalized probability histograms there is no missing group-size or row-
count factor:

```math
R_{\mu*\lambda}(\epsilon)
=\sum_{s,t}\mu(s)\lambda(t)K_p(st\epsilon)
=\sum_t\lambda(t)R_\mu(t\epsilon).
```

Probability averaging contracts `ell_infty`, and telescoping the factors
proves CR.4.  Every element of the projective Boolean group has order two.
Thus the `L`-fold power of
`(1-t)delta_e+t delta_a` is governed exactly by the odd-parity probability

```math
q_L(t)={1-(1-2t)^L\over2}.
```

Response is affine along this two-point segment, so the distance ratio is
exactly `q_L(t)/t`, not merely bounded by it.  Its derivative at zero is
`L`; the telescoping coefficient is locally sharp.

## 2. Occurrence trees, DAGs, and row counts

For each replica, multiplying independent samples from all leaf occurrences
has the convolution law represented by the subtree.  Replica banks are
independent, so each node marginal is an iid empirical measure and RC.1 plus
McDiarmid applies.  Node marginals may be dependent; the union bound in CR.11
does not require otherwise.

The word **occurrence** is essential.  If a semantic leaf appears twice in a
DAG expression, CR.2 requires two independently sampled occurrences after
tree expansion.  Reusing one coordinate bank gives diagonal products and is
not covered.  Likewise, the `T` in CR.11 counts the nodes of the declared
expanded occurrence tree being certified.  An adaptively chosen expansion,
or a DAG whose shared node bank is fed into both inputs, requires a union
bound over all predeclared possibilities or a different theorem.

For finite port systems the exact tensor rows are ordered pairs, including
all multiplicities.  Hence

```math
N(W\mathbin\otimes V)=N(W)N(V),
```

even when `W=V`; this is not the cardinality of a set-theoretic union or of
distinct row types.  Normalized histograms convolve, while the separately
stored product count recovers the unnormalized response `p N R`.  The draft's
row-count semantics are therefore correct.

## 3. Diagonal failure and the exact CR.14 supremum

Coordinatewise self-reuse sends every sample to `s^2=e`.  In contrast, the
uniform law on a nontrivial subgroup is convolution-idempotent, giving CR.13.
For the full uniform law, translation invariance makes its response the
constant

```math
c_p={\mathbb E|X_1+\cdots+X_p|\over p}.
```

The point response is `K_p`, lies in `[0,1]`, and attains one at the identity.
Moreover `c_2=c_3=1/2`, while Cauchy--Schwarz gives
`c_p<=1/sqrt(p)<=1/2` for `p>=4`.  Therefore every deviation is at most
`max(c_p,1-c_p)=1-c_p`, and the identity query attains it.  This proves the
exact equality

```math
d_p(u_p,\delta_e)=1-c_p,
```

as well as convergence to one.  The counterexample rules out diagonal sample
reuse, not all-pairs empirical convolution or every other self-convolution
algorithm.

## 4. Doeblin factor

Uniformly averaging a response difference gives zero because every translate
of `K_p` has the same uniform mean.  If
`lambda=alpha u_p+(1-alpha)lambda'`, its uniform part is therefore annihilated
and its remainder has operator norm at most `1-alpha` on the response metric.
This proves CR.17 with the exact advertised factor.  Iterating

```math
e_j\le\eta_j+(1-\alpha_j)e_{j-1}
```

gives CR.19 with the displayed indexing and CR.20 with geometric tail
`eta/alpha`.

## 5. Strengthened finite verification

The verifier passes after adding exact checks of the CR.14 supremum, the
response-level sharp parity ratio, and multiplication of integer tensor row
counts:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_port_convolution_reuse.py
```
