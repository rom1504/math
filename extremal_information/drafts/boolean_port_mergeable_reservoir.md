# A mergeable Boolean-port response coreset

**Status.** Rigorous randomized theorem.  It upgrades the static coreset in
RC.1 for disjoint-union composition: one shared min-priority sketch is an
exact associative state, so approximation errors do not add with tree depth.
The components are classical min-hash sampling and Rademacher contraction;
their combination gives the response-specific conclusion below.

## 1. Exact random carrier

Let `U` be a finite set of labelled rows and let each `u in U` carry a port
type `s(u) in G_p`.  Independently for every replica `ell<=k`, choose a
uniform random total order `pi_ell` of `U`.  For a nonempty subset `E subset
U`, define

```math
Z_\ell(E)=\mathop{\rm argmin}_{u\in E}\pi_\ell(u),
\qquad
\mathcal R_k(E)=
 \big((\pi_\ell(Z_\ell(E)),s(Z_\ell(E)))\big)_{\ell=1}^k.      \tag{MR.1}
```

The complete state also retains the cardinality,

```math
\mathcal C_k(E)=(|E|,\mathcal R_k(E)).              \tag{MR.1a}
```

For **disjoint** `E,F`, the merge operation adds cardinalities and takes the
lower-priority entry in each replica.  The reservoir component is an
associative, commutative, idempotent semilattice, while the count is an
additive commutative monoid.  Together they give the exact identity

```math
\boxed{\mathcal C_k(E\sqcup F)
=\big(|E|+|F|,\mathcal R_k(E)\wedge\mathcal R_k(F)\big).}     \tag{MR.2}
```

Thus `C_k` is a genuine randomized congruence for disjoint-row composition,
not a sequence of freshly chosen approximants.  The minimum reservoir alone
also composes under ordinary set union, but its idempotence cannot recover
the cardinality of overlapping multisets.

From the stored row types form the estimator

```math
\widehat R_E(\epsilon)
={1\over k}\sum_{\ell=1}^k
 {|s(Z_\ell(E))\cdot\epsilon|\over p}.             \tag{MR.3}
```

The true normalized response is

```math
R_E(\epsilon)={1\over |E|}\sum_{u\in E}
 {|s(u)\cdot\epsilon|\over p}.                    \tag{MR.4}
```

The count in `C_k(E)` therefore also recovers the corresponding unnormalized
response estimate

```math
\widehat L_E(\epsilon)=p|E|\widehat R_E(\epsilon). \tag{MR.4a}
```

### Theorem MR.1 (depth-independent mergeable response approximation)

For every fixed nonempty `E subset U`,

```math
\mathbb E\|\widehat R_E-R_E\|_\infty\le {4\over\sqrt k},      \tag{MR.5}
```

and, for every `t>0`,

```math
\mathbb P\left\{
 \|\widehat R_E-R_E\|_\infty>{4\over\sqrt k}+t
\right\}
\le e^{-2kt^2}.                                      \tag{MR.6}
```

Consequently, for any prescribed family `E_1,...,E_T` of aggregates that is
chosen independently of the random orders, with probability at least
`1-delta`, simultaneously for every `a<=T`,

```math
\boxed{
\|\widehat R_{E_a}-R_{E_a}\|_\infty
\le {4\over\sqrt k}
 +\sqrt{{\log(T/\delta)\over2k}}.}                 \tag{MR.7}
```

In particular, error at a root of an arbitrarily deep union tree depends on
`k`, not on its depth.  Certifying every one of `T` predetermined nodes costs
only the logarithmic simultaneous-confidence term in (MR.7).  For example,

```math
k\ge {1\over\eta^2}
 \max\{64,2\log(T/\delta)\}                        \tag{MR.8}
```

suffices for error at most `eta` at all nodes.

#### Proof

For a fixed `E`, the minimum of a uniform random order is uniform on `E`.
The `k` orders are independent, so

```math
Z_1(E),...,Z_k(E)\quad\hbox{are iid uniform on }E.  \tag{MR.9}
```

Theorem RC.1's symmetrization--contraction proof applies verbatim to this
empirical measure and gives (MR.5).

Let `D` denote the supremum error in (MR.5), viewed as a function of its `k`
iid winners.  Replacing one winner changes every empirical query mean by at
most `1/k`, and therefore changes `D` by at most `1/k`.  McDiarmid's bounded
difference inequality gives

```math
\mathbb P\{D-\mathbb ED>t\}
\le\exp\left(-{2t^2\over k(1/k)^2}\right)
=e^{-2kt^2},                                        \tag{MR.10}
```

proving (MR.6).  Apply (MR.6) marginally to each prescribed aggregate and
take a union bound; the sketches may be dependent across aggregates.  This
proves (MR.7).  Splitting the error budget equally between its two terms
gives (MR.8). `square`

## 2. Information and scope

With the random orders treated as public randomness, each aggregate stores
its cardinality, `k` projective row types, and `k` ranks.  Its description
uses at most

```math
\lceil\log_2(|U|+1)\rceil
+k\big((p-1)+\lceil\log_2|U|\rceil\big)            \tag{MR.11}
```

bits, up to delimiters.  The rank term is the price of exact future merging;
the response payload itself remains `O(p/eta^2)` bits at fixed confidence.

The theorem covers disjoint-union feature algebra.  It does not cover an
adaptive choice of aggregates based on the sketch, arbitrary max-plus
switching, or interacting quadratic child energies.  For adaptive queries
one must certify the complete possible aggregate family (or prove a separate
adaptive generalization).  MR.1 therefore resolves the static/dynamic gap
only for this particular exact composition law.

## 3. Finite diagnostic

The companion verifier enumerates random-order tuples on a small universe.
It checks exact uniformity and independence of the winners, as well as the
reservoir merge identity and associativity on every nonempty subset.  The
additive cardinality component is immediate:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_port_mergeable_reservoir.py
```
