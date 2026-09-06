# A finite-type bound for the one-row weave permanent

Date: 2026-09-06. Director derivation; full independent audit passed.
This is an explicit bound on the one-row quantity appearing in the
convergence agent's soft weave inequality. It does not evaluate the number
of Boolean rows with each profile, and is not an original cap certificate.

## 1. Exact finite table formula

Let a word of length d have l distinct types with multiplicities n_a>0,
sum_a n_a=d. Let K be a nonnegative l-by-l matrix, and let M be the
d-by-d matrix formed by repeating K_ab on the corresponding type blocks.
Then, exactly,

```math
\frac{\operatorname{per}M}{d!}
=\sum_{\substack{N_{ab}\in\mathbb Z_{\ge0}\\
                  \sum_bN_{ab}=n_a,\ \sum_aN_{ab}=n_b}}
 \frac{\prod_a(n_a!)^2}{d!\prod_{a,b}N_{ab}!}
                  \prod_{a,b}K_{ab}^{N_{ab}}.             (1)
```

To count permutations with table N, first partition the labeled columns
of each type b among the input row types (n_b!/prod_a N_ab! ways), then
bijection these assigned columns to the labeled rows of each type a
(n_a! ways). Multiplication gives (1). No independence of the entries of
a permutation is assumed.

Write p_a=n_a/d, let Pi(p,p) denote couplings with these two marginals,
and use natural logarithms. Set

```math
\Psi_K(p)=\sup_{\pi\in\Pi(p,p)}
 \left\{\sum_{a,b}\pi_{ab}\log K_{ab}
                    -D(\pi\Vert p\otimes p)\right\}.     (2)
```

Zero entries are permitted, with the usual extended-value convention;
a table using a zero entry contributes zero. Elementary factorial bounds
give the quantitative UPPER bound

```math
\frac{\operatorname{per}M}{d!}
\le e^{2l}(d+1)^{l^2+2l}\exp\{d\Psi_K(p)\}.             (3)
```

Indeed n!<=[e(n+1)](n/e)^n for n>=1, while n!>=(n/e)^n
(and 0!=1). Apply these to each summand of (1). Its logarithmic main term
is exactly d times the expression in (2), at pi=N/d. There are at most
(d+1)^{l²} tables. This proves (3) without an asymptotic limit or a lower
bound on the entries of K.

For fixed l and fixed strictly positive K, the corresponding normalized
logarithm converges to (2) whenever p converges with positive coordinates.
For completeness, approximate a maximizing coupling by integer tables
with the exact integer margins and entrywise O_l(1) error. One construction
rounds down all entries and fills the remaining nonnegative integer row
and column deficits by an integral bipartite flow; there are O(l²) remaining
units. Stirling's formula and continuity of entropy then give the lower
bound matching (3). This limiting assertion is not used when l grows.

## 2. Operational use on arbitrary spectral profiles

The actual soft kernel, for magnitudes u,v>=0, is

```math
K_t(u,v)=\tfrac12\big(e^{-t(u-v)^2/k}+e^{-t(u+v)^2/k}\big).
```

It is positive semidefinite: it is the inner product kernel of the
even projection of Gaussian-kernel feature vectors. That property is
used in the soft graph contraction BEFORE applying the present bound.

For an arbitrary one-row spectrum a and any declared finite bins I_b,
replace K_t(a_i,a_j) by

```math
\overline K_{bc}=\sup_{u\in I_b,v\in I_c}K_t(u,v).
```

Entrywise domination and nonnegativity of the permanent license (3) with
this bin matrix. The bin matrix need NOT itself be positive semidefinite.
Thus the square-root permanent appearing in the weave estimate is bounded by

```math
e^l(d+1)^{(l²+2l)/2}\exp\{d\Psi_{\overline K}(p)/2\}.    (4)
```

Diagonal-coordinate removal in that estimate is handled by applying (4)
to each of the at most m possible removed coordinates, then taking its
specified maximum. This is a deterministic bound, not a typical-spectrum
Gaussian replacement. A tail bin can be included, but its supremum may
be one and can make the estimate uninformative.

## 3. The remaining counting obligation is explicit

Suppose a family of length-k Boolean spin rows has at most
exp{m c(p)+o(m)} members with binned, diagonal-deleted spectral profile p,
uniformly over a fixed bin collection and the allowed target selectors.
Then the one-row permanent partition sum is bounded at exponential scale by

```math
\frac1m\log Z_T(t)
\le\sup_p\{c(p)+\tfrac12\Psi_{\overline K_t}(p)\}+o(1),  (5)
```

since d=m-1 and the number of bin profiles is polynomial in m. A uniform
negative bound on the right, stronger than the positive tilt penalty in
the soft weave inequality, would supply a genuine cap upper bound.

The trivial count c(p)<= (k/m)log 2 is insufficient near highly atomic
profiles. Exact bent/plateaued counts address some such profiles but do
not address every low-entropy spectrum or its perturbations. Equation (5)
isolates this gap at one-row scale; it does not assume the desired count.

The finite-type permanent variational formula is standard permutation
entropy/type counting, not claimed as a new general theory. Its use here
is the exact quantitative reduction from a candidate original-signing
construction to a smaller, testable spectral-profile obligation.
