# Independent audit: total-scale Walsh flux ceiling

Audit target: [`walsh_total_scale_flux_ceiling.md`](walsh_total_scale_flux_ceiling.md)
and
[`../experiments/verify_walsh_total_scale_flux_ceiling.py`](../experiments/verify_walsh_total_scale_flux_ceiling.py).

Verdict: **PASS**, after the scope and slot-count repairs now present in
Section 5.  The main diameter and packing theorems are rigorous.  One minor
wording clarification is recommended in Section 6, but it does not affect a
displayed theorem.

## 1. State-local diameter

For one block, symmetry and orthogonality of `C_a` give

```math
\left|{\sqrt n\over2}x^TC_ax\right|
\le {\sqrt n\over2}\|x\|_2^2
={1\over2}n^{3/2}.
```

Consequently changing the label in one slot changes its bounded onsite term
pointwise by at most `D n^(3/2)`.  The public term cancels before either
landscape is optimized, so neither its sign, size, density, nor arity matters.
The elementary inequality

```math
|\max f-\max g|\le\|f-g\|_\infty
```

then proves (TC.6).  The analogous assertions for minima and absolute maxima
are also one-Lipschitz.  Thus the factor `D/sqrt(k)` at total order `N=kn` is
correct.  In particular an arbitrarily large signed public connector cannot
amplify hidden onsite information.

The formulation also permits a query to use a public derived label: one must
simply count every actual derived or repeated onsite occurrence among the
`k` slots.  It does not permit the coefficient or a cross-block atom itself
to depend on the hidden state; that is exactly the advertised boundary.

## 2. Disjoint triangle cube

Flipping one flux coordinate keeps `u_i` fixed and changes precisely the
other two labels.  Hence (TC.10) has factor `2D`.  With `N=3hn`, the squared
normalized diameter and neighbour scales are respectively

```math
\left({2Dh n^{3/2}\over(3hn)^{3/2}}\right)^2
={4D^2\over27h},
\qquad
\left({2D n^{3/2}\over(3hn)^{3/2}}\right)^2
={4D^2\over27h^3}.
```

Thus coding can raise the possible normalized minimum distance only from
order `h^(-3/2)` to order `h^(-1/2)`; it cannot give a fixed total-scale gap.
This conclusion is independent of the public connector graph.

## 3. Unequal-cell packing theorem

After sorting the positive cell sizes, `i s_i<=N`, and therefore

```math
\sum_{i>r}s_i^{3/2}
\le N^{3/2}\sum_{i>r}i^{-3/2}
\le {2N^{3/2}\over\sqrt r}.
```

If a code has more than `2^r` words, two agree on its first `r` coordinates.
The hypothesis (TC.13) bounds their response distance by the displayed tail.
Taking `r=ceil(4L^2/epsilon^2)` proves (TC.14) under the draft's explicit
strict-packing convention.  The `L=0` case is indeed immediate.  If one
wants the more common convention `d>=epsilon N^(3/2)`, either note the strict
finite integral comparison when a tail is nonempty or replace `r` by
`floor(4L^2/epsilon^2)+1`.

This argument remains valid with arbitrary public interactions because
(TC.13) is a uniform pointwise assumption, not an independence assertion.

## 4. Shared-label count and the repaired scope

For `r` source labels there are at most

```math
2^{r(r+1)/2}
```

binary symmetric Gram matrices.  A relation kernel is a subspace of
`F_2^r`; padding a basis to `r` ordered vectors gives the valid crude count
`2^(r^2)`.  Hence (TC.16) is an upper bound, whether or not every formal pair
is realizable.  From `2^h<=2^((3r^2+r)/2)` one obtains exactly

```math
r\ge {\sqrt{1+24h}-1\over6}.
```

The original draft incorrectly combined this orbit count with the arbitrary
public baselines allowed in Section 1.  The repair is mathematically
necessary.  To see the obstruction, take one source label and two distinct
nonzero labels `a,b` having the same self-pairing.  They have the same
singleton `(G,R)`, while `C_a!=C_b`.  The two matrices have identical
diagonals, so some Boolean vector `z` has `H_a(z)!=H_b(z)`.  A public
coordinate-pinning penalty that is zero at `z` and sufficiently negative
elsewhere forces `z` to optimize and separates the two states.  Thus
arbitrary coordinate fields destroy the Witt-orbit quotient.

The repaired Corollary TC.3 now explicitly uses only Theorem 21.18's
coordinate-equivariant unrooted weighted Walsh-graph language and excludes
external poles, pins, and fields.  In that language `(G,R)` is sufficient,
including for synchronously derived labels.  It also correctly distinguishes
the `r` source labels from the `t` actual child slots.  Theorem TC.1 gives

```math
{\operatorname {diam}R\over(tn)^{3/2}}\le {D\over\sqrt t},
```

and the stated full-port hypothesis `t>=r` yields (TC.18).  Without that
slot hypothesis, a family could probe a large abstract tuple through only a
bounded number of derived ports, and the `h^(-1/4)` conclusion would not
follow.  The repaired text states this limitation explicitly.

For dimensional clarity, the phrase "total-scale response diameter" in
(TC.18) should be read as the diameter *divided by* `(tn)^(3/2)`.  Writing
that quotient explicitly would be a useful cosmetic edit.

## 5. Cross-block incidence escape

The added quantitative boundary (TC.19)--(TC.20) is correct.  If a flipped
bit replaces at most `d_i` atoms and each old or new atom is bounded by
`B n^(3/2)`, its pointwise effect is at most `2B d_i n^(3/2)`.  Requiring an
`epsilon(kn)^(3/2)` neighbour gap forces

```math
d_i\ge {\epsilon\over2B}k^{3/2}.
```

This rules out bounded-degree state-dependent bridge architectures but
deliberately does not rule out a nonlocal broadcast.  It is an incidence
ceiling, not a claim that every broadcast is a legitimate compressed state.

One wording clarification remains in item 4 of Section 6.  Scalar
coordinate-pinning baselines are allowed by TC.1 (and still cannot evade its
diameter bound), although they are excluded from TC.3's orbit count.  What is
outside TC.1 is a differently normalized rooted/vector output or a query
that omits most ports and normalizes only by its local order.  Stating this
distinction would avoid reading "precisely the query model excluded here" as
excluding the pinning baselines that Section 1 expressly allows.

## 6. Verifier

I reran

```text
./.venv/bin/python \
  extremal_information/experiments/verify_walsh_total_scale_flux_ceiling.py
```

and obtained

```text
Walsh total-scale flux ceiling checks passed: 101224
```

The repaired omitted-tail term `2/sqrt(199999)` is an upper bound for the
uncomputed sum beginning at `200000`; the previous `2/sqrt(200000)` was a
lower bound and therefore did not audit the desired direction.  The script
now accurately describes its mixture of exact arithmetic and numerical smoke
checks.  The Gaussian-binomial counts and equal-cell ratios are exact.  The
floating quadratic inversion is only a smoke check, but the displayed
quadratic solution is independently exact algebra and needs no computational
certification.

No remaining loophole was found involving signed or arbitrarily large public
connectors, synchronously derived labels, unequal cells, or total-scale
normalization, subject to the repaired language and full-port hypotheses of
TC.3.
