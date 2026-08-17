# Independent audit: approximate residual-shell law

**Audit status: PASS, with three non-substantive repairs recommended before
canonicalization.**

This audit was performed independently of the draft's self-audit.  I checked
the stated constants and query in ARS.1, the exact and approximate
rank-one/cycle claims in ARS.2--2a, the sharp `delta/4` example, small-shell
universality, the empty-core example, and overlap with the canonical theorem
set.

## 1. Verification performed

The repository verifier ran successfully:

```text
$ ./.venv/bin/python \
    extremal_information/experiments/verify_approximate_residual_shell_law.py
approximate residual-shell checks passed: 12283
```

I also independently generated 500 random `3 by 3` cyclic products of
entrywise-perturbed rank-one max-plus matrices, enumerated all simple cycles
to compute the product spectral radius, and checked

```math
\left|\rho(T_{v_1}\cdots T_{v_k})-
 \sum_s\max_j(p_{v_s}(j)+a_{v_{s+1}}(j))\right|
\leq\sum_s\epsilon_{v_s}.
```

All checks passed.

The checked-in verifier is useful but not a proof-complete test suite.  In
particular, its ARS.1 random test takes `eta=0`, and it does not directly test
the approximate spectral inequality ARS.15b or the `Delta` cycle
decomposition in ARS.15c--d.  Those omissions are covered by the proofs and
the independent random check above; they are not failures of the claims.

## 2. ARS.1: PASS

From

```math
T_v(k,j)=a_k+p_v(j)+E_{kj},\qquad |E_{kj}|\leq\epsilon,
```

right multiplication by the final block gives

```math
T_{uv}(i,j)=p_v(j)+
\max_k\{T_u(i,k)+a_k+E_{kj}\}.
```

If `b=max_k(T_u(i,k)+a_k)`, the last maximum differs from `b` by a
number in `[-epsilon,epsilon]`, uniformly in `j`.  Projective code error
`eta` therefore gives ARS.6 with `epsilon+eta`.

For the declared normalized terminal query

```math
\mathcal R_z(u)=\max_j(u_j+z_j)-\max_j u_j,
```

each maximum changes by at most `epsilon+eta`; hence the factor `2` in
ARS.7 is correct.  It cannot generically be removed: the terminal field may
select a coordinate with positive residual error while the unperturbed
normalizing maximum selects one with negative residual error.  The bound is
independent of the magnitude of `z`, as claimed.

The suffix-state count and right-congruence condition are also correct.
This theorem answers a scalar-normalized *rooted terminal* query and does not
control accumulated scalar reward.

## 3. ARS.2 and ARS.2a: PASS, with a scope clarification

For exact factors `a_e\otimes p_e`, multiplication gives

```math
(a_e\otimes p_e)(a_f\otimes p_f)
=\max_j(p_e(j)+a_f(j))+a_e\otimes p_f.
```

Iterating and closing the rank-one product at its spectral cycle proves
ARS.12 exactly.  A finite defect graph has depth-independent cyclic error if
and only if every repeatable directed cycle has zero defect sum.  On a
strongly connected component this is exactly the coboundary condition
`d(e,f)=psi(f)-psi(e)`; acyclic inter-component pieces contribute only a
bounded transient.

For approximate factors, max-plus multiplication and max-plus spectral
radius are each one-Lipschitz in entrywise sup norm.  Successive replacement
therefore costs at most `sum epsilon_v`, proving ARS.15b.  Every closed walk
decomposes into simple directed cycles, so the absolute defect is at most
`Delta` times its block length.  Division by block length `D` proves the
rate in ARS.15d.

**Recommended repair R1.**  State explicitly that the “legal block graph”
is the finite directed graph whose closed walks are exactly the cyclic block
words under consideration (the complete graph when all concatenations are
legal).  Define `Delta=0` if this graph has no directed cycle.  Without this
scope sentence, an arbitrary higher-order legal language need not be
captured by pairwise letter adjacency.  The theorem is correct for the
finite graph/regular presentation it actually uses.

Also make the bounded-transient sentence for a nonmultiple of `D` explicit:
there are only finitely many leftover products of length below `D`; adjoining
one to the approximate rank-one block product changes the reference endpoint
term by a bounded amount while the accumulated entrywise error remains the
sum of block errors.

## 4. ARS.3: PASS and the constant is sharp

For the displayed matrices, the compatibility table is indeed

```math
\delta\begin{pmatrix}0&1\\1&1\end{pmatrix}.
```

Thus a cyclic word of length `t` has response
`delta(t-N_AA)`.  The two profiles have projective distance `delta/2`, so
their optimal one-centre radius is `delta/4`.

A one-dynamic-state predictor may still assign letter-dependent tolls
`g_A,g_B`.  The cycles `A`, `B`, and `AB` force

```math
\max\{|g_A|,|g_B-\delta|,
|\tfrac12(g_A+g_B)-\delta|\}\geq\delta/4.
```

Equality holds at `(g_A,g_B)=(delta/4,5delta/4)`.  For a mixed cyclic word,
if `k` is its number of cyclic `A`-runs, the total defect is
`delta(k-t/4)` and `0<=k<=t/2`; constant words attain the other endpoint.
Therefore ARS.20 is exact, and repeating an offending cycle rules out a
depth-uniform absolute error.

**Recommended repair R2.**  Replace “one scalar per-letter state” by
“one dynamic state with a letter-dependent scalar toll.”  This avoids
confusion with a quotient compatibility table having literally one
state-pair entry.  The proof already uses the stronger, standard
one-state weighted-automaton interpretation.

## 5. ARS.4: PASS, with a sharper available constant

Positive scaling commutes with max-plus multiplication and spectral radius,
so ARS.22 is exact.  For products of length at least two, changing the first
and last edges of an optimal path changes its weight by at most `2alpha`,
which proves the displayed global entry-diameter bound.

In fact the residual-radius conclusion can be strengthened.  For a fixed
initial row, changing only the last edge shows that the row's terminal span
is at most `alpha`.  Centering each row therefore yields

```math
\operatorname{rad}(T_w;0)\leq\alpha/2
```

for every nonempty word, not merely `alpha`.  The weaker ARS.23 is valid and
the universality conclusion remains intact: below a constant multiple of the
shell scale, the scalar response can still contain an arbitrary scaled
all-finite weighted-automaton algebra.

**Recommended repair R3.**  Use the sharper `alpha/2` residual radius (or
say explicitly that `alpha` is a deliberately loose global-diameter bound).

## Repair resolution

Before canonicalization the draft was restricted to finite directed
legal-word presentations, “one state” was clarified to retain
letter-dependent tolls, and ARS.4 was sharpened to the rowwise
`alpha/2` radius.  The verifier now checks that sharper bound for every
nonempty tested word.  These changes implement R1--R3 without altering any
of the audited proofs.

## 6. Empty support core: PASS

For the displayed pair,

```math
pT_a=pT_b=p,
```

so every product has the finite left eigenprofile `p` with eigenvalue zero
and hence spectral radius zero.  The zero-threshold relations are
`{(0,0)}` and `{(1,0),(1,1)}` in zero-based coordinates.  The descending
common-core iteration first leaves `{0}` and then becomes empty.  This
correctly demonstrates that failure of this particular one-support witness
is not semantic scalar drift.

## 7. Novelty and canonical overlap

The draft should not present every component as new:

* ARS.2's exact cycle/coboundary statement is a specialization of canonical
  Theorem 17.1l and the residual-cycle part of Theorem 17.1u.
* The empty-core example is deliberately the mandatory counterexample from
  Theorem 17.1u.
* One-Lipschitz perturbation and cycle deletion in ARS.2a are elementary
  ingredients already used elsewhere in the dynamic theory.

The genuinely new combined content is nevertheless non-duplicative:

1. ARS.1 gives a depth-uniform *approximate* last-window theorem for the
   normalized terminal query without assuming contraction.
2. ARS.3 gives a sharp two-letter separation between zero projective
   contraction and positive scalar-response rate.
3. ARS.4 proves that a uniformly tiny one-profile residual shell can retain
   an arbitrary finite weighted response algebra.

Together these establish the claimed structural distinction: terminal
profile error is paid once, while unresolved scalar compatibility is paid on
repeatable cycles.  Subject to R1--R3, the draft is suitable for
canonicalization.
