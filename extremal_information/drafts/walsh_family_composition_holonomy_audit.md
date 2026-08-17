# Independent audit: Walsh-family composition and holonomy

Audited file:
[`walsh_family_composition_holonomy.md`](walsh_family_composition_holonomy.md)

Verdict: **ACCEPT with minor presentational repairs.**  The carrier identity,
gauge transport, commutation bit, bipartite saturating witness,
anticommuting spectral bound, and extensive path separation are all correct.
The exact verifier passes 39,212 scalable algebra checks and the complete
`n=4` wind tunnel.  No proof-level defect was found.

## 1. Child and carrier normalization

Write `U=R/sqrt(q)`.  In `(u,v)` order the switch used in the bounded-cap
family satisfies

```math
s_g(u,v)b(u,v)=(-1)^{g(v)}.
```

Thus, with `D_g` of order `q`,

```math
C_g=(I_q\otimes D_g)(R\otimes R)(I_q\otimes D_g)
=R\otimes(D_gRD_g),
```

which proves (WC.1).  Since `R^2=qI`, `C_g^2=nI`, so its operator norm is
`q=sqrt(n)`.  Moreover

```math
tr(C_g)=tr(R)tr(D_gRD_g)=0
```

for `m>=1`, because `tr(R)=sum_a(-1)^(a.a)=0`.  Deleting the diagonal
therefore changes the Boolean quadratic by the constant zero.  The Boolean
pole from the original construction attains `qn/2`, while the spectral
bound gives absolute value at most `qn/2`; (WC.2) is exact.

For the global symmetric block matrix, diagonal blocks are `C_gi` and edge
blocks are `W=R tensor R`.  Reordering `(i,u,v)` as `(u,i,v)` factors out
the first `R` and leaves exactly the carrier `K` in (WC.4).  Because a
quadratic form counts symmetric off-diagonal blocks twice, the prefactor
`1/2` recovers each bridge once.  Hence (WC.5) has no missing factor.

The proposition proves an exact `kq`-label-bit presentation (in addition to
the graph).  It does **not** prove minimality.  The sentence calling it the
“minimal place” where a positive theorem must act should be softened to
“an explicit reduced place” unless a separate minimality argument is added.

## 2. Gauge transport

Let `D_i=I_q tensor D_gi` on the order-`n` block and `F=W/q`.  The change
`z_i=D_i x_i` gives

```math
q^{-1}E
=\frac12\sum_i z_i^TFz_i
 +\sum_{ij}z_i^TD_iFD_jz_j,
```

so (WC.6)--(WC.7) have the right normalization.  Consecutive transports
cancel their shared diagonal:

```math
(D_iFD_j)(D_jFD_l)=D_iD_l,
```

and another factor restores `F`.  Induction proves the parity formula
(WC.8).  The draft correctly distinguishes this flat edge connection from
the relative on-site obstruction below.

## 3. The commutation bit

To make dimensions explicit, let

```math
M_a=diag((-1)^(a.v):v in F_2^m),
\qquad D_a=I_q\otimes M_a,
\qquad U=R/\sqrt q.
```

If `T_a` denotes translation by `a`, the Walsh identities are

```math
UM_a=T_aU,
\qquad
M_aU=UT_a,
\qquad
T_aM_a=(-1)^(a.a)M_aT_a.
```

Now

```math
F=U\otimes U,
\qquad
\widehat C_a=U\otimes(M_aUM_a).
```

Consequently

```math
F\widehat C_a=I\otimes(T_aM_a),
\qquad
\widehat C_aF=I\otimes(M_aT_a),
```

which proves (WC.10) with exactly the displayed sign.  Both factors are
symmetric orthogonal involutions.  The notation `D_a` should be declared as
the lifted `I_q tensor M_a` at first use, but the mathematics is correct.

For `m>=2`, a nonzero even-weight `a_0` and a nonzero odd-weight `a_1`
exist.  Their linear truth tables are both balanced.  In either constant
word every pairwise label correlation is `q`, and each bias is zero, so the
claimed collision of the proposed scalar summaries is genuine.

## 4. Even and odd composition bounds

For even `a_0`, let `s` be the Boolean `+1` eigenvector of
`widehat C_a0` and put `y=Fs`.  The transform construction makes `y`
Boolean.  Commutation gives

```math
\widehat C_{a_0}y=y.
```

Assign `s,y` to the two color classes of a bipartite `G`.  Each child term
is exactly `qn/2`, and every bridge is

```math
s^TW(Fs)=qn
```

(with the same value in the reverse orientation).  The separate spectral
bounds `x^TCx/2<=qn/2` and `x^TWy<=qn` give the matching global upper bound.
This proves (WC.14), including disconnected and isolated vertices.

For odd `a_1`, anticommutation cancels the cross term in the square:

```math
\mathcal M_{a_1}^2=(I_k+A(G)^2)\otimes I_n.
```

Since `A(G)` is real symmetric,

```math
\|\mathcal M_{a_1}\|=\sqrt{1+\|A(G)\|^2}.
```

A Boolean global vector has squared norm `kn`; multiplying the Rayleigh
bound by `q/2` gives exactly (WC.15).  Attainment is neither stated nor
needed.

For `P_k`, the even coefficient is `3k/2-1`, while the odd upper coefficient
is

```math
\frac k2\sqrt{1+4\cos^2(\pi/(k+1))}.
```

Their difference is positive: check `k=2` directly; for `k>=3`,
`3-2/k>sqrt(5)` while the square-root term is strictly below `sqrt(5)`.
The limit in (WC.17) follows.  There is a literal form-feed typo in
(WC.16), where `\frac` appears as `rac`; this should be repaired.

## 5. Verifier audit

Running

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_walsh_family_composition.py
```

produced:

```text
scalable Walsh carrier/holonomy checks passed: 39212
```

followed by the complete `n=4` state counts and eight pairwise-correlation
collision classes with distinct optima.  The scalable checks use integer
arithmetic for:

- `R^2=qI`, `W^2=nI`;
- commutation/anticommutation;
- Booleanity and child eigenvector identities;
- every even-parity saturating energy;
- the full squared-matrix identity in the odd case;
- every entry of the Kronecker carrier at `q=4,8`.

The call to `eigvalsh` is only a numerical regression for the already exact
squared-matrix identity, so it does not introduce a proof dependency.

## 6. Scope verdict

The result rigorously falsifies **bias plus pairwise truth-table overlap** as
a reusable state under repeated Walsh-bridge composition.  It does not show
that all `q` truth-table bits are necessary, does not lower-bound the full
cap-`1/2` class, and does not compute the carrier maximum.  The draft states
these limitations accurately.  Subject to the two wording/typographical
repairs above, promotion to the canonical theorem record is justified.
