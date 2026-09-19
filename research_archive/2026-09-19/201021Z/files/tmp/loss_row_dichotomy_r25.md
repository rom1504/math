# Wave 25 Route 2: loss--row dichotomy and retained-deficit collapse

## Status

The switching identities, fixed-size row moments, integer-cost extraction,
loss-Laplace coarea theorem, arbitrary-cut retained-deficit inequality, and
the `A_6,A_8,A_9` audits below are **Verified**.  The real-cost scalar game is
an **Abstract falsifier** of extraction without a cost gap; it is not a
signing.  The missing positive premise is an overlap lower bound for one
low-row-square exact ground, or, in the strictly weaker formulation found
below, for one arbitrary oriented cut.  No convergence theorem is proved.

The checks are in

```text
tmp/loss_row_dichotomy_r25.py
tmp/all_cut_effective_r25.py
```

and pass under `.venv/bin/python`.

## Main result first

The sharpest verified reduction from this route is not the old two-ground
parent-ground LP.  For **every** full oriented cut `d`, not necessarily a
parent ground, define

```math
\Delta(d)=q-\langle A,d\rangle,
\qquad
B_{n,m}=\left[(m/n)^{3/2}-p_2\right]q,
```

and

```math
\widehat\ell(S,d)
=Q(A[S])-c_A(S,d)-p_2\Delta(d)-B_{n,m}.
\tag{R25.0}
```

The exact centered decomposition and (10.759)--(10.760) show that it is
sufficient to find one cut with

```math
R_2(d)=O(n^{9/4-c}),\qquad
U_m\{\widehat\ell(S,d)\le O(n^{3/2-c})\}
\ge e^{-O(n^{3/4-c})}.
```

No nonnegative-row or exact-ground property is used in this deduction.  If a
mixed captured law is found instead, integrality of `R_2(d)` extracts one
such cut with only a one-unit cost enlargement and a polynomial coverage
loss.  Thus the old two-ground target is asymptotically unnecessary.  The
full derivation is in Section 4; Sections 1--3 record what remains valid in
the narrower parent-ground/Laplace route.

## 1. What switching gives, and what it still does not give

Fix an exact positive parent ground `g`, switch it to the all-one state, and
delete `T=[n]\setminus S`, where `|T|=k`.  Write

```math
r_i=\sum_{j\ne i}w_{ij},\qquad
R_T=\sum_{i\in T}r_i,\qquad
c_g=R_2(g)=\sum_i r_i^2.
```

The already identified exact identity is

```math
\ell(S,g)=R_T+b_T-d_T.
\tag{R25.1}
```

For a ground, `r_i\ge0` and `\sum_i r_i=q`.  Uniform fixed-size sampling
therefore gives the additional exact moment identities

```math
\boxed{
\mathbb E R_T=\frac{kq}{n},\qquad
\mathbb E R_T^2
=\frac{k}{n}c_g
+\frac{k(k-1)}{n(n-1)}(q^2-c_g).}
\tag{R25.2}
```

These are useful audits of the proposed row-square scale.  They do not close
the lower tail of (R25.1): `b_T-d_T` can correlate with `R_T`, and
`\mathbb E\ell` is exactly the unknown restriction excess.  Thus using a
bound on the mean loss here would be circular.  No hard row cap is introduced.

## 2. Integer row squares collapse the proposed two-ground dichotomy

The decisive elementary point is that every row-square cost is a
nonnegative integer.  The following statement does not require an optimizer
or a two-point support theorem.

### Integer extraction lemma

Let actions `g` have coverages `u_g\in[0,1]` and costs
`c_g\in\mathbb Z_{\ge0}`.  Let `\nu` be any prior with

```math
Z=\sum_g\nu_gu_g>0,
\qquad
\mu_g=\frac{\nu_gu_g}{Z},
\qquad
\sum_g\mu_gc_g\le C.
\tag{R25.3}
```

Put `B=\lfloor C\rfloor`.  Then

```math
\mu\{c\le B\}
\ge \frac{B+1-C}{B+1},
\tag{R25.4}
```

and consequently some action with `c_g\le B` obeys

```math
u_g\ge
Z\frac{B+1-C}{B+1}.
\tag{R25.5}
```

Indeed, if `p=\mu\{c\le B\}`, then integer costs give

```math
C\ge\mathbb E_\mu c\ge(1-p)(B+1).
```

Also, the captured overlap on the low set is `Zp`, while
`\sum_{c_g\le B}\nu_gu_g\le\max_{c_g\le B}u_g`.

The factor in (R25.5) can be arbitrarily small when a real `C` lies extremely
close below its next integer.  The uniform version needed here simply rounds
the allowed cost upward.  With `N=\lceil C\rceil`, some action satisfies

```math
\boxed{
c_g\le N,
\qquad
u_g\ge\frac{Z}{N+1}.}
\tag{R25.6}
```

This also covers `C<1`; for `C=0` it gives a cost-zero action with coverage at
least `Z`.  If `C` is already an integer, no budget enlargement occurs.

At the target scale `C=O(n^{9/4-c})`, the loss in log coverage is only
`O(\log n)`.  Hence

```math
Z\ge \exp\{-O(n^{3/4-c})\}
```

for any captured mixture already forces a single ground with the same
stretched-exponential coverage scale and low/moderate row square.  The
suggested alternative in which high-cost grounds rescue an exponentially
bad low-cost ground is therefore not a distinct asymptotic branch.

For comparison, if an extreme captured law uses
`c_-<C<c_+`, then exactly

```math
\mu_- =\frac{c_+-C}{c_+-c_-},\qquad
\mu_+=\frac{C-c_-}{c_+-c_-},\qquad
\frac1Z=\frac{\mu_-}{u_-}+\frac{\mu_+}{u_+}.
\tag{R25.7}
```

Without integrality this does admit an abstract wall.  Take `C=1`,
`(u_-,c_-)=(e^{-L},0)` and `(u_+,c_+)=(1,1+e^{-L})`.  The boundary mixture has

```math
Z=\frac{1+e^{-L}}2,
```

although the only below-budget action has exponentially small coverage.  If
the high cost is changed to the integer value `2`, then
`Z=2e^{-L}/(1+e^{-L})`.  Thus the near-budget continuous example explains
the apparent two-ground escape, while the actual integer cost lattice rules
it out up to a polynomial factor.

## 3. A rigorous loss-Laplace/coarea theorem

For the moment retain exact grounds and the nonnegative loss
`L_g(S)=\ell(S,g)`.  Define

```math
u_g(t)=U_m\{L_g\le t\},\qquad
v_g(\lambda)=\mathbb E_{U_m}e^{-\lambda L_g}.
```

For a prior `\nu` (possibly a subprobability after omitting a dummy action),
put

```math
K=\sum_g\nu_gv_g,
\qquad
K C_L=\sum_g\nu_gv_gc_g,
\qquad C_L\le C.
\tag{R25.8}
```

Because `L_g\ge0`, the layer-cake identities are exact:

```math
\boxed{
K=\lambda\int_0^\infty e^{-\lambda t}Z_\nu(t)\,dt,
\qquad
K C_L=\lambda\int_0^\infty e^{-\lambda t}N_\nu(t)\,dt,}
\tag{R25.9}
```

where

```math
Z_\nu(t)=\sum_g\nu_gu_g(t),
\qquad
N_\nu(t)=\sum_g\nu_gu_g(t)c_g.
```

Since `v_g\le1` and `\sum_g\nu_g\le1`, one has `0<K\le1`.

### Coarea recovery theorem

For every `a>2`, there is a threshold

```math
0\le t\le T_a:=\frac1\lambda\log\frac aK
\tag{R25.10}
```

such that

```math
\boxed{
Z_\nu(t)\ge K\left(1-\frac2a\right),
\qquad
\frac{N_\nu(t)}{Z_\nu(t)}\le aC.}
\tag{R25.11}
```

Proof: the part of the first integral above `T_a` is at most
`e^{-\lambda T_a}=K/a`.  Thresholds whose captured cost exceeds `aC`
contribute at most

```math
\frac1{aC}\lambda\int_0^\infty e^{-\lambda t}N_\nu(t)\,dt
\le\frac Ka
```

to the first integral.  The remaining thresholds in `[0,T_a]` therefore
carry weighted `Z`-mass at least `K(1-2/a)`.  The measure
`\lambda e^{-\lambda t}dt` of this interval is at most one, so one remaining
threshold has the claimed `Z` value.  The case `C=0` follows separately from
the vanishing of the nonnegative cost integral.

Combining (R25.11) with integer extraction gives one exact ground with

```math
\boxed{
c_g\le\lceil aC\rceil,\qquad
u_g(t)\ge
\frac{K(1-2/a)}{\lceil aC\rceil+1}.}
\tag{R25.12}
```

For instance `a=4` gives hard overlap at least `K/2`, captured cost at most
`4C`, and then individual coverage at least
`K/[2(\lceil4C\rceil+1)]`.  The earlier constants `K/4,4C` are valid but not
optimal; they result from discarding `K/2` rather than `K/4` in the upper
threshold tail.

The natural missing soft premise is the lower bound

```math
\boxed{
K_*(C;\lambda):=
\max_\nu\left\{
\sum_g\nu_gv_g:
\sum_g\nu_gv_g(c_g-C)\le0
\right\}
\ge e^{-O(n^{3/4-c})}.}
\tag{R25.13}
```

The dummy action is understood.  At
`\lambda\asymp n^{-3/4}` and `C=O(n^{9/4-c})`, (R25.10)--(R25.12) give

```math
t=O(n^{3/2-c}),\qquad
c_g=O(n^{9/4-c}),\qquad
-\log u_g(t)=O(n^{3/4-c}).
```

Thus (R25.13) closes the old exact-ground implementation without a loss-mean
bound or a hard row cap.  Nothing in the scalar switching identities proves
(R25.13); that is the explicit missing premise.

There is no universal monotonicity shortcut from a cost tilt.  For a base law
`w`, let

```math
\nu_\eta(g)=\frac{w_ge^{-\eta c_g}}{\sum_hw_he^{-\eta c_h}},
\qquad Z_\eta=\sum_g\nu_\eta(g)u_g.
```

Its captured mean cost is the negative logarithmic derivative of
`\sum_gw_gu_ge^{-\eta c_g}`, while exactly

```math
\frac{d}{d\eta}\log Z_\eta
=\mathbb E_{\nu_\eta}c-
\mathbb E_{\nu_\eta(\cdot\mid\mathrm{capture})}c
=-\frac{\operatorname{Cov}_{\nu_\eta}(u,c)}
        {\mathbb E_{\nu_\eta}u}.
\tag{R25.14}
```

The `A_9` uniform-ground covariance is positive at thresholds `0` and `4`
and negative at threshold `8`, so even its sign is not structurally fixed.

## 4. Stronger finding: retain parent deficit and the coefficient slack

The exact-ground restriction in (R25.13) is unnecessary for the analytic
information/mgf step.  Let `d` now be any full oriented cut and set

```math
\Delta(d)=q-\langle A,d\rangle\ge0,
\qquad
X(S,d)=c_A(S,d)-p_2\langle A,d\rangle,
\qquad
\ell(S,d)=Q(A[S])-c_A(S,d).
```

Then the following identity is exact for every selector and every cut:

```math
\boxed{
Q(A[S])-p_2q
=X(S,d)+\ell(S,d)-p_2\Delta(d).}
\tag{R25.15}
```

The domain-free slice mgf (10.759) is valid for these arbitrary cuts.  Its
proof uses the signed row square `R_2(d)=\sum_i r_i(d)^2`, the cross Frobenius
term, and the operator norm.  It never uses `r_i\ge0`, `\sum_i r_i=q`, or
exact parent minimality.  Moreover
`|\langle A,d\rangle|\le q` and switching preserves the operator norm for
every oriented cut.

Applying the optimized-reference identity (10.760) directly to an arbitrary
joint law `P` therefore retains the negative deficit term.  With
`p=m/n`, on the halved spectral domain,

```math
\boxed{
\begin{aligned}
V_{\rm ad}(A,m)-p_2q
\le{}&\mathbb E_P[\ell(S,D)-p_2\Delta(D)]
+\epsilon_{n,m}q\\
&+\frac{I(S;D)+D(P_S\Vert U_m)+\chi_{n,m}}{\lambda}
+\lambda\left[p^2\mathbb E_PR_2(D)+\frac{n^2}{2}\right].
\end{aligned}}
\tag{R25.16}
```

There is a second slack which (10.764) did not use.  Put

```math
B_{n,m}=\left[p^{3/2}-p_2\right]q\ge0,
\qquad
\widehat\ell(S,d)
=\ell(S,d)-p_2\Delta(d)-B_{n,m}.
\tag{R25.17}
```

Subtracting `B_{n,m}` from (R25.16) gives exactly

```math
V_{\rm ad}(A,m)-p^{3/2}q
\le \mathbb E_P\widehat\ell(S,D)
+\text{the same information/mgf costs}.
\tag{R25.18}
```

This verifies both the sign and normalization.  There is no hidden
exact-ground hypothesis.  For a parent ground `\Delta=0`, the favorable event
at residual tolerance `t` is

```math
\ell(S,g)\le B_{n,m}+t,
\tag{R25.19}
```

not the much stronger `\ell(S,g)\le t` used in (10.764).

Consequently the weakest current implementation is the following
single-cut statement.  Find one arbitrary oriented cut `d` such that

```math
\boxed{
\begin{aligned}
R_2(d)&=O(n^{9/4-c}),\\
U_m\{S:\widehat\ell(S,d)\le t\}
&\ge \exp\{-O(n^{3/4-c})\},\\
t&=O(n^{3/2-c}).
\end{aligned}}
\tag{R25.20}
```

Indeed, condition `U_m` on the displayed event and take the output cut to be
constant.  Then `I(S;D)=0`, the selector KL is exactly minus the logarithm of
the displayed coverage, and `\mathbb E\widehat\ell\le t`; (R25.18) gives the
desired restriction edge.

One may still formulate a captured-mass LP over all cuts, but (R25.6) shows
that any such mixed witness yields a single cut satisfying (R25.20), with at
most a one-unit row-cost enlargement and an `O(\log n)` loss in log coverage.
Thus mixtures are asymptotically unnecessary in both the old exact-ground
and the new all-cut formulations.

The loss `\widehat\ell` can be negative.  If a Laplace/coarea version is
wanted, it must use

```math
L_+(S,d)=[\widehat\ell(S,d)]_+.
```

For `t\ge0`, its sublevel event is exactly
`\{L_+\le t\}=\{\widehat\ell\le t\}`, and Section 3 applies.  Applying the
`[0,\infty)` layer cake directly to the signed `\widehat\ell` would be
incorrect.

## 5. Exact finite audits

The ground-only checker reproduces all one-deletion types in (10.77), checks
(R25.1)--(R25.2) for every `m\ge\lceil n/2\rceil`, and verifies the coarea
and integer-extraction conclusions for `\lambda=1/4,1/2`.

The arbitrary-cut checker enumerates every oriented projective cut and every
selector at the same sizes.  It checks (R25.15) in integer arithmetic, checks
the arbitrary-cut version of (10.759), and solves the cost LPs.  At one
deletion and residual threshold zero:

| signing and budget | exact-ground `\ell\le0` | all-cut `\ell-p_2\Delta\le0` | all-cut full slack `\widehat\ell\le0` |
|:---|---:|---:|---:|
| `A_6`, `C=30` | `5/6` | `5/6` | `5/6`, `(R_2,\Delta)=(30,0)` |
| `A_8`, `C=64` | `3/8` | `1/2` | `5/8`, `(R_2,\Delta)=(40,20)` |
| `A_9`, `C=80` | `1/9` | `8/21` | `4/9`, `(R_2,\Delta)=(40,24)` |
| `A_9`, `C=88` | `4/27` | `32/81` | `4/9`, `(R_2,\Delta)=(40,24)` |
| `A_9`, `C=96` | `2/9` | `16/39` | `4/9`, `(R_2,\Delta)=(40,24)` |
| `A_9`, `C=104` | `4/15` | `32/75` | `4/9`, `(R_2,\Delta)=(40,24)` |
| `A_9`, `C=112` | `1/3` | `4/9` | `4/9`, `(R_2,\Delta)=(40,24)` |

The `A_9` middle value uses two cuts, while the full-slack value is already
attained by the displayed single cut.  These finite examples do not prove an
asymptotic theorem, but they show that retaining both the parent deficit and
the coefficient slack is a strict, material relaxation rather than a formal
rewrite.

## Route conclusion

**Proved:** a captured low-row-square mixture cannot hide all of its coverage
in high-cost grounds at the relevant exponential scale; integer row squares
extract one favorable low/moderate-cost action.  A joint loss-Laplace lower
bound would rigorously recover the required hard threshold by coarea.

**Superseding verified reduction:** the entropy and row-square mgf argument
works for arbitrary cuts and retains both `-p_2\Delta` and the full leading
coefficient slack `B_{n,m}`.  The sharp current endpoint is therefore the
single-cut overlap statement (R25.20).

**Still open:** neither switching identities nor the coarea argument proves
the overlap lower bound in (R25.13) or (R25.20).  Proving (R25.20), or finding
an asymptotic signing family that falsifies it, is the next substantive task.
