# Wave 37 harmonic attack: dynamic selector erasure and matched mode transfer

## Status

- **Verified:** along the actual parent interpolation, the endpoint selector
  posterior gives an exact one-temperature-back conditional law: after
  conditioning on a selector which omits vertex `i`, the `i`th bit at time
  `t` has the parent conditional law at time `t-1`.
- **Verified:** this yields a dynamic version of the endpoint cancellation
  `J=C+I`, and convexity compares the ordinary binary cost from time `0` to
  the dynamic selector-erasure cost from time `t-1`.
- **Verified sufficient lemma:** actual adverse vertex migration is bounded
  by an explicit transient selector-erasure flux. Bounding the sum of these
  fluxes at `O(n^(1/2-2c))` is enough for the migration part of the harmonic
  package.
- **Open:** no exact-minimizer estimate of that flux is proved. Endpoint
  cost is only its `t=1` boundary and does not control it.
- **Numerical:** finite `A_4,A_6,A_8,A_9` audits show a large cancellation
  between omitted external-partition movement and exclusion-shape movement.
  Separately bounding those two pieces would be very wasteful.

The checker `tmp/matched_no_transient_r37_check.py` verifies every finite
identity, the dynamic cancellation, and the convexity inequality. Global
orientation, restoring correlation, and adjacent-selector Hellinger remain
separate throughout.

## 1. The omitted component at time `t` is the parent bit at time `t-1`

Write the matched external-partition weights as

```math
w_S(d)=e^{-F_{\beta,S}(d[S])},\qquad
U(d)=\mathbb E_{S\sim U_m}w_S(d),\qquad
f(d)=U(d)/\mathbb E_\nu U.
```

Let `pi_d(S)=w_S(d)/sum_T w_T(d)` be the endpoint selector posterior and

```math
\mu_t(d)=\frac{\nu(d)f(d)^t}{\mathbb E_\nu f^t},
\qquad P_t(S,d)=\mu_t(d)\pi_d(S).
\tag{R37.H1}
```

Fix a nonreference vertex coordinate `i`, a chart context
`e=D_{-i}`, and a selector `S` omitting `i`. The matched weight `w_S` is
constant across the `i`-edge. Therefore

```math
P_t(D_i=b\mid S,e)
\ \propto\ \nu(b\mid e)f(e,b)^t\frac{w_S(e)}{U(e,b)}
\ \propto\ \nu(b\mid e)f(e,b)^{t-1}.
```

The last expression is exactly the parent conditional at the previous
interpolation parameter:

```math
\boxed{
P_t(D_i=\cdot\mid S,e)=\mu_{t-1}(D_i=\cdot\mid e)
\quad(i\notin S).
}
\tag{R37.H2}
```

This holds for every real `t` for which the finite-space laws are defined;
in particular `t in [0,1]`. At `t=1` it reduces to the endpoint fact that
an omitted selector uses the base conditional `mu_0=nu`.

## 2. Dynamic selector cancellation

Put

```math
Q_{t,e}=\mu_t(D_i=\cdot\mid e),\qquad
q_{i,e}^{\leftarrow}(t)=D(Q_{t,e}\Vert Q_{t-1,e}).
\tag{R37.H3}
```

Use the previous-temperature bit law as a common reference and define the
lifted dynamic cost

```math
J_i^{\leftarrow}(t)
=\mathbb E_{P_t(S,e)}
D(P_t(D_i\mid S,e)\Vert Q_{t-1,e}).
```

The KL chain rule and (R37.H2) give exactly

```math
\boxed{
J_i^{\leftarrow}(t)
=\mathbb E_{e\sim(\mu_t)_{-i}}q_{i,e}^{\leftarrow}(t)
+I_{P_t}(S;D_i\mid D_{-i}),
}
\tag{R37.H4}
```

and every `S` omitting `i` contributes zero to the lifted cost. At `t=1`,
`Q_{t-1}=Q_0=nu_i`, so (R37.H4) is exactly (10.1010), with the first term
equal to `C_i`.

Thus selector information remains the subtractive cancellation at every
interpolation time, provided the reference is moved back by one unit. The
fixed-size identity `sum_i r_i(d)=n-m` also remains true under `P_t`, but, as
at the endpoint, it does not upper-bound the resulting logarithmic costs.

## 3. Convexity converts ordinary binary cost to dynamic erasure cost

For one context, `Q_s` is the one-dimensional exponential family obtained
by tilting `Q_0` with the two values of `g=log f`. Write

```math
H_e(s)=\log\mathbb E_{Q_{0,e}}e^{s g}.
```

The fixed-base binary cost in (10.980) and the one-step cost are

```math
k_e(t)=D(Q_{t,e}\Vert Q_{0,e})=tH'_e(t)-H_e(t),
```

```math
q_{i,e}^{\leftarrow}(t)
=D(Q_{t,e}\Vert Q_{t-1,e})
=H'_e(t)+H_e(t-1)-H_e(t).
```

Since `0=t(t-1)+(1-t)t`, convexity of `H_e` and `H_e(0)=0` give

```math
\boxed{
tq_{i,e}^{\leftarrow}(t)-k_e(t)
=tH_e(t-1)+(1-t)H_e(t)\ge0,
\qquad 0\le t\le1.
}
\tag{R37.H5}
```

This inequality is exact, including the single factor `t`. It is valid for
the orientation edge as a bare exponential-geodesic statement, but only
vertex edges have the selector-erasure interpretation (R37.H2)--(R37.H4).

## 4. A dynamic transient-flux sufficient lemma

For coordinate `i`, let `M_{t,i}(e)` be its context law. As in (10.1017),

```math
\operatorname{Cov}_{M_{t,i}}(A'_e,k_e)
=\sum_eM'_{t,i}(e)k_e(t).
```

Dropping contexts whose mass is increasing and using (R37.H5) proves

```math
\boxed{
\mathcal A_i
\le\mathfrak T_i
:=\int_0^1t\sum_e
q_{i,e}^{\leftarrow}(t)[-M'_{t,i}(e)]_+\,dt.
}
\tag{R37.H6}
```

Consequently the concrete matched no-transient-mode-transfer lemma

```math
\boxed{
\sum_{i\ne v_*}\mathfrak T_i=O(n^{1/2-2c})
}
\tag{R37.H7}
```

plus the separate orientation estimate
`A_0=O(n^(1/2-2c))` proves the adverse-migration part of (10.983).
Unlike (10.1019), (R37.H7) charges a context by its **current** dynamic
selector-erasure cost and suppresses early movement by the correct factor
`t`. Equation (R37.H4) says that this current cost still contains the full
selector-information cancellation.

This is an exact alternative majorant, not a uniformly smaller one:
`t q_e^{\leftarrow}(t)` and the endpoint cost `K_e` have no proved ordering.

This is a sharper sufficient target, not a proof. A falsifier for this
specific implementation is an unbounded exact-minimizer family with
`sum_i \mathfrak T_i=omega(n^(1/2-2c))`. Such a family would not by itself
falsify the harmonic route, because (R37.H6) may be loose. A harmonic
migration falsifier must instead make `sum_i A_i` itself exceed the target.

## 5. The exact matched external-partition remainder

The preceding cost has an explicit external-partition representation. Put

```math
B_i(e)=\mathbb E_{S\sim U_m}
[\mathbf1_{\{i\notin S\}}e^{-F_{\beta,S}(e)}],
\qquad r_i(d)=B_i(D_{-i})/U(d).
\tag{R37.H8}
```

Then `0<B_i(e)<=(n-m)/n`, and, up to the common normalization of `f`,
`f(d)r_i(d)=B_i(D_-i)`. Define

```math
\mathcal H_{i,e}(t)
=\log\mathbb E_{\nu(D_i\mid e)}r_i(e,D_i)^{-t}.
```

The conditional log partition, its score, and its binary cost split exactly
as

```math
\boxed{
A_{i,e}(t)
=t\log\frac{B_i(e)}{\mathbb E_\nu U}+\mathcal H_{i,e}(t),
\quad
A'_{i,e}(t)
=\log\frac{B_i(e)}{\mathbb E_\nu U}+\mathcal H'_{i,e}(t),
\quad
k_{i,e}(t)=t\mathcal H'_{i,e}(t)-\mathcal H_{i,e}(t).
}
\tag{R37.H9}
```

Therefore the vertex migration covariance is

```math
\boxed{
\operatorname{Cov}(A'_e,k_e)
=\operatorname{Cov}(\log B_i(e),k_e)
+\operatorname{Cov}(\mathcal H'_{i,e}(t),k_e).
}
\tag{R37.H10}
```

The first term is movement of the actual omitted-selector external
partition; the second is movement of the exclusion shape inside a context.
This is the exact minimizer-specific remainder absent from a purely endpoint
argument. A proof may use that every `B_i` is an average of matched external
partition functions of principal restrictions, not an arbitrary context
weight. No useful lower bound or covariance sign for `B_i` is currently
proved.

## 6. Finite audit: the two matched pieces nearly cancel

Gauss--Legendre quadrature of the proved split (R37.H10) gives the following
integrated signed vertex covariances `(total, omitted B, exclusion shape)`:

```text
A6, m=3, beta=1/2:  ( 0.000005609,  0.000333774, -0.000328165)
A8, m=5, beta=2:    (-0.00744106,   0.02980583,  -0.03724690)
A9, m=7, beta=2:    (-0.00146794,   0.08670755,  -0.08817549)
```

For `A_9,m=7,beta=2`, the actual adverse vertex integral is
`0.00285723`, whereas separately taking the adverse exclusion-shape part
would charge `0.08817549`, over thirty times larger. On all sampled
`A_4,A_6,A_8,A_9` nodes the omitted-partition term was nonnegative up to
quadrature/roundoff tolerance and the exclusion-shape term nonpositive, but
this is **numerical only** and no sign theorem is asserted.

The dynamic identities (R37.H2)--(R37.H5) were independently checked on
`A_6,m=3,beta=1/2` and `A_9,m=7,beta=2`, with maximum errors below
`7e-15`. The flux (R37.H6) was also checked coordinatewise. On the latter
example its two nonzero vertex contributions are about `0.00273519` each,
versus actual adverse contributions about `0.00142861` each. Thus the new
flux retains the correct finite scale here, although this is not asymptotic
evidence.

## Frontier

The endpoint-erasure theorem now extends exactly through the interpolation:
selector omission replaces the current bit law by the parent law one unit
back, and selector information remains the subtractive term. The precise
new sufficient quantity is the transient selector-erasure flux (R37.H7).
It is absent from endpoint data and is large in the redundant-mode wall.

A successful next theorem must exploit exact signing minimality and the
matched principal-restriction form of `B_i` to control the **combined**
movement in (R37.H10), or prove (R37.H7) directly. Bounding the two terms of
(R37.H10) separately is discouraged by the finite cancellation. Global
orientation has no omission identity; restoring and adjacent-selector
Hellinger remain separate. No convergence proof or actual-minimizer
falsifier is obtained.
