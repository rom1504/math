# Feasible center best response: a conditional route beyond all marked masks

Date: 2026-09-06. Author: resumed response track.
Status: the finite-matrix identity, degree-51 selection, and quantitative
Gaussian consequence below are proved. The required nonlinear conditional
transport statement is explicitly an additional hypothesis; it is being
investigated by the director and is not asserted proved here.

Subsequent same-session completion: the sufficient restricted test identities
have now been proved and independently reconstructed in
`resumed_bound_audit_restricted_channel_center_update_2026_09_06.md` and
`resumed_response_center_theorem_independent_audit_2026_09_06.md`.
The unconditional all-infinite-mask consequence is stated and proved in
`resumed_response_all_mask_uniform_escape_2026_09_06.md`. The conditional
form below is retained to identify exactly what this initial derivation used;
it does not assert an unqualified per-root conditional-mean theorem.

## 1. An exact feasible improvement mechanism requiring no zero strip

Let `B` be any symmetric hollow matrix. Suppose random vectors `F,H` obey
`0<=H<=1` and `|F|+H<=1` coordinatewise. Define

```
C=H circ sign(BF),
mu_plus=F+C, mu_minus=-F+C.
```

These are both in the cube exactly. For `E_B(x)=x^T Bx/2`,

```
[E_B(mu_plus)-E_B(mu_minus)]/2
 =F^T B C=sum_i H_i |(BF)_i|.                            (1)
```

Independent rounding preserves each hollow quadratic expectation. Hence
the absolute Boolean maximum `q(B)=max_signs |E_B(x)|` satisfies

```
q(B)>= E sum_i H_i |(BF)_i|.                              (2)
```

This optimizes the common center of the two feasible means while retaining
their opposite odd responses. It spends the existing marked amplitude `H`
instead of seeking slack where `W=0`. No local feasibility approximation,
clipping, or quadratic remainder is involved in (1)--(2).

For the old creation certificate,

```
W=UH, F=sign(W)(1-H), K=U*F,
J=E F W=E H K.
```

Here `U*` first projects onto first Gaussian chaos and then applies the
inverse creation isometry. The old center `SH` has energy `J`; the new
center in (1) can use the full local field `BF`.

## 2. Required transport hypothesis, with harmless drift retained

Consider any fixed finite old Gaussian construction and a coefficient-aligned
finite odd unmarked channel `Z_i`, of local limiting variance `v_i`. Let

```
s^2=sum_(r in D) |E[F h_r(G0)]|^2>0,
Z=sum_(r in D) (E[F h_r(G0)]/s) B h_r(BS),
v_i=(B R B)_ii,
R=sum_(r in D) (E[F h_r(G0)]^2/s^2) Q^(circ r), Q=B^2.
```

The existing Schur theorem gives `n^{-1}sum_i sqrt(v_i)>=1`. The missing
nonlinear projection statement would supply, in the required averaged local
limit, conditional means of the form

```
E[(BF)_i | X_i,S_i,Z_i]
 =T_i(X_i)+S_i K(X_i)+s Z_i,                             (3)
```

where `S_i` is uniform on signs and independent of the old Gaussian family
and `Z_i`; conditional on the old family, `Z_i~N(0,v_i)`. An arbitrary
old-field drift `T_i` is allowed. It cannot generally be omitted: already
`F=G0` gives `BF=B^2S`, which can retain nonzero old-field projections.

Only the part of (3) odd under simultaneous `(S_i,Z_i)->(-S_i,-Z_i)` is
actually needed. An equivalent sufficient test family is

```
E[M(X) psi_odd(Z) BF] = s E[M(X)] E[Z psi_odd(Z)],
E[S M(X) psi_even(Z) BF] = E[M(X) K(X)] E[psi_even(Z)],       (4)
```

with correct root-dependent variances and the averaged deterministic weights.
The existing **linear** `Z` identity is not enough to infer (3) or (4).

Indeed the desired center, when its conditional mean is used, decomposes as

```
sign(SK+sZ)
 =sign(sZ) 1{|sZ|>|K|}
  + S sign(K) 1{|sZ|<|K|}.
```

This explicitly displays why both the odd unmarked tests and the even
marked tests in (4) occur. Values on equality sets are immaterial.

## 3. Conditional quantitative Gaussian theorem

Assume (3), or the sufficient odd-part statement (4), has been established
with the approximation and integrability needed for (2). Conditional Jensen,
then symmetry and convexity in the drift `T_i`, give

```
E[H_i |(BF)_i|]
 >= E H E_(S,N)|T_i+SK+s sqrt(v_i)N|
 >= E H E_N|K+s sqrt(v_i)N|.
```

For fixed `H,K`, the function

```
Phi(t)=E H E_N|K+tN|,  t>=0,
```

is convex and nondecreasing. Jensen across root labels and the Schur
mean-standard-deviation inequality therefore give

```
liminf q(B)/n >= E H E_N|K+sN|.                            (5)
```

No bound on the largest `v_i` enters the positive Gaussian increment in (5).
A fixed operator bound can still enter the proof of transport; it does not
enter this conditional gain calculation.

Define, for `x>=0`,

```
g_s(x)=E|x+sN|-x
      =2[s phi(x/s)-x Phi(-x/s)]>0.
```

This function is positive, decreasing, and convex in `x`. Write `mu=EH` and
`p=EH^2`. Since `||K||_2^2<=||F||_2^2=1-2mu+p<=1-mu`, weighted Jensen
and Cauchy--Schwarz show

```
E H E|K+sN|-J
 >=E H g_s(|K|)
 >=mu g_s(E[H|K|]/mu)
 >=mu g_s(sqrt((1-mu)/mu)).                               (6)
```

Also `J^2<=p(1-2mu+p)<=mu(1-mu)`. Thus for `J>=j>0`,

```
mu>=m_j:=(1-sqrt(1-4j^2))/2,
gain >=m_j g_s(sqrt((1-m_j)/m_j))>0.                       (7)
```

The last step uses that `mu*g_s(sqrt((1-mu)/mu))` increases with `mu`:
both its positive prefactor and its decreasing-argument factor increase.

If an elementary positive lower expression is desired, put `a=x/s` and
`u=1/(a+1)`. Integrating only over `[a+u,a+2u]` gives

```
g_s(x)=2s integral_a^infinity (z-a)phi(z) dz
       >=2s u^2 phi(a+2u)>0.                             (8)
```

There is no claim that the tiny universal constant furnished below is sharp.

## 4. Independent improvement from degree 417 to degree 51

This section does not assume the unproved transport extension. For every
even mask `0<=H<=1`, let `F=sign(W)(1-H)`, `W=UH`, and `a=E G0 F`, with
`G0=U1`. The two-dimensional L2 Gram argument gives

```
a >= J-sqrt((1-p)(1-p-J^2/p)), p=EH^2.                    (9)
```

For completeness, write `mu=EH`; the exact covariance is
`Cov(G0,W)=mu`. Project `G0,F` orthogonally off `W`, apply Cauchy--Schwarz,
and then use `mu>=p` and `1-2mu+p<=1-p`. This proves (9) whenever `J>0`.

At `J>=43/100`, the right side of (9) is at least `141/1000`. Here is a
fully rational uniform check. With `j=.43` and `d=.289`, the square of the
radical is at most `d^2` exactly when

```
f(p):=p^3-2p^2+1.101379 p-.1849 <=0.
```

For `t=.39` and every `p in [0,1]`, Taylor's exact cubic identity gives

```
f(p)=f(t)+f'(t)(p-t)+(p+2t-2)(p-t)^2
    <=f(t)+f'(t)(p-t)-.22(p-t)^2
    <=f(t)+f'(t)^2/.88
     =-18965469/80000000000<0.
```

Thus `a>=.43-.289=.141` uniformly.

Now use the separating monomial from the earlier general-mask theorem with
`D=51`. Write `c_D=D!!` and `P_D(G)=G^D-c_D G`. It contains only Hermite
degrees `3,5,...,51`, and

```
|E F P_D(G0)| >=c_D[a-D^(-1/2)],
||P_D||_2^2=(2D-1)!!-(D!!)^2.
```

The double-factorial recurrence gives
`(2D-1)!!/(D!!)^2<=2^(D-1)`. The rational check

```
51*(14003/100000)^2>1
```

implies `1/sqrt(51)<.14003`, leaving margin `.00097`. Consequently

```
s_51:=sqrt(sum_(3<=r<=51, r odd) |E F h_r(G0)|^2)
 > 97/(100000*2^25)
 > 1/40000000000,

max_(3<=r<=51, r odd) |E F h_r(G0)|
 > 1/200000000000.                                      (10)
```

There are 25 degrees, so the last bound loses only `sqrt(25)=5`.
This improves both the fixed degree and the quantitative mass from the
previous degree-417 certificate. The coefficient-aligned finite mixture,
not the separating monomial itself, remains the transport direction.

## 5. What would follow, and why this does not justify iteration

If the nonlinear conditional transport statement were established uniformly
over fixed bounded-operator signings with the required limit order, (5)--(10)
would give a common strictly positive gain over **every** marked mask of
value at least `.43`. The gain would be independent of the fixed operator
bound. Selecting a target close to the supremum of all marked-mask values,
then using principal spectral deletion and letting its deleted fraction
tend to zero, would prove strict escape from that entire variational class.

This remains conditional on (3) or (4); neither a covariance identity nor
Gaussian local marginal independence alone establishes the needed conditional
mean. The zero-strip counterexample in
`resumed_response_high_value_zero_strip_2026_09_06.md` shows why the stronger
center-response mechanism would be consequential.

Even a successful first escape does **not** license indefinite iteration or
prove a sharp universal value. The new center is
`H sign(SK+sZ)`, not `SH` with a fresh independent own spin. Its spin is
correlated with the old fields and unmarked channel. A second iteration
therefore needs a closed transport theorem for this enlarged collection,
including its cross-root signed action. The creation-isometry identity
`W=UH` and the degree-51 argument cannot simply be reapplied after renaming
that dependent spin. No convergence result or matching upper mechanism is
claimed here.
