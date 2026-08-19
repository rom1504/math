# Exceptional bridge switches are a positive translation-convolution problem

Status: **proved exact all-order reduction and scalable switching-orbit
floor for actual own-scale pressure minimizers**.  The result treats the
exceptional aligned switch which is not covered by full-orbit arithmetic
averaging.  It shows that switching can rearrange the children’s positive
internal excess against a bipartite bridge kernel, but it can never cancel
the pure bipartite pressure itself.  It also identifies exactly what a
geometric-orbit/entropy proof can recover.

The theorem does not prove a sublinear cross-order defect.  Its obstruction
is conditional on the pure pressure of the chosen bridge templates being
too large relative to the child pressures; the presently known universal
constants do not activate that condition for optimal templates.

## 1. Notation and the rank-one switching group

For hollow sign matrices `A,D` of orders `m,n`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
\phi_A(t)=\log\mathbb E_x\cosh(tH_A(x)),
```

and define `H_D,phi_D` similarly.  For an `m` by `n` sign bridge `B`, let

```math
\psi_B(t)=\log\mathbb E_{p,q}e^{t p^{\mathsf T}Bq}.
                                                               \tag{1.1}
```

The group of rank-one sign words is

```math
\mathcal G_{m,n}
=\{pq^{\mathsf T}:p\in\{\pm1\}^m,q\in\{\pm1\}^n\},
\qquad |\mathcal G_{m,n}|=2^{m+n-1},                 \tag{1.2}
```

with entrywise multiplication.  If `g=pq^T`, write

```math
B^g=\operatorname {diag}(p)B\operatorname {diag}(q).
```

For a relative orientation `epsilon in {+-1}`, set

```math
L_{\epsilon,g}(B)=\log\mathbb E_{x,y}\cosh t\left{
H_A(x)+\epsilon H_D(y)+x^{\mathsf T}B^g y\right}.   \tag{1.3}
```

## 2. Exact positive convolution

For `Q in G_(m,n)`, choose any factorization `Q=x_Qy_Q^T` and define

```math
w_\epsilon(Q)
=\cosh t\{H_A(x_Q)+\epsilon H_D(y_Q)\},
\qquad
k_B(Q)=e^{t\langle B,Q\rangle}.                     \tag{2.1}
```

Both definitions are independent of the simultaneous gauge
`(x_Q,y_Q)->(-x_Q,-y_Q)`.  Let normalized group convolution be

```math
(f*h)(g)=\mathbb E_{Q\sim U_\mathcal G}f(gQ)h(Q).    \tag{2.2}
```

**Theorem 2.1 (exceptional-switch convolution identity).**  For every
finite order, every signing pair, every bridge, every orientation and every
switch,

```math
\boxed{
e^{L_{\epsilon,g}(B)}=(k_B*w_\epsilon)(g).}          \tag{2.3}
```

*Proof.*  Expand `cosh` with an independent `tau in {+-1}`:

```math
e^{L_{\epsilon,g}(B)}
=2^{-m-n-1}\sum_{\tau,x,y}
 e^{t\tau\{H_A(x)+\epsilon H_D(y)+x^TB^gy\}}.       \tag{2.4}
```

The map `(tau,x,y)->Q=tau xy^T` is uniform on `G_(m,n)` and has exactly
four preimages: two simultaneous gauges and two values of `tau`.  On a
fixed fibre the cross term is `t<B^g,Q>`, the gauge does not change either
internal energy, and averaging the two values of `tau` gives
`cosh t(H_A(x_Q)+epsilon H_D(y_Q))`.  Finally
`<B^g,Q>=<B,gQ>`, proving (2.3).  `square`

The identity separates a compulsory baseline from the part a switch may
rearrange.  Since `w_epsilon>=1`,

```math
e^{L_{\epsilon,g}(B)}
=e^{\psi_B(t)}+
 \mathbb E_Q k_B(gQ)\{w_\epsilon(Q)-1\}.             \tag{2.5}
```

The second term is nonnegative for every individual switch; it is not an
averaged assertion.

## 3. A floor for every exceptional switch

Equation (2.5) immediately proves the following.

**Corollary 3.1 (pure-bridge floor).**  For every `A,D,B,epsilon,g`,

```math
\boxed{L_{\epsilon,g}(B)\ge\psi_B(t).}               \tag{3.1}
```

More sharply,

```math
\boxed{
L_{\epsilon,g}(B)
\ge\psi_B(t)+\min_Q\log w_\epsilon(Q).}             \tag{3.2}
```

Indeed, convolution against `k_B` preserves constants and
`w_epsilon>=min w_epsilon`.  Thus for every collection `mathcal B` of
templates, of arbitrary cardinality,

```math
\boxed{
\min_{B\in\mathcal B,\epsilon,g}L_{\epsilon,g}(B)
\ge\min_{B\in\mathcal B}\psi_B(t).}                 \tag{3.3}
```

In particular (3.3) applies to every `exp(O(m+n))`-size union of switching
orbits.  Entropy in the switch catalogue cannot pay the pure bipartite
pressure; it can only anticorrelate the positive internal excess in (2.5).

There is also a global consequence not involving selected children.  If

```math
\Psi_{m,n}(t)=\min_B\psi_B(t),
```

then every order-`N=m+n` signing, split across these two shores, satisfies

```math
\boxed{P_N(\beta)\ge\Psi_{m,n}(\beta/\sqrt N).}      \tag{3.4}
```

For `S_r` a sum of `r` fair signs, row/column Jensen gives the explicit
all-order version

```math
\boxed{
P_N(\beta)\ge
\max\left\{
tm\,\mathbb E|S_n|-m\log2,
tn\,\mathbb E|S_m|-n\log2
\right\},\quad t={\beta\over\sqrt N}.}              \tag{3.5}
```

The scope is important: (3.4) is a lower bound, not the desired recurrence.

## 4. Exact implication for actual optimizing children

Let `A,D` now be exact own-scale minimizers for

```math
P_m(\beta)=\phi_A(\beta/\sqrt m),
\qquad
P_n(\beta)=\phi_D(\beta/\sqrt n),
```

and put

```math
t={\beta\over\sqrt N},\quad
\Delta_A=P_m(\beta)-\phi_A(t),\quad
\Delta_D=P_n(\beta)-\phi_D(t).                       \tag{4.1}
```

For any template family `mathcal B`, the direct orbit certificate is

```math
\mathcal R_{\mathcal B}(A,D)
=\min_{B\in\mathcal B,\epsilon,g}
 \{L_{\epsilon,g}(B)-P_m(\beta)-P_n(\beta)\}.       \tag{4.2}
```

Parent minimization and (2.3) give, with no surrogate law,

```math
\boxed{
E_{m,n}(\beta)
\le\mathcal R_{\mathcal B}(A,D)
=\min_{B,\epsilon,g}
 \{\log(k_B*w_\epsilon)(g)-P_m-P_n\}.}              \tag{4.3}
```

Consequently every proposed translation estimate has the required direct
arrow.  For example,

```math
\boxed{
\min_{B,\epsilon,g}\log(k_B*w_\epsilon)(g)
\le P_m+P_n+C_\beta N^{1-\delta}
\quad\Longrightarrow\quad
E_{m,n}(\beta)\le C_\beta N^{1-\delta}.}            \tag{4.4}
```

Conversely, (3.3) gives a scalable method-class obstruction:

```math
\boxed{
\min_{B\in\mathcal B}\psi_B(t)
\ge P_m+P_n+c_\beta N
\quad\Longrightarrow\quad
\mathcal R_{\mathcal B}(A,D)\ge c_\beta N.}         \tag{4.5}
```

Thus a switching-orbit architecture whose templates fail (4.5) cannot
give an `o(N)` recurrence, no matter how exceptional the selected switch
is.  Condition (4.5) is not known for pure-pressure-optimal templates and
the current universal constants do not imply it.

## 5. What orbit entropy sees

Put

```math
\overline w_\epsilon=\mathbb E_Qw_\epsilon(Q)
=e^{\phi_A(t)+\phi_D(t)}(1+\epsilon u_Au_D),          \tag{5.1}
```

where `u_A=E sinh(tH_A)/E cosh(tH_A)` and similarly for `D`.  Normalize

```math
a_\epsilon={w_\epsilon\over\overline w_\epsilon},
\qquad
b_B={k_B\over e^{\psi_B(t)}},
\qquad
r_{\epsilon,B}=a_\epsilon*b_B.                      \tag{5.2}
```

All three functions have uniform group mean one.  Let `Pi_(epsilon,B)` be
the probability law with density `r_(epsilon,B)` relative to uniform
switches.  Averaging the *logarithm* of (2.3) yields the exact quenched
identity

```math
\boxed{
\mathbb E_gL_{\epsilon,g}(B)
=\phi_A(t)+\phi_D(t)+\log(1+\epsilon u_Au_D)
 +\psi_B(t)-D(U_\mathcal G\Vert\Pi_{\epsilon,B}).}  \tag{5.3}
```

Hence

```math
\boxed{
E_{m,n}(\beta)
\le\min_{\epsilon,B}\{
\psi_B-\Delta_A-\Delta_D+\log(1+\epsilon u_Au_D)
-D(U\Vert\Pi_{\epsilon,B})\}.}                      \tag{5.4}
```

This strictly strengthens arithmetic orbit averaging.  The divergence is
the precise gain visible to a geometric-orbit/entropy proof.

In particular, for either chosen orientation and template, the exact
quantitative target is

```math
\boxed{
D(U\Vert\Pi_{\epsilon,B})
\ge\psi_B-\Delta_A-\Delta_D
 +\log(1+\epsilon u_Au_D)-\omega_N
\quad\Longrightarrow\quad
E_{m,n}(\beta)\le\omega_N.}                          \tag{5.4a}
```

Thus when the other displayed terms leave `cN+o(N)`, the orbit law must
carry reverse divergence at least `cN-o(N)`; a qualitative dependence or
collision statement is not enough.

The convolution channel gives two exact data-processing ceilings:

```math
\boxed{
D(U\Vert\Pi_{\epsilon,B})
\le\min\{D(U\Vert a_\epsilon U),
          D(U\Vert b_BU)\}
=\min\{D(U\Vert a_\epsilon U),\psi_B(t)\}.}         \tag{5.5}
```

The second equality uses `E_U log k_B=0`.  The first input divergence is
also explicit:

```math
\boxed{
D(U\Vert a_\epsilon U)
=\phi_A(t)+\phi_D(t)+\log(1+\epsilon u_Au_D)
 -\mathbb E_Q\log w_\epsilon(Q).}                  \tag{5.6}
```

Thus neither orbit cardinality nor generic data processing can force the
needed gain; it is determined by how the two *actual* functions in (5.2)
overlap under translation.

## 6. Why one exceptional switch is beyond fixed orbit averages

For completeness, define the exact negative-temperature soft minimum

```math
\mathsf S_\lambda(B,\epsilon)
=-{1\over\lambda}\log\mathbb E_g
 e^{-\lambda L_{\epsilon,g}(B)},\qquad\lambda>0.     \tag{6.1}
```

Since `|G|=2^(N-1)`, the elementary log-sum bounds are

```math
\boxed{
\min_gL_{\epsilon,g}(B)
\le\mathsf S_\lambda(B,\epsilon)
\le\min_gL_{\epsilon,g}(B)+{(N-1)\log2\over\lambda}.} \tag{6.2}
```

Therefore a fixed `lambda` has a worst-case `Theta(N)` resolution on this
orbit.  Resolving an `o(N)` exceptional-switch gain uniformly by this route
requires `lambda->infinity`; at that point (6.1) approaches the restricted
bridge minimization itself.  Equation (5.3) is the `lambda->0` geometric
endpoint and can entirely miss a single exponentially rare aligned switch.

This is a theorem about the resolution of the entropy surrogate, not a
claim that the actual optimizing-child orbit contains such a switch.

## 7. Verdict

The exceptional aligned-switch problem is now exactly the positive
translation problem (4.4).  The following are rigorously excluded:

1. using switch entropy to cancel the pure bridge pressure;
2. obtaining the geometric-orbit gain from data processing alone; and
3. uniformly resolving an `o(N)` minimum with a fixed negative-moment
   parameter.

What remains is narrower but still open: find an actual-child-dependent
template `B` for which a group translate anticorrelates `k_B` with the
positive internal excess `w_epsilon-1` strongly enough to prove (4.4), or
prove that no such template/translate exists.  The known exponential
conditional spread of actual children does not by itself decide this
translation geometry.
