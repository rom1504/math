# Gauge refresh versus semantic phase: an expander no-go theorem

Status: theorem draft for independent audit.  The first proposition gives a
small-description dense-sign realization of phase refresh which is not close,
in fixed coordinates, to one operator.  It also exposes why that example is
semantically trivial.  The main theorem shows that a nonconstant response
phase cannot admit the same kind of bounded-state expander refresh at small
transfer defect.

## 1. Gauge orbits give exact low-description refresh

Let `S_r` be self-adjoint operators on finite probability spaces and suppose
`L_r` is both an `L^2` and an `L^infinity` contraction with

```math
\|S_r-L_r^*S_(r+1)L_r\|_(2->2)\le eta_r.             \tag{ER.1}
```

For a finite phase set `G`, let `V_(r,g)` be signed coordinate permutations
and put

```math
T_(r,g)=V_(r,g)^*S_rV_(r,g).                          \tag{ER.2}
```

### Proposition ER.1 (gauge-orbit refresh)

For **every** Markov kernel `P_r` on `G`, the branch maps

```math
U_(r,g,h)=V_(r+1,h)^*L_rV_(r,g)                      \tag{ER.3}
```

give an operator-certified phase transfer with defect at most `eta_r` and
phase marginal `P_r`.  In particular one may choose the rank-one refresh
`P_r(g,h)=nu(h)` for any full-support law `nu`, independently of the phase.

Indeed, every branch has the same pullback

```math
U_(r,g,h)^*T_(r+1,h)U_(r,g,h)
=V_(r,g)^*L_r^*S_(r+1)L_rV_(r,g),                    \tag{ER.4}
```

so averaging does nothing and (ER.1) proves the claim.  Signed coordinate
permutations preserve the Boolean response set, hence all phase responses in
(ER.2) are already identical.  Thus this is a genuine low-description
realization of Theorem 31.1, but it refreshes only **gauge**, not semantic
response information.

This distinction is not the same as operator closeness in fixed coordinates.
For the order-four Walsh matrix `H`, let

```math
D=diag(1,-1,-1,1),\qquad H_r=H^(tensor r).            \tag{ER.5}
```

The two hollow normalized operators obtained from `H_r` and
`(D tensor I)H_r(D tensor I)` have operator distance exactly two for every
`r`: diagonal hollowing cancels in their difference and

```math
\left\|{H-DHD\over2}\right\|_(2->2)=2.              \tag{ER.6}
```

Consequently no single operator is within distance less than one of both
phases.  Nevertheless the phase needs one bit, the switch `D tensor I` has a
constant-size description, and (ER.3) gives arbitrary exact refresh (up to
the vanishing diagonal recovery defect).  The obstruction below therefore
cannot be an operator-diameter argument; it uses a response phase that is
actually observable.

## 2. Expander refresh forces either toll or memory

Write `Pi f=int f dpi` for averaging under a full-support probability `pi` on
a finite set `X`.  A kernel `P` is **`rho`-scrambling in `L^2(pi)`** if it
preserves `pi` and

```math
\|P-Pi\|_(L^2(pi)->L^2(pi))\le rho<1.                \tag{ER.7}
```

This includes reversible expanders with absolute second eigenvalue at most
`rho`, but reversibility is not needed.

### Theorem ER.2 (mixing-time lower bound for semantic refresh)

Let `P_j`, `j=r,...,r+t-1`, preserve the same law `pi` and satisfy (ER.7),
where `0<=rho<1`.  Let `g:X->[0,B]`, and suppose functions `f_j:X->R`
obey (ER.8) with `epsilon_j,omega_j>=0`.

```math
\|f_j-g\|_infinity\le omega_j,
\qquad f_j\le P_jf_(j+1)+epsilon_j.                  \tag{ER.8}
```

For every `x in X`, put

```math
D_x=g(x)-int g dpi,
\qquad delta_j=epsilon_j+omega_j+omega_(j+1).         \tag{ER.9}
```

Then for every `t>=1`,

```math
D_x\le {B rho^t\over sqrt(pi(x))}
       +sum_(j=r)^(r+t-1)delta_j.                    \tag{ER.10}
```

In particular, if `D_x=D>0`, `pi(x)>=kappa/S`, where `S=|X|`, `kappa>0`,
and `delta_j<=delta`, then for `0<rho<1` define

```math
t_*=
max\left\{1,
\left\lceil{\log(2B\sqrt{S/kappa}/D)\over
                   \log(1/rho)}\right\rceil\right\}. \tag{ER.11}
```

Then

```math
delta\ge {D\over2t_*}.                               \tag{ER.12}
```

Equivalently, with `lambda=log(1/rho)`, any `0<delta<D/2` certificate satisfies

```math
\log S
\ge \log kappa+{lambda D\over delta}
     -2lambda-2\log(2B/D).                           \tag{ER.13}
```

Thus a response phase with a fixed positive excess above its stationary
average has a quantitative dichotomy:

* a bounded-state uniformly scrambling quotient pays a fixed transfer toll;
* if the toll is `delta`, its phase description needs
  `Omega(lambda D/delta)` bits.

This is stronger than applying a one-step Doeblin coefficient when the
certificate is presented locally by a sparse expander: no minorization need
be listed.  The bound derives the relevant refresh window from the local
mixing certificate.

When `rho=0`, the one-step form of (ER.10) applies directly; the logarithmic
rearrangement (ER.11)--(ER.13) is not used.  When `delta=0<D/2`, sufficiently
long instances of (ER.10) are already contradictory.

#### Proof

The two estimates in (ER.8) imply

```math
g\le P_jg+delta_j.                                    \tag{ER.14}
```

Iterating positive kernels gives

```math
g\le P_rP_(r+1)...P_(r+t-1)g
       +sum_(j=r)^(r+t-1)delta_j.                    \tag{ER.15}
```

Each `P_j` fixes constants and contracts the mean-zero subspace by `rho`, so
their product contracts it by `rho^t`.  Evaluation at `x` is the inner
product in `L^2(pi)` against `1_x/pi(x)`.  Cauchy--Schwarz therefore gives

```math
|(P_r...P_(r+t-1)g)(x)-int g dpi|
\le rho^t\|g-int g dpi\|_(2,pi)
       \sqrt{1/pi(x)-1}
\le {B rho^t\over\sqrt{pi(x)}}.                     \tag{ER.16}
```

Equations (ER.15)--(ER.16) prove (ER.10).  The choice (ER.11) makes its first
term at most `D/2`, yielding (ER.12).  Since

```math
t_*\le 1+{\log(2B/D)+\tfrac12\log(S/kappa)\over
                 lambda},                            \tag{ER.17}
```

(ER.12) rearranges to (ER.13). `square`

The same proof works with a nonconstant sequence of scrambling factors by
replacing `rho^t` with `prod_j rho_j`.  More operationally, (ER.10) says that
**every** window which forgets its starting phase must accumulate at least
the response excess left after the mixing error.  It does not assume a
particular parametrization of the phase.

## 3. Mandatory Walsh test

Let `Phi(t)=2L(t)` be the order-four Walsh prefix response from Theorem 30.1.
The proved values and cap are

```math
Phi(1)=Phi(4)=1,
\qquad Phi(3)\ge c_*={89\over48\sqrt3},
\qquad 0\le Phi(t)\le2.                              \tag{ER.18}
```

Take any finite phase sample `X subset [1,4]` containing `1,3,4`.  Give `1`
and `4` mass `99/200` each and distribute the remaining mass `1/100` over
the other states with every atom at least `kappa/S` (the uniform distribution
on the remaining states permits `kappa=1/100`).  Then

```math
int Phi dpi\le1.01,
\qquad D_3:=Phi(3)-int Phi dpi
\ge D_*:={89\over48\sqrt3}-1.01
=0.06050362412... .                                  \tag{ER.19}
```

Suppose a finite-state pullback presentation on this sample has phase kernels
which preserve `pi`, have common `L^2(pi)` scrambling factor at most `rho`,
and whose recovered responses converge uniformly to `Phi|_X`.  Applying
Theorem ER.2 at `x=3`, with `B=2`, shows that its combined one-step operator
transfer and recovery toll `delta` must satisfy

```math
\log S
\ge \log kappa+{\log(1/rho)D_*\over delta}
 -2\log(1/rho)-2\log(4/D_*).                         \tag{ER.20}
```

For a half-scrambling presentation (`rho=1/2`) with
`delta<=C/sqrt(N)`, this becomes

```math
\log_2S\ge {D_*\over C}\sqrt N-O(1)
= {0.06050362412...\over C}\sqrt N-O(1).             \tag{ER.21}
```

Hence bounded, polynomial-state, and even `exp(o(sqrt(N)))` phase quotients
cannot erase the Walsh semantic phase at the natural `N^(-1/2)` operator
defect scale through a uniformly scrambling local kernel.  A raw operator
certificate is not the only possible source of quadratic description, but
an expander presentation does not magically reduce it to bounded state:
either it pays a fixed defect, weakens its mixing with scale, or stores at
least `Omega(sqrt(N))` phase bits.

This conclusion does **not** prove that every possible Walsh synchronization
certificate is large.  It isolates the precise escape routes: nonstationary
phase laws with no common reference measure, vanishing scrambling rate,
nonuniform recovery, or a transfer mechanism not expressible by the
one-sided Boolean pullback inequality.  Those are mathematical distinctions,
not variations hidden by the state count.

## 4. What is new and what is classical

The `L^2` mixing estimate is classical Markov-chain theory, and the gauge
orbit calculation is elementary.  The theorem-level contribution here is
their use as an information lower bound for operator-certified extremal
recovery: the observable excess `D`, local scrambling rate, transfer toll,
and number of reusable phase states enter one inequality.  It distinguishes
two families which fixed-coordinate operator diameter cannot distinguish:
Walsh gauge phases refresh with one bit, whereas the nonconstant Walsh scale
phase requires exponentially many states in inverse defect under the same
mixing architecture.
