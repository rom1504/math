# Finite square carriers: the exact channel theorem and its entropy ceiling

## 1. Question and verdict

The exact square-polynomial recovery theorem reduces every optimal row factor
to a nonnegative density of fixed literal Walsh degree and fixed `L^2` norm.
A natural next proposal is to put a finite net on that row carrier, average the
net points to obtain a common base measure, and use a singular-value or
strong-data-processing inequality for the resulting channel.  This audit
separates the part of that proposal which is true from the part which fails.

There are three conclusions.

1. Every finite carrier has an exact common-base channel.  Its contraction is
   exactly the Poincare gap of a weighted overlap graph, and tensor products
   attenuate each row-order channel by the corresponding singular values.
2. The full fixed-degree carrier has no uniform strict contraction: it
   contains disjoint square-polynomial densities.  Any accurate net contains
   a two-point subchannel whose contraction tends to one.
3. More decisively, even contraction tending to zero does not control the
   uniform minimum over the row-factor words.  A weak-bias family has top
   singular value `c/sqrt(n)` and bounded row divergence, yet a two-state
   rank-one pressure has an extensive exceptional response.  The same
   exceptional response survives, with fixed loss, inside the exact
   fixed-degree square carrier supplied by Theorem 37.39.

Thus a finite net plus a natural-channel SDPI is not a closure theorem.  It
controls typical `L^2` response, while the optimization asks for a rare
uniform extremum.  An additional optimizer-specific rare-event theorem is
necessary.

## 2. The finite common-base channel

Let `(Omega,U)` be finite and let

```math
\mathcal Q=\{q_a:a\in\mathcal A\},\qquad
q_a\ge0,\qquad E_Uq_a=1.
```

Fix a full-support prior `pi` on `A` and put

```math
\mu=\sum_a\pi_aq_aU.
```

The forward response operator and its adjoint are

```math
(SF)(a)=E_U[q_aF],
\qquad
(Th)(b)=\sum_a{\pi_aq_a(b)\over\mu(b)}h(a).
\tag{FC.1}
```

Thus `T h=E[h(A)|B]` under the joint law
`P(A=a,B=b)=pi_a q_a(b)U(b)`, and `S=T^*`.

Define symmetric overlap weights

```math
w_{ab}=\pi_a\pi_bE_U{q_aq_b\over\sum_c\pi_cq_c}
```

and the overlap gap

```math
\gamma(\mathcal Q,\pi)=
\inf_{E_\pi h=0}
{\frac12\sum_{a,b}w_{ab}(h_a-h_b)^2
 \over E_\pi h^2}.
\tag{FC.2}
```

### Theorem FC.1 (exact overlap-graph channel theorem)

On the mean-zero spaces,

```math
\boxed{\|S\|^2=\|T\|^2=1-\gamma(\mathcal Q,\pi).}
\tag{FC.3}
```

For `m` independent rows, let

```math
R(a_1,\ldots,a_m)
=E_{\otimes_iq_{a_i}}F.
```

If `F_S` denotes its orthogonal input component using a nonconstant
singular direction on exactly the row set `S`, then

```math
\boxed{
\|R_S\|_{L^2(\pi^{\otimes m})}
\le(1-\gamma)^{|S|/2}
\|F_S\|_{L^2(\mu^{\otimes m})}.}
\tag{FC.4}
```

In particular,

```math
\operatorname {Var}_{\pi^{\otimes m}}R
\le(1-\gamma)\operatorname {Var}_{\mu^{\otimes m}}F.
\tag{FC.5}
```

*Proof.*  Conditional variance gives, exactly,

```math
E_\pi h^2-\|Th\|_{L^2(\mu)}^2
=E\operatorname {Var}(h(A)\mid B)
=\frac12\sum_{a,b}w_{ab}(h_a-h_b)^2.
```

Taking the infimum proves (FC.3).  Singular-value decomposition followed by
tensorization proves (FC.4); summing the nonconstant components proves
(FC.5). `square`

The graph in (FC.2), rather than the literal polynomial degree, decides
strict contraction.  In particular, `gamma>0` exactly when the carrier
overlap graph is connected after deleting zero-prior vertices.

There is also an exact entropy price for a uniform extremum.  Suppose

```math
\log E_{\mu^{\otimes m}}
 e^{s(F-EF)}\le {s^2\sigma^2\over2}
\qquad(s\in\mathbb R).
\tag{FC.6}
```

Conditional Jensen gives the same moment bound for the random response
`R(A_1,...,A_m)`.  Hence, if `pi_min=min_a pi_a`,

```math
\boxed{
\operatorname {range}R
\le2\sqrt{2m\log(1/\pi_{\min})\,\sigma^2}.}
\tag{FC.7}
```

For a uniform alphabet of size `K`, the entropy charge is `m log K`.
Notice that the contraction coefficient does not appear in (FC.7): an
`L^2` singular-value estimate controls typical responses, not the least
likely word of mass `K^{-m}`.

## 3. The finite-net state size

Let

```math
\mathcal S_{n,d,K}=\{q\ge0:E_Uq=1,
 \deg q\le2d,\ \|q\|_2\le K\}.
```

This contains the exact square carrier of Theorem 37.39.  It lies in a
Euclidean coefficient space of dimension

```math
D_{n,2d}-1,
\qquad D_{n,2d}=\sum_{j=0}^{2d}{n\choose j}.
```

A standard volumetric argument gives an internal `L^2` epsilon-net with

```math
\boxed{
\log |\mathcal N_\epsilon|
\le D_{n,2d}\log(1+2K/\epsilon).}
\tag{FC.8}
```

Replacing a row density by its net representative changes a physical bridge
pressure by at most `beta epsilon`; on the fixed `L^2` ball it changes row
entropy by at most the uniform modulus `omega_K(epsilon)`.  Thus the net is a
valid finite discretization with total product-objective error

```math
m\{\beta\epsilon+\lambda^{-1}\omega_K(\epsilon)\}.
\tag{FC.9}
```

It is not yet a compression: (FC.8) is the size of a coefficient net, and
the product query alphabet has the corresponding `m`-fold entropy.

## 4. First obstruction: exact and near-exact disconnected fibres

Let `chi_S` be a nonconstant Walsh character with `|S|<=d`, and put

```math
q_+=1+\chi_S,
\qquad q_-=1-\chi_S.
\tag{FC.10}
```

These are exact members of the square carrier because

```math
q_\pm=\left({1\pm\chi_S\over\sqrt2}\right)^2.
```

They have `L^2` norm `sqrt(2)`, identical entropy `log 2`, and disjoint
supports.  For the uniform two-point prior, the channel reveals the query
sign exactly, so

```math
\gamma=0,\qquad \|S\|=\|T\|=1.                 \tag{FC.11}
```

This is robust under net approximation.  If probability densities `p_+`
and `p_-` satisfy

```math
\|p_\pm-q_\pm\|_2\le\epsilon,
```

then the natural binary-channel coefficient `kappa` obeys

```math
1-\kappa^2
=2E_U{p_+p_-\over p_++p_-}
\le2E_U\min(p_+,p_-)
\le2\epsilon.
\tag{FC.12}
```

The last inequality follows from
`TV(p_+,p_-)>=1-epsilon`.  Consequently every sufficiently accurate net of
the full carrier contains a two-point subchannel with

```math
\boxed{\kappa^2\ge1-2\epsilon.}                 \tag{FC.13}
```

No prior-uniform SDPI for the complete carrier can therefore have a fixed
gap.  Excluding (FC.10) would require an additional diffuse/optimizer-specific
hypothesis not present in the square-carrier theorem.

## 5. A weak-coordinate exponential packing

The disconnected example has strong individual coordinates.  The same net
problem persists inside the weak-coordinate, bounded-divergence class from
which the square carrier was obtained.

Fix `c>0`, set `a=c/sqrt(n)<1/2`, and, for
`y in {+-1}^n`, define

```math
f_y(b)=\prod_{j=1}^n(1+a y_jb_j).                \tag{FC.14}
```

These are product probability densities and

```math
D_2(f_yU\Vert U)=n\log(1+c^2/n)\le c^2,
\qquad
\operatorname {osc}_{b_j}\log f_y
=\log{1+a\over1-a}=O_c(n^{-1/2}).                \tag{FC.15}
```

If `h=d_H(y,y')`, then

```math
\|f_y-f_{y'}\|_2^2
=2(1+a^2)^n
 \left[1-\left({1-a^2\over1+a^2}\right)^h\right].
\tag{FC.16}
```

Choose a binary code `C_n` of relative distance `delta in (0,1/2)` and
cardinality `exp(r_delta n)`.  Equations (FC.15)--(FC.16) give a fixed
separation

```math
\|f_y-f_{y'}\|_2
\ge\Delta_{c,\delta}
:=\sqrt{2(1-e^{-2c^2\delta})}
\qquad(y\ne y',\ y,y'\in C_n).                 \tag{FC.17}
```

Apply Theorem 37.39 with the common constants in (FC.15), and fix `d` large
enough that its approximation error is below `Delta_(c,delta)/4`.  The
resulting exact degree-`2d` square carriers `q_(y,d)` form an
`exp(r_delta n)` packing at a fixed positive `L^2` scale.  Therefore

```math
\boxed{
\log N(\mathcal S_{n,d,K_1},L^2,\epsilon)
\ge r_\delta n
}
\tag{FC.18}
```

for one fixed `d,K_1,epsilon>0` and all large `n`.

Thus even the weak-coordinate source class forces a per-row alphabet with
`Omega(n)` bits under naive `L^2` discretization.  This does not by itself
lower-bound the smaller response quotient of an actual child; it shows that
coefficient-net entropy cannot be substituted for such a quotient.

## 6. Vanishing singular value but an extensive rare response

The family (FC.14) gives a sharper separation between typical channel
compression and uniform query minima.  Take the uniform prior on all
`2^n` directions `y`.  Its common base is exactly `U_n`, and the channel is
the product binary channel

```math
E_{f_y}\chi_S(B)=a^{|S|}\chi_S(y).               \tag{FC.19}
```

Hence its largest nonconstant singular value is

```math
\boxed{\kappa_n=a={c\over\sqrt n}\longrightarrow0.}
\tag{FC.20}
```

Now take `m` rows, let `N=m+n`, fix a sign matrix
`C=(C_ij)` (it may be rank one), put `t=beta/sqrt(N)`, and define

```math
L_C(B)=\log\left(2\cosh\left(t\sum_{i,j}C_{ij}B_{ij}\right)\right).
\tag{FC.21}
```

This is an exact two-state partition pressure.  It is invariant under global
bridge inversion and changing one bridge edge changes it by at most `2t`,
exactly the physical oscillation scale.

For a query word `Y=(y_1,...,y_m)`, set

```math
R(Y)=E_{\otimes_i f_{y_i}}L_C.
```

At the aligned word `Y=C`, Jensen gives

```math
R(C)\ge t a mn.
```

If `Y_0` is a word with `<C,Y_0>=0`, then the argument of the pressure has
mean zero and variance `mn(1-a^2)`, so

```math
R(Y_0)\le\log2+t\sqrt{mn(1-a^2)}.
```

Consequently

```math
\boxed{
\operatorname {range}_Y R(Y)
\ge {\beta c\,m\sqrt n\over\sqrt N}
     -\beta\sqrt{mn/N}-\log2.}                  \tag{FC.22}
```

For a balanced split this is `Theta(N)`, despite (FC.20).  On the other
hand, (FC.5) and bounded differences give

```math
\operatorname {Var}_{Y}R(Y)
\le {c^2\over n}\operatorname {Var}_{U^{\otimes m}}L_C
=O_{\beta,c}(1)                                  \tag{FC.23}
```

on balanced splits.  Thus the natural channel compresses the *typical*
response all the way to bounded variance while missing a linear exceptional
word among `2^{mn}` queries.

The obstruction survives inside a fixed exact square carrier, and the common
base and strict contraction can be preserved.  Apply Theorem 37.39 only to
`f_1` and choose a degree-`d` square approximation `q_1`.  Define the other
carriers equivariantly by

```math
q_y(b)=q_1(y\odot b).                             \tag{FC.24}
```

Then all `q_y` are degree-`2d` square densities with the same fixed `L^2`
norm and the same entropy, and

```math
\|q_y-f_y\|_2\le\epsilon_d
=O_c((d+1)^{-1/3}).                              \tag{FC.25}
```

Under the uniform prior on `y`, their common base is again exactly `U_n`.
The channel is group convolution, so its nonconstant singular values are
the absolute Walsh coefficients `|hat q_1(S)|`.  Parseval and (FC.25) give

```math
\boxed{
\kappa(q)^2=\max_{S\ne\varnothing}|\widehat q_1(S)|^2,
\qquad
\kappa(q)\le a+\epsilon_d.}                     \tag{FC.26}
```

Sequential row replacement changes (FC.21) by at most
`beta m epsilon_d`.  Choose one fixed `d` so this loss is less than one
quarter of the leading coefficient in (FC.22).  Then the response over the
degree-`2d`, fixed-`L^2` densities `q_y` still has linear range.  Since their
entropies agree exactly by (FC.24), the same is true of the
entropy-regularized response.

More strongly, given any declared `kappa_0>0`, first choose this fixed `d`
large enough that `epsilon_d<kappa_0/2` and then take `n` large enough that
`a<kappa_0/2`.  Equations (FC.22), (FC.25), and (FC.26) produce an exact
fixed-degree square-carrier channel with

```math
\boxed{
\kappa(q)\le\kappa_0
\quad\hbox{but}\quad
\operatorname {range}R_q\ge c_{\beta,c,\theta}N.}
\tag{FC.27}
```

Thus even an arbitrarily strong, dimension-independent SDPI on the exact
recovered carrier does not control its uniform `m`-fold extremum.

This is not an actual-minimizing-child counterexample: (FC.21) is a
two-state rank-one pressure and the factors (FC.14) are admissible weak-row
density models, not proved optimal factors of a contracted-temperature
child.  Its scope is exact and important: bounded row complexity, fixed
literal degree, a common base measure, and even a vanishing singular value
do not imply uniform product-phase closure.

## 7. Consequence for the smallest missing lemma

The finite-channel theorem is useful bookkeeping, but it does not close the
bounded-row-degree cross-row problem.  A valid next theorem must use an
optimizer-specific fact which excludes the rare alignment in (FC.22).  One
precise form would be a metric-entropy-sensitive tail bound for the *actual*
child pressure strong enough that

```math
\log N_{\rm eff}(\epsilon)\,\sigma_N^2(\epsilon)=o(N^2)
```

at the response scale being optimized, together with chaining across
`epsilon`; equivalently, an actual-child quotient whose effective query
entropy is `o(N)` per product rather than the `Theta(mn)` entropy in
(FC.22).  A singular-value estimate, even one with `kappa_n->0`, is not a
replacement for this rare-event statement.

This audit therefore rejects the proposed generic finite-net/channel route.
It does not falsify an optimizer-specific response quotient or
superconcentration theorem.
