# Regular-Hadamard recoupling in Walsh coordinates

## Status

This note is a coordinate dictionary and a structural audit.  It proves:

1. all row-sign, anchored-greedy, potential, terminal-margin, unmatched-core,
   and `kappa` quantities for the regular-Hadamard family have exact Walsh
   formulas;
2. the input `X=1` is the self-dual bent sequence `f`, while the originally
   exhibited bad terminal sequence is plateaued rather than bent; and
3. a static simultaneous diffuse-core configuration (whole-block row signs
   plus free-coordinate terminal stability) does **not** force the input or
   terminal to have bent, near-bent, or plateaued Walsh spectrum.  An
   explicit project-scale infinite family has full Walsh support on both
   sides, stays a fixed normalized distance from every flat spectrum, and
   has `kappa/m -> 1/40`.

The last result rejects a bent-specific inverse theorem based only on row
signs and terminal stability.  It does not show that the prescribed
field-initialized dynamics reaches these terminals, much less that they have
a large basin.  The counterfamily is tensor-structured, not a claim about a
generic random Walsh spectrum.

All energies use doubled normalization.  No asymptotic bound for the
original signing problem is improved here.

## 1. Walsh conventions and the switching dictionary

Fix `k>=2`, put

```math
V=\mathbb F_2^k\times\mathbb F_2^k,
\qquad s=2^k,
\qquad m=|V|=s^2,
```

and write

```math
\beta((a,b),(c,d))=a\mathbin\cdot d+c\mathbin\cdot b.
```

Let

```math
H_{u,v}=(-1)^{\beta(u,v)},
\qquad U={1\over s}H,
\qquad
\widehat h=Uh={1\over s}Hh.                         \tag{1}
```

Thus `U` is a symmetric orthogonal involution.  Inner products below are
counting inner products, so

```math
U^2=I,
\qquad \langle\widehat h,\widehat j\rangle
=\langle h,j\rangle.                                \tag{2}
```

When an unnormalized transform is more convenient, write
`\mathcal W h=Hh=s\widehat h`.

Set

```math
q(a,b)=a\mathbin\cdot b,
\qquad f(u)=(-1)^{q(u)}.
```

Polarization gives the regular symmetric Hadamard matrix

```math
K_{u,v}=f(u+v)=f(u)H_{u,v}f(v),\qquad K=D_fHD_f,
\qquad C=K-I.                                        \tag{3}
```

The quadratic Gauss-sum identity is

```math
Hf=sf,
\qquad\boxed{\widehat f=f}.                          \tag{4}
```

In the normalized convention, (4) says that `f` is self-dual bent.  In the
unnormalized convention it is `\mathcal Wf=sf`.

For a Boolean input `X`, introduce its bent-conjugated spin

```math
\boxed{g=fX}.                                        \tag{5}
```

Then every occurrence of `f` disappears from the switched signing:

```math
D_XCD_X=D_gHD_g-I.                                   \tag{6}
```

Moreover,

```math
(CX)(u)=f(u)\bigl(s\widehat g(u)-g(u)\bigr),
```

and the switched row field is

```math
\boxed{
\ell_g(u)=X(u)(CX)(u)
=s\,g(u)\widehat g(u)-1.}                           \tag{7}
```

With the zero-field convention used by the greedy audit,

```math
I_g=\{u:s g(u)\widehat g(u)-1\ge0\},
\qquad
J_g=\{u:s g(u)\widehat g(u)-1<0\}.                  \tag{8}
```

Thus the row-sign law is a coordinatewise threshold on the signed Walsh
response `g\widehat g`; no matrix multiplication remains hidden.

For `E\subset V`, abbreviate `g_E=g1_E`.  Its witnessed principal energy is

```math
\boxed{
\mathcal E_g(E)
=s\langle g_E,\widehat {g_E}\rangle-|E|.}            \tag{9}
```

For disjoint fixed and free shores `U_0,F`, the cross field on `F` is

```math
\boxed{
h_{U_0\to F}(u)=s\,g(u)\widehat {g_{U_0}}(u)
\quad(u\in F).}                                      \tag{10}
```

In particular,

```math
P=s\langle g_I,\widehat {g_I}\rangle-|I|,
\qquad
R=s\langle g_J,\widehat {g_J}\rangle-|J|,           \tag{11}
```

and the cross energy is
`s\langle g_J,\widehat {g_I}\rangle`.

## 2. Anchored recoupling and greedy ascent

This section treats both hard branches at once.  Let `U_0` be the fixed
shore, `F=V\setminus U_0` the free shore, and

```math
\sigma=\operatorname{sgn}\mathcal E_g(U_0).
```

In the hard orientation, `(U_0,F,\sigma)` is either `(I,J,+1)` or
`(J,I,-1)`.  Hence the free-shore weights

```math
\boxed{
a_u=-\sigma\ell_g(u)
=-\sigma\bigl(sg(u)\widehat g(u)-1\bigr)\ge0
\quad(u\in F)}                                      \tag{12}
```

are nonnegative (strict except for zero-field ties).

Let `r\in\{\pm1\}^F` be the free relative spin and let
`t\in\{\pm1\}` be the collapsed anchor.  Define

```math
z=g_Fr,
\qquad
\phi=tg_{U_0}+z.                                     \tag{13}
```

The complete sign-specific augmented objective is

```math
\boxed{
\sigma\left(r^T(D_gHD_g-I)[F]r+2t h^Tr\right)
=\sigma\left(
s\langle\phi,\widehat\phi\rangle
-s\langle g_{U_0},\widehat {g_{U_0}}\rangle-|F|
\right).}                                           \tag{14}
```

This identity retains the fixed/free cancellation before any absolute
value.  The field-aligned initialization becomes

```math
z^{(0)}(u)=\operatorname{sign}
\bigl(\sigma\widehat {g_{U_0}}(u)\bigr),             \tag{15}
```

with `z^(0)(u)=g(u)` under the prescribed zero tie.

For a free coordinate `u`, the exact flip gain is

```math
\boxed{
-4\sigma\bigl(s z(u)\widehat\phi(u)-1\bigr).}       \tag{16}
```

The anchor flip gain is

```math
\boxed{
-4\sigma t s\langle z,\widehat {g_{U_0}}\rangle.}   \tag{17}
```

Consequently a terminal state is characterized exactly by

```math
\sigma\bigl(s z(u)\widehat\phi(u)-1\bigr)\ge0
\quad(u\in F),
\qquad
\sigma t\langle z,\widehat {g_{U_0}}\rangle\ge0.    \tag{18}
```

At such a state the augmented terminal certificate is

```math
\boxed{
G_\sigma
=\sum_{u\in F}|s\widehat\phi(u)-z(u)|
+s\left|\langle z,\widehat {g_{U_0}}\rangle\right|.} \tag{19}
```

Equations (14)--(19) are the Walsh form of the anchored matrix, its cross
field, every greedy gain, and the terminal local-field identity.

## 3. Potential, terminal margins, unmatched mass, and `kappa`

Absorb a possible anchor sign into `r`, so that `t=1`.  There is then a
set `S\subset F` for which

```math
r=1-2 1_S,
\qquad
\phi=g-2g_S.                                         \tag{20}
```

The weighted-set potential has the particularly short Walsh form

```math
\boxed{
\Phi_{F,\sigma}(S)
=-\sigma s\langle g_S,\widehat {g-g_S}\rangle.}      \tag{21}
```

Indeed, the vertex part contributes
`-\sigma s\langle g_S,\widehat g\rangle+\sigma|S|`, while twice the effective
internal-edge sum contributes
`\sigma(s\langle g_S,\widehat g_S\rangle-|S|)`; the diagonal terms cancel exactly.

Put `T=F\setminus S`.  The removal margins on `S` and addition margins on
`T` are

```math
\boxed{
\begin{aligned}
p_u&=-\sigma\bigl(1+s g(u)\widehat\phi(u)\bigr)
&& (u\in S),\\
q_u&= \sigma\bigl(s g(u)\widehat\phi(u)-1\bigr)
&& (u\in T).
\end{aligned}}                                      \tag{22}
```

Singleton stability is exactly `p_u,q_u\ge0`.  With

```math
A_F=\sum_{u\in F}a_u,
```

the unmatched mass is

```math
\boxed{
u_{F,\sigma}
=a(T)-p(S)
=A_F-2\Phi_{F,\sigma}(S)
=A_F+2\sigma s
 \langle g_S,\widehat {g-g_S}\rangle.}              \tag{23}
```

The branch defect is

```math
\boxed{\delta_{F,\sigma}=2[u_{F,\sigma}]_+.}         \tag{24}
```

Finally, if `a_{(1)}\ge a_{(2)}\ge\cdots` are the weights of the points in `T` in
decreasing order, then the exact unmatched-core number is

```math
\boxed{
\kappa_{F,\sigma}
=\min\left\{j:\sum_{i=1}^j a_{(i)}
\ge u_{F,\sigma}\right\},}                           \tag{25}
```

with value zero when `u\le0`.  The joint quantity in the recoupling theorem
is `\kappa_*=\min(\kappa_{J,+},\kappa_{I,-})` in the hard branch.

Thus `\kappa` is also completely determined by the threshold deficits
`1-\sigma s g\widehat g` and one masked Walsh correlation.

## 4. The two-block signing and its bridge

For the full obstruction, let

```math
A=\begin{pmatrix}C&B\\B^T&-C\end{pmatrix},
\qquad
\mathcal B=D_fBD_f,
\qquad
\mathcal R={1\over s}\mathcal B.                    \tag{26}
```

For a full input `(X_+,X_-)`, put `g_+=fX_+`, `g_-=fX_-`, and introduce
the symmetric operator

```math
\mathcal T=
\begin{pmatrix}U&\mathcal R\\
\mathcal R^T&-U\end{pmatrix},
\qquad
J_0=\begin{pmatrix}I&0\\0&-I\end{pmatrix}.          \tag{27}
```

Then the full switched signing is exactly

```math
\boxed{
D_XAD_X=D_g(s\mathcal T-J_0)D_g.}                    \tag{28}
```

Write an index as `alpha=(epsilon,u)`, with
`tau_alpha=+1` in the first block and `tau_alpha=-1` in the second.  The
complete row-field formula is

```math
\boxed{
\ell_\alpha=s g_\alpha(\mathcal Tg)_\alpha
-\tau_\alpha.}                                      \tag{29}
```

Explicitly,

```math
\begin{aligned}
\ell_+(u)&=s g_+(u)\widehat g_+(u)-1
+g_+(u)(\mathcal B g_-)(u),\\
\ell_-(u)&=-s g_-(u)\widehat g_-(u)+1
+g_-(u)(\mathcal B^T g_+)(u).
\end{aligned}                                       \tag{30}
```

All preceding formulas persist after the following replacements:

```math
\begin{array}{c|c}
\text{one-block expression}&\text{two-block expression}\\ \hline
s\langle h,\widehat j\rangle&s\langle h,\mathcal Tj\rangle\\
-|E|&-\sum_{\alpha\in E}\tau_\alpha\\
s z(u)\widehat\phi(u)-1&
s z_\alpha(\mathcal T\phi)_\alpha-\tau_\alpha.
\end{array}                                         \tag{31}
```

For clarity, the resulting formulas are

```math
h_{U_0\to F}(\alpha)
=s g_\alpha(\mathcal T g_{U_0})_\alpha,              \tag{32}
```

```math
\sigma\left(
s\langle\phi,\mathcal T\phi\rangle
-s\langle g_{U_0},\mathcal T g_{U_0}\rangle
-\sum_{\alpha\in F}\tau_\alpha\right),             \tag{33}
```

```math
\begin{aligned}
\text{free gain}&=-4\sigma
 [s z_\alpha(\mathcal T\phi)_\alpha-\tau_\alpha],\\
\text{anchor gain}&=-4\sigma t s
 \langle z_F,\mathcal T g_{U_0}\rangle,              \tag{34}
\end{aligned}
```

and

```math
\boxed{
\Phi_{F,\sigma}(S)
=-\sigma s\langle g_S,\mathcal T(g-g_S)\rangle.}     \tag{35}
```

The weights and margins are

```math
\begin{aligned}
a_\alpha&=-\sigma
 [s g_\alpha(\mathcal Tg)_\alpha-\tau_\alpha],\\
p_\alpha&=-\sigma
 [s g_\alpha(\mathcal T\phi)_\alpha+\tau_\alpha]
&& (\alpha\in S),\\
q_\alpha&= \sigma
 [s g_\alpha(\mathcal T\phi)_\alpha-\tau_\alpha]
&& (\alpha\in T).                                   \tag{36}
\end{aligned}
```

Equations (23)--(25) then hold verbatim with (35)--(36).  This is the
requested bridge-aware Walsh dictionary: the only datum beyond `g_+`,
`\widehat g_+`, `g_-`, and `\widehat g_-` is the explicitly conjugated bridge
operator `\mathcal R`.

## 5. The known bad input and bad terminal

For `X_+=X_-=1`,

```math
g_+=g_-=f,
\qquad \widehat f=f.                                 \tag{37}
```

The bridge in the obstruction has `B1=B^T1=0`, and hence

```math
\mathcal Bf=\mathcal B^Tf=0.
```

Therefore

```math
\mathcal T(f,f)=(f,-f),
\qquad
\ell_+=s-1,
\qquad
\ell_-=1-s.                                          \tag{38}
```

The two shores are exactly the two blocks,

```math
P=m(s-1),
\qquad R=-m(s-1),                                    \tag{39}
```

and both anchored cross fields vanish.

Let `S_0` be the codimension-two subspace from the original obstruction,
put

```math
r_0=1-2 1_{S_0},
\qquad z_0=fr_0.                                     \tag{40}
```

The identity `K1_(S_0)=(s/2)r_0` is equivalent to

```math
\boxed{
\widehat {z_0}=f-z_0=2f1_{S_0}
=-2z_0 1_{S_0}.}                                    \tag{41}
```

Thus `z_0` is **not bent**: its normalized Walsh spectrum is zero off
`S_0` and has magnitude two on `S_0`.  It is a two-plateaued sequence.
The relative spin itself has

```math
\boxed{
\widehat {r_0}
=s\,1_{\{0\}}-{s\over2}1_{S_0^\perp},}              \tag{42}
```

so its Walsh support has size four and constant magnitude `s/2`; it too is
plateaued, not bent (apart from irrelevant smallest-order coincidences).

Equation (41) gives

```math
z_0(u)\widehat {z_0}(u)
=\begin{cases}0,&u\notin S_0,\\-2,&u\in S_0.
\end{cases}                                         \tag{43}
```

For the `-C` target, the strict local margins are therefore `1` off `S_0`
and `2s+1` on `S_0`.  Also,

```math
z_0^THz_0=-{1\over2}ms,
\qquad
-r_0^TCr_0=m\left({s\over2}+1\right).               \tag{44}
```

The unmatched-core quantities are now immediate from (12), (22), and
(25):

```math
\begin{aligned}
a_u&=s-1 &&(u\in F),\\
|S_0|&=m/4,
&p_u&=2s+1 &&(u\in S_0),\\
|T|&=3m/4,
&u&={m(s-4)\over4},\\
\delta&={m(s-4)\over2},
&\kappa&=\left\lceil{m(s-4)\over4(s-1)}\right\rceil.
                                                               \tag{45}
\end{aligned}
```

In particular, `kappa/m -> 1/4`.  Both anchored shores have this same
terminal point.

## 6. A useful flatness identity

There is an exact reason that the bad terminal is anti-flat rather than
approximately bent.  Use the unnormalized transform in this paragraph.
Consider the zero-cross two-block model with the same conjugated Boolean
input `g` in its `C` and `-C` blocks.  Suppose `g` is regularizing and an
anchored `-C` search terminates at an anti-regularizing Boolean vector `z`,
meaning

```math
d_u=g(u)(Hg)(u)\ge0,
\qquad
c_u=-z(u)(Hz)(u)\ge0.                                \tag{46}
```

Assume the row inequalities are strict enough to select the two whole
blocks (for this integral family, `d_u>=2` is sufficient).  Write

```math
E(v)=v^THv,
\qquad
B(g)=ms-E(g),
\qquad
A(z)=ms+E(z).                                        \tag{47}
```

Since `H^2=mI`,

```math
\sum_ud_u^2=\sum_uc_u^2=m^2.
```

Consequently

```math
\boxed{
B(g)={1\over2s}\sum_u(d_u-s)^2,
\qquad
A(z)={1\over2s}\sum_u(c_u-s)^2.}                    \tag{48}
```

The two whole blocks are then the row-sign shores, with witnessed magnitude
`E(g)-m`, and the recoupling defect of this terminal is exactly

```math
\boxed{
\Delta=[E(g)+E(z)-2m]_+
=[A(z)-B(g)-2m]_+.}                                  \tag{49}
```

For the self-dual bent input `g=f`, `d_u=s`, so `B(f)=0`.  A leading
defect therefore forces a leading **increase**, not decrease, in the
terminal Walsh-amplitude nonflatness `A(z)`.

The same conclusion is visible as an `L^1` deficit.  Since `c_u>=0`,

```math
A(z)=ms-\sum_uc_u.                                   \tag{50}
```

If `a=A(z)/(ms)<1`, interpolation gives

```math
\boxed{
\sum_u |(Hz)(u)|^4=\sum_uc_u^4
\ge {m^3\over(1-a)^2}.}                              \tag{51}
```

Indeed, `sum c^2=m^2` and `sum c=ms(1-a)`, and
`(sum c^2)^3 <= (sum c)^2 sum c^4`.  Bent flatness gives exactly `m^3`;
any fixed leading `a>0` forces a strictly larger fourth moment.  This
identity already points away from a near-bent stability theorem.

For (41), `A(z_0)=ms/2`, and (51) is sharp:
`sum |Hz_0|^4=4m^3`.

## 7. Full-support non-plateaued diffuse cores

The preceding warning can be made into an infinite exact counterfamily.
On `V_2=\mathbb F_2^2\times\mathbb F_2^2`, order coordinates by

```math
j=a_0+2a_1+4b_0+8b_1
```

and define the Boolean vector

```math
z_*=(1,1,1,-1,-1,-1,-1,-1,
     -1,-1,-1,-1,1,1,1,1).                           \tag{52}
```

A direct four-bit Walsh sum gives

```math
H_2z_*=(-2,-2,-2,14,2,2,2,2,
        2,2,2,2,-2,-2,-2,-2).                        \tag{53}
```

Hence

```math
z_*(u)(H_2z_*)(u)
=\begin{cases}-14,&u=3,\\-2,&u\ne3,
\end{cases}
\qquad z_*^TH_2z_*=-44.                              \tag{54}
```

This is a strict anti-regularizing vector, but its normalized Walsh
transform has full support and magnitudes

```math
{1\over2}\quad(15\text{ times}),
\qquad {7\over2}\quad(1\text{ time}).               \tag{55}
```

It is neither bent nor plateaued.

For every `k>=2`, split `V_k=V_2\oplus V_{k-2}` and put

```math
z_k=z_*\otimes f_{k-2},                               \tag{56}
```

where the second factor is the self-dual quadratic bent vector (and is
absent when `k=2`).  Tensor factorization of the symplectic Walsh transform
gives

```math
\widehat {z_k}=(U_2z_*)\otimes f_{k-2}.               \tag{57}
```

Thus (55) holds with the same proportions at every order.  Quantitatively,

```math
{1\over m}\sum_u\bigl(|\widehat {z_k}(u)|-1\bigr)^2
={5\over8},
\qquad
{1\over m}\sum_u|\widehat {z_k}(u)|^4
={151\over16}.                                       \tag{58}
```

The family is uniformly separated from all flat Walsh spectra.

Let `m=4^k`, `s=2^k`.  Equations (54)--(57) give

```math
z_k^TH_kz_k=-{11\over16}ms,                           \tag{59}
```

and every coordinate product `z_k(u)(H_kz_k)(u)` is either `-s/2` or
`-7s/2`.  Hence `r_k=f_kz_k` is a strict stable point for the same `-C`
anchored objective.  Since

```math
\langle f_2,z_*\rangle=2,
```

its flip set has size

```math
|S_k|={7m\over16},
\qquad |T_k|={9m\over16}.                             \tag{60}
```

At the self-dual input `g=f`, all weights are again `s-1`.  From (49) or
directly from (23),

```math
\boxed{
\begin{aligned}
\delta_k&=m\left[{5s\over16}-2\right]_+,\\
u_k&=m\left({5s\over32}-1\right),\\
\kappa_k&=\left\lceil
 {m((5/32)s-1)\over s-1}\right\rceil,
\qquad {\kappa_k\over m}\longrightarrow {5\over32}.
\end{aligned}}                                      \tag{61}
```

For `k>=3` this is a positive, linear static diffuse-core configuration.
Because the
two anchored objectives in the full construction are identical and have
zero cross fields, it occurs simultaneously on both shores.

This is an existence statement about strict terminal points.  It does not
claim that the prescribed field-initialized, best-improvement trajectory
reaches `r_k`; in particular it supplies no lower bound on its basin
probability.

One can also rule out the possibility that (55) is merely an accidental
two-amplitude substitute for plateauedness.  On `V_2`, let

```math
w(a,b)=-1\quad\text{if }b=0\text{ and }a\ne0,
\qquad w(a,b)=1\quad\text{otherwise}.                \tag{62}
```

Then `w(u)(H_2w)(u)` has values `10` once, `6` three times, and `2`
twelve times, all positive.  Explicitly,

```math
w=(1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1),
```

and

```math
H_2w=(10,-6,-6,-6,2,2,2,2,2,2,2,2,2,2,2,2),
```

so `w^TH_2w=52`.  Therefore

```math
z_*\otimes w^{\otimes t}                              \tag{63}
```

is strictly anti-regularizing, has full Walsh support, and has exactly
`(t+1)(t+2)` distinct Walsh magnitudes.  Its normalized energy is

```math
{z^THz\over ms}
=-{11\over16}\left({13\over16}\right)^t,             \tag{64}
```

so its unmatched-core density tends `1/2`, not zero.  Large diffuse
terminal cores can therefore have an unbounded number of Walsh levels.

## 8. A non-bent input with a simultaneous diffuse core

The family in (56) disproves terminal near-bentness, but its input is still
the bent vector `f`.  The following construction also makes the input
uniformly non-bent.

Use the positive regularizer `w` from (62) and set, on
`V_k=V_2\oplus V_{k-2}`,

```math
g_k=w\otimes f_{k-2},
\qquad
z_k=z_*\otimes f_{k-2}.                              \tag{65}
```

The input spectrum has full support and normalized magnitudes

```math
{1\over2}\quad(12/16\text{ of the coordinates}),
\quad
{3\over2}\quad(3/16),
\quad
{5\over2}\quad(1/16).                               \tag{66}
```

In particular,

```math
{1\over m}\sum_u\bigl(|\widehat {g_k}(u)|-1\bigr)^2
={3\over8},
\qquad
{1\over m}\sum_u|\widehat {g_k}(u)|^4
={55\over16}.                                       \tag{67}
```

It is neither bent, near-flat, nor plateaued.  Its unnormalized signed row
responses are

```math
g_k(u)(H_kg_k)(u)
\in\left\{{s\over2},{3s\over2},{5s\over2}\right\},  \tag{68}
```

so every row of `C` is favorable and every row of `-C` is unfavorable.
Also,

```math
E(g_k)={13\over16}ms,
\qquad
E(z_k)=-{11\over16}ms.                              \tag{69}
```

To embed this exact pair with zero cross field, let `B_0` be the bridge in
the original obstruction, let the physical input be

```math
x_k=f_kg_k,
```

and define

```math
B_k=D_{x_k}B_0D_{x_k},
\qquad
A_k=\begin{pmatrix}C&B_k\\B_k^T&-C\end{pmatrix}.    \tag{70}
```

Because `B_0 1=B_0^T1=0`,

```math
B_kx_k=B_k^Tx_k=0.                                  \tag{71}
```

The bridge remains a `+-1` matrix with exactly the original operator norm,
so the project-scale estimate for the full signing is unchanged.  At the
full input `(x_k,x_k)`, its two row-sign shores are the two whole blocks and
both anchored cross fields vanish.

The proposed relative terminal spin is

```math
r_k=g_kz_k=(wz_*)\otimes1.
```

Its flipped base coordinates are
`{1,2,4,5,6,7,8,9,10,11}`, so `|S|=10m/16` and
`|T|=6m/16`.  It is strict on every free coordinate by (54); the
collapsed anchor is neutral because its cross field is zero.  Equations
(49) and (69) give

```math
\boxed{
\delta_k=m\left[{s\over8}-2\right]_+,
\qquad
u_k=m\left({s\over16}-1\right).}                    \tag{72}
```

Within `T`, the `m/16` coordinates over base coordinate zero have the
largest weight `5s/2-1`, and that class alone pays `u_k` whenever the latter
is positive.  Hence, for `s>16`,

```math
\boxed{
\kappa_k
=\left\lceil {m(s-16)\over8(5s-2)}\right\rceil,
\qquad
{\kappa_k\over m}\longrightarrow {1\over40}.}       \tag{73}
```

The same free-coordinate terminal can be used in both identical anchored
objectives, so the two static cores are simultaneous.  This proves that even
the **input** of such a project-scale row-sign/terminal-stability
configuration need not be bent, plateaued, or quantitatively near-flat.

The bridge in (70) is explicitly conjugated to annihilate the selected
physical input.  Therefore (65)--(73) refute any uniform implication from
whole-block row signs plus terminal stability to Walsh flatness across the
project-scale class.  They do not refute an implication whose hypothesis
explicitly says that the state is the output of the prescribed initialized
trajectory, and they do not calculate a basin for the original fixed bridge
`B_0`.  Nor is
`g_k` a generic random spectrum: it is a fixed three-level tensor family.

## 9. Consequence for the research target

The exact conclusion is narrower than “bad cores are generic,” but it is
enough to reject the proposed bent-only mechanism:

- self-dual bentness describes the original exceptional input `X=1`, but
  (65)--(73) show it is not necessary for a static simultaneous stable-core
  configuration;
- the first bad terminal is plateaued;
- leading defect is measured by a Walsh `L^1` deficit and forces enhanced
  fourth moment, rather than flatness; and
- simultaneous linear static cores persist when both input and terminal
  have full-support spectra uniformly separated from flatness; terminal
  spectra may also have arbitrarily many Walsh magnitudes.

A separate exact computation now shows that the prescribed path from the
three-level base `g_k` reaches multilevel, full- or almost-full-support cores
with `kappa/m=0.0930862` at `k=9`; it shadows a different explicit tensor
terminal whose static core density tends `1/8`.  See
[`regular_hadamard_tensor_trajectory.md`](regular_hadamard_tensor_trajectory.md).
Thus the remaining issue is not whether the dynamics can reach a non-bent
core, but whether a nonnegligible set of uniform inputs does so, or whether
all such basins decay with `k`.  A useful external theorem would therefore
count basins or regularizing/anti-regularizing trajectories under the
quadratic Walsh automorphism group.  A stability theorem for bent or
plateaued functions alone cannot control `kappa_*`.
