# Gaussian conference rounding and the sharp local-search trap

Status: exact `1/pi` Gaussian rounding theorem and a sharp universal barrier
for deterministic one-flip local improvement.  This architecture does not
prove a `(1/2-o(1))N^(3/2)` conference lower bound.  A triangle flip escapes
the barrier and leaves `Omega(N)` improving additions at each subsequent
fixed radius; the uniform theorem continues only to radius
`Omega(sqrt(N))`.

## 1. Exact top-eigenspace Gaussian identity

Let `S` be any symmetric conference signing of even order `N`, and put

```math
q=N-1,\qquad \lambda=\sqrt q,\qquad
P_+={1\over2}(I+S/\lambda).                          \tag{GR1}
```

Draw a centered Gaussian `g` with covariance `P_+` and set
`x_i=sign(g_i)`.  Since `(P_+)_{ii}=1/2` and

```math
\operatorname{Corr}(g_i,g_j)={S_{ij}\over\lambda}, \tag{GR2}
```

the Gaussian arcsine identity gives

```math
\mathbb E[x_ix_j]
={2\over\pi}\arcsin(S_{ij}/\lambda).                \tag{GR3}
```

Every off-diagonal entry is a sign, so multiplication by `S_ij` removes its
sign inside the odd function `arcsin`.  Therefore

```math
\boxed{
\mathbb E H_S(x)
={N(N-1)\over\pi}
 \arcsin{1\over\sqrt{N-1}}.}                        \tag{GR4}
```

In particular,

```math
\operatorname{cap}(S)
\ge {N(N-1)\over\pi}
 \arcsin{1\over\sqrt{N-1}}
=\left({1\over\pi}+o(1)\right)N^{3/2}.             \tag{GR5}
```

This expectation is identical for every symmetric conference signing.  It is
the exact limit of what follows by selecting a sample no worse than the mean;
its asymptotic constant is `1/pi`, not `1/2`.

Sequentially flipping any coordinate with negative signed local field cannot
decrease the energy, so Gaussian rounding followed by coordinate ascent still
has the expectation lower bound (GR4).  The next sections show that no
pointwise improvement beyond it can be forced from local optimality.

## 2. Sharp energy floor for one-flip local maxima

For a Boolean spin `x`, define its signed local fields

```math
t_i=x_i(Sx)_i.                                      \tag{GR6}
```

Flipping coordinate `i` changes the energy by `-2t_i`.  Thus a strict local
maximum for positive energy has `t_i>0` for every `i`.  Each `t_i` is odd,
and conference orthogonality gives

```math
\sum_i t_i^2=\|Sx\|_2^2=Nq.                        \tag{GR7}
```

Also `1<=t_i<=q`.  On this interval,

```math
t_i^2-1\le(q+1)(t_i-1).                             \tag{GR8}
```

Summing (GR8), using `N=q+1` and (GR7), gives

```math
q^2-1
=\sum_i(t_i^2-1)
\le(q+1)\left(\sum_i t_i-N\right).
```

Consequently every positive one-flip local maximum satisfies

```math
\boxed{H_S(x)={1\over2}\sum_i t_i\ge q=N-1.}       \tag{GR9}
```

This linear bound is sharp for **every** conference signing.

## 3. Universal top-space rounding trap

Fix any vertex `r` and define

```math
x=e_r+Se_r.                                         \tag{GR10}
```

Because the diagonal of `S` is zero, (GR10) is Boolean: `x_r=1` and
`x_i=S_{ir}` off the root.  Using `S^2=qI`,

```math
Sx=Se_r+qe_r.                                       \tag{GR11}
```

Hence the signed local-field multiset is exactly

```math
\{t_i}=\{q,1,1,\ldots,1\},                        \tag{GR12}
```

and

```math
H_S(x)=q=N-1.                                       \tag{GR13}
```

Thus (GR10) attains equality in the sharp local-maximum bound (GR9), has no
improving coordinate flip, and is fixed by the synchronous majority update:

```math
\operatorname{sign}(Sx)=x.                         \tag{GR14}
```

More importantly, it is itself a top-eigenspace sign rounding.  Indeed

```math
g=P_+x={1\over2}(x+Sx/\lambda)                     \tag{GR15}
```

lies in the positive eigenspace, and

```math
x_i g_i={1\over2}(1+t_i/\lambda)>0                 \tag{GR16}
```

for every coordinate.  Therefore `sign(g)=x`.  Since all inequalities in
(GR16) are strict, a relatively open neighborhood of `g` in the top
eigenspace has the same sign pattern.  The Gaussian in Section 1 assigns this
cone positive probability.

Equations (GR10)--(GR16) are the sharp barrier:

> A top-eigenspace Gaussian sign pattern can be a strict coordinate-wise
> local maximum of energy only `N-1=o(N^(3/2))`.  One synchronous update, or
> any number of energy-improving single-coordinate updates, leaves it fixed.

This does not lower the expectation (GR4), because the trap cone may have
small probability.  It does rule out any per-sample or worst-cone theorem
claiming that top-space sign rounding plus local improvement must approach
the spectral constant `1/2`.

## 4. Exact flip radius: two is trapped, three always escapes

The root trap is robust to two simultaneous flips, but not to three once
`N>=10`.  This follows from an exact normalized-core reduction.

Switch `S` by the trap spin so that the trap becomes the all-one vector, and
retain the notation `T` for the switched conference matrix.  Its root row is
all positive.  Delete that row and column, obtaining the order-`q` normalized
core `C`.  Equations (GR11)--(GR12) imply

```math
C\mathbf1=0,
\qquad C^2=qI_q-J_q.                                \tag{GR17}
```

Any flip set is energy-equivalent, after taking its complement if necessary,
to a set `F` not containing the root.  Let `m=|F|`, and let `p(F)` be the
number of positive edges induced by `F` in `C`.  Its signed internal edge sum
is

```math
e_C(F)=2p(F)-{m\choose2}.                            \tag{GR18}
```

The signed cut from `F` to its complement is

```math
\operatorname{cut}_T(F,F^c)
=m-2e_C(F)=m^2-4p(F).                               \tag{GR19}
```

Flipping `F` negates exactly this cut, so the exact energy change is

```math
\boxed{
H_T(\mathbf1-2\mathbf1_F)-H_T(\mathbf1)
=2(4p(F)-m^2).}                                     \tag{GR20}
```

Therefore

```math
\boxed{F\text{ improves the root trap}
\iff p(F)>m^2/4.}                                   \tag{GR21}
```

For `m=1` this is impossible.  For `m=2`, even a positive edge has
`p(F)=1=m^2/4`, so it is a neutral move rather than an improvement.  Hence the
root construction is always a radius-two local maximum against strict
improvements, although positive-edge pairs give neutral moves.

The graph of positive entries in `C` is the conference graph with parameters

```math
\left(q,{q-1\over2},{q-5\over4},{q-1\over4}\right). \tag{GR22}
```

For `q>=9`, every positive edge lies in `(q-5)/4>=1` positive triangles.
A positive triangle has `m=3,p=3>9/4`; by (GR20), flipping it increases the
energy by exactly six.  No smaller strict improvement exists.  Thus

```math
\boxed{
N\ge10\quad\Longrightarrow\quad
\text{the smallest improving flip set has size exactly }3.} \tag{GR23}
```

At order six, the positive core graph has parameters `(5,2,0,1)` and is a
five-cycle.  It is triangle-free, so Mantel's theorem gives
`p(F)<=m^2/4` for every subset.  There is no improving flip set of any size;
the root trap is a global maximizer of energy five.  The trivial order-two
case behaves similarly.

Consequently the universal root construction is a scalable radius-two trap,
but **not** a scalable bounded-radius trap for any radius at least three.
Three-spin local search always escapes it at every nontrivial growing order.
This escape gains only six energy units, so (GR23) does not imply flow toward
the project scale; it precisely separates the one/two-flip obstruction from
a robust landing obstruction.

### 4.1 Exact post-triangle local fields

The first escaped state nevertheless has many improving directions.  Let
`F` be a positive triangle and let `rho` be the number of vertices completing
it to a positive four-clique.  For a vertex outside `F`, write `a` for its
number of positive neighbors in the triangle.  Its signed local field after
flipping `F` is

```math
t=7-4a.                                             \tag{GR24}
```

The root field is `q-6`, and each of the three flipped triangle vertices has
field three.  The conference-graph intersection numbers determine the full
outside distribution up to `rho`:

```math
\begin{array}{c|cccc}
a&0&1&2&3\\ \hline
\#\text{ vertices}
&{q-9\over4}-\rho&6+3\rho
&{3(q-9)\over4}-3\rho&\rho.
\end{array}                                        \tag{GR25}
```

The last two classes have fields `-1` and `-5`.  Since the first count in
(GR25) is nonnegative, `rho<=(q-9)/4`, and hence

```math
\boxed{
\#\{i:t_i<0\}
={3(q-9)\over4}-2\rho
\ge {q-9\over4}.}                                  \tag{GR26}
```

Thus the triangle escape creates `Omega(N)` simultaneously available
one-coordinate improvements.  A field `-1` flip gains two and a field `-5`
flip gains ten.  These directions are not independent: performing one flip
changes every other local field by `+/-2`, so (GR26) alone does not prove
`Omega(N)` successive improvements.

### 4.2 Majority dynamics and an `Omega(sqrt(N))` basin-radius theorem

There is an exact description at every later state.  For a flipped core set
`F` of size `m`, let `d_F(v)` be the number of positive core edges from `v`
to `F` (excluding `v` itself when `v in F`).  Directly from (GR19),

```math
t_r=q-2m,
\qquad
t_v=
\begin{cases}
2m+1-4d_F(v),&v\notin F,\\
4d_F(v)-2m+1,&v\in F.
\end{cases}                                        \tag{GR27}
```

Here `r` is the normalized root.  Therefore adding an outside vertex improves
precisely when it has a strict positive-edge majority into `F`; removing an
inside vertex improves precisely when it has less than half.  A nonempty
coordinate-wise local maximum must satisfy

```math
d_F(v)\ge\lceil m/2\rceil\quad(v\in F),
\qquad
d_F(v)\le\lfloor m/2\rfloor\quad(v\notin F),
\qquad
m\le(q-1)/2.                                       \tag{GR28}
```

The last condition is the root-coordinate condition.  It also ensures that
`m` is the projective Hamming distance from the root trap, rather than the
distance to the farther of the two globally opposite representatives.

Conditions (GR28) already exclude every bounded-size terminal basin as the
conference order grows.  Put `h=floor(m/2)`, let `M=binom(m,2)`, and write

```math
E=\sum_{v\notin F}d_F(v),
\qquad B=\sum_{v\notin F}{d_F(v)\choose2}.          \tag{GR29}
```

If the induced positive-edge count is `p`, the conference strongly regular
parameters give

```math
E=mk-2p,
\qquad
B=M\mu-p-\sum_{v\in F}{d_F(v)\choose2},             \tag{GR30}
```

where `k=(q-1)/2` and `mu=(q-1)/4`.  The outside condition in (GR28) implies

```math
B\le {h-1\over2}E.                                  \tag{GR31}
```

The same comparison is quantitative before terminality is assumed.  Let
`R(F)` count the outside vertices with `d_F(v)>=h+1`, hence the immediately
improving additions.  Since

```math
B-{h-1\over2}E
=\sum_{v\notin F}{d_F(v)(d_F(v)-h)\over2},
```

the nonviolating terms are nonpositive, while every term is at most
`m(m-h)/2`.  Substitution from (GR30) gives the exact identity

```math
B-{h-1\over2}E
={q-1\over4}\bigl(M-m(h-1)\bigr)
 +(h-2)p-\sum_{v\in F}{d_F(v)\choose2}.
```

For `m>=4` the coefficient of `p` is nonnegative.  Using
`sum_{v in F} binom(d_F(v),2)<=m binom(m-1,2)`, gives, for `m>=4`,

```math
R(F)\ge
\begin{cases}
\displaystyle {q-1\over2m}
-{2(m-1)(m-2)\over m},&m\text{ even},\\[6pt]
\displaystyle {q-1-2(m-1)(m-2)\over m+1},&m\text{ odd}.
\end{cases}                                        \tag{GR31a}
```

For `m=3`, an improving set must be a positive triangle and (GR26) supplies
the stronger exact count.  Consequently every fixed-size state reached from
the triangle has `Omega(q)` immediately improving additions; more uniformly,
the count is `Omega(q/m)` throughout `m=o(sqrt(q))`.  Choosing only such
additions gives a strictly energy-increasing trajectory out to
`Omega(sqrt(q))` flips.  This is a supply of linearly many *simultaneous*
moves at bounded radius, not linearly many successive moves.

If `F` is terminal and `m>=4`, then `R(F)=0`.  The lower bounds above
therefore give the necessary conditions

```math
q-1\le
\begin{cases}
4(m-1)(m-2),&m\text{ even},\\
2(m-1)(m-2),&m\text{ odd}.
\end{cases}                                        \tag{GR32}
```

For `q>=9`, sizes one and two are excluded directly.  At size three, the
inside condition forces a positive triangle, and (GR26) shows that terminality
is possible only when `q<=9`.  Combining this observation with (GR32), every
nonempty coordinate-wise local maximum at `q>=9` obeys

```math
\boxed{m\ge{3+\sqrt q\over2}.}                     \tag{GR33}
```

After the triangle move raises the energy above the root trap, monotone
coordinate ascent cannot return to the empty flipped set.  It must therefore
make enough net flips to reach Hamming distance `Omega(sqrt(N))` before it can
terminate.  This is stronger than a fixed-radius escape but much weaker than
`Omega(N)`, and it gives no project-scale energy: (GR20) can still be only
`O(1)` above the linear root energy for a majority-dense set.

Past (GR33), continuation depends on triple and higher intersection data of
the particular conference graph, which are not fixed by the conference
identities.  Maximizing (GR20) or classifying all sets satisfying (GR28) is
the original Boolean cap problem in the normalized-core language.  Thus this
is the natural stopping point for the local-ascent derivation: it proves an
unbounded basin radius, but not a new normalized cap constant or a monotone
bound beyond `1/pi`.

## 5. Why the arithmetic local-field deficit does not amplify rounding

On orders `N=4k^2+2`, the arithmetic identity from the two-fiber audit is

```math
2ks-H_S(x)
={1\over8k}\sum_i
(t_i-(2k-1))(t_i-(2k+1))                            \tag{GR34}
```

for positive energy.  It is an exact upper-deficit identity.  At the root
trap, all local fields are positive, so local ascent has no move, while the
multiset `{q,1,...,1}` gives a project-scale right side in (GR34).  Thus a
large arithmetic penalty does not imply the existence of a negative local
field or any coordinate improvement.

Taking expectations in (GR34) simply rewrites the gap between (GR4) and the
arithmetic upper bound; it supplies no amplification.  Any use of (GR34) must
control the **distribution** of the penalty across a randomized trajectory,
not merely its sum or the sign of the local fields.

## 6. Exact finite verification and research judgment

Exhaustive enumeration gives the following positive local maxima for the
cyclic conference certificates:

```text
order 6:  6 projective strict local maxima, all energy 5;
order 18: 120 projective strict local maxima,
          18 at energy 17 and 102 at energy 33.
```

The 18 energy-17 states are precisely the sharp scale predicted by the root
traps.  The reproducer also verifies the root construction at order 38.

The proved universal conference lower bound from this architecture remains
`(1/pi+o(1))N^(3/2)`.  A `(1/2-o(1))` theorem is not obtained.  The exact trap
shows that deterministic local improvement cannot bridge the gap uniformly.
A surviving probabilistic route would need a theorem that the low-energy
trap cones have sufficiently tiny Gaussian measure and that almost all of the
remaining mass flows to near-spectral basins.  Neither conference
orthogonality nor the local-field moment/deficit identities alone provides
such a basin-volume estimate.

## Reproduction

```bash
.venv/bin/python computations/audit_conference_gaussian_local_rounding.py \
  --output computations/results/conference_gaussian_local_rounding.json
```

The program verifies the root traps by exact integer arithmetic, exhausts all
projective spins at orders 6 and 18, and records the exact Gaussian expectation
formula and normalized constants.
