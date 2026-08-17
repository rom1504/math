# Collective cross-Gram packing and the hard-case response modulus

**Status.** Rigorous task-local lower theorem and continuity theorem.  This
note uses the collective Boolean quadratic metric suggested after CG.2.  It
shows simultaneously that the metric space has exponential Boolean-realized
packing entropy and that, at normalized total repeated-port budget, its
SA.3 spherical response is uniformly continuous.  The correct uniform
modulus is square-root, not linear, because of the trust-region hard case.

## 1. The collective metric

For two Gram--Rayleigh pairs on `p` labelled ports define

```math
d_q((G,R),(G',R'))
={1\over p^2}\max_{epsilon\in\{+-1\}^p,\ sigma\in\{+-1\}}
 \left|epsilon^T[(G-G')+sigma(R-R')]epsilon\right|.       \tag{CP.1}
```

Writing `K^+-=(G+-R)/2`, this is equivalently

```math
d_q={2\over p^2}\max_{epsilon}
 \max\left\{|epsilon^T(K^+-K'^+)epsilon|,
             |epsilon^T(K^--K'^-)epsilon|\right\}.       \tag{CP.2}
```

This is a pseudometric on raw pairs and a metric only after quotienting by
equality of all displayed quadratic responses (or on a suitable normalized
slice).  It retains collective quadratic queries, but forgets the particular
entrywise presentation of the two PSD sectors.

## 2. An exponential Boolean-realizable packing

Fix a regular symmetric Hadamard involution `J=mathcal H/sqrt(n)` and a
Boolean top eigenvector `w`, so

```math
Jw=w,
\qquad ||w||_2^2=n.
```

For `s in {+-1}^p`, take the labelled Boolean ports

```math
w_i^(s)=s_iw.                                            \tag{CP.3}
```

Their Gram--Rayleigh sectors are

```math
G_s=R_s=ss^T,
\qquad K_s^+=ss^T,
\qquad K_s^-=0.                                         \tag{CP.4}
```

Thus all ports have the same self data `G_ii=R_ii=1`.

### Theorem CP.1 (linear-rate collective packing)

If `s,t in {+-1}^p` have Hamming distance `h`, then

```math
d_q((G_s,R_s),(G_t,R_t))
={8h(p-h)\over p^2}
=2\left(1-\left({s^Tt\over p}\right)^2\right).          \tag{CP.5}
```

Consequently, for every fixed `0<eta<1/2`, there are

```math
|mathcal C_p|
\ge2^{(1-H_2(eta)-o(1))p}                              \tag{CP.6}
```

Boolean-realizable Gram--Rayleigh pairs with common self data and pairwise

```math
d_q\ge8eta(1-eta).                                      \tag{CP.7}
```

Here `H_2` is binary entropy.

#### Proof

The sign vectors `s` and `-s` define the same pair, so work in the projective
cube.  Gauge by `s` and split the coordinates into the agreement set `A`
and disagreement set `D`, of sizes `p-h` and `h`.  For
`u_i=epsilon_i s_i`,

```math
\begin{aligned}
(epsilon^Ts)^2-(epsilon^Tt)^2
&=\left(\sum_Au_i+\sum_Du_i\right)^2
 -\left(\sum_Au_i-\sum_Du_i\right)^2\\
&=4\left(\sum_Au_i\right)\left(\sum_Du_i\right).
\end{aligned}                                           \tag{CP.8}
```

Both factors can be aligned independently, so the largest absolute value is
`4h(p-h)`.  In CP.1 only the `sigma=+1` sector remains and contributes an
additional factor two, proving CP.5.

Let

```math
d_proj(s,t)=min\{d_H(s,t),p-d_H(s,t)\}.
```

A projective Hamming ball of radius `r-1` contains at most
`sum_(k<r) binom(p,k)` classes.  Greedy packing of the `2^(p-1)` projective
classes therefore gives a code of minimum projective distance
`r=ceil(eta p)` and size at least

```math
{2^{p-1}\over\sum_{k<r}\binom pk}
\ge2^{(1-H_2(eta)-o(1))p}.                         \tag{CP.9}
```

For every pair choose the representative distance `h<=p/2`.  Then
`h>=eta p`, and CP.5 gives CP.7. `square`

This is a packing of labelled contextual states.  If the endpoint signs are
subsequently maximized and forgotten, every rank-one state has the same
scalar maximum.  The lower bound concerns the response table available to
a future labelled query, not an unlabelled terminal scalar.

## 3. The normalized trust-region response

For one endpoint word `epsilon` put

```math
k^+-=epsilon^TK^+-epsilon.
```

In outer channel `sigma`, call `k^sigma` the dangerous sector and
`k^(-sigma)` the safe sector.  Set

```math
a={k^sigma\over p^2},
\qquad b={k^{-sigma}\over p^2},
\qquad kappa={pm\over r}.                              \tag{CP.10}
```

Both `a,b` lie in `[0,1]`.  The trust-region dual identity used in SA.19
remains valid for every `m >= 0`; the condition `2m>r` in that theorem was
needed only for the one-port anti-pin certificate.  Substituting
`t=2alpha-1` in SA.19 shows that the
normalized spherical response is exactly

```math
Psi_kappa(a,b)
=\inf_{t>0}\left\{
 {1+t\over2}+{kappa^2\over2}
 \left({a\over t}+{b\over t+2}\right)
 \right\}.                                             \tag{CP.11}
```

The budget `m asymp r/p` is precisely the regime `kappa=Theta(1)`.

### Theorem CP.2 (uniform square-root response continuity)

Let two PSD Gram--Rayleigh pairs have

```math
d_q((G,R),(G',R'))<=delta.                            \tag{CP.12}
```

At common `kappa`, their complete SA.3 spherical response tables obey

```math
\max_{epsilon,sigma}
|Psi_kappa(a_(epsilon,sigma),b_(epsilon,sigma))
 -Psi_kappa(a'_(epsilon,sigma),b'_(epsilon,sigma))|
\le kappa\sqrt{delta/2}+{kappa^2delta\over8}.          \tag{CP.13}
```

The same bound holds after maximizing over `(epsilon,sigma)`.  Hence for
bounded `kappa`, `d_q=o(1)` gives uniform `o(1)` error in the response
normalized by `rn`.

#### Proof

Equations CP.2 and CP.12 give, for every `epsilon`,

```math
|a-a'|<=delta/2,
\qquad |b-b'|<=delta/2.                              \tag{CP.14}
```

The function in CP.11 is monotone in each coordinate.  Suppose first that
only `a` is increased by at most `e`.  Take an arbitrarily accurate
minimizer `t` for the old pair.  If `t>=kappa sqrt(e)`, use the same `t`; the
new term is at most `kappa sqrt(e)/2`.  If
`t<kappa sqrt(e)`, evaluate the new pair at
`t'=kappa sqrt(e)`.  The linear cost increases by at most
`kappa sqrt(e)/2`, the old `a/t` and safe terms decrease, and the new
`e/t'` term costs at most `kappa sqrt(e)/2`.  Therefore

```math
0<=Psi_kappa(a+e,b)-Psi_kappa(a,b)<=kappa sqrt(e). \tag{CP.15}
```

Increasing only `b` by `e` costs at most

```math
{kappa^2e\over2(t+2)}<={kappa^2e\over4}.           \tag{CP.16}
```

For two unordered pairs, compare both to their coordinatewise maximum and
use monotonicity.  Equations CP.15--CP.16 give

```math
|Psi_kappa(a,b)-Psi_kappa(a',b')|
<=kappa sqrt(e)+kappa^2e/4.
```

Set `e=delta/2`.  Taking a maximum over a common finite channel set is
nonexpansive in sup norm, proving the last assertion. `square`

## 4. The hard case makes the exponent sharp

The square root in CP.13 cannot be replaced uniformly by a linear modulus.
For `kappa^2b<4`, the boundary value is

```math
Psi_kappa(0,b)={1\over2}+{kappa^2b\over4}.          \tag{CP.17}
```

Using `1/(t+2)>=1/2-t/4` in CP.11 and, for the reverse estimate, testing
`t=kappa sqrt(a)`, gives

```math
kappa\sqrt{a\left(1-{kappa^2b\over4}\right)}
\le Psi_kappa(a,b)-Psi_kappa(0,b)
\le kappa\sqrt a.                                  \tag{CP.18}
```

This hard case occurs inside the admissible Gram--Rayleigh cone even with
fixed self data.  Take `p=2`, `s=(1,1)`, `t=(1,-1)`, and `0<c<1`.  Compare

```math
\begin{array}{c|cc}
 &K^+&K^-\\ \hline
\mathcal A&css^T&(1-c)I\\
\mathcal B&ctt^T&(1-c)I.
\end{array}                                           \tag{CP.19}
```

Both pairs are PSD and have identical diagonal data
`G_ii=1`, `R_ii=2c-1`.  They are Euclidean-realizable using one common
involution: realize the rank-one factors in its positive eigenspace and the
orthogonal diagonal factors in its negative eigenspace.  Directly,

```math
d_q(\mathcal A,\mathcal B)=2c.                     \tag{CP.20}
```

At query `epsilon=t` in the positive outer channel, the two normalized
sector pairs are

```math
(a,b)=\left(0,{1-c\over2}\right),
\qquad
(a',b')=\left(c,{1-c\over2}\right).               \tag{CP.21}
```

For fixed `kappa<sqrt(8)` (more generally, with a fixed margin below
`kappa^2(1-c)=8`), CP.18 makes their response gap
`Theta(kappa sqrt(c))=Theta(kappa sqrt(d_q))`.  Thus the loss in CP.13 is a
real boundary phenomenon, not a defect of the proof.

## 5. Packing states are also contextually visible

For the Boolean rank-one states of CP.4, the positive-channel response has
`b=0` and hence the exact form

```math
Psi_kappa\left({(epsilon^Ts)^2\over p^2},0\right)
={1\over2}+{kappa|epsilon^Ts|\over p}.             \tag{CP.22}
```

If the projective Hamming distance between `s,t` is `h<=p/2`, choose
`epsilon=s`.  The two table entries differ by

```math
{kappa\over p}\left(p-|p-2h|\right)
={2kappa h\over p}.                                \tag{CP.23}
```

The code in CP.1 therefore gives an
`exp(Omega_eta(p))` packing of normalized SA.3 **response tables** at
separation at least `2kappa eta` whenever `kappa` is bounded below.  This is
a genuine collective-response lower bound in a Boolean-realizable class,
although, as noted above, optimizing away the query label collapses these
gauge-related states to one scalar value.

## 6. Consequence and scope

At total repeated-port budget `pm/r=Theta(1)`, the collective metric has the
right two properties for an approximate carrier:

1. `d_q`-accuracy `delta` controls the complete normalized spherical table
   to error `O(sqrt(delta))`;
2. fixed response accuracy still requires `exp(Omega(p))` states, already on
   Boolean rank-one top-eigenvector ports.

This does not yet give an `exp(O(p))` upper cover of the full PSD pair space,
nor does it transfer through the Boolean integrality gap of SA.3.  It does
show that the hard case is quantitatively manageable and identifies the
natural target: cover `(K^+,K^-)` in the collective quadratic metric at
constant accuracy, with exponential rather than quadratic-exponential
state count.
