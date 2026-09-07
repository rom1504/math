# Wave 54: orientation pairs and the exact remaining K-profile obstruction

## Status

- **Verified:** flipping only the orientation preserves the cross vector,
  `K_H`, and `b_H`, and gives an exact affine relation between the two local
  margins. Any saved population whose opposite orientations remain in an
  `O(T_n)` local band already proves the restriction recurrence by inverse
  Hanson--Wright; no K-profile is needed on that population.
- **Verified conditional theorem:** this gives a strictly weaker sufficient
  abundance target than (10.1285). The K-profile only has to be proved on
  pairs for which both orientations lie below the local band.
- **Verified:** the two orientations admit an exact interval-counting formula
  for far K-success. Failure has two explicit depth modes once the
  coefficient capacity reaches the far scale.
- **Verified, but only a repackaging of (10.1280)--(10.1282):** a
  one-dimensional Huber formula makes the active clipped head canonical. A
  high-energy K-failure must store its excess energy in that head. This is a
  sharper failure description, not a stronger deterministic lower bound.
- **Numerical finite diagnostic:** the complementary two-orientation-negative
  branch is large in the stored exact examples and can fail the full exact
  K-profile even with `b_H=0`. This is not asymptotic evidence.
- **Open:** prove saved abundance for the relaxed target below, or construct an
  exact-minimizer family that makes it superexponentially small on the
  `H=n^{3/4-c}` scale.

## 1. Exact orientation-pair algebra

Retain the notation of (10.1243)--(10.1244) and (10.1280)--(10.1285):

```math
h(S)=Q(A[S])-p^{3/2}q-t,\qquad
e=\sigma y^{\mathsf T}A[S]y,\qquad
g=(1-p_2)e-h,
```

and, when `g<0`,

```math
r=-g/p_2,\qquad u=A[T,S]y,\qquad
K_H(u)=\tfrac12\mathcal K_H(2\sigma u).
```

Let `iota` flip `sigma` and leave `(S,y)` fixed. It is a
measure-preserving involution of the uniform oriented local-state law.
The quantities `h`, `u`, `G`, `J`, `K_H(u)`, and `b_H` are invariant under
`iota`; the sign changes in `beta=2 sigma u` and `B=sigma A[T]` do not
change their relevant norms. Directly,

```math
\boxed{
g^\iota=-g-2h.
}
\tag{R54.1}
```

Thus a negative state with `g=-p_2r` has

```math
\boxed{
g^\iota=p_2r-2h.
}
\tag{R54.2}
```

There is an orientation-free form. Put

```math
a=\left|y^{\mathsf T}A[S]y\right|,
\qquad d=(1-p_2)a.
```

The two members of the pair have margins

```math
\boxed{
g_+=d-h,\qquad g_-=-d-h,\qquad
\max\{g_+,g_-\}=d-h.
}
\tag{R54.3}
```

In particular, for every `B>=0`,

```math
\boxed{
g_+<-B\ \hbox{and}\ g_-<-B
\quad\Longleftrightarrow\quad h-d>B.
}
\tag{R54.4}
```

The role of `h` is important: it depends on the selector but is independent
of both the local state `y` and the orientation.

## 2. Saved local-band mass already proves the recurrence

Fix a compact density window in `(1/2,1)`, `0<c<1/4`, and set

```math
H=\left\lceil n^{3/4-c}\right\rceil,\qquad
T_n=n^{3/2-c}.
```

The following implication is the inverse-Hanson--Wright step used in Wave
51, stated here in the form needed for orientation pairs.

**Verified local-band extraction.** Let `B_n=O(T_n)` and
`t=O(T_n)`. If a set `E_n` of oriented local states satisfies

```math
\nu_m(E_n)\ge e^{-C_0H},\qquad g\ge-B_n\quad\hbox{on }E_n,
```

then

```math
\boxed{
q_m\le p^{3/2}q_n+O(T_n).
}
\tag{R54.5}
```

Indeed `g>=-B_n` gives

```math
h\le(1-p_2)e+B_n.
```

Uniformly in `S`, Hanson--Wright and the exact-minimizer estimates

```math
\lVert A[S]\rVert_F\le n,\qquad
\lVert A[S]\rVert_{\rm op}\le\sqrt{2q_n}=O(n^{3/4})
```

show that, for a sufficiently large constant depending on `C_0`,

```math
\Pr_{\nu_m}\left\{
e>C\left(n\sqrt H+\sqrt{q_n}\,H\right)\right\}
<\tfrac12e^{-C_0H}.
```

The displayed threshold is `O(T_n)` exactly when `c<=1/4`. Hence some
state of `E_n` has `e=O(T_n)`, and for its selector

```math
Q(A[S])
=p^{3/2}q_n+t+h
\le p^{3/2}q_n+t+B_n+(1-p_2)e
=p^{3/2}q_n+O(T_n).
```

Since `q_m<=Q(A[S])`, (R54.5) follows. No completion tail or independence
claim is used.

## 3. A strictly weaker abundance target than (10.1285)

Choose any fixed sufficiently large constant `C_B`, put
`B_n=C_B T_n`, and let `omega_n` tend to infinity slowly. Define

```math
\mathcal R_n=
\left\{g<0,\ r\ge\omega_nT_n:
\quad g^\iota\ge-B_n
\quad\hbox{or}\quad
K_H(u)\ge r+b_H\right\}.
\tag{R54.6}
```

Then the following is a **Verified conditional convergence theorem**:

```math
\boxed{
\nu_m(\mathcal R_n)\ge e^{-C_0H}
\quad\Longrightarrow\quad
q_m\le p^{3/2}q_n+O(T_n).
}
\tag{R54.7}
```

To prove it, split `\mathcal R_n` disjointly according as
`g^\iota>=-B_n` or `g^\iota<-B_n`. One part has at least half the
displayed mass. On the first part, apply `iota`; its image has the same
mass and lies in `g>=-B_n`, so (R54.5) applies. On the second part,
membership in (R54.6) forces the K inequality, and (10.1276) gives
conditional completion
probability at least `(1/2)12^{-H}`, whence
`Z_t>=e^{-O(H)}` and the bare-tail extraction applies.

This is strictly weaker than (10.1285): states in the first branch need not
satisfy any scalar, truncated-head, or exact K-profile condition. For all
large `n`, a state in (R54.6) already has `g<-B_n` because
`p_2 omega_n T_n>B_n`. Hence only the complementary branch

```math
g^\iota<-B_n
```

has to pay the K-profile. By (R54.4), this is exactly the
**selector-driven pair region**

```math
\boxed{
h-(1-p_2)\left|y^{\mathsf T}A[S]y\right|>B_n,
}
\tag{R54.8}
```

where both orientations lie below the `O(T_n)` local band.

The new open abundance lemma is therefore the left side of (R54.7), rather
than (10.1285). It is still not proved.

## 4. Exact pair-window count and the two surviving failure modes

For one unoriented pair put

```math
r_{\min}=\frac{h-d}{p_2},\qquad
r_{\max}=\frac{h+d}{p_2},\qquad
\kappa=K_H(u)-b_H.
```

Here a negative value of `r_min` simply means that orientation has
nonnegative `g`. Let

```math
A=p_2\omega_nT_n,\qquad C=p_2\kappa.
```

The exact number of orientations in this pair that are both far and
K-qualified is

```math
\boxed{
N_{\rm pair}
=\mathbf1_{\{A\le h-d\le C\}}
+ \mathbf1_{\{A\le h+d\le C\}}.
}
\tag{R54.9}
```

If `bar nu_m` is the unoriented pair law, uniform orientation gives exactly

```math
\boxed{
\nu_m\{r\ge\omega_nT_n,\ K_H(u)\ge r+b_H\}
=\frac12\mathbb E_{\bar\nu_m}N_{\rm pair}.
}
\tag{R54.10}
```

Thus the orientation bookkeeping costs exactly a factor `1/2`, not an
exponential loss.

Suppose a pair has a far orientation (`h+d>=A`), lies outside the direct
branch (`h-d>B_n`), and has coefficient capacity reaching the far scale
(`C>=A`). If it nevertheless contributes zero to (R54.9), exactly one of
the following holds:

1. **selector-depth obstruction:** `h-d>C`, so even the shallower
   orientation is deeper than the K capacity;
2. **orientation-gap obstruction:** `h-d<A` and `h+d>C`, so the two
   discrete depths jump across the entire payable interval `[A,C]`.

If `C<A`, the separate obstruction is simply that the exact coefficient
profile does not reach the far scale. These alternatives are exhaustive;
they isolate depth correlation from coefficient concentration without an
independence assumption.

## 5. Canonical clipped-head certificate for a profile-bad state

For `b_i=|u_i|` and `theta>0`, define

```math
\phi_\theta(b)=
\begin{cases}
b^2/(2\theta),&b\le\theta,\\
b-\theta/2,&b\ge\theta.
\end{cases}
```

Lagrange duality applied to (10.1280) gives the exact
one-dimensional formula

```math
\boxed{
K_H(u)=
\min_{\theta\ge0}
\left\{\frac{H\theta}{2}+\sum_i\phi_\theta(b_i)\right\},
}
\tag{R54.11}
```

where the value at `theta=0` is the limiting `l_1` value when the support has
at most `H` coordinates. For support larger than `H`, the unique positive
minimizer satisfies

```math
\sum_i\min\{1,b_i^2/\theta^2\}=H.
```

Put `I={i:b_i>=theta}`, `s=|I|`, and
`E_tail=sum_{i notin I}b_i^2`. Then `s<H`,

```math
\theta^2=\frac{E_{\rm tail}}{H-s},
```

and (R54.11) becomes the exact active-head identity

```math
\boxed{
K_H(u)=
\lVert u_I\rVert_1+\sqrt{(H-s)E_{\rm tail}}.
}
\tag{R54.12}
```

Consequently, if `R=r+b_H`, `G>=R^2/H`, but the exact profile fails
`K_H(u)<R`, its canonical head obeys

```math
\boxed{
\begin{aligned}
\lVert u_I\rVert_1&<R,\\
E_{\rm tail}
&<\frac{(R-\lVert u_I\rVert_1)^2}{H-s},\\
\sum_{i\in I}u_i^2
&>G-\frac{(R-\lVert u_I\rVert_1)^2}{H-s}.
\end{aligned}
}
\tag{R54.13}
```

Thus every energy-only survivor that is K-bad has a canonical set of fewer
than `H` outside vertices carrying the required energy surplus. This avoids
an arbitrary union over exceptional coordinate sets. Equations
(R54.11)--(R54.13) are equivalent in strength to the clipped formula
(10.1281); they do not prove that such canonical heads have small saved
mass.

## 6. Stored finite audit

The checker `tmp/joint_far_profile_r54_check.py` exhausts all target-density
`m>n/2` cases for the stored exact `A6`, `A8`, `A9`, and one exact `A10`
sample. It verifies (R54.1)--(R54.3) on 53,878 negative oriented states and
checks (R54.9) in 133,600 pair/cutoff tests and (R54.11) against the exact
K-functional on 3,400 random vectors.

For the deliberately finite diagnostic normalization

```math
H=\min\{3,n-m\},\qquad t=0,\qquad b_H=0,
```

the 53,878 negative states split as follows:

```text
opposite orientation nonnegative: 12,238
both orientations negative:       41,640
K_H(u) >= r in the latter branch: 11,266
K_H(u) <  r in the latter branch: 30,374
```

For example, at `A10,m=6`, 8,060 negative states have both orientations
negative; 6,330 of those fail `K_H(u)>=r` at `H=3`. Some failures have
`K_H(u)=0`. This proves only that the selector-driven and clipped-profile
obstructions are genuine finite phenomena. The choices `H<=3` and `b_H=0`
do not model the asymptotic far scale and provide no population evidence for
or against (R54.7).

## 7. Exact remaining step

The sharp next target is to prove, for some slowly diverging `omega_n` and
fixed constants `C_B,C_0`,

```math
\nu_m\left\{
g<0,\ r\ge\omega_nT_n:
g^\iota\ge-C_BT_n
\ \hbox{or}\
K_H(u)\ge r+b_H
\right\}
\ge e^{-C_0H}.
\tag{R54.14}
```

By the theorem above, (R54.14) proves convergence. On its only genuinely
new branch, (R54.8) holds. A proof can now target either the exact pair
window (R54.9), or show that canonical heads satisfying (R54.13) cannot
carry all saved far mass. A family of exact minimizers for which the
left-hand mass in (R54.14) is `e^{-\omega(H)}` for every usable choice would
falsify this relaxed orientation/K implementation, but no such family is
known.
