# The prescribed descent from the non-bent tensor regularizer

This note isolates a concrete dynamical target suggested by the prescribed
least-index, best-improvement descent for the regular-Hadamard obstruction.
All transforms below are unnormalized.  The exact statements concern the
symplectic Walsh matrix; the trajectory table is certified finite
computation.

## 1. An exact tensor terminal which the trajectory shadows

Let

```math
V_k=V_2\oplus V_{k-2},\qquad m=4^k=16N,\qquad
s=2^k=4t,
```

where `N=4^{k-2}` and `t=2^{k-2}`.  Write `H_k=H_2\otimes H_{k-2}` and
let `f=f_{k-2}` be the self-dual quadratic bent vector, so
`H_{k-2}f=tf`.  On `V_2`, in the coordinate order
`j=a_0+2a_1+4b_0+8b_1`, put

```math
w=(1,-1,-1,-1,1^{12}),\qquad v=(-1,1^{15}).       \tag{1}
```

The structured regularizing input and comparison terminal are

```math
g_k=w\otimes f,\qquad v_k=v\otimes f.               \tag{2}
```

Direct Walsh summation gives

```math
w\,(H_2w)=(10,6,6,6,2^{12}),
\qquad
v\,(H_2v)=(-14,-2^{15}).                            \tag{3}
```

Consequently `g_k` is strictly regularizing and `v_k` is strictly
anti-regularizing.  Their Walsh energies are

```math
E(g_k)={13\over16}ms,
\qquad E(v_k)=-{11\over16}ms.                       \tag{4}
```

Relative to `g_k`, the vector `v_k` flips precisely base fibres
`0,1,2,3`.  Thus the selected set has size `m/4`; its complement consists
of the twelve fibres of weight `s/2-1`.  For the anchored `-C` branch the
exact unmatched quantities are therefore

```math
\boxed{
\begin{aligned}
\delta(v_k)&=m\left[{s\over8}-2\right]_+,\\
u(v_k)&=m\left[{s\over16}-1\right]_+,\\
\kappa(v_k)&=
 \left\lceil {m(s-16)\over8(s-2)}\right\rceil\quad(s>16),
& {\kappa(v_k)\over m}&\longrightarrow {1\over8}.
\end{aligned}}                                      \tag{5}
```

For `s<=16`, `kappa(v_k)=0`.

This is a second exact non-flat tensor terminal.  It has the same normalized
energy as the previously used `z_*\otimes f`, and it is the *actual*
prescribed endpoint at `k=2,3`.

## 2. Certified trajectory pattern

The implementation in
`computations/audit_regular_hadamard_walsh_basins.py` maintains the Walsh
field exactly after every flip.  Its tie rule is maximum gain, then least
physical coordinate.  Let `z_k` denote the resulting endpoint and let
`h_k=d_H(z_k,v_k)`.

| `k` | `m` | gross flips | net selected | `h_k/m` | `E(z_k)/(ms)` | `delta/m^(3/2)` | `kappa/m` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 16 | 4 | 4 | 0 | -0.687500 | 0 | 0 |
| 3 | 64 | 16 | 16 | 0 | -0.687500 | 0 | 0 |
| 4 | 256 | 62 | 58 | 0.0859375 | -0.734375 | 0 | 0 |
| 5 | 1,024 | 244 | 244 | 0.0429688 | -0.710938 | 0.0390625 | 0.0136719 |
| 6 | 4,096 | 880 | 880 | 0.0351563 | -0.669922 | 0.111328 | 0.0424805 |
| 7 | 16,384 | 4,030 | 3,790 | 0.0186768 | -0.686539 | 0.110336 | 0.0741577 |
| 8 | 65,536 | 16,178 | 15,404 | 0.0149536 | -0.685043 | 0.119644 | 0.0904541 |
| 9 | 262,144 | 65,830 | 61,530 | 0.0152817 | -0.685310 | 0.123284 | 0.0930862 |

For `k=6,7,8`, a fibre-by-fibre audit gives the following *net* flip
counts.  The order of the entries is the base coordinate `0,...,15`.

```text
k=6: (256,208,208,208,0,0,0,0,0,0,0,0,0,0,0,0)
k=7: (1024,922,922,922,0,0,0,0,0,0,0,0,0,0,0,0)
k=8: (4096,3770,3769,3769,0,0,0,0,0,0,0,0,0,0,0,0)
k=9: (16384,15049,15049,15048,0,0,0,0,0,0,0,0,0,0,0,0).
```

The central empirical pattern is thus much more specific than spectral
near-plateauedness: from `k=6` through `k=9`, fibre `0` is completely
flipped, fibres `4,...,15` are untouched, and all discrepancies from the
exact tensor terminal lie in the three medium-weight fibres.  The three
discrepancy masks agree exactly for `k=6,7`; at `k=8`, two agree and the
third differs at one coordinate.  At `k=9`, the three unflipped counts are
`1,335,1,335,1,336`, again synchronized to within one.  The terminal
spectra are nevertheless many-level and essentially full-support.

## 3. A rigorous near-`v` transfer lemma

Here is a sufficient dynamical statement considerably weaker than proving
that the trajectory reaches `v_k`.

**Lemma.**  Suppose a Boolean endpoint `z` differs from `v_k` at `h`
coordinates, all contained in base fibres `1,2,3`.  Put `eta=h/m`.  Then

```math
\delta(z)\ge ms\left({1\over8}-\eta-{2\over s}\right).             \tag{6}
```

If the right side is positive, the unmatched core obeys

```math
{\kappa(z)\over m}
\ge {1\over 3-2/s}
 \left({1\over8}-\eta-{2\over s}\right).             \tag{7}
```

In particular, it is enough to prove the fibre confinement and

```math
\limsup_{k\to\infty}{h_k\over m_k}< {1\over8};        \tag{8}
```

neither exact tensor recursion nor `h_k=o(m_k)` is required.  Under (8),
the prescribed unmatched core is linear.

**Proof.**  Let `P=v_k 1_D`, where `D` is the difference set, so
`z=v_k-2P`.  Every coordinate of `D` has signed Walsh response
`v_k(H_kv_k)=-s/2`.  Therefore

```math
E(z)-E(v_k)=2sh+4P^TH_kP.                            \tag{9}
```

Write the three fibre masks as `D_1,D_2,D_3` and put

```math
p(y)=f(y)\sum_{i=1}^3 1_{D_i}(y).
```

The `3` by `3` principal block of `H_2` on base coordinates `1,2,3`
is all ones.  Hence `P^TH_kP=p^TH_{k-2}p`.  Since
`H_{k-2}` has least eigenvalue `-t` and

```math
\lVert p\rVert_2^2
=\sum_y\left(\sum_i1_{D_i}(y)\right)^2\le3h,
```

we have `P^TH_kP>=-3th=-(3s/4)h`.  Equation (9) consequently gives
`E(z)>=E(v_k)-sh`.  Combining this with (4) and
`delta=[E(g_k)+E(z)-2m]_+` proves (6).

Fibre `0` remains selected, so every available outside weight is at most
`3s/2-1`.  Covering unpaid mass `delta/2` therefore requires at least
`\delta/(3s-2)` coordinates, which gives (7).

The current data have `eta` near `0.015`, almost an order of magnitude
below the permissive threshold `1/8` in (8).

## 4. Exact reduced fixed-point equation

The same fibre face has an exact lower-dimensional description.  For
`i=1,2,3`, define `h_i:V_{k-2}\to\{\pm1\}` by

```math
z_i=-h_i f,
```

so `h_i=+1` is a discrepancy from `v_k`.  Set

```math
y_i=h_if,
\qquad R_i=f\,H_{k-2}y_i,
\qquad R=R_1+R_2+R_3.                               \tag{10}
```

A direct block Walsh multiplication shows that the signed field in each
medium fibre is

```math
z_i(H_kz)_i=h_i(5t+R).                              \tag{11}
```

Thus terminal stability forces

```math
h_i(u)(5t+R(u))\le0\quad(i=1,2,3).                  \tag{12}
```

Away from the neutral level `R=-5t`, all three masks must therefore agree.
Their one-coordinate disagreement at `k=8` lies exactly on this neutral
level.  This explains the observed near-synchronization without invoking
bent or plateaued structure.

When the masks agree, write their common value as `h`, put
`D={h=+1}`, and set `p=f1_D` and `r=fH_{k-2}p`.  Since
`hf=-f+2p`, (12), together with stability of the twelve untouched fibres,
reduces to the scalar signed-Walsh threshold system

```math
\begin{aligned}
u\in D&\Longrightarrow r(u)\le -t/3,\\
u\notin D&\Longrightarrow -t/3\le r(u)\le t.
\end{aligned}                                       \tag{13}
```

The endpoint problem on this face is therefore a self-consistent threshold
set on the tail cube of order `m/16`, not a full `m`-coordinate
maximization.  This does not yet prove the trajectory enters or remains on
the face, but it gives a precise falsifiable recursion target.

## 5. Exact missing dynamical lemma

The strongest next target exposed by the computation is:

> For all sufficiently large `k`, the prescribed least-index trajectory
> from `g_k=w\otimes f_{k-2}` flips every site in base fibre `0`, flips no
> site in fibres `4,...,15` in its terminal parity, and leaves fewer than
> `(1/8-epsilon)m` unflipped sites across fibres `1,2,3`, for one fixed
> `epsilon>0`.

By (6)--(8), this statement alone proves a uniform linear diffuse core for
the *actual prescribed trajectory*.  The finite data through `k=9` are
consistent with it, and (10)--(13) supply the compressed tail-state in
which a tensor/lex recursion should be sought.  No exact recursion between
successive computed discrepancy sets was detected: their four tail slices
do not equal the preceding set or its complement.  Thus claiming literal
tensor self-similarity would currently overstate the evidence.

There is also a concrete warning against batching the flips.  If one first
flips all of fibre `0` and only then runs best improvement restricted to
fibres `1,2,3`, the computation for `k=4,...,8` reaches a different family:
it flips exactly `9N/4` medium coordinates and has
`E/(ms)=-0.6640625` at every tested order.  Those states are globally
one-flip stable, but they are not the prescribed endpoints (their Hamming
distances from the latter are `20,64,312,1020,3896`).  Any successful lex
recursion must therefore preserve the actual interleaving of fibre-0 and
medium-fibre flips.
