# Exact zero-Gaussian-variance hard-threshold counterexample

Date: 2026-09-06. Complete construction and proof, independently
reconstructed with the director. This is a sharp scope/limit-order
counterexample, not a new universal lower bound or a convergence result.

## 1. Statement

There is a sequence of ACTUAL symmetric hollow signings with uniformly
bounded normalized operator norm, and FIXED bounded functions f(g,y),
H(g,y), such that:

1. f is odd, H is even and nonnegative, and |f|+H=1.
2. Both functions are Gaussian-a.e.-continuous, and f is Gaussian-a.e.
   zero. Every Gaussian Hermite coefficient of f, including its residual
   variance tau^2, is therefore exactly zero.
3. On the actual old fields G=BS, D=S h_2(G), Y=BD, the feasible feedback
   C=H(G,Y) sign(Bf(G,Y)) satisfies

```math
\lim_n\frac{\mathbb E C^TBC}{2n}=\frac1{\pi\sqrt2}>0.       (1)
```

4. The normalized second moment of Bf tends to zero, and every fixed
   Lipschitz odd bounded feedback response in place of sign has energy
   tending to zero.

Thus the positive-residual-variance hypothesis in the hard-threshold
corollary cannot simply be omitted, even with Gaussian-a.e.-continuous
bounded input responses and actual dense signings. The smooth marked
energy theorem is not contradicted.

## 2. Squarefree principal compressions of actual Steiner signings

Use the exact binary-Steiner family already constructed in
`continued_director_steiner_feedback_test_2026_09_06.md`. Its orders N
grow geometrically, and its hollow sign matrix A_N satisfies

```math
A_N^2=c_N A_N+(N-1)I,\qquad
\gamma_N=c_N/\sqrt{N-1}\longrightarrow1/\sqrt2.
```

The normalized matrices B_N have uniformly bounded operator norm,
unit row squared norms, and Tr(B_N^3)/N=gamma_N.

For each sufficiently large N choose a squarefree integer m in
[N-h,N-1], where h=ceil(4 sqrt(N)), and put n=m+1. Such an m exists
elementarily. A nonsquarefree integer is divisible by d^2 for some
2<=d<=sqrt(N). The number of such integers in an interval of length h
is at most

```math
h\sum_{d=2}^{\infty}d^{-2}+\sqrt N
\le \frac{25}{36}h+\frac14h=\frac{17}{18}h<h.
```

Here sum_(d>=2) d^-2 <=1/4+1/9+integral_3^infinity x^-2 dx=25/36.
The intervals for successive sufficiently large Steiner orders are
disjoint, so the selected squarefree m's are distinct.

Take any n-vertex principal compression A of A_N, let J be the deleted
set of size ell=N-n=O(sqrt(N)), and set B=A/sqrt(m). This is still a
full hollow signing, and

```math
\|B\|_{op}\le\sqrt{\frac{N-1}{m}}\|B_N\|_{op}=O(1).
```

Writing Q=B^2, its diagonal entries equal one exactly. For i!=j in
the retained set I, the exact compression identity gives

```math
Q_{ij}=\frac{c_N A_{ij}-(A_{IJ}A_{JI})_{ij}}m.
```

Therefore max_(i!=j)|Q_ij|<= (c_N+ell)/m=O(n^-1/2).
The normalized cubic trace also survives:

```math
\frac1n\operatorname{Tr}(B^3)\longrightarrow1/\sqrt2.       (2)
```

For completeness, embed the compressed normalized matrix as P B_N P
in dimension N. Its difference from B_N has rank at most 2ell and
bounded operator norm, hence trace norm O(ell). Telescoping the cubes
shows their cubic traces differ by O(ell). The change in normalization
and the factor N/n both tend to one, proving (2).

## 3. A genuine fixed positive-variance zero-first test response

Let c=exp(3/2)/2 and choose

```math
u(g)=\frac{\sin g-c\sin(2g)}{1+c}.
```

This is a bounded continuous odd function, |u|<=1. For standard normal
N, E[N sin(tN)]=t exp(-t^2/2), so E[N u(N)]=0 exactly. Its variance
tau_u^2 is strictly positive and its Hermite expansion has only odd
degrees at least three.

Put

```math
R=\sum_{p\ge3\ {\rm odd}}u_p^2Q^{\circ p},\qquad T=BRB.
```

The off-diagonal bound from Section 2 and Tr(Q^2)=O(n) give

```math
\|R-\tau_u^2 I\|_F^2
\le \tau_u^4\big(\max_{i\ne j}|Q_{ij}|\big)^4
                    \sum_{i\ne j}Q_{ij}^2=O(n^{-1}).
```

Thus ||T-tau_u^2 Q||op=O(n^-1/2), and uniformly in i,
T_ii=tau_u^2+O(n^-1/2). Applying the independently audited zero-first
hard-feedback theorem to this GENUINE fixed positive-variance u gives

```math
\frac1{2n}\mathbb E[\operatorname{sign}(Bu(G))^T
                 B\operatorname{sign}(Bu(G))]
=\frac1{\pi n}\operatorname{Tr}
 (B D_{\sigma^{-1}}T D_{\sigma^{-1}})+o(1)
=\frac1{\pi n}\operatorname{Tr}(B^3)+o(1)
\longrightarrow\frac1{\pi\sqrt2}.                          (3)
```

There is no use of a zero-variance threshold theorem in this step.
The comparison variance for u is uniformly bounded below, and the
first-variance normalization in (3) is consequently harmless.

## 4. A fixed Thomae-type function encoding the actual discrete grids

For each chosen squarefree m let

```math
\mathcal D_m=\{k/\sqrt m:\ |k|\le m,\ k\equiv m\pmod2,
                                      \ k\ne0\}.
```

These finite nonzero grids are disjoint. Indeed an equality
k/sqrt(m)=k'/sqrt(m') with nonzero integers k,k' forces m=m', by the
parity of prime valuations in k^2 m'=k'^2 m for squarefree m,m'.

Define ONE function on the real line, once for the entire sequence:

```math
w(g)=\begin{cases}
 (m+1)^{-1},&g\in\mathcal D_m\text{ for a selected }m,\\
 0,&\text{otherwise}.
\end{cases}
```

It is even, bounded, and supported on a countable set. At every point
outside that set it is continuous with value zero: for each positive
accuracy only finitely many grids have weights exceeding that accuracy,
and their union is finite. Thus w is Gaussian-a.e.-continuous. The
point zero causes no exception and w(0)=0.

Set the FIXED old-frame functions

```math
f(g,y)=w(g)u(g),\qquad H(g,y)=1-|f(g,y)|.
```

They have all the regularity, parity and feasibility properties in
Section 1, and f is Gaussian-a.e. zero. They do not depend on n once
the selected infinite sequence and its grids have been defined.

For EVERY Boolean seed vector and every row at selected order n=m+1,
G_i is an integer k divided by sqrt(m), with the indicated parity and
range. If k!=0 it belongs to D_m; if k=0 then u(G_i)=0. Consequently
the following is an exact samplewise identity:

```math
f(G_i,Y_i)=\epsilon_m u(G_i),\qquad \epsilon_m=(m+1)^{-1}.
```

In particular Bf=epsilon_m Bu(G), and positive scaling disappears under
the hard threshold. The actual feedback differs uniformly by at most
epsilon_m from sign(Bu(G)) because H=1-epsilon_m|u(G)|. Bounded-op
energy continuity therefore transfers (3) to C, proving (1).

The endpoints f+C and -f+C lie in the cube exactly. Also
E||Bf||^2/n<=||B||op^2 epsilon_m^2. For a fixed odd Lipschitz psi,
|psi(t)|<=Lip(psi)|t|, so its feedback energy is O(epsilon_m^2).

## 5. Interpretation and limits of the counterexample

The discontinuity at zero amplifies a vanishing actual-law residual
that Gaussian L2 identifies with zero. This is an intentional example
of why taking n to infinity and then taking a hard-threshold limit
need not commute when tau=0. It is not an error in the smooth theorem,
which predicts vanishing energy here correctly.

The singular functions are fixed and regular in exactly the stated
Gaussian-a.e. sense; they would be excluded by a global continuity
assumption on f. In fact the bounded globally continuous tau=0 case is
trivial: f equals its Gaussian first linear projection almost everywhere,
hence everywhere by continuity and full Gaussian support; boundedness
then forces that linear projection and f itself to vanish. Thus no
unresolved globally continuous zero-variance case is asserted here.
The positive feedback value 1/(pi sqrt2) is below the preserved universal
bound and supplies no improvement to that bound.

## 6. Reproducible finite algebra checks

`computations/continued_audit_thomae_threshold_counterexample_2026_09_06.py`
checks squarefree grid disjointness using integer quadratic-irrational
representations, the exact integer compression identity, parity of every
actual field, and feasible endpoints. It also reports the positive-
variance u comparison separately from the identically zero Gaussian
variance of the fixed Thomae response. Its finite Monte Carlo checks
are corroboration only; the limit (1) is proved above.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_audit_thomae_threshold_counterexample_2026_09_06.py --samples 6000
```
