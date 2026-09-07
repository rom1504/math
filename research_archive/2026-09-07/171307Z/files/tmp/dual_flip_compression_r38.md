# Wave 38 S1: uniform-noise compression lies beyond the complement cap ceiling

## Outcome

The attempted high-cap compression theorem is rigorous as a conditional
statement, but its selector class is **empty under the theorem's own row
hypotheses** and it therefore gives no partial form of (10.1047).  Put

```math
B_S=A^{F_S}=-A+2P_SAP_S,\qquad M_S=Q(B_S),\qquad q=q_n.
```

For a row cutoff `R`, every feasible dual law in (10.1046) has only

```math
\exp\!\left\{O\!\left(\frac{nq}{t}+\log n\right)\right\}
```

total mass on selectors for which `M_S>=t`, provided

```math
t\ge K\frac{nq^2}{R}
```

for a universal constant `K`, in the project regime specified below.  At

```math
R\asymp n^{9/4-c},qquad q\asymp n^{3/2},qquad
t=K\frac{nq^2}{R}\asymp n^{7/4+c},
```

this is exactly

```math
\exp\{O(R/q)\}=\exp\{O(n^{3/4-c})\}.
```

The obstruction is even simpler than the spectral estimate: the structured
complement flip has the exact norm ceiling

```math
M_S=Q(-A+2P_SAP_S)
\le Q(A)+2Q(A[S])\le3q.
```

Indeed, the theorem's own hypothesis `R<=8nq` makes its required threshold
at least `16q`.  The high-cap selector class is therefore literally empty.
In particular, the unresolved part of (10.1047) is still the entire selector
family, even though it lies in the nominal range

```math
Q(A^{F_S})=O(nq^2/R)=O(n^{7/4+c}).
```

The conditional theorem below explains the failure sharply.  It noises a ground
of `B_S` almost halfway around the spin cube and would make high complement
cap pay simultaneously for row regularization and Hamming-sphere entropy.
But the triangle cap bound prevents a complement signing from having the
amount of cap needed to pay even the theorem's weakest admitted row budget.

The same audit also gives a sharp obstruction to the most direct weighted
first-moment argument.  A normalized dual can be exactly uniform on the
selector slice, in which case its mean complement matrix is merely
`(2p_2-1)A`, strictly below the minimality threshold.  Thus edgewise averaging
of the exact complement inequalities cannot control the remaining
moderate-cap selectors.

No convergence proof, nonvacuous asymptotic compression, or actual-minimizer
counterexample is obtained.

## 1. Exact dual consequences and the first-moment wall

Let

```math
L_S(d)=\langle B_S,d\rangle=2c_S(d)-E_d,
\qquad
I_d=\{S:L_S(d)\ge q\}.
```

Let `y_S>=0` be any feasible dual law,

```math
\sum_{S\in I_d}y_S\le1
\quad\text{for every }d\text{ with }R_2(d)\le R,
\qquad Y=\sum_Sy_S,
```

and, when `Y>0`, put `w_S=y_S/Y`.  Reversing the orientation sends
`L_S(d)` to `-L_S(d)` and leaves `R_2(d)` unchanged.  Consequently every
row-good physical cut satisfies the exact two-sided tail bound

```math
\boxed{
\Pr_{S\sim w}\{|L_S(d)|\ge q\}\le\frac2Y.
}
\tag{R38.1}
```

Because `|E_d|<=q`, `|c_S(d)|<=Q(A[S])<=q`, and hence `|L_S(d)|<=3q`,
(R38.1) also gives

```math
\boxed{
\mathbb E_w(|L_S(d)|-q)_+\le\frac{4q}{Y},
\qquad
|\mathbb E_wL_S(d)|\le q+\frac{2q}{Y}.
}
\tag{R38.2}
```

If `r_{ij}(w)=Pr_w(i,j in S)`, the averaged complement matrix is exactly

```math
\overline B_w
=\mathbb E_wB_S
=A\circ(2R_w-\mathbf1_{\ne}).
\tag{R38.3}
```

Thus (R38.2) says that its row-restricted quadratic cap is at most
`q+2q/Y`.  This is a genuine consequence of a large dual, but it cannot by
itself contradict minimality.

Indeed, let `U_m` be uniform and let

```math
D_R=\max_{d:R_2(d)\le R}|I_d|.
```

If `D_R>0`, the exactly uniform assignment `y_S=1/D_R` is dual feasible and
has mass

```math
Y=\binom nm/D_R=\{\max_{d:R_2(d)\le R}U_m(I_d)\}^{-1}.
\tag{R38.4}
```

If `D_R=0`, the dual is unbounded.  But under the normalized law associated
with (R38.4),

```math
\boxed{
\overline B_{U_m}=(2p_2-1)A,
\qquad Q(\overline B_{U_m})=(2p_2-1)q<q
}
\tag{R38.5}
```

whenever `1/2<=p_2<1`.  Under a persistent normalized gap, (10.1022) makes
the mass in (R38.4) `exp(Omega(n^(3/4)))`, so the large dual singled out by
the converse can be perfectly selector-regular while its first moment is a
strict contraction of `A`.  Therefore a proof based only on (R38.3), pair
marginals, or commuting the weighted sum with the selector-wise maximum is
falsified as a mechanism.  Some nonlinear use of the individual `B_S` is
necessary.

## 2. Noising one high-cap complement ground

Fix a selector `S`, abbreviate `B=B_S` and `M=Q(B)`, and choose an oriented
ground `(sigma,x)` with

```math
\sigma x^TBx=M.
```

For independent signs `z_i` with `E z_i=rho`, put `a=rho^2` and
`x'=x circ z`.  Conjugating by `D_x` gives the following exact identities:

```math
\boxed{
\mathbb E[\sigma(x')^TBx']=aM,
\qquad
\mathbb E R_2(x')=aR_2(x)+(1-a)n(n-1).
}
\tag{R38.6}
```

Let `C=sigma D_xBD_x`.  Writing `z=rho 1+xi`, all mixed linear-quadratic
covariances vanish, and independence gives the exact variance

```math
\boxed{
\operatorname{Var}(z^TCz)
=4a(1-a)\lVert C\mathbf1\rVert_2^2
+2(1-a)^2n(n-1).
}
\tag{R38.7}
```

The minimizer spectral estimate and the structured form of the complement
flip give

```math
\lVert A\rVert_{op}^2\le2q,
\qquad
\lVert B\rVert_{op}
\le\lVert A\rVert_{op}+2\lVert P_SAP_S\rVert_{op}
\le3\lVert A\rVert_{op}.
```

Therefore

```math
R_2(x)\le2nq,
\qquad
\lVert C\mathbf1\rVert_2^2=\lVert Bx\rVert_2^2\le18nq.
\tag{R38.8}
```

Choose

```math
a=\frac{4q}{M}.
\tag{R38.9}
```

For the high-cap range below, `a<1`.  Equations (R38.6)--(R38.8) imply

```math
\mathbb E[\sigma(x')^TBx']=4q,
\qquad
\mathbb E R_2(x')\le2anq+n^2,
\qquad
\operatorname{Var}(\sigma(x')^TBx')\le72anq+2n^2.
\tag{R38.10}
```

These formulas retain the exact complement-flip sign: the event
`sigma(x')^TBx'>=q` is precisely `S in I_(sigma,x')`.

## 3. A large typical sphere contains many row-good incidences

Here is a convenient statement with explicit constants.  Use the established
uniform signing lower bound

```math
q_n\ge\kappa n^{3/2}
\tag{R38.11}
```

for one absolute `kappa>0` and all sufficiently large `n`.  Suppose

```math
16n^2\le R\le\min\{8nq,q^2/4\}.
\tag{R38.12}
```

The project choice `R=Theta(n^(9/4-c))`, `0<c<1/4`, satisfies these
inequalities for all sufficiently large `n` by (R38.11).  Take `K=128`.  If

```math
M\ge t\ge K\frac{nq^2}{R},
\tag{R38.13}
```

then (R38.9) gives `a<=R/(32nq)`.  Hence

```math
\mathbb E R_2(x')\le R/8,
```

and Markov gives row-failure probability at most `1/8`.  Moreover

```math
\operatorname{Var}(\sigma(x')^TBx')
\le\frac{72}{32}R+2n^2\le\frac{19}{8}R.
```

The energy mean is `4q`, so Chebyshev and `R<=q^2/4` give the explicit bound

```math
\Pr\{\sigma(x')^TBx'<q\}
\le\frac{19R}{72q^2}\le\frac{19}{288}<\frac18.
\tag{R38.14}
```

Thus at least `3/4` of the product-noise law consists simultaneously of

```math
R_2(x')\le R,
\qquad \sigma(x')^TBx'\ge q.
\tag{R38.15}
```

It remains important not to use the largest atom of the biased product law;
that would lose `exp(O(n sqrt(a)))`.  Let `K_z` be the number of flipped
coordinates.  It is binomial with mean `n(1-sqrt(a))/2`.  Hoeffding with

```math
t_n=\sqrt{(n/2)\log16}
```

shows that it lies within `t_n` of that mean with probability at least
`7/8`.  Intersecting this typical event with (R38.15) has probability at
least `5/8`.  Conditional on `K_z=k`, the product law is exactly uniform on
the radius-`k` Hamming sphere.  Averaging over the typical radii therefore
gives at least one typical integer `k` for which a fraction
`c_0>=5/8` of the **uniform Hamming sphere**

```math
\{x circ z:|\{i:z_i=-1\}|=k\}
```

consists of row-good incidences for `S`.

For such a typical `k`, binary entropy (or the method-of-types lower bound)
gives the claimed size explicitly.  Put `u=k/n`.  Then

```math
|u-1/2|\le\sqrt a/2+t_n/n,
\qquad
\binom nk\ge\frac{\exp\{nh(u)\}}{n+1}.
```

For `a<=1/4` and all sufficiently large `n`, the elementary bound
`log 2-h(u)<=C(u-1/2)^2` applies.  Hence

```math
\boxed{
\binom nk
\ge2^n\exp\{-C(na+\log n)\}
\ge2^n\exp\left\{-C\left(\frac{nq}{t}+\log n\right)\right\}.
}
\tag{R38.16}
```

The dependence is on `na`, not `n sqrt(a)`, precisely because the radius is
conditioned.  This distinction is what recovers the project exponent.

## 4. Summing the spheres against the exact dual

Let

```math
\mathcal H_t=\{S:M_S\ge t\}.
```

For every `S in H_t`, choose the sphere furnished by Section 3 and let
`mu_S` be its uniform law.  It satisfies

```math
\mu_S\{d:R_2(d)\le R,\ S\in I_d\}\ge c_0,
```

and every atom of every `mu_S` is at most the reciprocal of the right side
of (R38.16).  To make all orientation factors explicit, the actual oriented
state space has

```math
2\cdot 2^{n-1}=2^n
```

states: two orientations and projective spins `[x]=\{x,-x\}`.  Instead sum
over the redundant space

```math
\widetilde{\mathcal D}=\{+1,-1\}\times\{+1,-1\}^n,
```

which has `2^(n+1)` elements and maps exactly two-to-one onto the actual
space.  Each redundant state inherits the same incidence and the same dual
constraint as its projective image.  Each `mu_S` is supported on one fixed
orientation and one full-spin Hamming sphere.  For every feasible dual `y`,
Tonelli and the dual constraints therefore give

```math
\begin{aligned}
c_0\sum_{S\in\mathcal H_t}y_S
&\le\sum_{\widetilde d\in\widetilde{\mathcal D}}
       \sum_{\substack{S\in\mathcal H_t\\
                       S\in I_{\pi(\widetilde d)}}}
       y_S\mu_S(\widetilde d)\\
&\le\left(\max_{S,\widetilde d}\mu_S(\widetilde d)\right)
       \sum_{\widetilde d\in\widetilde{\mathcal D}}
       \sum_{S\in I_{\pi(\widetilde d)}}y_S\\
&\le 2^{n+1}\max_{S,d}\mu_S(d).
\end{aligned}
```

Here `pi(tilde d)` is the projective image of the redundant state, and every
inner sum is at most one.  Using (R38.16) proves the promised rigorous dual
tail estimate

```math
\boxed{
\sum_{S:M_S\ge t}y_S
\le
\exp\left\{C\left(\frac{nq}{t}+\log n\right)\right\},
\qquad
t\ge K\frac{nq^2}{R}.
}
\tag{R38.17}
```

All constants are universal once the conditions in (R38.11)--(R38.13) are
imposed.  Finite LP duality also turns (R38.17) into the same
upper bound for the fractional cover restricted to `H_t`.

At the desired cutoff `R=Theta(n^(9/4-c))`, the known
`q=Theta(n^(3/2))` makes (R38.12) automatic and, at the threshold in
(R38.13),

```math
\frac{nq}{t}=O(R/q)=O(n^{3/4-c}).
\tag{R38.18}
```

Thus the conditional high-cap class would meet (10.1047).  However, the
exact triangle ceiling

```math
M_S\le q+2Q(A[S])\le3q
```

is smaller than the threshold in (R38.13): already `R<=8nq` and `K=128`
force `t>=16q`.  Hence this class is empty under the theorem's stated
hypotheses and (R38.17) supplies no compression.

## 5. Cap-ceiling obstruction and scope

Define

```math
M_0=K\frac{nq^2}{R}.
```

Formally, (R38.17) would leave the statement

```math
\boxed{
\sup_y\sum_{S:Q(A^{F_S})<M_0}y_S
\le\exp\{O(R/q)\},
}
\tag{R38.19}
```

where the supremum is over the same exact row-good dual constraints.  But
the cap ceiling puts every selector in this displayed range.  Therefore
(R38.19) is just the original dual target again, not a narrower lemma.

The argument uses no assumption on `p_2`.  The high-ratio requirement
`p_2>=1/2` enters only downstream, through (10.1045), when an incidence in
`I_d` is converted into `h_d(S)<=-B_(n,m)<=0`.  Thus no sign or regime has
been hidden.

The matrices `A^(F_S)` need not be near-minimizers on the `q` scale:
`M_0/q=Theta(n^(1/4+c))`.  Nevertheless, the conditional tail theorem
removes none of them.  Any continuation must use additional exact-minimizer
structure on the full family in (R38.19), not merely its pairwise mean,
because (R38.5) falsifies that first-moment mechanism.

The threshold `M_0` is also the precise limit of this uniform-noise argument
under the currently available worst-case row information.  Retaining a
constant-factor energy mean above `q` forces `a=Theta(q/M)`.  The only
uniform ground-row estimate is `R_2(x)<=2nq`, so (R38.6) can certify the row
budget only when

```math
a\,nq=O(R),
\qquad\text{equivalently}\qquad
M=\Omega(nq^2/R).
```

Below this scale, reducing `a` enough to certify the row cap drops the mean
energy below the incidence threshold.  Since the exact cap ceiling forces
**all** complement flips below this scale, this is a sharp obstruction to
the present uniform coordinate-noise proof.  More directly, its choice
`aM=4q` is impossible because `0<=a<=1` and `M<=3q`.  It is not a
counterexample to other minimizer-specific regularizations: a specially
low-row complement ground could still do better.

## 6. Exact computation audit

`tmp/dual_flip_compression_r38_check.py` enumerates three nonsymmetric finite
signings and verifies over Python `Fraction`, with no floating-point
arithmetic:

1. both expectations in (R38.6);
2. the full variance identity (R38.7); and
3. the exact fixed-Hamming-radius pair-correlation formulas underlying
   (R38.16).

It passes under `.venv/bin/python`.  No floating-point LP or numerical
asymptotic evidence is used in (R38.11)--(R38.19).
