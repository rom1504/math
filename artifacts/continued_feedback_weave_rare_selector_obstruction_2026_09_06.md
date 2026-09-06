# Rare selectors obstruct the Gaussian-optimal weave tilt

Date: 2026-09-06. Original one-row partition-sum LOWER bounds. These do
not obstruct a selector law conditioned to exclude the rare selectors,
and they do not prove an upper cap or convergence theorem.
The exact finite-spike proof, the nearby-selector/flip envelope, and the
repeated-fibre formulas have passed independent audit. Numerical profile
optimizations below remain uncertified.

Use the square-root permanent `L_t` and `Z_T(t)` of
`continued_convergence_restricted_weave_2026_09_06.md`, with a Walsh
Hadamard matrix of order m and a uniform k-selector, k/m=p. All logarithms
are natural; `h(p)=-p log p-(1-p)log(1-p)`.

## 1. An exact finite-spike obstruction

Fix a dyadic rational p=a/b in (0,1), where b is a fixed power of two.
Let m range over powers of two divisible by b. Choose a selector T0
which is the inverse image of an a-element subset of the b-point
quotient under a fixed linear quotient map. Its indicator has Fourier
support in the b-point annihilator U. With the all-plus spins on T0,
the length-m spectrum consequently has at most b nonzero coordinates.

For ANY such nonnegative spectrum, delete any one coordinate and put
d=m-1. Let r<=b be the number of remaining nonzero coordinates. In the
normalized permanent, retain only permutations fixing these r labeled
positions. The other d-r positions are zeros and can be permuted freely.
Since `K_t(v,v)=(1+exp(-4tv²/k))/2>=1/2` and `K_t(0,0)=1`, exactly

```math
\frac{\operatorname{per}K_t[a\setminus v]}{d!}
\ge 2^{-r}\frac{(d-r)!}{d!}\ge 2^{-b}d^{-b}.             (1)
```

This holds at every t>=0, independent of the size of the nonzero spikes.
The probability of this single selector is `1/binomial(m,pm)`. Therefore

```math
\boxed{\liminf_m\frac1m\log\mathbb E_T Z_T(t)\ge-h(p)}.  (2)
```

The maximum over diagonal removals causes no difficulty; (1) already
holds for every deletion. Counting more quotient selectors or signed
columns is unnecessary for the exponent.

The target strict negative exponent is thus impossible whenever

```math
t(1-\sqrt p)\ge h(p).                                   (3)
```

At p=15/16,

```text
h(p) = 0.2337916587064593
h(p)/(1-sqrt(p)) = 7.362551341912657
Gaussian-optimal t = sqrt(p)/(2(1-p)) = 7.745966692414834.
```

Thus the previously favorable Gaussian-optimal tilt is rigorously
falsified by actual rare selectors. This does not rule out smaller t:
the typical-Gaussian obstruction and (2) leave a nonempty window.
Since (1) is uniform in t, the same class also rules out every divergent
tilt sequence t_m->infinity for this unconditioned exact-dyadic setup:
the positive tilt penalty eventually exceeds h(p), without requiring
any uniform-in-t Gaussian approximation.

## 2. A nearby-selector lower envelope

The finite-spike argument extends to exponentially many nearby selectors.
Write q=1-p and fix 0<=alpha<=pq. Starting from T0, retain
`(p-alpha)m` of its pm points and add `alpha m` points of its complement.
The fraction of all k-selectors with this overlap has exponential rate

```math
-I_p(\alpha),\qquad
I_p(\alpha)=h(p)-p h(\alpha/p)-q h(\alpha/q).              (4)
```

For aligned-column spins, the spectral coefficients outside U have zero
mean and their variance, after division by k, tends to

```math
v_p(\alpha)=\frac{2\alpha-\alpha^2/(pq)}p.                (5)
```

Indeed the indicator variances in the two populations sum to
`p (alpha/p)(1-alpha/p)+q (alpha/q)(1-alpha/q)` times m.
For any frequency outside U its signs sum to zero separately in both
populations. Finite-population central limit theorems therefore give the
stated normal limit. This can also be proved directly by decomposing
each population into its two character-sign classes and applying the
hypergeometric central limit theorem.

For a pair of frequencies in different cosets of U, the corresponding
four character-sign classes have the balanced population sizes. Their
joint limit is independent. Only O(bm) ordered pairs lie in the same
coset, so bounded empirical spectral tests have variance tending to zero.
Fourth moments stay bounded at fixed alpha>0. Consequently the empirical
bulk magnitude profile converges with quadratic-tail control to
`|N(0,v_p(alpha))|`. Frequencies in U are exceptional and can all be fixed
in the permanent, at the polynomial cost in (1).

The same finite-type lower argument used for the Gaussian profile then
gives

```math
\liminf_m\frac1m\log\mathbb E_T Z_T(t)
\ge -I_p(\alpha)-\tfrac12F(t v_p(\alpha)),                (6)
```

where `F(s)=2s(1-rho)-log(1-rho²)/2` and
`2s=rho/(1-rho²)`. At alpha=0 use (1) directly; at alpha=pq this is
the aligned-column Gaussian-deletion profile of variance q.

One can additionally flip a proportion theta in [0,1/2] of retained
spins, fixing asymptotically that proportion in each population. There
are `exp{mp h(theta)+o(m)}` such spin choices per selector. The finite
population means are multiplied by `mu=1-2theta`; its bulk variance is
therefore `1-mu²(1-v_p(alpha))`. The same argument yields

```math
\liminf_m\frac1m\log\mathbb E_T Z_T(t)
\ge\sup_{\alpha,\theta}\left\{
 p h(\theta)-I_p(\alpha)
 -\tfrac12F\big(t[1-(1-2\theta)^2(1-v_p(\alpha))]\big)
 \right\}.                                             (7)
```

These are lower bounds, not upper evaluations of the full partition sum.
Rational parameter approximations handle integer population counts.

An uncertified grid optimization at p15/16 gives a best remaining
maximum lower exponent near -0.03745 at t6.16066, after adding the tilt
penalty. Thus this stronger family still does not falsify EVERY tilt.
The exact conclusion requiring no numerical optimization is (2)--(3).

## 3. A different scalable family does not currently obstruct the tilt

For comparison, make each spin row constant on b-point input fibres,
b a fixed power of two, and allow independent signs on the fibres hit
by T. A typical random selector hits a fraction `1-q^b` of them, giving
row-count rate `(1-q^b)log(2)/b`. The normalized bulk spectral profile is

```math
(1-1/b)\,|N(0,q)|+(1/b)\,|N(0,q+bp)|.                    (8)
```

To verify (8), in the zero internal-frequency class the random fibre
occupancy has second moment `bpq+b²p²`; in each nonzero class its signed
occupancy has second moment bpq. Divide by k=pm and use the independent
fibre signs. Distinct outer frequencies have asymptotically independent
normal limits except for O(m) same-outer-frequency pairs, and bounded
fourth moments again give empirical convergence. The exact count uses
only hit fibres, not all m/b fibre signs.

`computations/continued_feedback_weave_mixture_probe_2026_09_06.py`
evaluates the associated profile variational problem by Gaussian
quadrature and positive symmetric Sinkhorn scaling. The outputs are
UNCERTIFIED approximations, not upper enclosures. At p15/16 and the
Gaussian-optimal tilt, the tilted exponents for b=2,4,8,16,32,64 are
approximately

```text
-.26250, -.28814, -.24067, -.18049, -.13237, -.10267.
```

They do not presently reveal an obstruction. In contrast, Section 1
rigorously obstructs that same tilt through rare SELECTORS rather than
the typical-selector row families.
