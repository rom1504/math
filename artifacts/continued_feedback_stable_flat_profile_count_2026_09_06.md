# A stable count for nearly flat restricted Walsh spectra

Date: 2026-09-06. A quantitative supported-spin count using half-spectrum
signs and approximate pair sums. No divisibility, algebraic-degree bound,
or general half-entropy profile theorem is assumed.

Let m=2^d, H the unnormalized Walsh matrix, and k/m=p in (0,1). For a
ternary xi of weight k write h=H xi. Call it epsilon-flat if

```math
\frac1{km}\sum_u\big(|h_u|-\sqrt k\big)^2\le\epsilon^2.   (1)
```

Let N_T(epsilon) be the number of epsilon-flat spin rows supported
exactly on a uniform k-selector T. Put q=1-p and

```math
B(p)=\left(pq+\frac{p^2}4\right)\log2
 +\left(\frac{q^2}2+\frac{p^2}4\right)
       h\!\left(\frac{q^2/2}{q^2/2+p^2/4}\right).        (2)
```

For fixed p and epsilon with `e=16p epsilon²<5/6`, the theorem is

```math
\limsup_m\frac1m\log\mathbb E_T N_T(\epsilon)
\le c_{\rm flat}(p)+\frac12\{h(e)+e\log5\},              (3)
```

```math
c_{\rm flat}(p)=\frac12\log2+B(p)-h(p)
=(p+\tfrac12)\log2-\tfrac12 H(X_1+X_2),                 (4)
```

where X1,X2 are independent with probabilities
`P(X=0)=q`, `P(X=+1)=P(X=-1)=p/2`, and H denotes discrete entropy.
The same exponential bound holds after conditioning T on any event
whose probability tends to one.

## 1. Choose a direction with a short exact recovery description

For each nonzero a, partition the m inputs into m/2 unordered pairs
`{x,x+a}`. Let n00,n++,n--,n+-,n0+,n0- be their unordered symbol-type
counts for xi. Given all pair sums and these six counts, the number of
possible xi is at most

```math
2^{n_{0+}+n_{0-}+n_{+-}}
 \binom{n_{00}+n_{+-}}{n_{00}}.                           (5)
```

Indeed sums ±2 determine the pair. Sums ±1 require one orientation bit.
A zero-sum pair is either 00 or an oppositely signed pair: choose which
zero-sum pairs are 00, then choose the orientations of the others.
The logarithm of (5) is at most

```math
\mathcal B(n)=(n_{0+}+n_{0-}+n_{+-})\log2
 +(n_{00}+n_{+-})h(n_{00}/(n_{00}+n_{+-})).                (6)
```

This is a concave homogeneous function of the counts. Average over
a!=0. If the global symbol counts are N0,N+,N-, then the average pair
count for different symbols r,s is `N_r N_s/(m-1)`, and the same-symbol
count is `N_r(N_r-1)/(2(m-1))`: every unordered pair appears in exactly
one direction. Jensen gives a direction with B(n) at most B of these
averaged counts.

Here N0=qm, N++N-=pm, and `N+ N-<=p²m²/4`. In (6), for fixed n00 and
n0++n0-, the dependence on n+- is increasing, since its derivative is
`log2+log((n00+n+-)/n+-)>0`. The balanced split therefore maximizes the
averaged upper bound. It follows that some direction has

```math
\mathcal B(n)\le m B(p)+O_p(\log m).                     (7)
```

The O(log m) safely covers the elementary finite-denominator and entropy
endpoint corrections. One can encode the direction and six counts using
only O(log m) additional bits. No assumption of globally balanced xi
was needed.

## 2. Half-spectrum signs approximately determine the pair sums

For the chosen a, let U=a-perp, |U|=M=m/2. Choose representatives for
the input quotient by {0,a}, and let s_x=xi_x+xi_(x+a). The restricted
Walsh transform satisfies EXACTLY

```math
h|_U=H_M s,\qquad s=M^{-1}H_M^T(h|_U).                  (8)
```

Encode the signs z_u=sign(h_u) on U, assigning a fixed sign at zero.
These take M bits. Replace the encoded coefficients by `sqrt(k) z_u`
and let s' be their inverse transform in (8). By (1), orthogonality gives

```math
\|s-s'\|_2^2
=M^{-1}\|h|_U-\sqrt k z\|_2^2\le2\epsilon^2 k.          (9)
```

Round every s'_x to a nearest member of {-2,-1,0,1,2}. A wrong symbol
costs at least 1/4 in (9). There are consequently at most
`8epsilon² k=eM` erroneous pair sums. Encode their locations and their
correct values; the number of possibilities is at most

```math
\sum_{j\le\lfloor eM\rfloor}\binom Mj5^j
\le\exp\{M[h(e)+e\log5]+o(m)\}.                         (10)
```

The restriction e<5/6 keeps this below the maximum of the six-symbol
counting expression. A factor four could replace five after a specific
rounding convention; the looser five is sufficient here.

Now the exact pair sums are known. Recover xi using (5), with the chosen
direction satisfying (7). The resulting description is injective once
its direction, histogram, half-spectrum signs, correction data, and
recovery index are specified. Its total size is bounded by

```math
\exp\{m[\tfrac12\log2+B(p)
              +\tfrac12(h(e)+e\log5)]+o(m)\}.            (11)
```

Each xi belongs to exactly one selector. Divide (11) by binomial(m,k)
to obtain (3). The entropy identity in (4) follows by conditioning the
ordered iid pair (X1,X2) on its sum: its conditional entropy is exactly
2B(p). This also verifies every factor of two in the count.

## 3. Direct contribution to the one-row permanent bound

Normalize a_u=|h_u|/sqrt(k), so sum a_u²=m and
sum(a_u-1)²<=epsilon²m. For any retained d=m-1 coordinates and any
permutation pi, use

```math
\log K_t(a,b)=-t(a^2+b^2)+\log\cosh(2tab).
```

The function log cosh is 1-Lipschitz. Cauchy--Schwarz gives, uniformly
over permutations and the removed coordinate,

```math
\sum|a_u^2-1|\le2\epsilon m,
\qquad\sum|a_u a_{\pi(u)}-1|\le2\epsilon m.
```

Hence the logarithm of every permutation product differs from
`d log K_t(1,1)` by at most `8t epsilon m`. Averaging, taking the square
root, and maximizing over removals yields the deterministic upper bound

```math
\frac1m\log L_t(a)
\le\tfrac12\log[(1+e^{-4t})/2]+4t\epsilon+o(1).          (12)
```

Combining (3) and (12) proves the actual class-specific partition bound

```math
\limsup_m\frac1m\log\mathbb E_T
 \sum_{x\text{ epsilon-flat}}L_t(|H[T,:]^Tx|)
\le c_{\rm flat}(p)+\tfrac12\log[(1+e^{-4t})/2]
 +\tfrac12[h(e)+e\log5]+4t\epsilon.                      (13)
```

This is a stable, quantitative profile count and permanent bound. It
does not assume that epsilon-flat rows are near exact bent functions.

## 4. Quantitative limitation at retention 15/16

At p15/16 the constant in (4) is

```text
c_flat(p) = .3168851694355724,
```

compared with the trivial p log2=.6498254817749487. Despite this genuine
saving, the best limiting tilted upper expression from (13) remains
positive. Writing delta=1-sqrt(p), its unique minimizing tilt is

```math
t_*=\tfrac14\log(2/\delta-1)=1.03171853444778\ldots,
```

and the minimum at epsilon=0 is approximately .01107517452307. At the
tilts around 4--7 left open by the Gaussian and rare-selector LOWER
obstructions, the bound is weaker still (about .16084 at t6).

Thus (13) does not yet eliminate even the nearly flat class at a tilt
suitable for the complete one-row sum. A deeper robust recovery step,
or another count, is still necessary. Unlike the earlier sparse-span
module, this theorem controls a genuinely dense profile neighborhood.
