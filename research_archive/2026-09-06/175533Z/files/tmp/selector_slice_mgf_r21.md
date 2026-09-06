# Wave 21 memo: a defective Hanson--Wright bound on the fixed slice

## Status and normalization

The log-mgf theorem below is **proved self-containedly**.  It uses independent
Bernoulli selectors only as an auxiliary law and then conditions exactly on
their sum.  The conditioning loss is explicit and at most `log(n+1)`.  No
claim that variance alone implies an mgf bound is used.

Fix a full oriented cut `d=(sigma,x)` and put, on unordered edges,

```math
w_{ij}=\sigma a_{ij}x_ix_j,\qquad
T_S=\sum_{i<j}w_{ij}1_{\{i,j\in S\}},\qquad
W=\sum_{i<j}w_{ij}.
```

Thus the ordered selector payoff is `c_A(S,d)=2T_S` and the full payoff is
`<A,d>=2W`.  Write

```math
r_i=\sum_{j\ne i}w_{ij},\quad
R_2=\sum_i r_i^2,\quad R_\infty=\max_i|r_i|,
```

and let `\mathsf W=(w_{ij})` be the symmetric zero-diagonal matrix.  Finally,

```math
F_\times(\mathsf W)
=\max_{L\subset[n]}\|\mathsf W[L,L^c]\|_F^2.
```

For a signing, `F_x <= floor(n^2/4)` and
`||mathsf W||_op=||A||_op`.

## 1. Independent quadratic mgf

Let `xi_i` be independent Bernoulli(`p`) variables,
`eta_i=xi_i-p`, and `s=p(1-p)`.  Then

```math
T_\xi-p^2W
=p\sum_i r_i\eta_i+\sum_{i<j}w_{ij}\eta_i\eta_j
=:L+Q.
```

For every real `t` satisfying

```math
|t|<\min\left\{\frac3{2pR_\infty},
                 \frac1{\|\mathsf W\|_{op}}\right\},
```

with the usual infinite interpretation when a denominator is zero,

```math
\boxed{
\log\mathbb E e^{t(T_\xi-p^2W)}
\le
\frac{t^2p^2sR_2}{1-2|t|pR_\infty/3}
+\frac{t^2F_\times(\mathsf W)}
       {4(1-t^2\|\mathsf W\|_{op}^2)}.}
```

### Linear term

The elementary Bernstein mgf inequality

```math
\log\mathbb E e^{uX}
\le\frac{u^2\mathbb EX^2}{2(1-|u|b/3)}
\qquad(|X|\le b,\ \mathbb EX=0)
```

follows by expanding the exponential and using
`E|X|^k <= b^{k-2}EX^2` and `k! >= 2*3^{k-2}` for `k>=2`.
Applied independently to `X_i=pr_i eta_i`, it gives

```math
\log\mathbb E e^{uL}
\le\frac{u^2p^2sR_2}{2(1-|u|pR_\infty/3)}.
```

### Quadratic term

Choose a random vertex set `L_0` by independent fair coins and put

```math
C_{L_0}=\sum_{i\in L_0,j\notin L_0}w_{ij}\eta_i\eta_j.
```

Every unordered edge crosses with probability `1/2`, so
`Q=2 E_{L_0} C_{L_0}`.  Jensen gives

```math
e^{uQ}\le\mathbb E_{L_0}e^{2uC_{L_0}}.
```

For fixed `L_0`, write `B=mathsf W[L_0,L_0^c]`.  Conditioning first on
`eta_{L_0}` and using Hoeffding's lemma for the centered variables
`eta_j`, whose range has length one, gives

```math
\mathbb E_{\eta_{L_0^c}}e^{2uC_{L_0}}
\le \exp\left\{\frac{u^2}{2}\|B^T\eta_{L_0}\|_2^2\right\}.
```

For `alpha=u^2/2`, Gaussian linearization and Hoeffding once more give

```math
\begin{aligned}
\mathbb E_\eta e^{\alpha\|B^T\eta\|_2^2}
&=\mathbb E_g\mathbb E_\eta
  e^{\sqrt{2\alpha}\langle Bg,\eta\rangle}\\
&\le\mathbb E_g e^{(\alpha/4)\|Bg\|_2^2}\\
&=\prod_k(1-\alpha s_k(B)^2/2)^{-1/2}\\
&\le\exp\left\{
\frac{u^2\|B\|_F^2}
     {8(1-u^2\|B\|_{op}^2/4)}\right\}.
\end{aligned}
```

Here the last step uses `-log(1-z)<=z/(1-z_max)`.  Since a rectangular
compression satisfies `||B||_op<=||mathsf W||_op`, this proves

```math
\log\mathbb E e^{uQ}
\le
\frac{u^2F_\times(\mathsf W)}
     {8(1-u^2\|\mathsf W\|_{op}^2/4)}.
```

Cauchy--Schwarz,
`E exp(t(L+Q)) <= (E exp(2tL) E exp(2tQ))^{1/2}`, now gives the claimed
independent bound with exactly the displayed constants.

## 2. Conditioning to the uniform slice

Let `S` be uniform among the `m`-sets and set

```math
p=\frac mn,\qquad
p_2=\frac{(m)_2}{(n)_2},\qquad
\epsilon_{n,m}=p^2-p_2
=\frac{m(n-m)}{n^2(n-1)}.
```

For `K=sum_i xi_i`, the conditional law of `{i:xi_i=1}` given `K=m` is
exactly uniform on the slice.  Put

```math
\chi_{n,m}
=-\log\Pr\{\operatorname{Bin}(n,m/n)=m\}.
```

The value `m` is a mode of this binomial law, so
`chi_{n,m} <= log(n+1)` because its `n+1` point masses sum to one.
For every admissible real `t`, conditioning and then taking the absolute
value only in the deterministic centering correction gives

```math
\begin{aligned}
\log\mathbb E_{U_m}e^{t(T_S-p_2W)}
&\le \chi_{n,m}+|t|\epsilon_{n,m}|W|\\
&\quad+\frac{t^2p^2sR_2}{1-2|t|pR_\infty/3}
+\frac{t^2F_\times(\mathsf W)}
       {4(1-t^2\|\mathsf W\|_{op}^2)}.
\end{aligned}
```

Indeed, before taking absolute values the correction is exactly
`t(p^2-p_2)W`; its sign is not assumed.

For the ordered payoff and `lambda>0`, this is the fixed-slice theorem

```math
\boxed{
\begin{aligned}
\Lambda_d^U(\lambda)
&\le \chi_{n,m}
+\lambda\epsilon_{n,m}|\langle A,d\rangle|\\
&\quad+
\frac{4\lambda^2p^2sR_2}
     {1-4\lambda pR_\infty/3}
+\frac{\lambda^2F_\times(\mathsf W)}
     {1-4\lambda^2\|A\|_{op}^2},
\end{aligned}}
```

valid for

```math
0<\lambda<\min\left\{
\frac3{4pR_\infty},\frac1{2\|A\|_{op}}
\right\}.
```

The row-field coefficient `p^2s=p^3(1-p)` agrees to leading order, at
fixed `m/n`, with the exact slice coefficient `p_3-p_4` in (10.649).  The
theorem does not reproduce the exact finite variance: `chi` is the explicit
price of the simple conditioning proof.

## 3. Exact-minimizer corollary and the information budget

For an exact minimizer, (10.650) and (10.67) give

```math
R_2\le2(n-1)q_n,qquad
R_\infty\le n-1,qquad
\|A\|_{op}^2\le2q_n,qquad
|\langle A,d\rangle|\le q_n.
```

Restrict further to

```math
\lambda\le
\min\left\{\frac3{8pR_\infty},
             \frac1{2\sqrt2\|A\|_{op}}\right\}.
```

Both denominators are then at least `1/2`, and uniformly in `d`,

```math
\boxed{
\Lambda_d^U(\lambda)
\le\chi_{n,m}+\lambda\epsilon_{n,m}q_n
+\lambda^2\overline V_{n,m},
\qquad
\overline V_{n,m}
=16p^2s(n-1)q_n+\frac{n^2}{2}.}
```

Here the sharper cut-dependent coefficient is
`8p^2sR_2+2F_x <= overline V`; using the uniform ceiling in the optimizer
also guarantees that the proposed `lambda` is small enough even when a
particular cut has unusually small field energy.

Let

```math
J_{n,m}(\delta,\pi)
=R_\pi(\delta)+D_{KL}(\pi\Vert U_m)+\chi_{n,m}.
```

If

```math
\lambda_*=\sqrt{J_{n,m}/\overline V_{n,m}}
```

lies in the preceding domain, substitution in (10.671) proves

```math
\boxed{
V_{ad}(A,m)-p_2q_n
\le
\delta+\epsilon_{n,m}q_n
+2\sqrt{\overline V_{n,m}J_{n,m}(\delta,\pi)}.}
```

This explicitly optimizes the auxiliary parameter.  If `m/n` stays in a
fixed compact subinterval of `(0,1)` and `q_n=O(n^{3/2})`, then

```math
\overline V_{n,m}=O(n^{5/2}),\qquad
\epsilon_{n,m}q_n=O(n^{1/2}),\qquad
\chi_{n,m}=O(\log n).
```

Consequently:

* `delta=o(n^{3/2})` and
  `R_pi(delta)+D_KL(pi||U_m)=o(sqrt(n))` give an `o(n^{3/2})`
  revelation cost.  Here `lambda_*=o(1/n)`, so the required domain holds.
* Quantitatively, if for some `0<c<1/4`,

```math
\delta=O(n^{3/2-c}),\qquad
R_\pi(\delta)+D_{KL}(\pi\Vert U_m)=O(n^{1/2-2c}),
```

  then the right side is `O(n^{3/2-c})`.  The conditioning loss
  `sqrt(n^{5/2}log n)=n^{5/4}sqrt(log n)` fits this rate for every
  fixed `c<1/4`.

The latter statement, uniformly on any fixed-ratio comparison subwindow
whose ratios stay away from both zero and one, proves the corresponding
power-saving restriction edges.  Such fixed-ratio edges are enough for the
adaptive comparison framework when supplemented by an exact landing scheme;
the stronger all-`m` steering formulation, especially `n-m=o(n)`, is not
proved by this estimate and would need a separate near-top argument.  Wave 21
therefore closes the fixed-slice mgf half of the information interface at the
cost of the concrete remaining entropy target

```math
R_\pi(\delta)+D_{KL}(\pi\Vert U_m)=O(n^{1/2-2c}).
```

It does not prove that rate--distortion estimate.  The use of the uniform
reference in (10.671) means that a very concentrated nonuniform selector can
lose through `D_KL(pi||U_m)`; choosing `pi=U_m` removes that term but asks for
a low-information near-ground channel for the uniform slice.

## 4. Scope

The bound is deliberately called defective because of `chi_{n,m}` and
because its quadratic scale uses `F_x` and `||A||_op`, rather than the exact
degenerate variance in (10.649).  Those defects are power-saving at the
information scale above.  A direct slice martingale or slice log-Sobolev
argument could improve constants, but is no longer needed for the stated
rate target.

No pressure identity currently yields the remaining `O(n^{1/2-2c})`
shared-prior information.  Fixed-temperature child Gibbs channels have
`O(n)` distortion but their overlap across selectors is uncontrolled.  Thus
the next leading target is the rate--distortion/shared-prior half, not another
fixed-state variance calculation.

This theorem applies directly to the uniform-reference inequality (10.671).
It does not bound the law-adapted log-mgf in (10.670) for a nonuniform hidden
optimizer `pi_*`.  One must either use (10.671) and pay
`D_KL(pi_*||U_m)`, choose `pi=U_m`, or prove a separate comparison between
the nonuniform and uniform selector laws.
