# Wave 46A: noncentral complement slack and the arithmetic/geometric wall

## Status

- **Verified:** the exact fixed-layer mean and variance of the complement
  score, and the sharp one-sided second-moment survival bound below.
- **Verified:** this survival bound controls only the arithmetic incidence
  retention.  It cannot imply the geometric/logarithmic Harnack premise in
  (10.1159).
- **Verified conditional advance:** at an exact scalar optimizer, arithmetic
  retention can instead be combined with a centered fixed-slice row mgf.  If
  an `exp{-O(n^(3/4-c))}` incidence fraction has enough active slack, this
  gives the required `R_2=O(n^(9/4-c))` row bound without any geometric
  Harnack assumption.
- **Falsified as a literal finite mechanism:** on every scalar Pareto column
  of `A_8,m=6` and `A_9,m=7`, every nontrivial fixed block layer contains an
  empty column.  Hence the left side of the literal logarithmic Harnack
  premise is infinite at every `k=1,...,n-1`.  The project-price counts were
  independently reproduced.
- **Open:** no asymptotic exact-minimizer theorem supplies the required
  incidence-conditioned high slack, and no scalable exact-minimizer family
  falsifying it was found.  Thus the noncentral arithmetic/mgf route is
  narrowed to one explicit lemma, not closed.

## 1. Exact fixed-layer law

Fix an active pair `(d,S)`, put

```math
L=L_S(d)=q+\Gamma_S,
\qquad R_B=R_{B_S}(d)=\lVert B_Sx\rVert_2^2,
```

and let `U` be uniform among the `k`-subsets.  Write

```math
z=n-2k,
\qquad
\theta=\frac{z^2-n}{n(n-1)}
=1-\frac{4k(n-k)}{n(n-1)},
```

and

```math
\mu_4
=\frac{z^4-(6n-8)z^2+3n^2-6n}
       {n(n-1)(n-2)(n-3)}.
```

Gauging `B_S` by `d`, and separating equal, adjacent, and disjoint
unordered edge pairs, gives the exact identities

```math
\boxed{\mathbb E_U L_S(d^U)=\theta L,}
```

```math
\boxed{
V_S(k):=\operatorname {Var}_U L_S(d^U)
=(\mu_4-\theta^2)L^2
+4(\theta-\mu_4)R_B
+2(1-2\theta+\mu_4)n(n-1).}
\tag{R46.A1}
```

This includes the finite-population correction which is lost by replacing
the layer with independent biased signs.  For reference,

```math
\theta-\mu_4
=\frac{(n^2-z^2)(z^2-(n-2))}
       {n(n-1)(n-2)(n-3)}.
```

Also `mu_4-theta^2<=0` throughout the only mean-survival regime
`theta>=1/3`.  Indeed its numerator, over the positive common denominator,
is

```math
2(2n-3)(z^2-n^2)
\left(z^2-\frac{n(n-2)}{2n-3}\right).
```

Consequently, when `theta>=1/3`, the convenient rigorous upper bound is

```math
V_S(k)
\le 4(\theta-\mu_4)R_B
   +2(1-2\theta+\mu_4)n(n-1).
\tag{R46.A2}
```

For fixed `k/n -> kappa`, put `a=1-2kappa`.  Its leading form is

```math
V_S(k)
=4a^2(1-a^2)R_B+2(1-a^2)^2n^2
+O\!\left(\frac{L^2+R_B+n^2}{n}\right).
\tag{R46.A3}
```

Define the contracted mean surplus

```math
\delta_S(k)
=\theta(q+\Gamma_S)-q
=\theta\Gamma_S-(1-\theta)q.
\tag{R46.A4}
```

Cantelli's inequality now gives the exact moment-only conclusions

```math
\delta_S(k)>0
\quad\Longrightarrow\quad
\Pr_U\{L_S(d^U)\ge q\}
\ge\frac{\delta_S(k)^2}{V_S(k)+\delta_S(k)^2},
\tag{R46.A5}
```

and, with `g_S=q-theta L>0`,

```math
\Pr_U\{L_S(d^U)\ge q\}
\le\frac{V_S(k)}{V_S(k)+g_S^2}.
\tag{R46.A6}
```

These are sharp given only the displayed mean and variance.  Equation
(R46.A5), unlike an upper-tail estimate, is a genuine survival lower bound.

The slack cost is already severe.  Since

```math
|L_S(d)|\le Q(B_S)\le Q(A)+2Q(A[S])\le3q,
```

one has `Gamma_S<=2q`.  Thus (R46.A5) can be nontrivial only for
`theta>1/3`.  At a fixed ratio, if `g=Gamma_S/q` and
`tau=1-a^2=4kappa(1-kappa)`, positive contracted surplus is equivalent to

```math
\boxed{\tau<\frac{g}{1+g}.}
\tag{R46.A7}
```

The left side is exactly the asymptotic row-contraction factor in
(10.1159).  Hence for a given robust slack fraction `g`, the best possible
mean-survival optimization of that denominator is

```math
\frac1{1-\theta}\gtrsim\frac{1+g}{g}.
\tag{R46.A8}
```

In particular, any fixed macroscopic row contraction requires a fixed
positive fraction of the full `q` scale as active slack.  Taking `k/n` close
to zero reduces the slack tax, but loses the same factor in the row
contraction.

## 2. Arithmetic retention is not logarithmic Harnack

Let `alpha_d=U_m(I_d)` and condition `S` on `I_d`.  Averaging (R46.A5) gives

```math
\overline p_k(d)
:=\mathbb E_{S\mid I_d}
\left[
\mathbf1_{\{\delta_S(k)>0\}}
\frac{\delta_S(k)^2}{V_S(k)+\delta_S(k)^2}
\right],
```

and double counting proves only

```math
\boxed{
\mathbb E_{|U|=k}\alpha_{d^U}
\ge\alpha_d\,\overline p_k(d).}
\tag{R46.A9}
```

This is an arithmetic-mean statement.  Jensen has the wrong direction for
the desired estimate:

```math
\mathbb E_U[-\log\alpha_{d^U}]
\ge-\log\mathbb E_U\alpha_{d^U}.
```

There is no upper bound on the left side, and one empty orbit point makes it
infinite.  Therefore (R46.A9) cannot be inserted into (10.1159).  Any claim
that Cantelli survival proves the old Harnack premise is false.

## 3. A scalar-optimality replacement for geometric Harnack

There is, however, a rigorous way to use the arithmetic statement.  Let `d`
minimize

```math
F_\lambda(d)=-\log\alpha_d+\lambda R_2(d)
```

over nonempty columns.  Scalar optimality gives, including the trivial case
`alpha_(d^U)=0`,

```math
\frac{\alpha_{d^U}}{\alpha_d}
\le\exp\{\lambda(R_2(d^U)-R_2(d))\}.
```

Combining this pointwise inequality with (R46.A9), and writing
`D=n(n-1)`, yields

```math
\boxed{
\lambda(1-\theta)(R_2(d)-D)
\le \Psi_k(\lambda;d)+\log\frac1{\overline p_k(d)},}
\tag{R46.A10}
```

where

```math
\Psi_k(\lambda;d)
=\log\mathbb E_{|U|=k}
 \exp\{\lambda(R_2(d^U)-\mathbb E R_2(d^U))\}.
```

This is the correct arithmetic analogue of (10.1159).  It permits empty
columns and pays a centered row mgf instead of an average logarithm.

That mgf is usable.  Put

```math
C=\operatorname {diag}(x)A^2\operatorname {diag}(x),
\qquad C^\circ=C-(n-1)I,
\qquad u=C^\circ\mathbf1.
```

Let `p=k/n`, `a=1-2p`, and let
`chi_(n,k)=-log Pr{Bin(n,p)=k}<=log(n+1)`.  Center independent Bernoulli
flips and condition on their sum.  Hoeffding for the linear part and the
same centered quadratic-form mgf used in the proof of (10.690) give, for
universal `c,C>0` and
`0<lambda<=c/||C^circ||_op`,

```math
\boxed{
\begin{aligned}
\Psi_k(\lambda;d)
&\le \chi_{n,k}
+\lambda(a^2-\theta)|R_2(d)-D|\\
&\quad+
\frac{C\lambda^2
\{a^2\lVert u\rVert_2^2+\lVert C^\circ\rVert_F^2\}}
{1-C\lambda\lVert C^\circ\rVert_{\rm op}}.
\end{aligned}}
\tag{R46.A11}
```

For completeness, under the independent law the centered row is exactly

```math
-4a\sum_i u_i(\xi_i-p)
+8\sum_{i<j}C^\circ_{ij}(\xi_i-p)(\xi_j-p),
```

which proves the claimed linear/quadratic proxies.  The difference between
the independent and slice means is
`(a^2-theta)(R_2-D)`, and conditioning costs `chi_(n,k)`.

For an exact minimizer signing, the standard operator bound
`||A||_op^2<=2q` gives

```math
\begin{aligned}
\lVert C^\circ\rVert_{\rm op}&\le2q+n-1,\\
\lVert C^\circ\rVert_F^2&\le\operatorname {tr}A^4
\le2qn(n-1),\\
\lVert u\rVert_2^2&\le4qR_2(d)+2n(n-1)^2.
\end{aligned}
\tag{R46.A12}
```

The last line uses the biased first-chaos identity and bound

```math
\sum_i\{x_i(A^2x)_i\}^2=x^{\mathsf T}A^4x
\le\lVert A\rVert_{\rm op}^2x^{\mathsf T}A^2x
\le2qR_2(d),
```

followed by `||v-(n-1)1||^2<=2||v||^2+2n(n-1)^2`.

Fix `k/n -> kappa` with `0<kappa<(1-1/sqrt(3))/2`; the reflected range is
equivalent.  Let `q<=Q_0 n^(3/2)` for the verified universal constant
`Q_0`, and choose `lambda=b_kappa n^(-3/2)` with fixed `b_kappa>0` so small
that, for the universal constant in (R46.A11),

```math
C b_\kappa(2Q_0+1)<\frac12,
\qquad
\frac{4C b_\kappa Q_0(1-2\kappa)^2}
     {4\kappa(1-\kappa)}<\frac14.
\tag{R46.A12b}
```

The first condition is the quadratic-mgf domain and the second absorbs the
biased first-chaos `R_2` term.  These constants depend only on the fixed
layer ratio, not on `n`.
The `R_2` terms in (R46.A11) can then be absorbed into the left side of
(R46.A10).  The remaining quadratic term is `O(n^(1/2))`.  Thus

```math
\boxed{
R_2(d)
\le C_\kappa\left[
n^2+n^{3/2}\log\frac1{\overline p_k(d)}
\right].}
\tag{R46.A13}
```

In particular, for `0<c<1/4` and `H=n^(3/4-c)`,

```math
-\log\overline p_k(d)=O(H)
\quad\Longrightarrow\quad
R_2(d)=O(n^{9/4-c}).
\tag{R46.A14}
```

This is a minimizer-relevant project-scale sufficient lemma.  It still does
not prove the separate mass estimate `-log alpha_d=O(H)`.

There is a simpler exact hypothesis implying (R46.A14).  Since
`||B_S||_op<=3||A||_op`, (R46.A2) gives `V_S(k)=O(nq)=O(n^(5/2))` in this
regime.  A strictly positive `delta_S(k)` is a rational with denominator at
most `n(n-1)`, so its Cantelli factor is at least `n^(-7)` for all large
`n`.  Therefore it is enough to prove

```math
\boxed{
\Pr_{S\mid I_d}
\{\theta\Gamma_S>(1-\theta)q\}
\ge\exp\{-O(H)\}.}
\tag{R46.A15}
```

This high-slack incidence criterion, at a scalar optimizer which also has
the required mass, is the exact surviving noncentral target.

## 4. The central entropy wall extends to an `n^(7/8)` band

The conditioning argument of Wave 45 extends farther than the previously
recorded `O(sqrt(n))` band.  Fair iid flips make `d^U` uniform, so (10.1168)
and conditioning on `|U|=k` give

```math
\mathbb E_{|U|=k}\alpha_{d^U}
\le C\exp\left{-cn^{3/4}
-\log\Pr\{\operatorname {Bin}(n,1/2)=k\}\right}.
```

Equivalently, with `kappa=k/n`,

```math
\boxed{
\mathbb E_{|U|=k}[-\log\alpha_{d^U}]
\ge cn^{3/4}-nD(\kappa\Vert1/2)-O(\log n).}
\tag{R46.A16}
```

Here
`-log Pr{Bin(n,1/2)=k}=nD(kappa||1/2)+O(log n)`.  If
`|k-n/2|=o(n^(7/8))`, then

```math
nD(\kappa\Vert1/2)
=\frac{2(k-n/2)^2}{n}
+O\!\left(\frac{(k-n/2)^4}{n^3}\right)
=o(n^{3/4}).
```

Thus the `Omega(n^(3/4))` logarithmic wall holds throughout that band.  A
small constant multiple of `n^(7/8)` also works after shrinking the constant
using the `c` in (10.1168).  Any layer escaping this particular entropy
argument must have

```math
|1-2k/n|\not=o(n^{-1/8}).
```

The mean-survival criterion (R46.A7) is much more restrictive: because
`Gamma<=2q`, it requires `|1-2k/n|>1/sqrt(3)+o(1)`.

## 5. Exact finite obstruction

Exhaustive enumeration at the project price `lambda=n^(-3/2)` gives the
following empty-column counts for the displayed representative; the profiles
are symmetric under `k <-> n-k`.

| signing | active slacks | empty counts by `k=1,...,floor(n/2)` |
|:--|:--|:--|
| `A_5,m=4` | `0,0` | `1/5, 0/10` |
| `A_6,m=5` | `0,0,0` | `6/6, 0/15, 20/20` |
| `A_8,m=6` | `0,0` | `5/8, 6/28, 31/56, 16/70` |
| `A_9,m=7` | `4,4` | `4/9, 23/36, 48/84, 70/126` |

More strongly, the complete scalar Pareto types are

```text
A8: (R_2,alpha,slacks)=(40,1/14,(0,0)), (72,3/14,(0,...,0));
A9: (56,1/36,(4)), (72,1/18,(4,4)), (88,1/12,(4,4,4)).
```

Every one of these `A_8,A_9` Pareto columns has at least one empty column in
every nontrivial fixed layer.  Every scalar optimizer at positive price is a
Pareto column, so literal use of (10.1159) fails on every layer, not just the
central one.  In addition, (R46.A15) has zero left side on every listed
Pareto type for every nontrivial `k`: the active slack never pays even the
mean-contraction tax.

These are exact finite mechanism counterexamples, not a scalable
exact-minimizer obstruction.  They show that neither nonemptiness of the
full block orbit nor robust high slack follows formally from scalar
optimality and exact signing minimality.

## Verification

`tmp/noncentral_slack_r46_check.py` independently:

1. enumerates all oriented columns and scalar optimizers of
   `A_5,A_6,A_8,A_9`;
2. verifies (R46.A1) by exact rational enumeration for every active selector
   of a representative optimizer and every `k`;
3. verifies the Cantelli lower bound against the exact retained-pair count;
4. enumerates all scalar Pareto types and their robust-slack counts; and
5. checks the empty-column layer profiles quoted above.

It terminates with

```text
PASS noncentral_slack_r46_check
```
