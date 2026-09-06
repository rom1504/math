# Wave 49A: complete project-row non-strict coarea package

## Status

No asymptotic coarea theorem or convergence proof is obtained.  The main
advance is an exact reduction of the whole coarea clause to a cubic
high-harmonic inequality, equivalently a triple-ground retention theorem,
and a sharp one-quantity sufficient condition for the independent project-row
mass clause.  Exact finite calculations expose two constraints which a proof
must respect:

- prescribed near-core coarea is false on the exact minimizer `A9` and on all
  twelve sampled exact order-ten minimizers; and
- row truncation is not monotone for the normalized boundary.

The failures are finite, not a scalable family.  More nonlocal core scales
still work in every audited matrix, so the complete project-row route remains
open in the adaptive-scale form stated in Section 5.

All algebra and rational values are checked by
`tmp/coarea_package_r49_check.py`.  Its compact output is
`tmp/coarea_package_r49.out`.

## 1. Exact cubic and triple-retention reformulations

Fix a selector-independent base law on projective centers `z`, and let

```math
C=\{z:R_2(z)\le R_*\},\qquad
f_z(S)=\mathbf 1\{D_z(S)=0\},\qquad a_z=\mathbb E_S f_z(S).
```

The same formulas hold for a positive hard threshold after replacing the
zero-deficit indicator by its threshold indicator.  Let `K=K_ell` be the
Johnson common-core kernel, with eigenvalues `1=lambda_0 > lambda_1 >=
lambda_2 >= ...`, put `delta=1-lambda_1`, and set

```math
B_z=\langle f_z,(I-K)f_z\rangle,
\qquad
\Xi_z=\sum_{j\ge2}(\lambda_1-\lambda_j)
\lVert P_jf_z\rVert_2^2.
```

Exact spectral accounting gives

```math
B_z=\delta(a_z-a_z^2)+\Xi_z.
```

Therefore the complete non-strict coarea target is **exactly equivalent** to
the cubic high-harmonic inequality

```math
\boxed{
\mathbb E[\mathbf1_Ca_zB_z]\le\delta\mathbb E[\mathbf1_Ca_z^2]
\quad\Longleftrightarrow\quad
\mathbb E[\mathbf1_Ca_z\Xi_z]\le
\delta\mathbb E[\mathbf1_Ca_z^3].}
\tag{R49A.1}
```

There is a more concrete probabilistic form.  Conditional on a base center
`z`, sample `S_0` uniformly and independently sample a stationary Johnson
pair `(S_1,T)` with `T~K(S_1)`.  Since

```math
\mathbb P(C,f_z(S_0)=f_z(S_1)=1)
=\mathbb E[\mathbf1_Ca_z^2]
```

and

```math
\mathbb P(C,f_z(S_0)=f_z(S_1)=f_z(T)=1)
=\mathbb E[\mathbf1_Ca_z\langle f_z,Kf_z\rangle],
```

(R49A.1) is also **exactly equivalent** to

```math
\boxed{
\mathbb P\{f_z(T)=1\mid
C, f_z(S_0)=f_z(S_1)=1\}\ge\lambda_1.}
\tag{R49A.2}
```

Equivalently, for

```math
b_z(R)=\Pr\{f_z(S)=1\mid S\supset R\},
```

the needed incidence inequality is

```math
\boxed{
\mathbb E_{z,R}[\mathbf1_Ca_zb_z(R)^2]
\ge\lambda_1\mathbb E_z[\mathbf1_Ca_z^2].}
\tag{R49A.3}
```

Thus the missing statement is not an ordinary two-point expansion estimate.
One independent ground incidence size-biases the center, and a second ground
incidence must retain another ground after common-core resampling.

### The near-core demand is almost deterministic exchange stability

Write `ell=m-s` and `k=n-m`.  The common-core step stays at `T=S_1` with
exact probability

```math
h_s=\binom{k+s}{s}^{-1},
\qquad
\lambda_1=\frac{(m-s)k}{m(k+s)}.
```

Let `r_s` be the conditional retention probability in (R49A.2), additionally
conditioned on `T != S_1`.  Since the self step is automatically favorable,
(R49A.2) is exactly

```math
\boxed{r_s\ge\frac{\lambda_1-h_s}{1-h_s}.}
\tag{R49A.4}
```

For one replacement this becomes

```math
\boxed{r_1\ge1-\frac1m-\frac1k.}
\tag{R49A.5}
```

At fixed density and `s=o(n)`, (R49A.4) requires retention
`1-O(s/n)`.  This makes precise why the route can only work after the
double-incidence law has created dictator-scale exchange stability.

For `s=1`, the port identity (10.1211) gives an exact signing interpretation.
For a pair `R+v,R+w`, `f_z(R+v)=1` precisely when some orientation makes all
three nonnegative port regrets `o_v`, `g_v`, and `r_v` vanish; simultaneous
groundhood means this happens at both ports, with possibly different
orientations.  Hence (R49A.5) asks that, after the independent ground
incidence and one zero-regret port are imposed, all but an
`1/m+1/k` fraction of alternate ports also have a zero-regret triple.  The
identity itself gives no control of the alternate optimized response gap.

## 2. The mass clause has a sharper direct target

Let `S` have size `m`, put `T=S^c`, `k=|T|`, and let `y` be a child ground:

```math
Q_S=Q(A[S])=|y^{\mathsf T}A[S]y|,
\qquad \Delta=q_n-Q_S.
```

Define the internal and external ground-lift Grams

```math
G_{\rm int}=\lVert A[S]y\rVert_2^2,
\qquad
G_{\rm ext}=\lVert A[T,S]y\rVert_2^2.
```

For a uniform outside word `w`, the exact completion formulas are

```math
\mathbb E_wR_2(y,w)=G_{\rm int}+G_{\rm ext}+k(n-1),
\tag{R49A.6}
```

```math
\mathbb E_w\{(y,w)^{\mathsf T}A(y,w)\}^2
=Q_S^2+4G_{\rm ext}+2k(k-1)\le q_n^2.
\tag{R49A.7}
```

Consequently the sharpest simple sufficient mass lemma is only

```math
\boxed{
\text{some child ground }(S,y)\text{ has }
\lVert A[:,S]y\rVert_2^2
=G_{\rm int}+G_{\rm ext}=O(R_*).}
\tag{R49A.8}
```

Indeed `k(n-1)=O(n^2)=o(R_*)` at the project scale, so (R49A.6) supplies an
outside completion with `R_2=O(R_*)`.  That one incidence already makes
`D_C=E[1_Ca_z^2]` positive.  Principal shortfall is a sufficient way to
control `G_ext`, but is not logically necessary for mass.

If one wants a target involving only internal Gram and shortfall, (R49A.7)
gives the exact joint Pareto certificate

```math
\boxed{
\mathcal P(S,y):=
G_{\rm int}+\frac{\Delta(2q_n-\Delta)}4,}
\tag{R49A.9}
```

and

```math
\boxed{
\mathbb E_wR_2(y,w)
\le\mathcal P(S,y)+k(n-1)-\frac{k(k-1)}2.}
\tag{R49A.10}
```

Thus `min_(S,y ground) P(S,y)=O(R_*)` proves mass.  Because
`q_n <= 2q_n-Delta <= 2q_n`, this joint condition is, up to constants,
equivalent to the two former requirements
`G_int=O(R_*)` and `Delta=O(R_*/q_n)`; its value is exact bookkeeping and a
single Pareto objective, not an asymptotic exponent improvement.  The direct
total-Gram target (R49A.8) is genuinely weaker.

The remaining minimizer-specific mass lemma can therefore be stated sharply:

> For some fixed `c>0`, every target exact minimizer and every active
> fixed-density restriction window admit a child-ground incidence satisfying
> `||A[:,S]y||_2^2=O(n^(9/4-c))`.

Known operator control gives only `O(n^(5/2))`.  A scalable family for which
the minimum total Gram is `omega(n^(9/4-c))` would falsify this sufficient
total-Gram route, though direct box cancellation could still supply mass.

### Fixed-density scale audit

The shortfall bridge is much more exceptional at fixed density than the
near-diagonal examples suggest.  If `k=n-m`, (R49A.7) forces

```math
\boxed{
\Delta(2q_n-\Delta)=q_n^2-Q_S^2
\ge2k(k-1),
\qquad
\Delta\ge\frac{k(k-1)}{q_n}.}
\tag{R49A.11}
```

At fixed density and `q_n=Theta(n^(3/2))`, every child ground therefore has
`Delta=Omega(n^(1/2))`.  The Pareto bridge asks
`Delta=O(R_*/q_n)=O(n^(3/4-c))`.  Since `c<1/4`, these exponents leave a
nonempty window: the bridge is **not logically restricted to the diagonal**.
But it asks for a fixed-density principal submatrix retaining
`q_n-o(q_n)`.  The natural restriction scale would instead have a
`Theta(n^(3/2))` gap, which contributes `Theta(n^3)` to `P` and destroys the
project bound.  This natural-scale statement is a heuristic, not a theorem:
`A[S]` need not be an exact order-`m` minimizer and no known result excludes
an exceptional near-full principal norm.

In particular, all original mass examples use `n-m<=2`; their constant
shortfalls provide no fixed-density evidence for (R49A.9).  The right status
is: nonvacuous but wholly unproved and unusually rigid at fixed density.

There is an exact replacement which bypasses both shortfall and the uniform
completion average.  Define the partial-completion box discrepancy

```math
\boxed{
\mathcal V(S,y)=
\min_{w\in\{\pm1\}^{T}}
\lVert A[:,S]y+A[:,T]w\rVert_2^2.}
\tag{R49A.12}
```

Then zero-deficit project-row mass is **exactly equivalent** to

```math
\boxed{
\min_{S,\ y\text{ a ground of }A[S]}\mathcal V(S,y)\le R_*.
}
\tag{R49A.13}
```

The total-Gram lemma (R49A.8) proves this by random completion, but
(R49A.13) can also hold through cancellation between the fixed child field
and the outside columns.  Standard law-free vector balancing does not prove
the project scale (the all-positive signing remains a scalable obstruction);
a proof must use exact-minimizer structure.  Nevertheless (R49A.13), rather
than near-full principal norm, is the sharp mass obligation.

## 3. Exact finite phase diagram

The checker exhausts every projective center, selector, child ground, row
threshold, and core size for each displayed matrix.  `A6,A8,A9` are the
stored exact minimizers.  Each order-ten matrix is also an exact minimizer
(`Q=26`) produced at the proved cap, but the twelve-matrix collection is only
a deterministic MILP sample.

At the old diagnostic cap `2n(n-1)`, all positive-degree centers in every
audited matrix are included.  The full ground-lift coarea verdict is:

| matrix and selector size | core sizes satisfying (R49A.1) | failing cores |
|:--|:--|:--|
| `A6, m=5` | `ell=1,2,3,4` | none |
| `A8, m=6` | `ell=1,2,3,4,5` | none |
| `A9, m=7` | `ell=1,2,3,4,5` | `ell=6` |
| twelve sampled `A10, m=8` | `ell=1,2,3,4` in all 12 | `ell=5,6,7` in all 12 |

Two exact failed triples, recorded in the order
`(coarea ratio, cubic ratio, triple-retention ratio)`, are

```text
A9, ell=6:       (10154/9891, 10541/9489, 12925/13188)
A10 seed 0,ell=5:(11369/11200,16025/14504,10693/11200).
```

The first two ratios must be at most one and the third at least one.  These
are exact-minimizer counterexamples to a theorem asserted uniformly at each
prescribed core scale, but not to an existence-of-scale theorem: `A9` succeeds at `ell<=5`,
and every sampled `A10` succeeds at `ell<=4` with positive `kappa`.

Row truncation is genuinely nonmonotone.  On `A8,ell=3`, the criterion holds
at row cap `40`, fails at `56`, and holds again at `64`.  On `A9,ell=3`, it
alternates success/failure at caps `56,72,80,88`.  Moreover, the minimum-row
class fails for `A8` at `ell=4,5`, for `A9` at `ell=3,4,5,6`, and for the
sampled `A10` at every `ell>=2`.  Thus neither lowering nor enlarging the row
cap has a monotone favorable effect on the cubic inequality.

## 4. Finite mass diagnostics

The following minima are exhaustive within each displayed matrix.  `P` is
(R49A.9), `G=G_int+G_ext`, and `r_min` is the best actual completion row.

| matrix | best `P` profile `(Delta,G_int,G_ext)` | best total `G` | `r_min` |
|:--|:--|--:|--:|
| `A6,m=5` | `33 : (2,24,1)` | `25` | `30` |
| `A8,m=6` | `73 : (2,54,0)` | `46` (at `Delta=6`) | `40` |
| `A9,m=7` | `107 : (2,84,2)` | `70` (at `Delta=6`) | `56` |
| 12 sampled `A10,m=8` | `105:(2,80,0)` in 7; `113:(2,88,0)` in 5 | `80` or `88` | `74` |

These order-`n^2` values are positive finite evidence for (R49A.8), and the
zero external Gram in every sampled order-ten Pareto optimizer suggests that
external-field cancellation may be more useful than principal shortfall.
They do not imply an asymptotic bound.  The `A10` statement is exhaustive
inside each of twelve distinct normalized matrices, not exhaustive over
order-ten exact minimizers.

Lower-ratio exhaustive diagnostics make the scale distinction concrete:

| matrix | `m/n` | best `P` and its `Delta` | best total Gram and its `Delta` | best exact row and its `Delta` |
|:--|:--:|:--|:--|:--|
| `A8` | `5/8` | `96, 8` | `27, 12` | `32, 8` |
| `A9` | `5/9` | `136, 8` | `28, 16` | `16, 16` |
| sampled `A10` | `6/10` | `134, 4` | `38, 16` | `10, 16` |

For `A9` and the sampled `A10`, the best actual completion has an excellent
row despite a shortfall at least sixty percent of `q`.  This is exact
finite evidence that the Pareto shortfall certificate can be very wasteful,
and that direct total-Gram or box cancellation is the appropriate replacement
target.  It remains finite evidence only.

## 5. Research judgment and reduced target

The complete surviving package is now the following exact incidence theorem.
For every target exact minimizer, find a fixed-density scale `ell` with
`kappa=(lambda_1-lambda_2)/(1-lambda_1)>=kappa_0>0` such that:

1. some child ground satisfies the exact box target (R49A.13) (with the
   total-Gram bound (R49A.8) as a convenient sufficient input), hence `D_C>0`;
2. under the project-row-truncated triple experiment, the retention bound
   (R49A.2) holds.

Together with (10.1222), this is exactly sufficient for a constant-degree
project center and hence the bare restriction tail.  No strict coarea margin
is needed.

The audit gives no proof of item 2.  It does sharply reject three tempting
strengthenings: a fixed nearest-core schedule, a monotone-in-row-cap argument,
and an inference from mass or small total Gram alone.  A viable proof must
instead establish an adaptive nonlocal scale or a signing-specific exchange
theorem for zero port-regret profiles after double incidence bias.  A genuine
asymptotic falsifier would be an unbounded exact-minimizer family for which
either the exact box target (R49A.13) fails at the project scale or every
admissible positive-mass scale with `kappa>=kappa_0` has triple retention
strictly below `lambda_1`.
