# Global certificates for the actual-child row-product shadow

Status: **rigorous theorem note plus reproducible finite certificate**.  The
only landscapes used in the finite application are the bridge pressures of
the contracted-temperature minimizing children.  No conference or Paley
surrogate occurs.

The note has two purposes.  First, it gives a checkable strong-convexity
criterion under which the nonconvex row-product variational problem is in
fact globally solved by its coordinate Gibbs fixed point.  Second, when that
criterion fails at the target-reaching order-eight law, it obtains a global
lower certificate by pushing the actual escort through four deterministic
row features.  This turns the previously one-sided coordinate-descent value
into a genuine interval for the directed product projection.

## 1. Rectangle curvature and a global fixed-point certificate

Let `X_i` be finite sets, let `u_i` have full support on `X_i`, and put

```math
 {cal F}_\lambda(p_1,\ldots,p_m)
 =\mathbb E_{\otimes_i p_i}f
   +{1\over\lambda}\sum_iD(p_i\Vert u_i).                 \tag{GC.1}
```

For `i\ne j`, define the row-rectangle oscillation

```math
 C_{ij}=\sup\left|
 f(a,b,z)-f(a',b,z)-f(a,b',z)+f(a',b',z)
 \right|,                                                 \tag{GC.2}
```

where the supremum is over `a,a' in X_i`, `b,b' in X_j`, and
`z in product_(k notin {i,j})X_k`.  Set `C_ii=0`, so `C` is a symmetric
nonnegative matrix.

**Theorem GC.1 (rectangle-Hessian global certificate).**  Suppose
`r=tensor_i r_i` is a coordinate Gibbs fixed point:

```math
 {dr_i\over du_i}(a)
 \propto
 \exp\{-\lambda\mathbb E_{r_{-i}}f(a,X_{-i})\}.           \tag{GC.3}
```

If

```math
 \boxed{\lambda\lambda_{\max}(C)<4,}                    \tag{GC.4}
```

then `r` is the unique global minimizer of (GC.1).  More precisely, for
every product law `p`, with `d_i=TV(p_i,r_i)`,

```math
 {cal F}_\lambda(p)-{cal F}_\lambda(r)
 \ge {1\over2}d^T\left({4\over\lambda}I-C\right)d.       \tag{GC.5}
```

There is also an a posteriori version.  For an arbitrary full-support
product law `r`, define

```math
 \eta_i=\operatorname{osc}_{a\in X_i}\left[
  \mathbb E_{r_{-i}}f(a,X_{-i})
  +{1\over\lambda}\log{r_i(a)\over u_i(a)}
 \right].                                                \tag{GC.6}
```

If `H=(4/lambda)I-C` is positive definite, then

```math
 \boxed{
 {cal F}_\lambda(r)-{1\over2}\eta^TH^{-1}\eta
 \le \min_{p\ {\rm product}}{cal F}_\lambda(p)
 \le {cal F}_\lambda(r).}                              \tag{GC.7}
```

*Proof.*  Write `delta_i=p_i-r_i` and interpolate
`r_i(s)=r_i+s delta_i`.  If `g(s)=E_(tensor_i r_i(s)) f`, then

```math
 g''(s)=2\sum_{i<j}
 \mathbb E_{\delta_i\otimes\delta_j\otimes
             \bigotimes_{k\ne i,j}r_k(s)}f.             \tag{GC.8}
```

The Jordan decomposition
`delta_i=d_i(alpha_i-beta_i)` and the definition (GC.2) give

```math
 \left|\mathbb E_{\delta_i\otimes\delta_j\otimes\nu}f\right|
 \le C_{ij}d_id_j.                                      \tag{GC.9}
```

Taylor's formula therefore bounds the nonlinear energy remainder below by
`-sum_(i<j)C_ij d_i d_j`.  Expanding each entropy around `r_i` and using
(GC.3) cancels the linear energy term exactly.  Pinsker's sharp local-scale
bound `D(p_i||r_i)>=2d_i^2` proves (GC.5).  Under (GC.4) its right side is
strictly positive away from `r`.

Without (GC.3), the uncancelled linear term is at least
`-sum_i eta_i d_i`.  Hence

```math
 {cal F}_\lambda(p)-{cal F}_\lambda(r)
 \ge {1\over2}d^THd-\eta^Td
 \ge-{1\over2}\eta^TH^{-1}\eta,                        \tag{GC.10}
```

which proves (GC.7). `square`

The constant four is not a convention: it comes from the factor two in
Pinsker and the factor one-half in the symmetric quadratic form.  The
criterion is a genuine global certificate, not merely local Hessian
positivity.

## 2. Coarse row features give rigorous lower certificates

Let `q` be any full-support law on `product_i X_i`, and let
`phi_i:X_i -> Y_i` be deterministic row maps.  Write

```math
 Q=(\phi_1,\ldots,\phi_m)_\#q.                          \tag{GC.11}
```

**Theorem GC.2 (row-feature data-processing certificate).**

```math
 \boxed{
 \inf_{p=\otimes_i p_i}D(p\Vert q)
 \ge
 \inf_{P=\otimes_i P_i}D(P\Vert Q).}                   \tag{GC.12}
```

*Proof.*  The image of a row-product law under a rowwise deterministic map
is again row-product.  Apply KL data processing to each `p`, then minimize.
`square`

This lower bound is useful in exactly the direction required by the actual
escort identity

```math
 \mathcal I_\lambda^{\leftarrow}
 =\inf_{p\ {\rm row-product}}D(p\Vert q_\lambda)
 =\lambda(V_\lambda^{\rm row}-V_\lambda).              \tag{GC.13}
```

It is not forward total correlation.  Nor is it an unconstrained numerical
lower bound: after the feature maps are fixed, (GC.12) is an exact global
certificate.

For binary `Y_i={-1,1}`, put `s_i=E_P Y_i`.  If

```math
 \log Q(y)=\sum_{S\subseteq[m]}a_S\prod_{i\in S}y_i,    \tag{GC.14}
```

then its reverse product projection is the global minimum on `[-1,1]^m` of

```math
 -\sum_Sa_S\prod_{i\in S}s_i
 -\sum_i h\left({1+s_i\over2}\right),                  \tag{GC.15}
```

where `h` is binary entropy.  Thus GC.1 or elementary interval subdivision
can certify the coarse minimum without making any claim about coordinate
descent on the original row alphabets.

## 3. The target-reaching `N=8` actual child law

Take the balanced split `4+4`, `beta=4`, relative child orientation `-1`,
and

```math
 \lambda_*=5.382104195764755.                           \tag{GC.16}
```

Both children are the actual contracted-temperature minimizers selected by
complete signing/histogram enumeration.  On each four-bit bridge row `R_i`,
use the Walsh parities with masks

```text
 (2, 1, 4, 8).                                         (GC.17)
```

Exact integer pressure-signature counts show that the resulting sixteen
coarse atoms fall into four equality classes

```text
 {0,15}, {3,12}, {1,2,4,7,8,11,13,14}, {5,6,9,10}.     (GC.18)
```

At 80-digit outward interval precision their probabilities are respectively

```text
 .1951584116330018713321975922505871070...,
 .1733019534789941070628371479563939390...,
 .03279887792133212479506343585442076951...,
 .000172061601337761212355758187667937933....           (GC.19)
```

The equality classes in (GC.18) are checked at the integer-signature level,
before evaluating any exponential or logarithm.  Consequently (GC.14) has
exactly the form

```math
 \log Q(y)=a_0
 +J(y_1y_2+y_3y_4)
 +k(y_1+y_2)(y_3+y_4)
 -K y_1y_2y_3y_4,                                      \tag{GC.20}
```

with outward intervals centered at

```text
 -a_0 = 4.298927923054102131196931699355473...,
 J    = 1.743581562871090006699719064467845...,
 k    = 0.0148470159029785475578036352385368...,
 K    = 0.881566949161413464739541011709969....          (GC.21)
```

There is a useful exact dimension reduction.  Holding `s_3,s_4` fixed, the
coefficient of `-s_1s_2` in (GC.15) is
`J-Ks_3s_4>=J-K>0`, and the linear coefficients of `s_1,s_2` agree.  Every
minimizer is interior.  Subtracting its two first-order equations gives

```math
 \operatorname{atanh}s_1-\operatorname{atanh}s_2
 +(J-Ks_3s_4)(s_1-s_2)=0,                              \tag{GC.22}
```

so `s_1=s_2`.  Similarly `s_3=s_4`.  The `k>0` term and global sign symmetry
then reduce the global problem to `u,v in [0,1]`:

```math
 -a_0-J(u^2+v^2)-4kuv+Ku^2v^2-2h((1+u)/2)-2h((1+v)/2). \tag{GC.23}
```

Outward interval branch-and-bound on (GC.23), using only interval products
and exact lower bounds on `-h`, visits `752768` boxes and proves

```math
 \boxed{
 1.075
 \le \inf_{P\ {\rm product}}D(P\Vert Q)
 \le 1.075620.}                                        \tag{GC.24}
```

The upper endpoint is an evaluated feasible coarse product.  GC.2 therefore
gives the first nonzero global lower certificate for the original actual
escort's reverse row-product projection:

```math
 \boxed{1.075\le\mathcal I_{\lambda_*}^{\leftarrow}
                 \le4.506450.}                         \tag{GC.25}
```

The upper endpoint is the previously evaluated feasible full row-product
law.  Equivalently,

```math
 V_{\lambda_*}+{1.075\over\lambda_*}
 \le V_{\lambda_*}^{\rm row}
 \le12.751691.                                         \tag{GC.26}
```

An independent outward-interval recomputation encloses both the soft value
and the same-temperature child target at the displayed decimal
`lambda_*`; it gives directly

```math
 V_{\lambda_*}^{\rm row}-T\ge0.1997360067547\ldots .
```

Thus the numerically located threshold is not being treated as an exact
root.  In particular (GC.26) certifies a strictly positive target excess of
at least `0.19973`, or about `.02497 N`.  This is a finite theorem, not an
asymptotic claim.

## 4. Balanced `N=9` actual child laws at `beta=2,4`

For `lambda=1`, choose row-Walsh masks `(1,1,2,2)` at `beta=2` and
`(1,1,4,4)` at `beta=4`.  The resulting coarse log densities have only
degree two and four even Walsh terms.  Their rectangle matrices obey the
rigorous row-sum bounds

```text
 beta=2: max row sum(C) < .571,
 beta=4: max row sum(C) < 1.266.                       (GC.27)
```

Both are strictly below four.  The uniform coarse product is a coordinate
fixed point by complement symmetry, so GC.1 certifies it as the unique
global product minimizer.  Interval evaluation gives

```math
 \begin{array}{c|c}
 \beta&\inf_{P\ {\rm product}}D(P\Vert Q)\\ \hline
 2&0.0066376686168\ldots\\
 4&0.0273458654445\ldots
 \end{array}                                           \tag{GC.28}
```

and these are rigorous lower bounds on the full actual escort's directed
row dependence.

## 5. Sharp falsifier for direct strong convexity

The global theorem is nonvacuous on the same actual-child family at smaller
`beta` and `lambda`.  It does **not**, however, certify the full row problem
in the target regime.  Exact complete-cube rectangle computations give

```text
 case                              lambda_max(C)   lambda lambda_max(C)
 N=8, beta=4, lambda=5.382104       35.6652          191.95
 N=9, beta=2, lambda=1              15.0174           15.02
 N=9, beta=4, lambda=1              41.6182           41.62
```

against the required threshold four.  The failure factors are too large to
be repaired by rounding constants.  Thus the full-alphabet Dobrushin route
is decisively falsified at these actual children, while a deliberately
coarser exposed dependence witness still proves a positive global gap.

This distinction is the structural conclusion: a fixed point can be
globally certified in a compressed response channel even when the complete
row response map is far outside the contraction regime.

## 6. Scope

- GC.1 is a scalable sufficient condition, but computing every rectangle in
  (GC.2) can itself query the full landscape.  It is not by itself the
  low-information asymptotic lemma sought by the main campaign.
- GC.24--GC.28 concern exact finite actual-child laws.  They prove genuine
  directed dependence; they do not prove that this dependence stays linear
  in `N`.
- The Walsh maps were selected to expose dependence and are then frozen.
  Data processing makes post-selection harmless for the stated certificate.
- The result rules out the hypothesis that the order-eight target-reaching
  product excess is merely a coordinate-descent artifact.  It does not
  certify the full candidate as globally optimal.

The reproducible implementation and machine-readable intervals are in
[`../experiments/actual_child_row_product_certificate.py`](../experiments/actual_child_row_product_certificate.py)
and
[`../../computations/results/actual_child_row_product_certificate.json`](../../computations/results/actual_child_row_product_certificate.json).
