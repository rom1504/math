# Infinitesimal disorder geometry of the actual-child row shadow

Status: **task-local rigorous theorem note**.  This note computes the exact
small-disorder tangent of the row-product decomposition in (AC.24).  The
answer is the row Hoeffding/ANOVA decomposition of the *actual* parent
log-partition function.  It also identifies the first nonzero
bridge-amplitude coefficient of that ANOVA remainder with an explicit child
overlap tensor.  These are finite-system theorems.  No uniform passage from
`lambda -> 0` to fixed `lambda`, or from zero bridge amplitude to the physical
amplitude `t=beta/sqrt(N)`, is asserted.

The notation and normalization are those of
[`actual_child_negative_escort_structure.md`](actual_child_negative_escort_structure.md).
Fix one relative orientation and abbreviate its bridge pressure by `L(B)`.

## 1. Row ANOVA

Write the fair bridge law as

```math
U=U_1\otimes\cdots\otimes U_m,
\qquad R_i\in\{\pm1\}^n,
```

and put `mu=E_U L`.  The row Hoeffding decomposition is

```math
L=\mu+\sum_{\varnothing\ne S\subseteq[m]}L_S(R_S),   \tag{RA.1}
```

where every `L_S` is centered in each of its coordinates.  The summands are
orthogonal in `L^2(U)`.  Define

```math
L_{\rm add}=\sum_iL_{\{i\}},
\qquad L_{\rm cross}=\sum_{|S|\ge2}L_S,
```

and

```math
\sigma_{\rm add}^2=\|L_{\rm add}\|_2^2
=\sum_i\|L_{\{i\}}\|_2^2,
\qquad
\sigma_{\rm cross}^2=\|L_{\rm cross}\|_2^2.          \tag{RA.2}
```

Thus `Var_U(L)=sigma_add^2+sigma_cross^2`.

For `lambda>0`, recall

```math
V_\lambda=-{1\over\lambda}\log\mathbb E_Ue^{-\lambda L},
```

the row-product variational value `V_lambda^row` from (AC.15), its gain

```math
A_\lambda:=\mathbb E_UL-V_\lambda^{\rm row},          \tag{RA.3}
```

and the reverse information projection

```math
\mathcal I_\lambda^{\leftarrow}
=\inf_{p=\otimes_i p_i}D(p\Vert q_\lambda)
=\lambda(V_\lambda^{\rm row}-V_\lambda).             \tag{RA.4}
```

## 2. Exact infinitesimal split

**Theorem RA.1 (row-ANOVA tangent of the negative escort).**  For every
finite row-product space and every real landscape `L`, as `lambda downarrow 0`,

```math
\boxed{
\begin{aligned}
G_\lambda:=\mathbb E_UL-V_\lambda
 &= {\lambda\over2}
    (\sigma_{\rm add}^2+\sigma_{\rm cross}^2)
    -{\lambda^2\over6}\kappa_3(L)+O_L(\lambda^3),\\
A_\lambda
 &= {\lambda\over2}\sigma_{\rm add}^2
    -\lambda^2 C_{\rm row}(L)+O_L(\lambda^3),\\
\mathcal I_\lambda^{\leftarrow}
 &= {\lambda^2\over2}\sigma_{\rm cross}^2
    +\lambda^3\left\{C_{\rm row}(L)
                    -{\kappa_3(L)\over6}\right\}
    +O_L(\lambda^4).
\end{aligned}}                                      \tag{RA.5}
```

Here `kappa_3(L)=E_U[(L-mu)^3]` and the explicit second-order row term is

```math
C_{\rm row}(L)
={1\over6}\sum_i\mathbb E_U L_{\{i\}}^3
 +\sum_{i<j}\langle L_{\{i,j\}},
                         L_{\{i\}}L_{\{j\}}\rangle_{L^2(U)}.   \tag{RA.6}
```

Moreover, any global minimizing row product is unique for all sufficiently
small `lambda>0`, lies on the analytic critical-point branch through
`lambda=0`, and its row scores obey

```math
{dp_{\lambda,i}^*\over dU_i}
=1-\lambda L_{\{i\}}+O_L(\lambda^2).                 \tag{RA.7}
```

The ordinary row total correlation has the same quadratic tangent as the
oppositely directed projection:

```math
\boxed{
\operatorname{TC}(q_\lambda)
={\lambda^2\over2}\sigma_{\rm cross}^2+O_L(\lambda^3).}        \tag{RA.8}
```

In particular,

```math
\lim_{\lambda\downarrow0}{2A_\lambda\over\lambda}
=\sigma_{\rm add}^2,
\qquad
\lim_{\lambda\downarrow0}{2\mathcal I_\lambda^{\leftarrow}
                              \over\lambda^2}
=\lim_{\lambda\downarrow0}{2\operatorname{TC}(q_\lambda)
                              \over\lambda^2}
=\sigma_{\rm cross}^2.                               \tag{RA.9}
```

*Proof.*  The ordinary cumulant expansion gives

```math
V_\lambda
=\mu-{\lambda\over2}\operatorname{Var}_U(L)
     +{\lambda^2\over6}\kappa_3(L)+O_L(\lambda^3).   \tag{RA.10}
```

For completeness, multiply the row-product objective by `lambda` and write

```math
{dp_i\over dU_i}=1+\lambda h_i+\lambda^2k_i+O(\lambda^3),
\qquad \mathbb E h_i=\mathbb E k_i=0.                \tag{RA.11}
```

The coefficient of `lambda^2` in
`D(p||U)+lambda E_pL` is

```math
\sum_i\left({1\over2}\|h_i\|_2^2
             +\langle L_{\{i\}},h_i\rangle\right), \tag{RA.12}
```

so `h_i=-L_{\{i\}}`.  At this minimizer all terms involving `k_i`
cancel in the coefficient of `lambda^3`; that coefficient is exactly
(RA.6).  Hence

```math
V_\lambda^{\rm row}
=\mu-{\lambda\over2}\sigma_{\rm add}^2
     +\lambda^2C_{\rm row}(L)+O_L(\lambda^3).         \tag{RA.13}
```

Equations (RA.5) follow from (RA.4), (RA.10), and (RA.13).

To justify expansion at a global rather than merely local minimizer, note
that comparison with `U` gives `D(p_\lambda^*||U)=O_L(lambda)` for
`lambda>0`.
Therefore every global minimizer tends to `U`.  Near `U`, the Hessian of the
sum of marginal entropies is positive definite on the mean-zero tangent
space.  The analytic implicit-function theorem gives one analytic critical
point there, which contains all global minimizers for small `lambda>0`.
This also proves (RA.7).  Finally, the score of `q_lambda` at zero is
`-(L-mu)`, while the score of the product of its row marginals is
`-L_add`.  The standard quadratic expansion of relative entropy gives
(RA.8). `square`

**Interpretation.**  At infinitesimal disorder temperature the two terms in
the exact finite-`lambda` dichotomy (AC.24) are not mysterious.  Additive row
ANOVA is precisely the row-product gain, while all ANOVA orders involving at
least two rows are precisely the irreducible directed dependence, to leading
order.  The two KL directions agree only at this quadratic tangent; this does
not compare them at fixed `lambda`.

## 3. The actual-child overlap tensor

The preceding theorem is exact for the physical bridge pressure `L`.  We now
identify its cross-row tangent when the bridge interaction itself is turned
on from zero while the two actual children remain at their contracted raw
temperature `t=beta/sqrt(N)`.

Fix actual child minimizers `A,D` and an orientation `epsilon`.  Let
`mu_epsilon` denote the normalized augmented child-spin law in that sector at
zero bridge field.  It is the law of `(X,Y,tau)` proportional to

```math
\exp\{t\tau(H_A(X)+\epsilon H_D(Y))\}.                \tag{RA.14}
```

Set `Q_ij=tau X_iY_j`.  Global spin flip of either child makes the law of
`Q` centrally symmetric.  For a separate bridge amplitude `u`, define,
up to an irrelevant additive constant,

```math
L_u(B)=\log\mathbb E_{\mu_\epsilon}
             \exp\{u\langle Q,B\rangle\}.             \tag{RA.15}
```

The physical parent uses `u=t`.  Introduce the sector overlap tensor

```math
\Gamma_{ik;j\ell}^{(\epsilon)}
=\mathbb E_{\mu_\epsilon}[X_iX_kY_jY_\ell].           \tag{RA.16}
```

**Theorem RA.2 (overlap formula for infinitesimal cross-row mass).**  At
fixed finite `m,n,t,A,D,epsilon`,

```math
\boxed{
\lim_{u\to0}{\sigma_{\rm cross}^2(L_u)\over u^4}
=K_{\rm cross}(A,D,t,\epsilon)
:=\sum_{i<k}\sum_{j,\ell}
       \left(\Gamma_{ik;j\ell}^{(\epsilon)}\right)^2.}          \tag{RA.17}
```

The corresponding additive-row coefficient is

```math
\lim_{u\to0}{\sigma_{\rm add}^2(L_u)\over u^4}
=m\sum_{j<\ell}
 \left(\mathbb E_{\mu_\epsilon}Y_jY_\ell\right)^2.   \tag{RA.18}
```

In particular,

```math
K_{\rm cross}
\ge n\sum_{i<k}
 \left(\mathbb E_{\mu_\epsilon}X_iX_k\right)^2.       \tag{RA.19}
```

Writing `E_m=m(m-1)/2` and `E_n=n(n-1)/2`, a second useful scalar
consequence is

```math
K_{\rm cross}
\ge {\left(\mathbb E_{\mu_\epsilon}
                    H_A(X)H_D(Y)\right)^2\over E_mE_n}.          \tag{RA.19a}
```

(The ordered `j,l` sum in fact gives an extra factor two if only the
off-diagonal `j<l` terms are used; (RA.19a) records the safer weaker form.)

The tensor in (RA.16) is computable entirely from child correlations.  If
`s=tau`, then conditional on `s` the children are independent.  Thus, with
the exact sector weights `pi_s`,

```math
\Gamma_{ik;j\ell}^{(\epsilon)}
=\sum_{s=\pm1}\pi_s
 C_{A,s}(i,k)C_{D,\epsilon s}(j,\ell),                \tag{RA.20}
```

where `C_(A,s)(i,k)=E[X_iX_k|tau=s]`, and similarly for `D`.

*Proof.*  Central symmetry gives

```math
L_u(B)={u^2\over2}\mathbb E\langle Q,B\rangle^2+O(u^4).
                                                                  \tag{RA.21}
```

After deleting the constant diagonal terms, the quadratic Walsh polynomial
is

```math
u^2\sum_{e<f}\mathbb E[Q_eQ_f]B_eB_f.                \tag{RA.22}
```

Walsh monomials are orthonormal under `U`.  A monomial belongs to the
cross-row ANOVA space exactly when its two edges lie in different bridge
rows.  Taking those rows as `i<k` and their columns as arbitrary `j,l`
gives (RA.17).  The same-row monomials give (RA.18).  Restricting (RA.17)
to `j=l` gives (RA.19).  Applying Cauchy--Schwarz to the coefficients
`a_(ik)d_(jl)` over `i<k,j<l` gives (RA.19a).  Conditional independence
given `tau` proves (RA.20). `square`

This is a genuine statement about the actual optimized children, not a
conference surrogate.  It reduces the infinitesimal cross-row question to a
specific four-spin overlap tensor rather than to the full bridge landscape.

## 4. An exact physical-amplitude cross-response identity

There is also a characterization at the physical bridge amplitude, although
it uses parent Gibbs covariances rather than child overlaps alone.  For fair
independent replacement rows `R_i',R_k'`, let `B^(i)`, `B^(k)`, and
`B^(ik)` denote replacement of the indicated rows and put

```math
\square_{ik}L(B)
=L(B)-L(B^{(i)})-L(B^{(k)})+L(B^{(ik)}).              \tag{RA.21a}
```

**Proposition RA.3 (mixed row response counts cross ANOVA).**  One has

```math
{1\over4}\mathbb E(\square_{ik}L)^2
=\sum_{S\supseteq\{i,k\}}\|L_S\|_2^2,               \tag{RA.22a}
```

and hence, if

```math
\mathfrak J_2(L)
:={1\over4}\sum_{i<k}\mathbb E(\square_{ik}L)^2,
```

then

```math
\boxed{
\sigma_{\rm cross}^2
\le\mathfrak J_2(L)
\le {m\choose2}\sigma_{\rm cross}^2.}               \tag{RA.23a}
```

For the actual parent log-partition, each mixed response is exactly an
integrated Gibbs covariance.  Interpolate the two replacement rows by
`R_i'+s(R_i-R_i')` and `R_k'+v(R_k-R_k')`, and let `nu_(s,v)` be the parent
augmented Gibbs law at that real bridge.  With

```math
Z_i=\tau X_i\sum_j(R_{ij}-R'_{ij})Y_j,
\qquad
Z_k=\tau X_k\sum_j(R_{kj}-R'_{kj})Y_j,
```

one has

```math
\boxed{
\square_{ik}L
=t^2\int_0^1\!\int_0^1
 \operatorname{Cov}_{\nu_{s,v}}(Z_i,Z_k)\,ds\,dv.}   \tag{RA.24a}
```

*Proof.*  A row difference multiplies the squared norm of every ANOVA
component containing that row by two.  Applying two independent row
differences proves (RA.22a).  Summing over pairs weights `L_S` by
`binom(|S|,2)`, which lies between one and `binom(m,2)` for `|S|>=2`; this
proves (RA.23a).  Finally, first and mixed derivatives of a log-partition are
respectively a Gibbs mean and a Gibbs covariance.  The two-dimensional
fundamental theorem of calculus gives (RA.24a). `square`

This converts the fixed-amplitude question into a concrete response
criterion.  It is not yet a closure: without an upper bound on the ANOVA
order, `J_2` may overcount a high-order component by `Theta(m^2)`, and the
covariance in (RA.24a) is a parent quantity along an interpolation, not a
zero-bridge child statistic.

## 5. What this does and does not buy at the physical scale

Theorems RA.1 and RA.2 expose two independent non-uniformity barriers.

1. **Disorder-temperature barrier.**  The remainders in (RA.5) are
   finite-system Taylor remainders.  The elementary bridge flip bound controls
   `Var(L)=O(N)` but supplies no uniform third- or higher-cumulant estimate
   strong enough to continue (RA.5) to a fixed positive `lambda`.  In
   particular, an exponentially small lower tail can be invisible to the
   variance and dominate a fixed negative moment.

2. **Bridge-amplitude barrier.**  Formula (RA.17) is the derivative at
   `u=0`; the physical point is `u=t=beta/sqrt(N)`.  Although `u` tends to
   zero, the number of bridge coordinates is `Theta(N^2)`, so the `O(u^4)`
   remainder in (RA.21) is not uniform in `N`.  Higher connected child
   overlaps can contribute at leading order.  Replacing them by the pair
   tensor without a cluster/cumulant theorem would be an invalid
   high-temperature extrapolation.

There is also an optimizer-information mismatch.  The singleton child
minimality inequality (AC.33) controls signed, `tau`-odd correlations

```math
a_{ik}\,\mathbb E[\tau X_iX_k]\le\tanh t,             \tag{RA.23}
```

whereas (RA.16)--(RA.20) involve `tau`-even sector overlaps and products of
conditional correlations.  Neither (RA.23) nor its sum gives a positive
lower bound on `K_cross`.  The full edge-set family (AC.32) contains more
information, but using it to determine every tensor in (RA.20) would simply
move the missing structural problem into the complete flip-response family.

Consequently the rigorous infinitesimal conclusion is a **classification,
not a fixed-`lambda` closure**:

```math
\boxed{
\text{infinitesimal irreducible row dependence}
\quad\longleftrightarrow\quad
\text{cross-row ANOVA mass of the actual }L,
}                                                       \tag{RA.24}
```

and the first bridge-amplitude coefficient of that mass is the actual-child
overlap norm (RA.17).  A scalable use of this classification requires one of
the following genuinely new estimates:

- a uniform cumulant/cluster bound transporting (RA.17) to `u=t` and (RA.5)
  to the fixed target `lambda`; or
- a direct `Theta(N)` lower or `o(N)` upper bound on
  `sigma_cross^2(L_t)` from actual-child structure.

Without one of these, infinitesimal information geometry does not reset the
fixed-`lambda` smallest missing lemma in (AC.24).
