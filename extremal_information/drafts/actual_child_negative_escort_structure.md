# The actual child-induced negative escort: row regularity and directed dependence

Status: **task-local rigorous theorem note**.  All statements below concern the
bridge law induced by the actual contracted-temperature children.  The main
positive result is an exact product-shadow reduction: the optimal independent-
row approximation automatically has uniformly bounded row Renyi-two
complexity, while the value lost by imposing row independence is exactly a
directed relative-entropy distance.  This narrows the remaining question to a
single extensive-dependence alternative.  It does **not** prove that this
alternative occurs for minimizing children, nor that it is absent.

The normalization is that of
`artifacts/finite_temperature_reverse_kl_interface.md` and
`artifacts/two_temperature_bridge_audit.md`.

## 1. Exact setup and the orientation peel

Let `m+n=N`, let

```math
t={\beta\over\sqrt N},
```

and let `A,D` minimize the child pressures at this raw temperature (equivalently
at scaled temperatures `beta sqrt(m/N)` and `beta sqrt(n/N)`).  For a relative
orientation `epsilon in {+-1}` and bridge
`B in {+-1}^{m times n}`, put

```math
L_\epsilon(B)
=\log \overline Z_N(A,\epsilon D,B;t).                 \tag{AC.1}
```

On the joint uniform law `U_epsilon tensor U_B`, define the negative-disorder
escort

```math
{dq_\lambda\over d(U_\epsilon\otimes U_B)}(\epsilon,B)
={e^{-\lambda L_\epsilon(B)}
 \over \mathbb E_{U_\epsilon U_B}e^{-\lambda L}},
\qquad \lambda>0.                                      \tag{AC.2}
```

Conditionally on `epsilon`, this is the sector escort

```math
{dq_{\lambda,\epsilon}\over dU_B}(B)
={e^{-\lambda L_\epsilon(B)}
 \over \mathbb E_{U_B}e^{-\lambda L_\epsilon}}.        \tag{AC.3}
```

The orientation is one binary latent variable.  The KL chain rule gives

```math
D(q_\lambda\Vert U_\epsilon U_B)
=D(q_\lambda^\epsilon\Vert U_\epsilon)
 +\mathbb E_{q_\lambda^\epsilon}
   D(q_{\lambda,\epsilon}\Vert U_B),                   \tag{AC.4}
```

and the first term is at most `log 2`.  Thus an extensive bridge resource
cannot be hidden in the orientation.  Because a deterministic parent may
choose the better orientation for free, the structural results below are
stated sectorwise.  They hold uniformly in `epsilon,A,D`; in particular they
hold for the actual minimizing pair.

The exact forward output law `Pi` from the reverse-KL interface satisfies

```math
{d\Pi_\epsilon\over dU_B}(B)
=c_\epsilon e^{L_\epsilon(B)}.                          \tag{AC.5}
```

Consequently the actual negative escort is the inverse power escort

```math
{dq_{\lambda,\epsilon}\over dU_B}
= {p_\epsilon^{-\lambda}
    \over \mathbb E_{U_B}p_\epsilon^{-\lambda}},
\qquad p_\epsilon={d\Pi_\epsilon\over dU_B}.            \tag{AC.6}
```

This identity is useful conceptually: the child Gibbs variables give a
product-channel latent representation of `Pi`, not of its inverse escort.

## 2. A cube lemma

We use the following elementary bounded-oscillation statement.

**Lemma AC.1 (conditional escort regularity).**  Let `F` be a function on a
Boolean cube and suppose flipping coordinate `j` changes `F` by at most
`c_j`.  Let `q proportional exp(-lambda F) U`.  Fix any set of coordinates,
condition on arbitrary values of a disjoint set, and marginalize all remaining
coordinates.  If `S` is the retained set, then

```math
D_2(q_S\Vert U_S)
\le\sum_{j\in S}\log\!\left(
 1+\tanh^2(\lambda c_j/2)\right)
\le {\lambda^2\over4}\sum_{j\in S}c_j^2.                \tag{AC.7}
```

Moreover, in any sequential reveal order, the conditional KL increment of
coordinate `j` is at most

```math
\kappa(\lambda c_j/2),
\qquad
\kappa(a)=a\tanh a-\log\cosh a\le {a^2\over2}.          \tag{AC.8}
```

*Proof.*  After fixing and marginalizing coordinates, the unnormalized log
density `g` on `S` still has coordinate oscillation at most `lambda c_j`:
the pointwise integrands before and after a flip differ by a factor in
`[exp(-lambda c_j),exp(lambda c_j)]`.  Reveal the retained bits in any
order and write `a_j` for half their predictable conditional log odds.  Then
`|a_j|<=lambda c_j/2`, and the likelihood ratio with respect to fair bits is

```math
\ell=\prod_j{e^{a_jb_j}\over\cosh a_j}.
```

Now

```math
\ell^2
=\left[\prod_j{e^{2a_jb_j}\over\cosh(2a_j)}\right]
 \prod_j{\cosh(2a_j)\over\cosh^2a_j}.
```

The bracketed predictable product is a mean-one martingale under the fair
reveal law, while the last factor is at most
`prod_j(1+tanh^2(lambda c_j/2))`.  Taking expectations and logarithms proves
the first inequality in (AC.7); `log(1+tanh^2 u)<=u^2` proves the second.

For one revealed bit, its conditional log odds have magnitude at most
`lambda c_j`.  A Bernoulli law with log odds `2a` has KL from the fair bit
equal to `kappa(a)`, which is increasing in `|a|`.  Finally
`(a^2/2-kappa(a))'=a tanh^2(a)>=0`.  This proves (AC.8). `square`

## 3. Uniform conditional Renyi complexity of the actual bridge law

Flipping one bridge sign changes every parent Hamiltonian value by exactly
`2` in absolute value.  Comparing the augmented exponential sums term by
term therefore gives

```math
|L_\epsilon(B)-L_\epsilon(B^{(ij)})|\le2t.              \tag{AC.9}
```

Write the bridge as `m` rows `R_1,...,R_m`, each of length `n`, and set

```math
a=\lambda t={\lambda\beta\over\sqrt N}.
```

**Theorem AC.2 (actual-law row-filtration compactness).**  For every row
order, every prefix value, and every `i`,

```math
\boxed{
D_2\!\left(q_{\lambda,\epsilon}
 (R_i\mid R_{<i})\,\middle\Vert\,U_n\right)
\le n\log(1+\tanh^2(\lambda t))
\le\lambda^2t^2 n
=\lambda^2\beta^2{n\over N}.}                         \tag{AC.10}
```

The same bound holds if any collection of other bridge coordinates is fixed
and all unmentioned coordinates are marginalized.  In particular, for
comparable splits and fixed `beta,lambda`, every conditional row component
has a common bounded Renyi-two constant.  There is no escaping conditional
`D_2` mass in the natural row filtration of the actual negative escort.

For the whole sector,

```math
\begin{aligned}
D(q_{\lambda,\epsilon}\Vert U_B)
 &=\phi_\epsilon(\lambda)
   -\lambda\phi_\epsilon'(\lambda),
\qquad
\phi_\epsilon(\lambda)
 =-\log\mathbb E_{U_B}e^{-\lambda L_\epsilon},          \tag{AC.11a}\\
D(q_{\lambda,\epsilon}\Vert U_B)
 &\le mn\,\kappa(\lambda t)
 \le {\lambda^2t^2mn\over2},                            \tag{AC.11}\\
D_2(q_{\lambda,\epsilon}\Vert U_B)
 &\le mn\log(1+\tanh^2(\lambda t))
 \le\lambda^2t^2mn.                                    \tag{AC.12}
\end{aligned}
```

Consequently its row total correlation obeys

```math
\boxed{
\operatorname{TC}_\epsilon(R_1;\ldots;R_m)
=\sum_{i=1}^m I_{q_{\lambda,\epsilon}}(R_i;R_{<i})
=D(q_{\lambda,\epsilon}\Vert U_B)
 -\sum_{i=1}^mD((q_{\lambda,\epsilon})_{R_i}\Vert U_n).}
                                                               \tag{AC.13b}
```

and hence

```math
\boxed{
\operatorname{TC}_\epsilon(R_1;\ldots;R_m)
\le D(q_{\lambda,\epsilon}\Vert U_B)
\le {\lambda^2\beta^2mn\over2N}=O_{\beta,\lambda}(N).} \tag{AC.13}
```

In density form, if `m/N -> theta in (0,1)`, then

```math
\limsup {D(q_{\lambda,\epsilon}\Vert U_B)\over N},
\quad
\limsup {\operatorname{TC}_\epsilon\over N}
\le {\lambda^2\beta^2\over2}\theta(1-\theta).           \tag{AC.13a}
```

Also `L_epsilon(B)=L_epsilon(-B)` (send `x` to `-x` in the partition
sum), so `q_(lambda,epsilon)(B)=q_(lambda,epsilon)(-B)`.  Every individual
bridge bit is exactly unbiased.  Any positive row marginal KL in (AC.13)
therefore consists of within-row dependence rather than one-bit bias.

*Proof.*  Apply Lemma AC.1 with `c_j=2t`.  The KL bound follows by summing
the one-bit chain increments.  Total correlation is

```math
D\!\left(q_{\lambda,\epsilon}\middle\Vert
 \bigotimes_i(q_{\lambda,\epsilon})_{R_i}\right)
=D(q_{\lambda,\epsilon}\Vert U_B)
 -\sum_iD((q_{\lambda,\epsilon})_{R_i}\Vert U_n),
```

which proves (AC.13). `square`

Equation (AC.10) is the first strict resolution of one half of the previous
SML: **conditional component complexity is always tight** for this natural
filtration.  It uses the actual parent partition function, not a conference
surrogate.  It does not make the rows independent.

## 4. The optimal row-product shadow and an exact dependence charge

For a fixed orientation write

```math
V_\lambda
=-{1\over\lambda}\log\mathbb E_{U_B}e^{-\lambda L_\epsilon}.
                                                                  \tag{AC.14}
```

Let `P_row` be the set of all independent-row laws
`p=p_1 tensor ... tensor p_m` on the bridge, and define the best row-product
variational value

```math
V_\lambda^{\rm row}
=\min_{p\in\mathcal P_{\rm row}}
 \left\{\mathbb E_pL_\epsilon
       +{1\over\lambda}D(p\Vert U_B)\right\}.           \tag{AC.15}
```

**Theorem AC.3 (directed product-shadow identity).**  One has the exact
identity

```math
\boxed{
\lambda(V_\lambda^{\rm row}-V_\lambda)
=\inf_{p\in\mathcal P_{\rm row}}
 D(p\Vert q_{\lambda,\epsilon})
=:\mathcal I_\lambda^{\leftarrow}.}                    \tag{AC.16}
```

Every minimizing product shadow `p^*=tensor_i p_i^*` satisfies the
coordinate best-response equation

```math
{dp_i^*\over dU_n}(r)
\propto
\exp\{-\lambda\mathbb E_{p_{-i}^*}L_\epsilon(r,R_{-i})\},           \tag{AC.17}
```

and hence

```math
\boxed{D_2(p_i^*\Vert U_n)
\le n\log(1+\tanh^2(\lambda t))
\le\lambda^2t^2n
=\lambda^2\beta^2{n\over N}\quad(1\le i\le m).}      \tag{AC.18}
```

Thus the best independent-row competitor is automatically in the bounded
component-`D_2` class; bounded row complexity is a conclusion, not an
assumption.  The exact price of all dependence omitted by that class is the
reverse information projection `I_lambda^leftarrow`.

For comparison, the ordinary row total correlation is the opposite directed
projection,

```math
\operatorname{TC}_\epsilon
=\inf_{p\in\mathcal P_{\rm row}}
  D(q_{\lambda,\epsilon}\Vert p).                       \tag{AC.19}
```

Both directed resources are at most extensive:

```math
\begin{aligned}
\operatorname{TC}_\epsilon
&\le {\lambda^2t^2mn\over2},\\
\mathcal I_\lambda^{\leftarrow}
&\le D(U_B\Vert q_{\lambda,\epsilon})
=\lambda(\mathbb E_{U_B}L_\epsilon-V_\lambda)
\le {\lambda^2t^2mn\over2}.                            \tag{AC.20}
\end{aligned}
```

*Proof.*  The Gibbs identity gives, for every law `p`,

```math
D(p\Vert q_{\lambda,\epsilon})
=\lambda\left[
 \mathbb E_pL_\epsilon+{1\over\lambda}D(p\Vert U_B)-V_\lambda
 \right].                                               \tag{AC.21}
```

Minimize over row products to obtain (AC.16).  A minimizer exists by
compactness.  Holding all but row `i` fixed, strict convexity of entropy gives
(AC.17).  Its effective row potential retains the bit-flip bound `2t`, so
Lemma AC.1 proves (AC.18).

For the last inequality in (AC.20), bounded differences under `U_B` gives

```math
\log\mathbb E_{U_B}
e^{s(L_\epsilon-\mathbb E L_\epsilon)}
\le {s^2t^2mn\over2}.                                  \tag{AC.22}
```

Substitute `s=-lambda`. `square`

This identity gives a real, exhaustive dichotomy.  Put

```math
G_\lambda
=\mathbb E_{U_B}L_\epsilon-V_\lambda.                  \tag{AC.23}
```

Since `U_B` is a row-product candidate,

```math
\boxed{
G_\lambda
=\underbrace{\mathbb E_{U_B}L_\epsilon
              -V_\lambda^{\rm row}}_{\text{bounded-}D_2
                 \text{ row-product gain}}
+{1\over\lambda}
 \underbrace{\mathcal I_\lambda^{\leftarrow}}_
             {\text{irreducible directed row dependence}}.}       \tag{AC.24}
```

Both terms are nonnegative.  Therefore, if the actual child-induced escort
has a fixed linear gain `G_lambda>=eta N`, then at least one of the following
holds:

1. a row-product law whose every row obeys the uniform bound (AC.18) gains at
   least `eta N/2`;
2. `I_lambda^leftarrow >= lambda eta N/2`.

This is stronger than merely saying that some information resource must be
large: it identifies the product-side class and makes the dependence-side
quantity exactly equal to a variational pressure gap.

It is important not to replace `I_lambda^leftarrow` by total correlation.
They are opposite KL projections, and no dimension-free comparison follows
from the present hypotheses.  Proving that a target-reaching actual escort
has `TC=Omega(N)` would be a further theorem, not a consequence of (AC.16).

## 5. Effective coordinate support of any fixed linear escort gain

Let `d_j` be the one-bit conditional KL increments of
`q_(lambda,epsilon)` in any bridge reveal order, and, when their sum is
nonzero, define

```math
s_*(q)= {\left(\sum_j\sqrt{d_j}\right)^2
          \over \sum_jd_j}.                              \tag{AC.25}
```

By (AC.8), `d_j<=kappa(lambda t)`, and therefore

```math
s_*(q)\ge {D(q\Vert U_B)\over\kappa(\lambda t)}.        \tag{AC.26}
```

Concavity of
`phi(lambda)=-log E exp(-lambda L)` gives
`E_qL=phi'(lambda)<=phi(lambda)/lambda=V_lambda`.  Combining (AC.22) with
entropy duality yields

```math
\mathbb E_UL-\mathbb E_qL
\le\sqrt{2t^2mn\,D(q\Vert U_B)}.                        \tag{AC.27}
```

Hence, if `G_lambda>=eta N`, then

```math
D(q\Vert U_B)
\ge {\eta^2N^2\over2t^2mn},
\qquad
s_*(q)
\ge {\eta^2N^2
       \over2t^2mn\,\kappa(\lambda t)}.                 \tag{AC.28}
```

For comparable splits and fixed `beta,lambda`, these are respectively
`Omega(N)` and `Omega(N^2)`.  Thus a real phase is both linearly informative
and quadratically supported even though each conditional row has bounded
`D_2`.

This support conclusion is a general resource inequality.  It does not use
child optimality and should not be advertised as new optimizer structure.

## 6. The canonical child-spin latent does not productize the inverse escort

There is a useful exact falsifier of the most immediate latent-product hope.
Let `rho=tanh t`.  In the forward law `Pi_epsilon`, condition on the child
Gibbs variables.  The bridge bits are then independent with means
`rho q_(ij)`, where `q_(ij)=tau_1 x_i y_j`.

Take distinct rows `i,i'` and columns `j,j'`.  Conditional on `epsilon`,

```math
\begin{aligned}
\mathbb E_\Pi[B_{ij}B_{ij'}]
 &=\rho^2 c_{jj'},\\
\mathbb E_\Pi[B_{i'j}B_{i'j'}]
 &=\rho^2 c_{jj'},\\
\mathbb E_\Pi[B_{ij}B_{ij'}B_{i'j}B_{i'j'}]
 &=\rho^4,                                             \tag{AC.29}
\end{aligned}
```

where `c_(jj')=E[y_jy_(j')|epsilon]`.  At finite `t` the conditional child
law has full support, so `|c_(jj')|<1`.  If the two bridge rows were
independent, the last line would equal `rho^4 c_(jj')^2`, a contradiction.
Thus `Pi_epsilon`, and therefore its density `p_epsilon`, is not row-product
when `m,n>=2`.

Since `q_(lambda,epsilon) proportional p_epsilon^(-lambda)`, a row-product
inverse escort would force `p_epsilon` itself to be row-product.  Hence

```math
\boxed{q_{\lambda,\epsilon}\text{ is not row-product for }m,n\ge2,
\ t,\lambda>0.}                                        \tag{AC.30}
```

Nor does coupling `q_lambda` to the canonical child-spin posterior repair
this.  If `Z` denotes those latent spins, then

```math
q_\lambda(B\mid Z,\epsilon)
\propto \Pi(B\mid Z,\epsilon)
          p_\epsilon(B)^{-(\lambda+1)}.                 \tag{AC.31}
```

The first factor is product over bridge bits.  If (AC.31) were row-product,
the second factor, and hence `p_epsilon`, would be row-product, contradicting
(AC.29).  The obvious Gibbs latent therefore does not place the negative
escort inside the archived common-latent product class.

This is a finite exact obstruction, not an extensive lower bound:
(AC.29) alone does not imply that either directed dependence in
(AC.19)/(AC.16) is `Omega(N)`.

## 7. What actual child optimality adds

The results above use the actual child partition function but, except for
the selection of `A,D`, do not use their minimizing property.  The exact
optimizer-specific constraint currently available is as follows.  Let

```math
\nu_A(x,\tau)
={e^{t\tau H_A(x)}\over2^{m+1}\overline Z_m(A,t)}.
```

For every set `S` of child edges, flip precisely those signs.  Minimality of
`A` gives

```math
\boxed{
\mathbb E_{\nu_A}\exp\left{
-2t\tau\sum_{e\in S}a_ex_e\right}\ge1.}               \tag{AC.32}
```

In particular, for every edge `e={u,v}`,

```math
\boxed{a_e\,\mathbb E_{\nu_A}[\tau x_ux_v]\le\tanh t.} \tag{AC.33}
```

The same statements hold for `D`.

Indeed, the left side of (AC.32) is exactly the ratio of the flipped and
unflipped child partition functions.  For one edge it is
`cosh(2t)-a_e E[\tau x_ux_v]sinh(2t)>=1`, which is (AC.33).

Equations (AC.32)--(AC.33) are genuine actual-minimizer structure.  At
present they do not bound either term in (AC.24); claiming otherwise would
be circular.

## 8. The narrowed SML and evidentiary judgment

The previous question asked whether conditional component `D_2` is tight or
escapes.  Theorem AC.2 resolves it completely for the natural row filtration:
it is uniformly tight for every actual child pair.  The remaining issue is
not component density but whether row interaction buys a linear variational
advantage.

The new smallest missing lemma is:

> **Actual-child row-shadow lemma.**  For contracted-temperature minimizing
> children and every fixed `beta,lambda`, decide which term in (AC.24) can be
> linear.  Equivalently, prove either
> `E_U L-V_lambda^row=o(N)` (bounded-complexity independent rows cannot make
> the phase) and classify `I_lambda^leftarrow`, or prove
> `I_lambda^leftarrow=o(N)` and reduce the phase entirely to the explicit
> bounded-`D_2` row-product variational problem.

This is strictly narrower than “understand the child-induced law”: all
conditional Renyi complexity, KL scale, coordinate support, row-product
projection, and the exact interaction remainder have been identified.

It is not yet a Level-5-to-6 route.  The conference no-gain theorems close
the row-product term only for their conference pressure observable; applying
them verbatim to `L_epsilon` for arbitrary optimizing children would be a
surrogate-to-target overclaim.  A new use of the minimizer inequalities
(AC.32), or an actual-child bound on one term of (AC.24), is still required.
