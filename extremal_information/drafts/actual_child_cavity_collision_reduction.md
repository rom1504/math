# The actual-child row interaction as a collision--cavity partition

Status: **rigorous task-local reduction, sufficient synchronization theorem,
and scalable generic falsifier**.  The canonical interaction from
[`actual_child_cross_row_response_decomposition.md`](actual_child_cross_row_response_decomposition.md)
has an exact representation as the log partition function of iid row
factors acting on the two child Gibbs laws.  This integrates out the right
child explicitly, separates right-child posterior collision from left-child
cavity response without paying them separately, and gives a scalar
projective-synchronization condition implying `J=o(N)`.

The condition is not proved for optimizing children.  A centrally symmetric
rank-one channel has weak physical-scale coordinates and bounded row
Renyi-two complexity but a linear canonical cumulant.  Thus child optimality,
not rank one or weak noise alone, must drive any closure theorem.

## 1. Exact iid-row factor representation

Fix child orders `m,n`, internal amplitude `t`, channel amplitude `u`, and
orientation `epsilon`.  Use Proposition CR.0's notation and put

```math
 \nu_\epsilon(s,x,y)
 =\pi_s^{(\epsilon)}\mu_{A,s}(x)\mu_{D,\epsilon s}(y).
                                                               \tag{CC.1}
```

For a row word `b in {+-1}^n`, a latent state `(s,x,y)`, and
`rho=tanh u`, define

```math
 \begin{aligned}
 k_u(b\mid s,x_i,y)
  &:=\prod_{j=1}^n(1+\rho b_jsx_iy_j)
    ={e^{usx_i\langle b,y\rangle}\over(\cosh u)^n},\\
 z_u(b)&:=E_{\nu_\epsilon}k_u(b\mid s,X_i,Y)
    ={E_{s,Y}\cosh(u\langle b,Y\rangle)\over(\cosh u)^n},\\
 \ell_{u,b}(s,x_i,y)&:={k_u(b\mid s,x_i,y)\over z_u(b)}.
                                                               \tag{CC.2}
 \end{aligned}
```

The second line is independent of `i`, including when `m!=n`, because each
`X_i` is fair conditional on `s`.  It is precisely the row likelihood
`p_(row,u)` in CR.0.  Let

```math
 G_{A,D,u}(B)
 :=E_{\nu_\epsilon}\prod_{i=1}^m
                  \ell_{u,B_i}(s,X_i,y).             \tag{CC.3}
```

**Theorem CC.1 (collision--cavity representation).**  For every finite pair
of children,

```math
 \boxed{
 {p_u(B)\over\prod_i p_{i,u}(B_i)}=G_{A,D,u}(B),
 \qquad h_u(B)=\log G_{A,D,u}(B).}                   \tag{CC.4}
```

If

```math
 {dr_{\rm row,u}\over dU_n}(b)
 ={z_u(b)^{-\lambda}\over Z_{\rm row,u}},
 \qquad Z_{\rm row,u}=E_{U_n}z_u^{-\lambda},        \tag{CC.5}
```

then the rows `B_i` are iid under the canonical certificate and

```math
 \boxed{
 \mathcal J_u
 =\log E_{r_{\rm row,u}^{\otimes m}}G_{A,D,u}^{-\lambda}
  +\lambda E_{r_{\rm row,u}^{\otimes m}}\log G_{A,D,u}.}
                                                               \tag{CC.6}
```

Equivalently, `q_u/r_u` is `G^(-lambda)` divided by its expectation.

*Proof.*  Conditional on `(s,x,y)`, the channel rows are independent and
their likelihood ratios relative to fair rows are the first line of
(CC.2).  Averaging their product gives `p_u(B)`.  Each factor's marginal
average is `z_u(B_i)=p_(i,u)(B_i)`.  Division proves (CC.4).  Equations
(CC.5)--(CC.6) follow by substituting `p_u=G prod_i z_u` in the definitions
of `r_u`, `q_u`, and `D(r_u||q_u)`. `square`

This representation uses the two child Gibbs measures and a table of
`2^n` row factors, not the `2^(mn)` parent bridge-pressure table.  It does
not use the value of the target-order optimization.

## 2. Integrating the right child

Define its external-field log partition

```math
 \Lambda_{D,a}(h):=\log E_{Y\sim\mu_{D,a}}e^{\langle h,Y\rangle}.
                                                               \tag{CC.7}
```

All powers of `cosh u` cancel in (CC.3), giving the explicit formula

```math
 \boxed{
 G_{A,D,u}(B)
 ={\displaystyle
    \sum_{s=\pm1}\pi_s^{(\epsilon)}
      E_{X\sim\mu_{A,s}}
       \exp\{\Lambda_{D,\epsilon s}(usB^{\mathsf T}X)\}
   \over\displaystyle
    \prod_{i=1}^m\left[
      \sum_{s=\pm1}\pi_s^{(\epsilon)}
       E_{Y\sim\mu_{D,\epsilon s}}
        \cosh(u\langle B_i,Y\rangle)
    \right]}.}                                      \tag{CC.8}
```

Thus `J=o(N)` is a lower-tail concentration question for one random
left-child cavity free energy under iid rows from the explicit law (CC.5).
It is not a minimization over parent bridges.

There is a useful exact two-resource decomposition that retains their joint
cancellation.  On `v=(s,y)` put

```math
 \begin{aligned}
 \zeta(v)&=\pi_s^{(\epsilon)}\mu_{D,\epsilon s}(y),\\
 w_b(v)&={\cosh(u\langle b,y\rangle)\over(\cosh u)^n},
 &K_b(v)&={w_b(v)\over z_u(b)},\\
 a_b(v)&=s\tanh(u\langle b,y\rangle),\\
 F_{A,s}(a_1,\ldots,a_m)
   &:=E_{X\sim\mu_{A,s}}\prod_i(1+X_i a_i).
                                                               \tag{CC.9}
 \end{aligned}
```

Let

```math
 C_D(B)=E_\zeta\prod_iK_{B_i}(v),
 \qquad
 {d\zeta_B^\Delta\over d\zeta}(v)
 ={\prod_iK_{B_i}(v)\over C_D(B)}.                  \tag{CC.10}
```

Then

```math
 \boxed{
 G_{A,D,u}(B)
 =C_D(B)\,
   E_{v\sim\zeta_B^\Delta}
    F_{A,s}(a_{B_1}(v),\ldots,a_{B_m}(v)).}          \tag{CC.11}
```

Here `C_D` is the diagonal collision of the `m` right-channel row
posteriors conditional on the shared sector prior (the sector weights also
contain `Z_A^s`), while the second factor is the left-child cavity response
in that collided posterior.  Equation (CC.11) is an identity; bounding the
two logarithms separately would discard potentially leading cancellation
and is not asserted.

Global-flip symmetry gives the exact expansion

```math
 F_{A,s}(a)=
 \sum_{S\subseteq[m]}
 E_{\mu_{A,s}}X_S\prod_{i\in S}a_i,
 \qquad E_{\mu_{A,s}}X_S=0\quad(|S|\text{ odd}).     \tag{CC.12}
```

Thus ordinary pair overlap is not algebraically sufficient: without a new
synchronization theorem, the row environment can query every even child
correlation.  The optimizer inequality AC.33 controls one sector-weighted
pair response only and does not by itself bound (CC.12).

This response has two exact, named differential forms.  For
`v_i=atanh(a_i)` put

```math
 \mathscr L_{A,s}(v)
 :=\log E_{\mu_{A,s}}e^{\langle v,X\rangle}
   -\sum_i\log\cosh v_i
 =\log F_{A,s}(\tanh v).                             \tag{CC.12a}
```

If `m_(A,s)(v)=E_(mu_(A,s),v)X` is the magnetization after adding external
field `v`, then

```math
 \boxed{
 \mathscr L_{A,s}(v)
 =\int_0^1\langle v,
     m_{A,s}(rv)-\tanh(rv)\rangle\,dr.}              \tag{CC.12b}
```

Alternatively, let `nu_(q,v)` have density proportional to
`exp{qstH_A(x)+<v,x>}` relative to fair spins.  Direct interpolation in the
child coupling gives

```math
 \boxed{
 \mathscr L_{A,s}(v)
 =st\int_0^1\sum_{a<b}A_{ab}
  \{E_{\nu_{q,v}}X_aX_b-E_{\nu_{q,0}}X_aX_b\}\,dq.} \tag{CC.12c}
```

Both formulas subtract the independent-spin response exactly.  They reduce
the desired estimate to a left-child magnetization or pair-response
interpolation, but do not truncate it automatically.

Indeed, exact evaluation on the accessible field grid is generically not a
compression.  For any fixed `s,y` and `u>0`, choosing each row `b_i=y` or
`b_i=-y` realizes the Cartesian grid
`tanh(v_i) in {+-tanh(un)}`.  Since (CC.12) is multiaffine, its values on
this `2^m` grid recover every even correlation coefficient by Walsh
inversion.  Any useful theorem must bound a quotient of this response, not
retain its entire exact grid.

## 3. A strictly smaller sufficient observable

Let `mathsf K_u nu` denote the row-channel output likelihood

```math
 (\mathsf K_u\nu)(b)=E_{Z\sim\nu}
                  \prod_j(1+\rho b_jZ_j).           \tag{CC.13}
```

Given all output rows except `i`, let `nu_i^{B_{-i}}` be the forward
posterior law of the latent word `Z_i=sX_iY`.  Let `mu_row` be its
unconditioned law from CR.0.  Define the scalar projective response

```math
 \delta_i(B_{-i})
 :=\mathop{\rm osc}_{b\in\{\pm1\}^n}
   \log{(\mathsf K_u\nu_i^{B_{-i}})(b)
            \over(\mathsf K_u\mu_{\rm row})(b)},
 \qquad
 \Delta_u^2:=\sum_i\sup_{B_{-i}}\delta_i(B_{-i})^2. \tag{CC.14}
```

This is the Hilbert/projective diameter of one smoothed latent-row
posterior relative to its unconditioned output.  It is a scalar cavity
synchronization observable, not the table of parent bridge values.

**Theorem CC.2 (projective cavity synchronization).**  For every finite
actual-child channel,

```math
 \boxed{\mathcal J_u\le {\lambda^2\over8}\Delta_u^2.} \tag{CC.15}
```

In particular, along comparable splits,

```math
 \Delta_{\beta/\sqrt N}^2=o(N)
 \quad\Longrightarrow\quad
 \mathcal J_{\beta/\sqrt N}=o(N).                   \tag{CC.16}
```

*Proof.*  With `B_(-i)` fixed, Bayes' rule and conditional channel
independence give

```math
 {G(b,B_{-i})\over G(b',B_{-i})}
 ={(\mathsf K_u\nu_i^{B_{-i}})(b)/
       (\mathsf K_u\mu_{\rm row})(b)
   \over
   (\mathsf K_u\nu_i^{B_{-i}})(b')/
       (\mathsf K_u\mu_{\rm row})(b')}.             \tag{CC.17}
```

Thus the coordinate range of `log G` in row `i` is exactly
`delta_i(B_(-i))`, at most the corresponding term in `Delta_u`.  The rows
are independent under (CC.5).  The bounded-differences/Hoeffding martingale
lemma therefore gives

```math
 \log E_r e^{-\lambda(\log G-E_r\log G)}
 \le {\lambda^2\over8}
       \sum_i\sup_{B_{-i}}\delta_i(B_{-i})^2.
```

The left side is (CC.6). `square`

Knowing the single number `Delta_u` does not reconstruct any conditional
law.  However, verifying its supremum naively could still require an
exponential posterior-predictive table.  CC.2 is a noncircular sufficient
theorem, not yet an algorithmic compression and not a claim that actual
minimizers satisfy it.

The average version controls the initial interaction curvature even when
the uniform theorem is far too expensive.  Write the expectation below
under the iid canonical rows and define

```math
 \overline\Delta_u^2
 :=\sum_iE_{r_{-i}}\delta_i(B_{-i})^2.                \tag{CC.17a}
```

**Corollary CC.3 (average cavity sensitivity).**

```math
 \boxed{
 \operatorname{Var}_{r_u}(h_u)
 \le\sum_iE_{r_{-i}}
       \operatorname{Var}_{r_i}(h_u\mid B_{-i})
 \le {1\over4}\overline\Delta_u^2.}                 \tag{CC.17b}
```

*Proof.*  The first inequality is Efron--Stein for the independent rows.
For fixed other rows, (CC.17) says the range of `h_u` is
`delta_i(B_(-i))`; the variance of a variable in an interval of length
`delta` is at most `delta^2/4`. `square`

This controls the `s=0` endpoint of the interaction-curvature path IC.3.
It does not control the whole centered negative moment: an extensive
fluctuation can be created after tilting toward `q_s`.  A useful replacement
for the uniform CC.2 criterion would therefore propagate an averaged cavity
bound along that path.

## 4. A scalable rank-one falsifier

Weak physical-scale noise, central rank-one structure, and bounded row
Renyi complexity do not imply (CC.16).  Let `m=n=r`, let the planted bridge
word be `Q=sigma 1 1^T` for a fair sign `sigma`, and put

```math
 u_r={c\over\sqrt r},\qquad c>0.                     \tag{CC.18}
```

This prior can be represented using globally symmetric rank-one child
spins, but it is not asserted to be the finite-temperature law of optimized
children.  At the signing problem's physical scale `N=2r`, the notation in
(CC.18) corresponds to `c=beta/sqrt(2)`.  For `S_i=sum_jB_(ij)`,
`V_i=u_rS_i`, and
`g(v)=log cosh v`, direct channel summation gives

```math
 p_u(B)={\cosh(\sum_iV_i)\over(\cosh u_r)^{r^2}},
 \qquad
 z_u(B_i)={\cosh(V_i)\over(\cosh u_r)^r},
```

and hence

```math
 h_u(B)=g\left(\sum_iV_i\right)-\sum_i g(V_i).       \tag{CC.19}
```

Let `P_r` be the law of one `V_i` under `r_(row,u)`.  Relative to the fair
row law `U_r`,

```math
 {dP_r\over dU_r}(V)
 ={e^{-\lambda g(V)}\over Z_r},
 \qquad Z_r=E_{U_r}e^{-\lambda g(V)},               \tag{CC.20}
```

and put `kappa_r=D(P_r||U_r)>0` for `r>=2`.  Cancellation in (CC.6) yields the exact
identity

```math
 \begin{aligned}
 \mathcal J_{u_r}
 &=r\kappa_r
   +\lambda E_{P_r^{\otimes r}}
       g\left(\sum_iV_i\right)
   +\log E_{U_r^{\otimes r}}
       e^{-\lambda g(\sum_iV_i)}.                   \tag{CC.21}
 \end{aligned}
```

Indeed multiplying `P_r^(tensor r)` by
`exp(lambda sum_i g(V_i))` exactly cancels its product tilt.  Since
`0<=g(v)<=|v|`, Jensen and Cauchy--Schwarz imply

```math
 \boxed{
 \mathcal J_{u_r}\ge r\kappa_r-\lambda c\sqrt r.}   \tag{CC.22}
```

Under `U_r`, `V=c r^(-1/2)sum_jB_j` converges to `cZ` for standard Gaussian
`Z`.  The functions `e^(-lambda g)` and `g e^(-lambda g)` are bounded, so

```math
 \kappa_r\longrightarrow
 D\left(
 {\cosh(cZ)^{-\lambda}\over E\cosh(cZ)^{-\lambda}}
  d\gamma\ \middle\Vert\ d\gamma\right)
 =:\kappa(c,\lambda)>0.                              \tag{CC.23}
```

Consequently `J_(u_r)=Omega(r)`.  Yet AC.1 gives

```math
 D_2(P_r\Vert U_r)\le\lambda^2u_r^2r=\lambda^2c^2. \tag{CC.24}
```

This is a scalable falsifier to any generic implication from weak
coordinates plus tight row `D_2` to sublinear canonical interaction.  It
also shows why a left-only pair-overlap bound cannot suffice: one global
rank-one latent orbit already creates linear collision--cavity work.

There is a sharper obstruction to every fixed-order left-overlap
truncation.  Fix `k`, choose an even `ell>k`, and let `m` be divisible by
`ell`.  Partition the left spins into independent `ell`-blocks with law

```math
 \mu_\delta(x_I)=2^{-\ell}
       \left(1+\delta\prod_{i\in I}x_i\right),
 \qquad 0<|\delta|<1.                                \tag{CC.25}
```

Take a deterministic right word `y_0` and `Q=XY^T`.  Because `ell` is even,
the left law is globally flip invariant and `Q` is central.  All left
correlations of order at most `k` vanish exactly, just as for iid fair
spins.  Nevertheless, with

```math
 c_n(b)={\cosh(u\langle b,y_0\rangle)\over(\cosh u)^n},
 \qquad a_n(b)=\tanh(u\langle b,y_0\rangle),          \tag{CC.26}
```

one has the exact block factorization

```math
 {p_u(B)\over\prod_i p_{i,u}(B_i)}
 =\prod_{I}\left(1+\delta\prod_{i\in I}a_n(B_i)\right).  \tag{CC.27}
```

Under the canonical iid row law
`dr_n/dU_n proportional c_n^(-lambda)`, these block logarithms are iid, so

```math
 \mathcal J_u={m\over\ell}j_n,
 \quad
 j_n=\log E\exp\{-\lambda(g_n-Eg_n)\}>0,
 \quad
 g_n=\log\left(1+\delta\prod_{i=1}^{\ell}a_n(B_i)\right). \tag{CC.28}
```

If `m/N->theta in (0,1)` and `u=beta/sqrt(N)`, the tilted row CLT gives
`u<b,y_0> => alpha Z`, `alpha=beta sqrt(1-theta)`, under the Gaussian law
tilted by `cosh(alpha Z)^(-lambda)`.  Since `g_n` is uniformly bounded and
has a nonconstant limit, strict Jensen and bounded convergence yield

```math
 j_n\longrightarrow j_\infty>0,
 \qquad
 {\mathcal J_u\over N}\longrightarrow
 {\theta\over\ell}j_\infty>0.                       \tag{CC.29}
```

For `delta=0`, the correlations through order `k` are identical but
`G=1` and `J=0`.  Hence no fixed-order overlap summary controls the
physical-scale canonical interaction for general central rank-one priors.

## 5. Evidentiary conclusion and narrowed target

Equations (CC.8) and (CC.11) are a genuine algebraic compression from a
parent bridge table to two child Gibbs objects, but they do not by themselves
make the lower tail easy.  The full polynomial (CC.12) warns that simply
retaining all left responses would rebuild an exponentially rich child
landscape.  CC.2 instead identifies a one-number sufficient resource.

The optimizer-specific missing lemma is now:

> **Actual-child projective synchronization, or a weaker replacement.**  For
> contracted-temperature minimizing children at `u=beta/sqrt(N)`, prove
> `Delta_u^2=o(N)` in (CC.14), or prove directly that the iid-row
> collision--cavity partition (CC.8) has `o(N)` centered lower-tail
> log-MGF.

The first statement is strictly stronger than necessary but is not circular:
it concerns smoothed one-row posterior prediction and does not invoke a
target-order optimizer.  The rank-one example proves that some
optimizer-specific use of AC.32 or another rigidity theorem is indispensable.
