# Actual-child flip averaging: exact contraction and a radial-information ceiling

Status: **task-local rigorous theorem and ceiling note**.  This note starts
from (AC.32) for an actual pressure-minimizing child.  It derives the complete
inhomogeneous Bernoulli contraction identity, its radial and fixed-size
consequences, and an exact order-eight witness showing that the radial data do
not determine even two-replica overlap geometry or the response to a
one-vertex rank-one extension.  Thus the averaged constraints are genuine
optimizer structure, but their homogeneous part does not control the missing
row-product shadow in (AC.24).

The notation agrees with
`extremal_information/drafts/actual_child_negative_escort_structure.md`.

## 1. Child pressure and its multiaffine normalization

Let `E=binom([m],2)`, `K=|E|`, and let `A=(a_e)_(e in E)` be a hollow sign
matrix.  For an inhomogeneous coupling vector `s=(s_e)`, set

```math
 Z_A(s)
 =\mathbb E_{x,\tau}\exp\left\{
     \tau\sum_{e=\{u,v\}}s_ea_ex_ux_v\right\}.
 \tag{FC.1}
```

Write `Z_A(t)=Z_A(t 1)`, `psi_A(t)=log Z_A(t)`, and `rho=tanh t`.  Define

```math
 P_A(r)
 =\mathbb E_{x,\tau}\prod_{e\in E}
       (1+r_e\tau a_ex_e),
 \qquad x_e=x_ux_v.                                    \tag{FC.2}
```

Then, coordinatewise,

```math
 Z_A(s)=\prod_e\cosh(s_e)\,P_A(\tanh s).              \tag{FC.3}
```

Suppose throughout Sections 1--4 that `A` minimizes `Z_C(t)` over all
signings `C` at this fixed `t>0`.  With

```math
 \nu_A(x,\tau)
 ={\exp(t\tau H_A(x))\over 2^{m+1}Z_A(t)},
 \qquad Y_e=\tau a_ex_e,                               \tag{FC.4}
```

minimality is exactly (AC.32): for every `S subset E`,

```math
 \mathbb E_{\nu_A}\exp\{-2t\sum_{e\in S}Y_e\}
 ={Z_{A^S}(t)\over Z_A(t)}\ge1.                       \tag{FC.5}
```

Here `A^S` flips precisely the signs indexed by `S`.

## 2. The strongest Bernoulli average is an exact contraction principle

Choose each edge independently for `S`, with an edge-dependent probability
`p_e`.  Put

```math
 r_e=(1-2p_e)\rho,
 \qquad s_e=\operatorname{arctanh}(r_e).               \tag{FC.6}
```

The following identity is exact, before using minimality.

**Theorem FC.1 (inhomogeneous optimizer contraction).**  One has

```math
 \begin{aligned}
 \mathbb E_S{Z_{A^S}(t)\over Z_A(t)}
 &=\mathbb E_{\nu_A}\prod_{e\in E}
       \bigl[(1-p_e)+p_e e^{-2tY_e}\bigr]\\
 &=\left(\prod_e{\cosh t\over\cosh s_e}\right)
       {Z_A(s)\over Z_A(t)}
 ={P_A(r)\over P_A(\rho\mathbf1)}.                    \tag{FC.7}
 \end{aligned}
```

Consequently an actual minimizing child obeys

```math
 \boxed{
 {Z_A(s)\over\prod_e\cosh s_e}
 \ge {Z_A(t\mathbf1)\over(\cosh t)^K}
 \quad\text{for every }s\in[-t,t]^E.}                \tag{FC.8}
```

Moreover, (FC.8), the full inhomogeneous Bernoulli family in (FC.7), the
edge-set inequalities (FC.5), and exact minimization over all sign flips are
equivalent.

*Proof.*  For `y in {+-1}`,

```math
 e^{ty}\bigl[(1-p_e)+p_e e^{-2ty}\bigr]
 =(1-p_e)e^{ty}+p_e e^{-ty}
 ={\cosh t\over\cosh s_e}e^{s_ey},                   \tag{FC.9}
```

because `tanh s_e=(1-2p_e)tanh t`.  Multiplication and averaging prove
(FC.7), and averaging (FC.5) proves (FC.8).  Conversely, setting every
`p_e` to zero or one recovers every vertex `A^S`.  Equivalently, `P_A` is
multiaffine in `r`, so its minimum on the box `[-rho,rho]^E` occurs at a
vertex. `square`

This is the strongest direct consequence obtainable merely by averaging
(AC.32) with arbitrary independent edge probabilities: it has lost no
information at all, because its boundary values are the original
minimization problem.  It is genuinely optimizer-specific, but it is an exact
re-expression of minimization rather than a strict reduction.

## 3. Homogeneous noise, pressure tangent, and child entropy

For homogeneous `p_e=p`, put `q=1-2p` and

```math
 s=\operatorname{arctanh}(q\tanh t).
```

Then (FC.7) becomes the radial disorder-noise identity

```math
 \mathbb E_{S\sim\operatorname{Ber}(p)^E}
 {Z_{A^S}(t)\over Z_A(t)}
 ={P_A(q\rho\mathbf1)\over P_A(\rho\mathbf1)}\ge1.  \tag{FC.10}
```

Equivalently,

```math
 \boxed{
 \psi_A(t)-K\log\cosh t
 \le \psi_A(s)-K\log\cosh s,
 \qquad |s|\le t.}                                   \tag{FC.11}
```

Three explicit consequences are:

```math
 \begin{aligned}
 \psi_A(t)&\le K\log\cosh t,                         &&\tag{FC.12}\\
 0\le\psi_A'(t)&\le K\tanh t,                       &&\tag{FC.13}\\
 D(\nu_A\Vert U_{x,\tau})
 &=t\psi_A'(t)-\psi_A(t)
 \le Kt\tanh t.                                      &&\tag{FC.14}
 \end{aligned}
```

Here (FC.12) is the usual annealed-sign comparison (`p=1/2`), and the upper
tangent in (FC.13) follows by taking the left derivative of (FC.11) at `t`.
The lower tangent follows from evenness and convexity of `psi_A`.  Finally
`psi_A(t)>=0` by `cosh(tH)>=1`, which proves (FC.14).  Coordinatewise
differentiation of (FC.8) at its corner gives the stronger individual
tangents

```math
 a_e\mathbb E_{\nu_A}[\tau x_e]\le\tanh t,            \tag{FC.15}
```

which are exactly (AC.33).

At the contracted scale `t=beta/sqrt(N)` and comparable child sizes,
(FC.14) is only an `O_beta(N)` entropy-deficit bound.  It gives no `o(N)`
dependence estimate.

## 4. Fixed-size flips and the exact radial Fourier content

Expanding (FC.2), the only surviving edge sets are even-cardinality Eulerian
subgraphs:

```math
 P_A(r\mathbf1)
 =\sum_{\substack{F\subseteq E:\ \deg_F(v)\ \mathrm{even}\ \forall v,\\
                              |F|\ \mathrm{even}}}
       a_F r^{|F|}
 =\sum_{\ell=0}^K W_\ell(A)r^\ell,                   \tag{FC.16}
```

where `a_F=prod_(e in F)a_e`.  Let

```math
 k_k(\ell)
 ={1\over\binom Kk}\sum_j(-1)^j
   \binom\ell j\binom{K-\ell}{k-j}                  \tag{FC.17}
```

be the normalized Krawtchouk multiplier.  Uniformly averaging (FC.5) over
all `S` of size `k` gives

```math
 \boxed{
 {1\over\binom Kk}\sum_{|S|=k}{Z_{A^S}(t)\over Z_A(t)}
 ={\sum_\ell W_\ell(A)\rho^\ell k_k(\ell)
    \over P_A(\rho\mathbf1)}\ge1.}                  \tag{FC.18}
```

Equivalently, directly in the child Gibbs law, the left side is

```math
 \mathbb E_{\nu_A}{1\over\binom Kk}
 [z^k](1+ze^{-2t})^{(K+U)/2}
       (1+ze^{2t})^{(K-U)/2},
 \quad U=\tau H_A(x).                                 \tag{FC.19}
```

The Krawtchouk transform is invertible.  Hence, at any fixed `t>0`, the
complete list of *values* in (FC.18), not merely their lower bounds, is
equivalent to the radial coefficient list `(W_ell(A))`.  In turn

```math
 P_A(\tanh u\,\mathbf1)={Z_A(u)\over(\cosh u)^K},     \tag{FC.20}
```

so this is equivalent to the absolute-energy histogram of `H_A` (the
functions `cosh(hu)` for distinct `|h|` are linearly independent).

Thus homogeneous Bernoulli averaging and fixed-size averaging retain exactly
one radial object: the absolute-energy histogram.  They do determine scalar
quantities such as all pressure derivatives and the total Fourier norm of
the `x`-marginal Gibbs density

```math
 g_t(x)={\cosh(tH_A(x))\over Z_A(t)},
 \qquad
 \mathbb E_U g_t^2={Z_A(2t)+1\over2Z_A(t)^2}.          \tag{FC.21}
```

They do **not** determine how this Fourier mass is distributed among levels
or coordinates.

## 5. An exact actual-minimizer ceiling witness

The certified exhaustive order-eight classification in
`computations/results/m8_minimizer_orbits.json` contains two inequivalent
classes.  Representatives are

```text
 A0 =
  0  1  1  1  1  1  1  1
  1  0  1  1  1  1 -1 -1
  1  1  0  1 -1  1 -1  1
  1  1  1  0  1 -1 -1  1
  1  1 -1  1  0 -1  1 -1
  1  1  1 -1 -1  0  1 -1
  1 -1 -1 -1  1  1  0  1
  1 -1  1  1 -1 -1  1  0

 A1 =
  0  1  1  1  1  1  1  1
  1  0  1  1 -1 -1  1 -1
  1  1  0 -1  1 -1  1 -1
  1  1 -1  0  1 -1 -1  1
  1 -1  1  1  0 -1 -1  1
  1 -1 -1 -1 -1  0  1  1
  1  1  1 -1 -1  1  0  1
  1 -1 -1  1  1  1  1  0.
```

Exact enumeration gives the common projective absolute-energy histogram

| `|H|` | 0 | 2 | 4 | 6 | 8 | 10 |
|---:|---:|---:|---:|---:|---:|---:|
| count | 12 | 32 | 32 | 24 | 20 | 8 |

and both have cap `10=M_8`.  The exhaustive classification proves that every
other cap-10 signing lies in one of these two classes and has this same
histogram.  Every remaining signing has cap at least 12.  Therefore both
classes minimize `Z_A(t)` for all sufficiently large `t`; explicitly `t>=3`
suffices, because

```math
 {\cosh(12t)\over128}>\cosh(10t)\qquad(t\ge3).        \tag{FC.22}
```

The left side lower-bounds the pressure of a cap-at-least-12 signing, while
the right side upper-bounds that of either displayed cap-10 signing.

Nevertheless their nonradial Gibbs geometry differs.  Let

```math
 \mu_{A,t}(x)={\cosh(tH_A(x))\over2^8Z_A(t)},
 \qquad C_A(t)=\mathbb E_{\mu_{A,t}}xx^\top.           \tag{FC.23}
```

At zero temperature, exact ground-state enumeration gives

```math
 \begin{array}{c|c|c}
 &\operatorname{spec}C_A(\infty)&\operatorname{Tr}C_A(\infty)^2\\ \hline
 A_0& (1/2)^{\times6},(5/2)^{\times2}&14\\
 A_1& (1/2)^{\times4},(3/2)^{\times4}&10.
 \end{array}                                           \tag{FC.24}
```

In replica notation,

```math
 \operatorname{Tr}C_A(t)^2
 =\mathbb E_{x,x'\sim\mu_{A,t}}\langle x,x'\rangle^2. \tag{FC.25}
```

It is also `8+2 sum_(i<j) \hat g_t({i,j})^2`, so (FC.24) proves that the
same total norm (FC.21) can have different level-two Fourier mass and
different pair-overlap geometry.

There is an even more direct rank-one response separation.  Attach one new
spin with incident signs `b in {+-1}^8`, and define

```math
 K_A(b)=\max_x\left|H_A(x)+\sum_i b_ix_i\right|.       \tag{FC.26}
```

Exact enumeration over all 256 fields gives

| response | 12 | 14 | 16 | 18 |
|---:|---:|---:|---:|---:|
| number for `A0` | 24 | 112 | 104 | 16 |
| number for `A1` | 8 | 112 | 120 | 16 |

This is the smallest possible bridge to a new child: a one-column rank-one
channel.  Since the covariance limits in (FC.24) and the response limits in
(FC.26) differ, continuity and the zero-temperature Laplace principle imply
that, for all sufficiently large finite `t`, these two actual pressure
minimizers have different overlap geometry and different one-vertex bridge
pressure landscapes despite having identical pressure, entropy, every
homogeneous Bernoulli average, and every fixed-size flip average.

The exact finite calculations are independently reproducible with
`extremal_information/experiments/actual_child_radial_ceiling_witness.py`.

## 6. Consequences for the row-product shadow

The conclusion is deliberately limited but rigorous.

1. The full inhomogeneous contraction (FC.8) is equivalent to exact child
   minimization.  Invoking all of it without an additional theorem is not a
   strict reduction.
2. Its homogeneous and fixed-size radializations retain only the absolute
   energy histogram.  The actual-minimizer witness above proves that this
   information does not determine the child overlap law or even the simplest
   rank-one extension response.
3. Hence no derivation which uses only (FC.10)--(FC.21) can reconstruct the
   negative-disorder bridge law, establish rank-one channel resolvability, or
   identify either summand in the exact split (AC.24).  A universal upper
   bound could still exist, but it must use additional nonradial structure.
4. Tight autoregressive row `D_2` from (AC.10) does not supply that missing
   structure: it is not a common-latent product representation and does not
   transfer to a latent-iid no-gain theorem.  Moreover, the archived
   conference-child component estimate is specific to its projected
   conference observable; even an actual-child latent-iid representation
   would still require a new component no-gain theorem for `L_epsilon`.

The optimizer-specific theorem obtained here is therefore (FC.8), and the
sharp ceiling is (FC.24)--(FC.26).

## 7. Two primary-source mean-field checks

### 7.1 Augeri gives a legal `O(N)`, not `o(N)`, bound

Augeri's Theorem 2.6 in
[A transportation approach to the mean-field approximation](https://arxiv.org/abs/1903.08021)
states, for the uniform law on a `d`-dimensional Boolean cube, that the
bit-product mean-field error is at most a universal constant times the
Rademacher width `b(nabla f)` of the discrete-gradient set (using the harmonic
extension in Remark 2.2).

For the bridge potential `f(B)=-lambda L(B)`, `d=mn`, the natural continuous
gradient satisfies

```math
 \nabla f(B)=-\lambda t\,\mathbb E_B Q,
 \qquad Q\in\{\pm xy^\top:x\in\{\pm1\}^m,
                         y\in\{\pm1\}^n\}.            \tag{FC.27}
```

For a Rademacher matrix `Xi`, this gives

```math
 \mathbb E\sup_B\langle\nabla f(B),\Xi\rangle
 \le\lambda t\sqrt{mn}\,\mathbb E\|\Xi\|_{op}
 \le C\lambda t\sqrt{mn}(\sqrt m+\sqrt n).           \tag{FC.28}
```

The continuous diagonal Hessian obeys
`|partial_(ee)^2 f|<=lambda t^2`.  Therefore each discrete half-difference
differs from the continuous endpoint derivative by at most `lambda t^2`, and

```math
 b(\nabla^{\rm disc}f)
 \le C\lambda t\sqrt{mn}(\sqrt m+\sqrt n)
      +\lambda t^2mn.                                  \tag{FC.29}
```

For comparable splits and `t=beta/sqrt(N)`, this is `O_(beta,lambda)(N)`.
Bit-products are a subclass of row-products, so Augeri legally implies only

```math
 \mathcal I^{\leftarrow}_{\rm row}
 \le \mathcal I^{\leftarrow}_{\rm bit}=O(N),          \tag{FC.30}
```

not the `o(N)` needed to close (AC.24).  This estimate is universal and does
not use child minimality.

### 7.2 Continuous strong log-concavity does not transfer to the cube

The main theorem of Lacker--Mukherjee--Yeung,
[Mean field approximations via log-concavity](https://arxiv.org/abs/2206.01260),
concerns a `C^2`, strongly log-concave probability density on `R^d` (or a
strongly log-concave continuous product reference with a concave tilt).  Its
reverse-KL product gap is bounded using continuous conditional gradient
variances and cross Hessian squares.  It is not a theorem for atomic Boolean
laws.

There is a decisive extension falsifier.  Start with any smooth extension
`h` of an arbitrary Boolean potential, chosen with globally bounded Hessian
(a cutoff outside a neighborhood of the cube supplies one).  For `c` large,

```math
 h_c(b)=h(b)-c\|b\|_2^2                              \tag{FC.31}
```

is arbitrarily strongly concave on `R^d`, while on every Boolean vertex
`h_c(b)=h(b)-cd`.  Thus it induces **exactly the same discrete Gibbs law** as
`h`.  Indeed the Lacker--Mukherjee--Yeung continuous error can be driven down
by increasing the diagonal curvature, while the arbitrary discrete
reverse-KL product gap is unchanged.  The continuous measure has moved into
the interior and is not a realization of the Boolean law.

Consequently continuous strong log-concavity cannot control (AC.24) without
a new canonical-extension/rounding theorem which transfers both the
variational value and the product class.  No such theorem is supplied by the
paper.

## 8. Narrowed missing lemma

All consequences obtained by homogeneous or fixed-size averaging are now
classified and sharply falsified as channel-sufficient data.  The next lemma
must use the *nonradial* content of the exact contraction box without merely
restating all vertex inequalities:

> **Nonradial actual-child channel lemma.**  Extract from (FC.8), for
> contracted-temperature minimizing children, a statistic strictly smaller
> than the full sign-flip landscape which makes one summand in (AC.24)
> `o(N)` and supplies a separate no-gain/classification theorem for the
> remaining summand, or else prove that every such controlling statistic has
> extensive information complexity.

The order-eight witness closes radial data as a parameter-free exact
sufficient statistic.  It does not by itself exclude a new large-order,
small-`t=beta/sqrt(N)` rigidity theorem for actual optimizers.  Such a theorem
would be additional asymptotic structure, not a consequence of the averaged
identities above.  The full inhomogeneous route remains logically available,
but at present it is equivalent to exact minimization and therefore is not
yet a Level-5-to-6 mechanism.
