# Integrable frame roofs and directional conversion

**Status.**  Rigorous abstract conversion theorem for the exact rank-one
bridge channel, plus a sharp fixed-two-sided-frame specialization.  The
theorem identifies a sufficient form of frame synchronization which is
strictly stronger than pointwise low rank and strictly weaker than retaining
the complete posterior table: a subexponential **integrable upper roof** for
the bridge pressure.  Under that hypothesis, one additional one-row
observable decides the coherent-retuning branch.

The result does **not** prove that actual optimizing children possess such a
roof, nor that their one-row cost is extensive.  In particular, an
inverse-escort-only coreset for the posterior cavities does not automatically
satisfy the hypotheses below.  Those scope restrictions are part of the
statement.

## 1. Exact channel and the two product projections

Let `d=mn`, let `U` be fair measure on the `m by n` bridge cube, and let
`mu` be a centrally symmetric law on rank-one sign matrices

```math
Q=XY^{\mathsf T}\in\{-1,1\}^{m\times n},
\qquad \mu(Q)=\mu(-Q).
```

Fix `t>0` and define the bridge log likelihood, up to its irrelevant
`d log cosh(t)` constant, by

```math
Z(B)=\log E_\mu e^{t\langle B,Q\rangle}.          \tag{FF.1}
```

Central symmetry gives the useful pointwise normalization

```math
Z(B)=\log E_\mu\cosh(t\langle B,Q\rangle)\ge0.   \tag{FF.2}
```

For `lambda>0`, put

```math
{dq\over dU}={e^{-\lambda Z}\over z},
\qquad z=E_Ue^{-\lambda Z}.                       \tag{FF.3}
```

For row `i`, let `Q_i` denote the corresponding sign row and put

```math
Z_i(b)=\log E_\mu e^{t\langle b,Q_i\rangle},
\qquad
{dr_i\over dU_i}={e^{-\lambda Z_i}\over z_i},
\qquad z_i=E_{U_i}e^{-\lambda Z_i}.               \tag{FF.4}
```

The product `r=otimes_i r_i` is exactly the canonical erased-row inverse
escort.  Write

```math
D_i=D(r_i\Vert U_i),
\qquad
\mathcal J=D(r\Vert q),
\qquad
\mathcal I^{\leftarrow}
=\inf_{p=\otimes_i p_i}D(p\Vert q),
\qquad A=E_UZ.                                    \tag{FF.5}
```

Thus `I^leftarrow<=J`, and `J-I^leftarrow` is the coherent factor-retuning
advantage of the best row product over the canonical one.

## 2. An exact directional inequality

### Theorem FF.1 (fair roof plus row cost forces coherent retuning)

For every centrally symmetric finite rank-one channel above,

```math
\boxed{
\mathcal I^{\leftarrow}\le G\le \lambda A,
\qquad
\mathcal J-\mathcal I^{\leftarrow}
\ge \sum_{i=1}^mD_i-\lambda A-G
\ge \sum_{i=1}^mD_i-2\lambda A.}                 \tag{FF.6}
```

where

```math
G:=D(U\Vert q).
```

In particular, there is the exhaustive pressure/retuning alternative

```math
\boxed{
\max\{G,\mathcal J-\mathcal I^{\leftarrow}\}
\ge {1\over2}\left[\sum_iD_i-\lambda A\right]_+.} \tag{FF.6a}
```

More precisely, the identities and one-sided bounds behind (FF.6) are

```math
\begin{aligned}
G=D(U\Vert q)&=\lambda A+\log z\le\lambda A,\\
\mathcal J
&=\sum_iD_i+\lambda E_rZ+\log z
 \ge\sum_iD_i-\lambda A.                         \tag{FF.7}
\end{aligned}
```

*Proof.*  Since `Z>=0`, one has `z<=1`.  The first line of (FF.7) follows by
substituting (FF.3), and `U` is an admissible row product, so
`I^leftarrow<=G<=lambda A`.

From (FF.4),

```math
D_i=-\lambda E_{r_i}Z_i-\log z_i.                \tag{FF.8}
```

Substitution of the two densities into `D(r||q)` gives the exact second
line of (FF.7) before the inequality:

```math
\mathcal J
=\sum_iD_i+\lambda E_rZ+\log z.                 \tag{FF.9}
```

Now `E_rZ>=0` by (FF.2), while Jensen gives
`log z>=-lambda E_UZ=-lambda A`.  This proves the lower bound on `J`.
Subtracting the sharper upper bound `I^leftarrow<=G` proves the first lower
bound in (FF.6); `G<=lambda A` gives the second.  Put
`H=sum_iD_i-lambda A`.  If `G>=H/2`, the first branch of (FF.6a) holds; if
`G<H/2`, then (FF.6) gives
`J-I^leftarrow>=H-G>H/2`.  This proves (FF.6a).  `square`

The quantity `G` is exactly a pressure gain.  If `L=c+Z`, and

```math
V_\lambda=-{1\over\lambda}\log E_Ue^{-\lambda L},
```

then

```math
\boxed{G=\lambda(E_UL-V_\lambda).}                \tag{FF.6b}
```

Consequently, for a declared target `T` and tolerance `E>=0`, the pressure
branch is target-reaching whenever

```math
G\ge\lambda(E_UL-T-E),
```

because then `V_lambda<=T+E`.  Inequality (FF.6a) by itself gives an
extensive pressure gain, not the separate comparison of its constant with
an externally declared target.

The theorem is directional rather than merely a restatement of
`I^leftarrow<=J`: it replaces the complete product variational problem by
the fair pressure mean `A` and the sum of `m` one-row KL divergences, and
forces their excess into either pressure gain or coherent factor retuning.

## 3. Subexponential integrable frame roofs

The fair quantity `A` has a useful low-information upper certificate.  Let

```math
\mathcal C=\{C^1,\ldots,C^K\}\subset\mathbb R^{m\times n}
```

be centrally symmetric, with `||C^k||_F<=sqrt(d)`.  In the intended frame
application one also has `rank(C^k)<=R_N`, although rank is not needed for
the analytic inequality below.  Suppose the family is an
**integrable upper frame roof** in the following sense:

```math
\boxed{
Z(B)\le \varepsilon_NN
       +t\max_{1\le k\le K}\langle B,C^k\rangle
\quad\hbox{for every bridge }B.}                 \tag{FF.10}
```

This hypothesis concerns the scalar pressure itself.  It is stronger than
saying that, after computing the full posterior at `B`, its leading singular
vectors happen to lie near one of `K` frames.

### Corollary FF.2 (subexponential roof conversion)

Under (FF.10),

```math
\boxed{
A\le\varepsilon_NN+t\sqrt{2d\log K}.}            \tag{FF.11}
```

Consequently,

```math
\boxed{
\begin{aligned}
\mathcal I^{\leftarrow}
&\le\lambda\{\varepsilon_NN+t\sqrt{2d\log K}\},\\
\mathcal J-\mathcal I^{\leftarrow}
&\ge\sum_iD_i
-2\lambda\{\varepsilon_NN+t\sqrt{2d\log K}\}.
\end{aligned}}                                    \tag{FF.12}
```

In particular, if `mn=Theta(N^2)`, `t=beta/sqrt(N)`,

```math
\varepsilon_N=o(1),\qquad \log K=o(N),
\qquad \sum_iD_i\ge\eta N,                       \tag{FF.13}
```

then

```math
\boxed{
\mathcal I^{\leftarrow}=o(N),
\qquad
\mathcal J-\mathcal I^{\leftarrow}
\ge\eta N-o(N).}                                  \tag{FF.14}
```

*Proof.*  For fixed `C`, the Rademacher variable
`t<B,C>` is subgaussian with variance proxy `t^2||C||_F^2<=t^2d`.
The log-sum-exp bound, optimized over its exponential parameter, gives

```math
E_U\max_{k\le K}t\langle B,C^k\rangle
\le t\sqrt{2d\log K}.                            \tag{FF.15}
```

Averaging (FF.10) proves (FF.11), and Theorem FF.1 gives (FF.12).  At
physical scale,

```math
t\sqrt{d\log K}=O_\beta(\sqrt{N\log K})=o(N),    \tag{FF.16}
```

which proves (FF.14).  `square`

Thus the weakest synchronization hypothesis isolated by this argument has
two logically separate parts:

1. **integrability and response-image control:** a subexponential family of
   fixed frames upper-carries the scalar pressure, with fair mean error
   `o(N)`;
2. **direction:** the actual one-row inverse escorts have total KL cost
   `Omega(N)`.

The macroscopic posterior singular value supplies neither part by itself.
The second quantity is substantially cheaper than the complete bridge law,
but it remains a separate optimizer-specific observable.

Also, `log K=o(N)` controls the cardinality of the frame response image.  A
generic real `m by n` frame can still require `Theta(N^2)` coordinates, so
the condition does not by itself provide a mergeable low-bit description.

Also, `log K=o(N)` bounds the number of semantic response states, not the
cost of discovering or storing an arbitrary unstructured dictionary of
`K` matrices.  An operational quotient still needs a child-generated rule
for the low-rank frames.  Corollary FF.2 addresses directional conversion
after such a rule is available; it does not manufacture the rule.

## 4. Exact fixed two-sided frame

There is one case in which frame synchronization itself forces the complete
latent structure and hence permits a sharp conversion.

Let

```math
M(B)=E[Q\mid B]
= {E_\mu Q\prod_e(1+\rho B_eQ_e)
   \over E_\mu\prod_e(1+\rho B_eQ_e)},
\qquad \rho=\tanh t.                              \tag{FF.17}
```

### Proposition FF.3 (an exact one-dimensional posterior frame is antipodal)

Fix a rank-one sign matrix `C`.  If

```math
M(B)\in\operatorname {span}\{C\}
\quad\hbox{for every }B\in\{-1,1\}^{m\times n},  \tag{FF.18}
```

then

```math
\boxed{\operatorname {supp}\mu\subseteq\{C,-C\}.} \tag{FF.19}
```

Under central symmetry, `mu` is the fair law on `{C,-C}`.

*Proof.*  Let `Pi_perp` be orthogonal projection onto the complement of
`span{C}`.  The denominator in (FF.17) is positive.  Hence (FF.18) says that
the vector-valued multilinear polynomial

```math
F(B)=E_\mu\left[\Pi_\perp Q
       \prod_e(1+\rho B_eQ_e)\right]              \tag{FF.20}
```

vanishes on the entire bridge cube.  Its Walsh coefficient at a set `S` of
bridge coordinates is

```math
\widehat F(S)=\rho^{|S|}
 E_\mu\left[\Pi_\perp Q\prod_{e\in S}Q_e\right]. \tag{FF.21}
```

Since `rho>0`, every Fourier coefficient of the vector-valued signed
measure

```math
g(q)=\mu(q)\Pi_\perp q
```

on the sign cube is zero.  Walsh inversion therefore gives `g(q)=0` for
every `q`.  Thus every support atom lies in `span{C}`.  Both `q` and `C` have
all entries in `{+-1}`, so the only possible scalar multiples are `q=C`
and `q=-C`.  Central symmetry gives equal weights.  `square`

Gauge `C` to the all-one matrix and let `S(B)=sum_(ij)B_ij`.  In this case

```math
Z(B)=\log\cosh(tS(B)),
\qquad
M(B)=C\tanh(tS(B)).                               \tag{FF.22}
```

Thus the posterior frame is genuinely fixed on both sides and its leading
singular value is

```math
\sigma_1(M(B))=\sqrt d\,|\tanh(tS(B))|.           \tag{FF.23}
```

For a row sum `S_n=sum_(j=1)^nB_j`, define

```math
z_{n,N}=E\cosh(tS_n)^{-\lambda},
```

and

```math
D_{n,N}
=-\lambda {E[\cosh(tS_n)^{-\lambda}
                  \log\cosh(tS_n)]\over z_{n,N}}
 -\log z_{n,N}.                                   \tag{FF.24}
```

This is exactly every `D_i` in Theorem FF.1.  Since
`E_U log cosh(tS)<=t sqrt(d)`, (FF.6) gives the finite-order estimates

```math
\boxed{
\mathcal I^{\leftarrow}\le\lambda t\sqrt d,
\qquad
\mathcal J-\mathcal I^{\leftarrow}
\ge mD_{n,N}-2\lambda t\sqrt d.}                 \tag{FF.25}
```

If `n/N->alpha in (0,1)` and `t=beta/sqrt(N)`, then the central limit
theorem and bounded uniform integrability give

```math
D_{n,N}\longrightarrow d_{\alpha,\beta,\lambda}
:=-\lambda {E[\operatorname {sech}^{\lambda}(bG)
                  \log\cosh(bG)]
                 \over E\operatorname {sech}^{\lambda}(bG)}
  -\log E\operatorname {sech}^{\lambda}(bG)>0,   \tag{FF.26}
```

where `G` is standard Gaussian and `b=beta sqrt(alpha)`.  Positivity is
strict because (FF.26) is the KL divergence of the nonconstant density
proportional to `sech^lambda(bG)` relative to Gaussian measure.  Therefore

```math
\boxed{
\mathcal I^{\leftarrow}=O(\sqrt N),
\qquad
\mathcal J-\mathcal I^{\leftarrow}
\ge(1-\alpha)d_{\alpha,\beta,\lambda}N-o(N).}    \tag{FF.27}
```

This exact endpoint confirms the interpretation of (FF.12): a globally
synchronized phase has only sublinear reverse-product information, but
allowing every row to retune that common phase independently costs a positive
amount per row.

For completeness, along `q_{-a}` with any fixed `a>0`, a lattice local CLT
shows that `tS(B)` converges to the law on the real line with density
proportional to `sech^a(x)`.  Hence (FF.23) is macroscopic on the complete
negative path.  The corresponding cavity-overlap density converges to

```math
{\int\operatorname {sech}^a(x)\tanh^2(x)dx
 \over\int\operatorname {sech}^a(x)dx}
={1\over1+a},                                     \tag{FF.28}
```

using
`int sech^(a+2)=a/(a+1) int sech^a`.  Thus the integrated negative-path
overlap tends `log(1+lambda)/lambda>0` simultaneously with (FF.27).
The same limit holds for deleted cavities because Theorem PN.1 gives the
uniform comparison `|M_e-r_e|<=2rho=o(1)`.

## 5. What an approximate posterior coreset does and does not give

Theorem FF.1 is stable under a **fair one-sided pressure approximation**.
For example, if a coreset pressure `Z_tilde` has an upper roof of the form
(FF.10) and

```math
E_U(Z-\widetilde Z)_+\le\delta_NN,                \tag{FF.29}
```

then one simply replaces `epsilon_N` by `epsilon_N+delta_N` in
(FF.11)--(FF.14).

By contrast, each of the following is insufficient on its own:

- `E_q|Z-Z_tilde|=o(N)`;
- an `L^2(q)` approximation to the cavity matrix;
- a pointwise best-rank approximation whose singular frame is selected
  after evaluating the complete posterior;
- a finite list of frames with no curl-free pressure primitive.

The first failure is a genuine measure-change issue.  Since
`dq/dU` is proportional to `e^(-lambda Z)`, a set on which `Z` is large may
have exponentially negligible `q`-mass and macroscopic `U`-mass.  Therefore
escort error does not bound `A=E_UZ`.  The cavity statements additionally do
not integrate to a scalar roof unless a curl-free compatibility condition is
proved.

Nor does a coreset imply the row-cost hypothesis in (FF.13).  The fixed-
projective row-factor channel gives the sharp warning: all posterior means
are rank one and share one fixed right frame, with a macroscopic leading
singular value, while the negative escort is exactly a row product and

```math
\mathcal I^{\leftarrow}=\mathcal J=0.             \tag{FF.30}
```

Its missing left frame consists of independent row signs.  Thus a one-sided
fixed spectral frame, or pointwise rank one, is not the synchronization
needed by Corollary FF.2.  This example has an extensive `G` and therefore
falls into the pressure branch of (FF.6a); it is not a counterexample to the
exhaustive pressure/retuning alternative.

## 6. Exact boundary for the current SML

The directional part of posterior-frame synchronization can now be stated
without the full product oracle:

> Construct from the actual children an integrable frame roof satisfying
> (FF.10) with `log K=o(N)` and `epsilon_N=o(1)`, and prove that the explicit
> one-row cost `sum_iD(r_i||U_i)` is linear; or prove that one of these two
> requirements fails for actual optimizing children.

If both hold, (FF.14) gives the coherent-retuning alternative with an exact
linear constant.  This does not by itself give target reach: the latter still
requires comparing the actual soft value with the declared parent target.
The theorem therefore supplies a rigorous directional conversion for a
genuinely synchronized frame, while making explicit that neither adaptive
low rank nor positive overlap proves the synchronization or the row cost.
