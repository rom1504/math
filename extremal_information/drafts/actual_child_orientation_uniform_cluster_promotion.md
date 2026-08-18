# Orientation-uniform physical promotion from row Renyi two

Status: **rigorous theorem with two independent audits**.
This note removes the balanced-orientation premise from the actual-child
cluster promotion theorem. Central symmetry plus bounded row `D_2`, not
pointwise max-density, already gives the dimension-free vector-subgaussian
estimate needed by the off-block quadratic argument. Canonical actual-child
rows have that `D_2` bound in every orientation.

The result does not decide the balanced product phase. It removes a
separate obligation: a target-reaching orientation may be analyzed
directly, without proving that the bias-canceling orientation is itself
target-reaching.

The orientation-uniform Renyi-two constants and the balanced max-density
constants do not uniformly dominate one another.  Thus the result removes a
qualitative target-orientation obligation; it does not transfer the old
balanced Gram cutoff unchanged.

## 1. Renyi-two domination implies vector subgaussianity

Let `U_n` be uniform on the sign cube.

**Lemma OU.1 (symmetric Renyi-two subgaussian lemma).** Let `P` be a
centrally symmetric law on `{-1,1}^n` and suppose

```math
D_2(P\Vert U_n)\le C.
```

Then, for every `v in R^n`,

```math
\boxed{
E_Pe^{\langle v,R\rangle}
\le\exp\{2e^{C/2}\|v\|_2^2\}.}                 \tag{OU.1}
```

Equivalently, `P` is vector-subgaussian with variance proxy

```math
\sigma_C^2=4e^{C/2}.                               \tag{OU.2}
```

*Proof.* Put `f=dP/dU_n` and `a=\|v\|_2^2`. Central symmetry gives

```math
E_Pe^{\langle v,R\rangle}=E_P\cosh\langle v,R\rangle.
```

Because `E_Uf=1` and `E_Uf^2<=e^C`, Cauchy--Schwarz yields

```math
E_P\cosh\langle v,R\rangle
\le1+e^{C/2}
 \{E_U(\cosh\langle v,R\rangle-1)^2\}^{1/2}.  \tag{OU.3}
```

Now

```math
E_U\cosh(2\langle v,R\rangle)
=\prod_j\cosh(2v_j)\le e^{2a},
```

whereas `cosh x>=1+x^2/2` and
`E_U\langle v,R\rangle^2=a`. Hence

```math
E_U(\cosh\langle v,R\rangle-1)^2
\le{e^{2a}-1-2a\over2}\le a^2e^{2a}.            \tag{OU.4}
```

The last inequality is the Taylor remainder
`e^x-1-x<=x^2e^x/2` at `x=2a`. Therefore (OU.3) is at most
`1+e^(C/2)ae^a`. Since `ae^a<=e^(2a)-1` and
`1+K(e^x-1)<=e^(Kx)` for `K>=1`, this is at most
`exp{2e^(C/2)a}`. `square`

Central symmetry is essential for this proof: without it, the linear MGF
term need not vanish.

## 2. Uniform application to canonical actual-child rows

Fix arbitrary finite children, either relative orientation `epsilon`, a
choice of row direction with row width `n`, channel amplitude `u`, and
inverse exponent `lambda`. Let `r_(epsilon,u)` be the iid canonical
row-erased inverse escort from CR.3--CR.4.

The row output likelihood is even under `b mapsto -b`, by the global
right-child spin flip. Thus `r_(epsilon,u)` is centrally symmetric.
The exact bit-oscillation theorem CR.8, uniform in the children and
orientation, gives

```math
D_2(r_{\epsilon,u}\Vert U_n)
\le n\log\{1+\tanh^2(\lambda u)\}
\le\lambda^2u^2n.                                \tag{OU.5}
```

At `u=t=beta/sqrt(N)`, set

```math
C_2=\lambda^2t^2n\le\lambda^2\beta^2,
\qquad
\sigma_*^2=4e^{\lambda^2\beta^2/2}.              \tag{OU.6}
```

Lemma OU.1 shows that every canonical row block in either orientation and
either transpose direction has the common vector-subgaussian proxy
`sigma_*^2`.

The random-row-cut proof of Lemma SP.2 uses only independence between row
blocks and this linear-functional MGF bound. It therefore gives, for

```math
H=\sum_{i<k}R_i^{\mathsf T}M_{ik}R_k,
\qquad
V=\sum_{i<k}\|M_{ik}\|_F^2,
```

the orientation-uniform estimate

```math
\boxed{
\log E e^{\theta H}\le b_*\theta^2V
\quad\text{if}\quad
|\theta|\|M\|_{\rm op}\le a_*,}                \tag{OU.7}
```

where one may take

```math
a_*={1\over8\sqrt2e^{\lambda^2\beta^2/2}},
\qquad
b_*=64e^{\lambda^2\beta^2}.                      \tag{OU.8}
```

Indeed the determinant proof gives
`a_*=1/(2sqrt(2)sigma_*^2)` and `b_*=4sigma_*^4`.

## 3. Orientation-uniform actual-child cluster promotion

Use the exact zero-bridge child prior in a fixed orientation `epsilon`, and
let `K_epsilon` and `mathfrak C_(>=4)^epsilon(t)` be the sector--Gram
coefficient and absolute connected cross-row cluster tail in SP.1--SP.2,
computed in the chosen row direction. Retain SP.1's convergence hypothesis
for every full and row-restricted bridge word.

**Theorem OU.2 (target-orientation-uniform physical promotion).** Let
`m+n=N`, `t=beta/sqrt(N)`, and let the two children be the actual exact
contracted-temperature minimizers. Fix either relative orientation and
either row direction. If

```math
K_\epsilon\le\kappa N^2,
\qquad
\lambda\beta^2\sqrt{2\kappa}\le a_*,           \tag{OU.9}
```

then its exact canonical interaction cumulant obeys

```math
\boxed{
\mathcal J_t^\epsilon
\le b_*\lambda^2t^4K_\epsilon+2\lambda\mathfrak C_{\ge4}^\epsilon(t).}
\tag{OU.10}
```

Consequently a sublinear cluster tail implies
`\mathcal J_t^epsilon=o(N)` in that same orientation. Conversely, under
(OU.9),

```math
\mathcal J_t^\epsilon\ge\eta N
\quad\Longrightarrow\quad
\mathfrak C_{\ge4}^\epsilon(t)
\ge {\eta\over2\lambda}N-O_{\beta,\lambda,\kappa}(1).
\tag{OU.11}
```

*Proof.* Lemma SP.1 is orientationwise and gives
`h_t=t^2H_2+R_t+c_t` with
`osc R_t<=2mathfrak C_(>=4)^epsilon(t)`. For the canonical row product,
OU.5--OU.8 apply without any sector-bias premise. In (OU.7) take
`theta=-lambda t^2`, `V=K_epsilon`, and use
`\|M\|op<=sqrt(2K_epsilon)`. Centering `H_2` costs nothing because
independent centrally symmetric rows give `E H_2=0`. The same
oscillation argument as SP.15 proves (OU.10). Equations
`t^4=beta^4/N^2` and (OU.9) make the quadratic term `O(1)` and prove
(OU.11). `square`

The theorem needs no child minimality beyond identifying the actual child
pair; minimality remains essential elsewhere in the target recurrence.
This extra generality is not a surrogate substitution: the law,
orientation, `K_epsilon`, cluster tail, and `J` in OU.2 are exactly those
induced by the actual optimizing children.

## 4. SML consequence and limitation

Suppose a soft bridge reaches the desired target in some orientation
`epsilon_*`. Apply OU.2 directly to `epsilon_*`. There is no longer a need
to prove that the separate bias-canceling orientation of Theorem 37.32
reaches the same target, nor to pay an orientation-switching loss. Thus
target relevance is removed as an independent premise of the
physical-promotion branch.

What remains is still substantial. OU.2 neither proves (OU.9) nor cluster
tightness, and the all-order absolute tail is not known to have lower
operational complexity than the complete child Gibbs law. It also does not
decide whether a linear canonical interaction is carried by reverse-product
dependence or coherent product retuning. The strictly weaker SML is:

In whichever orientation actually reaches the target, decide from an
operationally lower-information child observable whether `J=o(N)`,
`I^leftarrow=Omega(N)`, or `J-I^leftarrow=Omega(N)`, with a directional
certificate in either extensive branch.

This is a genuine narrowing, not a Level-6 closure: one target-orientation
obligation disappears, while the three-way product phase remains open.

The proof and its scope were independently checked in
[`../audits/actual_child_orientation_uniform_cluster_promotion_adversarial_audit.md`](../audits/actual_child_orientation_uniform_cluster_promotion_adversarial_audit.md).
