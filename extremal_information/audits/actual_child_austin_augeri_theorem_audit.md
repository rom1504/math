# Actual-child bridge potential versus Austin--Augeri mean-field theorems

Status: **primary-source, theorem-level audit**.  This note maps the actual
optimized-child potential to the row- and bit-coordinate hypotheses in the
published theorems.  It is not a general literature survey.

**Conclusion.**  Neither theorem yields an `o(N)` row-product target excess.
Augeri's Bernoulli theorem gives only `O(N)`, exactly the already-known
scale.  Austin can use whole rows as coordinates, but its published
partition-function estimate would need an essentially singleton cover of
the full row-cavity gradient image at `o(1)` uniform accuracy.  Rank-one
latent support does not supply that cover.  The rank-one support relaxation
has Rademacher width `Theta(N)`, so the `O(N)` Augeri ceiling is sharp for
that method unless one proves a genuinely actual-child-specific collapse of
the gradient image.

Primary sources:

- Tim Austin, [*The structure of low-complexity Gibbs measures on product
  spaces*](https://arxiv.org/abs/1810.07278), Theorem A, Corollary A', and
  Proposition 5.1.
- Fanny Augeri, [*A transportation approach to the mean-field
  approximation*](https://arxiv.org/abs/1903.08021), Theorem 1.1 and its
  harmonic-extension remark.
- Fanny Augeri, [*Nonlinear large deviation bounds with applications to
  traces of Wigner matrices and cycles counts in Erdos--Renyi
  graphs*](https://arxiv.org/abs/1810.01558), Theorem 1.1 and Corollary 1.2.

## 1. Exact target and actual-child geometry

Fix actual pressure-minimizing children, one orientation, a comparable split
`m+n=N`, and

```math
d=mn,\qquad t={\beta\over\sqrt N},\qquad
f(B)=-\lambda L(B).
```

With `U` the fair bridge law and `q=e^fU/E_Ue^f`, Gibbs' identity gives

```math
\boxed{
\mathcal I_{\rm row}^{\leftarrow}
=\inf_{p=\otimes_i p_i}D(p\Vert q)
=\log E_Ue^f-
 \sup_{p=\otimes_i p_i}\{E_pf-D(p\Vert U)\}.}       \tag{MA.1}
```

Thus the desired statement is exactly an `o(N)` row-product mean-field
error.  Products of all `d` individual bits form a smaller variational class,
so a bit-product theorem is legal but potentially wasteful:

```math
\mathcal I_{\rm row}^{\leftarrow}
\le \mathcal I_{\rm bit}^{\leftarrow}.              \tag{MA.2}
```

Up to an irrelevant additive constant, the natural continuous extension is

```math
\widetilde f(b)
=-\lambda\log\sum_z w_z e^{t\langle b,Q_z\rangle},
\qquad Q_z=\tau xy^T.                                \tag{MA.3}
```

Here the weights `w_z` are the **actual child Gibbs weights**.  At finite
temperature the support consists of all rank-one sign words.  Directly,

```math
\nabla\widetilde f(b)=-\lambda tE_bQ,
\qquad
\nabla^2\widetilde f(b)
=-\lambda t^2\operatorname{Cov}_b(\operatorname{vec}Q).       \tag{MA.4}
```

In particular the continuous gradient lies in
`lambda t conv{+-xy^T}`.  This mapping uses the optimized children rather
than a conference surrogate; no property of their weights has yet been
discarded.

## 2. Augeri on the bit cube: a sharp `O(N)` method ceiling

Augeri's 2019 theorem says, for the harmonic extension of a function on a
Bernoulli cube,

```math
\mathcal I_{\rm bit}^{\leftarrow}
\le \kappa\, b(\mathcal G_d),
\qquad
b(V)=E_\varepsilon\sup_{v\in V}\langle\varepsilon,v\rangle,  \tag{MA.5}
```

where `G_d` is the image of the **discrete half-difference gradient**.  The
distinction from (MA.4) costs at leading scale here.  If `partial_e^d f` is
the half secant at a cube vertex, then the diagonal Hessian bound from
(MA.4) gives

```math
|\partial_e^df(B)-\partial_e\widetilde f(B)|
\le\lambda t^2.                                      \tag{MA.6}
```

Consequently, for a fair `m`-by-`n` Rademacher matrix `E`,

```math
\begin{aligned}
b(\mathcal G_d)
&\le \lambda tE\max_{x,y}|x^TEy|+\lambda t^2mn\\
&\le \lambda t\sqrt{mn}\,E\|E\|_{op}+\lambda t^2mn\\
&\le C\lambda t\sqrt{mn}(\sqrt m+\sqrt n)+\lambda t^2mn\\
&=O_{\lambda,\beta}(N).                              \tag{MA.7}
\end{aligned}
```

Combining (MA.2), (MA.5), and (MA.7) proves only

```math
\boxed{\mathcal I_{\rm row}^{\leftarrow}=O_{\lambda,\beta}(N).} \tag{MA.8}
```

This does not improve the elementary actual-law estimate
`I_row^leftarrow<=lambda^2t^2mn/2`.

The support relaxation in (MA.7) is sharp at the same scale.  For balanced
splits,

```math
E\max_{x,y}x^TEy
\ge E\sum_{j=1}^n\left|\sum_{i=1}^mE_{ij}\right|
\ge c n\sqrt m,
```

while the operator-norm upper bound has order
`sqrt(mn)(sqrt(m)+sqrt(n))`.  Hence

```math
b\big(\lambda t\{xy^T:x,y\}\big)=\Theta_{\lambda,\beta}(N). \tag{MA.9}
```

Equation (MA.9) is a **method-specific ceiling**, not a claim that the
actual gradient image fills its rank-one envelope.  It proves that rank-one
support alone cannot turn Augeri's theorem into `o(N)`.  The precise new
input that would do so is

```math
\boxed{b(\mathcal G_d)=o(N)
\quad\hbox{for the actual optimized-child gradient image}.}   \tag{MA.10}
```

This is an actual-child statement and is strictly narrower than controlling
the full bridge energy landscape, but it is not supplied by the imported
theorem.

Augeri's 2018 covering-number theorem is weaker at this scale.  Its Gaussian
width corollary has error `C d^(1/3)g(V)^(2/3)`.  Even granting
`g(V)=O(N)` and `d=Theta(N^2)`, it gives `O(N^(4/3))`, not `o(N)`.

Finally, Augeri's Bernoulli tilts are bit-product laws.  Treating a row as a
vector in `R^n` still produces only products within each row.  Representing
**every** row law as an affine tilt requires a one-hot alphabet embedding of
dimension `2^n-1` per row, at which point the mean-width/covering bounds are
exponential.  Her theorem therefore has no hidden row-coordinate shortcut.

## 3. Austin on whole rows: the hypothesis is the missing structure

Austin does allow the coordinates to be the `m` row alphabets
`K_i={+-1}^n`.  With a reference row `R_*`, his discrete row gradient is the
additively separable function

```math
\nabla_{\rm row}f(B,R)=\sum_i\partial_if(B,R_i).
```

For (MA.3), the exact component is

```math
\partial_if(B,R)
=-\lambda\log
 {E_{\pi_{i,B}}e^{t\langle R,Q_i\rangle}
  \over
  E_{\pi_{i,B}}e^{t\langle R_*,Q_i\rangle}},          \tag{MA.11}
```

where `pi_(i,B)` is the actual latent posterior after exposing every bridge
row except `i`.  Thus covering the row-gradient image means covering the
family of full row-cavity posterior response functions.  The fact that each
`Q` is rank one does not bound the number of such posteriors.

Austin's hypothesis is an explicit cover: for some `epsilon,delta`,

```math
\log\operatorname{cov}_{\delta m}
 \big(\operatorname{img}\nabla_{\rm row}f,\|\cdot\|_\infty\big)
\le\epsilon m.                                       \tag{MA.12}
```

Theorem A then gives a conditioning mixture whose average **forward** KL to
products is at most `(epsilon+delta)m`.  That statement alone is not a
single row product and does not bound (MA.1).  Austin's Proposition 5.1 does
bound the row-product variational gap:

```math
\mathcal I_{\rm row}^{\leftarrow}
\le(\epsilon+\delta)m
 +\sqrt{(\epsilon+\delta)/2}\;L_{\rm Ham}.            \tag{MA.13}
```

Using normalized within-row Hamming metrics, the outer product metric is
global normalized bridge Hamming distance.  The exact one-bit oscillation
gives the universal actual-child scale

```math
L_{\rm Ham}\le2\lambda t mn=\Theta_{\lambda,\beta}(N^{3/2}). \tag{MA.14}
```

Therefore the published estimate (MA.13), used with only the universally
available scale (MA.14), yields `o(N)` only if

```math
\epsilon_N+\delta_N=o(1/N).                           \tag{MA.15}
```

But then (MA.12) asks for a cover at uniform gradient error
`delta_Nm=o(1)` with logarithmic cardinality
`epsilon_Nm=o(1)`: eventually an essentially one-set cover.  Neither
rank-one latent support nor the already-proved bounded conditional row
Renyi-two complexity implies this collapse.  The latter controls each
conditional row density; Austin's hypothesis controls the entire family of
row-cavity response functions (MA.11).

Thus Austin supplies a valid conditional criterion, but verifying it at the
needed scale is at least the following new actual-child lemma:

```math
\boxed{
\operatorname{img}\nabla_{\rm row}f
\text{ has an }e^{o(1)}\text{-cover at }o(1)
\text{ uniform accuracy},}                           \tag{MA.16}
```

unless one first proves a substantially smaller actual Lipschitz scale or a
new replacement for the transport step.  The primary source contains no
such theorem.

## 4. Audited judgment

1. **No imported `o(N)` theorem.**  Augeri gives `O(N)` in bit coordinates;
   Austin's row theorem requires the unresolved cavity-gradient collapse and
   its structural mixture cannot be converted silently into reverse KL to a
   single product.
2. **Sharpest source-level ceiling.**  The rank-one atom envelope has width
   `Theta(N)` at physical scaling, so all support-only uses of Augeri stop at
   linear error.  Austin's published transport conversion likewise needs
   the near-singleton cover (MA.15)--(MA.16) at the universal row-change
   scale.
3. **Smallest theorem exposed by this audit.**  Prove either the actual
   discrete-gradient width collapse (MA.10), or a new row-cavity theorem
   bounding (MA.1) directly without the `N^(3/2)` transport conversion.
   Both explicitly concern optimized-child weights; neither follows from a
   conference example, generic row law, or raw spectral calculation.
