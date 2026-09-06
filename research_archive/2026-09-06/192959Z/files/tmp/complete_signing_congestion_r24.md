# Wave 24 Route 3: complete-signing competitor Fourier law and plaquette curvature

This note keeps the shared edge variables which were absent from the Wave 23
abstract star wall.  Combining *all* edge-perturbed partition sums gives an
exact Walsh/MacWilliams formula supported on even Eulerian subgraphs.  It
reduces the desired cavity estimate to stability of a highly cancelling
signed-cycle polynomial under vertex deletion.  Exact minimality supplies a
pointwise lower floor for all translates of that polynomial, but no upper
deletion ratio.

There is one genuinely favorable universal consequence of complete-signing
structure: every two-coordinate Gibbs square has fixed nonzero log-curvature.
A sharp four-variable optimization forces some coordinate affinity below a
fixed constant.  The resulting cavity reward is `O_beta(1)`, however, and
cannot approach the required `Theta(sqrt(r))` scale.

All identities and finite audits below are verified by
`complete_signing_congestion_r24.py`.  No asymptotic counterexample to
(10.617) is claimed.

## 1. Full-edge Walsh transform of every competitor partition sum

Let `B` be an order-`r` complete signing.  Use one representative of each
projective spin together with both orientations, so the state space
`Omega_r` has `2^r` elements.  Write

```math
s_{ij}(\omega)=\sigma b_{ij}x_ix_j,
\qquad
e_B(\omega)=2\sum_{i<j}s_{ij}(\omega),
\qquad
\Delta_\omega=q-e_B(\omega),
```

where at first `q` is any scalar baseline.  Put `E=E(K_r)`,
`N=|E|`, and for every perturbation `T\subseteq E` define

```math
D_T
=\sum_{\omega\in\Omega_r}
\exp\left\{-\beta\left[
\Delta_\omega+4\sum_{e\in T}s_e(\omega)
\right]\right\}.
\tag{R24.1}
```

For `A\subseteq E`, use the normalized Walsh transform

```math
\widehat D(A)
=2^{-N}\sum_{T\subseteq E}(-1)^{|A\cap T|}D_T.
```

Let `\mathcal E_{\rm even}` be the family of edge sets which have even degree at
every vertex and an even total number of edges.  For
`A\in\mathcal E_{\rm even}`
write `b_A=\prod_{e\in A}b_e`.  Then

```math
\boxed{
\widehat D(A)
=
\begin{cases}
2^r e^{-\beta q}
\cosh(2\beta)^{N-|A|}
\sinh(2\beta)^{|A|}b_A,
&A\in\mathcal E_{\rm even},\\
0,&A\notin\mathcal E_{\rm even}.
\end{cases}
}
\tag{R24.2}
```

### Proof and normalization

For a fixed state and edge, averaging the perturbation bit `t_e` gives

```math
\frac{1+e^{-4\beta s_e}}2
=e^{-2\beta s_e}\cosh(2\beta)
```

when the Walsh character omits `e`, and

```math
\frac{1-e^{-4\beta s_e}}2
=s_e e^{-2\beta s_e}\sinh(2\beta)
```

when it contains `e`.  Multiplying over all edges cancels the state energy:

```math
e^{-\beta\Delta_\omega}e^{-2\beta\sum_es_e(\omega)}
=e^{-\beta q}.
```

It remains to evaluate

```math
\sum_{\omega\in\Omega_r}\prod_{e\in A}s_e(\omega)
=b_A\sum_{\sigma,x}
\sigma^{|A|}\prod_vx_v^{\deg_A(v)}.
```

The orientation sum vanishes unless `|A|` is even.  The projective spin sum
vanishes unless every vertex degree is even; fixing one projective spin does
not change this condition because the degree parities sum to zero.  If both
conditions hold, every one of the `2^r` states contributes `b_A`.  This proves
(R24.2), including its factor `2^r`.

With `\rho=\tanh(2\beta)` and

```math
mathcal P_B(rho)
=\sum_{A\in\mathcal E_{\rm even}}b_A\rho^{|A|},
```

Walsh inversion is

```math
\boxed{
D_T
=2^r e^{-\beta q}\cosh(2\beta)^N
\sum_{A\in\mathcal E_{\rm even}}
b_A(-1)^{|A\cap T|}\rho^{|A|}.
}
\tag{R24.3}
```

The support `\mathcal E_{\rm even}` is exactly the dual of the augmented cut code.
Thus (R24.3) is the complete shared-edge version of all competitor
inequalities, not a product of unrelated star marginals.

The Fourier identity uses **every complete signing**, with no minimality
hypothesis and even with an arbitrary baseline `q`.  Exact global minimality
enters only afterward: if `B` is an exact minimizer and `q=q_r=Q(B)`, then
`Q(B^T)\ge q` for every `T`, so

```math
\boxed{D_T\ge1\qquad(T\subseteq E).}
\tag{R24.4}
```

It also gives `Delta_omega>=0` at the unperturbed signing.  Neither fact is
used in deriving (R24.2)--(R24.3).

Two useful consistency checks are

```math
2^{-N}\sum_TD_T
=2^r e^{-\beta q}\cosh(2\beta)^N
\tag{R24.5}
```

and the corresponding Parseval identity, whose right side is independent of
the cycle signs after squaring the coefficients in (R24.2).

## 2. The exact signed-cycle deletion ratio

Let `H=delta(i)` be the coordinate star, of size `m=r-1`.  Averaging (R24.3)
over all `T\subseteq H` kills precisely the Eulerian edge sets which use an
edge incident to `i`.  What remains is the same even-Eulerian polynomial of
the principal child:

```math
2^{-m}\sum_{T\subseteq\delta(i)}D_T
=2^r e^{-\beta q}\cosh(2\beta)^N
\mathcal P_{B[-i]}(\rho).
\tag{R24.6}
```

At `T=emptyset`, (R24.3) gives

```math
D_\beta(B)
=2^r e^{-\beta q}\cosh(2\beta)^N\mathcal P_B(\rho).
```

Comparing (R24.6) with (10.735) proves the exact complete-signing ratio

```math
\boxed{
\operatorname{BC}_i
=\frac{\mathcal P_{B[-i]}(\rho)}
       {\cosh(2\beta)^{r-1}\mathcal P_B(\rho)}.
}
\tag{R24.7}
```

Although the polynomials are written as signed sums, their displayed values
are positive because they are partition functions after the positive common
factor is removed.

For the threshold

```math
T_r=\alpha[r^{3/2}-(r-1)^{3/2}],
```

the precise signed-cycle statement sufficient for (10.617) at one step is

```math
\boxed{
\frac{\mathcal P_{B[-i]}(\rho)}{\mathcal P_B(\rho)}
\le
\cosh(2\beta)^{r-1}e^{-\beta(T_r-K)}
}
\tag{R24.8}
```

for some selectable coordinate `i`.  Thus the missing upper competitor bound
is exactly stability of the frustrated even-cycle cancellation under deleting
one vertex.

The pointwise floor (R24.4) does not give (R24.8).  It says that every edge-
character twist of the full polynomial in (R24.3) is at least a fixed positive
number.  Star averaging supplies a lower bound on the vertex-avoiding
polynomial.  A total-budget or Parseval upper bound loses a factor exponential
in the `Theta(r^2)` complementary edges.  Taking absolute values is worse:

```math
|\mathcal P_B(\rho)|
\le\sum_{A\in\mathcal E_{\rm even}}\rho^{|A|}
```

replaces the frustrated signing by the ferromagnetic one and destroys the
cycle cancellation that determines the cavity ratio.

The finite audits quantify this loss.  The ratios of the signed polynomial to
its coefficientwise unsigned majorant are

| matrix | `beta=0.5` | `beta=1` |
|---|---:|---:|
| `A_8` | `1.71e-7` | `1.94e-15` |
| `A_9` | `9.87e-10` | `3.57e-20` |

These finite numbers are not asymptotic lower bounds.  They are a precise
no-go for coefficientwise/absolute-cycle majorization: it discards many
orders of magnitude before any deletion comparison is made.

## 3. Complete-signing plaquette curvature gives a valid upper affinity

There is a different consequence which genuinely has the desired upper
direction.  Work on the fully lifted spin cube; the duplicate projective copy
does not change ratios.  Fix two coordinates `i ne j`, condition on the
orientation and all other spins, and write the square-root Gibbs weights on
the resulting two-dimensional face as positive numbers `a,b,c,d`, with
`a,d` on one diagonal.  Because the edge `ij` has coupling sign `+1` or `-1`,

```math
\left|\log\frac{ad}{bc}\right|=4\beta.
\tag{R24.9}
```

After exchanging the diagonals, put `ad=lambda bc`, where
`lambda=e^{4beta}>=1`.  The local contribution to the sum of the two
conditional affinities is

```math
L(a,b,c,d)
=\frac{2(ab+cd+ac+bd)}{a^2+b^2+c^2+d^2}
=\frac{2(a+d)(b+c)}{a^2+b^2+c^2+d^2}.
```

The sharp optimization is

```math
\boxed{
L(a,b,c,d)\le C_\beta,
\qquad
C_\beta=
\begin{cases}
2\operatorname{sech}(2\beta),&e^{4\beta}\le3,\\
\sqrt{\dfrac{2}{1-e^{-4\beta}}},&e^{4\beta}\ge3.
\end{cases}
}
\tag{R24.10}
```

### Sharp four-variable optimization

Set

```math
x=a+d,
\qquad y=b+c,
\qquad p=bc,
\qquad ad=\lambda p.
```

Then

```math
a^2+b^2+c^2+d^2
=x^2+y^2-2(\lambda+1)p,
```

while feasibility gives

```math
p\le\min\{y^2/4,x^2/(4\lambda)\}.
```

For fixed `x,y`, the ratio increases with `p`, so equality holds in the
smaller constraint.  With `t=x/y`, the branch `t>=sqrt(lambda)` is decreasing
and is maximized at the boundary, with value

```math
\frac{4\sqrt\lambda}{\lambda+1}
=2\operatorname{sech}(2\beta).
```

On `t<=sqrt(lambda)`, the ratio is

```math
\frac{2t}{1+[(\lambda-1)/(2\lambda)]t^2}.
```

Its unconstrained maximizer is
`t=\sqrt{2\lambda/(\lambda-1)}`.  This is feasible exactly when
`\lambda\ge3`, and its value is `\sqrt{2\lambda/(\lambda-1)}`.  Otherwise the
boundary value applies.  Both branches are attainable, proving sharpness of
(R24.10) as a local square inequality.

Partitioning the full Gibbs cube into these conditional squares and weighting
(R24.10) by their four-vertex masses gives, for every pair,

```math
\boxed{
\operatorname{BC}_i+\operatorname{BC}_j\le C_\beta.
}
\tag{R24.11}
```

Summing over all coordinate pairs yields

```math
\frac1r\sum_i\operatorname{BC}_i\le\frac{C_\beta}{2}.
```

Therefore every complete signing, with no minimality assumption, has

```math
\boxed{
\max_i\kappa_{\beta,i}
\ge
\frac1\beta\log\frac{2}{C_\beta}.
}
\tag{R24.12}
```

When `r` is even, the odd-field parity bound
`beta^{-1}log cosh(2beta)` may be stronger; (R24.12) is useful in particular
when zero fields are parity-allowed.  In all cases its right side depends only
on the fixed temperature.  Hence

```math
\frac{1}{\sqrt r}\frac1\beta\log\frac2{C_\beta}\longrightarrow0,
```

so plaquette curvature cannot establish the `Theta(sqrt(r))` reward in
(10.617).  Combining the pair inequalities linearly is already sharp at the
level of those inequalities: the feasible scalar profile
`\operatorname{BC}_i=C_\beta/2` for every `i` saturates all of them.  Any amplification must
use higher joint cancellation, not another summation of plaquettes.

## 4. Exact audits

The checker performs three independent tests.

1. On the exact order-five minimizer, it enumerates all `2^10` edge
   perturbations and all 32 oriented-projective states.  At both
   `beta=0.5,1`, the Walsh transform has exactly 32 nonzero coefficients,
   precisely the 32 even Eulerian edge sets predicted by (R24.2).  It also
   verifies (R24.4), (R24.5), and Parseval.
2. On `A_8` and `A_9`, it exhausts every orientation, coordinate pair, outside
   spin assignment, and conditional two-spin square.  Every cross ratio is
   `e^{plus_or_minus 4beta}` and every local value satisfies (R24.10).  The
   resulting global pair bounds pass.
3. It independently computes the signed cycle polynomial and the unsigned
   ferromagnetic majorant from primal partition sums, producing the
   cancellation table above.

Representative global values are:

| matrix | `beta` | max pair-affinity sum | `C_beta` | bound in (R24.12) | actual max reward |
|---|---:|---:|---:|---:|---:|
| `A_8` | 0.5 | 0.646749 | 1.520867 | 0.547734 | 2.257890 |
| `A_9` | 0.5 | 0.861473 | 1.520867 | 0.547734 | 2.665415 |
| `A_8` | 1.0 | 0.208461 | 1.427345 | 0.337331 | 2.261149 |
| `A_9` | 1.0 | 0.669428 | 1.427345 | 0.337331 | 1.760116 |

The complete-signing upper bound is valid but very slack on both exact
minimizers.

## 5. Scoped conclusion

- **Verified joint identity:** the full family of competitor partition sums
  has the exact even-Eulerian Walsh spectrum (R24.2).  This retains every
  shared edge and every block interaction.
- **Exact remaining target:** a square-root cavity reward is equivalent to
  the signed-cycle deletion estimate (R24.8).  Exact minimality presently
  supplies only the pointwise floor (R24.4), not this upper ratio.
- **Verified favorable but insufficient bound:** complete plaquette curvature
  proves (R24.11)--(R24.12), a strictly positive constant reward even when
  zero fields are allowed.  It is `o(sqrt(r))` at every fixed temperature.
- **Closed submechanisms:** coefficientwise absolute cycle bounds lose the
  essential frustration cancellation; total-budget/Parseval bounds lose the
  `Theta(r^2)` complementary-edge scale; and linear aggregation of sharp
  plaquette bounds cannot amplify their constant conclusion.
- **Still open:** a minimizer-specific theorem controlling how the signed
  even-cycle cancellation changes under deletion, or an equivalent genuine
  upper bound on a competitor-star partition sum.  No complete-signing
  asymptotic wall with all rewards `o(sqrt(r))` is constructed, so (10.617)
  remains open.
