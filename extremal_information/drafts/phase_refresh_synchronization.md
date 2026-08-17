# Phase refresh, Boolean pullbacks, and collapse of a scale carrier

Status: rigorous theorem draft.  This is a positive counterpart to the
nonconstant Walsh mantissa profile in Theorem 30.1.  It is an abstract
phase-synchronization theorem for dense quadratic hierarchies; it is not a
theorem about the minimizing values in the motivating signing problem.

## 1. Operator-certified phase transfer

Let `X` be a compact metric phase space.  At level `r` and phase `x in X`,
let `(Omega_(r,x),mu_(r,x))` be a finite probability space and let
`T_(r,x)` be a self-adjoint operator on its real `L^2` space.  Define the
same-spin response

```math
phi_r(x)=\sup_{\|f\|_\infty\le1}
          |\langle f,T_(r,x)f\rangle|.                 \tag{PR.1}
```

A **Boolean contraction** from the system at `(r,x)` to the system at
`(r+1,y)` is a linear map `U` satisfying

```math
\|U\|_(2\to2)\le1,
\qquad \|Uf\|_\infty\le\|f\|_\infty.                 \tag{PR.2}
```

Equal-fibre signed coordinate replication is the main example, and is in
fact an `L^2` isometry.

An operator-certified phase transfer at level `r` consists, for every
`x`, of a Borel probability kernel `Gamma_(r,x)` on pairs `(y,U)`, whose
phase marginal is a Markov kernel `P_r(x,dy)`, such that

```math
\left\|T_(r,x)-
 \int U^*T_(r+1,y)U\,dGamma_(r,x)(y,U)\right\|_(2\to2)
 \le epsilon_r.                                       \tag{PR.3}
```

The pulled-back operators are assumed strongly measurable and Bochner
integrable.  This is a matrix/operator assertion.  In a finite presentation it can be
checked from the branch weights, the coordinate maps, and an operator-norm
certificate; it does not ask for the maximizing Boolean spin or for the
full response landscape.

Write

```math
K_(r,l)=P_rP_(r+1)\cdots P_(r+l-1)                    \tag{PR.4}
```

for the `l`-step phase kernel, with the convention that kernels act on
functions from right to left.

## 2. Phase-refresh theorem

### Theorem PR.1 (Boolean pullback refresh forces phase collapse)

Assume:

1. `phi_r` converges uniformly on `X` to a continuous function `phi`;
2. (PR.3) holds at every level;
3. there are integers `l_r>=1`, numbers `alpha_r in (0,1]`, and one probability
   measure `nu` of full topological support on `X` such that

   ```math
   K_(r,l_r)(x,B)\ge alpha_r nu(B)                    \tag{PR.5}
   ```

   for every phase `x` and Borel set `B`;
4. with

   ```math
   omega_r=\|phi_r-phi\|_\infty,
   \qquad E_r=\sum_{j=r}^{r+l_r-1}epsilon_j,           \tag{PR.6}
   ```

   one has

   ```math
   {omega_r+omega_(r+l_r)+E_r\over alpha_r}\longrightarrow0.
                                                               \tag{PR.7}
   ```

Then `phi` is constant on `X`.

More quantitatively, if `M=max_X phi`, then

```math
M-\int_Xphi\,dnu
\le
\liminf_{r\to\infty}
{omega_r+omega_(r+l_r)+E_r\over alpha_r}.             \tag{PR.8}
```

If `X` is finite and `mu=min_x nu({x})>0`, then

```math
\mathop{osc}_X(phi)
\le {1\over mu}\liminf_{r\to\infty}
{omega_r+omega_(r+l_r)+E_r\over alpha_r}.             \tag{PR.9}
```

Thus a fixed Doeblin refresh fraction `alpha>0` and a bounded or growing
delay whose accumulated transfer defect tends to zero are sufficient.  The
theorem also allows the refresh fraction to vanish, but only when recovery
and transfer errors vanish faster than that fraction.

In particular, if `inf_r alpha_r>0`, then uniform phasewise recovery and
`E_r->0` imply collapse without any rate assumption on the recovery.

#### Proof

Fix `r,x` and `f` with `\|f\|_infty<=1`.  From (PR.2)--(PR.3),

```math
|\langle f,T_(r,x)f\rangle|
\le
\int |\langle Uf,T_(r+1,y)Uf\rangle|\,
          dGamma_(r,x)(y,U)+epsilon_r
\le (P_rphi_(r+1))(x)+epsilon_r.                      \tag{PR.10}
```

Taking the supremum over `f` gives the one-sided response inequality

```math
phi_r\le P_rphi_(r+1)+epsilon_r.                      \tag{PR.11}
```

Markov kernels preserve inequalities and constants, so iteration yields

```math
phi_r\le K_(r,l_r)phi_(r+l_r)+E_r.                   \tag{PR.12}
```

Choose `x_*` with `phi(x_*)=M`.  The minorization (PR.5) lets us write

```math
K_(r,l_r)(x_*,\cdot)
=alpha_rnu+(1-alpha_r)rho_r                           \tag{PR.13}
```

for a probability measure `rho_r` (the case `alpha_r=1` is immediate).
Using `phi<=M`, uniform recovery, and (PR.12),

```math
M-omega_r
\le alpha_r\int phi\,dnu+(1-alpha_r)M
    +omega_(r+l_r)+E_r.
```

Rearrangement proves the finite-`r` version of (PR.8), hence (PR.8).
Under (PR.7), `M=int phi dnu`.  Because `nu` has full support and `phi` is
continuous, a nonempty open set on which `phi<M` would make the integral
strictly smaller than `M`.  Hence `phi` is constant.  In the finite case,
if `m=min_X phi`, then

```math
M-\int phi\,dnu
=\sum_xnu({x})(M-phi(x))\ge mu(M-m),
```

which proves (PR.9). `square`

### Corollary PR.2 (balanced reorderings are enough)

Suppose that at each level one has finitely many signed equal-fibre
coordinate replications/reorderings `U_(r,x,j)` with rational weights
`w_(r,x,j)`.  If their phase marginals have refresh windows (PR.5) and

```math
\left\|T_(r,x)-
\sum_jw_(r,x,j)U_(r,x,j)^*T_(r+1,y_(r,x,j))
                   U_(r,x,j)\right\|_(2\to2)
\le epsilon_r,                                        \tag{PR.14}
```

with (PR.7), then all phases have the same limiting response.

This gives a finite spectral certificate for synchronization.  It is
strictly more structure than separate phasewise recovery, but it does not
assume pairwise equality of response carriers or even compare their
maximizers.  The common branch weights are essential: cancellation occurs
at the operator level before the absolute response is taken.

## 3. The condition is not disguised constancy

Let `H` be any symmetric regular Hadamard matrix and let `H_r=H^(tensor r)`.
For each phase in a fixed finite set, flip the signs on an arbitrary
phase-dependent perfect matching of `H_r`, and then hollow the matrix.
Write the result as `A_(r,x)`.  The different finite-level phase responses
are not required to agree.

Use the regular Boolean eigenvector of `H` to replicate every coarse spin
into the next tensor level.  The unperturbed normalized operator pulls back
exactly.  A matching flip has operator norm two, and diagonal deletion has
operator norm one.  Therefore, for *any* prescribed phase kernel `P`, the
average pullback defect is at most

```math
epsilon_r\le {3\over\sqrt{h^r}}+{3\over\sqrt{h^(r+1)}}. \tag{PR.15}
```

Taking a full-support fixed kernel makes `alpha_r` bounded below and gives
phase collapse by Theorem PR.1.  This is a nontrivial dense-sign class: the
phase perturbations may change `Theta(h^r)` edges at every level, and no
finite-level response equality was assumed.  Uniform recovery is explicit:
the normalized operator is within `3/sqrt(h^r)` of the unperturbed Hadamard
operator, so

```math
sup_x|phi_r(x)-1|\le {3\over\sqrt{h^r}}.               \tag{PR.15a}
```

This example also collapses directly from operator closeness; it demonstrates
that the certificate is realizable inside dense signings, not that PR.1 is
needed for this easiest class.

## 4. Scale-sharp falsifier: vanishing refresh can preserve a phase forever

Full support at every individual scale and a transfer defect tending to
zero do **not** suffice.  Take `X={0,1}`, one-dimensional Hilbert spaces,
and

```math
T_(r,0)=0,
\qquad T_(r,1)=I                                      \tag{PR.16}
```

at every level.  Thus `phi_r(0)=0`, `phi_r(1)=1` for all `r`.  Let `nu` be
uniform on `X`, `a_r=2^(-r)`, and

```math
P_r(x,\cdot)=(1-a_r)delta_x+a_rnu.                    \tag{PR.17}
```

Use the identity Boolean contraction on every branch.  The operator defect
in (PR.3) is exactly `epsilon_r=a_r/2`, and every `P_r(x,\cdot)` has full
support.  Nevertheless the limiting profile is nonconstant.  The ratio of
transfer error to newly mixed mass is exactly `1/2`, and the total future
refresh `sum_(j>=r)a_j` tends to zero.  This meets every tempting weakening
of PR.1 in which one asks only for pointwise full support and
`epsilon_r->0`; it fails precisely the error-versus-refresh condition
(PR.7).

The same example shows why a branch distribution allowed to depend on the
maximizing witness is dangerous.  Without one witness-independent phase
marginal, the apparent support of all branches says nothing about how much
extremal mass is actually transported.

## 5. Mandatory Walsh-prefix test

Apply the theorem to the phase family in Theorem 30.1 for the order-four
Walsh generator.  Use the operator response

```math
Phi(t)=2L(t).
```

The explicit certificate gives

```math
Phi(1)=Phi(4)=1,
\qquad Phi(3)\ge c_*:={89\over48\sqrt3}
                 =1.07050362412\ldots.                \tag{PR.18}
```

Coordinate compression of `H_(r+1)` and diagonal deletion give
`max_t Phi(t)<=2`.  Let `lambda` be normalized Lebesgue measure on `[1,4]`
and set

```math
nu_*={99\over200}(delta_1+delta_4)+{1\over100}lambda. \tag{PR.19}
```

This measure has full topological support, while

```math
\int Phi\,dnu_*
\le {99\over100}\cdot1+{1\over100}\cdot2
=1.01.                                                \tag{PR.20}
```

Consequently any proposed Boolean-pullback refresh windows for these Walsh
prefixes satisfying

```math
K_(r,l_r)(t,\cdot)\ge alpha_rnu_*(\cdot)              \tag{PR.21}
```

must obey the quantitative obstruction

```math
\liminf_{r\to\infty}
{omega_r+omega_(r+l_r)+E_r\over alpha_r}
\ge {89\over48\sqrt3}-1.01
=0.06050362412\ldots.                                 \tag{PR.22}
```

In particular no such inter-phase transfer can have accumulated operator
defect `o(alpha_r)`.  This is stronger than merely observing that the Walsh
profile is nonconstant: it identifies exactly which putative balanced
reordering or mixing certificate must pay a fixed response-scale defect.

## 6. What the theorem says about the theory

The phase carrier has two independent quantitative resources:

```text
semantic recovery error  +  transfer defect
---------------------------------------------
             refreshed phase mass
```

If this ratio vanishes over a forgetting window, a continuous phase cannot
survive.  If it remains positive, the two-state falsifier and the Walsh
hierarchy show that macroscopic phase information may persist.

The maximum-principle step is classical.  The nontrivial content for the
present program is the operator/Boolean-pullback interface (PR.3), which
turns a finite spectral certificate into the one-sided extremal inequality
(PR.11).  Unlike a Hausdorff comparison of full response carriers, it keeps
only an averaged normalized operator and a phase transition kernel.  This
is therefore a genuine candidate synchronization mechanism rather than a
restatement that the limiting response is already constant.

The certificate is smaller in its **query obligation** than comparison of
the full Boolean response carriers, but an unconstrained raw operator
certificate may still contain quadratically many coefficients.  Compression
of its description is a separate model-specific requirement.
