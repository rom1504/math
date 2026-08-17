# Regular-Hadamard amplification has an extremal limit

Status: rigorous task-local theorem draft.  This is a near-original
restricted model, not a statement about the minimum defining `M_n`.

## 1. Monotone amplification

Let `H` be a symmetric Hadamard matrix of order `h` and suppose it has a
Boolean eigenvector `u` with

```math
Hu=sqrt(h)u.                                             \tag{HA.1}
```

Thus `H^2=hI`, `u^Tu=h`, and `u^THu=h^(3/2)`.  Put
`H_r=H^(tensor r)` and `s_r=h^r`, with `H_0=(1)`.
The positive sign in (HA.1) is needed for the un-absolute upper functional.
If only `Hu=-sqrt(h)u` is available, the absolute statement below still
holds, but signed values reverse at each amplification. Since `Hu` is
integral, (HA.1) also forces `sqrt(h)` to be an integer.

For a fixed real symmetric `d by d` matrix `B`, define

```math
B_r=B tensor H_r,
\qquad N_r=ds_r,                                       \tag{HA.2}
```

and the upper and absolute normalized Boolean quadratics

```math
\begin{aligned}
q_r^+(B)&={1\over2N_r^(3/2)}
          \max_(x in {+-1}^N_r)x^TB_rx,\\
q_r^abs(B)&={1\over2N_r^(3/2)}
          \max_(x in {+-1}^N_r)|x^TB_rx|.
\end{aligned}                                           \tag{HA.3}
```

### Theorem HA.1 (monotone thermodynamic limit)

Both sequences in (HA.3) are nondecreasing and converge.  Quantitatively,

```math
q_r^+(B)<=q_(r+1)^+(B)<={||B||_(2->2)\over2sqrt(d)},
\qquad
q_r^abs(B)<=q_(r+1)^abs(B)<={||B||_(2->2)\over2sqrt(d)}. \tag{HA.4}
```

#### Proof

If `x` is any Boolean vector at level `r`, then `x tensor u` is Boolean and

```math
(x tensor u)^T(B_r tensor H)(x tensor u)
=(x^TB_rx)(u^THu)=h^(3/2)x^TB_rx.                     \tag{HA.5}
```

Since `N_(r+1)^(3/2)=h^(3/2)N_r^(3/2)`, every signed
Rayleigh value at level `r` reappears at level `r+1`.  This proves both
monotonicities.  On the other hand,

```math
||B_r||=||B||h^(r/2),
\qquad ||x||_2^2=dh^r,                                 \tag{HA.6}
```

so the Boolean Rayleigh bound gives the common upper bound in (HA.4).
Bounded monotone sequences converge. `square`

No subadditivity, probabilistic disorder, or optimizer consistency is used.
The regular Boolean eigenvector is the exact scale-preserving recovery map.

## 2. The limiting response set

Write a Boolean vector of length `ds_r` as `d` blocks
`x_1,...,x_d in {+-1}^(s_r)`.  Define

```math
K_r^(d)=conv\left\{
 \left(s_r^(-3/2)x_i^TH_rx_j\right)_(1<=i,j<=d):
 x_i in {+-1}^(s_r)
\right\}.                                               \tag{HA.7}
```

These are compact convex subsets of the symmetric cube
`[-1,1]^(d(d+1)/2)`.

### Theorem HA.2 (finite-dimensional amplification carrier)

The compact convex sets are nested,

```math
K_r^(d) subseteq K_(r+1)^(d),                          \tag{HA.8}
```

and converge in Hausdorff distance to the compact convex set

```math
K_infinity^(d)=closure(union_r K_r^(d)).                \tag{HA.9}
```

For every symmetric `B`,

```math
\lim_r q_r^+(B)
={1\over2d^(3/2)}\max_(K in K_infinity^(d))<B,K>,       \tag{HA.10}
```

with the analogous absolute support function for `q_r^abs`.  Convergence is
uniform for `B` in any fixed bounded coefficient set.  At entrywise
`l_infinity` accuracy `epsilon`, the limiting carrier has an external net of
size at most

```math
(1+2/epsilon)^(d(d+1)/2).                              \tag{HA.11}
```

Thus fixed outer dimension gives response complexity independent of the
amplification depth.

#### Proof

Amplify every block in (HA.7) as `x_i tensor u`.  Equation (HA.5), now with
two possibly different blocks, preserves every normalized matrix entry and
proves (HA.8).  The closure of the increasing union is compact.  Increasing
compact subsets of a fixed compact metric space converge in Hausdorff
distance to the closure of their union: cover the limit by finitely many
balls centered in the union and take the largest level of their centers.

Expanding the quadratic gives

```math
{x^T(B tensor H_r)x\over2(ds_r)^(3/2)}
={1\over2d^(3/2)}\sum_(i,j)B_(ij)K_(ij).                \tag{HA.12}
```

Linear optimization is unchanged by convexification, so Hausdorff
convergence proves (HA.10) and uniformity on bounded dual sets.  The cube
grid proves (HA.11). `square`

The set is a strict macroscopic state: it has fixed Euclidean dimension
while the Boolean landscape has `2^(dh^r)` configurations.  It is also
query-complete for every fixed-dimensional outer quadratic coefficient
matrix, not merely for one preselected scalar.

## 3. Walsh-graph and dense-signing corollaries

The order-four Walsh matrix is symmetric and regular.  For example,

```math
u=(1,1,1,-1),\qquad W_4u=2u.                           \tag{HA.13}
```

Fix a finite linear-label Walsh graph program at label dimension `m_0`:
the graph, scalar onsite weights, and labels supported on the first `m_0`
coordinates are fixed.  Extend every label by zero coordinates to
`m=m_0+r`.  After reordering old and new Walsh coordinates, its complete
quadratic matrix obeys

```math
M_m=M_(m_0) tensor W_4^(tensor r).                     \tag{HA.14}
```

Indeed, on `E_m=F_2^m direct-sum F_2^m`, regrouping old and new coordinates
gives

```math
W_(E_m)=W_(E_(m_0)) tensor W_(E_r),
\qquad D_((a,0^r))=D_a tensor I_(4^r).                 \tag{HA.14a}
```

Thus every zero-extended child `D_aW` and the common Walsh bridge acquire
the same `W_4^(tensor r)` factor. The identity need not survive newly
exposed ambient-coordinate labels, coordinate-dependent fields, or
microscopic coefficients that vary with `r`; those continuations are not
part of the corollary.

### Corollary HA.3 (a near-original structured limit)

For every fixed Walsh graph program, its normalized upper and absolute
Boolean quadratic maxima converge as `m->infinity`.  The limits are the
support functions of `K_infinity^(d)` with
`d=dim M_(m_0)`.

If the fixed outer template has only `+-1` entries and trace zero, then

```math
A_r=(B tensor W_4^(tensor r))^circ                     \tag{HA.15}
```

is a hollow `+-1` signing at every level, and its Boolean energy agrees with
the full quadratic because the trace vanishes.  Hence

```math
{Q(A_r)\over (d4^r)^(3/2)}                             \tag{HA.16}
```

converges. Complete signed Walsh block graphs give nontrivial examples of
this form.

This does not transfer a minimizer between arbitrary orders and does not
bound `M_n`.  It proves a genuine thermodynamic limit for a dense structured
signing hierarchy and identifies the finite-dimensional state responsible
for it. The natural next question is which non-tensor perturbations preserve
the nested response set up to a summable Hausdorff defect.
