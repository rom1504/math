# A cluster-tail promotion theorem for the actual-child sector--Gram state

Status: **rigorous conditional physical-scale theorem**.  This note uses only
the exact bridge law induced by two optimizing children.  It proves that the
sector--Gram tangent carrier controls the physical inverse-row cumulant once
one scalar, zero-bridge connected-cluster tail is sublinear.  Conversely, in
the same quantitative regime, a linear physical cumulant forces that cluster
tail to be linear.  The added hypothesis is strong (absolute cumulant
summability), but it is neither a conference surrogate nor a parent
external-field table.

## 1. Zero-bridge cumulants

Fix an orientation `epsilon` and let

```math
Q_{ij}=sX_iY_j,
\qquad
(s,X,Y)\sim\nu_\epsilon
```

be the exact zero-bridge actual-child prior from EO.5.  For an ordered edge
tuple `\boldsymbol e=(e_1,\ldots,e_k)`, write

```math
\kappa_\epsilon(\boldsymbol e)
=\operatorname {cum}_{\nu_\epsilon}
 (Q_{e_1},\ldots,Q_{e_k}).                         \tag{SP.1}
```

Let `rows(\boldsymbol e)` be the set of left endpoints represented in the
tuple.  At channel amplitude `u`, define the absolute connected cross-row
tail

```math
\boxed{
\mathfrak C_{\ge4}(u)
=\sum_{\substack{k\ge4\\ k\ {\rm even}}}{|u|^k\over k!}
  \sum_{|\operatorname {rows}(\boldsymbol e)|\ge2}
  |\kappa_\epsilon(\boldsymbol e)|.}               \tag{SP.2}
```

Whenever (SP.2) is finite, we additionally require that the corresponding
zero-field cumulant series converges to the log moment-generating function
at every real amplitude between `0` and `u`.  This is the usual absolute
cluster-expansion hypothesis; finiteness of a merely formal series is not
being used as a substitute for convergence.

The object in (SP.2) is child-only.  It is computed from replicas of the two
zero-bridge sector Gibbs measures.  It retains one nonnegative aggregate,
not the values of the parent pressure on the bridge cube.  Absolute values
make it deliberately stronger than the cancellation-sensitive condition
actually necessary.

**Lemma SP.1 (exact connected-cluster remainder).**  Let `h_u` be the exact
full-versus-row-erased interaction in CR.4 and let `H_2` be the quadratic
cross-row chaos in EO.11.  Under the convergence hypothesis above,

```math
h_u(B)=u^2H_2(B)+R_u(B)+c_u,
\qquad
\operatorname {osc}R_u\le2\mathfrak C_{\ge4}(u).   \tag{SP.3}
```

*Proof.*  The binary channel identity gives

```math
\log p_u(B)
=\log E_{\nu_\epsilon}
 \exp\left\{u\sum_eB_eQ_e\right\}
 -mn\log\cosh u.                                   \tag{SP.4}
```

Its absolutely convergent log-MGF expansion is the ordered cumulant series
from (SP.1).  The exact row likelihood has the same series restricted to
tuples contained in that row.  Hence subtracting all row log likelihoods
removes precisely the one-row tuples.  Central symmetry under
`X\mapsto-X` kills every odd cumulant.  The `k=2` cross-row term is exactly
`u^2H_2`; the triangle inequality bounds the absolute value of the
remaining nonconstant series by (SP.2), and therefore bounds its
oscillation by twice (SP.2). `square`

## 2. Optimizer extension support makes the row blocks subgaussian

Let `D` be the right child, of order `n`, and suppose it is an exact
minimizer of its contracted-temperature pressure.  Let `r_u` be the
canonical iid-row inverse escort at `u=t`, with inverse exponent `lambda`.
The erased-row identity EE.7 and optimizer bound EE.16 give

```math
{dr_{\rm row}\over dU_n}\le e^{C_{\rm row}},
\qquad
C_{\rm row}
=\lambda\{\delta_n(t)+2|\gamma_A(t)|\}.            \tag{SP.5}
```

At `t=beta/sqrt(N)`, `delta_n(t)<=beta^2/2`.  Thus a bounded opposite-child
sector bias makes `C_row=O(1)`.

**Lemma SP.2 (off-block quadratic MGF).**  Let `R_1,\ldots,R_m` be
independent centrally symmetric sign vectors and suppose
`d\mathcal L(R_i)/dU_n<=e^C`.  If

```math
H=\sum_{i<k}R_i^{\mathsf T}M_{ik}R_k,
\qquad
V=\sum_{i<k}\|M_{ik}\|_F^2,
```

then there are constants `a_C,b_C>0`, depending only on `C`, such that

```math
\boxed{
\log E e^{\theta H}\le b_C\theta^2V
\quad\hbox{whenever}\quad
|\theta|\,\|M\|_{\rm op}\le a_C.}                 \tag{SP.6}
```

Here `M` is the symmetric block matrix with off-diagonal blocks `M_(ik)`;
in particular `\|M\|_op<=\|M\|_F=\sqrt{2V}`.

*Proof.*  Uniform-cube Hoeffding concentration and the density domination
give

```math
\Pr\{|\langle v,R_i\rangle|>z\}
\le2e^C\exp\{-z^2/(2\|v\|_2^2)\}.                 \tag{SP.7}
```

Central symmetry and integration of this tail imply the joint subgaussian
MGF bound

```math
E e^{\langle v,R_i\rangle}
\le e^{\sigma_C^2\|v\|_2^2/2}                    \tag{SP.8}
```

for a finite `sigma_C` depending only on `C`.  Independence gives the same
bound for a concatenation of any collection of the row blocks.

Randomly bipartition the row indices.  Every off-block monomial is cut with
probability `1/2`; Jensen therefore decouples `e^(theta H)` into a constant
multiple of the MGF of a bilinear form between the two independent sides.
Conditioning on one side and applying (SP.8) turns this into the exponential
of a positive quadratic form on the other side.  The Gaussian identity

```math
e^{\alpha\|Tx\|_2^2}
=E_g e^{\sqrt{2\alpha}\langle g,Tx\rangle}         \tag{SP.9}
```

and (SP.8) bound that expectation by the corresponding Gaussian
determinant.  Expanding `-log det(I-S)` for `\|S\|op<1` yields
`b_C theta^2\|M\|_F^2`; since `\|M\|_F^2=2V`, constants can be absorbed as
in (SP.6).  This is the standard decoupled Hanson--Wright proof, included
here to note that independence is needed only between rows, not between
coordinates within a row. `square`

For `H_2`, the blocks are
`M_(ik)(j,l)=Gamma_(ik;jl)^epsilon`.  Therefore

```math
V=K_\epsilon,\qquad \|M\|_{op}\le\sqrt{2K_\epsilon}. \tag{SP.10}
```

## 3. Physical promotion and its converse

**Theorem SP.3 (sector--Gram/cluster-tail physical dichotomy).**  Let
`m+n=N`, `t=u=beta/sqrt(N)`, and let both children be exact minimizers at
the contracted temperature.  Fix `lambda,beta` and an orientation.  Assume

```math
|\gamma_A(t)|\le G,
\qquad
\lambda t^2\sqrt{2K_\epsilon}\le a_C,              \tag{SP.11}
```

where `C=lambda(beta^2/2+2G)` and `a_C` is from SP.2.  If the cluster series
in SP.1 converges at `t`, then the exact canonical interaction cumulant
satisfies

```math
\boxed{
\mathcal J_t
\le b_C\lambda^2t^4K_\epsilon
   +2\lambda\mathfrak C_{\ge4}(t).}                \tag{SP.12}
```

Consequently,

```math
K_\epsilon=O(N^2),
\qquad
\mathfrak C_{\ge4}(t)=o(N)                         \tag{SP.13}
```

(with the constant in the first bound satisfying (SP.11)) imply
`mathcal J_t=o(N)`; indeed the quadratic contribution is `O(1)`.
Conversely, under the same bias, convergence, and quadratic-size premises,

```math
\mathcal J_t\ge\eta N
\quad\Longrightarrow\quad
\boxed{
\mathfrak C_{\ge4}(t)
\ge {\eta\over2\lambda}N-O_{\beta,\lambda,G}(1).} \tag{SP.14}
```

*Proof.*  By (SP.5), every independent row block under `r_t` meets SP.2.
Apply SP.2 to `theta=-lambda t^2` and use (SP.10).  From (SP.3),

```math
\begin{aligned}
\mathcal J_t
&=\log E_{r_t}e^{-\lambda(h_t-E_{r_t}h_t)}\\
&\le\log E_{r_t}e^{-\lambda t^2(H_2-E_{r_t}H_2)}
    +\lambda\operatorname {osc}R_t,
\end{aligned}                                      \tag{SP.15}
```

which proves (SP.12).  Substitute `t^4=beta^4/N^2` for (SP.13), and
rearrange (SP.12) for (SP.14). `square`

The asymmetric statement used the right-child row decomposition.  One may
swap the children and retain whichever opposite-child sector bias is
smaller.

## 4. Information-footprint and frontier audit

The physical theorem uses only:

1. the four-coordinate sector--Gram states, which return `K_epsilon` by
   SQ.3;
2. the scalar sector bias and adjacent extension deficit controlling the
   actual canonical row law; and
3. the single absolute cluster mass `mathfrak C_(>=4)(t)`.

It does not retain `B -> log Z_parent(B)`, the inverse bridge escort, or a
table indexed by external fields.  The cluster mass can nevertheless be
hard to bound: it aggregates cumulants of all orders and all edge tuples.
Thus SP.3 has a strictly smaller **query footprint** than a bridge response
table, but it is not a finite-bit compression theorem: an exact real bound
on the all-order sum can itself be difficult to certify.  It is a
conditional closure criterion, not yet a closure theorem for actual
minimizers.

The theorem isolates a narrower physical question.  In the bounded-bias,
quadratic-tangent regime, actual optimizing children obey exactly one of:

```text
sublinear absolute connected cluster tail, hence no linear canonical gain;
or linear high-order connected cluster mass already at zero bridge.
```

The second branch is a named high-transport observable about the actual
children.  It is not escaping conditional row Renyi complexity, and it is
not visible in the sector--Gram tangent.  Proving which branch optimizing
children occupy is still open.

## 5. Audit of the tangent spectral corollary

The normalization in SQ.4 is correct.  The identity

```math
|D|+2G_D(a,b)=\operatorname {tr}(C_D^aC_D^b)\ge0
```

uses positive semidefiniteness of the two correlation matrices.  Together
with `|G_A(a,b)|<=g_A` it gives
`K<=g_A(|D|+2g_D)`.  If `K>=eta N^3` at comparable orders, this forces
`max(g_A,g_D)>=c_eta N^(3/2)`, and

```math
\lambda_{\max}(C_C^a)
\ge {\operatorname {tr}((C_C^a)^2)\over\operatorname {tr}C_C^a}
={|C|+2G_C(a,a)\over|C|}
\ge c'_\eta\sqrt N.
```

Thus the stated covariance-mode conclusion has the right power and no
missing factor.  The phrase “apply the full-correlation” in the proof is a
typographical fragment; the mathematics is unaffected.  As SQ.4 itself
notes, this is a tangent conclusion only.  SP.3 explains one rigorous way
to promote its harmless branch and names the high-order obstruction when
promotion fails.
