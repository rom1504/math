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

Whenever (SP.2) is finite, we additionally require that, for **every bridge
word** `B` and every row-restricted word obtained by making any set of rows
inactive, the corresponding zero-field cumulant series converges to the log
moment-generating function at every real amplitude between `0` and `u`.
This is the usual absolute cluster-expansion hypothesis; finiteness of a
merely formal series is not being used as a substitute for convergence.

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

## 2. Balanced optimizer extension support makes the row blocks subgaussian

Let `C,D` be exact contracted-temperature pressure minimizers.  Relabel them
so that `|gamma_C|<=|gamma_D|`, use `D` (of order `n`) as the row-base child,
and choose the relative orientation which cancels their bias signs.  Let
`r_u` be the resulting canonical iid-row inverse escort at `u=t`, with
inverse exponent `lambda`.  The sector-bias balancing theorem SB.3 gives

```math
{dr_{\rm row}\over dU_n}\le e^{C_{\rm row}},
\qquad
C_{\rm row}
=\lambda\{\delta_n(t)+\log2\}
\le\lambda\{\beta^2/2+\log2\}.                    \tag{SP.5}
```

This is uniform in both sector biases.  The row law is centrally symmetric,
because its forward likelihood is invariant under `b\mapsto-b` by the
global child-spin flip.

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

*Proof.*  Put `K_0=e^C`.  Central symmetry and density domination give the
dimension-free linear MGF bound

```math
\begin{aligned}
E e^{\langle v,R_i\rangle}
&=E\cosh\langle v,R_i\rangle\\
&\le1+K_0\{E_{U_n}\cosh\langle v,R\rangle-1\}\\
&\le1+K_0\{e^{\|v\|_2^2/2}-1\}
 \le e^{K_0\|v\|_2^2/2}.
\end{aligned}                                      \tag{SP.7}
```

The last step is Bernoulli's inequality; thus one may take
`sigma_C^2=K_0`.  Independence gives the same bound for a concatenation of
any collection of row blocks.

Let `H_delta` retain the pairs cut by a random bipartition of the row
indices.  Since `E_delta H_delta=H/2`, Jensen gives the exact decoupling

```math
e^{\theta H}\le E_\delta e^{2\theta H_\delta}.     \tag{SP.8}
```

For a fixed cut write `H_delta=R_S^TBR_T`.  Conditioning on one shore and
using (SP.7) turns its MGF into the exponential of a positive quadratic form
on the other shore.  The Gaussian identity

```math
e^{\alpha\|Tx\|_2^2}
=E_g e^{\sqrt{2\alpha}\langle g,Tx\rangle}         \tag{SP.9}
```

and a second application of (SP.7) give

```math
E e^{2\theta R_S^{\mathsf T}BR_T}
\le\det(I-4\sigma_C^4\theta^2B^{\mathsf T}B)^{-1/2}.
```

Here `\|B\|op<=\|M\|op` and `\|B\|F^2<=V`.  Using
`-\log(1-x)<=2x` for `x<=1/2`, one may take explicitly

```math
a_C={1\over2\sqrt2e^C},\qquad b_C=4e^{2C}.
```

This proves (SP.6) with no dimension-dependent prefactor.  Independence is
needed only between rows, not between coordinates within a row. `square`

For `H_2`, the blocks are
`M_(ik)(j,l)=Gamma_(ik;jl)^epsilon`.  Therefore

```math
V=K_\epsilon,\qquad \|M\|_{op}\le\sqrt{2K_\epsilon}. \tag{SP.10}
```

## 3. Physical promotion and its converse

**Theorem SP.3 (sector--Gram/cluster-tail physical dichotomy).**  Let
`m+n=N`, `t=u=beta/sqrt(N)`, and let both children be exact minimizers at
the contracted temperature.  Make the balanced row-direction and orientation
choice preceding SP.2.  Fix `lambda,beta`, put
`C=lambda(beta^2/2+log2)`, and assume

```math
K_\epsilon\le\kappa N^2,
\qquad \lambda\beta^2\sqrt{2\kappa}\le a_C,       \tag{SP.11}
```

where `a_C` is from SP.2.  If the cluster series in SP.1 converges at `t`
for every bridge word as specified above, then the exact canonical
interaction cumulant satisfies

```math
\boxed{
\mathcal J_t
\le b_C\lambda^2t^4K_\epsilon
   +2\lambda\mathfrak C_{\ge4}(t).}                \tag{SP.12}
```

Consequently,

```math
\mathfrak C_{\ge4}(t)=o(N)                         \tag{SP.13}
```

implies `\mathcal J_t=o(N)`; indeed the quadratic contribution is at most
the constant `b_C lambda^2 beta^4 kappa`.  Conversely, under the same
convergence and quadratic-size premises,

```math
\mathcal J_t\ge\eta N
\quad\Longrightarrow\quad
\boxed{
\mathfrak C_{\ge4}(t)
\ge {\eta\over2\lambda}N-O_{\beta,\lambda,\kappa}(1).} \tag{SP.14}
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

Equivalently, with

```math
\kappa_*={a_C^2\over2\lambda^2\beta^4},
```

a linear balanced canonical interaction forces either
`K_\epsilon>\kappa_*N^2`, failure of the required absolute cluster
expansion, or `\mathfrak C_{\ge4}(t)=\Omega(N)`.  This is an actual-child
structural trichotomy.  It still concerns the existential balanced
orientation and transpose; target relevance after bridge optimization is
not proved.

## 4. Information-footprint and frontier audit

The physical theorem uses only:

1. the four-coordinate sector--Gram states, which return `K_epsilon` by
   SQ.3;
2. the balanced orientation plus the adjacent extension deficit controlling
   the actual canonical row law; and
3. the single absolute cluster mass `mathfrak C_(>=4)(t)`.

It does not retain `B -> log Z_parent(B)`, the inverse bridge escort, or a
table indexed by external fields.  The cluster mass can nevertheless be
hard to bound: it aggregates cumulants of all orders and all edge tuples.
Thus SP.3 has a formally smaller **query footprint** than a bridge response
table, but no operational information advantage has been proved.  The exact
all-order scalar can require the complete high-order child Gibbs law, and
the absolute values may impose a stronger condition than the cancellation-
sensitive conclusion `J=o(N)`.  It is a conditional closure criterion, not
yet a strict lower-information reduction or a closure theorem for actual
minimizers.

The theorem isolates a narrower physical question.  Conditional on a linear
balanced canonical interaction along a subsequence, at least one of the
following nonexclusive alternatives occurs:

```text
a fixed normalized sector--Gram mass above the explicit threshold;
failure of absolute cluster expansion at physical amplitude;
or, if the canonical interaction is linear, linear high-order connected
cluster mass already at zero bridge.
```

The last branch is a named high-transport observable about the actual
children.  It is not escaping conditional row Renyi complexity and is not
visible in the sector--Gram tangent.  Proving which branch optimizing
children occupy, and pricing target relevance, are still open.

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
