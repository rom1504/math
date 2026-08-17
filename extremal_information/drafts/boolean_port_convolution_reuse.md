# Boolean-port convolution: tree sampling, self-reuse, and Doeblin forgetting

**Status.** Rigorous theorem suite with an exact finite verifier.  Projective
row histograms compose by group convolution under tensor product.  The
response metric is nonexpansive in each independent factor, but this alone
permits sharp linear error growth under semantic reuse.  Independent leaf
occurrences give a positive declared-tree compiler; diagonal reuse of one
sample bank fails macroscopically.  A uniform Doeblin component restores
geometric forgetting.

## 1. Response convolution algebra

Let

```math
G_p=\{+-1\}^p/\{s\sim-s\}
```

with coordinatewise multiplication, identity `e`, and uniform probability
`u_p`.  Write

```math
K_p(z)={1\over p}\left|\sum_{i=1}^pz_i\right|,
\qquad
R_\mu(\epsilon)=\sum_{s\in G_p}\mu(s)K_p(s\epsilon).          \tag{CR.1}
```

Thus `d_p(mu,nu)=||R_mu-R_nu||_infty`.  Tensoring two port systems multiplies
their projective row types, so their normalized histograms compose as
`mu*lambda`.

### Theorem CR.1 (semantic tensor law and sharp generic accumulation)

For probability measures `mu,lambda` on `G_p`,

```math
\boxed{R_{\mu*\lambda}=\lambda*R_\mu,}             \tag{CR.2}
```

where `(lambda*f)(epsilon)=sum_t lambda(t)f(t epsilon)`.  Hence

```math
d_p(\mu*\lambda,\nu*\lambda)\le d_p(\mu,\nu),     \tag{CR.3}
```

and, for arbitrary probability factors,

```math
\boxed{
d_p(\mu_1*\cdots*\mu_L,\nu_1*\cdots*\nu_L)
\le\sum_{i=1}^L d_p(\mu_i,\nu_i).}                \tag{CR.4}
```

The coefficient `L` for `L`-fold semantic reuse is locally sharp.  For any
nonidentity `a in G_p`, put

```math
\mu_t=(1-t)\delta_e+t\delta_a.
```

Since `a^2=e`, for every integer `L>=1`,

```math
\mu_t^{*L}=(1-q_L(t))\delta_e+q_L(t)\delta_a,
\qquad
q_L(t)={1-(1-2t)^L\over2},                         \tag{CR.5}
```

and therefore

```math
{d_p(\mu_t^{*L},\delta_e)
 \over d_p(\mu_t,\delta_e)}={q_L(t)\over t}
\longrightarrow L\qquad(t\downarrow0).            \tag{CR.6}
```

#### Proof

Associativity and commutativity give

```math
R_{\mu*\lambda}(\epsilon)
=\sum_{s,t}\mu(s)\lambda(t)K_p(st\epsilon)
=\sum_t\lambda(t)R_\mu(t\epsilon),                \tag{CR.7}
```

which is (CR.2).  Averaging translations contracts the supremum norm, proving
(CR.3), and a telescoping replacement of the `L` factors proves (CR.4).

The parity of `L` independent Bernoulli-`t` variables is odd with probability
`q_L(t)`, proving (CR.5).  Since response is affine in the measure,

```math
d_p((1-q)\delta_e+q\delta_a,\delta_e)
=q\,d_p(\delta_a,\delta_e).                       \tag{CR.8}
```

The nonzero common factor cancels, and `q_L(t)=Lt+O_L(t^2)`, proving (CR.6).
`square`

## 2. A positive declared-tree compiler

A **declared occurrence tree** is a finite rooted binary product tree.  Its
leaves are occurrences carrying measures on `G_p`; two occurrences may carry
the same semantic measure, but they remain separate leaves.

### Theorem CR.2 (iid occurrence-tree product sketch)

For every declared occurrence tree and every `k`, assign independently to
each leaf occurrence `v` and replica `ell<=k` a sample

```math
S_{v,\ell}\sim\mu_v.                               \tag{CR.9}
```

At an internal node, multiply its two child samples coordinatewise within
each replica.  Then at every node `z`, the `k` stored types are iid with law
equal to the convolution `mu_z` represented by that subtree.  Consequently,

```math
\mathbb E d_p(\widehat\mu_z,\mu_z)\le {4\over\sqrt k},         \tag{CR.10}
```

and for all `T` nodes simultaneously, with probability at least `1-delta`,

```math
\boxed{
\max_z d_p(\widehat\mu_z,\mu_z)
\le {4\over\sqrt k}
+\sqrt{{\log(T/\delta)\over2k}}.}                 \tag{CR.11}
```

Thus independent product composition has a depth-independent approximate
response carrier with `k=O(eta^{-2}(1+log(T/delta)))` samples per node.
However, expanding a DAG that reuses one semantic leaf requires independent
sample banks for its separate **occurrences**.  CR.2 does not turn one stored
sample bank into arbitrarily many independent copies.

For finite port systems, retain the total row count separately.  If the two
tensor factors have `N_1,N_2` rows, the child has `N_1N_2` rows, so counts
form an exact multiplicative monoid.  The normalized sample state estimates
`R`; the count then recovers the unnormalized response as `pN R`.  This
scalar bookkeeping does not affect (CR.10)--(CR.11).

#### Proof

For a fixed replica, products of independent leaf samples have the desired
convolution law by induction up the tree.  Different replicas use disjoint
leaf randomness and are independent.  Apply RC.1 at each node for (CR.10),
then its bounded-difference tail exactly as in MR.1 and union bound over the
`T` node marginals. `square`

## 3. A sharp diagonal-reuse failure

The distinction between an occurrence tree and a reused state is necessary.
Consider the tempting diagonal compiler

```math
(S_1,...,S_k)\odot(S_1,...,S_k)
=(S_1^2,...,S_k^2).                                \tag{CR.12}
```

Every output equals `e`.  Yet if `H={e,a}` is a nontrivial order-two subgroup
and `u_H=(delta_e+delta_a)/2`, then

```math
u_H*u_H=u_H,
\qquad
d_p(u_H,\delta_e)={1\over2}d_p(\delta_a,\delta_e). \tag{CR.13}
```

If the projective Hamming distance from `a` to `e` is `theta p`, this failure
is exactly `theta`.  Taking the full uniform law gives an even stronger
asymptotic diagnostic:

```math
u_p*u_p=u_p,
\qquad
d_p(u_p,\delta_e)
=1-{\mathbb E|X_1+\cdots+X_p|\over p}
\longrightarrow1.                                 \tag{CR.14}
```

Here the `X_i` are iid signs; the displayed equality holds for `p>=2`.
Indeed `R_(u_p)` is the constant `c_p=E|sum_iX_i|/p`, while the range of
`R_(delta_e)` contains `1`; the elementary bound `c_p<=1/2` for `p>=2`
(check `p=2,3`, then use Cauchy--Schwarz) shows that this upper endpoint gives
the largest deviation.
Thus a sample compiler correct on independent inputs can be wrong by a
leading constant when the *same randomness* is reused.  This falsifies
diagonal coordinatewise reuse, not all possible algorithms based on a sample
(for example, all-pairs empirical convolution is a different operation).

## 4. Doeblin contraction and geometric error control

Every response difference has zero uniform mean:

```math
\mathbb E_{\epsilon\sim u_p}
 (R_\mu(\epsilon)-R_\nu(\epsilon))=0,              \tag{CR.15}
```

because the uniform average of every translate of `K_p` is the same.

### Theorem CR.3 (uniform-component forgetting)

If a probability measure `lambda` has a uniform component

```math
\lambda\ge\alpha u_p
\quad\hbox{pointwise},\qquad0\le\alpha\le1,        \tag{CR.16}
```

then

```math
\boxed{
d_p(\mu*\lambda,\nu*\lambda)
\le(1-\alpha)d_p(\mu,\nu).}                       \tag{CR.17}
```

More generally, suppose

```math
\rho_j=\rho_{j-1}*\lambda_j,
\qquad
d_p(\widetilde\rho_j,
    \widetilde\rho_{j-1}*\lambda_j)\le\eta_j,
\qquad \lambda_j\ge\alpha_j u_p.                 \tag{CR.18}
```

Then

```math
d_p(\widetilde\rho_L,\rho_L)
\le d_p(\widetilde\rho_0,\rho_0)
       \prod_{j=1}^L(1-\alpha_j)
+\sum_{i=1}^L\eta_i\prod_{j=i+1}^L(1-\alpha_j).  \tag{CR.19}
```

In particular, if `alpha_j>=alpha>0` and `eta_j<=eta`, then

```math
d_p(\widetilde\rho_L,\rho_L)
\le(1-\alpha)^Ld_p(\widetilde\rho_0,\rho_0)
+{\eta\over\alpha}.                               \tag{CR.20}
```

#### Proof

For `alpha<1`, write `lambda=alpha u_p+(1-alpha)lambda'`; the case `alpha=1`
is immediate.  If
`h=R_mu-R_nu`, then (CR.15) gives `u_p*h=0`, while convolution by `lambda'`
is a supremum-norm contraction.  Equations (CR.2) and (CR.15) yield

```math
\|\lambda*h\|_\infty
=(1-\alpha)\|\lambda'*h\|_\infty
\le(1-\alpha)\|h\|_\infty,                       \tag{CR.21}
```

proving (CR.17).  Triangle inequality in (CR.18) gives the scalar recurrence

```math
e_j\le\eta_j+(1-\alpha_j)e_{j-1}.                 \tag{CR.22}
```

Iteration proves (CR.19), and the geometric series proves (CR.20). `square`

## 5. Interpretation

Tensor composition separates three dynamic regimes.

1. **Independent declared occurrences:** `k` iid samples close exactly under
   coordinatewise product, with no depth loss in their marginal response
   accuracy.
2. **Uncontrolled semantic reuse:** the response metric is only nonexpansive
   factor by factor, and the sharp local cost is the number of occurrences.
   Reusing the same sample coordinates can fail macroscopically even on a
   two-point subgroup.
3. **Mixing continuations:** a fixed uniform component annihilates the
   mean-zero response channel and converts additive local errors into a
   bounded geometric tail.

The new ingredient relative to static response entropy is therefore not a
larger response image but the **independence or forgetting budget attached to
reuse**.  The positive results contain strictly less information than the
full row histogram: `O(p/eta^2)` response bits per declared node, plus the
randomness needed for occurrence independence.  They do not supply a finite
carrier for arbitrary unbounded DAG reuse without mixing.

## 6. Verification

The companion exact verifier checks (CR.2), (CR.5), (CR.13), and the finite
form of (CR.17), together with occurrence-tree product marginals:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_port_convolution_reuse.py
```
