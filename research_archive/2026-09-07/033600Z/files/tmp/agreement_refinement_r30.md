# Wave 30 route 3: higher-order agreement and signature-constrained refinement

## Status

- **Verified analytic lemma:** once a batch of favorable completions has few
  coordinate signatures and its *coarsest signature coset* is row-good, one
  can split the signature classes into balanced blocks without losing the
  target row-square scale. Thus ``full balanced-refinement control'' is not a
  third independent hypothesis; up to constant/additive slack it follows from
  the necessary coarsest-coset bound.
- **Sharply reduced open hypothesis:** at the mesoscopic batch size dictated by
  the fractional-cover dual, one needs a non-negligible probability of jointly
  choosing favorable completions with both few signatures and a row-good
  coarsest signature coset. This is stated precisely in Section 3.
- **Literature audit:** existing direct-product/agreement results glue
  preselected local functions from an overlap-agreement probability. They do
  not supply that probability for child-ground lists, do not control the
  quadratic row statistic, and their usual average conclusion does not by
  itself imply the required signature entropy. No theorem located supplies
  the reduced hypothesis.
- **Mechanism counterexample:** even a perfectly globally glued list can have
  all `n` coordinate signatures. List agreement without a quantitative
  signature-entropy conclusion is therefore insufficient.

The exact algebra and the finite `A_9` slack check are independently exercised
by `tmp/agreement_refinement_r30_check.py`.

## 1. Signature-constrained hashing theorem

Let `A` be an `n` by `n` symmetric zero-diagonal sign matrix and let

```math
\mathcal F=\{x^{(0)},\ldots,x^{(q-1)}\}\subseteq\{\pm1\}^n.
```

Orient the words, put `D=diag(x^(0))`, and let
`C_1,...,C_J` be their relative-signature classes as in (10.856). Write

```math
u_j=AD\mathbf1_{C_j},\qquad
C_{\rm sig}(\mathcal F)
=\max_{\epsilon\in\{\pm1\}^J}
\left\|\sum_{j=1}^J\epsilon_ju_j\right\|_2^2.
\tag{R30.1}
```

This is exactly the maximum row square of the coarsest signature coset.
Every block coset containing `F` refines the signature classes and therefore
contains this coarsest coset. Consequently `C_sig` is a **necessary** lower
bound on the row cap of every containing refinement.

Fix an integer `b>=2`. For `N_j=|C_j|`, set

```math
r_j=\begin{cases}
1,&N_j\le b,\\
\lceil2N_j/b\rceil,&N_j>b,
\end{cases}
\qquad
K_0=\sum_jr_j\le J+\frac{2n}{b}.
\tag{R30.2}
```

Assume `b>=9 log(4K_0)`. Then there is a partition `P` refining every
`C_j`, with at most `K_0` nonempty blocks and maximum block size at most `b`,
such that, with

```math
\nu=\max\{\|A\|_{\rm op}^2,b(n-1)/2\},
\qquad u=\log(4(n+K_0)),
\tag{R30.3}
```

one has

```math
\boxed{
\max_{z\in\{\pm1\}^{K_0}}\|ADPz\|_2^2
\le
\left[
\sqrt{C_{\rm sig}(\mathcal F)}
+\sqrt{K_0}\left{
\sqrt{2\nu u}+\frac23\sqrt{n-1}\,u
\right}
\right]^2.}
\tag{R30.4}
```

Zero nominal columns may be retained in (R30.4); deleting them gives the
stated nonempty partition and does not change the maximum. A convenient
expanded version is

```math
\max_z\|ADPz\|_2^2
\le 2C_{\rm sig}(\mathcal F)
+8K_0\nu u+\frac{16}{9}K_0(n-1)u^2.
\tag{R30.5}
```

### Proof

For each large class independently hash every `i in C_j` uniformly into
`r_j` nominal bins; a small class is its sole bin. Put
`v_i=x_i^(0)Ae_i`. If `g(i)` is the bin of `i`, and `bar e_j` is the
uniform vector on the `r_j` bins of its class, then

```math
ADP=M_0+Z,
\qquad
M_0[:,(j,a)]=u_j/r_j,
\qquad
Z=\sum_i v_i(e_{g(i)}-\bar e_j)^{\mathsf T}.
\tag{R30.6}
```

The summands are independent and centered. Their two rectangular variance
matrices satisfy

```math
\begin{aligned}
\sum_i\mathbb E Z_iZ_i^{\mathsf T}
&=\sum_i(1-r_{j(i)}^{-1})v_iv_i^{\mathsf T}\preceq A^2,\\
\sum_i\mathbb E Z_i^{\mathsf T}Z_i
&=\bigoplus_jN_j(n-1)
\left(r_j^{-1}I-r_j^{-2}J\right).
\end{aligned}
\tag{R30.7}
```

For `r_j>=2`, `N_j/r_j<=b/2`; for `r_j=1` the corresponding covariance
block vanishes. Also `||Z_i||_op<=sqrt(n-1)`. The rectangular matrix
Bernstein inequality therefore gives, except on a set of probability at
most `(n+K_0)e^{-u}=1/4`,

```math
\|Z\|_{\rm op}
\le\sqrt{2\nu u}+\frac23\sqrt{n-1}\,u.
\tag{R30.8}
```

For a large class the load of each bin is binomial with mean
`mu_j=N_j/r_j in (b/3,b/2]`. Since `b>=2mu_j`, the standard Chernoff bound
gives probability at most `exp(-mu_j/3)<=exp(-b/9)` that a specified bin
has load above `b`. The union bound and the assumption on `b` make the
probability of any overflow at most `1/4`. Hence a realization satisfies
both load and operator bounds.

For any block-sign vector `z`,

```math
M_0z=\sum_j a_ju_j,
\qquad
a_j=r_j^{-1}\sum_{a=1}^{r_j}z_{j,a}\in[-1,1].
```

The norm of this linear map is convex on the cube, so its maximum occurs at
a vertex and is at most `sqrt(C_sig)`. Also
`||Zz||<=sqrt(K_0)||Z||_op`. This proves (R30.4); using
`(a+b)^2<=2a^2+2b^2` proves (R30.5).

### Target-scale consequence

For an exact minimizer, `||A||_op^2<=2q_n=O(n^(3/2))`. Fix
`0<c<1/4` and take

```math
b=\Theta(n^{1/4+c}\log n),\qquad
J=O(n^{3/4-c}/\log n).
\tag{R30.9}
```

Then `K_0=O(n^(3/4-c)/log n)`, `nu=O(n^(3/2))`, and the Bernstein part
of (R30.5) is `O(n^(9/4-c))`. Therefore

```math
C_{\rm sig}(\mathcal F)=O(n^{9/4-c})
\quad\Longrightarrow\quad
\text{some balanced containing refinement has row cap }
O(n^{9/4-c}).
\tag{R30.10}
```

Thus the row bound on the **full support of one selected balanced refinement**
is now proved with the constant slack allowed by the asymptotic target. The
route does not require every possible refinement to be row-good. What remains
is to construct the family with `J` and `C_sig` at these scales.

This does not contradict the cap-`80` `A_9` triangle (10.877). For the
two-signature witness choice in (R29.10), `C_sig=80`, while its best
max-size-two refinement has cap `88`. The theorem neither preserves the
same finite cap nor controls every refinement; it promises one refinement
with a universal additive/constant loss. Those are exactly the freedoms
available in (10.875).

## 2. Favorable completions turn the theorem into common-coset agreement

For an oriented full completion `(sigma,x)` of selector `S`, retain the
deficit notation

```math
\delta_S(\sigma,x)
=Q(A[S])-\sigma x_S^{\mathsf T}A[S]x_S.
```

The same expansion as (10.861) gives

```math
\sigma x^{\mathsf T}H_Sx
\ge Y_A(S)+B_{n,m}-\delta_S(\sigma,x).
\tag{R30.11}
```

Hence `delta_S<=B_(n,m)+t` implies that every block coset containing `x`
hits `S` at tolerance `t`. In particular, for a selector batch
`S_1,...,S_r`, suppose one can choose completions `(sigma_a,x^(a))` with

```math
\delta_{S_a}(\sigma_a,x^{(a)})\le B_{n,m}+t
\quad(1\le a\le r),
\tag{R30.12}
```

whose signature family has the bounds in (R30.9)--(R30.10). Then (R30.4)
constructs one eligible balanced row-good coset hitting every selector in
the batch.

This is a genuinely higher-order condition. It does not infer a common
coset from pairwise witnesses, so the pairwise-but-not-triple `A_9` wall is
fully respected.

## 3. Exact mesoscopic missing hypothesis

Put

```math
L=n^{3/4-c},\quad
k=\Theta(L/\log n),\quad
b=\Theta(n/k)=\Theta(n^{1/4+c}\log n).
```

There are at most `M<=(2k)^n` eligible labelled block cosets. For a selector
law `w`, let `J_r` be the event that `r` iid `w`-selectors share a hit
coset, and write `p_*=max_a w(H_a)`. Then exactly

```math
p_*^r\le\Pr_w(J_r)\le Mp_*^r.
\tag{R30.13}
```

Taking

```math
r\asymp\frac{n\log(2k)}{L}
=\Theta(n^{1/4+c}\log n)
\tag{R30.14}
```

makes `log(M)/r=O(L)`. Consequently the fractional target
`p_*>=exp(-O(L))` is exponent-equivalent to
`Pr_w(J_r)>=exp(-O(rL))`, uniformly for every adversarial `w`.

By Sections 1--2, the following more structured statement is sufficient:

> **Mesoscopic completion-agreement hypothesis (open).** For every selector
> law `w`, an iid batch of size (R30.14) has probability at least
> `exp(-O(rL))` of admitting favorable completions (R30.12) whose coordinate
> signature count is `J=O(L/log n)` and whose coarsest signature coset obeys
> `C_sig=O(n^(9/4-c))`.

This is a sharp sufficient higher-order target for the signature approach;
it is stronger than the bare common-coset event in (R30.13), not an equivalent
reformulation of it. The refinement portion is no longer missing.

The spanning-tree sufficient condition (10.876) makes its stringency
visible. For a batch of size `r~b`, a tree-variation proof of
`J=O(k)` would need total projective Hamming variation `O(k)`, hence average
tree-edge variation only

```math
O(k/r)=O(n^{1/2-2c}/(\log n)^2).
\tag{R30.15}
```

Local Johnson adjacency has variation two only when the selector vertices
themselves are adjacent and a coherent choice is already available; iid
vertices from an arbitrary `w` do not come with such a tree. Equation
(R30.15), together with the independent `C_sig` bound, is a concrete
falsification test for proposed agreement mechanisms.

## 4. Why current agreement/list-recovery theorems do not supply it

The closest directly checkable primary result is Theorem 1.1 in
[Dinur--Filmus--Harsha, *Agreement tests on graphs and hypergraphs*](https://arxiv.org/abs/1711.09426),
which restates the Dinur--Steurer direct-product theorem. For one prescribed
local function on each `k`-set, when `n>=Ck`, the tested pair has intersection
`t` with both `t>=alpha k` and `k-t>=beta k`, and pairwise overlap agreement
is at least `1-epsilon`, it yields one global function agreeing **exactly**
with a `1-O_(alpha,beta)(epsilon)` fraction of local functions. These
hypotheses do not match the present problem:

1. the theorem does not cover a dense-ratio window reaching `m/n` close to
   one; applicability even at `rho=1/2` would require its unspecified
   constant to satisfy `C<=2`;
2. child optimization gives a list of grounds/near-grounds, not a prescribed
   local function and not an overlap-agreement lower bound;
3. the needed result is uniform against every selector law `w`, whereas the
   theorem uses the uniform slice testing distribution; and
4. it has no `C_sig` or row-square conclusion.

The small-soundness result of
[Dinur--Steurer, *Direct Product Testing*](https://eccc.weizmann.ac.il/report/2013/179/)
reaches test success `delta>=exp(-k)`, but its stated conclusion is agreement
with a product function on local parts of the domain. Even if a suitable
version were transplanted to this slice, the research problem supplies no
test-success premise, and the theorem contains no quadratic row control.

[Gotlib--Kaufman, *List Agreement Expansion from Coboundary Expansion*](https://arxiv.org/abs/2210.15714)
explicitly requires extra global/coboundary structure beyond the natural
local list-agreement test. No such structure has been established for the
completion lists here. More importantly, merely gluing a list is
quantitatively insufficient for the signature budget.

Indeed, take `q=1+ceil(log_2 n)` full spins. Let the base spin be all `+`,
and let the other `q-1` spins record the binary digits of the coordinate
index. Give every selector the restrictions of this same global list. The
lists are perfectly consistent on every overlap and at every higher order,
but all `n` coordinate signatures are distinct. Thus

```math
J=n,
\tag{R30.16}
```

far above `k=Theta(L/log n)`. Any applicable list-recovery theorem must
control **coordinate-signature entropy**, not just the number of decoded
global functions or their overlap consistency, and it must additionally
establish the necessary coarsest-coset row bound (R30.1).

## 5. Bottom line

There is a rigorous positive reduction but no agreement theorem completing
the route:

```text
mesoscopic favorable-completion agreement
  + J = O(L/log n)
  + coarsest signature coset C_sig = O(n^(9/4-c))
        => constrained hashing (proved here)
        => one balanced row-good common coset for the batch
        => via (R30.13), the fractional-cover exponent.
```

The next useful attack should target the two genuinely minimizer-specific
batch statistics `J` and `C_sig`. A generic Johnson/list-agreement citation,
pairwise compatibility, or a bound on the number of decoded global words is
not enough.
