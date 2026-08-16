# Carrier capacity: when composition-created sets force response information

**Status.** The metric theorem and the rank-metric application are proved
below.  A finite rank-metric instance is checked by
[`verify_phase3_carrier_capacity.py`](../experiments/verify_phase3_carrier_capacity.py).
This is an abstraction of the multichannel holonomy packing, not a claim that
the underlying distance-function identity is new to metric geometry.

The multichannel theorem leaves a natural question: which part of the
`Dk` lower bound is special to binary Hamming space, and which part is a
general response law?  The exact answer has two ingredients:

1. composition creates a family of **carriers** with large Hausdorff metric
   entropy;
2. the cost of presenting a point of a carrier is smaller than the response
   scale.

Under precisely these hypotheses, carrier entropy and response entropy agree
up to the presentation radius.  Algebraic parameter dimension alone gives no
lower bound.  A rank-metric Cayley model below shows that the law is not
confined to Hamming geometry.

## 1. Weighted carrier responses

Let `(X,d)` be a finite metric space.  For every parameter `theta`, let
`C_theta subseteq X` be a nonempty carrier and let

```math
\pi_\theta:C_\theta\longrightarrow[0,p]
```

be a presentation cost.  Define the response

```math
F_\theta(x)
=\min_{c\in C_\theta}\{d(x,c)+\pi_\theta(c)\}.               \tag{CC.1}
```

Write `d_C(x)=min_{c in C}d(x,c)`, let `d_H` be Hausdorff distance
between nonempty subsets of `X`, and put

```math
d_{\rm resp}(\theta,\theta')
=\|F_\theta-F_{\theta'}\|_\infty.                            \tag{CC.2}
```

The carrier need not be recoverable as the zero set of `F_theta`: positive
presentation costs may hide all but a few carrier points.  The number `p`
measures exactly how much this can distort carrier geometry.

### Theorem CC.1 (carrier--response quasi-isometry)

For every pair `theta,theta'`,

```math
\boxed{
\left|d_{\rm resp}(\theta,\theta')
      -d_H(C_\theta,C_{\theta'})\right|\le p.}               \tag{CC.3}
```

In particular, when all presentation costs vanish, the map

```math
C\longmapsto d_C
```

is an isometric embedding of the Hausdorff hyperspace of nonempty subsets of
`X` into `ell_infinity(X)`.

#### Proof

For every `x`, choosing a nearest carrier point and using `0<=pi<=p` gives

```math
d_{C_\theta}(x)\le F_\theta(x)
\le d_{C_\theta}(x)+p.                         \tag{CC.4}
```

Thus `F_theta=d_{C_theta}+e_theta` for a function taking values in
`[0,p]`.  Consequently

```math
\|(F_\theta-F_{\theta'})
 -(d_{C_\theta}-d_{C_{\theta'}})\|_\infty\le p.             \tag{CC.5}
```

For any two nonempty subsets of a finite metric space,

```math
\|d_C-d_{C'}\|_\infty=d_H(C,C').                            \tag{CC.6}
```

The upper bound follows from the triangle inequality.  For the reverse
bound, evaluate at a point of `C` or `C'` attaining the larger directed
Hausdorff distance.  The reverse triangle inequality for norms in (CC.5)
now proves (CC.3). `square`

The constant is `p`, not `2p`, because the perturbations in (CC.4) are
one-sided and lie in the same interval.  For arbitrary two-sided
approximations `|F-d_C|<=p`, the corresponding universal constant is `2p`.

### Proposition CC.1a (query-mass exposure)

Let `mu` be any probability measure on `X`, and let `1<=s<infinity`.  Define

```math
d_{\mu,s}^{\rm car}(C,C')
=\|d_C-d_{C'}\|_{L^s(\mu)},
\qquad
d_{\mu,s}^{\rm resp}(\theta,\theta')
=\|F_\theta-F_{\theta'}\|_{L^s(\mu)}.          \tag{CC.6a}
```

Then

```math
\boxed{
|d_{\mu,s}^{\rm resp}(\theta,\theta')
 -d_{\mu,s}^{\rm car}(C_\theta,C_{\theta'})|
\le p.}                                      \tag{CC.6b}
```

The same assertion holds for a restricted uniform query set, with the sup
norm taken only over that set.  More quantitatively, suppose a full-query
Hausdorff witness `x_0` satisfies

```math
|d_C(x_0)-d_{C'}(x_0)|=\Delta.
```

For every `t>=0`,

```math
\|F_\theta-F_{\theta'}\|_{L^s(\mu)}
\ge
(\Delta-2t-p)_+\,\mu(B(x_0,t))^{1/s}.         \tag{CC.6c}
```

#### Proof

Equation (CC.5) is a pointwise error bound by `p`, so the reverse triangle
inequality in `L^s(mu)` proves (CC.6b).  Both distance-to-set functions are
one-Lipschitz.  Hence their difference changes by at most `2t` on
`B(x_0,t)`, and its absolute value there is at least `Delta-2t`.  Subtracting
the pointwise presentation error and integrating proves (CC.6c). `square`

This is the query-mass counterpart of the Hausdorff theorem.  It does not
turn an isolated witness into a diffuse lower bound for free: the mass of its
metric neighborhood is the exact extra charge.  In Hamming and rank-metric
carriers those balls can occupy an exponentially small fraction of the query
space.  Thus the full-rate theorems below concern uniform accuracy over all
endpoints; an average-query theorem additionally needs positive exposed
mass, exactly as in the earlier tropical-witness examples.

## 2. Exact capacity consequences

For a finite pseudometric space `(Y,r)`, write

* `Pack_r(s)` for the largest subset with all pairwise distances greater
  than `s`;
* `Cov_r(s)` for the smallest radius-`s` net whose centers lie in `Y`.

Both quantities below refer to the same indexed carrier family.

### Corollary CC.2 (packing and covering transfer)

For `s>p`,

```math
\operatorname{Pack}_{d_H}(s+p)
\le \operatorname{Pack}_{d_{\rm resp}}(s)
\le \operatorname{Pack}_{d_H}(s-p).                          \tag{CC.7}
```

where the last expression is interpreted using the indexed carriers when
different parameters have the same carrier.  For `s>=p`,

```math
\operatorname{Cov}_{d_{\rm resp}}(s)
\le \operatorname{Cov}_{d_H}(s-p),                           \tag{CC.8}
```

and

```math
\operatorname{Cov}_{d_H}(s)
\le \operatorname{Cov}_{d_{\rm resp}}(s-p).                  \tag{CC.9}
```

#### Proof

All statements follow by inserting `d_H-p<=d_resp<=d_H+p` into the
definitions. `square`

If parameters with one carrier can have different costs, `d_H` is only a
pseudometric on the parameter set.  This causes no problem in (CC.7), but it
explains why the upper packing comparison becomes informative only above the
presentation scale; below `p`, differently weighted copies of the same carrier
can have distinct responses.

Define `R_epsilon(mathcal F)` to be the base-two logarithm of the smallest
number of deterministic summary states from which a decoder can approximate
every `F in mathcal F` uniformly within `epsilon`.  Decoder outputs may be
arbitrary functions on `X`.

### Corollary CC.3 (extremal rate--distortion sandwich)

For `epsilon>p`,

```math
\boxed{
\log_2\operatorname{Pack}_{d_H}(2\epsilon+p)
\le R_\epsilon(\mathcal F)
\le\log_2\operatorname{Cov}_{d_H}(\epsilon-p).}              \tag{CC.10}
```

More directly, if `mathcal A` is a carrier family with pairwise Hausdorff
distance greater than `Delta`, then its responses are pairwise separated by
more than `Delta-p`.  Any deterministic summary at error
`epsilon<(Delta-p)/2` therefore needs at least `log_2|mathcal A|` bits.

If `Theta` is uniform on `mathcal A` and a randomized message permits
uniform reconstruction with probability at least `1-eta`, then

```math
I(\Theta;S)
\ge(1-\eta)\log_2|\mathcal A|-H_2(\eta).                     \tag{CC.11}
```

#### Proof

Profiles sharing a deterministic summary state are at response distance at
most `2epsilon`, which gives the packing lower bound using (CC.7).  A response
`epsilon`-net gives a valid encoder; (CC.8) gives the displayed upper bound.
For a `Delta`-separated family, nearest-profile decoding recovers `Theta`
whenever reconstruction succeeds.  Equation (CC.11) is Fano's inequality.
`square`

Consequently, for a sequence of models with macroscopic scale `a_n` and
presentation radii `p_n=o(a_n)`, normalized response packing and covering
rates are the same as Hausdorff packing and covering rates at every continuity
scale of those entropy functions.  This is the strongest general capacity
law supported by the multichannel proof: all further lower bounds must come
from actual Hausdorff entropy of the carrier class.

## 3. Recovery of the binary multichannel theorem

For the parallel-pair construction, `X=F_2^D` with Hamming distance,
`C_V=span(v_1,...,v_k)`, and

```math
\pi_V\left(\sum_j\alpha_jv_j\right)=2|\alpha|.
```

Thus `p=2k`, and (CC.3) is exactly the Grassmannian comparison in
[`phase3_multichannel_holonomy_packing.md`](phase3_multichannel_holonomy_packing.md).
Placing the carriers inside one good `[D,r,d]` linear code makes distinct
carriers Hausdorff-separated by at least `d`: if `c in C setminus C'`,
then every `c+c'` is a nonzero host-code word.  The Gaussian-binomial count
then gives `Theta(Dk)` response bits.  The Hamming proof is therefore an
instance of a general carrier-capacity law plus a model-specific source of
Hausdorff entropy.

## 4. A non-Hamming model: rank-metric shortcut landscapes

The same mechanism produces a full-rate response family in a metric whose
extreme geometry is genuinely different from Hamming space.

Let `q` be a prime power, let `E=F_{q^D}`, and let

```math
X=\operatorname{End}_{\mathbb F_q}(E)
```

with the translation-invariant rank metric

```math
d_{\rm rk}(A,B)=\operatorname{rank}_{\mathbb F_q}(A-B).      \tag{CC.12}
```

For `a in E`, write `M_a` for multiplication by `a`.  The set

```math
C_0=\{M_a:a\in E\}                                           \tag{CC.13}
```

is a `D`-dimensional `F_q`-linear subspace of `X`, and every nonzero member
has rank `D`.

For each `k`-dimensional subspace `U<=E`, choose an ordered basis
`u_1,...,u_k`.  Put `C_U={M_a:a in U}` and, for
`a=sum_i alpha_i u_i`, define

```math
\pi_U(M_a)=|\{i:\alpha_i\ne0\}|\le k.                       \tag{CC.14}
```

The associated response is

```math
F_U(A)=\min_{a\in U}
 \left\{\operatorname{rank}(A-M_a)+\pi_U(M_a)\right\}.      \tag{CC.15}
```

This is an actual Cayley shortest-path landscape.  Rank is the minimum
number of rank-one matrices summing to a residual, while `pi_U` is the
minimum number of the shortcut lines
`F_q^times M_{u_i}` needed to produce `M_a`.

### Theorem CC.4 (rank-metric carrier capacity)

For every `1<=k<D`, the profiles (CC.15), as `U` ranges over the
`k`-subspaces of `E`, form a family of size

```math
{D\brack k}_q\ge q^{k(D-k)}                                  \tag{CC.16}
```

with pairwise response distance at least

```math
D-k.                                                         \tag{CC.17}
```

In particular, for `k<=D/4`, any summary answering every rank-metric
endpoint query to error `epsilon D`, with fixed `epsilon<3/8`, needs at
least

```math
k(D-k)\log_2q\ge\frac34Dk\log_2q                             \tag{CC.18}
```

bits.

#### Proof

For distinct `U,U'`, choose `a in U setminus U'`.  For every `a' in U'`,
the difference `a-a'` is nonzero, so `M_{a-a'}` is invertible.  Hence

```math
d_{\rm rk}(M_a,C_{U'})=D.
```

The rank-metric diameter of `X` is `D`, and therefore
`d_H(C_U,C_{U'})=D`.  Theorem CC.1 with `p=k` proves (CC.17).  The count in
(CC.16) is the standard product formula for a `q`-Grassmannian; each factor
is at least `q^{D-k}`.  Finally, pairwise separation is at least `3D/4` when
`k<=D/4`, so radius-`epsilon D` decoding with `2epsilon<3/4` distinguishes
all profiles. `square`

This is a bona fide non-Hamming Cayley realization: its background moves are
low-rank updates and its carriers are subspaces of field-multiplication
operators.  The present packing proof, however, uses only that the
multiplication host is equilateral at full rank distance; it does not yet
exhibit a capacity phenomenon intrinsic to finer rank-metric ball geometry.
It still shows that the carrier theorem transfers outside the original
Hamming code landscape, while recording the exact scope of that transfer.

## 5. Why raw algebraic dimension is not a capacity law

It is tempting to promote the binary example to the assertion

```math
\text{response information}\asymp
(\dim\text{ cycle space})(\dim\text{ kernel}).               \tag{CC.19}
```

That statement is false without geometric hypotheses.

### Counterexample CC.5 (surjective holonomy collapse)

Let `K,W` be finite-dimensional vector spaces over `F_q`, with
`dim K>=dim W`, and let the parameter be a surjective linear map
`h:K->W`.  Give its image carrier zero presentation cost.  Although a raw
map has `(dim K)(dim W)` scalar coordinates and there are exponentially many
surjective maps,

```math
C_h=\operatorname{im}h=W,
\qquad F_h=d_W\equiv0.                                       \tag{CC.20}
```

Every response is identical.  The response sees the carrier, not the chosen
parametrization of that carrier.

### Counterexample CC.6 (presentation cost can erase a carrier)

Fix `c_0 in X`, let `A=diam(X)`, and assign `pi(c_0)=0` and
`pi(c)=A` to every other point of an arbitrary carrier containing `c_0`.
Then the triangle inequality gives

```math
d(x,c)+A\ge d(x,c_0),
```

so every such carrier has the identical response `F(x)=d(x,c_0)`.  Thus a
macroscopic presentation radius can erase arbitrarily rich carrier geometry.
The `-p` loss in CC.3 is indispensable.

Even with small presentation cost, a carrier class concentrated in a small
Hausdorff ball has small macroscopic response capacity.  Cycle rank and
ambient dimension say nothing about that separation.

## 6. Director-level assessment

Theorem CC.1 is the classical distance-to-set/Hausdorff identity plus a
one-sided bounded-penalty observation.  Presented in isolation, it is
elementary metric entropy, not a new mathematical theory.  The genuinely
useful project-level deduction is narrower:

> Composition forces extremal information precisely when it creates a
> Hausdorff-rich family of carriers and the cost of selecting a carrier point
> is below the target response scale.

This criterion explains both outcomes rather than just lower bounds:

* good host codes turn `k` mixed holonomies into `Theta(Dk)` information;
* synchronized or surjective maps can collapse many raw parameters to one
  carrier;
* large presentation costs hide carrier geometry; and
* once `p=o(a_n)`, the carrier Hausdorff entropy is the macroscopic response
  entropy.

The abstraction is therefore a rigorous diagnostic and a cross-model
transfer principle, but not by itself a deep replacement for metric entropy.
The next genuinely new theorem would have to predict the Hausdorff entropy of
composition-created carriers from local fragment data, or give structural
hypotheses under which those carriers synchronize into a smaller quotient.
