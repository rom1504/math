# Query sufficiency, closure, and deficiency: independent audit

Date: 2026-08-16.

Status: independent theory-incubator report.  This audits
`query_response_body.md` and `rate_distortion_report.md`; it does not alter any
main-state or final file.

## Executive verdict

The upper response roof is a sound and useful canonical quotient, but the
minimality claim needs an explicit quantifier.

* For the complete linear query family `Theta = R^d`, Theorem QR-A is correct:
  the response function and the concave roof determine one another.  The roof
  is then the coarsest deterministic exact summary, up to a bijective
  recoding.
* For a restricted declared query family `Theta`, the minimal exact statistic
  is the restricted answer function `V_H|Theta`, or equivalently the quotient
  by equality on `Theta`.  The full roof remains sufficient but is generally
  **not minimal**.
* QR-B (Minkowski addition/sup-convolution) and QR-C (one bilinear coupling)
  are correct.  QR-C has a particularly short proof by independent product
  mixtures.
* One-step exactness is not iterative closure.  A precise sufficient closure
  condition is available: every parent feature and cross potential must lie
  in the tensor product of the two child affine feature spaces.  In coordinate
  form they must be bi-affine in the child feature vectors.  Under this
  condition the entire parent roof is recoverable, not just its zero-field
  maximum.
* There is a genuine data-processing/composition inequality: full-roof error
  does not amplify under a common bi-affine composition; the parent error is
  at most the sum of the child errors.  A query-local version uses only the
  child directions induced by the other child's pure feature values.
* Ordinary total-variation Le Cam deficiency is too severe for approximate
  real-valued response reconstruction.  The natural notion here is a
  **metric (Wasserstein) deficiency on the answer-function space**, or its
  prior-averaged rate--distortion version.  It reduces exactly to expected
  query distortion when the ideal query experiment is deterministic.
* The full Boolean pinning collapse is absolute.  With `phi(x)=x`, the roof
  contains `H(x)` at every cube vertex.  Calling this a low-dimensional
  statistic would be misleading: it is a re-encoding of the full landscape,
  and the Ising hard family shows that even a structured quantized subclass
  retains `Theta(n^2)` counterfactual bits.
* QR-E is correct as a compactness and unrestricted-realization theorem.
  It does **not** prove constrained sign-matrix realizability, quantitative
  compression, or stability of exposed faces/optimizers.  Hausdorff or value
  convergence alone is insufficient for the last of these.

The resulting theory is therefore coherent only if every statement carries
three declarations: the allowed directions `Theta`, the feature algebra that
future compositions may use, and the distortion/loss class.  Omitting any of
the three turns a true local identity into a false global sufficiency claim.

## 1. Primary-source coordinates

The following are the relevant original sources, with the precise role each
plays here.

1. [Blackwell, *Equivalent Comparisons of Experiments*
   (1953)](https://doi.org/10.1214/aoms/1177729032) compares experiments by
   attainable risks and randomization/garbling.  It justifies speaking of an
   exact query experiment and its deterministic recodings, but it does not by
   itself make a landscape response map a sufficient statistic for the
   original landscape experiment.
2. [Le Cam, *Sufficiency and Approximate Sufficiency*
   (1964)](https://doi.org/10.1214/aoms/1177700372) quantifies comparison by
   randomizations and total variation.  Its operational dual is comparison
   over all bounded losses.  This is the correct classical reference, but the
   all-bounded-loss/TV metric is intentionally stronger than numerical value
   approximation.
3. [Fenchel, *On Conjugate Convex Functions*
   (1949)](https://doi.org/10.4153/CJM-1949-007-X) is the primary convex
   conjugacy source.  [Rådström's embedding theorem
   (1952)](https://doi.org/10.1090/S0002-9939-1952-0045938-2) is a primary
   source for linearizing Minkowski addition of compact convex sets using the
   Hausdorff geometry; [Hörmander's support-function treatment
   (1955)](https://archive.ymsc.tsinghua.edu.cn/pacm_paperurl/20170108203419298421913)
   gives the corresponding locally convex support embedding.  These are the
   convex facts behind QR-A, QR-B, and the response-body metric discussion.
4. [Kolmogorov and Tikhomirov, *epsilon-entropy and epsilon-capacity of sets in
   function spaces* (1959)](https://www.mathnet.ru/eng/rm7289) introduced the
   logarithmic covering and packing quantities.  The correct metric in this
   problem is the declared query pseudometric, not a norm on an undeclared raw
   parameterization.
5. In inverse optimization, [Burton and Toint
   (1992)](https://doi.org/10.1007/BF01585693) explicitly study recovering arc
   costs from shortest-path information, while [Ahuja and Orlin
   (2001)](https://doi.org/10.1287/opre.49.5.771.10607) formulate broad inverse
   linear-optimization problems.  The Walsh/pinning decoder in the
   rate--distortion draft is a particularly clean finite inverse-value-oracle
   result; neither inverse-optimization paper is needed as a black box.
6. The graph-sketching comparison must preserve its quantifiers.
   [Andoni--Krauthgamer--Woodruff
   (2014)](https://arxiv.org/abs/1403.7058) distinguish simultaneous
   `(for all)` cut approximation from one-fixed-query `(for each)` sketches.
   [Andoni--Chen--Krauthgamer--Qin--Woodruff--Zhang
   (2016)](https://doi.org/10.1145/2840728.2840753) prove quadratic lower
   bounds for general quadratic forms and sharply different bounds for PSD and
   Laplacian subclasses.  [Kapralov--Krachun
   (2019)](https://doi.org/10.1145/3313276.3316364) concern one-pass streaming
   approximation of one global Max-Cut value.  These are sound precedents,
   but none has the same oracle as exponentially many pinned optimum-value
   queries.

## 2. What “minimal exact statistic” can mean here

Fix a common finite state space, a common feature map
`phi: Omega -> R^d`, and a landscape class `H`.  For a declared set of fields
`Theta`, define

```math
A_\Theta(H)=(V_H(\theta))_{\theta\in\Theta},
\qquad
H\sim_\Theta G
\Longleftrightarrow A_\Theta(H)=A_\Theta(G).          \tag{2.1}
```

The quotient map `q_Theta: H -> H/~_Theta` is the coarsest deterministic
exact statistic in the following precise, elementary sense.

> **Proposition 2.1 (restricted-query minimality).**  A deterministic summary
> `S(H)` permits exact decoding of every query in `Theta` if and only if
> `A_Theta` factors through `S`.  Consequently `q_Theta`, or the equivalent
> object `A_Theta(H)`, factors through every exact summary.  It is unique up to
> a one-to-one recoding of its range.

This is the appropriate minimality theorem.  It is a quotient statement, not
a statement about the shortest bit representation of a real-valued function.

For `Theta=R^d`, extend the roof by `-infinity` outside `P_phi`.  QR-A says

```math
V_H(\theta)=\sup_u\{\widehat H_\phi(u)+\langle\theta,u\rangle\},
\qquad
\widehat H_\phi(u)=\inf_\theta
 \{V_H(\theta)-\langle\theta,u\rangle\}.             \tag{2.2}
```

The finite roof is proper, upper-semicontinuous, and concave, so the second
identity is precisely concave Fenchel--Moreau biconjugacy.  Thus the roof is a
canonical representative of `q_{R^d}` and QR-A is correct.

For a proper subset `Theta`, (2.2) generally cannot be inverted.  A two-state
counterexample suffices.  Let `Omega={-1,1}`, `phi(x)=x`, and declare only
`Theta={0}`.  The landscapes

```math
H(-1)=H(1)=0,
\qquad
G(-1)=0,\quad G(1)=-1                              \tag{2.3}
```

have the same answer `V_H(0)=V_G(0)=0` but different roofs.  Hence the full
roof is sufficient and strictly nonminimal for this declared query family.
A dense set of directions is sufficient to determine every finite convex
response by continuity, but an arbitrary bounded, finite, or structured query
set need not be determining.  The correct theorem should therefore say:

> The roof is minimal exactly for a query set that is determining for the
> chosen response class; `R^d` is the canonical such choice.  In general the
> minimal object is `V|Theta`.

There is also a statistical-language caveat.  Let the parameter be the
landscape `H`.  The identity experiment outputs `H`; the query experiment
outputs `A_Theta(H)`.  The latter is a deterministic garbling of the former.
For all fields, the roof is Blackwell equivalent to the **query experiment**.
It is Blackwell equivalent to the identity landscape experiment only when
the query map is injective.  Thus “query-sufficient quotient” is accurate;
unqualified “sufficient statistic for the landscape” is not.

## 3. Exact algebra audit

### 3.1 Additive composition

QR-B is correct under its stated common-field convention.  The lifted
generator set of the parent is the set sum of the two lifted generator sets,
and convex hull commutes with a finite set sum.  Therefore

```math
K_\oplus=K_1+K_2,
\qquad
\widehat H_\oplus=\widehat H_1\,\square\,\widehat H_2,
\qquad
V_\oplus=V_1+V_2.                                    \tag{3.1}
```

Here `square` denotes concave sup-convolution.  The only required warning is
semantic: if later fields can address the two children separately, the parent
feature must be `(phi_1,phi_2)`, not merely `phi_1+phi_2`.  Equation (3.1) is
not a license to forget which interface the next operation uses.

### 3.2 One bilinear coupling

QR-C is also correct.  Its cleanest proof exposes why it is only a one-step
statement.  For any `u,v`, choose distributions `lambda,mu` attaining the two
roofs.  Independence gives

```math
\widehat H_1(u)+\widehat H_2(v)+u^T Bv
=\mathbb E_{\lambda\otimes\mu}
 [H_1(X)+H_2(Y)+\phi_1(X)^TB\phi_2(Y)]
\le \max_{x,y}H_B(x,y).                               \tag{3.2}
```

Taking point masses gives the reverse inequality.  This proves QR-C without
any hidden minimax interchange.  The Jensen proof in the draft has the right
inequality direction but should explicitly keep the term `E H_2(Y)` when
passing from the convex response in `v` to the final pure maximum.

What QR-C does **not** establish is that a parent feature not expressible
through the child linear interfaces can be reconstructed.  The next section
gives the exact positive closure theorem and a sharp failure example.

## 4. New theorem: tensor-span closure under repeated composition

For child `i`, let

```math
\mathcal A_i
=\operatorname{span}\{1,\phi_{i1},\ldots,\phi_{id_i}\}
\subseteq\mathbb R^{\Omega_i}.                        \tag{4.1}
```

The tensor span `A_1 tensor A_2` consists of finite sums of products
`a(x)b(y)`.  In coordinates, its elements are exactly the restrictions of

```math
c+p^Tu+q^Tv+u^TBv                                    \tag{4.2}
```

to `u=phi_1(x), v=phi_2(y)`, modulo any affine relations in the realized
feature sets.

Let a parent be defined by

```math
H_P(x,y)=H_1(x)+H_2(y)+C(\phi_1(x),\phi_2(y)),
\qquad
\psi(x,y)=F(\phi_1(x),\phi_2(y)),                    \tag{4.3}
```

where the scalar `C` and every coordinate `F_k` have the bi-affine form
(4.2).  For `eta in R^r`, write

```math
C(u,v)+\langle\eta,F(u,v)\rangle
=c_\eta+p_\eta^Tu+q_\eta^Tv+u^TB_\eta v.            \tag{4.4}
```

> **Theorem 4.1 (exact bi-affine closure).**  The two child roofs determine
> the full parent response and hence the full parent roof.  Explicitly,
>
> ```math
> V_P^\psi(\eta)
> =c_\eta+
> \max_{u\in P_1,v\in P_2}
> \{\widehat H_1(u)+\widehat H_2(v)
>   +p_\eta^Tu+q_\eta^Tv+u^TB_\eta v\}.              \tag{4.5}
> ```
>
> Applying (2.2) to the right side reconstructs
> `widehat H_P^psi` at every parent feature value.

**Proof.**  For fixed `u,v`, take independent child mixtures attaining their
roofs.  Because (4.4) is separately affine, its value at the two means equals
its expectation under the product mixture.  The expression inside (4.5) is
therefore an average of the pure parent queried energies and cannot exceed
their maximum.  Point masses prove the opposite inequality.  Conjugacy then
recovers the parent roof.  `square`

This theorem is stronger than QR-C: it gives every future linear response of
the parent whenever the parent interface itself lies in the tensor span.

### 4.1 Query-local data processing

For a fixed parent query `eta`, define the induced child direction sets

```math
\Theta_1(\eta)
=\{p_\eta+B_\eta v:v\in\phi_2(\Omega_2)\},
\qquad
\Theta_2(\eta)
=\{q_\eta+B_\eta^Tu:u\in\phi_1(\Omega_1)\}.          \tag{4.6}
```

If `(H_1,H_2)` and `(G_1,G_2)` use the same feature maps and the same `C,F`,
then

```math
|V_{P(H)}^\psi(\eta)-V_{P(G)}^\psi(\eta)|
\le d_{\Theta_1(\eta)}(H_1,G_1)
   +d_{\Theta_2(\eta)}(H_2,G_2).                    \tag{4.7}
```

To prove it, write the parent maximum first as a maximum over `y` of
`H_2(y)+q_eta^Tv+V_{H_1}(p_eta+B_eta v)`, replace child 1, and then write the
remaining expression in the symmetric order to replace child 2.  The
elementary inequality `|max f-max g| <= ||f-g||_infinity` finishes the proof.
This is the promised query-distortion data-processing inequality.

There is also a full-roof version.  When two landscapes share the same
feature polytope `P`, conjugacy gives the exact isometry

```math
\sup_{\theta\in\mathbb R^d}|V_H(\theta)-V_G(\theta)|
=\sup_{u\in P}|\widehat H(u)-\widehat G(u)|.          \tag{4.8}
```

One inequality follows by taking suprema; the other follows by applying the
inverse formula (2.2).  Consequently Theorem 4.1 and (4.5) imply

```math
\|\widehat H_{P(H)}-\widehat H_{P(G)}\|_\infty
\le
\|\widehat H_1-\widehat G_1\|_\infty
+\|\widehat H_2-\widehat G_2\|_\infty.              \tag{4.9}
```

Thus a common exact composition is 1-Lipschitz for the sum metric.  By
induction, a composition tree amplifies leaf errors by at most their sum.
Internal approximation errors can be added in the same way.  No independence
assumption on the errors is involved.

### 4.2 A fixed-dimensional feature algebra

A particularly usable repeated theory takes one feature region `D subset
R^d`, a bi-affine map `F:D x D -> D`, and a bi-affine cross energy `C`.  Every
binary node uses (4.3).  Theorem 4.1 proves exact recursive closure for any
fixed parse tree while retaining feature dimension `d`.

If the state is also required to be independent of bracketing, the precise
additional identities are

```math
F(F(u,v),w)=F(u,F(v,w)),                              \tag{4.10}
```

and the energy two-cocycle identity

```math
C(u,v)+C(F(u,v),w)
=C(v,w)+C(u,F(v,w)).                                  \tag{4.11}
```

Equation (4.10) makes the feature composition associative; (4.11) says that
the recursively accumulated cross energy is associative.  In augmented
coordinates `(1,u)`, homogenization gives an associative bilinear algebra on
`span{(1,u):u in D}`; it extends to all of `R^(d+1)` only when `D` has full
affine span.  Iteration also requires the chosen feature domain to be closed
under `F`.  Commutativity, if desired, requires the corresponding symmetric
identities.  Scalar raw magnetization with
`F(u,v)=u+v` and pair energy `C(u,v)=beta uv` satisfies both (4.10) and
(4.11).  Normalized magnetization instead requires size-dependent weights or
a retained mass coordinate.

The full tensor lift

```math
\psi(x,y)=(1,\phi_1(x))\otimes(1,\phi_2(y))          \tag{4.12}
```

is always closed, but its dimension multiplies at every composition.  It is
an exact fallback, not a compression theorem.  A fixed-dimensional bi-affine
quotient is useful only when (4.10)--(4.11) and the intended future query
class genuinely descend to that quotient.

### 4.3 Why the tensor condition is substantive

Even a unary nonlinear new feature can resurrect information deliberately
discarded by a roof.  Let `Omega={-1,0,1}`, `phi(x)=x`, and

```math
H_a(-1)=H_a(1)=0,\qquad H_a(0)=-a,\qquad a>0.        \tag{4.13}
```

Every `H_a` has the identical roof `widehat H_a(u)=0` on `[-1,1]`.  Now
enlarge the future interface by `psi(x)=phi(x)^2` and query it with coefficient
`eta=-1`.  Then

```math
\max_x\{H_{1/2}(x)-\psi(x)\}=-1/2,
\qquad
\max_x\{H_2(x)-\psi(x)\}=-1.                         \tag{4.14}
```

The old roofs agree and the new answers differ.  The issue is not lack of a
clever composition proof: `u^2` is outside `span{1,u}`, so a nonlinear query
can expose an interior pure state lying below the old concave envelope.  By
taking the second child to be a singleton, (4.14) is also a binary-composition
counterexample.  This is the smallest possible warning against claiming
closure after an undeclared feature enlargement.

## 5. Approximate sufficiency: the natural deficiency

### 5.1 Why ordinary Le Cam deficiency is not query distortion

Let the ideal query experiment `Q_Theta` output the deterministic answer
function `A_Theta(H)`.  Classical one-sided Le Cam deficiency is

```math
\delta_{TV}(S,Q_\Theta)
=\inf_K\sup_H
 \|K P_H^S-\delta_{A_\Theta(H)}\|_{TV}.               \tag{5.1}
```

For a point mass, total variation satisfies
`TV(mu,delta_y)=1-mu({y})`.  A decoder that is always numerically very close
but never exactly equal therefore has deficiency one.  This is appropriate
for comparison over **all bounded losses**, including the discontinuous loss
`1{answer != y}`; it is not appropriate for additive response error.

The rate--distortion draft handles this correctly: it first obtains stable
exact decoding on a finite separated hard family and only then invokes
ordinary deficiency/Fano.  It should retain that separation.

### 5.2 Metric/Wasserstein query deficiency

Put a bounded metric on the answer space, for example

```math
\bar d_{\Theta,\tau}(f,g)
=\min\{1,\|f-g\|_{\infty,\Theta}/\tau\}.             \tag{5.2}
```

For bounded energies and bounded queries the cap can be omitted after a
fixed normalization.  Define

```math
\delta^W_{\Theta,\tau}(S,Q_\Theta)
:=\inf_K\sup_H
 W_1^{\bar d_{\Theta,\tau}}
   (K P_H^S,\delta_{A_\Theta(H)})
=\inf_K\sup_H
 \mathbb E_H\bar d_{\Theta,\tau}
   (\widehat A(Z),A_\Theta(H)).                       \tag{5.3}
```

The equality holds because the target law is a point mass.  This is the
natural approximate query-sufficiency notion: it is exactly minimax expected
uniform response distortion, expressed as an experiment deficiency.  Its
prior-averaged version is the operational distortion constraint in the
rate--distortion draft.  A tail version can instead minimize
`sup_H P{d_Theta>epsilon}`.

This metric deficiency has the right data-processing law.  If a map between
answer spaces is `L`-Lipschitz, pushforward contracts Wasserstein distance by
`L`.  In particular (4.7)--(4.9) show that bi-affine composition maps have
constant one for the sum metric, so two approximate child experiments yield
a parent deficiency at most the sum of their deficiencies.  Further garbling
of a sketch cannot improve the infimum in (5.3).

The terminology should remain explicit: (5.3) is a restricted-loss or
transport analogue of Le Cam deficiency, not classical TV deficiency.  Its
loss class consists of response-metric Lipschitz losses rather than every
bounded loss.

### 5.3 Metric entropy and the exact geometry

QR-D is correct for deterministic worst-case summaries.  After taking base-2
logarithms, its packing/covering sandwich is exactly the
Kolmogorov--Tikhomirov entropy/capacity comparison in the query pseudometric.
For randomized or prior-averaged summaries, the mutual-information
rate--distortion formulation is the right replacement; packing plus stable
decoding gives Fano-type converses.

Two geometric metrics should not be conflated.  For the full field family and
a common feature polytope, (4.8) identifies response distance with **vertical
sup-norm distance of roofs**.  For truncated hypographs `U_f,U_g` with the
same projection and lower base, Euclidean Hausdorff distance is instead

```math
d_H(U_f,U_g)
=\sup_{\theta\in\mathbb R^d}
 { |V_f(\theta)-V_g(\theta)|
   \over \sqrt{1+\|\theta\|_2^2} }.                  \tag{5.4}
```

Indeed, a unit support direction with positive last coordinate is a positive
multiple of `(theta,1)`; horizontal and downward directions depend only on
the common projection and lower base.  Thus ordinary Hausdorff distance is a
**weighted** all-direction response metric.  On a bounded `Theta`, the draft's
one-way Lipschitz estimate is correct, but it is only a directional support
metric and need not determine the full roof.

Finally, `K_epsilon` counts messages while treating the codebook/decoder as
shared and free.  This is standard communication complexity, but compactness
of the response-body class alone is not a constructive coding theorem and
does not control encoding time or the description length of a universally
stored codebook.

## 6. Audit of the finite rate--distortion theorem

The main finite Ising theorem is algebraically sound.

1. If `x` differs from the pinned state `u` in `k` coordinates, the field
   loses `2Mk`, while only `k(n-k)` quadratic edges change and their maximum
   possible gain is `2ak(n-k)`.  Hence `M>a(n-1)` pins `u` uniquely and gives
   equation (4.6) of that draft.
2. Subtracting `Mn` leaves `q_A(u)-c_A`.  The normalization `c_A` is the
   degree-zero Walsh coefficient, while the degree-two coefficient is exactly
   `a A_ij`.  Bessel's inequality and sign thresholding therefore turn mean
   squared query error into edgewise Hamming error exactly as claimed.
3. Entropy subadditivity and concavity of binary entropy give
   `I(A;Z)>=N[1-h_2(D)]`.  The argument permits randomized sketches and does
   not assume computationally efficient decoding.
4. Uniform additive error `<a` makes every degree-two coefficient error
   `<a`, so it recovers all edge signs.  The strict threshold and the `2a`
   packing separation are correct.
5. The bounded pair-coupling intervention also works: a nonconstant
   `z_i=u_i x_i` loses `2Lk(n-k)` in the rank-one coupling and can gain at
   most `2ak(n-k)` in the unknown energy.  For `L>a`, only `x=+/-u` remain.
6. The vertex-prize Max-Cut corollary is correct.  Pinning exposes the cut
   function, whose degree-two Walsh coefficient is `-B_ij/2`; uniform error
   `<1/4` identifies each unweighted edge.

The theorem's interpretation should nevertheless stay narrow.

* It is a simultaneous common-sketch result for an exponential query family,
  not a lower bound for one Max-Cut value or one preselected field.
* External-field pinning uses magnitude `Theta(an)`; bounded per-edge
  rank-one coupling queries are a different, richer interface.
* Additive error below a coupling quantum is a lossless inverse problem.
  Coarser or multiplicative distortion can have a completely different
  entropy.
* The graph and quadratic-sketch papers cited above are comparisons of models,
  not proofs of this theorem.  Their for-each/for-all, multiplicative, and
  streaming quantifiers must not be imported silently.
* The result establishes information cost, not computational hardness.

The Blackwell corollary is also correct on the finite hard family: the complete
answer vector and `A` deterministically decode one another.  Therefore the
query experiment and identity experiment are Blackwell equivalent there.
Approximate **value** fidelity is still not TV deficiency until the stable
Walsh decoder converts it to parameter error; the draft says this and should
continue to do so.

## 7. The full-landscape and pinning collapse

This is not a technical footnote; it is the controlling negative result.
For `Omega={-1,+1}^n` and `phi(x)=x`, every `x` is an extreme point of
`P_phi=[-1,1]^n`.  If a probability measure on the cube has mean `x`, it must
be the point mass at `x`.  Consequently

```math
\widehat H_\phi(x)=H(x)
\quad\text{for every cube vertex }x.                 \tag{7.1}
```

Thus the full roof is a lossless encoding of an arbitrary landscape table.
The fact that the feature vector has only `n` coordinates does not imply a
small state: a concave polyhedral function on the cube can carry independent
data at exponentially many vertices.  For dense quadratic Ising landscapes,
the model reduces the intrinsic degrees of freedom to `N=binom(n,2)`, and the
rate--distortion theorem proves that full pinning still recovers all `N` bits.

This collapse generalizes.  If `phi` is injective and every realized feature
point is extreme, the roof records the maximum energy in each feature fiber;
with singleton fibers it records every `H(x)`.  If the query set contains
fields strong enough to expose each such point relative to the allowed energy
range, even a finite bounded pinning family can operationally reconstruct the
same data.  Unbounded `R^d` directions merely make the statement independent
of a predeclared energy bound.

Accordingly, a theory that declares all labelled pinning fields has already
declared an inverse landscape oracle.  The roof remains mathematically
minimal for that oracle, but it supplies no nontrivial compression.  Escape
requires at least one genuine restriction:

* a smaller direction family;
* quotienting labels or symmetries;
* coarse distortion that merges response functions;
* a low-complexity model class or prior; or
* a closed low-dimensional feature algebra such as Section 4, together with
  a quantitative entropy theorem for its roofs.

Merely replacing the table by a convex roof is not an escape.

## 8. Verdict on QR-E

QR-E is correct as written, subject to fixing a norm for the Hausdorff metric
and using its dual in (QR.14).

* The hyperspace of nonempty compact subsets of a compact metric ambient
  space is compact in Hausdorff distance.
* Convexity and vertical downward closure are closed properties, so the stated
  family is a compact subspace.
* Given a finite `delta`-net `S subset U`, `conv(S)` lies in `U` and remains a
  `delta`-net.  Its downward extension is still contained in `U`.  Treating
  each `(u,t) in S` as a pure finite state produces exactly that downward
  extension as the truncated hypograph of the finite landscape roof.
* Support functions satisfy
  `|h_U(p)-h_V(p)|<=||p||_* d_H(U,V)`, which gives QR.14.

Three limitations are essential.

1. The realizing states have arbitrary real features and energies.  QR-E says
   nothing about complete sign matrices, graph constraints, integrality,
   exchangeability, or any other model-specific realization problem.
2. Compactness gives finite covers at each fixed accuracy but no useful
   quantitative entropy, algorithm, or polynomial state-size bound.
3. Hausdorff convergence controls support **values**, not exposed faces or
   optimizer selections without a uniform exposure gap.  For example the
   affine roofs `f_m(u)=(-1)^m u/m` on `[-1,1]` converge uniformly to zero,
   while their zero-field maximizers alternate between the two endpoints; the
   limit exposes the entire interval.  Thus QR-E does not establish item 4 of
   the draft's proposed axiom (“approximate the exposed faces”) by itself.

The right verdict is therefore: keep QR-E as a clean ambient compactness
lemma, but do not advertise it as recovery of constrained extremal objects or
of the interfaces needed by later optimizer-sensitive composition.

## 9. Recommended sufficiency axioms

The drafts can be made precise by adopting the following axioms/checks.

1. **Declared experiment.**  State `(phi,Theta,loss)` before naming a summary
   sufficient.  Exact equivalence means equality of `V|Theta`; approximate
   equivalence uses the induced response metric or an explicitly stated
   prior-average distortion.
2. **Minimal quotient.**  Call the roof minimal only when `Theta` determines
   it.  Otherwise call `V|Theta` the minimal exact quotient and the roof an
   optional sufficient extension.
3. **Closure certificate.**  At every binary composition, require the new
   interface and cross energy to lie in `A_1 tensor A_2`, or prove a
   model-specific synchronization theorem that makes the omitted observables
   functions of the retained ones.
4. **Repeated closure.**  For a fixed-dimensional, bracketing-independent
   algebra, verify the bi-affine form, associativity (4.10), and energy
   cocycle (4.11).  Feature dimension alone is not a complexity bound.
5. **Approximate data processing.**  Report the induced child query sets in
   (4.6) or use the full-roof norm, and propagate errors by (4.7)--(4.9).
6. **Deficiency matched to loss.**  Reserve “Le Cam deficiency” without a
   modifier for TV/all bounded losses.  Use “metric query deficiency” for
   (5.3), and say when stable decoding transfers it to classical deficiency.
7. **Pinning audit.**  Before claiming compression, test whether the declared
   fields expose every feature extreme point.  If they do, quantify the
   resulting landscape information rather than counting feature coordinates.
8. **Face stability and realization.**  Separate value/body compactness from
   exposed-face stability and from constrained finite realization.  Each
   needs its own hypothesis and theorem.

With these corrections, the roof is a rigorous response object and the
bi-affine algebra supplies a genuine compositional theory.  Without them, the
full-query version is exact but often tautologically lossless, while the
restricted-query version overclaims minimality and the iterated version can
silently resurrect discarded landscape data.
