# Upper response bodies: a candidate composable extremal object

Date: 2026-08-16.

Status: director draft for independent verification.  The statements below
are finite-dimensional convex-analysis claims, not a proposed solution of the
quadratic-signing convergence problem.

## 1. Why this object is being introduced

The obstruction atlas gives two exact lessons.

1. An unrooted energy histogram does not predict even one-vertex composition:
   the order-eight signings in `experiments/entropy_overlap_results.json` have
   the same complete energy histogram but caps `16` and `20` after the same
   all-negative one-vertex coupling.
2. Even the complete global energy--energy--overlap histogram does not predict
   a labeled block query.  The left/right Curie--Weiss pair in
   `overlap_entropy_report.md` has identical exact global pair data but
   left-block response `1/2` versus `1`.

The common defect is not lack of another scalar moment.  The summary omitted
the interface through which the future environment interrogates the
landscape.  The construction below therefore starts with a declared feature
map.  It is allowed to fail when the query interface is enlarged.

## 2. Definition after the obstructions

Let `Omega` be finite, let `H:Omega -> R`, and let

```math
\phi:\Omega\longrightarrow\mathbb R^d
```

be the feature visible to future linear couplings.  Put

```math
P_\phi=\operatorname{conv}\phi(\Omega).
```

For `u in P_phi`, define the **upper response roof**

```math
\widehat H_\phi(u)
=\max\left\{
 \sum_{x\in\Omega}\lambda_xH(x):
 \lambda\in\Delta(\Omega),\quad
 \sum_x\lambda_x\phi(x)=u
 \right\}.                                             \tag{QR.1}
```

Equivalently, it is the upper boundary of the lifted polytope

```math
K_\phi(H)=\operatorname{conv}
 \{(\phi(x),H(x)):x\in\Omega\}.                        \tag{QR.2}
```

The full polytope may contain irrelevant lower faces; the concave roof in
(QR.1) is the proposed exact state for this query class.

For `theta in R^d`, define the counterfactual optimum

```math
V_H(\theta)=\max_{x\in\Omega}
 \{H(x)+\langle\theta,\phi(x)\rangle\}.                \tag{QR.3}
```

## 3. Exact sufficiency and minimality

> **Theorem QR-A (response duality).**  For every finite landscape and feature
> map,
>
> ```math
> V_H(\theta)
> =\max_{u\in P_\phi}
>   \{\widehat H_\phi(u)+\langle\theta,u\rangle\},       \tag{QR.4}
> ```
>
> and conversely
>
> ```math
> \widehat H_\phi(u)
> =\inf_{\theta\in\mathbb R^d}
>   \{V_H(\theta)-\langle\theta,u\rangle\}.             \tag{QR.5}
> ```

Thus two landscapes give the same answer to every declared linear query if
and only if their upper response roofs agree.  In this exact sense, the roof
is the minimal quotient of a landscape sufficient for the query experiment.

**Proof.**  Maximizing a linear functional over (QR.2) may be done on its
generators, which gives (QR.4).  The function in (QR.1) is finite,
upper-semicontinuous, and concave on the compact polytope `P_phi`.
Concave Fenchel--Moreau duality applied to (QR.4) gives (QR.5).  Therefore the
response function and roof determine one another.  `square`

This is not the tautology “store the maximum.”  It determines all optimum
responses to a declared family of perturbations, and it discards every
feature--energy point lying strictly below the upper concave envelope.

## 4. Exact composition laws

### 4.1 Additive composition

Let `(H_i,phi_i)` be landscapes with features in the same vector space and
form

```math
H_\oplus(x,y)=H_1(x)+H_2(y),
\qquad
\phi_\oplus(x,y)=\phi_1(x)+\phi_2(y).
```

> **Theorem QR-B (sup-convolution law).**
>
> ```math
> \widehat H_{\oplus}(u)
> =\max_{u_1+u_2=u}
>   \{\widehat H_1(u_1)+\widehat H_2(u_2)\}.             \tag{QR.6}
> ```
>
> Equivalently, `K_phi(H_oplus)=K_phi1(H_1)+K_phi2(H_2)` under the lifted
> addition `(u,e)+(v,f)=(u+v,e+f)`, and
>
> ```math
> V_{H_\oplus}(\theta)=V_{H_1}(\theta)+V_{H_2}(\theta). \tag{QR.7}
> ```

**Proof.**  The generators of the lifted parent polytope are exactly all sums
of one generator from each child.  The convex hull of a set sum is the
Minkowski sum of the convex hulls.  Taking its upper boundary gives (QR.6),
and support functions turn Minkowski sum into addition, giving (QR.7).
`square`

### 4.2 One bilinear coupling

Let the feature spaces have dimensions `d_1,d_2`, let `B` be a fixed
`d_1` by `d_2` matrix, and put

```math
H_B(x,y)=H_1(x)+H_2(y)
          +\phi_1(x)^{\mathsf T}B\phi_2(y).             \tag{QR.8}
```

> **Theorem QR-C (one-step bilinear maximum).**
>
> ```math
> \max_{x,y}H_B(x,y)
> =\max_{\substack{u\in P_{\phi_1}\\v\in P_{\phi_2}}}
> \left\{
> \widehat H_1(u)+\widehat H_2(v)+u^{\mathsf T}Bv
> \right\}.                                            \tag{QR.9}
> ```

**Proof.**  The right side is at least the left by choosing point masses.
Conversely, choose distributions attaining the two roofs at `u,v`.  For
fixed `v`,

```math
\widehat H_1(u)+u^{\mathsf T}Bv
\le\max_x\{H_1(x)+\phi_1(x)^{\mathsf T}Bv\}.
```

The latter maximum is convex as a function of `v`.  Replace `v` by the mean
of a distribution attaining the second roof and apply Jensen:

```math
\max_x\{H_1(x)+\phi_1(x)^{\mathsf T}Bv\}
\le\mathbb E_y
 \max_x\{H_1(x)+\phi_1(x)^{\mathsf T}B\phi_2(y)\}.
```

After adding `E H_2(y)`, the result is at most the maximum over `(x,y)` in
(QR.8).  `square`

The theorem is useful but deliberately one-step.  To compute the complete
parent roof for a later labeled coupling, one generally needs the joint
feature correlations of mixtures of `(x,y)`, not merely their two marginal
means.  The left/right block counterexample is an exact witness: total-overlap
data is sufficient for total-overlap queries but is not closed under a
species-specific perturbation.  A closed iterative theory must either enlarge
the feature algebra or prove a synchronization theorem that makes the new
features functions of the old ones.

## 5. When the roof compresses, and when it cannot

### 5.1 A positive class

If `phi(x)` is scalar magnetization and `H(x)` depends only on that
magnetization, (QR.1) is the concave envelope of at most `n+1` numbers.
It answers every uniform-field query exactly and composes by (QR.6).  The
state has polynomial rather than exponential description size.  Curie--Weiss
and mean-field occupancy models fall in this class.

### 5.2 Full Boolean pinning retains the landscape

Take `Omega={-1,+1}^n` and `phi(x)=x`.  Every `x` is an extreme point of the
cube `P_phi=[-1,1]^n`.  The only probability distribution with mean `x` is
the point mass at `x`, so

```math
\widehat H_\phi(x)=H(x)
\qquad(x\in\{-1,+1\}^n).                               \tag{QR.10}
```

Therefore exact sufficiency for every linear pinning field retains the entire
landscape.  This is an information theorem, not a defect of the definition.
It explains why a full labeled bridge tends to reconstruct all Boolean
responses: the declared interface itself separates every cube state.

The distinction suggests that “how many features?” is the wrong first
question.  What matters is the facial complexity of `phi(Omega)` and the
metric entropy of its upper roof at the required accuracy.

## 6. Approximate rate--distortion is metric entropy of responses

Fix a class `H` of landscapes, a query set `Theta`, and define

```math
d_\Theta(H,G)
=\sup_{\theta\in\Theta}|V_H(\theta)-V_G(\theta)|.       \tag{QR.11}
```

Let `K_epsilon` be the smallest range size of a deterministic summary from
which one decoder estimates every response on `Theta` to uniform additive
error `epsilon`.

> **Theorem QR-D (packing/covering sandwich).**
>
> ```math
> \operatorname{Pack}(\mathcal H,d_\Theta,2\epsilon)
> \le K_\epsilon
> \le\operatorname{Cov}(\mathcal H,d_\Theta,\epsilon). \tag{QR.12}
> ```

Here `Pack(r)` means a set with pairwise distance strictly greater than `r`.

**Proof.**  Two landscapes sharing a summary state are each within
`epsilon` of the same decoded response, so their query distance is at most
`2epsilon`; this proves the packing bound.  An `epsilon`-net supplies a
summary by sending the index of a nearest center and using its exact response
as decoder; this proves the covering bound.  `square`

For a norm-bounded full set of support directions, the usual equality between
Hausdorff distance of convex bodies and uniform distance of their support
functions identifies (QR.11) with a directional Hausdorff metric on the upper
response bodies.  Restricted query classes deliberately induce only a
pseudometric.

The finite Ising theorem in `rate_distortion_report.md` now has a geometric
interpretation: under full pinning queries, the response roofs form a
`2a`-separated family of size `2^binom(n,2)` (hence an `r`-packing for every
`r<2a`), even after every landscape is shifted to have maximum zero.  The rate
is quadratic because the full cube interface in
(QR.10) is facially uncompressed, not because the scalar maximum is hard to
encode.

## 7. Compactness and unrestricted finite realization

Fix a compact convex feature region `C subset R^d`, an energy bound `B`, and
a norm on the finite-dimensional ambient space.  Replace
a roof by its truncated hypograph

```math
U_H=\{(u,t):u\in P_\phi,\ -B\le t\le\widehat H_\phi(u)\}
\subset C\times[-B,B].                                  \tag{QR.13}
```

Let `U(C,B)` be the family of all nonempty compact convex subsets `U` of the
ambient compact set that are vertically downward: if `(u,t) in U`, then
`(u,s) in U` for every `-B<=s<=t`.

> **Theorem QR-E (fixed-ambient compactness and realization).**  The family
> `U(C,B)` is compact in Hausdorff distance.  Every member of it is a
> Hausdorff limit of response bodies of finite landscapes whose features lie
> in `C` and whose energies lie in `[-B,B]`.  Moreover, for every bounded
> query set `Theta`,
>
> ```math
> \sup_{\theta\in\Theta}
> |h_U(\theta,1)-h_V(\theta,1)|
> \le
> \sup_{\theta\in\Theta}\|(\theta,1)\|_*\,d_H(U,V).     \tag{QR.14}
> ```

**Proof.**  The hyperspace of nonempty compact subsets of the compact ambient
set is compact.  Convexity and the downward property are closed under
Hausdorff limits, so `U(C,B)` is a closed subspace.  For realization, choose a
finite `delta`-net `S_delta subset U`.  Its convex hull lies in `U` and is
Hausdorff `delta`-close to `U`.  Treat each point `(u,t) in S_delta` as one
finite state with feature `u` and energy `t`, and take the downward extension
of its lifted convex hull.  This extension remains inside `U`, contains the
convex hull, and is therefore still `delta`-close to `U`.  Finally, support
functions are Lipschitz in Hausdorff distance, giving (QR.14).  `square`

This is a genuine compactness/recovery theorem for unrestricted finite
landscapes at fixed interface dimension.  It does **not** solve realization
inside a constrained class such as complete sign matrices: the finite states
used in the proof may have arbitrary features and energies.  The obstruction
has been localized to model-specific realizability rather than compactness of
the response object itself.  It also controls support values, not exposed
optimizer faces unless a uniform exposure gap is supplied.

## 8. Candidate axiom extracted from the theorem

A proposed extremal state should always be declared relative to a query
interface `(phi,Theta)`.  It passes the first composition test only if:

1. its decoded response agrees with the roof on `Theta`;
2. its state-size bound is expressed as metric entropy in `d_Theta`;
3. the intended composition operation maps the declared feature algebra into
   itself, or a theorem proves synchronization/closure; and
4. its finite realizations approximate both the roof and the exposed faces
   used by the next composition.

This is a candidate theory principle, not yet an axiom of a finished theory.
Theorems QR-A--QR-E are the rigorous content.  The most important falsifier is
an operation advertised as closed whose next response depends on a feature
correlation omitted by the state.
