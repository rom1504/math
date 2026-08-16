# Candidate definitions and surviving principles

These are working definitions forced by examples.  They are not axioms of an
established theory, and they should be renamed or discarded when a sharper
object appears.

## 1. Declare the experiment before the summary

A finite **landscape experiment** consists of

```math
(\Omega,H,\phi,\Theta),
```

where `Omega` is the state space, `H:Omega -> R` is a normalized energy,
`phi:Omega -> R^d` is the feature visible to the environment, and `Theta` is
the permitted set of linear interventions.  Its response is

```math
V_H(\theta)=
\max_{x\in\Omega}\{H(x)+\langle\theta,\phi(x)\rangle\}.
```

The normalization and the query family are mathematical data.  A summary
that preserves `V_H(0)` need not preserve `V_H(theta)`, and a state sufficient
for uniform fields need not be sufficient for a labeled bridge.

Two experiments are **exactly response-equivalent on `Theta`** when their
response functions agree on `Theta`.  For approximate questions use

```math
d_\Theta(H,G)=
\sup_{\theta\in\Theta}|V_H(\theta)-V_G(\theta)|.
```

This is generally a pseudometric: the theory should identify landscapes that
no declared query can distinguish.

## 2. The exact quotient for a linear feature interface

For `u` in `conv(phi(Omega))`, define the upper response roof

```math
\widehat H_\phi(u)=
\max\left\{
\mathbb E_\lambda H(X):
\mathbb E_\lambda\phi(X)=u
\right\}.
```

Equivalently, retain the upper faces of the lifted response body

```math
K_\phi(H)=
\operatorname{conv}\{(\phi(x),H(x)):x\in\Omega\}.
```

For the full query family `Theta=R^d`, convex duality proves that the response
function and the roof determine one another.  Hence the roof is the minimal
exact quotient for this experiment.  For a restricted `Theta`, only the
support values in those directions are operationally minimal; the whole roof
may retain unnecessary faces.

This restriction caveat is important.  “The response roof is minimal” is not
valid without specifying a determining query family.

## 3. Extremal rate--distortion

Let `H_n` be a declared class and let a summary `Z=S(H)` be decoded into a
response function.  The deterministic uniform information cost is

```math
R_n^{\rm det}(\epsilon;\Theta)
=\inf_S\left\{
\log_2|\operatorname{range}S|:
\sup_{H\in\mathcal H_n}d_\Theta(H,\widehat V_{S(H)})
\le\epsilon
\right\}.
```

Under a prior `Pi_n` and loss `ell(H,Z)`, define the Shannon version

```math
R_{\Pi_n}(D)=
\inf_{P_{Z|H},\widehat V:
\mathbb E\ell(H,Z)\le D} I(H;Z).
```

Candidate asymptotic rates must specify whether information is normalized by
`n`, by the number of interactions, or by another natural model size.  A
landscape with `Theta(n^2)` independent couplings and a Boolean optimizer with
`n` labels live at different rates.

For deterministic uniform distortion, covering and packing numbers in
`d_Theta` give the exact elementary bounds.  For ensembles, source
rate--distortion, Fano inequalities, and posterior entropy give lower bounds.

## 4. Composition must act on the state

For additive composition with a common feature space,

```math
H_\oplus(x,y)=H_1(x)+H_2(y),
\qquad
\phi_\oplus(x,y)=\phi_1(x)+\phi_2(y),
```

the roofs compose by sup-convolution.  For one bilinear coupling

```math
H_B(x,y)=H_1(x)+H_2(y)+\phi_1(x)^TB\phi_2(y),
```

the two child roofs determine the exact parent maximum.

That does **not** imply closure under iteration.  A later coupling may inspect
a joint feature of `(x,y)` not determined by the two marginal feature means.
A candidate state passes the composition test only if one of the following is
proved:

1. the feature class is closed under the operation;
2. the state is enlarged by a controlled tensor/feature algebra whose metric
   entropy remains sub-landscape scale; or
3. a synchronization theorem makes every newly created feature a function of
   the retained state.

Merely writing down the infinite hierarchy is not compression.

## 5. Support resolution is part of the topology

For a normalized sequence, define upper-tail complexity

```math
\Sigma_H^\uparrow(e)=
\liminf_{n\to\infty}{1\over n}
\log\#\{x:H_n(x)\ge e\}.
```

One must distinguish:

- an empty level (`log 0=-infinity`);
- a subexponential but nonempty level (rate zero); and
- a positive-rate cloud.

For homogeneous quadratic forms, a theorem in [`theorems.md`](theorems.md)
shows that the maximum is the closure of the positive-rate upper tail.  For a
general landscape this need not hold.  Any compactness notion must state
which edge resolution it preserves.

## 6. Rooted and unrooted information are different experiments

An energy histogram and global pair-overlap law average over the roots and
coordinate locations.  A restriction, puncture, inserted vertex, block
bridge, or external field singles out an apparatus.  Its response depends on
a rooted joint profile such as

```math
(H(x),\phi_{\rm root}(x))
```

or its upper response roof.

The code and Curie examples prove that exact unrooted pair data can coexist
with a leading rooted response gap.  Therefore rooted data should be added
only when the future query forces it, not hidden inside an ambiguous phrase
such as “overlap geometry.”

## 7. Compactness and realization are separate obligations

At fixed feature dimension and bounded energy, downward convex response
bodies are compact in Hausdorff distance.  Every such body is approximable by
an unrestricted finite landscape, and bounded response values are continuous.

This does not imply realization inside a constrained model such as complete
sign matrices, linear codes, or a fixed CSP.  A useful limiting theory needs
both:

1. compactness of its abstract extremal state; and
2. a model-specific recovery theorem at the required finite sizes.

Failure of the second is not repaired by strengthening the state until it
contains a finite optimizer.

## 8. Current survival table

| Candidate | Verdict | Reason |
|---|---|---|
| Scalar energy entropy | survives for recovering a homogeneous quadratic maximum; fails composition | Hamming noise thickens the edge, but labels are absent |
| Full global pair-overlap support | exact for every global two-replica query built from energy and total overlap; fails labeled queries | fixed-half example |
| Finite `k` global overlap hierarchy | rejected as a universal invariant | for every fixed `k`, parity-half-cube code pairs agree through `k` replicas but have a scalable covering-radius gap |
| Upper response roof for a fixed interface | survives | exact duality and one-step composition theorems |
| Full-spin response roof | rejected as compression | cube vertices recover every `H(x)` |
| Query metric entropy | survives as information definition | exact packing/covering theorem and quadratic-rate example |
| Bounded-dimensional universal state | rejected as a default assumption | scalable planted resonances evade fixed-arity tests |
| Growing feature algebra | tentative | useful only with a proved entropy bound and closure law |
| Deterministic Parisi-like overlap state | open | requires a deterministic synchronization mechanism, not vocabulary transfer |

## 9. Ten questions every candidate must answer

1. Which examples force the definition?
2. What is the declared normalization and query experiment?
3. Which previous failures become immediate in its language?
4. What does the state forget?
5. What does it retain exactly or approximately?
6. What is its packing/covering or description complexity?
7. What is its algebra under the intended composition?
8. In what topology is it compact, and is the queried extreme continuous?
9. Which finite constrained objects realize its limits?
10. What theorem about it is not the tautology “if the state converges, the
    maximum converges”?

The present roof answers these questions for fixed linear interfaces and
unrestricted finite landscapes.  Its bi-affine closure theorem also gives an
exact polynomial-state dynamic program for fixed-rank Curie--Weiss/Potts-type
mean-field ground states.  It does not yet answer them for growing bridge
interfaces in dense sign matrices.

## 10. A plausible next abstraction, not yet adopted

Repeated composition suggests a **query-generated feature algebra**: begin
with the observable seen by the environment, close it under the contractions
created by allowed couplings, and measure the metric entropy of only the
exposed response faces at each depth.  A synchronization theorem could make
this algebra collapse; without one it may grow to the entire state table.

This idea remains tentative until it yields either a bounded-complexity
closure theorem or a scalable lower bound proving that such closure requires
large information rate.
