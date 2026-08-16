# Heterogeneous mean field: exact roofs and depth-stable quantization

**Status.** Proof source for the benchmark theorem.  The finite checks are in
[`../experiments/verify_meanfield_response_state.py`](../experiments/verify_meanfield_response_state.py).

## 1. Declared experiment

A block is a finite multiset (A=\{h_i:i\in A\}\subset[-B,B]).  Its binary
local-field energy is

```math
H_A(x)=\sum_{i\in A}h_i x_i,\qquad x_i\in\{0,1\}.
```

Blocks compose by anonymous disjoint union.  A future may append another
block (C) and apply one scalar chemical potential (\lambda) to total
occupancy.  Thus

```math
{\cal R}_A(C,\lambda)
=\max_{x_A,x_C}\{H_A(x_A)+H_C(x_C)+\lambda(K_A+K_C)\}.
```

This is much poorer than an arbitrary lookup potential on occupancy.  Old
blocks are not separately labelled: if a future may address them with
different fields, the state below must retain the corresponding tuple.

Write the fields in decreasing order (a_1\ge\cdots\ge a_n), and put

```math
p_A(k)=\max_{K_A=k}H_A=\sum_{j=1}^k a_j,
\qquad
R_A(\lambda)=\max_k\{p_A(k)+\lambda k\}.
```

## 2. Exact quotient and composition

### Theorem 1 (one-parameter exposure, exact metric, and merge)

For the declared experiment:

1. (p_A) is discrete concave and its ordered slopes are precisely the
   sorted local fields.
2. Linear fields recover every occupancy fibre:

   ```math
   R_A(\lambda)=\sum_i(h_i+\lambda)_+,
   \qquad
   p_A(k)=\inf_{\lambda\in\mathbb R}
                    \{R_A(\lambda)-\lambda k\}.       \tag{2.1}
   ```

   At fixed mass and (h_i\in[-B,B]), (\lambda\in[-B,B]) suffices.
3. For blocks of the same mass,

   ```math
   \sup_{\lambda\in\mathbb R}|R_A(\lambda)-R_{A'}(\lambda)|
   =\max_k|p_A(k)-p_{A'}(k)|.                         \tag{2.2}
   ```

4. Anonymous union is max-plus convolution,

   ```math
   p_{A\sqcup C}(t)
   =\max_{k+\ell=t}\{p_A(k)+p_C(\ell)\}.             \tag{2.3}
   ```

   In slope coordinates this is sorted multiset union.
5. Up to injective recoding, (R_A), (p_A), and the slope multiset are the
   coarsest exact deterministic contextual state.  The state therefore
   closes associatively under repeated composition.

#### Proof

An exchange puts the (k) largest fields in every fixed-(k) optimizer, so
(p_A(k)-p_A(k-1)=a_k).  This proves concavity and the first formula in
(2.1).  If (1\le k<n), every
(-a_k\le\lambda\le-a_{k+1}) supports (k); the endpoint cases use
(\lambda\le-a_1) and (\lambda\ge-a_n).  This proves the biconjugacy in
(2.1), including ties.

Taking maxima proves the upper inequality in (2.2), and applying the inverse
formula proves the reverse inequality.  A size-(t) subset of a union uses
some (k) sites in the first block and (t-k) in the second, proving (2.3).
Finally, the empty future forces an exact state to determine (R_A), hence
(p_A) by (2.1), hence its slopes.  Conversely those slopes answer and
compose under every declared context.  `square`

This is not merely a conditional-table dynamic program: concavity proves
that a one-real-parameter query family exposes the entire exact state and
turns tropical convolution into ordinary histogram addition.

## 3. A depth-independent approximate state

Fix (\eta>0).  For (B>0), let

```math
M=1+\left\lceil{2B\over\eta}\right\rceil,
\qquad
\Delta={2B\over M-1}\le\eta,
```

and round every field to a nearest point of the common (M)-point grid in
([-B,B]).  Let (S_\eta(A)\) be the histogram of rounded fields.  When
(B=0), use the one-bin state.

### Theorem 2 (merge homomorphism and response rate)

For a mass-(n) block:

```math
S_\eta(A\sqcup C)=S_\eta(A)+S_\eta(C),              \tag{3.1}
```

```math
|\mathcal S_{n,\eta}|={n+M-1\choose M-1}\le(n+1)^M, \tag{3.2}
```

and the decoded profile and response satisfy

```math
\max_k|p_A(k)-\widetilde p_A(k)|\le{\eta n\over2},
\qquad
\sup_\lambda|R_A(\lambda)-\widetilde R_A(\lambda)|
\le{\eta n\over2}.                                  \tag{3.3}
```

On every merge tree of total leaf mass (N), the root error remains at most
(\eta N/2), independent of depth and bracketing.  The same bound holds
after adding any known aggregate interaction evaluated identically before
and after rounding.

#### Proof

The same sitewise quantizer is used at every leaf, so histograms add and no
site is rounded twice.  Weak compositions give (3.2).  Every selected set
of (k) sites changes energy by at most (\eta k/2); taking fixed-(k)
maxima proves the profile bound.  Alternatively,

```math
|(h+\lambda)_+-(Q(h)+\lambda)_+|\le|h-Q(h)|
```

proves the response bound for all (\lambda).  Summing the once-only
sitewise errors over all leaves proves the depth claim.  `square`

Thus, for fixed (B), choosing (\eta_N=B/\sqrt N) gives both
(O(B\sqrt N)) response error and (O(\sqrt N\log N)) state bits.  More
generally the proved presentation cost is

```math
O\bigl((1+B/\eta)\log(n+1)\bigr)\quad\hbox{bits}.    \tag{3.4}
```

On the exact grid subclass distinct histograms have distinct responses, so
(3.2) is also the exact-state lower bound; it remains necessary at error
less than (\Delta/2).  This is not a matching lower bound at macroscopic
error (\Theta(\eta n)).

For scale, a disjoint-triangular-perturbation construction gives the weaker
but genuinely macroscopic bound

```math
\log_2 K_{\varepsilon n}
\ge\Omega\left(\min\{n,\sqrt{B/\varepsilon}\}\right) \tag{3.5}
```

whenever the right side exceeds a fixed constant.  Use (q) disjoint
intervals of half-width (B/(2q)), replace (2s) copies of a centre by
(s) copies at its two endpoints, and take
(s=\lfloor n/(4q)\rfloor).  Each bit creates a response tent of height
at least (Bn/(16q^2)) on a disjoint (\lambda)-interval.  Choosing
(q^2\le B/(64\varepsilon)) separates all (2^q) responses by more than
(2\varepsilon n).  The gap between (3.4) and (3.5) is open.

## 4. A strict quotient for quadratic mean field

Now fix one coefficient (J\in\mathbb R) at every stage and set

```math
H_A^J(x)=\sum_i h_i x_i+J{K_A\choose2},
\qquad
q_A(k)=p_A(k)+J{k\choose2}.                          \tag{4.1}
```

The raw merge law contains a genuine interaction,

```math
q_{A\sqcup C}(t)
=\max_{k+\ell=t}\{q_A(k)+q_C(\ell)+Jk\ell\}.       \tag{4.2}
```

Let (\bar q_A) be the least concave majorant of the lifted points
`(k,q_A(k))`.  For a concave roof `f` on `[0,n]`, a concave roof `g` on
`[0,m]`, and `0<=t<=n+m`, define

```math
f\star_Jg
=\operatorname {cav}_t
  \max_{\substack{0\le u\le n,\ 0\le v\le m;\ u+v=t}}
  \{f(u)+g(v)+Juv\}.                                \tag{4.3}
```

### Theorem 3 (bilinear roof congruence)

Under repeated same-(J) anonymous block merges followed by a terminal
linear field:

```math
\bar q_{A\sqcup C}=\bar q_A\star_J\bar q_C.         \tag{4.4}
```

The operation is associative on realizable roofs, and ((n,\bar q_A)) is
the coarsest exact contextual state.  The roof can be strictly smaller than
the raw conditional profile.

#### Proof

At prescribed child means (u,v), take independent occupancy mixtures
attaining the two roofs.  Bilinearity gives expected cross energy (Juv).
Conversely, pure occupancies are allowed.  Hence for every terminal field
(\lambda), both sides of (4.4) have response

```math
\max_{u,v}\{\bar q_A(u)+\bar q_C(v)+Juv+\lambda(u+v)\}
=\max_{k,\ell}\{q_A(k)+q_C(\ell)+Jk\ell
                   +\lambda(k+\ell)\}.              \tag{4.5}
```

Linear biconjugacy identifies their roofs.  For three children either
bracketing has terminal response

```math
\max_{u,v,z}\{f(u)+g(v)+h(z)+J(uv+uz+vz)
                    +\lambda(u+v+z)\},              \tag{4.6}
```

so the operation is associative.  The empty future recovers the roof; (4.4)
proves sufficiency under every generated future.  `square`

For a strict example, if (0<a<\min\{B,J/2\}), the two mass-two blocks

```math
A=\{0,0\},\qquad A'=\{a,-a\}
```

have profiles ((0,0,J)) and ((0,a,J)), but the same endpoint-chord roof.
They remain indistinguishable under every declared same-(J) future.

The collapse threshold is sharp at each fixed mass.  For every
(h_i\in[-B,B]),

```math
p_A(k)-{k\over n}p_A(n)
\le {2B\,k(n-k)\over n}.                             \tag{4.7}
```

Therefore (J\ge4B/n) makes every point of (q_A) lie below its endpoint
chord, and the exact roof state at mass (n) is only the total field
(\sum_i h_i).  If (J<4B/n), taking (k) fields equal to (B) and the
remaining fields equal to (-B) violates the chord at that (k); thus no
smaller uniform threshold works.  The size-uniform sufficient threshold for
all (n\ge2) is (J\ge2B).  Equality can leave interior ties; uniqueness of
the endpoint optimizer is not claimed.

For (J\le0), the profile remains concave, so the roof retains the full
field histogram.  Positive intermediate (J) can create partial collapse.
Thus one context algebra exhibits both polynomial histogram growth and a
strict, sometimes one-dimensional, synchronized quotient.

## 5. Scope and normalization warnings

* In Ising variables (x_i=(1+\sigma_i)/2), fields, magnetization, and
  quadratic coefficients acquire the corresponding factors of two and
  additive baselines.  For physical minimization, negate the energy before
  using these max-response statements.
* (\lambda\in\mathbb R) is part of the exact query family.  Restriction to
  ([-B,B]) is valid for comparing equal-mass local-field blocks, but not
  automatically for changing masses or quadratic profiles.
* If children use a native Curie--Weiss coefficient (J/n) and a parent
  later uses (J/(n+m)), the curvature changes.  A roof formed at the child
  coefficient can have discarded points that become exposed.  Retaining the
  underlying histogram is safe; native-normalization roof closure is not
  proved.
* A non-biaffine cross interaction need not respect the scalar roof.  One
  must enlarge the feature vector or prove another model-specific
  congruence.
* The theorem controls optimal values, not optimizer identities.  Small
  perturbations can change an optimizer at a tie.

## 6. What the benchmark validates

The response framework predicts the classical sorted-field/transfer state
without assuming it: the allowed one-parameter queries expose a concave
profile, and composition forces max-plus convolution.  It also predicts two
facts not contained in the bare exact algorithm:

1. sitewise quantization is a merge homomorphism whose error is charged once
   per microscopic feature rather than once per composition layer; and
2. after a bilinear mean-field interaction destroys concavity, contextual
   equivalence discards the hidden fibres exactly, with a sharp regime in
   which only total field survives.

The lesson shared with separator and automaton benchmarks is not that all
states are boundary tables.  A useful state is a quotient congruence for the
declared future semigroup; its growth is controlled by the geometry and
algebra of the realizable response image.
