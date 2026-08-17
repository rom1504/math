# Independent audit of the local-field affine shell algebra

Date: 2026-08-17.

Audited draft: `local_field_affine_shell_algebra.md`.

## 0. Verdict

**REPAIR.**  The three mathematical conclusions are correct with their
displayed constants after one necessary orientation repair.  As written,
however, the opening instruction to “orient” a spin so that
`H_A(x)=Q(A)>0` is false for an arbitrary signing: replacing `x` by `-x`
does not change a quadratic energy.  For example, the all-negative signing
has its largest absolute value on the negative-energy side.  One must choose
an independent `rho in {+-1}` with

```math
rho H_A(x)=Q(A)
```

and define all local fields and switched edge signs using `rho A`.
Equivalently, globally negate `A` when necessary and explicitly say that the
oriented matrix is renamed `A`.  Without this repair, (LA.1), (LA.4), and
(LA.5) are false on the negative cap side.

After that repair:

| Item | Verdict | Finding |
|---|---|---|
| LA.1 local-field identities | **PASS** | Correct, including the less obvious upper bound `ell_i<=Q`. |
| LA.1 arbitrary-subset defect | **PASS** | The exact expansion and constants `4kQ/n+2k(k-1)` are correct. |
| Projective size and odd closure | **PASS** | Exactly `2^k` projective members for `k<=n-1`; odd products are XOR. |
| LA.2 all endpoints/ties | **PASS** | Both `t` and `t-2epsilon_j` are odd, so no tie is hidden.  The orbit can be described exactly. |
| LA.3 response normalization | **PASS** | No missing factor of two or factor of `m`; `m=0` is exact. |
| `O(k log n)` response-state claim | **PASS WITH SCOPE REPAIR** | Valid for the already declared frame/projective histogram; it is not a labelled description of the frame.  Storing the latter independently costs an additional `n`-bit ground-state gauge. |
| Level-5/composition interpretation | **PASS WITH QUALIFICATION** | A genuine arbitrary-order, one-shot Level-5 response law, but not a reusable composition law or a frontier arrow. |

The finite verifier
`experiments/verify_local_field_affine_shell_algebra.py` exhausts every
signing through `n=5` and checks random signings and their global negations
through `n=10`.  Its default run checked 1,138 signings, all affine subsets,
all majority endpoints, and full Boolean trust responses, and returned

```text
PASS: 1138 signings; all subsets/endpoints; n <= 10
```

## 1. LA.1 from first principles

Choose `rho in {+-1}` and `x in {+-1}^n` so that

```math
rho H_A(x)=Q,
```

and put

```math
\widetilde A=rho A,
\qquad
\ell_i=x_i(\widetilde A x)_i
       =rho x_i(Ax)_i.
```

Flipping coordinate `i` gives

```math
H_{\widetilde A}(x^{\{i\}})=Q-2\ell_i.
```

Positive-side one-spin optimality gives `ell_i>=0`.  The two-sided cap,
not one-spin optimality alone, gives

```math
Q-2\ell_i\ge -Q,
```

hence `ell_i<=Q`.  Finally,

```math
\sum_i\ell_i=x^T\widetilde A x=2Q.
```

Thus (LA.1) is correct after the `rho` repair.  Notice that the repository's
`nearmin_deterministic_inequalities.md`, (2.9), already contains this exact
local-field identity, with the same informal orientation shorthand.  The
identity itself is therefore not new to the repository.

Let

```math
s_{ij}=\widetilde a_{ij}x_ix_j=rho a_{ij}x_ix_j.
```

For every subset `S`, without a size or sign restriction,

```math
\begin{aligned}
H_{\widetilde A}(x^S)
 &=Q-2\sum_{i\in S}\ell_i
       +4\sum_{\{i,j\}\subseteq S}s_{ij}.       \tag{A.1}
\end{aligned}
```

Indeed, `sum_(i in S)ell_i` counts each internal edge twice and each
crossing edge once, while flipping `S` negates exactly the crossing edges.
Since `|s_ij|=1`, (A.1) yields

```math
Q-|H_A(x^S)|
=Q-|H_{\widetilde A}(x^S)|
\le 2\sum_{i\in S}\ell_i+4\binom{|S|}{2}.       \tag{A.2}
```

For the `k` smallest nonnegative local fields,

```math
\sum_{i\in I}\ell_i\le {k\over n}\sum_i\ell_i={2kQ\over n},
```

and every `S subseteq I` has no larger field sum.  Equation (A.2) therefore
gives exactly

```math
Q-|H_A(x^S)|\le {4kQ\over n}+2k(k-1).
```

No constant is missing.  The sentence `Q-|t|<=Q-t` does not require the
draft's extra condition `t<=Q`; it holds for every real `t`.

### Projective cardinality and product closure

The actual vectors `x^S`, `S subseteq I`, are distinct.  Projectively,
`x^S=+-x^T` is possible only when `S triangle T` is empty or all of `[n]`.
The second alternative is impossible because `S triangle T subseteq I`
and `|I|=k<=n-1`.  Hence the projective cardinality is exactly `2^k`.

For an odd number `2r+1` of factors,

```math
x^{S_1}\odot\cdots\odot x^{S_{2r+1}}
=x^{S_1\triangle\cdots\triangle S_{2r+1}}.
```

This proves actual, not merely projective, odd-product closure.  Moreover,
the odd products of the `k+1` ports in (LA.6) cover the entire affine coset:
use the single-flip ports indexed by `S` when `|S|` is odd, and include
`w_0` when `|S|` is even.

The asymptotics are also correct.  If `Q=O(n^(3/2))`, then

```math
\Delta_k=O(k\sqrt n+k^2).
```

Thus `k=o(n^(3/4))` gives `Delta_k=o(n^(3/2))`, and
`k=Theta(sqrt n)` gives `Delta_k=O(n)`.

## 2. LA.2 endpoint audit, including the tie-adjacent layers

For even `k`, `p=k+1` is odd.  Hence

```math
t=\sum_{j=0}^k\epsilon_j
```

is odd and nonzero.  Every `t-2epsilon_j` is also odd and nonzero.  This
checks both kinds of coordinates; there are no unexamined tie cases.

After multiplying the selector by the common projective sign `sgn(t)`, its
affine mask is exactly

```math
S(\epsilon)
=\{i_j:\operatorname{sgn}(t-2\epsilon_j)\ne\operatorname{sgn}(t)\}.
                                                               \tag{A.3}
```

The endpoint cases simplify to

```math
S(\epsilon)=
\begin{cases}
\varnothing, & |t|\ge3,\\
\{i_j:\epsilon_j=+1\}, & t=1,\\
\{i_j:\epsilon_j=-1\}, & t=-1.
\end{cases}                                                   \tag{A.4}
```

Consequently the selector orbit is not the whole `2^k` coset, but it is a
genuinely exponential subset:

```math
\{\varnothing\}\cup
\{S:|S|=k/2\}\cup
\{S:|S|=k/2+1\}.                              \tag{A.5}
```

This is fully consistent with the lemma, which claims containment rather
than equality.  It also makes the nonlinear content more transparent than
the parenthetical argument in the draft.

## 3. LA.3 normalization and presentation size

The response in (LA.8) has the repository's standard labelled normalization:

```math
\mathcal B_A(g)
=\max_y\{|H_A(y)|+g\mathbin\cdot y\}.
```

Here `H_A(y)=y^TAy/2`, so the leading cap term is `Q`, not `Q/2` or `2Q`.
For every `y,sigma`,

```math
sigma H_A(y)+m(W\epsilon)\mathbin\cdot y
\le Q+m\|W\epsilon\|_1.
```

Since all coordinates of `W epsilon` are nonzero, its sign selector pays
the full field norm.  LA.2 gives the same witness absolute quadratic energy
at least `Q-Delta_k`; optimizing `sigma` independently pays that energy with
the positive sign.  This proves (LA.9).  At `m=0`,
`mathcal B_A(0)=Q`, so the gap is exactly zero.  The proof in fact works for
every real `m>=0`, not only integer `m`.

The row histogram calculation is correct, and projectivization already
removes the row signs `x_i`; an actual switch of `A` need not be stored to
obtain it.  More explicitly,

```math
\|W\epsilon\|_1
=(n-k)|t|+\sum_{j=1}^k|t-2\epsilon_j|.            \tag{A.6}
```

Thus, conditional on the already declared port frame, `(n,k,Q)` and the
fixed star grammar even give an `O(log n)` numerical response presentation.
The stated `O(k log n)` is a safe sparse-histogram upper bound if the `k`
one-negative row types and their counts are explicitly listed by coordinate
index.

There is nevertheless an important state/accounting distinction:

* The projective response state of the declared pair `(A,W)` does not need
  the `n` signs of `x`.
* A labelled description that must reconstruct the actual columns of `W`
  independently of `A` must store a projective ground state (`n-1` bits in
  the worst case) and the labelled set `I` (`O(k log n)` bits).
* If `A` remains available and unlimited preprocessing may recompute a
  ground state, those `n` bits are construction advice rather than retained
  response state.

Accordingly, (LA.10) is correct only under its stated/implicit
already-declared-interface convention.  It must not be advertised as an
`O(sqrt n log n)`-bit labelled description of the bare signing together
with its interface.

## 4. The coding corollary

The proposed same-orientation packing corollary is valid.  A precise form is
as follows.

> **Corollary.**  Suppose `Q(A)<=C n^(3/2)`.  Fix `c>0` and, for all
> sufficiently large `n`, let `k` be the largest even integer at most
> `c sqrt n`.  There are `exp(Omega_c(sqrt n))` masks `S subseteq I` and
> one common `rho_0 in {+-1}` such that
> `rho_0H_A(x^S)>=Q-O_{C,c}(n)` for every selected `S`, and the corresponding
> signed cut words
> `z^S=(rho_0x_i^Sx_j^S)_(i<j)` have pairwise edge-Hamming distance
> `Theta_c(n^(3/2))`.

To prove it, color every mask `S subseteq I` by an energy orientation
`rho_S` satisfying `rho_SH_A(x^S)=|H_A(x^S)|` (choose either sign at zero).
One color class `F` has size at least `2^(k-1)`.  Greedily pack `F` by
deleting a Hamming ball of radius `r-1` after each chosen mask, with
`r=floor(k/4)`.  Since

```math
\sum_{j<k/4}\binom{k}{j}\le2^{h_2(1/4)k},
```

the code has at least

```math
2^{(1-h_2(1/4))k-1}=\exp(\Omega(k))
```

members and mask distance at least `r=Theta(k)`.  If two masks differ in
`d` coordinates, their spins have projective Hamming distance `d` once
`k<=n/2`, and their same-orientation signed cut words differ on exactly

```math
d(n-d)                                                   \tag{A.7}
```

edges.  With `d in [Theta(k),k]` and `k=Theta(sqrt n)`, (A.7) is
`Theta(n^(3/2))`.  This is edge-word distance; the spins themselves are only
`Theta(sqrt n)` apart.  The affine defect is

```math
\Delta_k\le(4cC+2c^2+o(1))n.
```

The parity phrase “choose even `k=floor(c sqrt n)`” should be replaced by
“choose the largest even `k<=floor(c sqrt n)`.”

This corollary is a large same-orientation near-top packing, but not a lower
bound on response-state size: the whole code is generated by the same
`k`-generator affine grammar.  Cardinality and generative information are
different resources.

## 5. Repository comparison, Level 5, and composition

The result is not wholly subsumed, but several ingredients are already in
the repository:

1. `nearmin_deterministic_inequalities.md`, (2.9), already records
   `0<=ell_i<=Q` and `sum ell_i=2Q`.
2. `boolean_port_product_algebra_closure.md`, (PC.13)--(PC.16), already
   identifies odd port products with an affine multiplicative coset, and
   PC.1 supplies the exact version of the selector/trust-response mechanism.
3. The generic cap-plus-field upper bound and selector-witness lower bound
   in LA.3 are direct approximate analogues of that existing response
   framework.

What appears genuinely new in the repository is the combined universal
observation: the `k` smallest oriented fields of an arbitrary cap state,
with the star frame `(x,x^{i_1},...,x^{i_k})`, generate an exponentially
large absolute near-top affine algebra whose majority selectors themselves
remain near top.  It should be described as a new synthesis/benchmark, not
as a new local-field identity or a new general odd-product theorem.

At the campaign's stated taxonomy this is honestly **Level 5 in the
one-shot necessary/response sense**: it holds for every bounded-cap exact
signing and therefore in particular for genuine minimizers and
near-minimizers.  It is stronger than a shell-cardinality statement because
it gives a uniform labelled response approximation on a growing declared
query language.

It does not prove the missing Level-5 synchronization-to-composition arrow:

* `W` is selected after seeing one ground state; it is not a response theorem
  for an arbitrary predeclared dense bridge.
* the affine grammars selected from two different children have no proved
  congruence or merge rule;
* a raw physicalization of the whole canonical frame has an all-aligned
  endpoint of norm `np-2k=Theta(n^(3/2))` when `p=Theta(sqrt n)`;
* one-step error cannot simply be accumulated through
  `Theta(sqrt n)` successive insertions.

The largest endpoint is not, by itself, a proof that every useful endpoint
is leading.  For `k>=4`, choose `t=3` with `epsilon_0=1`; then one can arrange

```math
\operatorname{sgn}(W\epsilon)=x,
\qquad \|W\epsilon\|_1=3n-4,
```

so a subleading exact anchor exists.  The separate draft
`local_affine_interface_composition.md` correctly uses this to obtain a
one-step anchored insertion, while also finding only the ordinary
`Theta(n sqrt s)` random-bridge error and no reusable recurrence.  Thus the
scope sentence about the *largest* endpoint is accurate, but it should not
be read as an impossibility theorem for one-step subleading anchoring.

The statement that the construction “survives the sparse-flip coherence
obstruction” is also correct only locally: it certifies its own nonlinear
selectors directly.  It does not show that an arbitrary FB.1 shell or an
arbitrary bridge has such a selector language.

## 6. Copy-ready repairs

### Required opening/orientation replacement

Replace the opening definition through (LA.1) by:

```text
Let A be a hollow symmetric signing of order n.  Choose
rho in {+-1} and x in {+-1}^n so that rho H_A(x)=Q(A)=Q>0,
and put \widetilde A=rho A and

    ell_i=x_i(\widetilde A x)_i=rho x_i(Ax)_i.

A one-spin flip at the positive maximizer of \widetilde A and the
two-sided cap give

    0<=ell_i<=Q,       sum_i ell_i=2Q.
```

Replace the definition before (LA.5) and (LA.5) itself by:

```text
If s_ij=rho a_ij x_i x_j, direct expansion gives

    rho H_A(x^S)
      =Q-2 sum_(i in S) ell_i
         +4 sum_({i,j} subseteq S) s_ij.
```

No later absolute-energy or response formula needs to change.

### Recommended endpoint clarification

After the definition of `t`, add:

```text
Because p is odd, both t and every t-2epsilon_j are odd and hence
nonzero.  Projectively orienting by sgn(t), the flip mask is empty for
|t|>=3; it is {i_j:epsilon_j=+1} for t=1 and
{i_j:epsilon_j=-1} for t=-1.  Thus the selector orbit consists exactly
of the empty mask and the masks of sizes k/2 and k/2+1.
```

### Recommended state-size qualification

Replace the state-presentation sentence after (LA.9) by:

```text
Equivalently,

    ||W epsilon||_1=(n-k)|t|+sum_j|t-2epsilon_j|.

Hence the response of the already declared frame has an O(k log n)-bit
sparse-histogram presentation (indeed the fixed star grammar makes the
numerical data smaller).  This count does not include a labelled
description of W itself: storing W independently of A requires the
projective ground-state gauge x and the labelled set I, costing
n-1+O(k log n) bits.  If A remains available, they may instead be
recomputed during interface construction.
```

### Recommended scope replacement

Replace the final Level-5 sentence by:

```text
Without such a merge law, LA.1--LA.3 form a genuine arbitrary-order
Level-5 one-shot response theorem for a ground-state-adapted declared
interface.  They do not prove the missing synchronization or cross-order
recurrence.  The local-field identity and abstract affine-product response
mechanism are already present in the repository; the universal low-field
star synthesis is the new benchmark.
```
