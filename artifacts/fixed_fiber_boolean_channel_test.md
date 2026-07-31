# Boolean-channel obstruction for fixed fiber algebras

Date: 2026-07-31. This is an agent-authored theorem-level family test. It
abstracts the exact balanced-Hadamard obstruction and applies to fixed
association-scheme, character-table, and design tensor algebras possessing a
Boolean extremal channel.

## 1. Fixed-generator algebra

Let `Q` be a fixed symmetric full sign matrix of order `m` with

```math
\operatorname{tr}Q=0.                                  \tag{FC1}
```

Let `G_1,...,G_h` be symmetric full sign matrices, where `G_j` has order
`k_j`, and suppose

```math
G_j^2=k_jI,qquad
G_jv_j=\sigma_j\sqrt{k_j},v_j,qquad
v_j\in\{+1,-1\}^{k_j},\quad\sigma_j\in\{+1,-1\}.      \tag{FC2}
```

Thus each generator is a symmetric Hadamard matrix and has a Boolean vector
in an extremal eigenspace. For a word `w=(i_1,...,i_s)` put

```math
R_w=G_{i_1}\mathbin\otimes\cdots\mathbin\otimes G_{i_s},
\qquad K_w=\prod_{a=1}^s k_{i_a},                      \tag{FC3}
```

and define the zero-diagonal signing

```math
S_w=Q\mathbin\otimes R_w
    -\operatorname{diag}(Q\mathbin\otimes R_w).        \tag{FC4}
```

The off-diagonal entries in (FC4) are signs. Its order is `N_w=mK_w`.
Because a Boolean quadratic form of a diagonal matrix is its trace, (FC1)
gives, for every Boolean `x`,

```math
H_{S_w}(x)
={1\over2}x^{\mathsf T}(Q\mathbin\otimes R_w)x.        \tag{FC5}
```

This is a constant-complexity algebraic state. It includes the common
construction in which a zero-diagonal macro signing `A` is completed by a
balanced diagonal `D`, so `Q=A+D`, and every fiber block is chosen from a
fixed Hadamard/design algebra.

## 2. Persistent-witness theorem

**Theorem (Boolean-channel persistence).** Suppose a word `w` has a Boolean
witness `z` satisfying

```math
|H_{S_w}(z)|\ge cN_w^{3/2}.                            \tag{FC6}
```

Then for every continuation word `a`,

```math
\operatorname{cap}(S_{wa})\ge cN_{wa}^{3/2}.           \tag{FC7}
```

**Proof.** Tensor the Boolean eigenvectors in (FC2) along `a` to obtain a
Boolean vector `y_a`. If `K_a` is the order of the continuation, then

```math
|y_a^{\mathsf T}R_ay_a|=K_a^{3/2}.                    \tag{FC8}
```

Use `z tensor y_a` in (FC5). Kronecker multiplication gives

```math
\begin{aligned}
|H_{S_{wa}}(z\mathbin\otimes y_a)|
 &=|H_{S_w}(z)|K_a^{3/2}\\
 &\ge c(N_wK_a)^{3/2}=cN_{wa}^{3/2}.
\end{aligned}                                         \tag{FC9}
```

This proves (FC7). No upper bound, optimizer computation, or asymptotic
argument is used.

### Landing consequence

Suppose the original problem has an all-order construction

```math
M_N\le(C_*+o(1))N^{3/2}.                               \tag{FC10}
```

If one finite witness has `c>C_*`, then every descendant in (FC7) has

```math
\operatorname{cap}(S_{wa})^{2/3}-M_{N_{wa}}^{2/3}
\ge\bigl(c^{2/3}-C_*^{2/3}-o(1)\bigr)N_{wa}.           \tag{FC11}
```

Thus a single finite witness proves a linear landing gap for an infinite
cofinal cone of the algebra. In this project `C_*=1/2` is available from the
square-field Paley construction.

## 3. Near-native restrictions do not remove the obstruction

Let `T` be any signing of order `N` with witness energy at least
`cN^(3/2)`, and delete an arbitrary set of `d` vertices. Restrict the same
witness. At most `dN` edges were incident with the deleted set, so

```math
|H_{T[K]}(z_K)|\ge cN^{3/2}-dN.                        \tag{FC12}
```

Consequently, if `d=o(sqrt(N))`, the restricted signing retains normalized
cap at least `c-o(1)`. Hence an order-filling rule whose native-order gaps
are `o(sqrt(N))` inherits every bad descendant cone. Larger deletions might
change the constant, but then require a new restriction theorem; algebraic
closure alone supplies no summable landing estimate.

## 4. Exact application to the order-14 seed

The report `phase2c_balanced_hadamard_lift_obstruction.md` takes the exact
order-14 conference minimizer `A`, a balanced diagonal `D`, `Q=A+D`, and the
single Sylvester generator `G=H_4`. The explicit word of length one has
order 56 and Boolean witness energy 220. Therefore

```math
c={220\over56^{3/2}}=0.524977439470833\ldots>\frac12. \tag{FC13}
```

Theorem (FC7) proves the obstruction simultaneously for every later
Sylvester continuation, while (FC11) gives the certified linear `b`-gap.
The verifier is
`computations/certify_balanced_hadamard_lift_obstruction.py`.

## 5. Scope and family-selection rule

The theorem is stronger than an isolated tensor example in three ways:

1. the seed witness can already be entangled across several fiber levels;
2. continuations may be arbitrary words in several fixed generators;
3. the conclusion rejects the entire descendant cone without computing any
   later cap.

It applies equally when the generators are presented as character tables or
rank-one Bose--Mesner idempotent channels: the only required mapping is the
exact full-sign Kronecker law (FC4) and the Boolean extremal character (FC2).
Calling the same matrices an association scheme or a code design does not
evade the theorem.

The theorem does **not** directly reject fiber algebras with no Boolean
extremal channel, non-Kronecker fusion rules, or a scale-dependent generator
whose state complexity grows. For symmetric Hadamard generators, however,
the first apparent escape is illusory. If `G` is symmetric Hadamard of order
`k`, then `vec(G)` is Boolean and the standard vectorization identity gives

```math
(G\mathbin\otimes G)\operatorname{vec}(G)
=\operatorname{vec}(G^3)
=k\operatorname{vec}(G).                            \tag{FC14}
```

Thus `G tensor G` always has a Boolean extremal eigenvector, even when `G`
itself does not. Grouping continuation letters in pairs brings every fixed
symmetric-Hadamard alphabet back under the persistence theorem. A missing
Boolean eigenvector can delay the channel by one level but cannot remove it.

The genuine remaining escapes are therefore a non-Kronecker fusion law, a
scale-dependent growing generator/state, or a correction that destroys the
persisted channel quantitatively. For any new fixed algebra satisfying
(FC1)--(FC4), the correct workflow is now
finite and falsifiable:

1. enumerate or optimize only modest prefix words;
2. stop immediately if one explicit witness crosses the best all-order
   constant `C_*`;
3. only if every tested prefix stays below `C_*`, attempt an entangled-spin
   upper theorem and an order-filling theorem.

This is a theorem-level obstruction to fixed fiber algebra as a generic
landing mechanism. It explains why preservation of separable seed energies,
spectral identities, and design closure are insufficient: a finite entangled
state, once found, is permanently carried by the Boolean extremal channel.
