# A two-port orientation holonomy for the BCX anti-pin metric

**Status.** Rigorous task-local falsifier.  This note shows that the BCX
one-port contextual metric is not an approximate congruence for general
dense exact-sign composition.  Two children can be `o(n^(3/2))`-close under
every one-port anti-pin query, while a fixed exact-sign bridge to a second
copy separates their composed caps by `Theta(n^(3/2))`.

This is not another planted-witness failure.  No optimizer is prescribed.
The separation is between exact global Boolean caps and follows from an
exact block spectral identity.

## 1. The orientation bit is invisible to one BCX port

Let `n=q^2`, `q=2^m`, and let `mathcal H` be the symmetric regular Walsh
matrix from BCX.0--BCX.1:

```math
\mathcal H^2=nI,
\qquad
\mathcal H\mathbf1=q\mathbf1,
\qquad
\operatorname{tr}\mathcal H=0.                    \tag{BH.1}
```

Put

```math
A=\mathcal H-\operatorname{diag}(\mathcal H).      \tag{BH.2}
```

Then `A` and `-A` are hollow complete signings, and on Boolean spins

```math
H_A(x)={1\over2}x^T\mathcal Hx,
\qquad
H_{-A}(x)=-{1\over2}x^T\mathcal Hx.                \tag{BH.3}
```

For an arbitrary query `t in {+-1}^n`, append `q` spins and use the BCX
anti-pin

```math
L_t(x,y)=(t\mathbin\cdot x)(\mathbf1\mathbin\cdot y),
\qquad
C=J_q-I_q.                                         \tag{BH.4}
```

Write

```math
F_\sigma(t)
=\max_{x,y}|\sigma H_A(x)+L_t(x,y)+H_C(y)|,
\qquad \sigma\in\{+-1\}.                           \tag{BH.5}
```

### Lemma BH.1 (uniform one-port orientation blindness)

For every query `t`,

```math
|F_+(t)-F_-(t)|
\le2Q(C)=q(q-1)<n.                                 \tag{BH.6}
```

Consequently the sup and projective distances between the two full one-port
response profiles are both `O(n)=o(n^(3/2))`.

#### Proof

First omit `H_C` and call the resulting cap `F_sigma^0(t)`.  The substitution
`y -> -y` gives

```math
\begin{aligned}
F_-^0(t)
&=\max_{x,y}|-H_A(x)+L_t(x,y)|\\
&=\max_{x,y}|-H_A(x)-L_t(x,y)|\\
&=\max_{x,y}|H_A(x)+L_t(x,y)|
=F_+^0(t).                                         \tag{BH.7}
\end{aligned}
```

The absolute-cap functional is one-Lipschitz in pointwise Hamiltonian norm,
so

```math
|F_\sigma(t)-F_\sigma^0(t)|\le Q(C)={q\choose2}.  \tag{BH.8}
```

The triangle inequality proves (BH.6).  A coordinatewise bound by
`2Q(C)` also bounds half the oscillation by `2Q(C)`. `square`

The proof did not use the Rayleigh code or any optimizer statement.  It
holds for every rank-one BCX query, including every query in any declared
code.  Without the exact-sign clique completion the two profiles are
identical, not merely close.

## 2. A fixed dense bridge exposes the forgotten orientation

Take a second order-`n` copy with child `A`.  Between the two shores put the
complete exact-sign bridge `mathcal H`.  The two possible orientation states
have energies

```math
\begin{aligned}
P_+(x,z)
 &=H_A(x)+H_A(z)+x^T\mathcal Hz,\\
P_-(x,z)
 &=H_A(x)-H_A(z)+x^T\mathcal Hz.                   \tag{BH.9}
\end{aligned}
```

Every off-diagonal coefficient in both order-`2n` parents is exactly a sign:
the two old blocks are `A` and `+-A`, while the bridge block is
`mathcal H`.

### Theorem BH.2 (two-port orientation holonomy)

The exact Boolean caps obey

```math
Q(P_+)=2qn,
\qquad
Q(P_-)\le\sqrt2\,qn.                               \tag{BH.10}
```

Therefore

```math
Q(P_+)-Q(P_-)
\ge(2-\sqrt2)n^{3/2}.                              \tag{BH.11}
```

In total-order units `N=2n`, this is the fixed normalized gap

```math
{Q(P_+)-Q(P_-)\over N^{3/2}}
\ge {2-\sqrt2\over2^{3/2}}>0.                     \tag{BH.12}
```

Both parents themselves have bounded cap `O(N^(3/2))`.

#### Proof

Because the diagonal trace in (BH.1) vanishes, the Boolean energies in
(BH.9) equal one half of the quadratic forms of the real symmetric block
matrices

```math
M_+=
\begin{pmatrix}\mathcal H&\mathcal H\\
                \mathcal H&\mathcal H\end{pmatrix},
\qquad
M_-=
\begin{pmatrix}\mathcal H&\mathcal H\\
                \mathcal H&-\mathcal H\end{pmatrix}.          \tag{BH.13}
```

Write

```math
T_+=\begin{pmatrix}1&1\\1&1\end{pmatrix},
\qquad
T_-=\begin{pmatrix}1&1\\1&-1\end{pmatrix}.       \tag{BH.14}
```

Then `M_+=T_+ tensor mathcal H` and
`M_-=T_- tensor mathcal H`.  Since

```math
||T_+||_{op}=2,
\qquad T_-^2=2I,
\qquad ||\mathcal H||_{op}=q,                     \tag{BH.15}
```

their operator norms are `2q` and `sqrt(2)q`.  Every Boolean vector on the
two shores has squared Euclidean norm `2n`, so

```math
Q(P_+)\le2qn,
\qquad Q(P_-)\le\sqrt2\,qn.                        \tag{BH.16}
```

At `x=z=mathbf1`, regularity in (BH.1) gives

```math
P_+(\mathbf1,\mathbf1)
={1\over2}qn+{1\over2}qn+qn=2qn.                  \tag{BH.17}
```

Thus the first upper bound is attained.  Equations (BH.11)--(BH.12) follow,
and (BH.16) proves the bounded-cap claim. `square`

The theorem does not need the exact value of `Q(P_-)`; its spectral upper
bound already yields a fixed cap gap.

The construction is gauge-covariant.  For every BCX switch `s`, put

```math
A_s=D_sAD_s,
\qquad \mathcal H_s=D_s\mathcal HD_s.              \tag{BH.17a}
```

The doubled family `{+-A_s:s in mathcal S}` has the same one-port blindness
pairwise between `A_s` and `-A_s`, while the continuation consisting of one
fixed `A_s` shore and bridge `mathcal H_s` exposes that orientation by the
same gap (BH.11).  Simultaneous gauge conjugation reduces the calculation to
BH.2.  Thus the falsifier persists over an exponential BCX switching code;
it is not tied to one isolated matrix.

## 3. Why this is composition-created holonomy

At one port, outer absolute value and inversion of the auxiliary shore erase
the sign of the child energy.  The exact-sign clique remembers it only at
the lower `O(n)` scale.  Thus the state

```math
[A]_{\rm one\ port}\simeq[-A]_{\rm one\ port}      \tag{BH.18}
```

is a valid vanishing-distortion quotient at the `n^(3/2)` scale.

After two systems are joined, the cross term fixes a **relative** orientation
between their quadratic channels.  The block coefficient is `T_+` or `T_-`,
whose norms differ by the fixed factor `2/sqrt(2)`.  This relative sign is
not a property of either marginal response profile.  It is created as an
observable compatibility bit at composition.

Equivalently, consider the two two-component states

```math
(A,A)
\quad\hbox{and}\quad
(A,-A).                                             \tag{BH.19}
```

Their first components agree, and their second components have vanishing
one-port BCX distance by BH.1.  The same public bridge `mathcal H` separates
the composed caps by BH.11.  This is true holonomy: no planted optimizer,
query-owned joint coefficient, or post hoc choice of bridge is used.

Equivalently again, fix the first `A` block and `mathcal H` bridge as one
order-`n` continuation acting on the second block.  That single continuation
distinguishes the two one-port-close inputs `A` and `-A`; the two-component
language is not essential to the logical counterexample.

## 4. A positive boundary for narrow staged continuations

The failure requires a second macroscopic shore.  It cannot be reproduced
by merely stacking a bounded number of BCX-width auxiliary ports onto one
child.

### Proposition BH.3 (orientation congruence below the internal-edge scale)

Let `A` be any hollow child signing on `n` spins.  Append `m` spins with an
arbitrary old--new matrix `B` and an arbitrary hollow auxiliary signing `C`:

```math
R_\sigma
=\max_{x,y}|\sigma H_A(x)+x^TBy+H_C(y)|.          \tag{BH.20}
```

Then

```math
|R_+-R_-|\le2Q(C)\le m(m-1).                      \tag{BH.21}
```

Hence `A` and `-A` remain congruent at the `n^(3/2)` scale under every such
continuation whenever `m=o(n^(3/4))`.  In particular any fixed collection
of BCX ports that can be flattened into one quadratic continuation of total
appended width `O(sqrt n)` cannot expose the orientation bit at leading
scale.  This statement does not cover adaptive re-encoding or unbounded
depth.

#### Proof

With `C` omitted, invert all auxiliary spins:

```math
\max_{x,y}|-H_A(x)+x^TBy|
=\max_{x,y}|H_A(x)+x^TBy|.                         \tag{BH.22}
```

Adding `H_C` changes either cap by at most `Q(C)`.  The triangle inequality
proves the first part of (BH.21), and every complete signing on `m` vertices
has `Q(C)<=binom(m,2)`. `square`

Thus BH.2 and BH.3 identify an actual transition.  Orientation is a reusable
quotient for sub-`n^(3/4)` auxiliary continuations, but not for composition
with a second order-`n` component carrying a spectrally flat dense bridge.

## 5. Scope and consequence

1. **What is falsified.**  Vanishing one-port BCX response distance is not a
   congruence for arbitrary dense block composition.  A theory retaining
   only that metric loses a relative orientation bit whose value is
   macroscopic after gluing.
2. **What survives.**  Proposition BH.3 is a positive congruence theorem for
   narrow disjoint continuations.  It includes arbitrary query-dependent
   old--new signs and arbitrary interaction among all appended spins.
3. **What is not claimed.**  The example does not concern the selected BCX
   Rayleigh code alone: `-A` need not lie in that switching subcode.  It
   concerns reuse of the metric as a state abstraction on the natural
   bounded-cap exact-sign class containing both orientations.
4. **No relevance overclaim.**  The two-port bridge is a structured regular
   Hadamard matrix.  The theorem proves failure of general congruence, not
   incompressibility of arbitrary dense bridges or a statement about exact
   minimizers.

The minimal repair to the state is one relative orientation bit for this
example.  More generally, a reusable carrier must retain how absolute-value
channels are coherently oriented before multiple systems share a bridge;
marginal cap profiles alone cannot reconstruct that compatibility.
