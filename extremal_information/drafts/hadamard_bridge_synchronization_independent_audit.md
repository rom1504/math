# Independent audit: Hadamard bridge synchronization

**Audited files**

- `drafts/hadamard_bridge_synchronization.md`;
- `experiments/verify_hadamard_bridge_synchronization.py`.

**Verdict:** PASS, with two presentation clarifications and one important
scope clarification.  HS.1--HS.4, including every displayed leading
constant, are correct under the symmetric Sylvester--Walsh convention used
by the verifier.  Hollowing the auxiliary Walsh block causes no hidden
diagonal error because its trace is zero.  HS.4 is a planted-witness
stability obstruction; it is not a proof that every global parent optimizer
lies off the full Hadamard relation set.

## 1. HS.1: pair and eliminated deficits

Put `z=W^Tx` and `v=z/sqrt(k)`.  Orthogonality gives

```math
||v||_2^2={1\over k}x^TWW^Tx=k,                    \tag{A.1}
```

and `||y||_2^2=k`.  Since
`x^TWy=sqrt(k) v dot y`, one has exactly

```math
{\sqrt k\over2}||y-v||_2^2
=\sqrt k(k-v\mathbin\cdot y)
=k^{3/2}-x^TWy.                                    \tag{A.2}
```

Coordinatewise maximization gives `max_y z dot y=||z||_1`.  Also

```math
\sum_a(|z_a|-\sqrt k)^2
=||z||_2^2+k^2-2\sqrt k||z||_1
=2k^2-2\sqrt k||z||_1,                             \tag{A.3}
```

so division by `2sqrt(k)` gives the middle identity in HS.4.  Finally

```math
\operatorname{dist}(v,\{+-1\}^k)^2
={1\over k}\sum_a(|z_a|-\sqrt k)^2,                \tag{A.4}
```

which gives its last identity.  The uniqueness statement is correct:
coordinatewise signs are unique exactly when no `z_a` vanishes; on the
Boolean regularizer code every coordinate equals `+-sqrt(k)`.

No quadratic-energy factor of two is missing here: the bridge occurs once
as `x^TWy` in the hollow parent energy.

## 2. HS.2: Lipschitz extension

Equation HS.5 is just (A.2) with `K(y)` added.  If `r=||y-v_x||`, then

```math
K(y)-{\sqrt k\over2}r^2
\le K(v_x)+Lr-{\sqrt k\over2}r^2
\le K(v_x)+{L^2\over2\sqrt k}.                    \tag{A.5}
```

For `K(v)=v^TCv/2` on the radius-`sqrt(k)` sphere,

```math
|K(y)-K(v)|
={1\over2}|(y-v)^TC(y+v)|
\le\sqrt k||C||_{op}||y-v||,                      \tag{A.6}
```

so HS.7 is exact.  The Frobenius observation is also correct for a hollow
dense signing:

```math
||C||_{op}\ge {||C||_F\over\sqrt k}=\sqrt{k-1}.    \tag{A.7}
```

### Presentation clarification 1

The sentence “the reverse inequality with zero error also holds” should be
read, or rewritten, as

```math
\max_y\{x^TWy+K(y)\}\ge k^{3/2}+K(v_x)
```

when `v_x` is Boolean.  It does not assert equality, because another
Boolean `y` could improve `K` enough to overcome part of the bridge penalty.

## 3. HS.3: Walsh convention and exact departure

Use the symmetric Walsh matrix

```math
W_{(a,b),(u,v)}=(-1)^{a\cdot u+b\cdot v}.           \tag{A.8}
```

Then

```math
\sum_{u,v}(-1)^{a\cdot u+b\cdot v+u\cdot v}
=q(-1)^{a\cdot b},                                 \tag{A.9}
```

so `Wx_0=qx_0`, exactly as claimed.  This convention is symmetric, hence
the draft's switch from `W^Tx` in HS.2 to `Wx` in HS.12 is harmless.

Writing `x_S=x_0-2x_0 1_S` gives

```math
{(Wx_S)_a\over q}
=x_0(a)\left(1-{2c_a\over q}\right).              \tag{A.10}
```

The assumption `d<q/2` makes every parenthesis strictly positive, so the
unique bridge sign remains `x_0`.  Symmetry and (A.9) give

```math
\sum_a c_a
=\sum_{i\in S}x_0(i)(W^Tx_0)_i=dq.                \tag{A.11}
```

Consequently both the `ell_1` optimum and the displayed planted cross term
equal `kq-2dq`, and the deficit is exactly `2dq=2d sqrt(k)`.  The positive
deficit for `d>0` indeed rules out Booleanity of `Wx_S/q`.

### Presentation clarification 2

The draft should state explicitly once that `W` is the symmetric
Sylvester--Walsh matrix.  For a nonsymmetric Hadamard matrix, HS.10 would
not justify the transpose step in HS.13 without a separate dual bent-vector
identity.

## 4. HS.4: prescribed cut and old-block cap

The construction of `A` is consistent.  Prescribe all `d(k-d)` crossing
edges by

```math
A_{ij}=-x_0(i)x_0(j),                               \tag{A.12}
```

then independently choose complete hollow filler signings on `S` and
`S^c`.  These three disjoint edge sets exhaust the old complete graph.
Flipping `S` leaves both fillers unchanged and changes every crossing
contribution from `-1` to `+1`, so

```math
H_A(x_S)-H_A(x_0)=2d(k-d).                          \tag{A.13}
```

For any old spin, the three parts separately give

```math
Q(A)
\le d(k-d)+2d^{3/2}+2(k-d)^{3/2}.                  \tag{A.14}
```

At `d=floor(q/4)` and `k=q^2`, the right side is below `3k^(3/2)` for all
sufficiently large `q`.  Thus HS.21 is sound; no compatibility between the
two probabilistic fillers is required.

The filler fact HS.14 has ample constants.  For one fixed spin, Hoeffding at
`2n^(3/2)` gives at most `2exp(-4n^2/(n-1))`, and union over `2^n` spins is
below one.  Switching or correlations are not used.

## 5. The hollow Walsh block and its diagonal

Let `W=H_2^{\otimes 2m}` in the convention of (A.8), and put

```math
C=W-\operatorname{diag}(W).                        \tag{A.15}
```

This is a hollow complete signing.  Since `tr H_2=0`,

```math
\operatorname{tr}W=(\operatorname{tr}H_2)^{2m}=0. \tag{A.16}
```

For every Boolean `y`, therefore,

```math
H_C(y)
={1\over2}y^TCy
={1\over2}y^TWy-{1\over2}\sum_iW_{ii}y_i^2
={1\over2}y^TWy.                                   \tag{A.17}
```

The spectral upper bound is

```math
|H_C(y)|\le {1\over2}\sqrt k\,||y||_2^2
={1\over2}k^{3/2},                                 \tag{A.18}
```

and equality is attained at `y=x_0` because `Wx_0=qx_0`.  Hence the exact
equality `Q(C)=k^(3/2)/2` is valid.

It is important that (A.17) is a Boolean identity.  Removing the diagonal
does change the quadratic extension away from the Boolean sphere, but HS.4
never invokes HS.2 with this `C`, so no continuous-extension error is hidden.

## 6. Parent normalization and gain

Because the Walsh matrix is symmetric, the block matrix in HS.16 is a
valid symmetric hollow signing; for a convention-independent display its
lower-left block could be written `W^T`.  Its Boolean energy is

```math
H_P(x,y)=H_A(x)+x^TWy+H_C(y).                       \tag{A.19}
```

The cross term obeys

```math
|x^TWy|\le||W||\,||x||_2||y||_2=k^{3/2}.           \tag{A.20}
```

Thus

```math
Q(P)\le3k^{3/2}+k^{3/2}+{1\over2}k^{3/2}
={9\over2}k^{3/2}.                                 \tag{A.21}
```

At the two compared states the `C` term cancels.  Combining (A.13) with
HS.11 gives

```math
H_P(x_S,x_0)-H_P(x_0,x_0)
=2d(k-d-q).                                        \tag{A.22}
```

For `q>=8`, `d=floor(q/4)>=q/8` and
`k-d-q>=k/2`; hence (A.22) is at least
`qk/8=k^(3/2)/8`.  Every constant in HS.15--HS.18 passes.

In total-order units `N=2k`, (A.21) is

```math
Q(P)\le {9\over2^{5/2}}N^{3/2},                   \tag{A.23}
```

so the claimed bounded-cap scope is genuine.

## 7. Scope correction

HS.4 proves that the **planted** exact relation state `(x_0,x_0)` is not
stable under arbitrary bounded-cap exact internal energies: the explicitly
constructed off-relation state beats it by a fixed `k^(3/2)` amount.

It does not prove either of the stronger statements

```math
\max_{x,y}H_P(x,y)>
\max_{x:\,Wx/\sqrt k\in\{+-1\}^k}
 H_P(x,Wx/\sqrt k),                                \tag{A.24}
```

or “every global optimizer is off relation.”  Another regularizer-code state
could still beat `(x_S,x_0)`.  The final scope paragraph of the draft already
acknowledges this, but the sentence immediately after HS.18 (“do not give
even one-layer robust synchronization”) should be tightened to “do not give
one-layer stability of a planted Boolean pullback witness.”  With that
wording the conclusion is fully supported.

This scope is still a scalable and useful no-go: it falsifies any proof step
which assumes that the bare Hadamard penalty preserves a selected pullback
witness uniformly over all cap-scale old landscapes.

## 8. Verifier audit

The supplied verifier runs successfully and checks:

- HS.3 and both forms of HS.4 on random Boolean pairs at orders
  `k=4,16,64`;
- the self-dual bent convention and exact optimizer sign after departures;
- the old-block cut gain and bridge loss;
- the finite old-block cap for the enumerably feasible case `k=16`.

It does not check:

- the abstract Lipschitz statement HS.2;
- `tr W=0`, the hollow block `C`, or `Q(C)=k^(3/2)/2`;
- the complete parent cap.

These omissions do not undermine the analytic proof, but adding the
following cheap assertions for `m=1,2` would improve normalization regression
coverage:

```python
assert sum(W[i][i] for i in range(k)) == 0
C = [[0 if i == j else W[i][j] for j in range(k)] for i in range(k)]
assert cap(C) == k * math.isqrt(k) // 2
```

The full parent cap need not be exhaustively enumerated; (A.21) is a direct
triangle inequality.

## 9. Final classification

| Item | Verdict | Action |
|---|---|---|
| HS.1 identities | PASS | none |
| HS.2 Lipschitz/error constant | PASS | clarify the lower-bound wording |
| HS.3 Walsh transform | PASS | state symmetric convention |
| HS.3 deficit and optimizer | PASS | none |
| HS.4 prescribed-cut construction | PASS | none |
| `Q(A),Q(C),Q(P)` constants | PASS | none |
| trace/diagonal handling | PASS | optionally add verifier assertions |
| synchronization conclusion | PASS after scope wording | say planted-witness stability |

The draft supplies a rigorous architecture-specific falsifier, not a full
Hadamard-relation no-go theorem.  No algebraic or normalization defect was
found.
