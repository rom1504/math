# Independent audit of the multiscale partition affine shell

Date: 2026-08-17.

Audited draft: `multiscale_partition_affine_shell.md`.

## 0. Verdict

**PASS, with scope qualifications.**  MP.1--MP.8 are correct with the
current constant `8`.  In fact, the crucial conclusion is genuinely
one-sided: after one global orientation, every member of the affine cube
has energy at least `(1-8/q)Q`.  No connectivity or bounded one-spin-jump
argument is needed.  The projective cardinality, edge-word packing scale,
shell-budget inversion, and one-sided response normalization all check out.

The finite verifier
`experiments/verify_multiscale_partition_affine_shell.py` exhausts every
hollow signing through order five and checks random signings through order
nine.  For every `q` it checks the two partition budgets, every subset of
every partition block, the selected oriented affine cube, projective
cardinality, every star-frame endpoint, and the full Boolean absolute and
one-sided responses.  Its default run returns

```text
PASS: 1146 signings through order 9
```

| Claim | Verdict | Audit finding |
|---|---|---|
| MP.1--MP.2 partition cap budget | **PASS** | Correct for every hollow symmetric real matrix.  The same partition inequality already appears in `artifacts/nested_restriction_paving.md`, equation (1). |
| MP.3 subset completion | **PASS** | Conditional random completion gives the sharper one-sided bound exactly as stated. |
| MP.4 oriented affine cube | **PASS** | The current direct proof is correct; `8Q/q` has no missing factor of two. |
| Projective size and odd-product closure | **PASS** | There are exactly `2^|I|` projective members, and odd products are mask XOR. |
| MP.6 edge-word packing | **PASS** | A constant-distance mask code gives exact edge distance `d(n-d)=Theta(n^2/q)`. |
| MP.6a shell-budget inversion | **PASS** | The stated range makes `q` admissible and the conservative constant `16` is valid. |
| MP.7 absolute response | **PASS** | Correct for every endpoint and every real `m>=0`. |
| MP.8 one-sided response | **PASS** | This is the strongest part: one common child orientation works for the whole language. |
| State-size interpretation | **PASS WITH SCOPE** | Sublinear for the already declared/gauged star language; a standalone labelled interface additionally needs the ground-state gauge and labelled block. |
| Physical/composition interpretation | **PASS WITH SCOPE** | It rigorously diagnoses the raw-frame and separately-paid compiler scales, but is not an impossibility theorem for joint cancellation or a global construction. |

## 1. Partition and subset budgets

For a hollow matrix, the uniform Boolean-spin mean of its quadratic energy
is zero.  Hence both one-sided caps are nonnegative.  On a partition
`J_1,...,J_q`, choose a positive optimizer within each block and independently
multiply the block spin by a Boolean sign.  The internal contribution is
fixed at `sum_a Q_+(D[J_a])`, and every cross edge has mean zero.  Some
block-sign choice therefore has full energy at least that sum.  The negative
argument is identical.  Thus

```math
\sum_a Q_+(D[J_a])\le Q_+(D)\le Q(D),
\qquad
\sum_a Q_-(D[J_a])\le Q_-(D)\le Q(D).             \tag{A.1}
```

This proves MP.1 and MP.2.  The slightly sharper right sides in (A.1) are
occasionally useful.  As a novelty matter, this exact fact is already the
first theorem of `artifacts/nested_restriction_paving.md`; it should not be
advertised as new by itself.

For `S subseteq J`, hold the spins in `S` at `+1` and make all other spins
in `J` independent uniform signs.  Every edge except those internal to `S`
averages to zero.  If

```math
P_S=\sum_{\{i,j\}\subseteq S}d_{ij},
```

then `P_S` is an average of energies in
`[-Q_-(D[J]),Q_+(D[J])]`.  Consequently

```math
-Q_-(D[J])\le P_S\le Q_+(D[J]).                  \tag{A.2}
```

This proves both assertions of MP.3.  It is also sharper than the earlier
polarization estimate `|P_S|<=3Q(D[J])/2`.

## 2. The oriented cube and its constant

Choose `rho,x` with `rho H_A(x)=Q(A)=Q` and switch to

```math
D=rho\,diag(x)A\,diag(x).
```

Then `1` is a positive optimizer of `D`.  Its row fields obey

```math
ell_i=\sum_jd_{ij}\ge0,
\qquad
\sum_iell_i=2Q.                                  \tag{A.3}
```

The first assertion follows by flipping one coordinate and the second from
`1^TD1=2Q`.  On a balanced `q`-block partition, put

```math
L_a=\sum_{i\in J_a}ell_i,
\qquad R_a=Q_-(D[J_a]).
```

Equations (A.1) and (A.3) give

```math
\sum_a(2L_a+4R_a)\le4Q+4Q=8Q.                   \tag{A.4}
```

Select a block `I` no larger than this average.  For every `S subseteq I`,
the exact flip identity is

```math
H_D(1^S)=Q-2\sum_{i\in S}ell_i+4P_S.             \tag{A.5}
```

Using nonnegativity of the fields and the *lower* side of (A.2),

```math
Q-H_D(1^S)
 =2\sum_{i\in S}ell_i-4P_S
 \le2L_a+4R_a
 \le {8Q\over q}.                               \tag{A.6}
```

This proves the oriented statement, not merely an absolute-energy
statement.  In particular, common orientation when `q>8` follows directly.
An earlier possible proof through cube connectivity and a `2(n-1)`
one-spin energy jump would have imposed an unnecessary extra hypothesis.

There is a harmless refinement.  If

```math
R=Q_-(D)
```

is the cap opposite the selected absolute ground-state side, then (A.4)
has right side `4Q+4R`.  Therefore the exact improved defect is

```math
{4(Q+R)\over q}\le {8Q\over q}.                 \tag{A.7}
```

The cleaner symmetric constant `8` is nevertheless valid and appropriate
for the principal statement.

Nothing in this derivation except later physical terminology uses
`|a_ij|=1`.  The affine-shell theorem and response formulas extend verbatim
to hollow symmetric real matrices.  The signing hypothesis is needed for
the intended cut-word and exact-sign-interface interpretation, not for
(A.6).

### Projective size and algebra

Every selected block has size between `floor(n/q)` and `ceil(n/q)`.  Since
`q>=2`, it omits a vertex.  Thus

```math
x^S=+-x^T
```

for `S,T subseteq I` can occur only when `S=T`; the alternative
`S triangle T=[n]` is impossible.  The projective cardinality is exactly
`2^|I|`.  For an odd number of masks,

```math
x^{S_1}\odot\cdots\odot x^{S_(2r+1)}
=x^{S_1\mathbin\triangle\cdots\mathbin\triangle S_(2r+1)},
```

so closure is actual, not merely projective.

## 3. Shell entropy and packing

Write `k=|I|`.  For `q<=n/2`, one has `k=Theta(n/q)`.  A standard greedy
packing of the `k`-cube supplies `exp(Omega(k))` masks with mutual mask
distance `d` between `delta k` and `k`, for a fixed `delta>0`.  MP.4 already
puts every mask in one energy orientation once `q>8`; no orientation
pigeonhole is needed.

Two spins whose masks differ in `d` coordinates have cut words differing
on exactly

```math
d(n-d)                                                   \tag{A.8}
```

edges.  Since `q->infinity` gives `d<=k=o(n)`, (A.8) is
`Theta(n^2/q)`.  Thus MP.6 is correct.  If one wants the number of words to
diverge, rather than merely the literal bound `exp(Omega(n/q))`, one should
also say `n/q->infinity`; the displayed example `q=log^2 n` has this
property.

For MP.6a, set `r=8Q/Delta`.  The assumed range gives `4<=r<=n/2`, so
`q=ceil(r)` is an admissible integer and `8Q/q<=Delta`.  Since `q<=2r`,

```math
\left\lfloor{n\over q}\right\rfloor
\ge\left\lfloor{nDelta\over16Q}\right\rfloor.
```

The stated inversion follows.  A mildly sharper but less memorable version
uses `q<=r+1` and returns

```math
\left\lfloor{nDelta\over8Q+Delta}\right\rfloor
```

coordinates.  No strengthening is required for correctness.

## 4. Response formulas

Let `k'=|I'|` be even, `p=k'+1`, and

```math
W=(x,(x^{\{i\}})_(i\in I')).
```

For `epsilon in {+-1}^p`, put `t=sum_j epsilon_j`.  The coordinate of
`W epsilon` is `x_i t` off `I'` and `x_i(t-2epsilon_i)` at the exceptional
row indexed by `i`.  Both values are odd and hence nonzero.  After one
global projective sign, `sgn(W epsilon)` is therefore `x^S` for some
`S subseteq I'`.

For the absolute response, cap plus Holder gives

```math
mathcal B_A(mW epsilon)\le Q+m||W epsilon||_1.
```

The selector pays the entire field norm and has absolute energy at least
`Q-8Q/q`, proving MP.7.  More strongly, (A.6) says that this same selector
has energy at least `Q-8Q/q` for `rho A`, independent of the endpoint.
This proves MP.8 with no hidden quadratic sign and with the displayed
normalization.  The argument works for every real `m>=0`, including `m=0`.

The state-size sentence needs its stated convention.  Conditional on the
already declared abstract star frame, the numerical norm is

```math
||W epsilon||_1
=(n-k')|t|+\sum_{i\in I'}|t-2epsilon_i|,          \tag{A.9}
```

so the response grammar itself needs only `(n,k',Q)` and is even smaller
than the safe `O(k' log n)` sparse-histogram count.  A standalone labelled
description of the physical frame additionally needs the projective
ground-state gauge (`n-1` bits) and the labelled set `I'`
(`O(k' log n)` bits).  Computing that gauge may require solving the child
ground-state problem.  This is a construction/preprocessing issue, not a
reconstruction of the full `2^n` energy landscape, but it must remain
visible in any algorithmic or gauge-free interpretation.

## 5. Physical and compositional scope

For the all-positive endpoint of the star frame,

```math
||W1||_1=n(k'+1)-2k'=Theta(n^2/q)
```

when `k'=Theta(n/q)`.  Thus the raw exact-sign frame has a macroscopically
unbalanced endpoint.  If `k'>=2`, an endpoint with `t=3` has exact norm
`3n-4` and selects `x`; when `k'=0`, the one ground-state port has norm
`n`.  These calculations are correct.

The cited microcanonical residual

```math
O(sqrt(nk'(n+k')))=O(n^(3/2)/sqrt(q))
```

is the separately proved bilinear compiler scale in
`local_affine_interface_composition.md` (and its independent audit), not a
new consequence of MP.3 alone.  Paying that residual separately at
`Theta(q)` steps indeed gives no summable recurrence.  This is a valid
ceiling for the raw-frame/scalar-compiler implementation.  It does not rule
out a construction that cancels residual channels jointly, chooses all
levels globally, or proves a new cross-level congruence; the draft correctly
leaves those possibilities open.

## 6. Repository novelty and research judgment

Three ingredients predate this draft:

1. the partition cap budget is equation (1) of
   `artifacts/nested_restriction_paving.md`;
2. the oriented local-field identities already occur in
   `nearmin_deterministic_inequalities.md`;
3. odd-product star frames and selector/trust-response bounds already occur
   in the PC/LA response algebra.

The new theorem is their nontrivial combination.  The one-sided block-cap
budget removes the `2k^2` internal-edge term in the earlier low-field cube.
For `k=n/q`, the earlier estimate was

```math
O(kQ/n+k^2),
```

whereas MP.4 gives `O(kQ/n)` even when `k` is nearly linear.  In the
important regime `q=o(sqrt n)`, this is a genuine asymptotic strengthening,
not a constant improvement.  It yields an orientation-pure affine cube of
dimension `Theta(n/log^2 n)` at normalized deficit `O(1/log^2 n)` for every
bounded-cap signing.  No earlier repository theorem located in this audit
has that conclusion.

Accordingly, MP.3/MP.8 are theorem-level progress and a strong universal
benchmark.  They are not near-minimizer-specific rigidity and do not yet
produce a physical cross-order state.  Their correct frontier role is:

```text
large orientation-pure vanishing-width affine response languages are
universal; the remaining obstruction is physical balance plus reusable
cross-level congruence, not shell existence or endpoint orientation.
```
