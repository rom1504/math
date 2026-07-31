# Independent audit of the low-scale Hadamard lift

Status: theorem independently verified.  One decimal coefficient in the
informal task statement was incorrect: the exact upper coefficient is
`207/8=25.875`, not `25.951`.

## 1. Lift and exact energy identity

Let `C` be an order-`m` signing, put `Q_C=C+I_m`, and let `H` be a symmetric
Hadamard matrix of order `k`.  Write `D=diag(H)` and define

```math
S_C(H)=C\otimes H+I_m\otimes(H-D)
      =Q_C\otimes H-I_m\otimes D.                    \tag{1}
```

Its diagonal is zero.  If two lifted vertices have different macro
coordinates, their edge is `c_ij h_ab`; if they have the same macro
coordinate and different micro coordinates, their edge is `h_ab`.  Every
off-diagonal entry is therefore a sign, so (1) is a signing.

For a Boolean vector `x=(x_1,...,x_m)` with each `x_i` in `{+/-1}^k`,

```math
H_{S_C}(x)
={1\over2}\sum_{i,j}(Q_C)_{ij}x_i^T Hx_j
 -{m\over2}tr(H).                                     \tag{2}
```

The last term follows from `x_(ia)^2=1`; it is not a norm error.  For the
Sylvester matrices used below, `k=4^s`, `H=H_4^(tensor s)`, and
`tr(H)=(tr H_4)^s=0`.

Put `U=H/sqrt(k)`.  Then `U` is a symmetric involution.  On writing
`x_i=sqrt(k)u_i`, with `||u_i||=1`, and decomposing
`u_i=p_i+n_i` into the positive and negative eigenspaces of `U`, (2) becomes

```math
{H_{S_C}(x)\over k^{3/2}}
={1\over2}\sum_{i,j}(Q_C)_{ij}
  (\langle p_i,p_j\rangle-\langle n_i,n_j\rangle).  \tag{3}
```

This checks the factor `1/2`, the sign of the negative eigenspace, and the
power `k^(3/2)`.

## 2. Rational dual certificate for the first base

For the base with code `5850642905`, let `Q=A+I_10` and set

```math
y={1\over200}(445,490,661,668,436,645,405,427,485,513). \tag{4}
```

Exact rational `LDL^T` factorization proves

```math
\operatorname{Diag}(y)-{Q\over2}\succ0,
\qquad \operatorname{Diag}(y)+{Q\over2}\succ0.       \tag{5}
```

The integer matrices actually checked are

```math
200(\operatorname{Diag}(y)\mathbin{\pm}Q/2)
=\operatorname{Diag}(200y)\mathbin{\pm}100Q.          \tag{6}
```

All exact `LDL^T` pivots are positive and are stored in the audit JSON.  This
is also why interpreting the listed numerators as fractions over `100`
would introduce an erroneous factor two.

Apply the minus certificate in (5) to the `p` terms of (3), and the plus
certificate to the negated `n` terms.  Since
`||p_i||^2+||n_i||^2=1`,

```math
H_{S_A}(x)/k^{3/2} <= sum_i y_i.
```

Swapping the two certificates gives the same bound for the negative of the
energy.  Finally,

```math
sum_i y_i={5175\over200}={207\over8},
```

and hence

```math
cap(S_A(H)) <= {207\over8} k^{3/2}.                   \tag{7}
```

This is an absolute-cap bound; it controls both energy orientations.

The independent audit also found a modest exact strengthening while retaining
the essential **common** dual vector.  Put

```math
y'={1\over100000}(222237,244660,330058,333605,217787,
322119,202501,213301,242265,256161).
```

Separate exact `LDL^T` factorizations prove both
`Diag(y')-Q/2` and `Diag(y')+Q/2` positive definite.  Since
`sum_i y'_i=1292347/50000=25.84694`, the same argument gives

```math
cap(S_A(H)) <= {1292347\over50000}k^{3/2}.            \tag{7a}
```

This is valid because the *same* `y'` passes both signed constraints.  Merely
optimizing unrelated vectors for the two signs would not be valid: the two
Hadamard Gram diagonals are complementary rather than separately equal to
one.

## 3. Exact Boolean witness for the second base

For

```text
H4 = [[1, 1, 1, 1],
      [1,-1, 1,-1],
      [1, 1,-1,-1],
      [1,-1,-1, 1]],
```

the Boolean vector `v4=(-1,-1,-1,1)` satisfies `H4 v4=2v4`.  Therefore
`v=v4^(tensor s)` is Boolean and

```math
Hv=sqrt(k)v,
qquad v^THv=k^{3/2}.                                  \tag{8}
```

The second order-ten base `B`, code `28771662001`, has an exactly enumerated
Boolean state `sigma` with `H_B(sigma)=21`.  Since

```math
{1\over2}sigma^T(B+I)sigma=21+{10\over2}=26,
```

substitution of `x=sigma tensor v` into (2) gives the exact positive witness

```math
cap(S_B(H)) >= 26 k^{3/2}.                            \tag{9}
```

At `k=4`, the independently reconstructed witness energies are 192 for the
first base state and 208 for the second, agreeing with
`(19+5)8` and `(21+5)8`.

Combining (7) and (9), for every `k=4^s`, `s>=1`,

```math
cap(S_B(H))-cap(S_A(H)) >= {1\over8}k^{3/2}.          \tag{10}
```

The lift order is `N=10k`, so the normalized separation is

```math
{1\over 8\,10^{3/2}}N^{3/2}
=0.003952847075... N^{3/2}.                           \tag{11}
```

The refined common certificate (7a) strengthens the exact gap to

```math
{7653\over50000}k^{3/2}
=0.004840182186...N^{3/2}.                            \tag{11a}
```

Both sequences are genuinely on the `N^(3/2)` scale.  For either fixed base,

```math
||S_C(H)||_op <= sqrt(k)||C+I||_op+1,
```

so the ordinary spectral estimate gives `cap(S_C)=O(N^(3/2))`.  Thus (10)
is not the earlier quadratic-cap twin obstruction: it falsifies fixed-profile
stability inside an `O(N^(3/2))` class.

## 4. Profile preservation

The required hypothesis is equality of the **oriented** switching/permutation
restriction histograms of the bases through order six.  Equality only after
global negation is not enough.  The two bases satisfy the stronger oriented
equality exactly at every order from two through six.

Here is a direct bijection proof.  Select `ell<=6` vertices of a lift.  Let
`I` be their set of distinct macro coordinates and, for each `i in I`, let
`R_i` be the nonempty subset of selected micro coordinates in that fiber.
For distinct selected vertices `(i,a),(j,b)`, including the case `i=j`, the
edge sign is

```math
(Q_C)_{ij}h_{ab}.                                     \tag{12}
```

Pair a base support `I` for `A` with a base support `J` for `B` in the same
oriented class.  Choose a witnessing macro permutation `pi:I->J` and
switches `s_i`, so

```math
(Q_B)_{pi(i),pi(j)}=s_i s_j(Q_A)_{ij}.                \tag{13}
```

Equation (13) also holds for `i=j`, because both diagonal entries are one.
Transport the micro assignment by `R_(pi(i))=R_i`, permute `(i,a)` to
`(pi(i),a)`, and switch every selected lifted vertex in fiber `i` by `s_i`.
The Hadamard factor in (12) is unchanged and (13) maps every edge exactly.

For each oriented base class, the numbers of supports agree.  All tuples of
micro subsets are transported bijectively, including tuples in which the
same micro coordinate occurs in different fibers.  Summing over support
sizes at most `ell` proves equality of the complete oriented lift profiles
through order six, and hence equality of `phi_6`, for every `k`.  No spectral
property of `H` is needed for this profile argument; symmetry is needed only
for the result to be an undirected signing.

## 5. Can the obstruction be moved toward normalized `1/2`?

The present normalized bracket is

```math
{207/8\over10^{3/2}}=0.8182393445...,
\qquad {26\over10^{3/2}}=0.8221921916....             \tag{14}
```

The refined common certificate improves the first displayed constant to
`(1292347/50000)/10^(3/2)=0.8173520094...`; it remains far above `1/2`.

Increasing the Hadamard order does not change these constants.  The micro
lift amplifies a base separation but does not improve it.  Progress toward
the conference-scale constant `1/2` therefore requires a better **macro
base collision**, not a larger `k`.

For an order-`m` base `C`, define the split-eigenspace elliptope certificate

```math
\Theta(C)=\min\left\{\sum_i y_i:
 \operatorname{Diag}(y)\mathbin{\pm}{C+I\over2}\succeq0\right\}. \tag{15}
```

and its positive Boolean witness value

```math
W(C)={m\over2}+max_sigma H_C(sigma).                  \tag{16}
```

The audited mechanism works for an oriented-profile collision `(A,B)`
whenever

```math
W(B)>\Theta(A).                                        \tag{17}
```

It then gives normalized upper/lower coefficients
`Theta(A)/m^(3/2)` and `W(B)/m^(3/2)`.  This is the correct falsifiable search
objective for this relaxation.  A difference of the base caps alone is
insufficient: the larger cap must cross the other base's split-eigenspace
vector-relaxation certificate.  The common dual vector in (15) is essential:
the positive and negative Hadamard Gram matrices have complementary
diagonals, not separate unit diagonals, so independently optimized signed
duals cannot simply be combined.

To put both lifted families below `(1/2+eta)N^(3/2)`, one additionally needs
certificates

```math
\Theta(A),\Theta(B)
 <= (1/2+eta)m^{3/2},                                 \tag{18}
```

or an equally strong Boolean upper bound.  Spectrally flat conference-like
bases are natural candidates because `||C+I||_op=(1+o(1))sqrt(m)` gives the
right coefficient.  Generic random signings have operator norm near
`2sqrt(m)` and do not meet this spectral surrogate.

There is no formal incompatibility between (17) and (18).  If collisions of
increasing base order satisfy both, with `eta` tending to zero, then for every
fixed accuracy one can choose a base and obtain an infinite equal-profile
pair at that accuracy.  However, no such collision is currently known, and
several cautions matter:

1. exact equality of oriented six-decks is stronger than equality of the
   original unoriented `phi_6`;
2. conference-like spectral flatness is much rarer than an unrestricted
   profile collision;
3. the useful certified gap is `W(B)-Theta(A)`, not
   `cap(B)-cap(A)`; and
4. if both normalized coefficients tend to `1/2`, their normalized gap must
   itself tend to zero.  Each fixed base can still yield a nonzero scalable
   gap, but there is no single uniform gap compatible with both constants
   converging to `1/2`.

The recommended computational extension is therefore to hash larger
spectrally flat or design-derived bases by their complete oriented profile,
then rank colliding pairs by (17), solving and rationalizing the two signed
SDP duals.  Conference switching classes, regular two-graphs, and
Hadamard/design-derived signings are more relevant than unrestricted random
matrices.  A solver should reject a pair before cap optimization whenever its
`Theta`/spectral surrogate is already too large.

This result already conclusively kills fixed `phi_6` as a universal state on
the correct `N^(3/2)` scale.  Moving it toward `1/2` would sharpen the
obstruction against conference-scale structured families, but would still
not by itself prove that an exact-minimizer profile has a bad canonical
representative, because the true asymptotic optimum is only known to lie
between `0.33649...` and `1/2`.

## Reproduction

```bash
.venv/bin/python computations/phase2b_hadamard_lift_theorem_audit.py \
  --output computations/results/phase2b_hadamard_lift_theorem_audit.json \
  > computations/logs/phase2b_hadamard_lift_theorem_audit.log
```

The audit uses a separate exact rational `LDL^T` implementation, reconstructs
both bases from their codes, enumerates their positive Boolean maxima, checks
the `H_4` eigenvector, builds both order-40 lifts, and verifies the witness
energies and every scaling factor.
