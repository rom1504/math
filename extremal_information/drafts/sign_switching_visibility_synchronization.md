# Switching visibility converts pole deficits into relative synchronization

**Status.** Rigorous exact-sign theorem, with an approximate edge-orbit
version.  The theorem identifies the extra hypothesis hidden between
"multiplicity-free symmetry" and a useful Gram-relative certificate: the
declared Boolean pole orbits must see every character with quantitative
weight.  The resulting certificate is sharp on the hollow order-16 PC.3
seed and exact on a non-Hadamard Cayley signing.

## 1. Signed switching actions and pole visibility

Let

```math
A=A^T\in\{0,\pm1\}^{n\times n},\qquad
A_{ii}=0,\qquad |A_{ij}|=1\quad(i\ne j),          \tag{SV.1}
```

and choose `r>=||A||_op`.  Put `T=A/r`.  Let
`Gamma=(Z/2Z)^s` act on `R^n` by commuting signed-permutation involutions
`P_1,...,P_s`.  Thus every `P_i` is orthogonal, `P_i^2=I`, and the `P_i`
commute.  Let `U` be a `Gamma`-invariant subspace on which the action is
multiplicity-free.  Since all characters are real,

```math
U=\bigoplus_{\chi\in X}U_\chi,
\qquad \dim U_\chi=1,                             \tag{SV.2}
```

for a set `X` of distinct characters.  Write `Pi_chi` for the orthogonal
projection onto `U_chi`.

Let `P` be a finite family of Boolean poles `z in {+-1}^n` which spans
`U` and is invariant projectively:

```math
P_gz\in\{\pm z':z'\in\mathcal P\}
\quad(g\in\Gamma,z\in\mathcal P).                \tag{SV.3}
```

Choose one representative `z_j` from every projective orbit, normalize
`u_j=z_j/sqrt(n)`, and define the orbit--character incidence weights

```math
p_{j\chi}=\|\Pi_\chi u_j\|_2^2,
\qquad
\nu=\min_{\chi\in X}\max_jp_{j\chi}.             \tag{SV.4}
```

The spanning assumption implies `nu>0`.  It need not be bounded away from
zero: multiplicity-free action by itself does not give quantitative
synchronization.  The number `nu` is computable from the finite signed
permutation action and the pole orbits; it contains no Boolean
maximization.

For a pole `z`, its individual upper-Rayleigh deficit is

```math
d(z)=1-{z^TAz\over rn}\ge0.                       \tag{SV.5}
```

For a generator, let `k_i` be the number of **unordered** off-diagonal
positions on which switching conjugacy fails:

```math
k_i=\#\{\{a,b\}:a<b,
 (P_i^TAP_i)_{ab}\ne A_{ab}\}.                   \tag{SV.6}
```

Thus `k_i=0` says that `P_i` is an exact switching automorphism of the
signed graph.  Notice that (SV.6) compares signed conjugates, not merely
the underlying vertex permutation.

## 2. The visibility theorem

On `U` define the compressed positive defect

```math
C=I_U-P_UT|_U\succeq0.                            \tag{SV.7}
```

For orbit `j`, let

```math
\bar d_j={1\over|\Gamma|}\sum_{g\in\Gamma}
 d(P_gz_j).                                       \tag{SV.8}
```

Projective invariance makes every term an individual deficit of a pole in
`P`.  Define the sharper visibility cost

```math
\delta_{\rm vis}
=\max_{\chi\in X}\min_{j:p_{j\chi}>0}
 {\bar d_j\over p_{j\chi}}.                      \tag{SV.9}
```

### Theorem SV.1 (sign-switching visibility synchronization)

Under (SV.1)--(SV.9),

```math
\boxed{
0\preceq C\preceq
\left(\delta_{\rm vis}
      +\sum_{i=1}^s{\sqrt{2k_i}\over r}\right)I_U.}          \tag{SV.10}
```

In particular, if every pole has `d(z)<=d`, then

```math
\boxed{
C\preceq
\left({d\over\nu}
      +\sum_{i=1}^s{\sqrt{2k_i}\over r}\right)I_U.}          \tag{SV.11}
```

For exact switching automorphisms the edge term vanishes.  If the declared
poles themselves form a character basis, `nu=1`, and separate Rayleigh
deficits give the relative certificate with no dimension loss.

More generally, let `V` be any (possibly redundant) pole presentation with
columns in `U`, and put

```math
G={V^TV\over n},\qquad
R={V^TAV\over rn},\qquad D=G-R.                  \tag{SV.12}
```

Then (SV.10) gives the Gram-relative statement

```math
\boxed{
D\preceq
\left(\delta_{\rm vis}
      +\sum_i{\sqrt{2k_i}\over r}\right)G.}       \tag{SV.13}
```

This remains meaningful when `G` is singular: it is a quadratic-form
inequality, and `ker G subseteq ker D` automatically.

#### Proof

Twirl the compressed defect:

```math
\overline C={1\over|\Gamma|}\sum_{g\in\Gamma}
 P_g^*CP_g.                                       \tag{SV.14}
```

Multiplicity freeness and character orthogonality make `bar C` diagonal
on (SV.2).  Write its nonnegative diagonal entries as `c_chi`.  Orbit
averaging gives the exact identity

```math
\bar d_j=\langle u_j,\overline C u_j\rangle
          =\sum_{\chi\in X}p_{j\chi}c_\chi.       \tag{SV.15}
```

Every summand is nonnegative.  For each `chi`, choose any orbit seeing it;
then `c_chi<=bar d_j/p_(jchi)`.  Minimizing over `j` and maximizing over
`chi` proves

```math
0\preceq\overline C\preceq\delta_{\rm vis}I_U.   \tag{SV.16}
```

If every deficit is at most `d`, (SV.4), (SV.8), and (SV.15) instead give
`c_chi<=d/nu`.

It remains to bound the nontwirled coherence.  Averaging the telescoping
commutator estimate over the `2^s` group words gives

```math
\|C-\overline C\|_{op}
\le {1\over2}\sum_{i=1}^s\|[T,P_i]\|_{op}.       \tag{SV.17}
```

Compression to the invariant space `U` cannot increase these norms.  Also

```math
P_i^TAP_i-A=P_i^T[A,P_i],                         \tag{SV.18}
```

so the two matrices have equal operator norm.  Every mismatched unordered
edge contributes two entries of magnitude two to (SV.18).  Consequently

```math
\|[T,P_i]\|_{op}
\le {\|P_i^TAP_i-A\|_F\over r}
={\sqrt{8k_i}\over r}.                           \tag{SV.19}
```

Equations (SV.17)--(SV.19) yield
`||C-bar C||<=sum_i sqrt(2k_i)/r`.  Combining this with (SV.16) proves
(SV.10)--(SV.11).  Finally,

```math
c^TDc={1\over n}\langle Vc,CVc\rangle,
\qquad c^TGc={1\over n}\|Vc\|_2^2,              \tag{SV.20}
```

which proves (SV.13). `square`

### Why this is more than commutator control

The commutator term in (SV.10) controls only off-character coherence.  The
new content is (SV.15): nonnegative character defects can be recovered
from **individual Boolean pole deficits** through a finite orbit incidence
matrix.  The visibility `nu` says exactly when that recovery has no growing
condition-number loss.  Without quantitative visibility, multiplicity
freeness can coexist with a pole orbit placing vanishing weight on a bad
character, and separate pole estimates do not give a dimension-free
relative certificate.

The factor `1/nu` is the optimal consequence of the nonnegative linear
system (SV.15): if one character has defect `c` and an orbit sees it with
weight `nu`, that orbit pays only `nu c`.  Extra structure can improve the
bound only by imposing additional relations among the character defects.

## 3. Sharp hollow PC.3 test

Index the order-16 coordinates by `x=(u,v) in F_2^2 times F_2^2`, put

```math
q(x)=u\cdot v,
\qquad H_{xy}=(-1)^{q(x)+x\cdot y+q(y)},          \tag{SV.21}
```

and obtain the hollow exact signing `A` by deleting the diagonal of `H`.
An exact characteristic-polynomial calculation gives

```math
\chi_A(\lambda)
=(\lambda^2-25)(\lambda^2-17)^4(\lambda^2-9)^3, \tag{SV.22}
```

so `||A||_op=5`.

Let `a,b,c` be the three PC.3 seed poles and `e=a odot b odot c`; explicit
vectors are recorded in the verifier.  Their span has dimension three,
with the sole displayed relation `a-b-c+e=0`.

For an even-weight translation `t`, let

```math
\ell_t(x)=q(x+t)+q(x)+q(t)+x\cdot t,
\qquad
(P_tf)(x)=(-1)^{\ell_t(x)}f(x+t).                \tag{SV.23}
```

The two involutions `P_3,P_5` commute and are exact switching
automorphisms of both `H` and `A`.  On the four projective poles,

```math
\begin{array}{c|rrrr}
 &a&b&c&e\\ \hline
P_3&-a&e&-c&b\\
P_5&-c&-e&-a&-b
\end{array}                                      \tag{SV.24}
```

so the two pole orbits are represented by `a` and `b`.  The action on
their three-dimensional span is multiplicity-free, with characters

```math
(+,-),\quad(-,+),\quad(-,-).
```

In that order, the two incidence rows are

```math
p_a=(0,1/2,1/2),
\qquad p_b=(1/2,1/2,0),                           \tag{SV.25}
```

and hence `nu=1/2`.

Every one of `a,b,c,e` has normalized Rayleigh value `4/5` against the
hollow signing, so `d=1/5`.  Theorem SV.1 proves

```math
\boxed{D\preceq {2\over5}G.}                     \tag{SV.26}
```

This is sharp: the three character defects are respectively

```math
{2\over5},\quad0,\quad{2\over5}.                 \tag{SV.27}
```

Thus the theorem recovers a nonzero, exact Gram-relative synchronization
constant for the **hollow exact signing**, not merely for the diagonal
Hadamard completion.  The example also shows why the visibility factor is
real: the four redundant Boolean poles do not themselves form the
orthogonal character basis.

## 4. Non-Hadamard Cayley test

Let `G=F_2^3`, set

```math
f=(0,-1,-1,-1,-1,-1,+1,+1),
\qquad A_{xy}=f(x+y),                             \tag{SV.28}
```

in lexicographic order.  This is a hollow symmetric exact signing.  Its
Walsh eigenvalues are

```math
(-3,1,-3,1,-3,1,5,1),                            \tag{SV.29}
```

so it is not Hadamard (`A^2` is not scalar) and `||A||_op=5`.  Translation
permutations form an exact `F_2^3` automorphism group.  The Boolean Walsh
characters `w_s(x)=(-1)^(s dot x)` are distinct one-dimensional character
spaces for that action.

Take poles `w_s` for `s in {1,3,5,6,7}`.  Each is projectively fixed by
translation, so `nu=1`.  Their relative defect spectrum is

```math
(4/5,4/5,4/5,0,4/5).                             \tag{SV.30}
```

Theorem SV.1 therefore gives `D<=(4/5)G`, with equality.  This confirms
that the mechanism is abelian switching diagonalization plus visibility,
not a special consequence of a Hadamard square identity.

## 5. Scope and stopping boundary

1. Multiplicity-free refers to the represented subspace `U`, not to the
   ambient coordinate representation.  Ambient character multiplicities
   are harmless.  Multiplicity inside `U` is not: twirling would leave a
   matrix block, and scalar orbit deficits need not control its top
   eigenvalue.
2. The pole family may be redundant and `G` singular.  What matters is
   projective orbit closure and spanning.  A pole list that is not closed
   under the declared action cannot use the orbit average (SV.15).
3. Approximate switching is useful at the dense `r~sqrt(n)` scale only
   when `sum_i sqrt(k_i)=o(sqrt(n))`, unless a stronger signed-discrepancy
   estimate replaces the Frobenius bound.  Merely changing `o(n^2)` edges
   is far too weak.
4. The theorem certifies one declared generated pole span.  It neither
   constructs such a span in arbitrary near-minimizers nor reconstructs
   the Boolean energy landscape.

## 6. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_sign_switching_visibility.py
```

The script checks the switching actions, exact orbit weights, PC.3
character polynomial and sharp relative spectrum, the mismatch constant,
and the non-Hadamard Cayley benchmark using exact integer/rational
arithmetic.
