# Exact invariant theory and the failure of fixed `phi_6`

Status: exact algebraic results, exact arithmetic verification, and a proved
scalable falsifier.  The proposed restricted alternatives at the end are
open assessments, not established routes.

## 1. What the 16 counts encode

Let `A` be a symmetric zero-diagonal sign matrix.  For a simple edge set
`F`, write

```math
chi_F(A)=\prod_{ij\in F}a_{ij}.
```

Under switching by signs `s_i`, this monomial is multiplied by
`prod_i s_i^(deg_F(i))`.  It is therefore switching invariant exactly when
every degree in `F` is even.  Global negation multiplies it by
`(-1)^|F|`, so invariance under global negation additionally requires an
even number of edges.  Finally, permutation invariance is obtained by
summing `chi_F` over the complete isomorphism orbit of `F`.

Exhaustive enumeration gives respectively 2, 4, and 10 such graph orbits on
4, 5, and 6 vertices.  These match exactly the numbers of
switching/permutation/global-negation signing classes.  Evaluating every
orbit sum on one representative of every signing class gives square
character transforms with exact determinants

| order | classes / graph orbits | determinant | rank over `Q` |
|---:|---:|---:|---:|
| 4 | 2 | `-4` | 2 |
| 5 | 4 | `512` | 4 |
| 6 | 10 | `-35184372088832` | 10 |

Thus the class histogram and the vector of all even-Eulerian signed orbit
moments are equivalent coordinates at each of these orders.  This is an
exact finite Fourier transform, not a heuristic correlation.

For reference, the transforms at orders four and five (rows are graph
orbits and columns are signing classes in the program's canonical order) are

```text
t=4:  [ 1,  1]
      [ 3, -1]

t=5:  [ 1,  1,  1,  1]
      [15,  3, -1, -5]
      [15, -3, -1,  5]
      [ 1, -1,  1, -1].
```

The nonempty graph orbits here are `C4`, `C4+I`, the bowtie (two triangles
sharing a vertex), and `K5`.  At order six the ten orbits have the following
descriptors; these descriptors distinguish all ten without relying on a
name for each graph.

| index | edges | degree sequence | component sizes | triangles | labeled orbit size |
|---:|---:|---|---|---:|---:|
| 0 | 0 | `0,0,0,0,0,0` | `1,1,1,1,1,1` | 0 | 1 |
| 1 | 4 | `2,2,2,2,0,0` | `4,1,1` | 0 | 45 |
| 2 | 6 | `4,2,2,2,2,0` | `5,1` | 2 | 90 |
| 3 | 6 | `2,2,2,2,2,2` | `6` | 0 | 60 |
| 4 | 6 | `2,2,2,2,2,2` | `3,3` | 2 | 10 |
| 5 | 8 | `4,4,2,2,2,2` | `6` | 0 | 15 |
| 6 | 8 | `4,4,2,2,2,2` | `6` | 2 | 180 |
| 7 | 10 | `4,4,4,4,4,0` | `5,1` | 10 | 6 |
| 8 | 10 | `4,4,4,4,2,2` | `6` | 6 | 90 |
| 9 | 12 | `4,4,4,4,4,4` | `6` | 8 | 15 |

There is also exact redundancy in the advertised 16 coordinates.  If `h_6`
is the ten-entry order-six histogram, double-counting pairs consisting of an
`r`-set and a containing six-set gives

```math
h_r={I_{r,6}h_6\over {n-r\choose 6-r}},\qquad r=4,5,              \tag{1}
```

where the integer incidence matrices `I_(r,6)` are recorded in the JSON
artifact.  Hence for every `n>=6`, the order-six histogram alone determines
the order-four and order-five histograms.  The 16-count state has only the
information of its ten order-six counts (which themselves have a fixed sum).

## 2. Exact spectral and Boolean-energy moments

Let `Z_G(A)` denote the signed orbit sum for graph `G`, over every appropriate
vertex subset, and put `E=binom(n,2)`.  Direct closed-walk and even-edge
multigraph enumeration gives

```math
tr(A^4)=n(n-1)(2n-3)+8Z_{C4},                                      \tag{2}
```

and

```math
tr(A^6)=2{n\choose2}+60{n\choose3}+120{n\choose4}
 +120Z_{C4}+48Z_{C4+I}+24Z_{bowtie}+12Z_{C6}.                      \tag{3}
```

For `H_A(x)=sum_(i<j) a_ij x_i x_j` and uniform Boolean `X`, the corresponding
energy identities are

```math
E H_A(X)^2=E,
```

```math
E H_A(X)^4=3E^2-2E+24Z_{C4},                                      \tag{4}
```

```math
E H_A(X)^6
=E+15E(E-1)+90{E\choose3}+(360E-960)Z_{C4}
 +720(Z_{bowtie}+Z_{C6}+Z_{2C3}).                                 \tag{5}
```

Equations (2)--(5) were independently checked by exact integer evaluation on
seeded random signings of orders 6, 7, 8, and 9.  They show precisely why the
profile was a good finite predictor: it fixes the first three even energy
moments and the first three even spectral moments.  They do not supply a tail
theorem.

Indeed, the strongest immediate spectral estimate is

```math
cap(A) <= {n\over2}||A||_op <= {n\over2}(tr(A^6))^{1/6}.            \tag{6}
```

Even at the optimal trace scale `tr(A^6)=Theta(n^4)`, (6) is only
`O(n^(5/3))`, rather than the required `O(n^(3/2))`.  More generally a known
`2q`-th uniform energy moment yields only the atom-mass bound

```math
cap(A) <= 2^{(n-1)/(2q)} (E|H_A(X)|^{2q})^{1/(2q)}.                 \tag{7}
```

Thus any generic moment-to-maximum argument needs `q=Omega(n)` even for a
constant-factor loss.  Taking local moments through only
`k=Theta(sqrt(n))` vertices leaves an `exp(Theta(sqrt(n)))` factor in (7).
A successful local-state theorem would need rigidity or a new tail mechanism,
not the finite moment conversion itself.

## 3. The exact collision kills universal profile control

The independently verified order-ten signings with root-gauged codes

```text
5850642905
28771662001
```

have equal complete oriented restriction-class histograms at every order at
most six, but have energy ranges `[-17,19]` and `[-17,21]`.  In particular
they have equal `phi_6` and caps 19 and 21.

Replace every base vertex by `L` positive twins.  Writing `z_i` for the sum
of the `L` spins in the `i`-th twin class and `Q=A+I`, the blowup energy is

```math
H_{T_L(A)}(z)={1\over2}z^TQz-5L.                                  \tag{8}
```

The oriented profile equality through order six is preserved: every small
restriction is specified by an oriented base restriction and an occupancy
vector, and twin replacement is a fixed operation on each oriented class.
Coordinatewise convexity of `z^TQz` puts its positive maximum at a box
vertex.  The supplied positive-definiteness certificates for `Q+4I` control
the negative extreme.  Consequently, for every integer `L>=3`,

```math
cap(T_L(A))=24L^2-5L,
\qquad cap(T_L(B))=26L^2-5L.                                      \tag{9}
```

The two order-`N=10L` signings have identical `phi_6` but cap gap
`2L^2=N^2/50`.  Therefore the proposed universal stability statement

```math
phi_6(A)=phi_6(B) => |cap(A)-cap(B)|=o(n^(3/2))
```

is false by much more than the project scale.  Fixed `phi_6`, and hence its
equivalent fixed collection of signed moments, is no longer a live universal
landing state.

## 4. Is a low-cap restriction still a useful target?

The blowups in (9) have quadratic cap, so they do not logically refute a
statement restricted to spectrally tame or otherwise certified low-cap
signings.  One exact noncircular candidate would be

```math
phi_6(A)=phi_6(B),
||A||_op,||B||_op <= C sqrt(n)
 => |cap(A)-cap(B)|=o(n^(3/2)).                                    \tag{10}
```

The added condition is polynomial-time checkable and excludes the positive-
twin obstruction.  But (10) is not presently a landing reduction.  Two new
bridges would still be required:

1. every exact minimizer, or at least one signing in its `phi_6` class, must
   satisfy the spectral condition uniformly; and
2. equality of the local state inside that condition must control the
   exponentially rare Boolean maximum.

Neither bridge follows from the fixed moments.  Restricting directly by
`cap(A)=O(n^(3/2))` is worse: it assumes essentially the property the family
was meant to certify, so it is circular as a tractable family definition.
Accordingly, a low-cap-only version remains logically open but is not a live
route without an independently checkable spike certificate that exact
minimizers are proved to possess.

## 5. How much growing restriction information gives a theorem?

There is a simple exact benchmark.  Deleting `r` vertices changes the cap by
at most the number of deleted incident edges,

```math
D(n,r)={n\choose2}-{n-r\choose2}=rn-{r(r+1)\over2}.                 \tag{11}
```

Hence, if two order-`n` signings have switching-equivalent principal
restrictions on some `n-r` vertices, then

```math
|cap(A)-cap(B)| <= 2D(n,r).                                        \tag{12}
```

For example, equality of the complete `(n-r)`-deck guarantees such a matched
restriction.  Formula (12) gives `o(n^(3/2))` only for
`r=o(sqrt(n))`; at the critical deletion radius `r=Theta(sqrt(n))` it gives
only an order-`n^(3/2)` error.  This state is also nearly the whole signing
and has exponentially many possible classes, so (12) is a scale benchmark,
not a simpler structured family.

The two direct growing-state ideas therefore bracket an obstruction:

- local signed moments on `Theta(sqrt(n))` vertices remain far too low-order
  for a generic maximum bound by (7);
- an almost-complete restriction with only `o(sqrt(n))` vertices omitted
  does give the desired error by (12), but retains essentially the full
  parent information.

A genuinely new state must sit between these extremes.  It must encode a
nonlocal Boolean-spike or bridge-response invariant, with growing but
substantially sub-full complexity, and it needs a quantitative update rule
under composition.  No such invariant is proved here.  The appropriate
classification is therefore:

1. universal fixed-profile landing: **falsified scalably**;
2. fixed profile plus a noncircular low-cap certificate: **open, with two
   missing bridges and no present landing theorem**;
3. near-complete restriction state: **proved landing bound but not simpler
   than full parent minimization**;
4. local state of size about `sqrt(n)`: **fixed-moment mechanism insufficient;
   a new nonlocal invariant is required**.

## Reproduction

```bash
.venv/bin/python computations/phase2_phi6_moment_transform.py \
  --output computations/results/phase2_phi6_moment_transform.json \
  > computations/logs/phase2_phi6_moment_transform.log
```

The JSON contains all three exact transform matrices, graph descriptors,
deck-incidence matrices, and the exact arithmetic formula checks.
