# Independent audit of the mesoscopic packing-inflation construction

Status: algebra and explicit compatible orders checked. This is a scalable
full-sign **bounded-cap** obstruction to no packing inflation, not an
obstruction for actual near-minimizers or selected optimal children.
The construction was proposed by the constructive agent.

## 1. Exact family and its symmetry requirement

Take b=2^r, s=64b, n=5bs. Let H_b,H_s be symmetric Sylvester Hadamards.
Within each of b mesoscopic blocks, place `A5 tensor J_s` on its five distinct
s-groups and fill each s-group by `F_s=H_s-diag(H_s)`. Fill the blocks between
distinct mesoscopic groups by

```math
(H_b-\operatorname{diag}H_b)\otimes H_{5s}.
```

Here H_(5s) must be **symmetric**; an unspecified Hadamard factor would not
guarantee that the assembled matrix is symmetric. Compatible explicit orders
are supplied by

```math
H_{5s}=H_{20}\otimes H_{16}\otimes H_b,
\qquad 5s=320b.
```

For a concrete symmetric H20, take the symmetric Paley conference C10 over
F9, C10²=9I, and set

```math
H_{20}=\begin{pmatrix}C_{10}+I&C_{10}-I\\
                     C_{10}-I&-C_{10}-I\end{pmatrix}.
```

It is symmetric, has all entries +-1, and squares to 20I. The independent
program `tmp/flatify_adversary_2026_09_07_symmetric_h20.py` builds F9 explicitly
as F3[a]/(a²+1), verifies the integer conference/Hadamard identities, and saves
both matrices in its JSON. Thus no conjectural Hadamard orders are needed.

The resulting A_n is an actual hollow symmetric full signing: its three
disjoint edge types are cross-mesoscopic, cross-s-group inside a mesoscopic
block, and inside an s-group.

## 2. Original cap bounds

The continuum-cube maximum of the hollow pentagon polynomial is its Boolean
cap 4. Therefore a mesoscopic `A5 tensor J_s` has cap at most 4s², because
the five group sums lie in [-s,s]. Its five diagonal fills cost at most
`5s(sqrt(s)+1)/2`. Across all mesoscopic blocks this is
`4bs²+O(bs^(3/2))`.

The cross-mesoscopic matrix has operator norm at most
`(sqrt(b)+1)sqrt(5s)`. Thus its quadratic cap is at most
`(1/2)(1+1/sqrt(b))n^(3/2)`. Combining the bounds gives

```math
{Q(A_n)\over n^{3/2}}
\le {32\over5\sqrt5}+{1\over2}+o(1).
```

For completeness this family really is far from near-minimal, not merely
poorly upper-bounded. Sylvester H_s has entry sum s and trace zero, so
`H_Fs(1)=s/2`. Choose a constant spin on each of the five s-groups according
to a positive pentagon maximizer, separately in every mesoscopic block.
Independently reversing whole mesoscopic blocks preserves their internal
energy and cancels cross-mesoscopic energy in expectation. Hence

```math
Q(A_n)\ge b(4s^2+5s/2),\qquad
\liminf {Q(A_n)\over n^{3/2}}\ge {32\over5\sqrt5}>2.
```

## 3. Both-polarity packing inflation

Inside each mesoscopic block use the five favorable pentagon triangles of
one polarity, replacing each selected seed vertex by its entire s-group.
Give each such ternary pattern weight 1/3. Each original coordinate belongs
to exactly three patterns, so its total packing load is one. Different
mesoscopic blocks have disjoint supports, so all these packings combine.

Each selected triangle contributes cross-group signed energy 3s². Its
within-group energy is sigma times the sum of the three F_s all-one energies.
Summing all five patterns with weight 1/3 yields exactly

```math
\kappa_\sigma(A_n)\ge5bs^2+\sigma\,5bs/2,
\qquad \sigma\in\{+1,-1\}.
```

Consequently both polarities satisfy

```math
\liminf {\kappa_\sigma(A_n)-Q(A_n)\over n^{3/2}}
\ge {8\over5\sqrt5}-{1\over2}> {1\over5}.
```

The final strict rational comparison follows by squaring positive numbers:
`8/(5sqrt5)>7/10` is equivalent to 6400>6125. This proves a fixed leading
gap for the specified full-sign bounded-cap class. The lower original-cap
bound in Section 2 explicitly prevents its misclassification as a sequence
of actual or asymptotic minimizers.

The generic obstruction does not establish packing inflation for any exact
large minimizer, and it does not decide the favorable flatification target.
