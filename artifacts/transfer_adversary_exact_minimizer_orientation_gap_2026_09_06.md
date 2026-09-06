# Exact-minimizer orientation balance: precise open target and a local trap

Date: 2026-09-06. The question whether EVERY exact minimizer satisfies
`|P-R|<=2`, or at least `|P-R|=O(N)`, is NOT resolved here.
The proved scalable construction below concerns strict local minima,
not global minima. Root and adversarial agents reconstructed the same
Hadamard trap independently during this check.

Conventions: `H_A(x)=sum_(i<j) A_ij x_i x_j`, `P=max H_A`,
`R=max(-H_A)`, `Q=max(P,R)`, and `M_N=min_A Q(A)`.

## 1. Exact coding and parity statements

Encode the negative edges of A by `a in F_2^m`, `m=binom(N,2)`.
Let C be the binary cut code of K_N and let j be the all-edge vector.
For a spin x, its disagreement edges form a cut c, and
`H_A(x)=m-2 wt(a+c)`. Therefore

```math
P=m-2d(a,C),\qquad R=m-2d(a,j+C).
```

For N>=3, `D=C union (j+C)` is the augmented binary linear cut code.
Exact Q-minimizers are precisely the deepest cosets of D, since
`Q(A)=m-2d(a,D)`. The proposed gap-two assertion is equivalent to
the distances of the two C-cosets inside EVERY deepest D-coset
differing by at most one. This is an exact reformulation, not an
imported covering-code theorem or an argument for the assertion.

At odd N, every cut has even size, so its signed cut sum is even.
Thus all energies of a fixed signing are congruent modulo four, and
`P+R` is divisible by four. In particular, when M_N is odd, possible
dominant gaps are 2,6,10,..., not 0,4,8,... . This explains part of
the finite observed pattern. At even N, flipping one spin changes
energy by two modulo four, so this stronger congruence is unavailable.

## 2. A finite neighborhood condition that really follows from global minimality

Orient an exact minimizer so `P=M_N>=R` and put `g=P-R`.
If F is a set of f edges and `A^F` reverses them, then pointwise

```math
R(A^F)\le R(A)+2f.
```

Whenever `2f<g`, this negative cap is strictly below M_N. Global
minimality consequently forces

```math
P(A^F)\ge M_N.                                         (1)
```

Thus a gap g forces one-sided P stability throughout the entire edge
Hamming ball of radius strictly less than g/2. This is stronger than
the one-edge witness condition, but is still only a necessary condition.
No uniform O(N) escape radius, and hence no O(N) global gap bound,
has been proved from it.

The continuous relaxation does not supply a missing stationarity law:
allowing coefficients in `[-1,1]` admits the zero matrix with cap zero.
The discrete vertex minimization is a different problem. Symmetry under
global edge reversal gives equally good matrices with reversed gaps,
not a single exactly minimizing signing with balanced oriented extrema.

## 3. Scalable strict one-edge local minima with gap N

Let `H_4=J_4-2I_4`, `W=H_4^(tensor k)`, `N=4^k`, `d=(-1)^k`,
and take k>=2. Define

```math
A=I_N-dW.                                              (2)
```

Since W has diagonal d, A is a symmetric hollow signing. The characters
of `(F_2^2)^k` form a Boolean orthogonal eigenbasis for W: on one
factor the trivial character has eigenvalue 2 and each nontrivial
character eigenvalue -2. Hence W has eigenvalues `+sqrt(N)` and
`-sqrt(N)`, each attained by Boolean vectors. It follows exactly that

```math
P(A)=N(\sqrt N+1)/2,\qquad
R(A)=N(\sqrt N-1)/2,\qquad P-R=N.                       (3)
```

Uniformly sample the Boolean eigencharacters in the eigenspace of W
with eigenvalue `-d sqrt(N)`. They are positive grounds of A. The
orthogonal projection and its rank show that their raw second moment is

```math
K=E xx^T=\frac{I-dW/\sqrt N}{1-1/\sqrt N}.
```

For every off-diagonal edge e={i,j}, this gives

```math
E[A_{ij}x_ix_j]=1/(\sqrt N-1)<1.                       (4)
```

The variable in (4) is a sign, so some positive ground has value -1
on that edge. Flipping the edge increases its energy by exactly two.
The reverse triangle bound changes Q by at most two. Consequently

```math
Q(A^e)=P(A^e)=Q(A)+2
```

for EVERY edge e. These matrices are strict one-edge local minima of
both Q and the dominant oriented cap P, despite gap N.

More generally, for f edges F the sum of their signed responses under
this same law has mean `f/(sqrt(N)-1)` and has parity f. If
`f<sqrt(N)-1`, some ground has that sum at most zero; for odd f it
is at most -1. Thus

```math
Q(A^F)\ge Q(A),
\quad\text{and for odd f,}\quad Q(A^F)\ge Q(A)+2.        (5)
```

These are local traps, not an exact-minimizer construction. Their
normalized caps tend to one half, strictly above the banked all-order
upper bound for the actual minimum. This family neither disproves
gap-two for global minimizers nor disproves their possible O(N) gap.

## 4. An exact order-sixteen escape and exhaustive finite checks

For k=2, A in (2) has `P=40`, `R=24`. Independently enumerating all
32768 projective spins, all 120 single-edge changes, and all 7140
two-edge changes verifies the preceding local statements. In the
natural tensor ordering, reverse precisely

```text
(0,5), (0,10), (0,15).
```

The resulting signing has `P=38`, `R=30`, `Q=38`. Hence the minimum
number of edge changes that can reduce either P or Q is exactly three,
equal here to `sqrt(N)-1`. In particular this matrix fails the stronger
necessary GLOBAL condition (1), which its gap 16 would require for
every f<8. Independently, the verified order-sixteen upper witness has
cap 30, so the local trap is not globally minimizing.

The same checker exhausts every root-normalized signing through N=8,
not only archived minimizer representatives. It finds no signing of
those sizes with `P-R>=4` that is a one-edge local minimum of P.
At N=8 it checks 2,097,152 signings, including 414,670 with such a gap.
The larger explicit trap demonstrates why that finite observation
cannot be extrapolated to a general local-balancing lemma.
The exact global minimizing gap sets in the complete N=3,...,8 screen
are respectively `{2},{0},{0},{0},{2},{0}`.

Source: `computations/transfer_adversary_orientation_local_vs_global_2026_09_06.py`.
Full matrices, edge lists and counts:
`computations/results/transfer_adversary_orientation_local_vs_global_2026_09_06.json`.
Replay with `--max-order 8`; all arithmetic and pass/fail comparisons
are integer. Root-gauge reduction is exhaustive because switching
preserves both oriented extrema and commutes with flipping a specified
edge. No solver status is used as an infeasibility certificate.

## 5. Unresolved boundary

No general proof or exact-minimizer counterexample for gap at most two
was obtained in this bounded check. A possible O(N) bound for actual
global minimizers is also still open. The planted-clique near-minimizer
construction and the present strict local minima do not answer either
global question. They explicitly identify two hypotheses that cannot
be substituted for exact global optimality.
