# Wave 15 memo: random balanced partitions

Status: the mean and variance identities below are exact and independently
checkable.  The random-principal-submatrix corollary uses the stated theorem
of Tropp/Rudelson--Vershynin; it is useful, but it also proves that this route
throws away the leading internal excess.  No uniform centered-mask theorem
that retains leading internal excess was proved.

## 1. Exact centered-mask reduction

Assume first that `s|n`, and let `P` be a uniform equipartition into blocks of
size `s`.  Write

```
D_P = direct_sum_{V in P} A[V],   C_P=A-D_P,
p=(s-1)/(n-1).
```

Every pair is in the same block with probability `p`, so, for every fixed
Boolean state `z`,

```
E_P z^T D_P z = p z^T A z.
```

Consequently, with `E_P=D_P-pA`,

```
C_P=(1-p)A-E_P,
Q(C_P) <= (1-p)Q(A)+Q(E_P).                 (1)
```

Thus `Q(E_P)=o(n^(3/2))` would prove the desired cross-only estimate even if
the uncentered internal norm is leading.  This is the non-tautological target.

By contrast, a partition satisfying

```
sum_V Q(A[V])=o(n^(3/2))
```

also satisfies

```
0 <= sum_V (Q(A[V])-q_s) <= sum_V Q(A[V])=o(n^(3/2)).
```

It therefore captures no leading internal excess in the common-mosaic
identity.  It proves `Q(C_P)<=q_n+o(n^(3/2))` only by making the quantity one
wanted to harvest subleading.

## 2. Exact fixed-state variance

Fix `z`, put `b_ij=A_ij z_i z_j` for `i<j`, and set

```
N=n(n-1)/2,  S=sum_{i<j} b_ij=z^T A z/2,
L=sum_i (sum_{j != i} b_ij)^2=||Az||_2^2.
```

For `I_e=1` when the endpoints of edge `e` share a block, the two distinct
indicator covariances are

```
ca = Cov(I_ij,I_ik)
   = -p(n-s)/((n-1)(n-2)),
cd = Cov(I_ij,I_kl)
   = 2p(n-s)/((n-1)(n-2)(n-3)).             (2)
```

Here the second line is for four distinct vertices.  Since

```
sum_adjacent_edge_pairs b_e b_f=(L-2N)/2,
sum_all_edge_pairs b_e b_f=(S^2-N)/2,
```

direct covariance expansion gives the exact formula

```
Var_P(z^T D_P z)
 =4{p(1-p)N+(ca-cd)(L-2N)+cd(S^2-N)}.       (3)
```

The apparently dangerous row-field term has a *negative* coefficient.  Put
`alpha=cd-ca=p(n-s)/((n-2)(n-3))`.  If `Q(A)=O(n^(3/2))`, then (3) implies,
uniformly in `z` and `1<s<=n/2`,

```
Var_P(z^T D_P z)
 <=4{pN+2 alpha N+cd Q(A)^2/4}=O(ns).        (4)
```

No operator norm is lost at the second-moment level.

## 3. Why (4) does not give the Boolean supremum

There are `2^(n-1)` projective Boolean states.  A genuinely subgaussian tail
with variance proxy `ns` would union-bound at deviation `O(n sqrt(s))`, which
is exactly what is wanted.  Available bounded-jump arguments do not give that
tail: exchanging two vertices changes `z^T D_P z` by `O(s)`.  Even granting
the optimistic Bernstein form

```
P(|T-ET|>=t) <= 2 exp[-c t^2/(ns+s t)],       (5)
```

at `t=n sqrt(s)` its exponent is only `Theta(n/sqrt(s))`, not `Theta(n)`.
Thus (5) cannot be union-bounded over Boolean states for growing `s`.

The same obstruction is the large-deviation/operator term in a
Potts/simplex Hanson--Wright formulation.  The Frobenius term is
`||A||_F^2/d = Theta(ns)` for `d=n/s` colors and is perfect; the second term
must also survive a `2^n` supremum.  A hypothetical dimension-free vector
Hanson--Wright bound would lead to the natural uniform scale

```
n sqrt(s) + s ||A||_op.
```

This is subleading at `s~sqrt(n)` because the verified global bound is
`||A||_op<=sqrt(2Q(A))=O(n^(3/4))`.  However the rescaled categorical simplex
row is not dimension-free subgaussian: its squared subgaussian parameter is
of order `d/log d` in the direction of a simplex vertex.  Therefore ordinary
vector Hanson--Wright cannot be quoted to obtain this bound; a
categorical/permutation-specific chaos theorem is still missing.

For perspective, the all-positive matrix (not a global minimizer) shows why
uniformity is real.  If the number of blocks is even, choose `z` constant on
each block with half the block signs positive and half negative.  Then

```
|z^T(D_P-pA)z| = n^2(s-1)/(n-1) = Theta(ns)
```

for every partition, despite the `O(ns)` variance for each state fixed before
the partition is drawn.

## 4. A rigorous but self-defeating random-restriction corollary

Tropp's random-principal-submatrix theorem (built on
Rudelson--Vershynin) says that for a Bernoulli coordinate projector `R` of
density `delta`, a zero-diagonal Hermitian matrix with entries bounded by one
satisfies

```
(E ||RAR||_op^r)^(1/r)
 <= C[log n + sqrt(delta*n*log n)+delta ||A||_op]       (6)
```

when one takes `r=2 log n`.  This is a direct specialization of Theorem 1.1
in J. Tropp, *Norms of random submatrices and sparse approximation*; the
underlying primary source is Rudelson--Vershynin, *Sampling from large
matrices: an approach through geometric functional analysis*,
arXiv:math/0503442.

A uniform fixed-size `s` subset is coupled to a Bernoulli subset of mean
`2s`: condition on having at least `s` vertices and choose a uniform
`s`-subset of it.  Principal-submatrix monotonicity and the constant lower
bound on this event transfer (6), up to an absolute constant, with
`delta=2s/n`.  Linearity of expectation over the (dependent) blocks of a
uniform equipartition then gives

```
E_P sum_V Q(A[V])
 <= C[n log n+n sqrt(s log n)+s ||A||_op].              (7)
```

Indeed `Q(A[V])<=s||A[V]||_op`, and each block marginal is a uniform
`s`-subset.  For a global minimizer, `Q(A)=q_n=O(n^(3/2))` and hence
`||A||_op=O(n^(3/4))`.  At `s~sqrt(n)`, (7) is
`O(n^(5/4)sqrt(log n))=o(n^(3/2))`; an equipartition with this bound exists.
The same argument works for almost-equal block sizes and a leftover block.

This rigorously proves the original small-total-internal-norm proposal at a
useful scale (up to `sqrt(log n)`), but Section 1 shows why it cannot capture
a leading internal excess.  It should be recorded as a no-go/tautology, not
as closure of the adaptive cross-mosaic route.  The surviving target is (1)
with centered cancellation while `sum_V Q(A[V])` is allowed to be leading.

## Checker

`tmp/adaptive_partition_r15_check.py` exhausts small equipartitions and
Boolean states.  It verifies (2), (3), the exact mean, inequality (1), and
the all-positive `Theta(ns)` adaptive-state witness.
