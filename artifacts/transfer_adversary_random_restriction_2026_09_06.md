# Random principal restriction loses cap on every actual near-minimizer

Date: 2026-09-06. Independently proved finite estimate and asymptotic
counterexample to a sampling transfer. This is NOT nonconvergence of `M_n`,
and does not rule out rare useful restrictions or any near-order transfer.

## Theorem

Let `A_n` be ANY deterministic sequence of hollow symmetric signings with
`Q(A_n)<=C n^(3/2)` for a fixed finite `C`. Put

```math
k_n=\left\lfloor\sqrt{\tfrac12\log_2 n}\right\rfloor.
```

If `I_n` is a uniform `k_n`-vertex subset, then

```math
\Pr\left\{
\frac{Q(A_n[I_n])}{k_n^{3/2}}
\ge \frac23\sqrt{\frac2\pi}-\varepsilon
\right\}\longrightarrow1
```

for every fixed `epsilon>0`. The explicit constant is approximately
`0.5319230405352436`, strictly larger than `1/2`.

In particular, take ACTUAL minimizing seeds `A_n`. Their normalized caps
are at most the proved `c_*+o(1)<0.499432220485404+o(1)`, but a random
principal restriction at the displayed slowly growing orders has a normalized
cap at least `0.5319230405352436-o(1)` with high probability. Thus this
operation has an actual-minimizer normalized loss at least
`(2/3)sqrt(2/pi)-c_*-o(1)`, approximately `0.0324908200498404`.

The same lower bound applies to the excess above the optimal child cap,
because the all-order upper theorem also applies at the growing orders `k_n`.

## 1. Boolean cap implies a quantitative cut bound

The maximum of the absolute multilinear quadratic polynomial on the entire
cube `[-1,1]^n` equals its maximum on the Boolean vertices. For Boolean
`x,y`, write `u=(x+y)/2`, `v=(x-y)/2`. Then

```math
|x^TAy|=|u^TAu-v^TAv|\le4Q(A).
```

By linearity in each argument the same bilinear bound holds for all vectors
in `[-1,1]^n`, in particular indicators of vertex subsets.

Make a step graphon `W_A` on `n` equal intervals with values
`W_A(i,j)=(1+A_ij)/2`. Its diagonal blocks have value `1/2`. If `delta`
denotes the unnormalized cut norm of `W_A-1/2`, then

```math
\delta:=\sup_{S,T\subseteq[0,1]}
\left|\int_{S\times T}(W_A-\tfrac12)\right|
\le\frac{2Q(A)}{n^2}.
```

Fractional intervals cause no problem: their normalized intersection sizes
are vectors in `[0,1]^n`, already covered by the bilinear bound.

## 2. Finite total-variation estimate for induced labeled signings

Write `r=binom(k,2)`. Sample `k` independent uniform points from `[0,1]` and
generate the `r` edges conditionally independently with their graphon
probabilities. For any particular labeled sign pattern `F`, telescope its
probability against the constant graphon `1/2`, replacing one edge factor
at a time. In a term belonging to edge `(i,j)`, fix all the other vertex
variables. All remaining factors involving `i` multiply to a function
`f(x_i)` in `[0,1]`, and all remaining factors involving `j` multiply to
`g(x_j)` in `[0,1]`. The rest is a factor in `[0,1]` independent of both.
Layer-cake integration therefore bounds this term by `delta`. This works
for nonedge factors as well, because their difference is `-(W_A-1/2)`.
Thus

```math
|\Pr_{W_A}(F)-2^{-r}|\le r\delta.
```

Sum over the `2^r` labeled patterns and divide by two. Coupling independent
vertex labels to a uniformly chosen injective tuple fails with probability
at most `r/n`. Conditional on distinct labels the graphon edges are
deterministically exactly the induced signing of `A`. Consequently

```math
d_{\rm TV}\bigl(\mathcal L(A[I]),\mathcal L(G_k^{\pm})\bigr)
\le 2^{r-1}r\delta+\frac rn
\le \frac{2^r rQ(A)}{n^2}+\frac rn.                 (1)
```

Here `G_k^pm` has independent uniform signs above the diagonal. Vertex labels
of the subset are ordered uniformly; the cap itself is permutation invariant.
No independence of the edges of `A` is assumed.

At `Q(A)<=C n^(3/2)` and the stated `k_n`, `r<=log_2 n/4` and (1) is at
most `C r n^(-1/4)+r/n`, which tends to zero. This uniform quantitative
bound is the reason the theorem applies to actual minimizers, not merely to
an explicitly pseudorandom algebraic family.

## 3. Elementary random-sign lower cap above one half

For an iid signing `G_k^pm`, assign spins sequentially. Set `x_1=1`; having
assigned `x_1,...,x_(i-1)`, choose `x_i` to make
`x_i sum_(j<i) A_ij x_j` nonnegative. The resulting energy is

```math
H=\sum_{i=2}^k Z_i,\qquad
Z_i=\left|\sum_{j<i}A_{ij}x_j\right|.
```

Conditionally on all earlier edges and spins, the fresh row signs multiplied
by the fixed earlier spins are independent uniform signs. Hence the `Z_i`
are INDEPENDENT and have the same laws as `|S_(i-1)|`, where `S_l` is a sum
of `l` independent Rademacher variables. The exact formulas

```math
\mathbb E|S_{2j}|=\frac{2j\binom{2j}{j}}{2^{2j}},\qquad
\mathbb E|S_{2j+1}|=\frac{(2j+1)\binom{2j}{j}}{2^{2j}}
```

and the central-binomial asymptotic give

```math
\mathbb E H=\left(\frac23\sqrt{\frac2\pi}+o(1)\right)k^{3/2},
\qquad \operatorname{Var}(H)\le\sum_{i=2}^k(i-1)\le\frac{k^2}{2}.
```

Chebyshev's inequality yields `H/k^(3/2) -> (2/3)sqrt(2/pi)` in probability.
Since `Q(G_k^pm)>=H`, this proves the iid lower cap. Combining it with (1)
proves the theorem. No Parisi formula or spin-glass universality is needed.

### Quantitative rarity of useful restrictions

The lower tail can be made exponential. For `Z=|S_l|`, changing one
Rademacher input changes `Z` by at most two. The martingale bounded-difference
moment-generating-function bound is
`E exp(lambda(Z-E Z)) <= exp(lambda^2 l/2)` for every real `lambda`.
Independence of the greedy increments and Chernoff optimization give

```math
\Pr\{H\le\mathbb EH-t\}
\le\exp\left(-\frac{t^2}{k(k-1)}\right).
```

Let `mu_k=sum_(l=1)^(k-1) E|S_l|`, using the exact binomial formulas above.
Combining this inequality with (1) yields, whenever `mu_k>c_0 k^(3/2)`,

```math
\Pr\{Q(A[I])\le c_0 k^{3/2}\}
\le \exp\left(-\frac{(\mu_k-c_0k^{3/2})^2}{k(k-1)}\right)
 +\frac{2^r rQ(A)}{n^2}+\frac rn.                 (2)
```

For every fixed `c_0<(2/3)sqrt(2/pi)` and the theorem's `k_n`, useful
restrictions therefore have probability at most
`exp(-[((2/3)sqrt(2/pi)-c_0)^2/4] k_n)+O_C(log n n^(-1/4))`
for all sufficiently large `n`. In particular, restrictions attaining the
known original upper interval are exponentially rare in their child order
within every bounded-cap seed. They are not proved nonexistent.

## What this falsifies, and what it does not

This theorem falsifies a lossless random restriction claim uniform over
actual optimizing seeds and widely separated scales. It also shows that
ordinary dense graph limit information, which sends every bounded-cap seed
to the constant graphon `1/2`, recovers high-cap random restrictions instead
of the extremal child cap.

It does NOT falsify the principal-restriction monotonicity `Q(A[I])<=Q(A)`,
whose normalization is different. It does NOT rule out the existence of
exceptional subsets with cap `M_k+o(k^(3/2))`, even though such subsets must
have probability tending to zero at these scales. It does not address a
fixed positive retention ratio, because `k_n/n ->0`, and it does not imply
that `M_n/n^(3/2)` fails to converge.

## Reproducible finite checks

`computations/transfer_adversary_random_restriction_2026_09_06.py` enumerates
all `32768` iid order-six signings and verifies the full joint increment law
equals the product law exactly, not merely its moments. It also computes an
exact induced-pattern total variation and exact cut norm for the archived
order-eight minimizing seed. These sanity checks pass; the latter small
finite TV bound is deliberately reported even though it is numerically
vacuous. The asymptotic conclusion rests on estimate (1), not that example.
