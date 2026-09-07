# Actual bounded-cap signings need not admit either small-rank ramp regime

Date: 2026-09-07. Status: proved deterministic counterexample to the two
rank/incoherence regimes in Sections 1--3 of
`decisive_audit_low_rank_midpoint_sign_surgery_2026_09_07.md` when imposed
on arbitrary bounded-cap parents. **Not a counterexample to midpoint
repair, or to a hypothesis restricted to actual width minimizers.**
Indeed this family has an explicit sparse repair, given in Section 4.

Write P=max H_A, R=max(-H_A), w=(P+R)/2, I=(P-R)/2. A ramp certificate
means an orthogonal projector Pi and e>=0 such that, for every Boolean x,

```math
H_A(x)\le w+{I\over n}x^T\Pi x+e.
```

## 1. Exact full-sign construction

Let F=H_2 tensor H_2, and let B be its diagonal switching by the vector
(-1,-1,1,-1). Explicitly,

```text
B = [[ 1, 1,-1, 1],
     [ 1,-1,-1,-1],
     [-1,-1,-1, 1],
     [ 1,-1, 1, 1]].
```

B is symmetric Hadamard and 1^T B 1=0. For s>=1 put

```math
b=4^s,\quad \ell=16b,\quad n=b\ell,
\quad H_b=F^{\otimes s},\quad H_\ell=B\otimes F^{\otimes(s+1)}.
```

Then both matrices are symmetric Hadamards, and 1^T H_l 1=0. Let D_b
be the diagonal matrix with the diagonal entries of H_b. Define

```math
A=(H_b-D_b)\otimes H_\ell+I_b\otimes(J_\ell-I_\ell).
```

A is hollow and symmetric. Within each of its b blocks all edges are +1;
between distinct blocks all entries are signs from the first tensor term.
There are no zero off-diagonal entries. Note l=4 sqrt(n).

Let C=(H_b-D_b) tensor H_l. Its operator norm is at most
sqrt(n)+sqrt(l). Every block-constant Boolean vector x=z tensor 1_l has
exactly zero C-energy, because 1^T H_l 1=0. Thus all 2^b such vectors have
the same exact energy

```math
E_0={n(\ell-1)\over2}.
```

The clique term has maximum E_0 and minimum at least -n/2. With
B_n=n(sqrt(n)+sqrt(l))/2, this gives

```math
E_0\le P\le E_0+B_n,\qquad R\le B_n+n/2,
```

and consequently

```math
Q(A)\le(5/2+o(1))n^{3/2},\quad
w\le(3/2+o(1))n^{3/2},\quad
(3/4-o(1))n^{3/2}\le I\le(5/4+o(1))n^{3/2}.
```

In particular I>0 eventually; no complementation is needed.

## 2. The ramp itself forces rank of order sqrt(n)

Suppose a ramp holds with e=o(n^{3/2}). Apply it to every block-constant
spin. The preceding endpoint bounds imply

```math
{x^T\Pi x\over n}\ge c_n,
\qquad
c_n={\ell-2-2\sqrt n-2\sqrt\ell-4e/n
       \over\ell-1+\sqrt n+\sqrt\ell}
\longrightarrow {2\over5}.
```

The numerator is positive eventually. This is a direct consequence of
the proportional RAMP, not the stronger near-complete high-set capture
condition previously ruled out.

Average over independent fair block signs z. Their covariance is
I_b tensor J_l=l P_block, where P_block is the projector onto block-
constant vectors. If r=rank(Pi), then

```math
c_n n\le\mathbb E_z[x^T\Pi x]
   =\ell\operatorname{Tr}(\Pi P_{\rm block})\le\ell r.
```

Therefore every such ramp satisfies

```math
\operatorname{rank}(\Pi)\ge(1/10-o(1))\sqrt n.
```

In particular the general surgery regime r=o(sqrt(n)) cannot be asserted
from bounded cap and full-sign structure alone.

## 3. The alternative incoherence regime also fails

Let kappa=max_{i!=j}|Pi_ij|. The same block-sign average gives

```math
c_n n\le\operatorname{Tr}\Pi+n(\ell-1)\kappa.
```

If r=o(n), it follows that

```math
\kappa\ge {c_n-r/n\over\ell-1}
          =(1/10-o(1))n^{-1/2}.
```

Thus there is no ramp with e=o(n^{3/2}), r=o(n), and simultaneously
kappa=o(n^{-1/2}). These arguments also apply to positive contractions
with the corresponding trace budget; they do not depend on a particular
choice of projector basis.

## 4. Essential scope: sparse midpoint repair is easy here

Change only the within-block clique edges back to the corresponding signs
of H_b tensor H_l, and hollow its diagonal. This changes at most
b*l(l-1)/2=Theta(n^{3/2}) undirected edges. The repaired matrix is a hollow
symmetric Hadamard tensor, with cap at most
(1/2+o(1))n^{3/2}. Meanwhile w>=P/2>=E_0/2=(1-o(1))n^{3/2}.
Thus this explicit repair more than attains the old midpoint width.
The family is demonstrably far from globally width-minimizing.

The block-constant projector has rank b=sqrt(n)/4 and off-diagonal
entries 1/l within blocks and zero outside. It is **sign-compatible**:
A_ij(P_block)_ij>=0. Consequently the director's newer compatible-sign
surgery theorem, which permits r=o(n) at bounded operator norm, is not
contradicted. This observation does not itself assert that P_block obeys
the exact ramp with the actual endpoints. Its point is that the new
compatibility regime deliberately avoids both necessary-condition failures
above.

The result refutes an unqualified small-rank/incoherence premise for all
bounded-cap signings. It says nothing against selecting a ramp for actual
width minimizers, or against a different actual-sign midpoint surgery.

Reproducibility: `computations/decisive_independent_ramp_block_counterexample_2026_09_07.py`
passes exact checks at n=256, including all 16 block-constant spins,
the Hadamard identities, sign compatibility of the block projector, and
the explicit repair (4032 changed edges). The finite zero-error ramp
capture fraction from Section 2 is exactly 14/87 at that order.
