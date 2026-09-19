# Stability-filtered selection: a repairable row-count obstruction

September19, uniform track. This is a bounded audit of the new selector
in `twisted_chiral_stable_permanent_2026_09_19.md`. It does not establish
the required average over child spin pairs or a native doubling theorem.

The solver-free companion
`computations/twisted_chiral_stability_row_bound_audit_2026_09_19.py`
checks 24 actual full-sign conditional profiles by direct integer fields,
280 rational scaled-permanent instances, and the Sylvester calibration
at orders 4,16,64. Every check passes. Its
results JSON has the same basename under `computations/results/`.

## 1. An actual full-sign high-energy profile with an impossible matching

Let m>=24, n=2m, and

```math
 A_m=\begin{pmatrix}J_m-I_m&J_m\\J_m&-J_m+I_m\end{pmatrix}.
```

Its quadratic cap is exactly Q(A_m)=m^2. Indeed its energy in the two
block sums u,v is `(u^2-v^2)/2+uv`; maximize first in u at an endpoint
and then in v, and similarly for the negative side.

Take the projective parent spins `x=(1,1), y=(1,-1)`, with blocks of
length m. Thus I_+ is the first block. Fix J_+ to be that same block.
Choose z negative on any h coordinates of the first block and positive
everywhere else, where `1<=h<=m-2`; put w=z on the first block and
w=-z on the second. This is an allowed conditional child pair in the
director's distribution, not an independently chosen artificial margin
matrix.

The parent-side row arrays are

```math
 alpha=(2m-1,1),       -gamma=(1,2m-1),
```

constant within their indicated blocks. For the columns, the exact
arrays and margins are as follows:

| column type | b_j | c_j | relevant-sector m_ij |
| --- | ---: | ---: | ---: |
| first block, z_j=+1 | -2h-1 | 2m-2h-1 | 2m-2h-2 |
| first block, z_j=-1 | 2h-1 | -2m+2h-1 | -2m+2h |
| second block | 2m-2h-1 | 2h-1 | 2m-2h |

Every nonforbidden margin is at least2. The plus-sector support has
all m rows adjacent to the **same** m-h columns. It has no perfect
matching. The minus-sector support is complete. Therefore the exact
conditional stability probability is zero for every matching energy,
and the exact generating permanent vanishes identically.

This example lies in the relevant high-energy tail, not in a tail
automatically discarded by the energy filter. Direct summation gives

```math
 E_0=4m^2-4mh+4h^2-2m >= 3m^2-2m.
```

Since the matching contributes at least-2m and m>=24,

```math
 E_0-2m >= 3m^2-4m > 2sqrt(2) Q(A_m).
```

Thus every matching energy would exceed the proposed threshold; only
the simultaneous local-stability compatibility rules these states out.

## 2. Quantitative loss of the row-only relaxation

Put `g(r)=(r!)^(1/r)`. All allowed edges in this example are strong, so
the exact coefficient-tail Brégman bound, its independent-edge variant,
and the optimized nonnegative-Chernoff bound all return

```math
 R_(m,h)={g(m-h)^m\over m!}>0,
```

instead of zero. The whole matching-energy polynomial belongs to the
tail, and the Chernoff expression is
`R_(m,h) exp(lambda(E0-L)) cosh(lambda)^(2m)`, minimized at lambda=0.
Consequently optimizing the Chernoff parameter does not repair this
example.

For h=1 the false-positive conditional mass satisfies exactly

```math
 R_(m,1)=g(m-1)/m -> e^(-1).
```

In particular, the row-count relaxation does not uniformly obtain even
an exponentially small stability factor in the high-energy regime.
There is also an extensive family of such impossible profiles: if
4 divides m and h=m/4, summing over the binom(m,m/4) choices gives

```math
 binom(m,m/4) R_(m,m/4)
 =Theta(m^(-1/3) 3^(m/4)),
```

by Stirling's formula, whereas the corresponding exact permanent sum
is zero. This is a leading exponential **unnormalized conditional**
overcount. The outer uniform z,J weights and the sum over parent spins
have not been evaluated here; this calculation does not imply that
the final averaged selector U_A(L) is at least one.

The child cap in this family is quadratic, not O(n^(3/2)). No
near-minimizer obstruction is claimed. Moreover this obstruction is
cheaply repairable: testing the support matrix for a perfect matching
zeros the conditional contribution immediately. It is not a barrier
to the exact permanent selector or to refinements retaining column
compatibility.

## 3. Audit of the director's scaled permanent refinement: PASS

Let K be an s-by-s nonnegative matrix, and suppose positive diagonal
matrices scale it to a doubly stochastic matrix
`B=diag(a) K diag(b)`. If, for an integer `1<=k<=s`, every entry of B
is at most1/k, then

```math
 per(K) <= [product_i a_i b_i]^(-1)
           (k!)^(s/k)/k^s.                         (1)
```

Proof: each row kB_i has entries in[0,1] summing to k, so belongs to the
convex hull of indicators of k-element subsets. Independently choose
one such indicator for each row, with mean kB_i, to produce a0-1
matrix U. Multilinearity in independent rows gives
`E per(U)=k^s per(B)`. Every row of U has degree exactly k; Brégman's
inequality therefore gives `per(U)<=(k!)^(s/k)` pointwise. Divide by
k^s and then undo the diagonal scaling. No independence of the original
bridge entries is asserted or needed.

For the usual capacity
`cap(K)=inf_(v>0) product_i(Kv)_i/product_j v_j`, the diagonal factor
in(1) is exactly cap(K). Indeed weighted AM--GM and the unit column
sums imply cap(B)>=1, while v=1 attains1. Capacity scales by
`product_i a_i b_i` under the displayed transformation.

Stirling gives

```math
 log per(K) <= log cap(K)-s+O((s/k) log(k+1)).
```

Thus k comparable to s/C gives logarithmic loss O(C log s), as stated
by the director. For k=s, B must be the uniform matrix and the bound
is exact; for k=1 it is the coarse bound per(K)<=cap(K).

Exact positive scaling is an explicit hypothesis. A matrix with a
perfect matching but edges belonging to no perfect matching need not
admit such a scaling before those irrelevant edges are removed. The
conditional theorem need not presume a scaling exists in every case.
The Hall-obstructed example above fails the support test outright.

Applied to K(exp(lambda)), this is a pointwise permanent/Chernoff
bound, not a coefficientwise Laurent-polynomial inequality: the
scaling may depend on lambda. Its improved column compatibility still
does not control the remaining average over child spin pairs, either
for low-cap seeds or near-minimizers. No uniform estimate of that
average has been proved in this audit.

## 4. Low-cap calibration: exact conditional tail probability one

The absence of a uniform per-quadruple stability penalty is not confined
to the preceding quadratic-cap family. Let H4 be the standard Sylvester
matrix, and take

```math
 x4=(1,1,1,-1),       y4=(-1,1,1,1),
 z4=(1,1,-1,1),       w4=(1,-1,1,1).
```

Direct multiplication gives `H4 x4=2x4`, `H4 y4=-2y4`,
`H4 z4=2w4`, and `H4 w4=2z4`. For n=4^k put H=H4 tensor-power k,
u=x4 tensor-power(k-1), and

```math
 x=x4 tensor u, y=y4 tensor u,
 z=z4 tensor u, w=w4 tensor u.
```

Thus Hx=sqrt(n)x, Hy=-sqrt(n)y, Hz=sqrt(n)w, and Hw=sqrt(n)z.
Both coordinate products x*y and z*w are balanced, so this quadruple
is admissible in the selector's conditional distribution. Put
`A=H-diag(H)` and `delta=diag(H)`. Since tr(H)=0, Boolean quadratic
energies of A and H coincide. The spectral bound and x's attained
eigenvalue prove exactly `Q(A)=n^(3/2)/2`.

Writing t'=z*w, the actual row arrays are

```math
 alpha_i=sqrt(n)-delta_i,
 -gamma_i=sqrt(n)+delta_i,
 b_j=c_j=sqrt(n)-delta_j t'_j.
```

Every margin therefore obeys

```math
 m_ij=2sqrt(n)-1-delta_j t'_j >= 2sqrt(n)-2 >= 2.
```

Both sector supports are complete and every edge is strong. Conditional
on this x,y,z,w, **every** sector-preserving permutation and **every**
matching is strictly locally stable. Moreover

```math
 E0=2n^(3/2)-sum_j delta_j t'_j,
 E0+sum_i u_i >= 2n^(3/2)-2n > 2sqrt(2)Q(A)
```

for n>=16. Thus the exact conditional probability of stability AND
energy above the target is one. The permanent formula, the row-count
bound, and the scaled refinement are all exact here; the sector-scaled
matrices are uniform doubly stochastic matrices. There is no Hall
obstruction or permanent-estimation loss to fix.

This falsifies a uniform small conditional bound based solely on
`Q(A)<=n^(3/2)/2`. It does **not** falsify the final average over
quadruples, since the frequency of these special quadruples has not
been estimated. It also does not concern selectable minimizing seeds:
the normalized cap1/2 is above the archived all-order near-minimizer
upper constant. The director independently reconstructed this example;
the companion script checks all fields, admissible sector sizes, and
one actual worst-matching parent at n=4,16,64 by integer arithmetic.
