# Wave 15 memo: natural amplifications of the `A_9` wall

Status: **verified exact identities/no-go results**, with the finite arithmetic
checked by `tmp/a9_amplification_r15.py`.  These results concern the order-nine
`3+6` wall in (10.550)--(10.551), in the `Q` normalization.

## 1. Seed data

Write the seed in its `3+6` form

```math
A_9=\begin{pmatrix}D_3&B\\B^{\mathsf T}&D_6\end{pmatrix},
\qquad
C_9=\begin{pmatrix}0&B\\B^{\mathsf T}&0\end{pmatrix}.
```

Exact enumeration reconfirms

```math
Q(A_9)=24,
\qquad
\min x^{\mathsf T}A_9x=-24,
\qquad
\max x^{\mathsf T}A_9x=24,
```

```math
\|B\|_{\infty\to1}=12,
\qquad Q(C_9)=24,
\qquad \mathbf 1^{\mathsf T}A_9\mathbf1=8.
```

Thus both signs of the seed ground energy really occur; this matters in the
tensor argument below.

## 2. Fixed-size replication cannot make selector excess leading

Take `t` copies of the seed `3+6` local blocks, with arbitrary signs between
different copies, and replace each order-three and order-six diagonal block by
a local minimizer.  If `C_t` denotes everything outside those `2t` blocks and
`K` the direct sum of all replacements, then

```math
Q(K)\le t(q_3+q_6)=16t
```

and the triangle/reverse-triangle inequality gives the exact uniform bound

```math
\boxed{
\left|Q(C_t+K)-Q(C_t)\right|\le16t.
}
```

Since the total order is `N=9t`,

```math
\frac{16t}{N^{3/2}}=\frac{16}{27\sqrt t}\longrightarrow0.
```

This does not say that the four-unit seed excess adds from copy to copy (it
need not).  It says something stronger for the asymptotic question: **every**
selector effect carried by fixed `3+6` blocks is only `O(N)`, independently of
how the copies are coupled.  A leading cross-only overshoot could still come
from the coupling, but it would not be an amplification of the finite selector
wall.

## 3. Regular vertex blow-ups lose global minimality

Replace each seed vertex by `k` clones and keep every inter-fibre sign constant.
Across the enlarged `3k+6k` split the cross block is `B\otimes J_k`.  Its
Boolean bilinear optimization is exactly the optimization over fibre
magnetizations in `[-k,k]`.  Bilinearity puts the extrema at the endpoints, so

```math
\boxed{
Q(C_9\otimes J_k)=k^2Q(C_9)=24k^2.
}
```

For any choice whatsoever of the two within-shore signings, flipping the
relative shore sign proves that the completed signing has norm at least this
cross norm.  It therefore has order-`N^2` norm rather than order-`N^{3/2}`
norm.

The exclusion is already finite for every nontrivial `k`.  For `k=2`, the
order-18 Paley conference matrix gives

```math
q_{18}\le18\sqrt{17}<96.
```

For `k=3`, restrict the order-30 Paley conference matrix to 27 vertices:

```math
q_{27}\le27\sqrt{29}<216.
```

For `k\ge4`, the random-sign upper bound gives

```math
q_{9k}
\le2\sqrt{(9k)(9k-1)(9k+2)\log2}
<24k^2.
```

The last inequality follows from
`4\sqrt k>\sqrt{\log2}(9+2/k)`, whose left side increases and right side
decreases and which already holds at `k=4`.

For reference, if all within-fibre edges are positive, the entire blow-up has
the exact norm

```math
\boxed{Q=33k^2-9k.}
```

Indeed its energy at fibre magnetization vector `m` is
`m^{\mathsf T}A_9m+\sum_i(m_i^2-k)`.  Multi-affinity gives
`|m^{\mathsf T}A_9m|\le24k^2`; the positive upper bound is attained by a
positive seed ground with every fibre aligned, while the negative bound is no
larger for `k\ge2` (and `k=1` is the seed itself).

## 4. Lexicographically coupling intact copies also loses minimality

Let `R` be any order-`t` signing and form

```math
L_t=R\otimes J_9+I_t\otimes A_9.
```

Choose the same all-one seed state in every copy, multiplied by an outer spin
`s`.  Since its magnetization is nine and its seed energy is eight,

```math
(s\otimes\mathbf1)^{\mathsf T}L_t(s\otimes\mathbf1)
=81s^{\mathsf T}Rs+8t.
```

Consequently

```math
\boxed{Q(L_t)\ge81Q(R)-8t\ge81q_t-8t.}
```

This construction is not globally minimizing for any `t\ge2`.  For
`2\le t\le9`, inserting the exact values

```math
(q_2,\ldots,q_9)=(2,6,8,8,10,18,20,24)
```

already makes the lower bound exceed the random-sign upper bound for
`q_{9t}`.  For `t\ge10`, finite Gaussian rounding and
`\arcsin u\ge u` give

```math
q_t\ge\frac{2}{\pi}t\sqrt{t-1}.
```

After bounding the order-`9t` random upper bound slightly upward and dividing
by `t^{3/2}`, it is enough that

```math
\frac{162}{\pi}\sqrt{1-\frac1t}-\frac8{\sqrt t}
>
54\sqrt{\log2}+\frac{12\sqrt{\log2}}t.
```

The difference is increasing in `t` and is greater than `0.43` at `t=10`.
Thus uniform inter-copy coupling has an unavoidable leading constant inflated
by the squared seed order; it cannot preserve global minimality.

## 5. Every diagonal-completion tensor square fails

There is a less naive tensor operation that must be checked separately.  For
an arbitrary diagonal sign matrix `D`, put

```math
P=A_9+D,
\qquad \delta=\operatorname{tr}D,
```

so `P` is a full symmetric `\{\pm1\}` matrix, and define the order-81 signing

```math
T_D=P\otimes P-\operatorname{diag}(P\otimes P).
```

This includes the usual star product `(A_9+I)^{\otimes2}-I`, but also all
balanced diagonal completions.  None is globally minimizing.

A symmetric Paley conference matrix of order 82 exists.  A principal
order-81 submatrix therefore proves

```math
q_{81}\le81\sqrt{81}=729.
```

Because `\delta` is odd, split into two cases.

If `|\delta|\ge5`, choose an `A_9` ground of energy
`24\operatorname{sgn}\delta` in both tensor factors.  The resulting product
state has exact `T_D` energy

```math
(24\operatorname{sgn}\delta+\delta)^2-\delta^2
=576+48|\delta|\ge816>729.
```

If `|\delta|\le3`, use the Boolean state `\operatorname{vec}(P)`.  Symmetry
gives

```math
\operatorname{vec}(P)^{\mathsf T}T_D\operatorname{vec}(P)
=\operatorname{tr}(P^4)-\delta^2.
```

Every diagonal entry of `P^2` is nine.  Every off-diagonal entry is the inner
product of two sign vectors of odd length nine, hence is a nonzero odd integer.
Therefore

```math
\operatorname{tr}(P^4)=\|P^2\|_F^2
\ge9\cdot9^2+9\cdot8=801,
```

and the displayed witness is at least `801-9=792>729`.  This proves the
no-go for **all 512 diagonal completions**, without trusting a heuristic
search.  Exact enumeration happens to sharpen the smallest witness supplied
by these two formulas to 896, but that sharpening is not needed.

## 6. Scope of the falsification

The fixed-block, regular-blow-up, uniform lexicographic, and same-seed star
tensor constructions cannot amplify the A9 wall while retaining global
minimality.  The conclusions do **not** cover a mixed product of the seed with
a growing conference/Hadamard-like factor, a nonuniform gadget substitution,
or a construction whose order-`3m` and order-`6m` local minimizers have a new
tensor-compatible description.  Those variants remain open.  In particular,
the six-state A9 certificate cannot simply be asserted to tensor: minimizer
sets at the enlarged block orders are not products of the order-three and
order-six minimizer sets.
