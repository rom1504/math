# Wave 19 memo: indexed child optimizers have an exact cross-puncture cocycle

## Status

The identities and inequalities below are **proved exactly**.  The `A_8` and
`A_9` audits are exhaustive integer computations over all projective child
spins and every optimal extension.  They are implemented in
`tmp/puncture_potential_r19.py`.  No asymptotic puncture-stability theorem is
proved.

## Indexed extension theorem

Use the one-copy normalization

```math
M(A)=\max_{q\in\mathcal C_n}\langle A,q\rangle,
\qquad
\langle A,q\rangle=\sum_{u<v}a_{uv}q_{uv},
```

where an oriented augmented-cut word has the form
$`q_{uv}=t x_u x_v`$.  Thus the ledger's doubled norm is $`Q(A)=2M(A)`$.
Let $`A`$ be an exact order-$`n`$ minimizer, write $`M=M_n`$, and put

```math
B_j=A[-j],\qquad
d_j=M-M(B_j),\qquad
e_j=M(B_j)-M_{n-1},\qquad
\delta_n=M_n-M_{n-1}.
```

Hence $`d_j+e_j=\delta_n`$.  For each $`i`$, choose an arbitrary oriented
child ground

```math
c_i=q^{t_i,x^{(i)}}\in\mathcal C_{n-1},
\qquad
\langle B_i,c_i\rangle=M(B_i)=M-d_i.
```

Its signed deleted-star field is

```math
h_i=t_i\sum_{j\ne i}a_{ij}x_j^{(i)}.
```

Choose the missing spin so that the restored star contributes $`|h_i|`$,
and call the resulting full word $`q_i`$.  Parent optimality gives

```math
\langle A,q_i\rangle=M-d_i+|h_i|\le M.
```

Consequently the extension slack

```math
s_i=d_i-|h_i|
```

obeys $`0\le s_i\le d_i`$, and

```math
\boxed{
\langle A,q_i\rangle=M-s_i,
\qquad
r_i(q_i)=d_i-s_i,
}
\tag{R19.1}
```

where

```math
r_j(q)=\sum_{k\ne j}a_{jk}q_{jk}
```

is the incident contribution at $`j`$.  The second equality in (R19.1) is
the promised normalization audit: the optimal choice of the missing spin
makes the incident contribution exactly $`|h_i|`$, with no factor of two.
In doubled normalization all quantities in (R19.1) are multiplied by two.

For every ordered pair define the **cross-puncture deficit**

```math
u_{ij}
=M(B_j)-\langle A[-j],q_i[-j]\rangle.
\tag{R19.2}
```

Since $`q_i[-j]`$ is an admissible oriented child word, $`u_{ij}\ge0`$.
Using

```math
\langle A[-j],q_i[-j]\rangle
=\langle A,q_i\rangle-r_j(q_i)
=M-s_i-r_j(q_i)
```

gives the pointwise identity

```math
\boxed{
u_{ij}=s_i+r_j(q_i)-d_j\ge0,
\qquad u_{ii}=0.
}
\tag{R19.3}
```

Thus every closest child word extends to a parent near-endpoint with its
own coordinate distinguished by zero cross-puncture deficit.  More
importantly, summing (R19.3) does not create a new inequality.  Every edge
of $`q_i`$ occurs in two incident fields, so

```math
\sum_jr_j(q_i)=2\langle A,q_i\rangle=2(M-s_i).
```

Therefore, for **every** $`i`$ and every choice of the child ground,

```math
\boxed{
\sum_{j=1}^n u_{ij}
=2M-\sum_{j=1}^n d_j+(n-2)s_i.
}
\tag{R19.4}
```

Writing $`E=\sum_j e_j`$ and using
$`\sum_jd_j=n\delta_n-E`$ yields the equivalent defect identity

```math
\boxed{
E
=n\delta_n-2M+
\sum_j u_{ij}-(n-2)s_i
\qquad\text{for every }i.
}
\tag{R19.5}
```

In particular,

```math
\boxed{
\overline e_n
=\delta_n-\frac{2M_n}{n}
+\frac1n\left(\sum_j u_{ij}-(n-2)s_i\right).
}
\tag{R19.6}
```

There is also a weighted form.  For any probability vector $`p`$ on the
deleted coordinate,

```math
\boxed{
\sum_jp_je_j
=\delta_n-s_i-\sum_jp_jr_j(q_i)+\sum_jp_ju_{ij}.
}
\tag{R19.7}
```

Indeed (R19.3) says pointwise that
$`r_j(q_i)-u_{ij}=d_j-s_i`$.  Thus a linear rule that credits the visible
field and debits its cross-puncture error collapses exactly to the original
deletion decrement; it is not a new selected-terminal potential.

## A pairwise constraint, and its scale

The cross-puncture matrix does retain nonlinear information.  Since
$`q_j[-j]`$ is a ground of $`B_j`$,

```math
u_{ij}=\langle A[-j],q_j[-j]-q_i[-j]\rangle.
```

For $`i\ne j`$, common edges cancel after symmetrization:

```math
\boxed{
\begin{aligned}
u_{ij}+u_{ji}
={}&\sum_{k\ne i,j}a_{ik}(q_{j,ik}-q_{i,ik})\\
&-\sum_{k\ne i,j}a_{jk}(q_{j,jk}-q_{i,jk}).
\end{aligned}
}
\tag{R19.8}
```

Both deficits are nonnegative and every word difference is in
$`\{0,\pm2\}`$, hence

```math
\boxed{0\le u_{ij}+u_{ji}\le4(n-2).}
\tag{R19.9}
```

Summing (R19.9) gives

```math
\sum_{i,j}u_{ij}\le2n(n-1)(n-2).
\tag{R19.10}
```

This universal bound is one square-root scale too weak for competitive
signings: (R19.4) naturally has total mass $`\Theta(n^{5/2})`$, while
(R19.10) allows $`\Theta(n^3)`$.  A useful refinement must therefore use
the cut/complement support, tails, or cycles of $`u`$, not just the
two-coordinate cancellation.

## Exact finite audits

For the balanced exact order-eight minimizer (10.445), exhaustive
enumeration gives

```math
M=10,\qquad M(B_i)=9,\qquad d_i=1,\qquad e_i=0
\quad(i=1,\ldots,8).
```

Every child has three projective optimal extensions.  For all 24 choices,
$`s_i=0`$, $`u_{ij}\in\{0,2,4\}`$, and

```math
\sum_j u_{ij}=12=2M-\sum_jd_j.
```

For the order-nine minimizer (10.298), exhaustive enumeration gives

```math
M=12,\qquad M(B_i)=12,\qquad d_i=0,\qquad e_i=2
\quad(i=1,\ldots,9).
```

Across all 28 projective optimal extensions, $`s_i=0`$,
$`u_{ij}\in\{0,2,4,6,8\}`$, and every row obeys

```math
\boxed{\sum_j u_{ij}=24=2M.}
\tag{R19.11}
```

There are zero off-diagonal entries as well as rows with maximal pair sum
$`u_{ij}+u_{ji}=16`$.  Thus isolated exact cross-compatibilities can coexist
with the full aggregate obstruction.  This passes both the no-child wall
and the flat two-level warning: the theorem does not infer descent from a
zero entry.

The verifier output is:

```text
A8 ... d=[1,...,1] s=[0] u_values=[0,2,4] row_sums=[12]
A9 ... d=[0,...,0] s=[0] u_values=[0,2,4,6,8] row_sums=[24]
```

## Consequence for the leading route

Equations (R19.4)--(R19.7) close a precise class of proposed arguments:
**linear summation or averaging of the indexed cross-puncture deficits and
extension slacks is an exact cocycle**.  It merely rewrites
$`\overline e_n`$ and cannot prove $`\overline e_n=O(1)`$ or construct the
potential (10.531).  This statement does not rule out nonlinear use of the
matrix $`(u_{ij})`$.

The additional lemma now needed can be stated sharply: exploit that each
$`q_i/q_j`$ is a cut or complemented cut to control a nonlinear feature of
$`u`$ strong enough that (R19.6) has a summable positive part.  Viable forms
include a tail bound coupling large $`u_{ij}`$ to negative credit, a
low-deficit directed cycle with a genuine multivertex exchange consequence,
or a deterministic landing rule for which the weighted quantity in
(R19.7) is small.  Any such result must use more than row totals and must
survive (R19.11).
