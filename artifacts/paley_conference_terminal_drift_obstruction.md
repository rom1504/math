# Square-field Paley conferences obstruct a unique-zero terminal drift reduction

## Status

**Verified scalable route obstruction.**  Square-field Paley conference
signings are strict one-edge local minima of the Boolean cap at infinitely
many orders.  Consequently, the proposed uniform terminal coset-graph drift
law with a unique zero would force that zero to be `1/2` and would prove the
full sharp lower constant, not merely existence of an unspecified limit.

This does not falsify the drift law.  It shows that the law is not a weaker
route to the convergence objective.

## Setup

Let `r` be an odd prime power, let `F=GF(r^2)`, and let `K=GF(r)` be its
subfield.  The symmetric Paley conference matrix `C` is indexed by
`F union {infinity}` and has order

```math
n=r^2+1.
```

Its entries are

```math
C_{infinity,z}=1,
\qquad
C_{z,w}=\chi(z-w)\quad(z\ne w),
\qquad
C_{z,z}=0,
```

where `chi` is the quadratic character of `F`.  It satisfies `C^2=r^2 I`.
Write

```math
Q(C)=\max_{x\in\{\pm1\}^n}{1\over2}|x^{\mathsf T}Cx|.
```

For a signing root `U`, put

```math
R(U)=d(U,\mathcal C_n^+),
\qquad
z(U)={\binom n2-2R(U)\over n^{3/2}},
```

and let `b(U)` be the number of single-edge flips that increase `R(U)` by
one, equivalently decrease `Q` by two.

## Theorem

For every such square-field Paley conference matrix,

```math
\boxed{Q(C)={rn\over2}.}                                      \tag{1}
```

Moreover, if `C^(e)` is obtained by flipping any one off-diagonal edge, then

```math
\boxed{Q(C^{(e)})=Q(C)+2.}                                   \tag{2}
```

Thus its augmented-cut-code coset satisfies

```math
\boxed{b(U_C)=0,
\qquad
z(U_C)={r\over2\sqrt{r^2+1}}\longrightarrow{1\over2}.}       \tag{3}
```

## Proof

Choose `t in F setminus K`.  If `f:K->{+-1}` satisfies
`sum_(b in K) f(b)=1`, define

```math
x_{infinity}=1,
\qquad
x_{a+bt}=f(b).
```

The standard subfield character sum gives `Cx=rx`.  Hence the spectral
identity `C^2=r^2I` gives both the lower and upper bounds in (1).

Explicitly, because every element of `K^*` is a square in `F`,

```math
\sum_{a\in K}\chi(a+ht)=
\begin{cases}r-1,&h=0,\\-1,&h\ne0.\end{cases}
```

At a finite coordinate in fibre `b`, the border term and these fibre sums
give

```math
1+(r-1)f(b)-\sum_{d\ne b}f(d)=rf(b),
```

while the infinity coordinate is
`r sum_(b in K)f(b)=r`.

We need the following strengthened form of the same construction: for every
edge `{i,j}`, there is a Boolean `r`-eigenvector `x` such that

```math
C_{ij}x_ix_j=-1.                                             \tag{4}
```

For a finite edge `{z,w}`, choose a square `u in F^*` such that
`u(z-w) notin K`.  Such a square exists because the square subgroup has
`(r^2-1)/2` elements, whereas at most `r-1` choices send `z-w` into `K^*`.
Multiplication by `u` fixes the Paley matrix because `chi(u)=1`.  Write

```math
uz=a+bt,
\qquad
uw=c+dt;
```

then `b!=d`.  Choose `f` with sum one and

```math
f(b)f(d)=-\chi(z-w).
```

This is possible for every odd `r>=3`: prescribe two positive values when
the required product is positive, or one positive and one negative value
when it is negative, and fill the remaining values so that there are
`(r+1)/2` positives.  Pulling the preceding Boolean eigenvector back through
the multiplication permutation gives (4).  For an edge `{infinity,z}`, choose
`f` with sum one and with the relevant fibre value equal to `-1`.

Now flip edge `{i,j}` and evaluate at its vector from (4).  In the one-copy
normalization the energy increases by exactly two:

```math
{1\over2}x^{\mathsf T}C^{(e)}x
={rn\over2}-2C_{ij}x_ix_j
={rn\over2}+2.
```

A single edge flip changes every Boolean energy by at most two, so this lower
bound is exact and proves (2).  Since

```math
R(U)={\binom n2-Q(U)\over2},
```

every edge flip is inward in the coset graph.  This proves (3).

## Consequence for the terminal drift proposal

Suppose there were a compact interval `I` containing `[0.33,0.51]`, a
continuous `beta:I->[0,1]` with unique zero `c`, and `epsilon_n->0` such that

```math
\left|{b(U)\over\binom n2}-\beta(z(U))\right|\le\epsilon_n   \tag{5}
```

for every sufficiently large order and every root with `z(U) in I`.
Applying (5) to (3) along square-field Paley orders gives

```math
\beta(1/2)=0,
```

so uniqueness forces `c=1/2`.  Applying (5) to deepest cosets, which also
have `b=0`, then gives

```math
\boxed{{M_n\over n^{3/2}}\longrightarrow{1\over2}.}          \tag{6}
```

The terminal statistic `(z,b/binom(n,2))` is genuinely smaller than a full
coset histogram, but its proposed uniform unique-zero theorem contains the
entire sharp-constant burden.  It therefore fails the project's criterion
that an execution lemma remove an obligation rather than replace it by an
equally strong or stronger one.

## Finite checks

The existing exact records independently verify the theorem at `r=3` and
`r=5`:

- at order `10`, `Q(C)=15` and all `45` edge flips have cap `17`;
- at order `26`, `Q(C)=65`; the saved `260` positive and negative projective
  extremizers collectively contain a witness of (4) for every one of the
  `325` edges.

The source and compact audit record are
`computations/audit_paley_edge_traps.py` and
`computations/results/paley_conference_edge_trap_audit.json`; they read the
previously certified exhaustive records
`computations/results/conference_order10_gf9.json` and
`computations/results/conference_order26_gf25.json`.
