# Independent audit: cross-Gram macroscopic dimension

Audited file:
[`cross_gram_macroscopic_dimension.md`](cross_gram_macroscopic_dimension.md)

**Verdict: PASS, with one Lipschitz-convention clarification.**  The PSD cube
lemma remains valid at singular centres, the constants in CG.2 are correct,
the tensor-Walsh cube is a genuine Boolean realization, and CG.25 follows
with the displayed constant from its explicitly assumed full-bit response
bound CG.23.  The result proves an `O(p)` law only for a full affine cube of
independently toggleable raw coordinates.  It does not prove `O(p)` metric
entropy for arbitrary nonlinear cross-Gram codebooks; the draft states this
limitation accurately.

## 1. CG.1 at a positive-definite centre

Let `M=K_0^{-1}` and

```math
T_\sigma=K_0^{-1/2}S_\sigma K_0^{-1/2}.
```

Because both sign words `sigma` and `-sigma` occur, positivity gives

```math
I+T_\sigma\succeq0,
\qquad I-T_\sigma\succeq0.
```

Thus `T_sigma` is a symmetric contraction and
`tr(T_sigma^2)<=p`.  Cyclicity of trace gives

```math
\operatorname{tr}T_\sigma^2
=\operatorname{tr}(MS_\sigma MS_\sigma).
```

Averaging the independent signs removes the terms indexed by distinct
edges.  Direct multiplication for `e={i,j}` gives

```math
\operatorname{tr}(MD_eMD_e)
=2(M_{ii}M_{jj}+M_{ij}^2).                           \tag{A.CG.1}
```

For a positive-definite correlation matrix, the Schur complement formula
indeed gives `M_ii>=1`.  Consequently

```math
p\ge2\sum_e\eta_e^2(M_{ii}M_{jj}+M_{ij}^2)
 \ge2\sum_e\eta_e^2,
```

which is exactly CG.3.  The factor `1/2` is sharp already for a disjoint
matching with centre `I` and amplitudes one.

## 2. Singular normalization is valid

The potentially delicate step is sound.  Put

```math
d_i=K_{0,ii}+\varepsilon,
\qquad D=\operatorname{diag}(d_i).
```

For every sign word,

```math
D^{-1/2}(K_\sigma+\varepsilon I)D^{-1/2}\succ0
```

has centre with diagonal one.  On edge `{i,j}` its amplitude is exactly

```math
\eta_e'=\frac{\eta_e}{\sqrt{d_id_j}}.               \tag{A.CG.2}
```

Since `K_(0,ii)<=1`, both `d_i,d_j<=1+epsilon`, and hence

```math
|\eta_e'|\ge\frac{|\eta_e|}{1+\varepsilon}.
```

Applying the positive-definite result yields

```math
\frac1{(1+\varepsilon)^2}\sum_e\eta_e^2\le p/2.
```

Letting `epsilon` decrease to zero proves the singular case.  There is no
division-by-zero assumption hidden here.  In particular, if a diagonal
entry of the original PSD centre is zero, the limiting inequality forces
every independently toggleable incident amplitude to vanish, as it should.

## 3. Constants in CG.2

The sector matrices

```math
K^\pm=(G\pm R)/2
```

are Gram matrices of the projections `(I+-J)w_i/2`, so they are PSD and
their diagonals lie in `[0,1]`.  Fixed one-port self data makes those
diagonals independent of the cube word.  On coordinate `a`, their affine
half-amplitudes are

```math
\eta_a^+=(g_a+r_a)/2,
\qquad
\eta_a^-=(g_a-r_a)/2.
```

CG.1 gives

```math
\sum_a(\eta_a^+)^2\le p/2,
\qquad
\sum_a(\eta_a^-)^2\le p/2.
```

Since

```math
(\eta_a^+)^2+(\eta_a^-)^2=(g_a^2+r_a^2)/2,
```

the displayed conclusion is precisely

```math
\sum_a(g_a^2+r_a^2)\le2p.
```

Thus neither a factor two nor a sector has been lost in CG.12--CG.14.  The
constant `2p` is not advertised as optimal; the theorem's claimed linear
order is what follows from the two separate PSD budgets.

One terminology point is worth keeping explicit.  In CG.11, changing a bit
from `sigma_a=-1` to `sigma_a=+1` changes the raw pair by

```math
2(g_a,r_a),                                           \tag{A.CG.3}
```

whereas `sqrt(g_a^2+r_a^2)` is its **half-amplitude**.  CG.13 uses the latter
convention consistently.  Calling it the full raw coordinate change would
lose a factor two.

## 4. Tensor-Walsh realization

Let `H_16 1=4 1` and `H_16 v_0=4v_0`, with `1` and `v_0` orthogonal Boolean
vectors.  Tensoring `j` choices gives `2^j=n^(1/4)` mutually orthogonal
Boolean vectors, all in the `+sqrt(n)` eigenspace of
`H_16^(tensor j)`.  Allocating two distinct words to every matching edge
therefore works for every even `p<=n^(1/4)`.

For bit zero the two ports are orthogonal; for bit one the second repeats
the first.  All cross terms between different allocated pairs vanish.  In
`+-1` bit coordinates the centre on each matching edge is `1/2` and
`g_a=r_a=1/2`.  Hence this is genuinely an affine cube with fixed diagonal
self data and `h=p/2` constant-amplitude coordinates.

The local response calculation also scales beyond the displayed base
check.  For any orthogonal Boolean `+sqrt(n)` eigenvectors `a,b`, every
signed sum `+-a+-b` has norm `sqrt(2n)`.  Therefore the repeated-pair cap is
exactly `5rn/2`, while the orthogonal-pair cap is at most
`(1/2+sqrt(2))rn`.  The difference is at least
`(2-sqrt(2))rn`, not an equality claim.  This verifies CG.20.

## 5. Total-order calculation and its exact scope

With `p` shores of width `r=sqrt(n)`, the vertex count is

```math
N=n+p\sqrt n.
```

Thus CG.22 is the exact conversion of an old-block gap `c n^(3/2)` into
total-order units.  In particular, `p<=n^(1/4)` leaves a nonvanishing ratio;
`p/sqrt(n)->infinity` kills it.

Now assume CG.23 literally: a **full flip** of coordinate `a` changes the
declared response by at most

```math
Lrn\sqrt{g_a^2+r_a^2}.
```

Separation by `epsilon N^(3/2)` forces CG.24.  Applying CG.12 to `h`
coordinates then gives

```math
h{\epsilon^2\over L^2}\left({N\over n}\right)^3
\le2p,
```

which is exactly CG.25.  The `O(p)` conclusion is therefore correct under
CG.23.

There is one convention repair to make in the explanatory sentence after
CG.23.  If a function of `(G,R)` is `L_0`-Lipschitz in the usual Euclidean
metric, (A.CG.3) gives the full-flip bound with `L=2L_0`, not `L=L_0`.
Equivalently, with the standard Lipschitz constant the last estimate reads

```math
h\le {8L_0^2\over\epsilon^2}
       p\left({n\over N}\right)^3.                  \tag{A.CG.4}
```

This is only a factor-four change in CG.25 and does not affect the linear
law.  CG.25 itself is correct because its `L` is defined directly by the
full-bit response bound CG.23.

More importantly, CG.25 is conditional, not a metric-entropy theorem for
the whole Gram pair.  It does not cover:

- a nonlinear code rather than a full affine bit cube;
- bits distributed coherently over many small entries;
- a high-arity decoder whose Lipschitz scale grows with `p`;
- a collective exact Boolean response not controlled by the spherical
  feature metric.

The draft explicitly preserves all four escape routes.  Thus the defensible
conclusion is “`O(p)` independently toggleable fixed-amplitude affine
coordinates,” not “every cross-Gram response quotient has `O(p)` bits.”

## 6. Verifier and final classification

Running

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_cross_gram_macroscopic_dimension.py
```

returns `PASS`.  It checks a finite PSD matching cube, four tensor-Walsh
states at order `256`, the exact order-`16` local caps `(160,96)`, and the
order-normalization ratios.  It does not numerically establish CG.1 or CG.2;
those are established by the algebra above.

| Item | Verdict | Repair |
|---|---|---|
| CG.1 positive-definite proof | PASS | none |
| CG.1 singular normalization | PASS | none |
| CG.2 factor and constants | PASS | none |
| Boolean tensor-Walsh cube | PASS | none |
| local SA.4 exposure | PASS | retain `>=` wording |
| CG.21--CG.22 order count | PASS | none |
| CG.25 from CG.23 | PASS | none |
| standard Lipschitz interpretation | PASS after clarification | use `L=2L_0` |
| unrestricted `O(p)` entropy law | not claimed/proved | retain current scope |

The draft's substantive conclusion survives the audit: a quadratic exact
table does not contain quadratically many independently toggleable
macroscopic affine coordinates.  Whether coherent small entries generate
larger collective response entropy remains genuinely open.
