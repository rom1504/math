# Per-signing projection barriers for the exposed-shell restriction target

**Status.**  Two rigorous proof-class obstructions.  They show that neither
a scalar parent-plus-residual decomposition nor concentration based only on
the operator norm forced by low cap can prove the exposed-layer
restriction-shadow statement (ERSR).  They do **not** falsify ERSR: an
argument using the joint incidence of restrictions across the uniform
exposed shell remains possible.

Throughout,

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|.
```

## 1. Every scalar coordinate projection has a leading residual

We first record an elementary complete-support lower bound.

### Lemma 1 (a universal cap lower bound)

If `B` is a hollow signing of order `t>=2`, then

```math
\boxed{Q(B)\ge c_0t^{3/2},\qquad c_0={1\over6}.}       \tag{PB.1}
```

**Proof.**  Partition the vertices into `I,J`, where
`a=|I|=floor(t/2)` and `b=|J|=ceil(t/2)`.  Choose the spins on `I`
independently and uniformly.  For `j in J`, put

```math
X_j=\sum_{i\in I}b_{ij}x_i.
```

The elementary interpolation inequality

```math
\|X_j\|_2\le \|X_j\|_1^{1/3}\|X_j\|_4^{2/3}
```

and

```math
\mathbb EX_j^2=a,
\qquad \mathbb EX_j^4=3a^2-2a\le3a^2
```

give `E|X_j|>=sqrt(a/3)`.  Hence some choice of the `I` spins has

```math
\sum_{j\in J}|X_j|\ge b\sqrt{a/3}.
```

Choose `x_j=sign(X_j)` on `J`.  Flipping all `J` spins preserves the two
internal block energies and reverses the cross energy.  Therefore one of
the two choices has absolute total energy at least the cross energy.  Since
`a>=t/3` and `b>=t/2`, this is at least `t^(3/2)/6`. `square`

Let `S subset [N]`, `|S|=m`, and let `A^S` be the `N` by `N` matrix equal
to `A` on `S times S` and zero elsewhere.  Thus `Q(A^S)=Q(A[S])`.

### Proposition 2 (scalar projection-width barrier)

For every hollow signing `A`, every `S`, and every real `lambda`,

```math
\boxed{
Q(A^S-\lambda A)
\ge c_0\max\left\{
 |1-\lambda|m^{3/2},
 |\lambda|(N-m)^{3/2}
\right\}.}                                           \tag{PB.2}
```

Consequently,

```math
\boxed{
\inf_{\lambda\in\mathbb R}Q(A^S-\lambda A)
\ge c_0,{m^{3/2}(N-m)^{3/2}
 \over m^{3/2}+(N-m)^{3/2}}.}                        \tag{PB.3}
```

In particular, uniformly for `m/N in [1/3,2/3]`, the right side is at
least

```math
{1\over36\sqrt3}N^{3/2}.                              \tag{PB.4}
```

**Proof.**  The principal restriction of the residual to `S` is
`(1-lambda)A[S]`; its restriction to `S^c` is `-lambda A[S^c]`.
The cap of a principal submatrix never exceeds the cap of the full matrix:
extend a maximizing spin randomly over the complementary coordinates and
use that the expected cross and complementary energies vanish.  Apply
Lemma 1 to both principal blocks.  Minimizing
`max{a|1-lambda|,b|lambda|}` with
`a=m^(3/2)` and `b=(N-m)^(3/2)` gives `ab/(a+b)`.  Finally,
`ab/(a+b)>=min(a,b)/2` proves (PB.4). `square`

The exact scalar decomposition is

```math
A^S=\lambda A+(A^S-\lambda A).                        \tag{PB.5}
```

Any proof which pays its two terms separately obtains at best

```math
Q(A[S])\le |\lambda|Q(A)+Q(A^S-\lambda A).           \tag{PB.6}
```

Proposition 2 says that the second term is never `o(N^(3/2))` on a
comparable restriction, regardless of how small `Q(A)` is.  Thus this
proof class cannot establish the cap error
`epsilon_N m^(3/2)=o(N^(3/2))` required by ERSR.  After applying
`u -> u^(2/3)`, a fixed leading cap payment is a `Theta(N)` payment in the
`b=Q^(2/3)` scale, so it cannot yield a sublinear cross-order defect or the
pressure recurrence deduced from ERSR.  Cancellation between the parent
and residual must be retained jointly.

## 2. Low cap does not force a usable spectral restriction theorem

The standard polarization/interpolation inequality gives, in the present
one-copy normalization,

```math
\|A\|_{op}^2\le4Q(A).                                \tag{PB.7}
```

Thus `Q(A)=O(N^(3/2))` forces only
`||A||_op=O(N^(3/4))`.  This exponent cannot be improved from even
asymptotic near-minimality.

### Proposition 3 (spectrally spiked exact-sign near-minimizers)

Let `M_N=min_A Q(A)`.  There are hollow signings `A_N` such that

```math
\boxed{
M_N\le Q(A_N)\le M_N+o(N^{3/2}),
\qquad
\|A_N\|_{op}\ge {N^{3/4}\over\sqrt{\log N}}(1-o(1)).} \tag{PB.8}
```

**Proof.**  Put `s=floor(sqrt(N)/log N)` and `r=N-s`.  On a set `U` of
size `r`, use an exact order-`r` minimizer.  Put arbitrary signs inside the
remaining set `T`, and put `+1` on every `T`--`U` edge.  The triangle
inequality gives

```math
Q(A_N)\le M_r+sr+{s\choose2}.
```

Principal deletion proves `M_r<=M_N`, while
`sr+binom(s,2)=o(N^(3/2))`; the opposite inequality
`M_N<=Q(A_N)` is definitional.  Finally, the rectangular compression
`P_TA_NP_U` is the all-one `s` by `r` matrix.  Its operator norm is
`sqrt(sr)`, and compression cannot increase operator norm. `square`

For completeness, fix a spin `x`, switch `A` by `x`, and reveal a
Bernoulli-`p` vertex set.  Writing `delta=p1+xi` gives

```math
H_{A[S]}(x_S)
=p^2H_A(x)+p\,\xi^T B1+{1\over2}\xi^TB\xi.           \tag{PB.9}
```

Subgaussian concentration for the linear term and Hanson--Wright for the
quadratic term imply

```math
\Pr\{|H_{A[S]}(x_S)-p^2H_A(x)|>u\}
\le C_p\exp\left[-c_p\min\left{
 {u^2\over N\|A\|_{op}^2},
 {u^2\over N^2},
 {u\over\|A\|_{op}}
\right\}\right].                                    \tag{PB.10}
```

Conditioning on `|S|=m=pN+O(1)` costs only a polynomial factor and does
not change the exponential scale.  At the ERSR deviation
`u=Theta(N^(3/2))`, the bound forced by (PB.7) has exponent only
`Theta(sqrt(N))`.  On the near-minimizer family in Proposition 3, the
minimum exponent in (PB.10) is at most
`N^(1/2+o(1))`.  This cannot be union-bounded over
`exp(Theta(N))` cut witnesses; degree-two hypercontractivity gives the same
subexponential obstruction.  Vertex bounded differences is weaker still:
one coordinate replacement can change the restricted cap by `Theta(N)`,
so its natural fluctuation scale is `Theta(N^(3/2))`.

## 3. Exact scope

These propositions prove that the following inputs are insufficient **for
the listed signing-by-signing scalar-projection and concentration proof
classes**:

1. low cap or even `o(N^(3/2))` additive near-minimality;
2. complete support and uniformly small individual edge influences;
3. a scalar conditional-expectation projection followed by a separately
   paid residual;
4. operator-norm/Hanson--Wright based only on (PB.7), hypercontractive, or
   bounded-difference concentration applied signing by signing.

They do not show that a typical restriction of a uniform exposed-shell
signing is bad.  The planted family in Proposition 3 may have negligible
mass in that shell, Proposition 2 leaves open cancellation visible only
after averaging jointly over signings and restrictions, and an
operator-restriction theorem using additional exposed-shell structure is
not ruled out.  ERSR can therefore still hold, but the audited generic
tools do not prove it; a live proof would need genuinely shell-level
restriction incidence, multiplicity, or another input absent from these
per-signing bounds.
