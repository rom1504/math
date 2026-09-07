# Wave 50C: the Bellman pressure is a switching-orbit circulation

## Status

- **Verified exact obstruction:** for every selector law, the canonical
  fractional vector `p_e=(1-w_e)/4` in `(10.1240)` is an explicit convex
  combination of cut-switching symmetries of the original signing.
- **Verified exact identity:** selector by selector, the full state-dependent
  Bellman margin is precisely the expected parent slack along that switching
  orbit.  Its positivity is therefore compatible with every integral
  realization having minimum margin zero.
- **Verified circulation:** a parent ground migrates bijectively with each
  switching outcome.  Because the switch is outside the sampled selector,
  its child restriction is unchanged.  Thus even a bad maximal block and its
  joint-gap budget are respected by the migrating witness.
- **Falsified, with scope:** no proof using only the canonical fractional
  point, the all-state identity, and an integral realization of those
  marginals can give descent.  This remains true if the law uses the complete
  maximal-selector family and retains all state-dependent budgets.
- **Open:** tight principal decomposition itself is not falsified.  A proof
  from the global all-maximal-bad premise would now need a genuinely
  non-switching certificate or an additional incompatibility not encoded in
  `(10.1240)`.

The exact checker and output are `tmp/global_all_bad_r50_check.py` and
`tmp/global_all_bad_r50.out`.

## 1. Conditional switching law and exact edge marginals

Fix an arbitrary signing `A`, an `m`-selector `S`, and write `T=S^c`.  For a
vertex set `U`, let

```math
\delta(U)=\{ij:\mathbf 1_U(i)\ne\mathbf 1_U(j)\}
```

be its cut.  Define a random edge-flip set `F_S` as follows:

1. with probability `1/2`, take `U=emptyset`;
2. with probability `1/2`, take `U` uniformly from all subsets of `T`;
3. set `F_S=delta(U)`.

For `e=ij`, there are three cases.  If `i,j in S`, then `e` is never
flipped.  If exactly one endpoint lies in `T`, its membership bit is fair in
the uniform-subset branch, so the conditional flip probability is
`(1/2)(1/2)=1/4`.  If both endpoints lie in `T`, their independent fair bits
differ with probability `1/2`, so the answer is again `1/4`.  Hence the exact
conditional marginal vector is

```math
\boxed{
\Pr\{e\in F_S\mid S\}
=\frac{1-\mathbf 1_{\{e\subset S\}}}{4}.}
\tag{R50C.1}
```

Now let `S` have any law `lambda` on maximal selectors and put

```math
w_e=\Pr_\lambda\{e\subset S\}.
```

Averaging `(R50C.1)` gives

```math
\boxed{
\Pr\{e\in F_S\}=\frac{1-w_e}{4}=p_e.}
\tag{R50C.2}
```

Thus the *exact* canonical fractional point of `(10.1240)` belongs to the
cut-switching polytope.  This is not an approximation or a generic
independent rounding.

## 2. All factors of two and four: the margin is orbit-averaged slack

Write an oriented state as `d=(sigma,x)` and use

```math
M_{d,ij}=a_{ij}\sigma x_ix_j,
\qquad E_A(d)=2\sum_eM_{d,e},
\qquad s_A(d)=q_n-E_A(d).
\tag{R50C.3}
```

Flipping an edge changes its contribution from `2M` to `-2M`, so for every
edge set `F`,

```math
E_{A^F}(d)=E_A(d)-4\sum_{e\in F}M_{d,e},
```

and consequently

```math
\boxed{
q_n-E_{A^F}(d)
=s_A(d)+4\sum_{e\in F}M_{d,e}.}
\tag{R50C.4}
```

For `F=delta(U)`, put `D_U=diag((-1)^{\mathbf1_U(i)})`.  Then

```math
A^{\delta(U)}=D_UAD_U,
```

so, writing `d^U=(sigma,D_Ux)`,

```math
q_n-E_{A^{\delta(U)}}(d)=q_n-E_A(d^U)=s_A(d^U).
\tag{R50C.5}
```

In the uniform-subset branch, the spins outside `S` are independently
randomized.  Every cross or outside-edge contribution has mean zero, while
the internal energy

```math
c_S(d[S])=2\sum_{e\subset S}M_{d,e}
```

is fixed.  Including the extra half-mass at `U=emptyset` proves the
**Verified identity**

```math
\boxed{
\begin{aligned}
\mathbb E\left[q_n-E_{A^{F_S}}(d)\mid S\right]
&=s_A(d)+4\sum_e
 \frac{1-\mathbf1_{\{e\subset S\}}}{4}M_{d,e}\\
&=\frac12\{s_A(d)+q_n-c_S(d[S])\}.
\end{aligned}}
\tag{R50C.6}
```

The middle expression also follows directly from
`E_A(d)=q_n-s_A(d)`; no Bellman inequality is being used.

For a maximal selector, recall

```math
c_S(d[S])=q_*+s_S(d[S])-j_S(d[S]).
```

Substitution into `(R50C.6)` gives

```math
\boxed{
\mathbb E\left[q_n-E_{A^{F_S}}(d)\mid S\right]
=\frac12\{q_n-q_*+s_A(d)-s_S(d[S])+j_S(d[S])\}.}
\tag{R50C.7}
```

Averaging `(R50C.7)` over `lambda` is exactly `(10.1240)`.  Therefore the
all-state Bellman pressure is not merely *representable* by a symmetry
mixture: selector by selector it is literally expected slack along that
mixture.

## 3. Explicit migrating-ground circulation

Every outcome in the preceding mixture is norm preserving:

```math
Q(A^{\delta(U)})=Q(D_UAD_U)=Q(A)=q_n.
\tag{R50C.8}
```

More concretely, fix any oriented parent ground `g=(sigma,x)` of `A`.  For
the outcome `U`, the state

```math
g_U=(sigma,D_Ux)
```

is a parent ground of `A^{delta(U)}`, since

```math
E_{D_UAD_U}(g_U)=E_A(g)=q_n.
\tag{R50C.9}
```

This is an exact witness circulation indexed by the switching group.  It
explains the strict Jensen gap

```math
\min_d\mathbb E_F\{q_n-E_{A^F}(d)\}>0,
\qquad
\mathbb E_F\min_d\{q_n-E_{A^F}(d)\}=0
\tag{R50C.10}
```

which `(10.1240)` can exhibit.

The circulation also respects the selector that created the pressure.  Every
sampled `U` lies in `S^c`, so

```math
(D_UAD_U)[S]=A[S],
\qquad (D_Ux)[S]=x[S].
\tag{R50C.11}
```

If `S` is bad, the migrated parent ground therefore retains exactly the same
nonground child restriction and the same child gap.  Mixing over the whole
maximal-selector family does not break this selector-by-selector witness
map.  In particular, the hypothesis `j_S>=4` increases the expected slack in
`(R50C.7)` but does not obstruct the circulation.

More globally, switching preserves `Q(A[R])` for every selector `R`, because
`(D_UAD_U)[R]` is switching-equivalent to `A[R]`.  It therefore preserves the
entire maximal-selector family and whether each member is good or bad.  A
hypothetical all-maximal-bad premise would hold on every outcome of this
circulation, not exclude its support.

## 4. Exact A8 check and scoped falsification

On the exact `A8,m=4` minimizer, use the equal law on the two bad maximal
selectors

```math
S_0=\{0,3,4,5\},\qquad S_1=\{1,2,6,7\}.
```

The resulting canonical coordinates are `p_e=1/8` inside one of the two
blocks and `p_e=1/4` across them.  Exhaustive exact evaluation recovers the
fractional minimum margin eight from Wave 49.  The new conditional
decomposition checks `32` switching-support outcomes (`16` for each
selector; the extra empty atom changes its weight but not the support).
Every outcome has
minimum integral margin zero, and the circulated ground retains child gap
at least four.  This is a finite illustration of the general proof, not the
proof's basis.

The conclusion must be scoped carefully.  Equations `(R50C.1)--(R50C.11)`
**falsify every attempt to round or anti-migrate the canonical point using
only its edge marginals, the full state-dependent identity, and the complete
selector law.**  The point has a universal, selector-compatible realization
entirely by norm-preserving symmetries, with explicit adaptive ground
witnesses.

This does **not** falsify tight decomposition: no exact minimizer with an
all-maximal-bad slice is produced.  Nor does it rule out a proof which uses
the all-bad premise to construct a genuinely noncanonical, non-switching
certificate or a new higher-order incompatibility absent from `(10.1240)`.
It does show that `(10.1240)` itself supplies no descent pressure after
quotienting by the switching symmetry, so the previously surviving
canonical global anti-migration implementation should be retired.
