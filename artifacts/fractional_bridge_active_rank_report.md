# Fractional bridge: near-cube-vertex reduction and its limit

This is a constructive-track proof report. It gives an exact rounding
implication and identifies the remaining obstruction; it does **not** claim
that the uniform lemma has been proved. The projected active-rank language
below is only a characterization of vertex fractionality, not a mechanism
that bounds it.

## 1. The automatically feasible fractional polytope

Fix child signings `A,B` of orders `m,n`, with caps `p,q`. Put

```math
T=(p^{2/3}+q^{2/3})^{3/2}.
```

For projective Boolean states `x,y`, write `h_A(x),h_B(y)` for their child
energies and set

```math
r_{x,y}=T-|h_A(x)+h_B(y)|.
```

Since `T>=p+q`, every radius is nonnegative. Define the central polytope

```math
\mathcal P_T(A,B)=\left\{C\in[-1,1]^{m\times n}:
 |\langle C,xy^{\mathsf T}\rangle|\le r_{x,y}
 \text{ for every }x,y\right\}.                    \tag{F1}
```

The zero bridge belongs to (F1). If `p,q>0`, then `T>p+q`, so zero is an
interior point and the polytope is full-dimensional. The exact block identity
shows that every `C in P_T` is a fractional bridge whose parent cap is at
most `T`.

This removes the feasibility obligation entirely at the fractional level;
the only question is how close a feasible point can be to a cube vertex.

## 2. Projected active rank is exactly vertex fractionality

Let `C` be a vertex of (F1), and let `F(C)` be the coordinate set on which
`|c_ij|<1`; write `f(C)=|F(C)|`. Project every state tensor `xy^T` whose
two-sided inequality in (F1) is active at `C` onto these fractional
coordinates, and let `R(C)` be the rank of the projected tensors.

**Proved.**

```math
f(C)=R(C).                                           \tag{F2}
```

Indeed, after the `mn-f(C)` integral coordinates are fixed by their active
box facets, a nonzero infinitesimal motion supported on `F(C)` would preserve
all active constraints unless the projected active state tensors span all of
`R^{F(C)}`. Thus `R(C)>=f(C)`; the reverse inequality is automatic because
the projected space has dimension `f(C)`.

For independent sign rounding `Z_ij` with `E Z_ij=c_ij`, define

```math
V(C)=\sum_{i,j}(1-c_{ij}^2).
```

Every integral coordinate contributes zero and every fractional coordinate
at most one, so (F2) gives

```math
V(C)\le f(C)=R(C).                                   \tag{F3}
```

This equality is **tautological as a research reduction**: vertexhood forces
the active normals to span every fractional direction. It does not make the
obligation `f(C)=o(mn)` easier. Using the unprojected rank would not fix this:
state inequalities can remain active at an integral cube vertex, so their
total rank need not vanish when `f=0`.

## 3. Uniform low-fractionality lemma sufficient for convergence

For each fixed state pair, the rounding error is a sum of independent,
mean-zero terms bounded by two and with total variance `V(C)`. Bernstein's
inequality and a union bound over at most `2^(m+n)` state pairs show that some
integral rounding has simultaneous error at most

```math
\Delta\le \sqrt{2V(C)L}+\frac43L,
\qquad L=(m+n)\log2+O(1).                            \tag{F4}
```

The constants in (F4) can be replaced by the precise constants in the main
fractional-rounding report; only the powers are used here.

Suppose `m,n` are comparable, `N=m+n`, and for at least one pair of exact
child minimizers at each pair of orders there is a vertex of (F1) satisfying,
uniformly,

```math
f(C)=O(N^{2-eta})                                    \tag{F5}
```

for some fixed `eta>0`. Equations (F3)--(F4) then give an integral parent of
cap

```math
T+O(N^{3/2-eta/2}+N).
```

Because `T=Theta(N^(3/2))` for comparable blocks, taking the `2/3` power gives

```math
M_N^{2/3}\le M_m^{2/3}+M_n^{2/3}
 +O(N^{1-eta/2}+N^{1/2}).                            \tag{F6}
```

The defect in (F6) is geometrically summable. Thus (F5), for balanced child
orders, is a precise direct route to convergence; it does not require a
structured family or a separate landing theorem.

## 4. Hidden-equivalence audit

An integral bridge feasible at cap `T` is itself a cube vertex of (F1), with
an empty fractional-coordinate space and hence `f=R=0`. Thus exact ideal
bridge feasibility implies (F5). The converse is
not required: (F5) permits `N^(2-eta)` genuinely fractional entries and only
produces an integral bridge after accepting the sublinear power defect in
(F6). It is therefore strictly weaker in its requested conclusion than exact
ideal bridge optimization.

Nor does (F5) merely rename fractional feasibility: (F1) is always feasible
at zero. It asks for a geometric property of at least one extreme point of a
known nonempty central polytope.

There remains a serious hidden-difficulty warning. Rank-one Boolean tensors
span the full `mn`-dimensional bridge space, so general polytope theory gives
only `f(C)<=mn`, which is useless. Equations (F2)--(F3) merely rename this
bound as projected active rank. Proving (F5) still needs an entropy, symmetry,
partial-coloring, or degeneracy theorem that forces a vertex near the cube
skeleton.
The earlier common-active-face theorem does not provide this: polynomial
moment control does not bound the rank of exponentially selected tight
states.

Accordingly, (F5) is a valid sufficient target but **not primary progress**
by itself. The only genuinely removed obligation is fractional feasibility,
which was trivial because the zero bridge works. Until a non-tautological
bound on `f(C)` is proved, this route is another precise statement of
near-integral bridge existence.

There is also no free bootstrap from a known integral bridge. Suppose `R` is
an integral bridge of parent cap `Q>T`. For every state its cross magnitude is
at most `Q-|h_A+h_B|`. Therefore `C=sR` is feasible in (F1) for

```math
s=\frac{T-(p+q)}{Q-(p+q)},                           \tag{F7}
```

because `(T-a)/(Q-a)` decreases with `a<=p+q`. Its variance is
`mn(1-s^2)`. Thus a power-saving variance obtained by this scaling already
requires `Q-T=o(N^(3/2))`, i.e. a power-saving approximate integral
composition theorem. Random rounding then generally halves the exponent.
Fractional scaling does not manufacture the first power saving; a useful
fractional argument must locate a near-cube feasible point by genuinely
different convex or algebraic structure.

## 5. Conference evidence and its exact limit

If `A,B` are complementary principal blocks of an integral conference parent
whose cap is at most `T`, its bridge belongs to (F1) and is an integral cube
vertex. Thus `f=0` (and the projected rank is zero) for those structured
pairs. This explains algebraically why internal conference composition is
easy. It makes no claim about the unprojected rank of other active state
constraints.

It does not transfer to true minimizers: replacing a structured child by a
lower-cap child changes every radius in (F1), and no conference identity
controls which state constraints become tight. Establishing that transfer is
again a landing theorem. Scaling an integral conference bridge also fails to
give (F5): if `C=sR` on `Theta(N^2)` entries with fixed `s<1`, then
`V=Theta(N^2)`.

## 6. A testable but not yet simpler formulation

Define

```math
\phi_T(A,B)=\min\{f(C):C\text{ is a vertex of }\mathcal P_T(A,B)\}. \tag{F8}
```

The uniform target is that at every balanced pair of orders some exact
minimizer representatives satisfy `phi_T(A,B)=O(N^(2-eta))`. A computational
test should optimize generic linear
objectives over (F1), record the number of fractional coordinates and the
projected rank of tight state tensors, and compare `phi/N^2` across increasing
orders. Merely adding more separator constraints does not test (F8): the LP
must expose an actual vertex and certify all omitted state inequalities by
exact separation.

**Current status: open, but noncircular and falsifiable.** To falsify the
existential uniform claim at a sequence of order pairs, one must show that
every choice of exact-minimizer representatives and every vertex of (F1) has
`f(C)>=c mn`. Conversely, any proof of (F5) removes the integer bridge and
structured-landing obligations at once.

## 7. A quantitative obstruction via restricted discrepancy

There is a proof-relevant way to falsify low fractionality without enumerating
LP vertices. For a set `E` of child-state pairs define

```math
\operatorname{disc}(E)=
 \min_{Z\in\{+1,-1\}^{m\times n}}
 \max_{(x,y)\in E}|x^{\mathsf T}Zy|,
\qquad
r(E)=\max_{(x,y)\in E}r_{x,y}.                       \tag{F9}
```

If `C in P_T` has variance `V`, the same simultaneous rounding argument,
restricted to `E`, proves

```math
\operatorname{disc}(E)
\le r(E)+\sqrt{2V L_E}+\frac43L_E,
\qquad L_E=\log(2|E|).                               \tag{F10}
```

Consequently, if along a balanced sequence `|E|<=2^(O(N))` and

```math
\operatorname{disc}(E)-r(E)\ge\gamma N^{3/2}        \tag{F11}
```

for a fixed `gamma>0`, then every feasible fractional bridge has
`V=Omega(N^2)`; in particular no vertex satisfies (F5). The additive
`O(N)` term in (F10) is negligible, and squaring the remaining inequality
gives the claim.

The natural choice of `E` is an orbit-compressed collection of aligned
near-extremal child pairs, because those pairs have the smallest radii. This
is materially different from adding unclassified separator states: (F11)
asks for a scalable lower bound on a restricted discrepancy norm. Conference
complement pairs cannot satisfy (F11), since their known integral bridge has
zero variance. For exact minimizers the sign and scale of (F11) are open.
