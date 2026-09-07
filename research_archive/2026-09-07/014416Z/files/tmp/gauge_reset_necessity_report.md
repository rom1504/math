# Gauge--reset necessity: a counterexample and a sharp selector converse

**Status.** Proof-level independent report. No repository terminology is
assumed beyond max-plus continuation and Hilbert's projective metric.

## Executive conclusion

There are two different converse questions, and they have different answers.

1. If a local defect is allowed to be an **arbitrary disturbance at every
   use** of a continuation, then on a finite coordinate-selector system there
   is a sharp gauge--reset converse: after subtracting an endpoint
   coboundary, depth-uniform `O(epsilon)` stability is equivalent to a
   bounded gap between rank-one (projectively zero) selector products.

2. If one compares two **fixed coherent families of exact continuations**,
   that converse is false. Exact semigroup relations can absorb the defect
   without either an entrywise kernel gauge or a small projective image. An
   all-finite `2 by 2` max-plus example consists of two nearby idempotent
   interval clamps. Their error is `O(delta)` at every depth, their kernel
   difference has nonzero rectangle defect, and every positive power retains
   projective diameter bounded away from zero.

Thus the existing gauge/reset mechanisms are complete for a natural robust
linearized model, but not for coherent approximate actions. A third mechanism
must be admitted there:

```math
\boxed{\text{bounded normal forms / exact semigroup absorption}.}
```

At the infinitesimal level the right holonomy is also *twisted* by the active
selector, rather than being an ordinary sum of labels. The exact one-cycle
criterion is proved below.

## 1. Why the statement needs a scale and a perturbation model

For an all-finite max-plus matrix `S`, the map

```math
(F_Su)_b=\max_a(u_a+S_{ab})
```

has finite projective image diameter. Consequently the image diameter of a
product is bounded by that of its last factor. Mere boundedness of product
diameters therefore has no rigidity content for a finite family of
all-finite matrices.

The meaningful assertion is a bound which vanishes with the one-step defect,
for example

```math
\sup_{w,x}d_H(F_wx,G_wx)\le C\epsilon,
\qquad
\sup_{\ell,x}d_H(F_\ell x,G_\ell x)\le\epsilon,       \tag{1.1}
```

with `C` independent of word length. Even (1.1) has two inequivalent
readings:

- **robust:** an adversary may choose a fresh error of size `epsilon` after
  every transition;
- **coherent:** `F_ell` and `G_ell` are fixed maps, so their errors satisfy
  every algebraic relation possessed by those maps.

The robust theorem below gives a converse. The coherent clamp example then
shows why it cannot be promoted to fixed exact actions.

## 2. A third mechanism: finite-semigroup absorption

### Theorem 2.1 (bounded normal forms give depth-uniform stability)

Let `A` be a finite alphabet generating a finite semigroup `S`. Let
`ell_S(s)` be the shortest length of a word over `A` representing `s`, and
put

```math
L=\max_{s\in S}\ell_S(s).
```

Suppose `F` and `G` are two exact actions of `S` by nonexpansive maps on a
metric space `X`. Equivalently, both generator families satisfy all the
relations of the same presented semigroup. If

```math
\sup_{x\in X}d(F_a x,G_a x)\le\epsilon
\quad(a\in A),                                      \tag{2.1}
```

then, for every word `w`,

```math
\boxed{\sup_x d(F_wx,G_wx)\le L\epsilon.}           \tag{2.2}
```

#### Proof

Let `s` be the element represented by `w`, and choose a word `v` of length at
most `L` representing `s`. Exact factorization through `S` gives
`F_w=F_v` and `G_w=G_v`. Replace the factors of `F_v` by those of `G_v` one
at a time. Nonexpansiveness of every suffix and (2.1) show that each
replacement costs at most `epsilon`; hence the total cost is at most
`|v|epsilon<=L epsilon`. `square`

The theorem is not a reset theorem: the maps can have large or unbounded
images. It is not additive cohomology either: a defining relation can be
nonlinear, such as idempotence. Its input is algebraic compatibility of the
two exact actions.

### Example 2.2 (all-finite idempotent max-plus counterexample)

Use the projective coordinate `z=u_2-u_1` on
`R^2/R1`. In this coordinate

```math
d_H(z,z')={1\over2}|z-z'|.
```

For `0<delta<1`, consider the all-finite matrices

```math
S_0=\begin{pmatrix}0&0\\-1&0\end{pmatrix},
\qquad
S_\delta=\begin{pmatrix}0&\delta\\-1&0\end{pmatrix}. \tag{2.3}
```

Direct calculation gives

```math
F_{S_0}(z)=\max(0,z)-\max(0,z-1)
          =\operatorname {clip}(z,0,1)=:P_0(z),
```

and

```math
F_{S_\delta}(z)=\max(\delta,z)-\max(0,z-1)
          =\operatorname {clip}(z,\delta,1)=:P_\delta(z). \tag{2.4}
```

Both are nonexpansive idempotents. Therefore, for every `t>=1`,

```math
P_0^t=P_0,
\qquad P_\delta^t=P_\delta,
\qquad
\sup_z d_H(P_0^tz,P_\delta^tz)={\delta\over2}.       \tag{2.5}
```

This is Theorem 2.1 for the two-element monoid presented by `p^2=p`.

It escapes both mechanisms currently recorded in the repository:

1. With `E=S_delta-S_0`,

   ```math
   E=\begin{pmatrix}0&\delta\\0&0\end{pmatrix},
   \qquad
   E_{11}+E_{22}-E_{12}-E_{21}=-\delta.             \tag{2.6}
   ```

   Hence `E` is not a row-plus-column potential and fails the exact
   rectangle/gauge test.

2. The projective images are the intervals `[0,1]` and `[delta,1]`, of
   Hilbert diameters `1/2` and `(1-delta)/2`. Every positive power is the
   same idempotent, so there is no word whose image diameter is `O(delta)`.

Thus entrywise gauge plus small **full-image** reset is not complete for
coherent fixed perturbations on the full projective domain.

There is an important semantic caveat. If “zero-holonomy recurrent
dynamics” is expanded to mean *any* exact finite semigroup relation on the
paired orbit space, then the example is absorbed by that phrase: after one
step both orbits are stationary. But that enlarged phrase is strictly more
general than the kernel-potential and additive-label criteria proved so far.
The finite-semigroup action, not an additive edge label, is the extra state.

### Contrast 2.3 (equally small local error can drift)

The nearby all-finite matrix

```math
\widehat S_\delta=
\begin{pmatrix}0&0\\-1+\delta&\delta\end{pmatrix}  \tag{2.7}
```

induces

```math
\widehat P_\delta(z)=\operatorname {clip}(z+\delta,0,1).
```

It is also within `delta/2` of `P_0` in uniform Hilbert distance, but at
`z=0`,

```math
P_0^t(0)=0,
\qquad
\widehat P_\delta^t(0)=\min(t\delta,1).             \tag{2.8}
```

So the error reaches order one at depth `Theta(1/delta)`. Here
`\widehat S_delta-S_0` is a row potential and has zero rectangle defect, but
successive source potentials are not matched by terminal potentials; the
adjacent-interface circulation is nonzero. This pair makes the logical
point sharply: neither one-step distance nor a within-factor rectangle test
distinguishes stability. The exact relations of the transition semigroup do.

## 3. A sharp robust gauge--reset converse for selector dynamics

The following is a complete theorem in the local linear model of a
max-plus selector cell.

Let

```math
V=R^r/R1,
\qquad \|[v]\|_H={1\over2}\operatorname {osc}(v).
```

For a function `sigma:[r]->[r]`, let

```math
(P_\sigma v)_j=v_{\sigma(j)}.                    \tag{3.1}
```

Then `P_sigma` is nonexpansive on `V`, and it is zero on `V` exactly when
`sigma` is constant. Coordinate-selector derivatives of a max-plus map have
this form on every tie-free cell.

Consider a declared factorial language of selector words and the disturbed
recursion

```math
e_t=P_t e_{t-1}+\eta_t,
\qquad e_0=0,
\qquad \|\eta_t\|_H\le\epsilon.                  \tag{3.2}
```

A contiguous word is a **tangent reset** if its selector product is zero on
`V` (equivalently, its composite selector is constant).

### Theorem 3.1 (syndetic tangent resets are necessary and sufficient)

The following are equivalent up to the explicit constants shown.

1. There is an `L` such that every allowed word of length `L` contains a
   tangent-reset factor.
2. Every allowed disturbed trajectory satisfies

   ```math
   \|e_T\|_H\le L\epsilon                              \tag{3.3}
   ```

   for all `T>=L` (with the obvious `T epsilon` bound before that).

Conversely, if there is a reset-free allowed word of length `T`, then there
are disturbances in (3.2) for which

```math
\boxed{
\|e_T\|_H\ge
\left\lfloor{T\over r(r-1)}\right\rfloor\epsilon.}  \tag{3.4}
```

In particular, a bound `\|e_T\|_H<=C epsilon` against all disturbances
forces reset-free words to have length less than `(C+1)r(r-1)`.

#### Proof

Unrolling (3.2) gives

```math
e_T=\sum_{s=1}^T P_TP_{T-1}\cdots P_{s+1}\eta_s. \tag{3.5}
```

If the final length-`L` window contains a zero product, every error inserted
before that factor is killed. At most `L` errors remain, and
nonexpansiveness gives (3.3).

Now let a length-`T` word have no zero-product factor. For every `s`, the
suffix product

```math
A_s=P_T\cdots P_{s+1}
```

(including the empty identity at `s=T`) is a nonconstant selector. Thus
there are output coordinates `(j_s,k_s)` which it sends to two distinct
input coordinates. Among the at most `r(r-1)` ordered pairs, one pair
`(j,k)` occurs for at least `floor(T/[r(r-1)])` values of `s`. For each such
`s`, choose `eta_s` with values `+epsilon` and `-epsilon` on the two selected
input coordinates, and zero elsewhere; put all other disturbances equal to
zero. Then `\|eta_s\|_H=epsilon`, and every chosen summand contributes
`2epsilon` with the same sign to `e_T(j)-e_T(k)`. Hence the Hilbert norm of
the sum is at least the number of chosen summands times `epsilon`, proving
(3.4). `square`

### Corollary 3.2 (gauge plus reset is complete in this robust model)

Suppose the actual error recursion has the form

```math
e_t=P_t e_{t-1}+h_t-P_t h_{t-1}+\eta_t,          \tag{3.6}
```

where `h_t` is an endpoint gauge and `eta_t` is an otherwise arbitrary
residual with `\|eta_t\|_H<=epsilon`. Setting
`\bar e_t=e_t-h_t` reduces (3.6) exactly to (3.2). Therefore bounded endpoint
gauges plus syndetic tangent resets give a depth-uniform bound; and in the
absence of syndetic tangent resets, arbitrary residuals force (3.4).

This is a genuine converse, but only after the perturbation quantifiers are
made adversarial. The idempotent clamp does not contradict it: its repeated
defects are a highly constrained coherent sequence, not arbitrary `eta_t`.

The theorem also reveals that a **full projective-image reset is stronger
than necessary**. What is needed for error robustness is a reset of the
transported error directions along every reachable selector path. A
nonlinear max-plus word may have several far-apart selector images while its
derivative is rank one on each reachable branch. Thus the next nonlinear
converse should use a selector-product automaton, not only the global image
diameter of the word.

## 4. Twisted holonomy on one recurrent selector cell

Ordinary additive labels are also insufficient once a continuation
transports error directions. The elementary exact replacement is the
following.

### Theorem 4.1 (affine-selector cycle criterion)

Let

```math
A(e)=P_\sigma e+b
```

act on `V=R^r/R1`. The projective iterates `A^k(e)` are bounded for every
`e` if and only if all directed cycles `C` of the functional graph of
`sigma` have the same mean

```math
\beta_C={1\over |C|}\sum_{j\in C}b_j.            \tag{4.1}
```

Equivalently, there are `p in R^r` and `beta in R` such that

```math
\boxed{b=p-P_\sigma p+\beta\boldsymbol1.}         \tag{4.2}
```

In that case

```math
A^k(e)=p+P_\sigma^k(e-p)+k\beta\boldsymbol1,     \tag{4.3}
```

so projective boundedness is immediate. If two cycles have different means,
the projective distance grows linearly at the difference of their means.

#### Proof

If (4.2) holds, (4.3) follows by induction and constants disappear in `V`.
Conversely, summing (4.2) around a cycle would force (4.1), so it remains to
construct `p`. If all cycle means equal `beta`, define `p` on each cycle by
the recursion

```math
p_j-p_{\sigma(j)}=b_j-\beta;
```

the cycle sum is zero, so this is consistent. Define `p` recursively on each
tree feeding a cycle. This proves (4.2). If cycles `C,C'` have distinct
means, coordinates in their basins accumulate respectively
`k beta_C+O(1)` and `k beta_C'+O(1)`, giving linear projective separation.
`square`

This is the minimal finite theorem behind a “zero-holonomy recurrent piece.”
It has two parts:

- transient trees are directionally reset into the recurrent cycles;
- on the recurrent permutation, the translation is a *twisted coboundary*.

The twist `P_sigma` matters. Ordinary sums of untransported local error
vectors can report false drift or miss real drift.

## 5. What is and is not proved about the desired converse

### Proved

1. The narrow claim

   ```math
   \text{coherent depth-uniform stability}
   \Longrightarrow
   \text{entrywise gauge or small full-image reset}
   ```

   is false even for one all-finite `2 by 2` max-plus generator.

2. For arbitrary disturbances on a fixed selector itinerary, gauge plus
   syndetic **tangent** reset is necessary and sufficient, quantitatively.

3. On one recurrent affine selector cell, boundedness is equivalent to a
   transported/twisted coboundary on its cycles.

4. Exact factorization through a finite semigroup supplies a third
   depth-uniform mechanism, with constant equal to the semigroup's bounded
   normal-form diameter.

### Not proved

No global converse is claimed for arbitrary all-finite max-plus maps with
ties and switching selector cells. A perturbation can change the active
cell, so a derivative argument must track paired cells, not one nominal
itinerary. Nor is every stable coherent pair claimed to factor through a
finite semigroup; near conjugacy and invariant graphs are other plausible
coherent mechanisms.

## 6. Recommended next theorem

The correct next target is narrower and testable:

> **Finite selector-skew-product decomposition.** Build the reachable paired
> selector automaton of two piecewise-affine max-plus actions. Prove that a
> uniform `O(epsilon)` bound follows from, and under adversarial residuals is
> equivalent to, the following finite certificate on every pumpable strongly
> connected component: (i) the affine translation satisfies the twisted
> cycle-mean criterion on the surviving selector coordinates; or (ii) every
> sufficiently long path contains a selector product that is rank one on the
> residual-error quotient. Treat exact finite-semigroup relations as zero
> residual, not as a metric reset.

This target is stronger than the existing additive graph-holonomy theorem,
weaker than a false global gauge/full-reset converse, and has a finite
falsifier: a pumpable paired selector cycle with unequal recurrent cycle
means.

For theory terminology, the evidence supports a three-way distinction:

1. **endpoint/twisted gauge** (cocycle cancellation);
2. **relative or tangent reset** (old error directions die, even if the full
   state image remains large);
3. **algebraic absorption** (bounded normal forms make coherent defects occur
   only boundedly many times).

Only the first two survive fully adversarial per-step error. The third is
real and useful precisely because actual continuation kernels are fixed
coherent objects rather than an adversarial error oracle.
