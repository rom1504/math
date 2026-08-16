# Adversarial audit of the projective Max-Cut boundary response draft

**Object audited:** `benchmark_maxcut_projective_response.md`, Propositions
MC.1--MC.2, Lemma MC.3, and Theorem MC.4.

**Verdict:** accept after two scope/rigor edits and one verifier repair.  The
pure nonnegative Max-Cut packing and its `Omega(2^w)`-bit conclusion at raw
additive error `epsilon*w` are correct.  That conclusion is deliberately
unnormalized: the constructed graphs have boundary sensitivity
`2^(w-2)` and total weight `w*2^(w-1)`.  It does not conflict with the
unit-sensitivity Lipschitz cover.

## 1. MC.1: exact state and pure-Max-Cut exposure

The conditional response

```math
h_G([\sigma])=\max_z \operatorname{Cut}_G(\sigma,z)
```

is projective because complementing all spins preserves every cut edge.
When private vertex sets are disjoint, conditioning on the common boundary
makes the two private optimizations independent, so

```math
\operatorname{MaxCut}(G\cup_B C)=
\max_{s\in X_w}\{h_G(s)+h_C(s)\}.
```

The pinning gadget is also correct.  Conditional on its anchor spin, a
direct edge for `tau_i=-1` earns one exactly when the desired inequality
holds.  A two-edge path for `tau_i=+1`, after maximizing its private middle
spin, earns two for equality and one for inequality.  Therefore

```math
p_\tau([\sigma])=c_\tau-
 \min\{d_H(\sigma,\tau),w-d_H(\sigma,\tau)\}.
```

Private copies add this profile.  More than
`max{osc(h_G),osc(h_G')}` copies make the target projective class the unique
maximizer for both profiles.  This exposes a coordinate attaining
`||h_G-h_G'||_infinity` and proves the isometry in MC.1.  All edges are
nonnegative and unit-weight; no field or signed factor is being smuggled in.

No correction is required in MC.1.  For clarity, its statement can say
`w>=1`, with the one-class case `w=1` understood as trivial.

## 2. MC.2: constant-weight spectral atom packing

Let `q=|S|`.  Greedy selection in the middle layer gives a constant-weight
code with relative Hamming distance at least `delta` and cardinality

```math
{2^{(1-H_2(\delta))q}\over q+1}
```

up to the stated integer rounding.  For incidence vectors `u,v`,

```math
\|P(u-v)\|_\infty
\ge q^{-1/2}\|P(u-v)\|_2
\ge q^{-1/2}s_{min}(P)\|u-v\|_2
\ge s_{min}(P)\sqrt\delta.
```

Thus MC.2 is correct.  Constant weight is not needed for the abstract lemma,
but it is essential in MC.4 because it cancels the padded constant profile.

For a completely finite statement, replace “up to an inessential rounding
change” by the following.  Put `d=ceil(delta*q)`.  Greedy selection gives

```math
|\mathcal F|\ge
{\binom q{\lfloor q/2\rfloor}\over
 \sum_{i=0}^{d-1}\binom qi},
\qquad d_H({\bf1}_U,{\bf1}_V)\ge d.
```

The displayed entropy bound then follows after replacing `delta` by
`delta+O(1/q)`.  This is a rigor edit, not a change to the asymptotic result.

## 3. Independent derivation of MC.13

Let `w=2m`, let `A` have size `2j`, and put

```math
a=j-1,\qquad b=m-j,\qquad a+b=m-1.
```

For independent signs write `S_n=sum_i X_i` and define the normalized full-
cube Fourier coefficient

```math
c_{2j}=\mathbb E\left[|S_{2m}|\prod_{i\in A}X_i\right].
```

With `D f(t)=(f(t+1)-f(t-1))/2`, conditioning on the signs in `A` gives
`c_(2j)=E[D^(2j)|T|]`, where `T` is a sum of `2b` signs.  On the even lattice,

```math
D^2|t|={\bf1}_{\{t=0\}}.
```

Moving the remaining even difference onto the law of `T`, or equivalently
extracting its generating-function coefficient, gives

```math
c_{2j}={1\over2^{2m-2}}
[u^{a+b}](1-u)^{2a}(1+u)^{2b}.                    \tag{A.1}
```

The central Chu--Vandermonde identity is

```math
[u^{a+b}](1-u)^{2a}(1+u)^{2b}
=(-1)^a{(2a)!(2b)!\over(a+b)!a!b!}.               \tag{A.2}
```

The projective sum is half the full-cube sum.  Since
`d_w=(2m-|S_(2m)|)/2` and every nontrivial character sums to zero, its
distance-matrix eigenvalue has magnitude

```math
2^{2m-2}|c_{2j}|
={(2j-2)!(2m-2j)!\over
  (m-1)!(j-1)!(m-j)!},
```

which is exactly the even line of MC.13.

For `w=2m-1`, condition a `2m`-spin sum on its last sign.  The remaining
sum `s` is odd and

```math
\tfrac12(|s+1|+|s-1|)=|s|.
```

The normalized even Fourier coefficients are therefore unchanged, while
the projective group has half as many elements.  Every nontrivial eigenvalue
is half its `2m`-spin counterpart, proving the odd line of MC.13.

This derivation verifies every sign-independent magnitude in MC.13 without
relying on the abbreviated Krawtchouk recurrence in the draft.  The recurrence
there is correct but should either cite (A.2) or include this coefficient
calculation; “Pascal's identity twice” is too compressed for the claimed
self-contained proof.

### Required scope correction

Lemma MC.3 must state `w>=2`.  At `w=1`, the projective space has one point,
`D_1=[0]`, so it is singular; there is no nontrivial even character.  All
uses in MC.4 have large `w`, and the finite verifier already starts at two,
so this correction has no downstream effect.

## 4. Least singular value

For even `w`, substitute `a=j-1,b=m-j` in MC.13.  The magnitude is

```math
L(a,b)={(2a)!(2b)!\over(a+b)!a!b!},
\qquad a+b=m-1.
```

Its consecutive ratio is

```math
{L(a+1,b-1)\over L(a,b)}={2a+1\over2b-1},
```

so the minimum occurs at the most balanced `a,b`.  If `m-1=2k` this gives
`binom(2k,k)`; if `m-1=2k+1` it gives `2 binom(2k,k)`.  Odd `w` divides these
values by two.  The trivial eigenvalue is the nonnegative row sum and has
absolute value at least every character sum, so it is not smaller.

The central-binomial estimate

```math
\binom{2k}{k}\ge {2^{2k}\over2k+1}
```

implies exactly

```math
s_{min}(D_w)\ge {2^{m-2}\over m}.
```

Thus MC.14--MC.15 are correct for `w>=2`.

Padding is legitimate.  The unpadded gadget has exactly `c_tau` unit edges;
adding `2w-c_tau` private disjoint unit edges changes its profile to
`2w-d_w`.  On nontrivial characters the padded matrix is `-D_w`; on the
constant character its eigenvalue is `2w*q-row_sum(D_w)>0`.  Hence no channel
is lost.

## 5. MC.4 and the information conclusion

Let `q=2^(w-1)` and `k=q/2`.  Every selected component contains exactly `k`
padded gadgets, so

```math
h_U=2wk\,{f1}-D_w{f1}_U.
```

The constant term cancels for the middle-layer code.  MC.2 applied directly
to `D_w` and MC.15 give

```math
\|h_U-h_V\|_\infty
\ge {\sqrt\delta\,2^{\lceil w/2\rceil-2}
       \over\lceil w/2\rceil}.
```

For fixed `epsilon,delta>0`, this eventually exceeds `2epsilon*w`.  The code
contains

```math
2^{(1-H_2(\delta)-o(1))2^{w-1}}
```

components.  MC.1 turns profile sup separation into future-response
separation.  Two components sharing one deterministic summary decoded to
error at most `epsilon*w` would be at distance at most `2epsilon*w`, a
contradiction.  The required message length is therefore
`Omega(2^(w-1))=Omega(2^w)` bits.  This proof is correct inside ordinary
unit-edge, nonnegative Max-Cut; the graph size is unrestricted.

## 6. Exact normalization reconciliation

The lower bound above is not a unit-sensitivity theorem.  Each padded gadget
has exactly `2w` unit edges and exactly one edge incident to each boundary
vertex.  A component has `k=q/2=2^(w-2)` gadgets.  Therefore

```math
\deg_B(b)=k=2^{w-2}\quad(b\in B),
\qquad
\text{total edge weight}=2wk=wq=w2^{w-1}.          \tag{A.3}
```

The raw error `epsilon*w` is only an `epsilon/q` fraction of total weight.
It is exponentially finer than the natural score scale of this family.

To impose weighted boundary degree one, multiply every edge weight by
`1/k=2/q`.  The total weight then becomes `2w`, but the proved separation
becomes only

```math
{2\over q}s_{min}(D_w)\sqrt\delta
\ge {\sqrt\delta\,2^{\lceil w/2\rceil-w}
       \over\lceil w/2\rceil},                    \tag{A.4}
```

which tends to zero and cannot support error `epsilon*w`.

Conversely, unit boundary degree makes every response one-Lipschitz in
boundary Hamming distance.  The coarse Lipschitz theorem in
`benchmark_maxcut_boundary_response.md` gives, modulo a constant, an
`epsilon*w` cover whose bit length is at most

```math
2^{(1-H_2(\epsilon/4)+o(1))w}=o(2^w)
```

for fixed `0<epsilon<1/4`.  In the rescaled projective family the total score
lies in an interval of length `O(w)`, so quantizing the additive offset costs
only `O(log(1/epsilon))` more bits.  Hence the two results are consistent:

- MC.4 is a valid unrestricted-size, unnormalized lower bound at a very fine
  raw error scale;
- the Lipschitz theorem is a normalized coarse-error upper bound and rules
  out extending MC.4 unchanged to unit terminal sensitivity.

The existing final caveat in MC.4 mentions polynomial size and bounded total
weight.  It should add the explicit quantities (A.3) and the rescaling
calculation (A.4), so the scale distinction is visible in the theorem itself.

## 7. MC.5 lookup universality

The later lookup theorem MC.5 is correct.  Its algebra can be checked without
an expressivity assumption.

For fixed boundary `s` and anchor `z`, put `x_i=s_i z` and
`y_a=(1+t_a z)/2`.  Then

```math
\sum_i a_i x_i-(w-1)=1-2d_H(a,x).
```

Exactly the oriented word `a=x` has positive coefficient, equal to one;
every other coefficient is at most minus one.  Since `lambda_a>=0` and the
`t_a` independently choose `y_a in {0,1}`, maximizing them returns
`lambda_x=F([s])`.  Symmetry of the lift removes the arbitrary orientation
of `z`.

The alleged cubic term is indeed pairwise:

```math
{1+t_a z\over2}\,a_i s_i z
={a_i\over2}(s_i z+t_a s_i),
```

and expanding the threshold gives MC.26 exactly.

For a signed pair term `Juv`, the two replacements have the advertised
constants.  A direct cut edge of weight `-2J` for `J<=0` gives
`-J+Juv`.  For `J>=0`, maximizing a fresh two-edge path of weights `2J`
gives `4J` when `u=v` and `2J` when `u=-v`, namely `3J+Juv`.

The offset audit is as follows for one antipodal pair `a,-a`, both of weight
`lambda`.

- Across the pair, the `s_i z` occurrences add replacement constant
  `2w lambda`.
- The `t_a s_i` occurrences add another `2w lambda`.
- The two negative `t_a z` terms add `(w-1)lambda`.
- The explicit constant in MC.26 is `-(w-1)lambda` for the pair, so removing
  it adds `(w-1)lambda`.

The net offset is therefore `(6w-2)lambda`, proving MC.22.  Summing over
projective classes gives at most `(6w-2)qW`; the isolated padding edge in
MC.23 has nonnegative weight and, being private and disconnected, adds
exactly that weight to every boundary fibre.  Hence the translated sup cube
is realized with a genuinely common offset.

MC.1 then makes its contextual metric exactly the cube sup metric.  Standard
grids give logarithmic packing and covering complexity

```math
Theta(q log(W/epsilon)),\qquad q=2^{w-1},
```

with the stated universal radius changes.  If `W=Theta(w)` and
`epsilon_0 w<=W/6`, this is `Theta(2^w)` bits.  Thus MC.24 and its corollary
are correct.

Two scope notes should accompany MC.5.

1. Literal per-occurrence replacement may create parallel direct edges;
   merging direct parallel edges after the replacement produces an ordinary
   weighted graph without changing the score.  Positive-coupling paths have
   fresh middle vertices.
2. The construction is still far from unit terminal sensitivity.  Each
   oriented word contributes boundary incident weight `2lambda_a` at every
   boundary coordinate, so

   ```math
   \deg_B(i)=4\sum_{t\in X_w}F(t),
   ```

   which can be `Theta(2^w W)`.  The common padding is private, but the lookup
   graph and its offset have exponential total scale.  MC.5 proves full
   response-shape universality for unrestricted weighted Max-Cut; it does not
   contradict the smaller coarse cover under unit sensitivity.

MC.5 also supersedes the weaker disclaimer in
`benchmark_maxcut_boundary_response.md`: the abstract response-cube packing
there is in fact realizable in pure Max-Cut after imposing global-flip
symmetry, provided unrestricted nonnegative weights and exponential-size
gadgets are allowed.

## 8. Finite-verifier correction and archive collision

The existing verifier correctly checks eigenvalue **magnitudes** through its
default widths.  Its comment that `distance_rows[0]` is the identity row and
that every character is one there is false: `projective_words` currently
lists `(1,-1,...,-1)` first.  Because the code takes `abs(lam)`, the tests
still pass.  Repair either by listing `(1,1,...,1)` first or by dividing the
computed character sum by `character(words[0],subset)`.

That repair has now been applied.  Exact enumeration through `w=11` verifies
MC.13--MC.15, and the lookup verifier through `w=4` verifies both the selected
table entry and the offset `(6w-2)sum F`.

There are three remaining editorial corrections in the audited draft.

- MC.3 must exclude `w=1`, as above.
- The Krawtchouk coefficient and the MC.4 family size both carry tag
  `MC.18`; renumber Section 4 onward.
- In MC.4, “boundary-normalized additive error” is misleading because the
  boundary load is `2^(w-2)`.  Use “width-scaled raw additive error” or simply
  “additive error `epsilon w`.”

After the blind derivation, comparison with the archive finds that the
conditional-separator response, pinning/sup isometry, and ambient response-
cube entropy collide with `phase2_feature_growth.md` and
`phase3_contextual_response_law.md`.  The two-reference/path Max-Cut exposure
gadget and the projective-distance spectral packing are model-specific.  No
archive passage inspected supplies the pure-Max-Cut `Omega(2^w)` packing in
MC.4.  Its genuine addition is precisely that restricted-language result,
with the normalization boundary in Section 6 above.
