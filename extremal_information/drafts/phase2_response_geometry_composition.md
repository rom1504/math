# Response geometry under composition: exact products, cancellation, and tropical rank

**Status.** Counterexample/axiomatizer report.  The Hilbert-space statements
and the robust tropical lemma are proved below.  The exact code-rank identity
is imported from Sheshadri and was independently reconstructed earlier in the
main ledger.  This report does not edit the theory surface files.

## 1. Director verdict

The response-separation polytope

```math
\Gamma(R)=\left\{\gamma\ge 0:
 \sum_{e:a_e\ne b_e}\gamma_e\le \|R_a-R_b\|^2
 \text{ for all }a,b\right\}
```

has an exact algebra for **orthogonally retained** channels.  It has no algebra
from the child polytopes alone when channels are added in the same response
space: their relative directions, which `Gamma` forgets, decide whether the
child displacements reinforce or cancel.  For universal max-plus boundary
gluing the failure is even stronger.  Every child entry has positive response
width before gluing, while after two or more factors the separation polytope
of the factor-to-parent map is exactly `{0}`.

This does not invalidate the posterior-width information theorem.  It gives
the correct information price for the response map to which it is applied.
It does show that `Gamma` is an information-converse certificate, not by
itself a reusable compositional state.  Under general composition it must be
recomputed from the parent response, or supplemented by joint displacement
geometry.

A stronger feature-algebra result survives in the tropical setting.  A
four-cell fooling-set gap gives a robust lower bound on approximate min-plus
factor rank.  Applied to the conditional distance table of a linear code, it
proves that every uniform approximation with error strictly below `1/2` has
at least `2^s` tropical channels, while the unperturbed table attains that
count.  This is a stable, structured
feature-growth theorem, rather than a statement about the universal matrix
cube.

## 2. Exact orthogonal product rule

The definition of `Gamma` makes sense on any finite binary parameter set; we
write the parameters as sign cubes only for consistency.  Let

```math
R:\{-1,+1\}^N\to\mathcal H,
\qquad
S:\{-1,+1\}^M\to\mathcal K
```

be maps into real Hilbert spaces.  For positive scalars `alpha,beta`, define

```math
T_{a,b}=\alpha R_a\oplus\beta S_b
\quad\text{in }\mathcal H\oplus\mathcal K.          \tag{RG.1}
```

For a set `C` of nonnegative vectors, write `c C={c gamma:gamma in C}`.

### Theorem RG.1 (orthogonal tensorization)

One has the exact identity

```math
\boxed{
\Gamma(T)=\alpha^2\Gamma(R)\times\beta^2\Gamma(S).}
                                                               \tag{RG.2}
```

Consequently,

```math
\boxed{
\kappa(T)=
\min\{\alpha^2\kappa(R),\beta^2\kappa(S)\}.}       \tag{RG.3}
```

The evident convention handles a zero scale: all weights on a response that
has been multiplied by zero must vanish.

#### Proof

Orthogonality gives, for every two parent parameters,

```math
\|T_{a,b}-T_{a',b'}\|^2
=\alpha^2\|R_a-R_{a'}\|^2
 +\beta^2\|S_b-S_{b'}\|^2.                          \tag{RG.4}
```

Thus a pair of scaled child certificates is a parent certificate.  Conversely,
if `(gamma,eta)` is a parent certificate, set `b=b'` in its defining
inequality.  This proves `gamma/alpha^2 in Gamma(R)`.  Setting `a=a'` proves
the analogous assertion for `eta`, hence (RG.2).

The quotient defining `kappa(T)` is bounded below by the smaller of the two
scaled child moduli, because its numerator and denominator are sums of the
two child quantities.  Varying only a minimizing child parameter attains
that smaller value.  This proves (RG.3). `square`

### Corollary RG.1 (positivity is not a scale statement)

Take `t` identical response maps and put the uniform-component `L^2` norm on
their orthogonal product:

```math
T^{(t)}_{a_1,\ldots,a_t}
={1\over\sqrt t}(R_{a_1}\oplus\cdots\oplus R_{a_t}).
```

Then

```math
\kappa(T^{(t)})={\kappa(R)\over t}.                 \tag{RG.5}
```

Every finite product is injective when `R` is, but its scalar modulus tends
to zero.  The posterior-width theorem correctly depends on the combination
`kappa times number of latent coordinates` and on the declared distortion
normalization; the bare assertion `kappa>0` carries no uniform content.

This exact product law is Hilbert product geometry.  Its value here is to
locate precisely the boundary of the response-separation object, not to claim
a new tensorization principle.

## 3. The smallest same-space cancellation counterexample

Orthogonality in Theorem RG.1 is essential.  Let the latent variables
`a,b` be signs and work in `R^2`.  Define

```math
R_a=a e_1,
\qquad S_b^\parallel=b e_1,
\qquad S_b^\perp=b e_2.                             \tag{RG.6}
```

All three one-bit response maps have the identical separation polytope

```math
\Gamma(R)=\Gamma(S^\parallel)=\Gamma(S^\perp)=[0,4].
                                                               \tag{RG.7}
```

Under same-space addition, however,

```math
T^\parallel_{a,b}=R_a+S_b^\parallel,
\qquad
T^\perp_{a,b}=R_a+S_b^\perp                       \tag{RG.8}
```

satisfy

```math
\boxed{
\Gamma(T^\parallel)=\{(0,0)\},
\qquad
\Gamma(T^\perp)=[0,4]^2.}                         \tag{RG.9}
```

Indeed, `T^parallel_(1,-1)=T^parallel_(-1,1)=0`.  The two parameters differ
in both bits, so nonnegativity and the separation constraint force both
coordinate weights to be zero.  In the orthogonal case the squared distance
is exactly four times Hamming distance, giving the displayed box.

Thus there is no rule that determines the parent separation polytope from
the two child separation polytopes under general addition.  The missing datum
is not another marginal modulus: it is the relative orientation of child
displacements.  Expanding a parent distance displays it explicitly:

```math
\|\Delta R+\Delta S\|^2
=\|\Delta R\|^2+\|\Delta S\|^2
 +2\langle\Delta R,\Delta S\rangle.                \tag{RG.10}
```

For reference, a sufficient joint-cancellation certificate is immediate.
If `gamma in Gamma(R)`, `eta in Gamma(S)`, and nonnegative vectors `c,d`
satisfy, for every four child parameters,

```math
2\langle R_a-R_{a'},S_b-S_{b'}\rangle
\ge
-\sum_{i:a_i\ne a'_i}c_i
-\sum_{j:b_j\ne b'_j}d_j,                          \tag{RG.11}
```

then the nonnegative coordinatewise differences
`(gamma-c,eta-d)` belong to `Gamma(R+S)`.  This is just (RG.10), but it states
the exact extra obligation: composition needs a joint cross-Gram lower bound.
Paying the two channels separately cannot see it.

The deterministic synchronization theorem addresses an analogous issue for
species profiles by imposing uniform no-crossing and scalar linkage.  It does
not automatically supply (RG.11): profile comparability and Hilbert
displacement orientation are different data.  Any claimed bridge between the
two must prove the cross-Gram bound rather than invoke the word
“synchronization.”

## 4. Complete collapse under universal max-plus gluing

Let `B` have cardinality `Q>=2`.  For `t>=1`, take binary kernels

```math
A^{(s)}\in\{0,1\}^{B\times B},\qquad 1\le s\le t,
```

and form their max-plus product

```math
F_t(\mathbf A)(i_0,i_t)
=\max_{i_1,\ldots,i_{t-1}\in B}
 \sum_{s=1}^t A^{(s)}(i_{s-1},i_s).                \tag{RG.12}
```

Regard the entries of all factors as `tQ^2` latent bits and flatten the
output kernel into a Hilbert space.  The normalization of the Hilbert norm is
irrelevant to the next result.

For one factor, with normalized Frobenius norm

```math
\|K\|_2^2={1\over Q^2}\sum_{a,b}K(a,b)^2,
```

the identity response `A -> A` has

```math
\Gamma(F_1)=[0,Q^{-2}]^{Q^2},
\qquad \kappa(F_1)=Q^{-2}.                          \tag{RG.13}
```

### Theorem RG.2 (max-plus posterior-width collapse)

For every `Q>=2` and every `t>=2`,

```math
\boxed{\Gamma(F_t)=\{0\},\qquad\kappa(F_t)=0.}     \tag{RG.14}
```

#### Proof

Equation (RG.13) follows because the squared distance between two binary
matrices is their Hamming distance divided by `Q^2`.  Single-coordinate
pairs impose the coordinatewise upper bounds, and those bounds imply every
remaining separation constraint.

Now fix any latent coordinate `(s,u,v)`.  Set every factor other than `s`
to the zero matrix.  In factor `s`, compare the all-one matrix `J` with
`J-E_(u,v)`.  These two input tuples differ in exactly the selected latent
bit.  Their max-plus products are nevertheless both the all-one matrix.  For
every fixed pair of external endpoints, a path can avoid the exceptional
transition `(u,v)`: at the first layer choose another outgoing state, at the
last layer choose another incoming state, and at an internal layer choose any
other transition.  Such a choice exists because `Q>=2`; the zero factors
place no restriction on the path.

The separation inequality for this colliding pair forces the weight of the
selected latent coordinate to be zero.  Since the coordinate was arbitrary,
every member of `Gamma(F_t)` is zero.  The formula for `kappa` follows.
`square`

This is a scalable counterexample, already present at two `2 by 2` factors.
It also answers the finite-modulus test decisively: a factor can have positive
`kappa`, yet one exact max-plus elimination can make the parent modulus zero.

There is no contradiction with the boundary-kernel composition theorem.
Two factor tuples with the same product kernel are equivalent for every
future serial endpoint context, so the eliminated bits are genuinely
irrelevant to that experiment.  The result instead says that child
posterior-width certificates do not propagate through tropical elimination.
The **kernel** is the reusable state; `Gamma` is a lower-bound certificate for
whichever kernel-valued response remains after elimination.

## 5. A robust tropical fooling-set theorem

The preceding collapse asks whether a structured class can nevertheless
force a large reusable kernel algebra.  Tropical factor rank measures this
different question.  For a finite real matrix `M`, let
`rank_(min,+)(M)` be the least `r` such that

```math
M(x,y)=\min_{1\le k\le r}\{u_k(x)+v_k(y)\}.         \tag{TR.1}
```

### Theorem TR.1 (robust tropical fooling set)

Suppose `M` has distinguished cells `(x_i,y_i)`, `1<=i<=r`, and `G>0`
such that for every `i!=j`,

```math
M(x_i,y_j)+M(x_j,y_i)
-M(x_i,y_i)-M(x_j,y_j)\ge G.                       \tag{TR.2}
```

If

```math
\|\widetilde M-M\|_\infty<G/4,                    \tag{TR.3}
```

then `rank_(min,+)(Mtilde)>=r`.

#### Proof and constant audit

Every rank-one term in a min-plus factorization majorizes the represented
matrix, and at every finite entry at least one term is tight.  If one term
were tight at distinguished cells `i` and `j`, separability would give

```math
\widetilde M(x_i,y_i)+\widetilde M(x_j,y_j)
\ge
\widetilde M(x_i,y_j)+\widetilde M(x_j,y_i).        \tag{TR.4}
```

The right side minus the left side is at least

```math
G-4\|\widetilde M-M\|_\infty>0,                   \tag{TR.5}
```

a contradiction.  Hence every distinguished cell needs a different tight
term.

All four errors are necessary in this argument, and the strict threshold is
globally sharp.  For

```math
M=\begin{pmatrix}0&1\\1&0\end{pmatrix},
```

one has `G=2`, whereas the constant matrix with every entry `1/2` is within
uniform error `G/4=1/2` and has min-plus factor rank one. `square`

This is a joint four-cell statement.  It does not separately bound row and
column channels, and it is not a mutual-information bound.

### Proposition TR.2 (a fooling set can have vanishing average mass)

Let `D_r` be the `r by r` matrix with zero diagonal and every off-diagonal
entry equal to one.  Then

```math
\operatorname{rank}_{\min,+}(D_r)=r,               \tag{TR.6}
```

but the all-one matrix has factor rank one and normalized entrywise
mean-square error

```math
{1\over r^2}\|D_r-\mathbf 1\|_F^2={1\over r}.      \tag{TR.7}
```

#### Proof

Theorem TR.1 at zero error gives the lower bound.  For the upper bound, use
one term for every `k`: let `u_k(i)=1_(i ne k)` and
`v_k(j)=1_(j ne k)`.  Their pointwise minimum is zero on the diagonal and
one off it.  The rank-one approximant differs only on the `r` diagonal
cells, proving (TR.7). `square`

Consequently, no lower bound under uniformly averaged entrywise loss can
depend only on the cardinality and gap of an unweighted tropical fooling set.
One also needs query mass, a distributional expansion property, or an
adversarial query that exposes the distinguished cells.  This is the tropical
version of the zero-entropy obstruction: exact and uniform response
complexity can live on a vanishing fraction of the query table.

## 6. Stable exponential feature growth for linear codes

Let `C<=F_2^m`, split the coordinates as `L disjoint_union R`, and define its
conditional distance table

```math
W(x_L,x_R)=d((x_L,x_R),C).
```

Put

```math
s=\dim C-\dim C_L-\dim C_R.                        \tag{TR.8}
```

Sheshadri's 2026 theorem proves

```math
\operatorname{rank}_{\min,+}(W)=2^s.               \tag{TR.9}
```

More precisely, its lower-bound proof selects one lifted codeword from every
class of `P_R(C)/C_R`.  The resulting `2^s by 2^s` block has diagonal zero
and every off-diagonal entry at least one.

### Corollary TR.3 (sub-half-error rank is exact)

For every `0<=epsilon<1/2`,

```math
\boxed{
\min_{\|\widetilde W-W\|_\infty\le\epsilon}
\operatorname{rank}_{\min,+}(\widetilde W)=2^s.}   \tag{TR.10}
```

#### Proof

For two different selected cells, the two cross entries are each at least
one and the two diagonal entries are zero.  Thus (TR.2) holds with `G=2`, so
Theorem TR.1 supplies the lower bound for every `epsilon<1/2`.  Taking
`Mtilde=W` and using (TR.9) supplies the upper bound. `square`

The exact theorem and transversal already occur in the main project ledger
and in Karthik Sheshadri,
[*Trellis State Complexity as an Exact Tropical Factorization Rank*](https://arxiv.org/abs/2607.23471),
Theorem 1 and Lemmas 2--3.  The source explicitly lists approximate
factorization as an open question (Remark 7(ii)); it does not state the
sub-half-error corollary.  A targeted search found general work on tropical
factor ranks and numerical approximate tropical factorization, but no source
for this precise robust fooling-set statement.  Because the proof is an
elementary stability argument, external novelty should nevertheless be
claimed cautiously.

The result settles only the lattice-scale regime of that open question.  If
distances are divided by block length `m`, its error threshold is
`1/(2m)`, not a macroscopic constant.  It says nothing about raw error at
least `1/2`, error growing with block length, relative error, or average
error.  Those regimes require new ideas.

## 7. What survives as theory

The experiments separate three objects which should not be conflated.

| Object | What it controls | Composition verdict |
|---|---|---|
| Boundary response kernel | all endpoint contexts | exact max-plus closure |
| Response-separation polytope | posterior variance/information for a fixed response embedding | exact for orthogonal products; no same-space or max-plus rule from child polytopes alone |
| Tropical factor rank | number of separable channels needed to represent a conditional kernel | stable `2^s` law for linear-code tables below error `1/2` |

The first and third together constitute a genuinely sub-landscape feature
algebra when the code state dimension `s` is small: `2^s` conditional channels
replace the full `2^|L| by 2^|R|` table, compose by min-plus elimination, and
are optimal even under arbitrary sub-half-error factorization.  When
`s=Theta(m)`, the same theorem is a structured exponential information
obstruction.

The response-separation polytope remains generative as the geometric input to
the sharp posterior-width information inequality.  It does **not** survive as
a self-contained feature algebra.  The collinear/orthogonal example proves
that this is an information loss in the object itself, not a missing proof.

The robust tropical corollary is the strongest result of this report.  It
turns an exact algebraic obstruction into a stable approximation theorem and
answers a published approximate-factorization question in a nontrivial
regime.  Its underlying crossing proof is classical in spirit, so the proper
claim is a generative theorem for this program, not a declaration of a new
general branch of tropical mathematics.

## 8. Minimal next theorem

The next useful result is not another scalar modulus.  It is a massive or
average-error tropical fooling-set theorem:

> If a conditional response table has a positive-measure family of
> distinguished cells with a uniform four-cell gap, quantify the minimum
> number of min-plus channels needed for small mean-square response error.

This would connect stable feature-algebra growth to the posterior-width
mutual-information theorem.  It is falsifiable on the transversal blocks:
if the distinguished block has vanishing mass under the declared query law,
average error may erase it, reproducing the zero-entropy obstruction.  Any
positive theorem must therefore include an explicit mass or query-exposure
hypothesis.
