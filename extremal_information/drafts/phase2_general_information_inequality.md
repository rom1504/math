# A posterior-width inequality for extremal response maps

**Status:** theorem-builder draft.  All mathematical claims below are proved
in the draft.  No claim has yet been promoted to `theorems.md`.

## 1. Result in one paragraph

There is a useful strengthening of the current pinned-query argument which is
not a packing bound and does not require the response to be linear.  Regard all
declared query answers of a landscape as one vector `R_a` in a Hilbert space,
where `a in {-1,+1}^N` is a latent combinatorial parameter.  If changing `r`
latent coordinates always moves the complete response vector by squared
Hilbert distance at least `kappa r`, then every transcript whose decoded
response has mean-square error `Delta` satisfies

```math
I(A;Z)\ge
N\left[1-g\left(\min\left\{{4\Delta\over\kappa N},1\right\}\right)\right],
\qquad
g(v)=h_2\left({1-\sqrt{1-v}\over2}\right).             \tag{1.1}
```

The same theorem has a weighted-coordinate form.  Its proof combines a
nonlinear posterior-width inequality for the response embedding with the
sharp entropy--posterior-variance curve of a binary variable.  The geometric
constant is computable directly from exposed response faces.  The theorem:

- strengthens the existing dense-Ising pinned-query rate from
  `N[1-h_2(D)]` to `N[1-g(D)]` and remains nonzero throughout `D<1`;
- applies without linearity to nearest-code/root-distance landscapes, giving
  an exponential information lower bound for approximating all rooted
  nearest-code queries at average squared error below `1/4`; and
- is sharp for scaled isometric hypercube response maps.

The new object suggested by the proof is the **posterior width of the response
embedding**, equivalently its squared inverse-Hamming modulus.  This is not a
new general information-theoretic principle: the information step is a sharp
binary rate--distortion calculation.  What is useful for this program is that
the response geometry supplies a checkable, nonlinear distortion measure
which applies to both quadratic optimization and code covering geometry.

## 2. Setup: a response is one Hilbert-space point

Let `mathcal Q` be a finite or measurable query space with probability measure
`mu`, and let

```math
R_a(q)=V_{H_a}(q)
```

be the zero-temperature optimum of landscape `H_a` under query `q`.  We regard
`R_a` as a vector in the real Hilbert space

```math
mathcal Y=L^2(mathcal Q,mu).
```

More generally, everything below holds for an arbitrary real Hilbert space.
The latent parameter is

```math
a=(a_1,\ldots,a_N)\in\{-1,+1\}^N,
```

and `A` denotes the uniform random parameter.  A possibly randomized summary
is a transcript `Z`; a decoder produces `Rhat_Z in mathcal Y`.  Its response
distortion is

```math
Delta=mathbb E\|R_A-\widehat R_Z\|_{mathcal Y}^2.       \tag{2.1}
```

All logarithms and mutual informations are in bits.

The response map carries the following geometry.

### Definition 2.1 (squared inverse-Hamming modulus)

```math
kappa(R)=
\min_{a\ne b}
{\|R_a-R_b\|_{mathcal Y}^2\over d_H(a,b)}.             \tag{2.2}
```

This is positive exactly when a finite response map is injective.  Its scale,
not just its positivity, matters.  It measures the least average squared
movement of all exposed query values per changed latent coordinate.

A more informative anisotropic object is the response-separation polytope

```math
Gamma(R)=\left\{gamma\in[0,infinity)^N:
 \sum_{e:a_e\ne b_e}gamma_e
 \le \|R_a-R_b\|_{mathcal Y}^2
 \text{ for every }a,b\right\}.                       \tag{2.3}
```

For a finite response table, `Gamma(R)` is an explicit linear program.  The
constant vector `(kappa(R),...,kappa(R))` belongs to it.

## 3. The binary entropy--variance lemma

Define, for `0<=v<=1`,

```math
g(v)=h_2\left({1-\sqrt{1-v}\over2}\right).             \tag{3.1}
```

### Lemma 3.1 (sharp entropy at fixed posterior variance)

Let `X in {-1,+1}` and put `w=mathbb E X` and
`v=Var(X)=1-w^2`.  Then

```math
H(X)=g(v).                                             \tag{3.2}
```

The function `g` is increasing and concave on `[0,1]`.  Consequently, if
`A` is uniform on `{-1,+1}^N`, `Z` is arbitrary, and

```math
v_e(z)=Var(A_e\mid Z=z),
\qquad
\bar v_e=mathbb E_Zv_e(Z),
```

then

```math
H(A\mid Z)
\le\sum_{e=1}^N g(\bar v_e)
\le N g\left({1\over N}\sum_e\bar v_e\right).        \tag{3.3}
```

#### Proof

The two probabilities of `X` are `(1+w)/2` and `(1-w)/2`; binary entropy is
symmetric, so (3.2) follows from `|w|=sqrt(1-v)`.

For concavity it is convenient to use natural entropy and divide by `log 2`
at the end.  Set `t=sqrt(1-v)`.  Direct differentiation gives

```math
g'(v)={operatorname{arctanh}t\over2t\log2},
```

and, for `0<t<1`,

```math
g''(v)=
-{t/(1-t^2)-operatorname{arctanh}t\over4t^3\log2}<0.  \tag{3.4}
```

The numerator is positive because its derivative is
`2t^2/(1-t^2)^2` and it vanishes at zero.  Endpoint values follow by
continuity.  Thus `g` is increasing and concave.

For each `z`, entropy subadditivity and (3.2) give

```math
H(A\mid Z=z)\le\sum_e g(v_e(z)).
```

Average over `z`, use concavity once in each coordinate, and then use it once
across coordinates.  This proves (3.3). `square`

This is the squared-error analogue of the binary entropy/error relations
behind Fano-type bounds.  The proof is included because the exact curve, not
only a Hamming-error surrogate, is important here.

## 4. Posterior width equals inverse-Hamming response geometry

For a probability distribution `pi` on the latent cube, write

```math
Var_pi(R)=mathbb E_pi\|R_A-mathbb E_piR_A\|^2,
\qquad
V_pi(A)=\sum_e Var_pi(A_e).                            \tag{4.1}
```

### Lemma 4.1 (nonlinear posterior-width lemma)

For every response map and every posterior `pi`, one has

```math
Var_pi(R)\ge {kappa(R)\over4}V_pi(A).                 \tag{4.2}
```

More generally, for every `gamma in Gamma(R)`,

```math
Var_pi(R)\ge {1\over4}\sum_e gamma_e Var_pi(A_e).     \tag{4.3}
```

The factor `1/4` in (4.2) is optimal.  Indeed,

```math
\inf_pi {Var_pi(R)\over V_pi(A)}={kappa(R)\over4},    \tag{4.4}
```

where the infimum is over posteriors with nonzero denominator.

#### Proof

Let `A,A'` be independent draws from `pi`.  The Hilbert variance identity and
the definition of `Gamma(R)` give

```math
Var_pi(R)
={1\over2}\mathbb E\|R_A-R_{A'}\|^2
\ge {1\over2}\sum_e gamma_e
       mathbb P\{A_e\ne A'_e\}.                       \tag{4.5}
```

If `p_e=mathbb P_pi(A_e=1)`, then

```math
mathbb P\{A_e\ne A'_e\}=2p_e(1-p_e)
={1\over2}Var_pi(A_e).
```

This proves (4.3) and hence (4.2).

For the reverse inequality in (4.4), choose a pair `a,b` minimizing (2.2) and
let `pi` put mass `1/2` at each.  Then

```math
Var_pi(R)={1\over4}\|R_a-R_b\|^2,
\qquad
V_pi(A)=d_H(a,b),
```

which gives equality. `square`

Thus `kappa/4` is not an arbitrarily chosen frame constant: it is exactly the
least amount of response variance forced by one unit of posterior coordinate
variance.  The statement is nonlinear and remains valid when optimizers and
exposed faces change with `a`.

## 5. General extremal information inequality

### Theorem 5.1 (weighted posterior-width rate bound)

Let `A` be uniform on `{-1,+1}^N`, let `Z` be any transcript, and let
`Rhat_Z` have distortion `Delta` as in (2.1).  For every
`gamma in Gamma(R)`,

```math
I(A;Z)\ge
N-
\max\left\{
 \sum_{e=1}^N g(v_e):
 0\le v_e\le1,
 \sum_e gamma_ev_e\le4\Delta
\right\}.                                             \tag{5.1}
```

In particular, if `kappa=kappa(R)>0`, then

```math
I(A;Z)\ge
N\left[1-g\left(
\min\left\{{4\Delta\over kappa N},1\right\}
\right)\right].                                      \tag{5.2}
```

The result permits randomized encoders and decoders and average, rather than
uniform, query error.

#### Proof

Condition on `Z=z`.  The conditional mean response minimizes squared Hilbert
loss.  Therefore

```math
Delta
\ge mathbb E_Z Var(R_A\mid Z)
\ge {1\over4}\sum_e gamma_e\bar v_e,                  \tag{5.3}
```

where the second inequality is Lemma 4.1 and
`bar v_e=mathbb E Var(A_e|Z)`.  Lemma 3.1 gives

```math
I(A;Z)=N-H(A\mid Z)
\ge N-\sum_e g(\bar v_e).                             \tag{5.4}
```

The vector `(bar v_e)` is feasible in (5.1), proving that bound.  For constant
weights `gamma_e=kappa`, Jensen's inequality maximizes the entropy term at
equal `v_e`, until the constraint reaches `v_e=1`.  This gives (5.2).
`square`

### Sharpness

Let `R_a=s a` in Euclidean `N`-space.  Then `kappa=4s^2`.  Pass each bit
independently through a binary symmetric channel with crossover `p<=1/2`, and
decode by the posterior mean.  One obtains

```math
I(A;Z)=N[1-h_2(p)],
\qquad
Delta=4s^2Np(1-p).
```

Since `g(4p(1-p))=h_2(p)`, equality holds in (5.2).  Therefore neither the
curve nor the factor four can be improved under only the stated geometric
hypothesis.

The weighted envelope is sharp as well.  Given nonnegative weights
`gamma_e`, take

```math
R_a=\left({\sqrt{\gamma_e}\over2}a_e\right)_{e=1}^N.
```

Then pairwise squared distance is exactly the sum of `gamma_e` over changed
coordinates.  Independent binary symmetric channels with posterior
variances `v_e` attain

```math
\Delta={1\over4}\sum_e\gamma_ev_e,
\qquad I(A;Z)=N-\sum_e g(v_e).
```

Thus no stronger weighted rate curve follows from membership in
`Gamma(R)` alone.

### Why this is more than response packing

A packing argument uses a threshold: below half the minimum separation it
forces exact identification, and above that threshold it often says nothing.
Theorem 5.1 instead charges every surviving posterior coordinate variance. It
gives a continuous, sharp rate curve under mean-square response distortion,
including regimes in which exact decoding has large error.  Its model-specific
content is the proof of a nonvanishing response-width modulus, not the
cardinality of a separated set.

## 6. Exposed-face and query-frame certificates

The constant in Theorem 5.1 can be certified from exposed optimizers.
Suppose queries have the form

```math
R_a(q)=\max_x\{H_a(x)+J_q(x)\},                       \tag{6.1}
```

and a state `x_q` is an optimizer for every latent `a`.  Then

```math
R_a(q)-R_b(q)=H_a(x_q)-H_b(x_q).                      \tag{6.2}
```

Consequently, any inequality

```math
mathbb E_q[H_a(x_q)-H_b(x_q)]^2
\ge\sum_{e:a_e\ne b_e}gamma_e                        \tag{6.3}
```

places `gamma` in `Gamma(R)`.  Strict uniqueness is unnecessary; only the
common exposed value is used.

An important special case is an affine exposed chart.  After removing a known
term and, if operationally justified, applying an orthogonal projection that
forgets declared nuisance directions, suppose

```math
R_a=r_0+\sum_e a_e psi_e,
\qquad psi_e\in mathcal Y.                            \tag{6.4}
```

Let `G_ef=<psi_e,psi_f>`.  Since `(a-b)/2` ranges over all nonzero vectors in
`{0,+1,-1}^N`,

```math
kappa(R)=
4\min_{0\ne t\in\{0,+1,-1\}^N}
{t^TGt\over |supp(t)|}.                               \tag{6.5}
```

This is a discrete restricted eigenvalue, not the ordinary least eigenvalue of
`G`.  It tests cancellation among all jointly evaluated channels before an
absolute value or information charge is applied.  If the `psi_e` are
orthogonal with squared norms `s_e^2`, then

```math
(4s_1^2,\ldots,4s_N^2)\in Gamma(R).                  \tag{6.6}
```

Conversely, a sign combination with small norm is an immediate obstruction to
a proposed coordinatewise information price.  This is exactly where hidden
feature-algebra cancellation enters the inequality.

Any contraction `P:mathcal Y -> mathcal Y_0` may be applied before computing
`Gamma`: project the decoded response as well, so its distortion cannot
increase.  A lower bound for this contracted task is therefore a valid lower
bound for the original task.  This does not assert that the contracted
response is sufficient for the original experiment.  Parameter-dependent
terms must not be called unobservable merely because a convenient projection
removes them.

## 7. Application I: dense quadratic Ising responses

Let `E_n` have `N=binom(n,2)` edges, fix `a>0`, and define as in the current
theory

```math
q_A(x)=a\sum_{i<j}A_ijx_ix_j,
\qquad
c_A=\max_xq_A(x),
\qquad
H_A=q_A-c_A.                                         \tag{7.1}
```

Every unqueried maximum is zero.  There are two equivalent exposing query
families.

1. Coordinate fields `J_u(x)=M<u,x>` with `M>a(n-1)` expose `u`.
2. Bounded-per-edge rank-one interactions

   ```math
   J_u(x)=L\sum_{i<j}u_iu_jx_ix_j,
   \qquad L>a,                                       \tag{7.2}
   ```

   expose the pair `{u,-u}`.  Indeed, after writing `y_i=u_ix_i`, a state with
   `k` negative `y_i` loses `2Lk(n-k)` in the query and can gain at most
   `2ak(n-k)` in the latent interaction.

For either family, after subtracting its known query constant,

```math
R_A(u)=q_A(u)-c_A.                                    \tag{7.3}
```

Use the uniform measure on `u in {-1,+1}^n`.  The Walsh characters
`u_i u_j` are orthonormal and orthogonal to constants, hence

```math
\|R_A-R_B\|_2^2
=4a^2d_H(A,B)+(c_A-c_B)^2
\ge4a^2d_H(A,B).                                     \tag{7.4}
```

Thus `kappa>=4a^2`.  Equality holds: switching one vertex changes exactly
`n-1` edge signs but only relabels the spin maximum, so `c_A=c_B` and (7.4)
has ratio exactly `4a^2`.  If a decoded query-response function has

```math
Delta=mathbb E_{A,Z,U}
[\widehat R_Z(U)-R_A(U)]^2,
```

then Theorem 5.1 gives

```math
I(A;Z)\ge
N\left[1-g\left(
\min\left\{{Delta\over a^2N},1\right\}
\right)\right].                                      \tag{7.5}
```

This strictly strengthens the existing `N[1-h_2(D)]` estimate on its stated
range `0<D<=1/2`, with `D=Delta/(a^2N)`, because

```math
{1-\sqrt{1-D}\over2}<D
```

whenever `D<3/4`.  The new bound moreover remains positive all the way to
`D<1`.
More directly, it is the correct binary squared-error curve whereas the old
proof first hard-decoded Walsh coefficients and paid Hamming error.

The rank-one version shows that fields of coordinate magnitude `Theta(n)` are
not essential.  It still uses a dense query of total interaction size
`Theta(n^2)`, so it does not settle the sparse- or low-total-budget query
problem.

## 8. Application II: nearest-code rooted responses

Let `Omega_m={-1,+1}^m`, let `K=2^m`, and fix an anchor `o in Omega_m`.
For each latent vector `a in {-1,+1}^{K-1}`, form the arbitrary (not
necessarily linear) code

```math
C_a={o}\cup\{u\ne o:a_u=+1\}.                        \tag{8.1}
```

There are `N=K-1` independent membership coordinates.  Give the code the
nearest-distance landscape

```math
H_{C_a}(x)=-d_H(x,C_a).                               \tag{8.2}
```

For a root `u`, query by `J_u(x)=M<u,x>` with any fixed `M>1/2`.  Hamming
distance to a nonempty set is one-Lipschitz, so if `d_H(x,u)=k`,

```math
H_{C_a}(x)-H_{C_a}(u)\le k,
\qquad
J_u(x)-J_u(u)=-2Mk.
```

Therefore `u` is the unique exposed state and

```math
R_a(u)=Mm-d_H(u,C_a).                                 \tag{8.3}
```

This is the complete rooted nearest-code profile.  If `C_a` and `C_b` differ
in `r=d_H(a,b)` membership coordinates, then at every root in their symmetric
difference one distance is zero and the other is at least one.  Under the
uniform root measure,

```math
\|R_a-R_b\|_2^2\ge {r\over K}.                        \tag{8.4}
```

Hence `kappa>=1/K`.  Equality holds: take one code to be `Omega_m` and the
other to omit a single nonanchor word.  Their distance profiles differ by one
only at the omitted word.  Thus

```math
kappa={1\over K}.                                     \tag{8.5}
```

For average squared root-response error `Delta`, Theorem 5.1 yields

```math
I(A;Z)\ge
(K-1)\left[1-g\left(
\min\left\{{4KDelta\over K-1},1\right\}
\right)\right].                                      \tag{8.6}
```

Consequences include:

- every fixed `Delta<1/4` has a positive exponential-in-`m` information
  price for all sufficiently large `m` (the exact finite threshold is
  `Delta<(K-1)/(4K)`);
- if `Delta=o(1)`, then `I(A;Z)=K-o(K)` bits; and
- uniform additive response error below `1/2` recovers every membership bit
  exactly by testing whether the decoded distance is below `1/2`.

This is genuinely nonlinear: distance to a code is not an affine function of
its membership vector.  The common theorem applies because the entire rooted
response map is inverse-Hamming separated.  It also quantifies why adding all
roots restores sufficiency only by essentially reconstructing an arbitrary
code.

## 9. Optional third application: counterfactual Max-Cut

Let the latent graph `B` be uniform on `{0,1}^N`, where
`N=binom(n,2)`, and

```math
C_B(u)=\sum_{i<j}B_ij{1-u_iu_j\over2}.
```

Add a vertex-prize query that gives reward `M` per vertex agreeing with `u`,
with `M>(n-1)/2`.  A one-vertex change loses `2M` in the query and can gain at
most `n-1` in the cut term, so it exposes `u` and the response is
`R_B(u)=Mn+C_B(u)`.  Give roots the uniform measure.  For two graphs,
orthogonality of the constant and degree-two Walsh terms gives

```math
\|R_B-R_D\|_2^2
={1\over4}\left(\sum_e(B_e-D_e)\right)^2
 +{1\over4}d_H(B,D)
\ge {1\over4}d_H(B,D).                               \tag{9.1}
```

Thus `kappa>=1/4` (in fact equality holds for `n>=3`).  For an arbitrary
transcript and decoded response define

```math
\Delta=\mathbb E_{B,Z,U}
[\widehat R_Z(U)-R_B(U)]^2.
```

Then

```math
I(B;Z)\ge
N\left[1-g\left(
\min\left\{{16Delta\over N},1\right\}
\right)\right].                                      \tag{9.2}
```

This gives the same continuous strengthening for counterfactual Max-Cut
responses, although it is algebraically close to the Ising application.

## 10. Boundary cases, hidden assumptions, and falsifiers

### Boundary cases

1. **`kappa=0`.**  Two latent landscapes are response-equivalent.  No theorem
   can charge all `N` bits from this query family.  For example, the scalar
   response `R_a=sum_e a_e` has many collisions.
2. **Vanishing scale.**  Injectivity alone is meaningless asymptotically.  If
   `kappa_N` is exponentially small, a normalized error can erase the whole
   latent cube.  Every application must state the scale of `kappa_N`.
3. **Large distortion.**  At `4Delta>=kappa N`, (5.2) correctly becomes zero:
   the posterior coordinate-variance budget permits the no-information
   posterior.
4. **Nonuniform prior.**  Uniform independent latent bits are used in
   `H(A)=N` and in the symmetric function `g`.  A nonuniform/product prior has
   an analogous but different entropy--variance allocation problem.  No such
   extension is claimed here.
5. **Query measure.**  Rescaling `mu` or query weights rescales both `Delta`
   and `kappa`.  They must use the same declared measure.
6. **Decoder randomness.**  It causes no issue: append all decoder randomness
   to `Z`.
7. **Projection.**  A contraction of responses is safe for proving a lower
   bound.  Silently deleting a parameter-dependent observable term is not.
8. **Exposure.**  Uniform exposure is only a certificate for `kappa`; Theorem
   5.1 itself does not assume common optimizers.  If optimizers change, compute
   the response distances directly.

### Exact finite falsifier

A claimed homogeneous constant `kappa_0` is false precisely when one pair
satisfies

```math
\|R_a-R_b\|_2^2<kappa_0d_H(a,b).                     \tag{10.1}
```

For a finite query table this is exhaustively decidable.  The weighted claim
`gamma in Gamma(R)` is a finite linear feasibility problem, and a violated
pair is a certificate.  The following minimal check verifies the sharp code
constant for all anchored codes through `m=3`:

```python
for m in range(1, 4):
    K = 1 << m
    tables = []
    for mask in range(1 << (K - 1)):
        C = {0} | {u for u in range(1, K)
                   if (mask >> (u - 1)) & 1}
        d = [min(bin(u ^ c).count("1") for c in C)
             for u in range(K)]
        tables.append((C, d))
    ratio = min(
        sum((x-y)**2 for x, y in zip(d, e)) / K / len(C ^ D)
        for i, (C, d) in enumerate(tables)
        for D, e in tables[i+1:]
    )
    assert ratio == 1 / K
```

The outputs are respectively `1/2`, `1/4`, and `1/8`.  The general proof in
(8.4)--(8.5) makes the finite pattern rigorous for every `m`.

The factor `1/4` in the posterior-width lemma has its own finite extremizer:
put posterior mass `1/2` on a pair attaining (2.2).  Any proposed improvement
under only pairwise response separation is therefore immediately false.

## 11. Relation to existing theory

- Shannon's 1959 rate--distortion theory supplies the general operational
  setting for minimizing mutual information at declared distortion
  ([primary source](https://ieeexplore.ieee.org/document/5311476)).
- Feder and Merhav study sharp relations between entropy and probability of
  error ([IEEE Transactions on Information Theory 40 (1994),
  259--266](https://doi.org/10.1109/18.272494)).  Lemma 3.1 instead uses the
  exact posterior variance of a binary coordinate; its short proof is given
  here.
- The language of inverse-Hamming response geometry is adjacent to nonlinear
  metric embeddings and metric cotype; see Mendel--Naor
  ([Annals of Mathematics 168 (2008),
  247--298](https://doi.org/10.4007/annals.2008.168.247)).  No embedding
  theorem from that paper is imported.  Here the image is the concrete
  zero-temperature response map and the direction of use is an information
  lower bound.

The theorem is therefore built from classical ingredients, but the
combination is not merely the definition of mutual information or a response
packing.  Its checkable mathematical obligation is:

```math
\boxed{
\text{prove a correctly scaled inverse-Hamming modulus for the exposed
response map.}
}
```

That obligation makes sense for nonlinear code-distance responses as well as
for an orthogonal Walsh response frame.

## 12. Director checkpoint

### What this newly formulates

The previous framework had two separate lower bounds: hard-decoding pinned
Walsh coefficients, and posterior sign polarization.  Theorem 5.1 identifies
the common missing object: the least response variance forced by posterior
uncertainty in latent coordinates.  Lemma 4.1 shows that this object is exactly
one quarter of a nonlinear inverse-Hamming modulus, and not an ad hoc moment.
This simultaneously:

- strengthens the quantitative dense-Ising result;
- handles nonlinear nearest-code landscapes with the same theorem;
- makes joint cancellation visible through the discrete restricted
  eigenvalue (6.5); and
- supplies an exact finite falsifier for proposed information lower bounds.

### Is this only repackaged information theory?

The entropy step is classical binary rate--distortion geometry and should be
presented as such.  The generative part is the response-width reduction and
the proof of its modulus in two unrelated extremal models.  This is a real
Level-3 theorem within the developing framework, but not by itself a new
field.  It becomes more than a useful lemma if response-width constants admit
composition, tensorization, or deterministic synchronization laws.

### Strongest next theorem suggested by this result

For a declared composition operation, determine whether response separation
weights tensorize or obey a data-processing law.  A precise first target is:

> Given child response maps with `gamma^(1) in Gamma(R_1)` and
> `gamma^(2) in Gamma(R_2)`, identify nontrivial hypotheses on a bi-affine
> coupling under which the parent response-separation polytope contains
> `(c gamma^(1), c gamma^(2))` with `c` bounded away from zero independently
> of composition depth.

A proof would turn the present one-shot information inequality into a feature-
algebra growth theorem.  A counterexample with `c` decaying exponentially
would quantify why a composition hierarchy reconstructs the landscape.
