# Rigorous results

This file contains proved statements only.  Exact finite computations belong
in [`examples.md`](examples.md).  Logarithms are natural except in explicitly
information-theoretic statements, where entropy is in bits.

## 1. Positive quadratic tail entropy reaches the maximum

Let `f:{-1,+1}^n -> R` be a homogeneous quadratic form and let

```math
M=\max_x f(x)>0.
```

### Theorem 1.1 (universal high-energy cloud)

For every `0<theta<1`,

```math
\liminf_{n\to\infty}{1\over n}
\log\#\{x:f_n(x)\ge\theta M_n\}
\ge
h\left({1-\sqrt\theta\over2}\right),                 \tag{1.1}
```

for every sequence of homogeneous quadratics `f_n` with `M_n=max f_n>0`,
where `h(p)=-p log p-(1-p)log(1-p)`.  If `M_n=0`, cube averaging forces the
quadratic to vanish identically, and every state lies at the edge.

More precisely, if

```math
0<\delta<{1-\sqrt\theta\over2},
```

a positive fraction of the Hamming sphere of radius `floor(delta n)` around
a maximizer is above `theta M_n`, for all sufficiently large `n`.

#### Proof

Fix a maximizer `x*`.  Flip a uniformly random `r`-subset of its coordinates,
writing the result as `X=x*z`, where exactly `r` entries of `z` are negative.
For `i!=j`,

```math
\mathbb E z_i z_j
=\lambda_{n,r}
:={(n-2r)^2-n\over n(n-1)}.
```

Homogeneity gives `E f_n(X)=lambda_(n,r) M_n`.  If
`p=P{f_n(X)>=theta M_n}`, then `f_n(X)<=M_n` implies

```math
\lambda_{n,r}M_n
\le pM_n+(1-p)\theta M_n,
```

so

```math
p\ge {\lambda_{n,r}-\theta\over1-\theta}.
```

For `r/n -> delta<(1-sqrt(theta))/2`, the right side stays positive.  The
sphere has size `binom(n,r)=exp((h(delta)+o(1))n)`.  Letting `delta` increase
to the stated endpoint proves (1.1). `square`

The same statement holds for `max_x|f_n(x)|` after replacing `f_n` by the
sign that attains the absolute cap.

### Corollary 1.2 (positive tail entropy determines the edge)

Let `a_n>0`, pass to a subsequence on which

```math
{\max_x f_n(x)\over a_n}\longrightarrow m,
```

and define

```math
s^\uparrow(t)=\liminf_n {1\over n}
\log\#\{x:f_n(x)\ge t a_n\}.
```

Then

```math
m=\sup\{t:s^\uparrow(t)>0\}.                          \tag{1.2}
```

For `t>m` the tail is eventually empty.  For every `t<m`, Theorem 1.1
supplies an exponentially large cloud above `t a_n`.  Thus resolved scalar
upper-tail entropy cannot coexist with a different quadratic normalized
maximum.  This theorem was previously proved in the main repository's
`entropy_energy_dichotomy.md`; its interpretation here is new, not the noise
argument itself.

## 2. Exact global pair data can miss labeled coupling response

Let `n=2m` with even `m`, and split the coordinates into fixed halves `L,R`.
Put

```math
Q_L(x)={ (\sum_{i\in L}x_i)^2-m\over m},
\qquad
Q_R(x)={ (\sum_{i\in R}x_i)^2-m\over m},
```

and `H_L=Q_L/(m-1)`, `H_R=Q_R/(m-1)`.  Both maxima equal one.

### Theorem 2.1 (block-location obstruction)

For every finite set of energy and overlap bins, indeed pointwise for every
exact triple, the counts

```math
#\{(x,y):(H_L(x),H_L(y),n^{-1}\langle x,y\rangle)
=(e_1,e_2,q)\}
```

and the corresponding counts for `H_R` are equal.  Nevertheless, with

```math
R_L(x,y)={1\over m}\sum_{i\in L}x_i y_i,
```

the fixed-left response satisfies

```math
{1\over2}\max_{R_L(x,y)=0}\{H_L(x)+H_L(y)\}
={m-2\over2(m-1)}\longrightarrow {1\over2},          \tag{2.1}
```

while the same expression for `H_R` is exactly one.

#### Proof

The coordinate permutation exchanging `L` and `R` maps `H_L` to `H_R` and
preserves global overlap, proving equality of every exact pair fiber.

For the response of `H_L`, set

```math
u={1\over m}\sum_{i\in L}x_i,
\qquad
v={1\over m}\sum_{i\in L}y_i.
```

Global sign flips allow `u,v>=0` while preserving the constraint
`R_L=0`.  The fraction of coordinates with `x_i=y_i=-1` is
`(1-u-v)/4`, hence `u+v<=1` and `u^2+v^2<=1`.  Therefore

```math
Q_L(x)+Q_L(y)=m(u^2+v^2)-2\le m-2.
```

Equality is attained with `x_L` constant and `y_L` balanced.  For `H_R`,
make both right blocks constant and choose left blocks of zero mutual overlap;
both energies attain one. `square`

The fixed block is part of the query.  The two landscapes become isomorphic
if the apparatus is relabeled with them.  The theorem says that unrooted
global pair data is not sufficient for an anchored experiment.

## 3. A tensor-stable coding obstruction

Let `Q_n={0,1}^n` with Hamming distance.  For a code `A subset Q_n`, write
`H_A(x)=1_A(x)` and

```math
I_A(j)=#\{(a,b)\in A^2:d(a,b)=j\}.
```

### Lemma 3.1 (inner distribution determines the ambient pair census)

If `|A|=s`, then for `epsilon,eta in {0,1}`,

```math
N_A^{11}(j)=I_A(j),
```

```math
N_A^{10}(j)=N_A^{01}(j)=s\binom nj-I_A(j),
```

and

```math
N_A^{00}(j)=(2^n-2s)\binom nj+I_A(j),                 \tag{3.1}
```

where `N_A^(epsilon eta)(j)` counts ambient ordered pairs with those two
membership energies and distance `j`.

#### Proof

Each fixed codeword has `binom(n,j)` ambient words at distance `j`; removing
the codewords gives the mixed cells.  Subtract all three remaining cells from
the `2^n binom(n,j)` ambient ordered pairs to obtain the last formula.
`square`

Set

```math
C=\{0000,0011,0101,0110\},
\qquad
D=\{0000,0011,0101,1001\}.
```

### Theorem 3.2 (same complete pair data, different rooted extreme)

For every `r>=1`, the membership landscapes of the Cartesian powers `C^r`
and `D^r` have identical exact ambient

```math
(H(x),H(y),d(x,y))
```

counts.  Their covering radii are respectively `2r` and `3r`.

For `lambda>4r`, the rooted response

```math
V_A(z;\lambda)=
\max_x\{\lambda 1_A(x)-d(x,z)\}
```

equals `lambda-d(z,A)`.  Hence its worst-root values differ by `r`.

#### Proof

Every two distinct words in either base code are at distance two, so both
ordered inner polynomials are `4+12z^2`.  Cartesian products multiply these
polynomials.  Lemma 3.1 gives equality of the full ambient pair census.

The code `C` is a fixed coordinate times the even-parity code on three bits,
so its radius is two.  The word `1110` is distance three from every word of
`D`, while no four-bit word is farther than three, so `rho(D)=3`.  Distance
to a Cartesian product adds blockwise, giving `rho(A^r)=r rho(A)`.

When `lambda>4r`, a codeword has query value at least `lambda-4r>0`, whereas
a noncodeword has value at most zero.  Thus the optimizer is a nearest
codeword and the displayed response formula follows. `square`

The missing minimal datum for this query is the rooted map `z -> d(z,A)`.
It is lossless for the code, because its zero set is `A`; the example does not
pretend that full labeled nearest-code response is a compression.

For `t>=1`, let `T_t(A)` be the complete ambient histogram of

```math
\left((1_A(x_i))_{i=1}^t,
      (d(x_i,x_j))_{1\le i<j\le t}\right),
\qquad (x_1,\ldots,x_t)\in Q_n^t.                    \tag{3.2}
```

### Theorem 3.3 (no fixed unrooted replica order controls covering radius)

For every fixed `k`, there are two binary codes `A_0,A_1` such that

```math
T_t(A_0)=T_t(A_1)
\qquad(1\le t\le k),                                 \tag{3.3}
```

but their covering radii differ.  Their Cartesian powers retain (3.3) and
have a fixed positive normalized radius gap.  Thus no fixed order of complete
unrooted membership/overlap data determines covering radius uniformly over
binary codes.

#### Construction and proof

Choose an odd `r=2s+1` with `M=r+1>=k` and put `N=2^(r-1)`.  For
`epsilon in {0,1}`, index the `N` coordinates by the parity half-cube

```math
P_\epsilon=\{v\in\mathbb F_2^r:|v|=\epsilon\pmod2\}.
```

Let `A_epsilon subset Q_N` consist of the zero word and the `r` coordinate
functions

```math
a_i(v)=v_i,
\qquad 1\le i\le r.                                  \tag{3.4}
```

Every proper selection of at most `r` distinct codewords has the same
embedded column-pattern multiset for the two parities.  Indeed, translate one
selected word to zero.  The remaining selected rows span a space `W` of
linear forms that does not contain the parity form `1 dot v`: if the zero row
was selected, at least one coordinate form is missing; if it was not, `W`
contains only even-weight coefficient vectors while `r` is odd.  Every fiber
of the selected-pattern map therefore meets `P_0` and `P_1` equally.  A
coordinate permutation, followed by the initial translation, identifies the
two selected configurations.

To count a cell of `T_t`, use inclusion--exclusion on the positions required
to be outside the code.  Every resulting term fixes some tuple of codewords
and counts ambient completions with prescribed distances.  For `t<=M`, either
the distinct fixed codewords form a proper subset and the preceding cube
isometry equates the completion counts, or all `M` codewords occur, leaving
no additional distinct ambient point and giving the same equilateral distance
matrix.  This proves (3.3) through `t=M`.

Define

```math
S_\epsilon=
\sum_{\substack{0\le w\le r\\w=\epsilon\ (2)}}
\binom rw\max\{w,M-w\}.                              \tag{3.5}
```

At a coordinate of weight `w`, the sum of a root's distances to all `M`
codewords contributes either `w` or `M-w`.  Hence every root has nearest-code
distance at most `S_epsilon/M`.  Put `p=s mod 2` and choose, at every
coordinate in `P_p`, the root bit opposite the strict majority of the `M`
codeword bits.  There is no tie because the tie weight `s+1` has the other
parity.  The disagreement sets are all parity-`p` subsets of the `M` codeword
labels of size greater than `M/2`, a permutation-invariant family.  Every
codeword is therefore at the same distance `S_p/M`, proving

```math
\rho(A_p)={S_p\over M}.                               \tag{3.6}
```

Pairing complementary binomial terms gives

```math
S_p-S_{1-p}=\binom{2s}{s}>0.                          \tag{3.7}
```

Consequently

```math
\rho(A_{1-p})
\le\left\lfloor{S_{1-p}\over M}\right\rfloor
<{S_p\over M}=\rho(A_p).
```

Finally, the `m`-fold Cartesian product combines membership vectors by
coordinatewise AND and distance vectors by addition.  Equality of every base
`T_t` therefore persists by convolution.  Covering radii add, so the product
pair of block length `mN` has normalized gap
`|rho(A_0)-rho(A_1)|/N>0`. `square`

## 4. Query-relative response geometry

Let `Omega` be finite, `H:Omega -> R`, and `phi:Omega -> R^d`.  Put

```math
P_\phi=\operatorname{conv}\phi(\Omega)
```

and define

```math
\widehat H_\phi(u)=
\max\left\{\sum_x\lambda_xH(x):
\lambda\in\Delta(\Omega),\ \sum_x\lambda_x\phi(x)=u\right\}. \tag{4.1}
```

### Theorem 4.1 (response duality and exact quotient)

For every field `theta in R^d`,

```math
V_H(\theta):=\max_x\{H(x)+\langle\theta,\phi(x)\rangle\}
=\max_{u\in P_\phi}
\{\widehat H_\phi(u)+\langle\theta,u\rangle\},       \tag{4.2}
```

and

```math
\widehat H_\phi(u)=
\inf_{\theta\in\mathbb R^d}
\{V_H(\theta)-\langle\theta,u\rangle\}.             \tag{4.3}
```

Thus the roof and the complete linear-response function determine each other.
For `Theta=R^d`, or any query set determining the roof, this is the coarsest
deterministic exact quotient up to a one-to-one recoding.  For a restricted
query set, the minimal quotient is only `V_H|Theta`; the full roof may be
strictly richer.

#### Proof

Equation (4.2) maximizes a linear functional over the lifted convex hull of
`(phi(x),H(x))`, so it may be evaluated on its generators.  The roof is a
closed concave function on the compact polytope `P_phi`; concave
Fenchel--Moreau biconjugacy gives (4.3).  Exact summaries are precisely maps
through which the answer function factors. `square`

### Theorem 4.2 (additive and one-step bilinear composition)

For

```math
H_\oplus(x,y)=H_1(x)+H_2(y),
\qquad
\phi_\oplus(x,y)=\phi_1(x)+\phi_2(y),
```

the parent roof is the sup-convolution

```math
\widehat H_\oplus(u)=
\max_{u_1+u_2=u}
\{\widehat H_1(u_1)+\widehat H_2(u_2)\}.             \tag{4.4}
```

For a fixed bilinear coupling `B`,

```math
\max_{x,y}\{H_1(x)+H_2(y)+\phi_1(x)^TB\phi_2(y)\}
```

equals

```math
\max_{u\in P_1,v\in P_2}
\{\widehat H_1(u)+\widehat H_2(v)+u^TBv\}.          \tag{4.5}
```

#### Proof

The lifted generator set of the additive parent is the set sum of the child
generator sets; convex hull turns this into Minkowski addition and its upper
boundary gives (4.4).  For (4.5), choose independent child distributions
attaining the two roofs at `u,v`.  The displayed expression is the expectation
of the pure coupled energy under their product and hence cannot exceed its
maximum.  Point masses prove the reverse inequality. `square`

### Theorem 4.3 (bi-affine closure and data processing)

Suppose a parent has energy and feature

```math
H_P(x,y)=H_1(x)+H_2(y)+C(\phi_1(x),\phi_2(y)),
```

```math
\psi(x,y)=F(\phi_1(x),\phi_2(y)),
```

where `C` and every coordinate of `F` are separately affine in the two child
features.  For each parent field `eta`, write

```math
C(u,v)+\langle\eta,F(u,v)\rangle
=c_\eta+p_\eta^Tu+q_\eta^Tv+u^TB_\eta v.            \tag{4.6}
```

Then the child roofs determine the complete parent response:

```math
V_P^\psi(\eta)=c_\eta+
\max_{u\in P_1,v\in P_2}
\{\widehat H_1(u)+\widehat H_2(v)
 +p_\eta^Tu+q_\eta^Tv+u^TB_\eta v\}.                \tag{4.7}
```

Consequently they determine the full parent roof by Theorem 4.1.

For two child pairs `H_i,G_i` with the same features and the same composition,
put

```math
\Theta_1(\eta)=\{p_\eta+B_\eta v:v\in\phi_2(\Omega_2)\},
```

```math
\Theta_2(\eta)=\{q_\eta+B_\eta^Tu:u\in\phi_1(\Omega_1)\}.
```

Then

```math
|V_{P(H)}^\psi(\eta)-V_{P(G)}^\psi(\eta)|
\le d_{\Theta_1(\eta)}(H_1,G_1)
   +d_{\Theta_2(\eta)}(H_2,G_2).                    \tag{4.8}
```

If full roofs on common feature polytopes are within `epsilon_1,epsilon_2` in
vertical sup norm, the parent roofs are within `epsilon_1+epsilon_2`.

#### Proof

Independent mixtures work exactly as in Theorem 4.2 because every term in
(4.6) is separately affine; this proves (4.7).  For (4.8), write the parent
maximum first as a maximum over pure `y` of a child-one response in direction
`p_eta+B_eta phi_2(y)`, replace `H_1` by `G_1`, and then reverse the order to
replace child two.  The inequality `|max f-max g|<=||f-g||_infinity` gives the
claim.

For common feature polytopes, conjugacy gives

```math
\sup_\theta|V_H(\theta)-V_G(\theta)|
=\sup_u|\widehat H(u)-\widehat G(u)|.
```

Apply (4.8) for every parent direction. `square`

If a common convex feature domain `D` is closed under a fixed bi-affine
operation `F`, that operation is associative on `D`, and the cross energy
obeys on `D`

```math
C(u,v)+C(F(u,v),w)=C(v,w)+C(u,F(v,w)),               \tag{4.9}
```

then repeated composition is independent of bracketing and leaf errors grow
by at most their sum.  The full tensor feature is always closed but its
dimension multiplies; closure is a compression theorem only when the feature
algebra and roof metric entropy remain controlled.

For normalized block features or energies, the corresponding size weights
must be included in `F,C` (or mass must be retained as a state coordinate);
the unweighted sum bound above must not be imported unchanged.

### Corollary 4.4 (exact polynomial-state mean-field composition)

Let each site have a finite spin alphabet, a feature
`phi(s) in Z^d`, and an arbitrary local energy `h_i(s)`.  For fixed `d` and a
fixed matrix `J`, consider

```math
E(s_1,\ldots,s_N)=
\sum_i h_i(s_i)+
\sum_{i<j}\phi(s_i)^TJ\phi(s_j).                     \tag{4.10}
```

Retain for each block the roof over its total feature
`u=sum_i phi(s_i)`.  Merging blocks uses

```math
F(u,v)=u+v,
\qquad C(u,v)=u^TJv.
```

These operations satisfy the hypotheses and the energy cocycle, so the roof
gives an exact, bracket-independent dynamic program for the ground-state
value.  A block of `N` sites has only polynomially many attainable total
features when `d` and the feature alphabet are fixed.  Thus this is genuine
polynomial-state extremal compression for finite-rank Curie--Weiss and
mean-field Potts-type models, despite an exponential spin landscape.

### Theorem 4.5 (fixed-ambient compactness and unrestricted recovery)

Fix a compact convex feature region `K subset R^d`, an energy bound `B`, and
a norm on `R^(d+1)`.  Consider all nonempty compact convex sets

```math
U\subset K\times[-B,B]
```

that are vertically downward: `(u,t) in U` and `-B<=s<=t` imply `(u,s) in U`.
This family is compact in Hausdorff distance.  Every member is a Hausdorff
limit of truncated response bodies of finite landscapes with features in `K`
and energies in `[-B,B]`.  Moreover,

```math
|h_U(\theta,1)-h_V(\theta,1)|
\le\|(\theta,1)\|_*d_H(U,V).                         \tag{4.11}
```

#### Proof

The hyperspace of nonempty compact subsets of a compact metric space is
compact.  Convexity and downward closure are preserved under Hausdorff
limits, so this family is closed.  For recovery, take a finite `delta`-net
`S subset U`.  Its convex hull lies in `U` and remains a `delta`-net; its
downward extension is still contained in `U`.  Treat each `(u,t) in S` as a
finite state with feature `u` and energy `t`.  Its truncated response body is
that downward extension and is `delta`-close to `U`.  The final inequality is
the standard Lipschitz bound for support functions. `square`

This is unrestricted fixed-interface recovery.  It gives neither constrained
sign-matrix realization nor stability of exposed optimizer faces.

### Proposition 4.6 (full Boolean pinning is lossless)

For `Omega={-1,+1}^n` and `phi(x)=x`,

```math
\widehat H_\phi(x)=H(x)
\qquad(x\in\{-1,+1\}^n).                             \tag{4.12}
```

Indeed, every Boolean vector is an extreme point of the cube, and the only
distribution on the cube with mean `x` is the point mass at `x`.  Therefore
the full response roof retains the complete landscape table despite having a
feature vector of dimension only `n`.

## 5. Extremal rate--distortion lower bounds

For a landscape class and declared query set `Theta`, let

```math
d_\Theta(H,G)=\sup_{\theta\in\Theta}|V_H(\theta)-V_G(\theta)|.
```

### Proposition 5.1 (metric-entropy sandwich)

If `K_epsilon` is the smallest number of messages in a deterministic summary
from which every response is decoded with uniform error at most `epsilon`,
then

```math
\operatorname{Pack}(\mathcal H,d_\Theta,2\epsilon)
\le K_\epsilon
\le\operatorname{Cov}(\mathcal H,d_\Theta,\epsilon), \tag{5.1}
```

where the packing uses pairwise distance strictly greater than `2epsilon`.

Two landscapes sharing a message are within `2epsilon` by the triangle
inequality; an `epsilon`-net supplies the decoder. `square`

Let `E_n` be the edge set of the complete graph, `N=binom(n,2)`, and fix
`a>0`.  For `A in {-1,+1}^N`, put

```math
q_A(x)=a\sum_{i<j}A_{ij}x_ix_j,
\qquad
c_A=\max_xq_A(x),
\qquad
H_A=q_A-c_A.                                         \tag{5.2}
```

Every landscape in this `2^N`-element family has maximum zero.  Query it at
the fields `h^u=M u`, where `M>a(n-1)`.

### Theorem 5.2 (quadratic pinned-query rate)

Let `A` be uniform, let `Z` be any possibly randomized transcript, and let a
decoder return `Vhat_Z(h^u)` for all `u`.  If

```math
{1\over a^2N}
\mathbb E_{A,Z,U}
[\widehat V_Z(h^U)-V_{H_A}(h^U)]^2\le D\le {1\over2},
```

then

```math
I(A;Z)\ge N[1-h_2(D)].                               \tag{5.3}
```

In particular, uniform additive error strictly below `a` recovers every edge
sign and requires at least `N` bits.

#### Proof

If `x` differs from `u` in `k` positions, its field loses `2Mk`; at most
`k(n-k)` interaction terms change and their total possible gain is
`2ak(n-k)`.  Thus `u` is the unique optimizer and

```math
V_{H_A}(h^u)=Mn+q_A(u)-c_A.                           \tag{5.4}
```

The degree-two Walsh coefficient of `V_{H_A}(h^u)-Mn` at character `u_i u_j`
is `aA_ij`; the unknown shift `c_A` is degree zero.  Bessel's inequality
therefore turns the mean-square response error into an edge-sign estimator
whose expected Hamming error fraction is at most `D`.  Entropy subadditivity
and concavity of binary entropy give

```math
I(A;Z)
\ge N-\sum_{i<j}h_2(P\{\widehat A_{ij}\ne A_{ij}\})
\ge N[1-h_2(D)].
```

Under uniform error below `a`, every Walsh coefficient error is below `a`, so
all signs are exact. `square`

The family is `2a`-separated in the query pseudometric, equivalently an
`r`-packing for every `r<2a`; it need not be a strict `2a`-packing.

### Theorem 5.3 (posterior sign-polarization price)

Let `A` be uniform on `N` independent signs, let `Z` be any transcript, and
put

```math
w_e(Z)=\mathbb E[A_e\mid Z],
\qquad
\mathcal V_Z=\sum_e(1-w_e(Z)^2).
```

Then

```math
I(A;Z)\ge
N\left[1-h_2\left({\mathbb E\mathcal V_Z\over2N}\right)\right]. \tag{5.5}
```

#### Proof

For a fixed posterior, entropy subadditivity gives the sum of the marginal
binary entropies.  Since `|w|>=w^2` and `h_2` increases on `[0,1/2]`,

```math
H(A\mid Z=z)
\le\sum_e h_2\left({1-|w_e(z)|\over2}\right)
\le N h_2\left({\mathcal V_z\over2N}\right).
```

Average and apply concavity once more, then subtract from `H(A)=N`. `square`

Thus a posterior barycenter with `E V_Z=o(N)` costs `N-o(N)` bits.  This
reinterprets the repository's sign-near recovery obstruction as an extremal
information lower bound rather than renaming it.

## 6. Application outside the motivating problem

Let `B in {0,1}^N` encode a graph and

```math
C_B(x)=\sum_{i<j}B_{ij}{1-x_ix_j\over2}
```

be its cut landscape.  For `M>(n-1)/2`, define the vertex-prize optimum

```math
W_B(u)=\max_x\{C_B(x)+M\langle u,x\rangle\}.
```

### Corollary 6.1 (counterfactual Max-Cut oracle is lossless)

The field pins `x=u`, so

```math
W_B(u)=Mn+C_B(u),
\qquad
\mathbb E_U[W_B(U)-Mn]U_iU_j=-{B_{ij}\over2}.        \tag{6.1}
```

Any one summary answering all such vertex-prize optima with uniform additive
error below `1/4` identifies every edge and needs `binom(n,2)` bits in the
worst case.

This is a simultaneous counterfactual-oracle theorem, not a space lower bound
for approximating one ordinary Max-Cut value.
