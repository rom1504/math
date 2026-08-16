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

## 7. Posterior width of an extremal response embedding

Let a uniform latent parameter `A` range over `{-1,+1}^N`, and let all
declared query answers be one vector

```math
R_A\in\mathcal Y,
```

where `mathcal Y` is a real Hilbert space (typically an `L^2` space of
queries).  Define

```math
\Gamma(R)=\left\{\gamma\in[0,\infty)^N:
 \sum_{e:a_e\ne b_e}\gamma_e
 \le\|R_a-R_b\|^2\quad\text{for every }a,b\right\}, \tag{7.1}
```

and

```math
\kappa(R)=\min_{a\ne b}
{\|R_a-R_b\|^2\over d_H(a,b)}.                     \tag{7.2}
```

Thus the constant vector with every coordinate `kappa(R)` belongs to
`Gamma(R)`.  Put, in bits,

```math
g(v)=h_2\left({1-\sqrt{1-v}\over2}\right),
\qquad 0\le v\le1.                                 \tag{7.3}
```

### Theorem 7.1 (posterior-width rate bound)

Let `Z` be any randomized transcript and let a decoder output
`Rhat_Z in mathcal Y`, with

```math
\Delta=\mathbb E\|R_A-\widehat R_Z\|^2.
```

For every `gamma in Gamma(R)`,

```math
I(A;Z)\ge N-
\max\left\{
 \sum_e g(v_e):
 0\le v_e\le1,
 \sum_e\gamma_ev_e\le4\Delta
\right\}.                                         \tag{7.4}
```

In particular,

```math
I(A;Z)\ge
N\left[1-g\left(
\min\left\{{4\Delta\over\kappa(R)N},1\right\}
\right)\right].                                   \tag{7.5}
```

The geometric constant has the exact posterior interpretation

```math
\inf_\pi
{\operatorname{Var}_\pi(R_A)\over
 \sum_e\operatorname{Var}_\pi(A_e)}
={\kappa(R)\over4},                                \tag{7.6}
```

where the infimum excludes a zero denominator.

#### Proof

For independent posterior draws `A,A'`, the Hilbert variance identity and
(7.1) give

```math
\operatorname{Var}_\pi(R_A)
={1\over2}\mathbb E\|R_A-R_{A'}\|^2
\ge {1\over4}\sum_e\gamma_e
       \operatorname{Var}_\pi(A_e).                \tag{7.7}
```

Taking a posterior supported equally on a pair attaining (7.2) proves
(7.6).  Conditional means minimize squared Hilbert loss, so (7.7), applied
given `Z`, implies

```math
4\Delta\ge\sum_e\gamma_e
  \mathbb E\operatorname{Var}(A_e\mid Z).           \tag{7.8}
```

For a sign of posterior mean `w`, its variance is `v=1-w^2` and its entropy
is exactly `g(v)`.  Direct differentiation shows that `g` is increasing and
concave.  Entropy subadditivity followed by Jensen therefore gives

```math
H(A\mid Z)\le
\sum_e g\bigl(\mathbb E\operatorname{Var}(A_e\mid Z)\bigr).
```

The actual variance vector is feasible in (7.4), proving it.  For constant
weights `gamma_e=kappa(R)`, a final Jensen step gives (7.5). `square`

The bound is sharp under only its geometric hypothesis.  For
`R_a=(sqrt(gamma_e)a_e/2)_e`, independent binary symmetric channels attain
equality in the weighted envelope.

### Corollary 7.2 (four exact response moduli)

The same theorem gives the following independently audited cases.

1. For the shifted quadratic family (5.2), with uniform pinning directions,

   ```math
   \kappa=4a^2.
   ```

   Orthogonality of degree-two Walsh characters gives the lower bound; a
   one-vertex switching changes `n-1` signs without changing the maximum and
   gives equality.  Consequently (7.5) strengthens Theorem 5.2 to the sharp
   curve `N[1-g(Delta/(a^2N))]`, nonzero for
   `Delta<a^2N`.
2. Let `K=2^m`, fix one anchor in the Hamming cube, and let the other `K-1`
   membership bits specify an arbitrary nonempty code.  The uniformly rooted
   nearest-code response has

   ```math
   \kappa={1\over K}.
   ```

   A changed membership root changes its distance by at least one; deleting
   one word from the full cube attains equality.  Thus any fixed per-root MSE
   below `1/4` costs `Theta(K)` bits for all sufficiently large `m`.
3. For a uniform binary `Q` by `Q` boundary kernel, endpoint-pinning queries
   return its entries and

   ```math
   \kappa={1\over Q^2},
   \qquad
   I(A;Z)\ge Q^2[1-g(\min\{4\Delta,1\})],           \tag{7.9}
   ```

   where `Delta` is uniform per-entry MSE.  The curve is sharp.  This is a
   Bayesian lossy extension of the deterministic `Q^2`-bit separator packing
   bound, not a claim about a restricted kernel semigroup.
4. For the uncentered Max-Cut landscape of a graph with independent edge
   bits, viewed under a uniform random spin query,

   ```math
   \kappa={1\over4}\qquad(n\ge3).
   ```

   Indeed, if the signed edge differences are `c_e in {-1,0,1}`, Walsh
   orthogonality gives squared response distance
   `((sum_e c_e)^2+|supp(c)|)/4`.  Dividing by the changed-edge count is at
   least `1/4`, with equality when one edge is added and another removed.
   Thus Theorem 7.1 also supplies a continuous information converse for
   counterfactual cut-value responses.  (At `n=2` the modulus is `1/2`.)

### Theorem 7.3 (orthogonal composition and its boundary)

For response maps `R,S` and positive scales `alpha,beta`,

```math
\Gamma(\alpha R\oplus\beta S)
=\alpha^2\Gamma(R)\times\beta^2\Gamma(S),          \tag{7.10}
```

and

```math
\kappa(\alpha R\oplus\beta S)
=\min\{\alpha^2\kappa(R),\beta^2\kappa(S)\}.      \tag{7.11}
```

#### Proof

Squared Hilbert distances add.  This proves sufficiency of the product
certificates.  Freeze one latent block at a time to prove necessity.  The
formula for `kappa` follows by taking a weighted average of the two child
ratios and then varying only the smaller child. `square`

Orthogonality is essential.  Same-space addition introduces the cross Gram
term `2<Delta R,Delta S>`; two one-bit children with identical `Gamma=[0,4]`
can either remain orthogonal or cancel at the parent.  Therefore `Gamma` is a
sharp information certificate for a fixed response embedding, not by itself
a reusable state under arbitrary composition.

## 8. A syndrome-rooted code feature algebra

Fix a labeled group `G=F_2^w`.  Let `H` be a full-row-rank binary parity-check
fragment with nonzero columns, and define

```math
C_H=\ker H,
\qquad
\lambda_H(s)=\min\{\operatorname{wt}(e):He=s\},
\qquad s\in G.                                     \tag{8.1}
```

For a future fragment `E` over the same interface, query

```math
\mathcal R_H(E)=\rho(\ker[H\ E]).                  \tag{8.2}
```

### Theorem 8.1 (exact syndrome response algebra)

Let `S_H` be the set of distinct nonzero column types of `H`.  Then:

```math
d(x,C_H)=\lambda_H(Hx),
\qquad
\rho(C_H)=\max_{s\in G}\lambda_H(s);               \tag{8.3}
```

```math
s\in S_H\quad\Longleftrightarrow\quad
\lambda_H(s)=1\qquad(s\ne0);                      \tag{8.4}
```

and concatenation obeys

```math
\lambda_{[H_1\ H_2]}(s)
=\min_{u\in G}
 \{\lambda_{H_1}(u)+\lambda_{H_2}(s+u)\}.          \tag{8.5}
```

Equivalently, `S_[H_1 H_2]=S_H1 union S_H2`.  For `w>=2`, `lambda_H` is,
up to injective recoding, the coarsest exact deterministic quotient for all
unrestricted future responses (8.2).

#### Proof

A correction sends `x` into `C_H` exactly when `He=Hx`, proving (8.3); full
row rank makes every syndrome occur.  A shortest correction never uses two
coordinates with the same column type, because deleting both preserves its
binary syndrome and lowers its weight.  Hence `lambda_H` depends only on
`S_H`, and its level-one set recovers `S_H`, proving (8.4).

Split a correction for `[H_1 H_2]` between the two blocks and condition on
the first syndrome.  This proves (8.5), associativity, and the support-union
law.

For operational minimality, for each nonzero `s` append the fixed fragment
`E_s` containing every nonzero column type except `s`.  It has `2^w-2`
columns and depends only on the declared apparatus.  The composite radius is

```math
\mathcal R_H(E_s)=
\begin{cases}
1,&s\in S_H,\\
2,&s\notin S_H.
\end{cases}                                        \tag{8.6}
```

In the second case write `s=u+(s+u)` with
`u notin {0,s}`.  Thus the complete future response recovers every support
bit, while (8.5) computes every future response from the profile. `square`

This quotient forgets duplicate columns and codeword geometry.  It is not
sufficient for named punctures, coordinate weights, finite-temperature
counts, or an incompatible syndrome labeling.

### Theorem 8.2 (worst-case response information)

For every `w>=2`, the deterministic message complexity of answering all
(8.2) with uniform additive error below `1/2` is `Theta(2^w)` bits.  More
precisely, it lies between

```math
2^w-1-w-\log_2(2^w-w)
```

and `2^w-1`.

#### Proof

The support indicator gives the upper bound.  Fix a basis of `G`, let
`N=2^w-1-w`, and add every `floor(N/2)`-subset of the nonbasis nonzero
types.  These full-rank fragments all have length `Theta(2^w)`.  The queries
(8.6) separate any two response vectors by one, so error below `1/2` requires
distinct messages.  Finally,

```math
\binom N{\lfloor N/2\rfloor}\ge {2^N\over N+1}.
```

Taking logarithms proves the lower bound. `square`

The statement is at unnormalized integer-lattice scale and does not imply an
exponential rate for `poly(w)`-length fragments or additive error
proportional to `w`.

### Theorem 8.3 (positive macroscopic syndrome-response rate)

Fix `L>=2`, put `w=Lq`, and decompose

```math
G=V_1\oplus\cdots\oplus V_q,
\qquad V_j\simeq\mathbb F_2^L.
```

In each block let `B_j` be a basis and `D_j=V_j\setminus{0}`.  For
`a in {0,1}^q`, let the column-type support of `H_a` use `B_j` when `a_j=0`
and `D_j` when `a_j=1`.  Repeating basis columns makes every `H_a` full rank
and of the same length `q(2^L-1)`.  For `P subset [q]`, append the valid
full-rank fragment `E_P` whose block support is `B_j` on `P` and `D_j` off
`P`.  Then

```math
\boxed{
\mathcal R_{H_a}(E_P)
=q+(L-1)|\{j\in P:a_j=0\}|.}                     \tag{8.7}
```

Consequently

```math
\max_P|\mathcal R_{H_a}(E_P)-\mathcal R_{H_b}(E_P)|
=(L-1)\max\{N_{01}(a,b),N_{10}(a,b)\}
\ge {L-1\over2}d_H(a,b).                         \tag{8.8}
```

If a deterministic summary answers every unrestricted future-fragment query
on this family with uniform error at most `eta`, then, for every `1<=d<=q`
such that

```math
2\eta<(L-1)\left\lceil{d\over2}\right\rceil,
```

its worst-case message length is at least

```math
q-\log_2\left(\sum_{i=0}^{d-1}\binom qi\right).   \tag{8.9}
```

In particular, for every fixed `epsilon<1/8`, there is
`c_epsilon>0` and an infinite sequence of widths for which error
`eta=epsilon*w` requires at least `c_epsilon*w` bits, even though all state
and query fragments have length `Theta_epsilon(w)`.  The latent vector `a`
itself is a matching `O(w)`-bit exact state for **all** future-fragment
responses on this family.  Thus its normalized response complexity is
`Theta_epsilon(w)` on this block source family.  The upper bound is not a
claim about arbitrary syndrome supports; it answers every unrestricted
future query only after the source has been restricted to the displayed
family.

#### Proof

Word length and covering radius add over the direct summands.  A basis block
has radius `L`, while a block containing every nonzero vector has radius one.
In `H_a` concatenated with `E_P`, a block remains a basis block exactly when
`j in P` and `a_j=0`, proving (8.7).  Maximizing the difference over `P`
selects either directed set difference between the zero coordinates of `a`
and `b`, proving (8.8).

A greedy binary Hamming packing of minimum distance `d` has size at least

```math
{2^q\over\sum_{i=0}^{d-1}\binom qi}.
```

Equation (8.8) and the decoder triangle inequality make distinct packing
points require distinct messages, proving (8.9).  For the asymptotic claim,
choose constants `delta<1/2` and `L` with

```math
\epsilon<{\delta(L-1)\over4L},
```

take `d=ceil(delta*q)`, and use the binary entropy estimate for the Hamming
ball.  Such constants exist for every `epsilon<1/8`; for example first take
`L>1/(1-8epsilon)` and then
`4epsilon*L/(L-1)<delta<1/2`.  Finally, `a` reconstructs the support of
`H_a`, which computes every response by (8.5). `square`

This theorem is macroscopic but deliberately restricted.  It proves only a
linear information rate on a block family; it neither proves
`Theta(2^w)` normalized complexity for arbitrary supports nor rules out a
subexponential approximate quotient of the full syndrome-support space.  The
strict forgetting property of the syndrome-support state itself is supplied
separately by Theorem 8.1 and the nonisometric duplicate-column examples.
The response formula and metric were exhaustively checked in five small cases
by
[`verify_phase2_normalized_code_rate_distortion.py`](experiments/verify_phase2_normalized_code_rate_distortion.py).

## 9. Finite deterministic synchronization

Let `Omega` be finite, let `E=binom(Omega,2)`, and let species overlaps
`R_s:E->[0,1]` have positive weights `lambda_s` summing to one.  Put

```math
q(e)=\sum_s\lambda_sR_s(e)                          \tag{9.1}
```

and define the cancellation defect

```math
\mathfrak c(e,f)=
\sum_s\lambda_s|R_s(e)-R_s(f)|-|q(e)-q(f)|.        \tag{9.2}
```

It is twice the smaller of the total upward and downward species movements.

### Proposition 9.1 (uniform cancellation implies synchronization)

If `sup_(e,f) mathfrak c(e,f)<=delta`, then for every species there is a
nondecreasing `1/lambda_s`-Lipschitz function `L_s` such that

```math
\max_e|R_s(e)-L_s(q(e))|
\le {\delta\over2\lambda_s}.                       \tag{9.3}
```

#### Proof

If `q(e)-q(f)=P-N`, where `P,N` are the weighted positive and negative
species movements, then `mathfrak c=2 min(P,N)`.  Hence

```math
R_s(e)-R_s(f)
\le{(q(e)-q(f))_++\delta/2\over\lambda_s}.         \tag{9.4}
```

The isotonic envelope

```math
L_s(p)=\inf_f\left\{R_s(f)+{(p-q(f))_+\over\lambda_s}\right\} \tag{9.5}
```

is nondecreasing and `1/lambda_s`-Lipschitz.  Equation (9.4), together with
the term `f=e`, sandwiches `L_s(q(e))` within the error in (9.3).  Clipping
to `[0,1]` does not increase it. `square`

A kernel `K` is `eta`-ultrametric when

```math
K(y,z)\ge\min\{K(x,y),K(x,z)\}-\eta                \tag{9.6}
```

for every triple of distinct states.  Two pair labels are adjacent if they
share a state.  Call `q` `(D,tau)`-monotone-linked if, whenever
`q(e)<=q(f)`, an adjacent path of length at most `D` joins them and its total
downward `q`-variation is at most `tau`.

### Theorem 9.2 (ultrametricity plus cross-root linkage)

Assume every `R_s` and every `R_s+R_t` for distinct species is
`eta`-ultrametric, and assume `q` is `(D,tau)`-monotone-linked.  Then the
functions in Proposition 9.1 can be chosen so that

```math
\boxed{
\max_e|R_s(e)-L_s(q(e))|
\le {\tau+3D\eta\over\lambda_s}.}                  \tag{9.7}
```

#### Proof

On two adjacent pair labels, opposite movements larger than `3eta` in two
species contradict ultrametricity of their sum: the third triangle edge is
within one `eta` of each lower endpoint, leaving the sum short by another
`eta`.  It follows that every adjacent step has cancellation defect at most
`6eta`.

Along a path, the triangle inequality adds these local defects.  A downward
increment is paid twice when comparing total variation with net `q` change.
Thus arbitrary endpoints have

```math
\mathfrak c(e,f)\le2\tau+6D\eta.                   \tag{9.8}
```

Apply Proposition 9.1 with this value of `delta`. `square`

### Corollary 9.3 (uniform zero-temperature control)

If a coupling potential satisfies

```math
|\Psi(r)-\Psi(r')|\le\sum_s\kappa_s|r_s-r'_s|,
```

then replacing every species profile by `(L_s(q))_s` changes

```math
\max_e\{G(e)+\Psi((R_s(e))_s)\}
```

by at most

```math
(\tau+3D\eta)\sum_s{\kappa_s\over\lambda_s}.       \tag{9.9}
```

This is a pointwise maximum estimate, not an averaged Gibbs statement.  The
cross-root linkage is substantive: PSD, weak exchangeability, ultrametricity
of every nonnegative species mixture, and even vanishing conditional
variance do not imply it; the matching counterexample is recorded in
[`examples.md`](examples.md).

### Corollary 9.4 (query-restricted exposed carrier)

Let `A subset E`, and let `H` be a subgraph of the pair-label line graph.
Assume every edge `ef` of `H` has

```math
\mathfrak c(e,f)\le\zeta,                          \tag{9.10}
```

and every ordered `e,f in A` with `q(e)<=q(f)` is joined in `H` by a path of
length at most `D` and total downward `q`-variation at most `tau`.  Then there
are nondecreasing `1/lambda_s`-Lipschitz functions `L_s` on the range `q(A)`
such that

```math
\max_{e\in A}|R_s(e)-L_s(q(e))|
\le {\tau+D\zeta/2\over\lambda_s}.                \tag{9.11}
```

Now let `G:E->R`, and let every allowed potential `Psi` have oscillation at
most `B` and satisfy

```math
|\Psi(r)-\Psi(r')|
\le\sum_s\kappa_s|r_s-r_s'|.
```

The carrier

```math
A_B(G)=\{e:G(e)\ge\max_fG(f)-B\}                  \tag{9.12}
```

contains every maximizer of `G(e)+Psi((R_s(e))_s)`.  If the preceding path
hypotheses hold for `A=A_B(G)`, replacing the species profiles by
`(L_s(q))_s` and maximizing only on that carrier changes every allowed optimum
by at most

```math
(\tau+D\zeta/2)\sum_s{\kappa_s\over\lambda_s}.    \tag{9.13}
```

#### Proof

Along an allowed path, summing the local defects and paying each downward
increment twice gives

```math
\mathfrak c(e,f)\le2\tau+D\zeta
\qquad(e,f\in A).
```

The isotonic-envelope proof of Proposition 9.1 restricted to `A` gives
(9.11).  A point outside (9.12) loses more than `B` in base score and can gain
at most `B` from the potential relative to a maximizer of `G`, so it cannot
win.  The pointwise Lipschitz bound on the common carrier then gives (9.13).
`square`

This is strictly weaker than global linkage: the rare-matching example can
have a two-label carrier satisfying (9.10)--(9.12) exactly while an arbitrarily
large unexposed fibre is disconnected.  The carrier is noncircular—it depends
only on `G` and the declared oscillation budget—but need not be succinct to
compute.  If the individual and pairwise-sum ultrametric inequalities hold
only on triangles traversed by `H`, the local no-crossing proof gives
`zeta=6eta`, recovering the error `tau+3D*eta` without global
ultrametricity.

## 10. Composition can amortize nonconvexity

Let nonempty compact component response sets `E_i subset R^p` compose by
Minkowski addition.  Put

```math
E=E_1+\cdots+E_n,
\qquad
K=\operatorname{conv}E
 =\operatorname{conv}E_1+\cdots+\operatorname{conv}E_n, \tag{10.1}
```

and let

```math
r=\dim\operatorname{span}\bigcup_i(E_i-E_i).        \tag{10.2}
```

Write the component diameters in decreasing order as
`Delta_1>=...>=Delta_n` in any fixed norm.

### Theorem 10.1 (Shapley--Folkman response bound)

```math
d_H(E,K)\le\sum_{i=1}^{\min(r,n)}\Delta_i.          \tag{10.3}
```

Consequently every `L`-Lipschitz aggregate query satisfies

```math
\left|\sup_{e\in E}\Psi(e)-\sup_{z\in K}\Psi(z)\right|
\le L\sum_{i=1}^{\min(r,n)}\Delta_i,               \tag{10.4}
```

and likewise for infima.

#### Proof

After translating one point from each component, all affine differences lie
in the `r`-dimensional space (10.2).  The Shapley--Folkman lemma represents
every `z in K` with all but at most `r` summands in the original `E_i`; only
the exceptional summands need lie in their convex hulls.  Replace each
exceptional summand by a point of its component.  The displacement is at most
the sum of their diameters, bounded by (10.3).  Since `E subset K`, this is
the Hausdorff bound.  Lipschitz continuity gives (10.4). `square`

The imported lemma originates in Starr's appendix to
[*Quasi-Equilibria in Markets with Non-Convex
Preferences*](https://doi.org/10.2307/1909201) (1969).  Equations
(10.3)--(10.4) are its response-theoretic specialization, not a claim of a
new convexity theorem.

### Corollary 10.2 (fixed-rank vector balancing)

For vectors `v_i` in a normed space, let

```math
S(V)=\left\{\sum_i\epsilon_iv_i:\epsilon_i\in\{-1,+1\}\right\},
\qquad
Z(V)=\sum_i[-v_i,v_i].                              \tag{10.5}
```

If `r=dim span{v_i}` and the norms `a_i=||v_i||` are decreasing, then

```math
d_H(S(V),Z(V))\le\sum_{i=1}^{\min(r,n)}a_i.         \tag{10.6}
```

Thus, uniformly in every target `t`,

```math
0\le
\min_{s\in S(V)}\|t-s\|-\operatorname{dist}(t,Z(V))
\le\sum_{i=1}^{\min(r,n)}a_i.                     \tag{10.7}
```

#### Proof

For `z=sum_i t_i v_i in Z(V)`, choose an extreme coefficient vector in the
fibre `{t in [-1,1]^n:sum_i t_i v_i=z}`.  At most `r` coordinates are
fractional.  Round each to its nearer sign; its displacement is at most
`||v_i||`.  This proves (10.6), and the triangle inequality gives (10.7).
`square`

For fixed `r` and bounded vector norms, the error is independent of the
number of composed components and therefore subextensive.  The state is
succinct for zonotope/generator or other controlled convex-body
representations; an arbitrary compact convex body can itself have high
description complexity.  Growing dimension is a genuine obstruction: two
vector families can have the same zonotope but an `ell_1` discrepancy gap
linear in `r`, as recorded in [`examples.md`](examples.md).

## 11. Robust tropical feature growth

The min-plus factorization rank of a finite matrix `M` is the least `k` such
that

```math
M(x,y)=\min_{t\le k}\{u_t(x)+v_t(y)\}.              \tag{11.1}
```

### Theorem 11.1 (robust tropical crossing bound)

Suppose `M` has distinguished cells `(x_i,y_i)`, `1<=i<=r`, and `G>0`
such that, for all `i!=j`,

```math
M(x_i,y_j)+M(x_j,y_i)
-M(x_i,y_i)-M(x_j,y_j)\ge G.                       \tag{11.2}
```

Every `Mtilde` with

```math
\|M-\widetilde M\|_\infty<G/4                     \tag{11.3}
```

has min-plus factorization rank at least `r`.

#### Proof

Every factor term majorizes the represented matrix, and one term is tight at
each distinguished cell.  If the same term were tight at cells `i,j`, its
separability would give

```math
\widetilde M(x_i,y_i)+\widetilde M(x_j,y_j)
\ge
\widetilde M(x_i,y_j)+\widetilde M(x_j,y_i).
```

The right side minus the left is at least
`G-4||M-Mtilde||_infinity>0`, a contradiction.  Distinct cells therefore
require distinct tight terms. `square`

The constant is sharp for this hypothesis: the matrix with diagonal zero and
off-diagonal one at order two is within `1/2=G/4` of the rank-one constant
matrix.

### Corollary 11.2 (stable code trellis channels at lattice scale)

For a binary linear code `C subset F_2^m` and a coordinate split `L sqcup R`,
put

```math
W(x_L,x_R)=d((x_L,x_R),C),
\qquad
s=\dim C-\dim C_L-\dim C_R.                        \tag{11.4}
```

Then, for every `0<=epsilon<1/2`,

```math
\min_{\|\widetilde W-W\|_\infty\le\epsilon}
\operatorname{rank}_{\min,+}(\widetilde W)=2^s.    \tag{11.5}
```

#### Proof

Sheshadri's exact theorem gives rank `2^s` and constructs a `2^s` square
submatrix with zero diagonal and all off-diagonal entries at least one.
Theorem 11.1 applies with `G=2`; taking `Mtilde=W` gives the matching upper
bound. `square`

The imported exact theorem is K. Sheshadri,
[*Trellis State Complexity as an Exact Tropical Factorization
Rank*](https://arxiv.org/abs/2607.23471) (2026).  The robust corollary is an
elementary project deduction; it was independently audited and is not stated
in that preprint.  Its scope is narrow but exact: after dividing distances by
block length `m`, the protected error is only `1/(2m)`.  It gives no
macroscopic or average-error lower bound.

## 12. Composition-created carriers and their information law

This section records the strongest theorem package from the third
investigation.  Detailed proofs and independent audits are in
[`phase3_mixed_circuit_hierarchy.md`](drafts/phase3_mixed_circuit_hierarchy.md),
[`phase3_multichannel_holonomy_packing.md`](drafts/phase3_multichannel_holonomy_packing.md),
[`phase3_carrier_capacity_law.md`](drafts/phase3_carrier_capacity_law.md), and
[`phase3_metric_quotient_synchronization.md`](drafts/phase3_metric_quotient_synchronization.md).

### Theorem 12.1 (mixed-relation exact sequence)

Let `G=W direct-sum Q`, with `W=F_2^D`, and let a labeled family of lifted
quotient columns be partitioned into fragments `E_1,...,E_m`.  Write

```math
U_j=\operatorname{span}q(E_j),
\qquad U=\sum_jU_j,
```

and let `Z` be the cycle space of all quotient columns while `Z_loc` is the
direct sum of the within-fragment cycle spaces.  Then

```math
\kappa:=\dim(Z/Z_{\rm loc})
=\sum_j\dim U_j-\dim U.                         \tag{12.1}
```

After every individual fragment gauge class is fixed, compatible global
kernel-fixing gauge classes form an affine space modeled on

```math
\operatorname{Hom}(Z/Z_{\rm loc},W).            \tag{12.2}
```

There are exactly `2^(D kappa)` labeled gluing classes.  For two fragments,
`kappa=dim(U_1 intersection U_2)`.

For a kernel coordinate basis `B`, column offsets `a_e`, their cycle
holonomy `h(R)=sum_(e in R)a_e` on `Z`, and target `(u,q_0)`, the associated
word response is exactly

```math
\ell(u,q_0)
=\min_{R:q_R=q_0}\{|R|+|u+a_R|\}.              \tag{12.3}
```

In particular, at the antipode `t=(1,...,1)`,

```math
\ell(t,0)
=D-\max_{R\in Z}\bigl(|h(R)|-|R|\bigr).        \tag{12.4}
```

#### Proof

A kernel-fixing shear adds one linear map `L(q_e)` to every offset and leaves
the sum on a quotient cycle unchanged.  Conversely, two offset assignments
with the same cycle sums differ by a map that vanishes on `ker q`, hence
factor through `im q` and differ by such a shear.  Gluing locally fixed
representatives therefore leaves exactly the homomorphisms on
`Z/Z_loc`; rank--nullity gives (12.1).  After a subset `R` of lifted columns
is chosen, its quotient must equal `q_0` and the coordinate basis supplies
the unique cheapest kernel correction, proving (12.3).  Boolean
complementation gives (12.4). `square`

The exact gauge count is not itself a response lower bound.  The next theorem
shows that a constant fraction can nevertheless remain operational.

### Theorem 12.2 (full macroscopic mixed-holonomy rate)

Let `q` be a prime power, `W=F_q^D`, `Q=F_q^k`, and use scalar-closed Hamming
alphabets.  For an independent tuple `V=(v_1,...,v_k)`, two individually
shear-trivial fragments have the unlabelled kernel-endpoint response

```math
F_V(u)=\min_{z\in F_q^k}
 \{2\operatorname{wt}(z)+\operatorname{wt}(u+Vz)\}.        \tag{12.5}
```

For all sufficiently large `D` and `1<=k<=D/32`, there is a family of at
least

```math
q^{3Dk/16}                                                   \tag{12.6}
```

such profiles whose pairwise uniform distances exceed `D/16`.  Therefore a
deterministic summary answering every endpoint to error `epsilon D`, for
fixed `epsilon<1/32`, needs at least

```math
{3\over16}Dk\log_2q                                        \tag{12.7}
```

bits.  Under a uniform prior, successful randomized reconstruction obeys the
corresponding Fano bound

```math
I(V;S)\ge(1-\eta){3\over16}Dk\log_2q-H_2(\eta).             \tag{12.8}
```

#### Proof

Let `C_V=im V`.  Dropping the coefficient toll and then representing a
nearest codeword with at most `k` basis columns gives

```math
d(u,C_V)\le F_V(u)\le d(u,C_V)+2k.             \tag{12.9}
```

A random `floor(D/4)`-dimensional linear code has minimum distance greater
than `D/8` for all large `D`.  Choose one basis for every `k`-subspace of
this host.  There are at least `q^(k(r-k))>=q^(3Dk/16)` choices.  For distinct
subspaces choose `c` in one but not the other.  The first response at `c` is
at most `2k`, while the second is greater than `D/8`, producing the stated
gap.  Packing and Fano give (12.7)--(12.8). `square`

No channel label is queried: the common query set is only `W`.  The theorem
holds over every finite field and applies equally to scalar-closed Cayley
metrics and systematic-code coset-leader profiles.

### Theorem 12.3 (presented-carrier capacity law)

Let `(X,d)` be finite.  For nonempty carriers `C_theta subseteq X` and costs
`pi_theta:C_theta->[0,p]`, define

```math
F_\theta(x)=\min_{c\in C_\theta}
 \{d(x,c)+\pi_\theta(c)\}.                    \tag{12.10}
```

Then

```math
\boxed{
\left|\|F_\theta-F_{\theta'}\|_\infty
      -d_H(C_\theta,C_{\theta'})\right|\le p.}              \tag{12.11}
```

For any probability measure `mu` on queries and `1<=s<infinity`, the same
statement holds with both functions measured in `L^s(mu)` and the carrier
metric replaced by

```math
\|d(\mathord\cdot,C_\theta)
  -d(\mathord\cdot,C_{\theta'})\|_{L^s(\mu)}.               \tag{12.12}
```

If a Hausdorff gap `Delta` is witnessed at `x_0`, then for every `r>=0`,

```math
\|F_\theta-F_{\theta'}\|_{L^s(\mu)}
\ge(\Delta-2r-p)_+\,\mu(B(x_0,r))^{1/s}.       \tag{12.13}
```

#### Proof

Write `F_theta=d_(C_theta)+e_theta`.  Nonnegative bounded presentation gives
`0<=e_theta<=p`, hence `||e_theta-e_theta'||<=p` in every displayed norm.
Distance-to-set functions embed the Hausdorff hyperspace isometrically in
the full sup norm, proving (12.11).  Their difference is two-Lipschitz, so a
gap `Delta` remains at least `Delta-2r` on the radius-`r` ball; subtract the
presentation error and integrate to obtain (12.13). `square`

Thus carrier Hausdorff packing and response packing agree up to the sharp
presentation-radius loss.  Equation (12.13) also states exactly why a
uniformly hard endpoint can disappear under a diffuse query law.

### Corollary 12.4 (a non-Hamming carrier family)

Let `E=F_(q^D)`, let `X=End_(F_q)(E)` carry rank distance, and let `M_a` be
multiplication by `a`.  For each `k`-subspace `U<=E`, choose a basis, let
`suppcost_U(a)` be the number of its nonzero basis coefficients, and put

```math
F_U(A)=\min_{a\in U}
 \{\operatorname{rank}(A-M_a)+\operatorname{suppcost}_U(a)\}.
                                                               \tag{12.14}
```

There are at least `q^(k(D-k))` such profiles, pairwise separated by at least
`D-k`.  For `k<=D/4` and `epsilon<3/8`, uniform error `epsilon D` therefore
requires at least `(3/4)Dk log_2q` bits.

#### Proof

Every nonzero multiplication map is invertible.  Distinct subspaces of the
multiplication host are consequently at full rank-metric Hausdorff distance
`D`; their basis-support presentation cost is at most `k`.  Apply Theorem
12.3 and the Gaussian-binomial count. `square`

This is a genuine rank-metric Cayley realization, although the lower proof
uses only its equilateral multiplication host rather than finer rank-ball
geometry.

### Theorem 12.5 (metric-quotient synchronization)

Let `varpi:(X,d_X)->(Y,d_Y)` be onto and one-Lipschitz.  Suppose every fibre
has diameter at most `a`, and that every quotient displacement lifts with
defect `b`: for all `x in X,y in Y`, some `z` over `y` satisfies

```math
d_X(x,z)\le d_Y(\varpi x,y)+b.                 \tag{12.15}
```

For every presented carrier `(C,alpha)` with `0<=alpha<=p`,

```math
\boxed{
0\le F_{C,\alpha}(x)-d_Y(\varpi x,\varpi C)
\le a+b+p.}                                    \tag{12.16}
```

Hence an `eta`-net for the projected carriers in Hausdorff distance gives an
`a+b+p+eta` response net of the same size.  Any fixed min-plus continuation
operator is nonexpansive in this uniform error.

If `X,Y` are abelian groups, `varpi` is a homomorphism, and carriers compose
by Minkowski sum, then projected carriers close exactly:

```math
\varpi(C+D)=\varpi C+\varpi D.                 \tag{12.17}
```

For translation-invariant metrics, infimal convolution of presentation costs
also gives the exact response update.  The maintained state is the projected
carrier plus a scalar presentation-radius certificate.

#### Proof

One-Lipschitzness gives the lower bound in (12.16).  Choose a nearest point
of `varpi C`, lift it within defect `b`, move within its fibre to the chosen
point of `C` at cost at most `a`, and then pay at most `p`.  Hausdorff nets
transfer through distance functions.  Min-plus nonexpansiveness follows by
preserving the order inequalities `g-epsilon<=f<=g+epsilon`; (12.17) is the
homomorphism identity. `square`

The theorem is a strict composable quotient only when four separately
verified conditions hold: `a+b+p` is subscale, projected carrier entropy is
smaller, the declared composition descends, and the presentation certificate
remains controlled.  Two-scale finite-field metrics give constant-state
macroscopic collapse, while rank-row projection and Hamming puncturing give
non-Hamming and coding small-error factors respectively.

## 13. A scale-rank response sandwich

For a translation-invariant metric on an `F_q`-space `W`, define

```math
s_W(\Delta)=\max\left\{
 \dim C:C\le W,
 \min_{c\in C\setminus\{0\}}\|c\|>\Delta
 \right\}.                                      \tag{13.1}
```

### Theorem 13.1 (separated rank versus synchronization rank)

For multichannel profiles

```math
F_V(u)=\min_{z\in F_q^k}
 \{2\operatorname{wt}(z)+\|u+Vz\|\},           \tag{13.2}
```

the following hold.

1. If `s=s_W(Delta)>=k` and `Delta>2k`, there are at least

   ```math
   q^{k(s-k)}
   ```

   distinct profiles separated by more than `Delta-2k`.  Uniform error
   `epsilon` with `2epsilon<Delta-2k` therefore costs at least
   `k(s-k)log_2q` bits.
2. If a linear `(a,b)` metric synchronization maps `W` onto an
   `r`-dimensional space, every profile is decoded from `varpi(im V)` to
   error `a+b+2k`, using at most

   ```math
   \sum_{j=0}^{\min\{r,k\}}{r\brack j}_q         \tag{13.3}
   ```

   states.
3. Every linear map to dimension `r` whose fibres have diameter at most `a`
   obeys the generalized Singleton inequality

   ```math
   \boxed{s_W(a)\le r.}                          \tag{13.4}
   ```

#### Proof

For (1), choose a dimension-`s` separated host and all its `k`-subspaces;
Theorem 12.3 transfers their Hausdorff separation through the `2k`
presentation toll.  For (2), apply Theorem 12.5 and count all possible
projected subspace dimensions.  For (3), a separated subspace intersects the
kernel trivially: a nonzero kernel member would be two points of one fibre at
distance greater than `a`.  It must therefore inject into the target.
`square`

Puncturing `d-1` Hamming coordinates turns (13.4) into the classical
Singleton bound.  Retaining rows of a matrix turns it into the rank-metric
Singleton bound.  On the two-scale carrier of Example 13, the complete curve
is

```math
s_W(\Delta)=
\begin{cases}
D,&0\le\Delta<1,\\
r,&1\le\Delta<L+1,\\
0,&\Delta\ge L+1,
\end{cases}                                      \tag{13.5}
```

so the lower and upper response exponents agree at every macroscopic scale
away from the presentation toll.

### Corollary 13.2 (intrinsic rank-metric response rate)

Let `E=F_(q^D)` and view `End_(F_q)(E)` as `D x D` matrices in rank metric.
For `1<=r<=D`, the Gabidulin space

```math
\mathcal G_r=
\left\{x\mapsto\sum_{i=0}^{r-1}a_ix^{q^i}:a_i\in E\right\}  \tag{13.6}
```

has dimension `rD` and minimum nonzero rank `D-r+1`.  Consequently, whenever
`2k<D-r+1`, its `k`-subspaces give at least

```math
q^{k(rD-k)}                                      \tag{13.7}
```

profiles separated by at least `D-r+1-2k`.

For `r=floor(D/2)` and `k<=D/16`, all sufficiently large `D` admit at least

```math
q^{kD^2/3}
```

profiles separated by more than `3D/8`.  Thus uniform error `epsilon D`, for
fixed `epsilon<3/16`, requires at least

```math
{1\over3}kD^2\log_2q                            \tag{13.8}
```

bits.

#### Proof

A nonzero linearized polynomial of `q`-degree at most `r-1` has at most
`q^(r-1)` roots, so its kernel dimension is at most `r-1` and its rank at
least `D-r+1`.  Its coefficient representation is unique, giving dimension
`rD`.  Apply Theorem 13.1.  With `r=floor(D/2)` and `k<=D/16`, the response
gap exceeds `3D/8` and `rD-k>=D^2/3` for large `D`. `square`

Unlike Corollary 12.4, this construction uses an MRD host and the sharp
rank-metric Singleton geometry.  Its information exponent
`Theta(D^2k log q)` is of the same order as the complete mixed-holonomy map.

### Theorem 13.3 (optimal quotient rank and code--anticode duality)

For `a>=0`, let

```math
A_W(a)=\max\{\dim K:K\le W,\ \operatorname{diam}K\le a\}.  \tag{13.9}
```

If `N=dim W`, the least dimension of a linear `(a,0)` metric
synchronization quotient of `W`, allowing the target to carry its induced
translation-invariant metric, is exactly

```math
\boxed{N-A_W(a).}                              \tag{13.10}
```

In particular,

```math
\boxed{s_W(a)+A_W(a)\le N.}                    \tag{13.11}
```

#### Proof

The kernel of any such quotient is an anticode of dimension `N-r`, giving
the lower bound on target dimension.  Conversely choose a maximizing
anticode `K`, pass to `W/K`, and use the quotient metric

```math
d(x+K,y+K)=\min_{k\in K}\|x-y+k\|.
```

The quotient map is one-Lipschitz, its fibre diameter is at most `a`, and
every quotient distance lifts with zero defect.  This proves (13.10), while
Theorem 13.1 gives (13.11). `square`

The code--anticode inequality is tight in the principal examples:

```math
\begin{array}{c|c|c}
\text{metric}&A_W(a)&s_W(a)\\ \hline
\text{Hamming on }F_q^D&\lfloor a\rfloor&\le D-\lfloor a\rfloor\\
\text{two-scale, }1\le a<L+1&D-r&r\\
\text{rank on }M_D(F_q),\ a\in\{0,\ldots,D\}&Da&D(D-a).
\end{array}                                      \tag{13.12}
```

For Hamming space, row-echelon pivots prove every `d`-dimensional subspace
has a word of weight at least `d`; equality in the separated-rank column
requires suitable MDS parameters and is not automatic.  For rank metric,
row-supported anticodes give `A_W(a)>=Da` and the Gabidulin host gives
`s_W(a)>=D(D-a)`; (13.11) forces equality.  The mismatch between lower and
upper response certificates is the scale-dependent code--anticode gap

```math
\gamma_W(a)=N-A_W(a)-s_W(a).                  \tag{13.13}
```

### Theorem 13.4 (binary Hamming has a leading duality gap)

For `W_D=F_2^D`, fixed `0<delta<1`, and
`a_D=floor(delta D)`,

```math
\liminf_{D\to\infty}{\gamma_{W_D}(a_D)\over D}
\ge H_2(\delta/2)-\delta>0.                    \tag{13.14}
```

#### Proof

The linear anticode dimension is `a_D`.  A code of minimum distance greater
than `a_D` has disjoint Hamming balls of radius `floor(a_D/2)`, so the sphere-
packing bound gives

```math
s_W(a_D)
\le D-\log_2\sum_{j\le\lfloor a_D/2\rfloor}{D\choose j}.
```

Divide by `D` and use the Hamming-ball entropy asymptotic.  Strict positivity
is `H_2(x)>2x` for `0<x<1/2`. `square`

Thus separated-host rank and optimal synchronization rank are not
asymptotically dual in general.  This does not yet determine the full
Grassmannian response entropy, because a large family of separated carrier
subspaces need not lie inside one common separated host.

## 14. Hamming Grassmannian response geometry

Let `W` be a `D`-dimensional `F_q`-space with a translation-invariant metric,
fix `C_0 in Gr_k(W)`, and give `W/C_0` its quotient norm.  Define

```math
L_{C_0}(\Delta)=\{\bar x:\|\bar x\|_{C_0}\le\Delta\},
```

```math
\Lambda_{C_0}(\ell,\Delta)=
\#\{U\le W/C_0:\dim U=\ell, U\subseteq L_{C_0}(\Delta)\}.   \tag{14.1}
```

### Theorem 14.1 (exact sparse-flat ball identity)

For directed Hausdorff distance,

```math
\#\{C\in\operatorname{Gr}_k(W):h^\to(C,C_0)\le\Delta\}
=\sum_{\ell=0}^{\min\{k,D-k\}}
 {k\brack\ell}_q q^{\ell^2}
 \Lambda_{C_0}(\ell,\Delta).                              \tag{14.2}
```

Moreover

```math
h^\to(C,C_0)=\rho_{C+C_0}(C_0),                            \tag{14.3}
```

so Hamming Grassmannian distance is the maximum of two relative covering
radii, rather than the usual injection distance.

#### Proof

The directed condition is exactly that the quotient image of `C` lies in
`L_(C_0)(Delta)`.  Put `ell=dim pi(C)`.  Choose
`J=C cap C_0` in `{k bracket ell}_q` ways and the quotient image in
`Lambda_(C_0)(ell,Delta)` ways.  Relative to a section, `C/J` is the graph of
a map from an `ell`-space to the `ell`-space `C_0/J`, giving `q^(ell^2)`
lifts.  This parametrization is unique.  For (14.3), write every member of
`C+C_0` as `c+c_0` and use translation invariance. `square`

Writing `lambda=log_q min(q^(D-k),V_W(Delta))`, ordered-basis counting in
the quotient leader ball and greedy deletion give

```math
\log_q\operatorname{Pack}(\operatorname{Gr}_k(W),d_H,\Delta)
\ge\log_q{D\brack k}_q
-\max_{0\le\ell\le\min\{k,D-k\}}
 \ell(k-\ell+\lambda)-O(\log D).              \tag{14.4}
```

For binary Hamming space, `k/D->kappa<=1/2` and
`Delta/D->delta<1/2`, this becomes

```math
\liminf {\log_2\operatorname{Pack}\over D^2}
\ge\kappa(1-\kappa)
-\max_{0\le\eta\le\kappa}
 \eta(\kappa-\eta+\min\{1-\kappa,H_2(\delta)\}).            \tag{14.5}
```

The proof is (14.2),
`Lambda(ell,Delta)<=|L(Delta)|^ell/|GL(ell,q)|`, and the Hamming-ball
entropy estimate.

### Theorem 14.2 (quotient leader geometry is not two-sided sufficient)

There are binary Hamming carriers with isometric quotient normed spaces at
every size but a linear gap in a Hausdorff query.  In one four-coordinate
block let

```math
C^{(2)}=\operatorname{span}(1100),\qquad
C^{(1)}=\operatorname{span}(1000).
```

Both quotients are the three-dimensional Hamming cube.  Nevertheless, for
`r` direct-sum blocks,

```math
d_H((C^{(2)})^{\oplus r},\{0\})=2r,
\qquad d_H((C^{(1)})^{\oplus r},\{0\})=r.       \tag{14.6}
```

Hence even the complete unlabelled quotient norm, and therefore every
sparse-flat spectrum (14.1), forgets macroscopic rooted lift information.

#### Proof

Quotienting the first two coordinates by `span(11)` leaves one parity-class
coordinate of leader cost one; quotienting by `span(10)` deletes one
coordinate.  Direct sums therefore give isometric `3r`-dimensional Hamming
quotients.  The largest kernel word has weight `2r` in the first family and
`r` in the second, proving (14.6). `square`

### Theorem 14.3 (line carriers recover unrestricted coding rate)

Let

```math
P_{D,k}(t)=\operatorname{Pack}
(\operatorname{Gr}_k(F_2^D),d_H,t)
```

with pairwise distance strictly greater than the integer `t`, and let
`A_2(D,d)` be the unrestricted binary coding number.  Then

```math
\boxed{A_2(D,t+1)-1\le P_{D,1}(t)\le A_2(D,t+1).}            \tag{14.7}
```

Consequently, for `t=floor(delta D)` and `0<delta<1/2`,

```math
1-H_2(\delta)
\le\liminf {\log_2P_{D,1}(t)\over D}
\le\limsup {\log_2P_{D,1}(t)\over D}
\le1-H_2(\delta/2).                                           \tag{14.8}
```

The puncturing quotient has `2^(D-t)` states, so its count exceeds the
maximum same-scale packing count by at least

```math
(H_2(\delta/2)-\delta-o(1))D.                                \tag{14.9}
```

More operationally, if `N_(D,1)(t)` is the radius-`t` covering number of all
binary lines, then

```math
N_{D,1}(t)
\le{(2^D-1)(\log(2^D-1)+1)\over
       \sum_{i=0}^t{D\choose i}}+1.                          \tag{14.9a}
```

Hence at `t=floor(delta D)`, puncturing uses at least

```math
(H_2(\delta)-\delta-o(1))D                                  \tag{14.9b}
```

more bits than an existing same-scale metric summary.

#### Proof

For `L_v=span(v)`, put `a=wt(v)`, `b=wt(w)`, and `c=wt(v+w)`.  Directly,

```math
d_H(L_v,L_w)=\max\{\min(a,c),\min(b,c)\}.                    \tag{14.10}
```

Thus separation greater than `t` means `c>t` and at least one of `a,b` is
greater than `t`.  Translate an optimal binary code to contain zero and
discard zero for the lower bound.  Conversely, representatives of a line
packing form a distance-`t+1` code and contain at most one word of weight at
most `t`; replace that word by zero if it exists.  This proves (14.7).
Gilbert and Hamming ball bounds prove (14.8), and puncturing `t` coordinates
proves the packing comparison (14.9).  Every radius-`t` line ball contains
at least the Hamming volume `V_D(t)`: for a heavy center it is one Hamming
ball, while for a light center it is a union of two.  Independently choose
each line center with probability `(log(2^D-1)+1)/V_D(t)` and then add every
uncovered line.  The expected size proves (14.9a); Hamming-ball entropy gives
(14.9b). `square`

### Theorem 14.4 (systematic-chart recoupling bound)

For every `1<=k<D`,

```math
P_{D,k}(t)
\le {D\choose k}A_{2^k}(D-k,t+1).                            \tag{14.11}
```

In particular, for `t<=D-k`,

```math
\log_2P_{D,k}(t)
\le k(D-k-t)+\log_2{D\choose k}.                             \tag{14.12}
```

#### Proof

Partition carriers by a chosen coordinate information set.  In one chart a
carrier has systematic generator `[I_k|X]`.  If matrices `X,Y` differ in `s`
columns, matching the same input row `u` gives

```math
d_H(C_X,C_Y)\le\max_u wt(u(X-Y))\le s.                       \tag{14.13}
```

The chart is therefore a `2^k`-ary code of distance at least `t+1`, proving
(14.11); Singleton gives (14.12).  When `2^k>=D-k`, Reed--Solomon evaluation
codes attain this column-code Singleton bound.  They need not attain the
Grassmannian bound because (14.13) discards same-input recoupling. `square`

### Theorem 14.5 (the injection/sum-weight route cannot beat a common host)

Let `C,C' in Gr_k(F_2^D)` have injection distance
`r=k-dim(C cap C')`.  If `d_H(C,C')<=t`, then

```math
2^r\le |(C+C')\cap B_D(t)|
\le\sum_{i=0}^{\min\{t,k+r\}}{k+r\choose i}.                \tag{14.14}
```

Hence, if the last sum is below `2^r`, an injection-distance-`r`
Grassmannian code is a Hausdorff-distance-`t` packing.  Greedy injection
packing gives one of size at least

```math
{2^{k(D-k)-(r-1)(D-r+1)}\over16r}.                           \tag{14.15}
```

However, put `k/D->kappa`, `r/D->rho`, and `t/D->delta`, with
`0<rho<=kappa<=1/2`.  Whenever this construction applies through

```math
(\kappa+\rho)
H_2\!\left({\delta\over\kappa+\rho}\right)<\rho,             \tag{14.16}
```

its exponent obeys

```math
\boxed{
\kappa(1-\kappa)-\rho(1-\rho)
\le\kappa(1-H_2(\delta)-\kappa).}                            \tag{14.17}
```

The right side is the elementary common-host Gilbert exponent.  Thus every
route that uses only a low-word count in `C+C'`, followed by injection
packing, is asymptotically unable to improve the host construction.

#### Proof

Hausdorff closeness supplies one distinct weight-`<=t` representative for
each of the `2^r` cosets of `C'` met by `C`.  An information set injects the
`(k+r)`-dimensional sum code into `F_2^(k+r)` without increasing weight,
proving (14.14).  The injection ball has size

```math
\sum_{s<r}2^{s^2}{k\brack s}_2{D-k\brack s}_2,
```

which is at most `16r 2^((r-1)(D-r+1))`, proving (14.15).

For (14.17), set `a=kappa+rho`, `p=delta/a`, and `h=H_2(p)`.  Condition
(14.16) is `ah<rho`.  The function

```math
x(1-x)-(a-x)H_2(ap)
```

is increasing for `x<=1/2`; at `x=ah` it equals

```math
a\{h(1-ah)-(1-h)H_2(ap)\}\ge0.
```

The last inequality follows because the braced expression decreases in
`a<=1` and vanishes at `a=1`.  Evaluation at `x=rho` gives
`rho(1-rho)>=kappa H_2(delta)`, equivalent to (14.17). `square`

Theorems 14.1--14.5 isolate the new middle layer.  One-sided leader-flat
counts are exact but not two-sided sufficient; at one channel the true
entropy is the nonlinear coding rate; at linearly many channels ordinary
column coding saturates the coarse quotient bound and cannot expose the
remaining recoupling loss.

## 15. A joint-channel response algebra

Let `(C_a,pi_a)`, `a in[q]`, be presented carriers in a finite metric space,
with `0<=pi_a<=p`, and let

```math
f_a(x)=\min_{c\in C_a}\{d(x,c)+\pi_a(c)\},
\qquad
r(a,b)=\sup_x\{f_a(x)-f_b(x)\}.                              \tag{15.1}
```

On Cartesian powers use the `ell_1` metric, product carriers, and additive
presentations.

### Theorem 15.1 (exact directed-response composition)

For words `a,b in[q]^m`, the product responses satisfy

```math
\boxed{
\|F_{\boldsymbol a}-F_{\boldsymbol b}\|_\infty
=\max\left\{
 \sum_i r(a_i,b_i),
 \sum_i r(b_i,a_i)
 \right\}.}                                                    \tag{15.2}
```

If every distinct local pair has both directed carrier distances at least
`d_0`, then

```math
r(a,b),r(b,a)\ge d_0-p.                                      \tag{15.3}
```

Consequently an outer `q`-ary code of relative distance `rho` produces
responses separated by at least

```math
(d_0-p)\rho m.                                                \tag{15.4}
```

For every `rho<1-1/q`, the family can have

```math
2^{(1-H_q(\rho)-o(1))m\log_2q}                               \tag{15.5}
```

members.  Thus `d_0>p` is a sufficient local condition for positive linear
response rate at arbitrary composition depth.

#### Proof

The product minimization factors, so
`F_boldsymbol_a(x)=sum_i f_(a_i)(x_i)`.  Suprema over the independent query
coordinates then add; applying this in both signs proves (15.2).  Since

```math
f_a(x)\ge d(x,C_a),\qquad f_b(x)\le d(x,C_b)+p,
```

and

```math
\sup_x\{d(x,C_a)-d(x,C_b)\}=h^\to(C_b,C_a),
```

both oriented inequalities (15.3) follow.  Equations (15.4)--(15.5) are the
`q`-ary Gilbert packing bound. `square`

The `q x q` directed table `r`, rather than the full local response
functions, is therefore an exact feature algebra for pairwise uniform
response distance under product composition.  The absolute value is taken
only after all channels of one orientation have added.  A coarse global
carrier estimate would give `d_0 rho m-pm`; (15.4) is larger by
`p(1-rho)m` because matching channels pay no presentation toll.

### Corollary 15.2 (two nontrivial model validations)

1. The seven nonzero words of the binary `[7,3,4]` simplex code, viewed as
   Hamming line carriers with nonzero access cost two, have `d_0=4,p=2`.
   A seven-letter outer code of relative distance `3/4` gives
   `2^((0.0573549...-o(1))m)` mixed-channel responses on `F_2^(7m)` separated
   by at least `3m/2`.
2. The seven nonzero multiplication maps of `F_8`, viewed as binary
   rank-metric line carriers with cost two, have `d_0=3,p=2`.  On
   block-diagonal products the same outer code gives the same number of
   responses separated by at least `3m/4`.

The first carrier family has bounded-weight words in every block, so it is
not contained in any common host of growing minimum distance.  The response
information is created by two-sided inter-carrier exposure, not by internal
separation in one code.

## 16. Benchmark laws: interface entropy and interacting continuation

This section records the theorem-level output of the benchmark campaign.
The standard conditional-profile and transfer-matrix identities are retained
in the benchmark drafts rather than restated here.  The results below are the
parts not supplied by the earlier universal-kernel theorem.

### Theorem 16.1 (pure-Max-Cut projective lookup and sharp response entropy)

Let `B=[w]` be an ordered boundary and

```math
X_w=\{+1,-1\}^w/\{\sigma\sim-\sigma\},
\qquad q=2^{w-1}.
```

For a nonnegatively weighted graph `G` with private vertices, define its
conditional cut profile

```math
h_G([\sigma])=max_z\operatorname {Cut}_G(\sigma,z). \tag{16.1}
```

Then the following hold.

1. If future contexts are all nonnegatively weighted Max-Cut attachments
   through `B`, their contextual metric is exactly

   ```math
   d_{\rm ctx}(G,G')=\|h_G-h_{G'}\|_\infty.        \tag{16.2}
   ```

   Hence `h_G` is the coarsest exact contextual state.
2. For every `F:X_w->[0,W]`, there is a pure weighted Max-Cut component of
   treewidth at most `w+1` with

   ```math
   h_G(s)=F(s)+(6w-2)\sum_{t\in X_w}F(t).           \tag{16.3}
   ```

   One additional private edge can pad every `F` in the cube to the common
   offset `C_(w,W)=(6w-2)qW`.
3. Consequently the translated cube of pure-Max-Cut responses has

   ```math
   \log\operatorname {Cov}_\epsilon,quad
   \log\operatorname {Pack}_\epsilon
   =\Theta\left(q\log {W\over\epsilon}\right)      \tag{16.4}
   ```

   up to universal changes of radius, for `0<epsilon<=W/6`.

#### Proof

Gluing is conditionally independent given the boundary, so every future
answer is `max_s(h_G(s)+h_C(s))`; this gives the upper bound in (16.2).  To
expose a target projective assignment `tau`, introduce a private anchor.  For
each coordinate prescribed opposite to the anchor use a direct unit cut
edge; for each coordinate prescribed equal use a two-edge unit path.  After
maximizing private spins the attachment profile is

```math
c_\tau-d_{\rm proj}(s,\tau).
```

Enough copies force the same target class for two fixed finite profiles, so
their continued-value difference exposes any chosen coordinate.  This proves
(16.2).

For the realization statement, symmetrically lift
`lambda_a=F([a])` to oriented words.  Introduce an anchor `z`, one private
spin `t_a` per oriented word, and put

```math
y_a={1+t_az\over2},
\qquad x_i=s_i z,
```

```math
E=\sum_a\lambda_a y_a
 \left(\sum_i a_ix_i-(w-1)\right).                 \tag{16.5}
```

For fixed `s,z`, exactly the term `a=(s_i z)_i` has positive activation,
equal to one; all others are nonpositive.  Maximizing the `t_a` therefore
returns `F([s])`.  Since `z^2=1`, (16.5) is pairwise:

```math
E=-{w-1\over2}\sum_a\lambda_a
 +{1\over2}\sum_{a,i}\lambda_aa_i(s_iz+t_as_i)
 -{w-1\over2}\sum_a\lambda_at_az.                 \tag{16.6}
```

A signed pair term `Juv` is implemented, up to a constant, by positive cut
edges.  For `J<=0`, a direct edge of weight `-2J` scores `-J+Juv`.  For
`J>=0`, a fresh two-edge path of weights `2J` scores `3J+Juv` after its
middle spin is maximized.  Keep the occurrences in (16.6) separate.  For the
oriented pair `a,-a`, the two boundary-anchor groups add `2w lambda_a`, the
two selector-boundary groups add another `2w lambda_a`, the selector-anchor
terms add `(w-1)lambda_a`, and removing the negative constant adds
`(w-1)lambda_a`.  The total offset per projective class is
`(6w-2)lambda_a`, proving (16.3).

The bags `B union {z}`, `B union {z,t_a}`, and three-vertex mediator bags
give width at most `w+1`.  A private edge supplies the nonnegative padding to
`C_(w,W)`.  Finally (16.2) reduces deterministic response covering and
packing to the ordinary sup cube, proving (16.4). `square`

The construction is exponential in `w` and has boundary load
`4 sum_tF(t)` at every boundary vertex.  Thus (16.4) is an exact
unrestricted-size response theorem, not a polynomial-size or
unit-boundary-sensitivity lower bound.

### Theorem 16.2 (metric-interface Lipschitz response entropy)

Let `(X,d)` be a nonempty finite metric space of diameter `D`, and let
`Lip_1(X)/R` be the one-Lipschitz functions modulo constants, with

```math
d_{\rm sh}([f],[g])=inf_c\|f-g-c\boldsymbol1\|_\infty
={1\over2}\operatorname {osc}(f-g).                \tag{16.7}
```

If `S` is an `r`-net of `X`, `eta>0`, and `r+eta/2<=delta`, then

```math
\operatorname {Cov}_\delta(Lip_1(X)/R)
\le\left({D\over\eta}+2\right)^{|S|}.              \tag{16.8}
```

Conversely, if `C subset X` is non-strictly `rho`-separated and `k=|C|`,
then, for every `0<rho'<rho`,

```math
d(c,c')\ge\rho\quad(c\ne c'),                      \tag{16.9}
```

then

```math
\operatorname {Pack}_{\rho'}(Lip_1(X)/R)
\ge {k\choose\lfloor k/2\rfloor}.                 \tag{16.10}
```

For the `w`-cube, write `V(w,r)=sum_(j<=r)binom(w,j)`.  One may take

```math
|S|\le\left\lceil{2^w\over V(w,r)}(w\log2+1)\right\rceil, \tag{16.11}
```

and, at separation `rho`, a greedy `C` of size at least

```math
\left\lfloor{2^w\over V(w,\lceil\rho\rceil-1)}\right\rfloor. \tag{16.12}
```

Therefore, for fixed `0<epsilon<1/4`,

```math
\log_2\operatorname {Cov}_{\epsilon w}
\le2^{(1-H_2(\epsilon)+o(1))w},
\qquad
\log_2\operatorname {Pack}_{\epsilon w}
\ge2^{(1-H_2(\epsilon)+o(1))w}.                   \tag{16.13}
```

In particular the actual radius-`epsilon w` covering number obeys

```math
2^{(1-H_2(2\epsilon)+o(1))w}
\le \log_2\operatorname {Cov}_{\epsilon w}
\le2^{(1-H_2(\epsilon)+o(1))w}.                  \tag{16.14}
```

#### Proof

Normalize `min_x f(x)=0`, store upward-rounded values
`q_s=eta ceil(f(s)/eta)` on `S`, and define

```math
u(x)=\min_s\{q_s+d(x,s)\},
\quad
\ell(x)=\max_s\{q_s-d(x,s)\},
\quad g={u+\ell\over2}.                            \tag{16.14a}
```

Both envelopes and their midpoint are one-Lipschitz.  The global Lipschitz
inequalities and a landmark within `r` give

```math
f\le u<f+\eta+2r,
\qquad f-2r\le\ell<f+\eta.
```

Thus `g-f` lies in `[-r,eta+r)`, so its oscillation is less than
`eta+2r` and its shape error is less than `r+eta/2<=delta`.  The normalized
sample range is contained in `[0,D]`, proving (16.8).

For each `U subset C` of cardinality `floor(k/2)`, prescribe value `rho` on
`U` and zero on `C setminus U`.  Condition (16.9) makes these partial data
one-Lipschitz, so they extend to `X` (and may be clipped to `[0,rho]`).  If
`U ne V`, the difference is `+rho` somewhere in `U setminus V` and `-rho`
somewhere in `V setminus U`; (16.7) makes their shape distance at least
`rho`, proving (16.10).

Random centers prove (16.11), while greedy deletion of radius-
`ceil(rho)-1` balls proves (16.12).  In (16.8) take `eta=o(w)` and
`r=(epsilon-o(1))w`; in (16.12) take `rho=epsilon w+o(w)`.  The standard
Hamming-ball asymptotics and
`binom(k,floor(k/2))>=2^k/(k+1)` give (16.13).  Finally, a radius-
`epsilon w` cover can contain at most one member of a set separated by more
than `2epsilon w`; applying (16.13) at that separation gives the lower half
of (16.14). `square`

This is an upper theorem for every effectively Lipschitz boundary-response
class.  The lookup constructions for pairwise Ising and Theorem 16.1 realize
the lower packing in those languages.  They do not impose a syntactic
unit-boundary-degree promise.

### Theorem 16.3 (exposed witnesses under min-plus continuation)

Let `Y,X` be finite, `K:X times Y->R`, and

```math
(T_Kf)(x)=\min_y\{K(x,y)+f(y)\},
\qquad r(f,g)=\max_y(f(y)-g(y)).                    \tag{16.15}
```

Define the exposure penalty

```math
e_K^f(y)=\min_x\{K(x,y)+f(y)-(T_Kf)(x)\}\ge0.      \tag{16.16}
```

Then

```math
\boxed{
\max_y\{f(y)-g(y)-e_K^f(y)\}
\le r(T_Kf,T_Kg)\le r(f,g).}                       \tag{16.17}
```

In particular, an exact maximizer of `f-g` which is exposed by `f` preserves
that directed response.  Along repeated continuations, the losses

```math
r(f,g)-\max_y\{f(y)-g(y)-e_K^f(y)\}                \tag{16.18}
```

add at most telescopically.  The analogous condition in the reverse
orientation controls the two-sided response separation.

There is a genuine interacting family on which (16.17) closes exactly.  Let
`(Y,d)` be a finite metric space, let `g` be a bijective self-isometry, and
for `lambda>=0` put

```math
D_{\lambda,g}(a,t)=\lambda d(t,g(a)).              \tag{16.19}
```

With min-plus kernel composition
`(K star L)(a,t)=min_u(K(a,u)+L(u,t))`, one has

```math
\boxed{D_{\lambda,g}\star D_{\mu,h}
=D_{\min\{\lambda,\mu\},h\circ g}.}               \tag{16.20}
```

Moreover,

```math
r(D_{\lambda,g}(a,\cdot),D_{\lambda,g}(b,\cdot))
=\lambda d(a,b).                                   \tag{16.21}
```

Thus a chain closes exactly on the bottleneck strength and isometry
holonomy:

```math
D_1\star\cdots\star D_T=D_{\lambda_*,G},
\qquad \lambda_*=\min_i\lambda_i,
\quad G=g_T\circ\cdots\circ g_1.                 \tag{16.22}
```

This closure is robust in the precise cumulative sense requested for an
interacting continuation.  If the actual `i`th kernel is within `eta_i` of
`D_(lambda_i,g_i)` in entrywise sup norm, then its full composite is within
`E=sum_i eta_i` of (16.22), and every directed row response is within `2E`
of `lambda_*d(a,b)`.  Hence `E=o(lambda_*d(a,b))` preserves the leading
directed scale for that declared pair through arbitrarily many steps.

#### Proof

For the upper inequality in (16.17), evaluate `T_Kf` at a minimizer used by
`T_Kg`.  For the lower inequality, take `x` attaining (16.16).  Then

```math
(T_Kf)(x)=K(x,y)+f(y)-e_K^f(y),
\qquad
(T_Kg)(x)\le K(x,y)+g(y).
```

Subtract and maximize over `y`.  Iteration gives (16.18).

For (16.20), move `t` by `h^(-1)`.  The triangle inequality gives

```math
\lambda d(u,g(a))+\mu d(u,h^{-1}(t))
\ge\min(\lambda,\mu)d(g(a),h^{-1}(t)).
```

Choosing `u=g(a)` or `u=h^(-1)(t)` attains the smaller endpoint cost.  The
reverse triangle inequality bounds (16.21) above, while `t=g(b)` attains the
bound.  Entrywise min-plus composition is one-Lipschitz in each factor, and
the directed maximum is two-Lipschitz in the two rows, proving the robust
claim. `square`

The exposure inequality is sharp throughout this family.  When a row of
strength `lambda` is continued by a link of strength `mu>=lambda`, every
intermediate point is exposed: choose its isometric image as the output and
use the triangle inequality.  If `mu<lambda`, only the row's centre is
exposed (for every other point, jumping first to the centre is strictly
cheaper).  The maximizing witness for the ordered pair `(a,b)` is the centre
of row `b`, so the weak link hides it and clips the response from
`lambda d(a,b)` to `mu d(a,b)`.

For the discrete metric on `[q]`, define the permutation-Potts reward kernel

```math
P_{J,\pi}(a,t)=-J{\bf1}\{t=\pi(a)\}
=D_{J,\pi}(a,t)-J.                                 \tag{16.23}
```

Then

```math
P_{J,\pi}\star P_{L,\rho}
=-\max(J,L)\boldsymbol1+P_{\min(J,L),\rho\circ\pi}, \tag{16.24}
```

and a chain has endpoint kernel

```math
-\left(\sum_iJ_i-\mu\right)\boldsymbol1+P_{\mu,\Pi},
\qquad \mu=\min_iJ_i,
\quad \Pi=\pi_T\circ\cdots\circ\pi_1.           \tag{16.25}
```

Under arbitrary labelled endpoint fields, distinct isometries at fixed
`mu>0` are separated by at least `mu delta_Y`, where `delta_Y` is the least
positive distance in `Y`; for the discrete metric this forces `q!` states,
or `Theta(q log q)` discrete bits, below error `mu/2`.  The directed row
table alone forgets the holonomy and retains only the bottleneck.  Thus the
declared future language determines whether the exact state is `(mu,Pi)` or
just `mu`; at `mu=0` the holonomy is operationally invisible.

Finally the signed binary Ising kernel has the exact affine identification

```math
K_J^{\rm Ising}(s,t)=-Jst
=|J|+P_{2|J|,\pi_J}(s,t),
\qquad \pi_J(s)=\operatorname {sgn}(J)s.           \tag{16.26}
```

Consequently its chain endpoint is

```math
-\left(\sum_i|J_i|-\nu\right)
-\operatorname {sgn}\left(\prod_iJ_i\right)\nu st,
\qquad \nu=\min_i|J_i|,                            \tag{16.27}
```

and both directed row gaps are `2nu`.

### Theorem 16.4 (normalized Max-Cut distance-shell characterization)

Let

```math
X_w=\{+1,-1\}^w/\{s\sim-s\},
```

and let `lambda_i>=0`.  The weighted projective Hamming pseudometric is

```math
d_\lambda([s],[t])=
\min\left\{
\sum_{i:s_i\ne t_i}\lambda_i,\,
\sum_{i:s_i=t_i}\lambda_i
\right\}.                                        \tag{16.28}
```

(It becomes a metric after zero-cost coordinates are quotiented.)  For a
nonnegatively weighted Max-Cut component `G` with boundary `[w]`, let
`ell_i(G)` be
the total weight of edges incident to boundary vertex `i`, and let `[h_G]`
denote its conditional cut profile modulo constants.  Then

```math
\boxed{
\{[h_G]:\ell_i(G)\le\lambda_i\ \hbox{for every }i\}
=\operatorname {Lip}_1(X_w,d_\lambda)/R.}         \tag{16.29}
```

If

```math
\Delta_i(f)=\max_s|f([s])-f([s^{(i)}])|,
```

the minimum load vector of the shape `[f]` is exactly

```math
\inf_{G:[h_G]=[f]}\ell_i(G)=\Delta_i(f)            \tag{16.30}
```

simultaneously in all coordinates.  In particular the minimum total exposed
load is `sum_i Delta_i(f)`.

The same normalized language remains fully exposing anisotropically: if both
tested components and all future attachments obey `ell_i<=lambda_i`, the
literal contextual metric is `||f-g||_infinity`, while after offset
calibration the shape metric is `osc(f-g)/2`.  Consequently, for the
unit-load class
`L_w=Lip_1(X_w,d_proj)/R` and fixed `0<epsilon<1/4`,

```math
\boxed{
2^{(1-H_2(2\epsilon)+o(1))w}
\le\log_2\operatorname {Cov}_{\epsilon w}(L_w)
\le2^{(1-H_2(\epsilon)+o(1))w}.}                 \tag{16.31}
```

Thus boundary sensitivity one still permits exponentially many response
bits at macroscopic distortion.

#### Proof

Flipping a set `D` of boundary spins while holding private spins fixed changes
the cut by at most `sum_(i in D)ell_i(G)`.  Comparing conditional maxima in
both directions, and then using global-flip invariance, proves the
`d_lambda`-Lipschitz necessity.

For sufficiency, translate a `d_lambda`-Lipschitz `f` to be nonnegative and
use Theorem 16.1 to compile `C+f([y])` at a private inner interface
`y=(y_i)`.  Join each true boundary spin `s_i` to `y_i` by a fresh two-edge
path whose two edges have weight `lambda_i`.  Maximizing the middle spin
gives `2lambda_i` when `s_i=y_i` and `lambda_i` otherwise.  The outer
response, apart from a common constant, is

```math
\max_y\left\{f([y])-\sum_{i:s_i\ne y_i}\lambda_i\right\}=f([s]). \tag{16.32}
```

The equality is the max-plus McShane identity: `y=s` attains the right side,
and Lipschitzness gives the reverse inequality.  Only the first edge of each
path meets the true boundary, proving (16.29).  Coordinate flipping gives
the lower half of (16.30), while telescoping shows that `f` is
`d_Delta`-Lipschitz and the same construction attains all coordinates at
once.

A weighted positive-edge pin has profile `C-d_lambda(x,t)` and exposed load
`lambda_i` in coordinate `i`.  The McShane identity exposes `f(t)` for every
`d_lambda`-Lipschitz `f`, proving the restricted contextual isometry.  The
upper half of (16.31) is Theorem 16.2 on the
projective cube, whose balls below radius `w/2` have volume `V(w,r)`.
For the lower half, take a projective code of distance
`2h>2epsilon w` and size at least
`2^(w-1)/V(w,2h-1)`.  On its codewords assign `+h` and `-h` according to
each half-size subset.  These labels extend one-Lipschitzly; two distinct
half-size subsets give both differences `+2h` and `-2h`, hence shape
distance at least `2h`.  Equation (16.29) realizes every extension at unit
load, and the standard Hamming-volume estimate proves (16.31). `square`

The private compiler in this proof has exponential size.  Therefore
(16.31) rules out boundary load alone as a compression promise; it does not
rule out a smaller response class under polynomial graph size or bounded
description complexity.  Unit load is not closed under parallel gluing,
because loads add, whereas serial distance shells are idempotent.

Abstractly, the shell is the max-plus distance projector

```math
(P_df)(x)=\max_y\{f(y)-d(x,y)\}.                  \tag{16.33}
```

Its fixed points are exactly the one-Lipschitz functions, and the distance
kernel is idempotent:

```math
\max_y\{-d(x,y)-d(y,z)\}=-d(x,z).                \tag{16.34}
```

Thus a private universal profile compiler followed by one resource-bounded
distance bridge realizes the full Lipschitz response ball at the bridge's
resource cost.  Equations (16.20) and (16.34) are the min-plus and max-plus
forms of the same metric closure mechanism.

### Theorem 16.5 (exact tropical resource distortion)

Let `(X,d)` be a nonempty finite metric space, `lambda>=0`, and `f:X->R`.
Define

```math
(P_\lambda f)(x)=\max_y\{f(y)-\lambda d(x,y)\},
\qquad
(Q_\lambda f)(x)=\min_y\{f(y)+\lambda d(x,y)\},   \tag{16.35}
```

and the directed Lipschitz defect

```math
\Delta_\lambda(f)=
\max_{x,y}\{f(y)-f(x)-\lambda d(x,y)\}.           \tag{16.36}
```

Then `P_lambda f` is the least `lambda`-Lipschitz majorant of `f`,
`Q_lambda f` is its greatest `lambda`-Lipschitz minorant, and

```math
P_\mu P_\lambda=P_{\min(\lambda,\mu)},
\qquad
Q_\mu Q_\lambda=Q_{\min(\lambda,\mu)}.            \tag{16.37}
```

The exact projective response distance to the admissible Lipschitz class is

```math
\boxed{
\inf_{h\in\operatorname {Lip}_\lambda}
d_{\rm sh}([f],[h])
={\Delta_\lambda(f)\over2}
=d_{\rm sh}([f],[P_\lambda f]).}                  \tag{16.38}
```

The exact literal sup distance is also

```math
\boxed{
\inf_{h\in\operatorname {Lip}_\lambda}\|f-h\|_\infty
={\Delta_\lambda(f)\over2},}                      \tag{16.39}
```

attained by `(P_lambda f+Q_lambda f)/2`.

For a chain of isometrically twisted distance shells with strengths
`lambda_i`, holonomy `G`, and `lambda_*=min_i lambda_i`, the response
distortion relative to the relabelled input is exactly one weakest-layer
defect:

```math
d_{\rm sh}\left(
[T_T\cdots T_1f],[f\circ G^{-1}]
\right)
={\Delta_{\lambda_*}(f)\over2}.                  \tag{16.40}
```

It is not paid once per layer.  If the `i`th shell has entrywise sup error
at most `eta_i`, the right side of (16.40) increases by at most
`sum_i eta_i`.

#### Proof

Every distance cone in (16.35) is `lambda`-Lipschitz, so both envelopes are
Lipschitz and `P_lambda f>=f>=Q_lambda f`.  If `h>=f` is Lipschitz, then

```math
h(x)\ge h(y)-\lambda d(x,y)
\ge f(y)-\lambda d(x,y),
```

so `h>=P_lambda f`; the minorant statement follows by sign reversal.
The triangle inequality, with equality by choosing either endpoint, proves
(16.37).

Directly from the definitions,

```math
\|P_\lambda f-f\|_\infty
=\|f-Q_\lambda f\|_\infty
=\Delta_\lambda(f).                              \tag{16.41}
```

At a global maximum of `f`, the first error is zero, so
`d_sh([f],[P_lambda f])=Delta_lambda(f)/2`.  Conversely, if `x,y` attain
(16.36), every Lipschitz `h` satisfies

```math
(f-h)(y)-(f-h)(x)\ge\Delta_\lambda(f).
```

This proves the lower bounds in (16.38)--(16.39).  The midpoint of the two
Lipschitz envelopes is Lipschitz, and its error from `f` is
`((P_lambda f-f)-(f-Q_lambda f))/2`, whose absolute value is at most
`Delta_lambda(f)/2`.  This proves (16.39).

Finally Theorem 16.3 identifies the whole twisted chain with the single
shell `P_(lambda_*)` after relabelling by `G`.  Equation (16.38) proves
(16.40), and uniform kernel nonexpansiveness gives the perturbation
statement. `square`

The envelope identities are classical tropical projection.  Their
extremal-information content is the exact operational formula (16.38) and
its compiler/bridge applications: in every language satisfying the
private-compiler and sensitivity hypotheses used in Theorem 16.4, a single
distance shell is a nearest resource-admissible response, and repeated
interacting composition pays only the weakest layer.

### Theorem 16.6 (anisotropic separator bottleneck lattice)

On the projective Boolean cube `X_w`, let `d_ell` be the weighted metric
(16.28), now with a vector `ell in R_+^w`, and let

```math
B_{\ell,g}(a,t)=-d_\ell(t,g(a)),                  \tag{16.42}
```

where `g` is a signed coordinate permutation (global sign is quotiented).
If `h_*ell` denotes the coordinate vector transported by `h`, then max-plus
composition obeys

```math
\boxed{
B_{\ell,g}\odot B_{m,h}
=B_{\,m\wedge h_*\ell,\,h\circ g},}              \tag{16.43}
```

where the wedge is coordinatewise minimum.  Thus arbitrary anisotropic
shell chains close on one load vector and one monomial holonomy.  In the
aligned case,

```math
B_\ell\odot B_m=B_{\ell\wedge m}.                 \tag{16.44}
```

At matched scalar precision this uses `w` resource values plus a signed
permutation, rather than an arbitrary `2^(w-1) by 2^(w-1)` transfer kernel.
For a profile transported through such a chain, its exact shape distortion
is `Delta_(ell_*)/2`, with `ell_*` the recursively transported wedge and
`Delta` defined as in (16.36) using `d_(ell_*)`.

#### Proof

First consider aligned shells.  For oriented words and fixed endpoint
orientation, the intermediate coordinates are independent:

```math
\min_y\left\{
\sum_i\ell_i{\bf1}\{x_i\ne y_i\}
+\sum_i m_i{\bf1}\{y_i\ne z_i\}
\right\}
=\sum_i\min(\ell_i,m_i){\bf1}\{x_i\ne z_i\}.     \tag{16.45}
```

Allowing the two independent global orientations of the intermediate
projective word leaves only their product at the endpoints, so minimizing
that product gives `d_(ell wedge m)([x],[z])`.  Negating proves (16.44).

For (16.43), change variables from the intermediate word `u` to `h(u)`.
The first distance becomes `d_(h_*ell)(h(u),h(g(a)))`, while the second has
weight vector `m`.  Equation (16.44) gives their coordinatewise wedge and
the endpoint map `h circ g`.  Induction proves chain closure.  The distortion
claim is Theorem 16.5 for the final weighted projective metric. `square`

This is a strict composable quotient for a nontrivial growing-interface
family.  It does not compress an arbitrary separator profile travelling
through the shell: the profile may still have exponential response
complexity.  It compresses the interacting *continuation algebra* and shows
that heterogeneous resource constraints combine by a transported bottleneck
lattice rather than by additive loss.

### Theorem 16.7 (sharp metric recognition and a long-depth obstruction)

Let `K` be a symmetric hollow real kernel on a `q`-point set, `q>=2`, with min-plus
square `K star K`, and put

```math
\tau(K)=\|K-K\star K\|_\infty,
\qquad
c_q=\max\left\{{1\over2},{q-2\over q}\right\}.    \tag{16.46}
```

Then:

1. `K star K=K` exactly when `K` is a pseudometric.
2. There is a pseudometric `d` with

   ```math
   \boxed{\|K-d\|_\infty\le c_q\tau(K).}          \tag{16.47}
   ```

   The coefficient `c_q` is the best possible universal coefficient for
   every `q`.
3. If a finite permutation group `Gamma` acts on the alphabet and

   ```math
   \omega_\Gamma(K)=
   \max_{\gamma,x,y}|K(\gamma x,\gamma y)-K(x,y)|,
   ```

   there is a `Gamma`-invariant pseudometric `d` satisfying

   ```math
   \|K-d\|_\infty
   \le\omega_\Gamma(K)+c_q\tau(K).                \tag{16.48}
   ```

Thus approximate idempotence recognizes one nearby metric shell
dimension-freely, and full-orbit invariance recognizes a nearby composable
metric-isometry family.

This one-step theorem cannot be iterated for free.  For `a>delta>0` on
`{0,...,q-1}`, define

```math
K_\delta(i,j)=
\begin{cases}
0,&i=j,\\
a|i-j|-\delta,&i\ne j.
\end{cases}                                      \tag{16.49}
```

It has `tau(K_delta)=delta` and lies within `delta` of the path metric
`a|i-j|`, but for `1<=T<=q-1`,

```math
K_\delta^{\star T}(i,j)=
\begin{cases}
0,&i=j,\\
a|i-j|-\delta\min\{T,|i-j|\},&i\ne j.
\end{cases}                                      \tag{16.50}
```

Consequently

```math
\|K_\delta^{\star T}-K_\delta\|_\infty
=(T-1)\delta,                                    \tag{16.51}
```

and the shape distance between row zero before and after composition is
`(T-1)delta/2`.  Taking `delta=c/q` and `T=q-1` makes the local defect and
distance to an exact metric vanish while leaving a fixed long-depth
response drift.

#### Proof

The choices of an endpoint in `K star K` show `K star K<=K`, and

```math
\tau(K)=\max_{i,j,k}\{K(i,j)-K(i,k)-K(k,j)\}.     \tag{16.52}
```

Thus zero defect is precisely the triangle inequality; diagonal triples and
symmetry also force nonnegativity.

For the quantitative repair, add `c_q tau` to every off-diagonal entry of
`K` and take the shortest-path pseudometric `d` of that complete weighted
graph.  Diagonal triples in (16.52) give `K(i,j)>=-tau/2`, so the shifted
edges are nonnegative.  The direct edge gives `d(i,j)<=K(i,j)+c_q tau`.
A shortest path is simple, with `ell<=q-1` edges, and repeated relaxed
triangles give

```math
K(i,j)\le\sum_{r=1}^{\ell}K(v_{r-1},v_r)
          +(\ell-1)\tau.
```

Its shifted length is at least
`K(i,j)+c_q tau+(\ell-1)(c_q-1)tau`, which is at least
`K(i,j)-c_q tau` precisely because `q c_q>=q-2`.  This proves (16.47).
For `q=2,3`, one off-diagonal entry `-tau/2` proves sharpness.  For `q>=4`,
the kernels (16.49) force any approximating pseudometric, by the endpoint
triangle along the line, to have error at least `(q-2)tau/q`.

For (16.48), average `K` over `Gamma`.  The average moves by at most
`omega_Gamma`, remains symmetric and hollow, and has triangle defect at most
`tau`; the shifted shortest-path repair is invariant.

Finally a `T`-factor path in (16.49) with `p` nonzero moves and total integer
variation `V` costs `aV-delta p`.  Since `V>=|i-j|` and `a>delta`, the
minimum uses `min(T,|i-j|)` monotone nonzero pieces, proving (16.50).
Comparing an adjacent endpoint with one at distance at least `T` gives
(16.51) and the shape claim. `square`

The final family isolates the information omitted by one-step recognition:
a per-transition toll is collected once per useful segment.  Approximate
tropical idempotence is therefore a stable local recognition criterion but
not a depth-stable quotient.  A long-composition theorem needs an exact
retraction, bounded useful path length, or a closed metric semilattice such
as Theorem 16.6.  Reading all triangle defects still requires the full
kernel, so (16.47) is a recognition theorem, not by itself a compression
algorithm.  No claim of external novelty is made for finite metric repair.

There is one sharp positive restoration.  Suppose `K>=0`, let `K_*` be its
min-plus shortest-path closure, and assume every pair has a shortest
`K`-path using at most `H` nonzero edges.  Then for every depth `T`,

```math
0\le K-K^{\star T}\le K-K_*
\le(H-1)\tau(K)\boldsymbol1.                     \tag{16.53}
```

Indeed, repeated relaxed triangles charge at most one `tau` per internal
vertex of a shortest path, while zero diagonal steps pad any shorter path.
If all positive off-diagonal costs lie in `[m,D]`, one may take
`H<=floor(D/m)`.  Thus bounded useful hop length is a concrete sufficient
mechanism for depth-uniform approximate composition; Theorem 16.7 shows that
some such global condition is necessary.

### Theorem 16.8 (shared-parameter tropical response entropy)

Let `X` have `q` queries.  For each `x`, let `A_x subset R^d` have at most
`r` elements, and consider

```math
f_\theta(x)=\max_{a\in A_x}\langle a,\theta\rangle,
\qquad \theta\in\Theta\subseteq R^d.             \tag{16.54}
```

Suppose all response shapes `[f_theta]` lie in the radius-`R` ball of
`R^X/R1` with norm `osc/2`.  Then, without any bound on the magnitude or
precision of `theta`,

```math
\boxed{
\log_2\operatorname {Cov}_\delta\{[f_\theta]\}
\le d\log_2(4(qr^2+1))
 +d\log_2\left(1+{2R\over\delta}\right).}        \tag{16.55}
```

If every witness vector is binary, all comparison normals belong to
`{-1,0,1}^d`; hence

```math
\log_2\operatorname {Cov}_\delta
=O\left(d^2+d\log(1+R/\delta)\right),             \tag{16.56}
```

independently of the number of queries and witnesses.

For a concrete consequence, let `C_(w,p,m)` be the response shapes of
unit-boundary-load pure weighted Max-Cut components with at most `p` private
vertices and `m` nonzero nonparallel edges.  Put

```math
H_*=\min\left\{2^{w+2p-2},{3^m-1\over2}\right\}.
```

Then

```math
\begin{aligned}
\log_2\operatorname {Cov}_\delta(C_{w,p,m})
\le{}&\log_2(m+1)+2m\log_2(w+p)\\
&+m\log_2(4(H_*+1))
+m\log_2\left(1+{w\over2\delta}\right).
\end{aligned}                                    \tag{16.57}
```

After boundary-disconnected private components are removed one may take
`p<=m`; therefore, at fixed macroscopic error,

```math
\log_2\operatorname {Cov}_{\epsilon w}(C_{w,m})
=O_\epsilon(m^2+m\log(w+m)).                     \tag{16.58}
```

This is semantic response compression for arbitrary real edge weights, not
finite-encoding cardinality.

It sharply separates the finite presentation from the universal unit-load
compiler of Theorem 16.4.  If `m=m(w)`-edge profiles form an
`epsilon w`-net of the entire unit-load Lipschitz response class, then

```math
\boxed{
\liminf_{w\to\infty}{\log_2m(w)\over w}
\ge {1-H_2(2\epsilon)\over2},
\qquad 0<\epsilon<1/4.}                          \tag{16.59}
```

In particular, no polynomial-size component family approximates every
unit-load response at fixed macroscopic accuracy, even with unlimited
weight precision.

#### Proof

For each query, changes of optimizer occur only on the comparison
hyperplanes

```math
\langle a-b,\theta\rangle=0
\qquad(a,b\in A_x).
```

There are at most `q binom(r,2)` such hyperplanes.  Their arrangement has at
most `[4(qr^2+1)]^d` relatively open faces, including lower-dimensional
tie faces.  On each face a fixed tie rule chooses one optimizer per query,
so the whole response is one linear map of `theta`.  Its image modulo
constants has dimension at most `d` and lies in a radius-`R` ball.  A
maximal separated subset and the usual volume comparison in that quotient
norm give a cover of size `(1+2R/delta)^d` on the face.  Union over faces
proves (16.55).  For binary witnesses, opposite nonzero vectors in
`{-1,0,1}^d` define the same hyperplane, leaving at most `(3^d-1)/2`;
this proves (16.56).

For a fixed Max-Cut topology with `e<=m` edges, each private cut gives an
incidence vector in `{0,1}^e`, and the conditional optimum is precisely
(16.54) with `theta` equal to the edge weights.  There are fewer than
`2^(w+2p-2)` listed comparisons and at most `(3^e-1)/2` distinct
comparison hyperplanes.  Unit boundary load makes every response
one-Lipschitz on the projective cube, so its shape norm is at most `w/4`.
There are at most `(m+1)(w+p)^(2m)` padded labelled topologies.  Multiplying
the face, volume, and topology covers proves (16.57).

Parallel edges combine, zero edges disappear, and every remaining connected
component with private vertices and a boundary vertex has at least as many
edges as private vertices.  Components disjoint from the boundary contribute
only a quotient constant.  Hence `p<=m` and (16.58) follows.

Finally let `C_(w,m)` be an `epsilon w`-net of the full unit-load class.  A
`tau w`-cover of `C_(w,m)` is an `(epsilon+tau)w`-cover of that class.
The balanced projective-code packing in Theorem 16.4 and (16.58) imply

```math
m^2+m\log(w+m)
\ge2^{(1-H_2(2\epsilon+2\tau)+o(1))w}.
```

First let `w` grow and then `tau` decrease to zero, proving (16.59).
`square`

The `m^2` term counts possible ternary-normal optimizer cells; whether its
robust macroscopic form can be reduced to `O(m log m)` or is sometimes
necessary is open.  A bare `B`-bit decoder admits no stronger general law
than `2^B` possible responses: it can hardwire any chosen packing.  The
nontrivial resource in (16.55) is a shared-parameter max-affine grammar.
The same proof applies to weighted finite-domain CSPs with Boolean-valued
factors and to acyclic max-plus path models with binary resource-incidence
vectors.  Its mechanism is classical hyperplane-arrangement geometry plus
finite-dimensional metric entropy; (16.58)--(16.59) are promoted as
semantic response-complexity consequences, not as an external novelty or a
computational lower bound.

### Theorem 16.9 (query entropy and exact tropical lumpability)

Let a `p`-state max-plus weighted automaton have forward vector `v`, and let
`H subset R^p` be the declared set of reachable suffix vectors.  Its suffix
response is

```math
F_v(h)=\max_i\{v_i+h_i\}.                         \tag{16.60}
```

For `v in [-B,B]^p`, let `R_B(H)` be the resulting class of functions on
`H`, and give `H` the projective metric

```math
d_{\rm pr}(h,k)=\inf_c\|h-k-c\boldsymbol1\|_\infty.
```

If `N_pr(delta;H)` is its covering number, then the external response
covering number satisfies

```math
\boxed{
\log_2\operatorname {Cov}^{\rm ext}_{\delta+\eta/2}(R_B(H))
\le N_{\rm pr}(\delta;H)\cdot
\log_2\left(\left\lceil{2B\over\eta}\right\rceil+2\right).}  \tag{16.61}
```

Conversely, suppose suffixes `h^(1),...,h^(k)` robustly expose distinct
coordinates `i_1,...,i_k`:

```math
h^{(j)}_{i_j}\ge
\max_{i\ne i_j}h^{(j)}_i+2B.                    \tag{16.62}
```

On the subfamily where those `k` forward coordinates vary in `[-B,B]`,
the residual map is an isometric sup-cube embedding.  Hence, for
`0<epsilon<B`,

```math
\log_2\operatorname {Cov}_\epsilon
=\Theta(k\log(B/\epsilon))                       \tag{16.63}
```

up to universal changes of covering radius.  In particular low affine
dimension of `H` alone does not imply compression: one affine line can
contain `p` suffixes satisfying (16.62).

There is an exact repeatedly composable quotient under a finite
synchronization hypothesis.  Partition the raw states into
`I_1,...,I_r`, choose gauges `c_i`, and assume every transition matrix
`T_ell` satisfies

```math
\max_{j\in I_b}\{T_\ell(i,j)+c_j\}-c_i
=S_\ell(a,b)
\quad(i\in I_a),                                  \tag{16.64}
```

independently of the microscopic `i`.  Suppose the terminal vector has the
form `beta_i=c_i+u_empty(a)` on `I_a`.  Then every suffix vector factors as

```math
h_y(i)=c_i+u_y(a),
\qquad
u_{\ell y}(a)=\max_b\{S_\ell(a,b)+u_y(b)\},       \tag{16.65}
```

and the nonlinear aggregates

```math
P_a=\max_{i\in I_a}\{v_i+c_i\}                    \tag{16.66}
```

are an exact state:

```math
F_v(h_y)=\max_a\{P_a+u_y(a)\},
\qquad
P'_b=\max_a\{P_a+S_\ell(a,b)\}.                  \tag{16.67}
```

Thus a `p`-coordinate presentation can quotient strictly to `r<p`
coordinates and remain closed for arbitrary concatenation depth.  If the
quotient suffixes robustly expose all `r` coordinates, the formal aggregate
response family on `[-B,B]^r` is coarsest there and has matching response
complexity

```math
\log_2\operatorname {Cov}^{\rm ext}_\epsilon
=\Theta(r\log(B/\epsilon)).                       \tag{16.68}
```

For a concrete fixed automaton this is a minimality statement only on its
reachable aggregate set.  The full rate (16.68) requires that set to contain
the displayed aggregate box, or a comparable product packing.

#### Proof

Choose a projective `delta`-net `k_1,...,k_N` in `H` and quantize each
`F_v(k_j)` at mesh `eta`.  For any query `h`, choose `j,c` with
`||h-k_j-c1||_infinity<=delta`.  The known query baseline `c` is restored by
the decoder, and the maximum in (16.60) is one-Lipschitz, giving total error
`delta+eta/2`.  These decoded vectors need not themselves be residuals,
which is why (16.61) is an external cover.  Choosing one class member from
each nonempty external ball gives an internal cover at twice the radius.
Under (16.62), query `h^(j)` returns
`v_(i_j)+h^(j)_(i_j)` throughout the box.  These queries read all `k`
coordinates, while (16.60) is one-Lipschitz in `v`, proving the isometry and
(16.63).  For example,

```math
h(t)_i=C(2it-i^2),\qquad t=1,...,p,\quad C\ge2B,
```

lies on one affine line and exposes coordinate `t`.

The terminal factorization in (16.65) is the base case.  If it holds for
`y`, then (16.64) gives

```math
h_{\ell y}(i)
=c_i+\max_b\{S_\ell(a,b)+u_y(b)\},
```

proving it inductively.  Grouping the maxima in (16.60) proves the first
identity in (16.67).  If `v'=v odot T_ell`, then

```math
\max_{j\in I_b}(v'_j+c_j)
=\max_a\max_{i\in I_a}\{v_i+c_i+S_\ell(a,b)\},
```

which proves the exact aggregate update.  Robust quotient pins give the
same cube isometry argument and (16.68). `square`

The exact residual quotient and derivative law are the classical weighted
Myhill--Nerode/Hankel viewpoint.  The benchmark-level additions are the
query-metric entropy sandwich, its affine-dimension falsifier, and a sharp
rate for a strict nonlinear aggregate that is compatible with every future
derivative.  This is a state other than a boundary-assignment table and a
second application of response exposure outside Max-Cut.

### Theorem 16.10 (gauges, repeatable holonomy, and tropical resets)

There are two finite mechanisms by which an approximate tropical state can
remain accurate at arbitrary depth: exact cohomological cancellation and
recurrent projective reset.  More precisely, the following hold.

**Gauged metric shells.**  On a finite metric space `(Y,d)`, let

```math
D_{\lambda,g}(a,b)=\lambda d(b,g(a))             \tag{16.69}
```

for a bijective isometry `g`.  If

```math
K_t(a,b)=D_{\lambda_t,g_t}(a,b)
 +\phi_{t-1}(a)-\phi_t(b)+c_t,                   \tag{16.70}
```

then

```math
\boxed{
K_1\star\cdots\star K_T
=D_{\lambda_*,G}+\phi_0\otimes\boldsymbol1
-\boldsymbol1\otimes\phi_T+\sum_tc_t,}         \tag{16.71}
```

where `lambda_*=min_t lambda_t` and
`G=g_T circ ... circ g_1`.  Its directed row table is

```math
r(K_{1:T}(a,\cdot),K_{1:T}(a',\cdot))
=\lambda_*d(a,a')+\phi_0(a)-\phi_0(a').          \tag{16.72}
```

Thus all internal perturbations cancel before the directed maximum, and the
uniform row-table error is at most `osc(phi_0)`, with no depth factor.
Projective distance between two rows is exactly `lambda_*d(a,a')`; the
projective error of the whole endpoint kernel, a different metric, is

```math
{\operatorname {osc}(\phi_0)+\operatorname {osc}(\phi_T)\over2}. \tag{16.73}
```

This gauge form has a finite exact recognition test.  With
`E_t=K_t-D_(lambda_t,g_t)`, simultaneous potentials in (16.70) exist if and
only if every within-factor rectangle and every adjacent-interface
circulation vanishes:

```math
\begin{aligned}
E_t(a,b)+E_t(a',b')-E_t(a,b')-E_t(a',b)&=0,\\
E_t(a,x)-E_t(a,x')+E_{t+1}(x,c)-E_{t+1}(x',c)&=0.
\end{aligned}                                    \tag{16.74}
```

This assertion is for complete real-valued interfaces.  On a sparse support,
rectangles must be replaced by all alternating support cycles.

**Repeatable holonomy.**  Let a finite directed graph have edge labels in a
normed real vector space.  The sums over all finite directed walks are
uniformly bounded if and only if every directed closed walk has label zero.
Equivalently, on each strongly connected component the labels are a vertex
coboundary.  A closed walk with holonomy `h ne 0` has `k`-fold holonomy
`kh`, so a nonzero local tolerance on a repeatable cycle cannot imply a
depth-uniform error.  For a regular declared language this criterion applies
to the reachable, co-reachable part of the automaton or product graph that
encodes the declaration.  For projective responses, labels are taken modulo
common scalar baselines.

**Max-plus reset dichotomy.**  For an all-finite real `r by s` matrix `S`,
put

```math
(F_Su)_b=\max_a\{u_a+S_{ab}\},
\qquad d_H([u],[v])={1\over2}\operatorname {osc}(u-v),          \tag{16.75}
```

and

```math
\Delta(S)={1\over2}\max_{b,c}
\left[\max_a(S_{ab}-S_{ac})-
      \min_a(S_{ab}-S_{ac})\right].              \tag{16.76}
```

On the full projective domain `R^r/R1`, the image diameter is `Delta(S)` and

```math
\boxed{
\operatorname {Lip}_H(F_S)=
\begin{cases}
0,&S_{ab}=\alpha_a+\beta_b,\\
1,&\text{otherwise}.
\end{cases}}                                     \tag{16.77}
```

Hence finite max-plus linear maps have no nontrivial global projective
contraction coefficient strictly between zero and one.  The conclusion is
false on an arbitrary restricted reachable subset.

Small image diameter nevertheless yields a different depth-uniform law.
Let `F_ell:X->X` be raw updates, `G_ell:Y->Y` nonexpansive quotient updates,
and suppose

```math
d(\pi F_\ell x,G_\ell\pi x)\le\epsilon.          \tag{16.78}
```

If a word `w` of length `m` has
`diam(G_w(Y))<=rho`, then every declared word `u=pwq` satisfies

```math
\boxed{
d(\pi F_u x,G_u\pi x)
\le\rho+(m+|q|)\epsilon.}                        \tag{16.79}
```

Consequently, if completed reset occurrences leave tails of length at most
`G` and the first completion is controlled, the error is at most
`rho+(m+G)epsilon` at every later depth.  In particular it is at most
`rho+2L epsilon` when `m,G<=L`.

#### Proof

The potentials `-phi_t` and `+phi_t` cancel at every internal minimization;
Theorem 16.3's shell identity then proves (16.71).  Terminal potentials and
scalar baselines cancel between rows, and a metric centre exposes (16.72).
For (16.74), a vanishing rectangle is exactly the additive-separability
criterion `E_t(a,b)=r_t(a)+s_t(b)`.  The second line says
`s_t(x)+r_(t+1)(x)` is constant, which is precisely the compatibility needed
to choose one shared `phi_t`.

For the graph claim, boundedness applied to repeated closed walks forces
zero holonomy.  Conversely, delete zero-label closed subwalks from any walk;
the remaining simple walk has at most `|V|-1` edges.  On a strongly connected
component, path integration from a root produces the stated potential.

Max-plus monotonicity and additive homogeneity prove nonexpansiveness.  Every
output-coordinate difference is bounded by the corresponding extrema in
(16.76), while making one input coordinate dominant realizes each row in the
limit, proving the diameter formula.  On a full-dimensional selector cell,
the derivative selects one input coordinate for each output.  Its projective
norm is zero if all outputs select the same input and one otherwise.  A
nonconstant continuous piecewise-affine projective map has a cell of the
second kind, proving (16.77); constant projective image is equivalent to the
displayed additive rank-one form.

Finally, iterate (16.78) through the reset word starting from the actual
encoded state, paying `m epsilon`.  The two quotient reset images are within
`rho` regardless of the pre-reset error.  Nonexpansiveness through the tail
adds `|q|epsilon`, proving (16.79). `square`

For a binary zero-temperature Ising link `S_J(s,t)=Jst`, formula (16.76)
gives `Delta(S_J)=2|J|`.  Every nonzero link has global coefficient one, but
a weak link resets all prior projective memory to diameter `2|J|`.  Thus a
weak link every `L` steps yields the uniform bound
`2|J|+2L epsilon`, rather than a depth times error.

The transition-toll family of Theorem 16.7 has rectangular defect
`2delta` on two adjacent indices and is not a gauge.  Its linear drift is
therefore the repeatable-holonomy side of this dichotomy.  The mechanisms in
this theorem are elementary max-plus algebra and graph cohomology; the
generative conclusion is the sharp separation between static response
compression and a state that can actually be reused indefinitely.
