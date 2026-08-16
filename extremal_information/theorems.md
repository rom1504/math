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
