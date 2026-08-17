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
=\Theta(k\log(1+B/\epsilon))                     \tag{16.63}
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
=\Theta(r\log(1+B/\epsilon)).                     \tag{16.68}
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

### Theorem 16.11 (landmarks and balanced query exposure)

Let `(X,d)` be a finite metric space of diameter `D`, and let
`F subset Lip_1(X)/R1`, with shape metric `d_sh=osc/2`.  Write `N_X(r)` for
the least size of an `r`-net.  Define the balanced exposure dimension
`E_gamma(F)` as the largest even `k` for which there are queries
`x_1,...,x_k`, thresholds `a_1,...,a_k`, and, for every balanced
`U subset [k]`, a response representative `f_U` satisfying

```math
f_U(x_i)\ge a_i+\gamma\ (i\in U),
\qquad
f_U(x_i)\le a_i-\gamma\ (i\notin U).             \tag{16.80}
```

Then, for `r,eta>0`,

```math
\boxed{
\log_2\operatorname {Cov}^{\rm ext}_{r+\eta/2}(F)
\le N_X(r)
\log_2\left(\left\lceil{D\over\eta}\right\rceil+1\right),}  \tag{16.81}
```

whereas, for `0<epsilon<gamma`,

```math
\boxed{
\log_2\operatorname {Cov}^{\rm ext}_{\epsilon}(F)
\ge
\log_2{E_\gamma(F)\choose E_\gamma(F)/2}.}       \tag{16.82}
```

For the full Lipschitz response language this lower certificate is computed
exactly from the query geometry:

```math
\boxed{
E_\gamma(Lip_1(X)/R1)
=2\left\lfloor{P_X(2\gamma)\over2}\right\rfloor,}             \tag{16.83}
```

where `P_X(s)` is the largest cardinality of a non-strictly `s`-separated
subset of `X`.  Thus the upper and lower objects are query landmarks and
query packings, rather than response covering numbers hidden in new notation.

There are two benchmark consequences.

1. For the unit-load pure-Max-Cut language on the projective `w`-cube,

   ```math
   E_{\epsilon w}
   \ge2^{(1-H_2(2\epsilon)+o(1))w},
   \qquad 0<\epsilon<1/4.                        \tag{16.84}
   ```

   Combining (16.82) with the shared-parameter bound in Theorem 16.8 gives

   ```math
   E_{\epsilon w}(C_{w,m})
   -\log_2(E_{\epsilon w}(C_{w,m})+1)
   \le O_\epsilon(m^2+m\log(w+m)).               \tag{16.85}
   ```

2. For max-plus residuals
   `g_v(h)=max_i(v_i+h_i)-max_i h_i`, `v in [-B,B]^p`, let `k_pin` be the
   number of coordinates robustly exposed with margin `2B`.  Ordinary fat
   dimension obeys

   ```math
   k_{\rm pin}\le\operatorname {fat}_\gamma\{g_v\}
   \le\min\{p,P_H(\gamma)\},
   \qquad0<\gamma<B.                             \tag{16.86}
   ```

   The exposed subfamily is an isometric sup cube, so for `0<epsilon<B`,

   ```math
   \max\{2,B/\epsilon\}^{k_{\rm pin}}
   \le\operatorname {Cov}^{\rm ext}_\epsilon
   \le(1+2B/\epsilon)^p                          \tag{16.87}
   ```

   up to universal changes of radius.  If all `p` coordinates are exposed,
   its response entropy is
   `Theta(p log(1+B/epsilon))`.  In particular, one affine line of suffixes
   can expose every coordinate; affine dimension alone is not the right
   complexity parameter.

#### Proof

Normalize `min f=0`, so `0<=f<=D`.  Store rounded values on an `r`-net and
take the midpoint of the upper and lower McShane envelopes.  Its error from
`f` has oscillation less than `2r+eta`, and each landmark has at most
`ceil(D/eta)+1` values, proving (16.81).  The decoder need not lie in `F`,
which accounts for the external cover; internal centres cost at most twice
the radius.

Two distinct balanced patterns differ in both orientations.  Their response
difference is at least `2gamma` at one query and at most `-2gamma` at
another, giving shape separation at least `2gamma` and proving (16.82).
Conversely, two queries that can be oppositely labelled in balanced patterns
must be at distance at least `2gamma`.  On any non-strictly
`2gamma`-separated query set, the values `+gamma` and `-gamma` assigned to a
balanced split are one-Lipschitz and extend by McShane.  This proves (16.83).

Theorem 16.4 identifies the Max-Cut response class with the full projective
Lipschitz ball, and a projective Hamming code proves (16.84).  Equations
(16.82) and (16.58) give (16.85).  For automata, `g_v` is two-Lipschitz in
the projective suffix metric, giving the packing term in (16.86).  Its
threshold subgraphs are complements of upper orthants in `R^p`, whose VC
dimension is `p`; hence the other upper bound.  Robust pins read the chosen
coordinates exactly.  Volume and corner packings of that cube, plus direct
coordinate quantization, prove (16.87). `square`

The one-scale exposure dimension does not determine entropy by itself: even
a one-dimensional response interval has constant fat dimension and
`Theta(log(1+B/epsilon))` precision cost.  Nor does this static theorem make
a cover reusable under composition.  The transition-toll family has small
static error and macroscopic depth drift; Theorem 16.10 supplies the
additional congruence/reset mechanisms.

### Theorem 16.12 (robust binary max-affine entropy lower bound)

Let `P=conv(V)` be a full-dimensional `0/1` polytope in `R^m`.  There is a
binary max-affine presentation with `m` shared real parameters and one query
for `P` plus one query for each facet, containing response shapes
`[f_(theta_F)]` such that

```math
\|[f_{\theta_F}]\|_{\rm sh}={1\over2},
\qquad
d_{\rm sh}([f_{\theta_F}],[f_{\theta_G}])=1
\quad(F\ne G).                                   \tag{16.88}
```

The parameters can be perturbed off every optimizer-tie hyperplane,
positively rescaled back to radius `1/2`, and kept pairwise more than
`1-2eta` apart for any `eta>0`.

There are full-dimensional `0/1` polytopes with at least

```math
\left({cm\over(\log m)^2}\right)^{m/2}           \tag{16.89}
```

facets for an absolute `c>0`.  Consequently, for every fixed
`delta<1/2`, some radius-`1/2` binary shared-parameter response class has

```math
\boxed{
\log\operatorname {Cov}^{\rm ext}_\delta
\ge {m\over2}(\log m-2\log\log m-O(1))
=\Omega(m\log m).}                               \tag{16.90}
```

Together with Theorem 16.8, the known general radius-bounded range is
therefore

```math
\Omega(m\log m)
\le\sup_{A_\bullet}
 \log\operatorname {Cov}_\delta
 \{[f_\theta]:\|[f_\theta]\|_{\rm sh}\le1\}
\le O(m^2).                                      \tag{16.91}
```

This is a statement for generic binary max-affine presentations, not yet for
unit-load Max-Cut components.

#### Proof

Use witness set `V` for a base query.  For every facet `F`, with vertex set
`V_F=V cap F`, use the deletion query `V setminus V_F`.  If `u_F` exposes
`F`, normalize it by the positive gap between `F` and the undeleted
vertices.  The base response and every other facet-deletion response equal
one common scalar `c_F`, while the `F`-deletion response equals `c_F-1`.
Thus

```math
f_{\theta_F}=c_F\boldsymbol1-e_F,                 \tag{16.92}
```

which proves (16.88).  Continuity and density off finitely many comparison
hyperplanes give the robust version.  The facet lower bound (16.89) is due
to Gatzouras--Giannopoulos--Markoulakis; substituting it in (16.88) proves
(16.90). `square`

The primary sources are
[Gatzouras--Giannopoulos--Markoulakis](https://arxiv.org/abs/math/0406125)
for (16.89) and
[Fleiner--Kaibel--Rote](https://doi.org/10.1006/eujc.1999.0326)
for the `exp(O(m log m))` face bound on one `0/1` polytope.  Hence the
single-polytope deletion mechanism has the `m log m` exponent.  The gap in
(16.91) concerns robust response shapes from common refinements of many
support polytopes: raw arrangement cells need not remain separated after
projective normalization, and intermediate exponents are not excluded.

### Theorem 16.13 (mean-field response quotient and depth-stable rate)

Let a binary block be a finite multiset of local scores
`A={h_i} subset [-B,B]`, anchored by `H_A(0)=0`.  The allowed contexts append
an anonymous block and apply one scalar field `lambda` to total occupancy.
If `a_1>=...>=a_n` are the sorted fields, define

```math
p_A(k)=sum_(j<=k)a_j,
\qquad
R_A(lambda)=max_k{p_A(k)+lambda k}.
```

Then `R_A`, the discrete-concave profile `p_A`, and its slope multiset are
equivalent coarsest exact contextual states.  More precisely,

```math
p_A(k)=inf_(lambda in R){R_A(lambda)-lambda k},      \tag{16.93}
```

```math
sup_lambda|R_A(lambda)-R_(A')(lambda)|
=max_k|p_A(k)-p_(A')(k)|                            \tag{16.94}
```

for equal masses, and anonymous union obeys

```math
p_(A sqcup C)(t)
=max_(k+l=t){p_A(k)+p_C(l)}.                        \tag{16.95}
```

In slope coordinates, (16.95) is sorted multiset union.

For `eta>0`, nearest rounding to one common grid of spacing at most `eta`
gives an exactly additive histogram with

```math
|S_(n,eta)|={n+M-1 choose M-1},
\qquad M=1+ceil(2B/eta),                             \tag{16.96}
```

when `B>0`; for `B=0` the state has one bin.  It uses
`O((1+B/eta)log(n+1))` bits and satisfies

```math
max_k|p_A(k)-p_tilde_A(k)|
vee sup_lambda|R_A(lambda)-R_tilde_A(lambda)|
<=eta n/2.                                          \tag{16.97}
```

On any merge tree of total mass `N`, the root error is at most `eta N/2`,
independently of depth: each microscopic field is rounded once.  Distinct
grid histograms are response-separated by at least the grid spacing, so
(16.96) is also the exact grid-state lower bound and remains necessary below
half-grid error.  At macroscopic error `epsilon n`, a disjoint-tent packing
gives the weaker general lower bound below, where `K_delta` is the minimum
number of deterministic summary states whose decoded responses have uniform
error at most `delta`:

```math
log_2 K_(epsilon n)
>=Omega(min{n,sqrt(B/epsilon)}),                     \tag{16.98}
```

for `B>0`, `epsilon>0`, and
`min{n,sqrt(B/epsilon)}` above a universal constant.  Equivalently, the
right side means `c min{n,sqrt(B/epsilon)}` for a universal `c>0` in that
range.  No matching macroscopic rate is claimed.

Now add the fixed pair score

```math
q_A(k)=p_A(k)+J {k choose 2}
```

at every composition stage, and let `bar q_A` be its piecewise-linear least
concave majorant on `[0,n]`.  For realizable roofs `f` on `[0,n]` and `g` on
`[0,m]`, and `0<=t<=n+m`, set

```math
f star_J g
=cav_t max_(substack{0<=u<=n,0<=v<=m; u+v=t})
  {f(u)+g(v)+Juv}.                                   \tag{16.99}
```

Then

```math
bar q_(A sqcup C)=bar q_A star_J bar q_C,            \tag{16.100}
```

`star_J` is associative on realizable roofs, and `(n,bar q_A)` is the
coarsest exact state for repeated same-`J` merges and a terminal linear
field.  This quotient can be strict: for
`0<a<min(B,J/2)`, the profiles `(0,0,J)` and `(0,a,J)` have the same
endpoint-chord roof.

The uniform fixed-mass collapse threshold is sharp.  Since

```math
p_A(k)-{k\over n}p_A(n)
<={2B k(n-k)\over n},                               \tag{16.101}
```

every `J>=4B/n` makes `bar q_A` the endpoint chord, so at mass `n` the exact
state is only `sum_i h_i`.  For `J<4B/n`, taking `k` fields equal to `B` and
the rest equal to `-B` violates that chord.  Thus `J>=2B` is the
size-uniform sufficient threshold for all masses at least two.  Endpoint
optimizer uniqueness at equality is not asserted.

#### Proof

Exchange puts the `k` largest fields in a fixed-occupancy optimizer, so the
slopes of `p_A` are the sorted fields.  Every `k` is supported by a scalar
field between `-a_k` and `-a_(k+1)`, proving (16.93).  Maxima are
nonexpansive for (16.94) in one direction; applying (16.93) gives the other.
Splitting a size-`t` subset between two blocks proves (16.95), while the
empty future proves contextual minimality.

The grid state adds coordinatewise.  Every selected set changes score by at
most `eta/2` per occupied site; maxima preserve that bound, and no merge
requantizes a site.  Stars and bars proves (16.96).  Distinct arithmetic-grid
histograms have different sorted slopes, and their first differing prefix
sum is separated by one grid spacing.  For (16.98), use disjoint triangular
response perturbations of height at least `Bn/(16q^2)` with
`q=Theta(min{n,sqrt(B/epsilon)})`.

For (16.100), take independent occupancy mixtures attaining the two child
roofs.  Separate affinity gives expected cross score `Juv`; pure occupancies
give the converse.  Both sides therefore have, for every `lambda`, response

```math
max_(u,v){bar q_A(u)+bar q_C(v)+Juv+lambda(u+v)}.
```

Linear biconjugacy identifies their roofs.  For three children this response
contains the symmetric term `J(uv+uz+vz)`, proving associativity.  Finally,
(16.101) follows by maximizing the deviation with `k` fields at `B` and the
rest at `-B`; subtracting the quadratic endpoint chord gives the sharp
threshold. `square`

The theorem uses a maximized score.  Physical minimum-energy conventions
reverse signs.  The full `lambda in R` family is essential unless mass and a
sufficient model-dependent range are separately declared.  Separately
addressable old blocks, changing `J/n` normalizations, and non-biaffine cross
terms need larger states; (16.100) does not cover them.

### Theorem 16.14 (algebraic absorption and the selector-reset converse)

Depth-uniform continuation has different converses for coherent fixed maps
and for fresh adversarial residuals.

First, let a finite alphabet generate a finite semigroup `S`, let

```math
L=max_(s in S) min{|w|:w represents s},              \tag{16.102}
```

and let `F,G` be two exact actions of `S` by nonexpansive maps on the same
metric space.  If corresponding generators are uniformly `epsilon`-close,
then every word, regardless of its written length, satisfies

```math
sup_x d(F_wx,G_wx)<=L epsilon.                       \tag{16.103}
```

This bounded-normal-form mechanism is not recognized by entrywise kernel
gauges or small full-image resets.  Indeed, in projective coordinate
`z=u_2-u_1`, with `d_H(z,z')=|z-z'|/2`, the all-finite max-plus matrices

```math
S_0=((0,0),(-1,0)),
\qquad S_delta=((0,delta),(-1,0))                   \tag{16.104}
```

induce the idempotents

```math
P_0(z)=clip(z,0,1),
\qquad P_delta(z)=clip(z,delta,1).                  \tag{16.105}
```

For every `t>=1`,

```math
sup_z d_H(P_0^t z,P_delta^t z)=delta/2.             \tag{16.106}
```

Nevertheless, `S_delta-S_0` has rectangular circulation `-delta`, and the
two positive-power image diameters are `1/2` and `(1-delta)/2`.  Thus the
pair is neither an entrywise endpoint gauge nor an `O(delta)` full-image
reset.  This refutes only that narrow dichotomy: the paired clamp orbit is
stationary after one transient and can be viewed as a finite zero-increment
recurrence.

For a sharp robust statement, let `r>=2`,
`V=R^r/R 1`, `||[v]||_H=osc(v)/2`, and for a selector
`sigma:[r]->[r]` set `(P_sigma v)_j=v_(sigma(j))`.  In a factorial declared
language, call a contiguous selector product a **tangent reset** when its
composite selector is constant.  For arbitrary disturbances

```math
e_t=P_t e_(t-1)+eta_t,
\qquad e_0=0,
\qquad ||eta_t||_H<=epsilon,                         \tag{16.107}
```

the following quantitative converse holds.

1. If every allowed word of length `L` contains a tangent-reset factor,
   then `||e_T||_H<=L epsilon` for every `T>=L`.
2. If there is a reset-free allowed word of length `T`, disturbances can be
   chosen so that

   ```math
   ||e_T||_H
   >=floor(T/[r(r-1)])epsilon.                       \tag{16.108}
   ```

Consequently, stability by `C epsilon` against all residuals forces every
reset-free word to have length less than `(C+1)r(r-1)`.  If the recursion
also contains an endpoint term

```math
h_t-P_t h_(t-1),                                    \tag{16.109}
```

subtracting `h_t` reduces it to (16.107).  Once the syndetic hypothesis
supplies a completed reset in the final length-`L` window, the bound is
`||h_T||_H+L epsilon`; before that, the transported initial endpoint term
must also be paid.

Finally, on one recurrent affine-selector cell, let
`A(e)=P_sigma e+b`.  Its projective iterates are bounded for every `e` if
and only if every directed cycle `C` of the functional graph of `sigma` has
one common `b`-mean.  Equivalently,

```math
b=p-P_sigma p+beta 1                                \tag{16.110}
```

for some `p,beta`, in which case

```math
A^k(e)=p+P_sigma^k(e-p)+k beta 1.                   \tag{16.111}
```

Different cycle means cause linear projective drift.  Thus the exact
holonomy is twisted by the active selector; ordinary untransported label
sums are not sufficient.

#### Proof

For (16.103), replace `w` by a representative of its semigroup element of
length at most `L`.  Both exact actions make the replacement without error;
hybridizing its factors costs at most `epsilon` each by nonexpansiveness.
Direct max-plus calculation proves (16.105), hence idempotence and
(16.106); the alternating rectangle of the kernel difference is `-delta`.

Unrolling (16.107) gives

```math
e_T=sum_(s=1)^T P_T...P_(s+1) eta_s.                \tag{16.112}
```

A reset in the last length-`L` window kills all earlier summands, leaving at
most `L`.  Conversely, every suffix of a reset-free word is a nonconstant
selector.  Choose for each suffix an ordered output pair sent to distinct
input coordinates.  One of the `r(r-1)` pairs recurs at least
`floor(T/[r(r-1)])` times.  At each selected time insert
`+epsilon,-epsilon` on the two distinct input coordinates selected by that
output pair.  All contributions have the same sign on the common final pair,
proving (16.108).

For (16.110), summing `p_j-p_(sigma(j))=b_j-beta` around a cycle gives the
common-mean necessity.  If the means agree, define `p` consistently around
each cycle and then recursively on its incoming trees.  Equation (16.111)
telescopes, while unequal cycle means create different linear coordinate
drifts. `square`

### Theorem 16.15 (arithmetic feature-algebra growth law)

Let `phi_1,...,phi_d` be bounded real response functions on a declared query
set.  At mass `n`, let

```math
H_(n,d)={c in N^d:sum_j c_j=n},
\qquad F_c=sum_j c_j phi_j,                          \tag{16.113}
```

and compose systems by histogram addition.  On
`V={z in R^d:sum_j z_j=0}`, put `Tz=sum_j z_j phi_j` and define

```math
Gamma_Phi=<phi_1-phi_d,...,phi_(d-1)-phi_d>_Z,
\qquad r_Z=rank_Z Gamma_Phi.                        \tag{16.114}
```

For projective responses, take this group and the norms below intrinsically
in the quotient by constant query functions.  On each fixed-mass slice,
for `c,c' in H_(n,d)`, contextual equivalence is exactly

```math
c sim c' iff T(c-c')=0,                             \tag{16.115}
```

it is a congruence for every future histogram addition, and the number
`N_n` of exact mass-`n` contextual states satisfies

```math
\boxed{N_n=Theta_Phi(n^(r_Z)).}                     \tag{16.116}
```

There is a robust version.  Let

```math
alpha=inf_(z in V,||z||_infty=1)||Tz||_infty,
\qquad L=max_j||phi_j||_infty,
\qquad sigma=inf_(0 ne z in V cap Z^d)||Tz||_infty. \tag{16.117}
```

If `d>=2` and `alpha>0`, the full histogram is the coarsest exact state,
`r_Z=d-1`, and, with
`s_epsilon=1+floor(2epsilon/alpha)`,

```math
(1+floor(n/[(d-1)s_epsilon]))^(d-1)
<=Cov^ext_epsilon{F_c:c in H_(n,d)}
<=(2+ceil(L(d-1)n/epsilon))^(d-1).                 \tag{16.118}
```

If `sigma>0`, then for every `0<epsilon<sigma/2`,

```math
\boxed{
Cov^ext_epsilon{F_c:c in H_(n,d)}
={n+d-1 choose d-1}.}                               \tag{16.119}
```

For the equally spaced heterogeneous mean-field atoms with `d>=2` and
`B>0`,

```math
gamma_j=-B+(j-1)Delta,
\qquad Delta={2B\over d-1},
\qquad phi_j(lambda)=(gamma_j+lambda)_+,
\quad lambda in [-B,B],                             \tag{16.120}
```

one has

```math
alpha>=Delta/4,
\qquad sigma=Delta                                  \tag{16.121}
```

in the literal anchored sup norm.  In the projective half-oscillation norm,
the corresponding safe constants are `alpha>=Delta/8` and
`sigma=Delta/2`.  Consequently the literal exact external response-cover
count below half-grid error is the binomial in (16.119), while (16.118)
gives polynomial two-sided response-rate bounds at coarser scales.

#### Proof

At fixed mass,

```math
F_c=n phi_d+sum_(j<d)c_j(phi_j-phi_d).              \tag{16.122}
```

This proves (16.115), and addition preserves it.  Identify the finitely
generated torsion-free group `Gamma_Phi` with `Z^(r_Z)`.  Its atom
generators have bounded integer coordinates, so every mass-`n` sum lies in
an `O_Phi(n)` box, proving the upper half of (16.116).  Choose `r_Z` of the
displayed differences independent over `Q`, vary each count up to
`floor(n/r_Z)`, and put all unused mass in type `d`.  The resulting sums are
distinct, proving the lower half; rank zero gives one state.

When `alpha>0`, `T` is injective on `V`.  For the lower bound in (16.118),
vary the first `d-1` counts in multiples of `s_epsilon`; distinct choices
are more than `2epsilon` apart in response norm.  For the upper bound,
externally round those counts on mesh `epsilon/[L(d-1)]` and adjust the last
coordinate to retain total mass.  Its response error is at most `epsilon`.
Distinct integer histograms are `sigma`-separated, which proves (16.119).

For (16.121), at the knots `lambda=-gamma_k`, write

```math
{Tz(-gamma_k)\over Delta}
=S_k=sum_(j>k)z_j(j-k).                             \tag{16.123}
```

The tail sums satisfy `sum_(j>=k+1)z_j=S_k-S_(k+1)`, hence
`||z||_infty<=4max_k|S_k|` and `alpha>=Delta/4`.  If `z` is a nonzero
integer vector, some `S_k` is a nonzero integer, giving norm at least
`Delta`.  Moving one site between adjacent bins attains equality, so
`sigma=Delta`. `square`

The arithmetic rank and robust conditioning are deliberately separate.  A
one-real-parameter query can assign rationally independent atom responses
and have large exact `r_Z` while `alpha=0`; unlimited-precision exact state
growth then need not represent a robust exposed direction.  Conversely, the
additive update is an exact depth-stable congruence, not merely a static
metric cover.

### Theorem 16.16 (finite suffix-product stability certificate)

Let all finite paths of a finite directed graph be legal selector
trajectories, with every edge labelled by a map `sigma:[r]->[r]`, `r>=2`.
Initialize the lift at `(v,emptyset)` for every permitted path-start vertex
(at every vertex when all graph paths are legal).  For a reset-free path
prefix ending at vertex `v`, let `S` be the set of composite selectors of all
its nonempty suffixes.  On appending an edge `sigma`, update

```math
S'={sigma} union {rho circ sigma:rho in S},          \tag{16.124}
```

and discard the transition if `S'` contains a constant map.  This defines a
finite lifted graph on at most

```math
|V_graph| 2^(r^r-r)                                 \tag{16.125}
```

states.

There are reset-free legal paths of arbitrarily large length if and only if
the reachable lifted graph contains a directed cycle.  If it is acyclic and
`H` is its maximum path length in edges, every legal word of length `H+1` contains a
tangent-reset factor.  Hence the arbitrary-residual recursion (16.107)
satisfies

```math
||e_T||_H<=(H+1)epsilon                             \tag{16.126}
```

once a final window of that length exists.  If a reachable lifted cycle
exists, it is an explicit instability certificate: pumping it gives
arbitrarily long reset-free words, and (16.108) produces errors growing
linearly with word length.

#### Proof

Induction shows that (16.124) stores exactly all products of suffixes ending
at the newest edge.  A new constant product is precisely a tangent-reset
factor ending there; every earlier factor was checked when its last edge was
appended.  Thus every reset-free legal path has a unique lift, and every
lifted path projects to one.  A finite directed graph has paths of unbounded length exactly
when a reachable directed cycle exists.  The acyclic height and the two
claims now follow from Theorem 16.14. `square`

This decides robust depth stability on a fixed tie-free selector language.
It does not yet decide switching perturbations which change the active cell;
that requires a paired-cell selector automaton.

### Corollary 16.16a (minimal kernel-partition reset state)

The suffix-set lift in Theorem 16.16 is correct but not minimal. If `rho` is
the whole-path selector product, store only its kernel partition `Pi`. Under
the composition convention in (16.124), appending `sigma` updates

```math
Pi'=sigma^(-1)Pi,
\qquad i~_(Pi')j iff sigma(i)~_Pi sigma(j).         \tag{16.126a}
```

The one-block partition is a reset sink. This recognizes exactly the
reset-free paths with at most

```math
|V_graph|(Bell(r)-1)+1                             \tag{16.126b}
```

states. At a fixed control vertex with the full transformation alphabet, the
`Bell(r)` partition states, including the sink, are the coarsest exact
future-reset quotient. This does not claim that distinct control vertices in
an arbitrary regular-language graph are themselves minimal.

Indeed, constant maps form a two-sided ideal. A word has a constant-product
factor if and only if its whole product is constant, and
`ker(rho circ sigma)=sigma^(-1)(ker rho)`. If two partitions differ, choose
`a,b` joined in exactly one and a continuation selector with image `{a,b}`;
it resets one state and not the other. A restricted language may admit a
further Myhill--Nerode quotient. For paired product selectors on `[r]^2`,

```math
ker(tau times sigma)=ker(tau) times ker(sigma),     \tag{16.126c}
```

so componentwise pullbacks on two channel partitions suffice rather than a
general partition of `r^2` coordinates. For the full independent product-
selector alphabet this gives `Bell(r)^2-1` accepting pairs plus one sink and
is worst-case minimal for the full reset query; restricted or diagonal-
observation languages can quotient further.

### Theorem 16.17 (depth-stable atomic type quantization)

Let an atom `z` contribute a vector `phi_z` in a normed response space, and
let a mass-`n` additive system have response

```math
F_(z_1,...,z_n)=sum_(i=1)^n phi_(z_i).              \tag{16.127}
```

Composition concatenates atom multisets.  Fix one root-scale external
`eta`-net of size `D` and one deterministic atom-to-centre quantizer, used at
every leaf of the target merge tree.  Store the resulting `D`-bin histogram.
This summary updates exactly by
histogram addition, has at most

```math
{n+D-1 choose D-1}                                  \tag{16.128}
```

mass-`n` states, and its decoded response satisfies

```math
||F-F_tilde||<=n eta.                               \tag{16.129}
```

The error depends on total mass, not merge depth or bracketing.  External
centres are allowed: the decoded vector need not be the response of a
physical atom system.

Writing `D(eta)` for the atom-response covering number, if

```math
eta_n -> 0,
\qquad D(eta_n)=o(n),                               \tag{16.130}
```

then both response error per atom and summary bits per atom vanish.  In
particular, if `D(eta)<=C eta^(-p)` for some `p>0`, choosing
`eta_n=n^(-a)` with `0<a<1/p` gives response error `n^(1-a)` and
`O(n^(ap)log n)` summary bits.

#### Proof

Replace each atom response by its assigned centre.  The triangle inequality
proves (16.129).  Concatenation adds the counts and never requantizes an
atom, so no depth factor appears.  Stars and bars gives (16.128), while

```math
log_2 {n+D-1 choose D-1}
<=(D-1)log_2(e(1+n/(D-1)))                          \tag{16.131}
```

for `D>=2` (the count is one for `D=1`).  Dividing by `n` and writing
`t=n/(D-1)` gives `log_2(e(1+t))/t -> 0`, which proves (16.130).  The
polynomial-cover specialization follows by substitution. `square`

The theorem holds intrinsically for projective responses when the net and
norm are taken in the quotient.  It is only an upper theorem: without the
arithmetic-rank, conditioning, or lattice-margin hypotheses of Theorem
16.15, different histograms may have identical or nearly identical total
response.  Heterogeneous local fields satisfy it with atom responses
`phi_h(lambda)=(h+lambda)_+` on the sufficient bounded query interval
`lambda in [-B,B]` at fixed known mass, and have one-dimensional atom-cover
growth.  Equivalently, for unrestricted `lambda`, subtract the common
baseline `(lambda)_+` and retain mass separately.

The common-net hypothesis is essential.  If children use coarser
mass-dependent nets and the parent later refines them, distinctions already
collapsed at a child cannot be reconstructed without additional refinement
data.

As a second model, let `v_i` lie in the Euclidean unit ball of fixed
dimension `p` and query the support of their signed-sum zonotope:

```math
F_V(theta)=max_(epsilon_i in {+-1})
 <theta,sum_i epsilon_i v_i>
=sum_i|<theta,v_i>|,
\qquad theta in S^(p-1).                             \tag{16.132}
```

Here `phi_v(theta)=|<theta,v>|` and
`||phi_v-phi_w||_infty<=||v-w||_2`.  A Euclidean `eta`-net has
`D(eta)<=(1+2/eta)^p`, so for every `0<a<1/p`, a composable type histogram
uses `O(n^(ap)log n)` bits and preserves every directional support query to
error at most `n^(1-a)`.  Thus Theorem 16.17 is not specific to mean-field
hinges.

The binomial upper bound and the `D=o(n)` threshold are sharp without
additional structure.  For query set `[D]`, take coordinate atoms

```math
phi_j(q)=1_(j=q).
```

Then `F_c(q)=c_q`, distinct mass-`n` histograms have sup distance at least
one, and for every `0<epsilon<1/2`,

```math
Cov^ext_epsilon={n+D-1 choose D-1}.                 \tag{16.133}
```

If `limsup D/n>0`, the logarithm of (16.133) is `Omega(n)` along an infinite
subsequence.  Hence no theorem based only on the atom-net cardinality can
guarantee vanishing bits per atom beyond the `D=o(n)` regime.

Indeed, the coordinate response identity proves the separation and stars
and bars proves (16.133).  Along a subsequence with `D>=c n`, monotonicity
in `D` and the standard binomial-entropy bound give
`log {n+floor(c n)-1 choose n}=Omega_c(n)`.

### Theorem 16.18 (tie-aware stochastic secants and exact response gain)

Let `r>=2` and let `P_1,...,P_T` be row-stochastic `r by r` matrices acting on
`V=R^r/R1` with `||[v]||_H=osc(v)/2`, and let

```math
e_t=P_te_(t-1)+eta_t,
\qquad ||eta_t||_H<=epsilon_t,
\qquad e_0=0.                                      \tag{16.134}
```

Write `R_(s,T)=P_T...P_(s+1)`, with `R_(T,T)=I`. Then

```math
\boxed{
\sup_(||eta_s||_H<=epsilon_s)||e_T||_H
=\max_(i,j)\sum_(s=1)^T epsilon_s
  TV(R_(s,T)(i,*),R_(s,T)(j,*)).}                  \tag{16.135}
```

If every length-`L` block occurring in a legal trajectory has Dobrushin
coefficient at most `rho<1`, then residuals of size at most `epsilon` satisfy

```math
||e_T||_H<={L\over1-rho}epsilon.                   \tag{16.136}
```

For a factorial language in a finite row-stochastic semigroup, uniform
fresh-residual gain holds if and only if identical-row products occur with
bounded gaps. If a reset-free legal word has length `T`, the semigroup has a
nonconsensus element; let `gamma` be its least positive row-pair total
variation over all semigroup elements. Then some residuals give

```math
||e_T||_H>={gamma T\over r(r-1)}epsilon.            \tag{16.137}
```

#### Proof

Unroll (16.134). The Hilbert norm is the maximum over ordered output pairs
of half their difference. For a zero-sum row difference `d`,

```math
sup_(||v||_H<=epsilon)d*v=epsilon||d||_1.
```

Residuals at different times optimize independently for the same terminal
pair; the outer half proves (16.135). Dobrushin coefficients are
submultiplicative, so grouping suffix lengths into blocks and summing a
geometric series proves (16.136).

In a finite semigroup every distinct row pair has separation at least
`gamma`. Along a reset-free word every terminal suffix is nonconsensus;
sum its Dobrushin coefficients and pigeonhole the maximizing pair among
`r(r-1)` choices to obtain (16.137). An identical-row factor is a two-sided
ideal under row-stochastic multiplication and kills every older residual,
proving the converse. `square`

This theorem applies across max-plus switches. For every all-finite max-plus
map `F_K` and `x,y`, there is a row-stochastic secant

```math
F_K(y)-F_K(x)=P_K[x,y](y-x),                        \tag{16.138}
```

whose rows are convex combinations of selectors active along the
corresponding line segment. Restrict each output maximum to the segment and
integrate its piecewise-constant slope. If trajectories use maps `F_t,G_t`,
decompose their difference into the secant of `G_t` and the same-input
residual `G_t(x)-F_t(x)`; (16.135)--(16.136) apply.

Fractional consensus is absent from the fixed-selector theorem. For
`P_alpha=((1-alpha,alpha),(alpha,1-alpha))`, the gain is at most
`epsilon/(1-|1-2alpha|)` although no finite power has identical rows unless
`alpha=1/2`. A tie alone does not provide this: secant weights may approach
zero. The lower formula permits fresh adversarial residuals and is not a
coherent-kernel converse. Recognizing which secant paths are dynamically
realizable remains a paired-cell reachability problem.

### Theorem 16.19 (observable cycles for exact selector presentations)

Let a finite directed multigraph with permitted starts carry affine-selector
edges

```math
e:q->q',
\qquad A_ez=P_(sigma_e)z+b_e,                     \tag{16.139}
```

and let terminal states declare ordered coordinate observations. Assume
every actual orbit is represented by a graph path (soundness), and initial
coordinate oscillations are at most `R_0`. Form the reverse witness graph

```math
(q',i,j)->(q,sigma_e(i),sigma_e(j))                \tag{16.140}
```

with weight `b_e(i)-b_e(j)`, retaining vertices on a path from a terminal
observation to an allowed start.

If the relevant graph has no positive cycle, every directed output is at
most `R_0+K`, where `K>=0` is the maximum weight of a simple relevant path,
including the empty path. If
every relevant path and cycle repetition is realized by an allowed finite
orbit (path realization), a positive cycle of weight `c` has fixed access,
cycle, and exit words `u,v,w` such that, for every `k`,
`sup_(allowed z_0)D_(uv^kw)(z_0)>=kc-C`. Hence absence of positive cycles is
necessary and sufficient for uniform directed upper boundedness in an exact
regular presentation. The maximizing finite-orbit seed may depend on `k`;
one common infinite orbit requires nested-cylinder realization.

Under the same path-realization assumption, if observations and relevance
are closed under pair reversal, two-sided projective boundedness is
equivalent to zero weight on every relevant cycle, or equivalently to a
vertex coboundary on every relevant strongly connected component.

#### Proof

One edge satisfies

```math
(A_ez)_i-(A_ez)_j
=z_(sigma_e(i))-z_(sigma_e(j))+b_e(i)-b_e(j).       \tag{16.141}
```

Iteration gives a bounded initial difference plus the reverse-path weight.
Deleting nonpositive cycles cannot lower that weight, so a maximum is
attained on a simple path. Under path realization, repeat a positive cycle;
(16.141) adds `c` each time. Pair reversal negates weights, so two-sided
boundedness forces zero cycle sums. The usual path-independence argument
identifies zero cycle sums with vertex potentials on an SCC. `square`

For paired channels with different selectors, diagonal error does not close.
The exact joint carrier is

```math
D_(ij)=y_i-x_j,
\qquad D'_(ij)=D_(tau(i),sigma(j))+t_i-s_j.         \tag{16.142}
```

It is an affine selector system on `r^2` coordinates, and Hilbert error is
observed by pairs `((i,i),(j,j))`. The witness graph has at most `|Q|r^4`
states and preserves joint cancellation before absolute values. This is
smaller than a trajectory language only when `Q` is itself a genuine
quotient. Applying the theorem also requires bounded initial oscillation of
the full cross carrier, not only zero initial diagonal error: if `x=y` is an
arbitrary projective input, off-diagonal `D_(ij)=x_i-x_j` may be unbounded
and differing selectors can expose it. A bounded-image prefix, reset, or
bounded initial projective domain supplies the required endpoint control.

That qualification is sharp. The all-finite map
`z->clip(z-delta,0,1)` has a slope-one cell with a nonempty one-step
self-intersection, yet every orbit exits it after finitely many steps. A
local face graph therefore has a spurious pumpable cycle. A finite invariant
cell partition is path-realizing when each map sends every whole cell into a
declared successor and has one selector-affine formula there; tie faces must
also refine by a common realizable tangent update. One-step feasibility,
including independent tie resolutions, is not enough. The remaining target
is a tropical lumpability theorem producing a small exact presentation.

## 17. Finite tropical lumpability and dynamic response memory

The results in this section give finite criteria and sharp boundaries for the
forward half of the target left by Theorem 16.19. They distinguish an exact symbolic quotient, an approximate
response simulator, and the information in a future response tree. Detailed
counterexample and benchmark arithmetic is recorded in
[`drafts/finite_tropical_lumpability.md`](drafts/finite_tropical_lumpability.md).

### Theorem 17.1 (contextual refinement and an invariant-arrangement certificate)

Let `E` be finite, let `F_e:X->X`, and let `P_0` be a finite observation
partition. Define

```math
x equiv_t y
iff
F_wx and F_wy lie in the same P_0 atom for every |w|<=t.          \tag{17.1}
```

There is a finite deterministic forward congruence refining `P_0` if and
only if the partitions `P_t` stabilize. The first stable partition is the
coarsest such congruence. If a congruence with `k` classes exists,
stabilization takes at most `k-|P_0|` strict refinements.

There is the following finite certificate. Let `X` be contained in a
`d`-dimensional affine space and let `H={h_1,...,h_m}` be oriented affine
forms whose sign partition refines `P_0`. On each nonempty sign atom

```math
C_s={x in X:sgn h_i(x)=s_i, 1<=i<=m},                            \tag{17.2}
```

assume every `F_e` agrees with one affine branch `g_(e,s)`. Suppose
also that, for every `i,e,s`,

```math
h_i circ g_(e,s)
in {0, nonzero constant}
 union {lambda h_j:lambda!=0, 1<=j<=m}.                          \tag{17.3}
```

Then the nonempty sign atoms form a sound, path-realizing forward
congruence. Their number is at most

```math
sum_(j=0)^d 2^j {m choose j}<=3^m.                               \tag{17.4}
```

Every quotient cycle is repeatable by every seed in its source atom. If the
branches are affine selectors, permitted starts have bounded coordinate
oscillation, and terminal observations are quotient-coordinate observations,
Theorem 16.19 decides observable drift on this quotient. The sign of
`lambda` in (17.3) is part of the certificate.

For selector-affine branches

```math
g_sigma(u)_a=u_(sigma(a))+b_a,                                   \tag{17.5}
```

start from finitely many oriented comparison forms
`u_p-u_q+theta`. Construct the finite ordered-pair graph

```math
(p,q)->(sigma(p),sigma(q))
```

with edge weight `b_p-b_q`, terminating when the pair collapses. The
unrestricted affine pullback closure of the starting comparisons is finite
if and only if every reachable directed cycle has weight zero. Hence this
zero-holonomy condition plus branch resolution supplies `H` in (17.3)
without enumerating orbit words.

#### Proof

The recursion

```math
P_(t+1)=P_t wedge bigwedge_(e in E)F_e^(-1)P_t                   \tag{17.6}
```

shows that a stable partition is a forward congruence. Conversely, every
forward congruence refining `P_0` refines every `P_t`. Therefore `|P_t|<=k`
when a `k`-class congruence exists, so every strictly increasing finite
refinement sequence stabilizes. The same observation proves coarseness.

On `C_s`, the branch in (17.3) is fixed. Every target sign is then a constant,
zero, or the sign of one source form, possibly reversed. Hence
`F_e(C_s)` lies in one sign atom. Induction realizes every labelled path and
cycle repetition. The hyperplane-arrangement face bound is (17.4).

Finally,

```math
(u_p-u_q+theta) circ g_sigma
=u_(sigma(p))-u_(sigma(q))+theta+b_p-b_q.                         \tag{17.7}
```

A pullback path changes the offset by its path weight. If every cycle has
zero weight, every repeated-vertex closed subwalk decomposes into reachable
directed cycles and has zero weight, so delete it without changing the
resulting form; only finitely many simple paths remain. A reachable cycle of nonzero
weight `c` produces distinct offsets `theta+kc`. `square`

The affine qualifications are necessary. Normal directions alone fail for
`F(x)=x+1`, because the negative half-line crosses its boundary. Unit
transport is also essential: `F(x)=2x` pulls `x-1` back to the infinitely
many boundaries `x=2^(-k)` despite zero additive translation. Closure under
formal but dynamically incompatible branches is sufficient, not necessary.
Tie-value paths require a fixed valid branch; tangent-selector paths need a
paired lift and a common perturbation cone. Raw observations transfer only
when they factor through any preliminary block quotient, or when fibre
oscillation is separately bounded. Raw path realization additionally requires
each declared quotient start atom to meet the image of the allowed raw start
set. In (17.4), `m` counts distinct full affine hyperplanes, including
different offsets of one normal; with `n` discrete control vertices the bound
is multiplied by `n`.

### Corollary 17.1a (observable refinement-growth sandwich)

Suppose the observation in Theorem 17.1 takes values in a finite metric set
whose distinct values are separated by at least `Delta`, and let
`N_t=|P_t|`. For every `0<=epsilon<Delta/2`,

```math
N_T<=C_T(epsilon)<=sum_(t=0)^T N_t,                              \tag{17.7a}
```

where `C_T` is the deterministic predictor complexity defined in
Theorem 17.3 below. The upper predictor is exact. Hence positive exponential
refinement growth forces linear predictive bits, while stabilization gives
one exact finite quotient for all depths.

Indeed, distinct `P_T` atoms have observation trees differing by at least
`Delta`, which proves the packing lower bound. For the upper bound, use states
`(t,[x]_(P_t))`, initialize at level `T`, and update

```math
(t,[x]_(P_t))->(t-1,[F_ex]_(P_(t-1))).                           \tag{17.7b}
```

This is well-defined by (17.1), and decoding the current observation is
exact. This corollary applies to selector labels only when those labels are
themselves declared, margin-separated responses; it does not charge a
numerically invisible tie refinement.

### Theorem 17.1b (compact rational unit-selector lumpability)

Let finitely many control fibres `X_q` be rational polytopes in
`R^r/R1`. Suppose every input map has finitely many closed rational
polyhedral branches covering its source fibre and, on each branch, is

```math
z|->P_sigma z+b,
\qquad (P_sigma z)_i=z_(sigma(i)),                               \tag{17.7c}
```

with rational `b`; overlapping tie branches agree on their overlap and map
into their declared target fibre. Let a finite rational affine sign coloring
(including observations, guards, and branch comparisons) be given, with all
normals perpendicular to `1`.

Restricted affine pullback saturation terminates effectively. Its nonempty
sign atoms form a finite exact forward congruence refining the coloring, and
every declared quotient path is a genuine raw trajectory. In particular,
every bounded rational piecewise-affine terminal response has, after rational
quantization at mesh `eta`, one finite path-realizing predictor with
depth-uniform response error controlled by `eta`.

Indeed, `P_sigma^Ta` is obtained by assigning the `r` labelled coefficients
of `a` to `r` bins and summing within bins, so one seed normal has at most
`r^r` reachable normals. Pullback offsets change by numbers `a.b` from a
finite rational set. All offsets therefore lie in a discrete lattice
`delta Z` after including the finitely many seed offsets. For fixed `q,a`, a
wall meets `X_q` only over an offset interval of length

```math
W_(q,a)=max_(X_q)a.z-min_(X_q)a.z,
```

so at most `1+floor(W_(q,a)/delta)` parallel walls survive (with the evident
one-wall convention for the zero subgroup). Rational linear programs decide
intersection and constant sign. The resulting finite pullback-closed
arrangement satisfies Theorem 17.1. If it has `m_q` walls in projective
dimension `r-1`, its total number of sign faces is at most

```math
sum_q sum_(j=0)^(r-1)2^j binom(m_q,j).                           \tag{17.7d}
```

This is a symbolic/terminal-response theorem. An accumulated real reward
still needs a compatible scalar cocycle; a repeatable reward-cycle defect can
grow linearly. Compactness and discreteness are independent: noncompact
rational translation realizes arbitrarily long transient words without an
infinite orbit, while irrational circle rotation gives dense pullbacks of a
single observation cut on a compact carrier.

### Theorem 17.1c (finite repeatability test for one selector cycle)

Let the return along one fixed legal branch cycle be

```math
Az=P_sigma z+b
```

on `R^r/R1`, and let its closed legal-return domain `D` be projectively
compact. Let `p` be the least common multiple of the directed-cycle lengths
of `sigma`. The following are equivalent:

1. for every `k`, a seed realizes `k` consecutive traversals;
2. one seed remains in `D` under every iterate of `A`;
3. some `y in D` and `lambda in R` satisfy

   ```math
   y,Ay,...,A^(p-1)y in D,
   \qquad A^py=y+lambda1.                                      \tag{17.7e}
   ```

For rational polyhedral `D` and rational selector data, item 3 is one finite
rational LP (not necessarily polynomial-size, since `p` can be exponential).
To prove the equivalence, nested compactness gives 1 implies 2. After the
at-most-`r-1` transient of the functional graph of `sigma`, projective
boundedness forces all selector cycles to have the same mean increment;
hence an appropriate iterate of the infinite seed satisfies (17.7e).
Additive homogeneity makes its first `p` iterates repeat modulo constants,
proving 3 implies 2.

Consequently a nonzero projective observable holonomy on such a cycle proves
that it is not indefinitely repeatable. Pumpable nonzero drift must use a
noncompact legal domain or a setting where arbitrarily long paths cannot be
compactly recoupled to one infinite trajectory.

### Proposition 17.1d (compact equicontinuity does not imply finite prediction)

Let `F_alpha(x)=x+alpha mod 1` be an irrational rotation of the circle and
observe `h(x)=x` in geodesic distance. Then

```math
C_infinity(epsilon)=+infinity
\qquad(0<=epsilon<1/2),                                          \tag{17.7f}
```

while one constant state is valid at error `1/2`. The same impossibility
holds for exact prediction of any nontrivial finite arc coloring, even if
errors are ignored inside a fixed sufficiently thin neighborhood of its cut
points.

Indeed, the output path of any finite autonomous predictor is eventually
periodic. If its eventual period is `p`, then at every fixed phase the true
subsequence advances by `palpha` and is dense. Its supremal distance from the
one decoder value assigned to that phase is `1/2`. For an arc coloring, each
phase subsequence visits robust interiors of two colors. Taking a badly
approximable irrational `alpha` shows that both equicontinuity and
Diophantine regularity are insufficient. In interval coordinates the map is
a compact two-chart unit-selector PWA system; precisely its translations
generate a nondiscrete offset group.

### Theorem 17.1e (finite affine-germ reward dichotomy)

In the rational compact setting of Theorem 17.1b, refine to its exact
sign-atom presentation. For every realizable path from atom `s_0` to atom
`s`, record the projective affine germ

```math
G_w(z)=P_(sigma_w)z+b_w.
```

Only finitely many germs occur. Therefore vertices `(s_0,s,G_w)`, with
composition transitions, form a finite path-realizing germ lift.

Let every raw transition `e` have a rational affine reward `r_e(z)`. On a
lifted edge `v--e-->v'`, define the affine seed response

```math
ell_(v,e)(z)=r_e(G_vz),
\qquad z in C_(s_0),                                            \tag{17.7g}
```

and take its class in `Aff(C_(s_0))/R`, since a scalar finite-state edge toll
may pay any constant. The following are equivalent:

1. scalar tolls on the finite germ lift approximate all cumulative rewards
   with one depth-independent uniform error;
2. every reachable directed cycle in the germ lift has zero total label in
   `Aff(C_(s_0))/R`.

The criterion is effective. Choose an anchor `z_*` in every starting atom
and toll `g_(v,e)=ell_(v,e)(z_*)`. Under item 2, deleting directed cycles
does not change the residual function, so every residual equals that of a
simple lifted path. A valid sharp finite certificate is

```math
B=max_(s_0,z,p simple)
 |sum_(e in p)(ell_e(z)-ell_e(z_*))|,                            \tag{17.7h}
```

computed by finitely many rational LPs. If a lifted cycle has nonconstant
label `L`, choose seeds `z,z'` with `L(z) ne L(z')`. Returning to the same
germ returns both raw states exactly, so `k` repetitions separate their
residuals by `k(L(z)-L(z'))`. No scalar toll machine on this finite control
can have bounded error.

For finiteness, selector parts have at most `r^r` values. Projective
translation differences belong to one rational lattice. A realized germ
satisfies `b_w=G_wz-P_(sigma_w)z`; compact source and target fibres bound the
right side. Hence only finitely many lattice offsets occur.

Exact endpoint telescoping is strictly stronger than bounded error. Affine
potentials `H_s(z)=p_s.z+d_s` and tolls `g_e` exist exactly when

```math
a_e=p_s-P_e^Tp_t,
\qquad c_e-g_e=d_s-d_t-p_t.b_e                                 \tag{17.7i}
```

on every edge. These incidence equations impose equality on transient
coterminal paths as well as recurrent cycles; on lower-dimensional fibres,
the coefficient equality is interpreted modulo affine forms vanishing on the
source fibre.

### Proposition 17.1f (periodic data miss noninvertible mergers)

On `X=[0,1]^2`, the rational unit selectors

```math
A(x,y)=(0,0),
\qquad B(x,y)=(0,x)
```

with rewards `r_A(x,y)=x`, `r_B=0` have zero reward on every periodic orbit
and uniformly bounded cumulative reward. Nevertheless no scalar tolls and
no potential satisfy `r_e-c_e=V-V circ e`. The fixed zero loops force both
tolls to vanish; the `A` equation gives `V(x,y)=V(0,0)+x`, while `B` followed
by the `A` equation gives `V(x,y)=V(0,x)=V(0,0)`. Thus ordinary periodic data
do not characterize exact cohomology for noninvertible selector systems.
The germ lift correctly records the discrepancy only on the transient
coterminal paths `A` and `BA`; because it is not recurrent, it costs a bounded
endpoint error rather than extensive drift.

### Corollary 17.1g (cycle-response packing forces dynamic memory)

Let one input word `w` of length `ell` return every point of a set `Y` to
itself, and let its accumulated reward be `R_w(y)`. Put

```math
m_w(y)={R_w(y)\over ell}.                                       \tag{17.7j}
```

Suppose an `S`-state deterministic simulator, initialized from `y`, predicts
the accumulated reward after every `k` repetitions with

```math
|R_(w^k)(y)-Rhat_(w^k)(y)|
<=epsilon k ell+o_y(k).                                         \tag{17.7k}
```

Then

```math
|S|>=Pack_(>2epsilon)(m_w(Y)).                                  \tag{17.7l}
```

Indeed, two seeds with one initial simulator state produce the same predicted
sequence under `w^k`, whereas their true totals differ by
`kell|m_w(y)-m_w(y')|`. Divide the two error bounds by `kell` and pass to the
limit. Thus a nonconstant affine cycle mean on a connected atom rules out
every finite simulator with bounded absolute error, and its one-dimensional
metric entropy lower-bounds memory at positive asymptotic distortion. This
turns recurrent holonomy into exposed response information rather than only
a yes/no drift certificate.

### Theorem 17.1h (cycle-response rate--distortion law)

Let `Gamma=(V,E)` be a finite path-realizing graph reached from one initial
control, let `Y` be a compact hidden-seed space, and let each lifted edge have
a continuous reward `ell_e:Y->R`. The legal edge path is supplied to an
open-loop weighted transducer, every walk from the initial control is legal
for every seed, and `R_p(y)=sum_(e in p)ell_e(y)`. All seed dependence of a
simulator passes through its finite initial state. Define

```math
d_circ(y,y')=
 max_(reachable directed simple cycles C)
 {1\over|C|}|sum_(e in C)(ell_e(y)-ell_e(y'))|,                  \tag{17.7m}
```

with value zero when no reachable cycle exists. Let `C_epsilon` be the
minimum total number of states in a deterministic simulator
for which one finite `B` gives

```math
|R_p(y)-Rhat_p(y)|<=epsilon|p|+B                                \tag{17.7n}
```

on every legal path. Then

```math
Pack_(>2epsilon)(Y,d_circ)
<=C_epsilon
<=|V|Cov^int_epsilon(Y,d_circ).                                 \tag{17.7o}
```

If the visible graph control is externally supplied and not charged, the
factor `|V|` is removed from the upper bound.

For the lower bound, two seeds with one initial state have identical
predictions along a common access path followed by repetitions of any cycle.
Subtract (17.7n), divide by the repeated length, and pass to the limit. For
the upper bound, retain one internal net center `c` and emit `ell_e(c)` on
edge `e`. Delete directed simple cycles from an arbitrary walk. Each deleted
cycle costs at most `epsilon` times its length, while the remaining simple
path contributes at most

```math
(|V|-1)max_(e,y,c)|ell_e(y)-ell_e(c)|.                           \tag{17.7p}
```

This proves the bound.

Writing

```math
M=max_e sup_(y,y' in Y)|ell_e(y)-ell_e(y')|,
```

the same loop erasure gives the stronger pairwise estimate

```math
|R_p(y)-R_p(y')|
<=|p|d_circ(y,y')+(|V|-1)M.                                    \tag{17.7q}
```

Thus approximate compositional complexity is not a formal product of static
entropy and graph memory in this class. It is the metric entropy of the
**recurrent cycle-response image**, with the finite control congruence carried
separately. Transient merging diamonds affect only `B`; recurrent holonomy is
the information amplified to macroscopic scale.

This is a theorem for cumulative path reward, not simultaneous prediction of
every individual edge reward. Opposite edge labels can cancel on a cycle;
per-edge queries require their own stronger response metric. Unreachable
components are omitted throughout. A syntactic control cycle which does not
return the affine germ is likewise not a cycle in this theorem.

### Theorem 17.1i (finite invariant-grid shadowing)

Let finitely many deterministic continuous PWA maps act from one whole
compact convex rational projective polytope `X_q` to another. Suppose every
affine branch has a rational unit-selector linear part and rational
translation. Then for every
`eta>0` there are finite internal nets `C_(q,eta) subset X_q` such that every
map sends its source net exactly into its target net. Moreover,

```math
sum_q|C_(q,eta)|
<=sum_q O_q(eta^(-dim X_q)),                                    \tag{17.7r}
```

and a raw state initialized within Hilbert distance `eta` of a net point is
shadowed by that point under every switching word with the same error
`eta`. Thus every `L`-Lipschitz terminal projective response has one finite
path-realizing simulator of error `Leta` at arbitrary depth, without a
contraction assumption.

To prove this, use the common ambient coordinate-difference lattice in a
fixed coordinate gauge and clear all vertex and projective-translation
denominators by `D`. On mesh `h=1/(DN)`, the regauged selector matrices are
integral, so selector maps preserve the ambient grid. Their intersection
with each invariant polytope is therefore forward invariant. A rational
triangulation and denominator-`N` barycentric rounding
give an `O_q(1/N)` internal net and the lattice-point count (17.7r), including
on lower-dimensional rational faces.

Each selector derivative contracts Hilbert norm. Continuity across the
finite PWA cells and subdivision of the segment between two points therefore
make the whole map nonexpansive. Exact evolution of the chosen grid point
introduces no fresh residual, proving depth-uniform shadowing.

This is not generic finite quantization. Independently chosen intrinsic
lattices on lower-dimensional fibres need not be compatible. State-dependent
edge enabling must first be refined into whole-source controls. Irrational
translations preserve no such grid, discontinuous selector switches need not
be nonexpansive, and accumulated affine rewards still require the
recurrent-germ criterion of Theorem 17.1e. Even the rational contraction
`x|->x/2` on `[0,1]` has no finite forward-invariant internal net except
`{0}`: every positive point has an infinite halving orbit. Unit/integer
lattice preservation, not rationality alone, is essential.

In particular, if a finite family of max-plus maps has a common eigenprofile

```math
F_e(v)=v+lambda_e1,
```

then one projective state with per-letter toll `lambda_e` shadows every word.
Any initial projective radius `R` from `[v]` remains at most `R` forever. This
is an implicit exposed carrier certified generator by generator, rather than
an enumeration of active germs.

If transition rewards `ell_e` are `L_e`-Lipschitz, the same invariant-grid
simulator predicts accumulated reward, not only terminal response:

```math
|R_w(x)-R_w(c)|
<=eta sum_(e in w)L_e
<=eta L|w|.                                                     \tag{17.7r'}
```

Thus per-step distortion `epsilon` costs at most
`sum_q O_q((L/epsilon)^(dim X_q))` states. This does not give bounded absolute
error; it gives a depth-uniform *rate*. Theorem 17.1j shows that exponential
dependence on growing projective dimension cannot be removed in general.

### Theorem 17.1j (fixed-alphabet exponential exposed-germ memory)

For every `r>=2`, put `k=floor(r/2)` and let

```math
X_r={[z] in R^r/R1:osc(z)<=1},
qquad
Y_r={1_A:A subseteq [r], |A|=k}.                               \tag{17.7s}
```

There is a three-letter rational max-plus selector system on `X_r` whose
exact cumulative-response complexity on the orbit `Y_r` is

```math
C_epsilon=
\begin{cases}
\binom r{\lfloor r/2\rfloor},&0<=epsilon<1/2,\\
1,&epsilon>=1/2.
\end{cases}                                                     \tag{17.7t}
```

Here complexity means the least total number of deterministic predictor
states satisfying `|R_w-Rhat_w|<=epsilon|w|+B` for every word, from one fixed
raw initial seed, with arbitrary finite `B`.

For a permutation `sigma`, take the all-finite matrix

```math
K^sigma_(ij)=0 if j=sigma(i), and -2 otherwise.                 \tag{17.7u}
```

On `X_r` its intended selector is uniquely maximizing with margin at least
one, so it realizes `P_sigma` globally. Two letters use the cycle
`(1 2 ... r)` and transposition `(1 2)`, both with zero reward. They generate
`S_r`. The third letter is the identity selector with projectively invariant
reward

```math
rho([z])=z_1-r^(-1)sum_j z_j.                                  \tag{17.7v}
```

The orbit of one weight-`k` indicator is all of `Y_r`. For distinct
`1_A,1_B`, choose `j in A triangle B` and a common suffix permutation pulling
coordinate `j` to coordinate one. The repeatable probe loop then has reward
rates differing by exactly one. If histories reaching `A` and `B` ended in
one predictor state, their predicted future probe totals could differ only
by a fixed prefix constant. Repeating the probe `n` times and dividing the
two error inequalities by `n` would give `1<=2epsilon`. Thus all
`binom(r,k)` orbit states are necessary below error `1/2`. Storing the current
subset is exact. At error `1/2`, one state emits `1/2-k/r` on the probe and
zero otherwise, proving (17.7t).

Consequently a fixed input alphabet, compactness, rationality, globally
unique active selectors, exact isometry, and zero translation holonomy do
not imply a subexponential exposed carrier:

```math
log_2 C_epsilon
>=r-(1/2)log_2 r-O(1).                                         \tag{17.7w}
```

The full germ group has `r!` elements, but the exact response quotient is the
orbit `S_r/(S_k times S_(r-k))`; stabilizer germs are correctly discarded.
Thus the lower bound charges exposed response information, not syntactic
germ enumeration.

### Proposition 17.1k (one selector has subexponential recurrent germ count)

Let one globally active affine unit selector

```math
A[z]=[P_sigma z+b]
```

preserve a nonempty projectively compact set on `r` coordinates. If `k_0`
is the longest tail in the functional graph of `sigma`, and `p` is the least
common multiple of all its directed-cycle lengths, then

```math
[A^(n+p)]=[A^n]
\qquad(n>=k_0).                                                 \tag{17.7x}
```

Hence it has at most `k_0+p<=r-1+p` iterate germs and at most `p` recurrent
germs. In particular

```math
p<=g(r)=exp((1+o(1))sqrt(r log r)),                             \tag{17.7y}
```

where `g` is Landau's maximal permutation order.

Indeed, after `k_0` steps every coordinate path lies on a selector cycle.
The difference between iterates `n+p` and `n` is `p` times that cycle's mean
translation. Two unequal cycle means would make a coordinate difference
grow linearly along one orbit, contradicting projective compactness. All
means are therefore equal, and the difference is scalar. Theorem 17.1j
shows the sharp structural boundary: two switched permutation selectors,
not one, already expose exponentially many recurrent classes. This does not
settle one genuinely PWA selector generator whose active selector changes
with the cell; Theorem 17.1o below does.

### Theorem 17.1l (cycle LP for a proposed reward congruence)

Let `Q` be a finite deterministic reward system with alphabet `E`, transition
`delta_e`, and edge reward `r(q,e)`. Let `pi:Q->S` be an input congruence, so
`pi(delta_e q)` depends only on `pi(q)`. For scalar quotient tolls `g(s,e)`,
write `R_w(q)` and `Rhat_w(pi(q))` for the two cumulative rewards. Define

```math
D(pi)=inf_g inf{epsilon>=0: exists B<infinity such that
 |R_w(q)-Rhat_w(pi(q))|<=epsilon|w|+B for all q,w}.              \tag{17.7z}
```

Then, with `C(Q)` the finite set of directed simple cycles in the raw
transition graph,

```math
D(pi)=inf_g max_(C in C(Q))
 {1\over|C|}|sum_((q,e) in C)(r(q,e)-g(pi(q),e))|.              \tag{17.7aa}
```

Thus the best asymptotic distortion of a fixed finite congruence is one
finite linear program. In particular, it has a depth-independent absolute
error iff tolls can make every raw cycle defect zero. Necessity follows by
repeating a raw cycle. For sufficiency, erase directed cycles from an
arbitrary path. Each erased cycle costs at most `epsilon` per edge, while the
remaining simple path has fewer than `|Q|` edges and contributes one uniform
constant. Taking the infimum over tolls proves (17.7aa).

Feasible bounded-error congruences need not have a unique coarsest member and
are not closed under join. On states `{I,A,B}`, let every `a` go to `A`, every
`b` go to `B`, and let the only nonzero reward be `r(A,b)=1`. Both

```math
{I,A}|{B}
\quad\hbox{and}\quad
{A}|{I,B}                                                       \tag{17.7ab}
```

are bounded-error congruences: toll the recurrent `A --b--> B` edge by one
and all other recurrent edges by zero; the transient state costs at most
one. Their join is the one-block partition. Constant words force its two
letter tolls to be zero, whereas `(ba)^n` from `A` has reward `n`. The join
therefore fails. Moreover every pair of raw starting states has same-word
reward difference at most one, since the first letter synchronizes them.
Pairwise asymptotic response distance can consequently vanish while every
one-state compositional quotient fails: the missing datum is cycle incidence
created by the quotient, not another static metric.

### Theorem 17.1m (finite projective-semigroup cycle realization)

Let two finite alphabets of all-finite max-plus matrices `A_a,B_a` have
finite projective product semigroups `P_A,P_B`. Normalize a matrix `P` by
subtracting `nu(P)=max_(ij)P_(ij)`, and write

```math
P A_a=kappa_A(P,a)+tau_A(P,a)                                  \tag{17.7ac}
```

with `tau_A` normalized; define the analogous `B` transition. Synchronize
the two projective Cayley graphs under the same input letter and label an
edge by `kappa_A-kappa_B`. On the reachable synchronized graph `G`, let
`mu_max(G),mu_min(G)` be its extreme directed-cycle means. Then

```math
sup_(w ne empty){rho(A_w)-rho(B_w)\over|w|}=mu_max(G),
\qquad
inf_(w ne empty){rho(A_w)-rho(B_w)\over|w|}=mu_min(G).           \tag{17.7ad}
```

In particular the aligned-word recurrent response distance is the larger
absolute value of these two means.

For one direction, a carrier cycle labelled by `v` gives
`P A_v=kappa+P`. Every row of the all-finite `P` is a finite left eigenvector
of the irreducible all-finite `A_v`, so `kappa=rho(A_v)`; similarly for `B`.
Conversely, fix any word `w` and follow the synchronized sequence of
projective powers `([A_(w^j)],[B_(w^j)])`. Finiteness gives `i>=1,p>=1` for
which the pair at powers `i` and `i+p` agrees. The intervening closed walk is
labelled by `w^p`. Its two finite-row eigenvector equations make its tolls
`rho(A_(w^p))=p rho(A_w)` and `rho(B_(w^p))=p rho(B_w)`, so its mean is exactly
`(rho(A_w)-rho(B_w))/|w|`. No claim is made that the first copy of `w` closes
or that its max-entry normalization toll already equals its spectral radius.

The theorem gives a genuine all-word path-realizing carrier and reduces its
response comparison to two finite maximum-cycle-mean computations. It is
output-sensitive rather than automatically compressive: the projective
semigroups may be exponential or infinite. A fixed-letter critical graph is
not a substitute. For one-state systems with rewards
`A_a=1,A_b=0` and `B_a=0,B_b=1`, the two letterwise envelopes have the same
critical graph, but aligned word responses differ by `+1` on `a` and `-1`
on `b`.

If one all-finite alphabet is pairwise max-plus commuting, its exact carrier
collapses further. The common-eigenvector theorem of
[Katz--Schneider--Sergeev](https://arxiv.org/abs/1005.1424) supplies a finite
`h` with

```math
A_a h=rho(A_a)+h
```

for every letter. Therefore

```math
rho(A_w)=sum_a |w|_a rho(A_a).                                 \tag{17.7ae}
```

The recurrent response factors through the Parikh vector, and a scalar
per-letter toll replaces the projective Cayley state. This is a checkable
algebraic synchronization mechanism, not an assumption about common active
faces.

### Theorem 17.1n (dominating quotient with a coherent path lift)

Let `T_e` be finite max-plus matrices on raw states `I`, let `S_e` be finite
max-plus matrices on coarse states `J`, and let `pi:I->J` be onto. Suppose
nonnegative numbers `eta_e^+,eta_e^-` satisfy

```math
T_e(i,j)<=S_e(pi(i),pi(j))+eta_e^+                 for all i,j,
```

and, for every `i in pi^(-1)(a)` and every coarse target `b`, there is some
`j in pi^(-1)(b)` such that

```math
T_e(i,j)>=S_e(a,b)-eta_e^-.                                    \tag{17.7af}
```

For every word `w=e_1...e_t`,

```math
T_w(i,j)<=S_w(pi(i),pi(j))+sum_s eta_(e_s)^+,
max_(j in pi^(-1)(b))T_w(i,j)
 >=S_w(a,b)-sum_s eta_(e_s)^-  (i in pi^(-1)(a)).               \tag{17.7ag}
```

Consequently

```math
-sum_s eta_(e_s)^-
<=rho(T_w)-rho(S_w)
<=sum_s eta_(e_s)^+.                                            \tag{17.7ah}
```

The proof maps every raw path to its coarse path for the upper bound. For the
lower bound, choose a maximizing coarse path and lift it inductively, allowing
the raw representative of each coarse state to depend on the preceding path.
This proves (17.7ag). A raw cycle projects to a coarse closed path. Conversely,
repeat a critical coarse cycle and lift it forever. At returns to one coarse
phase, finiteness of its raw fibre makes two raw representatives repeat; the
intervening lift is a genuine raw cycle with the claimed average. This proves
(17.7ah). At zero defect, the generator inequalities certify
`rho(T_w)=rho(S_w)` for **every** aligned word without enumerating products or
active cells. The lift relation is the path-realization datum; domination is
the response-approximation datum.

This applies beyond max-plus notation to any finite additive path dynamic
program. It explains the corrected four-to-two weighted-automaton benchmark:
the even representatives form an exact section and every microscopic edge is
dominated by its corrected block edge. A quotient which only matches
letterwise critical graphs lacks the coherent lift required by (17.7af). A
fixed coherent section is a convenient stronger certificate, not a necessary
hypothesis: it satisfies (17.7af) by choosing `j=iota(b)` at every step.
For a transitive permutation action, the gap kernel with weight zero on
`i->g_ei` and `-C` elsewhere has an exact one-state relational lift, choosing
`j=g_ei`; a fixed section would require a nonempty invariant singleton and
therefore generally fails. The two-state flip is the antiferromagnetic Ising
parity mechanism. Allowing the witness to move is essential.

### Theorem 17.1o (one lattice-PWA generator has exponential exposed memory)

For `m>=1`, put `r=2m` and use dual-rail coordinates
`(u_0,...,u_(m-1),v_0,...,v_(m-1))` on

```math
X={ [u,v] in R^(2m)/R1:osc(u,v)<=1}.
```

There is one continuous rational PWA unit-selector map `F:X->X`, described by
`O(m)` shared min/max gates, with a genuine orbit of period `2^m`. Define

```math
u'_0=v_0,\qquad v'_0=u_0,
```

and, for `i>=1`,

```math
c_i=min_(j<i)u_j,\qquad d_i=max_(j<i)v_j,

u'_i=max{min(u_i,d_i),min(v_i,c_i)},
v'_i=max{min(v_i,d_i),min(u_i,c_i)}.                            \tag{17.7ai}
```

Every component is a continuous additively homogeneous lattice polynomial.
On each cell it selects one input coordinate, and every output lies between
the input minimum and maximum, so `F` is a rational unit-selector PWA map
preserving `X`.

For `b in {0,...,2^m-1}`, let `y(b)` have

```math
u_i=bit_i(b),\qquad v_i=1-bit_i(b).
```

Then `c_i` is exactly the carry into bit `i` when adding one and
`d_i=1-c_i`; (17.7ai) gives

```math
F(y(b))=y(b+1 mod 2^m).                                       \tag{17.7aj}
```

Add one identity input `p` with bounded projective reward

```math
rho([u,v])={u_(m-1)-v_(m-1)\over2}.                            \tag{17.7ak}
```

The probe word around the counter orbit is one half `-1/2` and one half
`+1/2`, with primitive cyclic period `2^m`. Thus for every two phases a
common power of `F` makes their probe rewards differ by one. Repeating `p`
and using the same merged-state argument as Theorem 17.1j proves

```math
C_epsilon=
\begin{cases}
2^m,&0<=epsilon<1/2,\\
1,&epsilon>=1/2.
\end{cases}                                                     \tag{17.7al}
```

The upper predictor stores the counter phase; above threshold it emits zero.
This is path-realized exponential memory from one polynomial-description PWA
generator plus one non-evolving scalar query, not an enumeration of `2^m`
active cells.

The query is essential. With only the letter `F` and any bounded one-step
reward, phase trajectories are rotations of one finite cycle. Their
cumulative totals differ by a bounded amount, and one state emitting the
cycle mean has bounded absolute error. Path complexity becomes response
information only when the declared future can freeze and repeatedly expose a
phase.

This also separates lattice-PWA from max-plus-linear dynamics. For one
all-finite max-plus matrix, the classical cyclicity theorem makes projective
powers ultimately periodic with the cyclicity of the critical graph, at most
the Landau lcm bound `g(r)` (and at most `r` for one critical component); see
[Sergeev--Schneider](https://arxiv.org/abs/0912.2534). The min/max carry
circuit in (17.7ai) lies outside coordinatewise-convex max-plus linearity.

The exponential base can be arbitrarily close to two. Fix an even block
length `ell`, let `U` be the `q=binom(ell,ell/2)` constant-weight binary
words, and order them cyclically with words omitting coordinate one followed
by words containing it. For `u^a in U`, put

```math
h_a(z)=min_(j:u^a_j=1)z_j,
\qquad
S(z)_j=max_(a:u^(a+1)_j=1)h_a(z).                              \tag{17.7am0}
```

Equal support sizes give `h_a(u^b)=1_(a=b)`, hence `S(u^a)=u^(a+1)`.
The detector `h_(q-1)` and the maximum over coordinates outside `u^(q-1)`
are complementary carry rails on `U`. More explicitly, for blocks
`z^(0),...,z^(m-1)`, let

```math
c_i=min_(j<i)h_(q-1)(z^(j)),
\qquad
d_i=max_(j<i)max_(k:u^(q-1)_k=0)z^(j)_k,
```

with `c_0=1,d_0=0`, and set

```math
z'^(i)_k=max{min(z^(i)_k,d_i),min(S(z^(i))_k,c_i)}.             \tag{17.7am0a}
```

On `U^m`, `c_i` is one exactly when every lower digit is `q-1`, and
`d_i=1-c_i`.  Thus (17.7am0a) leaves digit `i` fixed unless there is a carry
and otherwise applies `S`; it is one `O_ell(m)`-gate lattice-PWA base-`q`
counter of period `q^m`. The most-significant coordinate-one probe again has
one negative half and one positive half. Globally it is the projectively
invariant functional

```math
p(z)=z_1-{1\over ell}\sum_(j=1)^ell z_j                         \tag{17.7am0b}
```

on the most-significant block, equal to `-1/2,+1/2` on the
constant-weight digit orbit. The exact sub-half-error complexity is `q^m`.

Since

```math
binom(ell,ell/2)^(1/ell)->2,                                   \tag{17.7am1}
```

for every fixed `c<2` there is a linear-description one-generator family on
`r=m ell` coordinates requiring more than `c^r` states. This strengthening
uses a fixed local table once `c` is chosen; it does not enumerate the global
orbit.

Allowing the local table to grow still keeps the presentation polynomial and
makes the exponent asymptotically optimal. Choose `m=2^ell`. Computing all
local `h_a` once per block uses `O(q ell)` shared lattice gates, so the whole
map has

```math
O(m q ell)=O(r^2/ell^(3/2))
```

gates. Since `r=ell 2^ell` and
`log_2q=ell-O(log ell)`, its exposed orbit has

```math
q^m=2^(r-O(r log ell/ell))=2^(r-o(r))                           \tag{17.7am2}
```

states. Thus even polynomial circuit description does not prevent
essentially maximal Boolean dynamic response information.

### Theorem 17.1p (minimum approximate reward congruence is NP-complete)

Given a rational finite deterministic reward system, rational `epsilon`, and
`k`, deciding whether an input congruence `pi` with at most `k` classes has
`D(pi)<=epsilon` in Theorem 17.1l is NP-complete. This already holds for

```math
k=3,\qquad epsilon=1/2,
```

rewards in `{-1,0,1}`, and identity transitions.

For identity transitions every partition is a congruence, every raw state is
fixed, and independent midpoint optimization of each block/letter toll gives

```math
D(pi)={1\over2}max_(B,e)
 \left(max_(q in B)r(q,e)-min_(q in B)r(q,e)\right).             \tag{17.7am}
```

Given a graph `G=(V,E)`, use states `V` and one letter `e_(uv)` for each
edge, with reward `+1` at `u`, `-1` at `v`, and zero elsewhere. A block has
diameter at most one in every reward coordinate iff it contains no edge.
Thus `D(pi)<=1/2` iff every block is independent, and the minimum state count
is exactly the chromatic number of `G`. Three-colorability proves hardness.

The problem is in NP. Guess `pi`, check the congruence, and introduce quotient
tolls `g` and two vertex potentials. The polynomial system

```math
r(q,e)-g(pi(q),e)-epsilon
 <=p^+(delta_e q)-p^+(q),

-r(q,e)+g(pi(q),e)-epsilon
 <=p^-(delta_e q)-p^-(q)                                      \tag{17.7an}
```

is feasible exactly when the positive and negative defect graphs have no
cycle of mean above `epsilon`. This is equivalent to `D(pi)<=epsilon` and is
rational LP feasibility.

Hardness persists under maximal dynamic forgetting. Replace each graph edge
`uv` by letters `U_(uv),V_(uv)` resetting every state respectively to `u,v`,
with reward tables

```math
r(q,U_(uv))=1_(q=u)-1_(q=v),
\qquad
r(q,V_(uv))=1_(q=v)-1_(q=u).                                  \tag{17.7an1}
```

Independent blocks admit edgewise defect at most `1/2`. If `u,v` merge, the
two endpoint self-loops force the two corresponding quotient tolls to be at
least `1/2`. The raw two-cycle `v -> u -> v` has reward `-2` but quotient
toll at least one, hence absolute mean defect at least `3/2`. Thus the optimum
is again `chi(G)` although every generator is a rank-one reset.

This is an algorithmic obstruction, not an information lower bound on every
instance. With identity dynamics and one scalar reward coordinate, sorting
and greedily grouping intervals of diameter `2epsilon` is optimal. The
hardness is created by a growing family of jointly declared future queries.

### Proposition 17.1q (exact word responses need not admit a small path lift)

There is a useful exact relational characterization behind the separation.
Let every entry of `T_e` lie in `{0,-C}` and let `R_e` be its zero-edge
relation.  For a word `w`, write `R_w` for relational composition in word
order.  Then

```math
rho(T_w)=0
\quad\Longleftrightarrow\quad
R_w\hbox{ contains a directed cycle}.                           \tag{17.7ao0}
```

Indeed, the zero entries of `T_w` are exactly `R_w`; since all entries are
nonpositive, spectral radius zero is equivalent to a zero-weight directed
cycle.  By contrast, an exact Theorem 17.1n lift to the scalar zero system
exists exactly when every `R_e` is left-total.  If (17.7ao0) fails for a
word, every simple cycle of `T_w` contains a negative edge, so

```math
rho(T_w)<=-C/|I|,
\qquad rho(T_(w^k))=k rho(T_w).                                 \tag{17.7ao1}
```

Thus failure has a finite, pumpable linear-drift witness.  What can fail is
the converse: periodic completeness of every word relation need not supply
one locally total relation.

For every `r>=2` and `C>0`, let raw states and letters both be `[r]`, and set

```math
T_e(i,j)=0 if i=e, and -C otherwise                             \tag{17.7ao}
```

for every target `j`. Then `rho(T_w)=0` for every nonempty word `w`: the
closed path following the successive letter names has weight zero, and every
edge has weight at most zero. Thus the one-state alphabet `S_e=[0]` has
exactly the same aligned-word spectral responses.

Nevertheless every Theorem 17.1n certificate of local errors satisfying

```math
eta_e^++eta_e^-<C\qquad(e in[r])                                \tag{17.7ap}
```

has at least `r` coarse states. Indeed, for a fibre map `pi` define

```math
m_e(i,b)=max_(j in pi^(-1)(b))T_e(i,j).
```

The upper and lower path-lift inequalities put `m_e(i,b)` within
`[S_e(pi(i),b)-eta_e^-,S_e(pi(i),b)+eta_e^+]`. If distinct `i,i'` share one
fibre, letter `e=i` gives `m_e(i,b)=0` and `m_e(i',b)=-C` for every `b`,
contradicting (17.7ap).

The obstruction survives every common diagonal gauge. After replacing
`T_e(i,j)` by `T_e(i,j)+h_j-h_i`, the two relevant block-row gaps for letters
`i` and `i'` are

```math
|C+h_(i')-h_i|\quad\hbox{and}\quad|C-h_(i')+h_i|,
```

whose maximum is at least `C`.

The alphabet can be fixed while the gap becomes exponential. Let
`I_m={0,1}^m`, use letters `e in {0,1}`, and put a zero edge from `s` to `t`
exactly when

```math
s_1=e,\qquad (t_1,...,t_(m-1))=(s_2,...,s_m),                  \tag{17.7aq}
```

with `t_m` free; give every other edge weight `-C`. Every word has a closed
zero lift given by the length-`m` windows of its periodic extension, so the
scalar zero system is again exact. But the length-`m` word naming `s` has a
zero row only at `s`. Iterating a relational certificate shows that merging
two states requires `m max_e(eta_e^++eta_e^-)>=C`. Hence exact path lifting
needs all `2^m` states; taking `C=m` gives the same conclusion below a fixed
unit local defect. Every raw state is genuinely exposed by its own word.

This same family also rules out contraction or unique criticality as the
missing synchronization hypothesis. If `w` has length `L>=2m` and

```math
s_w=(w_1,...,w_m),
```

then direct window propagation gives

```math
T_w(s,t)=
\begin{cases}
0,&s=s_w,\\
-C,&s\ne s_w
\end{cases}
\qquad(t\in I_m).                                               \tag{17.7aq1}
```

The final `m` appended bits can be chosen freely. From a wrong initial
window, one `-C` edge resets the window and the remaining at least `m` steps
both synchronize it to the word and choose the endpoint. Hence every such
word matrix has max-plus projective rank one and Hilbert contraction
coefficient zero. Moreover its unique critical node is `s_w`: repeating
`w` forces the unique length-`m` window of `w^infinity` at the word boundary.
Every competing simple cycle has mean at most `-C/2`, a uniform exposed
margin.

Nevertheless, over the blocked alphabet of all length-`L` words, any
Theorem 17.1n quotient merging `s\ne t` has block defect at least `C`: choose
a word beginning with `s` and compare the two constant rows in (17.7aq1).
The scalar model attains this threshold with `(eta^+,eta^-)=(0,C)`. Thus even
exact one-state spectra, a unique uniformly exposed critical node for every
product, and one-block contraction zero coexist with `2^m` path-lift states.
If a letterwise certificate has defect at most `delta`, its induced block
defect is at most `L delta`; at the natural horizon `L=2m,C=m`, every
`delta<1/2` certificate is therefore injective. Longer blocking can amortize
this finite synchronization toll, so the obstruction is to fixed local
accuracy, not to vanishing error per original step as `L/m->infinity`.

Thus Theorem 17.1n is a strict structural certificate, not a converse to
all-word spectral agreement. Spectral optimization has quantifiers
`for every word, there exists a critical raw cycle`; a uniform path lift asks
one relation to continue from every represented raw state. Exchanging these
quantifiers can cost the entire microscopic state space even when the scalar
response function is already exact. Contraction forgets the *incoming
vector* after the future word is known; it does not choose one microscopic
seed coherently before that word is revealed.

### Theorem 17.1r (finite survival lumpability for scalar tropical response)

Let `T_e` be finite real max-plus matrices on an `r`-point set `I`. Fix
scalar tolls `lambda_e` and errors `theta_e>=0`, assume

```math
T_e(i,j)<=lambda_e,
```

and let `R_e` contain the **good edges**

```math
(i,j) in R_e
\quad\Longleftrightarrow\quad
T_e(i,j)>=lambda_e-theta_e.                                    \tag{17.7ar}
```

A finite survival carrier consists of a finite deterministic input system
`(Q,q_0,delta)` and nonempty subsets `K_q subseteq I` such that

```math
K_(delta(q,e)) subseteq Delta_e(K_q),
\qquad
Delta_e(S)={j:exists i in S\ (i,j) in R_e}.                     \tag{17.7as}
```

Then every word `w=e_1...e_t` satisfies

```math
sum_s(lambda_(e_s)-theta_(e_s))
 <=rho(T_w)<=sum_s lambda_(e_s).                                \tag{17.7at}
```

The carrier is anticipatory rather than rowwise: after a finite future is
declared, choose any terminal point in its nonempty `K`-set and use (17.7as)
backwards to obtain a genuine good-edge path. Apply this to `w^r`. The
resulting `r` macro-edges contain a repeated raw vertex, hence a cycle of
`T_w` whose every macro-edge has weight at least
`sum(lambda_e-theta_e)`. Domination gives the upper bound. Thus symbolic path
existence, scalar response approximation, and finite control state are
separate and explicit; no active-cell or relation-semigroup enumeration is
part of the certificate.

At zero error this criterion is complete. Put
`R_e={(i,j):T_e(i,j)=lambda_e}`. The following are equivalent:

1. `rho(T_w)=sum_(e in w)lambda_e` for every nonempty word;
2. every `R_w` contains a directed cycle;
3. `R_u` is nonempty for every word `u`;
4. the subset trajectory `Delta_u(I)` never reaches the empty set;
5. a finite survival carrier exists.

Only `3=>2` is not immediate: `R_(w^r)` is nonempty, so the `r`-vertex graph
`R_w` has a path of `r` edges and therefore a directed cycle. The canonical
carrier consists of the reachable nonempty endpoint subsets

```math
Q={Delta_u(I):u in E^*},
\qquad |Q|<=2^r-1,                                             \tag{17.7au}
```

with `K_S=S`; smaller certificates may use proper inclusions in (17.7as).
This is exponentially smaller than the `2^(r^2)` relation semigroup in the
worst ambient count, and it retains strictly less information than all raw
rows.

The powerset scale is nevertheless sharp for *monitoring* survival, even
with two letters. On `I=Z/rZ`, let `P` be cyclic permutation and let `D` be
the partial identity undefined only at zero. Their subset actions are

```math
P(S)=S+1,
\qquad D(S)=S\setminus\{0\}.                                   \tag{17.7au1}
```

Conjugating `D` by powers of `P` deletes any chosen point, so all `2^r`
subsets are reachable from `I`. If `S\ne T`, preserve one point of
`S triangle T` and delete all others; the suffix survives from exactly one
of `S,T`. Hence the exact alive/dead monitor has `2^r` Myhill--Nerode states.
This is a verification lower bound, not a claim that every positive instance
needs an exponential survival certificate.

Failure has a finite pumpable witness. A shortest word `u` with
`Delta_u(I)=emptyset` has

```math
|u|<=2^r-1.                                                     \tag{17.7av}
```

If every nongood edge has the uniform gap
`T_e(i,j)<=lambda_e-C`, then

```math
T_(u^k)(i,j)<=k sum_(e in u)lambda_e-kC
\quad\hbox{for all }i,j,k.                                     \tag{17.7aw}
```

Indeed a shortest subset trajectory cannot repeat a nonempty subset, and
every `u`-block path uses a nongood edge. Repetition is therefore a genuine
linear response-drift witness.

The hierarchy is strict. An exact scalar Theorem 17.1n path lift requires
every `R_e` to be left-total. Survival lumpability only requires nonmortality
of the endpoint subset dynamics. In the binary de Bruijn family,
`Delta_e(I)=I` for both letters, so `Q={I}` is already an exact one-state
survival carrier, while every exact rowwise path lift has `2^m` states. This
is a finite tropical lumpability theorem tailored to the unrooted response
query, not a weakened simulation of the entire landscape.

### Theorem 17.1s (anticipatory-support lumpability over a coarse system)

Let `T_e` be raw max-plus matrices on a finite nonempty set `I`, let `S_e` be
coarse max-plus matrices on a finite nonempty set `J`, and let `pi:I->J` be
onto. All errors and gauge coordinates below are finite real numbers. A
common raw gauge replaces

```math
T_e(i,j)\quad\hbox{by}\quad
That_e(i,j)=T_e(i,j)+h_i-h_j
```

without changing any word spectral radius. Suppose

```math
That_e(i,j)<=S_e(pi(i),pi(j))+alpha_e.                          \tag{17.7ax}
```

For each `a in J`, choose a finite nonempty family `K_a` of nonempty subsets
of `pi^(-1)(a)`. For every `K in K_a` and every finite coarse edge
`S_e(a,b)`, choose a successor support

```math
K'=sigma(K,e,b) in K_b
```

and a shortfall `d(K,e,b)>=0` such that

```math
K' subseteq {j in pi^(-1)(b):exists i in K,
              That_e(i,j)>=S_e(a,b)-d(K,e,b)}.                 \tag{17.7ay}
```

Assume there are a potential `psi` on all selected supports and numbers
`beta_e>=0` with

```math
d(K,e,b)<=beta_e+psi(K')-psi(K).                               \tag{17.7az}
```

Then every word `w=e_1...e_t` with finite `rho(S_w)` satisfies

```math
-sum_s beta_(e_s)
 <=rho(T_w)-rho(S_w)
<=sum_s alpha_(e_s).                                          \tag{17.7ba}
```

If `rho(S_w)=-infinity`, (17.7ax) forces `rho(T_w)=-infinity` as
well.

For the lower bound, take a critical cycle of the product `S_w`, unfold each
of its macro-edges into a maximizing generator-level coarse path, and repeat
the resulting periodic path. Evolve the finite support state using `sigma`.
At returns to one coarse phase, some support `K` repeats. Composing (17.7ay)
around the intervening repetitions gives `K subseteq Delta(K)`. The induced
finite relation on `K` has positive indegree at every vertex and hence a
directed cycle. Its raw weight is the critical coarse weight minus the sum of
the shortfalls. The potential in (17.7az) telescopes on the closed support
path, leaving at most `sum beta_e` per copy of `w`. Since
`rho(T_(w^p))=p rho(T_w)`, this proves the lower bound. Projection of every
raw cycle through (17.7ax) proves the upper bound.

The certificate separates four resources:

1. `sigma` gives symbolic coarse-path realization;
2. (17.7ay) gives backward-surjective endpoint realization;
3. (17.7az) controls response error modulo a finite cocycle;
4. `sum_a |K_a|` is the support-state information complexity.

All conditions are generator-checkable; (17.7az) is a finite system of
difference constraints, equivalently a cycle-mean feasibility problem. A
canonical verifier, when all coarse paths have a near-tight support lift, has
at most

```math
sum_(a in J)(2^|pi^(-1)(a)|-1)                                 \tag{17.7bb}
```

support states. This certificate is stronger than the nonlocal necessary and
sufficient statement “every word has some critical coarse cycle with a tight
raw cyclic lift.” Source-total rowwise lifting and target-surjective support
lifting are otherwise incomparable. The toggle and width-two examples below
make support lifting smaller, while the deterministic de Bruijn shift in
Section 18 has a one-state rowwise lift but needs exponentially many support
states.

Strictness already occurs in a two-state coarse system. Fix `r>=2`, use
alphabet `E=[r]`, let every letter have coarse weight zero on `0<->1` and
weight `-C` on the two fixed-state edges, take raw states
`I={0,1} times [r]`, and put

```math
T_e((a,k),(b,l))=
\begin{cases}
0,&b=1-a\hbox{ and }k=e,\\
-C,&\hbox{otherwise}.
\end{cases}                                                     \tag{17.7bc}
```

The two supports `K_a={a} times [r]` lift both types of coarse edge exactly.
For the zero toggle edge, the single source `(a,e)` reaches every target in
the opposite support; every fixed-state raw edge has the coarse weight
`-C`. Theorem 17.1s therefore gives exact equality of every word spectrum
with only two support states. But a
row `(a,k)` with `k\ne e` has no tight `e`-successor, so no corresponding
rowwise lift through this displayed `pi,S` exists.

There is also a pumpable flat-toll converse. Suppose

```math
S_e(a,b)=lambda_e\quad\hbox{on an allowed relation, and }-infinity
\quad\hbox{otherwise},
```

where `lambda_e` is finite and `gamma>0`. Call an allowed projected raw edge
tight when its weight is exactly `lambda_e`, and assume every other finite
raw edge is at most `lambda_e-gamma`. For a word
whose coarse relation has a cycle,

```math
rho(T_w)=sum_(e in w)lambda_e
```

iff its tight raw word relation has a directed cycle. If it has none, every
simple raw cycle pays at least one `gamma` defect, so

```math
rho(T_w)<=sum_(e in w)lambda_e-gamma/|I|.                       \tag{17.7bd}
```

Repeating `w` pumps this fixed loss linearly. Failure to find one particular
support certificate is not itself a falsifier; absence of a tight word cycle
is the observable obstruction.

### Corollary 17.1t (width-two Ising has a strict anticipatory quotient)

Let raw boundary states be `s=(s_1,s_2) in {+-1}^2`, take letters
`{a,b,c}`, and fix `C>0`. Define `T_e(s,t)=g_e(s)+h_e(t)` by

```math
g_a=C(s_2-1)/2,       &h_a=C(t_1-1)/2,
g_b=C(s_1s_2-1)/2,    &h_b=C(t_1-1)/2,
g_c=-C(1+s_1s_2)/2,   &h_c=-C(1+t_1)/2.                       \tag{17.7be}
```

These are width-two Ising fields and pair couplings. Their zero relations
are `D_e times K_(tau(e))`, where

```math
D_a={s:s_2=1},\quad D_b={s:s_1s_2=1},\quad
D_c={s:s_1s_2=-1},

K_+={t:t_1=1},\quad K_-={t:t_1=-1},\qquad
tau(a)=tau(b)=+,\quad tau(c)=-.                                \tag{17.7bf}
```

Every `K_q` meets every `D_e`, hence
`Delta_e(K_q)=K_(tau(e))`. The two supports `K_+,K_-` form an exact
anticipatory carrier and `rho(T_w)=0` for every word. They are minimal among
support carriers: a one-state support stable under both `a` and `c` would
have to lie in the disjoint sets `K_+` and `K_-`.

By contrast, every exact rowwise path-lift partition retains all four raw
states. Their source-membership signatures in `(D_a,D_b,D_c)` are

```math
110,\quad001,\quad101,\quad010.                                \tag{17.7bg}
```

Thus any merged pair has a zero versus negative block-row gap for some
letter. A common diagonal gauge cannot repair the merge: the four source
reward vectors

```math
(0,0,-C),\ (-C,-C,0),\ (0,-C,0),\ (-C,0,-C)
```

have no pair whose difference is constant across all letters. The
response-derived support state is strictly smaller than the usual transfer
state for a genuine switching Ising strip.

The quotient also carries order-sensitive weights without returning to four
states. For `C>4`, replace only

```math
T_a(s,t)\quad\hbox{by}\quad Ttilde_a(s,t)=g_a(s)+h_a(t)-s_1t_1. \tag{17.7bh}
```

In a cyclic word, the zero baseline constraints select the unique boundary
state in `K_(tau(e_(k-1))) intersect D_(e_k)` at phase `k`. If another
periodic path differs at `d` boundary positions, its baseline loses at least
`Cd`; changed endpoints can improve the added `a`-bonds by at most `4d`.
The selected path is therefore uniquely optimal and

```math
rho(Ttilde_w)=2N_(ca)^cyc(w)-N_a(w).                            \tag{17.7bi}
```

The same two-state carrier updates `q'=tau(e)` and uses toll

```math
r(q,a)=-q,\qquad r(q,b)=r(q,c)=0.                              \tag{17.7bj}
```

More explicitly, with `pi(s)=s_1=q`, its exact dominating coarse matrices
are

```math
S_a(q,q')=C(q'-1)/2-qq',
\qquad S_b(q,q')=C(q'-1)/2,
\qquad S_c(q,q')=-C(1+q')/2.                                  \tag{17.7bj1}
```

The supports `K_q=pi^(-1)(q)` lift these matrices with zero defect in
Theorem 17.1s. The `C>4` comparison above shows that their critical cycles
use precisely `q'=tau(e)`, on which (17.7bj1) reduces to the displayed toll.

It is minimal: `aabccb` and `abbcac` have identical Parikh counts but
responses `-2` and `+2`, so no one-state per-letter toll can answer both.
Letters `b,c` restrict a gauge-compatible raw merge to two opposite-spin
pairs, and letter `a` separates each of those by a nonzero row difference
when `C>2`; hence exact rowwise lifting still needs four states. The three
distinct resource counts are

```math
dim(response)=1,\qquad C_(anticipatory)=2,
\qquad C_(forward\ path)=4.                                    \tag{17.7bk}
```

### Theorem 17.1u (bounded-delay residual-core lumpability)

Let decorated letters `c` carry finite real max-plus blocks
`T_c in R^(I_s times I_t)` between finite typed fibres.  Suppose that for
some fixed `D>=1`, every legal `D`-letter product has max-plus row rank one:

```math
T_v(i,j)=A_v(i)+p_v(j),\qquad max_jp_v(j)=0.             \tag{17.7bl}
```

Let `Q` be the finite set of distinct terminal phase/profile pairs.  For
every enabled letter there are unique maps

```math
delta:Q times E->Q,\qquad kappa:Q times E->R
```

such that

```math
p_q tensor T_c=kappa(q,c)1+p_(delta(q,c)).              \tag{17.7bm}
```

This quotient has three exact, distinct consequences.

1. **Residual-cycle response.**  For every context cycle `C`, with word `w`,

   ```math
   rho(T_w)=sum_((q,c) in C)kappa(q,c).                 \tag{17.7bn}
   ```

2. **Accumulated-error dichotomy.**  If every entry of `T_c` is at most
   `lambda_c`, put `d(q,c)=lambda_c-kappa(q,c)>=0`.  Given budgets
   `beta_c>=0`, a potential satisfying

   ```math
   d(q,c)<=beta_c+psi(delta(q,c))-psi(q)                \tag{17.7bo}
   ```

   exists if and only if every residual-context cycle has
   `sum_C(d-beta)<=0`.  A violating cycle is an observable pumpable witness:
   its `k`th repetition exceeds the declared deficit budget by
   `k sum_C(d-beta)`.

3. **Greatest fixed-context support core.**  Define the scalar-threshold
   relation

   ```math
   R_(q,c)(K)=\{j:exists i in K,\ T_c(i,j)>=kappa(q,c)\}. \tag{17.7bp}
   ```

   Starting from allowed endpoint sets `K_q^0=A_q`, iterate

   ```math
   K_q^(r+1)=K_q^r intersection
      intersection_((s,c):delta(s,c)=q)R_(s,c)(K_s^r). \tag{17.7bq}
   ```

   It stabilizes after at most `sum_q|A_q|` deletions.  Its limit is the
   greatest family with one support per residual context satisfying

   ```math
   K_q subseteq A_q,
   \qquad K_(delta(q,c)) subseteq R_(q,c)(K_q).         \tag{17.7br}
   ```

   Hence all required limit supports are nonempty exactly when such a
   one-support-per-context backward-surjective lift exists.

Short prefixes of length below `D` contribute only the finite transient
`sum_(ell<D)|E|^ell`; the recurrent quotient has at most `|E|^D` profiles.

#### Proof

If `v=c_1...c_D` and `v'` is the length-`D` suffix of `vc`, factor
`T_(vc)=T_(c_1)T_(v'c)`.  Row-rank one of the suffix product forces the
normalized row profile of `p_vT_c` to be `p_(v')`.  Equality of two source
profiles makes their one-letter images equal, proving (17.7bm) is well
defined.

Around a context cycle, iteration gives

```math
p_qT_w=K1+p_q,
\qquad K=sum_Ckappa.                                   \tag{17.7bs}
```

A finite max-plus left eigenprofile has eigenvalue equal to the maximum cycle
mean: profile telescoping bounds every cycle above by `K`, while selecting
one maximizing predecessor at each target yields a tight directed cycle.
This proves (17.7bn).  The potential equivalence is the finite difference-
constraints criterion; (17.7bn) makes every violating cycle and all its
repetitions exact response witnesses.

For the core, the right-hand side of (17.7bq) is monotone.  Every family
satisfying (17.7br) lies in each descending iterate, and the stabilized
family itself satisfies (17.7br).  Backward selection through (17.7br)
realizes every finite context path by raw edges of weight at least their
`kappa` tolls. `square`

The threshold relation in (17.7bp) is not the automatically target-surjective
row-argmax relation

```math
\{(i,j):p_q(i)+T_c(i,j)
       =kappa(q,c)+p_(delta(q,c))(j)\}.                \tag{17.7bt}
```

The profile terms prevent either relation from containing the other.  The
core is therefore a substantive local witness certificate, but failure of
this particular core is not a scalar-response falsifier.  For example, with
`p=(0,-1)` and

```math
T_a=\begin{pmatrix}0&-1\\-2&-3\end{pmatrix},\qquad
T_b=\begin{pmatrix}-2&-3\\1&0\end{pmatrix},           \tag{17.7bu}
```

both generators have row rank one and `pT_a=pT_b=p`, so every word has
spectral radius zero.  But their zero-threshold relations are respectively
`{(1,1)}` and `{(2,1),(2,2)}`, whose greatest one-context core is empty.
Allowing multiple support states or using the eigenprofile directly may
still answer the scalar query.

For arbitrary rooted terminal fields, the projective residual metric is
exactly

```math
inf_(a in R)sup_z
|max_i(p_i+z_i)-max_i(p_i'+z_i)-a|
={1\over2}osc(p-p').                                  \tag{17.7bv}
```

Coordinate pins attain both endpoint differences.  Thus the distinct
profiles are minimal exact rooted residual states under complete terminal
probes, but not necessarily under the weaker unrooted spectral query.

The finite row-residual quotient is weighted-subset/determinization algebra.
The fixed-delay row-rank-one hypothesis is a strong checkable synchronizer,
not a necessary condition and not equivalent to the classical twins
property.  The genuinely additional object here is the greatest-core test
for when those residual classes also admit one economical backward extremal-
witness presentation.  The corrected independent audit and source boundary
are recorded in
[`drafts/weighted_residual_core_independent_audit.md`](drafts/weighted_residual_core_independent_audit.md).

### Theorem 17.2 (approximate block lumpability with depth-uniform error)

Use the column max-plus convention

```math
(F_Kx)_i=max_j(K_(ij)+x_j).
```

Let `Pi={I_1,...,I_r}`, choose a gauge `c`, and set

```math
Lambda_c(x)_a=max_(i in I_a)(x_i+c_i).                            \tag{17.8}
```

Suppose a finite alphabet of raw matrices `K_e` and quotient matrices `S_e`,
all with finite entries,
satisfies

```math
|max_(i in I_a)(K_e(i,j)+c_i-c_j)-S_e(a,b)|<=epsilon
quad(j in I_b).                                                   \tag{17.9}
```

Then

```math
||Lambda_c(F_(K_e)x)-F_(S_e)(Lambda_cx)||_H<=epsilon.             \tag{17.10}
```

Assume all actual aggregate and quotient states of interest lie in a
forward-invariant set `Y` contained in a Hilbert ball of radius `R`, every
quotient map is defined and nonexpansive on `Y`, and every contiguous legal
length-`L` quotient composition is `rho`-contractive on `Y`, `rho<1`. For every
`h>0` there is a deterministic finite simulator with at most

```math
(1+2R/h)^(r-1)                                                    \tag{17.11}
```

states and aggregate shadow error, uniformly over all depths and legal switch
words, at most

```math
h+{L(epsilon+h) over1-rho}.                                      \tag{17.12}
```

If, uniformly on all raw states under consideration, the raw response in its
declared norm is within `kappa` of an `L_O`-Lipschitz function of
`Lambda_cx`, its response error is at most

```math
kappa+L_O[h+{L(epsilon+h) over1-rho}].                            \tag{17.13}
```

At `epsilon=0`, before quantization,
`Lambda_cF_(K_e)=F_(S_e)Lambda_c` is an exact semiconjugacy. This
certificate uses the block coefficients and an `(r-1)`-dimensional response
carrier, rather than an active-cell orbit language.

#### Proof

Regrouping maxima gives

```math
Lambda_c(F_(K_e)x)_a
=max_b max_(j in I_b){x_j+c_j+L^e_(aj)},                          \tag{17.14}
```

where `L^e_(aj)` is the block coefficient in (17.9). Replacing it by
`S_e(a,b)` changes every output coordinate by at most `epsilon`, proving
(17.10). Take an internal `h`-net `C subset Y` of `Y` and round after each
quotient update. Explicitly, for a deterministic nearest-net map `Q`, set
`c_0=Q(Lambda_cx_0)` and
`c_t=Q(F_(S_(e_t))c_(t-1))`. Each step injects residual at most
`epsilon+h`. Grouping
transported residuals into length-`L` blocks gives the geometric sum in
(17.12), plus the initial net error `h`. The standard volumetric covering
bound in dimension `r-1` gives (17.11); the response decoder proves (17.13).
`square`

The finite input-word path is exact, while the state trajectory is a metric
shadow, not an exact tie-selector itinerary. Without contraction the same
construction gives `h+T(epsilon+h)`, and that linear loss can be sharp.
The count (17.11) excludes any separate automaton needed to recognize a
restricted legal language. The contraction is on `Y`: globally, a
nonconstant finite max-plus map has Hilbert coefficient one, so strict global
contraction is a projective reset.

The weighted-automaton verifier gives a sharp finite example. A four-state
two-letter automaton has a coarsest strong two-block quotient and 4,802 exact
aggregation checks. Increasing one microscopic self-loop by `delta>0` makes
strong refinement split all four states. The old two-block quotient has
one-step defect exactly `delta`, and a repeatable maximizing self-loop has
error exactly `n delta` at depth `n`. Thus the noncontractive linear estimate
is attained at every horizon, while exact state cardinality jumps
discontinuously.

The cycle-response audit finds a strictly smaller bounded-error repair. If
the perturbed entry is `T_A(2,2)=1+delta`, change only the old quotient toll
`S_A(1,1)` from `1` to `1+delta`. Even microscopic representatives `0,2`
then give an exact section of both quotient matrices, every raw edge is at
most its corrected block edge, and an odd initial representative loses at
most `2+delta` on its first transition. Therefore the corrected two-state
quotient has uniform all-word endpoint error at most `2+delta`, attained by
state `3` under sufficiently many `A` steps, although exact strong refinement
still has four states. The old quotient's `n delta` error was a repairable
cycle toll, not proof that four runtime states were needed.

Only exposed/maximizing cycles may enter this conclusion. The unperturbed raw
cycle `0->1->0` under `A` has mean discrepancy `-3/2` from its quotient loop
but is never maximizing; charging every syntactic microscopic cycle would
falsely reject an exact response quotient.

### Corollary 17.2a (intrinsic contraction from a weighted control graph)

Let a finite legal-control graph have invariant metric fibres `Y_q`. Give
each edge `e:u->v` a restricted Lipschitz coefficient `lambda_e<=1`, and put

```math
beta_n=max_(legal |p|=n) prod_(e in p)lambda_e,
\qquad beta_0=1.                                                 \tag{17.14a}
```

For exact and approximate trajectories with edge defect at most `eta`,

```math
D_n<=beta_nD_0+eta sum_(m=0)^(n-1)beta_m.                        \tag{17.14b}
```

This is sharp from the local coefficients alone. If `rho=beta_L<1`, then

```math
sup_n D_n
<=D_0+eta{sum_(r=0)^(L-1)beta_r\over1-rho}
<=D_0+{Leta\over1-rho}.                                         \tag{17.14c}
```

Equivalently, some `beta_L<1` precisely when every directed cycle has
coefficient product below one, or when the subgraph of `lambda_e=1` edges is
acyclic. If

```math
gamma=max_(cycles C)(prod_(e in C)lambda_e)^(1/|C|)<1,
```

then for every `mu in (gamma,1)` logarithmic difference constraints produce
positive vertex weights `a_q` with

```math
lambda_e{a_v\over a_u}<=mu,
```

and hence `beta_n<=Kmu^n`, `K=max a/min a`. Thus the block-contraction
hypothesis of Theorem 17.2 is certified by a maximum-cycle-mean computation,
not by enumerating products.

The proof of (17.14b) transports each fresh defect through its exact suffix;
submultiplicativity gives (17.14c). The cycle equivalences are the standard
finite weighted-graph potential criterion. A directed `L`-cycle with `L-1`
unit edges and one coefficient `theta<1`, with aligned scalar defects,
attains `Leta/(1-theta)`. At `theta=1` its error grows linearly. Hence both
the contraction denominator and the failure at unit cycle mean are sharp.

For a max-plus map on a convex full-dimensional projective polytope, its
restricted Hilbert coefficient is either zero or one. On a tie-free selector
cell a nonconsensus row-selector has Dobrushin coefficient one; if no such
cell meets the interior, continuity makes the projective map constant. A
zero edge is checked intrinsically by a dominance cone: input coordinate `k`
wins every row on `Y` whenever

```math
x_k-x_j>=max_i(S_(ij)-S_(ik))
\quad(x in Y,\ j).                                               \tag{17.14d}
```

Deleting these reset edges and checking that the remainder is acyclic is
therefore an exact graph certificate for full-dimensional max-plus fibres.
Fractional coefficients remain meaningful on restricted/thin carriers and
for stochastic secants.

### Theorem 17.2b (finite Bellman envelope: potential or pumpable lasso)

On a finite reachable control graph, suppose a nonnegative comparison radius
obeys on edge `e:u->v`

```math
r^+<=eta_e+lambda_e r,
\qquad eta_e,lambda_e>=0,                                      \tag{17.14e}
```

with initial radii represented by source edges. For a path `p`, write its
ordered affine comparison map as

```math
F_p(x)=Lambda_p x+A_p.
```

The worst envelope

```math
B_v=sup_(p:source->v)F_p(0)                                    \tag{17.14f}
```

is the least extended fixed point of

```math
B_v=max_(e:u->v)(eta_e+lambda_eB_u).                            \tag{17.14g}
```

The following are equivalent:

1. every `B_v` is finite;
2. there is a finite Bellman potential `P>=0` satisfying
   `P_v>=eta_e+lambda_eP_u` on every reachable edge;
3. there is no reachable pumpable lasso `pC`.

If `x=F_p(0)` is the entry radius of a closed walk `C`, the exact pump
criterion is

```math
Lambda_C>=1
\quad\hbox{and}\quad
A_C+(Lambda_C-1)x>0.                                            \tag{17.14h}
```

For `Lambda_C>1` repetition grows exponentially; for `Lambda_C=1,A_C>0`
it grows linearly. If `Lambda_C<1`, a positive one-cycle defect converges to
a finite fixed radius and is not a pump. Path iteration proves (17.14g);
a finite supersolution bounds every path, while a divergent path in a finite
graph yields a repeated closed walk satisfying (17.14h). Conversely,

```math
F_C^k(x)=
\begin{cases}
Lambda_C^kx+A_C( Lambda_C^k-1)/(Lambda_C-1),&Lambda_C ne1,\\
x+kA_C,&Lambda_C=1,
\end{cases}                                                     \tag{17.14i}
```

proves pumping. The least finite potential is obtained by the linear
inequalities in item 2. Positive log-gain cycles are found by maximum cycle
mean; zero-log cycles carrying positive defect are the critical linear case.

This is an iff for the nonnegative comparison envelope, not for a fixed
signed or vector error. A forced two-cycle `z|->z+1`, `z|->z-1` is bounded,
whereas replacing the second map by `z|->z+1` drifts linearly; both have the
same local absolute data `lambda=eta=1`. Thus an envelope pump is a realized
drift witness only under an alignment/attainability certificate. This is
exactly why joint cancellation cannot be inferred from separately paid local
channels.

### Theorem 17.3 (response packing, memory gain, and predictive compactness)

Let `E` be finite and let `R:X->Z` be `L_R`-Lipschitz, `L_R>0`, and define

```math
d_T(x,y)=max_(|w|<=T)d_Z(R(F_wx),R(F_wy)).                        \tag{17.15}
```

An `epsilon`-predictor consists of a finite set `S`, an encoder `q:X->S`,
deterministic input updates `delta_e`, and a decoder `g:S->Z` such that

```math
d_Z(R(F_wx),g(delta_wq(x)))<=epsilon                             \tag{17.16}
```

for all `|w|<=T`. Let `C_T(epsilon)` be its minimum state count. Define

```math
G_T=sup_(t<=T,e_1...e_t)
[Lip(F_(e_t)...F_(e_1))
 +sum_(s=1)^t Lip(F_(e_t)...F_(e_(s+1)))].                       \tag{17.17}
```

Then

```math
Pack_(>2epsilon)(X,d_T)
<=C_T(epsilon)
<=Cov^int_(epsilon/(L_RG_T))(X).                                \tag{17.18}
```

For nonexpansive maps, `G_T<=T+1`. If every contiguous length-`L` block
occurring in a legal trajectory is also `rho`-contractive, then
`G_infinity<=L/(1-rho)`.

If `Z` is compact and `C_infinity(epsilon)` denotes the minimum state count
of one predictor valid for every finite word, with value `+infinity` if none
exists, then

```math
\boxed{C_infinity(epsilon)=sup_(T>=0)C_T(epsilon).}               \tag{17.19}
```

#### Proof

Two points with one encoded state receive the same prediction after every
word, so their contextual distance is at most `2epsilon`, proving the lower
bound. For the upper bound, encode by a nearest center of an internal
`eta`-net and quantize after every exact update. The initial and fresh errors
are each at most `eta`; (17.17) transports their sum. Setting
`eta=epsilon/(L_RG_T)` proves (17.18). The two gain estimates follow by
summing unit suffix bounds or a block-geometric series.

For (17.19), choose a strictly increasing sequence of horizons with predictors
on a common `S`-state set, padding unused states with self-loops and a fixed
decoder value. Pass to a subsequence on which their finitely many transition
maps are identical, then use compactness of `Z^S` to make every decoder value
converge. On this final subsequence, for each `x` choose a state occurring
infinitely often among its horizon-dependent initial encodings. Along that
occurrence subsequence,
(17.16) passes to the limit for every fixed finite word. This gives one
infinite-depth `S`-state predictor. The reverse inequality is restriction.
`square`

This predictor need not be an exact quotient: the encoder is not asserted to
satisfy `qF_e=delta_eq`. Theorem 17.1 is the stronger exact-congruence law.

### Corollary 17.3a (behavioral recoupling of every finite predictor)

Let `(S,q,delta,g)` be an infinite-depth `epsilon`-predictor. Define

```math
d_S(s,t)=sup_(w in E^*)d_Z(g(delta_ws),g(delta_wt)).              \tag{17.19a}
```

Then `d_S` is a pseudometric, `g` is one-Lipschitz, every transition is
nonexpansive, and

```math
d_S(q(F_ex),delta_eq(x))<=2epsilon.                              \tag{17.19b}
```

Moreover, `d_S` is determined by the finite synchronized pair graph:

```math
d_S(s,t)=max{d_Z(g(u),g(v)):(u,v)
               reachable from (s,t)}.                            \tag{17.19c}
```

After deleting states unreachable from `q(X)`, quotient `S` by zero
behavioral distance. This is the minimal Moore realization of this
predictor's exact decoded response trees (not necessarily the smallest
`epsilon`-predictor for the physical system). If `gamma` is the least positive
distance between distinct quotient states, with `gamma=+infinity` for a
one-state quotient, then `2epsilon<gamma` makes the induced encoder an exact
forward semiconjugacy.

Indeed, the empty word proves the decoder bound, and prepending one letter
proves transition nonexpansiveness. For every future word `w`, the predictions
from `q(F_ex)` and `delta_eq(x)` both approximate `R(F_wF_ex)` within
`epsilon`, proving (17.19b). Formula (17.19c) is the finite word orbit of the
state pair. Zero distance is therefore a transition congruence, and the
strict margin forces the defect in (17.19b) to be zero. The pair search is
effective whenever the finitely many decoder distances are effectively
represented. `square`

The constant is sharp. A one-point zero-response system has an error-one
predictor which outputs `-1` initially and then transitions to a fixed state
outputting `+1`; its behavioral defect is exactly two. Thus equality at the
margin need not yield exact semiconjugacy.

### Theorem 17.3b (behavioral entropy times dynamic suffix gain)

Let `(S,q,delta,g)` be any infinite-depth `epsilon`-predictor. Remove
unreachable states and quotient by the zero classes of (17.19a). Choose an
internal `eta`-net `C` in its finite behavioral metric and a retraction
`Q:S->C`, and define

```math
qhat=Qq,
\qquad deltahat_e(c)=Qdelta_e(c),
\qquad ghat=g|_C.                                                \tag{17.19d}
```

Let `G_T` be (17.17) for the finite transition maps `delta_e` in `d_S`.
Then the rounded predictor has error at most

```math
epsilon+eta G_T                                                 \tag{17.19e}
```

through horizon `T`. If every legal length-`L` transition block is
`rho`-contracting, it is valid at every depth with error

```math
epsilon+{Leta\over1-rho}                                       \tag{17.19f}
```

and uses at most `Cov^int_eta(S,d_S)` states. The proof compares the rounded
trajectory with the original predictor trajectory; the initial and fresh
retraction residuals are transported by suffix maps, while `g` is
one-Lipschitz. No physical semiconjugacy is needed. Corollary 17.3a's strict
gap condition supplies the separate conclusion that re-encoding the evolved
physical state agrees exactly with its predicted behavioral class.

The contraction factor is an essential multiplier, not decoration. Let

```math
S_n={0,...,n},
\qquad delta(i)=min(i+1,n),
\qquad g(i)=alpha i,                                             \tag{17.19g}
```

and take this exact predictor as the physical system. Its behavioral metric
is `alpha|i-j|`. For `eta=kalpha`, `1<=2k<=n`,

```math
Cov^int_eta(S_n)=ceil{n+1\over2k+1},
\qquad C_infinity(eta)=n-2k+1.                                  \tag{17.19h}
```

For the lower bound, any repeated predictor state along the orbit from zero
is eventually periodic. Every decoder value on its cycle must be within
`eta` of `nalpha`, while its first occurrence at time `i` must be within
`eta` of `ialpha`; thus `(n-i)alpha<=2eta`. The first `n-2k+1` states are
distinct. Keeping the earlier transient exactly and merging the last
`2k+1` states into one midpoint-decoded sink proves equality. Taking
`k` proportional to `n` gives a constant static cover but linear reusable
memory. This finite grid is invariant under the compact two-piece PWA map
`F(x)=min(x+alpha,nalpha)`; its forgetting time grows like `n`.

### Theorem 17.4 (contraction-weighted context-tree entropy)

For `h:X->Z`, put

```math
S_h(r,epsilon)=sup_(A subseteq X,diam A<=r)
Cov^ext_epsilon(h(A)).                                           \tag{17.20}
```

Suppose `|E|=q`, `diam X<=D`, and every `F_e` is `rho`-Lipschitz,
`0<rho<1`. For the response tree

```math
Phi_T(x)=(h(F_wx))_(|w|<=T),                                    \tag{17.21}
```

one has

```math
log Cov^ext_epsilon(Phi_T(X))
<=sum_(k=0)^T q^k log S_h(D rho^k,epsilon).                       \tag{17.22}
```

If `epsilon>0` and `h:X->R^p` is `L`-Lipschitz in sup norm, with finite
`L,D` and `q>=1`, the right side yields

```math
Cov^ext_epsilon(Phi_T(X))
<=prod_(k=0)^T
 max{1,ceil(LD rho^k/(2epsilon))}^(p q^k).                        \tag{17.23}
```

Every factor and exponent is sharp over finite affine systems. In
particular, for `A=LD/(2epsilon)>1`, `q>=2`, and horizons beyond the mixing
scale,

```math
log Cov^ext_epsilon(Phi_T(X))
=Theta_(p,q,rho)(A^(log q/log(1/rho)))                            \tag{17.24}
```

in a matching example. For `q=1`, the sharp order is
`Theta_(p,rho)((log A)^2)`.

#### Proof

Each word image at level `k` has diameter at most `D rho^k`; independently
cover its response image and multiply over the `q^k` words, proving (17.22).
A diameter-`r` response lies in a `p`-box of side at most `Lr`, proving
(17.23).

For sharpness, for each `T` put
`X=prod_(|u|<=T)[0,D]^p`, let `h(x)=Lx_empty`, and set

```math
(F_ax)_u=cases(rho x_(au),&|u|<T;
               0,&|u|=T).                                       \tag{17.25}
```

The word `a_1...a_k` reads `Lrho^k x_(a_1...a_k)`, so the response image is

```math
prod_(k=0)^T[0,LD rho^k]^(p q^k).                                \tag{17.26}
```

Its exact sup-covering number is (17.23): equal subintervals cover, while
equally spaced points including endpoints give the matching strict packing.
Only `Theta(log A/log(1/rho))` levels are active. The deepest active level
and a geometric sum give (17.24); with `q=1`, summing the linearly decreasing
logarithms gives quadratic order. `square`

This sharpness uses a horizon-dependent family of finite affine systems whose
dimension grows with the active context tree; it is a distribution-free
sharpness statement, not a fixed-dimensional asymptotic. Every predictor
gives an external response-tree center, so
`Cov^ext_epsilon(Phi_T(X))<=C_T(epsilon)`. The reverse is false in general because
external table centers need not admit shift-consistent transitions.

### Theorem 17.5 (bounded static response, exponential dynamic memory)

On `X=[0,1]`, observe `h(x)=x`, use

```math
E_0(x)=x/3,
\qquad E_1(x)=(x+2)/3,
```

and the continuous piecewise-affine decoder

```math
R(x)=cases(3x,&0<=x<=1/3;
           2-3x,&1/3<=x<=2/3;
           3x-2,&2/3<=x<=1).                                    \tag{17.27}
```

Although the numerical one-step response vector
`(h,h circ E_0,h circ E_1,h circ R)` has an `epsilon`-cover of
size at most `ceil(3/(2epsilon))+1`, after `t` encoders every predictor which
answers the decoder continuations to error `epsilon<1/6` needs at least
`2^t` states.

Also, for `m>=3`, pairwise and indeed every proper-subfamily compatibility of
max-plus tie faces is insufficient. For cyclic coordinates modulo constants, choose
`sum_i c_i=1` and set

```math
phi_i(x)=max{x_(i+1),x_i+c_i}.
```

Every proper subfamily of the `m` tie faces intersects, but

```math
min_x max_i|x_(i+1)-x_i-c_i|=1/m.                                \tag{17.28}
```

#### Proof

The decoder satisfies `RE_b=id`. Two different encoder words first differ in
one reversed ternary digit. After stripping their common leading digits with
`R`, one state lies in `[0,1/3]` and the other in `[2/3,1]`; some future
response differs by at least `1/3`. A common predictor state would put every
future pair within `2epsilon`, proving the packing. The one-step vector is
`3`-Lipschitz, giving its scalar-grid cover.

For the tie claim, write `r_i=x_(i+1)-x_i-c_i`. Deleting one cycle equation
leaves a forest and hence a soluble proper subsystem. But
`sum_i r_i=-1`, so `max_i|r_i|>=1/m`. Setting every `r_i=-1/m` is integrable
around the cycle and attains equality. `square`

### Theorem 17.6 (strict-strip response isometry and restricted lumpability)

For a width-`w` Ising strip prefix, let `f(x)` be the maximum prefix energy
conditioned on boundary spins `x in {-1,1}^w`. Under arbitrary real legal
strip continuations,

```math
D_ctx(f,g):=sup_C|Opt_f(C)-Opt_g(C)|=||f-g||_infinity.            \tag{17.29}
```

Hence the exact contextual state is the full `2^w`-entry boundary table, and
one-time sup approximation has the same error under every future depth.
Modulo constants the metric is `osc(f-g)/2`.

For a fixed finite rational column alphabet whose coefficients and seed lie
in `eta Z`, normalized profiles lie in a finite set. If

```math
B=max_c max_(y,z)[A_c(y)-A_c(z)
 +2sum_(i:y_i!=z_i)|J_i^c|],                                    \tag{17.30}
```

then, with `q=2^w` and `K=floor(B/eta)`, there are at most

```math
(K+1)^q-K^q                                                     \tag{17.31}
```

normalized profiles. Weighted partition refinement by `(successor,
baseline toll)` gives the coarsest exact restricted quotient.

#### Proof

Max-plus nonexpansiveness gives the upper half of (17.29). To expose
coordinate `a`, append one column with horizontal couplings `M>0`, fields
`La_i` with `L>M`, and no vertical interactions. The new boundary is uniquely
`a`; choosing `M` larger than half every competing profile gap also makes
the old boundary uniquely `a`. The continued optimum difference is
`f(a)-g(a)`, proving the reverse inequality.

For (17.31), the difference of two maxima lies between the minimum and
maximum pointwise kernel differences. Thus every normalized successor has
span at most (17.30), lies on the `eta`-lattice, and has maximum zero.
Counting all `(K+1)^q` grid vectors and excluding the `K^q` vectors with no
zero gives (17.31). Finite weighted Myhill--Nerode refinement proves
the last statement. `square`

The exact verifier exhibits a strict width-three response chart with
determinant `-1024`, a closed eight-dimensional sup-cube of radius `1/2`, and
seven projective degrees. A rational width-two two-letter alphabet has seven
reachable normalized profiles and a strict two-state weighted quotient.

Repeatedly re-quantizing need not inherit the one-time isometry. In the
embedded one-spin model with coupling `-K`, field `0<s<K`, and grid mesh
`Delta>4s`, every approximate projective state rounds to flat and pays toll
`K+s`; the exact optimum after `n` steps is

```math
nK+(n mod 2)s.                                                    \tag{17.32}
```

The absolute error is `ns` for even `n` and `(n-1)s` for odd `n`. Both exact
and approximate projective controls return after two steps, but their reward
cycle differs by `2s`. Thus projective lumpability and scalar-cocycle
compatibility are separate necessities.

There is also a positive switching approximation not requiring a common
optimizer cell. For one column with horizontal couplings `J_i`,

```math
||[T_cf]-[T_cg]||_H<=2sum_i|J_i|                                 \tag{17.33}
```

for arbitrary incoming profiles. If every `L` consecutive columns contain
one with `2sum_i|J_i|<=delta`, and each normalized update is rounded with
Hilbert error at most `eta`, then after the first such column the projective
error is at most `delta+Leta` at every depth. Indeed, the difference of the
two output maxima has oscillation at most `4sum_i|J_i|`; a weak column resets
all prior error, while subsequent nonexpansive updates add at most one fresh
`eta` each. This is a restricted-image reset, not a strict global Lipschitz
coefficient.

### Theorem 17.7 (terminal residuals, scalar compatibility, and small-shell universality)

Let finite max-plus matrices act on rows.  For a row `u` and terminal field
`z`, define the normalized terminal response

```math
R_z(u)=max_j(u_j+z_j)-max_j u_j.                    \tag{17.34}
```

Fix `D`.  Suppose every legal length-`D` product `T_v` has a factorization

```math
|T_v(i,j)-a_v(i)-p_v(j)|<=epsilon,                  \tag{17.35}
```

and the profiles `[p_v]` have an `eta`-net `{[r_c]}` in projective
sup distance.  If a word `w` ends in `v`, then for every raw initial state
`i` and terminal field `z`,

```math
|R_z(T_w(i,.))-R_z(r_(chi(v)))|
<=2(epsilon+eta).                                    \tag{17.36}
```

The bound is independent of word depth and of `z`.  The suffix is always a
finite reusable state; its smaller code label is reusable when the label is
a right congruence under suffix shift.

This terminal theorem does not control accumulated scalar response.  If
exact rank-one generators are `T_e=a_e tensor p_e`, their directed
compatibility table

```math
varphi(e,f)=max_j(p_e(j)+a_f(j))                     \tag{17.37}
```

obeys, on every legal cyclic word,

```math
rho(T_(e_1)...T_(e_t))=
sum_(s=1)^t varphi(e_s,e_(s+1)),\qquad e_(t+1)=e_1. \tag{17.38}
```

Approximating length-`D` blocks entrywise by rank-one factors with errors
`epsilon_v`, and replacing `varphi` by a quotient table, gives asymptotic
spectral error per original letter at most

```math
{max_v epsilon_v+Delta\over D},                      \tag{17.39}
```

where `Delta` is the largest absolute mean defect on a directed simple
cycle of the finite legal-block graph.  At zero factor error, bounded
absolute error is possible exactly when every cycle defect vanishes, or
equivalently when the defect is a potential on each strongly connected
component.

Both distinctions are sharp.  For every `delta>0`, the two exact rank-one
matrices

```math
A=\begin{pmatrix}0&0\\0&0\end{pmatrix},\qquad
B=\begin{pmatrix}delta&2delta\\0&delta\end{pmatrix} \tag{17.40}
```

have zero projective contraction and profiles with optimal one-centre
radius `delta/4`, yet every one-dynamic-state predictor with
letter-dependent tolls has optimal scalar response distortion exactly
`delta/4` per letter.  More generally, scale any all-finite weighted-
automaton alphabet with entries in `[-1,0]` by `alpha`.  Every nonempty
product then lies in a one-profile row shell of radius `alpha/2`, while its
entire spectral response algebra is the original algebra scaled by `alpha`.

#### Proof

Write `w=uv`.  Substituting (17.35) after the prefix maximum leaves one
scalar plus `p_v(j)` and an error in `[-epsilon,epsilon]`; projective
replacement costs `eta`, and the two maxima in (17.34) give (17.36).
Rank-one multiplication gives

```math
(a_e tensor p_e)(a_f tensor p_f)
=varphi(e,f)+a_e tensor p_f,                          \tag{17.41}
```

which proves (17.38).  Entrywise perturbations and spectral radius are
one-Lipschitz; simple-cycle deletion proves (17.39) and the potential
criterion.

For (17.40), the compatibility table is
`delta [[0,1],[1,1]]`.  The cycles `A,B,AB` force error at least
`delta/4`; tolls `(delta/4,5delta/4)` attain it because a cyclic word with
`k` runs of `A` has defect `delta(k-|w|/4)` and
`0<=k<=|w|/2`.  Finally, changing only the last edge of an optimal path in a
scaled `[-1,0]` alphabet changes any fixed product row by at most `alpha`.
Rowwise centering gives radius `alpha/2`, while positive scaling commutes
with max-plus multiplication and spectral radius. `square`

Thus terminal response error is paid once, while unresolved compatibility
is paid once per repeatable cycle.  Small projective radius, even combined
with perfect one-step forgetting, cannot bound cumulative-response memory.
The full proof, support-core guardrail, and independent audit are in
[`drafts/approximate_residual_shell_law.md`](drafts/approximate_residual_shell_law.md)
and its adjacent audit files.

## 18. Bridge interfaces, synchronization, and proof-memory separation

The results in this section were proved during the bridge-hierarchy campaign.
They separate algebraic bridge rank, exposed response information, graph
sparsity, deterministic synchronization, and certificate memory.

### Theorem 18.1 (bridge-query isometry)

Let `X,Y` be finite, let `B:X times Y -> R`, and for `h:X->R` define

```math
(P_Bh)(y)=max_(x in X){h(x)+B(x,y)}.                           \tag{18.1}
```

With arbitrary future weights `g:Y->R`, put

```math
Opt_B(h,g)=max_(x,y){h(x)+B(x,y)+g(y)}.
```

Then

```math
sup_g|Opt_B(h,g)-Opt_B(h',g)|
=||P_Bh-P_Bh'||_infinity.                                     \tag{18.2}
```

Modulo additive constants, the quotient distance is

```math
inf_c||P_Bh-P_Bh'-c1||_infinity
={1\over2}osc(P_Bh-P_Bh').                                    \tag{18.3}
```

#### Proof

The maximum is one-Lipschitz in the sup norm, proving one inequality. Choose
a coordinate attaining the largest signed response difference and make every
other future weight sufficiently negative. Both optimizations are then pinned
to that coordinate, proving the reverse inequality. Best uniform
approximation of a finite vector by a constant is half its oscillation.
`square`

Thus the exact state of a declared bridge is its *realizable response image*,
not rank, edge count, or symmetry in isolation. A rank factorization restricts
this table to a finite-dimensional field set; a width-`w` separator gives a
`2^w` table; a symmetry quotient gives one coordinate per relevant orbit.

### Theorem 18.2 (exact finite-rank roof algebra)

A finite `r`-featured landscape is `(X,H,phi)` with
`phi:X->R^r`. Its linear-field response and upper concave roof are

```math
V_L(t)=max_x{H(x)+<t,phi(x)>},

Hbar_L(u)=max\left\{sum_xp_xH(x):p in Delta_X,
                         sum_xp_xphi(x)=u\right\}.             \tag{18.4}
```

For two such landscapes define

```math
(H star G)(x,y)=H(x)+G(y)+<phi(x),psi(y)>,
\qquad (phi star psi)(x,y)=phi(x)+psi(y).                      \tag{18.5}
```

Then

```math
V_(L star K)(t)=max_(u,v)
{Hbar_L(u)+Gbar_K(v)+<u,v>+<t,u+v>}.                           \tag{18.6}
```

On lifted points use

```math
(u,h)o(v,k)=(u+v,h+k+<u,v>)                                   \tag{18.7}
```

and then take the upper concave hull. This roof operation is associative; for
`m` landscapes its energy is

```math
sum_iH_i(x_i)+sum_(i<j)<phi_i(x_i),phi_j(x_j)>.                \tag{18.8}
```

The final concavification is essential. If each child has features `+-1` and
height zero, the fixed-total raw profile is `z^2/4`, while the parent upper
roof is identically one on `[-2,2]`.

#### Proof

The expression in (18.6) is separately affine in the two convexifying
probability vectors, so its global maximum is attained on original states.
Bi-affinity gives

```math
conv(A o B)=conv(conv(A) o conv(B)).                            \tag{18.9}
```

Both parenthesizations of three lifted points have energy
`h+k+l+<u,v>+<u,w>+<v,w>`, proving associativity. `square`

If `R=UV^T` has rank at most `r`, take `phi(x)=U^Tx` and
`psi(y)=V^Ty`. Thus arbitrary internal landscapes coupled through this fixed
port have an exact `r`-dimensional semantic interface. The minimal interface
is query-relative: on a fixed future signature set `Theta` it is only
`V_L|_Theta`; under all singleton fields it is the complete roof by concave
Fenchel biconjugacy.

### Theorem 18.3 (fixed-rank approximation and the growing-rank boundary)

Assume `||phi(x)||<=P` and allowed fields satisfy `||t||_*<=Q`. Quantize the
feature ball by a map with error at most `eta`, retain only nonempty buckets,
and store

```math
w(c)=max_(x:Q_eta(phi(x))=c)H(x).
```

Then

```math
Vtilde_L(t)=max_c{w(c)+<c,t>},
\qquad |V_L(t)-Vtilde_L(t)|<=Q eta,                            \tag{18.10}
```

using at most

```math
(1+2P/eta)^r                                                   \tag{18.11}
```

buckets. Quantizing heights at mesh `zeta` adds at most `zeta` error. If two
bridge sides of radii `P,Q` are quantized at `eta_1,eta_2`, the cross-term
error is at most

```math
Q eta_1+(P+eta_1)eta_2.                                       \tag{18.12}
```

This exponential dependence on rank is unavoidable at fixed error. For every
fixed `0<rho<1/2`, there is a class with unit feature radius containing

```math
2^(2^(Omega_rho(r)))
```

response functions separated by a fixed positive constant. Hence a uniform
fixed-error summary needs `2^(Omega(r))` bits.

#### Proof

Replacing one affine term by its bucket center changes it by at most `Qeta`,
which proves (18.10); a volumetric net proves (18.11), and expansion of the
bilinear term proves (18.12).

For the lower bound choose a code `C subset {-1,1}^r` of relative distance
`rho` and exponential size, put `p_x=x/sqrt(r)`, and, for
`b in {0,1}^C`, define

```math
H_b(x)=a b_x\ (x in C),\qquad H_b(x)=-1\ (x notin C),
\qquad 0<a<2rho.                                               \tag{18.13}
```

At the actual Boolean field `q_c=c/sqrt(r)`, the state `c` is uniquely
optimal and

```math
V_b(q_c)=1+a b_c.                                              \tag{18.14}
```

The `2^|C|` responses are pairwise `a`-separated, so error below `a/2` needs
at least `|C|=2^(Omega(r))` bits. `square`

At bounded fixed scale, (18.11) gives polynomial state size when
`r=O(log n)`, subpolynomial size when `r=o(log n)`, and reaches arbitrary
`n`-spin table scale only around `r=Theta(n)`. Exact complexity behaves very
differently: already at rank one, setting

```math
p_x={sum_i2^(i-1)x_i\over2^n-1},\qquad H(x)=-p_x^2             \tag{18.15}
```

makes all `2^n` configurations uniquely exposed by the fields `2p_x`. Their
exposure margins shrink exponentially, reconciling this with fixed-error
compression.

### Theorem 18.4 (bounded bridge degree does not imply compression)

For every fixed `0<delta<1/2`, there are codes

```math
C_n subset {-1,1}^n,
\qquad |C_n|>=2^((1-h_2(delta)-o(1))n),                        \tag{18.16}
```

and landscapes indexed by `sigma in {0,1}^C_n` such that the degree-one
matching bridge `B(x,y)=<x,y>` has

```math
||P_Bh_sigma-P_Bh_tau||_infinity=delta n
\quad(sigma!=tau).                                             \tag{18.17}
```

Consequently every absolute-score summary with error below `delta n/2`
needs `2^|C_n|` states, or `|C_n|` bits. The same bit exponent holds
projectively.

#### Proof

Take a greedy code of distance at least `delta n` and put
`h_sigma(c)=delta n sigma(c)`. At query `c`, another codeword scores at most
`n-delta n`, while `c` scores `n+delta n sigma(c)`. Thus

```math
(P_Bh_sigma)(c)=n+delta n sigma(c).                            \tag{18.18}
```

For the projective claim restrict to constant-weight labels; two distinct
labels have both signs in their response difference, so (18.3) gives distance
`delta n`. `square`

The bridge graph has treewidth one but an extensive live interface. Sparse
dynamic programming compresses through small live separators, not through
edge degree alone. The internal landscapes here are arbitrary, so this is a
theory falsifier rather than a lower bound for quadratic signings.

### Theorem 18.5 (anticipatory support is certificate complexity)

Let `E` have `q>=2` letters, `I=E^m`, `m>=1`, and `C>0`. Define

```math
F_e(s_1...s_m)=s_2...s_me,

T_e(s,t)=cases(0,&t=F_e(s); -C,&t!=F_e(s)).                   \tag{18.19}
```

Against the scalar coarse system `S_e=0`, every word has spectral response
zero. Nevertheless every exact Theorem-17.1s anticipatory-support carrier has
at least `q^m` states, and this is attained.

More generally, if a carrier of `N` states has uniform certificate toll
`beta`, then

```math
ell beta<C,\quad ell<=m\quad\Longrightarrow\quad N>=q^ell.   \tag{18.20}
```

Conversely, for `0<=L<m`, suffix-cylinder supports give

```math
N_L={q^(L+1)-1\over q-1},\qquad beta_L={C\over L+1}.           \tag{18.21}
```

Thus for `1<=N<q^m`, the optimal support-certificate toll obeys

```math
beta_N=Theta\left({C\over1+log_qN}\right).                    \tag{18.22}
```

#### Proof

Every word map has a fixed point obtained from the periodic extension of the
word, proving zero spectral response. Exact upper domination forces any
allowed gauge to be constant, since `h(s)<=h(F_e(s))` on the strongly
connected de Bruijn graph.

Start from a carrier support of maximal potential. If `beta=0`, its successor
must have zero shortfall and the same maximal potential. After any word
`u in E^m`, the deterministic image of a nonempty support is the singleton
`{u}`. The `q^m` words therefore reach distinct support states.

For (18.20), the same maximal-potential start gives total shortfall at most
`ell beta<C`; every step is consequently a zero edge, and the `q^ell` word
images lie in disjoint suffix cylinders. For (18.21), retain all suffix
cylinders of depths at most `L`, follow a zero edge below depth `L`, and reset
to the full support at cost `C`. The potential

```math
psi(K_u)=-{|u|C\over L+1}
```

makes every ordinary edge and reset meet the toll inequality exactly.
`square`

The semantic response state is one point, and a one-state rowwise path lift
also exists because every raw state has the successor `F_e(s)`. Hence the
exponential quantity is proof/certificate memory, not intrinsic response
information. Together with the free-tail de Bruijn and width-two Ising
examples, this proves that source-total rowwise lifts and target-surjective
anticipatory supports are incomparable.

For a fixed support architecture, its optimal toll is exactly a finite
mean-payoff game: the controller chooses successor supports, edge costs are
their least backward-surjective shortfalls, and a potential exists iff every
selected cycle has nonpositive adjusted cost. Positional strategies suffice.

### Theorem 18.6 (signed-balance synchronization and its limit)

For blocks `x^a in {-1,1}^n`, let `k_a` be the number of plus spins and let
the internal energy be an arbitrary `h_a(k_a)`. Couple block pairs by

```math
R_ab=alpha_ab I+beta_ab J.                                    \tag{18.23}
```

If every `alpha_ab>=0`, the exact optimum is

```math
max_(k_1,...,k_m)
\left\{sum_ah_a(k_a)+sum_(a<b)[
alpha_ab(n-2|k_a-k_b|)
+beta_ab(2k_a-n)(2k_b-n)]\right\}.                            \tag{18.24}
```

More generally, (18.24) holds after vertex sign gauges iff the signed graph
of nonzero `alpha_ab` has positive sign product around every cycle. This is
equivalent to `sgn(alpha_ab)=epsilon_aepsilon_b` for vertex signs
`epsilon_a`.

An isolated unbalanced unit-sign cycle of length `ell`, with even `n` and all
blocks pinned to zero magnetization, has true identity-channel optimum

```math
(ell-2)n,                                                      \tag{18.25}
```

although the sum of separately optimized edge responses is `ell n`.

#### Proof

For plus sets `P_a`,

```math
(x^a)^Tx^b=n-2|P_a\triangle P_b|
<=n-2|k_a-k_b|.                                                \tag{18.26}
```

Nested initial segments `P_a={1,...,k_a}` attain equality for every pair
simultaneously. The `J` term depends only on magnetizations. Signed cycle
balance is exactly the path-independence condition for the vertex gauge.

On an unbalanced cycle, realized edge products multiply to `+1` at every
coordinate while the desired signs multiply to `-1`, so one edge is
unsatisfied and the coordinate reward is at most `ell-2`. Pair one maximizing
assignment with its global negative on equal halves of the coordinates to
make every block balanced and attain the bound. `square`

The condition therefore characterizes when one common representative section
attains all edgewise conditioned optima. The failure is an extensive
holonomy obstruction to the separable pair-potential algebra, not to every
possible joint response table.

There is a genuine thermodynamic consequence. Fix `m` and a signed-balanced
graph, gauge to coefficients `a_ab=|alpha_ab|`, scale the dense coefficient
as `b_ab/n`, and assume uniformly

```math
n^(-1)h_(a,n)(x)=f_a(n^(-1)sum_i x_i)+o(1)                    \tag{18.27}
```

for continuous `f_a`. Then the normalized optimum converges to

```math
max_(u in [-1,1]^m)\left\{
sum_af_a(u_a)+sum_(a<b)[a_ab(1-|u_a-u_b|)+b_abu_au_b]
\right\}.                                                     \tag{18.28}
```

Indeed (18.24) reduces the finite problem to the parity grid in the compact
cube, and the objectives converge uniformly to (18.28).

The matrices `alpha I+beta J` are dense when `beta!=0` and full rank when
`alpha!=0` and `alpha+n beta!=0`. Thus deterministic synchronization can
beat algebraic rank in a nontrivial full-rank dense class. The state grows
with the number of blocks; this is an exact factor algebra, not a bounded
universal state for arbitrary dense couplings.

### Theorem 18.7 (scale-sensitive spectral bridge reduction)

For arbitrary Boolean state sets `X subseteq {+-1}^p`,
`Y subseteq {+-1}^q`, arbitrary internal landscapes `H,K`, and real bridges
`R,S`,

```math
\left|\max_{x,y}\{H(x)+K(y)+x^TRy\}
      -\max_{x,y}\{H(x)+K(y)+x^TSy\}\right|
\le\sqrt{pq}\,||R-S||_{2\to2}.                          \tag{18.29}
```

This remains true after adding any shared state-specific future. For a best
rank-`r` truncation `R_r`, the right side is
`sqrt(pq) sigma_(r+1)(R)`. The constant is sharp on pinned Boolean singular
vectors.

At balanced size `n` and scale `n^(3/2)`, define

```math
r_\epsilon(R)=#{j:\sigma_j(R)>\epsilon\sqrt n}.         \tag{18.30}
```

Then the interaction is uniformly within `epsilon n^(3/2)` of one factoring
through `r_epsilon(R)` singular features. If `||R||<=Csqrt(n)`, continuations
are declared to depend only on these retained features, and their total dual
field radius is bounded by a fixed number of opposite ports, a feature net
with

```math
\left(1+O_C(1/\epsilon)\right)^r                         \tag{18.31}
```

cells and one additive baseline gives a further
`epsilon n^(3/2)` response approximation. Spectral-tail and quantization
budgets add.

For a graph of bridges `R_e`, rank-`r_e` replacements `S_e` satisfying
`||R_e-S_e||<=delta_e sqrt(n)` change every global objective pointwise by at
most

```math
n^{3/2}\sum_e\delta_e.                                  \tag{18.32}
```

The truncated system has an exact local roof factorization of dimension

```math
d_v\le\sum_{e\ni v}r_e                                  \tag{18.33}
```

at vertex `v`. This local presentation composes associatively while every
exposed port is retained. It is not a bounded global quotient: eliminating a
region creates a joint boundary factor controlled by separator size or
treewidth, which may be linear even at bounded degree.

#### Proof

For every Boolean pair,

```math
|x^T(R-S)y|\le||x||_2||R-S||_{2\to2}||y||_2,            \tag{18.34}
```

and maxima are one-Lipschitz under a pointwise perturbation. Eckart--Young
proves (18.30). Factor
`R_r=U Sigma V^T` and use features
`Sigma^(1/2)U^Tx`, `Sigma^(1/2)V^Ty`; their radii are at most
`sqrt(C)n^(3/4)`. Theorem 18.3, with mesh proportional to
`epsilon n^(3/4)`, gives (18.31) for the bounded retained-feature query
language. Summing (18.34) over physical edges proves (18.32), and
concatenating incident singular features proves (18.33). `square`

The conclusion is scale-sensitive. A bounded-row-and-column-degree bridge is
subleading at `n^(3/2)`, while an iid dense sign bridge has extensive
numerical rank for thresholds inside its spectral bulk. The general
`2^(Omega(r))` lower bound of Theorem 18.3 applies to arbitrary featured
roofs; optimality for every Boolean SVD port is not claimed.

### Theorem 18.8 (bounded-operator sign bridges have extensive visible rank)

Let `R in {+-1}^{n by n}` satisfy
`||R||_(2->2)<=C sqrt(n)`, with `C>=1`, and put

```math
r_epsilon(R)=#{j:sigma_j(R)>epsilon sqrt(n)},
\qquad 0<=epsilon<1.                                 \tag{18.35}
```

Then

```math
r_epsilon(R)>= n{1-epsilon^2\over C^2-epsilon^2}.    \tag{18.36}
```

Also every rank-`r` dense sign matrix obeys

```math
||R||_(2->2)>={n\over sqrt r}.                       \tag{18.37}
```

#### Proof

The sign alphabet gives

```math
n^2=||R||_F^2=sum_j sigma_j(R)^2.                    \tag{18.38}
```

The `r_epsilon` visible terms are at most `C^2n`; all others are at most
`epsilon^2n`.  Rearranging
`n^2<=r_epsilon C^2n+(n-r_epsilon)epsilon^2n` proves (18.36).
If only `r` singular values are nonzero, (18.38) is at most
`r||R||^2`, proving (18.37). `square`

Consequently the SVD/operator certificate of Theorem 18.7 cannot interpolate
from fixed rank to a bounded-operator dense sign bridge through a
subextensive interface at fixed `n^(3/2)` accuracy: it must retain a linear
number of singular features.  Hadamard matrices attain (18.36) at `C=1`,
and the all-ones matrix attains (18.37).

This is a scoped stable-rank barrier, not a contextual-information lower
bound.  Full-rank permutation-invariant bridges can still compress by
Theorem 18.6, while the rank-one all-ones matrix compresses by magnetization.
Only synchronization or another nonlinear congruence can beat the spectral
interface in a bounded-operator sign family.  The audited benchmark proof is
in [`drafts/bounded_operator_rank_barrier.md`](drafts/bounded_operator_rank_barrier.md).

## 19. Extremal cut-norm replacement

### Theorem 19.1 (uniform labeled replacement)

Let `A,B` be real matrices on the same finite labeled interface and define the
unnormalized cut norm by

```math
||A-B||_square=max_(S,T)|(A-B)(S,T)|.
```

For `q` labels, `J in R^(q times q)`, and an arbitrary conditional future
`F:[q]^V->R`, put

```math
M(A;J,F)=max_sigma\left\{
F(sigma)+sum_(u,v)A_(uv)J_(sigma_u,sigma_v)\right\}.
```

Then

```math
|M(A;J,F)-M(B;J,F)|
<=||J||_1||A-B||_square.                                      \tag{19.1}
```

#### Proof

For fixed `sigma`, its label fibres `V_i` partition the interface, and the
energy difference is

```math
sum_(i,j)J_(ij)(A-B)(V_i,V_j).
```

This is bounded by the right side of (19.1). Apply the pointwise bound to an
optimizer in each direction. `square`

The future may pin a single rare labeling, so this is genuinely an extremal
replacement theorem. Frieze--Kannan weak regularity gives, for bounded dense
matrices, block representatives with

```math
k<=2^(O(epsilon^(-2))),
\qquad ||A-B||_square<=epsilon n^2.                            \tag{19.2}
```

Therefore a finite block/occupancy state preserves every fixed-label future
response at the dense `n^2` scale. Independent block rounding adds only
`O(n^(3/2))` cut error for fixed block data, supplying finite realizations.

The scale qualification is essential. Preservation at a declared scale
`L_n` needs cut error `o(L_n)`. At `L_n=n^(3/2)`, the generic weak-regularity
bound requires `epsilon=o(n^(-1/2))` and may have exponentially many blocks.
Thus Theorem 19.1 does not by itself compress the motivating signing scale.

The criterion has an explicit converse witness up to a universal constant.
For symmetric `D=A-B`, an Alon--Naor cut-norm rectangle `S,T` can be encoded
by four labels `(1_S,1_T)` and the pair reward

```math
J_((a,b),(c,d))=(ad+bc)/2.                                   \tag{19.3}
```

A sufficiently strong labeled unary future pins this configuration, making
the optimized response gap exactly `|D(S,T)|`. Hence a residual of order
`L_n` is an observable, scalable falsifier of replacement at that scale.

## 20. When hidden witness phase becomes response information

Let finite max-plus generators satisfy `T_e(i,j)<=lambda_e`.  Choose good
relations `R_e`, tolerances `theta_e>=0`, and a gap `C>0` so that

```math
(i,j) in R_e
  ==>T_e(i,j)>=lambda_e-theta_e,

(i,j) notin R_e
  ==>T_e(i,j)<=lambda_e-C.                                  \tag{20.1}
```

For a root set `A`, let `A_u=R_u(A)`, put
`Theta(u)=sum_s theta_(e_s)`, and define

```math
V_A(u;f)=max_(i in A,j in I)\{T_u(i,j)+f(j)\}
          -sum_slambda_(e_s).                               \tag{20.2}
```

### Theorem 20.1 (gap filtering and orbit-profile isometry)

If `A_u` is nonempty and `Theta(u)+osc(f)<C`, then

```math
max_(j in A_u)f(j)-Theta(u)
<=V_A(u;f)<=max_(j in A_u)f(j).                             \tag{20.3}
```

In particular, at `theta_e=0` equality holds.  If a finite group `G` acts by
permutations on `I`, `r:I->R`, and

```math
Phi_r(B)(g)=max_(i in B)r(gi),                              \tag{20.4}
```

then, for `k osc(r)<C`,

```math
sup_(g in G)|V_A(u;kr_g)-V_A(v;kr_g)|
=k||Phi_r(A_u)-Phi_r(A_v)||_infinity.                       \tag{20.5}
```

With nonzero tolerances, the two sides differ by at most
`Theta(u)+Theta(v)` under the corresponding gap condition.

#### Proof

An all-good path has weight between the declared toll minus `Theta(u)` and
the toll.  A path using a bad edge loses at least `C`, more than the total
good-path uncertainty plus `osc(f)`, so it cannot optimize.  This gives
(20.3); apply it to every orbit probe `kr_g` to get (20.5). `square`

Now let phases `q` carry nonempty backward supports `K_q`, forward envelopes
`B_q`, and a deterministic update `delta`, satisfying

```math
K_(delta(q,e)) subseteq R_e(K_q),

K_q subseteq B_q,
\qquad R_e(B_q) subseteq B_(delta(q,e)),                    \tag{20.6}
```

and suppose every declared orbit probe is saturated:

```math
max_(i in K_q)r_g(i)=max_(i in B_q)r_g(i).                  \tag{20.7}
```

### Theorem 20.2 (observable support sandwich)

From root `A=K_(q_0)`, every history `u` obeys

```math
K_(delta(q_0,u)) subseteq A_u subseteq B_(delta(q_0,u)),

Phi_r(A_u)=Phi_r(K_(delta(q_0,u))).                          \tag{20.8}
```

Thus `M` reachable phases with pairwise `gamma`-separated profiles give an
intrinsic `M`-state rooted-response packing at error below `kgamma/2`, under
the gap hypothesis of Theorem 20.1.

#### Proof

The first and second inclusions propagate respectively by the two halves of
(20.6).  Monotonicity of maxima and (20.7) force equality of the profiles;
(20.5) and the merged-state triangle inequality prove the packing. `square`

This is not arbitrary lookup exposure.  In the deterministic de Bruijn shift
on `I=E^m`, cyclic coordinate rotations and the one-symbol readout
`r(s)=s_1/(q-1)` separate all `q^m` singleton phases.  For nontight cost
`-C`, every integer `k<C` gives response separation at least `k/(q-1)`.
With hard (`-infinity`) leakage the separation may be amplified indefinitely.

The finite-gap qualification is necessary.  For length-`m` de Bruijn
histories, their rooted endpoint vectors differ by at most `C` in sup norm.
Every common max-plus continuation and terminal maximum is nonexpansive, so
no common future can separate them by more than `C`.  Anticipatory support
memory becomes semantic information only through a phase-preserving filter,
hard reachability, bounded probe horizon, or an equivalent observable
sandwich.

## 21. Positive and negative endpoints of the bridge hierarchy

### Theorem 21.1 (dense sign bridges carry target-scale response bits)

There are constants `a,gamma,L>0` such that, for every sufficiently large
`n`, some `B in {-1,1}^(n times n)` with
`||B||_(2->2)=O(sqrt n)` has `N>=exp(gamma n)` state--query pairs
`(x_c,y_c)`.  For each `sigma in {0,1}^N`, define

```math
h_sigma(x_c)=a n^(3/2)sigma_c,
\qquad
h_sigma(x)=-Ln^(3/2)\quad(x notin \{x_c\}).                 \tag{21.1}
```

Then

```math
(P_Bh_sigma)(y_c)=D_c+a n^(3/2)sigma_c,                     \tag{21.2}
```

with `D_c` independent of `sigma`.  The `2^N` response functions are pairwise
exactly `a n^(3/2)` apart.  Uniform error below half this amount needs
`2^(exp(Omega(n)))` states, or `exp(Omega(n))` bits.

Even the linear landscapes

```math
h_c(x)=-x^TBy_c                                             \tag{21.3}
```

contain an `exp(gamma n)`-element projective response packing at separation
`Omega(n^(3/2))`, and therefore need `Omega(n)` bits.

#### Proof

Take a random sign matrix and independent random sign queries, and set
`x_c=sign(By_c)`.  Standard Rademacher concentration and a sufficiently small
`gamma` give simultaneously

```math
||B||_(2->2)<=C_0sqrt n,
\quad ||By_c||_1>=d_0n^(3/2),

x_d^TBy_c<=d_1n^(3/2)\quad(d!=c),                           \tag{21.4}
```

for constants `d_1<d_0`.  Conditional on `B,y_d`, the last tail is Hoeffding
with coefficient norm at most `C_0n`; a union bound over `N^2` pairs is still
exponentially small.  Choose `a<(d_0-d_1)/2` and `L` larger than the
`infinity->1` bridge bound.  This proves (21.2), while sup-norm
nonexpansiveness proves exact separation.

After also requiring the random queries to have linear pairwise Hamming
distance, the same rowwise concentration gives
`||B(y_c-y_d)||_1>=cn^(3/2)`.  But

```math
(P_Bh_c)(y)=||B(y-y_c)||_1.                                 \tag{21.5}
```

At `y_c,y_d`, a pairwise response difference takes both signs with this
magnitude, proving the projective linear-child packing. `square`

Thus neither density nor spectral control compresses arbitrary child
landscapes at the motivating scale.  The theorem does not address the much
smaller class of fixed-magnitude quadratic near-minimizers.

### Theorem 21.2 (multitype common-section optimization)

Partition a common coordinate set into types `C_c`, `|C_c|=n_c`.  Block `a`
has an internal energy depending only on the typewise plus counts `k_(a,c)`.
Couple blocks by

```math
R_ab=sum_c alpha_(ab,c)D_c
     +sum_(c,d)beta_(ab,cd)1_(C_c)1_(C_d)^T,
\qquad alpha_(ab,c)>=0.                                     \tag{21.6}
```

Then the exact joint optimum is

```math
max_(k_(a,c))\left\{
sum_a h_a(k_a)+sum_(a<b)\left[
sum_c alpha_(ab,c)(n_c-2|k_(a,c)-k_(b,c)|)
+sum_(c,d)beta_(ab,cd)s_(a,c)s_(b,d)
\right]\right\}.                                           \tag{21.7}
```

One nested representative inside each type attains every pairwise bound at
once.  The reduced joint search uses at most `mKlog_2(n+1)` label bits for
`m` blocks and `K` types, even though the matrices can be dense and full
rank.

For signed nonzero `alpha_(ab,c)` and `n_c>=2`, a common section exists
uniformly for all count assignments iff each type's signed block graph is
balanced.  An isolated unbalanced unit cycle at even `n_c` and zero
type-magnetization loses exactly `2n_c` relative to separately optimized pair
responses.

#### Proof

Inside one type,

```math
sum_(i in C_c)x_i^ax_i^b
=n_c-2|P_(a,c)\triangle P_(b,c)|
<=n_c-2|k_(a,c)-k_(b,c)|.                                  \tag{21.8}
```

Common initial segments attain equality for every pair.  Signed balance is
the vertex-gauge criterion.  For necessity on an unbalanced cycle, gauge a
spanning path and assign count one or `n_c-1` according to the path gauge;
path-edge optimality forces one fixed singleton or its complement, while the
closing edge has the incompatible sign.  The quantitative cycle loss follows
coordinatewise as in Theorem 18.6. `square`

This is a reduced **jointly reoptimizable** search, not an arbitrary-future
state for a fragment with a frozen microscopic alignment.  Under fixed
`m,K`, convergent type proportions, `beta=b/n`, and uniformly convergent
gauged internal energies, (21.7) yields a finite-dimensional compact
variational limit by uniform convergence on the type-magnetization grids.

### Theorem 21.3 (vertex-cover bridge width is worst-case sharp)

Let a bipartite bridge `R` have left variables `A union U`, right variables
`B union V`, where `A union B` is a vertex cover of size `k`; hence
`R_(UV)=0`.  For an arbitrary left landscape `F`, define its conditioned
table

```math
Q_F(a,beta)=a^TR_(AB)beta
+max_u\{F(a,u)+u^TR_(UB)beta\}.                              \tag{21.9}
```

It has `2^k` entries and determines the exact outgoing bridge response:

```math
(P_RF)(beta,v)=max_a\{Q_F(a,beta)+a^TR_(AV)v\}.              \tag{21.10}
```

Thus equality of this envelope is the coarsest exact state for arbitrary
right futures; the raw table is a universal sufficient representation.

This dependence is worst-case sharp.  For a `k`-edge matching of weight `W`,
let normalized child tables range over `[-D,0]^({-1,1}^k)` with `2W>D`.
Then

```math
(P_RQ)(y)=kW+Q(y).                                           \tag{21.11}
```

Every coordinate is exposed.  The exact projective lattice with mesh `eta`,
`M=D/eta`, has

```math
(M+1)^(2^k)-M^(2^k)                                         \tag{21.12}
```

classes, and the worst-case response rate is

```math
Theta(2^k log(1+D/epsilon))                                  \tag{21.13}
```

bits.

#### Proof

Because no edge joins `U` to `V`, conditioning on `(a,beta)` separates the
left maximization and gives (21.10).  Arbitrary future weights compare only
the displayed envelope and can pin each right spin assignment.  In the
matching instance, changing one spin from `x=y` loses `2W>D`, proving
(21.11).  Sup-norm response isometry, grid packing, and coordinate
quantization give (21.12)--(21.13). `square`

The matching has degree and treewidth one but minimum vertex cover `k`.
Sparse compression is controlled by simultaneous live interface size, not
local degree.  Fixed numerical degeneracies of `R` may quotient the table
further; `2^k` is the sharp support-only law.

### Theorem 21.4 (complete sign quadratics retain extensive bridge information)

There are universal `gamma,g>0` such that, for every sufficiently large
`n`, some sign bridge `B in {-1,1}^(n times n)` with
`||B||_(2->2)=O(sqrt n)` and complete sign quadratics

```math
H_c(x)=\sum_{i<j}z_(c,i)z_(c,j)x_ix_j
      ={(x^Tz_c)^2-n\over2}                                \tag{21.14}
```

have `N>=exp(gamma n)` pairwise projectively separated responses:

```math
{1\over2}\operatorname{osc}(P_BH_c-P_BH_d)
\ge g n^{3/2}\qquad(c\ne d).                              \tag{21.15}
```

Consequently every state that answers arbitrary later continuations on the
full class of complete sign quadratics to error smaller than
`g n^(3/2)/2` requires `exp(gamma n)` states, or `Omega(n)` bits, even after
one additive calibration per child.

#### Proof

The deterministic ingredient is pole locking.  If
`||h||_infinity<n/2`, then

```math
\max_x\{((x^Tz)^2-n)/2+h^Tx\}
={n^2-n\over2}+|z^Th|.                                    \tag{21.16}
```

Indeed, choose the nearer pole `sz`, at Hamming distance `r<=n/2` from
`x`.  The quadratic loss is `2r(n-r)>=nr`, while the field can recover less
than `2r(n/2)=nr`.

Take a random sign bridge and independent random sign queries `y_c`, put
`h_c=By_c` and `z_c=sign(h_c)`.  For a sufficiently small fixed `gamma`,
standard Rademacher tails and a union bound over `exp(2gamma n)` queries and
their ordered pairs give simultaneously

```math
||B||_(2->2)<=C_0sqrt n,
\quad ||h_c||_1>=d_0n^{3/2},
\quad ||h_c||_infinity<n/2,
\quad |z_d^Th_c|<={d_0\over2}n^{3/2}\ (c\ne d).             \tag{21.17}
```

For the last bound, condition on `B,y_d`; then `z_d` is fixed and
`z_d^TBy_c` is a Rademacher sum whose coefficient norm is at most `C_0n`.
Equation (21.16) makes the diagonal response exceed every off-diagonal
response in its query coordinate by `(d_0/2)n^(3/2)`.  Reversing `c,d`
gives both signs, so half the oscillation has the same lower bound. `square`

This theorem is a genuine quadratic strengthening of Theorem 21.1, but its
scope is essential: (21.14) has cap `(n^2-n)/2` and projective spread
`Theta(n^2)`.  It does not address the bounded-cap or near-minimizing
quadratics relevant to the motivating signing problem.

### Theorem 21.5 (sparse surrogates and the coefficientwise ceiling)

Let `m=binom(n,2)`, fix any bridge `B`, and let `0<epsilon<=1`. For
`n>=64/epsilon^2`, every complete sign quadratic has a sparse weighted
surrogate with coefficients in

```math
\{-1/q,0,1/q\},\qquad q=1-\epsilon^2/2,                   \tag{21.18}
```

which is uniformly within `epsilon n^(3/2)` before and after every bridge
or later maximization. One universal list of at most `m+2` support masks
suffices for all `2^m` inputs, and the summary uses at most

```math
\left\lceil(1-\epsilon^2/4)m\right\rceil
+\left\lceil\log_2(m+2)\right\rceil                       \tag{21.19}
```

bits. Thus a strict constant fraction of the exact coefficient information
can be discarded at the target scale. The decoder's surrogate is
`2`-bounded and sparse, not itself a sign quadratic.

If code centers must remain complete sign quadratics, set
`E=epsilon n^(3/2)`, `r=floor(E/2)<=m/2`. There is an internal codebook of
size at most

```math
K\le\left\lceil {2^m\over\sum_{j=0}^r{m\choose j}}
                  (m\log2+1)\right\rceil                  \tag{21.20}
```

whose responses approximate every complete sign-quadratic response within
`E`.  For fixed `epsilon>0`, its bit count is

```math
m-\log_2\sum_{j=0}^r{m\choose j}+O(\log m)
\le m-{\epsilon\over4}n^{3/2}\log_2n+O_epsilon(n^{3/2}).   \tag{21.21}
```

The coefficient-Hamming architecture is optimal up to `O(log n)` bits:
the sphere-covering inequality gives the reverse bound for covers by
radius-`r` coefficient balls.  This is not a response-entropy lower bound,
because Boolean maximization may identify Hamming-distant coefficients.

For the larger class `a_ij in [-1,1]`, there is a simultaneous unbiased
rounding to a grid with `O(1/epsilon)` levels such that

```math
\sup_x|H_A(x)-H_(\widehat A)(x)|
\le\epsilon n^{3/2}.                                      \tag{21.22}
```

Hence `O(n^2 log(1+1/epsilon))` bits suffice uniformly over `B` and every
later continuation.

#### Proof

For the sparse surrogate, keep each coefficient independently with
probability `q` and divide retained coefficients by `q`. For any fixed input
signing and spin, the error is a sum of centered variables bounded by one
with total variance at most `epsilon^2n^2/2`. Bernstein and a union bound over
the `2^n` spins show that one mask is uniformly accurate with probability
greater than `3/4`. Chernoff shows with probability greater than `3/4` that
it erases at least `(epsilon^2/4)m` coefficients. Thus a fixed coefficient
signing has a good mask with probability at least `1/2`. Sampling `m+2`
independent masks and union-bounding over all `2^m` coefficient signings
leaves positive probability that the one list covers every input. This gives
(21.19); max-type response operators are sup-norm nonexpansive.

The Hamming cube has a radius-`r` covering of size (21.20) by the
probabilistic covering argument.  Changing `r` signs changes every energy by
at most `2r`, and maximization is sup-norm nonexpansive.  The binomial lower
bound `binom(m,r)>=(m/r)^r` gives (21.21), while the sphere-covering bound
proves the architecture-specific converse.

For (21.22), round each coefficient independently and without bias to one of
its two adjacent grid points of spacing `delta<=epsilon`.  For fixed `x`,
Hoeffding bounds the rounding error by
`delta sqrt(m(n+2)log(2)/2)` except with probability below `2^(-n-1)`.
A union bound over all `2^n` spins leaves positive probability and the
displayed quantity is below `epsilon n^(3/2)`.  Sup-norm nonexpansiveness
again transfers the estimate through the bridge and all futures. `square`

Together, Theorems 21.4--21.5 leave a real information gap:

```math
Omega(n)\le R_{\rm sign\ quadratic}
             (epsilon n^{3/2})\le O(n^2).                  \tag{21.23}
```

Closing it requires optimizer geometry shared by Hamming-distant
coefficients; coefficientwise Lipschitz control, switching, and generic cut
regularity cannot do so.

### Theorem 21.6 (exact cap-`1/2` children still carry growing response state)

On the subsequence `n=2^(2m)`, there is one dense sign bridge `B_n` with
`||B_n||_(2->2)=sqrt n` and `N>=exp(c sqrt n)` hollow symmetric signings
`A_1,...,A_N` such that

```math
Q(A_c)=\max_x|H_(A_c)(x)|={1\over2}n^{3/2}                \tag{21.24}
```

and

```math
{1\over2}\operatorname{osc}(P_(B_n)H_(A_c)-P_(B_n)H_(A_d))
\ge {1\over8}n^{3/2}\qquad(c\ne d).                       \tag{21.25}
```

Hence the cap-`1/2` class needs `exp(Omega(sqrt n))` states, or
`Omega(sqrt n)` bits, for error below `n^(3/2)/16` under arbitrary
coordinate-pinning continuations. The explicit family itself has a matching
`O(sqrt n)`-bit description.

#### Construction and proof

Put `q=2^m`, `n=q^2`, and index coordinates by
`(u,v) in F_2^m times F_2^m`. Let `W` be the Walsh matrix and put

```math
b(u,v)=(-1)^(u dot v),
\qquad mathcal H=D_bWD_b.                                  \tag{21.26}
```

Then `W^2=nI`, `Wb=qb`, and

```math
mathcal H^2=nI,
\quad mathcal H1=q1,
\quad tr(mathcal H)=0.                                     \tag{21.27}
```

Thus `A=mathcal H-diag(mathcal H)` is a hollow sign matrix and, on Boolean
vectors, `H_A(x)=x^Tmathcal Hx/2`. The spectral bound and the all-one
eigenvector prove `Q(A)=qn/2`.

For each Boolean table `g:F_2^m->F_2`, define

```math
s_g(u,v)=(-1)^(u dot v+g(v)),
\qquad y_g=q^(-1)Ws_g.                                    \tag{21.28}
```

Direct Walsh summation gives

```math
(Ws_g)(a,b)=q(-1)^(a dot b+g(a)),                          \tag{21.29}
```

so `y_g` is a sign vector and `Wy_g=qs_g`. A random code supplies at least
`exp(q/32)` tables with pairwise bias

```math
|S(g,h)|=\left|\sum_v(-1)^(g(v)+h(v))\right|\le q/2.      \tag{21.30}
```

Indeed, each pair violates (21.30) with probability at most `2e^(-q/8)`,
and the pair union bound succeeds. Switch the child by
`A_g=D_(s_g)AD_(s_g)` and use the common bridge `B=W`.

For `w=s_g odot s_h`, the exact rooted Rayleigh coordinate is

```math
w^Tmathcal Hw=qS(g,h)^2,
\qquad rho={w^Tmathcal Hw\over qn}\le1/4.                 \tag{21.31}
```

At query `y_h`, changing variables gives

```math
(P_WH_(A_g))(y_h)
=\max_u\{u^Tmathcal Hu/2+q w^Tu\}.                        \tag{21.32}
```

For the diagonal `g=h`, `u=1` gives `3qn/2`. Off the diagonal, set
`K=2qI-mathcal H`. Completing the square and using
`K^(-1)=(2qI+mathcal H)/(3q^2)` gives, even on the containing Euclidean
sphere,

```math
u^Tmathcal Hu/2+qw^Tu
\le qn+{2qn+w^Tmathcal Hw\over6}
\le {11\over8}qn.                                         \tag{21.33}
```

Thus child `g` beats child `h` at `y_g` by `qn/8`, and the reverse holds at
`y_h`; (21.25) follows. `square`

This theorem closes the bounded-cap qualitative gap left by Theorem 21.4.
It does not prove a linear information rate, apply at every order, or say
that exact minimizers below the conference-scale cap share this complexity.
It is a response lower bound for the declared arbitrary-pinning future
family, not for futures restricted only to hollow sign-quadratic fragments.

### Theorem 21.6a (explicit permutation-Walsh rate)

On the same subsequence `n=q^2=2^(2m)`, there are the deterministic Walsh
bridge `W` and a family of exact-cap switched children of size

```math
\exp(\Omega(q\log q))=\exp(\Omega(\sqrt n\log n))        \tag{21.33a}
```

whose projective response distance is at least `n^(3/2)/8`. The family has a
matching `Theta(sqrt(n)log n)`-bit description.

#### Proof

For every permutation `pi:F_2^m->F_2^m`, use the Maiorana--McFarland switch

```math
s_\pi(u,v)=(-1)^(u dot \pi(v)).                          \tag{21.33b}
```

It is bent, since direct summation gives

```math
(Ws_\pi)(a,b)=q(-1)^(b dot \pi^{-1}(a)).                 \tag{21.33c}
```

Thus `y_pi=q^(-1)Ws_pi` is Boolean. Switch the regular child in Theorem 21.6
by `s_pi`. For a pair `pi,sigma`, put

```math
\tau(v)=v+\pi(v)+\sigma(v),
\qquad
\rho(\pi,\sigma)
=E_(x,y)(-1)^(x dot y+\tau(x) dot \tau(y)).              \tag{21.33d}
```

Exact Walsh summation identifies `rho` with the normalized rooted Rayleigh
coordinate `w^Tmathcal Hw/(qn)`, where `w=s_pi odot s_sigma`.

It remains to find many pairs with `rho<=1/4`. First take `tau` to be a
uniform random function. With

```math
f_a(y)=(-1)^(a dot \tau(y)),
```

one has

```math
\rho(\tau)=q^{-1}\sum_x\widehat f_(\tau(x))(x).          \tag{21.33e}
```

If this exceeds `1/4`, more than `q/7` summands exceed `1/8`. Parseval allows
fewer than `64` such frequencies for one value of `tau(x)`, so the good
image spans at least `r=m-8` independent output characters. For fixed
independent characters and witnesses, the resulting `rq` Fourier signs are
independent Rademachers. Hoeffding and a union bound give

```math
Pr\{\rho(\tau)>1/4\}
\le q^{2r}\exp(-rq/128).                                \tag{21.33f}
```

For independent random functions `pi,sigma`, the values of `tau` are iid
uniform. Conditioning both functions to be permutations costs at most
`(q^q/q!)^2<=e^(2q)`. Turan's reciprocal-density bound on the resulting bad-
pair graph therefore supplies a permutation code of size at least

```math
\exp(rq/128-2r\log q-2q)=\exp(\Omega(q\log q)).          \tag{21.33g}
```

For every good pair, the same resolvent completion as in (21.33) puts the
cross response at most `11qn/8`, while the matched response is `3qn/2`.
Evaluating in both query directions gives projective gap `qn/8`. Listing the
permutation costs `log(q!)=Theta(q log q)` bits and reconstructs the child,
which proves the matching family-level upper. `square`

This explicit theorem is dominated in rate by Corollary 21.8a, whose bridge
and code are probabilistic. Its additional structural content is the
approximate-inner-product-isometry tail (21.33f).

### Lemma 21.7 (switching response equals weighted near-top deficit)

Let `P=max_u H_A(u)` and fix a top state `u_*`. For a field `h`, define

```math
Delta_A(h)=P+||h||_1-\max_u\{H_A(u)+h^Tu\}.                \tag{21.34}
```

Then the exact variational identity is

```math
Delta_A(h)=\min_u\left\{P-H_A(u)
+2\sum_(i:u_i ne sign(h_i))|h_i|\right\}.                 \tag{21.35}
```

For every bridge query `y`, put
`s_y=u_* odot sign(By)` and `A^(s_y)=D_(s_y)AD_(s_y)`. Its matched response
attains the joint roof:

```math
(P_BH_(A^(s_y)))(y)=P+||By||_1.                            \tag{21.36}
```

For two queries `y,z`, their switched children obey

```math
d_proj(A^(s_y),A^(s_z))
\ge {1\over2}\{Delta_A(s_z odot By)+Delta_A(s_y odot Bz)\}.\tag{21.37}
```

#### Proof

Changing variables under switching turns the response into
`max_u(H_A(u)+(s odot By)^Tu)`. Expanding a linear score as its `l_1` roof
minus twice the weighted Hamming disagreement proves (21.35). For `s_y`,
the field is `u_* odot |By|`, so `u_*` maximizes both terms and gives
(21.36). At query `y`, the response difference between its matched child and
the `z` child is the first deficit in (21.37); at query `z` it is the
negative of the second. Half the oscillation gives (21.37). `square`

Thus an `Omega(n)`-bit cap-bounded packing reduces to a concrete geometric
statement: exponentially many linked fields must avoid the weighted Hamming
neighborhoods of the near-top set. Exact enumeration supports this at orders
through fourteen but is not used as asymptotic evidence.

### Theorem 21.8 (near-top entropy amplifies to contextual information rate)

Let `H:{-1,1}^n->R`, let `P=max_uH(u)=H(u_*)`, and define the switching
orbit by `H^s(x)=H(s odot x)`. Fix `d_0,kappa>0`. If

```math
#\{u:P-H(u)<d_0n^{3/2}\}
\le \exp\{(\log2-\kappa)n\},                              \tag{21.38}
```

then there are constants `C,gamma,d>0`, depending only on `d_0,kappa`,
such that, for all sufficiently large `n`, one sign bridge `B` obeys

```math
||B||_(2->2)\le C\sqrt n                                \tag{21.39}
```

and at least `exp(gamma n)` switched responses satisfy

```math
d_proj(P_BH^s,P_BH^t)\ge d n^{3/2}\qquad(s\ne t).       \tag{21.40}
```

In particular, let an encoder map these switched children to `K` states and
let one common decoder map `(state,y)` to a predicted response. Uniform
error `epsilon n^(3/2)` at every Boolean query, with `epsilon<d/2`, forces

```math
K\ge\exp(\gamma n).                                     \tag{21.41}
```

Thus positive near-top entropy deficit implies positive contextual
response-information rate. Since the switch itself is an `n`-bit label,
the constructed switching family has `Theta(n)` response bits at this scale,
up to constants.

#### Proof

For every row tolerance `eta>0`, there are `a,rho>0` such that, uniformly
for `|y^Tz|<=rho n`, a random sign row `R`, and either `t in {+-1}`,

```math
Pr\{|R^Ty|\ge a\sqrt n,
 sign(R^Ty)sign(R^Tz)\ne t\}\ge {1\over2}-\eta.          \tag{21.42}
```

Indeed, the normalized pair has covariance
`Sigma_theta=((1,theta),(theta,1))`, where `theta=y^Tz/n`. For its independent
summands `X_i`,

```math
\sum_i E||\Sigma_\theta^{-1/2}X_i||_2^3
\le {2^{3/2}\over(1-\rho)^{3/2}\sqrt n}.                \tag{21.43}
```

The convex-set Lyapunov bound of Bentkus therefore gives a uniform
`O_rho(n^(-1/2))` bivariate Gaussian approximation. The event in (21.42) is
the disjoint union of two convex quadrant half-strips. The Gaussian arcsine
law makes either sign occur with probability at least
`1/2-arcsin(rho)/pi`, and deleting `|R^Ty|<a sqrt(n)` costs `Pr{|G|<a}`.
Choosing `rho,a` small and then `n` large proves (21.42).

Choose `alpha=1/2-eta` sufficiently close to `1/2`, and then
`delta in (0,alpha)` sufficiently small, that

```math
D(\delta||\alpha)>\log2-{\kappa\over2}.                 \tag{21.44}
```

Shrink `a` so `d=2a delta<d_0`. A random Boolean code gives
`Y subset {-1,1}^n`, `|Y|>=exp(gamma n)`, with
`|y^Tz|<=rho n` for distinct words, where
`gamma<min(rho^2/8,kappa/8)`.

Choose `B` with iid sign entries and put

```math
s_y=u_* odot sign(By).                                  \tag{21.45}
```

For a fixed ordered pair `y!=z` and a spin in the set (21.38), the row
events (21.42), with the target sign supplied by that spin, are independent
and have probability at least `alpha`. The binomial lower-tail bound makes
the probability of fewer than `delta n` successful rows at most
`exp[-D(delta||alpha)n]`. Union-bounding over the near-top set and all
ordered query pairs leaves probability at most `exp(-kappa n/4)`. The
standard iid-sign norm tail intersects this event with (21.39).

Every successful row contributes at least `a sqrt(n)` of weighted mismatch.
The exact identity (21.35) therefore gives a cross deficit at least
`d n^(3/2)`: outside the near-top set the energy deficit gives it, and inside
the set twice the weighted mismatch gives it. Both ordered directions hold.
Equation (21.37) now proves (21.40). If two children shared a decoded state,
their common sup-error approximation would put their projective distance at
most `2 epsilon n^(3/2)`, proving (21.41). `square`

The probability inputs are the Rademacher specialization of
[Rudelson--Vershynin's Hanson--Wright theorem](https://doi.org/10.1214/ECP.v18-2865)
and [Bentkus's independent-summand Lyapunov bound](https://doi.org/10.1137/S0040585X97981123).
The latter is applied after the explicit covariance normalization (21.43),
not as an iid black box.

### Corollary 21.8a (exact cap `1/2` has linear response rate)

For every sufficiently large `n=2^(2m)`, there are one sign bridge of
operator norm `O(sqrt n)` and `exp(gamma n)` hollow signings, all switchings
of one regular-Walsh child, such that

```math
Q(A_s)={1\over2}n^{3/2},
\qquad
d_proj(P_BH_(A_s),P_BH_(A_t))\ge d n^{3/2}.             \tag{21.46}
```

Hence the arbitrary-coordinate-pinning response complexity of this exact-cap
family is `Theta(n)` bits at fixed target-scale error.

For the Walsh matrix `mathcal H` in Theorem 21.6, choose `d_0=1/8`.
Membership in the near-top set implies
`u^Tmathcal Hu>(3/4)n^(3/2)`. Since

```math
E u^Tmathcal Hu=0,
\quad ||mathcal H||_F=n,
\quad ||mathcal H||_(2->2)=\sqrt n,
```

Hanson--Wright gives (21.38) with an absolute `kappa>0`; Theorem 21.8
applies. This strictly strengthens the explicit `Omega(sqrt n)` subfamily in
Theorem 21.6, while the latter retains the advantage of a deterministic
Walsh bridge and an explicit matching small state.

More generally, if symmetric `K_n` obeys

```math
||K_n||_F\le Ln,
\quad ||K_n||_(2->2)\le L\sqrt n,
\quad
\max_u{u^TK_nu\over2}-E_U{U^TK_nU\over2}
\ge p n^{3/2},                                           \tag{21.47}
```

then its switching orbit has positive contextual response rate, with
constants depending only on `p,L`. Thus the amplification law is not tied
to Walsh algebra.

### Theorem 21.9 (exact Walsh graph carrier)

Let `q=2^m`, `n=q^2`, and let `R` be the order-`q` Walsh matrix. For a truth
table `g:F_2^m->F_2`, put

```math
C_g=R\otimes(D_gRD_g),
\qquad H_g(x)={1\over2}x^TC_gx.                          \tag{21.48}
```

This is exactly the hollow switched-child energy from Theorem 21.6, since
`tr(C_g)=0`, and `Q(C_g)=n^(3/2)/2`. If a graph `G` joins `k` such blocks by
the common bridge `W=R tensor R`, define the `kq`-square carrier

```math
(K_(G,g))_(ii)=D_(g_i)RD_(g_i),
\qquad
(K_(G,g))_(ij)=cases(R,&ij in E(G);0,&otherwise).        \tag{21.49}
```

After reordering `(i,u,v)` to `(u,i,v)`, the complete quadratic matrix is
exactly

```math
R\otimes K_(G,g).                                       \tag{21.50}
```

Thus graph composition has an exact coefficient-level presentation using
the graph and `kq=k sqrt(n)` truth-table bits. This is a strict description
reduction from a general `kn`-variable quadratic, but it is not yet an
extremal-response quotient: maximizing (21.50) still ranges over all
`kq^2` Boolean spins.

#### Proof

Writing the original switch as the bent factor times `(-1)^{g(v)}` cancels
the regularizing diagonal and gives (21.48). Diagonal blocks factor as
`R tensor (D_gRD_g)`, while every edge block is `R tensor R`; reordering
factors out the first `R` and proves (21.50). `square`

### Theorem 21.10 (composition exposes a Walsh commutation cocycle)

For `a in F_2^m`, let `g_a(v)=a dot v`, put `F=W/q`, and let
`widehat C_a=C_(g_a)/q`. Then

```math
F\widehat C_a=(-1)^(a dot a)\widehat C_aF.              \tag{21.51}
```

Choose nonzero `a_0,a_1` with `a_0 dot a_0=0` and
`a_1 dot a_1=1`. Their truth tables are balanced. For either constant word
on `k` blocks, every within-word pair correlation is `q` and every bias is
zero. Nevertheless, on every bipartite graph `G`,

```math
\max_XE_(G,a_0)(X)
=\left({k\over2}+|E(G)|\right)n^{3/2},                  \tag{21.52}
```

whereas

```math
\max_XE_(G,a_1)(X)
\le {k\over2}\sqrt{1+||A(G)||_(2->2)^2}\;n^{3/2}.      \tag{21.53}
```

For the path `P_k`, the missed amount is at least

```math
\left[{3k\over2}-1
-{k\over2}\sqrt{1+4\cos^2{\pi\over k+1}}\right]n^{3/2},\tag{21.54}
```

whose coefficient divided by `k` tends to `(3-sqrt5)/2`.

#### Proof

Walsh modulation and translation give (21.51). Both `F` and
`widehat C_a` are symmetric orthogonal involutions. In the even case they
commute. A Boolean `+1` child eigenvector `s` and its Boolean Walsh dual
`Fs` are both `+1` child eigenvectors; assigning them to the two color
classes saturates every child and every bridge, proving (21.52).

In the odd case they anticommute. For

```math
mathcal M_a=I_k tensor \widehat C_a+A(G) tensor F
```

the cross terms cancel, so

```math
mathcal M_a^2=(I_k+A(G)^2) tensor I_n.
```

Hence `||mathcal M_a||=sqrt(1+||A(G)||^2)`. A Boolean global vector has
squared Euclidean norm `kn`, and the Rayleigh bound, with the prefactor
`q/2`, proves (21.53). The path spectrum gives (21.54). `square`

This is a scalable composition-created-information counterexample. Biases
and all pairwise truth-table overlaps are not a reusable state, even though
they suffice for the earlier one-bridge scalar certificate. The missing
observable is a relative commutation cocycle between the flat edge transport
and the on-site involution. The result does not prove that the full
`sqrt(n)`-bit truth table is necessary under repeated composition.

### Theorem 21.11 (relative involution composition law)

Let `C,F` be real symmetric `N by N` involutions and, for a graph `G` on `k`
vertices, set

```math
\mathcal M_G=I_k\otimes C+A_G\otimes F,
\qquad E_G(X)={\lambda\over2}X^T\mathcal M_GX.           \tag{21.55}
```

Write `rho=||A_G||`, `eta_+=||CF+FC||`, and
`eta_-=||CF-FC||`. Then

```math
\max_{X\in\{\pm1\}^{kN}}E_G(X)
\le {\lambda kN\over2}
     \sqrt{1+\rho^2+\rho\eta_+}.                        \tag{21.56}
```

If `G` is bipartite with classes `L,R` and there is a Boolean `s` with
`Cs=s` and `Fs` Boolean, then

```math
\max_XE_G(X)\ge\lambda N\left(
 {k\over2}+|E(G)|-{\eta_-^2\over4}
 \sum_{H\in\operatorname{cc}(G)}
 \min\{|L\cap H|,|R\cap H|\}\right).                  \tag{21.57}
```

In particular, exact commutation gives the termwise optimum

```math
\max_XE_G(X)=\lambda N\left({k\over2}+|E(G)|\right).    \tag{21.58}
```

If `G` is `r`-regular bipartite with `r>0`, comparison of such a commuting
pair with a pair satisfying `||CF+FC||<=eta<2` gives the extensive gap

```math
{\lambda kN\over2}
\left(1+r-\sqrt{1+r^2+r\eta}\right)>0.                 \tag{21.59}
```

#### Proof

The exact square is

```math
\mathcal M_G^2=(I_k+A_G^2)\otimes I_N
 +A_G\otimes(CF+FC).                                    \tag{21.60}
```

The two summands commute, so in fact

```math
||\mathcal M_G||^2
=\max_{\alpha\in\operatorname{spec}(A_G),
       \sigma\in\operatorname{spec}(CF+FC)}
 (1+\alpha^2+\alpha\sigma).                            \tag{21.61}
```

The triangle bound in (21.61), followed by the Boolean Rayleigh bound
`||X||_2^2=kN`, proves (21.56). For the lower bound put `t=Fs`, assign `s`
and `t` to opposite color classes, and orient each connected component so
that `t` occurs on the smaller side. Every bridge is saturated. Moreover

```math
N-t^TCt={1\over2}||Ct-t||_2^2
       ={1\over2}||[C,F]s||_2^2.                        \tag{21.62}
```

Thus every `t`-child loses at most `lambda eta_-^2N/4`, proving (21.57).
When the commutator vanishes, all child and edge terms are saturated and
their separate norm bounds prove equality. Finally a regular bipartite
adjacency has eigenvalues `+-r`; pairing the sign with an extremal eigenvalue
of `CF+FC` makes (21.56) exact at the operator level and yields (21.59).
`square`

The theorem is a robust certificate, not a complete response state. Its
lower half also depends essentially on the transported pole `Fs` remaining
Boolean. Generic perturbations need not preserve that property.

### Theorem 21.12 (binary Gram data alone is not a Walsh extremal quotient)

Let `m>=3`, `q=2^m`, `n=q^2`, `F=W/q`, and

```math
\widehat C_a=D_aFD_a,
\qquad D_a(u,v)=(-1)^{a\cdot v}.                         \tag{21.63}
```

For the path on three blocks, give a label word `c=(c_1,c_2,c_3)` the energy

```math
E_c(X)={q\over2}\sum_i x_i^T\widehat C_{c_i}x_i
       +q x_1^TFx_2+q x_2^TFx_3.                        \tag{21.64}
```

There are words `c^-` and `c^+` with identical binary Gram matrices
`(c_i dot c_j)`, but

```math
\max E_{c^+}={7\over2}n^{3/2},
\qquad
\max E_{c^-}\le {3\sqrt3\over2}n^{3/2}.                \tag{21.65}
```

#### Proof

Take

```math
a=(1,1,1,0,\ldots,0),\qquad b=(0,0,1,0,\ldots,0),
```

and `c^-=(a,a,a)`, `c^+=(a,b,a)`. All pairwise binary inner products are
one. The Boolean sign of

```math
Q(u,v)=u_0u_1+v_0v_1+u_2v_2+v_0+v_1+v_2
       +\sum_{j\ge3}u_jv_j                              \tag{21.66}
```

is bent. Calling it `x`, direct two-bit Walsh factorization gives

```math
\widehat C_ax=x,\qquad \widehat C_bx=-x,
\qquad Fx\in\{\pm1\}^n,
\qquad \widehat C_b(Fx)=Fx.                             \tag{21.67}
```

Hence `(x,Fx,x)` saturates all five terms for `c^+`. For `c^-`, `F` and
`widehat C_a` anticommute, so the square of the normalized three-block
operator is `(I_3+A(P_3)^2) tensor I_n`; its norm is `sqrt(3)`, proving the
second bound. `square`

This rejects **Gram alone**. It does not reject truth-table overlaps: for
linear labels their overlap is `q 1_(a=b)`, which already separates these
two words.

### Theorem 21.13 (the characteristic root is response-visible)

Put `omega=(1,...,1) in F_2^m`. For an ordered label tuple define

```math
G(\mathbf a)=(a_i\cdot a_j)_{ij},\qquad
\mathcal R(\mathbf a)=\{c:\sum_i c_i a_i=0\},
\qquad
\mathcal R_\omega(\mathbf a)=\{c:\sum_i c_i a_i=\omega\}. \tag{21.68}
```

If `m>=3` is odd, use the singleton labels `omega` and `e_j`; if `m>=4` is
even, use `omega` and `e_1+e_2`. In either case the two labels have the same
`G` and `mathcal R` but different `mathcal R_omega`. Through the rooted Walsh
future from Theorem 21.6 their projective responses are separated by at
least `n^(3/2)/6`.

#### Proof

In the odd case both labels have norm one; in the even case both are nonzero
isotropic. They have no nonzero singleton relation, while only the first
equals `omega`. With `s_c(u,v)=(-1)^(u dot v+c dot v)` and
`y_c=q^(-1)Ws_c`, the matched response is `3n^(3/2)/2`. Distinct linear
truth tables have zero correlation, hence the rooted Rayleigh coordinate in
(21.31) is zero. The resolvent calculation (21.33), now with `rho=0`, bounds
the crossed response by `4n^(3/2)/3`. Reversing the queries gives the
projective gap. `square`

### Theorem 21.14 (all-dimensional rooted relation-form quotient)

For every `m>=1`, two ordered linear-label tuples have identical triples

```math
(G(\mathbf a),\mathcal R(\mathbf a),
  \mathcal R_\omega(\mathbf a))                         \tag{21.69}
```

if and only if one is carried to the other by a common
`O in O(m,2)`. Consequently (21.69) is an exact `O(k^2)`-bit quotient,
independent of `m`, for the Boolean maximum on every Walsh-bridged graph of
`k` blocks. It also preserves rooted responses after the same coordinate
relabeling.

#### Proof

Equality of relation kernels makes `a_i -> a_i'` a well-defined linear
isomorphism of their spans, and equality of `G` makes it an isometry. The
rooted coset says exactly that membership of `omega` agrees and that the
isometry fixes `omega` whenever it is present.

For odd `m`, `omega dot omega=1` and
`H=omega^perp` is nondegenerate alternating. Every vector decomposes
uniquely as

```math
u=(u\cdot u)\omega+h(u),\qquad h(u)\in H.               \tag{21.70}
```

The rooted condition makes `h(u)->h(phi(u))` well defined on the projected
span; Gram preservation makes it symplectic. The symplectic Witt extension
lemma extends it to `S in Sp(H)`. Then
`O(c omega+h)=c omega+Sh` is the required orthogonal extension.

For even `m`, choose

```math
V=\langle e,\omega\rangle\perp W,
\qquad e\cdot e=e\cdot\omega=1,                         \tag{21.70a}
```

with `W` nondegenerate alternating. Every orthogonal map has, and every
choice of the displayed parameters defines, the form

```math
\begin{aligned}
T(\omega)&=\omega,\\
T(w)&=Sw+(t\cdot Sw)\omega,\\
T(e)&=e+t+c\omega,
\end{aligned}                                           \tag{21.70b}
```

where `S in Sp(W)`, `t in W`, and `c in F_2`. If `omega` is absent from a
label span, first adjoin it on both sides; the rooted condition makes this a
well-defined isometric extension. The induced partial isometry on `W`
extends by the symplectic Witt lemma. Nondegeneracy of `W` then chooses `t`
to realize the remaining linear `omega`-coefficient, and `c` matches one
odd vector if present. Formula (21.70b) therefore extends the original span
isometry. The audited all-parity draft records the explicit substitution.

Conversely, every orthogonal map fixes the characteristic vector `omega`,
since `x dot x=omega dot x` characterizes it, and therefore preserves
(21.69).

Applying `O` simultaneously to both Walsh coordinates is a permutation of
the Boolean cube, preserves `W`, and conjugates every `C_(a_i)` to
`C_(a_i')`. This proves the extremal and equivariant-response claims.
`square`

This is a strict structured quotient, but not an independently composable
one. Gluing two tuples creates cross-Gram values and cross-relations not
determined by their isolated states. Those missing fibres are a concrete
form of composition-created information.

### Theorem 21.15 (rooted bilinear amalgamation)

For a label tuple `a=(a_1,...,a_k)`, let

```math
\alpha:F_2^k\to F_2^m,
\quad \alpha(c)=\sum_i c_i a_i,
\quad R_a=\ker\alpha,
\quad U_a=F_2^k/R_a.                                    \tag{21.71}
```

For two tuples `a,b`, their isolated states from (21.69) become exactly
composable after supplying three relative objects:

```math
\begin{aligned}
\kappa([c],[d])&=(\alpha c)\cdot(\beta d),\\
J_{ab}&=\ker\{U_a\oplus U_b\to F_2^m,
                  (u,v)\mapsto\bar\alpha u+\bar\beta v\},\\
Z_{ab}^{\times}&=\{(u,v):\bar\alpha u+\bar\beta v=\omega\}.
\end{aligned}                                           \tag{21.72}
```

Indeed the concatenated rooted state is reconstructed by

```math
\begin{aligned}
G_{a\sqcup b}((c,d),(c',d'))
 &=G_a(c,c')+G_b(d,d')\\
 &\quad+\kappa([c],[d'])+\kappa([c'],[d]),\\
R_{a\sqcup b}
 &=\{(c,d):([c],[d])\in J_{ab}\},\\
Z_{a\sqcup b}
 &=\{(c,d):([c],[d])\in Z_{ab}^{\times}\}.             \tag{21.73}
\end{aligned}
```

Conversely, given the two marked isolated states, their combined rooted
state determines all three objects in (21.72). Thus (21.72) is the minimal
relative datum for a lossless **orbit-complete** carrier. Composition is
associative on presented actual spans:

```math
(U_a\oplus U_b)/J_{ab}
\simeq\operatorname{span}(\mathbf a,\mathbf b),         \tag{21.74}
```

and every parenthesization pulls back the same form, kernel, and inverse
image of `omega` from the map `direct-sum_i U_i -> F_2^m`.

#### Proof

The concatenated coefficient map is
`gamma(c,d)=alpha(c)+beta(d)`. Expanding `gamma(c,d) dot gamma(c',d')`
proves the Gram formula; its zero and `omega` fibres are exactly the last two
formulas. Theorem 21.14 turns equality of the reconstructed state into a
global Walsh-coordinate conjugacy. Conversely, the off-diagonal Gram block
descends to `kappa`, while the images of the combined relation and root
fibres in `U_a direct-sum U_b` are `J_ab` and `Z_ab^times`. Associativity is
the first isomorphism theorem applied to the one accumulated coefficient
map. `square`

The theorem classifies actual realizations. It does not claim that every
abstract amalgamation datum embeds in a prescribed ambient dimension, nor
that every bit is distinguished by one fixed scalar graph maximum.

### Theorem 21.16 (quadratic cross-memory and a ternary obstruction)

The three resources in (21.72) are independent, and orbit-complete gluing
can require quadratically many bits.

1. If odd `m>=2(r+s)+1`, there are fixed isolated states, fixed `J=0`, and
   empty root fibre realizing all `2^(rs)` cross forms. Hence `rs` cross bits
   are necessary.
2. If odd `m>=2r+1`, there are fixed isolated states, zero cross form, and
   empty root fibre realizing `|GL(r,2)|` different intersection
   correspondences. Hence

```math
\log_2|GL(r,2)|=r^2+O(1)                                \tag{21.75}
```

   intersection bits are sometimes necessary.
3. The combined root fibre is not determined by the other two resources.
4. Complete pairwise amalgamation data does not determine a multi-piece
   composite: three singleton pieces can agree in every singleton and pair
   state while differing by a ternary relation.

At least one intersection bit is semantically visible at the full leading
scale: the two three-block Walsh path systems of Theorem 21.12 have identical
isolated piece states, cross form, and root fibre when the endpoints are one
piece and the middle is the other, but their maxima differ by at least

```math
{7-3\sqrt3\over2}n^{3/2}.                               \tag{21.76}
```

#### Proof

In `H=omega^perp`, choose a symplectic basis `(p_i,q_i)`. Fix
`a_i=p_i` and, for any binary `r by s` matrix `K`, take

```math
b_j^K=p_{r+j}+\sum_iK_{ij}q_i.                           \tag{21.77}
```

Private `p_(r+j)` coordinates make the second tuple independent and
disjoint from the first; both internal Gram forms vanish, while
`a_i dot b_j^K=K_ij`. This proves the first claim.

For the second, keep `a_i=p_i` and let
`b_j^P=sum_i P_ij p_i` for `P in GL(r,2)`. Every internal and cross Gram
form vanishes, but

```math
J_{ab}^P=\{(Pv,v):v\in F_2^r\}.                          \tag{21.78}
```

The standard product formula for `|GL(r,2)|` gives (21.75). For root
independence at odd `m>=5`, use
`a=e_1+e_2`, `b^+=omega+a`, and `b^-=e_3`: all unrooted relative data agree,
but only `a+b^+` equals `omega`.

Finally choose totally isotropic independent `p_1,p_2,p_3`. The triples

```math
(p_1,p_2,p_1+p_2),\qquad(p_1,p_2,p_3)                  \tag{21.79}
```

agree on every singleton and pair amalgam but only the first has the ternary
relation summing to zero. The path example is exactly the coincidence
collision already proved in Theorem 21.12. `square`

Thus the accumulated presented span, rather than a collection of pairwise
edge labels, is the dynamic memory that makes future coincidences meaningful.
Its `O(t^2)`-bit state is still strictly smaller than `mt` raw label bits
when a length-`t` word satisfies `t=o(m)`.

### Theorem 21.17 (unrooted Walsh spectra forget the root fibre)

Let `m>=1`, let `F` be the normalized order-`2^m` Walsh involution, and put

```math
J_a=D_aFD_a.                                             \tag{21.80}
```

If two ordered linear-label tuples have the same binary Gram matrix and
relation kernel, then for every real weighted graph on the marked tuple, the
coefficient carriers `K_(G,a)` from Theorem 21.9 have exactly the same
spectrum with multiplicity. Their characteristic-root fibres need not agree.

#### Proof

Write modulation and translation operators as `M_a,T_a`. Then

```math
J_a=M_aT_aF,
\qquad FM_aF=T_a,
\qquad T_bM_a=(-1)^(a\cdot b)M_aT_b.                    \tag{21.81}
```

Induction reduces every word to

```math
J_(a_1)\cdots J_(a_l)
=(-1)^theta M_sT_sF^(l\bmod2),
\qquad s=\sum_ja_j,                                    \tag{21.82}
```

where `theta` is a quadratic polynomial in the pairwise Gram values. Hence

```math
\operatorname{tr}(J_(a_1)\cdots J_(a_l))
=\begin{cases}
(-1)^theta2^m,&l\text{ even and }s=0,\\
0,&\text{otherwise}.
\end{cases}                                             \tag{21.83}
```

For odd `l`, the diagonal character sum is proportional to
`sum_x(-1)^(x dot x)=sum_x(-1)^(omega dot x)=0`.

Expand every power trace of `K_(G,a)` over closed block walks. Formula
(21.83) says that each summand depends only on whether one coefficient mask
lies in the relation kernel and on a phase determined by the Gram matrix.
All power traces, hence the finite symmetric-matrix spectra, agree. `square`

In particular, the two children in Theorem 21.13 have identical spectra in
every weighted graph made from synchronized copies of the compared child,
yet one canonical rooted Boolean future separates their projective responses
by `n^(3/2)/6`. An arbitrary common appended label is not covered: its cross
pairing can make the complete tuple Gram states differ.

Thus complete spectral information can be strictly smaller than rooted
extremal information even in a family with an exact orbit quotient.  The
next theorem resolves the corresponding semantic question: the entire
unrooted graph landscape also forgets this root, for a larger symmetry
reason not visible in the spectral proof alone.

### Theorem 21.18 (unrooted Walsh graphs forget the label-space root)

Let `V=F_2^m`, `E=V direct-sum V`, `n=2^(2m)`, and let `W_E` be the
order-`n` Walsh matrix.  Embed a linear label by

```math
iota(a)=(0,a),
\qquad C_a=D_(iota(a))W_ED_(iota(a)).                  \tag{21.84}
```

Suppose two ordered source tuples `a,b` have the same binary Gram matrix and
relation kernel.  Then one coordinate permutation `P` satisfies

```math
PW_EP^T=W_E,
\qquad PC_(a_i)P^T=C_(b_i)\quad\hbox{for every }i.     \tag{21.85}
```

More generally, a continuation may use any synchronously derived labels
`a[c]=sum_i c_i a_i`, arbitrary real scalar onsite weights, and arbitrary
real symmetric edge weights.  Applying the same `P` in every block
conjugates the **entire** two Boolean landscapes pointwise.  Hence their
upper and absolute maxima, minima, histograms, and optimizer multiplicities
all agree.  The pair

```math
(G(a),R(a))                                             \tag{21.86}
```

is therefore an exact unrooted weighted-graph response state; no
label-space characteristic-root fibre is required.

#### Proof

Equality of relation kernels and Gram forms makes

```math
(0,a_i)\longmapsto(0,b_i)                              \tag{21.87}
```

a well-defined isometry between the presented subspaces of `E`.  The
characteristic vector of `E` is

```math
Omega_E=(omega,omega).                                 \tag{21.88}
```

It lies in neither subspace, since their first coordinate is zero.  Adjoin
it on both sides and fix it.  This is still an isometry because

```math
Omega_E dot(0,a)=a dot a=(0,a) dot(0,a).               \tag{21.89}
```

The characteristic-rooted Witt extension lemma extends the partial map to
`O in O(2m,2)`.  The permutation `(Pf)(z)=f(O^(-1)z)` preserves `W_E` and
carries every modulation `D_(iota(a[c]))` to
`D_(iota(b[c]))`.  This proves simultaneous conjugacy of every displayed
child and bridge matrix, and hence the pointwise landscape identity.
`square`

This theorem strictly sharpens Theorem 21.17.  For odd `m>=3`, the constant
tuples `omega^k` and `e_1^k` have equal `(G,R)`; for even `m>=4`, use
`omega^k` and `(e_1+e_2)^k`.  Their label-space root fibres differ, yet all
unrooted weighted graph landscapes agree.  A fixed external pole or
coordinate-dependent field is outside the theorem: the canonical rooted
future of Theorem 21.13 breaks the enlarged ambient symmetry and separates
the singletons by `n^(3/2)/6`.

### Theorem 21.19 (query-local Walsh amalgamation and a sharp local rate)

Present two marked pieces by maps `alpha:F_2^r->V` and
`beta:F_2^s->V`, and let `U_a,U_b` be their quotient image spaces.  Their
relative rooted data are the cross form, coincidence relation, and root
fibre

```math
\begin{aligned}
kappa(u,v)&=u dot v,\\
J&=\{(u,v):u=v\},\\
Z^times&=\{(u,v):u+v=omega\}.
\end{aligned}                                           \tag{21.90}
```

Fix a finite family `Theta` of graph queries.  A connected support is the
active-vertex set of one component of the graph of nonzero bridge weights;
rooted terms must be componentwise.  Let `C_Theta` be the inclusion-maximal
supports, and let `P_C subset U_a,Q_C subset U_b` be the spans generated by
the marked labels in `C`.  Then the isolated local rooted states together
with

```math
\left(
kappa|_(P_C times Q_C),
J cap(P_C direct-sum Q_C),
Z^times cap(P_C direct-sum Q_C)
\right)_(C in C_Theta)                                 \tag{21.91}
```

answer every query in `Theta` exactly.  Among states reconstructing each
local rooted orbit triple, (21.91) is coarsest.  On a purely unrooted
component, Theorem 21.18 deletes the last coordinate exactly: `(G,R)` and
the relative pair `(kappa,J)` suffice.

If `w=max_C|C|` and `L=sum_C|C|`, a direct relative presentation uses

```math
O\left(sum_C|C|^2\right)=O(wL)                         \tag{21.92}
```

bits.  Thus bounded `w` and `L=O(r+s)` give a linear-size exact carrier,
instead of the unrestricted quadratic orbit state.

#### Proof

On one support the concatenated coefficient map is
`gamma_C(c,d)=alpha(c)+beta(d)`.  Its pulled-back Gram form uses only the
restricted cross form, while its zero and root fibres are precisely the
preimages of the two restricted sets in (21.91).  Conversely the local
concatenated rooted state recovers all three restrictions.  The rooted orbit
theorem supplies an independent coordinate conjugacy on every connected
component; disconnected Boolean maxima split, so no global conjugacy is
needed.  The unrooted deletion of `Z^times` is Theorem 21.18.  A cross form,
a row-reduced basis for `J`, and an optional affine root representative use
`O(|C|^2)` bits per maximal support. `square`

The linear order in (21.92) can be necessary even for ordinary unrooted
scalar maxima.  For every `h`, use disjoint three-coordinate chunks with

```math
a_i=(1,1,1),\qquad c_i=(0,0,1),                        \tag{21.93}
```

leave at least one ambient coordinate unused, put two copies of every `a_i`
in the first piece, and choose the second-piece label independently as
`b_i=a_i` or `c_i`.  All isolated states, the complete cross form, and every
root fibre are fixed; only `J` records the `h` choices.  The `i`th unrooted
three-block path query compares `(a_i,a_i,a_i)` with `(a_i,c_i,a_i)`.
Theorem 21.12 gives the gap

```math
Delta n^(3/2),\qquad Delta={7-3sqrt3\over2}.            \tag{21.94}
```

Hence error below `Delta n^(3/2)/2` requires `2^h` response states, or `h`
bits.  Here `n` is one Walsh-block order.  This is a semantic lower bound on
the surviving unrooted coincidence memory, not merely an orbit count.

### Theorem 21.20 (interaction mass controls approximate compatibility memory)

For a real weighted linear-label Walsh graph, write

```math
E_(G,a)(x)=sum_vh_vH_(a_v)(x_v)
 +sum_(uv in E(G))w_(uv)x_u^TWx_v,                    \tag{21.95}
```

where `H_a(x)=x^TC_ax/2` and one Walsh block has order `n`.  For a public
partition `P` of the active vertices put

```math
d_G(P)=sum_(uv:[u]_P ne[v]_P)|w_(uv)|.                 \tag{21.96}
```

Deleting the cross-part edges and decoding each retained component from its
exact unrooted state `(G(a|_C),R(a|_C))` approximates the upper optimum with
error at most

```math
d_G(P)n^(3/2).                                         \tag{21.97}
```

The same error holds for the absolute optimum when the component decoder is

```math
max\left\{sum_Cmax E_C,-sum_Cmin E_C\right\}.          \tag{21.98}
```

For a finite declared query family, use one public partition per query and
store only inclusion-maximal parts.  A simultaneous carrier uses

```math
O\left(sum_(C\ maximal)|C|^2\right)                    \tag{21.99}
```

bits.  In particular a unit path on `t` blocks has an
`O(t/eta)`-bit carrier at error `eta t n^(3/2)`.  Conversely, for a complete
graph whose edge magnitudes are at least `c`, every partition with deleted
mass at most `eta t^2` satisfies

```math
sum_C|C|^2>=\left(1-{2eta\over c}\right)t^2.           \tag{21.100}
```

The last inequality is a ceiling on this deletion architecture, not a
semantic lower bound for every possible code.

#### Proof

The Walsh bridge has operator norm `sqrt n`, so
`|x^TWy|<=n^(3/2)`.  The pointwise difference between the original and
truncated energies is therefore at most (21.97); maxima, minima, and absolute
maxima are one-Lipschitz under a uniform perturbation.  Theorem 21.18 gives
the exact component state, and disconnected upper optima add.  A maximal
state restricts to every contained part; Gram forms and row-reduced relation
kernels use quadratic bits, proving (21.99).  Cut a path every
`ceil(1/eta)+1` vertices.  On `K_t`, the number of cross-part edges is exactly
`(t^2-sum_C|C|^2)/2`, which proves (21.100). `square`

### Theorem 21.21 (an off-diagonal Gram flux is scalar-visible)

For every `m>=5`, let

```math
\begin{aligned}
a^0&=(e_1+e_2,e_3+e_4,e_1+e_2+e_3+e_4),\\
a^1&=(e_1+e_2,e_1+e_3,e_2+e_3).
\end{aligned}                                          \tag{21.101}
```

Both tuples have zero self-pairings, relation kernel `{000,111}`, and empty
characteristic-root fibre.  Their off-diagonal Gram values are respectively
all zero and all one.  For the ordinary unrooted, unweighted triangle query,

```math
max E_(a^0)={9\over2}n^(3/2),                           \tag{21.102}
```

whereas

```math
max E_(a^1)
 <={3(1+sqrt17)\over4}n^(3/2).                         \tag{21.103}
```

Thus one Gram/triangle-flux bit is exposed with gap

```math
{3(5-sqrt17)\over4}n^(3/2).                            \tag{21.104}
```

#### Proof

For even labels `a,b` and `c=a+b`, the reduced Walsh involutions satisfy

```math
J_aJ_bJ_c=(-1)^(a dot b)F,                             \tag{21.105}
```

and all four involutions commute.  In the zero-flux tuple, the self-dual
chirp `(-1)^(u dot v)` and a characteristic-rooted Witt transport give one
Boolean vector fixed by `F` and all three children.  Putting it in every
block saturates the three child and three bridge terms, proving (21.102).

On a joint eigenspace in the unit-flux case, multiply the three child signs
by the Walsh sign.  Their product is `-1`, and the normalized three-block
operator becomes, up to a scalar sign,

```math
diag(s_1,s_2,s_3)+A(K_3),\qquad s_1s_2s_3=-1.          \tag{21.106}
```

Its norm is at most `(1+sqrt17)/2`: the all-negative case has norm two, and
the one-negative case has eigenvalues `0,(1+-sqrt17)/2`.  A Boolean
three-block vector has squared norm `3n`; the prefactor `q/2`, with
`qn=n^(3/2)`, gives (21.103). `square`

Theorem 21.18 proves `(Gram,relations)` sufficient.  The present theorem is
a genuine partial converse: after relations and all rooted information are
fixed, scalar unrooted maxima cannot discard every off-diagonal Gram bit.
It does not prove entrywise recovery of a general Gram matrix; the exposed
quantity is naturally a relation-cycle flux.

### Theorem 21.22 (connected packing of independent Walsh fluxes)

For every `h>=1` there are `2^h` marked tuples of `3h` linear Walsh
children, indexed by `sigma in {0,1}^h`, with every self-pairing, the
complete relation kernel, the characteristic-root fibre, and every
cross-gadget Gram entry fixed.  Only one binary off-diagonal Gram flux in
each of `h` disjoint relation triangles varies.  There are `h` public scalar
unrooted queries, all using the same connected support on all `3h` children,
whose response vectors obey

```math
||R(sigma)-R(tau)||_infinity
\ge \Delta_* n^(3/2)\quad(sigma\ne tau),               \tag{21.107}
```

where one Walsh child has order `n=2^(2(4h+1))` and

```math
\Delta_*={3(5-\sqrt {17})\over4}-{9\over100}
=0.5676707807\ldots.                                   \tag{21.108}
```

The common support may be chosen either with maximum degree four or as the
complete graph.  Consequently every summary answering these queries to
uniform error below `Delta_* n^(3/2)/2` has at least `2^h` states.

#### Proof

In the `i`th disjoint four-coordinate chunk take

```math
u_i=e_(4i-3)+e_(4i-2),
```

and choose either

```math
v_i^0=e_(4i-1)+e_(4i)
\quad\hbox{or}\quad
v_i^1=e_(4i-3)+e_(4i-1),                               \tag{21.109}
```

with third label `u_i+v_i^(sigma_i)`.  All labels are even.  Each local
relation kernel is generated by `111`; different chunks are orthogonal and
linearly disjoint; and one unused ambient coordinate makes the root fibre
empty.  The three local off-diagonal pairings are all zero for `sigma_i=0`
and all one for `sigma_i=1`.

The triangle query in Theorem 21.21 has favorable value
`M_0 n^(3/2)` and unfavorable upper bound `M_1 n^(3/2)`, where

```math
M_0={9\over2},\qquad M_1={3(1+\sqrt {17})\over4}.       \tag{21.110}
```

Connect all marked blocks by any nonnegative public Walsh-bridge graph of
total weight `B`.  Every connector contributes at most `B n^(3/2)`.  In the
zero-flux state the common Boolean witness from Theorem 21.21 may be put in
every block, so it attains every nonnegative connector ceiling
simultaneously.  Thus connective padding adds the same public ceiling to the
favorable value and cannot erode the local gap.

To make every nontarget child active, add all other local triangles with
weight `gamma=1/(100(h-1))` when `h>=2`.  Their total absolute contribution
is at most `9n^(3/2)/200`; comparing two states pays this twice.  The
remaining gap is (21.108).  A Hamiltonian path through the blocks, together
with the local triangles, has maximum degree four; all cross-gadget edges
give the complete-support version.  A differing flux coordinate selects a
public query that realizes (21.107), and the packing conclusion follows.
`square`

The normalization matters.  The full query has `N=3hn` Boolean variables,
so its separation in units of `N^(3/2)` is
`Delta_*/(3h)^(3/2)`.  The theorem proves additive scalar visibility of `h`
independent cycle fluxes at fixed one-port accuracy, not an extensive
free-energy-density packing for the whole connected graph.  The complete
construction and independent audit are in
[`drafts/walsh_connected_flux_packing.md`](drafts/walsh_connected_flux_packing.md)
and
[`drafts/walsh_connected_flux_packing_independent_audit.md`](drafts/walsh_connected_flux_packing_independent_audit.md).

### Theorem 21.23 (state-local information vanishes at total scale)

Partition `N=kn` Boolean variables into `k` blocks.  Let a hidden state enter
a scalar query only through bounded onsite Walsh children:

```math
E_a^theta(X)=S_theta(X)+sum_(v=1)^k d_(theta v)H_(a_v)(x_v),
\qquad |d_(theta v)|<=D,                               \tag{21.111}
```

where `S_theta` is arbitrary and state-independent and
`|H_a(x)|<=n^(3/2)/2`.  Then every public query family has response metric

```math
d_Theta(a,b)
<=D\#\{v:a_v!=b_v\}n^(3/2)
<= {D\over sqrt k}N^(3/2).                            \tag{21.112}
```

The same holds for minima and absolute maxima.  Thus no two such states have
a fixed positive total-scale separation as `k->infinity`, regardless of how
dense, signed, or large the public bridge term is.

More generally, suppose hidden bit `i` affects only a disjoint cell of `s_i`
variables and

```math
||E_sigma^theta-E_tau^theta||_infinity
<=Lsum_(i:sigma_i!=tau_i)s_i^(3/2).                   \tag{21.113}
```

For every fixed `epsilon>0`, a response code whose distinct pairs have
distance greater than `epsilon N^(3/2)` satisfies

```math
log_2|C|<=ceil(4L^2/epsilon^2).                        \tag{21.114}
```

#### Proof

The public term cancels pointwise.  One changed onsite child changes the
landscape by at most `D n^(3/2)`, and `max`, `min`, and `max|.|` are all
one-Lipschitz in uniform norm.  This proves (21.112).

For (21.114), order `s_1>=s_2>=...`.  Since `i s_i<=N`,

```math
sum_(i>r)s_i^(3/2)
<=N^(3/2)sum_(i>r)i^(-3/2)
<=2N^(3/2)/sqrt r.                                    \tag{21.115}
```

More than `2^r` codewords contain two agreeing on the first `r` bits.  With
`r=ceil(4L^2/epsilon^2)`, their distance is at most the forbidden threshold.
`square`

For the `h`-flux cube in Theorem 21.22, `k=3h` and changing one flux changes
only two onsite children.  Therefore

```math
d_Theta(sigma,tau)<=2D d_H(sigma,tau)n^(3/2),         \tag{21.116}
```

the whole cube has total-scale diameter `O(h^(-1/2))`, and its minimum
pairwise separation is at most `2D/(3h)^(3/2)` in units of `N^(3/2)`.
Theorem 21.22 attains the correct full-cube order.  Error-correcting subcodes
cannot make the diameter constant.

Within the coordinate-equivariant unrooted Walsh-graph query language of
Theorem 21.18, `r` source labels have at most
`2^((3r^2+r)/2)` possible `(Gram,relation)` states.  Encoding `h` scalar-
distinct states requires

```math
r>={sqrt(1+24h)-1\over6}.                              \tag{21.117}
```

If those labels occupy `t>=r` child slots, (21.112) gives normalized diameter
at most `D/sqrt(t)=O(h^(-1/4))`.  This orbit-count corollary excludes
coordinate pins, which break the ambient symmetry; the diameter theorem
(21.112) itself allows them.

Finally, if flipping one hidden bit can alter at most `d_i` unit-weight atoms,
each bounded by `B n^(3/2)`, a total-scale neighboring gap
`epsilon(kn)^(3/2)` requires

```math
d_i>={epsilon\over2B}k^(3/2).                          \tag{21.118}
```

Therefore a successful total-scale escape must broadcast hidden
compatibility into genuinely state-dependent cross-block coefficients,
use unbounded normalized weights, or change the output scale/language.
Public connector padding alone cannot amplify local state.  The repaired
proof and audit are in
[`drafts/walsh_total_scale_flux_ceiling.md`](drafts/walsh_total_scale_flux_ceiling.md).

### Theorem 21.24 (alternating-form broadcast attains positive total-scale rate)

For every `r>=2`, put `V=F_2^r`,

```math
h={r(r-1)\over2},\qquad k=64r^2.                     \tag{21.119}
```

There are one public label list `P=(p_1,...,p_k) in V^k` and one public
hollow signing `A` such that, for all alternating forms `B` on `V`,

```math
A_B(i,j)=A(i,j)(-1)^(B(p_i,p_j))                    \tag{21.120}
```

defines an exact hollow signing with

```math
||A_B||_(2->2)<=8sqrt(k),\qquad Q(A_B)<=4k^(3/2).    \tag{21.121}
```

Declare the same-support additive contexts `-H_T`, one for every alternating
form `T`, where `H_B(s)=s^TA_Bs/2`.  Their absolute responses satisfy

```math
R_B(B)=0,
\qquad {sqrt2\over32}k^(3/2)
       <=R_T(B)=Q(A_B-A_T)<=8k^(3/2)\quad(B!=T).     \tag{21.122}
```

Consequently these `2^h` bounded-cap sign quadratics have pairwise response
separation `(sqrt2/32)k^(3/2)`.  Uniform error below half that gap requires
at least

```math
h>={k\over256}                                      \tag{21.123}
```

bits.  This is a positive contextual information rate at the total
extremal scale.

#### Proof

Choose the `p_i` independently.  Every nonzero alternating form has rank at
least two, so `Pr(B(p_i,p_j)=1)>=3/8`.  McDiarmid gives probability at most
`exp(-k/128)` that fewer than `binom(k,2)/4` pairs are detected.  Since
`2^h exp(-k/128)<1`, one list simultaneously detects every nonzero form on
at least that many pairs.

Now choose the upper entries of `A` independently.  For fixed `B` and a
unit vector, the quadratic form is a Rademacher sum with squared coefficient
norm at most two.  A `1/4` sphere net and a union bound over all `2^h` forms
give

```math
2exp\{klog9+hlog2-4k\}<1,                            \tag{21.124}
```

so one `A` satisfies (21.121) simultaneously.

For `B!=T`, the difference has entries in `{0,+-2}` and at least `k^2/16`
nonzero unordered entries.  Some vertex bipartition cuts at least `k^2/32`
of them.  Optimizing one side after randomizing the other and using the sharp
`p=1` Khintchine inequality gives a cross value at least
`(sqrt2/32)k^(3/2)`.  Flipping all spins on one side shows the absolute full
quadratic value is at least this large.  The upper bound is
`k||A_B-A_T||/2<=8k^(3/2)`.  The predeclared coordinate `T=B` separates any
pair of response vectors, and `h>=r^2/4=k/256`. `square`

This theorem realizes the escape left open by Theorem 21.23: one hidden
compatibility coordinate is broadcast through a dense set of state-dependent
interaction atoms.  The edge phases form a linear alternating-form
evaluation code and have only `h=Theta(k)` hidden bits, not independent
quadratic edge data.

The scope is essential.  The random shared base `A` is `Theta(k^2)` bits of
nonuniform public advice.  Child and context are exact signings separately,
but their same-edge overlay has coefficients in `{0,+-2}` and is not an
exact-sign parent or a disjoint appended future.  Thus this is a structured
contextual incompressibility theorem, not closure of the original signing
class.  Fixed-order regular-Hadamard block lifts preserve the normalized gap
and positive rate; growing the inner order makes the hidden rate
`Theta(1/n)`.  The complete audited proof is in
[`drafts/state_dependent_gram_broadcast.md`](drafts/state_dependent_gram_broadcast.md).

### Theorem 21.25 (bounded fan-in has a sharp broadcast law)

Let `z in {0,1}^h` and suppose a fixed child presentation has bounded atoms

```math
H_z(x)=P(x)+\sum_(e=1)^E c_e(z)\phi_e(x),
\qquad |\phi_e(x)|\le1,                              \tag{21.125}
```

where `c_e` depends on at most `t` coordinates of `z` and has oscillation
`omega_e`.  For any fixed public continuation language and any response
functional that is one-Lipschitz in uniform norm, let `d` be the resulting
contextual pseudometric.  Then, at every `z`,

```math
\sum_(i=1)^h d(z,z+e_i)
\le\sum_(e=1)^E |I_e|\omega_e
\le t\sum_(e=1)^E\omega_e.                          \tag{21.126}
```

Thus, for distinct quadratic pair atoms on `N` spins with `|c_e|<=B`, if
every neighbouring hidden state is separated by `epsilon N^(3/2)`, then

```math
h<{Bt\over\epsilon}\sqrt N.                         \tag{21.127}
```

This scale is optimal for unrestricted exact sign quadratics.  There is an
absolute `c>0` such that, for all sufficiently large `N` and
`1<=t<=cN`, one can construct `h=t floor(sqrt N)`-bit exact signings and a
fixed negative-clone continuation language such that every edge reads at
most `t` hidden bits, every neighbour has response gap
`Omega(N^(3/2))`, and an `exp(Omega(t sqrt N))` subfamily is pairwise
separated by `Omega(N^(3/2))`.

Finally, in the alternating-form broadcast of Theorem 21.24, every choice
of linear coordinates has

```math
{1\over E}\sum_e |I_e|\ge {h\over4},
\qquad \max_e|I_e|\ge {h\over4},                    \tag{21.128}
```

and each hidden coordinate changes at least `E/4=Theta(N^2)` edge phases.
The flat construction is therefore genuinely high-fan-in and high-influence,
not a disguised collection of independently paid local channels.

#### Proof

For a fixed future, uniform-norm nonexpansiveness cancels that future and
gives

```math
d(z,z+e_i)
\le\sum_(e:i\in I_e)|c_e(z)-c_e(z+e_i)|.            \tag{21.129}
```

Summing over `i` proves (21.126), and `E<=N(N-1)/2` gives (21.127).

For sharpness, partition almost all edges into `g=floor(sqrt N)` equal
cells.  In each cell, use a lookup among `2^t` public switching vectors.
A probabilistic library makes every pair of switches disagree on a fixed
fraction of every cell.  Changing one hidden bit then flips
`Theta(N^(3/2))` coherently signed coefficients.  A positive-rate outer
code over the `2^t` symbols makes two selected words differ on a fixed
fraction of all `Theta(N^2)` edges.  If a symmetric hollow difference has
`m` unordered entries of magnitude two, a random vertex bipartition and the
sharp `p=1` Khintchine inequality give

```math
\max_x\left|\sum_(u<v)D_(uv)x_ux_v\right|
\ge {m\over\sqrt(2N)},                              \tag{21.130}
```

which supplies the pairwise total-scale gap.

For (21.128), write edge evaluation in any basis of alternating forms.  A
nonzero basis vector is detected on at least `E/4` edges by the sampler in
Theorem 21.24.  Double-counting coordinate--edge incidences proves both
bounds. `square`

Neighbour separation alone forces only a two-colouring of the hidden cube;
the outer code is what turns incidence into a genuine information packing.
The sharpness family need not be spectrally flat, uses an exponential public
spin library, and uses same-support additive overlays.  The theorem therefore
classifies the bounded-atom incidence resource, not low public-description
or exact-sign disjoint composition.  Full details and audit are in
[`drafts/bounded_fanin_broadcast_law.md`](drafts/bounded_fanin_broadcast_law.md)
and
[`drafts/bounded_fanin_broadcast_law_independent_audit.md`](drafts/bounded_fanin_broadcast_law_independent_audit.md).

### Theorem 21.26 (short-seed flat Gram broadcast)

For every `r>=2`, put

```math
h={r(r-1)\over2},\qquad k=256r^2.                    \tag{21.131}
```

There is a uniformly computable public description of `O(k log k)` bits
defining `2^h` exact hollow sign children `A_B`, indexed by alternating
forms `B` on `F_2^r`, such that

```math
||A_B||_(2->2)<=8sqrt(k),\qquad Q(A_B)<=4k^(3/2).    \tag{21.132}
```

For the predeclared same-support contexts `-H_T`, with
`H_B(x)=x^TA_Bx/2`,

```math
R_B(B)=0,
\qquad {sqrt2\over16}k^(3/2)
 <=R_T(B)<=8k^(3/2)\quad(B!=T).                      \tag{21.133}
```

Consequently uniform response error below
`(sqrt2/32)k^(3/2)` requires at least

```math
h>={k\over1024}                                      \tag{21.134}
```

bits, even though the shared base no longer stores a quadratic edge table.

#### Proof

First choose a `1/8`-biased multiset `S subset F_2^r` of size `256r` and
repeat every indexed occurrence `r` times to form the public list
`P=(p_1,...,p_k)`.  Hoeffding and a union bound over the nonzero linear
characters prove existence; storing `S` costs exactly `k` bits.  For every
nonzero alternating form `B`, Fourier expansion of its radical and the bias
bound give

```math
#\{i<j:B(p_i,p_j)=1\}>={k^2\over8}.                 \tag{21.135}
```

Thus the alternating-form evaluation words form an injective distance code.

Let `E=binom(k,2)`, `d=ceil(log_2E)`, and index the edges by distinct points
of `F_(2^d)`.  The traces of a uniformly random polynomial of degree below
`6k` give `6k`-wise-independent edge signs from a seed of

```math
6k\lceil\log_2E\rceil=O(k\log k)                   \tag{21.136}
```

bits.  For a fixed form and a fixed unit vector, the `6k`-th moment agrees
with the fully independent Rademacher moment.  Taking `m=3k` gives

```math
Pr\{|z^T(A\odot\chi_B)z|>4sqrt(k)\}
 <=(3/8)^(3k).                                      \tag{21.137}
```

A `1/4` sphere net and a union bound over all `2^h` forms have failure
probability

```math
9^k2^h(3/8)^(3k)<e^(-0.74k).                        \tag{21.138}
```

One seed therefore gives (21.132) simultaneously.  If `B!=T`, (21.135)
makes `A_B-A_T` nonzero with magnitude two on at least `k^2/8` unordered
edges.  The bipartition--Khintchine bound (21.130) yields the lower half of
(21.133), and the operator bounds yield its upper half.  The coordinate
`T=B` in the already-declared response vector separates every pair; (21.134)
follows from `h/k=(r-1)/(512r)`. `square`

The sample space is uniformly computable by exhaustive finite search, and a
random short seed succeeds with exponentially high probability.  No
deterministic polynomial-time certification is proved.  The public-bit count
excludes the `h` hidden bits and the exponential same-support query language.
Most importantly, the overlay remains weighted in `{0,+-2}` rather than an
appended disjoint exact-sign composition.  The remaining cost is therefore
closure/realization, not quadratic public randomness.  Full details and
audit are in
[`drafts/short_seed_gram_broadcast.md`](drafts/short_seed_gram_broadcast.md)
and
[`drafts/short_seed_gram_broadcast_independent_audit.md`](drafts/short_seed_gram_broadcast_independent_audit.md).

### Theorem 21.27 (exact sparse compilation and independent-star barriers)

Let `x in {+-1}^k`, `E=binom(k,2)`, and

```math
H_T(x)=\sum_(i<j)T_(ij)x_ix_j,
\qquad C_T(x)=E-H_T(x).                              \tag{21.139}
```

There is an exact disjoint unit-edge realization of `C_T`: introducing one
spin `y_(ij)` per old edge and setting

```math
G_T(x,y)=\sum_(i<j)y_(ij)(x_i-T_(ij)x_j)             \tag{21.140}
```

gives `max_y G_T(x,y)=C_T(x)`.  This compiler has `E` auxiliaries.

More generally, suppose `C_T` is exactly realized by independent sign
stars,

```math
C_T(x)=c+\sum_(a=1)^m
 \left|\sum_(i\in S_a)\sigma_(ai)x_i\right|,
\qquad d_a=|S_a|,\quad \Delta=\max_a d_a.            \tag{21.141}
```

Then

```math
\sum_a d_a^(3/2)\ge k(k-1),
\qquad
m\ge {k(k-1)\over\Delta^(3/2)},
\qquad
\sum_a d_a\ge {k(k-1)\over\sqrt\Delta}.            \tag{21.142}
```

Endpoint-local stars attain the sharp count `m=E`.  If every star instead
has full support `[k]`, then

```math
m\gamma_k\ge\max_xH_T(x),
\qquad
\gamma_d=2^(-(d-2)){d-2\choose\lfloor(d-1)/2\rfloor}.
                                                               \tag{21.143}
```

Consequently, if separate full-support compilers with at most `m` stars are
available for both `T` and `-T`, then

```math
m\ge {Q(T)\over\gamma_k}=\Omega(k^2).               \tag{21.144}
```

Thus neither bounded-fan-in nor fully dense independent selectors give a
linear-order, scale-preserving exact compiler.

#### Proof

For one edge, pointwise optimization gives

```math
\max_(y=+-1)y(x_i-T_(ij)x_j)
=|x_i-T_(ij)x_j|=1-T_(ij)x_ix_j,                    \tag{21.145}
```

proving the sparse identity.  The pair Fourier coefficient of a star of
size `d` is

```math
\widehat{\left|\sum_(i\in S)\sigma_ix_i\right|}({i,j})
=\gamma_d\sigma_i\sigma_j.                         \tag{21.146}
```

Conditioning on the other `d-2` spins proves this formula, and the central
binomial estimate gives

```math
0<\gamma_d\le {1\over\sqrt(d-1)},
\qquad
\gamma_d{d\choose2}\le {1\over2}d^(3/2).           \tag{21.147}
```

The target pair coefficients all have modulus one.  Summing their absolute
values and applying the triangle inequality proves (21.142).

For full-support stars, antipodally symmetrize their empirical sign law
`mu`.  The function `|sum_i x_i|` has a nonzero Fourier coefficient at every
even level.  Since `C_T` has no levels four and above, all corresponding
moments of `mu` vanish, while its pair moments are
`-T_(ij)/(m gamma_k)`.  Fourier inversion is therefore exact:

```math
\mu(\sigma)=2^(-k)\left(1-{H_T(\sigma)\over m\gamma_k}\right). \tag{21.148}
```

Nonnegativity proves (21.143).  Apply it separately to `T` and `-T`, use
`gamma_k=Theta(k^(-1/2))`, and use the universal
`Q(T)=Omega(k^(3/2))` signing bound to obtain (21.144). `square`

The two-orientation statement uses two separate compilers; it does not turn
absolute response into one future.  Independent-star factorization is the
essential scope.  Complete details, including the middle-slice Krawtchouk
calculation proving nonvanishing of every even Fourier level, are in
[`drafts/exact_disjoint_star_compiler_barrier.md`](drafts/exact_disjoint_star_compiler_barrier.md)
and its
[`independent audit`](drafts/exact_disjoint_star_compiler_independent_audit.md).

### Theorem 21.28 (interacting selector and bounded-cap barriers)

Let `k` be even and let

```math
C_+(x)={k^2-(\sum_i x_i)^2\over2}.                  \tag{21.149}
```

Consider an arbitrary correlated selector envelope

```math
F(x)=\max_(y\in{+-1}^m)\{c(y)+b(y)\mathbin\cdot x\},
\qquad |b_i(y)|\le d_i,                              \tag{21.150}
```

and suppose `||F-C_+||_infinity<=eta<k^2/4`.  Put

```math
a={k^2\over2}-2\eta,
\qquad D^2=\sum_i d_i^2.                             \tag{21.151}
```

Then

```math
\sum_i d_i\ge a,
\qquad
m\log2\ge {a^2\over2D^2}-\log(k+1).                \tag{21.152}
```

For a complete unit-sign old--new block, `d_i=m`; hence every
`o(k^2)`-accurate compiler satisfies

```math
m\ge(8\log2)^(-1/3)k-o(k)=0.565\ldots k-o(k).       \tag{21.153}
```

There is also a query-relative form.  Suppose

```math
F(x)=\max_(q\in[K])\{c_q+b_q\mathbin\cdot x\},
\qquad ||b_q||_2\le D,
\qquad ||F-f||_infinity\le\eta.                     \tag{21.154}
```

If some antipodal pair `p,-p` and set `X` obey

```math
f(x)\ge\max\{f(p),f(-p)\}+a+2\eta
\quad(x\in X),                                     \tag{21.155}
```

then

```math
\log K\ge\log|X|-k\log2+{a^2\over2D^2}.           \tag{21.156}
```

Finally, if `b(y)=By` and the selector is induced by a complete quadratic
parent `P` on `N` spins, then

```math
\mathop{osc}_xF\le2||B||_(infinity->1)\le2Q(P).     \tag{21.157}
```

Thus no parent with `N=O(k)` and `Q(P)=O(N^(3/2))` uniformly approximates
the all-positive cut shell to `o(k^2)`.

#### Proof

At every balanced `x`, choose an affine piece active in (21.150).  Comparing
it with the envelope at the two constant spins and adding the inequalities
gives `x dot b(y)>=a`.  Each selector covers at most

```math
2^k\exp\{-a^2/(2D^2)\}                              \tag{21.158}
```

cube points by Hoeffding.  The `2^m` selectors cover the balanced slice,
whose size is at least `2^k/(k+1)`, proving (21.152).  The same two-antipode
comparison on `X` proves (21.156).

For (21.157), maxima of affine functions give the first Lipschitz bound.
Flipping the whole auxiliary shore reverses `x^TBy` and preserves both
internal quadratic energies, so one of the two parent energies has absolute
value at least `|x^TBy|`; maximize over the shores. `square`

Unlike Theorem 21.27, the selector may have arbitrary auxiliary--auxiliary
interaction.  Its universal target is nevertheless the high-oscillation cut
shell.  Theorem 21.28 does not rule out contextual cancellation on a common
background or a linear-order compiler tailored to the flat Gram family.  It
identifies the extra input needed there: exponentially many antipodally
exposed configurations at a `Theta(k^(3/2))` gap, or a different lower-bound
mechanism.  The complete audited proof and finite verifier are in the same
draft and audit cited above.

### Theorem 21.29 (exact disjoint coordinate compiler)

For a hollow complete signing `A` of order `k`, write

```math
H_A(x)={1\over2}x^TAx.
```

Given `u in {+-1}^k`, append `k` spins `y`, use the complete sign bridge
`R_u=u1^T`, and put the positive clique `J-I` on the new block.  The
one-sided complete-parent response satisfies the exact identity

```math
F_u(A)=\max_(x,y)\{H_A(x)+x^TR_uy+H_(J-I)(y)\}
={3k^2-k\over2}+H_A(u).                             \tag{21.159}
```

Consequently

```math
\sup_u|F_u(A)-F_u(A')|
=\max_u|H_A(u)-H_(A')(u)|
=Q(A-A').                                           \tag{21.160}
```

Thus the full same-support quadratic response metric embeds isometrically
into disjoint complete exact-sign contexts on `N=2k` vertices.  Applied to
the short-seed family of Theorem 21.26, it gives `2^h` exact sign children
with

```math
h\ge {N\over2048},
\qquad
d_(response)(B,T)\ge {N^(3/2)\over32}.              \tag{21.161}
```

#### Proof

Put `a=x dot u` and `b=y dot 1`.  The appended energy is

```math
ab+{b^2-k\over2}.                                   \tag{21.162}
```

Global sign changes let an optimizer have `a,b>=0`.  Convexity in `b`
then puts `b=k`.  If `x` has projective Hamming distance `d<=k/2` from `u`,
the lock loses `2kd` from its planted value.  Switching those `d` old spins
changes only `d(k-d)` child edges, so the child can gain at most
`2d(k-d)`.  The net change is at most `-2d^2`.  This proves (21.159), and
subtraction proves (21.160).  The constants in (21.161) follow from
Theorem 21.26 and `N=2k`. `square`

The response `F_u` is one-sided, and its bridge depends on the declared
coordinate `u`.  Most importantly, (21.159) contains a common
`Theta(N^2)` calibration.  It cancels in response differences but is not a
low-cap signing construction.  Exact-sign/disjoint **metric closure** is
therefore solved; bounded-cap closure remains open.

### Theorem 21.30 (universal equality-lock ceilings)

Let a complete signing `L` act on paired blocks
`X=(x_1,...,x_k)` and `Y=(y_1,...,y_k)`.  Every duplicate `(u,u)` is a
global maximizer if and only if some `s in {+-1}^k` satisfies

```math
L_(x_i,y_i)=1,
\qquad
\begin{pmatrix}
L_(x_i,x_j)&L_(x_i,y_j)\\
L_(y_i,x_j)&L_(y_i,y_j)
\end{pmatrix}
=-s_is_j\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.    \tag{21.163}
```

Equivalently, with `d_i=(x_i-y_i)/2`,

```math
H_L(x,y)=k-2\left(\sum_i s_id_i\right)^2.           \tag{21.164}
```

Hence a universal one-layer equality lock never isolates the duplicate
relation: for `k>=2` it has nonduplicate two-mismatch ground states.

Two further quantitative ceilings hold.

1. If `q` signed balance squares have no nonzero ternary mismatch in their
   common kernel, then

   ```math
   q\ge {k\log2\over\log(2k+1)}.                    \tag{21.165}
   ```

2. Let a symmetric full sign matrix `W` have a positive
   `lambda`-eigenspace containing at least `2^(alpha k)` Boolean vectors.
   Then `lambda<=sqrt(k/alpha)`.  In a universal repeated-lock architecture
   with `s` size-`k` blocks and cross block `W`, robustness even to every
   exact hollow sign child of norm at most `9sqrt k` requires

   ```math
   s\ge1+{k-1\over\lambda}
    \ge1+(k-1)\sqrt{\alpha/k}.                       \tag{21.166}
   ```

   It therefore uses `Omega(k^(3/2))` vertices at positive Boolean-code
   rate.

#### Proof

The one-spin optimality inequalities at every duplicate force each matching
edge to be positive and every `2 by 2` pair block to have zero row and
column sums.  Write it as `a_(ij)[[1,-1],[-1,1]]`.  A positive triangle
product would allow a three-pair mismatch increasing the energy, so every
triangle product is negative.  The signed complete graph is antibalanced,
`a_(ij)=-s_i s_j`, and direct substitution gives (21.164).  Two opposite
signed mismatches give the tie.

For (21.165), the `2^k` subset sums of the `k` column words must all be
different, while only `(2k+1)^q` integer sums are possible.  For (21.166),
the Boolean intersection of a `d`-dimensional subspace has at most `2^d`
points, so the eigenspace dimension is at least `alpha k`.  Frobenius mass
gives `alpha k lambda^2<=k^2`.  Finally prescribe one child row to have
signed field `-(k-1)` and complete the remaining signing with operator norm
at most `9sqrt k`.  The repeated codeword's local field is
`(s-1)lambda-(k-1)`, which must be nonnegative. `square`

Theorem 21.30 excludes universal duplicate locks and universal repeated
eigenspace locks, not a child-specific alternating-form regularizer.  The
exact proof, Hadamard pullback obligation, verifier, and independent audit
are in
[`drafts/algebraic_exact_sign_locking.md`](drafts/algebraic_exact_sign_locking.md)
and
[`drafts/algebraic_exact_sign_locking_independent_audit.md`](drafts/algebraic_exact_sign_locking_independent_audit.md).

### Theorem 21.31 (universal coordinate pins have quadratic cap)

Fix an old state `u in {+-1}^k`.  Let a single appended quadratic future
have old--new matrix `B`, arbitrary auxiliary quadratic energy `H_C`, and
effective old landscape

```math
g(x)=\max_y\{x^TBy+H_C(y)\}.                        \tag{21.167}
```

Suppose this same future makes `u` projectively `eta`-optimal for every
hollow complete sign child `A`:

```math
H_A(x)+g(x)\le H_A(u)+g(u)+\eta                     \tag{21.168}
```

for every `A,x`.  Put `d_*=floor(k/2)`.  Then

```math
\mathop{osc}g\ge2d_*(k-d_*)-\eta,
\qquad
||B||_(infinity->1)\ge d_*(k-d_*)-\eta/2.           \tag{21.169}
```

Every completed parent
`P_A(x,y)=H_A(x)+x^TBy+H_C(y)` therefore satisfies

```math
Q(P_A)\ge d_*(k-d_*)-\eta/2.                        \tag{21.170}
```

In particular, an exact universal coordinate pin with `N=O(k)` total
vertices has cap `Omega(k^2)`, not `O(N^(3/2))`.

#### Proof

For any `x` at projective Hamming distance `d` from `u`, choose the crossing
edges of `A` so that

```math
H_A(x)-H_A(u)=2d(k-d).                               \tag{21.171}
```

Equation (21.168) gives `g(u)-g(x)>=2d(k-d)-eta`; take `d=d_*`.
Maxima of affine pieces give `osc g<=2||B||_(infinity->1)`.  Finally,
flipping the entire auxiliary shore changes the sign of `x^TBy` and leaves
both internal quadratic energies fixed, so one of the two parent energies
has absolute value at least `|x^TBy|`. `square`

This theorem explains the quadratic calibration in Theorem 21.29 without
claiming that every possible metric embedding uses a common pinned witness.
The audited proof and verifier are in
[`drafts/universal_pin_cap_barrier.md`](drafts/universal_pin_cap_barrier.md)
and
[`drafts/flat_gram_universal_pin_independent_audit.md`](drafts/flat_gram_universal_pin_independent_audit.md).

### Theorem 21.32 (signed quadratic contrasts have full exposed entropy)

Every nonzero homogeneous quadratic Boolean polynomial `P` obeys

```math
Pr\{P>0\}\ge {1\over324},
\qquad
Pr\{P<0\}\ge {1\over324}.                          \tag{21.172}
```

Consequently, for every nonzero hollow symmetric `D`, one may orient
`f=sH_D` and choose antipodal `p,-p` and `X subset {+-1}^k` so that

```math
|X|\ge {2^k\over324},
\qquad
f(p)=f(-p)=-Q(D),
\qquad
f(x)\ge0\quad(x\in X).                              \tag{21.173}
```

If an affine selector envelope with `K` pieces, slope norm at most `D_0`,
and uniform error `eta<Q(D)/2` realizes that signed orientation, then

```math
\log K\ge-\log324+
 { (Q(D)-2\eta)^2\over2D_0^2}.                      \tag{21.174}
```

For every distinct pair in the short-seed flat Gram family of Theorem
21.26, `Q(D)>=c_0k^(3/2)` with `c_0=sqrt2/16`.  If the selectors are
`y in {+-1}^m` with cross block `C`, (21.174) implies

```math
m\ge\left({1\over256\log2}\right)^(1/3)
 k^(2/3)(1-o(1))                                    \tag{21.175}
```

for an arbitrary complete sign cross block.  Under the extra flatness
`||C||_(2->2)<=Lsqrt k`, it implies

```math
m\ge {1\over16L\sqrt{\log2}}k-o(k).                 \tag{21.176}
```

#### Proof

Degree-two hypercontractivity gives `||P||_4<=3||P||_2`.  Interpolation
between `L^1` and `L^4` yields `||P||_1>=||P||_2/9`.  Since `EP=0`, each
of `EP_+` and `EP_-` is at least `||P||_2/18`; Cauchy--Schwarz proves
(21.172).  Orient the larger absolute extremum to be the minimum and take
`X={f>=0}` to obtain (21.173).  The antipodal selector-covering inequality
of Theorem 21.28 gives (21.174).  Finally use

```math
D_0\le||C||_(2->2)\sqrt m,
\qquad ||C||_(2->2)\le\sqrt{km}                    \tag{21.177}
```

and `K<=2^m` to obtain (21.175)--(21.176). `square`

The signed orientation is essential.  For `k=s^2`, partition the variables
into `s` blocks of size `s` and put coefficient two inside each block and
zero between blocks.  Then

```math
||D||_(2->2)<2\sqrt k,
\qquad Q(D)=k^(3/2)-k,                              \tag{21.178}
```

but for every fixed `alpha>0`,

```math
Pr\{|H_D|\ge\alpha k^(3/2)\}
\le\exp\{-\alpha k/4+O(\sqrt k)\}.                 \tag{21.179}
```

Thus flatness plus a large absolute maximum does not give the constant-
deficit exposed bulk needed by (21.174).  Moreover, when the old child
remains outside the selector maximum, the inequality applies to the residual
seen by the selectors, not automatically to the child--query difference.
The exact scope, proof of (21.179), verifier, and independent audit are in
[`drafts/flat_gram_exposed_entropy.md`](drafts/flat_gram_exposed_entropy.md)
and the audit cited after Theorem 21.31.

### Theorem 21.33 (quadratic-character locks lose the leading scale)

Call `phi:{+-1}^k->{+-1}^k` quadratic-character preserving if every
pair product `phi_a phi_b` is, up to sign, a degree-two input character.
For `k>=5`, every such map has the rigid form

```math
\phi(x)=g(x)DPx,                                    \tag{21.180}
```

where `g(x)` is an arbitrary common Boolean gauge, `D` is a sign diagonal,
and `P` is a permutation.

For every fixed complete sign bridge `R`, define the defect between its free
new-shore optimum and the best globally gauged intended copy by

```math
\Delta_(R,\phi)(x)=||R^Tx||_1-|x^TR\phi(x)|.        \tag{21.181}
```

Then

```math
\max_x\Delta_(R,\phi)(x)
\ge k\sqrt{k/2}-\sqrt3 k.                           \tag{21.182}
```

Bounded coordinate replication cannot remove this loss.  If
`phi_a(x)=g(x)s_ax_(pi(a))`, where `pi:[m]->[k]` is onto and every fibre has
size at most `L`, then for every `k by m` complete sign bridge

```math
\max_x\{||R^Tx||_1-|x^TR\phi(x)|\}
\ge m\sqrt{k/2}-\sqrt{m^2+2kLm}.                   \tag{21.183}
```

Thus `m=Theta(k)` and `L=O(1)` still incur an
`Omega(k^(3/2))` worst defect.

Finally, in the bare one-layer pin architecture, suppose a target `u` is a
global maximizer of `H_A(x)+||R^Tx||_1` for every exact hollow sign child
with operator norm at most `9sqrt k`.  Then necessarily

```math
R=us^T,
\qquad ||R||_(infinity->1)=k^2,                    \tag{21.184}
```

so every complete parent containing that bridge has cap at least `k^2`.

#### Proof

Fix one output coordinate.  The supports of its pair characters with all
other outputs are distinct pairwise-intersecting two-subsets.  A family of
at least four such sets is a star; symmetric difference reconstructs all
remaining pair supports, and the signs satisfy a one-cocycle identity.  This
gives (21.180).

The gauge disappears inside the absolute value in (21.181), leaving
`|x^TCx|` for a sign matrix `C`.  Under uniform `x`, Khintchine gives

```math
E||R^Tx||_1\ge k\sqrt{k/2},                         \tag{21.185}
```

while Fourier orthogonality gives

```math
E(x^TCx)^2
=(tr C)^2+\sum_(i<j)(C_(ij)+C_(ji))^2<3k^2.        \tag{21.186}
```

Averaging proves (21.182).  With replication, collapse each fibre into a
column of `C`; then `||C||_F^2<=kLm` and `|tr C|<=m`, proving (21.183).

For the final claim, prescribe for each coordinate a bounded-operator sign
child whose signed local field at `u` is `-(k-1)`.  Comparing `u` with its
one-bit flip forces

```math
||R^Tu||_1-||R^Tu^(i)||_1\ge2(k-1).                 \tag{21.187}
```

Columnwise, every summand is in `{-2,0,2}`.  The inequalities for all `i`
force every column to agree with `u` up to one column sign, hence (21.184).
The shore-flip argument gives the cap. `square`

The theorem closes dominant equality-lock proofs for the **full** quadratic
character algebra.  It does not cover query-dependent coordinate pins,
which succeed in Theorem 21.29, or a nonlinear pullback proved only for the
narrower alternating-form language.  The full audited proof and verifier
are in
[`drafts/exact_sign_disjoint_compiler.md`](drafts/exact_sign_disjoint_compiler.md)
and
[`drafts/exact_sign_disjoint_compiler_independent_audit.md`](drafts/exact_sign_disjoint_compiler_independent_audit.md).

### Theorem 21.34 (bounded cap prices witness diversity, not fibre size)

Bounded cap does not force child-dependent optimizer switching.  For every
`u in {+-1}^k`, append one spin with incident signs `u_i`, so the effective
future is

```math
g_u(x)=|u\mathbin\cdot x|.                          \tag{21.188}
```

There are at least

```math
2^{ {k\choose2}-k}                                 \tag{21.189}
```

distinct complete sign children `A` such that

```math
Q(A)\le2k^(3/2),
\qquad Q(P_A)\le2k^(3/2)+k,                        \tag{21.190}
```

and the same old states `+-u` optimize `H_A+g_u` for every child.

The response-sensitive replacement is quantitative.  Let `mathcal F` be a
family with `Q(A)<=C_0k^(3/2)`, pairwise separated over a declared future
language by

```math
\sup_g|R_A(g)-R_(A')(g)|\ge\epsilon k^(3/2).        \tag{21.191}
```

Suppose one dictionary `U subset {+-1}^k` answers every child--future pair
within `tau k^(3/2)`:

```math
0\le R_A(g)-\max_(u\in U)\{H_A(u)+g(u)\}
\le\tau k^(3/2),
\qquad 2\tau<\epsilon.                              \tag{21.192}
```

Then

```math
|\mathcal F|
\le\left(1+\left\lceil{2C_0\over\epsilon-2\tau}
                 \right\rceil\right)^{|U|}.        \tag{21.193}
```

In particular an `exp(alpha k)` bounded-cap response packing requires
`|U|=Omega(k)`.  If the language has futures `g_1,...,g_q` and `U_j` is the
set of exact old witnesses used under `g_j`, then

```math
\sum_(j=1)^q|U_j|=\Omega(k).                        \tag{21.194}
```

Thus `q=O(1)` forces extensive child-dependent witness switching in some
context, while one common optimizer per context forces `q=Omega(k)`.

#### Proof

Vertex switching acts freely on complete signings modulo the global sign,
giving `2^{binom(k,2)-k+1}` orbits.  Every orbit has a representative whose
one-sided quadratic maximum is attained at `u`.  A random signing has
`Q(A)<=2k^(3/2)` with probability above one half by Hoeffding and a union
bound over the cube; the event is switching invariant.  Select one
representative from each bounded orbit.  Both `H_A` and (21.188) are
maximized at `u`, and the appended linear term has magnitude at most `k`,
proving (21.189)--(21.190).

For the second claim, truncated maxima are nonexpansive:

```math
|R_A^U(g)-R_(A')^U(g)|
\le\max_(u\in U)|H_A(u)-H_(A')(u)|.                 \tag{21.195}
```

Equations (21.191)--(21.192) therefore make the evaluation vectors
`(H_A(u))_(u in U)` an `ell_infinity` packing at separation
`(epsilon-2tau)k^(3/2)`.  Each coordinate lies in an interval of length
`2C_0k^(3/2)`; bin it to obtain (21.193).  Taking the union of the exact
witness supports gives (21.194). `square`

For an uncalibrated hollow quadratic parent, parent cap implies the child-cap
hypothesis by averaging over all auxiliary spins.  Projective and signed
absolute variants obey the same law after anchoring one evaluation
coordinate or adjoining the sign channel.  Full details and the independent
audit are in
[`drafts/bounded_cap_optimizer_switching.md`](drafts/bounded_cap_optimizer_switching.md)
and
[`drafts/bounded_cap_optimizer_switching_independent_audit.md`](drafts/bounded_cap_optimizer_switching_independent_audit.md).

### Theorem 21.35 (Hadamard synchronization has an exact deficit but no universal planted-witness stability)

Let `W in {+-1}^{k times k}` satisfy `W^TW=kI`, and put

```math
v_x={W^Tx\over\sqrt k}\qquad(x\in\{+-1\}^k).
```

Then the bridge deficit is exactly

```math
k^(3/2)-x^TWy
={\sqrt k\over2}||y-v_x||_2^2,                    \tag{21.196}
```

and, after optimizing the new shore,

```math
k^(3/2)-||W^Tx||_1
={1\over2\sqrt k}\sum_a
 (|(W^Tx)_a|-\sqrt k)^2.                          \tag{21.197}
```

Consequently, if a future energy `K` extends to the radius-`sqrt(k)` sphere
with Euclidean Lipschitz constant `L`, then

```math
\max_y\{x^TWy+K(y)\}
\le k^(3/2)+K(v_x)+{L^2\over2\sqrt k}.             \tag{21.198}
```

This is a useful positive synchronization criterion, but a Boolean cap bound
does not imply it at a subleading scale.  For all sufficiently large
`k=2^(2m)` there are complete hollow signings `A,C` and the symmetric
Sylvester--Walsh matrix `W` such that

```math
Q(A)\le3k^(3/2),\qquad Q(C)={1\over2}k^(3/2),      \tag{21.199}
```

the complete exact-sign parent

```math
P=\begin{pmatrix}A&W\\W&C\end{pmatrix}
```

satisfies `Q(P)<=9k^(3/2)/2`, yet a state off one planted Boolean pullback
relation beats its planted relation state by at least `k^(3/2)/8`.

#### Proof

Equation (21.196) follows by expanding `||y-v_x||_2^2`; coordinatewise
optimization and `||W^Tx||_2^2=k^2` give (21.197).  Completing the square in
the distance variable proves (21.198).

For the obstruction, write `q=sqrt(k)` and take the self-dual bent vector

```math
x_0(u,v)=(-1)^(u\mathbin\cdot v),\qquad Wx_0=qx_0.
```

If `x_S` flips `d<q/2` coordinates of `x_0`, its bridge optimizer remains
`x_0` and its exact bridge deficit is `2dq`.  Put
`d=floor(q/4)`, prescribe the old edges across `(S,S^c)` as
`-x_0(i)x_0(j)`, and fill both principal shores by signings of cap at most
twice their order to the power `3/2`.  Then

```math
H_A(x_S)-H_A(x_0)=2d(k-d),
```

while `Q(A)<=3k^(3/2)`.  Take `C` to be the hollow part of `W`.  Its trace
is zero, so on Boolean states `H_C(y)=y^TWy/2` and
`Q(C)=k^(3/2)/2`.  The two compared parent states use the same new spin;
hence

```math
H_P(x_S,x_0)-H_P(x_0,x_0)
=2d(k-d-q)\ge k^(3/2)/8.                           \tag{21.200}
```

Triangle inequality gives the stated parent cap. `square`

The last assertion is deliberately about **planted-witness stability**.  It
does not assert that every global optimizer lies outside the entire Hadamard
relation set.  The complete proof, verifier, and independent audit are in
[`drafts/hadamard_bridge_synchronization.md`](drafts/hadamard_bridge_synchronization.md),
[`experiments/verify_hadamard_bridge_synchronization.py`](experiments/verify_hadamard_bridge_synchronization.py),
and
[`drafts/hadamard_bridge_synchronization_independent_audit.md`](drafts/hadamard_bridge_synchronization_independent_audit.md).

### Theorem 21.36 (a bounded-cap exact-sign contextual metric compiler)

There is an absolute `gamma>0` and, for every sufficiently large
`n=q^2` with `q` a power of two, a family `mathcal S subset {+-1}^n` of
size at least `exp(gamma n)` with the following properties.

Let `mathcal H` be the regular symmetric Walsh signing satisfying

```math
\mathcal H^2=nI,\qquad \mathcal H\mathbf1=q\mathbf1,
\qquad \operatorname{tr}\mathcal H=0,             \tag{21.201}
```

put `A=mathcal H-diag(mathcal H)` and
`A_s=D_sAD_s`.  Every child is an exact hollow signing with

```math
Q(A_s)={1\over2}n^(3/2).                           \tag{21.202}
```

For every declared query `t in mathcal S`, append `q=sqrt(n)` spins, use
the exact cross block `t1_q^T`, and put a public positive clique on the new
shore.  If `F_s(t)` is the absolute cap of this order-`N=n+sqrt(n)` parent,
then

```math
F_s(s)={3\over2}n^(3/2)+{q\choose2},               \tag{21.203}
```

whereas, for `s!=t`,

```math
F_s(t)\le{11\over8}n^(3/2)+{q\choose2}.           \tag{21.204}
```

All compiled parents therefore have cap `O(N^(3/2))`.  If

```math
d_0(s,t)=Q(A_s-A_t),
```

which is the projective contextual metric of the same-support language
`r -> Q(A_s-A_r)`, and

```math
d_C(s,t)={1\over2}\operatorname{osc}_{r\in\mathcal S}
               (F_s(r)-F_t(r)),                   \tag{21.205}
```

then

```math
{1\over8}d_0(s,t)\le d_C(s,t)\le d_0(s,t).        \tag{21.206}
```

Thus exact signs, disjoint edge ownership, linear total order, and bounded
parent cap are compatible with constant-distortion compilation of a
linear-rate contextual metric.  The optimizer is allowed to depend jointly
on the child and query.

#### Proof

A Rademacher Hanson--Wright bound gives

```math
\Pr\{|w^T\mathcal Hw|>qn/4\}\le2e^(-c n).
```

A union bound selects `exp(gamma n)` switches such that every distinct pair
`s,t`, with `w=s odot t`, obeys

```math
|w^T\mathcal Hw|\le qn/4.                          \tag{21.207}
```

After switching the old spin by `s`, omit the new clique and split the
absolute cap into its two signs.  The two channels are

```math
R_\pm(w)=\max_u\left\{
 \mathord\pm{1\over2}u^T\mathcal Hu+qw^Tu\right\}. \tag{21.208}
```

Because `mathcal H^2=q^2I`, completing squares with
`2qI minus-or-plus mathcal H` gives

```math
R_+(w)\le qn\left(1+{2+\rho(w)\over6}\right),
\quad
R_-(w)\le qn\left(1+{2-\rho(w)\over6}\right),     \tag{21.209}
```

where `rho(w)=w^T mathcal H w/(qn)`.  Equation (21.207) gives the
`11/8` off-diagonal ceiling.  On the diagonal, `u=1` and aligned new spins
attain the sum of the three separate caps, proving (21.203).

At query `s`, the response difference between child `s` and any other child
is at least `n^(3/2)/8`; at query `t` its sign reverses.  Hence
`d_C>=n^(3/2)/8`.  The cap is one-Lipschitz in the child Hamiltonian, so
`d_C<=d_0`, while the triangle inequality gives `d_0<=n^(3/2)`.  This proves
(21.206). `square`

The `sqrt(n)` query shore is order-optimal within repeated rank-one query
bridges.  With only `m` repeated columns, changing a query changes the
projective profile by at most `2nm` when the auxiliary block is public, and
at most `2nm+m(m-1)` if that block is query-owned.  Thus a fixed
`n^(3/2)` response gap forces `m=Omega(sqrt(n))`.

This is a one-layer metric embedding, not yet a reusable compositional
congruence or a coordinatewise reconstruction of the old response table.
Its scope is one linear-rate switching subcode of a regular-Hadamard family,
not arbitrary near-minimizers.  The complete proof, verifier, and audit are
in
[`drafts/bounded_cap_contextual_metric_compiler.md`](drafts/bounded_cap_contextual_metric_compiler.md),
[`experiments/verify_bounded_cap_contextual_metric_compiler.py`](experiments/verify_bounded_cap_contextual_metric_compiler.py),
and
[`drafts/bounded_cap_contextual_metric_compiler_independent_audit.md`](drafts/bounded_cap_contextual_metric_compiler_independent_audit.md).

### Theorem 21.37 (spectral anti-pins have a finite-port Gram certificate)

Let `H` be a symmetric sign matrix satisfying

```math
H\mathbf1=r\mathbf1,qquad ||H||\le r,qquad
\operatorname{tr}H=0,                              \tag{21.210}
```

and put `A=H-diag(H)`, `A_s=D_sAD_s`.  Then
`Q(A_s)=rn/2`.  Fix an integer `m` with `2m>r`, and for
`w in {+-1}^n` define

```math
\Psi_\sigma(w)={m\over2n}
 w^T(2mI-\sigma H)^{-1}w,qquad\sigma\in\{+-1\}.   \tag{21.211}
```

If a switch code `mathcal S` obeys, for distinct `s,t`,

```math
\max_\sigma\Psi_\sigma(s\odot t)
\le {r\over2m}-\delta,                             \tag{21.212}
```

then the repeated-column query `t1_m^T`, completed by a public positive
clique, gives a bounded-cap exact-sign compiler whose projective response
metric satisfies

```math
{\delta m\over r}Q(A_s-A_t)
\le d_C(s,t)\le Q(A_s-A_t).                        \tag{21.213}
```

If additionally `H^2=r^2I` and `m=r`, the two features collapse to the one
Rayleigh coordinate `rho(w)=w^THw/(rn)`:

```math
\Psi_+(w)={2+\rho(w)\over6},\qquad
\Psi_-(w)={2-\rho(w)\over6}.                       \tag{21.214}
```

Thus every two-sided Rayleigh code `|rho|<=theta<1` works with
`delta=(1-theta)/6`; the Walsh realization is not essential.

For a fixed collection of `ell` repeated ports `w_1,...,w_ell`, the exact
trust-region spherical relaxation depends only on the two Gram matrices

```math
G_(ab)={w_a^Tw_b\over n},
\qquad R_(ab)={w_a^THw_b\over rn}.                 \tag{21.215}
```

Indeed, for each sign vector `epsilon` and channel `sigma` it equals

```math
rn\inf_(\alpha>1/2)\left\{
\alpha+{(m/r)^2(2\alpha\epsilon^TG\epsilon
                  +\sigma\epsilon^TR\epsilon)
             \over2(4\alpha^2-1)}\right\}.        \tag{21.216}
```

At fixed `ell`, this is an `O(ell^2 log n)`-bit one-layer certificate; any
public exact-sign completion on the `ell m` auxiliary spins adds only
`O_ell(n)` when `m=Theta(sqrt(n))`.

#### Proof

After switching the old spin, the two absolute channels are

```math
\max_u\left\{{\sigma\over2}u^THu+mw^Tu\right\}.
```

Completing the square with `2mI-sigma H` bounds this by
`mn(1+Psi_sigma(w))`.  The diagonal query attains
`rn/2+mn` before the public completion, so (21.212) gives opposite gaps
`delta mn` at the two diagonal coordinates.  Lipschitzness in the old child
and `Q(A_s-A_t)<=rn` prove (21.213).  The inverse formula

```math
(2rI\mp H)^{-1}={2rI\pm H\over3r^2}
```

proves (21.214).

For several ports, endpoint optimization produces
`v=sum_a epsilon_a w_a`.  The exact Lagrange dual for maximizing a quadratic
plus a linear form on `||u||_2^2=n`, followed by

```math
(2\alpha rI-\sigma H)^{-1}
={2\alpha rI+\sigma H\over r^2(4\alpha^2-1)},
```

gives (21.216). `square`

The Gram state is a spherical certificate, not an exact Boolean-cap
reconstruction: it may retain a fixed integrality gap.  It is also not yet a
congruence.  For tensor powers of one order-16 regular Walsh matrix there are
orthogonal Boolean `+sqrt(n)` eigenvectors `1,v` with identical separate
one-port states, but

```math
\mathcal Q(\mathbf1,\mathbf1)-\mathcal Q(\mathbf1,v)
\ge(2-\sqrt2)n^(3/2),                               \tag{21.217}
```

where `mathcal Q` is the two-port cap before the lower-order completion.
The full Gram pair sees the newly created cross entries; the product of the
separate one-port states does not.  Full details and audit are in
[`drafts/spectral_antipin_feature_algebra.md`](drafts/spectral_antipin_feature_algebra.md),
[`experiments/verify_spectral_antipin_feature_algebra.py`](experiments/verify_spectral_antipin_feature_algebra.py),
and
[`drafts/spectral_antipin_feature_algebra_independent_audit.md`](drafts/spectral_antipin_feature_algebra_independent_audit.md).

### Theorem 21.38 (one-port orientation blindness and macroscopic holonomy)

Let `n=q^2` and let `mathcal H` be a regular symmetric Walsh signing with
`mathcal H^2=nI`, `mathcal H1=q1`, and trace zero.  Put
`A=mathcal H-diag(mathcal H)`.  For **every** rank-one anti-pin query `t`
using `q` appended spins and a public positive clique, the cap responses of
`A` and `-A` satisfy

```math
|F_A(t)-F_(-A)(t)|\le q(q-1)<n.                    \tag{21.218}
```

Nevertheless, join either orientation to one fixed `A` shore through the
public exact-sign bridge `mathcal H`.  The two order-`2n` parents

```math
P_+(x,z)=H_A(x)+H_A(z)+x^T\mathcal Hz,
\quad
P_-(x,z)=H_A(x)-H_A(z)+x^T\mathcal Hz              \tag{21.219}
```

obey

```math
Q(P_+)=2n^(3/2),qquad Q(P_-)\le\sqrt2n^(3/2).      \tag{21.220}
```

Thus a relative orientation bit invisible to every one-port query becomes
observable after dense composition with gap
`(2-sqrt(2))n^(3/2)`.

There is a matching positive narrow-continuation law.  For any old child
`A`, arbitrary old--new sign block `B`, and arbitrary complete auxiliary
signing `C` on `m` spins, let

```math
R_\sigma=\max_(x,y)|\sigma H_A(x)+x^TBy+H_C(y)|.
```

Then

```math
|R_+-R_-|\le2Q(C)\le m(m-1).                      \tag{21.221}
```

Hence orientation is a valid vanishing-distortion quotient for every
quadratic continuation flattenable onto total auxiliary width
`m=o(n^(3/4))`, but not for a second macroscopic shore.

#### Proof

With `C` omitted, replacing `y` by `-y` and then negating the expression
makes the two orientation caps identical.  Adding `C` changes either cap by
at most `Q(C)`, proving (21.218) and (21.221).

The real symmetric block matrices of (21.219) are

```math
\begin{pmatrix}1&1\\1&1\end{pmatrix}\otimes\mathcal H,
\qquad
\begin{pmatrix}1&1\\1&-1\end{pmatrix}\otimes\mathcal H.
```

Their operator norms are `2q` and `sqrt(2)q`.  The spectral bound on a
Boolean vector of squared norm `2n` gives (21.220), and the first bound is
attained at `(1,1)`. `square`

The result is a global-cap congruence falsifier, not a planted-witness
failure.  Simultaneous switching makes it gauge-covariant over the
exponential BCX code.  It does not cover adaptive re-encoding or
unbounded-depth narrow composition.  The draft, verifier, and audit are in
[`drafts/bcx_two_port_holonomy.md`](drafts/bcx_two_port_holonomy.md),
[`experiments/verify_bcx_two_port_holonomy.py`](experiments/verify_bcx_two_port_holonomy.py),
and
[`drafts/bcx_two_port_holonomy_independent_audit.md`](drafts/bcx_two_port_holonomy_independent_audit.md).

### Theorem 21.39 (the sharp interface scale for orientation visibility)

Let `H` be any even landscape on `n` Boolean spins.  For an arbitrary sign
bridge `B` to `m` new spins and an arbitrary even future `K`, define

```math
R_\sigma=\max_(x,y)|\sigma H(x)+x^TBy+K(y)|.
```

Then

```math
|R_+-R_-|\le2||K||_\infty.                        \tag{21.222}
```

Consequently, if `K=H_C` and `Q(C)<=K_0m^alpha`, an
`epsilon n^(3/2)` orientation gap requires

```math
m\ge\left({\epsilon\over2K_0}\right)^{1/\alpha}
       n^{3/(2\alpha)}.                            \tag{21.223}
```

For unrestricted complete sign continuations `alpha=2`, the exponent
`3/4` is sharp even when the final parent must have natural cap.  More
precisely, along regular-Hadamard orders there are complete exact-sign
parents on

```math
N=n+\lfloor n^(3/4)\rfloor
```

spins such that

```math
Q(P_+)\ge(2-o(1))n^(3/2),
\qquad
Q(P_-)\le(5/4+o(1))n^(3/2).                       \tag{21.224}
```

Both caps are `O(N^(3/2))`.  If instead every future shore is itself
bounded-cap, `Q(C)=O(m^(3/2))`, (21.223) forces `m=Omega(n)`.

#### Proof

With `K` omitted, invert `y` and negate the expression to identify the two
orientation caps.  Adding `K` perturbs either by at most
`||K||_infinity`, proving (21.222)--(21.223).

For sharpness, write `S=n^(3/2)`, put
`m=floor(n^(3/4))`, `a=m/n`, and choose an exact sign bridge

```math
B=aJ+E,
\qquad ||E||=O(\sqrt n).                           \tag{21.225}
```

Such a bridge exists by taking independent signs of mean `a` and applying a
rectangular subgaussian norm bound.  Use a regular symmetric Hadamard old
block `mathcal H` and a positive clique on the new shore.  At
`x=y=1`, the positive orientation has energy

```math
S/2+3m^2/2-o(S)=(2-o(1))S.                        \tag{21.226}
```

For the negative orientation set
`p=(1^Tx)/n`, `s=(1^Ty)/m`, and `lambda=m^2/S<=1`.  The two outer absolute
channels, after division by `S`, have the scalar upper envelopes

```math
{1\over2}-p^2+\lambda ps+{\lambda\over2}s^2
\le {5\over4},
\qquad
{1\over2}+\lambda(-ps-s^2/2)\le1.                \tag{21.227}
```

The error bridge contributes at most
`O(sqrt(n)sqrt(nm))=O(n^(11/8))=o(S)`.  This proves
(21.224); termwise caps also give the final-parent bound. `square`

Thus approximate contextual equivalence must be indexed by a continuation
**cap budget**, not width alone.  The construction spends quadratic cap on
the sublinear clique; it is not a universal pin or a statement about
near-minimizers.  The full proof, exact finite checks, and audit are in
[`drafts/orientation_visibility_threshold.md`](drafts/orientation_visibility_threshold.md),
[`experiments/verify_orientation_visibility_threshold.py`](experiments/verify_orientation_visibility_threshold.py),
and
[`drafts/orientation_visibility_threshold_independent_audit.md`](drafts/orientation_visibility_threshold_independent_audit.md).

### Theorem 21.40 (PSD geometry limits macroscopic cross-Gram dimension)

Let `K_0` be a positive semidefinite `p by p` matrix with diagonal in
`[0,1]`.  For distinct off-diagonal pairs `e_a={i_a,j_a}`, put
`D_e=e_ie_j^T+e_je_i^T`.  If

```math
K_\sigma=K_0+\sum_(a=1)^h\sigma_a\eta_aD_(e_a)
\succeq0
\qquad(\sigma\in\{+-1\}^h),                       \tag{21.228}
```

then

```math
\sum_a\eta_a^2\le {p\over2}.                     \tag{21.229}
```

Consequently, let `J` be a symmetric involution and let `p` port vectors
have squared norm `n`.  Define

```math
G_(ij)={w_i^Tw_j\over n},
\qquad R_(ij)={w_i^TJw_j\over n}.                 \tag{21.230}
```

Suppose every word in a full affine coordinate cube is realizable with fixed
self data and

```math
G_\sigma=G_0+\sum_a\sigma_ag_aD_(e_a),
\qquad
R_\sigma=R_0+\sum_a\sigma_ar_aD_(e_a).            \tag{21.231}
```

Then

```math
\sum_a(g_a^2+r_a^2)\le2p.                         \tag{21.232}
```

In particular, coordinatewise independently toggleable features of raw
amplitude at least `delta` have number at most `2p/delta^2`.

This order is sharp for Boolean ports.  Tensor powers of the order-16
regular Walsh matrix contain `n^(1/4)` mutually orthogonal Boolean top
eigenvectors; pairing duplicates versus distinct vectors realizes
`floor(p/2)` independent constant-amplitude Gram--Rayleigh bits with the
same one-port self state.

If, additionally, a fixed-arity declared response changes by at most
`L n^(3/2)sqrt(g_a^2+r_a^2)` when bit `a` changes, and the materialized total
order is `N`, then `epsilon N^(3/2)`-visible independent bits obey

```math
h\le {2L^2\over\epsilon^2}
       p\left({n\over N}\right)^3.                \tag{21.233}
```

#### Proof

For positive definite `K_0` with diagonal one, conjugate the perturbation by
`K_0^(-1/2)`.  Since both sign words occur, every resulting symmetric
matrix `T_sigma` has spectrum in `[-1,1]`, so
`tr(T_sigma^2)<=p`.  Averaging over signs gives

```math
p\ge2\sum_a\eta_a^2
  (M_(i_ai_a)M_(j_aj_a)+M_(i_aj_a)^2),
\qquad M=K_0^{-1}.                                 \tag{21.234}
```

Every diagonal entry of the inverse of a correlation matrix is at least
one, proving (21.229).  Adding `epsilon I`, normalizing the diagonal, and
letting `epsilon` vanish handles singular `K_0`.

The two spectral sectors

```math
K^+=(G+R)/2,\qquad K^-=(G-R)/2                   \tag{21.235}
```

are Gram matrices of the projections `(I+-J)/2`.  Apply (21.229) to both;
their perturbation amplitudes are `(g_a+-r_a)/2`, which proves (21.232).
The tensor construction supplies the lower example.  Finally a visible bit
must have raw amplitude at least
`(epsilon/L)(N/n)^(3/2)`; substitution in (21.232) gives (21.233). `square`

The theorem concerns full coordinatewise affine cubes.  It does not bound
dense affine directions, nonlinear codes, or collective metric entropy, and
order normalization alone does not prove that all local pair gaps survive
simultaneous port activation.  Thus the exact `p^2` Gram table is not an
information lower bound, but an `O(p)` approximate sufficient state is still
open.  The full proof, verifier, and audit are in
[`drafts/cross_gram_macroscopic_dimension.md`](drafts/cross_gram_macroscopic_dimension.md),
[`experiments/verify_cross_gram_macroscopic_dimension.py`](experiments/verify_cross_gram_macroscopic_dimension.py),
and
[`drafts/cross_gram_macroscopic_dimension_independent_audit.md`](drafts/cross_gram_macroscopic_dimension_independent_audit.md).

### Theorem 21.41 (the exact orientation--cycle gluing carrier)

Fix a regular symmetric Hadamard matrix `H` of order `n`, and let a finite
graph `G=(V,E)` carry onsite orientations `sigma_i in {+-1}` and bridge
signs `b_ij in {+-1}`.  Put

```math
T_(ii)=\sigma_i,
\qquad T_(ij)=b_(ij)\ (ij\in E),
\qquad T_(ij)=0\ (ij\notin E).                    \tag{21.236}
```

The corresponding complete block energy is

```math
\mathcal E_T(X)={1\over2}X^T(T\otimes H)X.        \tag{21.237}
```

Its entire Boolean energy multiset up to global sign, hence its absolute
cap, factors through the projective switching quotient

```math
T\sim\epsilon DTD,
\qquad \epsilon\in\{+-1\},\quad D^2=I.           \tag{21.238}
```

If `G` has `e` edges and `c` connected components, this labelled quotient
has exactly

```math
2^(e+c-1)                                           \tag{21.239}
```

classes.  A spanning forest presents its `e+c-1` coordinates as `|V|-1`
relative onsite orientations and `e-|V|+c` modified fundamental-cycle
products.

The quotient has an exact gluing law.  Fix the carrier of each of `s`
connected pieces and join them by `r` cross edges whose piece graph is
connected.  Every fibre of joined carriers has exactly

```math
2^r,                                                \tag{21.240}
```

with intrinsic dimension split

```math
r=(s-1)+(r-s+1)                                    \tag{21.241}
```

between relative marginal antipodes and new cross-cycle holonomies.  After
choosing marginal gauge sections and a cross-edge forest, these bits and the
piece carriers reconstruct the joined carrier.  Repeated gluing is
associative.

Both kinds of compatibility information can affect the Boolean cap at
leading order.  One bridge gives the orientation gap
`(2-sqrt(2))n^(3/2)` of Theorem 21.38.  On a triangle, positive versus
negative bridge product gives

```math
Q_H(J_3)={9\over2}n^(3/2),
\qquad
Q_H(T_\mathrm{unbal})\le3n^(3/2),                 \tag{21.242}
```

so one cycle bit changes the cap by at least `3n^(3/2)/2`.

#### Proof

Blockwise spin inversion sends `T` to `DTD`; the outer absolute value
identifies `T` with `-T`.  Vertex switching has effective rank
`|V|-c` on edge signs, while the global antipode acts freely because all
onsite signs are nonzero.  Quotienting the `|V|+e` coefficient bits proves
(21.239), and forest normalization gives the displayed coordinates.

The connected joined graph has `sum_a e_a+r` carrier bits while the fixed
connected marginals have `sum_a e_a`; restriction is onto, proving
(21.240).  Relative choices of marginal antipode give `s-1` bits, a cross
forest normalizes `s-1` edge signs, and each remaining edge closes one
independent cycle.  This proves (21.241) and exact reconstruction.

For the unbalanced triangle, the scalar block matrix has operator norm two;
the Boolean vector has squared norm `3n`, giving the second bound in
(21.242).  The balanced regular pole saturates every term and gives the first.
`square`

This is minimal as a coefficient-conjugacy carrier, not a claim that one
bare scalar cap is injective on every quotient class.  Spectrum is also too
coarse: an explicit order-four pair of cospectral block matrices has exact
caps `32` and `34`.  The theorem is a strict state only for the common-factor
orientation family, not arbitrary dense bridges.  The full proof and two
independent finite checks are in
[`drafts/regular_hadamard_orientation_carrier.md`](drafts/regular_hadamard_orientation_carrier.md),
[`experiments/verify_regular_hadamard_orientation_carrier.py`](experiments/verify_regular_hadamard_orientation_carrier.py),
[`drafts/regular_hadamard_orientation_carrier_independent_audit.md`](drafts/regular_hadamard_orientation_carrier_independent_audit.md),
and
[`experiments/verify_regular_hadamard_orientation_carrier_independent_audit.py`](experiments/verify_regular_hadamard_orientation_carrier_independent_audit.py).

### Theorem 21.42 (collective Gram states have linear response entropy)

For two positive spectral sectors on `p` labelled ports, put

```math
\mathcal K_p=\{(K^+,K^-):K^\pm\succeq0,
                         \operatorname{tr}K^\pm\le p\},       \tag{21.243}
```

write `G=K^++K^-`, `R=K^+-K^-`, and define the collective query
pseudometric

```math
d_q((G,R),(G',R'))={1\over p^2}
 \max_{\epsilon\in\{+-1\}^p,\sigma\in\{+-1\}}
 \left|\epsilon^T(\Delta G+\sigma\Delta R)\epsilon\right|.
                                                               \tag{21.244}
```

For every `0<eta<=1`,

```math
\log\operatorname{Cov}(\eta,\mathcal K_p,d_q)
\le2p\left\lceil{4\over\eta}\right\rceil
       \log\left(1+{16\over\eta}\right).             \tag{21.245}
```

The centres may be chosen PSD of rank at most `ceil(4/eta)` in each
sector.  Conversely, for every fixed `0<theta<1/2`, Boolean ports in one
regular-Hadamard top eigenspace contain

```math
2^{(1-H_2(\theta)-o(1))p}                           \tag{21.246}
```

states separated by at least `8theta(1-theta)` in `d_q`.  Thus fixed-scale
metric entropy is `Theta_eta(p)`, despite the exact table having
`Theta(p^2)` entries.

This metric controls a genuine collective response.  Put

```math
\kappa={pm\over r},\qquad
a={\epsilon^TK^\sigma\epsilon\over p^2},\qquad
b={\epsilon^TK^{-\sigma}\epsilon\over p^2}.          \tag{21.247}
```

The normalized spherical trust response in that channel is

```math
\Psi_\kappa(a,b)=\inf_{t>0}\left\{
 {1+t\over2}+{\kappa^2\over2}
 \left({a\over t}+{b\over t+2}\right)\right\}.       \tag{21.248}
```

If two states have `d_q<=eta`, their entire labelled trust-response tables,
and hence their maxima, differ by at most

```math
\kappa\sqrt{\eta/2}+{\kappa^2\eta\over8}.            \tag{21.249}
```

The square-root exponent is sharp at the boundary where the dangerous
sector has zero mass.  If all relevant trust minimizers instead obey
`t>=tau>0`, the right side improves to

```math
{\kappa^2\eta\over4}
\left({1\over\tau}+{1\over\tau+2}\right).            \tag{21.250}
```

At total repeated-port mass `pm<=r`, an error-`epsilon` spherical response
carrier therefore has

```math
\log N_\epsilon
=O\left({p\over\epsilon^2}\log{1\over\epsilon}\right). \tag{21.251}
```

This dependence on `p` cannot be removed: the Boolean rank-one family in
(21.246) has labelled response-table separation `2kappa theta` when
`kappa` is bounded below.

#### Proof

For one sector, discard eigenvalues at most `eta p/4`.  The discarded
operator norm contributes at most `eta/4` in the normalized Boolean
quadratic metric, while the retained rank is at most `4/eta`.  Factor the
retained part as `BB^T`, `||B||_F<=sqrt(p)`, and cover its Euclidean factor
ball at radius `eta sqrt(p)/8`.  The identity

```math
||BB^T-CC^T||_{op}
\le(||B||_{op}+||C||_{op})||B-C||_{op}              \tag{21.252}
```

gives the other `eta/4`; taking two sector nets proves (21.245).

For the lower family, take ports `w_i=s_iw`, where `w` is a Boolean top
eigenvector and `s in {+-1}^p`.  Then `K_s^+=ss^T`, `K_s^-=0`, and for
projective Hamming distance `h`,

```math
d_q(s,t)={8h(p-h)\over p^2}.                        \tag{21.253}
```

Greedy projective Hamming packing proves (21.246).

Formula (21.248) is the exact one-dimensional trust dual.  Coordinatewise
monotonicity and comparison at `t'=kappa sqrt(e)` show that changing the
dangerous coordinate by `e` costs at most `kappa sqrt(e)`, while changing
the safe coordinate costs at most `kappa^2e/4`.  Definition (21.244) gives
`e=eta/2`, proving (21.249).  Restricting the comparison to `t>=tau` proves
(21.250).  Finally use (21.245) at radius `epsilon^2`; the rank-one packing
has exact response `1/2+kappa|epsilon^Ts|/p`. `square`

The upper theorem is for the relaxed PSD/spherical carrier.  It neither
rounds it to exact Boolean old spins nor supplies physical low-rank centres
inside every realizable subset.  At the original one-port anti-pin scaling
`m=r`, one has `kappa=p`, so fixed response accuracy requires much finer
metric resolution; the linear fixed-accuracy statement applies at bounded
total port mass, not at arbitrary simultaneous amplification.  Full proofs,
checks, and independent audits are in
[`drafts/cross_gram_response_metric_entropy.md`](drafts/cross_gram_response_metric_entropy.md),
[`drafts/cross_gram_response_metric_entropy_independent_audit.md`](drafts/cross_gram_response_metric_entropy_independent_audit.md),
[`drafts/collective_cross_gram_packing_and_response_modulus.md`](drafts/collective_cross_gram_packing_and_response_modulus.md),
[`drafts/collective_cross_gram_packing_and_response_modulus_independent_audit.md`](drafts/collective_cross_gram_packing_and_response_modulus_independent_audit.md),
[`experiments/verify_cross_gram_response_metric_entropy.py`](experiments/verify_cross_gram_response_metric_entropy.py),
[`experiments/verify_collective_cross_gram_packing.py`](experiments/verify_collective_cross_gram_packing.py),
and
[`experiments/verify_collective_cross_gram_packing_independent_audit.py`](experiments/verify_collective_cross_gram_packing_independent_audit.py).

### Theorem 21.43 (PSD gluing has a contraction fibre of bounded fixed-scale cost)

Let `K_i=Y_iY_i^T\succeq0`, where each `Y_i` has full column rank, and
consider a two-piece joined sector

```math
K=\begin{pmatrix}K_1&C\\C^T&K_2\end{pmatrix}.       \tag{21.254}
```

Then `K\succeq0` if and only if

```math
\boxed{C=Y_1WY_2^T,\qquad ||W||_{op}\le1.}          \tag{21.255}
```

For fixed factor frames `W` is unique.  Under `Y_i\mapsto Y_iO_i`, it
transforms as `W\mapsto O_1^TWO_2`; the presentation, not `W` alone, is
quotiented by this gauge.  For `s` pieces the associative compatibility
object is one block correlation operator

```math
\Omega\succeq0,qquad \Omega_{ii}=I,qquad
K_{ij}=Y_i\Omega_{ij}Y_j^T.                         \tag{21.256}
```

Pairwise contraction tests alone are insufficient when `s>=3`.

The exact fibre admits a uniform approximate compression.  Suppose
`tr K_i<=p_i`, put `p=p_1+p_2`, and retain in each marginal only eigenvalues
larger than `tau p_i`.  If `K^h` is the resulting joined compression, then

```math
q_p(K,K^h):={1\over p^2}\max_{\epsilon\in\{+-1\}^p}
 |\epsilon^T(K-K^h)\epsilon|
\le\sqrt\tau+{\tau\over2}.                          \tag{21.257}
```

For both spectral sectors this gives `d_q<=2sqrt(tau)+tau`.  The square-root
loss is sharp for this marginal truncation architecture.  At target joined
error `eta`, taking `tau=Theta(eta^2)` leaves rank `O(eta^(-2))` per
marginal sector.  Approximating its factors in Frobenius norm by
`delta sqrt(p_i)` and its contraction in operator norm by `zeta` adds

```math
d_q\le4\delta+\zeta.                                \tag{21.258}
```

In particular the compatibility data add only `O_eta(1)` parameters once
the marginal factor states are fixed.

More intrinsically, fix arbitrary Gram--Rayleigh marginals on `p_1,p_2`
ports and metrize their two-sector compatibility fibre by the restriction
of `d_q`.  There is an absolute constant `C_0` such that

```math
\log\operatorname{Cov}_\eta(\mathfrak F,d_q)
\le C_0\eta^{-4}\log(C_0/\eta),                    \tag{21.259}
```

independently of `p_1+p_2`.  This is not a one-state theorem: for all small
`eta` there are Boolean-realizable fixed marginals with a compatibility
packing of size

```math
\exp\{c_0\eta^{-1}\log(1/\eta)\}.                  \tag{21.260}
```

#### Proof

Positivity of (21.254), tested on `(x,ty)` and optimized in `t`, gives

```math
|x^TCy|^2\le(x^TK_1x)(y^TK_2y).                    \tag{21.261}
```

This defines a contractive bilinear form on the two factor support spaces,
which is represented uniquely by `W`; conversely
`[[I,W],[W^T,I]]` is PSD exactly when `||W||op<=1`.  A common Gram
realization proves (21.256).

For the truncation, write `P_i,Q_i` for the high and low projections.  The
diagonal tails contribute at most `tau(p_1^2+p_2^2)`.  Contractivity in
(21.255) bounds each cross-tail bilinear form by
`sqrt(tau)p_1p_2`; remembering both symmetric cross blocks gives

```math
q_p(K,K^h)
\le{\tau(p_1^2+p_2^2)+4\sqrt\tau p_1p_2\over p^2}
\le\sqrt\tau+{\tau\over2}.                          \tag{21.262}
```

The rank bound follows from trace.  Expanding
`Y_1WY_2^T-Z_1VZ_2^T` proves (21.258), and Euclidean nets of the fixed-rank
factor and contraction balls give the stated parameter bound.

For (21.259), truncate each fixed marginal at eigenvalue
`Theta(eta^2(p_1+p_2))`; the cross-tail error is `O(eta)`, the retained
support ranks are `O(eta^(-2))`, and a Frobenius net of their contraction
balls proves the bound.  For (21.260), take `r=Theta(eta^(-1))` orthogonal
Boolean top eigenvectors, duplicate each on both shores, and vary the
relative frame by a permutation.  A relative-Hamming-distance permutation
code has `exp(Omega(r log r))` words, while

```math
||P_\pi-P_\sigma||_{\infty\to1}
\ge {4\over3}d_H(\pi,\sigma),                       \tag{21.263}
```

which gives fixed `eta` separation after normalization. `square`

Unordered interchange of the two spectral sectors supplies the relative
antipode bits in Theorem 21.41, but its cycle holonomy exposes a sharp
boundary of the PSD theory.  Rank-one edge contractions have the same local
switching action, yet one global PSD Gram realization forces
`C_(ij)=q_iq_j` and hence positive product around every cycle.  A negative
orientation-cycle product is coefficient-side interaction information, not
a PSD compatibility coordinate.  Thus sector ordering embeds in the Gram
carrier, whereas nontrivial coefficient holonomy must be retained by a
separate dynamic state.

The theorem controls the spherical carrier only at bounded total port mass
`mp/r=O(1)`.  It does not solve the exact Boolean old-spin integrality gap,
and at anti-pin scaling `m=r` the required metric accuracy deteriorates with
`p`.  Full proofs, checks, and an independent audit are in
[`drafts/psd_compatibility_fibre.md`](drafts/psd_compatibility_fibre.md),
[`drafts/psd_compatibility_fibre_independent_audit.md`](drafts/psd_compatibility_fibre_independent_audit.md),
[`drafts/psd_gluing_compatibility_entropy.md`](drafts/psd_gluing_compatibility_entropy.md),
[`experiments/verify_psd_compatibility_fibre.py`](experiments/verify_psd_compatibility_fibre.py),
[`experiments/verify_psd_compatibility_fibre_independent_audit.py`](experiments/verify_psd_compatibility_fibre_independent_audit.py),
and
[`experiments/verify_psd_gluing_compatibility_entropy.py`](experiments/verify_psd_gluing_compatibility_entropy.py).

### Theorem 21.44 (bounded port mass does not close the Boolean trust gap)

Let `H` be a symmetric entrywise sign matrix satisfying

```math
H^2=r^2I,\qquad \operatorname{tr}H=0,               \tag{21.264}
```

and suppose it has orthogonal Boolean top eigenvectors `a,b`.  Put
`A=H-diag(H)` and append two disjoint `m`-spin shores with repeated old--new
columns `a` and `b`.  Before filling the auxiliary block, the exact Boolean
and spherical trust responses are respectively

```math
\mathcal B_m=\max_{x\in\{+-1\}^n}
 \left\{|H_A(x)|+m|a^Tx|+m|b^Tx|\right\},           \tag{21.265}
```

```math
\mathcal S_m=\max_{||u||_2^2=n}
 \left\{\left|{1\over2}u^THu\right|
          +m|a^Tu|+m|b^Tu|\right\}.                \tag{21.266}
```

For every integer `m>=0`, these values are exactly

```math
\boxed{\mathcal B_m={rn\over2}+mn,\qquad
       \mathcal S_m={rn\over2}+\sqrt2mn.}           \tag{21.267}
```

Hence, with `p=2` and `m=r/2`, the total port mass is `mp/r=1` but

```math
{\mathcal S_m-\mathcal B_m\over rn}
={\sqrt2-1\over2}.                                  \tag{21.268}
```

This is scalable through tensor powers of the regular order-16 Walsh
matrix.  Their total parent order is `N=n+r`, so the same limiting gap holds
in `N^(3/2)` units.  Completing all auxiliary pairs by any exact signing
changes the difference by at most `2Q(C)=O(r^2)=O(n)` and cannot remove the
leading gap.

#### Proof

Orthogonal Boolean vectors agree on half the coordinates and disagree on
half, so for every endpoint signs

```math
||\epsilon_1a+\epsilon_2b||_1=n,qquad
||\epsilon_1a+\epsilon_2b||_2=\sqrt{2n}.            \tag{21.269}
```

Cube duality bounds the Boolean field by `mn`; the spectral bound gives
`rn/2`, and `x=a` attains both.  Euclidean duality bounds the spherical
field by `sqrt(2)mn`; `u=(a+b)/sqrt2` lies in the top eigenspace and attains
both.  This proves (21.267).  Tensoring the two orthogonal Boolean top poles
of the order-16 example proves scalability.  Finally the auxiliary energy
is uniformly bounded by `Q(C)` in each optimization, proving the completion
claim. `square`

This theorem falsifies uniform rounding of the **spherical trust value** at
bounded total mass, even in rank two.  It does not prove that `(G,R)` is
information-insufficient for every conceivable Boolean-response functional;
that stronger claim would require an equal-Gram collision with separated
Boolean responses.  Full details and independent checks are in
[`drafts/regular_hadamard_boolean_spherical_gap.md`](drafts/regular_hadamard_boolean_spherical_gap.md),
[`drafts/regular_hadamard_boolean_spherical_gap_independent_audit.md`](drafts/regular_hadamard_boolean_spherical_gap_independent_audit.md),
[`experiments/verify_regular_hadamard_boolean_spherical_gap.py`](experiments/verify_regular_hadamard_boolean_spherical_gap.py),
and
[`experiments/verify_regular_hadamard_boolean_spherical_gap_independent_audit.py`](experiments/verify_regular_hadamard_boolean_spherical_gap_independent_audit.py).

### Theorem 21.45 (one exposed flat optimizer suffices for Boolean recovery)

Let `H` be symmetric with `H^2=r^2I`, let `w_1,...,w_p` be Boolean ports
of common width `m`, and maximize the channel family

```math
F_{\sigma,\epsilon}(u)
={\sigma\over2}u^THu
 +m\left(\sum_{i=1}^p\epsilon_iw_i\right)^Tu.       \tag{21.270}
```

Write `mathcal S` and `mathcal B` for its maxima on the sphere
`||u||^2=n` and the Boolean cube.  If one sphere point in one channel obeys

```math
F_{\sigma,\epsilon}(u)\ge\mathcal S-\xi rn,
\qquad
\phi(u):=1-{||u||_1\over n}\le\varphi,              \tag{21.271}
```

then, with `c=mp/r`,

```math
\boxed{
0\le\mathcal S-\mathcal B
\le rn\left[\xi+(1+c)\sqrt{2\varphi}\right].}      \tag{21.272}
```

The condition must concern a near-global exposed channel; locating that
channel is a separate response-state problem.  It is nevertheless a strict
recovery certificate: after exposure, it stores one `l_1` statistic and
uses coordinatewise sign rounding, not a Boolean response table.

The stronger demand that every direction in a `d`-dimensional trust span be
uniformly flat is impossible for `d>=2`.  Every such subspace contains
`u`, `||u||^2=n`, with

```math
{||u||_1\over n}
\le\gamma_d
:=\sqrt d\,{\Gamma(d/2)\over
                 \sqrt\pi\Gamma((d+1)/2)}<1.        \tag{21.273}
```

Yet exposed flatness occurs in a nontrivial scalable family.  For every
even `d`, put `q=2^d`, `n=q^2`, and regularize the Walsh matrix on
`F_2^d\times F_2^d` by the self-dual bent vector
`y_0(x,z)=(-1)^(x\cdot z)`.  A self-dual linear subspace `L` of size `q`
produces two distinct Boolean top poles

```math
a=\mathbf1,
\qquad b=\mathbf1-2\mathbf1_L,
\qquad {a^Tb\over n}=1-{2\over q}.                  \tag{21.274}
```

For two ports of width `m=q/2`, total port mass is one and

```math
{\mathcal S-\mathcal B\over rn}
=\sqrt{1-{1\over q}}-\left(1-{1\over q}\right)
=O(n^{-1/2}).                                       \tag{21.275}
```

Thus rank two can exhibit either the fixed gap of Theorem 21.44 or
asymptotically exact Boolean recovery; the distinguishing mechanism is
exposed pole synchronization, not rank.

#### Proof

For `x=sgn(u)`,

```math
{||x-u||_2^2\over n}=2\phi(u).                      \tag{21.276}
```

The quadratic rounding loss is at most
`rn sqrt(2phi)`, and the port field costs at most
`mpn sqrt(2phi)=crn sqrt(2phi)`, proving (21.272).

For (21.273), choose a uniform direction `theta` on the unit sphere of the
subspace.  Rotational invariance and Cauchy--Schwarz give

```math
\mathbb E||\sqrt nQ\theta||_1\le n\gamma_d,         \tag{21.277}
```

and strictness follows from `E|theta_1|<sqrt(E theta_1^2)`.

For the construction, take the graph
`L={(x,Mx)}` of a fixed-point-free symmetric coordinate pairing.  Then
`L=L^perp`, `y_0|_L=1`, and Walsh inversion gives
`W1_L=q1_L`.  Hence `y_0-2 1_L` is another Boolean top eigenvector.
For two Boolean top poles of correlation `rho>=0`, cube and sphere duality
give exactly

```math
\mathcal B={rn\over2}+mn(1+\rho),\qquad
\mathcal S={rn\over2}+mn\sqrt{2(1+\rho)}.           \tag{21.278}
```

Substitution proves (21.275). `square`

An arbitrary exact-sign completion on the `q` auxiliary vertices changes
the gap by only `O(q^2)=O(n)`, so the little-oh recovery survives, although
the coefficient in (21.275) need not.  Full proof and audits are in
[`drafts/exposed_boolean_synchronization.md`](drafts/exposed_boolean_synchronization.md),
[`drafts/exposed_boolean_synchronization_independent_audit.md`](drafts/exposed_boolean_synchronization_independent_audit.md),
[`experiments/verify_exposed_boolean_synchronization.py`](experiments/verify_exposed_boolean_synchronization.py),
and
[`experiments/verify_exposed_boolean_synchronization_independent_audit.py`](experiments/verify_exposed_boolean_synchronization_independent_audit.py).

### Theorem 21.46 (common-pole synchronization is a multiplicative certificate)

In the setting of Theorem 21.45, suppose `x_0` is a Boolean top pole,
`Hx_0=rx_0`, and define

```math
\delta(x_0;W)
=1-{1\over pn}\sum_{i=1}^p|w_i^Tx_0|.              \tag{21.279}
```

Then

```math
\boxed{0\le\mathcal S-\mathcal B
       \le {mp\over r}\delta(x_0;W)rn.}             \tag{21.280}
```

An arbitrary exact-sign completion `C` on the `pm` auxiliary spins adds at
most `2Q(C)` to this gap.  Thus at bounded total mass, `delta=o(1)` is a
strict one-witness Boolean recovery certificate.

It has an exact composition algebra.  For two tensor systems, use the
product pole and the Cartesian multiset of tensor-product ports.  If their
deficits are `delta_1,delta_2`, then

```math
\boxed{1-\delta_{12}=(1-\delta_1)(1-\delta_2).}     \tag{21.281}
```

After `L` factors,

```math
1-\delta_{[L]}=\prod_{j=1}^L(1-\delta_j),
\qquad
\delta_{[L]}\le\sum_j\delta_j.                     \tag{21.282}
```

The actual recovery condition is
`c_[L]delta_[L]=o(1)`, including the composed total port mass; deficit alone
is sufficient only when `c_[L]=O(1)`.

#### Proof

The sphere is bounded by `rn/2+mpn`.  In the Boolean response, take
`x=x_0`, the top quadratic channel, and align every endpoint independently;
the value is at least

```math
{rn\over2}+m\sum_i|w_i^Tx_0|
={rn\over2}+mpn(1-\delta),                          \tag{21.283}
```

which proves (21.280).  Tensor inner products multiply, so their normalized
average absolute correlations multiply, proving (21.281)--(21.282).
`square`

The scalar deficit is an exact compositional **certificate observable**, not
a complete contextual state: explicit order-four families have equal
`delta=1/2` and Boolean responses `8` and `10`.  Full proof and independent
audit are in
[`drafts/common_pole_synchronization_algebra.md`](drafts/common_pole_synchronization_algebra.md),
[`experiments/verify_common_pole_synchronization.py`](experiments/verify_common_pole_synchronization.py),
[`drafts/exposed_flatness_common_pole_independent_audit.md`](drafts/exposed_flatness_common_pole_independent_audit.md),
and
[`experiments/verify_flatness_common_pole_independent_audit.py`](experiments/verify_flatness_common_pole_independent_audit.py).

### Theorem 21.47 (active-eigenspace recovery has sphere-covering complexity)

Let `V` be a `d`-dimensional subspace of the `+r` eigenspace of a symmetric
`H` with `||H||<=r`, and let `P=P_V`.  For Boolean `x`, put

```math
q(x)={x^THx\over rn},
\qquad u(x)={Px\over\sqrt n}.                       \tag{21.284}
```

For a unit `v in V` and field strength `beta>0`, the spherical response is
`S_beta(v)=1/2+beta`, while

```math
B_\beta(v)=\max_x\{q(x)/2+\beta\langle v,u(x)\rangle\}. \tag{21.285}
```

Let `Gamma_beta=sup_v(S_beta(v)-B_beta(v))` and
`W_delta={x:1-q(x)<=delta}`.  If `Gamma_beta<=eta`, then every direction has
a witness in `W_(2eta)` with support at least `1-eta/beta`; consequently,
for `d>=2` and `0<eta/beta<=1/2`,

```math
\boxed{|W_{2\eta}|
\ge\left({\beta\over2\eta}\right)^{(d-1)/2}.}      \tag{21.286}
```

Conversely, if a finite library `L subset W_delta` has support at least
`1-epsilon` in every unit direction, then

```math
\Gamma_\beta\le\delta/2+\beta\epsilon.             \tag{21.287}
```

These orders match.  If pointwise witnesses achieve support deficit
`epsilon<=1/4`, one can retain only

```math
|L|\le(C/\sqrt\epsilon)^{d-1}                       \tag{21.288}
```

of them and obtain support deficit `4epsilon` uniformly.  With a spectral
gap `gamma` below `V`, a defect-`delta` angular cover of deficit `alpha`
certifies

```math
\Gamma_\beta
\le\delta/2+\beta(\alpha+\delta/\gamma).            \tag{21.289}
```

If `V` is the complete strict top eigenspace, it contains at most `2^d`
Boolean eigenvectors.  Hence exact-eigen witnesses alone cannot achieve
uniform support deficit below

```math
2^{-1-2d/(d-1)}\ge1/32.                             \tag{21.290}
```

#### Proof

For a Boolean optimizer the gap decomposes exactly into two nonnegative
terms,

```math
S_\beta(v)-B_\beta(v)
={1-q(x)\over2}+\beta(1-\langle v,u(x)\rangle).     \tag{21.291}
```

This proves necessity and sufficiency.  Unit directions served with support
`1-epsilon` form spherical caps of angular radius
`arccos(1-epsilon)`; their surface mass is at most
`(2epsilon)^((d-1)/2)`, proving (21.286).  A
`sqrt(epsilon)` sphere net and one pointwise witness per net point give
(21.288).  The spectral gap yields
`1-||u(x)||^2<=delta/gamma`, proving (21.289).
Finally choose `d` coordinate functionals injective on `V`; they inject
`V intersect {+-1}^n` into `{+-1}^d`, and combine with the cap bound to get
(21.290). `square`

This is a continuum fixed-strength query theorem.  A finite `p`-port trust
table exposes at most `2^p` directions with direction-dependent strengths,
and the two quadratic eigenspaces must be treated separately.  Thus
(21.286) prices explicit witness libraries; it is not an unconditional bit
lower bound for every algebraically generated carrier.  Full proof and audit
are in
[`drafts/boolean_active_eigenspace_synchronization.md`](drafts/boolean_active_eigenspace_synchronization.md),
[`drafts/boolean_active_eigenspace_synchronization_independent_audit.md`](drafts/boolean_active_eigenspace_synchronization_independent_audit.md),
[`experiments/verify_boolean_active_eigenspace_synchronization.py`](experiments/verify_boolean_active_eigenspace_synchronization.py),
and
[`experiments/verify_boolean_active_eigenspace_synchronization_independent_audit.py`](experiments/verify_boolean_active_eigenspace_synchronization_independent_audit.py).

### Theorem 21.48 (exposed flatness has an exact hierarchical composition law)

Partition a sphere vector `u in R^N` into blocks of relative sizes
`lambda_i=n_i/N`.  Put

```math
\rho_i={||u_i||_2\over\sqrt{n_i}},
\qquad
\phi_i=1-{||u_i/\rho_i||_1\over n_i},               \tag{21.292}
```

with zero contribution for a zero block.  Since
`sum_i lambda_i rho_i^2=1`, its global flatness satisfies the exact chain
rule

```math
\boxed{
\phi(u)
={1\over2}\sum_i\lambda_i(\rho_i-1)^2
 +\sum_i\lambda_i\rho_i\phi_i.}                   \tag{21.293}
```

On a partition tree, if `A_v` is the local amplitude-allocation defect and
`omega_v=(N_v/N)R_v` is its transported RMS weight, then

```math
\phi(u)=\sum_{v\ {m internal}}\omega_vA_v
       +\sum_{\ell\ {m leaf}}\omega_\ell\phi_\ell,
\qquad
\sum_{v\ {m at\ depth}\ d}\omega_v\le1.         \tag{21.294}
```

Consequently, if every node at depth `j` has `A_v<=a_j` and every leaf has
flatness at most `b`, then a spherical vector exposed within
`xi Lambda N` for `Lambda>0` and

```math
F(u)={1\over2}u^TMu+h^Tu,
\qquad ||M||\le\Lambda,
\qquad \kappa={||h||\over\Lambda\sqrt N},           \tag{21.295}
```

obeys

```math
S-B\le\Lambda N\left[
\xi+(1+\kappa)\sqrt{2\left(b+\sum_j a_j\right)}
\right].                                            \tag{21.296}
```

The summability condition is real.  On equal binary splits take relative
RMS amplitudes `rho_+-=sqrt(1+-delta)` and
`s=(rho_++rho_-)/2<1`.  Every scalar leaf is Boolean, but after depth `D`,

```math
{||u_D||_1\over N}=s^D,
\qquad \phi(u_D)=1-s^D.                             \tag{21.297}
```

For the rank-one linear landscape `F_D(x)=u_D^Tx`, the sphere and cube
maxima are exactly `N` and `Ns^D`; the normalized gap tends to one.

#### Proof

Additivity of `l_1` gives
`phi=1-sum lambda_i rho_i+sum lambda_i rho_i phi_i`, while

```math
1-\sum_i\lambda_i\rho_i
={1\over2}\sum_i\lambda_i(\rho_i-1)^2,             \tag{21.298}
```

proving (21.293).  Iteration telescopes the weights in (21.294), and
Cauchy--Schwarz bounds their sum at every depth.  Apply sign rounding and
(21.276) to prove (21.296).  Multiplication of the binary amplitudes gives
`||u_D||_1/N=s^D`, proving the pumpable example. `square`

The additive pair `(||u||_2^2,||u||_1)` is an exact state for a **supplied
exposed witness**, not an optimizer-selection state for a switching response
roof.  The theorem nevertheless identifies a new dynamic obstruction:
perfect local Booleanity can be destroyed solely by repeated amplitude
allocation.  Full proof and audit are in
[`drafts/exposed_flatness_composition_law.md`](drafts/exposed_flatness_composition_law.md),
[`experiments/verify_exposed_flatness_composition.py`](experiments/verify_exposed_flatness_composition.py),
[`drafts/exposed_flatness_common_pole_independent_audit.md`](drafts/exposed_flatness_common_pole_independent_audit.md),
and
[`experiments/verify_flatness_common_pole_independent_audit.py`](experiments/verify_flatness_common_pole_independent_audit.py).

### Theorem 21.49 (generic delocalization does not imply Boolean flatness)

Two independent benchmarks rule out a generic exposed-flatness theorem.

First, specialize Bourgade--Yau eigenvector universality to a real symmetric
Rademacher generalized Wigner matrix with independent entries **including
the diagonal**, normalized by `1/sqrt(N)`.  For every eigenvector index in
their stated edge/central-bulk set, if `v_k` has norm one and
`u_k=sqrt(N)v_k`, then

```math
{||u_k||_1\over N}\longrightarrow\sqrt{2/\pi}
\quad\hbox{in }L^2.                                 \tag{21.299}
```

Thus its flatness deficit tends to
`1-sqrt(2/pi)`, not zero.  No hollow-matrix extension is asserted.

Second, for `n=2m` the explicit hollow signing

```math
A_m=\begin{pmatrix}
J_m-I_m&J_m\\J_m&-(J_m-I_m)
\end{pmatrix}                                       \tag{21.300}
```

has a simple top eigenvalue

```math
r_m=\sqrt{(m-1)^2+m^2}.                             \tag{21.301}
```

Its norm-`sqrt(n)` top eigenvector is constant on the two blocks with ratio

```math
t_m={r_m-(m-1)\over m},                             \tag{21.302}
```

and satisfies

```math
{||u_m||_1\over n}
={1+t_m\over\sqrt{2(1+t_m^2)}}
\longrightarrow\cos(\pi/8)<1.                     \tag{21.303}
```

#### Proof

For the random model, Bourgade--Yau gives uniform one- and two-coordinate
Gaussian moment convergence.  Averaging those coordinates and using uniform
integrability yields convergence of the first and second moments of the
normalized `l_1` norm to `sqrt(2/pi)` and `2/pi`, proving (21.299).

For (21.300), the block-constant quotient is
`[[m-1,m],[m,-(m-1)]]`; the orthogonal complements have eigenvalues `-1`
and `+1`.  Diagonalizing the quotient gives (21.301)--(21.302), and direct
normalization gives (21.303). `square`

The deterministic family has quadratic cap of order `n^2`, so it does not
falsify a near-minimizer-specific rigidity theorem.  It does prove that
ordinary delocalization, isotropy, Kashin/Dvoretzky-type Euclidean sections,
vector balancing, and inverse Littlewood--Offord structure cannot by
themselves supply `phi=o(1)`: Gaussian coordinates sit at the wrong equality
constant.  Any positive route must force near-cube equality, not merely
spread mass.  Exact checks and the source audit are in
[`drafts/hypercube_flatness_literature_scout.md`](drafts/hypercube_flatness_literature_scout.md)
and
[`experiments/verify_hypercube_flatness_scout.py`](experiments/verify_hypercube_flatness_scout.py).

### Theorem 21.50 (Boolean ports have an exact Fourier algebra and linear fixed-scale entropy)

For `W in {+-1}^(n times p)`, define its labelled endpoint response

```math
L_W(epsilon)=||Wepsilon||_1
```

and let `mu_W` be the histogram of its rows in the projective group
`G_p={+-1}^p/{+-1}`.  With `g_p(z)=|sum_i z_i|`, one has the convolution

```math
L_W(epsilon)=sum_(s in G_p)mu_W(s)g_p(sepsilon).   \tag{21.304}
```

Characters of `G_p` are the even subsets `S subseteq[p]`.  If
`P=2ceil(p/2)` and `|S|=2k>=2`, then

```math
hat g_p(emptyset)={(P-1)!!\over(P-2)!!},
\qquad
hat g_p(S)=(-1)^(k-1)
{(2k-3)!!(P-2k-1)!!\over(P-2)!!}.                 \tag{21.305}
```

Every multiplier is nonzero.  Hence two systems have the same complete
labelled response table if and only if they have the same projective row
histogram.  The exact quotient on `n` rows has

```math
{n+2^(p-1)-1\choose2^(p-1)-1}                    \tag{21.306}
```

states, and it composes by histogram addition.  The same state answers all
real weighted queries `||Wa||_1`, so it is the exact zonotope-support
carrier, not merely a Boolean-endpoint code.

Its fixed-scale response image is far smaller.  For probability measures
on `G_p`, set

```math
d_p(mu,nu)=max_epsilon
\left|E_mu{|s dot epsilon|\over p}
      -E_nu{|s dot epsilon|\over p}\right|.        \tag{21.307}
```

Every `mu` has an equally weighted `eta`-coreset on at most

```math
k=ceil(16/eta^2)                                   \tag{21.308}
```

row types.  Therefore

```math
log Cov_p(eta)<=k(p-1)log2.                        \tag{21.309}
```

Conversely, point masses satisfy the exact identity

```math
d_p(delta_s,delta_t)={2d_H^proj(s,t)\over p},      \tag{21.310}
```

so projective Hamming codes give

```math
\boxed{log Cov_p(eta)=Theta_eta(p)}                \tag{21.311}
```

for every fixed `0<eta<1/2`.

#### Proof

Fourier transformation of (21.304) multiplies `hat mu_W(S)` by the
nonzero coefficient (21.305), proving invertibility and the weak-composition
count.  Formula (21.305) follows by adjoining an unused sign when `p` is odd
and iterating the even-order binomial recurrence

```math
{hat g_P(2k)\over hat g_(P-2)(2k)}
={P-2k-1\over P-2}.                                \tag{21.312}
```

For the coreset, sample `S_1,...,S_k` from `mu`.  Symmetrization and scalar
Rademacher contraction give

```math
E d_p(mu,mu_k)
<= {4\over kp}E_sigma sum_(i=1)^p
       \left|sum_(j=1)^ksigma_j(S_j)_i\right|
<= {4\over sqrt k}.                               \tag{21.313}
```

Some sample therefore proves (21.308)--(21.309).  Reverse triangle
inequality gives the upper side of (21.310), while querying `epsilon=s`
gives equality; projective Hamming packing completes the proof. `square`

Gram data are exactly the constant and degree-two Fourier truncation.  They
therefore determine the full pure-linear table for `p<=3`; at `p=4`, the
uniform row distribution and doubled even-parity distribution have equal
Gram matrices but joint supports `3n/2` and `2n`.  Exact feature dimension is
exponential in `p`, whereas fixed-error response information is linear.
The coreset is static: resparsifying after each merge can accumulate fresh
error and is not an all-depth congruence.  Full proofs and audits are in
[`drafts/boolean_port_fourier_feature_algebra.md`](drafts/boolean_port_fourier_feature_algebra.md),
[`drafts/boolean_port_fourier_feature_algebra_independent_audit.md`](drafts/boolean_port_fourier_feature_algebra_independent_audit.md),
[`drafts/boolean_port_dimension_free_coreset.md`](drafts/boolean_port_dimension_free_coreset.md),
[`drafts/four_port_gram_boolean_collision.md`](drafts/four_port_gram_boolean_collision.md),
[`experiments/verify_boolean_port_fourier_feature_algebra.py`](experiments/verify_boolean_port_fourier_feature_algebra.py),
[`experiments/verify_boolean_port_dimension_free_coreset.py`](experiments/verify_boolean_port_dimension_free_coreset.py),
and
[`experiments/verify_four_port_gram_boolean_collision.py`](experiments/verify_four_port_gram_boolean_collision.py).

### Theorem 21.51 (equal Gram--Rayleigh states can hide a leading Boolean response)

There is a regular symmetric Hadamard matrix `H_0` of order 16 and two
four-port tuples `W^A,W^B` of Boolean `+4` eigenvectors such that

```math
G^A=G^B=R^A=R^B=
\begin{pmatrix}
1&1/2&0&0\\1/2&1&0&0\\0&0&1&0\\0&0&0&1
\end{pmatrix},                                    \tag{21.314}
```

but their joint Boolean field supports are respectively `32` and `28`.
This finite collision amplifies without diluting the scale.  Put

```math
H_j=H_0^(tensor j),\quad n_j=16^j,\quad r_j=4^j,
\qquad W_j=W_0\otimes1_(16^(j-1)),                \tag{21.315}
```

and give every port width `m_j=r_j`.  Then the two states remain exactly
equal, total port mass is `pm_j/r_j=4`, and

```math
\boxed{
\mathcal B_(r_j)(H_j;W_j^A)
-\mathcal B_(r_j)(H_j;W_j^B)
\ge {r_jn_j\over8}={n_j^(3/2)\over8}.}            \tag{21.316}
```

Filling all pairs among the `4r_j` auxiliary vertices by the same arbitrary
exact signing changes the difference by at most `O(r_j^2)=O(n_j)`.  Hence
no decoder receiving only `(G,R,n,r,m,p)`, even together with the common
`H_j` and public completion, can have uniform Boolean error `o(rn)` on this
class.

#### Proof

At the seed, endpoint enumeration gives supports `32` and `28`.  One Boolean
word on the high side has child energy `24` and field support `32`.  Tensor
it with the common all-one pole.  Writing `N=16^(j-1)`, this gives the lower
bound `38r_jN`; the spectral child bound plus the low-side support gives the
upper bound `36r_jN`.  Their difference is `2r_jN=r_jn_j/8`.  Gram and
Rayleigh products scale identically under the tensor lift, and cap
Lipschitzness proves the completion statement. `square`

There is also a scalable one-port collision outside the involutive class.
Two order-16 zero-Rayleigh ports have responses `64` and `78`.  Tensoring
the common child with `J_k`, the port with `1_k`, and the width by `k` gives

```math
\mathcal B_(km)(H\otimes J_k;w\otimes1_k)
=k^2\max_(y\in\Gamma_k^n)
 \left\{|y^THy|/2+m|w^Ty|\right\},                \tag{21.317}
```

and independent rounding bounds the box maximum by the seed cube maximum
plus `sum_i|H_ii|/2`.  The two equal singleton `(G,R)` states consequently
remain separated by at least `(3/32)rn`.  A strict regular-Hadamard common
factor does not admit this magnetization quotient: an orthogonal top mode
with zero fibre magnetization can retain child energy `rn/2`.

Theorem 21.51 upgrades the earlier spherical integrality gap to an
information collision.  Pairwise geometry is not a complete Boolean state;
the first scalable witness is a fourth-order row-pattern channel.  Full
proofs and independent exact derivations are in
[`drafts/regular_hadamard_equal_gram_rayleigh_collision.md`](drafts/regular_hadamard_equal_gram_rayleigh_collision.md),
[`drafts/three_port_gram_closure.md`](drafts/three_port_gram_closure.md),
[`drafts/rank_one_common_factor_amplification.md`](drafts/rank_one_common_factor_amplification.md),
[`experiments/verify_regular_hadamard_equal_gram_rayleigh_collision.py`](experiments/verify_regular_hadamard_equal_gram_rayleigh_collision.py),
[`experiments/verify_three_port_gram_closure.py`](experiments/verify_three_port_gram_closure.py),
and
[`experiments/verify_rank_one_common_factor_amplification.py`](experiments/verify_rank_one_common_factor_amplification.py).

### Theorem 21.52 (odd product closure gives a growing exact Boolean carrier)

Let `||H||op<=r`, let `w_1,...,w_p` be Boolean `+r` eigenvectors, and fix
an antipodally odd tie-broken majority selector `tau`.  If every coordinate
product

```math
w_S=\bigodot_(i\in S)w_i
```

appearing with nonzero coefficient in the Boolean Fourier expansion of
`tau` is also a `+r` eigenvector, then every labelled trust channel is exact:

```math
\boxed{
B_epsilon={rn\over2}+m||sum_i epsilon_iw_i||_1,
\qquad
S_epsilon={rn\over2}+m\sqrt n||sum_i epsilon_iw_i||_2.}
                                                               \tag{21.318}
```

The witness is the coordinatewise selector itself.  Its Fourier expansion
is a linear combination of positive top poles, so it pays the quadratic and
all port channels jointly before absolute values.  Consequently the
projective histogram of Theorem 21.50 is the minimal exact labelled Boolean
carrier on this closed class.

The selector-independent full condition has an intrinsic form.  Fixing
`w_1` and

```math
\mathcal C=<w_1\odot w_i:2<=i<=p>,                \tag{21.319}
```

all odd port products form exactly the affine multiplicative coset
`w_1 mathcal C`.  Thus full closure means that this coset lies in the
Boolean positive top eigenset.  It is preserved by block concatenation and
tensor product; on tensor ports the exact state law is projective-group
convolution

```math
\mu_(W\boxtimes V)=\mu_W*_(G_p)\mu_V.             \tag{21.320}
```

This mechanism has a dense growing-arity realization.  The product-closed
triple in the order-16 regular Hadamard seed supplies an affine coset of four
positive poles.  Tensor `j` seeds and take one base pole plus the `2j`
factor generators as ports.  Then

```math
n_j=16^j,\quad r_j=4^j,\quad
p_j=2j+1={1\over2}\log_2n_j+1,                    \tag{21.321}
```

every odd product is a positive pole, and the exact histogram carrier has

```math
2^(p_j-1)=\sqrt{n_j}quad\hbox{bins},\qquad
O(\sqrt{n_j}\log n_j)\quad\hbox{bits}.            \tag{21.322}
```

At width `m_j=floor(r_j/p_j)`, total port mass is at most one.  Arbitrary
public exact-sign completion on the at most `r_j` new vertices costs only
`O(n_j)=o(r_jn_j)`.  This is therefore a dense, growing-interface,
exact-sign benchmark with a strict sub-landscape composable state.

For `p<=2` closure is automatic for top ports.  At `p=3` it asks only that
the triple product be top and yields an explicit Gram formula.  The
four-port collision in Theorem 21.51 fails the product condition, proving
that it is genuine synchronization rather than disguised pairwise data.
Full proofs and independent audit are in
[`drafts/boolean_port_product_algebra_closure.md`](drafts/boolean_port_product_algebra_closure.md),
[`drafts/boolean_port_product_algebra_closure_independent_audit.md`](drafts/boolean_port_product_algebra_closure_independent_audit.md),
[`drafts/three_port_gram_closure.md`](drafts/three_port_gram_closure.md),
and
[`experiments/verify_boolean_port_product_algebra_closure.py`](experiments/verify_boolean_port_product_algebra_closure.py).

### Theorem 21.53 (Boolean port states have exact merge and tensor reuse laws)

The fixed-scale coreset in Theorem 21.50 has two genuine dynamic extensions.

First, give every labelled row occurrence `u` and each of `k` replicas an
independent public random priority.  An aggregate stores its cardinality and
the minimum-priority row type in each replica.  For disjoint union,
coordinatewise priority minimum and cardinality addition form an exact
associative carrier.  At every one of `T` aggregates declared independently
of the priorities, with probability at least `1-delta`, its normalized
response error is simultaneously at most

```math
{4\over\sqrt k}+\sqrt{{\log(T/delta)\over2k}}.     \tag{21.323}
```

Thus union-tree depth does not add approximation error; only the logarithm
of the simultaneously certified node family enters.  Adaptive subsets and
reconvergent overlap without relabelling are outside this theorem.

Second, normalized tensor histograms compose by probability convolution on
`G_p`.  If `R_mu(epsilon)=E_mu|s dot epsilon|/p`, then

```math
\boxed{R_(mu*lambda)=lambda*R_mu},                 \tag{21.324}
```

so

```math
d_p(mu*lambda,nu*lambda)<=d_p(mu,nu),
```

```math
d_p(mu_1*cdots*mu_L,nu_1*cdots*nu_L)
<=sum_(i=1)^Ld_p(mu_i,nu_i).                       \tag{21.325}
```

The coefficient `L` is locally sharp.  For a nontrivial order-two element
`a` and `mu_t=(1-t)delta_e+t delta_a`,

```math
mu_t^(*L)=(1-q_L)delta_e+q_Ldelta_a,
\qquad q_L={1-(1-2t)^L\over2},                    \tag{21.326}
```

and the ratio of the `L`-fold to one-step response error tends to `L` as
`t` tends to zero.

Independent leaf **occurrences** nevertheless admit a depth-independent
compiler: draw `k` independent types at every leaf occurrence and multiply
corresponding replicas up the declared tensor tree.  Every node then contains
`k` iid samples from its exact convolution law and obeys (21.323).  Reusing
the same sample bank in both inputs is not valid: coordinatewise squaring
returns the identity, whereas the uniform law on `{e,a}` is convolution
idempotent and has fixed positive response distance from that identity.

There is also an exact forgetting mechanism.  Every response difference has
zero uniform mean.  If `lambda>=alpha u_p` pointwise, where `u_p` is uniform
on `G_p`, then

```math
\boxed{d_p(mu*lambda,nu*lambda)
       <=(1-alpha)d_p(mu,nu).}                    \tag{21.327}
```

If local approximation errors are `eta_j` and successive contexts have
uniform masses `alpha_j`, the final error is bounded by

```math
e_0prod_(j=1)^L(1-alpha_j)
+sum_(i=1)^Leta_i prod_(j=i+1)^L(1-alpha_j).       \tag{21.328}
```

This theorem makes the static/dynamic law exact in one nontrivial algebra:
linear fixed-scale response information is reusable under declared
independence, accumulates sharply under uncontrolled semantic reuse, and is
forgotten geometrically under a Doeblin component.  Full proofs and checks
are in
[`drafts/boolean_port_mergeable_reservoir.md`](drafts/boolean_port_mergeable_reservoir.md),
[`drafts/boolean_port_mergeable_reservoir_audit.md`](drafts/boolean_port_mergeable_reservoir_audit.md),
[`drafts/boolean_port_convolution_reuse.md`](drafts/boolean_port_convolution_reuse.md),
[`drafts/boolean_port_convolution_reuse_audit.md`](drafts/boolean_port_convolution_reuse_audit.md),
[`drafts/boolean_port_convolution_reuse_independent_audit.md`](drafts/boolean_port_convolution_reuse_independent_audit.md),
[`experiments/verify_boolean_port_mergeable_reservoir.py`](experiments/verify_boolean_port_mergeable_reservoir.py),
and
[`experiments/verify_boolean_port_convolution_reuse.py`](experiments/verify_boolean_port_convolution_reuse.py).

### Theorem 21.54 (robust product synchronization and its dynamic boundary)

Let `tau:{+-1}^p->{+-1}` be an antipodally odd majority selector and let
`Z` contain the Boolean products in the nonzero Fourier support of `tau`.
For a symmetric child `H` with `||H||_op<=r`, set

```math
G={Z^TZ\over n},\qquad R={Z^THZ\over rn},\qquad D=G-R\succeq0. \tag{21.329}
```

If `a^epsilon_S=hat(tau)(S)prod_(i in S)epsilon_i`, define the intrinsic
joint defect

```math
Delta_tau=\max_epsilon (a^epsilon)^TD a^epsilon.  \tag{21.330}
```

Parseval gives `||a^epsilon||_2=1`, while `Za^epsilon` is the one Boolean
selector that simultaneously maximizes every port field.  Therefore the
labelled Boolean trust response obeys

```math
0\le {rn\over2}+m||Wepsilon||_1-B_epsilon
\le {rn\over2}(a^epsilon)^TD a^epsilon
\le {rn\over2}Delta_tau.                          \tag{21.331}
```

This pays all Fourier product channels jointly.  The operator certificate
`delta=||D||_op` is sufficient, but can be much stronger than the semantic
defect.  Individual product deficits `d_S=D_SS` instead give the explicit
but potentially expensive bound

```math
(a^epsilon)^TD a^epsilon
\le\left(\sum_S|hat(tau)(S)|\sqrt{d_S}\right)^2. \tag{21.332}
```

For corresponding-port tensor products, `G_12=G_1 circ G_2` and
`R_12=R_1 circ R_2`.  Row characters translate the selector orbit and the
Schur product theorem controls the remaining cross term, giving the
dimension-free semantic law

```math
\boxed{Delta_(tau,12)\le Delta_(tau,1)+Delta_(tau,2).}         \tag{21.333}
```

The operator certificate also satisfies

```math
||D_12||_op\le
\min\{delta_1+kappa_1delta_2,\ kappa_2delta_1+delta_2\},       \tag{21.334}
```

where `kappa_i=max_S||H_i z_(i,S)||/(r_i sqrt(n_i))<=1`.
The indefinite Rayleigh Schur multiplier is contractive because its entries
factor as `<u_S,T u_T>` with `||T||<=1`; entrywise boundedness alone would
not suffice.

If a projective histogram approximation has response error `eta` and
`c=mp/r`, the resulting decoder has uniform error

```math
rn(Delta_tau/2+c eta),                            \tag{21.335}
```

using `O(p/eta^2)` fixed-scale bits.  With `pm=O(r)=o(n)`, arbitrary public
exact-sign completion of the auxiliary shore adds only `O(r^2)=o(rn)`.

This state is deliberately an approximate certificate, not an exact tensor
quotient.  Already for three-port majority there are two admissible
contractions with the same projective histogram and the same **entire**
endpoint defect table, all equal to `1/2`, whose self-tensors have endpoint
defects `48/64` and `47/64`.  Off-query Schur coherence is invisible in one
step but becomes observable after composition.  Thus (21.333) is a reusable
upper-error law; exact updating requires richer defect geometry.

### Theorem 21.55 (relative synchronization survives growing Cartesian pole algebras)

The growing affine-coset construction uses a different tensor operation.
If `V_i` lists every pole in factor `i`, its all-tuples pole list has

```math
G=\bigotimes_iG_i,\qquad R=\bigotimes_iR_i.       \tag{21.336}
```

Define the Gram-relative defect on the represented pole span by

```math
delta_i=\sup_{c:c^TG_ic>0}{c^T(G_i-R_i)c\over c^TG_ic}.       \tag{21.337}
```

For symmetric contraction factors,

```math
\boxed{delta_(1\boxtimes\cdots\boxtimes L)
       \le\sum_i delta_i.}                       \tag{21.338}
```

If every compressed factor is positive semidefinite and
`0<=D_i<=delta_iG_i` with `delta_i<=1`, the sharper sharp law is

```math
\boxed{0\le D\le
\left(1-\prod_i(1-delta_i)\right)G.}             \tag{21.339}
```

No inverse of `G_i` is needed; the inequalities automatically descend to
its support.  This exactly covers the PC.3 hierarchy: its `4^j` odd product
poles are the Cartesian tensor list of the four seed affine-coset poles.
Consequently factorwise defects with `sum_i delta_i=o(1)` retain
`o(rn)` Boolean trust error even though the port arity grows and the raw
operator norm `||D||` can acquire a spurious factor `2^j`.  Exact PC.3 is
the case `delta_i=0`.

The theorem is not a generic growing-arity result.  It assumes the declared
pole span factors, and repeated fixed positive defect tends to a leading
loss.  For a factored construction its certificate is constant-size per
factor; without that presentation the final relative matrix may require
`Theta(q^2)` parameters.

Full proofs, independent audits, sharp counterexamples, and exact checks are
in
[`drafts/robust_boolean_product_synchronization.md`](drafts/robust_boolean_product_synchronization.md),
[`drafts/robust_product_synchronization_independent_audit.md`](drafts/robust_product_synchronization_independent_audit.md),
[`experiments/verify_robust_boolean_product_synchronization.py`](experiments/verify_robust_boolean_product_synchronization.py),
and
[`experiments/verify_robust_product_synchronization.py`](experiments/verify_robust_product_synchronization.py).

### Theorem 21.56 (marginal near-top poles do not control coherent synchronization)

For every odd `p`, there is a generated odd-character pole algebra of
dimension

```math
q=2^{p-1}={\sqrt N\over2},\qquad N=2^{2p},        \tag{21.340}
```

and two symmetric PSD contraction children with `G=I_q` such that their
complete degree-one generator blocks, every individual active-product
Rayleigh deficit, and their trace and average defects agree.  Nevertheless
one has relative defect `1`, while the other has relative defect `O(1/p)`.

Explicitly, if `a` is the Fourier vector of odd majority, remove its
singleton part and normalize the remaining vector `u`.  Then

```math
D_(coh)=uu^T,\qquad D_(diag)=diag(u_S^2).          \tag{21.341}
```

The majority selector exposes

```math
a^TD_(coh)a\longrightarrow1-{2\over\pi},
\qquad a^TD_(diag)a=O(1/p).                       \tag{21.342}
```

Taking spectral scale `r=sqrt N`, this is a leading
`Theta(N^(3/2))` separation between prescribed joint-selector losses even
though every separately inspected pole is asymptotically top-Rayleigh.  The
construction is a weighted contraction and does not assert separated
Boolean optima.

There is a matching positive structural criterion.  On an orthogonal
character pole span, twirl the compressed child `A` by its finite abelian
generator group.  If every diagonal deficit is at most `d`, then

```math
\boxed{D\preceq(d+||A-\overline A||_op)G.}        \tag{21.343}
```

If the group has commuting involutory generators,

```math
||A-\overline A||_op
\le {1\over2}\sum_i||[A,rho(g_i)]||_op.           \tag{21.344}
```

Thus exact or approximate equivariance converts marginal near-top data into
the relative synchronization needed by Theorem 21.55.  Without coherence,
separate product estimates cannot do so.  The full proof and 70,333 exact
checks are in
[`drafts/gram_relative_coherence_blindness.md`](drafts/gram_relative_coherence_blindness.md),
[`drafts/gram_relative_coherence_blindness_audit.md`](drafts/gram_relative_coherence_blindness_audit.md),
and
[`experiments/verify_gram_relative_coherence_blindness.py`](experiments/verify_gram_relative_coherence_blindness.py).

### Theorem 21.57 (factorized first-moment phase law and exact-sign obstruction)

Let a Cartesian pole interface have independent factor row blocks
`X_t in {+-1}^{q_t}`, total port arity `p_L=1+sum_tq_t`, and means
`mu_(t,j)=E X_(t,j)`.  Put

```math
V_L=\sum_tq_t^2,
\qquad
theta_L={1+\sum_(t,j)|mu_(t,j)|\over p_L}.         \tag{21.345}
```

Then its exact maximum linear support satisfies

```math
\boxed{
\left|{1\over n_Lp_L}\max_epsilon||W_Lepsilon||_1-theta_L\right|
\le {\sqrt{V_L}\over p_L}.}                      \tag{21.346}
```

This follows from Jensen and the variance bound
`Var(sum_j epsilon_jX_(t,j))<=q_t^2`; no independence within a factor is
needed.  Thus bounded factor arity collapses an exponentially large exact
histogram to one scalar phase for the max-over-label query.

If the represented tensor pole span has relative synchronization defect
`e_L`, and `c_L=m_Lp_L/r_L`, every labelled Boolean trust channel obeys

```math
\left|{B_L(epsilon)\over r_Ln_L}
-\left({1\over2}+c_L|d_L(epsilon)|\right)\right|
\le {e_L\over2}+c_L{\sqrt{V_L}\over p_L},         \tag{21.347}
```

where `d_L` is the signed normalized first-moment phase.  Consequently,

```math
{\max_epsilon B_L(epsilon)\over r_Ln_L}
={1\over2}+c_Ltheta_L
+O\left(e_L+c_L{\sqrt{V_L}\over p_L}\right).     \tag{21.348}
```

If the child factors are entrywise signs, have trace-zero tensor product,
`r_L/sqrt(n_L)->rho`, and the auxiliary exact-sign completion has squared
size `o(r_Ln_L)`, convergence of `c_L`, `theta_L`, and vanishing `e_L`
therefore gives a completed-parent limit

```math
{Q(P_L)\over|P_L|^{3/2}}
\longrightarrow rho(1/2+c theta).                \tag{21.349}
```

For the repeated order-16 PC.3 seed, `theta_L->1/4`, `c_L->1`, and
`rho=1`, proving the exact-sign limit `3/4` from Theorem PL.2.

The phase hypothesis is necessary.  In every factor of the same exact
regular-Hadamard child choose either a relative pole of mean `1/2` or one of
mean zero.  Both choices have exact pole closure and zero relative defect.
Arrange their types in alternating blocks, each dominating the whole prior
prefix.  The completed exact-sign parent sequence then has

```math
\limsup_L{Q(P_L)\over|P_L|^{3/2}}=1,
\qquad
\liminf_L{Q(P_L)\over|P_L|^{3/2}}={1\over2}.      \tag{21.350}
```

Thus exact synchronization and exact tensor composition do not themselves
force a thermodynamic limit; the empirical factor phase is an independent
state.  This is a nonconvergent structured signing sequence, not a statement
about the minimizing values `M_n`.  Proof, independent audit, and 44 finite
checks are in
[`drafts/factorized_port_phase_law.md`](drafts/factorized_port_phase_law.md),
[`drafts/factorized_port_phase_law_audit.md`](drafts/factorized_port_phase_law_audit.md),
and
[`experiments/verify_factorized_port_phase_law.py`](experiments/verify_factorized_port_phase_law.py).

### Theorem 21.58 (exact signs permit a fixed-arity coherence gap)

There is a symmetric regular Hadamard matrix `H_0` of order `16` and five
Boolean ports `w_1,...,w_5` for which every one of the `16` odd product
poles `z_S=product_(i in S)w_i` satisfies

```math
{z_S^TH_0z_S\over4\cdot16}={52\over64},
\qquad d_S={3\over16},                             \tag{21.351}
```

while the all-positive majority selector `x` satisfies

```math
{x^TH_0x\over4\cdot16}={32\over64},
\qquad d_(maj)={1\over2}.                         \tag{21.352}
```

Thus the coherent selector loss exceeds every separately inspected pole
loss by `5/16`.  This persists exactly under

```math
H_j=H_0^{\otimes j},\qquad
w_(i,j)=w_i\otimes\mathbf1_(16^{j-1}),            \tag{21.353}
```

at orders `n_j=16^j` and scale `sqrt(n_j)n_j`.  Since
`tr(H_j)=0`, deleting the diagonal produces a hollow exact signing without
changing any Boolean quadratic energy.  The exact displayed deficits use
the diagonal-completed Hadamard roof `sqrt(n_j)`.  For the hollow child
`A_j=H_j-diag(H_j)`, one may take the valid contraction scale
`r'_j=sqrt(n_j)+1`; then the two deficits tend respectively to `3/16` and
`1/2`, and their excess tends to `5/16`.  The result is therefore an
exact-sign, external-`n^(3/2)` certificate obstruction, not just a weighted-
contraction example, but the hollow and completed spectral roofs must not be
identified.

The qualification is essential: arity is fixed, the marginal deficit is
the nonvanishing constant `3/16`, and (21.352) concerns the prescribed
same-field selector rather than a separated full Boolean trust optimum.
At the order-16 three-port PC.3 seed the four projective selector witnesses
are exactly its four projective odd-product poles, so no such separation is
possible for any quadratic child on that seed.  At PC.3 tensor depth two
the analogous two projective sets intersect in only one of sixteen points;
the rigidity is genuinely seed-specific.  The proof and 98,467 exact checks
are in
[`drafts/exact_sign_product_coherence_gap.md`](drafts/exact_sign_product_coherence_gap.md)
and
[`experiments/verify_exact_sign_product_coherence_gap.py`](experiments/verify_exact_sign_product_coherence_gap.py).

### Theorem 21.59 (switching visibility converts marginal deficits into synchronization)

Let a group `Gamma=(Z/2Z)^s` act by commuting signed-permutation
involutions on a subspace `U`, and suppose its representation on `U` is
multiplicity-free.  Let a projectively `Gamma`-invariant Boolean pole family
span `U`.  For one representative `z_j` of each projective orbit, put

```math
p_(j,chi)=||Pi_chi(z_j/sqrt n)||_2^2,
\qquad
nu=\min_chi\max_jp_(j,chi)>0.                    \tag{21.354}
```

For a hollow signing `A`, choose `r>=||A||_op`, let every declared pole have
Rayleigh deficit at most `d`, and let `k_i` count unordered edges on which
the `i`th signed switching conjugate differs from `A`.  If `G` and `R` are
the pole Gram and normalized Rayleigh matrices and `D=G-R`, then

```math
\boxed{
0\preceq D\preceq
\left({d\over nu}+\sum_i{\sqrt{2k_i}\over r}\right)G.}
                                                               \tag{21.355}
```

The sharper nonuniform coefficient is

```math
delta_(vis)=\max_chi\min_(j:p_(j,chi)>0)
 {\bar d_j\over p_(j,chi)},                     \tag{21.356}
```

in place of `d/nu`, where `bar d_j` is the orbit-average deficit.  The proof
twirls the positive compressed defect.  Multiplicity freeness diagonalizes
the twirl, positivity lets each visible orbit bound a character defect, and

```math
||C-\bar C||_op
\le{1\over2}\sum_i||[A/r,P_i]||_op
\le\sum_i{\sqrt{2k_i}\over r}.                  \tag{21.357}
```

This identifies visibility—not multiplicity freeness alone—as the missing
quantitative datum.  On the hollow PC.3 order-16 signing,
`r=5,nu=1/2,d=1/5`, and the sharp relative spectrum is
`(0,2/5,2/5)`.  A non-Hadamard order-eight Cayley signing gives an exact
`nu=1` character-diagonal test.  Under `j`-fold PC.3 tensoring, however, the
two orbit-incidence rows tensor and

```math
\nu_j=2^{-j}=N_j^{-1/4},\qquad N_j=16^j.         \tag{21.358}
```

Thus this declared symmetry proves vanishing relative defect from uniform
marginal deficit only at the stronger scale `d_j=o(N_j^(-1/4))`; visibility
itself can be lost under growth.  The theorem neither constructs a visible
action in arbitrary near-minimizers nor makes a generic `o(n^2)` edit budget
small at `r=Theta(sqrt n)`.  Full proof, audit, and 3,363 exact checks are in
[`drafts/sign_switching_visibility_synchronization.md`](drafts/sign_switching_visibility_synchronization.md),
[`drafts/sign_switching_visibility_synchronization_audit.md`](drafts/sign_switching_visibility_synchronization_audit.md),
and
[`experiments/verify_sign_switching_visibility.py`](experiments/verify_sign_switching_visibility.py).

### Theorem 21.60 (a gapped conjugation twirl forces collective coherence)

Let `T` be a real symmetric contraction, let `Z` be a finite nonempty set of
Boolean vectors, and put

```math
B_z=M_zTM_z,
\qquad
K={1\over|Z|}\sum_(z\in Z)B_z,
\qquad M_z=diag(z).                               \tag{21.359}
```

Suppose

```math
spec(K)\subseteq[-1,1-gamma]\cup\{1\}.           \tag{21.360}
```

Then every Boolean `x in span_R(Z)` satisfies the operator inequality

```math
\boxed{
I-M_xTM_x\preceq {2\over gamma}(I-K).}           \tag{21.361}
```

Indeed the eigenvalue-one space of `K` is the intersection of the kernels of
the positive defects `I-B_z`, equivalently the common eigenvalue-one spaces
of the `B_z`.  Linearity puts this common kernel inside the eigenvalue-one
space of `M_xTM_x`; the gap of `K` then compares the two positive operators.
Hence,
for `||s||_2^2=n`, the selector deficit is at most

```math
{2\over gamma|Z|}\sum_(z\in Z)
 \left(1-{s^TM_zTM_zs\over n}\right).            \tag{21.362}
```

For the order-16 PC.3 seed, the conjugation twirl over its four active poles
obeys the exact finite identity

```math
K_0^3=K_0,
\qquad mult_(K_0)(1,0,-1)=(3,10,3).              \tag{21.363}
```

At tensor depth `j`, the `4^j` active odd-product poles give
`K_j=K_0^(tensor j)`, so `gamma=1` at every depth.  Since every majority
selector lies in their span, define

```math
d_z(s)=1-{s^TM_zT_jM_zs\over n_j}
```

and define `d_epsilon(s)` analogously with the endpoint selector in place of
`z`.  Every diagonal switching `s` of the PC.3 child then satisfies

```math
\boxed{
d_epsilon(s)\le {2\over4^j}\sum_zd_z(s)
                 \le2\max_zd_z(s).}             \tag{21.364}
```

Thus vanishing marginal product deficits force every selector deficit to
vanish throughout this exact-sign perturbation class.  Diagonal deletion
preserves all Boolean energies because the completed signing has trace zero.
For the hollow child the scale `r'_j=sqrt(n_j)+1` is a valid contraction
roof; its deficits are `(1-c_j)+c_jd` with
`c_j=sqrt(n_j)/(sqrt(n_j)+1)`, so (21.364) and the same vanishing implication
continue to hold.  The inherited `sqrt(n_j)n_j` scale alone is a Boolean
energy roof, not necessarily the hollow operator norm.
The theorem does not cover genuinely non-switching sign children or compare
full Boolean trust optima.  Proof and exact order-256/4096 stress tests are
in
[`drafts/pc3_diagonal_switching_coherence.md`](drafts/pc3_diagonal_switching_coherence.md)
and
[`experiments/verify_pc3_diagonal_switching_coherence.py`](experiments/verify_pc3_diagonal_switching_coherence.py).

### Theorem 21.61 (the fixed five-port seed cannot be diluted by monomial amplification)

For the five-port seed of Theorem 21.58, every odd seed product has
normalized Rayleigh value `13/16`, while its majority selector has value
`1/2`.  First let `p` be odd and allow arbitrary monomial seed factors,
port-dependent Boolean auxiliary factors, and arbitrary weighted direct
mixtures of blocks.  If the majority output agrees projectively with the
prescribed selector `x` on every concrete seed row in every block, the odd
active-product deficits obey

```math
\boxed{
{1\over2^{p-1}}\sum_(B\text{ odd})d_B\ge {3\over32},
\qquad \max_(B\text{ odd})d_B\ge {3\over32}.}    \tag{21.365}
```

Indeed the positive-roof seed monomials are exactly the subgroup
`K={0,5,9,12,17,20,24,29}` of `F_2^5`; two concrete seed rows agree on every
`K`-monomial but `x` separates them.  Thus each block's exponent map is
nontrivial modulo `K`.  At least half of the uniformly sampled odd subsets
map outside `K`, where every product loses at least `3/16`.  Averaging over
blocks leaves one common labelled channel exposed.

There is a sharper quantitative law for the common-factor subclass.  Let
`p` be odd, let each lifted port be a signed odd seed monomial times one
common auxiliary Boolean factor `u`, and suppose the lifted majority witness
remains projectively `x tensor u`.  If `s=<u,Su>` for a self-adjoint
auxiliary contraction `S`, then the full active Fourier pole has an odd seed
product on the first factor and

```math
d_(sel)=1-{s\over2},
\qquad
d_(full)=1-{13s\over16},
\qquad
\boxed{d_(full)\ge {3\over8}d_(sel).}            \tag{21.366}
```

The constant is sharp at `s=1`, and the stable form is

```math
d_(full)\ge {3\over8}d_(sel)-{11\over8}eta.      \tag{21.367}
```

Repeated selector composition is worse: for `L` seed copies, an active
full-leaf pole has value `(13/16)^L` and hence deficit tending to one.
Exact-sign tensor decoration and trace-zero hollowing satisfy the theorem;
the declared completed roof and the hollow operator roof differ by at most
one.

Thus arbitrary monomial tensor decoration and direct-mixture dilution, odd
replication or cancelling pairs, and normalized-operator-`o(1)` repairs
cannot turn the fixed seed into a vanishing-marginal obstruction.  A
nonmonomial/nonlocal lift, a leading cross-block completion, or a genuinely
growing seed remains outside scope.  Full trust optima are not compared.
The proof, independent audit, and 140,457 checks are in
[`drafts/fixed_seed_amplification_no_go.md`](drafts/fixed_seed_amplification_no_go.md),
[`drafts/fixed_seed_amplification_no_go_audit.md`](drafts/fixed_seed_amplification_no_go_audit.md),
and
[`experiments/verify_fixed_seed_amplification_no_go.py`](experiments/verify_fixed_seed_amplification_no_go.py).

### Theorem 21.62 (sparse exact-sign flips create vanishing marginal loss and constant coherent loss)

Let `H` be a symmetric entrywise-sign matrix of square order `N`, with
`H^2=NI` and `tr(H)=0`.  Put `r=sqrt N`.  Let `Z` be Boolean positive
`r`-eigenvectors and let another Boolean positive eigenvector `x` satisfy

```math
alpha_N=\max_(z\in Z){|x^Tz|\over N}=o(1),
\qquad \log|Z|=o(N^{3/2}).                       \tag{21.368}
```

For every fixed `0<kappa<1`, there is a symmetric exact-sign matrix `H'`
with `r'=||H'||_op=(1+o(1))sqrt N` such that

```math
0\le1-{z^TH'z\over r'N}
 \le kappa alpha_N^2+o(1)\quad(z\in Z),          \tag{21.369}
```

uniformly, while

```math
{x^TH'x\over r'N}=1-kappa+o(1).                 \tag{21.370}
```

Construct `H'` by independently flipping, with probability
`q=kappa/sqrt N`, every off-diagonal edge positive for `x`.  The expected
energy loss of a pole `z` is exactly

```math
q\{rN+(x^Tz)^2-N\},                             \tag{21.371}
```

whereas `x` has `[N(N-1)+rN]/4` eligible edges.  Scalar Bernstein is uniform
under (21.368), matrix Bernstein contributes only
`O(N^(1/4)sqrt(log N)+log N)`, and

```math
E H'=(1-q)H-qxx^T+q(diag(H)+I).                 \tag{21.372}
```

The first two terms commute and have norm at most `sqrt N+o(sqrt N)`;
Frobenius norm gives the reverse bound.  Hollowing preserves every Boolean
energy and changes the operator roof by at most one.

The diffuseness hypothesis holds in the PC.3 tensor tower.  For depth `j`,
take endpoint signs

```math
epsilon_0=1,
\qquad epsilon_(2t+1)=(-1)^t,
\qquad epsilon_(2t+2)=1.                        \tag{21.373}
```

Then its majority selector `x_epsilon` and `4^j` active poles satisfy

```math
\boxed{
\max_(z\in Z_j){|x_epsilon^Tz|\over16^j}
 =O(j^{-1/2}).}                                  \tag{21.374}
```

After cancelling the common base pole, rows are iid under the three-atom law
`(X,Y)=(1,1),(1,-1),(-1,1)` with probabilities `1/4,1/2,1/4`.  The selector
is the sign of `1+sum_t((-1)^tX_t+Y_t)` and every pole is a product of local
characters in `{1,X,Y,XY}`.  Berry--Esseen gives uniform
`O(m^(-1/2))` anti-concentration after arbitrary factor deletion; each
nonconstant character has mean modulus at most `1/2`.  Eliminating its
factors gives a geometrically weighted error series and proves (21.374).

Combining (21.369)--(21.374) yields hollow exact signings with every declared
active-product deficit `o(1)` but one prescribed same-field selector deficit
`kappa+o(1)`.  This is the sought scalable exact-sign **certificate**
obstruction.  It does not separate full Boolean trust optima: another spin
may repair the response at this stage; Theorem 21.64 subsequently rules out
that repair for one leading labelled field.  Proof, independent audits, and exact finite
checks through depths two and three plus the factor law through depth ten are
in
[`drafts/pc3_diagonal_switching_coherence.md`](drafts/pc3_diagonal_switching_coherence.md),
[`drafts/pc3_sparse_flip_coherence_audit.md`](drafts/pc3_sparse_flip_coherence_audit.md),
and
[`experiments/verify_pc3_diagonal_switching_coherence.py`](experiments/verify_pc3_diagonal_switching_coherence.py).

### Corollary 21.63 (the diffuse PC.3 field reaches leading scale on a sublinear shore)

Let `W_j` have the `p_j=2j+1` repeated PC.3 ports as columns, let `epsilon_j`
be the periodic endpoint in (21.373), put `h_j=W_jepsilon_j`, and set
`N_j=16^j`.  Then

```math
\boxed{
 ||h_j||_1=
 \sqrt{7\over2\mathop{\rm pi}}\,N_j\sqrt j+O(N_j).}
                                                               \tag{21.375}
```

In particular, if

```math
m_j=\left\lfloor\lambda\sqrt{N_j/j}\right\rfloor
```

for fixed `lambda>0`, repeating every port `m_j` times gives

```math
m_j||h_j||_1=
 \left(\lambda\sqrt{7\over2\mathop{\rm pi}}+o(1)\right)
 N_j^{3/2},
\qquad
m_jp_j=O(\sqrt{N_jj})=o(N_j).                  \tag{21.376}
```

Indeed the uniform PC.3 row law turns `||h_j||_1/N_j` into
`E|1+sum_(t<j)A_t|`, with the independent variables from Theorem 21.62.
Their means are `(-1)^t/2`, and their variances are `3/4` for even `t` and
`11/4` for odd `t`.  Pairing consecutive increments gives iid bounded
centred variables of variance `7/2`.  Wasserstein Berry--Esseen applied to
the one-Lipschitz absolute-value function gives the Gaussian first absolute
moment with `O(1)` error; the initial constant and possible odd leftover
increment are bounded.  This proves (21.375).

Thus the sparse-flip certificate can be queried at the full `N_j^(3/2)`
scale by only a sublinear number `s_j=m_jp_j` of port vertices.  Any exact-
sign completion internal to that shore costs at most
`O(s_j^2)=O(N_jj)=o(N_j^(3/2))`.  The remaining issue is genuinely
contextual: `m_jh_j` is the labelled channel in which the replicas use the
prescribed endpoint, whereas in a full parent the new spins and the child
spin optimize jointly and may repair it.

### Theorem 21.64 (a labelled full Boolean response retains a fixed sparse-flip gap)

For a hollow child `A`, define its labelled trust response by

```math
\mathcal B_A(g)=
 \max_{y\in\{\mathord\pm1\}^N,\ \sigma\in\{\mathord\pm1\}}
 \left\{{\sigma\over2}y^TAy+y^Tg\right\}.      \tag{21.377}
```

Fix `0<kappa<1`.  Let `A_j` be the unflipped hollow PC.3 Hadamard child.
There is a choice of sparse-flip children `A'_j` from the construction of
Theorem 21.62 for its periodic selector such that the following holds.  Set

```math
rho=\sqrt{2/\mathop{\rm pi}},\qquad
s=\sqrt{1-rho^2},\qquad b={\lambda\sqrt7\over2}.
```

If

```math
0<b<\min\left\{{2kappa(rho-s)\over rho^2},
                  2(2-kappa)rho\right\},       \tag{21.378}
```

then, for `m_j=floor(lambda sqrt(N_j/j))`,

```math
\boxed{
 \mathcal B_{A_j}(m_jh_j)-\mathcal B_{A'_j}(m_jh_j)
 \ge(\delta+o(1))N_j^{3/2},}
\quad
\delta=\min\left\{
b(rho-s)-{b^2rho^2\over2kappa},
brho-{b^2\over2(2-kappa)}\right\}>0.            \tag{21.379}
```

To prove this, put `T_j=H_j/sqrt(N_j)`, `e_j=x_j/sqrt(N_j)`, and
`f_j=h_j/||h_j||_2`.  The exact row law and Corollary 21.63 give

```math
T_je_j=e_j,\qquad T_jf_j=f_j,\qquad
e_j^Tf_j\longrightarrow rho,\qquad
||h_j||_2^2=N_j(7j/4+O(1)).                     \tag{21.380}
```

On the simultaneous concentration event from Theorem 21.62,

```math
{H'_j\over\sqrt{N_j}}=T_j-kappa e_je_j^T+o_{op}(1). \tag{21.381}
```

Hollowing contributes only `o_op(1)`.  Relax `y=sqrt(N_j)u` to the unit
sphere and decompose `u` into the two eigenspaces of the involution `T_j`.
Writing `t=e_j^Tu`, the negative component only subtracts from the quadratic
term, while (21.380) gives

```math
f_j^Tu\le rho|t|+s+o(1).
```

Since `m_j||h_j||_2/N_j -> b`, the positive quadratic channel is at most

```math
{1\over2}+bs+max_{a\ge0}\{brho a-kappa a^2/2\}+o(1)
\le {1\over2}+bs+{b^2rho^2\over2kappa}+o(1).    \tag{21.382}
```

In the negative channel, if `a` is the norm of the positive-eigenspace
component, the quadratic and linear terms are at most
`1/2-(1-kappa/2)a^2+ba+o(1)`.  Hence that channel is at most

```math
{1\over2}+{b^2\over2(2-kappa)}+o(1).            \tag{21.383}
```

For the unflipped child, the Boolean pole `x_j=sgn(h_j)` gives the lower
bound `1/2+brho+o(1)`.  Subtraction proves (21.379).  For example,
`kappa=1/2,lambda=1/10` gives `delta>0.0146`.

This upgrades the sparse-flip construction from a prescribed-spin
certificate failure to a genuine full Boolean response separation for one
fixed **labelled linear-field context**.  It is not yet an unconstrained
exact-sign parent separation: if the repeated port vertices are free, their
endpoint word optimizes too.  Proof and normalization details are in
[`drafts/pc3_diagonal_switching_coherence.md`](drafts/pc3_diagonal_switching_coherence.md).
An independent audit is in
[`drafts/pc3_labelled_response_gap_audit.md`](drafts/pc3_labelled_response_gap_audit.md).

### Corollary 21.65 (the direct free-shore realization has super-target cap)

The repeated-port implementation of Theorem 21.64 cannot simply leave its
shore spins free.  For the all-positive endpoint `eta`, the exact PC.3 row
law gives

```math
{||W_jeta||_1\over N_j}
=\mathbb E\left|1+\sum_{t<j}(X_t+Y_t)\right|
\ge1+{j\over2}.                                  \tag{21.384}
```

Consequently `m_j=floor(lambda sqrt(N_j/j))` copies of every port expose
field reward at least

```math
(\lambda/2+o(1))N_j^{3/2}\sqrt j.               \tag{21.385}
```

This dominates an `O(N_j^(3/2))` child energy and the entire internal energy
`O(N_jj)` of the `O(sqrt(N_jj))`-vertex shore.  Thus every exact-sign
completion of this direct unconstrained-shore construction has cap
`Omega(N_j^(3/2)sqrt j)`.  A useful unlabelled realization must balance or
restrict endpoint words, or use a different interface; subleading shore
calibration cannot suffice.

### Theorem 21.66 (rowwise microcanonical compilation balances every endpoint)

Let `g in Z^N`, let `s>=2`, and suppose `|g_i|<=s` and
`g_i=s (mod 2)` for every row.  There is an exact-sign matrix
`B in {+-1}^{N times s}` and a target endpoint `eta_* in {+-1}^s` such that

```math
B eta_*=g,
\qquad
\max_eta\left\|B eta-{<eta,eta_*>\over s}g\right\|_1
\le\sqrt{2s}\,N+C s^{3/2}\sqrt N,              \tag{21.386a}
```

and consequently

```math
||B||_(infinity to 1):=\max_eta||B eta||_1
\le||g||_1+\sqrt{2s}\,N+C s^{3/2}\sqrt N.      \tag{21.386}
```

Here `C` is absolute.  To prove this, switch columns so `eta_*=1` and sample
row `i` uniformly among sign vectors of sum `g_i`.  For a fixed endpoint
`eta`, with `c=1^Teta`, exchangeability and sampling without replacement give

```math
\mathbb E(b_i^Teta)={g_ic\over s},
\qquad \mathop{\rm Var}(b_i^Teta)\le2s.          \tag{21.387}
```

Thus, writing `a=c/s`,

```math
\sum_i\mathbb E|b_i^Teta-a g_i|\le\sqrt{2s}\,N.
                                                               \tag{21.388}
```

The `N` absolute centred row responses are independent and bounded by `2s`.
Hoeffding at deviation `C s sqrt(Ns)` followed by a union bound over the
`2^s` endpoints proves (21.386a), and the triangle inequality gives
(21.386).

For the PC.3 target `g_j=m_jh_j`, the parity condition is automatic because
`h_j` is a sum of `2j+1` signs in every row.  With
`s_j=m_j(2j+1)=Theta(sqrt(N_jj))`,

```math
N_j\sqrt{s_j}+s_j^{3/2}\sqrt{N_j}=o(N_j^{3/2}). \tag{21.389}
```

Hence a balanced exact-sign cross block can realize the labelled field with
global cross cap `||g_j||_1+o(N_j^(3/2))`.  This repairs the endpoint-bias
failure of the repeated-column lift.  The stronger affine estimate
(21.386a), not the scalar cap bound alone, is what Theorem 21.67 uses to
control every competing endpoint and preserve the response gap.

### Theorem 21.67 (balanced compilation realizes the labelled gap in unconstrained exact-sign parents)

Use `A_j,A'_j,g_j=m_jh_j` and parameters satisfying Theorem 21.64, and let
`B_j` be a microcanonical compiler from Theorem 21.66 with
`s_j=Theta(sqrt(N_jj))`.  Complete the new shore by any common hollow
exact-sign matrix `C_j`, and form

```math
P_j=\begin{pmatrix}A_j&B_j\\B_j^T&C_j\end{pmatrix},
\qquad
P'_j=\begin{pmatrix}A'_j&B_j\\B_j^T&C_j\end{pmatrix}.          \tag{21.390}
```

Then these are unconstrained hollow exact signings of common order
`N_j+s_j=(1+o(1))N_j`, and

```math
\boxed{Q(P_j)-Q(P'_j)\ge(\delta+o(1))N_j^{3/2},}              \tag{21.391}
```

with `delta>0` from (21.379).

Indeed, for every free shore endpoint `eta`, put
`a=<eta,eta_*>/s_j`.  The strengthened microcanonical bound (21.386a)
and Lipschitzness of Boolean trust response in field `l_1` give

```math
\mathcal B_{A'_j}(B_jeta)
\le\mathcal B_{A'_j}(a g_j)+o(N_j^{3/2}).        \tag{21.392}
```

Trust response is even in its field.  Repeating the two spherical estimates
of Theorem 21.64 with strength `|a|b` shows that their upper bounds are
increasing for `0<=|a|<=1`; hence (21.392) is at most
`\mathcal B_{A_j}(g_j)-(delta+o(1))N_j^{3/2}` uniformly in `eta`.
Precisely, changing `y` to `tau y` in the cross term gives, for either child
`D`,

```math
Q\begin{pmatrix}D&B_j\\B_j^T&C_j\end{pmatrix}
=\max_{eta,tau}\left[
 \max_y\left\{{tau\over2}y^TDy+y^TB_jeta\right\}
 +{tau\over2}eta^TC_jeta\right].               \tag{21.393}
```

Thus the flipped parent is at most
`max_eta \mathcal B_{A'_j}(B_jeta)+Q(C_j)`, whereas the target endpoint gives
the unflipped parent at least `\mathcal B_{A_j}(g_j)-Q(C_j)`.  Since
`2Q(C_j)=O(s_j^2)=o(N_j^{3/2})`, (21.391) follows.

This is a physical, all-spins-free contextual collision at the target dense
quadratic scale.  It does **not** compare minima over signings, select these
children as near-minimizers, or provide a cross-order recurrence for `M_n`.
It proves that a strict balanced response compiler can preserve coherent
information invisible to every active product pole.

## 22. Finite-port response dimension

For a message `m in R^q`, write `R_m(g)=max_j(m_j+g_j)`.  On projective
classes define

```math
d_proj([m],[m'])
=inf_(c in R)sup_(g in R^q)|R_m(g)-R_(m')(g)-c|.             \tag{22.1}
```

### Theorem 22.1 (polyhedral response rate)

If `S subset R^q/R1` is a fixed compact finite ordinary polyhedral complex
of topological dimension `d`, arbitrary coordinate futures make `S` the
coarsest exact projective response state and

```math
d_proj([m],[m'])={1\over2}osc(m-m').                         \tag{22.2}
```

As `epsilon` tends to zero,

```math
Cov(S,d_proj,epsilon)=Theta_S(epsilon^(-d)).                 \tag{22.3}
```

The same exponent governs optimal response codebooks up to the usual
packing/covering factor two.  A bounded scalar baseline adds one scalar
precision term; an unbounded baseline has no finite absolute codebook.

#### Proof

Coordinate-pinning futures recover every coordinate difference.  Optimizing
one future-independent scalar calibration gives the midrange and (22.2).
A finite `d`-dimensional polyhedral complex has an
`O_S(epsilon^(-d))` cover, while a relatively open ball in a maximal cell
gives the matching packing. `square`

### Theorem 22.2 (all-finite tropical-kernel carrier)

For finite `K in R^(p times q)`, let

```math
(U_Ku)_j=max_i\{u_i+K_(ij)\},
\qquad C_K=\{[U_Ku]:u in R^p\}.                              \tag{22.4}
```

Then `C_K` is a compact finite polyhedral complex of ordinary dimension
`d_K<=min(p-1,q-1)`, and hence has response rate (22.3) with `d=d_K`.
At a tree node with child messages `b_c1+r_c in R^p`, unary score `a`, and
outgoing kernel `K`, its exact state update is

```math
u=a+sum_c r_c,
\quad w=U_Ku,

b_out=sum_cb_c+max_jw_j,
\quad r_out=w-(max_jw_j)1.                                  \tag{22.5}
```

Thus the response image is an exact semantic congruence.  If all microscopic
data lie in one lattice `eta Z`, the normalized carriers are finite and
exactly invariant.  Rounding `F` arbitrary factors first gives a carrier for
the rounded model and total conditional-optimum error at most `Feta/2`; it
does not give depth-independent error at fixed mesh.

#### Proof

On each selector region `U_K` is affine.  The projective image of that
polyhedron is a polyhedron, and

```math
min_i(K_(ij)-K_(ik))
<=(U_Ku)_j-(U_Ku)_k
<=max_i(K_(ij)-K_(ik))                                      \tag{22.6}
```

makes every normalized image bounded.  A finite common refinement is the
classical tropical type complex.  Formula (22.5) is max-sum elimination;
the lattice and rounding claims follow assignmentwise. `square`

The geometry in Theorem 22.2 is classical tropical convexity.  The response
metric, rate statement, and separation between static image and reusable
arithmetic congruence are the deductions used here.

For a ferromagnetic `q`-state Potts edge of strength `K>0`, direct calculation
gives the exact image

```math
C_K=\{r in [-K,0]^q:max_jr_j=0\},                            \tag{22.7}
```

so `d_K=q-1` and its sharp contextual codebook has
`Theta_q((K/epsilon)^(q-1))` states.  This independently derived benchmark
also gives the binary clipped-gap update.  Bounded Viterbi survivor families
and bounded parity tables are instances of Theorem 22.1, but not literal
all-finite compact-kernel images because their structural infeasibility uses
infinite semiring entries.

## 23. Contraction turns response scale into forgetting time

Let a deterministic discounted-control block of duration `h` act on a
terminal value vector by

```math
(U_(K,h)v)_i=max_j\{K_(ij)+lambda^h v_j\},
\qquad 0<lambda<1.                                          \tag{23.1}
```

Blocks compose by the semidirect max-plus rule

```math
(K odot_h L)_(ik)=max_j\{K_(ij)+lambda^hL_(jk)\}.           \tag{23.2}
```

### Theorem 23.1 (discounted response isometry and rate)

If arbitrary deterministic routing continuations and initial-state probes
are declared, then the exact depth-`h` response metrics are

```math
d_h(v,w)=lambda^h||v-w||_infinity,

dbar_h([v],[w])={lambda^h\over2}osc(v-w).                    \tag{23.3}
```

Suppose immediate rewards lie in `[-R,R]` and a horizon-`H` terminal reward
is zero. The realizable value class is exactly the cube

```math
[-B_H,B_H]^n,
\qquad B_H=R{1-lambda^H\over1-lambda}.                       \tag{23.4}
```

Consequently its optimal depth-`h`, error-`epsilon` absolute response
codebook satisfies

```math
log_2M=Theta\left(
n log\left(1+{B_Hlambda^h\over epsilon}\right)\right),      \tag{23.5}
```

with exponent `n-1` projectively. The constants are absolute in the
nontrivial resolution regime; explicit coordinate-grid packing and covering
bounds govern the one-state endpoint.

#### Proof

The maximum inequality contracts sup distance by `lambda^h`. A hard
zero-reward route to a coordinate attaining the sup difference gives
equality, and scalar calibration gives the projective formula. Every value
lies in (23.4). Conversely, a statewise self-loop reward

```math
r_i=v_i{1-lambda\over1-lambda^H}                             \tag{23.6}
```

realizes every point of that cube without exceeding `R`. Grids of mesh
`epsilon/lambda^h` give matching covers and separated packings. `square`

The infinite-horizon radius is `R/(1-lambda)`, so old information becomes
`epsilon`-irrelevant after `Theta(log(B/epsilon)/(-log lambda))` steps. If a
fresh local error `delta` enters at every Bellman update, however,

```math
e_t<=lambda^te_0+delta{1-lambda^t\over1-lambda},             \tag{23.7}
```

and the denominator is attained by self-loops. This benchmark quantifies the
dynamic factor missing from static response entropy: contraction changes the
scale at which the same response image must be covered.

## 24. Critical branching creates an asymptotic rare-event state

The next result is an imported benchmark, not a theorem about deterministic
finite ports.  Consider a boundary-case branching random walk satisfying

```math
E\sum_(|u|=1)e^{-V(u)}=1,
\qquad E\sum_(|u|=1)V(u)e^{-V(u)}=0,                        \tag{24.1}
```

together with Madaule's nonlattice, second-moment, and logarithmic
integrability hypotheses.  Write

```math
W_n=\sum_(|u|=n)e^{-V(u)},
\qquad Z_n=\sum_(|u|=n)V(u)e^{-V(u)}.                       \tag{24.2}
```

### Imported Theorem 24.1 (derivative-mass response collapse)

Conditional on survival, `W_n` tends to zero, `Z_n` tends to a positive
`Z_infinity`, and the extremal process centered by `(3/2)log n` converges,
conditionally on `Z_infinity`, to a Cox decorated Poisson process whose
cluster-leader intensity is

```math
lambda Z_infinity e^x dx.                                  \tag{24.3}
```

For a fixed decoration law and every nonnegative compactly supported test
function `f`, its conditional Laplace response is

```math
-\log E[e^{-<E_Z,f>}\mid Z]
=Z A_D(f).                                                  \tag{24.4}
```

Thus the entire unmarked limiting Laplace-functional query family has one
realized scalar state `Z`.  It composes under a cut at generation `r` by the
exact limiting smoothing transform

```math
Z_infinity=\sum_(|u|=r)e^{-V(u)}Z_infinity^(u).             \tag{24.5}
```

At finite depth the exact identity instead is

```math
Z_(r+m)=\sum_(|u|=r)e^{-V(u)}
       (Z_m^(u)+V(u)W_m^(u)),                               \tag{24.6}
```

so the prelimit state is at least the pair `(W,Z)`.  The scalar collapse is
created by critical centering and the vanishing of `W`, not by an exact
finite-depth quotient.

#### Source and deduction

Madaule's [Theorem 1.1](https://arxiv.org/abs/1107.2543) proves joint
convergence of the shifted extremal process and `Z_n`, with an independent
decorated exponential Poisson limit.
Undoing its random `log Z_infinity` shift gives (24.3), and the Poisson
Laplace functional gives (24.4).  Expanding the definition of the derivative
martingale below a generation-`r` cut gives (24.6); sending `m` to infinity
gives (24.5).  The corresponding branching-Brownian decorated process is due
to [Aidekon--Berestycki--Brunet--Shi](https://arxiv.org/abs/1104.3738).

This state is query-relative.  If future queries retain a genealogy label,
two branch-mass allocations with the same total `Z` are distinguishable; the
state must then be the derivative-mass measure over labels.  The benchmark
therefore supplies a genuine orthogonal mechanism—renormalized rare-event
universality—while also validating the rule that the allowed future query
class determines what can be compressed.

## 25. A finite-dictionary response sparsification law

Let `X` be finite, `|X|=L`, and fix public features
`phi_e:X->[-1,1]`, `1<=e<=m`. For sign coefficients `a`, put

```math
H_a(x)=\sum_ea_e\phi_e(x),
\qquad V_Phi=\max_x\sum_e\phi_e(x)^2.                     \tag{25.1}
```

### Theorem 25.1 (universal Bernoulli-mask response code)

Fix `0<p<=1/2`, put `q=1-p`, and let `E>0`. Define

```math
Delta=2L\exp\left\{-{E^2\over2(pV_Phi/q+E/3)}\right\}
      +\exp(-pm/8).                                       \tag{25.2}
```

If `Delta<=delta<1`, then there is one public list of at most

```math
N=\left\lceil{(m+1)\log2\over\log(1/delta)}\right\rceil  \tag{25.3}
```

masks, each retaining at most `floor((1-p/2)m)` features, such that every
`a in {-1,1}^m` has a listed mask `S` satisfying

```math
\sup_x\left|H_a(x)-{1\over q}\sum_(e in S)a_e\phi_e(x)\right|
\le E.                                                     \tag{25.4}
```

Thus the absolute response code uses at most

```math
\lfloor(1-p/2)m\rfloor+\lceil\log_2N\rceil                \tag{25.5}
```

bits. The same error holds after every shared max-type future
`T_KH(y)=max_x(H(x)+K(x,y))`, without a depth factor.

#### Proof

Retain each feature independently with probability `q` and importance-weight
it by `1/q`. For fixed `a,x`, the error summands are centered, bounded by one,
and have total variance at most `pV_Phi/q`. Bernstein followed by a union
bound over `X` gives the first term of (25.2). Chernoff says that fewer than
`pm/2` erasures has probability at most the second term.

For fixed `a`, a mask is therefore bad with probability at most `delta`.
For `N` independent masks the expected number of coefficient signings missed
by the whole list is at most `2^m delta^N<=1/2`; hence one deterministic list
covers all signings. The decoder stores its mask index and the retained
signs. Finally, max-type response operators are sup-norm nonexpansive. `square`

A convenient regime follows by setting `t=log(8L)` and

```math
p=\min\{1/2,E^2/(8V_Phi t)\}.                              \tag{25.6}
```

If `E>=4t/3` and `pm>=8log4`, then `delta=1/2`, `N=m+1`, and

```math
b\le m-\min\left\{{m\over4},
 {mE^2\over16V_Phi\log(8L)}\right\}
 +\lceil\log_2(m+1)\rceil.                                \tag{25.7}
```

This law applies to every public Boolean monomial or bounded-CSP dictionary,
and to code/coset correlation landscapes `H_a(c)=a dot c`. In the latter,
uniform correlation error `E` gives nearest-code-distance error `E/2`.
These applications share no special quadratic algebra: the common resource
is a finite exposed row set with controlled feature variance.

The theorem is a strict multi-model upper law, not an internal-family or
algorithmic result. Its sparse weighted centers may leave the original model
class; a list can be existential; and duplicating the approximated landscape
several times in one future incurs the corresponding multiplicity. It
therefore separates **one-shot semantic response compression** from an
invariant compositional quotient.

## 26. Discrete adversarial thermodynamics through a lower spectral radius

Fix a spin alphabet `[q]`, a finite disorder alphabet `D`, and bounded local
rewards `h_d:[q] times [q]->R`. For a disorder word
`w=(d_1,...,d_n)`, put

```math
H_w(\sigma)=\sum_(j=1)^n h_(d_j)(\sigma_(j-1),\sigma_j),
\qquad
Z_(n,beta)(w)=\sum_\sigma e^(beta H_w(\sigma)),            \tag{26.1}
```

and retain the discrete adversarial choice

```math
F_(n,beta)=\min_(w in D^n)\log Z_(n,beta)(w),
\qquad
G_n=\min_(w in D^n)\max_\sigma H_w(\sigma).              \tag{26.2}
```

### Theorem 26.1 (adversarial finite-range thermodynamic limit)

Let

```math
T_d(i,j)=e^(beta h_d(j,i)),
\qquad mathcal T_beta=\{T_d:d in D\}.
```

Then, for every fixed `beta>0`,

```math
p(\beta):=\lim_(n->infinity){F_(n,beta)\over n}
=\log\check\rho(mathcal T_beta),                          \tag{26.3}
```

where

```math
\check\rho(mathcal T)
=\lim_n\min_(d_1,...,d_n)
 ||T_(d_n)\cdots T_(d_1)||_1^(1/n).                      \tag{26.4}
```

Moreover the zero-temperature adversarial density exists and satisfies

```math
g:=\lim_n{G_n\over n}
=\lim_(beta->infinity){p(beta)\over beta},
\qquad
0\le {p(beta)\over beta}-g\le{\log q\over beta}.         \tag{26.5}
```

No distribution or convex hull on the disorder alphabet is introduced.

#### Proof

For every nonnegative `q by q` matrix `P`,

```math
||P||_1\le 1^TP1\le q||P||_1.                            \tag{26.6}
```

The log of the minimum norm of a length-`n` product is subadditive, since
minimizing products can be concatenated and the matrix norm is
submultiplicative. Fekete and (26.6) prove (26.3).

There are `q^(n+1)` spin paths, so for every fixed disorder word

```math
\max_\sigma H_w(\sigma)
\le\beta^{-1}\log Z_(n,beta)(w)
\le\max_\sigma H_w(\sigma)+{(n+1)\log q\over\beta}.      \tag{26.7}
```

Take the minimum over the same word set. Because the middle term divided by
`n` converges for each `beta`, the limsup--liminf gap of `G_n/n` is at most
`log(q)/beta`. Sending `beta` to infinity proves convergence, and then the
same sandwich proves (26.5). `square`

### Theorem 26.2 (contractive projective cavity and finite mean-pressure state)

Assume, for a fixed `beta`, that

```math
0<a\le T_d(i,j)\le b<infinity                              \tag{26.8}
```

uniformly. On the probability simplex define

```math
tau_d(p)={T_dp\over||T_dp||_1},
\qquad r_d(p)=\log||T_dp||_1,
\qquad
(mathcal Vf)(p)=\min_d\{r_d(p)+f(tau_d(p))\}.             \tag{26.9}
```

There are a scalar `lambda` and a continuous potential `u` such that

```math
mathcal Vu=u+lambda,
\qquad
lambda=\log\check\rho(mathcal T_beta),
\qquad
Lip_(d_H)(u)\le{1\over1-kappa},                          \tag{26.10}
```

where

```math
kappa\le{b-a\over b+a}<1.                                \tag{26.11}
```

Furthermore, a Hilbert-metric `delta`-net with

```math
O_(q,a,b)(delta^(-(q-1)))                                \tag{26.12}
```

states has minimum cycle mean `lambda_delta` satisfying

```math
|lambda_delta-lambda|\le{delta\over1-kappa}.             \tag{26.13}
```

This approximates asymptotic mean pressure. It does not claim a
depth-independent additive approximation to every unnormalized finite-
horizon response.

#### Proof

Every update lands in the invariant compact interior set

```math
X=\{p in Delta_(q-1):a/(qb)\le p_i\le b/(qa)\}.
```

[Birkhoff's projective contraction theorem](https://doi.org/10.1090/S0002-9947-1957-0087058-6)
and the cross-ratio bound give (26.11), while
`|r_d(p)-r_d(p')|<=d_H(p,p')`. Consequently

```math
Lip(mathcal Vf)\le1+kappa Lip(f).                         \tag{26.14}
```

On normalized `1/(1-kappa)`-Lipschitz functions, Arzela--Ascoli and Schauder
give the additive eigenfunction. Sup-norm nonexpansiveness yields

```math
||mathcal V^n0-(u+n lambda)||_infinity\le||u||_infinity.
```

For fixed interior `p`, `||Pp||_1` and `||P||_1` differ by a factor
independent of product length, identifying `lambda` with (26.4).

Round each `tau_d(p)` at a net point to `p'`. With
`L=1/(1-kappa)`, every rounded edge satisfies

```math
r_d(p)+u(p')\ge lambda+u(p)-L delta.                     \tag{26.15}
```

Every cycle therefore has mean at least `lambda-L delta`. Following a
minimizing symbol from every net point gives the reverse inequality on its
eventual cycle, proving (26.13). `square`

### Proposition 26.3 (why the local transfer mechanism stops at a dense cut)

For every `n by n` sign matrix `B`,

```math
\max_(x,y in {+-1}^n)x^TBy\ge {n^{3/2}\over\sqrt3}.      \tag{26.16}
```

Indeed, for uniform `x`, every coordinate `S_j=(B^Tx)_j` has
`ES_j^2=n` and `ES_j^4<=3n^2`. Interpolation gives
`E|S_j|>=sqrt(n/3)`. Some `x` therefore has
`||B^Tx||_1>=n^(3/2)/sqrt3`, and optimizing `y` proves the claim.

Thus the fixed-width lower-spectral mechanism glues across a bounded
interface, while a balanced dense quadratic split has a provably leading
interface response and a standard transfer dimension `2^n`. This falsifies
the naive local-transfer import, not every possible nonlocal multiplicative
norm or quotient.

## 27. Microcanonical hypograph compactness

Let `K` be a compact metric descriptor space and let bounded-above usc
profiles take values in `[-infinity,M]`. Write `f_n ->^h f` for hypograph
convergence: every convergent input sequence obeys the limsup inequality and
every limit point has a recovery sequence. Closed downward hypographs in the
compactified cylinder `K times [-infinity,M]` form a compact hyperspace.

### Theorem 27.1 (compact sup-convolution and finite count recovery)

Let `K_1,K_2,K` be compact, `m:K_1 times K_2->K` continuous, and
`f_n ->^h f`, `g_n ->^h g`. Define

```math
(f\star_mg)(z)=\max_{m(x,y)=z}\{f(x)+g(y)\},            \tag{27.1}
```

with value `-infinity` on an empty fibre. Then

```math
f_n\star_mg_n\xrightarrow{h}f\star_mg.                 \tag{27.2}
```

Every continuous tilt `V` has convergent maximum value; cluster points of
finite maximizers are limit maximizers, and every finite-valued limit
maximizer has an asymptotically maximizing recovery sequence.

Suppose now `A_n,B_n` are finite multiplicity profiles at speed `a_n`, their
normalized log profiles have the preceding limits, and

```math
C_n(z)=\sum_{m(x,y)=z}A_n(x)B_n(y).                     \tag{27.3}
```

If `D_n(z)` is the number of positive decompositions and

```math
\sup_z\log\max(1,D_n(z))=o(a_n),                        \tag{27.4}
```

then the normalized log-count profile of `C_n` has limit (27.1). If the
number of occupied descriptor fibres is also `exp(o(a_n))`, a recovered
fibre lying `Delta` below the maximum has uniform probability
`exp(-a_n Delta+o(a_n))`.

#### Proof

For the upper hypograph bound, take maximizing decompositions and pass to a
compact subsequence; continuity of `m` and the two limsup inequalities give
(27.2). For recovery, take a maximizing limit decomposition and combine the
two component recovery sequences. Adding continuous `V` preserves these
arguments.

On the positive effective domain, largest-summand bounds give

```math
0\le {\log C_n(z)\over a_n}
-\max_{m(x,y)=z}\left({\log A_n(x)\over a_n}
                      +{\log B_n(y)\over a_n}\right)
\le {\log\max(1,D_n(z))\over a_n}.                      \tag{27.5}
```

Both sides are `-infinity` off that domain. Uniformity in (27.4) transfers
hypograph convergence. The same largest-fibre estimate over a
subexponential descriptor image proves the probability statement. `square`

This is the compact sign-dual form of standard epi/Gamma stability under
infimal convolution; see
[Rockafellar--Wets, *Variational Analysis*, Section 7](https://sites.math.washington.edu/~rtr/papers/rtr169-VarAnalysis-RockWets.pdf).
The useful project-specific point is its query and
speed boundary. A full usc microcanonical hypograph retains nonconcave
finite-rate branches that a bounded family of linear-temperature pressures
can convexify away.

### Proposition 27.2 (bounded-temperature pressure can miss a rare maximum)

Fix `B>0` and `0<delta<1/B`. Let one landscape have `ceil(e^n)` states at
energy density zero and one state at density `delta`; let another have only
the bulk. Uniformly for `|beta|<=B`, their normalized log partition functions
differ by `o(1)`, while their normalized maxima differ by `delta`.

Indeed the difference is

```math
{1\over n}\log\left(1+{e^{\beta n\delta}\over\lceil e^n\rceil}\right)
\le {1\over n}\log(1+e^{-n(1-B\delta)+o(n)}).           \tag{27.6}
```

The speed-`n` hypograph sees the isolated energy branch. It still identifies
one maximal state with `e^(sqrt(n))` maximal states, so it does not determine
extremal spacings, Cox decorations, or other subexponential structure.

The deterministic mean-field Blume--Emery--Griffiths model is a classical
benchmark: its occupation entropy is `-sum_j L_j log L_j`, its energy is
`L_++L_- -K(L_+-L_-)^2`, and its nonconcave microcanonical branch is lost by
canonical Legendre data
([Ellis--Touchette--Turkington](https://doi.org/10.1016/j.physa.2003.11.028)).
The present theorem is a speed-sensitive contextual
roof under continuous descriptor tilts, not a wholly separate algebra and
not a structured recovery theorem for sign quadratics.

### Theorem 27.3 (finite lexicographic count algebra)

Fix `L` scales

```math
a_(1,n)>>a_(2,n)>>cdots>>a_(L,n)->infinity,
\qquad a_(j+1,n)/a_(j,n)->0.                           \tag{27.7}
```

For a nonnegative sequence `w_n`, write `nu(w)=u in R^L` when

```math
\log w_n=\sum_(j=1)^L a_(j,n)u_j+o(a_(L,n)),          \tag{27.8}
```

and use `-infinity` for an eventually zero sequence. On
`T_L=R^L union {-infinity}`, let addition be lexicographic maximum and
multiplication be coordinatewise addition. Then, whenever the valuations
exist,

```math
nu(uv)=nu(u)+nu(v),\qquad nu(u+v)=max_lex{nu(u),nu(v)}. \tag{27.9}
```

Consequently, for fixed finite descriptor alphabets and a fixed map
`m:K_1 times K_2->K`, the convolution

```math
C_n(z)=\sum_(m(x,y)=z)A_n(x)B_n(y)                    \tag{27.10}
```

has exact asymptotic profile

```math
nu(C_n(z))=
\max_(lex,m(x,y)=z){nu(A_n(x))+nu(B_n(y))}.            \tag{27.11}
```

Every finite `T_L` profile has an all-order integer realization up to one
common additive shift in its leading coordinate. For growing fibres, the
same rule holds uniformly if the term remainders are uniform and
`log |I_n|=o(a_(L,n))`. This branching condition is worst-case sharp: if
`exp(ca_(L,n))` terms tie, the last coordinate increases by `c`.

#### Proof

Products add expansions. If `u>_lex v` and `j` is their first unequal
coordinate, then

```math
\log w_n-\log z_n
=a_(j,n)(u_j-v_j)+o(a_(j,n))->+infinity,
```

so the larger term determines the sum; equal valuations cost only `log 2`.
Finite iteration proves (27.11). For recovery, add a constant `C` making
every leading coefficient positive and set

```math
A_n(q)=floor\exp\left(a_(1,n)(C+u_1(q))
             +\sum_(j=2)^L a_(j,n)u_j(q)\right).      \tag{27.12}
```

Flooring costs `o(1)`. Finally

```math
0\le\log\sum_(i in I_n)w_(n,i)-\max_i\log w_(n,i)
\le\log|I_n|,                                         \tag{27.13}
```

which proves the growing-fibre statement and its equal-term sharpness.
`square`

### Proposition 27.4 (pointwise multiscale profiles miss saddle mass)

At scales `(n,log n)`, fix `p in (0,1)` and an integer subsequence with `pn`
integral. Stirling's formula gives

```math
nu {n choose pn}=(h(p),-1/2).                          \tag{27.14}
```

In the Vandermonde convolution

```math
{2n choose 2pn}=\sum_k {n choose k}{n choose 2pn-k},  \tag{27.15}
```

the pointwise largest summand has vector `(2h(p),-1)`, while the left side
has `(2h(p),-1/2)`. The missing half-coordinate is the mass of
`Theta(sqrt(n))` near-saddle decompositions. Thus a bare pointwise finite-
speed roof does not compose once the decomposition fibre is large at the
smallest retained speed; a tangent-density or equivalent decoration is
required even under a smooth strictly concave leading profile.

## 28. Regular-Hadamard amplification limits

### Theorem 28.1 (monotone Boolean amplification)

Let `H` be a symmetric Hadamard matrix of order `h` and suppose that some
`u in {+-1}^h` satisfies `Hu=sqrt(h)u`. For a fixed real symmetric
`d by d` matrix `B`, put

```math
B_r=B tensor H^(tensor r),\qquad N_r=dh^r,
```

and define

```math
q_r^+(B)={1\over2N_r^(3/2)}\max_(x in {+-1}^N_r)x^TB_rx,
\qquad
q_r^abs(B)={1\over2N_r^(3/2)}\max_x|x^TB_rx|.          \tag{28.1}
```

Both sequences are nondecreasing and converge, with

```math
q_r^+(B),q_r^abs(B)
\le {||B||_(2->2)\over2sqrt(d)}.                       \tag{28.2}
```

#### Proof

The Boolean lift `x mapsto x tensor u` obeys

```math
(x tensor u)^T(B_r tensor H)(x tensor u)
=h^(3/2)x^TB_rx.                                      \tag{28.3}
```

The denominator is multiplied by the same factor, so every signed normalized
value at level `r` reappears at level `r+1`. This proves both monotonicities.
Moreover `||B_r||=||B||h^(r/2)` and `||x||_2^2=dh^r`, which give (28.2).
Bounded monotone sequences converge. The positive eigenvalue is needed for
`q^+`; eigenvalue `-sqrt(h)` still gives monotonicity of the absolute
functional. `square`

### Theorem 28.2 (a finite-dimensional limiting response set)

Write `s_r=h^r` and split a Boolean vector into `d` blocks. Let

```math
K_r^(d)=conv\left\{
 \left(s_r^(-3/2)x_i^TH^(tensor r)x_j\right)_(i<=j):
 x_i in {+-1}^(s_r)
\right\}.                                             \tag{28.4}
```

Then

```math
K_r^(d) subseteq K_(r+1)^(d) subseteq[-1,1]^(d(d+1)/2) \tag{28.5}
```

and the sets converge in Hausdorff distance to

```math
K_infinity^(d)=closure(union_r K_r^(d)).               \tag{28.6}
```

For every symmetric `B`,

```math
lim_r q_r^+(B)
={1\over2d^(3/2)}\max_(K in K_infinity^(d))<B,K>,      \tag{28.7}
```

with the analogous absolute support formula. Convergence is uniform on
bounded coefficient sets. An external entrywise-`epsilon` net has at most

```math
(1+2/epsilon)^(d(d+1)/2)                              \tag{28.8}
```

points, independently of amplification depth.

#### Proof

Applying `x_i mapsto x_i tensor u` to every block preserves each normalized
entry in (28.4), proving nesting. Increasing compact subsets of a fixed
compact cube converge in Hausdorff distance to the closure of their union.
Expanding the quadratic gives

```math
{x^T(B tensor H^(tensor r))x\over2(ds_r)^(3/2)}
={1\over2d^(3/2)}\sum_(i,j)B_(ij)
 {x_i^TH^(tensor r)x_j\over s_r^(3/2)}.                \tag{28.9}
```

Linear optimization is unchanged by convexification. Hausdorff convergence
gives (28.7), its absolute analogue, and uniform convergence on bounded dual
sets. A grid in the ambient cube gives (28.8). `square`

### Corollary 28.3 (a dense structured signing hierarchy has a limit)

The order-four Walsh matrix has a Boolean vector `u=(1,1,1,-1)` with
`W_4u=2u`. Hence every fixed linear-label Walsh graph program whose labels
are extended only by zero coordinates and whose scalar onsite and bridge
coefficients remain fixed has convergent normalized upper and absolute
Boolean maxima. After regrouping coordinates,

```math
M_(m_0+r)=M_(m_0) tensor W_4^(tensor r).               \tag{28.10}
```

If a fixed outer template `B` is symmetric and entrywise `+-1`, hollowing
`B tensor W_4^(tensor r)` gives a valid dense signing `A_r` of order `d4^r`.
Diagonal removal changes the normalized absolute quadratic by at most
`1/(2sqrt(d4^r))` (and changes it exactly by zero when `tr(B)=0`). Thus

```math
{Q(A_r)\over(d4^r)^(3/2)}
```

converges. This is an exact tensor-hierarchy theorem on geometric orders. It
does not transfer near-minimizers to arbitrary orders and therefore does not
imply convergence of the motivating sequence.

## 29. Summable Boolean recovery

For nonempty compact sets, write

```math
e(K,L)=\sup_(x in K)dist(x,L)
```

for directed Hausdorff excess.

### Theorem 29.1 (one-sided recovery forces a response limit)

Let `K_r` be nonempty compact subsets of one compact metric space. If

```math
e(K_r,K_(r+1))\le epsilon_r,
\qquad \sum_r epsilon_r<infinity,                     \tag{29.1}
```

then `K_r` converges in Hausdorff distance to a compact set `K_infinity`.
Moreover

```math
e(K_r,K_infinity)\le\sum_(s>=r)epsilon_s.             \tag{29.2}
```

No reverse bound tending to zero can depend only on this tail: exact forward
inclusions may introduce new points arbitrarily late.

#### Proof

Starting from each `x_r in K_r`, choose a forward chain with successive
distances at most `epsilon_s`. Its limit lies within the tail in (29.2).
Let `E` be the set of all such limits and set `K_infinity=closure(E)`. This
proves (29.2). A finite net of `K_infinity` with centers in `E` has all of
its forward chains represented in every sufficiently late `K_r`, proving
the reverse excess tends to zero. `square`

Now let

```math
kappa_r(x,y)={x^TA_ry\over a_r},\qquad x,y in {+-1}^(n_r),
```

be uniformly bounded, and let `K_r^(d)` be the convex image of all fixed-`d`
cross-response matrices `(kappa_r(x_i,x_j))_(i<=j)`. A common Boolean map
`L_r:{+-1}^(n_r)->{+-1}^(n_(r+1))` has pair distortion

```math
alpha_r=\sup_(x,y)|kappa_(r+1)(L_rx,L_ry)-kappa_r(x,y)|. \tag{29.3}
```

Applying the same map to every port gives

```math
e_infinity(K_r^(d),K_(r+1)^(d))\le alpha_r.           \tag{29.4}
```

Hence `sum alpha_r<infinity` makes every fixed-port carrier converge and
makes all bounded linear and absolute support queries converge uniformly.
A backward Boolean map with distortion `beta_r` strengthens this to

```math
d_H(K_r^(d),K_infinity^(d))
\le\sum_(s>=r)max(alpha_s,beta_s).                    \tag{29.5}
```

### Theorem 29.2 (a checkable matrix and amplification certificate)

If `T_r in {0,+-1}^(n_(r+1) times n_r)` has one nonzero entry per row, then
`L_rx=T_rx` is Boolean and its exact distortion is

```math
alpha_r=left\|T_r^T{A_(r+1)\over a_(r+1)}T_r
                 -{A_r\over a_r}\right\|_(infinity->1). \tag{29.6}
```

Entrywise `l_1` and `n_r||.||_(2->2)` bound this norm. In particular, take
`a_r=n_r^(3/2)` and suppose

```math
n_(r+1)=h_rn_r,
\qquad A_(r+1)=A_r tensor H_r+E_r.                    \tag{29.7}
```

For Boolean `u_r`, put

```math
rho_r={u_r^TH_ru_r\over h_r^(3/2)},\quad
M_r={||A_r||_(2->2)\over sqrt(n_r)},\quad
e_r={||E_r||_(2->2)\over sqrt(n_(r+1))}.              \tag{29.8}
```

The replication `x mapsto x tensor u_r` satisfies

```math
alpha_r\le M_r|rho_r-1|+e_r.                         \tag{29.9}
```

Thus bounded `M_r` and summability of the right side force convergence of
every fixed-port response carrier and outer quadratic support query. A
checkable sufficient bound for `M_r` is

```math
||H_r||<= (1+sigma_r)sqrt(h_r),
\qquad \sum_r(sigma_r+e_r)<infinity.                 \tag{29.10}
```

#### Proof

Compressing (29.7) by `T_r=I tensor u_r` gives

```math
T_r^T{A_(r+1)\over n_(r+1)^(3/2)}T_r
-{A_r\over n_r^(3/2)}
=(rho_r-1){A_r\over n_r^(3/2)}
 +{T_r^TE_rT_r\over n_(r+1)^(3/2)}.                  \tag{29.11}
```

The two Boolean bilinear norms are at most `M_r|rho_r-1|` and `e_r`.
Also `M_(r+1)<=(1+sigma_r)M_r+e_r`, proving the uniform bound and the
claim by Theorem 29.1. `square`

### Corollary 29.3 (a non-tensor dense-sign hierarchy)

Fix a nontrivial symmetric regular Hadamard `H` with a positive Boolean
eigenvector. Starting from any symmetric full sign matrix `C_0`, form
`C_r tensor H` and flip every sign on a fixed-point-free perfect matching to
obtain `C_(r+1)`. The perturbation is a direct sum of signed
`[[0,2],[2,0]]` blocks, so its operator norm is exactly two and

```math
e_r={2\over sqrt(N_(r+1))}.
```

These errors are summable. Therefore the fixed-port carriers converge, and
the hollow dense signings `C_r^circ` satisfy

```math
{Q(C_r^circ)\over N_r^(3/2)}\longrightarrow q
```

for some `q`. This changes `Theta(N_r)` undirected edges at every level and
is genuinely non-tensor.

The finite-total-drift hypothesis has a sharp abstract boundary. There are
growing Boolean kernels with distortions `alpha_r->0` and
`sum alpha_r^2<infinity` whose scalar carriers oscillate: take
`n_r=2^r`, `c_r=sin(log(r+2))`, and
`A_r=c_r sqrt(n_r)I`. Hence vanishing or square-summable error alone cannot
replace (29.1) without an additional cancellation mechanism.

## 30. Automatic tensor-prefix phase laws

### Theorem 30.1 (continuous phase law for regular-Hadamard prefixes)

Let `H` be a symmetric Hadamard matrix of order `h>1`, with top-left entry
one and a Boolean vector `u` satisfying `Hu=sqrt(h)u`.  Use the lexicographic
ordering

```math
H_r=H^(tensor r),\qquad H_(r+1)=H tensor H_r.          \tag{30.1}
```

The `H_r` are compatible leading principal blocks of one infinite symmetric
sign matrix.  Let `A_n` be its hollow leading `n by n` block, put `R_r=h^r`,
and define

```math
F_r(t)={Q(A_(floor(tR_r)))\over R_r^(3/2)},
\qquad 1\le t\le h.                                   \tag{30.2}
```

Then `F_r` converges uniformly to a continuous nondecreasing function `F`.
Moreover

```math
G_r(t)={Q(A_(floor(th^r)))\over floor(th^r)^(3/2)}
\longrightarrow L(t):={F(t)\over t^(3/2)}              \tag{30.3}
```

uniformly on `[1,h]`.  If

```math
r(n)=floor(log_h n),\qquad t_n={n\over h^(r(n))},       \tag{30.4}
```

then

```math
{Q(A_n)\over n^(3/2)}-L(t_n)\longrightarrow0.          \tag{30.5}
```

Thus this coherent all-order sequence converges precisely when its continuous
scale-phase profile `L` is constant.  Always `L(1)=L(h)=1/2`.

#### Proof

For a base-`h` rational `t=p/h^k`, with integer
`h^k<=p<=h^(k+1)`, the relevant prefix factors exactly as

```math
C_(p h^(r-k))=B_(p,k) tensor H_(r-k),                 \tag{30.6}
```

where `B_(p,k)` is the leading `p`-square block of `H_(k+1)`.  The outer
template is fixed.  The Boolean lift `x mapsto x tensor u` and Theorem 28.1,
with the vanishing diagonal-removal error, make `F_r(t)` converge on this
dense set.

For `1<=t<=s<=h`, write `n=floor(tR_r)`, `m=floor(sR_r)`, and `d=m-n`.
The block decomposition of `A_m` over `A_n` has cross block `C` and new
diagonal block `D` satisfying

```math
||C||<=sqrt(hR_r),\qquad ||D||<=sqrt(hR_r)+1.          \tag{30.7}
```

Indeed both unhollowed blocks are coordinate compressions of `H_(r+1)`.
For Boolean block spins, their new half-quadratic energy is at most

```math
sqrt(hR_rnd)+{d\over2}(sqrt(hR_r)+1).                  \tag{30.8}
```

Principal deletion gives `Q(A_n)<=Q(A_m)`: average over independent unbiased
missing spins and select an extension with the desired sign.  Hence

```math
0\le F_r(s)-F_r(t)
\le sqrt(h)sqrt{(n/R_r)(d/R_r)}
 +{sqrt(h)\over2}{d\over R_r}+{d\over2R_r^(3/2)}.     \tag{30.9}
```

This is a common `O_h(sqrt(s-t)+s-t)+o_r(1)` modulus.  Convergence on a
finite base-`h` rational net makes the functions uniformly Cauchy; the
retained modulus makes the limit continuous.  Dividing by the uniformly
positive factor `(floor(th^r)/h^r)^(3/2)` proves (30.3).  Choosing
`t=t_n` gives (30.5) with no floor error.  At geometric orders, the spectral
upper bound and the Boolean eigenvector give normalized value
`1/2+O(h^(-r/2))`, proving the endpoint identities. `square`

### Corollary 30.2 (one Walsh hierarchy genuinely does not converge)

For

```math
H=\begin{pmatrix}
1&1&1&1\\1&-1&1&-1\\1&1&-1&-1\\1&-1&-1&1
\end{pmatrix},
\qquad u=(1,1,1,-1)^T,                                 \tag{30.10}
```

let `A_n` be the coherent hollow prefix signing from Theorem 30.1.  Then

```math
{Q(A_(4^r))\over(4^r)^(3/2)}={1\over2},               \tag{30.11}
```

whereas

```math
{Q(A_(3*4^r))\over(3*4^r)^(3/2)}
\ge {89\over96sqrt3}=0.535251812061\ldots\quad(r>=2). \tag{30.12}
```

Consequently the one explicit all-order sequence `Q(A_n)/n^(3/2)` does not
converge, and its complete asymptotic behavior is the nonconstant continuous
mantissa law `L(t_n)` from Theorem 30.1.

#### Proof

At order `4^r`, `H_r` has operator norm `2^r`, Boolean eigenvector
`u^(tensor r)`, and trace zero, giving (30.11).  Let `B` be the leading
three-square block of `H`.  At order 48, the explicit Boolean vector

```text
(+ + + -  + + - +  + - + +  - + + -)
(+ - - -  + + + -  + - + +  + + - -)
(+ + + +  - + - +  - + + -  - - + -)
```

has quadratic value `356` against `B tensor H_2`.  Tensoring it with
`u^(tensor s)` multiplies the value by `8^s`; trace zero makes hollowing
exact.  Division by `2(48)^(3/2)8^s` gives (30.12).  Theorem 30.1 supplies
the full phase law. `square`

This is a near-original scalable obstruction, but not nonconvergence of the
minimizing values `M_n`: the displayed prefixes are one fixed family and are
not known to be near-minimizers.  What it disproves is the general inference
that convergence on every fixed amplification phase, even with exact
all-order realization, synchronizes the phases automatically.  Detailed
proofs, audits, and the exact certificate are in
[`drafts/automatic_tensor_prefix_phase.md`](drafts/automatic_tensor_prefix_phase.md),
[`drafts/walsh_prefix_nonconvergence.md`](drafts/walsh_prefix_nonconvergence.md),
and their adjacent audit files.

### Theorem 30.3 (scale-phase averaging law and its information boundary)

Let `h>=2` be an integer, let `(q_n)` be bounded, and suppose a continuous
`L:[1,h]->R` satisfies

```math
sup_(h^r<=n<h^(r+1))|q_n-L(n/h^r)|->0.               \tag{30.13}
```

Then every continuous `psi` has the unique logarithmic empirical limit

```math
{1\over\log N}\sum_(n<=N){\psi(q_n)\over n}
->{1\over\log h}\int_1^h{\psi(L(t))\over t}dt.       \tag{30.14}
```

Positive power biases behave differently.  For `alpha>0`, `s in [1,h]`,
and `N_R=floor(sh^R)`,

```math
{\sum_(n<=N_R)n^(alpha-1)q_n\over
 \sum_(n<=N_R)n^(alpha-1)}
->C_alpha(s),                                        \tag{30.15}
```

where

```math
C_alpha(s)={alpha\over s^alpha}\left{
 {1\over h^alpha-1}\int_1^h t^(alpha-1)L(t)dt
 +\int_1^s t^(alpha-1)L(t)dt\right}.                \tag{30.16}
```

The map `L mapsto C_alpha` is injective, with

```math
L(s)=C_alpha(s)+{s\over alpha}C_alpha'(s).            \tag{30.17}
```

Thus `C_alpha` is constant exactly when `L` is constant.  Logarithmic
averaging retains only the pushforward of `dt/(t log h)` by `L`: for
`u=log_h t`, the distinct profiles `cos(2pi u)` and `cos(4pi u)` have the
same logarithmic law, whereas every positive-power response phase
distinguishes them.

#### Proof

On a complete scale block, harmonic and power-weighted Riemann sums give,
respectively,

```math
\sum_(n=h^r)^(h^(r+1)-1){\psi(q_n)\over n}
->\int_1^h{\psi(L(t))\over t}dt,                     \tag{30.18}
```

and

```math
\sum_(n=h^r)^(h^(r+1)-1)n^(alpha-1)q_n
=h^(alpha r)\int_1^h t^(alpha-1)L(t)dt+o(h^(alpha r)).\tag{30.19}
```

Cesaro averaging the `log_hN+O(1)` bounded harmonic blocks proves (30.14).
Geometric summation of (30.19), adding the final partial block, and using
`sum_(n<=N)n^(alpha-1)=N^alpha/alpha+o(N^alpha)` prove
(30.15)--(30.16).  Differentiation yields
`C_alpha'=(alpha/s)(L-C_alpha)`, hence (30.17).  Finally `dt/(t log h)=du`,
and both displayed cosines have the arcsine pushforward under uniform `u`.
`square`

For Corollary 30.2, the logarithmic response law exists and is nontrivial,
while every positive power-weighted mean retains a nonconstant subsequential
phase.  Averaging the scale phase is therefore a declared-query quotient,
not phase synchronization.  The proof and numerical wind tunnel are in
[`drafts/phase_averaging_laws.md`](drafts/phase_averaging_laws.md).

## 31. Phase refresh and Boolean pullback synchronization

### Theorem 31.1 (refresh dominates transfer defect)

Let `X` be a compact metric phase space.  At level `r` and phase `x`, let
`T_(r,x)` be a self-adjoint operator on a finite probability space and put

```math
phi_r(x)=sup_(||f||_infinity<=1)|<f,T_(r,x)f>|.        \tag{31.1}
```

Assume `phi_r` converges uniformly to a continuous `phi`.  Suppose there is
a Borel phase kernel `P_r` and, over every branch `x->y`, an `L^2` contraction
`U` which is also an `L^infinity` contraction, such that the strongly
measurable pulled-back operators obey

```math
\left\|T_(r,x)-\int U^*T_(r+1,y)U\,dGamma_(r,x)(y,U)
\right\|_(2->2)\le epsilon_r,                         \tag{31.2}
```

where `Gamma` has phase marginal `P_r`.  For delay windows `ell_r`, write

```math
K_r=P_r...P_(r+ell_r-1),\qquad
E_r=sum_(j=r)^(r+ell_r-1)epsilon_j.                   \tag{31.3}
```

If one full-support probability measure `nu` and numbers `alpha_r>0` satisfy

```math
K_r(x,B)>=alpha_rnu(B)                                \tag{31.4}
```

for every `x,B`, and, with `omega_r=||phi_r-phi||_infinity`,

```math
{omega_r+omega_(r+ell_r)+E_r\over alpha_r}->0,        \tag{31.5}
```

then `phi` is constant.  Quantitatively,

```math
max_Xphi-\int phi\,dnu
\le liminf_r{omega_r+omega_(r+ell_r)+E_r\over alpha_r}. \tag{31.6}
```

If `X` is finite and `mu=min_xnu(x)>0`, the oscillation of `phi` is at most
the right side of (31.6) divided by `mu`.

#### Proof

For every Boolean-bounded `f`, (31.2) and the two contraction properties give

```math
|<f,T_(r,x)f>|
<=\int |<Uf,T_(r+1,y)Uf>|dGamma+epsilon_r.
```

Taking the supremum and iterating yields

```math
phi_r<=K_rphi_(r+ell_r)+E_r.                          \tag{31.7}
```

At a maximizer `x_*` of `phi`, decompose
`K_r(x_*,.)=alpha_rnu+(1-alpha_r)rho_r`.  Uniform recovery and (31.7) give

```math
max phi-omega_r
<=alpha_r\int phi dnu+(1-alpha_r)max phi
  +omega_(r+ell_r)+E_r.
```

This is (31.6).  Under (31.5), the full-support average equals the maximum;
continuity forces `phi` to be constant.  In the finite case, the average
deficit is at least `mu osc(phi)`. `square`

Equal-fibre signed coordinate replications and reorderings are Boolean
contractions.  Thus a finite convex combination of their operator pullbacks,
together with a Doeblin refresh window, is a checkable phase-collapse
certificate.  It compares neither maximizing spins nor full response
carriers; its description size is nevertheless model-dependent because a
raw operator certificate may still have quadratically many coefficients.

The error/refresh ratio is scale-sharp.  On `X={0,1}`, take scalar operators
zero and one, and

```math
P_r(x,.)=(1-2^(-r))delta_x+2^(-r)nu.                  \tag{31.8}
```

Identity pullbacks have defect `2^(-r-1)` and every one-step kernel has full
support, yet the response profile stays `(0,1)`.  Thus full support and
vanishing defect without (31.5) do not synchronize phases.

### Corollary 31.2 (the Walsh phase prices every refresh certificate)

For the order-four Walsh prefix profile in Corollary 30.2, use the operator
normalization `Phi=2L`.  Let

```math
nu_*={99\over200}(delta_1+delta_4)+{1\over100}lambda, \tag{31.9}
```

where `lambda` is normalized Lebesgue measure on `[1,4]`.  Every Boolean-
pullback refresh satisfying `K_r(t,.)>=alpha_rnu_*` must obey

```math
liminf_r{omega_r+omega_(r+ell_r)+E_r\over alpha_r}
\ge {89\over48sqrt3}-1.01
=0.06050362412\ldots.                                 \tag{31.10}
```

Indeed `Phi(1)=Phi(4)=1`, `Phi(3)>=89/(48sqrt3)`, and the coordinate-
compression bound gives `Phi<=2`; hence `int Phi dnu_*<=1.01`, and (31.6)
applies.  This identifies the quantitative price any proposed balanced
reordering of the bad Walsh hierarchy must pay.  The independently audited
full statement is in
[`drafts/phase_refresh_synchronization.md`](drafts/phase_refresh_synchronization.md).

### Theorem 31.3 (semantic expander refresh forces toll or memory)

Let finite-state Markov kernels `P_j` preserve one full-support law `pi` and
satisfy

```math
||P_j-Pi||_(L^2(pi)->L^2(pi))<=rho<1,                \tag{31.11}
```

where `Pi f=int f dpi`.  Let `g:X->[0,B]`, and suppose

```math
||f_j-g||_infinity<=omega_j,
\qquad f_j<=P_jf_(j+1)+epsilon_j,                    \tag{31.12}
```

with nonnegative errors.  If
`delta_j=epsilon_j+omega_j+omega_(j+1)`, then for every `x` and every
window of length `t`,

```math
g(x)-int g dpi
<= {B rho^t\over sqrt(pi(x))}+sum_(j=r)^(r+t-1)delta_j.\tag{31.13}
```

In particular, suppose `D=g(x)-int g dpi>0`, `|X|=S`,
`pi(x)>=kappa/S`, `0<rho<1`, and `delta_j<=delta<D/2`.  Then

```math
log S>=log kappa+{log(1/rho)D\over delta}
       -2log(1/rho)-2log(2B/D).                      \tag{31.14}
```

Thus a bounded-state uniformly scrambling semantic quotient pays a fixed
transfer toll; at toll `delta`, its description needs
`Omega(log(1/rho)D/delta)` bits.

#### Proof

The two inequalities in (31.12) give `g<=P_jg+delta_j`.  Iterate positive
kernels.  Their product contracts the mean-zero subspace by `rho^t`.
Point evaluation at `x` has centered `L^2(pi)` norm
`sqrt(1/pi(x)-1)`, proving (31.13).  Choose

```math
t=ceil{log(2B sqrt(S/kappa)/D)\over log(1/rho)}       \tag{31.15}
```

so the mixing term is at most `D/2`; then `delta>=D/(2t)`, which rearranges
to (31.14). `square`

The semantic qualifier is essential.  If
`T_(r,g)=V_(r,g)^*S_rV_(r,g)` is merely a finite signed-permutation gauge
orbit, the branch maps

```math
U_(r,g,h)=V_(r+1,h)^*L_rV_(r,g)                     \tag{31.16}
```

have a pullback independent of `h`, so **any** phase kernel is exactly
realizable while all Boolean responses are already identical.  For the
order-four Walsh matrix, one such one-bit orbit has fixed-coordinate
operator diameter exactly two.  Operator diameter alone therefore does not
measure semantic phase memory.

For the nonconstant Walsh prefix response `Phi=2L`, put mass `99/200` at
each endpoint and the remaining `1/100` on a finite phase sample containing
`3`.  Then

```math
D>=D_*={89\over48sqrt3}-1.01=0.06050362412... .      \tag{31.17}
```

A half-scrambling certificate with `delta<=C/sqrt(N)` must consequently
have

```math
log_2S>={D_*\over C}sqrt(N)-O(1).                    \tag{31.18}
```

Hence bounded, polynomial, and `exp(o(sqrt N))` stationary expander phase
quotients cannot synchronize that semantic phase at the natural
`N^(-1/2)` transfer scale.  Nonstationary laws, vanishing scrambling, or a
non-pullback mechanism remain outside the theorem.  The independently
audited proof is in
[`drafts/expander_phase_refresh_complexity.md`](drafts/expander_phase_refresh_complexity.md).

## 32. Tangent mass and a finite Gaussian response semigroup

The lattice-Laplace estimate below is classical.  Its role here is to repair
the precise response coordinate missed by the pointwise multiscale roof in
Theorem 27.3; the finite-parameter closure in Theorem 32.2 is the new theory
benchmark.

### Theorem 32.1 (Morse tangent-mass convolution)

Let positive arrays on `nD_A intersect Z^d` and `nD_B intersect Z^d` obey

```math
A_n(k)=Theta(n^alpha e^(nf(k/n))),\qquad
B_n(k)=Theta(n^beta e^(ng(k/n)))                     \tag{32.1}
```

uniformly on the full feasible fibre.  At an admissible grid query `z`, let

```math
F_z(x)=f(x)+g(z-x)
```

have a unique interior maximizer `x_z`, with uniform two-sided quadratic
exposure there and a strict global gap.  Then

```math
sum_kA_n(k)B_n(nz-k)
=Theta\left(n^(alpha+beta+d/2)e^(nF_z(x_z))\right).    \tag{32.2}
```

If the input amplitudes are locally

```math
A_n(k)=n^alpha e^(nf(k/n))(a(k/n)+o(1))
```

and similarly for `B_n`, and
`J_z=-D^2F_z(x_z)` is positive definite, then on the unit-density standard
lattice the leading amplitude is

```math
{(2pi)^(d/2)a(x_z)b(z-x_z)\over sqrt(det J_z)}.        \tag{32.3}
```

A tangent lattice of covolume `v` adds the factor `v^(-1)`.  The same result
holds when the precise input asymptotics are local, provided a polynomial
global envelope and an exponential off-saddle gap are available.

#### Proof

After extracting the prefactors and `e^(nF_z(x_z))`, the upper quadratic
bound gives a shifted Gaussian lattice sum `O(n^(d/2))`.  A radius-`sqrt n`
ball contains `Omega(n^(d/2))` feasible lattice points, each with a constant-
factor lower weight.  This proves (32.2).  At `k=nx_z+sqrt(n)y`, Taylor
expansion and the rescaled lattice Riemann sum give

```math
n^(-d/2)sum_k e^(n(F_z(k/n)-F_z(x_z)))
->\int_(R^d)e^(-y^TJ_zy/2)dy,
```

which is (32.3).  The localization variant discards only an exponentially
smaller total mass. `square`

This recovers the missing `+1/2` logarithmic exponent in Vandermonde
convolution and, in `q-1` tangent dimensions, the usual multinomial
`-(q-1)/2` exponent.  Boundary summands are handled by the global entropy
upper bound, not by falsely extending interior Stirling asymptotics to the
simplex boundary.

### Theorem 32.2 (finite-parameter Gaussian tangent semigroup)

For positive-definite `P`, define on `Z^d`

```math
G_n^(c,mu,P,alpha,a)(k)
=n^alpha a\exp\{nc-(k-nmu)^TP(k-nmu)/(2n)\}.          \tag{32.4}
```

Write `Sigma=P^(-1)` and

```math
m=a(2pi)^(d/2)/sqrt(det P).                           \tag{32.5}
```

On compact parameter sets with precision eigenvalues bounded away from zero
and infinity, convolution is uniformly represented by the associative,
commutative parameter law

```math
(c,mu,Sigma,alpha,m) star (c',mu',Sigma',alpha',m')
=\left(c+c',mu+mu',Sigma+Sigma',
       alpha+alpha'+d/2,mm'\right).                   \tag{32.6}
```

More precisely, at compact scaled queries,

```math
(G_n^theta*G_n^phi)(ell)
=(1+o(1))G_n^(theta star phi)(ell).                   \tag{32.7}
```

This is a strict finite-dimensional composable state of dimension
`d(d+1)/2+d+3`, independent of the lattice support size.  It has finite
integer realizations on every large order: truncate to one common linear-size
box, multiply by `e^(nC)`, and floor.  The recovered convolution still obeys
(32.7), with `C` added to each leading `c` coordinate.

#### Proof

Completing the square gives output covariance `Sigma+Sigma'`, mean
`mu+mu'`, and saddle Hessian `P+P'`.  Theorem 32.1 gives raw amplitude

```math
a_out={(2pi)^(d/2)aa'\over sqrt(det(P+P'))}.           \tag{32.8}
```

Since

```math
det(P+P')=det(P)det(P')det(P^(-1)+P'^(-1)),            \tag{32.9}
```

(32.8) is exactly `m_out=mm'`; all other coordinates update additively.
This proves (32.6)--(32.7) and associativity.

For integer recovery, compactness places every relevant saddle a fixed
scaled distance inside common boxes.  Choose `C` so both factors exceed
`e^(delta n)` on one saddle neighbourhood.  There

```math
0<=XY-floor(X)floor(Y)<=X+Y,
```

so relative rounding loss is exponentially small.  Outside the
neighbourhood, including truncated tails, the whole unrounded product mass
is exponentially smaller by the quadratic gap. `square`

A quartic saddle has tangent mass `Theta(n^(3/4))`, and a flat exposed
interval has `Theta(n)` mass.  Thus `d/2` is not a universal correction and
arbitrary Morse amplitude fields are not finite states.  The theorem is a
genuine closed response algebra only on the declared Gaussian class.  The
complete repaired proof and independent audits are in
[`drafts/morse_tangent_mass_composition.md`](drafts/morse_tangent_mass_composition.md).

### Theorem 32.3 (power roofs close, but their tangent family is Gaussian-rigid)

Fix `p>1` and put `I_(p,a)(x)=a|x|^p`.  The leading roofs form the exact
infimal-convolution semigroup

```math
inf_x\{a|x|^p+b|z-x|^p\}
=\left(a^(-1/(p-1))+b^(-1/(p-1))\right)^(-(p-1))|z|^p.\tag{32.10}
```

For the lattice arrays

```math
A_n^a(k)=exp\{-a|k|^p/n^(p-1)\},                    \tag{32.11}
```

the central convolution has tangent mass `Theta(n^(1-1/p))`, whereas every
fixed nonzero macroscopic output has a nondegenerate saddle and tangent mass
`Theta(n^(1/2))`.  More decisively, the normalized tangent densities

```math
g_(p,a)(x)=Z_(p,a)^(-1)e^(-a|x|^p)                  \tag{32.12}
```

satisfy

```math
g_(p,a)*g_(p,a)=g_(p,b)                              \tag{32.13}
```

for some `b>0` if and only if `p=2`, in which case `b=a/2`.

#### Proof

Strict convexity and homogeneity give (32.10).  At output zero, rescale
`k=n^(1-1/p)u`; at a nonzero output the unique optimal split has positive
finite Hessian, giving the two tangent exponents by lattice Riemann sum and
Laplace's method.

If (32.13) holds and `X,X_1,X_2` have density `g_(p,a)`, then
`X_1+X_2` is a scale copy of `X`.  Finite nonzero variance forces the scale
to be `sqrt2`; iteration makes every normalized dyadic iid sum have the law
of `X`.  The central limit theorem forces `X` to be Gaussian, hence `p=2`.
The Gaussian converse is exact.  Finally, with `s_n=n^(1-1/p)`, the locally
uniform tangent limit

```math
{1\over s_n}(A_n^a*A_n^a)(floor(s_n u))
->\int_R e^(-a|v|^p)e^(-a|u-v|^p)dv                 \tag{32.14}
```

follows by dominated Riemann sums, so the continuous rigidity is inherited
by the discrete carrier. `square`

The leading rate coefficient therefore closes for every power while the
next response object does not: for `p!=2`, repeated composition generates
new convolution shapes.  The Gaussian semigroup in Theorem 32.2 is
exceptional, not the first member of a fixed-power finite hierarchy.  This
does not exclude every non-Gaussian finite semigroup; it closes the proposed
power-exponential extension.  Full details and audit are in
[`drafts/non_gaussian_tangent_closure.md`](drafts/non_gaussian_tangent_closure.md).

## 33. Contracting fibres and persistent reward cohomology

### Theorem 33.1 (three-channel response decomposition)

Let `G=(Q,E)` be a finite directed control graph.  Vertex `q` carries a
finite probability space `(X_q,pi_q)`, and edge `e:q->q'` carries a Markov
operator satisfying

```math
P_e1=1,qquad pi_qP_e=pi_(q'),qquad
||P_eh||_(2,pi_q)<=rho||h||_(2,pi_(q'))              \tag{33.1}
```

for every centred `h`, with one `rho<1`.  Give edge `e` a reward residual
`a_e` and vertex `q` a terminal residual `u_q`.  Write

```math
m_e=pi_qa_e,quad b_e=a_e-m_e1,quad
bar u_q=pi_qu_q,quad v_q=u_q-bar u_q1,              \tag{33.2}
```

and let `B=max_e||b_e||_2`, `R=max_q||v_q||_2`.  For a path
`p=e_1...e_t`, its expected accumulated response residual `D_p` satisfies

```math
\left\|D_p-\left(bar u_(q_t)+\sum_(s=1)^t m_(e_s)\right)1\right\|_2
<=rho^tR+B{1-rho^t\over1-rho}.                       \tag{33.3}
```

The constant and centred parts are orthogonal.  If

```math
chi_G(m)=max_C{1\over|C|}\left|\sum_(e\in C)m_e\right|, \tag{33.4}
```

where `C` ranges over directed simple cycles in the reachable recurrent
graph, then

```math
limsup_(t->infinity){1\over t}
 sup_(|p|=t)||D_p||_2=chi_G(m).                       \tag{33.5}
```

Moreover, `chi_G(m)` is exactly the maximum absolute pairing of `m` with a
normalized stationary edge flow.  It vanishes precisely when `m` is a
vertex gradient on every strongly connected component.  Hence a uniformly
bounded all-depth response approximation exists exactly when the recurrent
scalar reward cochain is potential; strict mixing cannot pay nonzero scalar
holonomy.

For a fixed strongly connected `G`, put `r_G=|E|-|Q|+1`.  Let `N_U(eta)`
and `N_B(eta)` cover the terminal and centred-edge dictionaries in maximum
`L^2` distance, and let scalar cohomology classes lie in the radius-`L` ball
for `chi_G`.  Then the fixed certified carrier, modulo exact endpoint gauges,
has a composable codebook of size

```math
N_U(epsilon)N_B((1-rho)epsilon)
\left(1+{2L\over epsilon}\right)^(r_G)               \tag{33.6}
```

whose error on every length-`t` path is

```math
epsilon t+C_Gepsilon.                                \tag{33.7}
```

The cycle exponent is sharp up to graph-dependent constants: cochain-centre
codes must cover the cycle-norm ball, and arbitrary encoders must distinguish
a strict `2epsilon`-packing.  Thus response entropy is sampled at the
forgetting scale `(1-rho)epsilon`, while persistent dynamic memory costs the
cycle rank independently.

#### Proof

Constants and centred functions are invariant orthogonal channels under
(33.1).  Expanding the path response gives the exact identity

```math
D_p=\left(bar u_(q_t)+\sum_sm_(e_s)\right)1+P_pv_(q_t)
 +\sum_(s=1)^tP_(e_1)...P_(e_(s-1))b_(e_s).          \tag{33.8}
```

Geometric summation proves (33.3).  Delete directed cycles from a path.  Its
remaining simple part has at most `|Q|-1` edges, whereas each deleted cycle
has mean at most `chi_G(m)`; repeating a maximizing cycle proves the reverse
asymptotic bound.  Circulation decomposition proves the stationary-flow
dual, and the usual path-integral argument identifies zero cycle cochains
with gradients on a strongly connected component.

The quotient by gradients has dimension `r_G`.  Cover its radius-`L` ball
by `(1+2L/epsilon)^(r_G)` cycle-norm balls; a bounded linear section and loop
erasure turn each such error into `epsilon t+O_G(epsilon)`.  Recenter any
external centred-dictionary codeword orthogonally.  Terminal error costs
`epsilon` once, while centred reward error `(1-rho)epsilon` sums to at most
`epsilon`.  Repeated exposed cycles give the lower packing law. `square`

Every factor is independently sharp: a two-state eigenmode realizes
`rho^tR`, repeated centred reward realizes `B(1-rho^t)/(1-rho)`, and
singleton fibres on a directed cycle realize `t chi_G(m)`.  The theorem
assumes a fixed common-law contracting carrier; it neither discovers one for
arbitrary hidden dynamics nor counts unbounded endpoint potentials for free.

### Theorem 33.2 (nonlinear secants and the variance toll)

Theorem 33.1 remains valid for switching nonlinear homogeneous maps along
one declared visible path if every realized stochastic secant transports the
same full-support law, contracts its centred `L^2` space by `rho`, and every
same-input approximation defect has a control-visible mean `m_e` and centred
norm at most `B`.  Optimizer switches and ties then create no fourth channel.

Conversely, let a Markov operator preserve `pi` and satisfy
`||P-Pi||_(2->2)<=rho<1`.  If

```math
||f-g||_infinity<=omega,qquad f<=Pf+epsilon1,        \tag{33.9}
```

then

```math
epsilon>=
{(1-rho)^2(\sqrt(Var_pi(g))-omega)_+^2
 \over osc(g)+2omega}.                               \tag{33.10}
```

Thus a response phase with positive stationary variance must survive in the
recovery radius or be paid as fresh transfer toll.

#### Proof

For the nonlinear claim, subtract the two trajectories.  The stochastic
secant propagates the old difference and the same-input defect supplies a
new reward, so (33.8) applies verbatim.

For (33.10), put `h=f-Pf`.  It has mean zero and lies in
`[-osc(f),epsilon]`, whence `||h||_2^2<=osc(f)epsilon`.  With
`v=f-pi f`,

```math
(1-rho)||v||_2<=||(I-P)v||_2=||h||_2.                \tag{33.11}
```

Finally `||v||_2>=sqrt(Var_pi(g))-omega` and
`osc(f)<=osc(g)+2omega`. `square`

For all-finite max-plus maps stochastic secants exist, but a common
transported law and strict centred contraction generally do not.  These are
generator-level falsifiers, not synonyms for the desired quotient.  The
proofs, two independent audits, and 8,099 finite checks are in
[`drafts/contracting_fibre_cocycle_decomposition.md`](drafts/contracting_fibre_cocycle_decomposition.md),
[`drafts/contracting_fibre_cocycle_independent_audit.md`](drafts/contracting_fibre_cocycle_independent_audit.md),
[`drafts/contracting_fibre_cocycle_second_audit.md`](drafts/contracting_fibre_cocycle_second_audit.md),
and
[`experiments/verify_contracting_fibre_cocycle.py`](experiments/verify_contracting_fibre_cocycle.py).

## 34. Dynamic broadcast incidence

### Theorem 34.1 (only recurrent scalar incidence survives at positive rate)

Fix one strongly connected common-law carrier satisfying Theorem 33.1 with
centred contraction `rho<1`.  Let `z in {0,1}^h` index terminal and edge-
reward dictionaries on this carrier.  For a bit flip `z->z+e_i`, write

```math
Delta u_q^i=Delta bar u_q^i1+Delta v_q^i,
\qquad
Delta a_e^i=Delta m_e^i1+Delta b_e^i,                \tag{34.1}
```

and define

```math
U_i=max_q|Delta bar u_q^i|,\quad
R_i=max_q||Delta v_q^i||_2,\quad
B_i=max_e||Delta b_e^i||_2,
```

```math
M_i=max_e|Delta m_e^i|,\qquad
C_i=chi_G(Delta m^i).                                \tag{34.2}
```

If `d_L(z,z+e_i)` is the largest `L^2` response difference over visible
paths of length `L`, then

```math
d_L(z,z+e_i)
<=U_i+rho^LR_i+{1-rho^L\over1-rho}B_i
  +LC_i+(|Q|-1)M_i,                                  \tag{34.3}
```

and, with the usual vacuous convention when there is no recurrent cycle,

```math
limsup_(L->infinity){d_L(z,z+e_i)\over L}=C_i.       \tag{34.4}
```

Now suppose the dictionaries have bounded-atom presentations

```math
u_q^z=\sum_alpha c_alpha(z)phi_(alpha,q),\qquad
a_e^z=\sum_beta d_beta(z)psi_(beta,e),               \tag{34.5}
```

where every atom has maximum component `L^2` norm at most one, coefficient
oscillation `omega`, and hidden dependency set `I`.  Put

```math
J_U=\sum_alpha|I_alpha|omega_alpha,
\qquad J_A=\sum_beta|I_beta|omega_beta.              \tag{34.6}
```

Then

```math
\sum_(i=1)^h d_L(z,z+e_i)
<=2J_U+
\left(1+{1-rho^L\over1-rho}+L+|Q|-1\right)J_A,      \tag{34.7}
```

while the exact rate law sharpens to

```math
\sum_(i=1)^h\limsup_(L->infinity)
 {d_L(z,z+e_i)\over L}
<=J_A.                                               \tag{34.8}
```

Thus `h` neighbouring hidden states each separated by persistent rate
`epsilon` require `h epsilon<=J_A`.  In particular, `E_0` scalar atoms of
oscillation at most `2B` and fan-in at most `t` support at most

```math
h<={2BtE_0\over epsilon}                             \tag{34.9}
```

such rate-visible coordinates.

#### Proof

Apply (33.3) to the difference dictionaries.  Loop erasure bounds the scalar
path sum by `LC_i+(|Q|-1)M_i`, and repeating a maximizing cycle proves the
rate equality.  Orthogonal mean/centred projection and the triangle
inequality charge each terminal atom once to `U_i,R_i` and each edge atom
once to `B_i,M_i`.  Since `C_i<=M_i`, summing bit--atom incidences proves
(34.7)--(34.9). `square`

The three scales are separately sharp by two-state eigenmodes and singleton
cycle fibres.  This theorem does not discover the common-law carrier.  It
does prove a conversion law once one exists: centred broadcast is forgotten,
whereas positive-rate hidden information must be represented in recurrent
cohomology or freshly rebroadcast.  The full proof and audit are in
[`drafts/dynamic_broadcast_incidence_law.md`](drafts/dynamic_broadcast_incidence_law.md)
and
[`drafts/dynamic_broadcast_incidence_law_independent_audit.md`](drafts/dynamic_broadcast_incidence_law_independent_audit.md).

## 35. Extremal witness transversals and orbit-query rates

### Theorem 35.1 (query complexity is reciprocal extremal mass)

Let a finite group `G` act on itself by translation, let `f:G->R`, and put

```math
Q=\max_x|f(x)|,
\qquad
W_\alpha=\{x:|f(x)|\ge\alpha Q\},
\qquad p_\alpha={|W_\alpha|\over|G|}.               \tag{35.1}
```

A coordinate library `X subset G` sees an `alpha`-extremal witness in every
translate `f_s(x)=f(sx)` if and only if it meets every translate of
`W_alpha`.  Its minimum size `L_alpha` satisfies

```math
{1\over p_\alpha}\le L_\alpha
\le\left\lceil{\log|G|+1\over p_\alpha}\right\rceil.
                                                               \tag{35.2}
```

This bound is exponentially sharp.  Let `k=s^2`, partition the Boolean cube
into `s` blocks of size `s`, and put coefficient two within each block and
zero between blocks.  The resulting quadratic `f=H_D` has

```math
||D||_(2->2)<2\sqrt k,
\qquad Q=k^(3/2)-k,                                 \tag{35.3}
```

but, for every fixed `0<alpha<1`,

```math
p_\alpha\le\exp\{-\alpha k/4+O(\sqrt k)\}.         \tag{35.4}
```

Hence every coordinate-pin library for its full switching orbit has size

```math
\exp\{\alpha k/4-O(\sqrt k)\}.                     \tag{35.5}
```

Moreover `D=A-A'` for two exact hollow signings with
`||A||+||A'||=O(sqrt k)`, and switching preserves these bounds.

#### Proof

Meeting all translates is equivalent to the product-set identity
`G=W_alpha X`, which gives the lower bound in (35.2).  Independent uniform
sampling misses a fixed translate with probability at most
`exp(-p_alpha|X|)`; union bound over `G` gives the upper bound.

Writing `M_b` for the block magnetizations gives

```math
H_D(x)=\sum_(b=1)^s(M_b^2-s).                       \tag{35.6}
```

For `lambda=1/4`, Gaussian integration and
`cosh t<=exp(t^2/2)` bound
`E exp(lambda M_b^2/s)` uniformly.  Chernoff over the independent blocks
proves (35.4).  To realize exact flat signings, prescribe opposite signs on
within-block edges and one common random sign completion between blocks;
a standard net bound gives completion norm `O(s)`. `square`

On the Hamming-distance landscape, (35.2) is exactly the sphere-covering
law: an `alpha`-extremal library is, after antipodal complementation, a
binary covering code of radius `(1-alpha)k`.  Thus the same state law applies
outside quadratic forms.

### Theorem 35.2 (product-orbit query rate is the Cramer rate)

Let `G_0` be a finite abelian group of order `q`, let nonconstant
`f:G_0->R`, and define

```math
F_n(x)=\sum_(i=1)^n f(x_i).                         \tag{35.7}
```

For `Ef<a<max f`, let `L_n(a)` be the smallest coordinate library which,
for every translate `s`, contains an `x` with `F_n(s+x)>=an`.  Define

```math
\Lambda_f(\theta)=
 \log\left({1\over q}\sum_(z\in G_0)e^(\theta f(z))\right),
\qquad
I_f(a)=\sup_(\theta\ge0)\{\theta a-\Lambda_f(\theta)\}.
                                                               \tag{35.8}
```

Then

```math
\boxed{\lim_(n->infinity){1\over n}\log L_n(a)=I_f(a).}
                                                               \tag{35.9}
```

For heterogeneous products with asymptotic type proportion `lambda`, the
state composes by

```math
\Lambda_\lambda
=\lambda\Lambda_f+(1-\lambda)\Lambda_g,
\qquad
I_\lambda=\Lambda_\lambda^*.                       \tag{35.10}
```

The corresponding minimum query libraries have exponential rate
`I_lambda(a)`.

#### Proof

Theorem 35.1 sandwiches `L_n(a)` between the reciprocal upper-tail
probability and that quantity times `n log q+1`.  Chernoff gives one
large-deviation inequality.  For the reverse inequality, group equal values
of `f`; a type `nu` has probability
`exp(-nD(nu||mu)+O(log n))`.  Minimization subject to `E_nu f>=a` and
finite-dimensional entropy duality give exactly (35.8).  The polynomial
factor disappears on the logarithmic scale.  The heterogeneous statement
uses the same type proof and additivity of log moments. `square`

The theorem is an exact fusion of contextual response complexity and
classical large deviations, not a claim that one-point energy histograms
survive arbitrary interaction.  Cross-block overlap-dependent couplings are
its sharp falsifier.  Full proofs, computations, and independent audit are
in
[`drafts/extremal_witness_transversals.md`](drafts/extremal_witness_transversals.md),
[`drafts/orbit_query_large_deviations.md`](drafts/orbit_query_large_deviations.md),
and
[`drafts/extremal_witness_ld_independent_audit.md`](drafts/extremal_witness_ld_independent_audit.md).

## 36. Near-minimizer frontier

### Theorem 36.1 (diffuse exact-sign selectors carry a growing physical message)

At PC.3 tensor order `N_j=16^j`, there are a common family of hollow
exact-sign children `(A(u))_(u in U_j)` and a common externally labelled
bank of balanced exact-sign compiler contexts such that

```math
|U_j|\ge e^{j/10000}=N_j^{1/(10000\log16)},                  \tag{36.1}
```

every sparse-flip selector is `O(j^(-1/2))`-diffuse against the entire
active-product library, and, uniformly for `u!=v`, query `u` gives the
oriented all-spins-free parent gap

```math
Q(P_j(v\mid u))-Q(P_j(u\mid u))
\ge(\delta_*+o(1))N_j^{3/2},
\qquad \delta_*>0.012.                                      \tag{36.2}
```

The endpoint code is obtained by balancing the first PC.3 coordinate in
every factor and choosing exponentially many tensor-depth words with
pairwise centred-field covariance at most `j/10`.  Their Boolean selectors
then have pairwise overlap at most

```math
{2\over\pi}\arcsin(2/15)+O(j^{-1/2}).                       \tag{36.3}
```

Sparse flips converge to the corresponding rank-one penalties, so each
query damages its own child but not any cross child by a fixed leading
amount.  A rowwise microcanonical block realizes the query exactly and puts
every free endpoint within `o(N_j^(3/2))` of the balanced scalar query
segment, preserving the gap after all spins optimize.

This is a polynomial **one-hot state packing**, carrying
`log_2|U_j|=Theta(log N_j)` bits.  It is not a `2^{|U_j|}` edit cube and does
not place all query shores in one matrix.  For the native-selector model,

```math
N_j^{1/(10000\log16)}\le |U_j|\le\sqrt {N_j}.                \tag{36.4}
```

The result closes Phase I qualitatively: the physical collision is growing
information, not one isolated bit.  It remains a structured Level-4 result,
not a statement about minimizers.  Proof, scope ceilings, and independent
audit are in
[`drafts/multiselector_sparse_flip_rate.md`](drafts/multiselector_sparse_flip_rate.md)
and
[`drafts/phase1_fractional_balance_independent_audit.md`](drafts/phase1_fractional_balance_independent_audit.md).

### Theorem 36.2 (near-minimality forces a thick fractionally balanced shell)

Let `E=binom(n,2)`, let `a` be a hollow signing, and let

```math
Z_n^+=\{(\sigma x_ix_j)_(i<j):\sigma,x_i\in\{+-1\}\},
\qquad
S_u(a)=\{z:\langle a,z\rangle\ge Q(a)-un^{3/2}\}.            \tag{36.5}
```

There is an absolute `C` such that, if

```math
Q(a)\le M_n+\epsilon n^{3/2},
\quad 0<\kappa<1/4,
\quad
\epsilon+C(\sqrt\kappa n^{-1/4}+n^{-1/2})<\kappa/2,          \tag{36.6}
```

then a probability measure `mu` on `S_(2\kappa)(a)` satisfies

```math
{1\over E}\sum_e
 \left(\mathbb E_\mu[a_ez_e]\right)_+
\le
{\epsilon+C(\sqrt\kappa n^{-1/4}+n^{-1/2})
 \over\kappa(1-1/n)}.                                       \tag{36.7}
```

When `Q(a)>2\kappa n^(3/2)`, the normalized `l_1` norm of the complete
edge-mean vector is at most twice the right side.  Thus fixed small
`epsilon`, followed by `n->infinity` and `kappa=sqrt(epsilon)`, gives shell
width and total edge bias `O(sqrt(epsilon))`.  For exact minimizers,
`kappa=n^(-1/6)` gives both quantities `O(n^(-1/6))`.

Finite minimax identifies the left side of (36.7) as the best common-law
obstruction to a weighted sparse flip.  Independently realize that flip at
rate `kappa/sqrt n`; Bernstein plus a union bound over the `2^n` signed cuts
costs only

```math
O(\sqrt\kappa n^{5/4}+n).                                   \tag{36.8}
```

If (36.7) failed, the realized exact signing would have cap below `M_n`.
This proves the result without knowing a minimizer or its active set.

There is also a cardinality consequence.  If eventually
`M_n>=c_0n^(3/2)` with `c_0>0`, every exact minimizer and every integer
sequence `n<=r=o(n^(3/2))` obey

```math
|S_(2r/n^{3/2})(a)|\ge\exp(c_1r/n)                           \tag{36.9}
```

for `c_1=c_1(c_0)>0`.  Otherwise one random `r`-edge flip would have
positive correlation with every member of the shell, while exact
minimality produces a shell member with nonpositive correlation.  In
particular a shell of vanishing normalized width `2/log n` contains
`exp(Omega(sqrt n/log n))` signed cuts.

Equation (36.7) is the first arbitrary-order Level-5 collective law in this
program.  It is only a first-marginal equilibrium, not a response carrier;
(36.9) is cardinality, not a contextual packing or a lower bound on
generative description length.  The proof and independent audit are in
[`drafts/near_minimizer_fractional_balance.md`](drafts/near_minimizer_fractional_balance.md)
and
[`drafts/phase1_fractional_balance_independent_audit.md`](drafts/phase1_fractional_balance_independent_audit.md).

### Theorem 36.3 (every near-minimizer halo has linear pinned-response entropy)

Let

```math
R_a(g)=\max_x\{H_a(x)+g\mathbin\cdot x\}.
```

There are absolute `c,C,delta>0` such that, whenever
`C/sqrt(n)<=kappa<=1/10`, every exact minimizer `a` has exact-sign
perturbations `(b^u)_(u in U)` with

```math
|U|\ge e^{cn},
\qquad
Q(b^u)\le M_n+2\kappa n^{3/2},                                \tag{36.10}
```

and the common amplitude-`n` query bank satisfies

```math
R_(b^u)(nu)-R_(b^v)(nu)
\ge\delta\kappa n^{3/2}qquad(u\ne v).                       \tag{36.11}
```

Choose exponentially many query spins with
`|u dot v|<=n/2`.  For each `u`, independently flip a
`kappa/sqrt(n)` fraction of the edges on which `a_eu_iu_j=-1`.  The exact
expected oriented gap is

```math
{\kappa\over2\sqrt n}\left(n^2-(u\mathbin\cdot v)^2\right), \tag{36.12}
```

because the base energy cancels between the target and cross channels.
Bernstein concentration closes simultaneously over `e^{2cn}` ordered
pairs, while every perturbation changes at most `kappa n^(3/2)` edges.
Finally the exact pinning identity

```math
R_b(nu)=n^2+H_b(u)                                              \tag{36.13}
```

follows because changing `d` spins can gain at most `2d(n-d)` internally
and loses `2nd` in the field.

Thus any uniform decoder below half the gap in (36.11) needs `Omega(n)`
bits even on a genuine near-minimizer class.  This is an unrestricted-query
negative theorem, not yet a low-cap physical obstruction: realizing full
pinning in an all-spins-free parent can introduce a common
`Theta(n^2)` baseline.  It narrows the remaining incompressibility lemma to
balanced continuations whose whole cap stays `O(n^(3/2))`.  Proof and audit
are in
[`drafts/nearmin_pinned_response_packing.md`](drafts/nearmin_pinned_response_packing.md)
and
[`drafts/nearmin_level5_theorem_independent_audit.md`](drafts/nearmin_level5_theorem_independent_audit.md).

### Theorem 36.4 (Grothendieck peeling envelope and its sharp near-minimizer limitation)

For every hollow symmetric real matrix `A`,

```math
\|A\|_(\infty\to1)\le4Q(A),
\qquad
\boxed{\|A\|_*\le4K_GQ(A)}.                                  \tag{36.14}
```

Indeed the bilinear form polarizes into two quadratic forms on disjoint
principal supports, and symmetric Grothendieck factorization writes
`A=DTD`, with `tr(D^2)=1` and
`||T||<=K_G||A||_(infinity->1)`.  Consequently, for every `L>=1`, fewer
than `n/L` vertices can be removed so that the surviving principal matrix
obeys

```math
\|A[V\setminus S_L]\|_(2\to2)
\le {4K_GLQ(A)\over n}.                                      \tag{36.15}
```

The trace-norm envelope is edit robust: changing `r` unordered signs changes
`||A||_*` by at most `4r`.  It is nevertheless too weak for target-scale
response compression.  Sublinear spectral rank is obtained only above a
diverging multiple of `sqrt n`, whose discarded Boolean error exceeds the
target scale.

Moreover this limitation cannot be repaired by a polynomially sublinear
fixed-root peel.  Around any exact minimizer, overwrite `m` disjoint
`k`-vertex principal blocks by all-positive cliques.  The resulting signing
`B` satisfies

```math
Q(B)\le M_n+mk(k-1),                                           \tag{36.16}
```

while every deletion set producing operator norm at most `R` has size at
least

```math
m(k-1-R)_+.                                                    \tag{36.17}
```

With `k=Lsqrt n`, `m=sqrt n/L^3`, and a subpolynomial
`L->infinity`, the excess is `O(n^(3/2)/L)` but an `O(sqrt n)` core needs
`n^(1-o(1))` deleted vertices.  Thus every claimed
`O(n^(1-delta))` peel is false.  A qualitative `o(n)` peel remains open and
would still lack a composable response theorem.

Proof, primary-source mapping, and independent audit are in
[`drafts/nearmin_spectral_harmonic_report.md`](drafts/nearmin_spectral_harmonic_report.md)
and
[`drafts/nearmin_level5_theorem_independent_audit.md`](drafts/nearmin_level5_theorem_independent_audit.md).

### Theorem 36.5 (deep-hole shell balance forces global edge geometry)

Use the augmented-cut-code dictionary

```math
M_n=N-2\rho(C_n),\qquad N={n\choose2}.
```

If a shell law from Theorem 36.2 has normalized signed edge barycentre
`N^(-1)sum_e|m_e|<=delta`, two independent shell words `Z,Z'` satisfy

```math
0\le \mathbb E R(Z,Z')={1\over N}\sum_em_e^2\le\delta,
\qquad
R(Z,Z')=1-{2d_H(Z,Z')\over N}.                   \tag{36.18}
```

Consequently the shell has Hamming diameter at least
`(1-delta)N/2`.  For an exact minimizer with the explicit
`kappa=n^(-1/6)` law, this is

```math
\left(1-O(n^{-1/6})\right){N\over2}              \tag{36.19}
```

inside excess code radius `n^(4/3)`.

The cardinality theorem also has a metric strengthening.  If its excess
radius `r` obeys `r/(n log n)->infinity`, the same shell contains

```math
\exp(\Omega(r/n))
```

codewords at mutual edge-Hamming distance `Omega(r/log n)`.  This is raw
carrier geometry, not contextual response separation.

The usual uniformly-packed-code shortcut is unavailable.  The augmented
cut code has external distance `N/2-O(n)`, while its covering radius is at
most `N/2-O(n^(3/2))`; hence the equality needed for uniform packing fails
by `Omega(n^(3/2))`.  Deep cosets already have nonconstant leader
multiplicities `3,4,7` at order seven.  Proof, literature mapping, and
independent audit are in
[`drafts/nearmin_deep_hole_shell_report.md`](drafts/nearmin_deep_hole_shell_report.md)
and
[`drafts/nearmin_deep_hole_shell_audit.md`](drafts/nearmin_deep_hole_shell_audit.md).

### Theorem 36.6 (two sharp limits of shell-to-spectral synchronization)

First, an exact-sign spectral-roof selector lemma is too strong to be a
strict convergence reduction.  If exact signings `A,B` obey

```math
d_\square(A,B)\le h,
```

and one Boolean selector of `B` has defect `beta` relative to a roof
`r>=||B||_(2->2)`, then

```math
\boxed{Q(A)\ge {1-\beta\over2}n\sqrt{n-1}-h.}     \tag{36.20}
```

Thus `beta=o(1)` and `h=o(n^(3/2))` for every exact minimizer would already
prove the conjectural limit `1/2`, not merely convergence to an unknown
constant.  The Frobenius identity for exact signings is the obstruction;
weighted cores evade this particular implication but lack an all-order
realization theorem.

Second, optimal first-marginal shell balance does not supply a coherent
spectral side.  At Walsh orders `N=2^(2m)` there is a hollow signing of cap
`N^(3/2)/2` with an exact active augmented-cut law whose normalized edge
bias is the minimum possible `sqrt(N)/(N-1)`, yet a
`Theta(log N)`-port majority selector has spectral defect
`1-O(1/log N)` (and defect exactly one on an infinite subsequence).  A
different orientation-pure sublaw has exact odd-product closure and defect
`O(N^(-1/2))`.  Hence the missing datum is coherent product orientation,
not another first-moment estimate.  This Walsh falsifier is Level 4 because
its cap is not known to be `o(1)`-near-minimal.

The two audited conclusions respectively falsify exact-sign spectral-roof
coercivity as a demonstrably weaker SML and the generic arrow
“FB.1 balance implies selector coercivity.”  Details are in
[`drafts/selector_roof_assumption_gap.md`](drafts/selector_roof_assumption_gap.md),
[`drafts/selector_roof_assumption_gap_audit.md`](drafts/selector_roof_assumption_gap_audit.md),
and
[`drafts/shell_to_selector_coercivity.md`](drafts/shell_to_selector_coercivity.md).

### Theorem 36.7 (universal multiscale orientation-pure affine shells)

Let `A` be any hollow signing of order `n`, let `Q=Q(A)`, and let
`2<=q<=n`.  There are an absolute ground state `x`, an orientation
`rho in {+-1}`, and a set `I` with `|I|>=floor(n/q)` such that

```math
\boxed{
\rho H_A(x^S)\ge\left(1-{8\over q}\right)Q
\quad\hbox{for every }S\subseteq I.}             \tag{36.21}
```

The `2^|I|` spins are projectively distinct and closed under every odd
coordinatewise product.  The proof combines three exact budgets.  After
orienting and switching at `x`, the nonnegative local fields sum to `2Q`;
on every vertex partition, the negative one-sided caps of the principal
blocks sum to at most `Q`; and a subset-edge sum in one block is bounded
below by minus that block's negative cap.  Averaging
`2` times local-field mass plus `4` times negative block cap gives the
constant eight.

Consequently every bounded-cap sequence and every `q_n->infinity` has an
orientation-pure affine cube of dimension `floor(n/q_n)` in normalized
shell width `O(1/q_n)`.  Equivalently, a shell of absolute budget `Delta`
in the stated nontrivial range contains an affine cube of dimension at
least

```math
\left\lfloor {n\Delta\over16Q}\right\rfloor.      \tag{36.22}
```

For `Q=O(n^(3/2))`, this is
`log|S_Delta|=Omega(Delta/sqrt n)`.  A constant-distance mask subcode gives
`exp(Omega(n/q_n))` same-orientation cut words at mutual edge distance
`Theta(n^2/q_n)`.

There is also a growing one-sided response law.  An even subset of `I`
defines the star frame consisting of `x` and its single-coordinate flips.
Every majority endpoint selector lies back in the affine cube, and

```math
0\le Q+m||W\epsilon||_1
-\max_y\{H_(\rho A)(y)+m(W\epsilon)\mathbin\cdot y\}
\le {8Q\over q}                                      \tag{36.23}
```

for every `m>=0` and endpoint `epsilon`.  Conditional on the declared
gauge, this nearly linear language has an `O((n/q)log n)`-bit star
presentation.

This is a new universal Level-5 response benchmark, not near-minimizer
rigidity.  It proves that vanishing-width entropy, odd-product closure, and
common orientation are too generic to identify the minimizers.  Its raw
physical frame has an endpoint of order `n^2/q`; balanced scalar compilation
retains the ordinary `n^(3/2)/sqrt q` random-bridge error and has no
cross-level congruence.  Proof and independent verification are in
[`drafts/multiscale_partition_affine_shell.md`](drafts/multiscale_partition_affine_shell.md),
[`drafts/multiscale_partition_affine_shell_audit.md`](drafts/multiscale_partition_affine_shell_audit.md),
and
[`experiments/verify_multiscale_partition_affine_shell.py`](experiments/verify_multiscale_partition_affine_shell.py).

### Theorem 36.8 (multi-anchor compilation and its projection ceiling)

The favorable coordinate support in Theorem 36.7 can be chosen commonly
for any finite family.  On one equitable `q`-cell partition and for
signings `A^(1),...,A^(t)`, one cell of size at least `floor(n/q)` has,
for every child, one-sided affine defect at most

```math
{8\sum_(u=1)^tQ(A^(u))\over q}.                 \tag{36.24}
```

Thus coordinate support is compatible for `t=o(q)` bounded-cap children;
their gauges and mixed-field responses need not be.

There is also a genuinely multi-anchor exact-sign compiler.  Partition an
`s`-vertex shore into `d` nonempty cells, let `P` average inside each cell,
and prescribe every old-row sum in every cell.  Some exact sign bridge `B`
realizes all those sums and obeys

```math
\max_(\eta\in\{+-1\}^s)||B\eta-BP\eta||_1
\le C\sqrt{n(s-d)(n+s)}.                         \tag{36.25}
```

All block-constant endpoints are therefore preserved jointly rather than
as separately paid scalar channels.

The natural projection architecture nevertheless has a sharp target-scale
ceiling.  For every rank-`d` orthogonal projection `P` and every exact sign
bridge,

```math
\boxed{
\max_\eta||B(I-P)\eta||_1
\ge {ns-||BP||_F^2\over\sqrt{2s}}.}              \tag{36.26}
```

If `P` has a Boolean orthonormal anchor basis whose bridge images each have
`l_1` norm at most `Ln`, then

```math
\max_\eta||B(I-P)\eta||_1
\ge n\sqrt{s/2}\left(1-{dL^2n\over s^2}\right). \tag{36.27}
```

At fixed shore ratio, a sublinear number of balanced linear anchors hence
leaves an `Omega(n^(3/2))` channel whenever that channel is paid separately
by trust-response Lipschitzness.  This is not a lower bound on the actual
parent cap and not an information-theoretic impossibility theorem: joint
child--bridge cancellation before absolute values, or a nonlinear
cross-level congruence, remains outside its scope.

Proof and audit are in
[`drafts/multiscale_partition_composition_audit.md`](drafts/multiscale_partition_composition_audit.md)
and
[`drafts/multiscale_partition_composition_independent_audit.md`](drafts/multiscale_partition_composition_independent_audit.md).

### Theorem 36.9 (every exact minimizer has a balanced affine atlas)

Let `A` be a hollow signing of cap `Q`.  For a signed cut
`z(sigma,x)` of deficit

```math
d=Q-\langle a,z(sigma,x)\rangle,
```

switch and orient at `(sigma,x)`, and let `ell_i` be the resulting local
fields.  Their total negative mass satisfies

```math
\boxed{\sum_i(-\ell_i)_+\le2\sqrt{Qd}.}             \tag{36.28}
```

Consequently, for every `2<=q<=n`, the original atom lies in an affine
chart of dimension `floor(n/q)` on which every atom retains the same
orientation and has deficit at most

```math
d+{8Q+4\sqrt{Qd}\over q}.                           \tag{36.29}
```

This thickening is performed in place; it does not apply a local-ascent map
and hence does not move the centre's edge word.

Suppose a shell law has deficit at most `d_0` and normalized signed edge
bias at most `delta`.  For every `K,q`, there is an existential atlas of at
most `K` such charts whose uniform chart-and-mask law has bias at most

```math
\delta+K^{-1/2}+{2\lfloor n/q\rfloor\over n-1},      \tag{36.30}
```

and whose atoms have deficit at most

```math
D(q,d_0)=d_0+{8Q+4\sqrt{Qd_0}\over q}.              \tag{36.31}
```

Each chart also gives a declared one-block response language.  If `W_r`
is the star frame formed from its centre and an even set of singleton
flips, then for every endpoint `epsilon` and `m>=0`,

```math
0\le Q+m\|W_r\epsilon\|_1
-\max_y\{\sigma_rH_A(y)+m(W_r\epsilon)\mathbin\cdot y\}
\le D(q,d_0).                                         \tag{36.32}
```

The same estimate holds for the common absolute trust response after
maximizing over the quadratic sign.

Applying Theorem 36.2 to an exact minimizer with

```math
\kappa=n^{-1/6},\qquad K=\lceil n^{1/3}\rceil,
\qquad q=\lceil n^{1/6}\rceil
```

gives, at every sufficiently large order, at most `O(n^(1/3))` charts of
dimension `Theta(n^(5/6))`, normalized atlas bias `O(n^(-1/6))`, and
common-oriented shell and declared-response error `O(n^(4/3))`.  A direct
labelled presentation costs `O(n^(4/3))` bits while generating
`exp(Theta(n^(5/6)))` endpoints per chart.

This is a strict, Level-5, signing-dependent **designed-interface
certificate**.  It is not an efficient encoder, a common exogenous query
quotient, a mixed-chart response theorem, or a reusable cross-order
congruence.  Walsh mixed orientations falsify the generic cross-chart
implication, while Theorem 36.8 obstructs separately paying the omitted
linear channel at fixed shore ratio.  Proof, finite verification, and
independent audit are in
[`drafts/nearmin_balanced_affine_atlas.md`](drafts/nearmin_balanced_affine_atlas.md),
[`experiments/verify_nearmin_balanced_affine_atlas.py`](experiments/verify_nearmin_balanced_affine_atlas.py),
and
[`drafts/nearmin_balanced_affine_atlas_audit.md`](drafts/nearmin_balanced_affine_atlas_audit.md).

### Theorem 36.10 (the universal thin-shell entropy floor)

Let `A` be a hollow signing with `Q(A)<=Cn^(3/2)`, orient a ground state,
and fix `c>0`.  At least

```math
{1\over2}{\lfloor n/2\rfloor\choose\lfloor c\sqrt n\rfloor}
=2^{(c/2+o(1))\sqrt n\log_2n}                     \tag{36.33}
```

projectively distinct spins have one common energy orientation and deficit
at most

```math
(8cC+2c^2+o(1))n.                                 \tag{36.34}
```

Indeed at least `n/2` ground-state local fields are at most `4Q/n`; flip
every `floor(c sqrt n)`-subset of that pool and use the exact flip identity.
Unlike Theorem 36.7, this family is not asserted to be one affine algebra.

The scale is attained already at exact cap by a structured bounded-cap
family.  If `N=2^d` with even `d`, the hollow Sylvester Walsh signing has
cap `N^(3/2)/2` and

```math
\log_2|\mathcal G_d^+|
\ge(1/2+o(1))\sqrt N\log_2N.                    \tag{36.35}
```

projective positive maximizers.  This follows by injecting every bent
function in `d-2` variables into a self-dual bent function in `d` variables
and applying Haugland's bent-function lower bound.  It applies to the
unflipped hollow PC.3 child tower under diagonal conjugacy, not to its
sparse-flipped or completed descendants.  No matching Walsh upper count is
known.

Two stability facts delimit the conclusion.  If `B` is obtained from an
exact minimizer `A` by `r` edge flips, then

```math
\mathcal S_B(\Delta)\subseteq\mathcal S_A(\Delta+2r),
\qquad
\mathcal S_A(\Delta)\subseteq\mathcal S_B(\Delta+4r).       \tag{36.36}
```

Thus black-box `o(n^(3/2))` halo proximity does not preserve an `O(n)`
shell.  Conversely, a ground-state vertex set `J` with total local-field
mass `L(J)` and induced subset discrepancy

```math
P(J)=\max_{S\subseteq J}
 \left|\sum_{\{i,j\}\subseteq S}a_{ij}x_ix_j\right|
```

generates at least `2^(|J|-1)` common-oriented shell states of deficit at
most `2L(J)+4P(J)`.  A mesoscopic `J` with `L(J)+P(J)=O(n)` is therefore a
precise falsifier for any proposed smaller thin-shell upper bound.

This theorem concerns witness cardinality, not contextual response entropy
or generative description length.  Proof, primary-source mapping, and
independent audit are in
[`drafts/nearmin_thin_shell_entropy_falsifier.md`](drafts/nearmin_thin_shell_entropy_falsifier.md)
and
[`drafts/nearmin_thin_shell_entropy_falsifier_audit.md`](drafts/nearmin_thin_shell_entropy_falsifier_audit.md).

### Theorem 36.11 (projective shell packing has a low-cap physical compiler)

Let `A` be a bounded-cap signing and suppose its positive deficit-`d` shell
contains signed cuts `z^u=sigma_uc(u)`, `u in U`, with

```math
{ |\langle z^u,z^v\rangle|\over {n\choose2}}
\le1-\gamma\qquad(u\ne v).                         \tag{36.37}
```

For fixed `alpha,lambda,gamma>0`, independently flip an
`alpha/sqrt n` fraction of the edges on which `A` disagrees with each
`z^u`, producing exact-sign children `b^u`.  Query child `v` with the
all-spins-free rank-one shore

```math
B^u=u\mathbf1_h^T,qquad h=\lfloor\lambda\sqrt n\rfloor.
```

One simultaneous realization satisfies, for every `u ne v`,

```math
Q(P^{u|u})-Q(P^{v|u})
\ge\left({\gamma\over4}\min\{\alpha,\lambda\}-o(1)\right)n^{3/2}
 -(1-\alpha/\sqrt n)d.                              \tag{36.38}
```

The matrices are hollow exact signings of order `n+O(sqrt n)` and their
whole caps are `O(n^(3/2))`.  The proof keeps the sparse-flip and field
channels joint until after maximization.  Its key spherical inequality is

```math
{\alpha\over2}(V\mathbin\cdot Y)^2
+\lambda|U\mathbin\cdot Y|
\le {\alpha\over2}+\lambda
-{1-|U\mathbin\cdot V|^2\over4}\min\{\alpha,\lambda\}.     \tag{36.39}
```

Thus a `K`-point vanishing-width projective shell packing yields a common
low-cap physical query bank carrying `log_2K` bits.  This conclusion is
conditional on the projective packing, not on signed first marginals.

There is also an unconditional ordered consequence.  Take one ground state
`u` and any Boolean `v` with `|u dot v|<=1`.  Only the target `u` must be
near-top, so every bounded-cap signing has two halo children and one common
free-shore context with

```math
Q(P^{u|u})-Q(P^{v|u})
\ge(\alpha/4-o(1))n^{3/2}.                          \tag{36.40}
```

For an exact-minimizer base this is one genuine Level-5 physical hidden bit.
It persists in a vanishing halo: if `alpha_n->0` and
`alpha_n sqrt n->infinity`, both children lie in
`\mathcal N_n(2alpha_n)` and the gap is

```math
(\alpha_n/4-o(\alpha_n))n^{3/2}.                   \tag{36.41}
```

This does not amplify to a growing state family without additional near-top
directions.

Finally, fractional edge balance has an exact alternative.  If a positive
`d`-shell law has normalized edge bias at most `delta` but no pair satisfying
(36.37), then relative to any support atom it lies in two projective caps
whose orientation bit has bias at most `delta+gamma`.  Positivity forces

```math
\gamma>{Q(A)-d\over {n\choose2}}.                  \tag{36.42}
```

Hence Theorem 36.2 alone reaches only the `Theta(n^(-1/2))` separation
scale: either it gives such a pair or it permits two nearly antipodal caps.
At that scale (36.38) has only `O(n)` gap.  Ruling out this two-cap geometry
at a fixed scale, or exploiting it with a different compiler, is the precise
remaining structural obligation.

Proof, independent audit, and the overflow-safe finite verifier are in
[`drafts/nearmin_absolute_overlap_physical_compiler.md`](drafts/nearmin_absolute_overlap_physical_compiler.md),
[`drafts/nearmin_absolute_overlap_physical_compiler_independent_audit.md`](drafts/nearmin_absolute_overlap_physical_compiler_independent_audit.md),
and
[`experiments/verify_nearmin_absolute_overlap_compiler.py`](experiments/verify_nearmin_absolute_overlap_compiler.py).

### Theorem 36.12 (vanishing near-minimizers can have one-cap shells)

Let `A` be a bounded-cap signing, let `z_0` be an oriented ground cut, and
put

```math
D=\{e:a_e(z_0)_e=-1\},\qquad |D|={{n\choose2}-Q(A)\over2}.
```

For every sequence `n=o(r)` and `r=o(n^2)`, some `r`-set `F subset D`
has the following property.  Flip precisely `F`, obtaining `B`.  Then

```math
\boxed{Q(B)=Q(A)+2r}                                      \tag{36.43}
```

and every positive deficit-`Delta` shell atom obeys

```math
{d_E(z,z_0)\over {n\choose2}}
\le {32n\over r}+{\Delta\over r}
       +{\Delta\over2{n\choose2}}.                         \tag{36.44}
```

The cap identity is exact:

```math
Q(A)+2r-\langle B,z\rangle
=Q(A)-\langle A,z\rangle+4|F\cap\{e:z_e\ne(z_0)_e\}|.      \tag{36.45}
```

The new ingredient beyond the archived planted-face identity is a uniform
hypergeometric hitting set over all `2^n` augmented cuts.  It converts the
frozen planted edges into the global radius bound (36.44).

Consequently, for every prescribed `Delta_n=o(n^(3/2))`, there is a
certified sequence

```math
Q(B_n)=M_n+o(n^{3/2})                                      \tag{36.46}
```

whose entire positive `Delta_n` shell lies in one vanishing projective cap
and has only `exp(o(n))` projective support.  One may take

```math
r_n=\left\lceil\sqrt{n^{3/2}\max\{n,\Delta_n\}}\right\rceil.
```

This is a scalable Level-5 falsifier for structural claims quantified over
**all** vanishing near-minimizers.  It does not falsify Theorem 36.11's
selected packing lemma when that lemma is restricted to exact minimizers,
nor a shell whose width is comparable to or larger than the planted excess
`2r`.  In particular it does not contradict Theorem 36.2, which deliberately
chooses a shell wider than the near-minimality slack.

Proof and independent audit are in
[`drafts/nearmin_planted_ground_cap_falsifier.md`](drafts/nearmin_planted_ground_cap_falsifier.md)
and
[`drafts/nearmin_planted_ground_cap_falsifier_independent_audit.md`](drafts/nearmin_planted_ground_cap_falsifier_independent_audit.md).

### Theorem 36.13 (deep-hole projective floor and the cut-specific two-cap frontier)

Let `C <= F_2^N` be a binary linear code and let `y` be a deep hole of
covering radius `rho`.  Its radius-`rho+1` error supports cover all `N`
coordinates, and hence

```math
\left|\{c\in C:d(y,c)\le\rho+1\}\right|
\ge\left\lceil {N\over\rho+1}\right\rceil.          \tag{36.47}
```

For the augmented cut code, `rho=(N-M_n)/2`.  Once `M_n>2`, (36.47)
gives at least three distinct projective words in every exact minimizer's
positive deficit-two shell.  Their unconditional separation is only the
projective minimum distance `n-1`, or AO overlap parameter `4/n`.  If three
of their error supports already cover all coordinates, each pair is
separated projectively by at least `M_n-2`, still only normalized scale
`Theta(n^(-1/2))`.

Generic deepest-coset theory cannot improve this to fixed scale.  The exact
antipodal codes

```math
C_m=\operatorname{Rep}_{m^4-m^3}\oplus\mathbb F_2^{m^3}
```

have length `N_m=m^4`, dimension `Theta(N_m^(3/4))`, covering radius
`N_m/2-Theta(N_m^(3/4))`, and a deep hole whose complete
`o(N_m^(3/4))` shell has projective diameter `o(N_m)`.  This shell also
satisfies every simultaneous finite-flip majority-cover certificate.  The
sphere-covering inequality

```math
D\ge {2t^2\over k\log2}                              \tag{36.48}
```

for a `[D,k]` factor of radius `D/2-t` shows why this particular separable
collapse cannot simply be reduced to the cut code's `Theta(sqrt N)`
dimension.

There is nevertheless one exact cut-specific conclusion.  If two positive
shell words of an exact minimizer have deficits at most `2s`, actual mutual
distance `N-D`, and

```math
D<N-\lfloor n^2/4\rfloor,
```

then their agreement coordinates form one complete bipartite cut
`delta(S)` and

```math
D=|S|(n-|S|)\ge M_n-2s,
\qquad
\sum_{e\in\delta(S)}a_e\ge M_n-2s.                 \tag{36.49}
```

Thus a genuinely collapsed `D=o(N)`, `s=o(n^(3/2))` two-cap obstruction
must concentrate a leading bias on a shore with
`Omega(sqrt n)<=|S|=o(n)`.  This is anatomy, not a contradiction.  The
remaining fixed-scale three-point statement is genuinely cut-specific and
open; if proved, Theorem 36.11 compiles it into a three-state low-cap
physical response packing, not convergence.

Proof, finite checks, and the repaired independent audit are in
[`drafts/deep_hole_projective_packing_frontier.md`](drafts/deep_hole_projective_packing_frontier.md)
and
[`drafts/deep_hole_projective_packing_independent_audit.md`](drafts/deep_hole_projective_packing_independent_audit.md).

### Theorem 36.14 (projective shell covers give local-field response roofs)

For a signing `A`, let `S_G^+(A)` be its positive augmented-cut shell of
deficit `G`.  Suppose `L` shell centres projectively cover this shell at
edge radius `R<E-floor(n^2/4)`, and put

```math
k_R=\max\{0\le d\le\lfloor n/2\rfloor:d(n-d)\le R\}.
```

At each centre `z^r=sigma_r c(u^r)`, retain only its baseline and oriented
local fields

```math
h_r=\sigma_rH_A(u^r),
\qquad
\ell_{r,i}=\sigma_ru_i^r(Au^r)_i.                 \tag{36.50}
```

These `O(Ln log n)` bits define an explicitly sortable atlas roof for the
absolute trust response

```math
\mathcal B_A(g)=\max_{\sigma,x}\{\sigma H_A(x)+g\mathbin\cdot x\}.
```

Uniformly for every `||g||_1<=G`, its error is

```math
\boxed{|\widehat{\mathcal B}(g)-\mathcal B_A(g)|
       \le2k_R(k_R-1).}                            \tag{36.51}
```

The proof uses the exact identity

```math
\sigma_rH_A((u^r)^S)
=h_r-2\sum_{i\in S}\ell_{r,i}
 +4\sum_{\{i,j\}\subseteq S}\sigma_ra_{ij}u_i^ru_j^r, \tag{36.52}
```

and drops only the final internal-edge term.  The shell width is exactly
`G`: a response optimizer has oriented quadratic value at least `Q(A)-G`.

If `R=floor(gamma E)`, then `k_R<=gamma(n-1)`.  Thus the atlas error is
`o(n^(3/2))` for `gamma=o(n^(-1/4))`.  An inclusion-maximal
`R`-separated shell net gives at the same time the opposite packing branch,
but Theorem 36.11's physical compiler resolves that branch only for
`gamma>>n^(-1/4)`.  The critical scale remains open.  Same-radius covering
number and maximal-net size are not interchangeable; an `R/2` cover does
bound the latter.

Two qualifications are essential.  First, when `G=o(n^(3/2))`, the scalar
bound `Q<=B_A(g)<=Q+G` already gives target-scale compression, so the atlas
is substantively finer only when `k_R^2=o(G)` or `G` is macroscopic.
Second, this is a static one-block representation: it provides neither an
efficient encoder nor a cross-order update law.  It sharpens the multiscale
frontier but proves no near-minimizer arrow by itself.

Proof, independent audit, and the exhaustive finite verifier are in
[`drafts/nearmin_projective_shell_roof.md`](drafts/nearmin_projective_shell_roof.md),
[`drafts/nearmin_projective_shell_roof_independent_audit.md`](drafts/nearmin_projective_shell_roof_independent_audit.md),
and
[`experiments/verify_nearmin_projective_shell_roof.py`](experiments/verify_nearmin_projective_shell_roof.py).

### Theorem 36.15 (all sublinear principal restrictions remain near-minimal)

Let `n=m+k`.  A random exact-sign rectangular bridge between minimizers of
orders `m` and `k` gives

```math
M_n\le M_m+M_k+\sqrt{2(\log2)mkn},                 \tag{36.53}
```

and a random hollow signing gives

```math
M_k\le\sqrt{(\log2)(k^3-k)}                       \tag{36.54}
```

for `k>=2` (with `M_1=0`).  If an order-`n` signing `A` satisfies
`Q(A)<=M_n+eta`, then **every** principal set `U` of size `m` satisfies

```math
\boxed{
Q(A[U])-M_m
\le eta+M_k+\sqrt{2(\log2)mkn}.}                  \tag{36.55}
```

Indeed, averaging the omitted spins gives `Q(A[U])<=Q(A)`, and (36.53)
then supplies (36.55).  Consequently, if `A_n` is a vanishing
near-minimizer, `k_n=o(n)`, and `U_n` is any sequence of `n-k_n` principal
sets, then

```math
Q(A_n[U_n])-M_{n-k_n}=o(n^{3/2}).                 \tag{36.56}
```

This extends the direct edge-count deletion window from `o(sqrt n)` to all
`o(n)`, but only for the scalar cap.  It transports neither ground states
nor shell/response geometry.  In particular, deleting the exceptional shore
in Theorem 36.13 preserves near-minimality of the large core while erasing
the leading bipartite response channel.

The random-bridge inequality and the consequence that relatively dense
liminf-realizing orders would imply convergence were already in the archive.
The newly isolated content is the elementary uniform principal-heredity
synthesis, not a new composition method and not a fixed-ratio frontier
advance.  Proof and independent audit are in
[`drafts/nearmin_macroscopic_principal_heredity.md`](drafts/nearmin_macroscopic_principal_heredity.md)
and
[`drafts/nearmin_macroscopic_principal_heredity_independent_audit.md`](drafts/nearmin_macroscopic_principal_heredity_independent_audit.md).

### Theorem 36.16 (localized exact flips force a third energy-scale direction)

Let `A` be an exact minimizer, `E=binom(n,2)`, and `M=M_n`.  Suppose
positive shell words `z_1,...,z_k`, each of deficit at most `2s`, have a
common-correct reservoir

```math
R_k=\{e:a_e(z_i)_e=+1\text{ for every }i\},
\qquad p_k=|R_k|.                                  \tag{36.57}
```

Fix `0<theta<1/2`.  If `1<=r<=p_k`, `2r<M`, `s+r<M`, and

```math
2^n\exp\{-2(1/2-\theta)^2r\}<1,                  \tag{36.58}
```

then the same deficit-`2max{s,r}` shell contains a further positive word
`z_(k+1)` satisfying, simultaneously for all old anchors,

```math
\boxed{
d_{\rm P}(z_{k+1},z_i)
\ge\min\{\theta p_k,M-s-r\}.}                    \tag{36.59}
```

Choose an `r`-set uniformly inside `R_k`.  A hypergeometric union bound over
all `2^n` augmented cuts makes every response which is negative on half the
sample negative on more than `theta p_k` reservoir edges.  Flip the sample.
Exact minimality supplies such a response with deficit at most `2r`.  Its
actual distance from every anchor is therefore large, while positivity of
the two words gives the complementary bound `M-s-r`.  Every new witness is
obtained from a fresh perturbation of the same `A`, so iterating (36.59)
does not accumulate shell width.

In the collapsed branch of Theorem 36.13, the first pair has
`p_2>=M-2s`.  Taking

```math
\theta={1\over2}-{1\over\log n},
\qquad r=O(n\log^2n)=o(M),                         \tag{36.60}
```

therefore forces a third word of `o(n^(3/2))` deficit with

```math
d_{\rm P}(z_3,z_i)\ge(1/2-o(1))M_n
\qquad(i=1,2).                                    \tag{36.61}
```

More generally, for fixed `beta>0` the same argument either builds any
prescribed number of mutually `(beta/4-o(1))M`-separated words or reaches an
explicit higher-order failure `|R_k|<beta M`.  Pairwise reservoir mass does
not prevent this failure, even for actual cut words in a bounded-cap exact
signing (Example 159).

This is a theorem at the energy scale, not fixed projective scale:
`M/E=Theta(n^(-1/2))`.  The complementary term in (36.59) is exactly the
archived two-positive-word inequality AO.20.  Thus the displayed
one-reservoir entropy proof reaches its certifiable ceiling; a fixed-scale
packing needs a cut-specific higher-order non-recycling fact or a new way to
control the opposite signed lift.  Proofs and independent audits are in
[`drafts/mesoscopic_bipartite_core_third_witness.md`](drafts/mesoscopic_bipartite_core_third_witness.md),
[`drafts/mesoscopic_bipartite_core_third_witness_independent_audit.md`](drafts/mesoscopic_bipartite_core_third_witness_independent_audit.md),
[`drafts/mesoscopic_core_iteration_frontier.md`](drafts/mesoscopic_core_iteration_frontier.md),
and
[`drafts/mesoscopic_core_iteration_frontier_independent_audit.md`](drafts/mesoscopic_core_iteration_frontier_independent_audit.md).

### Theorem 36.17 (low dimension forces deficit-scale shell diffusion)

Let `C<=F_2^N` be an antipodal `[N,k]` code, let `y` be a deep hole, and
write

```math
\rho(C)=N/2-t,
\qquad
\mathcal L_r(y)=\{c\in C:d(y,c)\le\rho+r\}.
```

If `1<=r<t`, `h=ceil(r/2)`, `m=N/2+t`, and

```math
D=\operatorname{diam}_{\rm P}\mathcal L_r(y)<2(t-r),
```

then the exact cap-cover entropy inequality

```math
\boxed{
|\mathcal L_r(y)|\binom rh(D/m)^h\ge1}            \tag{36.62}
```

holds.  Consequently, whenever

```math
k=o(t),\qquad t=o(N),
```

every integer sequence `k=o(r)=o(t)` obeys

```math
\operatorname{diam}_{\rm P}\mathcal L_r(y)
\ge(2-o(1))t.                                     \tag{36.63}
```

For the augmented cut-code scaling `N=Theta(n^2)`, `k=Theta(n)`, and
`t=Theta(n^(3/2))`, one may take `r=Theta(n^(5/4))`.  Thus low dimension
alone forces thin-shell spread of order `M_n`, but not a fixed fraction of
`N`.

The proof flips every `r`-subset of the complement of one nearest error
support.  Covering-radius optimality forces a shell leader which agrees with
at least half of the flipped coordinates.  Below the threshold `2(t-r)`, an
oppositely oriented Hamming lift cannot enter, so one family of at most
`D`-sets must half-cover every `r`-set.  Counting those incidences gives
(36.62).  The proof stops sharply at `2(t-r)`: beyond that point an opposite
lift can answer many local queries at once.

Two further audited barriers exclude common collapsed constructions.  A
genuine direct product whose radius-one alternatives have aggregate
activation mass at least `N/3` already has a fixed-scale pair in an
`o(t)` shell.  If a projectively collapsed shell contains an affine
subspace of dimension `k-o(k)`, its direction subcode is supported on
`o(N)` coordinates; a split-discrepancy bound then contradicts
`t=Theta(sqrt(Nk))`.  These do not rule out a nonlinear,
affine-subspace-evasive two-cap shell.

This is a generic mesoscopic theorem and a construction barrier, not the
selected fixed-scale `L_projective`.  Its complete proof and independent
finite/proof audit are in
[`drafts/generic_low_dimension_localization_barriers.md`](drafts/generic_low_dimension_localization_barriers.md)
and
[`drafts/generic_low_dimension_localization_barriers_independent_audit.md`](drafts/generic_low_dimension_localization_barriers_independent_audit.md).

### Theorem 36.18 (an exact ground links three shell directions, with a sharp four-anchor ceiling)

Let `A` be an order-`n` signing, let `M=Q(A)`, and let
`z_0,z_1,z_2` be positively oriented augmented-cut words with

```math
\langle a,z_i\rangle=M-d_i.
```

Their common-correct reservoir satisfies the exact cap consequence

```math
\boxed{
|R(z_0,z_1,z_2)|
\ge {M\over2}-{d_0+d_1+d_2\over4}.}              \tag{36.64}
```

Indeed, after gauging by `z_0`, the signed mass of the cell on which the
three words agree is

```math
{1\over4}\big(
 \langle a,z_0\rangle+
 \langle a,z_1\rangle+
 \langle a,z_2\rangle+
 \langle a,z_0z_1z_2\rangle\big),               \tag{36.65}
```

and the odd product is another augmented cut, so its response is at least
`-M`.  The analogous pair bound is

```math
|R(z_0,z_1)|\ge M-{d_0+d_1\over2}.               \tag{36.66}
```

Combining (36.64)--(36.66), Theorems 36.16--36.17, and two fresh localized
flip arguments proves a uniform consequence for exact minimizers.  There is
`D_n=o(M_n)` such that the positive `D_n`-deficit shell of **every** exact
order-`n` minimizer contains four words `g,u,v,w` with

```math
\boxed{
d_{\rm P}(z,z')\ge(1/4-o(1))M_n
\quad(z\ne z').}                                  \tag{36.67}
```

One may use code-shell excess `u_n=Theta(n^(5/4))`, localized-flip sample
size `r_n=O(n log^2 n)`, and `D_n=2max{u_n,r_n}`.  Deficits do not
accumulate: each new word is obtained from a fresh perturbation of the same
exact minimizer.

The corresponding four-anchor identity has a genuine orientation ceiling.
For positive words `z_0,...,z_3`, put `h=z_0z_1z_2z_3`.  The signed mass of
their common cell is

```math
{1\over8}\sum_{i=0}^3
 \big(\langle a,z_i\rangle+\langle a,hz_i\rangle\big). \tag{36.68}
```

Thus a joint bound

```math
\sum_i\langle a,hz_i\rangle\ge-(4-eta)M
```

gives an `eta M/8-o(M)` common-correct reservoir when the deficits are
`o(M)`.  Cap boundedness alone gives only zero at exact ground energy, and
this is sharp: an explicit order-five exact minimizer has four positive
ground words with `h=-mathbf 1`, all pairwise projective distances equal to
`M_5=4`, and empty fourwise reservoir.  Negative holonomy therefore kills
the literal common intersection pointwise, although it does not classify
all possible four-anchor collapse.

This theorem removes the former generic Example 159 as a depth-three
obstruction for genuine thin-shell words.  It remains an energy-scale
result: `M_n/binom(n,2)=Theta(n^(-1/2))`, so it does not prove the selected
fixed-scale `L_projective`.  The proof, exact order-five witness, and
independent audit are in
[`drafts/exact_ground_endpoint_linkage.md`](drafts/exact_ground_endpoint_linkage.md)
and
[`drafts/exact_ground_endpoint_linkage_independent_audit.md`](drafts/exact_ground_endpoint_linkage_independent_audit.md).

### Theorem 36.19 (fractional reservoirs force a growing energy-scale packing)

For every fixed `K` there is a finite constant `C_K` with the following
finite-phase trimming property.  If a multiset of sign patterns
`r_e in {+-1}^K` satisfies

```math
\sum_e r_{i,e}\ge m>0\qquad(1\le i\le K),
```

then there are weights `0<=w_e<=1` such that

```math
m\le W:=\sum_e w_e\le C_Km,
\qquad
\sum_e w_er_{i,e}\ge m\quad(1\le i\le K).       \tag{36.69}
```

One explicit bound is

```math
C_K\le\max\left\{1,
 \max_{2\le p\le K}{p^{1+p/2}\over2^{p-1}}
 \right\}.                                      \tag{36.70}
```

This follows by applying Minkowski--Weyl to the box-capacity vector and
trimming it to a convex combination of vertices.  A vertex uses at most
`K` sign types; Cramer's rule, sign-determinant divisibility, and Hadamard's
inequality give (36.70).  The dependence on `K` cannot be removed from this
abstract lemma: the sharp universal constant is at least `K` for odd `K`
and at least `K-1` for even `K`.

The trimming lemma has a stronger exact-minimizer consequence.  Let
`z_1,...,z_K` be positive responses of one exact minimizer, each of deficit
at most `2s`, and put `m=M_n-2s`.  Independently include edge `e` with
probability `rw_e/W`, where `r=Theta(C_K^2n)`.  A uniform Bernstein bound
over all augmented cuts produces one actual edge set whose sampled response
approximates every weighted response.  After flipping that set, exact
minimality selects a new positive response `z_(K+1)` of deficit at most
`4r` and

```math
\boxed{
d_{\rm P}(z_{K+1},z_i)
\ge\min\left\{{m\over4},M_n-s-2r\right\}
\quad(1\le i\le K).}                              \tag{36.71}
```

The separation constant does **not** deteriorate with `C_K`: the latter
only prices the number of physical edge flips needed to approximate all
responses before the adaptive maximizer is chosen.  Independent Bernoulli
inclusion is a no-replacement construction, so no collision error is
hidden.

Consequently, for every fixed `0<alpha<1/2`, every exact order-`n`
minimizer has a positive shell of deficit

```math
D_n=n^{1+\alpha+o(1)}=o(M_n)                     \tag{36.72}
```

containing

```math
L_n=\left\lfloor{\alpha\log n\over\log\log n}\right\rfloor
\longrightarrow\infty                            \tag{36.73}
```

words with pairwise projective distance

```math
\boxed{(1/4-o(1))M_n.}                            \tag{36.74}
```

Each extension uses a fresh perturbation of the original exact minimizer,
so the shell deficit does not accumulate.  This crosses the literal
four-anchor holonomy ceiling and completes the previously missing growing
packing step at the **energy scale**.  It still does not prove fixed-edge-
scale `L_projective`, because `M_n/binom(n,2)->0`, nor does it by itself
activate the existing low-cap compiler at leading `n^(3/2)` scale.

Proof and independent audit are in
[`drafts/fractional_reservoir_localized_flip.md`](drafts/fractional_reservoir_localized_flip.md)
and
[`drafts/fractional_reservoir_localized_flip_independent_audit.md`](drafts/fractional_reservoir_localized_flip_independent_audit.md).
The exact-arithmetic basic-support verifier and its output are
[`experiments/finite_phase_fractional_reservoir_constants.py`](experiments/finite_phase_fractional_reservoir_constants.py)
and
[`experiments/finite_phase_fractional_reservoir_constants.json`](experiments/finite_phase_fractional_reservoir_constants.json).

## 37. Presented rare-event states

### Theorem 37.1 (generic finite-rank spikes survive bulk normalization)

Let `B_N` be deterministic real symmetric matrices with uniformly bounded
operator norm, empirical spectral laws converging to a compactly supported
probability measure `mu`, and

```math
\lambda_{\max}(B_N)\longrightarrow
b:=\max\operatorname{supp}\mu.
```

For a fixed `theta>0`, let `u_N` be Haar-uniform and independent of `B_N`.
After freezing an almost-surely good sequence of directions,

```math
\lambda_{\max}(B_N+\theta u_Nu_N^{\mathsf T})
\longrightarrow \mathcal R(\mu,\theta),           \tag{37.1}
```

where, with `G_mu(z)=int (z-t)^(-1)dmu(t)` for `z>b`,

```math
\mathcal R(\mu,\theta)=
\begin{cases}
\rho,&\theta G_\mu(b+)>1,
       \quad G_\mu(\rho)=\theta^{-1},\ \rho>b,\\
b,&\theta G_\mu(b+)\le1.
\end{cases}                                      \tag{37.2}
```

At finite `N` the exact carrier is the rooted spectral measure

```math
\nu_{N,u}=\sum_i|\langle u,v_i\rangle|^2
                  \delta_{\lambda_i(B_N)}.        \tag{37.3}
```

Haar isotropy synchronizes this rooted carrier to the unrooted bulk law.
For a fixed jointly Haar-generic frame and fixed positive spike multiset
`Theta`, the finite-rank secular equation therefore yields an asymptotically
exact **presented** multiset-union update for all upper outliers.  The update
is not uniform over directions selected adversarially after the presentation
is frozen.

The spike is invisible to every fixed collection of normalized trace
moments and to the limiting weak empirical law.  In the supercritical case,
nevertheless,

```math
\lambda_{\max}(B_N+\theta u_Nu_N^{\mathsf T})
-\lambda_{\max}(B_N)\longrightarrow\rho-b>0.      \tag{37.4}
```

Thus empirical mass `1/N` can carry order-one extremal response and must be
retained at a separate normalization.  This is a strict rare-event
composition theorem for a spherical quadratic benchmark, but only at Level
3 here: `mu` may be infinite-dimensional, the state is presentation-
dependent, and correlated spikes require relative Gram geometry.  It has no
current Boolean near-minimizer implication.

The proof, primary citations, scope qualifications, and independent audit
are in
[`drafts/rare_event_spectral_spike_state.md`](drafts/rare_event_spectral_spike_state.md)
and
[`drafts/rare_event_spectral_spike_state_independent_audit.md`](drafts/rare_event_spectral_spike_state_independent_audit.md).

### Theorem 37.2 (rare-phase activation has a sharp adaptive replica exponent)

Let a landscape on `Omega_n` be zero off a marked set `P_n` and have reward
`E_n` on `P_n`, where

```math
q_n={|P_n|\over|\Omega_n|}=\exp\{-nI+o(n)\},
\qquad E_n=n\delta+o(n),                            \tag{37.5}
```

with `I,delta>0`.  At inverse temperatures `0<=beta<=B<I/delta`, set

```math
a_B=I-B\delta>0.
```

Then the complete normalized pressure curve differs from the zero
landscape by at most

```math
{1\over n}\exp\{-na_B+o(n)\},                     \tag{37.6}
```

and this exponential rate is attained at `beta=B`.  More strongly, an
adaptive observer choosing each temperature from the preceding Gibbs
samples has transcript total variation at most

```math
K_n\exp\{-na_B+o(n)\}                              \tag{37.7}
```

after `K_n` samples.  Every strict smaller sample exponent is therefore
blind, including all sampled multi-overlap statistics, while every strict
larger exponent is sufficient: sampling only at `B` and counting visits to
`P_n` distinguishes the landscapes consistently.  The exact critical
conditions are `K_nd_n(B)->0` for blindness and
`K_np_n(B)->infinity`, `q_n/p_n(B)->0` for phase-count recovery; the
subexponential critical window is not universal.

The normalized maxima nevertheless differ by `delta+o(1)`.  This also gives
a precise robust minimax statement: if an adversarial `L^infinity[0,B]`
pressure oracle has radius at least half the pressure-curve separation and
`K_nd_n(B)->0`, the minimax absolute error in the normalized maximum is at
least `delta/2-o(1)`.

For the declared marked-phase pressure and mass queries, `(q_n,E_n)` is an
exact presented carrier.  Under hard conjunctive product its extensive
rarity cost `-log q` and reward add.  The same identities hold over a
background with carrier `(\beta\mapsto q_K(\beta),E)`.  This state is a strict
quotient for code-membership and hard-CSP bonuses, but not for labelled
geometric futures or ordinary additive Hamiltonian composition, which
create intermediate phases.

The entropy--reward balance and phase-count product specialize earlier
Theorems 27.2--27.3.  The new operational content is the exact adaptive
transcript bound and sharp replica-budget exponent.  It is a scoped Level-3
theorem/no-go, not a new adversarial-statistical-mechanics architecture and
not a Boolean near-minimizer implication.  Proof and post-repair independent
audit are in
[`drafts/thermal_activation_replica_budget.md`](drafts/thermal_activation_replica_budget.md)
and
[`drafts/thermal_activation_replica_budget_independent_audit.md`](drafts/thermal_activation_replica_budget_independent_audit.md).
