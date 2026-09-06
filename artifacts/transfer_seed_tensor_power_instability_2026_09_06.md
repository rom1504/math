# Complete even self-tensor classification for full symmetric sign seeds

Date: 2026-09-06. Seed-transfer track. The main theorem and its application
to the actual strict-upper weave seeds are proved. The root and the
adversarial researcher independently reconstructed the elementary proof.

The even-power normalized cap tends to one half for a fixed symmetric
Hadamard seed, and to infinity for every other fixed full symmetric sign
seed of order at least two. This is a scalable obstruction for ACTUAL sign constructions, not merely
an obstruction to a Gaussian upper certificate. It does not concern
padding by an external Hadamard, one-time tensoring, or all-order
convergence of the original minima.

## 1. Explicit strict vector witness from a nonorthogonal signing

Let `B` be a full symmetric sign matrix of fixed order `d>=2`; its
diagonal entries are signs as well. Set `G_0=B^TB`. Define

```math
u_i=\frac{B_{i,:}}{\sqrt d},\qquad
v_j=\frac{(G_0)_{:,j}}{\|(G_0)_{:,j}\|},\qquad
V(B)=\frac1{\sqrt d}\sum_j\|(G_0)_{:,j}\|.
```

All vectors are unit. The denominator defining `v_j` is nonzero since
`(G_0)_jj=d`. Direct summation gives

```math
\sum_{i,j}B_{ij}\langle u_i,v_j\rangle=V(B).
```

Each column of `G_0` has its diagonal entry `d`, so its norm is at least
`d`. Equality in every column is equivalent to `G_0=dI`. Consequently

```math
V(B)\ge d^{3/2},\qquad
V(B)>d^{3/2}\ \Longleftrightarrow\ B\text{ is not Hadamard}. (1)
```

Here 'Hadamard' means precisely `B^TB=dI`; no Boolean eigenvector
assumption is included in the term.

There is an explicit Gram-mass lower bound. Let
`Omega(B)=sum_(i!=j)(G_0)_ij^2`. Since
`||(G_0)_:j||<=sqrt(d)||B||op<=d sqrt(d)`,

```math
\frac{V(B)}{d^{3/2}}-1
\ge \frac{\Omega(B)}{d^3(\sqrt d+1)}.                      (2)
```

Indeed rationalize each difference `||(G_0)_:j||-d` and sum. The exact
quantity `V(B)` is usually much sharper than (2), but any nonzero
off-diagonal Gram mass supplies a strictly positive exponential base.

## 2. Tensor even powers make every Gaussian rounding correction positive

Write `c_ij=<u_i,v_j>`. For an EVEN positive integer `r`, tensor these
unit vectors over `r` copies and round them with a common standard
Gaussian vector `g`:

```math
f_{i_1,\ldots,i_r}=\operatorname{sign}
       \langle g,u_{i_1}\otimes\cdots\otimes u_{i_r}\rangle,
\quad
h_{j_1,\ldots,j_r}=\operatorname{sign}
       \langle g,v_{j_1}\otimes\cdots\otimes v_{j_r}\rangle.
```

The elementary bivariate Gaussian sign identity is
`E sign(U)sign(V)=(2/pi)arcsin(corr(U,V))`. Expand
`arcsin(z)=sum_(k>=0) a_k z^(2k+1)`, where
`a_k=binom(2k,k)/(4^k(2k+1))>0` and `a_0=1`. The series is absolutely
convergent even at `z=+-1`. Tensor factorization gives exactly

```math
\mathbb E f^T B^{\otimes r}h
=\frac2\pi\sum_{k\ge0}a_k
       \left[\sum_{i,j}B_{ij}c_{ij}^{2k+1}\right]^r
\ge\frac2\pi V(B)^r.                                    (3)
```

Every bracket is real and `r` is even, so all omitted terms are
nonnegative. This is the reason to keep the EVEN-power qualification.
No quantum XOR multiplicativity theorem or Grothendieck inequality is
needed. Equation (3) implies an actual Boolean bilinear witness by
expectation.

## 3. Conversion to the original same-spin hollow cap

Let `N=d^r`, `C=B^(otimes r)`, and `A=C-diag(C)`. Then `A` is a genuine
symmetric hollow sign matrix. For any symmetric hollow matrix,
`beta(A)=max_(x,y)|x^TAy|<=4Q(A)`, by polarization with the disjoint
vectors `(x+y)/2` and `(x-y)/2` in the unit cube. Also
`beta(C)<=beta(A)+N`. Combining with (3) gives

```math
\boxed{\quad
\frac{Q(\operatorname{hollow}(B^{\otimes r}))}{d^{3r/2}}
\ge\frac1{2\pi}\left(\frac{V(B)}{d^{3/2}}\right)^r
                  -\frac1{4d^{r/2}},\qquad r\text{ even}.
\quad}                                                     (4)
```

Therefore the normalized caps of the even tensor powers tend to
INFINITY for every fixed non-Hadamard full sign seed. A signing which
is good at one order cannot simply be self-tensored indefinitely.

This is not a statement about the sequence `M_N`: at each tensor order,
other signings may have much smaller caps. It also does not say that
every seed with a finite sub-half cap is non-Hadamard; `H_2` is already
a counterexample to such an assertion.

### The Hadamard branch gives the exact other alternative

If `B` is symmetric Hadamard, the Boolean vector `vec(B)` satisfies

```math
(B\otimes B)\operatorname{vec}(B)
=\operatorname{vec}(B^3)=d\operatorname{vec}(B).
```

For every even `r`, tensoring this witness `r/2` times gives a Boolean
eigenvector of `B^(otimes r)` at its positive spectral endpoint
`d^(r/2)`. Its full cap is exactly `d^(3r/2)/2`; hollowing changes a
quadratic value by at most `d^r/2`. Therefore

```math
\lim_{r\to\infty,\ r\ {m even}}
\frac{Q(\operatorname{hollow}(B^{\otimes r}))}{d^{3r/2}}
=\frac12.
```

Together with (4), this is the COMPLETE fixed-seed even-power
classification. No fixed full symmetric sign seed yields a strict
sub-half asymptotic coefficient by literal self-tensoring. This includes
Hadamard seeds with a strict sub-half cap at their initial finite order.

## 4. Every diagonal completion of the strict-upper restricted-weave seeds

The actual strict-upper construction starts with a symmetric Hadamard
`W` of order `M=m^2`, so `W^2=m^2 I`. It keeps `N=mk` coordinates, with
`p_m=N/M=k/m -> p=31/32`, to produce a full sign principal matrix `K`.
Its hollowing is `A=K-diag(K)`.

Partition `W` according to retained and removed coordinates:

```math
W=\begin{pmatrix}K&L\\L^T&J\end{pmatrix},\qquad
K^2+LL^T=m^2 I_N.
```

Since `rank(LL^T)<=M-N`, at least `2N-M` singular values of `K` equal
`m`. Now let `B=A+D` be ANY sign-diagonal completion. Then
`B=K+E`, where `E` is diagonal and `||E||op<=2`. Thus

```math
\|B\|_{op}\le m+2,
\quad\text{and at least }2N-M\text{ singular values of }B
\text{ are at least }m-2.                                 (5)
```

In particular every completion is non-Hadamard once
`m-2>sqrt(N)=m sqrt(p_m)`. This eventually holds because `p<1`.
The conclusion is proved for this principal-restriction family; it is
not inferred merely from having a sub-half Boolean cap.

## 5. A UNIFORM exponential base for that actual family

Let `G=B^TB/N`, so `G_jj=1` and `tr G=N`. Put

```math
\alpha=\frac{2N-M}{N}=2-\frac1{p_m},\qquad
a_-=(m-2)^2/N,\qquad a_+=(m+2)^2/N.
```

For sufficiently large orders, `a_->1` and `0<alpha<1`. Equation (5)
says that a fraction at least `alpha` of the eigenvalues of `G` are
at least `a_-`, while all eigenvalues are at most `a_+` and their
average is one. Cauchy--Schwarz on that portion and its complement
(or minimizing their two average squares at fixed total mean) gives

```math
\frac1N\operatorname{tr}G^2-1
\ge\frac\alpha{1-\alpha}(a_--1)^2.                         (6)
```

The normalized explicit vector value is
`V(B)/N^(3/2)=N^(-1)sum_j sqrt((G^2)_jj)`. Each `(G^2)_jj` lies
between one and `a_+`, since `G^2<=a_+G` and `G_jj=1`. Rationalizing
the square roots in this average and applying (6) gives

```math
\frac{V(B)}{N^{3/2}}-1
\ge\frac{\alpha(a_--1)^2}
             {(1-\alpha)(\sqrt{a_+}+1)}.
```

Uniformly over every signing and diagonal completion in the restricted
family, the right side tends to

```math
\kappa(p)=\frac{(1-p)(2p-1)}{p^2(1+1/\sqrt p)}.
```

At `p=31/32`, this is approximately `0.0154848541950232`, in particular
strictly greater than `0.0154`. For a simple exact lower check, replace
`1/sqrt(p)` by the rational upper bound `51/50`; this still gives
`kappa(p)>77/5000` by rational arithmetic.

Consequently every sufficiently large seed from the strict-upper family,
with EVERY sign-diagonal completion, obeys
`V(B)/N^(3/2)>5077/5000`. Formula (4) at the single EVEN power `r=76`
then yields

```math
\liminf_{N\to\infty}
\frac{Q(\operatorname{hollow}(B_N^{\otimes76}))}{N^{114}}
\ge\frac1{2\pi}\left(\frac{5077}{5000}\right)^{76}>0.508.
```

The final strict comparison can be checked with integers by using
`pi<22/7`: `(7/44)(5077/5000)^76>508/1000`.
The original seed family has normalized caps at most
`0.499432220485404+o(1)`, even after completing its diagonal. Thus this
fixed tensor power produces a genuine, uniformly separated cap increase
for actual good signings. No optimizing-seed assumption is involved.

## 6. A fixed positive random retention does not repair unbounded self-powers

Fix the non-Hadamard seed `B` and `0<p<=1`, and let the even tensor
power `r` tend to infinity. Choose a deterministic maximizing-sign witness
`x_r` for its hollowing `A_r`, of energy at least the right side of (4)
before normalization. Let `T_r` be a uniform random subset of exact size
`floor(p d^r)`.

For a fixed witness on `N` vertices, with `pi_j=(k)_j/(N)_j`, the
exact-size restriction variance is

```math
\operatorname{Var}\left(\sum_{i<j}A_{ij}x_ix_j1_{i,j\in T}\right)
=(\pi_2-2\pi_3+\pi_4)\sum_{i<j}A_{ij}^2
 +(\pi_3-\pi_4)\|Ax\|^2
 +(\pi_4-\pi_2^2)P_A(x)^2
\le N^2/2+N\|A\|_{op}^2.                                  (7)
```

This identity follows by classifying pairs of edges by their intersection
size. The final coefficient is nonpositive; the other two are inclusion/
exclusion probabilities bounded by one. Its mean is `pi_2 P_A(x)`.

Here `||A_r||op<=||B||op^r+1`. Dividing (7) by the square of the positive
inherited witness energy gives exponentially vanishing terms bounded by
constant multiples of

```math
\left(\frac{d\|B\|_{op}^2}{V(B)^2}\right)^r
 +\left(\frac{d^2}{V(B)^2}\right)^r.
```

Both bases are strictly below one: `||B||op<=d` and
`V(B)^2>d^3`. Consequently the inherited witness is relatively
concentrated, even though its normalized energy is diverging. Thus

```math
\frac{Q((A_r)_{T_r})}{|T_r|^{3/2}}
\ge\left(\frac{\sqrt p}{2\pi}-o_{\mathbb P}(1)\right)
          \left(\frac{V(B)}{d^{3/2}}\right)^r
\longrightarrow\infty\quad\text{in probability}.           (8)
```

So a fixed positive random retention AFTER unbounded literal self-tensoring
does not repair this instability. This conclusion is about typical
selectors; specially chosen exceptional restrictions and retentions
tending to zero are not covered.

## 7. Scope and next alternatives

Literal self-tensoring is ruled out, and a constant small diagonal change
cannot repair it for the existing strict-upper seeds. External Hadamard
padding, thin or optimized restrictions chosen AFTER a tensor, switching/gating, and a
different seed-transfer construction remain separate questions.
Section6 treats a fixed positive random retention after unbounded
self-powers. No blanket claim about arbitrarily thin or specially
optimized restrictions is made by this theorem.

## 8. Direct same-spin rounding, without bilinear polarization

There is also a direct same-spin version of the mechanism in Section2.
Keep its unit vectors and set `kappa=max_i |<u_i,v_i>|<=1`. At even
power `r`, write `u_I=otimes_l u_(i_l)`, `v_I=otimes_l v_(i_l)`, and

```math
w_I=\frac{u_I+v_I}{\sqrt{2(1+\kappa^r)}}.
```

Every `w_I` has norm at most one. Add a private orthogonal coordinate
to each one to make it unit, then round all of them with a common
standard Gaussian. For distinct `I,J`, the resulting correlation is
that of the unextended vectors. Write
`a_ij=<u_i,u_j>`, `b_ij=<v_i,v_j>`, and `c_ij=<u_i,v_j>`.
Their correlation before the private extension is

```math
\langle w_I,w_J\rangle
=\frac{\prod_l a_{i_lj_l}+\prod_l b_{i_lj_l}
       +\prod_l c_{i_lj_l}+\prod_l c_{j_li_l}}
       {2(1+\kappa^r)}.
```

Expand the arcsine and then each power by the multinomial theorem.
Each monomial's FULL contraction with `B^(otimes r)` is a nonnegative
coefficient times a real scalar raised to even power `r`. Thus all such
contractions are nonnegative. The two cross terms in the linear arcsine
term contribute `V(B)^r/(1+kappa^r)`. Removing the diagonal from the
full contraction costs at most `pi d^r/2` before the factor `1/pi`.
Hence the direct same-spin witness gives

```math
Q(\operatorname{hollow}(B^{\otimes r}))
\ge\frac{V(B)^r}{\pi(1+\kappa^r)}-\frac{d^r}{2},
\qquad r\text{ even}.                                    (9)
```

Absolute convergence of the arcsine series on `[-1,1]` justifies the
termwise contraction. The private coordinates affect only the diagonal,
which is absent from the objective. If `kappa<1`, this asymptotically
improves the leading prefactor in (4) by a factor of two; the exponential
growth mechanism is the same. Equation (9) is useful conceptually
because the witness is same-spin from its construction, without passing
through a bilinear norm or an absolute PSD majorant.
