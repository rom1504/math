# Uniform first-marked local rate and cap-only energy transport

2026-09-06. Director/fresh-agent derivation and archive audit.

The QUALITATIVE local Gaussian law is not new: the independently audited
rooted odd-degree tree theorem in
`fresh_limit_second_rooted_tree_2026_09_05.md`, Sections 1 and 4, already
proves it uniformly at every root under the weaker assumption
`beta(A)=o(n^2)`. The elementary spectral bound used below is also archived
in `fresh_limit_rooted_response_2026_09_05.md`, Section 2.

This note supplies an explicit smooth-test rate with a self-contained
finite-cubic argument, audits the alternative forced-root core argument,
and records the cap-only ENERGY-scale consequence for source approximation.
It does not prove a full untruncated feedback covariance or Gaussian law
for BF, BC, QS, or QD.

## 1. Normalization and the archived sharp elementary bootstrap

Let A be an n-by-n hollow symmetric signing, m=n-1, and

```
B=A/sqrt(m),  Q=B^2,  L=||B||op,
Lambda(B)=max_(x Boolean) |x^T B x|/(2n)<=C,
beta(B)=max_(x,y Boolean) |x^T B y|<=4 C n.
```

For an eigenpair `Bv=lambda v`, multilinearity on the two cubes gives

```
beta(B) >= ||B(v/||v||infty)||1
         = |lambda| ||v||1/||v||infty.
```

The eigen-equation at a largest coordinate gives
`|lambda| ||v||infty <= maxentry(B) ||v||1`. Hence

```
L^2 <= maxentry(B) beta(B) <=4 C n/sqrt(m).           (1)
```

Also `1<=L<=sqrt(m)` and every row of B has squared norm one. No complex
interpolation is required for (1).

For reference, the proposed complex-interpolation alternative is valid
but loses a factor two. The complex infinity-to-one norm is at most
`2 beta_real(B)`: rotate the dual pairing to be real and expand its real
part as two real cube pairings. Complex Riesz--Thorin between (infinity,1)
and (1,infinity), at parameter one half, then gives
`L^2<=2 beta(B) maxentry(B)`. The endpoint theorem was checked in
[MIT's own interpolation lecture, Theorem 22.1 and proof](https://ocw.mit.edu/courses/res-18-015-topics-in-fourier-analysis-spring-2024/mitres_18_015_s24_lec22.pdf).
We use the stronger elementary (1), not that imported alternative.

## 2. Exact old fields and the cubic kernel cuts

On independent Boolean signs S keep

```
G=BS,  D_i=S_i h2(G_i),  Y=BD.
```

Y_i is a pure multilinear cubic. Its symmetric ordered kernel is

```
k_i(a,b,c) = 1_{a,b,c distinct}/(3 sqrt(2))
 * [B_ia B_ab B_ac+B_ib B_ab B_bc+B_ic B_ac B_bc],

Y_i=sum_(a,b,c) k_i(a,b,c) S_a S_b S_c.
```

The factor is exact: summing the six orders of a triple gives its
coefficient `sqrt(2)` times the bracketed sum. Let kappa_i be the operator
norm of any 1-versus-2 flattening of this symmetric kernel. Then

```
kappa_i <=sqrt(2) L/sqrt(m).                         (2)
```

Here is the complete factorization. First omit distinctness and consider
`F_i(a,b,c)=B_ia B_ab B_ac`. Its a-versus-(b,c) flattening is
`diag(B_i) K`, where `K_(a,(b,c))=B_ab B_ac` has Gram matrix `Q circ Q`.
Schur multiplication by the correlation matrix Q is operator-contractive,
so `||K||op<=L`, while `||diag(B_i)||op<=1/sqrt(m)`.
Its b-versus-(a,c) flattening is B times the disjoint-row lift
`J_(a',(a,c))=1_{a'=a}B_ia B_ac`, of norm at most `1/sqrt(m)`.
The third cut is identical after interchange.

The a=b and a=c diagonals already vanish in this unsymmetrized term.
The only remaining exclusion is b!=c. On a-versus-(b,c) it is a column
projection. On either other cut its complement is a rectangular
block-diagonal pinching, of norm at most one, so imposing the inequality
costs at most two. Symmetrizing the three terms with the displayed factor
proves (2). No growing-L constant is hidden in this argument.

The exact covariance identity from the archive is

```
Cov(D)=(1-3/m)I+(2/m)Q,
v_i:=E Y_i^2=1-3/m+(2/m)(B^4)_ii.                  (3)
```

Thus `v_i<=3`, `v_i>=1-1/m`, and
`|v_i-1|<=(3+2L^2)/m`. G has variance one and is exactly orthogonal to Y.

## 3. Boolean-to-Gaussian replacement with an explicit constant

For `V_i=(G_i,Y_i)`, the exact coefficient bound gives

```
sum_a I_a(V_i)<=10,
max_a I_a(V_i)<=10/m.                              (4)
```

Indeed a cubic coefficient is at most `3sqrt(2)/m^(3/2)`, so
`Inf_a(Y_i)<=9(n-2)/m^2`; its total influence is `3v_i<=9`.

Replace the Boolean seeds one at a time by independent standard Gaussians.
At a replacement coordinate, write `V_i=V0+X_a d`; the vector d has degree
at most two and is independent of X_a. Its squared L2 norm is I_a(V_i),
unchanged by any preceding replacements. The elementary fourth-moment
induction works on every Boolean/Gaussian hybrid: in the splitting
`P=P0+X P1`, both laws have mean zero, variance one, zero third moment,
and fourth moment at most three. Therefore

```
(E||d||^4)^(1/4)<=3 sqrt(I_a(V_i)),
E||d||^3<=27 I_a(V_i)^(3/2).
```

For a C3 test phi on R squared let M3 bound the operator norm of its third
derivative. Taylor expansion through order two has matching first two
terms for the two seed laws. Summing the third-order remainder gives

```
|E phi(V_i^Boolean)-E phi(V_i^Gaussian)|
 <=(27/6)[1+E|N|^3] M3 sum_a I_a(V_i)^(3/2)
 <=400 M3/sqrt(m).                                 (5)
```

The constant 400 exceeds
`(27/6)(1+2sqrt(2/pi))10sqrt(10)`. This is a direct finite replacement,
not an appeal to an invariance theorem with unspecified parameters.

## 4. Gaussian cubic brackets and a direct interpolation comparison

Use the same kernels on independent Gaussian seeds. In conventional Wiener
integral notation `G_i=I1(b_i)`, `Y_i=I3(k_i)`, with
`||b_i||2=1` and `6||k_i||HS^2=v_i`. Put

```
Gamma_ab=(1/d_a) <D V_a,D V_b>,  d_G=1, d_Y=3.
```

Gaussian integration by parts in homogeneous polynomials gives
`E[V_a f(V)]=sum_b E[Gamma_ab partial_b f(V)]` exactly. In particular
`E Gamma=diag(1,v_i)`. The product of two second Gaussian chaoses gives

```
Gamma_YY-v_i
 =3 I4(sym(k_i contract_1 k_i))
    +12 I2(k_i contract_2 k_i),
Gamma_YG=I2(k_i contract_1 b_i),
Gamma_GY=3 Gamma_YG,  Gamma_GG=1.
```

These identities follow by grouping zero, one, or two cross contractions
in the product of the two degree-two derivatives. Orthogonality of
different Gaussian chaoses then gives

```
Var(Gamma_YY)
 <=216||k_i contract_1 k_i||HS^2
    +288||k_i contract_2 k_i||HS^2 <=252 kappa_i^2,
Var(Gamma_YG)<=2 kappa_i^2,
Var(Gamma_GY)<=18 kappa_i^2.                         (6)
```

For the last bound on Gamma_YY, both contraction norms squared are bounded
by `kappa_i^2 ||k_i||HS^2<=kappa_i^2/2`: they are the two Gram matrices of
the same flattening. This checks every dimension factor and normalization.

Let Z be an independent Gaussian pair of covariance diag(1,v_i).
Differentiate `E phi(sqrt(t)V+sqrt(1-t)Z)` and use the exact integration
by parts for each vector. The derivative equals

```
(1/2) sum_ab E[(Gamma_ab-Cov(V)_ab)
                    partial_ab phi(sqrt(t)V+sqrt(1-t)Z)].
```

If M2 bounds the Hessian operator norm, integrating and using (6) bounds
the error by
`(M2/2)(sqrt(252)+4sqrt(2)) kappa_i`.
Changing the Gaussian variance v_i to one costs at most
`M2 |v_i-1|/2`, by the same elementary interpolation.
Since `L<=sqrt(m)`, equations (2)--(3) bound their sum by
`32 M2 L/sqrt(m)`.

Combining with (5) proves the uniform estimate

```
sup_i |E phi(G_i,Y_i)-E phi(N0,N1)|
 <=400 M3/sqrt(m)+32 M2 L/sqrt(m)
 <=400 M3/sqrt(m)+64 sqrt(C) M2 n^(1/2)/m^(3/4),     (7)
```

where N0,N1 are independent standard Gaussians. The test may be taken
bounded C3; its first derivative and value need not enter the displayed
constant. The rate is uniform over all original signings with Lambda<=C.

## 5. Rectangles, hard source responses, and uniform row L2

Mollify an axis-aligned rectangle from above and below at width eta.
The tests have Hessian and third-derivative bounds O(eta^-2), O(eta^-3).
The independent Gaussian boundary-strip mass is O(eta), uniformly in all
four endpoints, including infinite endpoints. The resulting sandwich in
(7) gives

```
sup_i sup_(rectangles R)
 |P((G_i,Y_i) in R)-gamma_2(R)|
 <=O_C(eta+n^-1/4 eta^-2+n^-1/2 eta^-3).
```

Taking `eta=n^-1/12` proves `O_C(n^-1/12)` rectangle error. Open, closed,
and half-open conventions have the same estimate by the strip bound.
This is a local statement only; it is not a small-ball estimate for BF.

For every fixed polynomial p and bounded Gaussian-a.e.-continuous f,
the uniform local law and uniform fixed moments imply

```
sup_i |E|f(G_i,Y_i)-p(G_i,Y_i)|^2
                 -E_gamma|f-p|^2| ->0.             (8)
```

The moment input requires no operator cap: G has variance one and degree
one, Y has variance at most three and degree three. Iterating the elementary
fourth-moment estimate on their fixed powers bounds every fixed moment,
uniformly over all rows and signings. This supplies uniform integrability
for the polynomial growth in (8). Approximation degree remains FIXED
before n tends to infinity; (8) is not uniform over increasing degrees.

## 6. Cap-only transport at energy scale, not covariance scale

The Hilbert-space form of Grothendieck gives, for arbitrary random vectors
X,Z on the SAME probability space,

```
|E X^T B Z|/n
 <=4 KG C [sup_i ||X_i||2][sup_j ||Z_j||2].           (9)
```

No independence, centering, or bounded operator norm is used. In
particular if `F_i=f(G_i,Y_i)`, `P_i=p(G_i,Y_i)` and
`||f-p||L2(gamma)<=delta`, then (8)--(9) imply

```
E||B(F-P)||1/n <=4 KG C [delta+o(1)].               (10)
```

To see (10), choose the random Hilbert vectors
`X_i=sign((B(F-P))_i)` in (9). Their individual L2 norms are at most one.
This is genuine full-parent transport, with no discarded-source term.
The uniform-in-root conclusion (8), not a merely averaged L2 conclusion,
is essential to this inference.

There is also a useful two-source version. Suppose |f|,|h|<=1 and choose
fixed polynomials p,q with Gaussian L2 errors at most delta. Set
`F=f(G,Y), H=h(G,Y), P=p(G,Y), H_p=q(G,Y)`. Then

```
|n^-1 E sum_i H_i |(BF)_i|
       -n^-1 E sum_i H_p,i |(BP)_i||
 <=4 KG C delta(2+delta)+o(1).                      (11)
```

For the first term of the difference, apply (9) to
`X=(H-H_p)sign(BF)`, Z=F. For the second, bound its absolute value by
`E sum_i |H_p,i| |(B(F-P))_i|` and apply (9) to
`X=|H_p|sign(B(F-P))`, Z=F-P. The row norms are respectively bounded by
delta, 1, 1+delta, and delta, up to o(1), proving (11).

Equations (10)--(11) do NOT assert that `E||B(F-P)||2^2/n` is small,
nor normalized nuclear covariance convergence. Replacing the hard returned
sign in `C=H sign(BF)` still requires a separate small-ball/joint-comparison
argument for BF. No such unrestricted statement follows merely from the
old pair's uniform Gaussian marginal.

## 7. Audit of the alternative seed-core proof and normalization

Let R be a diagonal-majorant core deleting at most epsilon n vertices.
To study a designated output root i, force it into the core by replacing
R with `R union {i}`. Adding one signing row and column raises the operator
norm by at most sqrt(n). Thus the child, at fixed epsilon, has a fixed
normalized operator bound, uniformly in the designated root.

For m_R=|R|, put `r=(m_R-1)/(n-1)`, and let child fields use the signing
`A_R/sqrt(m_R-1)`. On any retained output root, the EXACT seed Fourier
projections are

```
E[G_i|S_R]=sqrt(r) G_i^child,
E[Y_i|S_R]=r^(3/2) Y_i^child.                       (12)
```

The restricted quadratic center is r, not one: the projected inner field
is `S_j[(sum_(k in R) B_jk S_k)^2-r]/sqrt(2)`.
The influence bounds give, if t seeds were removed,

```
E|G_i-E[G_i|S_R]|^2 <=t/(n-1),
E|Y_i-E[Y_i|S_R]|^2 <=9t(n-2)/(n-1)^2.              (13)
```

Use the archived fixed-operator, uniform-root local law on the child at
fixed epsilon; (12)--(13) cost O(sqrt(epsilon)) in a bounded-Lipschitz
comparison, and the Gaussian rescaling costs O(epsilon). Then send
epsilon to zero. This proves the same qualitative uniform local law.
Without forcing the output root into the core, the argument proves only
an averaged law, which would be insufficient for (9)--(10).

The seed projection is a proof device for the local query law. It does
not substitute child optimization caps for the original parent or incur
a principal retention factor in the energy transport conclusion.

## 8. Scope and archive status

Already archived: uniform qualitative old Gaussian law, rooted-tree
moment mechanism, exact D covariance, and the O(n^1/4) spectral bootstrap.
Recorded here: the explicit smooth/rectangle rate and its uniform-row
source-approximation consequence at cap-only energy scale, together with
the checked seed-projection normalization.

No unrestricted feedback ascent, improved decimal lower bound, or
convergence of the original minima is claimed.

Independent adversarial full audit reported PASS on 2026-09-06: symmetric
cubic normalization/cuts, Gaussian bracket factors and both numerical
comparison constants, rectangle rate, uniform-row L2 passage, cap-only
energy transport, and forced-root normalization were reconstructed.
