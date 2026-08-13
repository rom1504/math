# Symmetric Grothendieck/Krivine audit for quadratic signings

## Status

**New result: a scalable route falsifier, not a bound improvement.** A broad,
precisely defined class of *canonical same-map Gaussian/Krivine roundings* has
an exact leading ceiling on symmetric conference inputs. The class allows
the Gaussian dimension, tensor degrees, Boolean partition, mixtures, and all
coefficients to depend on the order. Its undoubled expected-energy constant
is at most

```math
\frac1\pi=0.318309886\ldots,
```

strictly below the already proved project constant

```math
\frac{c_*}{2}=0.336493364431\ldots .
```

Thus the recent Grothendieck work changes the search in a precise negative
way: its asymmetric left/right correlation engineering should **not** be
pursued by simply using the same Gaussian map at both endpoints of the
canonical conference Gram representation. The theorem below remains uniform
for growing dimension and growing Hermite degree, so it is stronger than the
observation that a fixed cubic term is lower order.

The theorem does **not** cover arbitrary matrix-dependent SDP embeddings,
field--plus--spin rounding, dependent/Onsager probes, asymmetric partitions
followed by a new same-switch coupling, or nonlocal postselection. It is not
a theorem about every conceivable method called “symmetric Grothendieck.”
No rigorous lower bound for the project is improved here.

The two primary sources were downloaded as PDFs and source archives and
inspected in full. Their interval-arithmetic programs were not independently
rerun, so their computer-assisted certificates are reported as paper results
rather than independently certified project results.

## 1. Normalizations and the exact polarization loss

For a symmetric zero-diagonal sign matrix \(A\), distinguish

```math
M(A)=\max_{z\in\{\pm1\}^n}
\left|\sum_{i<j}a_{ij}z_iz_j\right|,
\qquad
Q(A)=\max_z|z^{\mathsf T}Az|=2M(A),
```

and

```math
B(A)=\max_{x,y\in\{\pm1\}^n}|x^{\mathsf T}Ay|.
```

Then

```math
\boxed{Q(A)\le B(A)\le2Q(A)=4M(A).}             \tag{1}
```

The first inequality uses \(x=y\). For the second, set
\(u=(x+y)/2\) and \(v=(x-y)/2\). Symmetry gives

```math
x^{\mathsf T}Ay=u^{\mathsf T}Au-v^{\mathsf T}Av. \tag{2}
```

Both \(u\) and \(v\) lie in \(\{-1,0,1\}^n\). If \(w\) is any such ternary
vector, fill each zero coordinate independently by a uniform sign. Because
\(A\) has zero diagonal, the expected quadratic energy of the completed cube
point is \(w^{\mathsf T}Aw\), so

```math
|w^{\mathsf T}Aw|\le Q(A).
```

Applying this twice in (2) proves (1). Consequently,
\(2M(A)\le B(A)\le4M(A)\) in the user-supplied normalization is correct;
writing that \(M\) as the ledger's doubled \(Q\) would introduce a factor-two
error.

A direct Rademacher calculation already shows the limitation of the bilinear
detour. For uniform \(x\), choose
\(y_i=\operatorname{sgn}(Ax)_i\). Every row sum has the law of a sum of
\(n-1\) independent signs, hence

```math
B(A)\ge \mathbb E_x\|Ax\|_1
=\left(\sqrt{\frac2\pi}+o(1)\right)n^{3/2}.
```

Equation (1) gives only

```math
M(A)\ge\left(\frac1{\sqrt{8\pi}}-o(1)\right)n^{3/2}
=\left(0.1994711402\ldots-o(1)\right)n^{3/2}.    \tag{3}
```

This is below \(0.3364933644\ldots\). Merely substituting the new numerical
upper bound for \(K_G\) into a bilinear SDP comparison cannot remove the
polarization loss and is still weaker.

## 2. What the August 2026 Grothendieck papers prove

For odd measurable signs \(f,g:\mathbb R^d\to\{\pm1\}\), let \(X,Y_t\)
be standard Gaussian vectors with cross-covariance \(tI_d\), and use the
paper's normalization

```math
H_{f,g}(t)=\frac\pi2\mathbb E[f(X)g(Y_t)]
=b_1t+b_3t^3+b_5t^5+\cdots .                    \tag{4}
```

The lower-bound paper proves the affine Hermite strip

```math
\boxed{b_3\ge2b_1-\frac{11}{6}.}                \tag{5}
```

With \(P_r\) denoting Gaussian Hermite projection,

```math
b_1=\frac\pi2\langle P_1f,P_1g\rangle,
\qquad
b_3=\frac\pi2\langle P_3f,P_3g\rangle.          \tag{6}
```

The proof writes \(h=(f+g)/2\), \(k=(f-g)/2\). These are disjoint ternary
functions with \(h^2+k^2=1\), and reduces (5) to

```math
2\|P_1h\|_2^2-\|P_3h\|_2^2+\|P_3k\|_2^2
\le\frac{11}{6}\nu^2,
\qquad \nu=\sqrt{\frac2\pi}.                   \tag{7}
```

After rotating \(P_1h\) to one Gaussian coordinate, its agreement term is
controlled by a one-dimensional rearrangement involving

```math
\Delta(x)=x+2(1+x)\log\frac{1+x}{2},
\qquad
\|P_3h\|_2^2\ge\frac{\nu^2}{6}\Delta_+(x)^2.   \tag{8}
```

The disagreement term is controlled by a sharp degree-three fibre inequality
for \(u:\mathbb R\to\{-1,0,1\}\). The remaining scalar cases yield (7).
Affineness in \((b_1,b_3)\) is crucial: (5) survives convex mixtures and
coefficientwise limits.

If

```math
H^{-1}(z)=a_1z+a_3z^3+\cdots,
```

then \(a_1=1/b_1\) and \(a_3=-b_3/b_1^4\). If \(b_1\le11/12\), the linear
inverse term already bounds the inverse-majorant radius by \(11/12\). If
\(b_1>11/12\), (5) forces a positive \(b_3\), and the linear plus cubic
absolute inverse coefficients enforce the same strict barrier. A reworking
of the Naor--Regev optimality construction supplies limiting mixed Krivine
schemes whose admissible radii approach \(\pi/(2K_G)\). Hence

```math
\frac{\pi}{2K_G}\le\frac{11}{12},
\qquad
K_G\ge\frac{6\pi}{11}.                           \tag{9}
```

The companion case-study paper records the broader candidate family

```math
b_3\ge(1+\lambda)b_1-\left(\lambda+\frac56\right),
\qquad \lambda\ge\frac12,                       \tag{10}
```

whose member \(\lambda=1\) is (5). The tangent case \(\lambda=1/2\) would
give \(K_G\ge9\pi/16\), but remains open there. Likewise, the stronger
intermediate values \(27\pi/49\) and \(51\pi/92\) are described as
machine-verified research claims awaiting human verification; they are not
the companion mathematical paper's theorem.

On the upper-bound side, the new construction engineers signed cubic and
quintic correlation coefficients and uses different functions or
preprocessings on the two sides of a bilinear form. That asymmetry is not
cosmetic: it permits negative cross-Hermite products. A same-map
self-correlation has squared Hermite coefficients instead.

## 3. New theorem: the canonical same-map ceiling

### 3.1 The class being ruled out

For each order \(n\), allow an arbitrary dimension \(d_n\). For each
coordinate \(j\le d_n\), allow an odd tensor/Schoenberg self-kernel

```math
\rho_{j,n}(t)=\sum_{r\text{ odd}}\alpha_{j,r,n}t^r,
\qquad
\alpha_{j,r,n}\ge0,
\qquad
\sum_r\alpha_{j,r,n}\le1.                       \tag{11}
```

This includes direct sums of normalized odd tensor powers; the coefficients
are squared feature norms. From unit vectors \(v_i\), make independent
Gaussian channels \(G_i=(G_{i1},\ldots,G_{id_n})\) satisfying

```math
\mathbb E[G_{ij}G_{kj}]
=\rho_{j,n}(\langle v_i,v_k\rangle),
\qquad
\mathbb E[G_{ij}G_{k\ell}]=0\quad(j\ne\ell).    \tag{12}
```

Finally choose any odd measurable common partition

```math
f_n:\mathbb R^{d_n}\longrightarrow[-1,1]
```

and output at vertex \(i\) according to \(f_n(G_i)\), with independent local
randomization if the value is not already a sign. Convex mixtures of these
schemes are also allowed. Every ingredient may depend on \(n\). Call this
the **canonical same-map tensor/Krivine class**.

The adjective “canonical” matters below: on a conference matrix we feed the
scheme the natural Gram representation \(I+C/\sqrt{n-1}\). The theorem does
not optimize over all matrix-dependent Gram representations.

### 3.2 Coefficient theorem

Let \(K_n(t)\) be the output self-correlation when the input unit vectors have
inner product \(t\). Then

```math
K_n(t)=\sum_{r\text{ odd}}w_{r,n}t^r,
\qquad
w_{r,n}\ge0,
\qquad
\sum_rw_{r,n}\le1,                              \tag{13}
```

and, uniformly in \(n\),

```math
\boxed{w_{1,n}\le\frac2\pi.}                    \tag{14}
```

Consequently, for every \(0\le t\le1\),

```math
\boxed{
K_n(t)\le\frac2\pi t+\left(1-\frac2\pi\right)t^3.}
                                                               \tag{15}
```

#### Proof

Expand \(f_n\) in the orthonormal multivariate Hermite basis:

```math
f_n(g)=\sum_{\beta\in\mathbb N^{d_n}}
\widehat f_n(\beta)h_\beta(g).
```

Multivariate Mehler gives

```math
K_n(t)=\sum_\beta \widehat f_n(\beta)^2
       \prod_{j=1}^{d_n}\rho_{j,n}(t)^{\beta_j}. \tag{16}
```

All coefficients on the right are nonnegative. Oddness of \(f_n\) means only
odd \(|\beta|\) occur, and oddness of every \(\rho_{j,n}\) makes every
surviving monomial odd. Parseval and (11) give

```math
\sum_rw_{r,n}=K_n(1)\le\mathbb E f_n(G)^2\le1.
```

The coefficient of \(t\) can only come from \(\beta=e_j\) and the linear
term of \(\rho_{j,n}\), so

```math
w_{1,n}
=\sum_j\widehat f_n(e_j)^2\alpha_{j,1,n}
\le\|P_1f_n\|_2^2.                              \tag{17}
```

Let \(c=\mathbb E[Gf_n(G)]\); then
\(\|P_1f_n\|_2^2=\|c\|_2^2\). If \(c\ne0\), put \(u=c/\|c\|\). Since
\(|f_n|\le1\),

```math
\|c\|_2
=\mathbb E[(u\mathbin{\cdot}G)f_n(G)]
\le\mathbb E|u\mathbin{\cdot}G|
=\sqrt{\frac2\pi}.                              \tag{18}
```

This proves (14). For \(0\le t\le1\), every remaining odd degree is at least
three, so

```math
K_n(t)\le w_{1,n}t+(1-w_{1,n})t^3.
```

The right side is increasing in \(w_{1,n}\), giving (15). Convex mixtures
preserve every inequality. Notice that no fixed-dimension limit and no
Taylor remainder estimate was used. This is why the theorem remains valid
for \(d_n\to\infty\), unbounded tensor degree, and \(n\)-dependent maps.

### 3.3 Conference corollary

Let \(C\) be a symmetric conference matrix of order \(n\): it is symmetric,
zero diagonal, has off-diagonal signs, and

```math
C^2=(n-1)I.
```

Then

```math
G_C=I+\frac{C}{\sqrt{n-1}}                       \tag{19}
```

is positive semidefinite with diagonal one, hence is a Gram matrix of unit
vectors. Set \(\varepsilon=(n-1)^{-1/2}\). For a scheme in the class above,
oddness gives

```math
\mathbb E[x_ix_j]=K_n(\varepsilon c_{ij})
=c_{ij}K_n(\varepsilon).
```

Therefore

```math
\begin{aligned}
\mathbb E[x^{\mathsf T}Cx]
&=n(n-1)K_n(\varepsilon)\\
&\le\frac2\pi n\sqrt{n-1}
+\left(1-\frac2\pi\right)\frac{n}{\sqrt{n-1}}.
                                                               \tag{20}
\end{aligned}
```

The remainder is \(O(\sqrt n)\), not merely \(O(n)\). In the project's
undoubled normalization, the leading expected-energy guarantee is at most

```math
\frac1\pi n^{3/2}+o(n^{3/2}),                   \tag{21}
```

which is below \(c_*n^{3/2}/2\). Infinite Paley families of symmetric
conference orders make the obstruction scalable rather than an isolated
finite example.

Equation (20) is a ceiling on the *rounding certificate* furnished by this
class, not an upper bound on the actual value \(Q(C)\). A rare output can
have larger energy than its mean, and an arbitrary optimizer of another SDP
Gram matrix is outside the statement. Any use of (20) must retain these two
qualifications.

## 4. Why the new affine strip is subleading here

For the canonical conference Gram, the direct edge correlation is

```math
\varepsilon=(n-1)^{-1/2}.
```

If a fixed analytic odd response is

```math
K(t)=k_1t+k_3t^3+k_5t^5+\cdots,
```

its degree contributions after summing \(n(n-1)\) directed edges have scales

```math
n^2\varepsilon=\Theta(n^{3/2}),\qquad
n^2\varepsilon^3=\Theta(n^{1/2}),\qquad
n^2\varepsilon^5=\Theta(n^{-1/2}).               \tag{22}
```

Thus the paper's restriction on \(b_3\) does not by itself constrain the
leading project constant. In the Grothendieck proof it becomes powerful only
after a global optimality transfer identifies the inverse-majorant radius
with \(\pi/(2K_G)\). No corresponding theorem currently says that canonical
same-map schemes are asymptotically optimal for the symmetric same-switch
Seidel minimax problem.

The new upper construction gets leverage by allowing signed correlation
coefficients from different functions \(f,g\). In the same-map case the
Hermite products become squares, which is exactly the positivity behind
(13). Keeping the asymmetric maps returns to independent left/right signs;
recovering one common switch invokes (1), unless a genuinely new coupling
mechanism is proved.

The only visible ways around the scale argument are:

1. an \(n\)-dependent singular response outside the nonnegative same-map
   coefficient class;
2. a matrix-dependent or nonlocal selector whose performance is not captured
   by the pairwise mean response;
3. dependent cavity/Onsager probes whose endpoint state is not an iid
   Gaussian noise channel; or
4. a new symmetric optimality/duality transfer that changes the relevant
   global functional.

## 5. Audit of the existing \(0.336493\ldots\) proof

The current project theorem uses doubled energy and proves

```math
\liminf_{n\to\infty}\frac{Q(A)}{n^{3/2}}
\ge c_*=0.672986728863\ldots,
```

hence

```math
\liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\ge\frac{c_*}{2}=0.336493364431\ldots .          \tag{23}
```

Its local rule is field--plus--spin rounding,

```math
Y_i^\sigma
=\operatorname{sgn}\!\left(
\sigma\frac{(A\xi)_i}{\sqrt{n-1}}+t\xi_i+\tau Z_i
\right),
```

and the direct-edge coefficient after smoothing and Lindeberg replacement is

```math
4\phi_{1+\tau^2}(t)
\left(2\Phi\!\left(\frac{t}{\sqrt{1+\tau^2}}\right)-1\right).
                                                               \tag{24}
```

Sending \(\tau\downarrow0\) and optimizing at
\(t_*=0.876902\ldots\) gives \(c_*\).

The leading-constant audit is:

- \(Q=2M\) creates the final factor \(1/2\), but this is an identity, not an
  analytic loss.
- Comparing two orientations with the quadratic range uses
  \(\max x^{\mathsf T}Ax-\min x^{\mathsf T}Ax\le2Q(A)\). This is a possible
  centering loss for any asymmetric orientation method, but all factors are
  already present in (24) and (23).
- The spectral bootstrap, smoothing, Lindeberg replacement, and
  \(\tau\downarrow0\) contribute \(o(n^{3/2})\), not a further fixed
  constant.
- Within the entire one-probe coordinatewise Boolean class, \(c_*\) is sharp
  (ledger Section 3.13).
- Within every fixed number of independent field--spin replicas, it remains
  sharp (ledger Section 3.20 and the independent check below).

There is therefore no unidentified numerical slack in the final
\(0.336493\ldots\) conversion that the improved value of \(K_G\) can repair
by substitution.

## 6. Independent check of the multireplica ceiling

Let

```math
G\sim N(0,I_d),\qquad
S\sim\operatorname{Unif}\{\pm1\}^d,
```

independently, and let \(f(G,S)\in[-1,1]\). Put

```math
a=\mathbb E[Gf],\qquad b=\mathbb E[Sf].
```

The first-order direct-edge response of a fixed-\(d\) independent-replica
rule is \(2a\mathbin{\cdot}b\). The exact dimension-free inequality is

```math
\boxed{
2|a\mathbin{\cdot}b|
\le\max_{t\ge0}4\phi(t)(2\Phi(t)-1)=c_*.}        \tag{25}
```

This was already proved in ledger Section 3.20; the present audit checked it
independently. A short proof is as follows. Define

```math
m(s)=\mathbb E_G f(G,s),\qquad
r^2=\mathbb E_Sm(S)^2.
```

Walsh Bessel gives \(\|b\|\le r\). In the direction of \(a\), Gaussian
centroid rearrangement at fixed mean gives

```math
\|a\|\le\mathbb E_S J(|m(S)|),
\qquad
J(u)=2\phi\!\left(\Phi^{-1}\!\left(\frac{1+u}{2}\right)\right).
```

The function \(v\mapsto J(\sqrt v)\) is concave on \([0,1]\). In the scalar
parameter \(z=\Phi^{-1}((1+\sqrt v)/2)\), the sign of its second derivative
is exactly the elementary inequality

```math
2\Phi(z)-1\ge2z\phi(z),
```

whose difference has derivative \(2z^2\phi(z)\ge0\). Jensen therefore gives
\(\|a\|\le J(r)\), and

```math
2|a\mathbin{\cdot}b|\le2rJ(r)\le c_*.
```

Equality ignores all but one replica and uses
\(f(G,S)=\operatorname{sgn}(G_1+t_*S_1)\). This corroborates the existing
result; it is not a new frontier theorem. It also explains why the
field--plus--spin construction lies outside the same-map theorem in Section
3: the local state contains a Rademacher companion with an oriented
cross-correlation to the opposite endpoint, not just copies of one Gaussian
self-correlation channel.

## 7. Endpoint-operator check and the remaining obstruction

For one *idealized isolated edge*, the conditional-expectation channel
between a normalized Gaussian coordinate \(p\) and its opposite endpoint
spin \(q\) has the four eigenvalues

```math
1,\qquad \rho,\qquad-\rho,\qquad\rho^2
```

on the basis \(1,p,q,pq\), with
\(\rho=\varepsilon+O(\varepsilon^5)\) in the conference-scale normalization.
Tensoring independent replicas would make every orientation-odd term beyond
first level \(O(\rho^3)\). That algebra is correct under the isolated iid
channel assumptions, and then (25) is already the stronger sharp
optimization.

It does not close the actual project law. After deletion of one direct edge,
the row fields retain the matrix-dependent correlation profile

```math
q_{ij}=\frac{(A^2)_{ij}}{n-1}.
```

Dependent second-step probes share residuals and backtracking terms and are
not tensor powers of isolated binary-symmetric endpoint channels. A local
operator norm of \(O(\rho^3)\) therefore does not control the sum over the
full signing-dependent correlation kernel. If the number of states grows,
coefficient mass can also compensate a small eigenvalue unless a uniform
norm theorem is supplied.

The exact missing object is a **transport/closure inequality for the
dependent endpoint law**, uniform over the admissible row-correlation
profiles of exact signings. Merely diagonalizing the isolated direct-edge
channel proves a statement in the already-falsified independent-replica
class.

## 8. Research judgment and concrete escape targets

The 2026 Grothendieck result contributes two useful lessons here:

1. affine low-chaos inequalities become globally decisive only when a
   verified optimality transfer makes their coefficients decisive; and
2. tuning inside a response class should stop once a structural coefficient
   obstruction explains the failures.

Equations (13)--(21) give that structural explanation for the canonical
same-map Gaussian/Krivine route. It should be marked inactive. The ordinary
bilinear insertion and fixed independent replicas should also remain
inactive.

The closest genuinely new theorem that would reopen this lower-bound track is
one of:

- **Dependent-channel Hermite inequality.** Define a finite-stage cavity or
  Onsager state including its complete backtracking data, and prove a uniform
  direct-response certificate exceeding \(c_*\) for every signing with
  \(Q(A)=O(n^{3/2})\), with error \(o(n^{-1/2})\) per edge after summation.
- **Same-switch optimality transfer.** Prove that a specified symmetric
  response class is asymptotically optimal for the Seidel quadratic minimax,
  analogously to the Naor--Regev transfer for rectangular bilinear forms.
  Only then could an affine Hermite strip become a project bound.
- **Asymmetric-to-common coupling.** Turn two signed-coefficient left/right
  Krivine partitions into one Boolean vector without the factor-two
  polarization loss, and verify the conversion uniformly for symmetric
  zero-diagonal sign matrices.

The dependent-channel inequality is the most concrete extension of the
current proof architecture. It is also where the present argument has real
unresolved structure rather than a hidden copy of the one-probe optimization.
No theorem in the inspected papers supplies this bridge.

## Sources and audit provenance

- Rahul Saha, Alan Li, Anton Xue, Swarat Chaudhuri, Adam Klivans, Pravesh K.
  Kothari, and Raghu Meka, *New Lower and Upper Bounds for the Grothendieck
  Constant*, [arXiv:2608.11158](https://arxiv.org/abs/2608.11158).
- Alan Li, Rahul Saha, Anton Xue, Swarat Chaudhuri, Adam Klivans, Pravesh K.
  Kothari, and Raghu Meka, *Long-Horizon AI Research for Grothendieck
  Constant: A Case Study in Human--AI Mathematical Collaboration*,
  [arXiv:2608.11195](https://arxiv.org/abs/2608.11195).

The downloaded PDFs, source archives, and extracted text used for this audit
are under /home/math/quadra/tmp/grothendieck_2608/, in compliance with the
repository's no-/tmp rule. Equations (1), (13)--(21), (25), and the mapping
claims were independently derived or checked here. The primary paper's
one-dimensional interval certificate and upper-bound quadrature certificate
were inspected but not rerun.
