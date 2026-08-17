# Exact cap `1/2` does not compress switching response: a linear-rate theorem

Status: rigorous theorem.  The proof combines the exact weighted-deficit
identity BC.1 with three standard quantitative inputs (Rademacher
Hanson--Wright, a uniform two-dimensional Berry--Esseen theorem, and the
spectral norm tail for an iid sign matrix).  All scale and union-bound
inequalities are displayed below.  This is a response-information theorem,
not a convergence result for the original signing problem.

## 1. Statement

For a hollow symmetric signing `A`, write

```math
H_A(x)=\frac12x^TAx,
\qquad
Q(A)=\max_x|H_A(x)|,
```

and, for a sign bridge `B`, define its Boolean response

```math
(P_BH_A)(y)=\max_x\{H_A(x)+x^TBy\}.
```

The projective response distance is

```math
d_{\rm proj}(R,S)
=\frac12\operatorname{osc}_{y\in\{-1,1\}^n}(R(y)-S(y)).
```

### Theorem BCL.1 (bounded-cap linear response rate)

There are absolute constants

```math
c,\gamma,d>0
```

such that, for every sufficiently large

```math
n=2^{2m},
```

there are a dense sign bridge `B_n in {+-1}^{n times n}` satisfying

```math
\|B_n\|_{2\to2}\le c\sqrt n,                            \tag{BCL.1}
```

and at least `exp(gamma n)` hollow symmetric signings `A_1,...,A_N`
such that

```math
Q(A_j)=\frac12n^{3/2}                                   \tag{BCL.2}
```

and

```math
d_{\rm proj}(P_{B_n}H_{A_j},P_{B_n}H_{A_k})
\ge d n^{3/2}
\qquad(j\ne k).                                         \tag{BCL.3}
```

Consequently, let an encoder map every constructed child to one of `K`
states, and let one common decoder map `(state,y)` to a predicted response.
If, for some fixed `epsilon_resp<d/2`,

```math
\sup_y|\widehat R(\operatorname{Enc}(A_j),y)
          -(P_{B_n}H_{A_j})(y)|
\le \epsilon_resp n^{3/2}\qquad(1\le j\le N),            \tag{BCL.3a}
```

then two children cannot share a state: otherwise their projective distance
would be at most `2 epsilon_resp n^(3/2)`.  Thus `K>=N`, and the summary
needs

```math
\exp(\Omega(n))\text{ states},
\qquad \Omega(n)\text{ bits}.                           \tag{BCL.4}
```

All children are switchings of one regular-Walsh signing.  Hence isolated
cap, spectrum, and every switching-invariant statistic are identical across
the family.  The linear information is a genuinely contextual gauge exposed
by the bridge.

The scope of (BCL.4) is the declared response problem `y -> P_BH_A(y)`.
It also lower-bounds any richer continuation class that can pin an external
Boolean query `y` (for example by allowing sufficiently strong coordinate
fields).  It is not automatically a lower bound for a context class that
forbids pinning fields or restricts every appended coefficient to a fixed
small alphabet.

### Corollary BCL.1a (sharp switching-orbit response rate)

For the known regular-Walsh base child and the bridge furnished by Theorem
BCL.1, the worst-case fixed-length response-description complexity of its
entire switching orbit, at every fixed error below `d n^(3/2)/2`, is

```math
\Theta(n)\text{ bits}.                                  \tag{BCL.4a}
```

The theorem supplies the lower bound on a subfamily.  For the upper bound,
encode the switch `s` modulo the immaterial global sign; this takes `n-1`
bits and identifies `A^s` exactly.  A decoder knowing the public base child
and bridge can then evaluate its response (with no claim of computational
efficiency).  Thus the theorem determines an information rate, not only a
one-sided obstruction.

The bridge argument is not specific to quadratic forms.  It proves the
following general response-amplification statement.

### Theorem BCL.0 (extremal-entropy-to-contextual-rate amplifier)

Fix `kappa,d_0>0`.  There are constants `C,gamma,d>0`, depending only on
`kappa,d_0`, with the following property.  Let

```math
H:\{-1,1\}^n\longrightarrow\mathbb R,
\qquad P=\max_u H(u),
```

and suppose

```math
\#\{u:P-H(u)<d_0n^{3/2}\}
\le \exp((\log2-\kappa)n).                              \tag{BCL.0a}
```

For `s in {+-1}^n`, let `H^s(x)=H(s odot x)`.  Then, for all sufficiently
large `n`, there are a sign bridge `B` with
`||B||_(2->2)<=C sqrt(n)` and a set `S` of at least `exp(gamma n)` switches
such that

```math
d_{\rm proj}(P_BH^s,P_BH^t)\ge d n^{3/2}
\qquad(s\ne t\text{ in }S).                             \tag{BCL.0b}
```

The proof is Sections 4--7 below verbatim.  The only facts about `H` used
there are a chosen maximizer `u_*`, the entropy gap (BCL.0a), and the exact
weighted-deficit identity (BCL.27), which holds for every function on the
Boolean cube.  Sections 2--3 verify (BCL.0a) with absolute constants for the
regular-Walsh quadratic landscape, while switching preserves its exact cap;
Theorem BCL.1 follows.

This is a genuine implication from an unrooted extremal statistic to rooted
future complexity:

```math
\text{positive near-top entropy deficit}
\quad\Longrightarrow\quad
\text{positive contextual response-information rate}.  \tag{BCL.0c}
```

It does not say that extremal entropy is a sufficient compositional state.
On the contrary, it shows that even a one-number entropy gap can certify
that the hidden switching root must become linearly observable under a
suitable low-operator-norm interaction.

### Corollary BCL.0a (bounded-operator quadratic landscapes)

Fix `p,L>0`.  Let `K_n` be symmetric and put

```math
H_n(u)=\frac12u^TK_nu,
\qquad \bar H_n=\mathbb E_U H_n(U).
```

If

```math
\|K_n\|_F\le Ln,
\qquad
\|K_n\|_{2\to2}\le L\sqrt n,
\qquad
\max_uH_n(u)-\bar H_n\ge p n^{3/2},                    \tag{BCL.0d}
```

then the switched family of `H_n` has a positive contextual response-
information rate in the sense of (BCL.0b), with constants depending only
on `p,L`.

Indeed, take `d_0=p/2`.  A spin within `d_0 n^(3/2)` of the top has
`H_n(U)-bar H_n >=p n^(3/2)/2`; Hanson--Wright and (BCL.0d) give
(BCL.0a) with an absolute positive `kappa(p,L)`.  Theorem BCL.0 applies.
This includes bounded-operator SK-type deterministic realizations with a
macroscopic ground-state excess, as well as regular Hadamard/conference
quadratics.  No randomness of the child is required.

## 2. The regular-Walsh child

Put `q=2^m`, so `n=q^2`, and let `W` be the Sylvester matrix indexed by
`F_2^{2m}`:

```math
W_{a,u}=(-1)^{a\cdot u}.
```

There is a Boolean vector `b` with `Wb=qb`; explicitly, after writing a
coordinate as `(u,v) in F_2^m times F_2^m`, take

```math
b(u,v)=(-1)^{u\cdot v}.
```

Define

```math
\mathcal H=D_bWD_b,
\qquad
A=\mathcal H-\operatorname{diag}(\mathcal H).           \tag{BCL.5}
```

Then

```math
\mathcal H\mathbf1=q\mathbf1,
\qquad
\mathcal H^2=nI,
\qquad
\operatorname{tr}\mathcal H=0.                         \tag{BCL.6}
```

Thus `A` is a hollow symmetric sign matrix and, for Boolean `u`,

```math
H_A(u)=\frac12u^T\mathcal Hu,
\qquad
P:=\max_uH_A(u)=Q(A)=\frac12qn.                         \tag{BCL.7}
```

We use the top state `u_*=1`.  Every switched child

```math
A^s=D_sAD_s
```

has the same exact cap.

## 3. Only exponentially few spins are near the top

Fix once and for all

```math
d_0=\frac18,
\qquad
T_0=\{u:P-H_A(u)<d_0qn\}.                               \tag{BCL.8}
```

### Lemma BCL.2 (near-top entropy gap)

There is an absolute `kappa>0` such that, for all sufficiently large `n`,

```math
|T_0|\le\exp((\log2-\kappa)n).                          \tag{BCL.9}
```

#### Proof

For uniform Rademacher `U`, membership in `T_0` implies

```math
U^T\mathcal HU>(1-2d_0)qn=\frac34n^{3/2}.              \tag{BCL.10}
```

The quadratic form has mean `tr(mathcal H)=0`, while

```math
\|\mathcal H\|_F=n,
\qquad
\|\mathcal H\|_{2\to2}=\sqrt n.
```

The Rademacher Hanson--Wright inequality therefore bounds the probability
of (BCL.10) by

```math
2\exp\left[-c_{HW}
 \min\left\{\frac{(3n^{3/2}/4)^2}{n^2},
                    \frac{3n^{3/2}/4}{\sqrt n}\right\}
 \right]
\le e^{-\kappa n}                                      \tag{BCL.11}
```

after reducing an absolute `kappa` and taking `n` large.  Multiplication by
`2^n` proves (BCL.9). `square`

The fixed superset `T_0`, rather than a deficit-dependent estimate, is used
below.  This avoids any circular dependence between the Hanson--Wright
constant and the final, much smaller response gap `d`.

## 4. A uniform two-query row lemma

For `r in R`, let

```math
D(\delta\|\alpha)
=\delta\log\frac\delta\alpha
 +(1-\delta)\log\frac{1-\delta}{1-\alpha}              \tag{BCL.12}
```

denote binary relative entropy.

### Lemma BCL.3 (uniform weighted sign disagreement)

For every `epsilon_row>0` there are `a,rho>0` and `n_0` such that the following
holds.  If `n>=n_0`, `y,z in {+-1}^n` obey

```math
|y^Tz|\le\rho n,
```

`R` is a uniform sign vector, and `t in {+-1}`, then, with

```math
S=R^Ty,
\qquad T=R^Tz,
```

one has

```math
\Pr\{|S|\ge a\sqrt n,\ \operatorname{sign}(S)
       \operatorname{sign}(T)\ne t\}
\ge\frac12-\epsilon_row.                               \tag{BCL.13}
```

Any fixed convention at zero may be used.

#### Proof

Put `xi_i=R_i y_i`, and split the coordinates according to
`y_i z_i=+1` and `-1`.  For independent Rademacher sums `U,V` of lengths

```math
\frac{n+y^Tz}{2},
\qquad
\frac{n-y^Tz}{2},
```

respectively,

```math
(S,T)=(U+V,U-V).                                       \tag{BCL.14}
```

The covariance of `(S/sqrt(n),T/sqrt(n))` has unit diagonal and
off-diagonal `theta=y^Tz/n`.  More explicitly, with

```math
X_i=n^{-1/2}\xi_i(1,c_i),
\qquad c_i=y_iz_i,
\qquad
\Sigma_\theta=
\begin{pmatrix}1&\theta\\ \theta&1\end{pmatrix},
```

one has `sum_i Cov(X_i)=Sigma_theta`.  If `|theta|<=rho<1`, then

```math
\sum_i\mathbb E\|\Sigma_\theta^{-1/2}X_i\|_2^3
\le {2^{3/2}\over(1-\rho)^{3/2}\sqrt n}.                \tag{BCL.14a}
```

Bentkus's Lyapunov bound for independent, not necessarily identically
distributed vectors therefore compares `(S/sqrt(n),T/sqrt(n))` with its
centered Gaussian limit, uniformly over convex Borel sets, with error
`C_rho/sqrt(n)`.  The event in (BCL.13) is the disjoint union of two convex
quadrant half-strips, so its total comparison error is at most
`2C_rho/sqrt(n)`.  See V. Bentkus, *A Lyapunov-type bound in R^d*, Theory of
Probability and Its Applications 49 (2005), 311--323,
<https://doi.org/10.1137/S0040585X97981123>.

For the limiting Gaussian pair, the arcsine law gives

```math
\Pr\{\operatorname{sign}(S)\operatorname{sign}(T)=+1\}
=\frac12+\frac{\arcsin\theta}{\pi}.                    \tag{BCL.15}
```

Thus either requested sign has probability at least
`1/2-arcsin(rho)/pi`.  Removing the strip `|S|<a sqrt(n)` costs at most
`Pr{|G|<a}` in the limit, uniformly in `theta` because the `S` marginal is
standard normal.  First choose `rho` and `a` small, then `n` large, so

```math
{\arcsin\rho\over\pi}+\Pr\{|G|<a\}
 +{2C_\rho\over\sqrt n}<\epsilon_row.
```

This also handles either convention at zero: the Gaussian boundaries have
zero mass and the strong `|S|` threshold excludes a zero first coordinate.
`square`

The lemma deliberately controls both target signs.  This is what permits a
union bound over an adversarial near-top spin rather than an averaged target.

## 5. Code and constants

Take `kappa` from Lemma BCL.2.  Choose `epsilon_row>0` so small that, with

```math
\alpha=\frac12-\epsilon_row,
```

```math
D(0\|\alpha)=-\log(1-\alpha)
>\log2-\frac\kappa4.                                   \tag{BCL.16}
```

Choose `delta in (0,alpha)` sufficiently small that

```math
D(\delta\|\alpha)>\log2-\frac\kappa2.                 \tag{BCL.17}
```

Apply Lemma BCL.3, shrinking `a` if necessary so that

```math
d:=2a\delta<d_0.                                        \tag{BCL.18}
```

Finally choose

```math
0<\gamma<\min\left\{\frac{\rho^2}{8},\frac\kappa8\right\}.
                                                                    \tag{BCL.19}
```

A standard random-code argument now gives a fixed set

```math
Y\subset\{-1,1\}^n,
\qquad |Y|\ge e^{\gamma n},                             \tag{BCL.20}
```

with

```math
|y^Tz|\le\rho n\qquad(y\ne z).                         \tag{BCL.21}
```

Indeed, for two random codewords Hoeffding gives probability at most
`2e^{-rho^2 n/2}` of violating (BCL.21), and the union bound over fewer than
`e^{2gamma n}` pairs tends to zero by (BCL.19).

## 6. A random bridge avoids every weighted near-top neighborhood

Choose `B` with independent uniform sign entries.  For `y in Y`, define the
query-linked switching

```math
s_y=u_*\odot\operatorname{sign}(By)
=\operatorname{sign}(By).                              \tag{BCL.22}
```

Fix an ordered pair `y!=z` and a spin `u in T_0`.  For row `i` of `B`, put

```math
S_i=(By)_i,
\qquad T_i=(Bz)_i,
\qquad t_i=u_i(u_*)_i=u_i.
```

The rows are independent.  By Lemma BCL.3, each row has probability at
least `alpha` of satisfying

```math
|S_i|\ge a\sqrt n,
\qquad
\operatorname{sign}(S_i)\operatorname{sign}(T_i)\ne t_i.
                                                                    \tag{BCL.23}
```

The second condition says exactly that `u_i` disagrees with the sign of the
cross field

```math
h=s_z\odot By.
```

The lower Chernoff bound, valid also for independent non-identically
distributed indicators with success probabilities at least `alpha`, gives

```math
\Pr\left\{
 \sum_{i:u_i\ne\operatorname{sign}(h_i)}|h_i|
 <\frac d2n^{3/2}
 \right\}
\le e^{-D(\delta\|\alpha)n}.                            \tag{BCL.24}
```

Here (BCL.18) is used: `delta n` successful rows already contribute
`delta n a sqrt(n)=d n^(3/2)/2`.

Union-bound (BCL.24) over all `u in T_0` and all ordered distinct pairs in
`Y`.  By (BCL.9), (BCL.17), and (BCL.19), the failure probability is at
most

```math
\begin{aligned}
\exp\{[(\log2-\kappa)+2\gamma
       -D(\delta\|\alpha)]n\}
&\le e^{-\kappa n/4}.                                  \tag{BCL.25}
\end{aligned}
```

For an iid sign matrix, the standard net argument also gives

```math
\Pr\{\|B\|_{2\to2}>C\sqrt n\}\le2e^{-c_Bn}            \tag{BCL.26}
```

for absolute `C,c_B>0`.  Hence, for large `n`, there is a deterministic
bridge satisfying (BCL.1) and all the weighted-neighborhood inequalities in
(BCL.24).

## 7. BC.1 converts the neighborhood estimate into response packing

For any field `h`, define its top deficit

```math
\Delta_A(h)
=P+\|h\|_1-\max_v\{H_A(v)+h^Tv\}.
```

The exact identity BC.1 is

```math
\Delta_A(h)
=\min_v\left\{P-H_A(v)
 +2\sum_{i:v_i\ne\operatorname{sign}(h_i)}|h_i|\right\}.
                                                                    \tag{BCL.27}
```

Put `eta=d n^(3/2)`.  If `v notin T_0`, then, because `d<d_0`, its first
term in (BCL.27) is already at least `eta`.  If `v in T_0`, the simultaneous
event furnished by (BCL.24)--(BCL.25) makes the second term at least `eta`.
Therefore, for every ordered distinct `y,z in Y`,

```math
\Delta_A(s_z\odot By)\ge\eta.                          \tag{BCL.28}
```

Let `R_y=P_BH_{A^{s_y}}`.  The query-linked exposure identity gives

```math
R_y(y)=P+\|By\|_1,
```

while BC.1 gives

```math
d_{\rm proj}(R_y,R_z)
\ge\frac12\left[
 \Delta_A(s_z\odot By)+\Delta_A(s_y\odot Bz)
 \right]
\ge\eta.                                               \tag{BCL.29}
```

In particular the switchings in (BCL.22) cannot collide.  There are
`|Y|>=e^(gamma n)` distinct cap-`1/2` children, and (BCL.29) proves Theorem
BCL.1.

## 8. Interpretation

This closes the qualitative bounded-cap packing question.

1. A cap bound at the natural `n^(3/2)` scale does not imply sublinear
   future-response information.  Even the exact regular-Hadamard cap `1/2`
   class contains a fixed-bridge packing with a positive information rate.
2. The mechanism is not coefficient decoding.  All children are gauge
   switchings of one signing.  What the bridge observes is weighted distance
   from an exponentially sparse near-top set.
3. The proof separates bulk entropy from rare extremal geometry: a
   Hanson--Wright entropy gap of only `kappa n` is amplified by a row event
   whose optimal exponent approaches `log 2`.
4. The theorem is stronger than the explicit
   `Omega(sqrt(n) log n)` Maiorana--McFarland code, but less explicit: the
   bridge and query code are obtained probabilistically.
5. It still does not address a response class in which bridge and child are
   quotiented jointly by switching, nor does it prove that arbitrary dense
   composition requires the same linear rate.  It also does not distinguish
   exact minimizers below the conference/regular-Hadamard cap.

The next discriminator is no longer whether bounded-cap children can carry
linear contextual information; they can.  It is whether a **joint
gauge-covariant carrier** can quotient this exposed switching phase while
remaining strictly smaller than the full bridge/child response landscape.

## 9. Standard probability inputs

The exact forms used above are available in the following primary sources.

- M. Rudelson and R. Vershynin,
  ["Hanson--Wright inequality and sub-gaussian concentration"](https://doi.org/10.1214/ECP.v18-2865),
  *Electronic Communications in Probability* **18** (2013), 1--9.  Their
  Theorem 1.1 applies directly to the Rademacher vector in (BCL.10).
- V. Bentkus,
  ["A Lyapunov-type bound in `R^d`"](https://doi.org/10.1137/S0040585X97981123),
  *Theory of Probability and Its Applications* **49** (2005), 311--323.
  This version permits independent, non-identically distributed summands
  and controls every convex Borel set; after covariance normalization it
  gives the uniform `O_rho(n^(-1/2))` error used in Lemma BCL.3.

The iid-sign spectral-norm estimate (BCL.26) also follows directly from the
standard one-quarter-net proof: for fixed unit `x,y`, `y^TBx` is
subgaussian with absolute scale one; a union bound over two nets of size at
most `9^n` gives `||B||<=C sqrt(n)` outside `2e^{-c_Bn}`.  No asymptotic
random-matrix limit is used.
