# Polynomial BH applied to reduced quadratic powers

Status: **proved; independently reconstructed by root and the paper verifier**, 2026-09-16.
This note uses the normalization
`H_A(x)=sum_{i<j} a_ij x_i x_j`, `Q(A)=||H_A||_infinity`, and
`M_n=min_A Q(A)` for hollow symmetric full sign matrices. It does not use
the doubled `Q` notation in the older `high_order_flat_chaos.md` artifact.

## 1. Primary theorem and proof reconstruction

The complete primary paper
[Ivanisvili, arXiv:2609.12427v1](https://arxiv.org/html/2609.12427v1)
was read before consulting the two requested prior artifacts.
Its claims are a dimension-free bound `B_m <= K m^27` and the weighted
bound

```math
\mathcal W(f):=\left(\sum_{r=1}^m r^{-10}
  \|\widehat f\vert_{|S|=r}\|_{2r/(r+1)}^2\right)^{1/2}
\le K m^{22}\|f\|_\infty.
```

The proof mechanism checks as follows. Work at fixed dimension and take
the supremum `M` of this weighted functional over the bounded degree-`m`
unit ball. Levels `r<=11` cost `O(m^22)` by iterated Markov coefficient
extraction and the finite fixed-degree homogeneous BH constants. Randomly
partition coordinates into two equiprobable colors; for `r>=12`, the
event `r/6<=d,e<=5r/6` captures probability at least
`p=1-2 exp(-8/3)`. Comparing the norm over degree blocks costs
`kappa=13^(1/24)`. Mixed-norm interpolation yields
`b_(d,e)<=X_(d,e)^(d/r)Y_(d,e)^(e/r)`. Inserting
`rho_d=sqrt((d-1)/(d+1))` and its counterpart costs at most `3`.
The weight conversion contributes
`[(d/r)^(d/r)(e/r)^(e/r)]^5 <= (2/3)^5`.
After squaring and AM--GM, each of the two damped weighted sums is at most
`M^2`: Minkowski collects **all** degrees in the other color, and
hypercontractivity with `rho_d^2=q_d-1` produces the same functional on
coordinate restrictions. Thus the high-level part is at most
`p^(-13/24) kappa 3 sqrt(2) (2/3)^5 M < 3M/4`.
Closing the fixed-dimensional contraction gives `O(m^22)`;
removing `r^-5` gives `O(m^27)`. The argument never assumes a priori
dimension-uniform finiteness and never treats a homogeneous projection as
an infinity contraction. The stated level-flat influence corollary follows
from `S_r(f)^2=binom(n,r)^(1/r) v_r >= (n/r)v_r` and Cauchy--Schwarz.

The reconstruction found no gap in these steps. The exponent `27` is
not a sharp degree-two constant.

## 2. Exact coefficient object and cancellation

Write the reduced Walsh power as

```math
H_A(x)^k=\sum_{S\subseteq[n]}c_{k,S}(A)x^S.
```

Reduction uses `x_i^2=1`, so multiplication is symmetric-difference
convolution, not ordinary homogeneous polynomial multiplication:

```math
c_{k,S}(A)=\sum_{e_1\triangle\cdots\triangle e_k=S}
                    a_{e_1}\cdots a_{e_k},
\qquad
\sum_S c_{k,S}(A)^2=\mathbb E H_A^{2k}.
```

Only even boundaries with `|S|<=min(2k,n)` occur. At the top level,

```math
c_{k,S}=k!\operatorname{haf}(A[S]),\qquad |S|=2k.
```

Each hafnian is an odd integer, but the parity floor alone yields only
order `n`, rather than `n^(3/2)`, for linear powers. The squared boundary
coefficients are nonnegative **after** all signed contributions have been
collected; the summands inside a boundary coefficient need not have one
sign. This is the same exact cancellation obligation recorded in the
older `high_order_flat_chaos.md` note, not a new positivity mechanism.

For independent computational checks, the first three powers satisfy

```math
c_{2,\varnothing}=N,\quad
c_{2,\{i,j\}}=2(A^2)_{ij},\quad
c_{2,S}=2\operatorname{haf}(A[S])\ (|S|=4),
\qquad N=\binom n2;
```

```math
c_{3,\varnothing}=\operatorname{tr}(A^3),\qquad
c_{3,\{i,j\}}=6(A^3)_{ij}+(3N-12n+16)a_{ij},
```

and for a four-set `S`,

```math
c_{3,S}=6\left[
\sum_{\{i,j\}\subset S}a_{S\setminus\{i,j\}}(A^2)_{ij}
-2\sum_{v\in S}\prod_{w\in S\setminus\{v\}}a_{vw}
\right].
```

Here `a_T` for a two-set means its edge sign. The last identity counts
length-two paths plus a disjoint edge and three-edge stars; the correlation
sum counts each internal star three times, explaining the subtraction.

## 3. The exact boundary-entropy gain

For the direct degree-`2k` BH application, put

```math
p_k=\frac{4k}{2k+1},\qquad
\mathcal F_{n,k}(A)=\left(\sum_S|c_{k,S}(A)|^{p_k}\right)^{1/(kp_k)},
\qquad L_{n,k}(A)=\|H_A\|_{2k}.
```

Let `D_(n,k)` be the number of even subsets through size `2k`. Ordinary
finite-dimensional norm comparison and Parseval give the exact sandwich

```math
L_{n,k}(A)\le \mathcal F_{n,k}(A)
\le D_{n,k}^{1/(4k^2)}L_{n,k}(A).                 \tag{3.1}
```

There is also an exact entropy identity. Set
`pi_S=c_(k,S)^2/sum_T c_(k,T)^2` and `alpha=p_k/2`.
For Renyi entropy `H_alpha(pi)=log(sum pi^alpha)/(1-alpha)`,

```math
\log\frac{\mathcal F_{n,k}(A)}{L_{n,k}(A)}
=\frac{H_\alpha(\pi)}{4k^2}.                    \tag{3.2}
```

Thus the sole gain over the ordinary moment is boundary entropy, bounded
by `log D_(n,k)`. If `2k<=n/2`, then
`D_(n,k)<=sum_(j<=2k) binom(n,j)<=(en/(2k))^(2k)`. In particular,

```math
\mathcal F_{n,k}(A)=(1+o(1))L_{n,k}(A)
\quad\hbox{uniformly in }A\quad\hbox{if }k/\log n\longrightarrow\infty.
                                                               \tag{3.3}
```

For `k>=n/4`, the bound `D_(n,k)<=2^(n-1)` gives the same conclusion.
At `k=alpha n`, the ratio is `1+O_alpha(1/n)`. Using the valid degree bound
`m=min(2k,2 floor(n/2))` and exponent `2m/(m+1)` only replaces the entropy
denominator in (3.2) by `2mk`; it does not change these conclusions.
This `m` need not be the actual reduced degree: additional cancellations
can lower it. Choosing any smaller positive valid degree bound works too.

The weighted theorem does not avoid the obstruction. Since

```math
\|c^{=r}\|_{2r/(r+1)}
\le \binom nr^{1/(2r)}\|c^{=r}\|_2
\le \sqrt{en/r}\,\|c^{=r}\|_2,
```

we have, for every signing and every power,

```math
\mathcal W(H_A^k)^2
\le en\sum_{r\ge1}r^{-11}\|c_k^{=r}\|_2^2
\le en\,\mathbb E H_A^{2k},
\qquad
\mathcal W(H_A^k)^{1/k}\le(en)^{1/(2k)}L_{n,k}(A). \tag{3.4}
```

Consequently any leading-order lower bound produced by the weighted
functional at `k>>log n` is already a leading-order lower bound on the
ordinary all-boundary moment, up to a factor tending to one.

## 4. Stronger no-go: every sublinear power is asymptotically uninformative

**Theorem (derived here).** For every integer sequence
`2<=k_n=o(n)`,

```math
\frac{1}{n^{3/2}}\min_A\mathcal F_{n,k_n}(A)\longrightarrow0,
\qquad
\frac{1}{n^{3/2}}\min_A\mathcal W(H_A^{k_n})^{1/k_n}
\longrightarrow0.                                           \tag{4.1}
```

The conclusion also holds for the global coefficient norm with any
positive valid bound on the reduced degree. A single explicit all-order signing family
witnesses both limits. Thus applying either of the paper's inequalities
to a sublinear power `k>=2` cannot give even a positive universal
`n^(3/2)`-scale lower constant. This is not a restriction on all possible
uses of the paper, only on these two reduced-power applications.
The witnesses need not be actual minimum-cap matrices, and no assertion
is made that actual minimizers have their moment law. Existence of these
admissible controls is enough to refute a proposed **universal** coefficient
lower bound at that scale.

In fact, the estimates below prove the stronger adaptive-power statement:
for every `K_n=o(n)` with `K_n>=2`,

```math
\frac{1}{n^{3/2}}\min_A\max_{2\le k\le K_n}
 \max\{\mathcal F_{n,k}(A),\mathcal W(H_A^k)^{1/k}\}
\longrightarrow0.                                           \tag{4.1a}
```

Thus choosing the sublinear power after seeing the signing does not
remove this obstruction.

### 4.1 A shorter minimax proof by averaging coefficient signs

For independent uniform edge signs and each fixed spin configuration,
`H_A(x)` has the law of a sum of `N=binom(n,2)` independent signs.
Therefore

```math
\mathbb E_A\mathbb E_x H_A(x)^{2k}
=\mathbb E\left(\sum_{e=1}^N\varepsilon_e\right)^{2k}
\le(2k-1)!!\,N^k\le(2kN)^k.
```

The first inequality follows coefficientwise from
`(2r)!>=2^r r!` in the product of the even power series of cosh; odd
moments vanish. Equivalently, a sum of independent signs has its even
moments bounded by those of the Gaussian with the same variance.
Some admissible signing consequently satisfies
`L_(n,k)<=sqrt(2kN)`. Combining this with the coefficient norm comparisons
gives the sharper minimax bounds

```math
\min_A\mathcal F_{n,k}(A),\quad
\min_A\mathcal W(H_A^k)^{1/k}
\ \le\ (en)^{1/(2k)}\sqrt{2kN}.                  \tag{4.2a}
```

The same one-variable maximization used below yields, uniformly for
`2<=k<=K<=n`,

```math
\frac{\min_A\mathcal F_{n,k}(A)}{n^{3/2}},\quad
\frac{\min_A\mathcal W(H_A^k)^{1/k}}{n^{3/2}}
\le C\left(n^{-1/4}+\sqrt{K/n}\right).            \tag{4.2b}
```

This proves (4.1) immediately. The averaged signing may depend on `k`;
this is allowed under the minimization and does not describe the moment
law of any actual minimum-cap family. The explicit construction below
is stronger in a different direction: one predetermined all-order family
works simultaneously for every sublinear power.

One can also make the averaging argument simultaneous, proving (4.1a)
without spectral input. Put

```math
Z_A=\sum_{k=2}^{\infty}2^{-k}
 \frac{\mathbb E_x H_A(x)^{2k}}{(2k-1)!!\,N^k}.
```

Tonelli and the same moment comparison give `E_A Z_A<=1/2`.
There is therefore a signing with `Z_A<=1/2`, and for that single signing
every moment satisfies

```math
\mathbb E_x H_A(x)^{2k}\le2^k(2k-1)!!N^k,
\qquad \|H_A\|_{2k}\le2\sqrt{kN}\quad(k\ge2).
```

The harmless weakened factor `2^k` avoids tracking the extra half.
The uniform maximization leading to (4.2b) now proves (4.1a), with
`C(n^-1/4+sqrt(K/n))`. This existence proof chooses one signing for each
`n` that works at **all** powers simultaneously, but does not provide the
explicit matrix supplied by the Sylvester construction.

The averaging proof works unchanged for **any** fixed Walsh support
`Lambda` of `N` characters of degrees at most `d`, with independent
coefficient signs. The standard coefficient exponent is then
`q_(dk)`, and the bound `#supp(f^k)<=(n+1)^(dk)` gives the same
`(en)^(1/(2k)) sqrt(2kN)` estimate. Fixed coefficient magnitudes are also
allowed, replacing `N` by their squared sum. Thus the obstruction is not
specific to complete quadratic support.

### 4.2 An elementary moment bound for a fixed signing

For any real hollow symmetric `A`, let
`N_A=sum_(i<j) a_ij^2` and `L=||A||_op`. Then, with an absolute constant,

```math
\|H_A\|_{2k}\le C\left(\sqrt{kN_A}+kL\right).       \tag{4.2}
```

Here is a direct proof, so (4.2) does not depend on importing a
Hanson--Wright theorem. Split coordinates randomly into `I,J` and put
`H_cross=x_I^T B x_J`. Each edge is retained with probability one half,
so Jensen gives

```math
\mathbb E_x e^{tH_A(x)}
\le\mathbb E_{I,J}\mathbb E_{x_I,x_J}e^{2t x_I^T Bx_J}.
```

Condition on `x_J` and use `cosh u<=exp(u^2/2)`. Linearize the resulting
square with a standard Gaussian `g`, and use the same cosh bound for
`x_J`. This gives

```math
\mathbb E e^{2t x_I^T Bx_J}
\le\mathbb E_{x_J}e^{2t^2\|Bx_J\|_2^2}
=\mathbb E_g\mathbb E_{x_J}e^{2t g^T Bx_J}
\le\mathbb E_g e^{2t^2\|B^Tg\|_2^2}
=\det(I-4t^2BB^T)^{-1/2}.
```

If `|t|<=1/(sqrt(8)L)`, every eigenvalue `z` of `4t^2BB^T` is at most
one half. The estimate `-log(1-z)<=2z` then gives

```math
\log\mathbb E e^{tH_A}\le4t^2\|B\|_F^2\le4t^2 N_A.
```

The bound is uniform in the partition, hence remains true after averaging.
It also holds for negative `t`. Chernoff optimization yields
`Pr(|H_A|>=u)<=2 exp[-c min(u^2/N_A,u/L)]` for an absolute `c>0`.
Integrating this tail proves (4.2).

### 4.3 Explicit all-order witnesses

Let `h` be the least power of two at least `n`, and let `S_h` be the
symmetric Sylvester matrix, recursively `S_(2h)=[[S_h,S_h],[S_h,-S_h]]`.
It has sign entries and `S_h^2=hI`. Take its leading `n` by `n` principal
submatrix and set the diagonal to zero, obtaining a full signing `A_n`.
Since principal compression is an operator-norm contraction,

```math
\|A_n\|_{\rm op}\le\sqrt h+1\le\sqrt{2n}+1.
```

Equation (4.2) therefore gives

```math
L_{n,k}(A_n)\le C(n\sqrt k+\sqrt n\,k).           \tag{4.3}
```

Both coefficient functionals in (4.1) are bounded above by
`(en)^(1/(2k))` times (4.3), up to an immaterial absolute constant.
For the global functional, use the elementary support bound
`D_(n,k)<=(n+1)^(2k)`; for actual reduced degree `m`, use
`D<=(n+1)^m`. The weighted case is (3.4).

Uniformly for `2<=k<=K<=n`, elementary one-variable maximization gives

```math
\frac{\max\{\mathcal F_{n,k}(A_n),
                 \mathcal W(H_{A_n}^k)^{1/k}\}}{n^{3/2}}
\le C'\left(n^{-1/4}+\sqrt{K/n}+K/n\right).       \tag{4.4}
```

Indeed `sqrt(k) exp(log(en)/(2k))` decreases and then increases, with
its unique minimum at `k=log(en)`; its maximum over an interval occurs
at an endpoint. The analogous function `k exp(log(en)/(2k))` has its
minimum at `k=log(en)/2`. Separate intervals below and above those
minima to obtain (4.4). Taking `K=o(n)` proves (4.1).

### 4.4 Stronger actual-class statement: all bounded-cap signings through o(n^(3/4))

Status: **independent proof-verifier audit PASS**, including all endpoint
and actual-degree cases.

The preceding existential controls do not identify the moment behavior
of actual minimizers. A separate cap-to-spectrum argument now gives a
uniform upper bound for **every** bounded-cap signing, and hence for
every exact minimizer.

**Theorem.** Fix `C>0` and `s>=0`. There is `C_C<infinity` such that
for all `n>=2`, all integers `K>=2`, and every full signing with
`Q(A)<=C n^(3/2)`,

```math
\frac1{n^{3/2}}\max_{2\le k\le K}
 \max\{\mathcal F_{n,k}(A),\mathcal W_s(H_A^k)^{1/k}\}
\le C_C\left(n^{-1/4}+\sqrt{K/n}+K/n^{3/4}\right).          \tag{4.5}
```

The constant can be chosen uniformly over `s>=0`. In particular, if
`K_n=o(n^(3/4))`, all these adaptive power certificates are
`o(n^(3/2))`, uniformly over the entire bounded-cap class. This includes
**all choices of exact minimizers**: the explicit Sylvester controls
already imply `M_n<=C_0 n^(3/2)` with an absolute `C_0`. No
bounded-operator assumption or conjecture on minimizer distributions is
being imposed. This is an upper moment/certificate theorem, not a
determination of their moment law.

For completeness, the independently audited cap-to-spectrum estimate in
`artifacts/bh_2026_09_16_nonlinear_filter_audit.md`, Section 9, is

```math
\|A\|_{\rm op}\le\sqrt{8Q(A)}\le\sqrt{8C}\,n^{3/4}.
```

Its proof uses real polarization to give the bilinear cube norm
`B(A)<=4Q(A)`, the safe complex endpoint bound `2B(A)`, and three-lines
interpolation against `max|a_ij|<=1`. Combining this with the elementary
moment estimate (4.2), rather than any distributional assumption, yields

```math
\frac{\|H_A\|_{2k}}{n^{3/2}}
\le C_C\left(\sqrt{k/n}+k/n^{3/4}\right)\qquad(k\ge1).      \tag{4.6}
```

The support bound and the weighted comparison in Section 3 give, for
both quantities on the left of (4.5), the pointwise upper bound

```math
(en)^{1/(2k)}\|H_A\|_{2k}.
```

Here is an explicit endpoint lemma that removes any hidden dependence
on `K`. For `alpha>0`, `b>1`, and `K>=2`,

```math
\max_{2\le k\le K}k^\alpha b^{1/(2k)}
\le\max\{2^\alpha b^{1/4},e^\alpha K^\alpha\}.             \tag{4.7}
```

Indeed, the logarithmic derivative has the sign of
`alpha k-(log b)/2`, so the function decreases and then increases.
If `K` is below the turning point, the maximum is at two. Otherwise the
other possible endpoint obeys `b^(1/(2K))<=e^alpha`. Apply (4.7) with
`b=en` and `alpha=1/2,1` to (4.6). The two small-endpoint terms are
`O(n^(-1/4))` and `O(n^(-1/2))`; the latter is absorbed by the former.
This proves (4.5), uniformly in every integer `K>=2`.

**Reduced-degree caveat.** The displayed `F_(n,k)` uses the valid upper
degree bound `2k`, not a claim that reduction preserves that degree.
The same estimate also holds if one uses a smaller valid positive
degree bound `m<=2k`, including the actual nonzero degree. In fact,

```math
\|\widehat{H_A^k}\|_{q_m}^{1/k}
\le D^{1/(2mk)}\|H_A\|_{2k}
\le(n+1)^{1/(2k)}\|H_A\|_{2k},
```

because `D<=sum_(r<=m)binom(n,r)<=(n+1)^m`. A constant reduced power
can simply be assigned the valid degree bound one. Thus cancellations
that lower the actual degree do not create an escape from (4.5).

The exponent `3/4` is sharp for this **whole bounded-cap class**, though
not asserted sharp for minimizers. The planted-clique family in
`artifacts/bh_2026_09_16_planted_tail_sharpness.md` satisfies
`Q(A)<n^(3/2)` and, along `n=16^t`,
`Pr(H_A/n^(3/2)>=a)>=c_a2^(-n^(3/4))` for each fixed `0<a<1/2`.
For `k_n=floor(tau n^(3/4))`, `tau>0`, this implies

```math
\liminf_{n=16^t\to\infty}\frac{\mathcal F_{n,k_n}(A)}{n^{3/2}}
\ge\liminf_{n=16^t\to\infty}\frac{\|H_A\|_{2k_n}}{n^{3/2}}
\ge a\,2^{-1/(2\tau)}>0.
```

Thus the little-o conclusion cannot hold uniformly over bounded-cap
signings at power order `n^(3/4)`. This does not improve the original
minimum or establish convergence of its normalized values.

## 5. Comparison with the old subexponential constants

For the standard power application, BH says

```math
Q(A)\ge B_{2k}^{-1/k}\mathcal F_{n,k}(A).
```

The new factor is `exp[-O(log k/k)]`. The old bound from
Defant--Mastylo--Perez gives `exp[-O(sqrt(log k/k))]`.
**Both tend to one whenever `k` tends to infinity.** A smaller vanishing
factor is a finite-size improvement, not a new asymptotic constant.
For bounded `k>=2`, (4.1) shows that the entire coefficient lower-bound
route is lower order. For `k=1`, the problem is the sharp degree-two
constant, which `K m^27` does not determine or improve numerically.

At linear powers, (3.1) makes the coefficient functional an asymptotically
identical moment. At superlinear powers, the finite-cube squeeze

```math
\|H_A\|_{2k}\le Q(A)
\le2^{(n-1)/(2k)}\|H_A\|_{2k}
```

makes the moment itself asymptotically identical to the original supremum.
At `k=alpha n`, this latter factor remains `2^(1/(2alpha))`, even though
the BH constant and boundary-entropy factor both tend to one.

In fact, at linear powers neither degree-uniform BH estimate is needed:
the dimension-dependent elementary estimate already gives

```math
\|c_k\|_{p_k}\le 2^{(n-1)/(4k)}\|c_k\|_2
\le2^{(n-1)/(4k)}Q(A)^k.
```

The prefactor is `O_alpha(1)` at `k=alpha n`, smaller than the polynomial
degree bound before taking the root, and its `k`th root is `1+O_alpha(1/n)`.
For the weighted theorem, (3.4) similarly gives the elementary bound
`W(H_A^k)<=sqrt(en) Q(A)^k`. Thus the informative regime already has a
stronger finite-dimensional bound than either old or new BH growth bound.

The same pointwise sandwiches pass through `min_A`. Hence convergence of
the global power-coefficient minima in the superlinear regime is equivalent
to the original convergence problem, not a new composition theorem. In the
linear regime the missing lower estimate is exactly the signed
all-boundary moment estimate from the older artifact.

The weighted functional also becomes equivalent at superlinear powers,
but this requires a separate lower bound, not the one-sided estimate
(3.4). Fixed degree-two BH gives `Q(A)>=c n^(3/2)` uniformly, and hence
`tau_n=sqrt(2N)/Q(A)=O(n^-1/2)`. Markov's inequality gives probability at
least one half to `|H_A|<=sqrt(2N)`. An absolute maximizer has projective
probability at least `2^(-(n-1))`. Using independent spins `X,X'`,

```math
\operatorname{Var}(H_A^k)
=\frac12\mathbb E\big(H_A(X)^k-H_A(X')^k\big)^2
\ge2^{-n}\left(Q(A)^k-(2N)^{k/2}\right)^2.
```

This holds for all sufficiently large `n`, when `Q(A)>sqrt(2N)`; the
absolute-difference lower bound works for both odd and even `k`.
Since `W(H_A^k)>=n^-5 sqrt(Var(H_A^k))`, we obtain

```math
n^{-5/k}2^{-n/(2k)}Q(A)(1-\tau_n^k)^{1/k}
\le\mathcal W(H_A^k)^{1/k}
\le(en)^{1/(2k)}Q(A).
```

For `k/n->infinity` both ratios tend to one uniformly in `A`; thus the
weighted minima have the same superlinear convergence equivalence too.

## 6. Reproducible finite checks

Run

```sh
.venv/bin/python computations/bh_2026_09_16_power_checks.py \
  --output computations/results/bh_2026_09_16_power_checks.json
```

The script checks convolution support, Parseval, all top-layer hafnians for
`k=1,2,3`, and the displayed lower-level second/third-power formulas using
exact integer arithmetic on 36 random/Sylvester matrices. It additionally
loads the repository's actual stored minimum-cap witnesses of orders
`6,...,14`, rechecks their caps `(5,9,10,12,13,17,18,20,21)` by exhaustive
spin evaluation, and runs the same identities and coefficient-norm checks.
The sources and matrix hashes are recorded in the JSON. The witnesses'
global optimality is inherited from the repository, not reproved here.

For the explicit Paley symmetric conference signings at orders `6` and
`14`, the script verifies `A^2=(n-1)I` in exact arithmetic. At order `6`,
every four-vertex hafnian has absolute value `1`, although it is the sum
of three signs. At order `14`, the sum of squared six-vertex hafnians is
`27027`, compared with the random-sign expectation `45045`; on eight
vertices the corresponding totals are `125307` and `315315`. These are
finite demonstrations of signed cancellation, not asymptotic hafnian laws.

The same controls show cancellation in the complete Eulerian expansion.
At order `14` and `k=14`, the coefficient-independent family in which `k`
distinct edges each occur twice contributes approximately `868950.48`
times the **entire** moment `E H^(2k)`. The script stores both integers
exactly, along with the negative difference made up by all other diagrams.
This replays, rather than replaces, the older scalable conference
cancellation argument.

All tests passed. The reconstructed primary contraction factor is
`0.6742048108733384`, safely below `3/4`. Floating-point norm checks are
labeled numerical; the asymptotic conclusions rest on the proofs above.
The script also exhaustively averages all coefficient signings at orders
four and five, verifying the independent-sign even-moment identity and
the simultaneous weighted-moment selection argument with exact fractions.

## 7. Research status

No improved universal numerical bound, recurrence, or convergence result
has been obtained. The genuinely sharper route diagnosis is (4.1):
sublinear powers are excluded altogether, while informative powers lose
their coefficient advantage and recover the old all-boundary moment
problem. The new polynomial growth theorem is valid and quantitatively
strong, but its reduced-power applications do not remove that obligation.
Section 4.4 additionally excludes adaptive `o(n^(3/4))` powers for
every actual exact minimizer, not only for the existential controls.

This note does not rule out other uses of the weighted contraction,
coefficient flattening with additional structure, or a new theorem on the
signed convolution orbit. Those would need separate input and auditing.

## 8. Related completed refinements

The director's nonlinear-filter analysis is preserved in
`artifacts/bh_2026_09_16_nonlinear_filter_audit.md`. Its coefficient
Cauchy bound, weighted exponential sum, full coefficient-norm estimate,
regular-Hadamard false-cap controls, and rare-tail extension have received
independent power-agent checks. They are broader filter-specific negative
statements, not claims about exact minimizer moment laws.

The discriminating question of a hypothetical **linear** weighted BH
growth bound is treated separately in
`artifacts/bh_2026_09_16_linear_weighted_filter_analysis.md`. Its uniform
all-sign expansion shows that `o(sqrt(n))`-degree filters have exactly
the same normalized leading test as the quadratic form itself. An
all-degree linear constant is constrained by the order-six witness;
an asymptotic-only even-class slope has different quantifiers and is not
ruled out by that fixed finite witness.

Finally, high-degree exponential filters have an elementary
thermodynamic comparison that does not need polynomial BH. For any
Walsh function `f` with positive degree bound `m`,

```math
\|f\|_2\le\|\widehat f\|_{2m/(m+1)}
\le2^{n/(2m)}\|f\|_2.
```

Thus normalized logarithms differ by at most `log(2)/(2m)`. Taking the
valid degree bound `m=n` for `exp(beta H/sqrt(n))` or its cosh analogue,
the full coefficient norm has exactly the ordinary `L_2`/free-energy
rate. The dimension-dependent prefactor is only `sqrt(2)`, already
stronger than the degree-polynomial bound here. This does not alter the
minimization-over-disorder quantifier of the finite-temperature problem.

## 9. Appendix: coefficient flattening by permutation lifts

This bounded follow-up isolates exact scalar and Hilbert norm costs.
The entire appendix, including both normalization identities and the
finite replay, received an independent proof-verifier audit: **PASS**.
These results do not rule out transformations using new
structure outside these hypotheses.

### 9.1 Scalar signed averaging exposes the missing dispersion estimate

Write `c_S=fhat(S)`, `D_r=binom(n,r)`, and
`V_r=sum_(|S|=r)|c_S|^2`. For

```math
Tf=\sum_{\pi\in S_n}w_\pi f\circ\pi,\qquad
\sum_\pi|w_\pi|\le1,
```

both the supremum norm and every level's coefficient absolute sum are
contracted. If the output has constant coefficient magnitude on level
`r`, then necessarily

```math
V_r(Tf)\le\frac{\|c^{=r}\|_1^2}{D_r}
=\eta_rV_r(f),\qquad
\eta_r=\frac{\|c^{=r}\|_1^2}{D_r\|c^{=r}\|_2^2}\le1.       \tag{9.1}
```

For a nonzero level, equality `eta_r=1` holds exactly when the input
magnitudes were already flat. Thus lossless scalar flattening by such
averages is impossible for a genuinely nonflat level. Near-lossless
flattening would require a new near-flatness/dispersion estimate; the
averaging operation does not supply it. This is a necessary bound, not
an assertion that every allowed output or bound is attainable.

The ordinary permutation average has coefficient
`D_r^(-1) sum_(|S|=r)c_S` on that level. In particular it annihilates
`H_A` if the total edge-sign sum is zero. More generally, if a signed
orbit average happens to be permutation-invariant, it equals
`(sum_pi w_pi)` times the ordinary average. Hence signed averaging cannot
produce a nonzero symmetric output from a zero-average orbit.

Noise multiplies the entire level by `rho^r`, so for `rho>0` it leaves
`eta_r` and all within-level magnitude ratios unchanged. Averaging
independent random coordinate resampling gives this same noise operator.
An individual restriction can change coefficients, but does not
automatically preserve full-level flatness on the smaller cube.

### 9.2 The Hilbert orbit does flatten, but scalar BH fails on this very class

For the normalized Hilbert orbit

```math
G(x)=\big(f(\pi x)/\sqrt{n!}\big)_{\pi\in S_n},
```

permutation transitivity gives

```math
\|G\|_{\infty;H}\le\|f\|_\infty,\qquad
\|\widehat G(S)\|_H^2=V_r/D_r,\qquad
\operatorname{Inf}_i(G)=\frac1n\sum_{r\ge1}rV_r.             \tag{9.2}
```

Thus coefficient-vector lengths really are level-flat. Nevertheless,
take `f=x_1...x_m`, with `m` fixed and `D=binom(n,m)`. This orbit has
pointwise Hilbert norm one, variance one, and every influence `m/n`.
Its level-`m` coefficient vectors are orthogonal, each of norm
`D^(-1/2)`, and hence

```math
\left(\sum_{|S|=m}\|\widehat G(S)\|_H^{q_m}\right)^{1/q_m}
=D^{1/(2m)}\longrightarrow\infty.
```

Consequently neither scalar BH nor the paper's level-flat influence
corollary extends dimension-freely to Hilbert-valued functions, **even
when restricted to permutation orbit lifts**. The scalar-field step,
not the orbit identities, is the failure.

For this same monomial orbit, any Hilbert functional `u` satisfies the
sharp scalar-projection bound

```math
\operatorname{Var}\langle u,G\rangle\le\|u\|_H^2/D.        \tag{9.3}
```

Restoring scalar variance one therefore requires `||u||>=sqrt(D)`.
More generally, any scalar degree-`m` homogeneous polynomial with all
coefficient magnitudes `D^(-1/2)` necessarily has
`||g||infinity>=D^(1/(2m))/B_m>=sqrt(n/m)/B_m`. Thus at **fixed degree**,
preserving both full flatness and order-one variance imposes an
unavoidable order-`sqrt(n)` scalar supremum cost, regardless of the
construction. This is a norm obstruction, not a new quadratic bound;
the statement does not claim the same loss uniformly in growing degree.

The elementary auxiliary-variable scalarization
`F(x,y)=D^(-1/2) sum_(|S|=m)y_S x^S` has exactly
`||F||infinity=sqrt(D)`, since the independent `y_S` can align all terms.
It also is not full-level-flat on the enlarged `(n+D)`-cube: only `D`
of the degree-`m+1` supports occur. More elaborate encodings would need
their own degree, supremum, and complete-support analysis.

### 9.3 Exact scalar-projection cost for full row-balanced quadratic signings

There is a target-class version, not just a sparse-monomial example.
Suppose a full hollow symmetric sign matrix has every row sum zero,
and put `D=binom(n,2)`. Such signings exist for every `n=4t+1`: on the
cyclic vertex set, assign `+1` at distances `1,...,t` and `-1` at
distances `t+1,...,2t`. For the orbit of its `H_A`, the Gram matrix of
the edge coefficient vectors has entries

```math
K_{e,f}=\begin{cases}
1,&e=f,\\
-1/(n-2),&|e\cap f|=1,\\
2/((n-2)(n-3)),&e\cap f=\varnothing.
\end{cases}
```

The adjacent entry follows by squaring each zero row sum; the disjoint
entry follows from the zero total edge sum. If `B` is the unsigned
vertex-edge incidence matrix, direct summation then gives

```math
K=\frac{n-1}{n-3}P_{\ker B},\qquad
\dim\ker B=D-n.
```

Therefore the exact minimum Hilbert functional norm needed to produce
scalar coefficient vector `b` in `ker B` is

```math
\min\{\|u\|_H:\widehat{\langle u,G\rangle}=b\}
=\sqrt{\frac{n-3}{n-1}}\,\|b\|_2.                         \tag{9.4}
```

In particular a norm-one functional retains variance at most
`(n-1)/(n-3)`, although the orbit variance is `D`. Producing a full
flat sign vector requires norm `sqrt(D-n)`, and this is attained for
the original coefficient vector. Crucially, this is a cost in the
**generic Hilbert-Cauchy cap guarantee** `||<u,G>||infinity<=||u||Q(A)`.
It is not a lower bound on the actual scalar cap: recovering the
original polynomial plainly recovers its actual cap `Q(A)`.

In fact this same full-sign orbit gives a stronger vector-AA falsifier.
Its **actual Hilbert supremum**, denoted `R`, is exactly

```math
R^2=\frac{(n-1)^2(n+1)}{2(n-2)}.                           \tag{9.5}
```

To verify this, put `s=sum_i x_i`, `v_(ij)=x_i x_j`. Since
`BB^*=(n-2)I+J` and `Bv=sx-1`,

```math
\|P_{\ker B}v\|_2^2
=D-s^2-\frac n{n-2}
 +\frac{(s^2-n)^2}{2(n-2)(n-1)}.
```

This is convex in `s^2` on `[1,n^2]`; the endpoint `n^2` gives zero,
and the endpoint `1` gives its maximum. Multiplication by
`(n-1)/(n-3)` proves (9.5). Thus `G/R` has Hilbert supremum one and

```math
\operatorname{Var}(G/R)=\frac{n(n-2)}{n^2-1}\longrightarrow1,
\qquad
\operatorname{Inf}_i(G/R)=\frac2n\operatorname{Var}(G/R)
\longrightarrow0.
```

The scalar level-flat influence conclusion therefore fails already on
**permutation orbits of full row-balanced quadratic signings**. Also,
`R=Theta(n)` while the scalar cap is at least order `n^(3/2)`: the
Hilbert lift has lost a factor of at least order `sqrt(n)` in its actual
supremum. The weaker upper bound `R<=Q(A)` conceals this distortion.

### 9.4 The power hypothesis must be checked, including the order-six exception

The primary Corollary 1.3 requires equal magnitudes on **every** support
of each level, not just on nonzero coefficients. Powers do not generally
meet it. For the stored cap-nine order-seven witness, exact Walsh
coefficients of `H_7^2` have absolute histograms

```text
level 2: 2 occurs 18 times; 6 occurs 3 times;
level 4: 2 occurs 30 times; 6 occurs 5 times.
```

Both levels have `eta_r=27/35`. Thus (9.1) forces every level-flat
signed permutation average to lose at least `8/35` of the nonconstant
variance of this square. No positive noise parameter fixes flatness.

The stored order-six witness is a genuine exception, not a
counterexample: writing `X=x_1...x_6`, exact coefficients give
`H_6^2=15+2XH_6`. Since `X^2=1`, every power lies in the span of
`1,H_6,XH_6,X`, which occupy separate flat levels `0,2,4,6`. Therefore
all of its powers are level-flat. This exceptional identity cannot be
assumed for arbitrary signings or dimensions.

Finally, normalizing a power orbit by its **scalar** cap is different
from normalization by its actual Hilbert cap. For the orbit of
`f=(H_A/Q(A))^k`, write `V=Var(G)>0`, `I=Inf_i(G)`, and let `m` be
its actual degree. Evenness and `|H_A/Q(A)|<=1` give

```math
I\ge2V/n,\qquad V\le\mathbb E|H_A/Q(A)|^{2k}\le N/Q(A)^2,
\qquad \frac I{V^2}\ge\frac{2Q(A)^2}{nN}.
```

Equality holds at `k=1`, where `m=2`. Consequently, for every `a>=0`,

```math
\inf_{A,\ k\ge1:\,V>0}m^a\frac I{V^2}
=\frac{2^{a+1}M_n^2}{nN}.                                 \tag{9.6}
```

So a restricted influence theorem for these `Q(A)`-normalized power
orbits can be true, but its best degree-weighted constant exactly
repackages the original quadratic minimum and is already attained at
the first power. This avoids overextending the vector counterexample
to a different normalization and hypothesis class.

The exact checks are reproducible with
`.venv/bin/python computations/bh_2026_09_16_orbit_flattening_checks.py`.
No improved original-value bound or convergence statement is claimed.
