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
