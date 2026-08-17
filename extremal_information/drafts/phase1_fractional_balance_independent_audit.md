# Independent audit: multiselector packing and fractional balance

**Scope.**  I examined only
`multiselector_sparse_flip_rate.md`,
`near_minimizer_fractional_balance.md`, and the cited canonical statements
Theorem 21.62--Theorem 21.67.  I did not use the other linked PC.3 drafts or
verifier.  This audit treats signed cuts as vectors, rather than counting
their redundant `(sigma,x)` representations.

**Overall verdict.**  The four requested central claims pass:

| Claim | Verdict | Qualification |
|---|---|---|
| Theorem MS.3 | **PASS** | It is a polynomial-in-`N_j` one-hot state packing, not a `k_j`-bit cube. |
| Corollary MS.4 | **PASS** | It is a contextual packing with one externally selected compiler block per query, not one parent containing all query shores. |
| Theorem FB.1 | **PASS** | The minimax identity, perturbation normalization, union bound, and shell constants are correct. |
| Corollary FB.2 | **PASS** | The empirical first-marginal estimate is correct with no union bound over edges. |
| Theorem FB.3 | **PASS with statement cleanup** | The Serfling exponent and exact-flip repair are correct; `c_0`, integrality, and asymptotic quantifiers should be made explicit. |

There are two false or overbroad subsidiary sentences, neither of which
invalidates a theorem: the remaining-field mean bound in Lemma MS.1 has the
wrong additive constant, and FB.3 cardinality does not by itself rule out a
bounded *generative description* of the shell.

## 1. Common normalizations

For a hollow signing `a` and its symmetric matrix `A`,

```math
 {\sigma\over2}x^TAx=\sum_{i<j}a_{ij}(\sigma x_ix_j)
 =\langle a,z\rangle.
```

Thus the conventions in (MS.1), (FB.1), and the cited definition (21.377)
agree.  Flipping one unordered edge changes this normalized energy for a
fixed signed cut by **two**, not four.  Four is the change in the unhalved
quantity `x^TAx`; consequently (MS.37) is valid but non-sharp under the
draft's `Q` normalization.  The sharpened version is

```math
 |F_ell|\ge {c\over2}N_j^{3/2},
 \qquad k\le {1+o(1)\over c}\sqrt {N_j}.
```

The parent metric normalization is also consistent.  Since
`s_j=Theta(sqrt(N_j j))=o(N_j)`,

```math
 {N_j^{3/2}\over(N_j+s_j)^{3/2}}=1+o(1).
```

Accordingly the raw parent gap in (MS.32) is exactly what is needed for
(MS.3), with a uniform `o(1)`.

Hollowing causes no hidden energy shift.  The sparse-flip construction does
not alter the diagonal of the trace-zero base matrix, so `tr(H(u))=0`.
Putting `A(u)=H(u)-diag(H(u))` therefore preserves every Boolean quadratic
energy.  It changes the scaled operator approximation by at most
`1/sqrt(N_j)`.  The matrices in (MS.2) are consequently hollow exact
signings: the child and shore are hollow, and every cross-block entry is
exactly `+-1`.

## 2. Audit of the multiselector draft

### 2.1 Metric and counterexample: PASS

The metric (MS.3) compares one common state family under one common labelled
query bank.  For a pair `u ne v`, selecting label `u` makes the oriented gap
(MS.4) imply the absolute pairwise gap in (MS.3).  The definition does not
claim that the different blocks `B(ell)` coexist in one physical matrix.

The collapse example (MS.8)--(MS.9) is normalized correctly.  Swapping one
`+1` and one `-1` in `alpha` changes the field by at most four, so uniform
anti-concentration gives selector disagreement probability `O(j^(-1/2))`.
Hence

```math
 e^Te'=1-O(j^{-1/2}),
 \qquad
 \|ee^T-e'e'^T\|_{op}
   =\sqrt{1-(e^Te')^2}=O(j^{-1/4}).
```

The response Lipschitz bound is

```math
 |\mathcal B_A(g)-\mathcal B_{A'}(g)|
 \le {N_j\over2}\|A-A'\|_{op},
```

so the difference is `o(N_j^(3/2))`.  The same bound applies after a common
compiler because only the child block changes.

### 2.2 Lemma MS.1 and moment formulas: PASS after one local correction

The moment calculation (MS.12) is exact:

```math
 E A_t={alpha_t\over2},
 \qquad Var(A_t)={7\over4}-alpha_t beta_t.
```

The variance is uniformly between `3/4` and `11/4`, and the local
nonconstant character means have modulus at most `1/2`.  These are the
uniform ingredients needed by the factor-deletion proof of (MS.11).

One sentence in that proof is false as written.  After deleting `k` local
factors, the remaining score has mean

```math
 1+{1\over2}\sum_{t\notin K}alpha_t,
```

so (MS.10) gives the safe bound

```math
 \left|1+{1\over2}\sum_{t\notin K}alpha_t\right|
 \le {k+3\over2},
```

not `(k+1)/2`.  Already at `k=0` and `sum alpha_t=0`, the mean equals one.
This additive correction does not change the geometric factor-deletion
argument or the `O(j^(-1/2))` conclusion.

A copy-ready replacement for the affected sentence is:

> After deleting any `k` local factors, the remaining score has variance
> comparable to `j-k`, while (MS.10) bounds the modulus of its mean by
> `1+(k+1)/2=(k+3)/2`.

The formulas (MS.14)--(MS.15) pass.  Uniform Wasserstein Berry--Esseen gives

```math
 E|L|=sqrt(2/pi)sqrt(V)+O(1),
 \qquad E L^2=V+O(1),
```

because the field mean is bounded.  Since `x=sgn(h)`, their ratio is
`sqrt(2/pi)+O(j^(-1/2))`, uniformly in the balanced endpoint class.

For standalone readability, (MS.23) should define `d_z(A)` explicitly, or
replace it by the normalized deficit from (21.369).  This is a notation
omission, not a normalization failure.

### 2.3 Lemma MS.2: PASS

The covariance kernel (MS.16) follows directly from
`EX=1/2`, `EY=0`, and `EXY=-1/2`.  In (MS.19), splitting at `tj/2` and using
the coefficient ranges stated in the draft gives the advertised even-`j`
bound

```math
 P\{|C|>tj\}\le4e^{-t^2j/25}.
```

For a sample of
`K=exp(t^2j/100)` endpoints, the even-`j` pair union bound is at most

```math
 {K^2\over2}\,4e^{-t^2j/25}
 \le2e^{-t^2j/50}=o(1).
```

With the stated odd-`j` denominator `30`, it is at most
`2e^(-t^2j/75)=o(1)`.  Thus the union bound is genuinely uniform over all
pairs.  Repeated samples are automatically excluded on the good event,
since a repeated endpoint has `C=V>=3j/4>tj` for `t<1/2`.

The correlation conversion also passes.  Both variances are at least
`3j/4`, so `|Corr|<=4t/3`, while the bivariate covariance has least
eigenvalue at least `(3/4-t)j`.  Fixed-dimensional Berry--Esseen is therefore
uniform, and the bounded field means move quadrant thresholds by only
`O(j^(-1/2))`.  This proves (MS.18).

### 2.4 Theorem MS.3: PASS

Choosing the nearest even `m_u` is harmless: its step-size error in
`m_u||h_u||_2/N_j` is `O(sqrt(j/N_j))=o(1)`.  The bounds in (MS.13) make
`m_u=Theta(sqrt(N_j/j))` uniform in `u`.

The cited sparse-flip concentration can be selected separately for each
state.  No probability union bound across `U_j` is needed for existence.
Moreover, all relevant estimates have a common deterministic error bound:
diffuseness is `O(j^(-1/2))`, the active set has common size `4^j`, and the
matrix concentration parameters do not depend on the endpoint signs.
Thus the `o(1)` in (MS.24), (MS.28), and (MS.29) can be taken uniformly over
the polynomial-size state family.

For the cross state, using `x_u` and the positive trust channel gives

```math
 {\mathcal B_{A(v)}(g(u))\over N_j^{3/2}}
 \ge {1\over2}-{kappa\over2}(e_u^Te_v)^2+b rho+o(1),
```

because both `e_u` and the normalized field lie in the positive eigenspace,
and `m_u x_u^Th_u/N_j^(3/2)=b rho+o(1)`.  The target-state spherical bound is
the uniform extension of (21.382)--(21.383), so subtraction gives exactly

```math
 delta_*=delta-{kappa\over2}mu_{1/10}^2.
```

Numerically, for the declared parameters,

```text
delta = 0.0146650563...,
mu_(1/10) = 0.0851361740...,
delta_* = 0.0128530143... .
```

Thus the displayed weaker claims `delta>0.0146` and `delta_*>0.012` are
correct.

The rate conversion is exact, apart from an immaterial integer floor:

```math
 e^{j/10000}=N_j^{1/(10000\log16)}.
```

### 2.5 Corollary MS.4: PASS

The microcanonical compiler hypotheses hold for every query.  The field
coordinates satisfy

```math
 |g_i(u)|\le m_u(2j+1)\le s_j,
```

and `h_i(u)` is odd while `m_u` and `s_j` are even, so
`g_i(u)=s_j (mod 2)`.  A common even
`s_j=Theta(sqrt(N_jj))` can be chosen because the constants in (MS.22) are
uniform.

The affine compiler error has the scale

```math
 N_j\sqrt{s_j}+s_j^{3/2}\sqrt{N_j}
 =O\left(N_j^{5/4}(j^{1/4}+j^{3/4})\right)
 =o(N_j^{3/2}).
```

For a free endpoint `eta`, putting
`a=<eta,eta_*>/s_j` gives `|a|<=1`.  Trust response is even in its field and
is `l_1`-Lipschitz.  Both spherical upper bounds are increasing in `|a|b`,
so the target parent is bounded by (MS.33) uniformly over every free shore
endpoint.  For the cross parent, selecting `eta_*`, the positive witness
`x_u`, and the appropriate trust sign proves (MS.34).

There is no missing shore factor: the upper bound pays `+Q(C_j)` and the
lower bound pays `-Q(C_j)`, hence subtraction pays `2Q(C_j)`.  Since

```math
 2Q(C_j)=O(s_j^2)=O(N_jj)=o(N_j^{3/2}),
```

(MS.32) follows.  Dividing by `(N_j+s_j)^(3/2)` proves the metric claim.

### 2.6 What rate was and was not proved

The construction has `k_j` states and `k_j` externally labelled one-hot
queries.  Its message set therefore carries

```math
 \log_2 k_j=Theta(j)=Theta(\log N_j)
```

bits.  It does **not** carry `k_j` independently writable bits and does not
produce `2^(k_j)` states.  Calling `k_j` an independent-bit capacity would
overstate the result by an exponential.

The native-selector ceiling `k_j<=4^j=sqrt(N_j)` is correct only for the
stated model assigning a distinct native endpoint query to each one-hot
state.  It is not an upper bound for arbitrary fields, arbitrary query
banks, or independent-bit cubes.

The disjoint-layer ceiling (MS.37) is valid, with the factor-two sharpening
noted above.  The rank ceiling (MS.40) also passes in its explicitly
conditional model: on the selector span `T_j=I`, the operator roof implies
`lambda_max(sum e_ell e_ell^T)<=2/kappa+o(1)`; trace and rank then give
`k<= (2/kappa+o(1))3^j`.  Neither ceiling is a general nonlinear-cube
impossibility theorem.

## 3. Audit of fractional balance

### 3.1 Theorem FB.1: PASS

The minimax identity (FB.7) has the correct order and sign.  Writing a pure
minimum as a minimum over `mu in Delta(S)`, finite bilinear minimax gives

```math
 \max_w\min_mu {1\over E}\sum_e w_e E_mu[a_ez_e]
 =\min_mu\max_w {1\over E}\sum_e w_em_e(mu).
```

The coordinatewise maximum over `0<=w_e<=1` is
`sum_e(m_e)_+`, exactly `E V(mu)`.  A maximizing `w` therefore satisfies
`sum_e w_ea_ez_e>=Ev(S)` for every `z in S`.

For the independent perturbation, the fixed-query variance is at most
`4qE=O(kappa n^(3/2))`, and each centred summand is bounded by two.
Bernstein at logarithmic cost `log|Z_n^+|=O(n)` gives

```math
 C\left(\sqrt{kappa}\,n^{5/4}+n\right),
```

uniformly over all signed cuts.  This is precisely (FB.9) after division by
`n^(3/2)`; the union bound does not depend on the chosen maximizing `w`.

The factor in (FB.10) is exact:

```math
 2qE={2kappa\over\sqrt n}{n(n-1)\over2}
 =kappa(1-1/n)n^{3/2}.
```

For a shell member this yields the decrement in (FB.10).  For a non-shell
member,

```math
 E\langle a^F,z\rangle
 <Q(a)-2kappa n^{3/2}+2qE
 =Q(a)-kappa(1+1/n)n^{3/2},
```

which is slightly stronger than (FB.11).  If `v(S)` strictly exceeds the
right side of (FB.5), the realized shell energies are strictly below
`Q(a)-epsilon n^(3/2)`.  The non-shell energies are too because
`epsilon+eta_n(kappa)<kappa/2` implies
`kappa-eta_n(kappa)>epsilon`.  Hence
`Q(a^F)<M_n`, the required strict contradiction.

Take a minimizing measure in (FB.7).  It has `V(mu)=v(S)`, so it is the same
measure promised in (FB.5).  On the shell,

```math
 P-N=E_mu\langle a,z\rangle
 \ge Q(a)-2kappa n^{3/2}>0,
```

giving `N<=P` and the factor two in (FB.6).  No absolute-value factor is
missing.

The two parameter consequences are correct.  For fixed small `epsilon`,
`kappa=sqrt(epsilon)` gives bias `O(sqrt(epsilon))` after `n -> infinity`.
To invoke (FB.6), the text should explicitly say that the known positive
lower bound on `M_n/n^(3/2)` makes
`Q(a)>2sqrt(epsilon)n^(3/2)` for sufficiently small `epsilon`.  For an exact
minimizer and `kappa=n^(-1/6)`,

```math
 eta_n(kappa)=O(n^{-1/3}),
 \qquad {eta_n(kappa)\over kappa}=O(n^{-1/6}),
```

so both the shell width and normalized total bias are `O(n^(-1/6))`.

### 3.2 Corollary FB.2: PASS

For `X_{r,e}=a_e(z_r)_e in {+-1}` with mean `m_e`,

```math
 E\left|{1\over K}\sum_rX_{r,e}-m_e\right|
 \le\sqrt{{Var(X_{1,e})\over K}}\le K^{-1/2}.
```

Triangle inequality, averaging over edges, and then expectation prove
(FB.13).  One only needs existence of a sample with small *edge average*, so
there is correctly no union bound over `E`.  Taking `K>=delta^(-2)` makes
the additional error at most `delta` (for `delta>0`).

### 3.3 Theorem FB.3: PASS with precise repairs

Let the rigorous lower constant be positive and choose
`0<c_0` below it.  Uniformly for a cut in the deficit-`2r` shell,

```math
 {1\over E}\sum_eg_e(z)
 ={\langle a,z\rangle\over E}
 \ge {c\over\sqrt n}
```

for some `c=c(c_0)>0`, because `r=o(n^(3/2))`.  For an `r`-element sample
without replacement from a `+-1` population of mean at least
`c/sqrt(n)`, Hoeffding--Serfling gives

```math
 P\left\{\sum_{e\in F}g_e(z)\le0\right\}
 \le \exp\left(-c_{tail}{r\over n}\right).
```

This exponent is correct.  Choose the theorem constant
`0<c_1<c_tail`.  If the shell size were below `exp(c_1r/n)`, its bad-event
union probability would be at most
`exp(-(c_tail-c_1)r/n)<1`, since `r>=n`.  Hence one `F` has positive sum for
every shell member.

The exact-flip repair is also correct.  Minimality gives a maximizer `z_*`
of `a^F` satisfying

```math
 d_a(z_*)+2\sum_{e\in F}a_e(z_*)_e\le0.
```

Since the sum is at least `-r`, `d_a(z_*)<=2r`, so `z_*` is in the shell;
the same inequality forces its sum on `F` to be nonpositive, contradicting
the selected `F`.

The following statement repairs should be made:

1. Require `0<c_0`, not merely “any constant below” the lower bound.
2. State that `r=r_n` is an integer sequence with
   `n<=r_n` and `r_n/n^(3/2)->0`, or replace the `o` notation by an explicit
   uniform range.  Sampling an `r`-element set otherwise lacks an integrality
   declaration.
3. In the `L_n` corollary take
   `r=floor(n^(3/2)/L_n)` and, if necessary, decrease `c_1` by an absolute
   factor.  This preserves both the stated energy threshold and the exponent
   `Omega(sqrt(n)/L_n)`.
4. Replace “rules out any bounded witness-list description” by “rules out a
   bounded list that explicitly enumerates the whole shell.”  Exponential
   cardinality alone does not rule out a bounded-size algebraic, symmetric,
   or generative description; the preceding sentence in the draft correctly
   acknowledges this distinction.

Here is a precise replacement theorem statement with an explicit admissible
constant:

> Fix `c_0>0` such that `M_n>=c_0 n^(3/2)` for all sufficiently large `n`.
> Let `r_n` be an integer sequence satisfying `n<=r_n` and
> `r_n/n^(3/2)->0`.  Then, for all sufficiently large `n`, every exact
> minimizer `a` satisfies
> `|S_(2r_n/n^(3/2))(a)|>=exp((c_0^2/4)r_n/n)`.

Indeed, eventually `2r_n<=(c_0/2)n^(3/2)`, so every shell member has
population mean at least `c_0/sqrt(n)`.  Hoeffding--Serfling for a `+-1`
population then gives the explicit tail

```math
 P\left\{\sum_{e\in F}g_e(z)\le0\right\}
 \le \exp\left(-{c_0^2r_n\over2n}\right).
```

The union bound closes with `c_1=c_0^2/4`.  For the corollary, the copy-ready
choice `r_n=floor(n^(3/2)/L_n)` yields, for all sufficiently large `n`,

```math
 |\{z:\langle a,z\rangle\ge M_n-2n^{3/2}/L_n\}|
 \ge \exp\left({c_0^2\over8}{\sqrt n\over L_n}\right),
```

where the factor `1/8` only absorbs the floor.

With these repairs, setting `r=n^(3/2)/L_n` gives normalized shell width
`2/L_n` and logarithmic cardinality `Omega(sqrt(n)/L_n)`.  In particular,
`L_n=log n` proves the stated `exp(n^(1/2-o(1)))` count.  This is a count of
near-top signed cuts, not a response-metric packing and not an independent-
bit lower bound.

## 4. Final scope verdict

The audited results establish two different phenomena and should not be
conflated:

* MS.3--MS.4 give a common-family, externally queried one-hot packing with a
  polynomial number of states and only logarithmically many message bits.
* FB.1--FB.2 give a low-bias probability mixture and a finite empirical
  witness list; they do not encode the full response landscape.
* FB.3 gives a stretched-exponential *cardinality* lower bound for an exact
  minimizer's near-top shell; it supplies neither pairwise response
  separation nor independently decodable coordinates.

No audited theorem proves an independent-bit capacity, a `2^(k_j)` exact
edit cube, or simultaneous realization of all MS query shores in one
physical parent.
