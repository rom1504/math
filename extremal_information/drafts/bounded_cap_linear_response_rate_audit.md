# Adversarial audit of the bounded-cap linear response-rate theorem

Audit target:
[`bounded_cap_linear_response_rate.md`](bounded_cap_linear_response_rate.md).

## Verdict

**PROMOTE, after the two proof-hardening edits stated in Section 8.**  I found
no false implication, normalization error, circular dependence, or missing
probabilistic intersection.  In particular, the proposed argument really does
give one common dense sign bridge and an `exp(gamma n)` projective response
packing inside the exact cap-`1/2` switching orbit.

The two requested edits are local rather than mathematical changes:

1. expand the invocation of multivariate Berry--Esseen into the covariance
   normalization and two-convex-set estimate given below, with a precise
   citation to Bentkus's independent, non-identically-distributed theorem;
2. formalize the packing-to-summary sentence by declaring an encoder, a common
   decoder, and either pointwise sup error or projective error.

There is also one duplicated sentence in the row-lemma proof.  None of these
repairs changes a constant or conclusion.

## 1. Child construction and normalization

Write `n=q^2=2^(2m)`.  For

```math
b(u,v)=(-1)^{u\cdot v},
```

direct summation gives `Wb=qb`.  Hence

```math
\mathcal H=D_bWD_b,
\qquad \mathcal H\mathbf 1=q\mathbf 1,
\qquad \mathcal H^2=nI.
```

Moreover

```math
\operatorname{tr}\mathcal H
=\operatorname{tr}W
=\sum_{a\in\mathbb F_2^{2m}}(-1)^{a\cdot a}=0
```

for `m>=1`.  Thus, if `A=mathcal H-diag(mathcal H)`, then for every Boolean
`u`,

```math
\frac12u^TAu
=\frac12u^T\mathcal Hu-\frac12\operatorname{tr}\mathcal H
=\frac12u^T\mathcal Hu.
```

The spectral bound and the all-one witness give

```math
\left|\frac12u^T\mathcal Hu\right|
\le\frac12\sqrt n\,\|u\|_2^2
=\frac12n^{3/2},
\qquad
\frac12\mathbf1^T\mathcal H\mathbf1
=\frac12n^{3/2}.
```

Therefore both the positive cap `P` and the absolute cap `Q(A)` equal
`n^(3/2)/2`.  The use of the **positive** cap in BC.1 is compatible with the
statement about the **absolute** cap.  Switching by `D_s` preserves
hollowness, signs, spectrum, and the exact absolute cap.

## 2. Hanson--Wright near-top count

Membership in `T_0` is equivalent to

```math
H_A(u)>\left(\frac12-\frac18\right)n^{3/2}
=\frac38n^{3/2},
```

and hence implies `u^T mathcal H u>(3/4)n^(3/2)`.  For uniform Rademacher
`U`,

```math
\mathbb E U^T\mathcal HU=\operatorname{tr}\mathcal H=0,
\qquad
\|\mathcal H\|_F=n,
\qquad
\|\mathcal H\|_{op}=\sqrt n.
```

Rudelson--Vershynin's Hanson--Wright theorem therefore gives

```math
\Pr\{U^T\mathcal HU>(3/4)n^{3/2}\}
\le 2\exp[-c_{HW}(9/16)n].
```

For example, any fixed `kappa<9c_HW/16` works after increasing `n_0` to
absorb the factor `2`.  Multiplication by `2^n` proves the stated entropy
gap.  The diagonal causes no problem: its contribution is the deterministic
quantity `tr(mathcal H)=0`.

The source theorem applies exactly in the required form: M. Rudelson and
R. Vershynin, *Hanson--Wright inequality and sub-gaussian concentration*,
Electronic Communications in Probability 18 (2013), Theorem 1.1,
<https://doi.org/10.1214/ECP.v18-2865>.

The proof correctly fixes `d_0=1/8` before choosing the eventual response
gap.  There is no hidden dependence of `kappa` on the later parameters.

## 3. Uniform row lemma

The row reduction is correct.  With `xi_i=R_i y_i` and
`c_i=y_i z_i`,

```math
(S,T)=\sum_i\xi_i(1,c_i),
\qquad
\operatorname{Cov}\frac{(S,T)}{\sqrt n}
=\Sigma_\theta
=\begin{pmatrix}1&\theta\\\theta&1\end{pmatrix},
\qquad
\theta=\frac{y^Tz}{n}.
```

Here is the missing explicit Berry--Esseen calculation.  Put

```math
X_i=n^{-1/2}\xi_i(1,c_i).
```

If `|theta|<=rho<1`, the least eigenvalue of `Sigma_theta` is at least
`1-rho`, and

```math
\sum_{i=1}^n
\mathbb E\|\Sigma_\theta^{-1/2}X_i\|_2^3
\le
\frac{2^{3/2}}{(1-\rho)^{3/2}\sqrt n}.
```

Bentkus's Lyapunov bound for convex Borel sets therefore compares `(S,T)`
with the centered Gaussian having covariance `n Sigma_theta`, uniformly in
the pair `y,z`, with error `C_rho/sqrt(n)` per convex set.  The relevant event
is the disjoint union of two convex quadrant/half-strip sets.  Thus its total
comparison error is at most `2C_rho/sqrt(n)`.  This also handles either fixed
zero convention: the corresponding strict or non-strict halfspaces are
convex Borel sets, while their Gaussian boundaries have zero mass.

For the Gaussian pair,

```math
\min_{\sigma\in\{-1,1\}}
\Pr\{\operatorname{sign}(G_1)\operatorname{sign}(G_2)=\sigma\}
=\frac12-\frac{|\arcsin\theta|}{\pi}.
```

Consequently the probability in BCL.3 is at least

```math
\frac12-\frac{\arcsin\rho}{\pi}
-\Pr\{|G|<a\}
-\frac{2C_\rho}{\sqrt n}.
```

Choosing `rho`, then `a`, then `n_0` proves the lemma for both requested
target signs.  This is exactly the uniformity required later when `t_i=u_i`
is adversarial.

The applicable primary result is V. Bentkus, *A Lyapunov-type bound in
R^d*, Theory of Probability and Its Applications 49 (2005), 311--323,
<https://doi.org/10.1137/S0040585X97981123>.  It explicitly permits
independent non-identically distributed summands and gives the convex-set
bound in terms of
`sum_i E|C^{-1}X_i|^3`.

## 4. Code and exponent choices

All parameter choices are feasible:

- as the row tolerance tends to zero,
  `alpha=1/2-epsilon` tends to `1/2` and
  `D(0||alpha)` tends upward to `log 2`, so BCL.16 can be met;
- `D(delta||alpha)` tends to `D(0||alpha)` as `delta` decreases to zero, so
  BCL.17 can then be met;
- decreasing `a` preserves the row lemma and allows
  `d=2a delta<d_0`;
- positive `gamma<min(rho^2/8,kappa/8)` exists.

A random code of `ceil(exp(gamma n))` independent Boolean words has a bad
pair with probability at most

```math
2\exp\{(2\gamma-\rho^2/2)n+o(n)\}=o(1).
```

The stricter choice `gamma<rho^2/8` leaves ample slack.  Duplicate and
antipodal words are automatically counted as bad pairs, so passing this
test really gives a set with the asserted projective correlation bound.

## 5. Weighted mismatch probability and union bounds

For fixed ordered `y!=z` and `u in T_0`, the row indicators are independent
and each has success probability at least `alpha`.  They therefore
stochastically dominate `Bin(n,alpha)`.  Since `delta<alpha`,

```math
\Pr\{N_{\rm good}<\delta n\}
\le \exp[-nD(\delta\|\alpha)].
```

Every successful row contributes at least `a sqrt(n)` to the weighted
mismatch, so a weighted sum strictly below
`a delta n^(3/2)=d n^(3/2)/2` implies this lower-tail event.  Nonintegrality
of `delta n` changes nothing.

The combined exponent is

```math
(\log2-\kappa)+2\gamma-D(\delta\|\alpha)
<2\gamma-\frac\kappa2
<-\frac\kappa4.
```

Thus the displayed union bound is valid for every `u in T_0` and every
ordered query pair.  A standard iid-subgaussian matrix norm tail gives
`Pr{||B||op>C sqrt(n)}<=2e^{-c_B n}`.  The union of that bad event with the
weighted-neighborhood bad event has probability below one for all large
`n`; independence between those two events is neither claimed nor needed.
Hence one deterministic **common** bridge satisfies both properties.

## 6. Switching, BC.1, and projective separation

For `s_z=sign(Bz)` and row sums `S_i=(By)_i`, `T_i=(Bz)_i`,

```math
\operatorname{sign}(s_{z,i}S_i)
=\operatorname{sign}(T_i)\operatorname{sign}(S_i),
\qquad |s_{z,i}S_i|=|S_i|.
```

Therefore the row event is exactly a weighted mismatch for the cross field
`s_z odot By`; it is not an independently paid surrogate channel.

For a spin outside `T_0`, the energy deficit is at least
`d_0 n^(3/2)>d n^(3/2)`.  For a spin in `T_0`, the twice-weighted mismatch
is at least `d n^(3/2)`.  BC.1 therefore proves each directed deficit.
The union bound was over ordered pairs, so it also proves the reversed
deficit needed in BCL.29.

If `R_y=P_BH_{A^{s_y}}`, direct change of variables gives

```math
(R_y-R_z)(y)=\Delta_A(s_z\odot By),
\qquad
(R_y-R_z)(z)=-\Delta_A(s_y\odot Bz).
```

Thus

```math
d_{\rm proj}(R_y,R_z)
=\frac12\operatorname{osc}(R_y-R_z)
\ge\frac12(\Delta_{z\to y}+\Delta_{y\to z})
\ge d n^{3/2}.
```

The factor `1/2` and both directed terms are correct.  If two switchings
collided (including the unavoidable identification `s~-s`), their response
functions would coincide, contradicting this positive separation.  Hence
the construction contains the asserted number of distinct children.

## 7. Information consequence and scope

The response scope is correctly limited to the fixed public bridge `B_n`
and Boolean external queries.  To make BCL.4 formally self-contained, state
the following standard packing consequence.

Let an encoder map each child to one of `K` states, and let a common decoder
map `(state,y)` to a predicted response.  If

```math
\sup_y|\widehat R(A,y)-P_BH_A(y)|
\le\varepsilon n^{3/2}
```

for every child in the constructed family, then two children sharing a state
would have projective distance at most `2 epsilon n^(3/2)`.  Therefore
`epsilon<d/2` forces `K>=exp(gamma n)`, or at least
`gamma n/log 2` bits in a worst-case fixed-length representation.  The same
argument works if approximation is declared directly in `d_proj`.

The theorem does **not** establish such a lower bound for a future class
that cannot implement or pin the Boolean queries.  It also does not lower
bound a representation allowed to quotient child and bridge jointly under
switching.  Those limitations are already stated accurately in the target.

## 8. Required repairs and nonissues

Before canonical promotion:

1. Replace the abbreviated Berry--Esseen paragraph by the normalized
   Bentkus calculation in Section 3, and cite the 2005 independent-summand
   theorem.  Remove the duplicated sentence about the standard-normal
   marginal.
2. Replace “a summary which answers all queries” by the explicit
   encoder/decoder statement in Section 7.  Rename the row-lemma tolerance
   (for example `epsilon_row`) so it is not confused with the later response
   approximation tolerance.

No repair is needed for the following possible concerns:

- positive versus absolute cap: the regular-Walsh child has both equal;
- lattice zeros: convex Borel halfspaces, or an `O(n^(-1/2))`
  anti-concentration term, handle either sign convention;
- dependence of switchings on the random bridge: the proof union-bounds the
  underlying row events before defining the children, so this is legitimate;
- simultaneous spectral norm: a final union bound gives a common bridge;
- hidden full-parent optimization: BC.1 is used only as an exact response
  identity, and no target-order optimum or unknown minimizer enters the
  construction;
- reconstruction of coefficients: the lower bound concerns contextual
  response states, not coefficient descriptions.

The exact regression
[`verify_bounded_cap_linear_response.py`](../experiments/verify_bounded_cap_linear_response.py)
also passes all `29,625` deterministic checks.  Those checks validate the
Walsh identities, BC.1, row reduction, and query-linked projective identity;
they are appropriately not presented as evidence for the asymptotic
concentration steps.
