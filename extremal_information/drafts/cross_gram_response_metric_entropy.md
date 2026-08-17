# Linear metric entropy for cross-Gram response states

**Status.** Rigorous task-local covering theorem and spherical-response
corollary.  The theorem concerns the relaxed PSD Gram--Rayleigh carrier.  It
does not prove an approximate quotient for the exact Boolean cap.

The exact `(G,R)` table has `Theta(p^2)` entries.  Nevertheless, in the
quadratic query metric relevant to all signed port sums, its fixed-accuracy
covering entropy is only linear in `p`.  The hard trust-region boundary loses
Lipschitz continuity but retains a sharp square-root modulus.

## 1. State space and query metric

Let

```math
\mathcal K_p={(K^+,K^-):K^\pm\succeq0,
                         \ \operatorname{tr}K^\pm\le p}.       \tag{GE.1}
```

Write

```math
G=K^++K^-,
\qquad R=K^+-K^-.                                    \tag{GE.2}
```

For two states define

```math
d_q((G,R),(G',R'))
={1\over p^2}
 \max_{\epsilon\in\{+-1\}^p,\ s\in\{+-1\}}
 \left|\epsilon^T(\Delta G+s\Delta R)\epsilon\right|. \tag{GE.3}
```

Equivalently, if

```math
q_p(K,L)={1\over p^2}\max_\epsilon
                   |\epsilon^T(K-L)\epsilon|,
```

then

```math
d_q=2\max\{q_p(K^+,K^{+'}),q_p(K^-,K^{-'} )\}.       \tag{GE.4}
```

Actual port states satisfy the stronger identities
`K^+ + K^- = G`, `diag(G)=1`, and
`tr(K^+)+tr(K^-)=p`; hence they form a subset of (GE.1).

## 2. Spectral truncation and a low-rank factor net

### Theorem GE.1 (linear cross-Gram metric entropy)

For every `0<eta<=1` and every `p>=1`, the relaxed state space has a
`d_q`-cover of radius `eta` with

```math
\boxed{
\log \operatorname{Cov}(\eta,\mathcal K_p,d_q)
\le
2p\left\lceil{4\over\eta}\right\rceil
 \log\left(1+{16\over\eta}\right).}                 \tag{GE.5}
```

The centres may be chosen as pairs of PSD matrices, each of rank at most
`ceil(4/eta)` and trace at most `p`.  In particular,

```math
\log \operatorname{Cov}(\eta,\mathcal K_p,d_q)
=O_\eta(p),                                           \tag{GE.6}
```

not `O(p^2)`.

#### Proof

It suffices to cover one sector in `q_p` to radius `eta/2`.  Let
`K\succeq0`, `tr K<=p`, and split spectrally

```math
K=K_{\rm hi}+K_{\rm lo},
```

where `K_hi` retains eigenvalues strictly greater than `eta p/4`.  Then

```math
\operatorname{rank}K_{\rm hi}
\le\left\lfloor{4\over\eta}\right\rfloor,
\qquad
\|K_{\rm lo}\|_{op}\le {\eta p\over4}.              \tag{GE.7}
```

Every Boolean query has squared norm `p`, so

```math
q_p(K_{\rm lo},0)\le\eta/4.                          \tag{GE.8}
```

Put `r_eta=ceil(4/eta)` and factor, after zero padding,

```math
K_{\rm hi}=BB^T,
\qquad B\in\mathbb R^{p\times r_\eta},
\qquad \|B\|_F^2=\operatorname{tr}K_{\rm hi}\le p. \tag{GE.9}
```

The Euclidean ball of radius `sqrt(p)` in dimension `p r_eta` has an
internal Frobenius net of radius `eta sqrt(p)/8` with at most

```math
\left(1+{16\over\eta}\right)^{p r_\eta}              \tag{GE.10}
```

points.  If `C` is the selected net point, both factors have operator norm
at most `sqrt(p)`, and

```math
\begin{aligned}
\|BB^T-CC^T\|_{op}
&\le(\|B\|_{op}+\|C\|_{op})\|B-C\|_{op}\\
&\le {\eta p\over4}.                                 \tag{GE.11}
\end{aligned}
```

Therefore `q_p(K_hi,CC^T)<=eta/4`.  Combining this with (GE.8) covers one
sector to radius `eta/2`.  Each centre `CC^T` is PSD, has rank at most
`r_eta`, and has trace `||C||_F^2<=p`.

Take the Cartesian product of the two sector nets.  By (GE.4) it is a
`d_q`-cover of radius `eta`; squaring (GE.10) proves (GE.5). `square`

No diagonal information, Boolean realization, or fixed-self-data hypothesis
was used.  The result is thus an upper bound for the larger relaxed space.
If centres must belong to a smaller realizable subset, choose one realizable
representative from every occupied net cell; this at most doubles the
covering radius.

## 3. Exact trust-coordinate formula

Let `H^2=r^2I`, and let `w_1,...,w_p` be port vectors of squared norm `n`.
For a port width `m`, put

```math
\mu={m\over r}.                                       \tag{GE.12}
```

For a signed endpoint word `epsilon` and child channel `sigma`, set

```math
g=\epsilon^TG\epsilon,
\qquad h=\epsilon^TR\epsilon,
\qquad
a=g+\sigma h\ge0,
\qquad b=g-\sigma h\ge0.                             \tag{GE.13}
```

After `t=2alpha-1`, the SA.3 spherical trust value divided by `rn` is

```math
\mathcal V_\mu(a,b)
={1\over2}+\inf_{t>0}
\left\{
 {t\over2}+{\mu^2a\over4t}
 +{\mu^2b\over4(t+2)}
\right\}.                                             \tag{GE.14}
```

The boundary value at `t=0` is understood by a limit when `a=0`.  The full
collective spherical response is

```math
\mathcal S_\mu(G,R)
=\max_{\epsilon,\sigma}\mathcal V_\mu(a,b).           \tag{GE.15}
```

Formula (GE.14) is valid for every `m>=0`.  The inherited SA.3 condition
`2m>r` is needed for the original one-port anti-pin compiler, not for the
trust-region dual identity itself.  This distinction matters below because
the compressible scaling has `m` of order `r/p`.

## 4. Global hard-edge modulus

### Lemma GE.2 (sharp trust-boundary continuity)

For `A,B,A',B'>=0`, define

```math
F(A,B)={1\over2}+\inf_{t>0}
\left\{{t\over2}+{A\over4t}+{B\over4(t+2)}\right\}.   \tag{GE.16}
```

Then

```math
|F(A,B)-F(A',B')|
\le
\sqrt{{|A-A'|\over2}}+{|B-B'|\over8}.                \tag{GE.17}
```

The square-root term is sharp: `F(A,0)-F(0,0)=sqrt(A/2)`.

#### Proof

Changing `B` at fixed `t` changes the objective by at most
`|B-B'|/[4(t+2)]<=|B-B'|/8`; taking infima preserves that bound.

Assume next that `A'=A+delta`.  For any `t>0`, replace `t` in the primed
objective by `t+s`.  The `B` term decreases, while

```math
{A+\delta\over4(t+s)}-{A\over4t}
\le{\delta\over4s}.                                  \tag{GE.18}
```

Hence

```math
F(A+\delta,B)-F(A,B)
\le {s\over2}+{\delta\over4s}.
```

Choosing `s=sqrt(delta/2)` gives `sqrt(delta/2)`.  Approximation and passage
to `t downarrow0` cover the hard case `A=0`; interchanging `A,A'` gives the
absolute bound.  Finally direct minimization of `t/2+A/(4t)` proves
sharpness. `square`

### Theorem GE.3 (when the query metric controls SA.3 response)

Let two Gram--Rayleigh states obey `d_q<=eta`, and put

```math
c=\mu p={mp\over r}.                                  \tag{GE.19}
```

Then their normalized collective spherical responses satisfy

```math
\boxed{
|\mathcal S_\mu(G,R)-\mathcal S_\mu(G',R')|
\le c\sqrt{\eta/2}+{c^2\eta\over8}.}                 \tag{GE.20}
```

#### Proof

For each common channel `(epsilon,sigma)`, definition (GE.3) gives

```math
|a-a'|\le\eta p^2,
\qquad |b-b'|\le\eta p^2.                            \tag{GE.21}
```

Apply Lemma GE.2 with `(A,B)=mu^2(a,b)`.  The resulting channelwise bound is
uniform, and taking the maximum over the same finite channel set preserves
it. `square`

Thus the exact amplification parameter is the **total repeated-port mass**
`c=mp/r`.  Fixed `d_q` accuracy controls the normalized response uniformly
when `mp=O(r)`.  The theorem gives no dimension-free response control when
`mp/r` diverges.

The square-root loss cannot be removed uniformly on the relaxed PSD cone.
Let `u=p^(-1/2)1`, compare the zero state to

```math
K^+={\eta p\over2}uu^T,
\qquad K^-=0.                                         \tag{GE.22}
```

Their query distance is `eta`, and the all-positive hard channel has
`mu^2a=c^2eta`, `b=0`.  Its response difference is exactly
`c sqrt(eta/2)`.

## 5. Linear continuity under a trust margin

The hard coordinate is the mass in the eigensector aligned with `sigma`.
Write `A=mu^2a`, `B=mu^2b`.  The objective in (GE.16) is convex, and its
minimizer lies at or beyond `t=tau>0` exactly when

```math
{A\over\tau^2}+{B\over(\tau+2)^2}\ge2.               \tag{GE.23}
```

Indeed, (GE.23) says that the derivative at `tau` is nonpositive.
A simple sufficient condition is `A>=2tau^2`.

### Proposition GE.4 (trust-margin Lipschitz law)

Suppose every relevant channel in both states has a trust minimizer
`t>=tau>0`.  Then

```math
|\mathcal S_\mu(G,R)-\mathcal S_\mu(G',R')|
\le {c^2\eta\over4}
 \left({1\over\tau}+{1\over\tau+2}\right).           \tag{GE.24}
```

In the original variable `alpha>=1/2+gamma`, this becomes

```math
{c^2\eta\over8}
\left({1\over\gamma}+{1\over1+\gamma}\right).        \tag{GE.25}
```

To compare the two infima, use the minimizer of either state in the other
objective.  On `t>=tau`, the coefficients of `A` and `B` are at most
`1/(4tau)` and `1/[4(tau+2)]`; (GE.21) proves (GE.24).  The reverse
comparison uses the other state's minimizer.

The hypothesis may be restricted to channels capable of attaining the two
outer maxima, but the all-channel statement is the clean checkable version.

## 6. Response-entropy corollary at fixed total port mass

### Corollary GE.5 (linear approximate spherical state)

Assume `mp<=r`, so `c<=1`.  For every `0<epsilon<=1`, the collective
spherical responses admit error-`epsilon` representatives with

```math
\log \operatorname{Cov}_{\rm response}(\epsilon)
\le
2p\left\lceil{4\over\epsilon^2}\right\rceil
\log\left(1+{16\over\epsilon^2}\right)
=O\left({p\over\epsilon^2}\log{1\over\epsilon}\right). \tag{GE.26}
```

Indeed, use the `d_q` cover at radius `eta=epsilon^2`.  By (GE.20), its
response error is at most

```math
{\epsilon\over\sqrt2}+{\epsilon^2\over8}<\epsilon.
```

The natural integral realization is `m=floor(r/p)` for `p<=r`.  Then the
total auxiliary width is at most `r=sqrt(n)`, so an arbitrary exact-sign
completion on all auxiliary vertices costs only `O(r^2)=O(n)`, negligible
against the `rn=n^(3/2)` scale.

Under a uniform trust margin, (GE.24) permits a `d_q` radius proportional to
`epsilon` instead of `epsilon^2`, improving the dependence on accuracy while
retaining linear dependence on `p`.

## 7. Scope and research consequence

1. **What is proved.**  The entire relaxed Gram--Rayleigh pair has
   fixed-scale query-metric entropy `exp(O_eta(p))`.  This is stronger than
   the affine-coordinate cube bound CG.2: no codebook or affine structure is
   assumed.
2. **Why it controls a collective query.**  At total port mass `mp=O(r)`,
   the same metric controls the maximum over all `2^(p+1)` spherical trust
   channels with a sharp square-root modulus.  Exponentially many channels
   cost nothing further because the metric is already uniform over them.
3. **What is not proved.**  The spherical relaxation may have a fixed
   leading integrality gap from the exact Boolean old-spin cap.  GE.5 is not
   an exact-sign Boolean response quotient and does not repair the BCX
   compositional congruence.
4. **Why the original SA scaling differs.**  With `m=r`, one has `c=p`; the
   modulus (GE.20) amplifies state error by `p`.  The fixed-radius entropy
   theorem alone then does not give an `exp(O(p))` fixed-response-error
   quotient.  Compression comes from sharing a fixed total auxiliary mass,
   not from the Gram table in isolation.
5. **Next theorem.**  Determine whether Boolean trust values obey an
   analogous modulus under a verifiable rounding/stability hypothesis, or
   construct a bounded-integrality-gap family which separates the spherical
   quotient from every exact Boolean response state.

## 8. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_cross_gram_response_metric_entropy.py
```

The verifier checks the metric/sector identity, spectral truncation ranks
and tail bounds, the exact `t`-coordinate formula, the hard-edge equality,
and randomized instances of (GE.17), (GE.20), and the trust-margin criterion.
