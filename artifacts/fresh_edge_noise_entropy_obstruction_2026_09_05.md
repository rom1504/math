# Edge noise, microcanonical entropy, and a sharp annealed obstruction

Date: 2026-09-05. Status: proved sufficient criterion and proved obstruction
to that criterion. Neither result determines the typical noisy optimum.

## 1. Normalization and exact one-state rate

Let `A_n` be symmetric hollow sign matrices and write

\[
 e_n(x)=\frac{H_{A_n}(x)}{n^{3/2}},\qquad
 H_A(x)=\sum_{i<j}a_{ij}x_ix_j.
\]

Independently multiply each undirected edge by a sign `xi_ij` having
`Pr(xi_ij=-1)=p`, where `0<p≤1/2` is fixed. Set `kappa=1-2p` and denote
the resulting sign matrix by `A_n'`. For each deterministic spin `x`,

\[
 H_{A_n'}(x)=\kappa H_{A_n}(x)+Z_{n,x},\qquad
 Z_{n,x}=\sum_{i<j}a_{ij}x_ix_j(\xi_{ij}-\kappa).
\]

For every fixed real `lambda`, uniformly in the underlying signing and spin,

\[
 \frac1n\log E\exp\left(\lambda Z_{n,x}/\sqrt n\right)
   =p(1-p)\lambda^2+o(1).                                      \tag{1}
\]

Indeed each summand is bounded by two, has mean zero and variance
`4p(1-p)`. Its logarithmic moment generating function at `lambda/sqrt(n)`
is its quadratic term plus `O_p,lambda(n^{-3/2})`; there are
`n(n-1)/2` terms. Therefore, for fixed `t>0`,

\[
 \Pr\{Z_{n,x}\ge t n^{3/2}\}
   =\exp\left[-n\frac{t^2}{4p(1-p)}+o(n)\right],               \tag{2}
\]

and likewise for the lower tail. The upper bound is Chernoff. For the
matching lower bound, exponentially tilt by a parameter whose limiting
mean is `t+epsilon`. Under that tilt the normalized variable has variance
`O(1/n)`, so a fixed interval around `t+epsilon` has probability tending
to one. Undo the tilt and let `epsilon` tend to zero. All bounds remain
uniform when `t` ranges in a fixed compact subset of `(0,infinity)`.

Thus the Gaussian-looking exponent is the exact speed-`n` moderate
deviation rate, not merely a variance heuristic.

## 2. An explicit sufficient microcanonical criterion

Suppose `|e_n(x)|≤E+o(1)` uniformly. Let `S:[0,E]→[0,log 2]` be an
upper-semicontinuous upper entropy profile: for each closed interval `I`,

\[
 \limsup_n\frac1n\log\#\{x:|e_n(x)|\in I\}
 \le\sup_{e\in I}S(e).
\]

If a proposed constant `c` satisfies `c>kappa E` and

\[
 \boxed{\qquad
 \sup_{0\le e\le E}
 \left[S(e)-\frac{(c-\kappa e)^2}{4p(1-p)}\right]<0,
 \qquad}                                                     \tag{3}
\]

then with probability tending to one,

\[
 \max_x |H_{A_n'}(x)|<c n^{3/2}.
\]

To prove this, partition the compact energy range into finitely many small
bins and apply (2), a union bound within each bin, and then a union bound
over bins and the two tail signs. The strict margin makes fixed sufficiently
fine bins legitimate. Equivalently, the profile-based threshold is

\[
 c_{\rm ann}(p)=\sup_e\bigl[\kappa e+
                  2\sqrt{p(1-p)S(e)}\bigr].                  \tag{4}
\]

Strictly exceeding this threshold gives (3), subject to the endpoint
condition. An exact entropy profile can improve a global `2^n` union bound,
but the next theorem shows why it cannot cross `1/2` when the starting
family itself approaches `1/2`.

## 3. Unavoidable entropy from a single near-extremizer

Assume `max_x |H_{A_n}(x)|=O(n^{3/2})`, and suppose some spins `x_n` satisfy
`e_n(x_n)→e_0>0`. For every fixed `0<delta≤1/2`, putting `r=1-2delta`,

\[
 \boxed{\quad
 \lim_{eta\downarrow0}\liminf_n\frac1n
 \log\#\{x:|e_n(x)-r^2e_0|<eta\}\ge h(delta),
 \quad}                                                      \tag{5}
\]

where `h(delta)=-delta log(delta)-(1-delta)log(1-delta)`.

Proof. Switch so that `x_n` is the all-one spin, and sample independent
coordinates `z_i` with mean `r`. Write `epsilon_i=z_i-r` and
`v=1-r²`. The degree-one and degree-two parts are orthogonal, giving exactly

\[
 E H_A(z)=r^2 H_A(1),\qquad
 \operatorname{Var}H_A(z)
   =r^2v\|A1\|_2^2+v^2\frac{n(n-1)}2.                       \tag{6}
\]

The low-cap hypothesis implies `||A||op²=O(n^{3/2})`. One elementary proof
uses `||A||op²≤beta(A)≤4 max_x|H_A(x)|`, where
`beta(A)=max_{x,y in {±1}^n}|x^TAy|`. The first inequality follows by bounding
the maximum absolute row sum of `A²` using a row of `A` in the cube; the
second follows by polarization on disjoint supports. Hence the variance
in (6) is `O(n^{5/2})=o(n³)`. For a spectral-cap family it is even `O(n²)`.

Consequently the sampled energy converges in probability to `r²e_0`.
Simultaneously the empirical fraction of negative coordinates converges to
`delta`. A typical sample has point probability
`exp[-n h(delta)+o(n)]`. Intersecting the two typical events, which has
probability tending to one, proves the cardinality bound (5). No count of
bent functions or exact eigenvectors is needed.

## 4. Why no exact entropy union bound crosses one half

Now take `e_0=1/2`, as holds for any asymptotically saturated spectral-cap
family, after reversing the whole signing if needed. In (5), set
`delta=p`. Then `r=kappa`, so states with energy near

\[
 e_p=\kappa^2/2
\]

have entropy at least `h(p)`. At the proposed noisy upper constant `c=1/2`,
their tail rate is

\[
 I_p=\frac{(1-\kappa^3)^2}{16p(1-p)}
     =\frac{(1-\kappa)(1+\kappa+\kappa^2)^2}{4(1+\kappa)}.
\]

For `0≤kappa≤1`,

\[
 I_p\le\frac9{16}(1-\kappa^2)
        =\frac94p(1-p),                                    \tag{7}
\]

because `(1+kappa+kappa²)/(1+kappa)≤3/2`. Meanwhile binary entropy obeys

\[
 h(p)\ge4\log(2)p(1-p).                                    \tag{8}
\]

For completeness, expand in `kappa=1-2p`:

\[
 h(p)=\sum_{j\ge1}\frac{1-\kappa^{2j}}{2j(2j-1)}
      \ge(1-\kappa^2)\sum_{j\ge1}\frac1{2j(2j-1)}
      =(1-\kappa^2)\log2.
\]

Since `4 log2>9/4`, (7)--(8) imply the strict, explicit deficit

\[
 \boxed{\quad
 h(p)-I_p\ge(4\log2-9/4)p(1-p)>0.
 \quad}                                                     \tag{9}
\]

Thus the actual entropy profile, not merely a coarse bound on it, violates
(3) at `c=1/2` for every fixed flip probability. It certainly violates the
criterion for any smaller `c`. In fact (2) and (5) show that the expected
number of spins whose noisy energy exceeds `(1/2)n^{3/2}` grows
exponentially. This is an annealed statement; high correlation between the
events can prevent it from deciding the typical maximum.

If `p=p_n→0`, a fixed near-maximizer itself has noisy energy converging in
probability to `1/2`, by its expectation and variance. Such vanishing noise
therefore cannot yield a fixed positive improvement below `1/2` with
probability tending to one either. Probabilities above one half are reduced
to those below one half by globally negating the noisy matrix.

## 5. Consequence and remaining escape

Independent edge noise plus any *one-state* microcanonical first-moment
certificate cannot improve the original `1/2` upper bound starting from a
family saturating it. Hamming-band entropy already precludes the proposal;
more detailed counts of exact Walsh extremizers cannot remove this obstacle.

This is not a proof that an appropriately chosen noisy signing never
improves the upper bound. A successful construction would have to exploit
joint overlap information, atypical noise, or a mechanism beyond the
one-state entropy union bound. If a starting family is already bounded away
from `1/2`, it already provides the desired upper improvement and the present
obstruction does not apply in its saturated form.
