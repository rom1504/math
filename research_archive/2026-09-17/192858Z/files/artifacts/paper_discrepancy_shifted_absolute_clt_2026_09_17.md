# Uniform shifted absolute-value CLT at order 1/q

2026-09-17 campaign. A self-contained discrepancy-track proof prompted by
the Bernoulli track's proposed `3/q` bound for the sparse sign compiler.
The numerical test is independent evidence, not part of the proof.

## 1. Finite theorem

Let ε₁,…,ε_q be independent fair signs, let Z_q=q^(-1/2)Σ ε_i,
and let G be a standard Gaussian. For every integer q≥1,

$$
\boxed{\quad
\sup_{t\in\mathbb R}
\left|\mathbb E|Z_q-t|-\mathbb E|G-t|\right|
\le \frac{C_0}{q}<\frac{3}{2q},\qquad
C_0=\sqrt{\frac2\pi}
 \left(1+\frac{\sqrt2}{4}+\frac4{\pi^2}\right)
=1.403349785\ldots .\quad}                         \tag{1}
$$

In particular the proposed constant 3 is valid, with room to spare.
Equivalently, for unnormalized sums and arbitrary deterministic shifts a,

$$
\sup_a\left|\mathbb E|\textstyle\sum_i\varepsilon_i-a|
       -\mathbb E|\sqrt q G-a|\right|\le C_0/\sqrt q. \tag{2}
$$

### Proof

The elementary absolute-value Fourier identity is

$$
|x|=\frac2\pi\int_0^\infty\frac{1-\cos(ux)}{u^2}\,du.
$$

Applying it to the two symmetric random variables and subtracting gives,
uniformly over the shift t,

$$
\left|\mathbb E|Z_q-t|-\mathbb E|G-t|\right|
\le\frac2{\pi\sqrt q}\int_0^\infty
 \frac{|\cos(v)^q-e^{-qv^2/2}|}{v^2}\,dv.             \tag{3}
$$

The substitution is u=√q v. The difference is integrable both at zero
and at infinity; alternatively subtract truncated Fourier integrals and
pass to the limit. The omitted shift factor has absolute value at most 1.

Assume first q≥2. On 0≤v≤π/2, write a=exp(−v²/2), b=cos v.
Then 0≤b≤a and

$$
0\le a-b\le v^4/8.
$$

Indeed cos v≥1−v²/2, exp(−z)≤1−z+z²/2 for z≥0, and
log cos v≤−v²/2. Consequently the central-lobe integral is at most

$$
\int_0^{\pi/2}\frac{a^q-b^q}{v^2}\,dv
\le\frac q8\int_0^\infty v^2e^{-(q-1)v^2/2}\,dv
=\frac q8\sqrt{\frac\pi2}(q-1)^{-3/2}.              \tag{4}
$$

After multiplying by 2/(π√q), this contributes at most

$$
\frac1q\frac{\sqrt{2/\pi}}8
 \left(\frac q{q-1}\right)^{3/2}
\le\frac1q\sqrt{\frac2\pi}\frac{\sqrt2}{4}.          \tag{5}
$$

On the remaining half-line, triangle inequality separates the Gaussian
tail and the cosine tail. The former is bounded by

$$
\int_{\pi/2}^\infty\frac{e^{-qv^2/2}}{v^2}\,dv
\le\frac4{\pi^2}\sqrt{\frac\pi{2q}},                 \tag{6}
$$

so its contribution to (3) is at most
`q^(-1) sqrt(2/pi) 4/pi²`.

For the cosine tail, partition into the intervals
I_k=[(k−1/2)π,(k+1/2)π], k≥1. On I_k,

$$
|\cos v|^q\le e^{-q(v-k\pi)^2/2},\qquad
v^{-2}\le[\pi(k-1/2)]^{-2}.
$$

Therefore

$$
\int_{\pi/2}^\infty\frac{|\cos v|^q}{v^2}\,dv
\le\sqrt{\frac{2\pi}q}
 \sum_{k\ge1}\frac1{\pi^2(k-1/2)^2}
=\frac12\sqrt{\frac{2\pi}q}.                        \tag{7}
$$

Its contribution to (3) is `sqrt(2/pi)/q`. Equations (5)–(7)
give exactly C₀/q, for both even and odd q.

Finally, for q=1 couple ε=sign G. Every shifted absolute value is
1-Lipschitz, and hence its expectation discrepancy is at most

$$
\mathbb E\big||G|-1\big|
\le\sqrt{\mathbb E(|G|-1)^2}
=\sqrt{2-2\sqrt{2/\pi}}<1<C_0.
$$

This proves (1), and scaling gives (2). No limit theorem, smooth test
function approximation, or numerical certification is used. □

## 2. Exhaustive finite-shift stress test

The script
[shifted_abs_clt.py](../computations/paper_discrepancy_2026_09_17_shifted_abs_clt.py)
checks every q=1,…,1000. For each q, it examines every nonnegative atom
and every stationary point in every inter-atom gap. This is an exhaustive
list of possible finite extrema, since the derivative of the difference
is 2(F_q(t)−Φ(t)), and the binomial CDF is constant on each gap.

For numerical stability it evaluates the stop-loss difference from upper
tail probabilities and upper tail first moments. At a stationary point,
the shift-times-tail terms cancel exactly. It checked 319,809 extrema.
The largest observed q-scaled discrepancy was

$$
0.20442252894471957\quad(q=3,\ t=0).
$$

The largest observed even-q value was 0.19944617514877638 at q=1000,
t=0. Arithmetic is floating point, so these observations do not certify
the optimal constant or that the maximum is always attained at zero.
The analytic theorem above needs neither assertion.
