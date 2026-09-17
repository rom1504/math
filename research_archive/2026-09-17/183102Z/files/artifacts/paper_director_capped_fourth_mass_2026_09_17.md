# Capped fourth-mass Laplace envelope and the critical-frame regime

2026-09-17. **Verified: self-contained derivation and independent BH-track
reconstruction, including sharpness and critical-scale normalization.** This targeted
fourth paper was added because it addresses the precise fourth-order
quantity in the new common-Gibbs universality bound:
[Gao--Qian, Fourth-Moment Geometry of Rademacher Sums,
2608.17802v1](https://arxiv.org/html/2608.17802v1).
Only its Laplace-transform mechanism is reconstructed here; the paper's
other moment and stability theorems are NOT imported without audit.
External novelty of the extension below has not been established.

## 1. A self-contained reconstruction and capped extension

Let S=sum_i a_i epsilon_i, with independent fair signs. Put

```math
v=\sum_i a_i^2,\quad w=\sum_i a_i^4,\quad
B\ge\max_i|a_i|>0,\quad m=\lfloor w/B^4\rfloor,\quad
r=w-mB^4\in[0,B^4).
```

Then for every real t,

```math
\boxed{\log\mathbb Ee^{tS}\le
 \frac{t^2}{2}(v-mB^2-\sqrt r)
 +m\log\cosh(tB)+\log\cosh(t r^{1/4}).}              \tag{1}
```

The Gaussian variance v-mB^2-sqrt(r) is nonnegative for every feasible
coefficient vector. The weaker but convenient consequence is

```math
\boxed{\log\mathbb Ee^{tS}\le\frac{t^2v}{2}
 -\frac{w}{B^4}\left[\frac{t^2B^2}{2}-\log\cosh(tB)\right].} \tag{2}
```

**Proof.** For u>=0 set
f_t(u)=t^2u/2-log cosh(t sqrt(u)), and g_t(z)=f_t(sqrt(z)).
The function g_t is nonnegative, vanishes at zero, and is concave.
To check the only nontrivial assertion, assume t>0, put u=sqrt(z)
and b=t sqrt(u), and differentiate:

```math
g_t''(z)=\frac{t^2}{16u^3}
 \left[3\frac{\tanh b}{b}-\operatorname{sech}^2 b-2\right]\le0.
```

Indeed L(b)=2b-3 tanh b+b sech^2 b satisfies L(0)=0 and
L'(b)=2 tanh b[tanh b-b sech^2 b]>=0, since
sinh b cosh b>=b. The displayed bracket is -L(b)/b.
Evenness treats negative t, and t=0 is immediate.

At fixed sum w, a concave function's sum over coordinates in [0,B^4]
is minimized by packing the coordinates at the endpoints, with at
most one remainder: m coordinates B^4 and one r. Therefore

```math
\sum_i g_t(a_i^4)\ge m g_t(B^4)+g_t(r).
```

Subtract from t^2v/2 to prove (1). The same endpoint-packing statement
for sqrt(z), another concave function, gives
sum_i a_i^2>=mB^2+sqrt(r). Finally concavity gives
g_t(z)>=z g_t(B^4)/B^4 on [0,B^4], proving (2).

Without an active coefficient cap, take B>=w^(1/4). Formula (1)
reduces exactly to the paper's one-Bernoulli-spike plus Gaussian
Laplace envelope. Thus the Laplace specialization has a short direct
proof independent of the paper's more general moment-extremizer argument.

## 2. Sharpness and limit order

For fixed feasible v,w,B, the right side of (1) is the MGF of m signs
with coefficient B, one with coefficient r^(1/4), and an independent
Gaussian of variance a=v-mB^2-sqrt(r). It is the supremum in the
closure of finite sign sums with these exact v,w and the same cap B.

If a=0 the displayed finite sum attains it. If a>0, choose one positive
packed squared coefficient u0 (sqrt(r) when r>0, otherwise B^2).
For large integer L solve

```math
u_L+a_L=u_0+a,\qquad u_L^2+a_L^2/L=u_0^2,
```

using the solution u_L<=u0 tending to u0. Replace that coefficient
by sqrt(u_L) and add L coefficients sqrt(a_L/L). All second and fourth
masses are EXACT, all coefficients obey B for large L, and their MGFs
converge to the asserted packed-sign/Gaussian MGF. This limit concerns
the number of summands, not the order of an original optimizing matrix.
The zero vector is trivial and excluded from the sharpness construction.

## 3. Why the cap matters at the new transfer boundary

For a balanced Hadamard-mode bridge with n old and q new vertices,
assume k divides both n and q to display the normalization cleanly.
For each pair of Boolean words x,y, its independent-sign expansion has

```math
V_{x,y}=qn,\qquad |c_b(x,y)|\le k,
\qquad W_4(x,y)=\frac qk\sum_{a,t}|(h_a^Tx_{:,t})|^4.
```

In particular W4 is independent of y. With psi(z)=z^2/2-log cosh z,
(2) gives the exact physical-sign MGF bound

```math
\log\mathbb E e^{t x^TCy}
 \le\frac{t^2qn}{2}-\frac{W_4(x)}{k^4}\psi(tk).       \tag{3}
```

Suppose q/n->epsilon, k/sqrt(n)->lambda>0 and
W4(x)/n^3>=omega. At t=beta/sqrt(n), (3) has normalized exponent

```math
\frac1n\log\mathbb Ee^{\beta x^TCy/\sqrt n}
 \le\frac{\epsilon\beta^2}{2}
 -\frac{\omega}{\lambda^4}\psi(\beta\lambda)+o(1).    \tag{4}
```

Thus a macroscopic fourth-mass query receives a macroscopic negative
pressure correction. The uncapped one-spike envelope permits a spurious
coefficient of size W4^(1/4), and only subtracts O(sqrt(n)) from the
raw exponent when W4 is order n^3; it misses this O(n) correction.
The physical cap k is essential information, not a cosmetic parameter.

Conversely the independently audited common-Gibbs theorem, combined
with querywise large-coefficient truncation, transfers ANY declared
query sector with max W4=o(n^3) and exp(O(n)) queries to its Gaussian
model at o(n^(3/2)) error. The explicit bound is
8 W4^(1/4)[log(2K)]^(3/4), without a coefficient-cap hypothesis.
This DOES give zero-temperature transfer on diffuse sectors at
k=Theta(sqrt(n)). It supersedes an earlier draft's smoothing caveat;
see paper_symmetric_frame_universality_2026_09_17.md, Section1.1.

The combination retains the exact correction (3) on coherent sectors
and uses Gaussian transfer on diffuse sectors. A favorable full-parent bound still needs
actual old-child deficits or sector counts. This file does not assume
them and does not change the original signing interval.

## 4. Explicit transfer boundary to be investigated

Pointwise Laplace improvement is not by itself an upper bound for the
maximum over all words. One must pay their joint entropy, retain both
polarities, and keep the old child energies as offsets. The correct
annealed identity uses W4(x) inside the old child's partition sum;
replacing it by a favorable value selected separately for each extremum
would be invalid. The next step is a rigorous all-query certificate
with the exact losses displayed, followed by a test on actual sign
children. No claim of convergence is made from the scalar envelope.
