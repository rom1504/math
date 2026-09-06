# Generic original-signing gain from genuine fractional slack

Date: 2026-09-06. Status: the director, convergence agent, and bound-audit
agent independently reconstructed
conditional rounding, the smooth Lindeberg estimate, the Gaussian-excess
bound, and the exact quadratic gain.
This applies to arbitrary bounded-op hollow signings. It also explains
why it does not directly improve the currently banked pure policy.

## 1. A conditional Gaussian field requiring no tree or Haar theorem

Let B=A/sqrt(n-1) be a symmetric hollow signing with ||B||op<=L. Let
u in [-1,1]^n be any mean vector, possibly computed from B and previous
randomness. Condition on that information and independently round to
Boolean v with E[v_i|u]=u_i. Put

    d=1-||u||^2/n,
    a=Bu,
    epsilon=v-u.

For every row i, exactly

    Var((B epsilon)_i |u)
       = [n d-(1-u_i^2)]/(n-1)=d+O(1/n),        (1)

uniformly in i. Summands are independent, centered, and bounded by
2/sqrt(n-1). In fact an elementary smooth Lindeberg replacement gives

    E[|(Bv)_i| |u]=Psi(a_i,d)+O(n^(-1/6))       (2)

uniformly in i and ALL u, including zero variance. Here is a direct
proof avoiding a variance-dependent Berry--Esseen invocation. Smooth
|x| to sqrt(x^2+epsilon^2), with uniform function error at most epsilon
and third derivative bounded by C/epsilon^2. Replace each summand by an
independent Gaussian with the same variance. The sum of their third
absolute moments, original and Gaussian, is O(n^(-1/2)): use
E|v_j-u_j|^3<=2(1-u_j^2), and the maximal matrix coefficient.
Taylor expansion through second order therefore gives error
2epsilon+C epsilon^(-2)n^(-1/2), uniformly in the shift a_i. Take
epsilon=n^(-1/6). Finally the exact variance (1) differs from d by
O(1/n), so the Gaussian absolute-value expectations differ by O(n^(-1/2))
using |sqrt(s)-sqrt(t)|<=sqrt(|s-t|). No lower variance bound is needed
for (2); it is needed only for a fixed positive gain in Section 2.

Because B is HOLLOW, v_i is independent of (Bv)_i conditional on u.
Consequently

    E[v_i(Bv)_i |u]=u_i a_i.

There is no analogous claim here when the same rounding randomness
has already been fed back into u: conditioning precedes this NEW round.

## 2. Uniform gain for a mean with fixed positive slack

Assume 0<d_0<=d<=1, and define the exact positive constant

    g_0=Gamma_sqrt(d_0)(L sqrt(1-d_0))>0.

The expected actual stability gap of v is therefore

    g_n=(1/n) sum_i E[|(Bv)_i|-v_i(Bv)_i |u]
        =(1/n) sum_i [Psi(a_i,d)-u_i a_i]+o(1).

Since |u_i|<=1 and ||a||^2/n<=L^2(1-d), the convex-in-a_i^2 Gaussian
excess bound gives

    g_n >= Gamma_sqrt(d)(L sqrt(1-d))-o(1)
         >= g_0-o(1).                          (3)

Gamma_sigma(z)=E|z+sigma N|-z for z>=0 increases with sigma and decreases
with z. This is a lower bound, not an asserted exact distribution of
an adaptively generated old field.

Put w=sign(Bv), with any Boolean choice at zero, and choose the FIXED
eta=g_0/(4L). The feasible mean m=(1-eta)v+eta w obeys

    [Q_B(m)-Q_B(v)]/n
       >=eta (1/n)sum_i(|(Bv)_i|-v_i(Bv)_i)-2L eta^2.

Taking conditional expectations gives gain at least

    g_0^2/(8L)-o(1).                            (4)

Initial independent rounding has E Q_B(v)=Q_B(u) exactly because B is
hollow. Final independent Boolean rounding of m also preserves its
expected energy exactly. Thus (4) is a genuine original-class signing
improvement over any fractional mean having deterministic average slack
at least d_0. This improves the Boolean spin vector for the SAME fixed
coefficient matrix A; it does not edit that matrix's signs. If u is random,
the same assertion holds conditionally
when that slack condition holds; no such condition is inferred merely
from positive expected slack on an unspecified event.

This argument uses one new query Bv. It needs neither involution
structure nor the restricted forest-response transport identities.

## 3. Why the banked policy does not supply this slack

The banked theta=1 optimized two-dimensional rectangle policy has

    F in {0,+1,-1}, H=1-|F| in {0,1},
    C=H sign(BF).

F and C have disjoint supports and |F|+|C|=1 if the zero-field sign is
chosen Boolean. Therefore BOTH endpoints +/-F+C are Boolean and d=0
EXACTLY. Choosing sign(0)=0 merely creates slack on zero fields; no
positive average lower bound on that event is banked. Finite Lipschitz
approximants create slack which can vanish with approximation tolerance,
not a uniform d_0. Earlier relaxed mixtures were purified by an even
gate precisely to license their larger residual variance; the purified
terminal endpoints are likewise saturated.

One can deliberately damp a Boolean vector x to u=(1-rho)x. Then
d=2rho-rho^2 and its positive normalized energy e is reduced to
(1-d)e: the energy cost is d e. The generic Gaussian-only certificate
(4) cannot pay this cost at the current frontier. Indeed

    Gamma_sqrt(d)(L sqrt(1-d))<=sqrt(2/pi) sqrt(d),
    [Gamma_sqrt(d)(L sqrt(1-d))]^2/(8L)<=d/(4pi L).

Since L>=1 and .433322...>1/(4pi), this guaranteed increment is smaller
than the damping cost for every d>0. For d tending to zero the excess
is actually exponentially small when L stays positive.

This is a limitation of this conservative rounding certificate, not a
theorem that deliberate damping or a more detailed field calculation
cannot help. Extra pre-existing instability, nonuniform coordinate
damping, or new mixed transport might improve the estimate. The current
pure endpoint, however, does not itself furnish the fixed fractional
slack needed by (3)--(4).
