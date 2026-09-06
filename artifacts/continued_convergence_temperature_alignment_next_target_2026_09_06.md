# The exact information slack behind temperature alignment

Date: 2026-09-06. This note isolates one remaining sufficient inequality
for `B E_t<=E_t`. It proves the information identity, the optimal-channel
moment restriction, and the Gaussian-residual case. It does not prove
the general inequality or a new bound on the original minimax problem.

The latent envelope, fixed-temperature supersolutions, and unconditional
Gaussian-boundary Bellman limit are in the companion convergence artifacts.
Write `g=g_t`, and `E(nu)=sup_L[g(E Var(X|L))-I(X;L)]`.

## 1. An exact two-label identity

Take a Bellman pair `(A,B)` with both signed marginals `nu`, after the
always-safe global reversal and input-swap symmetrization. Set
`U=(A+B)/sqrt(2)`, `V=(A-B)/sqrt(2)`. Choose child channels `M|U` and
`N|V`, sampled independently conditional on `(U,V)`, and put `L=(M,N)`.
Let

`V_+=E Var(U|M)`, `V_-=E Var(V|N)`,
`V_0=[E Var(A|L)+E Var(B|L)]/2`.

Choose an independent fair selector `R`, take `X=A` or `B` according to
it, and include `R` in the parent latent label. The resulting source is
again `nu`, and its information is exactly

`I(X;L,R)`
` =[I(U;M)+I(V;N)+I(A;B)-S]/2`,                            (1)

where

`S=I(M;N)+I(A;B|M,N)>=0`.                                 (2)

Indeed local conditional independence gives
`I(U,V;M,N)=I(U;M)+I(V;N)-I(M;N)`, while
`I(A;L)+I(B;L)=I(A,B;L)+I(A;B)-I(A;B|L)`.
Orthogonal invariance identifies the two joint informations. Also

`V_0=[E Var(U|L)+E Var(V|L)]/2 <=(V_++V_-)/2`.              (3)

Consequently the actual constructed parent channel proves

`E(nu)>=g(V_0)-[I(U;M)+I(V;N)+I(A;B)]/2+S/2`.             (4)

The child-channel contribution to `B E` would be bounded by this if

`S >= g(V_+)+g(V_-)-2g(V_0)`.                              (5)

The stronger but simpler sufficient test replaces `V_0` on the right by
`(V_++V_-)/2`, because `g` is decreasing. Establishing (5) for suitable
optimal or almost-optimal child channels, or showing that a refinement
of the constructed parent label always compensates for its failure,
would prove temperature alignment. A failure of (5) for a particular
pair of channels alone would not refute `B E<=E`: another parent channel
could give the necessary larger value.

## 2. Why arbitrary channels do not satisfy this test

Take a sparse ternary source
`nu_epsilon=(1-epsilon)delta_0+(epsilon/2)delta_(+1/sqrt(epsilon))`
` +(epsilon/2)delta_(-1/sqrt(epsilon))`, let `A=B`, and use constant
child labels. Then

`S=H(nu_epsilon)->0`, `V_+=2`, `V_-=0`, `V_0=1`.

For fixed `t>0`, strict convexity of `g` gives
`g(2)+g(0)-2g(1)>0`. Hence (5) fails for all sufficiently small
`epsilon`. In this example the child channel on the sparse nonzero
output is not an optimal envelope channel: encoding the rare large
component greatly improves its value. This falsifies only an unrestricted
variance-information shortcut, not the optimal-channel target.

## 3. An exact restriction on optimal quadratic rate-distortion posteriors

For a finite source, let a reproduction distribution `q` minimize

`J_lambda(nu)=inf_q -integral log[(k_lambda q)(x)] dnu(x)`,
`k_lambda(x,y)=exp[-lambda(x-y)^2]`.

The reproduction range may be restricted to the convex hull of the
source, so a minimizer exists. Its directional optimality condition is

`Z(y)=integral k_lambda(x,y)/(k_lambda q)(x) dnu(x)<=1`

for every real `y`, with equality at `q`-almost every active reproduction
point. Although the search range can be restricted to the convex hull,
the same inequality holds outside it: projection onto that hull only
increases the kernel pointwise.

At an active point `y`, the posterior source law is
`dnu_y(x)=k_lambda(x,y)dnu(x)/(k_lambda q)(x)`. Comparing `Z(y+h)` to
`Z(y)=1` gives the exact tilted-moment inequality

`E_(nu_y) exp[2lambda h(X-y)] <= exp(lambda h^2)`.

It follows that `E[X|y]=y` and

`log E[exp(s(X-y))|y] <= s^2/(4lambda)`                    (6)

for every real `s`. Thus every active posterior is subGaussian with
variance proxy `1/(2lambda)`. This is stronger than the second-derivative
necessary condition alone.

For a differentiable interior envelope optimizer, the Gaussian supporting
line relation additionally gives

`lambda=-g_t'(V)=t(1-rho)`,
`2lambda V=(t-lambda)/(2t-lambda)<=1/2`.                    (7)

The proxy in (6) can therefore be much larger than the actual average
residual variance. Equation (6) must not be silently replaced by a sharp
variance-proxy bound. The unresolved question is whether this exact
optimal-channel structure, together with (1), pays the Jensen gap in
(5), perhaps after a useful parent-label refinement.

## 4. A sharp all-temperature Gaussian-residual comparison

The Gaussian potential satisfies

`v^2 g_t''(v)=rho^2/[2(1+rho^2)]<=1/4`,
`2tv=rho/(1-rho^2)`.

Thus `g_t(v)+(1/4)log v` is concave on `v>0`. For all `a,b>0`,

`g_t(a)+g_t(b)-2g_t((a+b)/2)`
` <=(1/4)log[(a+b)^2/(4ab)]`.                              (8)

Suppose the chosen valid child channels have conditionally independent
Gaussian residuals of fixed variances `a,b`: conditional on `(M,N)`,
`U` and `V` are independent Gaussians whose means depend only on their
own labels `M,N`, respectively, and whose variances are `a,b`. Then
`V_+=a`, `V_-=b`, `V_0=(a+b)/2`, while the rotated residuals give

`I(A;B|M,N)=(1/2)log[(a+b)^2/(4ab)]`.

By (2), this is a lower bound on `S`, and it is twice the upper bound
in (8). Therefore (5) holds throughout this class, with a factor-of-two
margin whenever the two variances differ. In particular this settles
the proposed slack test for independent Gaussian residuals at every
temperature and variance ratio, not merely asymptotically.

The non-Gaussian optimal-posterior extension is the precise remaining
target. None of these identities asserts that a general subGaussian
conditional law has the same variance-ratio mutual-information lower
bound as a Gaussian law.

For the specific restricted-weave application, a universal theorem on
all finite-second-moment laws is stronger than necessary. Starting from
the ternary source, every finite Bellman tree has finite symmetric
lattice-valued states, with the lattice rescaled by `sqrt(2)` at each
level. Alignment on this reachable state class suffices for its finite
trees and the stopped-tree upper argument. An alignment counterexample
outside that class would need a separate reachability or approximation
argument before it could rule out the live construction.
