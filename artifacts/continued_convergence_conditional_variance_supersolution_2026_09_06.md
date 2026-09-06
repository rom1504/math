# A conditional-variance envelope is a Bellman supersolution

Date: 2026-09-06. Define

`T_t(nu)=sup_L { E_L g_t(Var(X|L))-I(X;L) }`,

where `X` has a finite symmetric law `nu`. The Gaussian
potential is the same `g_t` as in the other recursive-pairing artifacts.
This differs from `E_t`: the Gaussian potential is averaged *after*
taking the classwise conditional variances, not evaluated at their
average. The director proposed this envelope; the sequential-label
factorization below proves its supersolution property.

## 1. Scalar properties

The function `g=g_t` is decreasing and convex, with `g(0)=0`. In addition
it is concave as a function of log variance:

`d/ds g(exp(s))=v g'(v)=-rho/[2(1+rho)]`,
`2tv=rho/(1-rho^2)`.

The last quantity decreases with `v`, so `s -> g(exp(s))` is concave.

For a positive definite two-by-two covariance matrix `Sigma`, let
`b=Sigma_22`, `a=det(Sigma)/b`, and let its eigenvalues be `lambda_1,
lambda_2`. The pair `(a,b)` has the same product as the eigenvalue pair
and is less spread on the log scale: each of `a,b` lies between the two
eigenvalues. Hence

`g(a)+g(b)>=g(lambda_1)+g(lambda_2)`.                       (1)

If `(u,v)` are the diagonal entries after any orthogonal rotation,
ordinary convexity and eigenvalue majorization give

`g(lambda_1)+g(lambda_2)>=g(u)+g(v)`.                       (2)

Singular cases follow by continuity. Directly, when `b=0`, the second
coordinate is deterministic and one uses `a=Sigma_11`; when `b>0` but
the determinant is zero, `a=0`.

## 2. Sequential parent labels

Take an admissible Bellman pair `(A,B)` whose two signed marginals are
both `nu`, using the safe reversal/swap symmetrization. Set
`U=(A+B)/sqrt(2)`, `V=(A-B)/sqrt(2)`. Choose child channels `M|U` and
`N|V`, independently conditional on `(U,V)`, and put `L=(M,N)`.

Use `(L,B)` as a parent label for `A`, and `L` as a parent label for
`B`. Their information costs satisfy the exact identity and bound

`I(A;L,B)+I(B;L)=I(A,B;L)+I(A;B)`
` <=I(U;M)+I(V;N)+I(A;B)`,                                (3)

because the difference in the second line is `I(M;N)>=0`.

Fix a value of `L`, and let `Sigma` be the conditional covariance of
`(A,B)`. Conditional linear regression gives

`E_(B|L) Var(A|B,L)<=det(Sigma)/Sigma_22`

when `Sigma_22>0`. Jensen and the fact that `g` decreases therefore show

`E_(B|L) g(Var(A|B,L))+g(Var(B|L))`
` >=g(det(Sigma)/Sigma_22)+g(Sigma_22)`
` >=g(Var(U|L))+g(Var(V|L))`,                              (4)

where the last line combines (1)--(2). The stated singular conventions
make the same inequality valid in degenerate cases.

Finally, refining a label increases the averaged Gaussian reward.
For example total variance and Jensen give

`E_(N|M) g(Var(U|M,N))`
` >=g(E_(N|M) Var(U|M,N))>=g(Var(U|M))`.                    (5)

Average (4), apply (5) and its counterpart for `V`, and subtract (3).
The two parent channels together have envelope value at least

`E_M g(Var(U|M))+E_N g(Var(V|N))`
` -I(U;M)-I(V;N)-I(A;B)`.

Since each parent marginal is `nu`, their sum is at most `2 T_t(nu)`.
Taking the supremum over the child channels and then the Bellman pair
proves

`B T_t<=T_t`.                                              (6)

For a finite source, the envelope can be optimized over finite labels:
its posterior-simplex formulation has a continuous reward on a compact
simplex, and a finite posterior mixture suffices by Caratheodory's
theorem. Thus every label above, including `B`, can be finite. The
finite-source statement covers every state reachable from the ternary
root and is sufficient for the cap consequence below. An extension to
arbitrary unbounded source laws is not needed or claimed here.

## 3. Terminal and deep-limit comparisons

The trivial label gives `G<=T_t`, where `G(nu)=g_t(m_2(nu))`, and Jensen
gives `E_t<=T_t`. Also `T_t<=Phi_t`. Indeed, conditionally self-couple
given a latent label to obtain

`F_t(nu)<=I(X;L)+E_L F_t(nu_|L)`.

This follows by bounding the self-coupling information by
`I(X;L)+I(X;Y|L)`; the conditional costs add. Thus

`Phi_t(nu)>=E_L Phi_t(nu_|L)-I(X;L)/2`
` >=E_L g_t(Var(X|L))-I(X;L)`.

Taking the supremum proves the comparison. The already proved
Gaussian-boundary replacement theorem and (6) now imply

`lim_r B^r Phi_t(nu)=lim_r B^r G(nu)<=T_t(nu)`.             (7)

Consequently a rigorously negative value of
`p log2+T_t(nu_p)+t(1-sqrt(p))` is sufficient for some finite-depth
recursive-weave upper certificate. An independent envelope computation
is still required to establish that negative value; this proof alone
does not assert a numerical minimax bound or convergence of the original
normalized minimax sequence.
