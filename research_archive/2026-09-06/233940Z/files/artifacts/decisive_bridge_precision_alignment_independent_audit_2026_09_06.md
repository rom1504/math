# Independent audit of the director's precision-side alignment proof

Date: 2026-09-06. Status: PASS. The proof was supplied by the director
after the classwise-envelope equality investigation; the argument below
is my independent reconstruction and adversarial check. It proves the
previously open smaller-envelope supersolution, not original minimax
convergence.

## 1. Definitions and scalar sign check

```math
J_\lambda(X)=\inf_L\{I(X;L)+\lambda\mathbb E\operatorname{Var}(X\mid L)\},
\quad
E_t(X)=\sup_{0<\lambda\le t}\{c_t(\lambda)-J_\lambda(X)\},
```

where `c_t(lambda)=log(lambda(2t-lambda)/t²)/4`.
The two optimizations are both suprema after negating J, so this is
exactly `sup_L[g_t(E Var(X|L))-I(X;L)]`, not the classwise envelope T.

For lambda in (0,t],

```math
{d\over d\log\lambda}c_t(\lambda)
={t-\lambda\over2(2t-\lambda)},\qquad
{d^2\over d(\log\lambda)^2}c_t(\lambda)
=-{t\lambda\over2(2t-\lambda)^2}<0.
```

Thus c_t is strictly concave in log precision. This is the relevant
curvature on the PRECISION side, distinct from the previously used
concavity of g_t in log variance.

## 2. Weighted reconstruction and information

Take any finite pair (A,B) and set U=(A+B)/sqrt(2), V=(A-B)/sqrt(2).
Fix positive child precisions lambda_1,lambda_2<=t and arbitrary finite
child channels M|U,N|V independently given (U,V). Put
`u_hat=E[U|M]`, `v_hat=E[V|N]`,
`a_hat=(u_hat+v_hat)/sqrt(2)`, `b_hat=(u_hat-v_hat)/sqrt(2)`.

Define

```math
a={\lambda_1+\lambda_2\over2},\quad
h={2\lambda_1\lambda_2\over\lambda_1+\lambda_2},\quad
d={\lambda_1-\lambda_2\over\lambda_1+\lambda_2}.
```

The exact pointwise square identity is

```math
\lambda_1(U-\widehat u)^2+\lambda_2(V-\widehat v)^2
=a[A-\widehat a+d(B-\widehat b)]^2+h(B-\widehat b)^2.
```

The A reconstruction `a_hat-d(B-b_hat)` is measurable from (B,M,N).
The B reconstruction b_hat is measurable from (M,N). Replacing each
by the corresponding conditional mean lowers its own positive-weighted
mean squared error. Hence the weighted sum of the parent conditional
MSEs is at most the displayed child weighted MSE.

With L=(M,N), the exact information identity is

```math
I(A;B,L)+I(B;L)
=I(A;B)+I(U;M)+I(V;N)-I(M;N).
```

Thus the parent feasible objectives yield

```math
J_a(A)+J_h(B)
\le J_{\lambda_1}(U)+J_{\lambda_2}(V)+I(A;B).
```

For the last display take arbitrarily accurate child optimizers; no
claim of simultaneous or special-structure optimizers is required.
No posterior Gaussianity, equal distortion, equal temperature, or
independence of the marginalized labels M,N is assumed.

## 3. The smaller envelope really is a supersolution

The two parent precisions a,h lie between lambda_1 and lambda_2, have
the same product as those two values, and belong to (0,t]. Their logs
are therefore less spread with the same sum. Section 1 gives

```math
c_t(a)+c_t(h)\ge c_t(\lambda_1)+c_t(\lambda_2).
```

Combining with Section 2 proves the heterogeneous-parent inequality

```math
E_t(U)+E_t(V)-I(A;B)\le E_t(A)+E_t(B).
```

In the original Bellman problem, the safe simultaneous-reversal and
input-interchange averaging preserves each child absolute law while
decreasing relative entropy. It makes BOTH actual signed input marginals
equal to the symmetric source nu. Its cost is then I(A;B). Dividing the
last inequality by two and optimizing the pair proves

```math
\mathcal B E_t\le E_t.
```

The former apparent obstacle was requiring one common parent precision
inside a single preselected channel. That is unnecessary: the two parent
sources have the same law, and each can use a different admissible branch
when bounded by E_t(nu).

Precision endpoints cause no issue. Every selected lambda is positive,
so both a,h are positive even if source covariances are singular. If an
outer or channel optimum is not attained, epsilon-optimal positive
precisions/channels suffice. The lambda=t endpoint is admissible.

## 4. Consequences with audited dependencies

The archived conditional-copy theorem gives `B E_t>=E_t`, so E_t is now
a genuine fixed point. Combining the already audited Gaussian-boundary
replacement with `E_t<=H_t<=T_t` gives

```math
\lim_r\mathcal B^r\Phi_t=H_t=E_t.
```

The independently proved classwise equality criterion remains consistent:
`BT=T iff E=T iff H=T`. The explicit sources with E<T therefore show
that the larger classwise supersolution was genuinely loose, while the
smaller average-variance envelope now exactly describes the binary
Gaussian-boundary Bellman certificate limit. The actual one-row annealed
pressure is only bounded above by this certificate; no equality is claimed.

I independently replayed the old exact E-root certificate, obtaining the
identical rational `a_E=19678127864847/800000000000000` at p=31/32,t=4.
The replay is in
`computations/decisive_bridge_supersolution_ternary_E_replay_2026_09_06.json`.
Under the newly established alignment theorem this certificate can now
be passed through the existing stopped-tree and weave chain. The same
all-order cap formula becomes `1/2-a_E/(8sqrt(31/32))`. The director owns
the final cap statement and its independent all-order verification.

Neither this characterization of the particular recursive ensemble nor
its improved original-signing upper bound identifies actual minimizing
signings with the ensemble. The original convergence question remains
separate.
