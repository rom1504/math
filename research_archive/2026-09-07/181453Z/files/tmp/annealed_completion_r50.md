# Wave 50A: a sharp two-moment lower envelope for the annealed completion tail

## Status

- **Verified:** the completion CDF in (10.1232) has an explicit two-branch
  lower bound using only the conditional mean, variance, and the parent cap.
  The positive-margin branch is one-sided Cantelli.  The negative-margin
  branch is a support-and-second-moment inequality, equivalently a
  near-saturation condition for conditional Parseval.
- **Verified:** averaging this envelope gives a concrete sufficient scalar
  lemma for the Wave 49 annealed target.  It separates two genuinely
  different populations: moderately high local energy and unusually large
  completion variance.
- **Verified scope audit:** the positive-margin population contains the old
  local-deficit population and has a strictly wider numerical cutoff whenever
  $q-Q+B+t>0$.  It is not algebraically the principal recurrence, but it is
  still a local-energy subroute and cannot represent the full completion
  mechanism.
- **Numerical:** exhaustive calculations inside `A6`, `A8`, `A9`, and one
  deterministic exact order-ten minimizer find a positive two-branch bound in
  every tested pair.  At high ratios most actual favorable incidences instead
  come from the negative-margin population, while the two-moment reverse
  branch captures only a small fraction of them.
- **Open:** no asymptotic lower bound on the averaged two-moment envelope, or
  on the exact annealed mass, is proved.

## 1. Conditional completion coordinates

Fix a selector $S$, put $T=S^c$, $k=|T|$, and abbreviate

```math
q=Q(A),\qquad Q=Q(A[S]),\qquad
p_2=\frac{m(m-1)}{n(n-1)},\qquad
Y=Q-p^{3/2}q.
```

Let $a=(\sigma,y)$ be an oriented local projective state and write

```math
e=\sigma y^{\mathsf T}A[S]y.
```

For a uniform outside completion $w$, its full oriented energy can be written

```math
E_{S,a}(w)=e+Z_{S,a}(w),\qquad \mathbb E_w Z_{S,a}=0.
```

The linear and quadratic Walsh levels in the outside spins are orthogonal, so
the conditional variance is the exact quantity

```math
\boxed{
V_{S,a}:=\mathbb E_wZ_{S,a}^2
=4\lVert A[T,S]y\rVert_2^2+2k(k-1).}
\tag{R50A.1}
```

The orientation does not affect $V$.  Since every full oriented energy lies
in $[-q,q]$, conditional Parseval gives

```math
\boxed{e^2+V_{S,a}\le q^2.}
\tag{R50A.2}
```

Define the conditional Parseval slack

```math
D_{S,a}=q^2-e^2-V_{S,a}\ge0.
\tag{R50A.3}
```

For tolerance $t$, put

```math
h=Y-t,\qquad g=(1-p_2)e-h.
\tag{R50A.4}
```

Directly rearranging the favorable inequality, or (10.1232), gives the exact
normal form

```math
\boxed{
F_{S,a}\!\left(q+\frac{B+t-\delta_S(a)}{p_2}\right)
=\Pr_w\left\{Z_{S,a}\le\frac{g}{p_2}\right\}.}
\tag{R50A.5}
```

Thus $g$ is precisely the threshold measured from the conditional mean.  It
is important that (R50A.5) still uses all outside completions; no complement
cut or robust cylinder has been imposed.

## 2. The two sharp moment branches

When $g>0$, one-sided Cantelli applied to the mean-zero variable $Z=Z_{S,a}$
gives

```math
\boxed{
\Pr\left\{Z\le\frac g{p_2}\right\}
\ge C_+(e,V):=
\frac{g^2}{g^2+p_2^2V}.}
\tag{R50A.6}
```

Using the strict event $Z<g/p_2$ makes (R50A.6) insensitive to atoms at the
threshold.  Cantelli is sharp given only the first two moments.

There is also a rigorous branch when $g\le0$.  Set

```math
r=-\frac g{p_2}\ge0,\qquad
a=q+e,\qquad b=q-e,
```

so $Z\in[-a,b]$ and the favorable event contains $\{Z<-r\}$.  Suppose
$0\le r<a$.  The quadratic

```math
\phi(z)=(z+r)(z-b)
```

is nonpositive on $[-r,b]$ and lies between zero and
$(a-r)(a+b)$ on $[-a,-r)$.  Since

```math
\mathbb E\phi(Z)=V-rb,
```

we obtain the **Verified reverse bound**

```math
\boxed{
\Pr\{Z\le-r\}\ge\Pr\{Z<-r\}
\ge C_-(e,V):=
\frac{(V-rb)_+}{(a-r)(a+b)}.}
\tag{R50A.7}
```

This is the sharp lower bound from support, mean, and variance, apart from the
irrelevant choice of whether an atom exactly at $-r$ is counted.  Indeed a
three-point law on $\{-a,-r,b\}$ attains its strict-tail value whenever
$rb\le V\le ab$.

The reverse branch is more transparent in Parseval-slack coordinates.  With
$D=D_{S,a}$, direct cancellation gives

```math
\boxed{
C_-(e,V)=
\frac{\big[(q-e)(p_2q+e-h)-p_2D\big]_+}
{2q(p_2q+e-h)},}
\tag{R50A.8}
```

under

```math
(1-p_2)e\le h<p_2q+e.
```

The right inequality is exactly $r<q+e$.  Therefore this is not a generic
``large variance helps'' assertion: it asks the conditional second moment to
come sufficiently close to its cap.  At exact saturation $D=0$, the bound is
$(q-e)/(2q)$.

Define the two-moment envelope

```math
\mathcal C_t(S,a)=
\begin{cases}
C_+(e,V),&g>0,\\
C_-(e,V),&g\le0\text{ and }r=-g/p_2\in[0,q+e),\\
0,&\text{otherwise}.
\end{cases}
\tag{R50A.9}
```

The interval condition is empty when $q+e=0$; hence every zero-denominator
case, as well as every case with $r\ge q+e$, is assigned the zero branch.
Equations (R50A.6)--(R50A.7) prove the averaged **Verified bound**

```math
\boxed{
Z_t\ge \mathcal H_t(A,m):=
\mathbb E_{S\sim U_m,\ a}\mathcal C_t(S,a).}
\tag{R50A.10}
```

For a cleaner final statement, (R50A.9)'s middle condition can equivalently
be written

```math
g\le0,\qquad h<p_2q+e.
```

## 3. Exact asymptotic sufficient statements

Fix a compact density window and $0<c<1/4$, and put
$L_n=n^{3/4-c}$.  The following is a direct sufficient lemma for the annealed
route:

```math
\boxed{
t=O(n^{3/2-c}),\qquad
\mathcal H_t(A,m)\ge \exp\{-O(L_n)\}.}
\tag{R50A.11}
```

Indeed (R50A.10) gives $-\log Z_t=O(L_n)$, and
(10.1228)--(10.1230) then prove the bare target (10.795).

Two population versions make clear what would establish (R50A.11).

1. If a set $\mathcal E_+$ of $(S,a)$ has mass $\mu_+$ and obeys

   ```math
   g\ge s>0,\qquad V\le Ks^2,
   ```

   then

   ```math
   \mathcal H_t\ge\frac{\mu_+}{1+p_2^2K}.
   \tag{R50A.12}
   ```

   Hence
   $-\log\mu_++\log(1+K)=O(L_n)$ suffices.  In particular, because
   $V\le q^2=O(n^3)$, saved mass of states with any inverse-polynomial
   positive margin is enough up to a polynomial loss.

2. If a set $\mathcal E_-$ has mass $\mu_-$ and obeys the reverse-branch
   conditions together with $C_-(e,V)\ge\eta$, then

   ```math
   \mathcal H_t\ge\mu_-\eta.
   \tag{R50A.13}
   ```

   Thus $-\log\mu_- -\log\eta=O(L_n)$ suffices.  By (R50A.8), this asks
   for saved mass of local states with a quantitatively small conditional
   Parseval slack relative to their threshold geometry.

These clauses are independently testable.  Failure of both clauses does not
falsify the exact CDF target, because higher conditional moments may still
put mass below a negative threshold.  An unbounded exact-minimizer family
with $-\log\mathcal H_t=\omega(L_n)$ would falsify this two-moment
implementation; only $-\log Z_t=\omega(L_n)$ would falsify the annealed
implementation itself.

## 4. Relation to the retired local-deficit condition

Write the local deficit as $\delta=Q-e$ and retain
$B=(p^{3/2}-p_2)q$.  Formula (R50A.4) becomes

```math
\boxed{
g=B+t+p_2(q-Q)-(1-p_2)\delta.}
\tag{R50A.14}
```

The old all-completions condition is $\delta\le B+t$.  It implies

```math
g\ge p_2(q-Q+B+t)>0
```

when the displayed quantity is positive.  Conversely, the new
positive-margin cutoff is

```math
\boxed{
\delta<\frac{B+t+p_2(q-Q)}{1-p_2}.}
\tag{R50A.15}
```

Its excess over the old cutoff is

```math
\frac{p_2(q-Q+B+t)}{1-p_2}.
```

Thus $g>0$ is a strict relaxation of the old local-deficit population, not a
renaming of it.  It is also not algebraically the principal recurrence.  If
$Y>t$, it asks for the high local-energy layer

```math
e>\frac{Y-t}{1-p_2}
```

to have saved entropy; it does not upper-bound $Q$.  However, it remains a
local-energy-only subroute.  A proof by a robust coordinate cylinder would
re-enter the circular wall (10.1233), and the negative-margin branch shows
that discarding completion deviations would lose real incidences.

## 5. Finite falsification audit

`tmp/annealed_completion_r50_check.py` exhausts every selector, oriented
local state, and outside completion in `A6`, `A8`, `A9`, and one
deterministically generated exact order-ten minimizer.  It checks (R50A.1),
(R50A.2), both forms (R50A.7)--(R50A.8), and the pointwise lower bound.
All percentages below are **Numerical finite evidence**, not asymptotic
claims.

Across the sixteen tested pairs at $t=0$:

- $\mathcal H_0/Z_0$ ranges from $0.06968$ to $0.63887$;
- the two-moment bound is strictly positive in every case;
- at `A9,m=8`, $Z_0=0.2526042$, the positive branch is $0.0131798$,
  the combined bound is $0.0176024$, and the combined bound is only $6.97\%$
  of the exact mass;
- in that same case, $88.66\%$ of exact favorable incidence comes from
  $g\le0$, while the $g>0$ local-state mass is $0.045139$ versus only
  $0.012153$ for the old local-deficit population;
- at the sampled `A10,m=8`, $73.66\%$ of exact incidence comes from
  $g\le0$, and the combined bound captures $10.84\%$ of $Z_0$.

Therefore the finite search does not falsify (R50A.11), but it does falsify
the heuristic that the positive-margin/Cantelli population explains nearly
all of (10.1232).  The reverse moment branch is real but quantitatively modest
at the hardest audited high-ratio cases.  An asymptotic proof must either
establish saved mass in one of (R50A.12)--(R50A.13), or use information beyond
the first two conditional moments.
