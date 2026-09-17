# Wave 51A: exact completion balances and collapse of the two-moment support

## Outcome and scope

- **Verified:** for every selector, the local energy, conditional completion
  variance, Parseval slack, and the raw reverse-branch numerator obey exact
  averages.  In particular both the mean margin and the mean raw reverse
  pressure have the wrong sign when the principal excess is positive.
- **Verified:** on any compact density window inside `(1/2,1)`, every state
  on which the Wave 50 two-moment envelope is positive has
  `g >= -O(n)`.  The apparently different near-Parseval branch therefore
  cannot obtain mass from a genuinely macroscopic negative margin.
- **Verified:** if the averaged two-moment envelope has the saved mass
  `exp{-O(n^(3/4-c))}`, then a direct conditional Hanson--Wright inversion
  already gives one principal submatrix with
  `Q(A[S]) <= p^(3/2)q_n+O(n^(3/2-c))`.  Hence the scalar recurrence, and
  convergence through the established exact-landing chain, hold before the
  annealed row extraction is used.
- **Verified scope audit:** the same conclusion applies to the stronger
  degree-two Bonami event `g >= -alpha p_2 sqrt(V)` for every fixed `alpha`.
  Its extra negative-margin width is only `O(n^(5/4))`.
- **Not proved:** none of these results forces saved mass for the envelope,
  nor do they bound the full completion incidence `Z_t` from below.  Higher
  conditional information can still make the true CDF large well outside the
  two-moment support.  Thus the result retires the two-moment envelope as an
  independently weaker mechanism; it does not falsify the full annealed route.

Throughout, write

```math
p=\frac mn,\qquad
p_2=\frac{m(m-1)}{n(n-1)},\qquad a=1-p_2,
```

and keep `p` in a fixed compact interval `[p_0,p_1]` contained in `(1/2,1)`.
No replacement of `p_2` by `p^2` is made in an identity below.  Let `A` be an
exact order-`n` minimizer, `q=Q(A)=q_n`, let `S` have size `m`, and put

```math
Q=Q(A[S]),\qquad h=Q-p^{3/2}q-t,
\qquad g=ae-h.
```

Here `t >= 0`, as in the restriction tolerance, and the local oriented state
is `a_0=(sigma,y)` with
`e=sigma y^T A[S]y`.  To avoid collision with `a=1-p_2`, the state itself is
denoted `a_0` in this memo.  Put `T=S^c`, `k=n-m`, and retain

```math
V=4\lVert A[T,S]y\rVert_2^2+2k(k-1),
\qquad D=q^2-e^2-V.
```

## 1. Exact fixed-selector balance laws

For a uniform oriented projective local state, orientation symmetry and
Walsh orthogonality give

```math
\boxed{
\begin{aligned}
\mathbb E_{a_0}e&=0,\\
\mathbb E_{a_0}e^2&=2m(m-1),\\
\mathbb E_{a_0}V&=4mk+2k(k-1)=2k(n+m-1),\\
\mathbb E_{a_0}D&=q^2-2n(n-1).
\end{aligned}}
\tag{R51A.1}
```

Indeed `E_y yy^T=I`, so
`E_y ||A[T,S]y||^2=||A[T,S]||_F^2=mk`, while the local quadratic energy has
second moment twice the number of ordered off-diagonal entries.  The last
line uses

```math
m(m-1)+k(n+m-1)=n(n-1).
```

There is a second exact cancellation which is important for the proposed
near-Parseval pressure:

```math
p_2\,\mathbb E_{a_0}V
=(1-p_2)\,\mathbb E_{a_0}e^2.
\tag{R51A.2}
```

This follows from

```math
1-p_2=\frac{k(n+m-1)}{n(n-1)}.
```

The untruncated numerator in the reverse formula (10.1245) is

```math
N=(q-e)(p_2q+e-h)-p_2D
  =p_2V+(q-e)g.
\tag{R51A.3}
```

Using (R51A.1)--(R51A.2) gives the **exact fixed-selector identities**

```math
\boxed{
\mathbb E_{a_0}g=-h,
\qquad
\mathbb E_{a_0}N=-qh.}
\tag{R51A.4}
```

Thus when `h>0`, neither the margin nor the untruncated reverse numerator has
positive average pressure.  Positive parts can exist, as the finite examples
show, but no first-moment averaging of these quantities forces them.
Also, (R51A.1) says that the mean conditional Parseval slack is the same for
every selector and is of order `q^2`: it carries no favorable selector signal.

## 2. Positive reverse pressure forces a near-zero margin

First note the elementary but crucial principal cap

```math
\boxed{Q(A[S])\le q.}\tag{R51A.5}
```

For any local spin, extend it by independent uniform outside spins.  The mean
full energy equals the local energy, while every full energy lies in
`[-q,q]`; maximize over local spins.

Let

```math
\beta_0=\min_{p\in[p_0,p_1]}(p^{3/2}-p^2)>0.
```

Since `p_2 <= p^2`, the buffer

```math
B=(p^{3/2}-p_2)q
```

obeys `B >= beta_0 q` exactly.  Consider a state in the reverse branch, so
`g<=0`, and put `r=-g/p_2`.  From (R51A.5),

```math
ae\le h\le(1-p^{3/2})q-t,
```

and hence

```math
\boxed{
q-e\ge\frac{(p^{3/2}-p_2)q+t}{1-p_2}
=\frac{B+t}{a}\ge\beta_0q.}
\tag{R51A.6}
```

Strict positivity of the reverse bound is exactly

```math
V-r(q-e)>0.
```

The exact-minimizer spectral estimate and `||y||_2^2=m` give

```math
\boxed{
V\le4\lVert A\rVert_{op}^2m+2n^2
\le8qm+2n^2.}
\tag{R51A.7}
```

Using the established uniform exact-minimizer bounds
`q=Theta(n^(3/2))`, (R51A.6)--(R51A.7) imply

```math
0\le r<\frac{V}{q-e}
\le\frac{8qn+2n^2}{\beta_0q}=O(n).
\tag{R51A.8}
```

Therefore `g=-p_2r >= -O(n)` on every positive reverse state.  The positive
Cantelli branch already has `g>0`.  If `C_t(S,a_0)` denotes the full Wave 50
two-moment envelope, we have proved the uniform support theorem

```math
\boxed{
C_t(S,a_0)>0\quad\Longrightarrow\quad g\ge-C_0n,}
\tag{R51A.9}
```

where `C_0` depends only on the fixed density window and the universal
two-sided constants for `q_n/n^(3/2)`.

This is stronger than calling the reverse branch "near Parseval": at fixed
density, its positive numerator forces its negative threshold to be only
linear in `n`, even though `q` is of order `n^(3/2)`.

## 3. Saved two-moment mass already implies the principal recurrence

Fix `0<c<1/4` and set

```math
L_n=n^{3/4-c},\qquad 0\le t=O(n^{3/2-c}).
```

Suppose the proposed Wave 50 sufficient estimate holds:

```math
\mathcal H_t(A,m)=\mathbb E_{S,a_0}C_t(S,a_0)
\ge e^{-C L_n}.
\tag{R51A.10}
```

Since `0<=C_t<=1`, (R51A.9) implies

```math
\Pr_{S,a_0}\{g\ge-C_0n\}\ge e^{-CL_n}.
```

There is consequently one selector `S` for which the same local event has
mass at least `e^{-CL_n}`.  On that selector,

```math
e\ge u:=\frac{h-C_0n}{1-p_2}.
\tag{R51A.11}
```

If `u<=0`, then `h=O(n)` and the desired conclusion is immediate.  If
`u>0`, orientation symmetry and conditional Hanson--Wright for the local
zero-diagonal signing give

```math
\begin{aligned}
e^{-CL_n}
&\le\Pr_{a_0}\{e\ge u\}\\
&=\frac12\Pr_y\{|y^TA[S]y|\ge u\}\\
&\le C_1\exp\left[-c_1\min\left\{
\frac{u^2}{n^2},\frac{u}{n^{3/4}}
\right\}\right].
\end{aligned}
\tag{R51A.12}
```

Here `||A[S]||_F^2=m(m-1)<=n^2` and
`||A[S]||_op<=||A||_op<=sqrt(2q)=O(n^(3/4))`; these are the only matrix
bounds used.  Inverting (R51A.12),

```math
u=O\!\left( n\sqrt{L_n}+n^{3/4}L_n\right)
=O(n^{3/2-c}),
\tag{R51A.13}
```

because `c<1/4` makes
`n sqrt(L_n)=n^(11/8-c/2)<=n^(3/2-c)`.  Equations
(R51A.11)--(R51A.13) yield

```math
h=Q(A[S])-p^{3/2}q-t=O(n^{3/2-c}).
```

Finally `A[S]` is an order-`m` signing, so

```math
\boxed{
q_m\le Q(A[S])
\le p^{3/2}q_n+O(n^{3/2-c}).}
\tag{R51A.14}
```

This is the scalar restriction recurrence (10.1024), which already proves
convergence with the established ratio-window uniformity and exact landing.
Thus a saved lower bound for `H_t` cannot be powered by a broad population of
high-`h`, negative-margin states.  It is a stronger certificate than the
principal recurrence it was intended to prove, not an independent route
around that recurrence.

## 4. Audit of the degree-two Bonami enlargement

For completeness, a centered Rademacher polynomial `Z` of degree at most two
and variance `V` satisfies

```math
\lVert Z\rVert_4\le3\lVert Z\rVert_2,
\qquad
\mathbb E|Z|\ge
\frac{(\mathbb EZ^2)^{3/2}}{(\mathbb EZ^4)^{1/2}}
\ge\frac{\sqrt V}{9}.
```

The first inequality is Bonami hypercontractivity and the second follows by
`L^1`--`L^2`--`L^4` interpolation.  Centering gives
`E Z_-=E|Z|/2>=sqrt(V)/18`.  Splitting at `r>=0` and using Cauchy--Schwarz,

```math
\mathbb EZ_-
\le r+\sqrt{V\Pr\{Z<-r\}},
```

so

```math
\boxed{
\Pr\{Z\le-r\}\ge
\left(\frac1{18}-\frac r{\sqrt V}\right)_+^2.}
\tag{R51A.15}
```

(When `V=0`, interpret only the deterministic threshold directly.)  Hence
`g>=0` contributes at least `1/324`, and

```math
g\ge-\frac{p_2\sqrt V}{36}
```

contributes at least `1/1296`.  This is a valid improvement over the Wave 50
moment envelope.  However, (R51A.7) gives

```math
\sqrt V=O(n^{5/4}).
```

Thus saved mass of the Bonami event implies saved mass of
`g>=-O(n^(5/4))`.  Repeating (R51A.11)--(R51A.13), with `n^(5/4)` in place
of `n`, again proves (R51A.14), because
`n^(5/4)=o(n^(3/2-c))` for every `c<1/4`.  The enlarged event is therefore
also principal-recurrence-strong at the target entropy scale.

## 5. Falsification criteria and surviving target

The finite checker exhausts all selectors and local states in the stored
`A_6`, `A_8`, `A_9`, and one deterministic exact `A_10` minimizer.  It
reproduces (R51A.1)--(R51A.4), verifies (R51A.6) for every positive reverse
state, and verifies the support bound algebra.  These computations are an
audit of exact formulas, not asymptotic evidence.

The conclusions have precise scope:

1. Any proposed proof that forces `H_t` by averaging `g`, `D`, or the raw
   numerator `N` must overcome the exact wrong-sign balances (R51A.4), not
   merely reweight their first moments.
2. A putative family with positive reverse pressure at
   `-g=omega(n)` would falsify (R51A.9), but cannot exist under the stated
   fixed-density, exact-minimizer, and nonnegative-tolerance assumptions by
   (R51A.6)--(R51A.8).
3. An unbounded exact-minimizer family with
   `-log H_t=omega(n^(3/4-c))` would still falsify only the two-moment
   implementation.
4. The full annealed target remains
   `Z_t>=exp{-O(n^(3/4-c))}`.  It may use the detailed linear-plus-quadratic
   completion law when `g` is substantially negative; neither the support
   collapse nor Hanson--Wright above rules that out.
