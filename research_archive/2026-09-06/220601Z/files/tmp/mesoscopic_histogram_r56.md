# Wave 56 memo: mesoscopic full-histogram phase diagram

Status: analytic derivation independently checked against the exact finite
factorial formulas and the log-domain checker
`tmp/mesoscopic_histogram_r56_check.py`.  Fix `p=m/n` in a compact subset of
`(1/2,1)`, put `alpha=ell/n`, and assume

```math
H=o(\ell),\qquad \ell=o(n).
```

Write `a=(1-p)/p`, `N=binom(n,m)`, and use the common-core kernel and
normalizations of (10.1266)--(10.1270), (10.1292).

## 1. Uniform kernel asymptotics and extraction cost

Factorial cancellation gives exactly

```math
\rho_\ell={\binom m\ell\over\binom n\ell},\qquad
\lambda _1={\ell(n-m)\over m(n-\ell)},
```

```math
\lambda _2={\ell(\ell-1)(n-m)(n-m-1)
 \over m(m-1)(n-\ell)(n-\ell-1)}.
```

Uniformly for `alpha=o(1)` (and with harmless integer-rounding errors),

```math
-\log\rho_\ell
=n\{H(\alpha)-pH(\alpha/p)\}+O(\log n)
=\ell\log(1/p)+{1-p\over2p}{\ell^2\over n}
 +O\!\left({\ell\over n}+{\ell^3\over n^2}+\log n\right),
```

```math
\lambda _1={a\alpha\over1-\alpha},\qquad
\lambda _2=\lambda _1^2\{1-\ell^{-1}+O(n^{-1})\}
=a^2{\ell^2\over n^2}\{1+O(\alpha+\ell^{-1})\}.
```

There is also the useful exact extraction-denominator identity

```math
D_\ell=(1-\lambda _1)+n(\lambda _1-\lambda _2)
=(1-\lambda _1)
\left\{1+{\ell(n-m)(n-2)\over(m-1)(n-\ell-1)}\right\}.
```

Consequently `D_ell=(1+o(1))a ell`.  The denominator costs only a
polynomial and does not change an `exp{-O(H)}` excess.  On the other hand,
`-log rho_ell=omega(H)`, so the direct term in (10.1267) is below the target
scale.  This is a genuinely spectral regime, unlike `ell=O(H)`.

The diagonal is negligible:

```math
h_\ell={1\over\binom{n-\ell}{m-\ell}}
       ={1\over N\rho_\ell}=e^{-\Theta(n)}=o(\lambda _2).
```

## 2. Exact transition law

For a fixed base selector `S`, let

```math
A_j=\binom mj\binom{n-m}{m-j}.
```

Then `A_j/N` is the ordinary intersection law `J_0=|S cap T|`, while

```math
A_jK_\ell(j),\qquad K_\ell(j)={\binom j\ell\over d_\ell},
```

is exactly the probability mass function of

```math
J_\ell-\ell\sim\operatorname{Hypergeom}(n-\ell,m-\ell,m-\ell).
```

Thus

```math
\mu_\ell=\mathbb E J_\ell
=\ell+{(m-\ell)^2\over n-\ell},
```

```math
\operatorname{Var}J_\ell
={(m-\ell)^2(n-m)^2\over(n-\ell)^2(n-\ell-1)}.
```

Relative to the ordinary mean `p^2 n`,

```math
\mu_\ell-p^2n={(1-p)^2\ell\over1-\alpha},\qquad
\sigma_\ell=(1+o(1))p(1-p)\sqrt n.
```

The standardized transition displacement is therefore

```math
\tau_\ell={\mu_\ell-p^2n\over p(1-p)\sqrt n}
=(1+o(1))a{\ell\over\sqrt n}\longrightarrow\infty.
```

## 3. Sharp full-histogram partner capacity

Let `C_ell(r)` be (10.1330), with the self-loop omitted, and let `r_ell^*`
be the least number of highest-intersection partners whose exact total
`K_ell` mass is at least `lambda_2-h_ell` (or `lambda_2-h_ell+epsilon`,
where `epsilon=exp{-O(H)}`; this has the same asymptotics).  If `q_ell`
satisfies `barPhi(q_ell)=lambda_2`, then

```math
q_\ell=2\sqrt{\log(n/\ell)}
 +O\!\left({\log\log(n/\ell)+1\over\sqrt{\log(n/\ell)}}\right)
```

and the boundary shell is

```math
s_\ell^*=\mu_\ell+\sigma_\ell q_\ell+o(\sqrt n/q_\ell).
```

In logarithmic precision, writing

```math
\beta_*=\alpha+{(p-\alpha)^2\over1-\alpha},
```

the shell entropy gives

```math
\log r_\ell^*
=nV(p,\beta_*)
+O\!\left({\ell\over\sqrt n}\sqrt{\log(n/\ell)}
          +\log n\right),
```

```math
\boxed{\log{N\over r_\ell^*}
={a^2\ell^2\over2n}\{1+o(1)\}.}
```

The last formula remains uniform throughout `ell=o(n)` because the cubic
term in the entropy expansion is only `O(ell^3/n^2)=o(ell^2/n)`, and the
polynomial-tail quantile correction is lower order since `ell >> H >>
sqrt(n log n)`.

The sign and size of that correction can also be resolved.  Put

```math
I_*=n\{H(p)-V(p,\beta_*)\},\qquad
\tau=a\ell/\sqrt n.
```

Then the upper transition quantile gives

```math
\log{N\over r_\ell^*}
=I_*+\tau q_\ell
+O(q_\ell^2+\alpha\tau q_\ell+\log n).
```

Thus the complete histogram permits an additional
`exp{tau q_ell+o(tau q_ell)}` thinning beyond the transition-mode Johnson
ball.  This is subleading relative to `I_*` here, but is superpolynomial.

This is a necessary per-base capacity statement.  If a favorable family has
`r=N exp{-L}` members, then every base has off-diagonal load at most
`C_ell(r-1)`.  Hence it cannot reach `lambda_2+epsilon` unless

```math
{a^2\ell^2\over2n}\ge L-o(L).
```

For a migration-scale entropy loss `L=Theta(H)`, the count boundary is

```math
\boxed{\ell_{\rm count}\asymp\sqrt{nH}
=n^{7/8-c/2}.}
```

If `H << ell << sqrt(nH)` and the actual entropy loss is `Theta(H)`, even
placing every available partner in the best possible intersection shells
cannot reach the polynomial baseline `lambda_2`; this is a law-free
feasibility wall at that count.  A mere lower bound `r >= N exp{-C H}` does
not exclude the possibility that the actual family is much denser; rather,
`ell >> sqrt(nH)` is what makes that worst-case guaranteed count sufficient.
This is only count-feasibility, not a signing theorem or a simultaneous
neighborhood construction.

## 4. Exact loss of one-threshold compression

The strongest possible one-threshold rectangle, even using the entire
slice, is

```math
R_\ell:=\max_{\ell\le s\le m-1}
\mathcal B_{n,m}(s){\binom s\ell\over d_\ell}.
```

The exact ratio

```math
{A_{j+1}\over A_j}
={(m-j)^2\over(j+1)(n-2m+j+1)}
```

and a local limit theorem at the `J_ell` mode show that the optimizer is
`s=mu_ell+O(n/ell)`.  There

```math
{\mathcal B_{n,m}(s)\over A_s}
=(1+o(1)){p^2n\over\ell},\qquad
A_sK_\ell(s)=(1+o(1)){1\over\sqrt{2\pi}p(1-p)\sqrt n}.
```

Therefore

```math
\boxed{R_\ell=(1+o(1)){p\over\sqrt{2\pi}(1-p)}
{\sqrt n\over\ell}.}
```

This extends the linear-core `Theta(n^{-1/2})` wall and gives the sharper
comparison

```math
\boxed{{R_\ell\over\lambda _2}
=(1+o(1)){1\over\sqrt{2\pi}a^3}{n^{5/2}\over\ell^3}.}
```

Thus one-threshold compression is capacity-impossible for
`ell >> n^(5/6)`.  At `ell=kappa n^(5/6)`, its exact constant threshold is

```math
\kappa\le(2\pi)^{-1/6}{p\over1-p}.
```

Below `n^(5/6)` this particular capacity wall disappears, but the family
must still contain enough partners.  Combining both requirements gives the
only mesoscopic one-threshold feasibility window when its actual entropy
loss is `Theta(H)`:

```math
\boxed{\sqrt{nH}\ \lesssim\ \ell\ \lesssim\ n^{5/6}.}
```

At exponent level this window is nonempty exactly when `c>=1/12` (with
constants decisive at equality).  For `c<1/12`, no mesoscopic
one-threshold choice can simultaneously overcome the saved-degree entropy
loss and retain the `lambda_2` baseline.  The full histogram has no
`n^(5/6)` upper wall: it sums the transition mass that every single
rectangle discards.

There is a sharper comparison below the upper wall.  If

```math
U_\ell={R_\ell\over\lambda _2}
\sim{1\over\sqrt{2\pi}a^3}{n^{5/2}\over\ell^3}\longrightarrow\infty,
```

then the highest threshold which can still meet the baseline is

```math
s_{\rm rect}=\mu_\ell+\sigma_\ell u_{\rm rect}+o(\sqrt n),
\qquad
u_{\rm rect}=\sqrt{2\log U_\ell}+o(\sqrt{\log U_\ell}).
```

Its necessary partner count obeys

```math
\log{N\over r_{\rm rect}}
=I_*+\tau u_{\rm rect}+O(u_{\rm rect}^2+\log n),
```

whereas the complete histogram has `u=q_ell`.  Since
`q_ell>u_rect` for every power scale between `sqrt(n)` and `n^(5/6)`, the
full histogram tolerates fewer partners by the factor

```math
\exp\{\tau(q_\ell-u_{\rm rect})+o(\tau\sqrt{\log n})\}.
```

This is a genuine, superpolynomial advantage even in the range where the
one-threshold route is not absolutely impossible; it is still subexponential
on the `H` scale at `ell asymp sqrt(nH)`.

More geometrically, a threshold rectangle keeps only one boundary atom
times an ordinary-tail Mills factor.  Near the transition this is
`Theta(sqrt(n)/ell)`, whereas the exact histogram integrates all relevant
transition levels.  To collect merely `lambda_2`, its upper-tail window has
effective width `Theta(sqrt(n)/(q_ell))` levels; compressing those weights to
one boundary value is precisely the avoidable loss.

For a power scale `ell=n^gamma` and entropy loss `L=Theta(H)`, with
`h=3/4-c`, the complete phase diagram is therefore:

- `gamma<=h`: this is the direct-core range already identified in
  (10.1270), not the present mesoscopic regime;
- `h<gamma<(1+h)/2`: direct extraction is too small and a
  `Theta(H)`-loss family has insufficient full-histogram capacity;
- `gamma>(1+h)/2=7/8-c/2`: the full histogram is count-feasible;
- one-threshold compression is additionally feasible only when
  `gamma<5/6` (subject to the signing-specific cluster theorem).

Thus if `c>1/12`, the interval
`7/8-c/2<gamma<5/6` is a genuinely open one-threshold mesoscopic window.
If `c<1/12`, the count boundary lies above `5/6`: no power scale meets both
one-threshold requirements, although the full histogram becomes
count-feasible above `7/8-c/2`.  At either equality, constants and the
displayed logarithmic corrections decide.

## 5. Research judgment

The first genuinely plausible mesoscopic full-histogram scale is
`ell~sqrt(nH)`, not merely `ell/H -> infinity`.  Existing migration
(10.1323)--(10.1326) produces a near-parent witness for each chosen excess
block but proves neither `N exp{-O(H)}` distinct child selectors nor bounded
witness multiplicity.  Therefore it does not yet meet the capacity theorem.
An aggregate migration theorem that really yields saved selector degree
could plausibly feed the full histogram at or above `sqrt(nH)`.

For `c>1/12`, a narrow one-threshold experiment between `sqrt(nH)` and
`n^(5/6)` is not ruled out and is genuinely different from the fixed-linear
failure.  For `c<1/12`, one-threshold compression is globally infeasible at
every mesoscopic scale; only the complete histogram can remain active.  In
all cases, count-feasibility alone does not supply the signing-specific
correlation needed to make a positive fraction of biased bases realize the
optimal shell placement.
