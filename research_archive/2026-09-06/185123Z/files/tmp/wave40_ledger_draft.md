### 10.93 Wave 40: flat complement exchange, curvature-defect walls, and exact favorable soft trees

Wave 40 begins after the first blank-slate abstraction audit.  It attacks the
anchored-conflict clause of (10.1089), the scalar harmonic defect (10.1096),
and the terminal soft-potential fallback.  Complement ground maps are exactly
cyclically monotone and admit a four-corner exchange cone, but exact finite
minimizers show that the cone has no Hamming modulus.  Harmonic curvature has
an exact third-derivative and two-replica calculus; literal scalar monotonicity
fails, and an abstract fixed-size erasure mixture shows that endpoint cost,
parent entropy, and the matched exclusion identities alone cannot bound the
isotonic defect.  Favorable-fiber distance is an exact flip-compatible soft
potential, but its dynamic program is precisely the old hard terminal forest
and center/scaffold problem.

All displayed algebraic identities, finite exhaustive claims, and checker
outputs below are **Verified**, including independent audits.  Asymptotic
normal-fan curvature, signing-specific harmonic control, and the project-scale
soft-scaffold event are **Open**.  Finite quadrature values are **Numerical**.
Nothing in this wave proves convergence or gives an asymptotic falsifier.

#### 10.93.1 Complement ground maps are cyclically monotone but not metrically coherent

For an `m`-selector `S`, let `B_S=A^{F_S}` flip every edge outside `E(S)`.
Thus, for an oriented cut `d`,

```math
H_{B_S}(d)=2c_S(d)-E_d.
```

Let `P(d)` be any selector-independent penalty and choose

```math
d_S\in\mathop{\rm argmax}_d\{H_{B_S}(d)-P(d)\}.
```

Pairwise optimality makes both the penalty and parent energy cancel:

```math
\boxed{
\langle B_S-B_T,d_S-d_T\rangle\ge0,
\qquad
c_S(d_S-d_T)-c_T(d_S-d_T)\ge0.
}
\tag{10.1103}
```

For a cycle `S_1,...,S_k,S_{k+1}=S_1`, the stronger cyclic form is

```math
\boxed{
\sum_{j=1}^k c_{S_j}(d_{S_j}-d_{S_{j+1}})\ge0.
}
\tag{10.1104}
```

Optimization errors `zeta_S` and `zeta_T` change the second right side in
(10.1103) to at worst `-(zeta_S+zeta_T)/2`.  Equations
(10.1103)--(10.1104) remain exact after a common row, Hamming-to-center, or
other price.  They organize all such canonical response maps as a normal-fan
subgradient, but supply no strong monotonicity.

The missing modulus is visible in a separate four-corner identity.  Thin by
at most a constant factor to a common response orientation `sigma`, put two
projective words `x,w` in the common anchor gauge, and let `U` be their global
disagreement set.  Put

```math
H=U\cap((S\cap T)\setminus\{v\}),
\qquad R=U\setminus H,
```

and use signed parent shores `C_x(H),C_w(H)`.  Then, exactly,

```math
\boxed{
C_x(H)+C_w(H)=2J,
\quad
J=\sum_{i\in H,\ j\in R}\sigma a_{ij}x_i x_j,
}
\tag{10.1105}
```

and simultaneous crossover of `H` gives

```math
E_A(x^H)+E_A(w^H)=E_A(x)+E_A(w)-8J.
```

If `E_A(x)=q_n-\Delta_x` and `E_A(w)=q_n-\Delta_w`, the parent cap yields only

```math
\boxed{
-\frac{\Delta_x+\Delta_w}{8}
\le J\le
\frac{q_n}{2}-\frac{\Delta_x+\Delta_w}{8}.
}
\tag{10.1106}
```

The difference `P=(C_x(H)-C_w(H))/2` has the analogous centered interval
obtained from the two individual cap inequalities.  Repeating the calculation
in `B_S,B_T`, with complement cap slacks `\Gamma_x,\Gamma_w` and
`\Xi=(C_{B_S,x}(H)+C_{B_T,w}(H))/2`, gives

```math
\boxed{
\Xi\ge-\frac{\Gamma_x+\Gamma_w}{8}.
}
\tag{10.1107}
```

None of (10.1103)--(10.1107) contains `|H|`.  This is a real finite wall, not
just missing proof technique.  The exact order-five minimizer

```math
A_5=\begin{pmatrix}
0&-1&-1&-1&-1\\
-1&0&-1&-1&1\\
-1&-1&0&1&-1\\
-1&-1&1&0&1\\
-1&1&-1&1&0
\end{pmatrix},
\qquad Q(A_5)=q_5=8,
```

has `S=0123,T=0124`, a common anchor and same-orientation labels which are
simultaneously exact parent grounds and exact `B_S,B_T` grounds.  Both have
row square `24`, both parent deficits and complement cap slacks vanish, yet
they disagree on both nonanchor overlap vertices.  Exhaustion shows no such
proper `m=n-1` zero-gap pair at orders three or four.  On the exact `A_9`
minimizer there are exact `m=7` complement maximizers of cap `28` and row
`112` with nonzero overlap disagreement but equality in (10.1103).

These examples defeat a Hamming modulus based only on the displayed scalar
gaps, row scale, or cyclic monotonicity.  They do not rule out selecting a
different coherent representative from the same fibers.  The exact positive
successor is an averaged minimizer-specific **normal-fan curvature** or
cross-overlap selection theorem: a common penalty must preserve complement
incidence and the row budget while its Bregman gap controls corrected anchored
disagreement.  Cyclic monotonicity without such a gap does not prove the
conflict clause of (10.1089).

#### 10.93.2 Harmonic curvature has exact differential calculus but needs signing structure

For coordinate `i`, context `e`, and the two conditional score values, put
`delta=g_+-g_-` and let `p_s` be the time-`s` conditional probability of the
upper endpoint.  Direct binary differentiation gives

```math
\boxed{
v_{i,s}(e)=p_s(1-p_s)\delta^2,
\qquad
\partial_s v_{i,s}(e)=p_s(1-p_s)(1-2p_s)\delta^3.
}
\tag{10.1108}
```

The second expression is the third conditional central moment `A'''_{i,e}`.
Since the conditional law of `D` given `g` is time-independent,

```math
\boxed{
\partial_s m_s(u)
=\mathbb E_\nu\left[\sum_i A'''_{i,D_{-i}}(s)\mid g=u\right].
}
\tag{10.1109}
```

There is also the exact two-replica representation

```math
v_{i,s}(e)
=\frac12\mathbb E[(g(B)-g(B'))^2\mid e],
\tag{10.1110}
```

where `B,B'` are independent from the time-`s` binary conditional.  In the
matched restriction endpoint, its jump remains the exact exclusion log ratio
`delta=log(r_i(x)/r_i(y))`.  These formulas give no sign to (10.1109).

Indeed literal isotonicity is **Numerically false** on most or all quadrature
nodes for the audited `A_6,A_8,A_9` cases.  The weighted defect can still be
small: for `A_9,beta=2,m=7`, `K^2=0.000997862...` versus endpoint vertex cost
`0.050868191...`; for `A_8,beta=2,m=4`, the corresponding values are
`0.002631271...` and `0.596918558...`.  These finite ratios are evidence only.

There is a sharp abstract wall to promoting the numerical comparison to a
generic theorem.  Start with two chart bits, uniform base law, and score

```math
g_L(00)=0,
\qquad g_L(01)=g_L(10)=g_L(11)=L.
```

Writing

```math
v_s=L^2\frac{e^{sL}}{(1+e^{sL})^2},
```

one has statewise curvature `(2v_s,v_s,v_s,0)` and hence the antitone scalar
regression

```math
m_s(0)=2v_s,
\qquad m_s(L)=\frac23v_s.
```

If `w_s` is the low-score mass under `(lambda_s+lambda_1)/2`, exact two-point
isotonic regression gives

```math
\boxed{
\overline{R}_s=w_s(1-w_s)\left(\frac43v_s\right)^2.
}
\tag{10.1111}
```

The substitution `u=sL` proves

```math
\boxed{
K_L^2\sim c_*^2L,
\qquad
C_V(1)\longrightarrow\frac23\log2,
\qquad
H\longrightarrow\log\frac43,
}
\tag{10.1112}
```

with `c_*=0.1234111478...>0`.  This wall even has a positive fixed-size
erasure lift.  Add an inert third bit, fix `0<\varepsilon<1/2`, and define

```math
\boxed{
r_1=r_2=\varepsilon e^{-g_L},
\qquad r_3=1-2\varepsilon e^{-g_L},
\qquad b_i=e^{g_L}r_i.
}
\tag{10.1113}
```

Here `\sum_i r_i=1`, `b_1=b_2=\varepsilon`, and
`b_3=e^{g_L}-2\varepsilon`; every `b_i` is positive and independent of its omitted
coordinate, `sum_i b_i=e^(g_L)`, and all matched jump identities hold.  The
same example has

```math
\frac{\operatorname{Var}_{\mu_s}(g_L)}{\mathbb E_{\mu_s}V_s}
=\frac{3(1+e^{sL})}{2(1+3e^{sL})}\le\frac34,
```

bounded time-weighted scalar Dirichlet integral, and actual adverse signed
migration tending only to `0.05133578656...`.  Thus the isotonic-information
majorant itself can lose a growing factor.  No bound whose right side stays
bounded with endpoint cost, parent entropy, restoring, and scalar Dirichlet
cost, even after imposing the generic matched-erasure identities, can prove
(10.1096).  The construction is not a quadratic-signing Gibbs endpoint: its
positive component likelihoods are arbitrary and its uniform base law is not
supplied by a nonzero-temperature exact-minimizer signing.  The harmonic route
therefore survives, but its next theorem must use the actual
outside-completion polynomial/Gibbs form or exact-minimizer restrictions, or
control signed migration more directly.

#### 10.93.3 Favorable distance is an exact soft terminal potential, but no new estimate

Make every constrained terminal a leaf by permitting free full-projective
Steiner nodes.  For a selector `S` and its nonempty actual favorable fiber
`F_S`, define

```math
\pi_S(y)=d_{\rm pr}(y,F_S).
```

Changing one terminal label changes the optimal completion-tree length by at
most its projective Hamming change.  Nearest-fiber projection therefore proves
the exact infimal-convolution identity

```math
\boxed{
\min_{\mathbf y}\left\{D_T(z;\mathbf y)
+\sum_v d_{\rm pr}(y^v,F_v)\right\}
=\min_{\mathbf f\in\prod_vF_v}D_T(z;\mathbf f).
}
\tag{10.1114}
```

The left domain is fully flip-closed, and a minimizer may be chosen favorable.
For its edge disagreement set `U_e`, (10.1099) becomes

```math
\boxed{
|U_e|\le\sum_{v\in H_e}
d_{\rm pr}((y^v)^{U_e\cap S_v},F_v).
}
\tag{10.1115}
```

The apparent ancestor reuse is removed globally, not by summing (10.1115):
projecting all leaf labels changes the tree length by at most the **single**
sum of their fiber distances.  Equivalently, the min-plus descendant message

```math
M_u(a)=\min_x\left\{d_{\rm pr}(a,x)+\pi_u(x)
+\sum_{w:\,w\text{ child of }u}M_w(x)\right\}
\tag{10.1116}
```

computes both sides of (10.1114) exactly.

For any unconstrained scaffold `bar y`, the theorem yields

```math
\boxed{
\min_{\mathbf f\in\prod_vF_v}D_T(z;\mathbf f)
\le D_T(z;\bar{\mathbf y})
+\sum_v d_{\rm pr}(\bar y^v,F_v).
}
\tag{10.1117}
```

The common-root scaffold is exactly the total-cost center star (10.926), and
minimizing over all scaffolds returns the original hard forest target
(10.923).  Thus the soft potential repairs flip admissibility and charging,
but does not supply the missing project estimate.  The exact successor is an
`O(k_0)` bound on the right side of (10.1117), per group on the required
`exp{-O(rL_0)}` event, from a canonical minimizer-specific scaffold.

Generic log-sum-exp smoothing also does not create it.  Eliminating one full
projective state costs at most `tau(n-1)log2`; over `O(s)` states its error is
`O(tau n s)`.  At project parameters an `O(k_0)` conclusion therefore forces

```math
\boxed{
\tau=O\left(
\frac{n^{-1/2-2c_0+\eta}}{(\log n)^2}
\right)
}
\tag{10.1118}
```

without a new conditional-entropy theorem, essentially returning to min-plus.
The smallest boundary example is an exact order-four minimizer with
`S=012,T=023`: its actual `t=0` favorable fibers are singleton projective
words `+++` and `+--`, their completion cylinders have distance one, and a
row-`20` root realizes that cost.  Both local singleton log normalizers are
zero while the global compatibility partition is `e^(-lambda)`.  This
falsifies cancellation of boundary compatibility from local normalizers, not
the correct dynamic program, which retains the nonzero root message.

#### 10.93.4 Updated frontier

Wave 40 leaves the strategic leader unchanged but closes several tempting
shortcuts:

- the low-information, low-row, low-conflict kernel (10.1089) remains the
  strongest structured implementation of the bare tail.  Penalized canonical
  complement grounds are cyclically monotone, but exact minimizing signings
  show zero cap slack, low row, four-corner exchange, and even a flat normal
  cone do not price overlap disagreement.  Continue only with a genuinely
  quantitative normal-fan-curvature/cross-overlap selection input, or with a
  direct exceptional-center/nonlinear rare-tail theorem.  The finite walls do
  not falsify existence of a coherent joint selection;

- the harmonic scalar target (10.1096) survives, but literal monotonicity is
  retired.  The exact derivative and two-replica formulas do not control its
  sign, and a positive fixed-size erasure mixture has unbounded `K^2` with
  bounded endpoint cost and entropy.  Any next proof must exploit the actual
  quadratic outside-completion/Gibbs form or exact minimizer, beyond the
  generic matched exclusion identities;

- favorable distance is the correct exact soft terminal potential.  It solves
  flip closure and ancestor accounting without loss, but is equivalent to the
  original favorable forest and center/scaffold targets.  Generic finite
  temperature has an entropy error at the wrong scale.  Terminal min-cut
  remains a fallback rather than a new leader;

- no convergence proof or asymptotic falsifier has appeared, the rigorous
  interval is unchanged, and the Wave 40 blank-slate judgment remains valid.
  The next mandatory `STEERING.md` refresh and blank-slate audit are still due
  at the Wave 44 boundary unless a decisive result intervenes.

Wave 41 should compare three independent attacks: an entropy-regularized
complement kernel with an explicit quantitative cross-overlap price; a
signing-specific exclusion of the harmonic erasure wall using the actual
outside-completion polynomial; and a direct nonlinear exceptional-center or
other rare-tail theorem.  Do not retry unpriced cyclic monotonicity, generic
endpoint-cost control, or log-sum-exp terminal smoothing.
