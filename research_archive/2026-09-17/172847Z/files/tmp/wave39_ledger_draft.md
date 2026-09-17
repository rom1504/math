### 10.92 Wave 39: information-to-Gram transport, scalar curvature transport, and terminal cut compression

Wave 39 attacks the three interfaces selected at the mandatory strategic
boundary.  It removes internal favorable Gram as an independent obstruction:
low selector--certificate information and low parent row already transport to
the **full** selected-column Gram.  The harmonic route gains an exact
one-dimensional Bregman-curvature transport and a scalar isotonic-defect
bootstrap.  Exact terminal-min-cut screening gains a no-Johnson subtree
Euler calculus, but favorability and signed laminar charging remain open.  A separate row-price audit gives a clean
pointwise complement-witness inequality but confirms that pointwise row
control is not cross-selector congestion.  The two analytic reductions were
independently reconstructed; all finite checkers pass, and numerical LP or
quadrature values below are labelled numerical.

#### 10.92.1 Selector--certificate information transports parent row to full Gram

Let `P(S,D)` be any joint law of an `m`-selector and a full projective word
`z^D`.  Put

```math
\mathsf H=D(P_S\Vert U_m)+I(S;D),
\qquad
\overline R=\mathbb E R_2(z^D).
```

Then every exact minimizer satisfies the uniform transport inequality

```math
\boxed{
\begin{aligned}
\mathbb E\lVert A[:,S]z^D_S\rVert_2^2
\le{}&p^2\overline R+p(1-p)n(n-1)\\
&+C\sqrt{q_n(n^2+\overline R)
                 (\mathsf H+\log(n+1))}\\
&+Cq_n(\mathsf H+\log(n+1)).
\end{aligned}}
\tag{10.1087}
```

To prove it for a fixed `z`, set `G_z=D_zA^2D_z` and let `xi` have iid
`Ber(p)` coordinates.  Then

```math
K_\xi(z)=\xi^TG_z\xi,
\qquad
\mathbb EK_\xi=p^2R_2(z)+p(1-p)n(n-1).
```

The centered Bernoulli expansion splits into a linear part along
`G_z 1` and a centered quadratic part.  Bounded-variable Hanson--Wright,
the linear subgaussian mgf, and Cauchy--Schwarz combine because

```math
\boxed{
\lVert G_z\rVert_{op}\le2q_n,
\quad
\lVert G_z\rVert_F^2\le2q_n n(n-1),
\quad
\lVert G_z\mathbf1\rVert_2^2\le2q_nR_2(z).
}
\tag{10.1088}
```

Conditioning on `sum xi_i=m` gives `U_m` and costs at most `log(n+1)`
because `m` is a binomial mode.  Entropy duality yields the fixed-word
estimate.  Applying it under `P(S|D=d)`, using

```math
\mathbb E_DD(P_{S|D}\Vert U_m)
=D(P_S\Vert U_m)+I(S;D),
```

and Cauchy--Schwarz proves (10.1087).  At the project scales

```math
\overline R=O(n^{9/4-c}),
\qquad
\mathsf H=O(n^{3/4-c}),
```

both nonlinear error terms are `O(n^(9/4-c))`.  Thus the whole column Gram,
and hence its internal part, has the required scale.

This simplifies the complement/agreement implementation more than the
completion-Parseval bridge (10.1078).  In the high-ratio window `p_2>=1/2`,
support on complement incidences makes every restricted label favorable by
(10.1045), regardless of the parent deficit.  For an anchored family of
relative density `beta`, the exact sufficient kernel package is now

```math
\boxed{
\begin{aligned}
\log\beta^{-1}+I(S;D)&=O(n^{3/4-c}),\\
\mathbb ER_2(D)&=O(n^{9/4-c}),\\
\mathcal C_v(\pi)&\le C_{\rm conflict}(n).
\end{aligned}}
\tag{10.1089}
```

Here `C_v(pi)` is the expected corrected disagreement of two independent
selector--certificate draws after the common projective anchor is fixed.
Sampling one certificate per selector and applying the probabilistic method
to the three normalized nonnegative costs derandomizes (10.1089).  Equation
(10.1087) supplies (10.1057), and the decoder (10.1056) supplies the global
row cut.  The deficit and cap-persistence clause in (10.1080) is unnecessary
for this direct full-Gram route.

Fractional covers pay the information term exactly.  If weights `w_d` cover
the family, `Z_S=sum_(d:S in I_d)w_d>=1`, `W=sum_dw_d`, and

```math
\pi_S(d)=\frac{w_d\mathbf1_{\{S\in I_d\}}}{Z_S},
\qquad r(d)=\frac{w_d}{W},
```

then

```math
\boxed{
\mathbb E_SD(\pi_S\Vert r)
=\mathbb E_S\log\frac W{Z_S}
=I(S;D)+D(P_D\Vert r),
\qquad I(S;D)\le\log W.
}
\tag{10.1090}
```

What remains is therefore an affordable, low-parent-row complement kernel
with project-scale anchored conflict.  On the entire slice with pointwise
row cap, its weight bound already implies (10.1047), so (10.1087) is not a
circular proof of the fractional cover.  Its genuine gain is on a merely
affordable anchored agreement family, where internal Gram no longer needs a
separate theorem.  The checker verifies the exact mean and matrix identities
on `A_6,A_8,A_9` and the cover-information identity on `A_6,m=5`.

#### 10.92.2 Harmonic migration is one-dimensional Bregman-curvature transport

For a nonreference vertex coordinate `i` and context `e`, put

```math
A_{i,e}(t)=\log\mathbb E_{\nu(D_i\mid e)}e^{tg},
\qquad
v_{i,s}(e)=A''_{i,e}(s),
\qquad
V_s(d)=\sum_{i\in V_*}v_{i,s}(d_{-i}).
```

The binary conditional KL has the curvature representation

```math
k_{i,e}(t)=tA'_{i,e}(t)-A_{i,e}(t)
=\int_0^t s,v_{i,s}(e)\,ds.
```

Fubini at the endpoint, conditional variance at contemporaneous time, and
(10.1082) give the exact transport identity

```math
\boxed{
\int_0^1\operatorname{Cov}_{\mu_t}(g,L_t)\,dt
=\int_0^1s\left{
\mathbb E_{\mu_1}V_s-\mathbb E_{\mu_s}V_s
\right}\,ds.
}
\tag{10.1091}
```

This transport is exactly scalar.  Since

```math
\frac{d\mu_1}{d\mu_s}(d)
=\exp\{(1-s)g(d)-\psi(1)+\psi(s)\},
```

the conditional law of the state given `g` is independent of time.  With

```math
m_s(u)=\mathbb E[V_s(D)\mid g(D)=u]
```

and `lambda_t` the law of `g` under `mu_t`, one has

```math
\boxed{
\mathbb E_{\mu_t}V_s=\mathbb E_{\lambda_t}m_s,
\qquad
I(J;D\mid g(D))=0,
\qquad I(J;D)=I(J;g(D)),
}
\tag{10.1092}
```

where the balanced label `J` selects time `s` or `1`.  The endpoint scalar
law is an increasing-likelihood-ratio tilt of the time-`s` law, so monotone
quantile coupling gives an equivalent signed one-dimensional formula.  No
pointwise monotonicity of `m_s` is asserted.

Let `bar lambda_s=(lambda_s+lambda_1)/2`, let
`eta_s(u)=P(J=1|g=u)`, and define

```math
\overline{\mathcal R}_s
=\inf_{\phi\ {m nondecreasing}}
\mathbb E_{\overline\lambda_s}[m_s-\phi]^2.
```

The posterior `eta_s` is nondecreasing, and directly

```math
\mathbb E_{\mu_1}V_s-\mathbb E_{\mu_s}V_s
=4\operatorname{Cov}_{\overline\lambda_s}(m_s,\eta_s).
```

Comonotonicity, Cauchy--Schwarz, and binary Pinsker prove

```math
\boxed{
\left[\mathbb E_{\mu_s}V_s-\mathbb E_{\mu_1}V_s\right]_+
\le\sqrt{8\mathcal I_s\overline{\mathcal R}_s},
\qquad \mathcal I_s=I(J;D).
}
\tag{10.1093}
```

Jensen--Shannon information is at most one quarter of the Jeffreys
divergence.  For the exponential family this yields

```math
\boxed{
\mathcal I_s
\le\frac{1-s}{4}\{\psi'(1)-\psi'(s)\}
\le\frac{1-s}{4s}\mathscr H,
\qquad
\mathscr H=\operatorname{Ent}_\nu(f).
}
\tag{10.1094}
```

Consequently, with

```math
\mathcal K
=\int_0^1\sqrt{s(1-s)\overline{\mathcal R}_s}\,ds,
```

the signed migration obeys

```math
\boxed{
\left[-\int_0^1\operatorname{Cov}_{\mu_t}(g,L_t)\,dt\right]_+
\le\sqrt{2\mathscr H}\,\mathcal K.
}
\tag{10.1095}
```

Put `a_n=n^(1/2-2c)`.  If endpoint vertex cost and orientation cost are
`O(a_n)`, restoring gives
`Var_(mu_t)(g)<=C_R E_t(g)`, and

```math
\boxed{
\left(
\int_0^1\sqrt{s(1-s)\overline{\mathcal R}_s}\,ds
\right)^2=O(a_n),
}
\tag{10.1096}
```

then (10.1082), (10.1095), and restoring give, for
`h=sqrt(mathscr H)`,

```math
h^2\le O(a_n)+O(\mathcal K)h,
```

and hence `mathscr H=O(a_n)`.  This is a closed, noncircular conditional
bootstrap.  It does not prove (10.1096), endpoint cost, restoring,
orientation, or adjacent-selector Hellinger control.  An independent audit
reconstructed every sign and factor.  Quadrature on the four finite
minimizers checks the identities below `3e-15`; on `A_9,m=7,beta=2`, the
signed migration is `-0.001467944`, the isotonic-information bound is
`0.005469408`, and `K^2=0.000997854` (numerical only).

#### 10.92.3 Terminal min-cuts admit entry-free subtree flips but not laminar charging

Retain a rooted completion tree with optimal full Hamming completions `x^v`.
For an edge `e=(p,c)`, let `H_e` be its descendant component and put

```math
U_e=\{i:x_i^p\ne x_i^c\}.
```

Flip every coordinate in `U_e` at every node in `H_e`, and flip
`U_e cap S_v` in each affected partial label.  Internal descendant edges have
both endpoints flipped, all other edges are unchanged, and `e` loses every
disagreement.  Therefore

```math
\boxed{D_T(y^{(e)})\le D_T(y)-|U_e|.}
\tag{10.1097}
```

This pays no Johnson-entry distance.  Moreover the root and the edge sets
`U_e` reconstruct the completion exactly, with
`sum_e|U_e|=D_T`.  For a fixed `s`-node tree, the number of codes with cost at
most `k_0` is

```math
\sum_{j\le k_0}\binom{ns}{j}
\le\exp\left\{O\left(k_0\log\frac{ens}{k_0}\right)\right\}
=\exp\{O(L_0)\}
```

at the project scales.  This confirms that exact free-coordinate screening,
unlike (10.953), is combinatorially affordable.

There is an exact Euler interface.  Give terminal `v` an orientation sector
`epsilon_v`, deficit

```math
\delta_v(y)=Q(A[S_v])-\epsilon_v y^TA[S_v]y,
```

and signed shore `C_v(F;y)`.  Since

```math
\delta_v(y^F)-\delta_v(y)=4C_v(F;y),
```

a minimizer over a descendant-flip-closed domain of

```math
D_T(y)+\gamma\sum_{v\ne0}\delta_v(y)
```

obeys, edge by edge,

```math
\boxed{
|U_e|\le4\gamma\sum_{v\in H_e}
C_v(U_e\cap S_v;y),
\qquad
D_T\le4\gamma\sum_e\sum_{v\in H_e}
C_v(U_e\cap S_v;y).
}
\tag{10.1098}
```

Only each aggregate descendant sum is forced positive; individual shores
may have either sign.  More generally, for additive soft potentials `pi_v`,
global minimality gives

```math
\boxed{
|U_e|\le\sum_{v\in H_e}
\{\pi_v((y^v)^{U_e\cap S_v})-\pi_v(y^v)\}.
}
\tag{10.1099}
```

The exact reduced lemma is now: on the required
`exp{-O(rL_0)}` event, find in every group a low-row root, tree, additive
potentials, and penalized minimizer whose labels remain in their actual
favorable fibers and whose total laminar potential increment in (10.1099) is
`O(k_0)`.  Then `D_T=O(k_0)` and (10.923)--(10.924) close the forest input.
This remains open for two precise reasons.  Hard favorable fibers are not
closed under descendant flips, and a terminal is reused in every ancestor
shore; the signed shore terms are evaluated at the same word and do not
telescope as successive energy differences.

The cancellation wall can be made with complete sign coefficients, although
not with a known exact order-minimizer.  Partition `n=gb`, take a
`g`-vertex signing `K` with operator norm `O(sqrt g)`, put

```math
H=\begin{pmatrix}R&-R\\-R&R\end{pmatrix},
\qquad H\mathbf1=0,
\qquad \lVert H\rVert_{op}=O(\sqrt b),
```

use `J_b-I_b` on diagonal blocks and `K_(ac)H` off diagonal.  If
`kappa=||K||op||H||op/b<=1-2/b`, decomposition into block means and
orthogonal parts gives

```math
\boxed{
Q(A)=n(b-1),
\qquad
R_2(x)=n(b-1)^2
}
\tag{10.1100}
```

for every block-uniform ground.  Taking `b=Theta(sqrt n)` gives the correct
`Q=Theta(n^(3/2))` and `R_2=O(n^2)` scales.  Two such grounds can be
projectively `Theta(n)` apart across a zero signed shore.  Thus arbitrary
exact grounds are not cut-compressed by deficits, shore totals, cap scale,
and row scale alone.  The construction is not known to minimize `Q` globally
and its fibers also contain compatible common grounds, so it does not
falsify a joint actual-minimizer selection theorem.  The checker verifies
(10.1097), the factor four, sparse reconstruction, and an exact order-16
instance with `Q=112`, ground row `784`, zero shore, and distance `8`.

#### 10.92.4 A common row price gives pointwise witnesses, not congestion

For a selector `S`, maximize over oriented cuts

```math
\Phi_{S,\alpha}(d)=\alpha L_S(d)-R_2(d),
\qquad \alpha\ge0.
```

Write

```math
T_i=x_i(A^2x)_i-(n-1),
\qquad v_i=\sigma x_i(B_Sx)_i.
```

The one-vertex increments are `-4T_i` for row and `-4v_i` for flipped
energy.  Global optimality gives `T_i<=alpha v_i`; summing proves

```math
\boxed{
R_2(d)\le n(n-1)+\alpha L_S(d)
\le n(n-1)+3\alpha q_n.
}
\tag{10.1101}
```

Let `alpha_S^*` be the first price at which some optimizer has
`L_S(d)>=q_n`.  Optimizer energy is nondecreasing in `alpha`.  Therefore

```math
\max_{|S|=m}\alpha_S^*=O(n^{3/4-c})
\tag{10.1102}
```

would produce a project-row complement witness for every selector.  It
would not control how many distinct cuts are selected or the fractional
cover weight.  Exact enumeration gives maximum critical prices `0,14/3,8`
on `A_6,m=5`, `A_8,m=6`, `A_9,m=7`.  At those common prices the numerical
fractional cover values are respectively `2,16,36`; the unrestricted `A_9`
value was `13.75`.  Thus common row pricing can worsen collision even at
finite order.  These data neither falsify (10.1102) nor the project cover.

#### 10.92.5 Updated frontier

Wave 39 preserves adaptive optimized principal restriction and the bare
arbitrary-cut tail (10.795) as the leading framework, but materially
simplifies one structured implementation:

- favorable internal Gram is no longer an independent obstruction for a
  complement-certificate family.  The selector--certificate transport
  (10.1087) shows that project information and parent row control the entire
  column Gram.  In the high-ratio window complement support already supplies
  favorability, so neither certificate deficit nor cap persistence is needed
  for this route.  The exact remaining structured package is the affordable
  low-information, low-parent-row, low-conflict kernel (10.1089).  On the
  full slice this is at least as hard as the fractional cover (10.1047); on
  an affordable anchored family it is a genuine reduction of the earlier
  Gram package;

- the harmonic alternative gains the exact scalar curvature transport
  (10.1091)--(10.1092).  Its new migration-side sufficient target is the
  time-weighted scalar isotonic defect (10.1096), which combines with
  restoring through a closed quadratic bootstrap.  No estimate of that
  defect is known, and endpoint cost, restoring, orientation, and
  adjacent-selector Hellinger control remain separate;

- exact terminal min-cut screening really avoids the Johnson-entry wall:
  descendant flips delete coordinate separators at unit cost and give the
  laminar Euler inequality (10.1099).  It does not yet give likely cut
  compression because favorable fibers are not flip-closed and signed
  laminar increments do not telescope.  The complete-signing block wall
  rules out a theorem about arbitrary prescribed grounds, but not a
  canonical joint selection for actual minimizers.  Retain this route only
  for a favorable soft-potential theorem or a new canonical-selection input;

- common row pricing gives the exact pointwise bound (10.1101), but finite
  data show that price-optimal witnesses may have worse fractional
  congestion than unrestricted complement witnesses.  A slope theorem alone
  is not a collision theorem;

- no convergence proof or asymptotic falsifier has appeared.  Constant
  shortfall remains below the sharper routes, and the rigorous interval is
  unchanged.

Before selecting Wave 40, apply the newly requested blank-slate abstraction
audit from the original problem, update the research instructions, and
refresh `STEERING.md` with the independently generated candidates and their
ledger comparison.  The mathematical evidence favors attacking (10.1089)
and (10.1096), but the audit must be allowed to identify a genuinely stronger
composition, interpolation, compactness, exchange, or other convergence
mechanism before the Wave 40 ranking is fixed.
