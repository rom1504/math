### 10.110 Wave 57: sector cancellation, Johnson-correlated reuse, and the row-incidence dual

Wave 57 audited the three bridges isolated in Wave 56.  The active-face
escape has an exact sector-disagreement interpretation, but an exact finite
minimizer defeats every pointwise escape-to-favorability implication using
only sector, PSD, and full-energy identities.  At the mesoscopic core, the
missing quantity is now an exact size-biased same-witness retention
probability; balanced incidence partitions falsify degree, distinctness,
bounded multiplicity, and ordinary dependent-random-choice packages.  For
the selected prior, KL projection reduces the row tradeoff to a two-dimensional
column convex hull and an exact scalar separation certificate.  High row in
the correct energy sector forces migration, but neither row descent nor
favorable-fibre stability follows.

The sector, mesoscopic, and row-incidence checkers were rerun independently
and passed.  Algebraic identities and LP statements below are **Verified**.
The `A_8,A_9` statements are **Exact finite** obstructions; the balanced
partition and hub constructions are **Abstract scalable** incidence
obstructions, not exact-minimizer sequences.  No convergence or new
restriction estimate is proved.

#### 10.110.1 Escape is positive-edge disagreement, while favorability is signed cancellation

For an oriented child ground `(sigma,y)` on `S` and an oriented parent state
`(tau,x)`, put

~~~math
c_{ij}=\sigma a_{ij}y_iy_j,\qquad z_i=x_iy_i,
\qquad\kappa=\tau\sigma,
~~~

and define the relative sector-disagreement set

~~~math
D_S=\{ij\subset S:\ \kappa z_iz_j=-1\}.
~~~

It is a cut when `kappa=+1` and an uncut when `kappa=-1`.  If `E_S` is any
legal child-positive excess block, direct edge accounting proves the
**Verified unified sector identities**

~~~math
\boxed{
\delta_S=4\sum_{e\in D_S}c_e,
\qquad
S_{E_S}=|E_S|-2|E_S\cap D_S|.
}
\tag{10.1358}
~~~

Thus active-face escape `S_{E_S}<=|E_S|/2` means
`|E_S cap D_S|>=|E_S|/4` in either sector.  Writing `P_D,N_D` for the
numbers of child-positive and child-negative edges in `D_S`, bare
favorability is exactly

~~~math
\boxed{
4(P_D-N_D)\le B_{n,m}+p_2\Delta+C T_n.
}
\tag{10.1359}
~~~

Escape controls `P_D` from below but says nothing about the compensating
`N_D`.  The canonical second perturbation also has the wrong sign.  If the
internal edges in `D_S` are flipped, the same parent state has energy

~~~math
\boxed{
E_{A^{D_S}}(\tau,x)=q-\Delta+\delta_S.
}
\tag{10.1360}
~~~

Large local deficit therefore supports, rather than contradicts, the lower
cap required by exact signing minimality.

This pointwise gap is real on the stored exact order-nine minimizer.  At
`m=7`, one exact parent ground has energy `24`, row `96`, and full edge sum
`12=q/2`.  With canonical child grounds and excess blocks it escapes `28/36`
selectors, is zero-slack favorable on `5/36`, and the two sets are disjoint.
The escaped `(kappa,delta)` counts are

~~~text
(-1,4): 8,   (-1,12): 3,   (+1,8): 14,   (+1,16): 3.
~~~

Here `B_(9,7)=56sqrt(7)/9-14=2.46245...`, so every escaped incidence has
strictly positive `hat ell=delta-B`.  More strongly, `15/36>1/3` selectors
are unfavorable escapes for every legal choice of child ground and excess
block.  This example satisfies the total-energy identity and has a rank-one
PSD orientation matrix.  It is an **Exact finite falsification of pointwise
or identity-only bridges**, not an asymptotic `T_n`-scale falsifier.

The precise aggregate successor is to prove, for the particular common law
`mu` in (10.1341),

~~~math
\boxed{
(U_m\otimes\mu)\left\{
S_{E_S}\le |E_S|/2,
\ 4(P_D-N_D)\le B_{n,m}+p_2\Delta+C T_n
\right\}\ge e^{-C'H},
}
\tag{10.1361}
~~~

together with captured row `O(n^(9/4-c))` on the same event.  This would
feed the selected-prior extraction, but it requires genuinely joint
selector--state cancellation or a common perturbation; sector splitting and
marginal balance cannot prove it.

#### 10.110.2 Mesoscopic excess is exactly Johnson-correlated witness reuse

For a controlled-row center `z`, write `F_z` for its favorable selector
family and `r_z=|F_z|`.  For `ell<=s<m`, let `Gamma_s(S)` be the `D_s`
distinct selectors `T` with `|S cap T|>=s`.  If

~~~math
E_s(z)=|\{(S,T)\in F_z^2:S\ne T,\ |S\cap T|\ge s\}|,
~~~

then the size-biased probability that an `s`-overlap move retains the same
center is

~~~math
\boxed{
\Pi_s=
\frac{\sum_z r_zE_s(z)}{D_s\sum_zr_z^2}.
}
\tag{10.1362}
~~~

Diagonal separation and the minimum weight of a threshold pair prove the
**Verified exact threshold-retention lemma**

~~~math
\boxed{
P_\ell\ge h_\ell+R_{\ell,s}\Pi_s,
\qquad
R_{\ell,s}=\frac{D_s\binom s\ell}{d_\ell}.
}
\tag{10.1363}
~~~

Consequently
`Pi_s>=(lambda_2-h_ell+epsilon)/R_(ell,s)` implies
`P_ell-lambda_2>=epsilon`.  Equivalently, with weighted codegree
`c_r(S,T)=sum_(z:S,T in F_z)r_z`, the numerator of `Pi_s` is the total
`Gamma_s`-edge codegree.  This is the exact shared-witness quantity which
ordinary incidence degrees omit.

At the clean choice `c=1/8`, `H=n^(5/8)` and
`ell=kappa n^(13/16)`, the optimal threshold has

~~~math
\boxed{
\Pi_{s_*}\ge
(1+o(1))\sqrt{2\pi}\,a^3\kappa^3n^{-1/16}
+e^{-(C+o(1))H},
\qquad a=\frac{1-p}{p},
}
\tag{10.1364}
~~~

as the exact one-threshold sufficient target.  The polynomial term is the
degree-two baseline; only the final saved margin is needed.

There is a sharp generic obstruction.  Partition the selector slice
uniformly into fixed blocks of sizes `r_1,...,r_B`, and use each block as one
center fibre.  Exact pair counting gives

~~~math
\boxed{
\mathbb E P_\ell=h_\ell+(1-h_\ell)q_{\rm part},
\qquad
\mathbb E\Pi_s=q_{\rm part},
\qquad
q_{\rm part}=
\frac{\sum_i r_i^2(r_i-1)}{(N-1)\sum_i r_i^2}.
}
\tag{10.1365}
~~~

Taking `B=ceil(e^H)` balanced blocks gives
`q_part=e^{-(1+o(1))H}<<lambda_2`.  Hence some deterministic partition has

~~~math
\boxed{P_\ell-\lambda_2<-\lambda_2/2.}
\tag{10.1366}
~~~

Nevertheless every selector has one witness, every center has
`(1+o(1))Ne^{-H}` distinct selectors, and every selector pair has codegree
at most one.  Labeled copies preserve the failure at any fixed left degree
and bounded pair multiplicity.  This **falsifies degree, distinctness,
bounded-multiplicity, and ordinary DRC implications**; it does not falsify a
signing-specific theorem forcing the Johnson-correlated reuse (10.1364).

The 2025 Khot--Minzer--Moshkovitz--Safra small-set-expansion theorem for the
Johnson graph supports the core-structure heuristic but does not supply this
claim: its stated regime takes the vertex-set size small relative to the
ambient dimension after sequential limits, treats fixed expansion defect,
and provides no uniform additive `e^{-Theta(H)}` control on the present
fixed-density diagonal.  The applicable input remains the exact core-load
identity (10.1266), not an unstated quantitative extension of that theorem.

#### 10.110.3 The row bridge is a two-dimensional column convex hull

For the favorable incidence `F_t`, put `u_d=U_m(F_t^d)`,
`h_d=log(1/u_d)`, and let `U_d` be uniform on the column fibre.  Every
supported joint law with output marginal `pi` obeys the **Verified KL
projection identity**

~~~math
\boxed{
K(P)=I(S;D)+D(P_S\Vert U_m)
=\sum_d\pi_d\{h_d+D(P_{S|d}\Vert U_d)\}.
}
\tag{10.1367}
~~~

Thus the lower attainable information--row region is exactly the convex hull
of the column points `(h_d,R_d)`.  The least mean row at information budget
`H_0` and its exact scalar dual are

~~~math
\boxed{
\mathcal R(H_0)
=\min_{\pi:\,\mathbb E_\pi h\le H_0}\mathbb E_\pi R
=\sup_{\lambda\ge0}
\left[\min_d\{R_d+\lambda h_d\}-\lambda H_0\right].
}
\tag{10.1368}
~~~

An optimizer uses at most two columns.  Failure at row budget `R_0` is
equivalent to a nonnegative separating price with
`min_d{alpha(h_d-H_0)+beta(R_d-R_0)}>0`; this, rather than the verbal
statement that low-information columns have high row, is the exact
obstruction.

The corresponding selected-prior LP reduces to

~~~math
\boxed{
Z_{\rm sel}(R_0)^{-1}
=\min_{\pi:\,\mathbb E_\pi R\le R_0}
\mathbb E_\pi e^{h_D}.
}
\tag{10.1369}
~~~

For a fixed reference prior `nu`, row tilting has the exact Gibbs dual

~~~math
\boxed{
-\log\sum_d\nu_du_de^{-\gamma R_d}
=\inf_{P:P(F_t)=1}
\{K(P)+D(P_D\Vert\nu)+\gamma\mathbb E_PR_D\}.
}
\tag{10.1370}
~~~

Optimizing `nu` collapses to `min_d(h_d+gamma R_d)`.  Arbitrary-prior KL
therefore never controls row without a proved reference moment or an
explicit row constraint.

High row does produce one exact-minimality response in the correct sector.
For signed parent fields `r_i=sum_(j ne i)s_ij`, if `Delta(d)<=q`, then

~~~math
R_d\le2(n-1)\sum_i(r_i)_+.
~~~

If `R_d>=2(n-1)L` with `L asymp n^(5/4-c)`, a positive-edge block of size
between `L/2` and `L` exists.  The localized migration theorem gives

~~~math
\boxed{
\Delta(\omega)=O(n^{9/8-c/2})=o(T_n),
\qquad
|\{e\in E:s_e(\omega)=-1\}|>|E|/4.
}
\tag{10.1371}
~~~

It does not bound `R(omega)` or preserve the original fibre.

Both omissions are sharp for current information.  On exact `A_8,m=5`, all
eight maximum bare-fibre columns have

~~~math
\boxed{
u_t=40/56,\qquad R=64,\qquad\Delta=2q=40,
\qquad\sum_i(r_i)_+=0.
}
\tag{10.1372}
~~~

They lie in the opposite sector; reversing orientation activates the fields
but reduces coverage.  Separately, an abstract family avoiding `k` hubs has

~~~math
\boxed{
-\log\frac{\binom{n-k}m}{\binom nm}
=k\log\frac1{1-p}+O(k^2/n)=o(H),
\qquad
R\asymp kn^2\gg n^{9/4-c}
}
\tag{10.1373}
~~~

when `k=omega_n n^(1/4-c)` with a suitable slowly diverging `omega_n`.
The favorable selectors can avoid every edge supporting the high row, so
entropy alone cannot make the migration response fibre-stable.  This is an
**Abstract incidence wall**, not an exact-minimizer construction.

The exact surviving row target is a fibre-stable row descent at the
separating price: a row-high favorable column `d` must have a response
`omega` satisfying

~~~math
h_\omega+\lambda R_\omega<h_d+\lambda R_d,
~~~

or, equivalently for the present purpose, a saved overlap of favorable
fibres together with a controlled-row response.  Neither follows from
one-block migration.

#### Updated frontier

1. **The selected-prior lemma remains the sharp leader, but its obstruction
   is now exact.**  The feasible information--row pairs form the convex hull
   of `(h_d,R_d)`; a failure is a scalar separating price.  The next row
   theorem must descend that scalar objective or construct the two-column
   selected prior directly.  High row by itself supplies migration, not
   descent.

2. **Common active-face escape is not a pointwise favorable bridge.**  Bare
   favorability requires child-negative cancellation inside the relative
   sector-disagreement set.  Exact `A_9` defeats sector, PSD, total-energy,
   and canonical-second-flip arguments simultaneously.  Only the aggregate
   cancellation-and-row statement (10.1361) survives.

3. **Mesoscopic one-threshold work now has one diagnostic.**  Any proposed
   cluster mechanism must be evaluated by the size-biased same-center
   Johnson retention `Pi_(s_*)` in (10.1364).  Target-size fibres, distinct
   selector labels, bounded witness multiplicity, and ordinary DRC are all
   insufficient.  The full histogram can combine shells but still needs
   Johnson-correlated same-witness reuse.

4. **Wave 58 should test mechanisms, not restate interfaces.**  The strongest
   independent targets are: (i) a fibre-stable row-descending migration at
   the LP separating price; (ii) a signing-specific source of the exact
   Johnson codegree (10.1364), possibly from aggregate sector cancellation;
   and (iii) direct construction of the optimal two-column prior in
   (10.1369).  A separate route should revisit exact-minimizer exchange only
   if it controls one of these quantities.

No decisive leading-route change occurs, so `STEERING.md` is not refreshed
at this boundary.  The next mandatory refresh remains Wave 61 and must
include a blank-slate abstraction audit.
