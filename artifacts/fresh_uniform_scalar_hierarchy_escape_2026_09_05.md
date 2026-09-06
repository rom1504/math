# Uniform strict escape from the entire marked scalar-mask hierarchy

Date: 2026-09-05. This is a variational incompleteness theorem, not a
convergence theorem. It strengthens escape from one specified certificate
to a strict gap above the supremum of all old scalar central-mask
certificates.

## 1. The variational family and claim

Use the countable Gaussian marked-tree space and its creation isometry
`U`. Let `G0=U1`. For every unit first-chaos Gaussian variable `V` and
every `alpha>0`, define

`H=1{|V|<=alpha}`, `W=U H`, `F=sign(W)(1-H)`,

`J(V,alpha)=E F W=E|W|(1-H)`.

Let `C_scalar` be the supremum of these values over the stated family.
This allows arbitrary unit first-chaos `V`, not only solutions obtained
by a scalar or finitely anchored contraction iteration. All such values
have finite smooth Boolean realizations by the ordered approximation
argument already proved.

The banked 21-anchor construction implies `C_scalar>0.4306`. The theorem
proved below is

`liminf_n M_n/n^(3/2) > C_scalar`.                            (1)

No explicit numerical size for this strict gap is asserted. In
particular it does not show that the enlarged hierarchy of marked and
unmarked responses is complete for the Boolean problem.

## 2. A compact mask-mass range near the scalar supremum

Put `p=E H=2Phi(alpha)-1`. Isometry gives `Var W=p`, and rearrangement
of the folded-normal variable `|W|` gives

`J(V,alpha) <= R(p):=2 sqrt(p) phi(alpha(p))`,

where `alpha(p)=Phi^(-1)((1+p)/2)`. Also `J<=sqrt(p)` by Cauchy--Schwarz.
Thus `J>=0.43` implies `p>=0.43^2>0.18`.

For clarity, the upper mass bound uses only generous exact rational
inequalities. The previously audited Gaussian interval functions give

`Phi(22/25)<13/16<Phi(9/10)`, `phi(22/25)<271/1000`,

and `(791/1000)^2>5/8`. Hence

`22/25<alpha(5/8)<9/10`,

`R(5/8)<2*(791/1000)*(271/1000)=0.428722<0.43`.

These two Gaussian bounds were independently rerun with the exact
`Fraction` interval functions from `fresh_limit_rooted_lower_certificate.py`:

`Phi(.88) in [.810570345223287870723974356154005569074784419623356383015321,`
`             .810570345223287870723974356154005569074784419623356383015344]`,

`phi(.88) in [.270863971798338004780850036241922084730592234050324090273384,`
`             .270863971798338004780850036241922084730592234050324090273409]`.

For every `p>=5/8`,

`R'(p)=phi(alpha)/sqrt(p)-alpha sqrt(p)<0`,

since `p alpha >=(5/8)(22/25)>.271>=phi(alpha)`. Therefore every
certificate with `J>=.43` lies in the fixed compact mass interval

`I=[9/50,5/8]`.                                               (2)

## 3. A uniformly positive edge-linear response

Write `rho=Cov(G0,V)` and `w=Cov(V,W)`. Because `G0=U1` and `W=UH`,

`Cov(G0,W)=p`.

The first Hermite coefficient `a1=E[G0 F]` has the uniform lower bound

`a1 >= sqrt(2p/pi)-p`.                                       (3)

Indeed the unmasked sign has
`E[G0 sign(W)]=sqrt(2p/pi)`. On the masked event,

`E[G0^2 H]=p-2rho^2 alpha phi(alpha)<=p`,

so Cauchy--Schwarz gives
`E[G0 sign(W) H] <= E|G0|H <=p`. Subtract to obtain (3).

On the entire interval (2), the right side of (3) exceeds `1/200`.
One elementary check uses `pi<22/7`, so it is bounded below by
`sqrt(7p/11)-p`. This is concave in `p`, and both endpoint values exceed
`1/200`: at `p=5/8`, use `sqrt(35/88)>63/100`; the other endpoint is
much larger. Consequently

`a1>=1/200` for every certificate under consideration.         (4)

This argument needs no sign assumption on `rho` or `w`.

## 4. Uniform nondegeneracy from the tree-height obstruction

Replace `V` by `-V` if necessary so that `w>=0`; the central mask, `W`,
and certificate value do not change. For each mass `p`, let

`g_p(v)=1{|v|<=alpha(p)}/sqrt(p)`.

Its Hermite coefficient squares define the kernel

`K_p(q)=E[g_p(N)g_p(N')]`, `Cov(N,N')=q`.

For every `0<p<1`, there exists `q_p<1` with `K_p(q_p)<q_p`. For example,
the threshold has infinite Gaussian Dirichlet energy, so
`K_p'(1-)=infinity`; equivalently its near-one noise defect has a positive
square-root leading term. This is the supercritical central-indicator
case of the already proved tree-height barrier.

For a fixed `q<1`, `K_p(q)` is continuous in `p`. The height-barrier
lower bound therefore supplies a strictly positive residual bound on a
neighborhood of each `p` in the compact interval (2). Taking a finite
subcover gives a single `delta>0` such that every unit first-chaos `V`
at every mass in (2) obeys

`||V-U g_p(V)||_2>=delta`.

We can decrease `delta` so that `delta<1`. Since `U g_p(V)=W/sqrt(p)`,

`w/sqrt(p)<=1-delta^2/2`,

`Var(V|W)=1-w^2/p >=delta^2-delta^4/4>0`.                     (5)

This is uniform over all feasible `V`, with no bounded-depth assumption.

## 5. Compactness supplies a finite nonzero Hermite direction

Consider the closed parameter set of all triples `(p,rho,w)` satisfying
`p in I`, `w>=0`, (5), and positive semidefiniteness of

`Cov(G0,V,W)=[[1,rho,p],[rho,1,w],[p,w,p]]`.

It is compact. The pair `(V,W)` is uniformly nondegenerate, although the
full triple is allowed to be singular. On this parameter set, define

`F_theta=sign(W)1{|V|>alpha(p)}`,

`b_r(theta)=E[F_theta h_r(G0)]`.

Every `b_r` is continuous. One way to check this even at singular full
triples is to keep `G0` fixed and represent the conditional Gaussian pair
`(V,W)|G0` by the positive square root of its two-dimensional covariance.
The square root varies continuously. The boundaries `W=0` and
`|V|=alpha(p)` have zero probability, so dominated convergence and fixed
Gaussian moments prove continuity of the displayed coefficients.

For every parameter in the set, some odd `r>=3` has `b_r(theta)!=0`.
Otherwise Hermite completeness and oddness imply

`E[F_theta|G0=g]=b_1(theta) g` almost everywhere.

Its left side is bounded by one, whereas (4) gives `b_1(theta)>=1/200`;
this is impossible. The open sets `{theta:|b_r(theta)|>0}` therefore
cover the compact parameter set. A finite subcover supplies a fixed
finite set of odd degrees and a constant `b>0` such that

`max_(3<=r<=R, r odd) |b_r(theta)|>=b`                         (6)

uniformly. This is a genuine finite-degree selection argument, not a
claim that an arbitrary bounded response has a uniformly large cubic.

## 6. Uniform gain over the whole near-optimal family

Apply the all-odd weighted projection theorem from
`fresh_all_odd_weighted_projection_2026_09_05.md`, choosing for each
target a degree and sign supplied by (6). The odd variance normalization
has mean normalized coefficient at least `1/2` and local Gaussian
variance at most one, uniformly in the chosen finite degree set.

Every Gaussian constant in the normalized-gain argument is now uniform
over the parameter set: the mask mass is compactly bounded; (5) supplies
a positive joint density on a fixed rectangle outside all central masks;
(6) supplies a uniformly nonzero Hermite coefficient; and Gaussian
threshold boundaries permit uniform smoothing estimates. One can for
example choose the `V` cutoff inside `[1,2]` and its negative, since
all `alpha(p)<.9` on (2), together with a `W` cutoff of width `t`.

The finite smooth approximation index may depend on the chosen target
and on `L`; that causes no loss because it is chosen before `n` grows.
Exactly the argument in
`fresh_normalized_unmarked_gain_and_limit_order_2026_09_05.md` therefore
gives constants `c>0,L0<infinity`, independent of the target, such that
for each fixed `L>=L0`,

`liminf_n M_n(L)/n^(3/2) >= C_scalar+c/log L`.                (7)

For precision, first choose a scalar certificate within an arbitrarily
small fixed fraction of `c/log L` of its supremum, choose one finite
smooth realization with another such small loss, and only then take the
matrix-size limit. One can decrease `c` once to absorb these losses.
The same positive constant works for every sufficiently large `L`.

Principal spectral deletion now converts (7) to the unrestricted
problem. On a low-cap minimizing sequence, retain at least
`(1-epsilon)n` vertices with normalized operator bound
`L=C_1/epsilon`. The parent's Boolean optimum dominates that of the
principal submatrix. Hence

`liminf_n M_n/n^(3/2)`

` >= (1-epsilon)^(3/2) [C_scalar+c/log(C_1/epsilon)]`.

For sufficiently small fixed `epsilon`, this exceeds `C_scalar`, proving
(1). The use of a supremum introduced no exchange of a growing tree
degree with the matrix limit: compactness fixed a finite degree set
before the target and finite approximation were chosen.

## 7. What the result does not settle

The scalar marked hierarchy is rigorously incomplete for the universal
Boolean lower bound. The theorem does not identify its supremum
numerically, nor does it identify the best enlarged marked/unmarked
variational family. In particular there is still no theorem that finite
local Gaussian observables determine the global Boolean maximum, no
matching upper construction, and no proof that `M_n/n^(3/2)` converges.

Independent audit record: root checked the full theorem, literature
independently checked Sections 2--3 and all rational endpoint bounds,
and algebra independently checked Sections 4--6, including the enlarged
compact parameter set and finite Hermite-degree selection. The upper
bound `alpha(5/8)<.9` used for the fixed cutoff was separately verified
by the same exact `Fraction` Gaussian interval code.
