# Variance-normalized local gain and the finite-approximation order

Date: 2026-09-05. This derives an asymptotic strict improvement from a
fixed Gaussian hierarchical certificate whose final odd response has a
nonzero cubic edge coefficient. All matrix limits below use a fixed finite
construction. The eventual numerical size of the strict improvement is
not evaluated here.

## 1. Inputs and target Gaussian certificate

Use the old Gaussian tree isometry `U` from even Gaussian responses to
the first Gaussian chaos. Let `G0=U1` be its edge coordinate. Suppose `V`
is a unit first-chaos Gaussian and, for a fixed `alpha>0`, set

`H=1{|V|<=alpha}`, `W=U H`,

`J*=E |W|(1-H)`.

Assume the limiting covariance triple `(V,W,G0)` satisfies

`Var W>0`, `Var(V|W)>0`,

`b*=E[sign(W)(1-H) h3(G0)] != 0`.                             (1)

These are strict finite-dimensional conditions on that triple. The last
condition requires an actual verification for the selected certificate;
it is not inferred from genericity.

Two proved matrix inputs are used:

1. For every fixed finite old construction and bounded deterministic
   `d_i`, the weighted projection identity is

   `n^-1 E[d circ M(X) circ Z]^T B F(X)`
   ` = E M * E[F h3(G0)] * n^-1 sum_i d_i v_i + o(1)`,

   where `Z=B h3(BS)` and `v_i=(B (B^2)^(circ3) B)_ii`.
   See `fresh_full_weighted_unmarked_projection_2026_09_05.md`.

2. The variance inequality is `n^-1 sum_i sqrt(v_i)>=1`.
   See `fresh_unmarked_variance_normalization_2026_09_05.md`.

Both hold for actual hollow signings. The projection theorem assumes
`||B||op<=L` for fixed `L`. Also the uniform single-root joint limit is
`(X_i,Z_i) -> (gamma,sqrt(v_i)N0)`, with `N0` independent of the old
Gaussian family; no convergence of the deterministic numbers `v_i` is
needed.

## 2. Finite smooth approximations with an exact energy identity

Choose finite-coordinate unit first-chaos variables `V_j -> V` in `L2`.
Choose smooth even functions `H_j(v)` between zero and one, equal to one
on `|v|<=alpha-delta_j`, and equal to zero on
`|v|>=alpha+delta_j`, where `delta_j->0`. Then

`H_j(V_j) -> H` in Gaussian `L2`.

Choose a finite coordinate projection `P_j` containing every coordinate
of `V_j` and `G0`, enlarged enough that

`W_j=P_j U H_j(V_j) -> W` in `L2`.

The omitted first-chaos tail `U H_j(V_j)-W_j` is independent of all
coordinates in `P_j` and has mean zero. Therefore, for every response
`F` measurable in those finitely many coordinates,

`E F U H_j(V_j) = E F W_j`.                                   (2)

This is the exact old paired Gaussian energy identity for the finite
construction, not just an approximation. In particular

`J_j=E |W_j|(1-H_j(V_j)) -> J*`,

`b_j=E[sign(W_j)(1-H_j(V_j))h3(G0)] -> b*`.

The finite matrix response will use the smooth `H_j`, not a hard mask.
Thus no discontinuous indicator is being approximated in Gaussian
`W^{1,4}`. Derivative bounds can depend on `j`: it is fixed before the
matrix size tends to infinity.

The covariance triples `(V_j,W_j,G0)` lie eventually in a fixed compact
neighborhood of their limit, in which `Var W_j` and `Var(V_j|W_j)` are
bounded away from zero. All elementary Gaussian estimates below are
uniform in that neighborhood.

## 3. Uniform localized slack and cubic coefficient

Fix a smooth odd function `psi` with `|psi|<=1`,
`psi(s)=s` for `|s|<=1/2`, and `psi(s)=sign(s)` for `|s|>=1`, with
`psi(s)` having the sign of `s`. For small `t>0`, set

`F_(j,t)=psi(W_j/t)(1-H_j(V_j))`.

The old energy loss is uniformly quadratic:

`E F_(j,t) W_j >= J_j-C t^2`,                                 (3)

where we enlarge `C` so that `C>=1`. Indeed the difference is supported
on `|W_j|<=t`, bounded by `2|W_j|`, and the density of `W_j` is uniformly
bounded.

Choose a fixed smooth even cutoff of `V_j` supported on
`|V_j|>=alpha+1`, positive on a further fixed open interval. Multiply it
by a smooth even cutoff of `W_j/t`, supported on `|W_j|<=t/4`, and equal
to one on a smaller central interval. Denote the product by `M_(j,t)`.
For all large `j` and all sufficiently small `t`,

`0<=M_(j,t)<=1`,

`kappa t <= E M_(j,t) <= C_M t`,                               (4)

with positive constants independent of `j,t,L`. This follows from a
uniform positive lower bound on the joint Gaussian density on a fixed
compact rectangle with `W_j=0`, and a uniform upper bound on the
`W_j` density. On the support of `M_(j,t)`,

`H_j(V_j)=0`, `|F_(j,t)|<=1/4`.                               (5)

Moreover, after making the upper bound on `t` smaller and `j` larger,

`|E[F_(j,t)h3(G0)]| >= b0 >0`,                               (6)

with its sign equal to that of `b*`. To check uniformity, first use
`b_j->b*`, and then Cauchy--Schwarz on the set `|W_j|<=t` to bound the
change in this coefficient by `O(sqrt(t))`, uniformly in the covariance
neighborhood. No matrix derivative constant occurs in (3)--(6).

## 4. Normalize local variances before clipping

Fix `eta=1/4` and put

`d_i=1/sqrt(max(v_i,eta))`, `Ztilde_i=sign(b*) d_i Z_i`.

Then `0<d_i<=2`, and the local limiting variance of `Ztilde_i` is at
most one. Also

`n^-1 sum_i d_i v_i >= n^-1 sum_i sqrt(v_i)-sqrt(eta) >=1/2`.  (7)

For the finite matrix realization of `F=F_(j,t)` and `M=M_(j,t)`, the
weighted projection theorem and (4),(6),(7) imply

`liminf n^-1 E[M circ Ztilde]^T B F >= a t`,

where `a=kappa b0/2>0` is independent of `j,t,L` in the stated ranges.

Let `clip_R(z)=max(-R,min(z,R))` and define

`D_i=M(X_i) clip_R(Ztilde_i)`.

The loss from clipping is controlled by the full matrix norm and the
local Gaussian independence:

`limsup |n^-1 E[M circ (Ztilde-clip_R(Ztilde))]^T B F|`

` <= L sqrt(C_M t tau_R)`,                                    (8)

where

`tau_R=E[(|N|-R)_+^2] <= E[N^2 1{|N|>R}]`
`                         <= C_0(R+1) exp(-R^2/2)`.

Here `N` is a standard normal. A centered Gaussian of variance at most
one has a no larger clipping second moment. The deterministic variances
may vary with the root, but this common bound is uniform. Fixed local
joint convergence and polynomial moment bounds justify (8) before any
parameter is changed.

## 5. Exact feasible paired perturbation

Set `epsilon=1/(2R)` and use the two cube means

`mu_+=F(X)+S circ H_j(V_j(X))+epsilon D`,

`mu_-=-F(X)+S circ H_j(V_j(X))+epsilon D`.

Off the support of `M`, feasibility follows from
`|F|+H_j<=1`. On its support, (5) gives `H_j=0`, `|F|<=1/4`, and
`epsilon |D|<=1/2`. Thus both means lie in `[-1,1]^n` exactly.

For `H_B(x)=x^T Bx/2`, symmetry gives the exact identity

`[E H_B(mu_+)-E H_B(mu_-)]/(2n)`

` = n^-1 E F^T B[S circ H_j] + epsilon n^-1 E F^T B D`.        (9)

There is no quadratic perturbation remainder in (9). Independent
rounding of each mean to signs preserves the expected hollow quadratic
energy. Therefore the Boolean absolute optimum divided by `n` dominates
the right side of (9), and its full range divided by `n` dominates twice
that expression.

By the finite old hierarchical energy identity (2), (3), and (8), this
gives the asymptotic bound

`J_j-Ct^2+(1/(2R))[a t-L sqrt(C_M t tau_R)]`.                  (10)

## 6. A gain of order 1/log L, with limits in the safe order

For all sufficiently large fixed `L`, choose

`R=3 sqrt(log L)`, `t=a/(8 C R)`.

Then `t` lies in the uniform range of Section 3. Since
`tau_R<=C_0(R+1)L^(-9/2)`, for every sufficiently large `L` one has

`L sqrt(C_M t tau_R) <= a t/2`.

Expression (10) is consequently at least

`J_j-C t^2+a t/(4R) = J_j+a^2/(64 C R^2)`.                  (11)

Now, after `L,R,t` are fixed, take a finite approximation index `j`
large enough that all the uniform conditions hold and

`J_j >= J*-a^2/(128 C R^2)`.

Finally take the matrix size `n` to infinity at this fixed construction.
Equations (9)--(11) prove that every sequence of bounded-operator
signings with normalized operator norm at most `L` satisfies

`liminf Q(A)/n^(3/2) >= J*+c/log L`,                          (12)

where `c=a^2/(1152 C)>0`, for all sufficiently large fixed `L`.
The negligible distinction between `sqrt(n-1)n` and `n^(3/2)` does not
change (12).

The order is: choose `L`, then `R,t`, then one finite smooth Gaussian
construction, and last let `n` grow. Matrix error constants may depend
arbitrarily on the chosen finite construction. The positive constants
in (12) do not, because they came only from a fixed Gaussian covariance
neighborhood and the nonzero coefficient in (1).

## 7. Returning to unrestricted signings

Principal deletion suffices; no random refill is needed. Take an
unrestricted minimizing sequence with `Q(A)<=C n^(3/2)` for fixed `C`.
The all-order upper bound supplies such a `C`. The proved diagonal
Grothendieck regularization gives a principal submatrix of order
`r>=(1-epsilon)n` such that

`||A_sub||op <= (4 K C/epsilon) sqrt(n)`,

where `K` is the fixed real Grothendieck constant used in
`fresh_range_and_spectral_regularization_2026_09_05.md`. For
`0<epsilon<=1/2` and large `n`, its normalized operator norm is at most
`L=C_1/epsilon`, where `C_1` depends only on `C,K`.

The original optimum dominates that of its principal submatrix exactly:
fix the submatrix spins, average the omitted spins independently, and
all omitted energy terms vanish. Taking a maximizing absolute submatrix
energy gives `Q(A)>=Q(A_sub)`.

Apply (12) at the fixed value `L=C_1/epsilon`. Since `r` grows with `n`,

`liminf_n M_n/n^(3/2)`

` >= (1-epsilon)^(3/2) [J*+c/log(C_1/epsilon)]`.

For sufficiently small fixed `epsilon`, the positive term of order
`1/log(1/epsilon)` dominates the multiplicative loss of order `epsilon`.
Thus, under the explicitly verified Gaussian conditions (1),

`liminf_n M_n/n^(3/2) > J*`.                                  (13)

This is a strict lower-bound improvement, not a convergence theorem or
an upper construction. Its numerical increment is not certified by
this argument alone.

## 8. Application to the banked 21-anchor certificate

The exact instance is supplied by
`fresh_finite_anchor_fixed_point_2026_09_05.md`, with
`alpha=361/500`, edge covariance `rho=2027/2500`, and the remaining
20 anchor coefficients recorded in its exact certificate. Its value
`J*` is enclosed by

`0.430658179405528602724053804634711026327173238336325190455581`

and

`0.430658179405528602724053804634711026327173238336325190455840`.

The exact calculation in
`fresh_finite_anchor_cubic_coefficient_2026_09_05.md` verifies

`-0.063914331766535141564927573349386587491619191945279997110215`
` <= b* <=`
`-0.063914331766535141564927573349386587491619191945279997110143`.

It also proves `Var(V|W)>0.07385942754623718`; `Var W` is enclosed
strictly above `0.5297`. Thus every hypothesis in (1) holds for this
specific target, and (13) applies unconditionally to it:

`liminf_n M_n/n^(3/2) > J*`.

In particular the new unmarked construction strictly exceeds the entire
mathematically specified 21-anchor certificate, not merely a chosen
finite approximation to it. No explicit positive decimal increment
above `J*` is asserted. Such an increment would additionally require
bookkeeping the constants `kappa,b0,C,C_M,L`; existence of a positive
increment is already proved by Sections 3--7.

Independent audit record: the root, algebra, and literature agents each
checked the full normalized-gain argument, including the safe order of
limits, exact paired energy factor, clipping estimate, and principal
deletion. The weighted projection theorem and the nonzero cubic closed
formula also received independent derivations and audits.
