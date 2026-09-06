### 10.111 Wave 58: relative-fibre transfer and scalar migration walls

This wave tested the three mechanism-level targets from Wave 57.  The main
positive result is an exact concentration theorem for transferring the bare
favourable fibre between two oriented cuts.  Its controlling statistic is the
row of their relative cut/uncut mask, and it works at precisely the project
entropy speed.  It also shows why the present high-row migration theorem does
not suffice: at a genuinely high endpoint the exponent loses the entire row
factor.  The proposed two-column construction collapses, up to harmless
constant slack, to the existing one-column target.  Exact finite minimizers
separately block migration-only scalar descent and monotone one-spin descent.
No convergence or new restriction recurrence is proved.

All identities and LP statements below are **Verified**.  The relative
formula and variance were checked independently by exhaustive fixed-slice
enumeration.  The `A_8` statements are **Exact finite mechanism walls**, not
asymptotic counterexamples.

#### 10.111.1 Bare fibres transfer through the relative masked row

Let `d=(tau,x)` and `omega=(upsilon,z)` be oriented global cuts.  Gauge at
`d`, put `kappa=tau upsilon`, `y_i=x_i z_i`, and define

~~~math
s_{ij}=\tau a_{ij}x_ix_j,\qquad
\mathcal D=\{ij:\kappa y_iy_j=-1\},\qquad
b_i=\sum_{j:ij\in\mathcal D}s_{ij},\qquad
W=\sum_{ij\in\mathcal D}s_{ij}.
~~~

The set `mathcal D` is a cut when `kappa=1` and an uncut when `kappa=-1`.
Writing `r_i(d)=sum_j s_(ij)`, direct switching gives

~~~math
\boxed{
\begin{aligned}
\Delta(\omega)-\Delta(d)&=4W,\\
r(\omega)&=r(d)-2b,\\
R(\omega)-R(d)&=4\lVert b\rVert_2^2-4\langle r(d),b\rangle.
\end{aligned}}
\tag{10.1374}
~~~

The row identity is the orientation-unified form of (10.914), not a new row
mechanism.  The new object comes from subtracting the full bare losses.  For
`S sim U_m`, set `W_S=sum_(ij in mathcal D, i,j in S)s_(ij)`.  Cancellation
of the principal norm and use of the retained parent deficit give

~~~math
\boxed{
\widehat\ell(S,\omega)-\widehat\ell(S,d)
=4\{W_S-p_2W\}.}
\tag{10.1375}
~~~

This variable is exactly centred.  Equal, adjacent, and disjoint edge-pair
counting yields

~~~math
\boxed{
\begin{aligned}
\operatorname{Var}(W_S-p_2W)
={}&(p_2-2p_3+p_4)|\mathcal D|\\
&+(p_3-p_4)\lVert b\rVert_2^2
+(p_4-p_2^2)W^2.
\end{aligned}}
\tag{10.1376}
~~~

The last coefficient is nonpositive.  More sharply, the masked matrix
`G_(ij)=s_(ij)1_(ij in mathcal D)` satisfies

~~~math
G=\tfrac12(S_d-\kappa\operatorname{diag}(y)S_d
\operatorname{diag}(y)),\qquad
G\mathbf1=b,\qquad \lVert G\rVert_{\rm op}\le\lVert A\rVert_{\rm op}.
~~~

The Bernoulli-conditioning/Hanson--Wright proof of (10.1124), applied to
`G`, therefore proves

~~~math
\boxed{
\Pr\{|\widehat\ell(S,\omega)-\widehat\ell(S,d)|\ge a\}
\le C\sqrt n\exp\left[-c\min\left\{
\frac{(a-Cq_n/n)_+^2}{n^2+\lVert b\rVert_2^2},
\frac{(a-Cq_n/n)_+}{\lVert A\rVert_{\rm op}}
\right\}\right].}
\tag{10.1377}
~~~

Here the slice-centering correction is `O(q_n/n)` because the total energy
of `G` is half an energy difference of two cuts.  Consequently

~~~math
u_{t+a}(\omega)\ge u_t(d)-\varepsilon_a(d,\omega),
\tag{10.1378}
~~~

where `epsilon_a` is the right side of (10.1377).  With

~~~math
H=n^{3/4-c},\qquad T_n=n^{3/2-c},\qquad R_*=n^{9/4-c},
~~~

the exact project-scale consequence is

~~~math
\boxed{
\lVert b\rVert_2^2\le C_bR_*,\quad u_t(d)\ge e^{-C_0H}
\quad\Longrightarrow\quad
u_{t+KT_n}(\omega)\ge\tfrac12u_t(d)
}
\tag{10.1379}
~~~

for a sufficiently large fixed `K=K(C_0,C_b)`.  Indeed the two exponents
are `K^2H/C_b` and `KH`.  Bounds `R(d),R(omega)=O(R_*)` suffice since
`||b||<=(sqrt(R(d))+sqrt(R(omega)))/2`, but the exact hypothesis is the
relative row.

This also quantifies the obstruction.  At the generic endpoint scale
`||b||^2=O(n^(5/2))`, the Gaussian exponent is only

~~~math
\boxed{T_n^2/n^{5/2}=n^{1/2-2c}=o(H).}
\tag{10.1380}
~~~

Thus ordinary concentration cannot preserve an `e^{-Theta(H)}` fibre during
a genuine high-to-low row migration.  A response with controlled endpoint
row, a sharper relative-row susceptibility, or a multistage entropy--row
tradeoff is still required.  Widening by `O(T_n)` is also essential: at a
fixed threshold a fibre may lie on its boundary.

#### 10.111.2 Two-column optimization is not a separate asymptotic route

Put `a_d=e^(h_d)=1/u_d`.  The selected-prior LP and its exact scalar dual are

~~~math
\boxed{
Z_{\rm sel}(R_0)^{-1}
=\min_{\pi:\,\mathbb E_\pi R\le R_0}\mathbb E_\pi a_D
=\sup_{\gamma\ge0}
\left[\min_d\{a_d+\gamma R_d\}-\gamma R_0\right].}
\tag{10.1381}
~~~

Its supporting hull is therefore `(R,e^h)`, whereas the information--row LP
(10.1368) uses `(R,h)`.  These scalarizations cannot be interchanged without
a spread bound.  Nevertheless, for
`U_*(L)=max_(d:R_d<=L)u_d`, Markov and conditional averaging give, for every
fixed `b>1`,

~~~math
\boxed{
U_*(R_0)\le Z_{\rm sel}(R_0)
\le\frac b{b-1}U_*(bR_0).}
\tag{10.1382}
~~~

This is the row-only specialization of (10.1353).  Hence any two-column
selected prior at incidence `e^{-O(H)}` and row `O(R_*)` already rounds to
one cut with the same exponent and row scale.  Strict caps can benefit from
mixing, but constructing the optimal pair is not an independent asymptotic
proof obligation.

The retained-deficit switching identities do not create local descent.  On
the exact `A_8,m=5,t=0` bare incidence, one column has

~~~math
(|F|,R,\Delta)=(32,72,20).
~~~

Its orientation mate is identical and all eight one-spin neighbours lower
row to `64` while remaining nonempty: four have `(|F|,Delta)=(19,8)` and
four have `(15,32)`.  At price `1/20` the best-neighbour increments are

~~~math
\boxed{
\Delta(h+R/20)=\log(32/19)-2/5>0,
\qquad
\Delta(e^h+R/20)=56/19-7/4-2/5>0.}
\tag{10.1383}
~~~

Thus it is a strict local minimum modulo orientation for both scalarizations,
although a two-spin block does descend and a distant low-row column is much
better.  This is an **Exact finite wall to monotone orientation/one-spin
descent**, not to global block switching or an orbit-Harnack theorem.

There is a complementary migration wall on exact `A_8,m=6`.  A correct-sector
column globally minimizes `h+10^(-3)R` and has

~~~math
(u,R,\Delta)=(18/28,72,4).
~~~

A positive two-edge block has an exact-ground response reversing half its
edges and lowering row to `64`, but that response has `u=11/28` and

~~~math
\boxed{\Delta(h+10^{-3}R)=\log(18/11)-0.008>0.}
\tag{10.1384}
~~~

Another escaping response keeps row `72`.  This grants the qualitative
migration conclusion directly but shows that reversed block mass alone
controls neither the complete signed relative row nor fibre surprise.  The
block is too small for the asymptotic rounding premise of (10.1323), so the
example is a scoped algebraic wall rather than an asymptotic falsifier.

#### 10.111.3 Local child-sector consistency still permits Johnson-random cancellation

The Wave 57 balanced-partition wall used unrelated center labels.  A stronger
local obstruction survives after the common-state Fubini step.  If one common
parent label escapes on every selector but receives the required signed
cancellation exactly on a family `F subset Omega=binom([n],m)` of size `r`,
then a uniformly random choice of `F` has

~~~math
\boxed{
\mathbb E P_\ell(F)=h_\ell+(1-h_\ell)\frac{r-1}{N-1},
\qquad
\mathbb E\Pi_s(F)=\frac{r-1}{N-1}.}
\tag{10.1385}
~~~

Thus at `r/N=e^(-Theta(H))` some deterministic cancellation family has full
histogram value `e^(-Theta(H))+e^(-Theta(n))`, far below the polynomial
`lambda_2`, and its optimal-threshold retention is exponentially below the
baseline in (10.1364).  This does not obstruct the direct one-column theorem
if the common parent already has controlled row; it says only that the
mesoscopic spectral certificate is not forced by common escape and arbitrary
cancellation labels.

The labels can preserve more signing structure than arbitrary incidence.
There is an exact order-six child signing `C` with cap `10`, exposed all-one
ground, fields `(5,1,1,1,1,1)`, and legal five-edge positive star `E`.  Two
actual relative cuts have

~~~math
\boxed{
(\sum_{e\in D}c_e,\ S_E,\ \delta)=(0,1,0),
\qquad (1,-1,4).}
\tag{10.1386}
~~~

Both escape the same block, while any bare local threshold in `[0,4)` accepts
only the first.  Assigning the two records according to an arbitrary family
`F` therefore preserves exact child optimality, a genuine cut sector, the
identities (10.1358), legality of the Wave 56 block, and a constant child cap.
At `(n,m,ell)=(11,6,3)` an explicit eight-selector family gives

~~~math
\boxed{
P_3=\frac3{112}<\lambda_2(3)=\frac1{14},
\qquad
\Pi_4=\frac1{120}<\frac1{12}
=\frac{\lambda_2-h_3}{R_{3,4}}.}
\tag{10.1387}
~~~

This is a **Local/abstract sector-compatible wall**.  Its selector-wise
child patterns do not glue to restrictions of one global signing, and it
controls neither parent row nor global deficit nor exact-minimizer status.
It therefore leaves precisely a global restriction-consistency theorem open.

Elementary adjacent-selector exchange is too local to supply that theorem.
For adjacent selectors `S,T`, extending a child ground across the exchanged
vertex and comparing a fixed parent's induced energies gives

~~~math
\boxed{|Q(A[S])-Q(A[T])|\le2(m-1),\qquad
|\delta_S(d)-\delta_T(d)|\le6(m-1).}
\tag{10.1388}
~~~

Target slack transports this bound through only `O(T_n/n)=O(n^(1/2-c))`
exchanges, whereas a mesoscopic `K_ell` partner at `ell=o(n)` is at Johnson
distance `Theta(n)`.  Any positive codegree result must use overlapping-edge
consistency, parent row/global energy, or many-selector exact minimality, not
scalar exchange Lipschitzness.

#### Updated frontier

1. **The selected low-row bare column remains the leader, and two-column
   optimization is no longer a separate target.**  The exact selected-prior
   hull lives in `(R,e^h)`, but constant-slack rounding (10.1382) turns every
   successful mixture into the one column already required by (10.795).

2. **Bare-fibre stability now has an exact diagnostic.**  The relative
   cut/uncut mask controls the entire change through (10.1375), and relative
   row `O(R_*)` gives entropy-speed transfer after `O(T_n)` threshold
   enlargement.  A genuinely high-to-low migration loses the factor
   `R/R_*` in its exponent.  The next row theorem must control the complete
   signed relative row, exploit an entropy--row tradeoff over global block
   orbits, or produce the low-row column directly.

3. **Local escape/cancellation data do not force Johnson reuse.**  The new
   template preserves child grounds, legal excess blocks, genuine sectors,
   and the exact deficit identity, while allowing a Johnson-random
   cancellation family.  It does not preserve global restriction glue.
   Hence a surviving overlap theorem must be genuinely global and
   minimizer-specific; another local star, degree, or exchange argument is
   insufficient.

4. **Wave 59 should test global mechanisms.**  The strongest independent
   targets are: (i) a global block-orbit or entropy--row inequality beyond
   monotone one-spin descent; (ii) a many-selector consistency theorem that
   couples the active-face cancellation labels; and (iii) a direct
   low-row/cancellation consequence of exact minimality, bypassing
   mesoscopic spectral reuse if possible.

No decisive leading-route change occurs, so `STEERING.md` is not refreshed
at this boundary.  The next mandatory refresh remains Wave 61 and must
include a blank-slate abstraction audit.
