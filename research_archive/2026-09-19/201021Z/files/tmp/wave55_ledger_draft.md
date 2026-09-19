### 10.108 Wave 55: joint completion localization and successor structure

Wave 55 began with the orientation-paired clipped-profile route, but a direct
joint Hanson--Wright audit found a decisive obstruction.  Under the uniform
local-state/outside-completion law, the entire completion increment is one
quadratic form with parent-controlled operator norm.  Its genuinely far tail
is too small on the target entropy scale.  Consequently every saved annealed
completion incidence already contains saved local-band mass and proves the
restriction recurrence without completion.  A sharper cross-Gram calculation
separately falsifies the proposed far K-profile abundance lemmas.

Exact-minimizer block attenuation turns every sufficiently large hidden child-
field spike into migration to a subproject-deficit parent state, but also
isolates the still-missing child-fibre/common-witness compatibility theorem.
On the overlap side, nonnegative common-core mixtures are only convex
repackagings, and every one-threshold compression at a linear core is now
falsified, including its formerly tuned entropy surface by a sharp
square-root shortfall.  The full multilevel histogram remains live.

All quadratic-form identities, norm estimates, tail comparisons, replacement
and histogram theorems, and logical implications below are **Verified**.
Stored-instance statements are **Numerical exact finite diagnostics**.  All
four independent checkers were rerun and reproduced their saved outputs.  The
rigorous convergence interval is unchanged.

#### 10.108.1 The full uniform-completion incidence localizes before coefficient profiles

Fix a compact density window, (0<c<1/4), and retain

~~~math
H=\lceil n^{3/4-c}\rceil,\qquad T_n=n^{3/2-c},\qquad
g=(1-p_2)e-h.
~~~

For a selector (S), write (T=S^c), let (P_S) be coordinate
projection, and define

~~~math
C_S=A-P_SAP_S=
\begin{pmatrix}0&A[S,T]\\A[T,S]&A[T]\end{pmatrix}.
~~~

If (a=(\sigma,y)) is a uniform oriented local projective state, (w) is
a uniform outside completion, and (x=(y,w)), then the centered completion
increment in (10.1244) is exactly

~~~math
Z(w)=\sigma x^{\mathsf T}C_Sx.
~~~

The joint local/completion law gives the uniform full projective law for this
globally sign-invariant quadratic form.  Moreover, for an exact minimizer,

~~~math
\boxed{
\operatorname{tr}C_S=0,\qquad
\lVert C_S\rVert_F^2=2m(n-m)+(n-m)(n-m-1)<n^2,\qquad
\lVert C_S\rVert_{\rm op}\le2\sqrt{2q_n}.
}
\tag{10.1314}
~~~

The last inequality uses
(\lVert P_SAP_S\rVert_{\rm op}\le\lVert A\rVert_{\rm op}) and the exact-
minimizer estimate (\lVert A\rVert_{\rm op}^2\le2q_n).  Hanson--Wright
therefore proves, uniformly in (S),

~~~math
\boxed{
\Pr_{a,w}\{Z(w)\le-u\}
\le2\exp\left[-c_0\min\left\{
\frac{u^2}{n^2},\frac{u}{\sqrt{q_n}}
\right\}\right].
}
\tag{10.1315}
~~~

Since (q_n=O(n^{3/2})), at (u=a_nT_n), (a_n\ge1), the two
exponents are respectively
(\Omega(a_n^2n^{1-2c})) and (\Omega(a_nH)).  Their ratio is
(\Omega(a_nn^{1/4-c})\to\infty), so

~~~math
\boxed{
\Pr_{a,w}\{Z(w)\le-a_nT_n\}\le2e^{-c_1a_nH}.
}
\tag{10.1316}
~~~

Let (Z_t) denote the annealed incidence in (10.1232).  The exact normal
form (10.1244) says

~~~math
Z_t=\mathbb E_{S,a}\Pr_w\{Z(w)\le g/p_2\}.
~~~

For every fixed (D>0), split this incidence according to (g).  On
(g<-p_2DT_n), a favorable completion requires (Z(w)<-DT_n), and
(10.1316) gives the **Verified annealed localization inequality**

~~~math
\boxed{
Z_t\le\nu_m\{g\ge-p_2DT_n\}+2e^{-c_1DH}.
}
\tag{10.1317}
~~~

If (Z_t\ge e^{-C_0H}), choose a fixed (D) with
(c_1D>C_0+1).  Then the local-band event in (10.1317) has mass at least
(e^{-C_0H}/2).  Conditioning on one selector and applying the established
inverse Hanson--Wright estimate to its local quadratic form gives
(h=O(T_n)), hence, for (t=O(T_n)),

~~~math
\boxed{
Z_t\ge e^{-O(H)}
\quad\Longrightarrow\quad
\nu_m\{g\ge-O(T_n)\}\ge e^{-O(H)}
\quad\Longrightarrow\quad
q_m\le p^{3/2}q_n+O(T_n).
}
\tag{10.1318}
~~~

Thus the full uniform arbitrary-cut annealed tail is **Retired as an
independent nonlocal implementation**.  Its lower bound remains a logically
valid sufficient statement, but any such bound already pays for the local
recurrence before completion is used.  This does not retire the bare fixed-
cut selector tail (10.795): selecting one global cut destroys the product
local-state law used in (10.1315).  It also does not retire box/coarea
extraction of such a selected low-row cut.

#### 10.108.2 Cross energy sharply falsifies far K-profile abundance

There is a stronger local-state-only estimate for the K branch.  Fix (S),
put (B=A[T,S]), (R=B^{\mathsf T}B), and
(G=\lVert By\rVert_2^2=y^{\mathsf T}Ry).  Uniform projective (y) has
the full Rademacher law for (G), and

~~~math
\boxed{
\mathbb E G=\operatorname{tr}R=m(n-m),\qquad
\lVert R\rVert_{\rm op}\le2q_n,\qquad
\lVert R\rVert_F^2\le2q_nm(n-m).
}
\tag{10.1319}
~~~

Here positivity gives
(\operatorname{tr}(R^2)\le\lVert R\rVert_{\rm op}\operatorname{tr}R).
Hanson--Wright and the exact dual upper bound
(K_H(u)\le\sqrt H\lVert u\rVert_2) now show, for every (a_n\ge1),

~~~math
\boxed{
\nu_m\{g<0, r\ge a_nT_n, K_H(u)\ge r+b_H\}
\le2e^{-c_2a_n^2H}.
}
\tag{10.1320}
~~~

Indeed K success forces
(G\ge r^2/H\ge a_n^2T_n^2/H).  The two centered Hanson--Wright
exponents are then at least constant multiples of
(a_n^4n^{1-2c}) and (a_n^2H), with the first larger by a diverging
factor.

Taking (a_n=\omega_n\to\infty) makes (10.1320)
(e^{-\omega(H)}), uniformly in the selector, orientation, exact minimizer,
and compact density window.  This **Falsifies** the proposed lower bounds
(10.1277) and (10.1285), and makes the K branch of (10.1297) negligible.
Any (e^{-C_0H}) mass in the orientation-relaxed union must come from its
opposite-orientation local-band branch and hence already proves recurrence.

The same estimate at fixed (a_n=D), with (D) chosen after (C_0),
shows more generally that every (e^{-O(H)})-abundant K-success population
can be truncated to (r=O(T_n)).  Thus the exact conditional completion
theorem (10.1276) remains correct, but the K condition is **Retired as any
independent entropy-(H) implementation under the uniform local-state law**.
Canonical clipped heads cannot rescue it: the Euclidean energy necessary
condition already removes the required far population before head exchange.

#### 10.108.3 Exact minimality turns a hidden field spike into witness migration

The low-cross box route gains a minimizer-specific theorem, but not yet the
box witness.  Let (A) be an exact minimizer, (q=Q(A)), and work on the
oriented parent state space with

~~~math
E_A(\omega)=\tau x^{\mathsf T}Ax,\qquad
\Delta_A(\omega)=q-E_A(\omega),\qquad
s_e(\omega)=\tau a_ex_ix_j.
~~~

For an edge set $\mathcal E$, put $e=|\mathcal E|$, let $v$ be the
number of incident vertices, and define

~~~math
S_{\mathcal E}(\omega)=\sum_{e\in\mathcal E}s_e(\omega),\qquad
\eta_{\mathcal E}=\sqrt{8ev\log2}.
~~~

The **Verified attenuated-block replacement theorem** is: if
(0<\lambda\le1), (0<t\le4e), and every parent state satisfies

~~~math
\Delta_A(\omega)\le t
\quad\Longrightarrow\quad
S_{\mathcal E}(\omega)\ge\lambda e,
~~~

then some integral replacement (A'), changing only signs in
(\mathcal E), obeys

~~~math
\boxed{Q(A')\le q-\frac{\lambda t}{2}+\eta_{\mathcal E}.}
\tag{10.1321}
~~~

Indeed first attenuate those coefficients to
(z_e=(1-t/(4e))a_e).  The oriented energy becomes

~~~math
E_z(\omega)=q-\Delta_A(\omega)
-\frac{t}{2e}S_{\mathcal E}(\omega)
\le q-\frac{\lambda t}{2}
\tag{10.1322}
~~~

both on the (t)-near layer and, using
(S_{\mathcal E}\ge-e), off it.  Independently sign-rounding the attenuated
edges has error (2\sum_{e\in\mathcal E}(\zeta_e-z_e)\tau x_ix_j).
There are at most (2^v) oriented edge patterns.  Hoeffding plus a union
bound gives a rounding with uniform additive error at most
(\eta_{\mathcal E}); a limiting finite-choice argument removes an
arbitrarily small strict slack.  Exact minimality therefore proves the
**Verified localized migration theorem**

~~~math
\boxed{
\frac{2\eta_{\mathcal E}}{\lambda}<t\le4e
\quad\Longrightarrow\quad
\exists\omega:\ 
\Delta_A(\omega)\le t,\quad
S_{\mathcal E}(\omega)<\lambda e.
}
\tag{10.1323}
~~~

Apply this to a child ground (y) on (S).  Orient its internal fields as

~~~math
r_i=\sigma y_i(A[S]y)_i\ge0,\qquad I=\sum_{i\in S}r_i^2.
~~~

For any (U\subseteq S), choose (r_i) child-positive incident edges at
each $i\in U$, and take their union $\mathcal E$.  Every edge is counted
at most twice, so, with (R_U=\sum_{i\in U}r_i),

~~~math
\boxed{
\frac{R_U}{2}\le e\le R_U,\qquad v\le m.
}
\tag{10.1324}
~~~

If (r_{(1)}\ge\cdots\ge r_{(m)}), (R_h=\sum_{j\le h}r_{(j)}), and
(H_m=\sum_{j\le m}j^{-1}), the Lorentz-prefix estimate

~~~math
\boxed{
\max_h\frac{R_h}{\sqrt h}\ge\sqrt{\frac I{H_m}}
}
\tag{10.1325}
~~~

follows from
(r_{(h)}\le h^{-1}R_h\le h^{-1/2}\max_jR_j/\sqrt j).
Choose a maximizing prefix in (10.1324).  If
(I\gg n^{9/4-c}), then

~~~math
e\gg n,\qquad
\eta_{\mathcal E}/e=o(1),\qquad
\eta_{\mathcal E}\le\sqrt{8nq_n\log2}=O(n^{5/4})=o(T_n).
~~~

Taking (t=3\eta_{\mathcal E}/\lambda) in (10.1323) gives the
**Verified hidden-spike migration consequence**

~~~math
\boxed{
\Delta_A(\omega)=o(T_n),\qquad
|\{e\in\mathcal E:s_e(\omega)=-1\}|
>\frac{1-\lambda}{2}e.
}
\tag{10.1326}
~~~

This uses exact discrete minimality and removes the positive-star baseline:
the block size is controlled by field excess (r_i).  It does not prove box
cancellation.  The quantities $X,d,\mathsf H_T$ constrain the child word
$y$, whereas (10.1326) may be witnessed by a different parent near-ground
state.  The exact remaining target is therefore aggregate child-fibre/parent-
witness compatibility, not pointwise curvature exclusion.

Stored exact (A_9,A_{10}) hidden spikes with
(I=120) and $X=d=\mathsf H_T=0$ already have escaping parent grounds
which reverse respectively as many as (7/15) and (9/17) of the
extracted edges.  These are **Numerical exact finite diagnostics** showing
that pointwise common-witness compatibility is false.  Conversely, the
order-22 scalable generic hidden-spike template is visibly nonminimal:
replacing its internal hub block by a pair duplication of the stored
(A_8) lowers the enumerated cap from (170) to (118).  This finite
replacement is **Numerical exact** and is not asserted uniformly for the
scalable family.

#### 10.108.4 The full weighted histogram survives a sharp square-root threshold wall

The overlap route was audited without collapsing the intersection histogram
to one level.  On the (m)-slice put

~~~math
d_\ell=\binom m\ell\binom{n-\ell}{m-\ell},\qquad
K_\ell(S,T)=\frac{\binom{|S\cap T|}{\ell}}{d_\ell}.
~~~

Every nonzero binomial polynomial
(Q(j)=\sum_{\ell=0}^{m-1}c_\ell\binom j\ell), (c_\ell\ge0), after
row normalization is exactly a convex common-core mixture:

~~~math
\boxed{
K_Q=\sum_\ell w_\ell K_\ell,\qquad
w_\ell=\frac{c_\ell d_\ell}{\sum_a c_ad_a}.
}
\tag{10.1327}
~~~

Consequently its harmonic eigenvalues, self-loop, and weighted excess mix
linearly:

~~~math
\boxed{
\Lambda_j(Q)=\sum_\ell w_\ell\lambda_j(\ell),\quad
h_Q=\sum_\ell w_\ell h_\ell,\quad
P_Q-\Lambda_2(Q)=\sum_\ell w_\ell
\{P_\ell-\lambda_2(\ell)\}.
}
\tag{10.1328}
~~~

The extraction denominator
(D_\ell=(1-\lambda_1(\ell))+n(\lambda_1(\ell)-\lambda_2(\ell)))
also mixes.  Hence, whenever the mixture numerator is positive,

~~~math
\boxed{
\frac{P_Q-\Lambda_2(Q)}{D_Q}
\le\max_{\ell:w_\ell>0}
\frac{P_\ell-\lambda_2(\ell)}{D_\ell}.
}
\tag{10.1329}
~~~

Thus nonnegative core mixtures cannot outperform their best component; they
can only help identify one.

There is also an exact law-free partner-capacity formula.  For fixed (S),
the shell (j=|S\cap T|) has
(A_j=\binom mj\binom{n-m}{m-j}) points, and every (K_Q(S,T)) is
nondecreasing in (j).  Among (r) distinct partners, the maximum weighted
load is obtained by filling the largest-intersection shells first.  If
(\sum_{j>s}A_j\le r\le\sum_{j\ge s}A_j), it equals

~~~math
\boxed{
\mathcal C_Q(r)=
\sum_{j>s}A_jK_Q(j)
+\left(r-\sum_{j>s}A_j\right)K_Q(s),
}
\tag{10.1330}
~~~

with the self-loop (j=m) omitted.

For (m/n\to p), (\ell/n\to\alpha), at a fixed linear core, the
(K_\ell)-transition intersection (J) satisfies

~~~math
J-\ell\sim\operatorname{Hypergeom}(n-\ell,m-\ell,m-\ell),\qquad
\frac{\mathbb EJ}{n}\to
\beta_*=\alpha+\frac{(p-\alpha)^2}{1-\alpha}.
\tag{10.1331}
~~~

If

~~~math
v_*=pH\!\left(\frac{p-\beta_*}{p}\right)
+(1-p)H\!\left(\frac{p-\beta_*}{1-p}\right),
~~~

uniform Stirling bounds and hypergeometric concentration prove the sharp
exponential capacity boundary

~~~math
\boxed{
r\le e^{n(v_*-\eta)}\Longrightarrow
\mathcal C_\ell(r)=e^{-\Omega_\eta(n)},\qquad
\mathcal C_\ell(e^{n(v_*+\eta)})=1-o(1).
}
\tag{10.1332}
~~~

Since
(\lambda_2(\ell)\to\theta^2>0),
(\theta=\alpha(1-p)/(p(1-\alpha))), any positive fraction of the
required off-diagonal load needs (e^{nv_*+o(n)}) partners for a positive
fraction of the biased base incidences.  This is a necessary capacity
theorem, not yet a signing theorem.

Most decisively, let

~~~math
\mathcal B_{n,m}(s)=\sum_{j=s}^{m-1}A_j.
~~~

Uniformly when ((p,\alpha)) ranges over a compact subset of
(0<\alpha<p<1),

~~~math
\boxed{
\sup_{\ell\le s\le m-1}
\mathcal B_{n,m}(s)\frac{\binom s\ell}{d_\ell}
=O(n^{-1/2}).
}
\tag{10.1333}
~~~

Above a fixed point between the ordinary Johnson mean (p^2n) and
(\beta_*n), shell sizes decrease geometrically, so the left side is at
most a constant times the largest atom of the hypergeometric law in
(10.1331), which is (O(n^{-1/2})).  Below that point the entropy comparison
is uniformly strict.  The formerly tuned equality surface (10.1307) is
exactly (s/n=\beta_*): matching exponential rates still loses the local-
limit factor (\sqrt n).

But

~~~math
\frac{b_\ell-\binom m\ell}{d_\ell}
=\lambda_2(\ell)-h_\ell\longrightarrow\theta^2>0.
~~~

Therefore the one-threshold condition (10.1294) is **Falsified for every
fixed linear parameter choice**, including the tuned surface: away from it
the shortfall is exponential, and on it the shortfall is
(\Theta(\sqrt n)).  The same uniform estimate shows that no nonnegative
mixture supported on compact linear-core ranges repairs one-threshold
compression.

The full multilevel identity (10.1292) remains live.  It sums
(\Theta(\sqrt n)) transition-typical intersection levels, precisely the
mass discarded by replacing all pair weights by a single boundary weight.
A successful theorem must retain that histogram and use exact-minimizer and
controlled-row structure to create the capacity in (10.1332).

#### 10.108.5 Canonical-head exchange has exact drift formulas but cannot evade localization

For completeness, the selected head attack also closes algebraically.  Let
(I\subset T) be the canonical clipped head, (s=|I|<H), and put
(a_i=\operatorname{sgn}(u_i)), (z_I=\sigma a_I).  For an (s)-set
(R\subset S), write (K=S\setminus R) and choose (R) to minimize the
oriented deleted energy

~~~math
D_R=\sigma\{y^{\mathsf T}A[S]y-y_K^{\mathsf T}A[K]y_K\}.
~~~

Set (S'=K\cup I), (T'=(T\setminus I)\cup R), and
(y'=(y_K,z_I)).  Direct expansion gives

~~~math
\boxed{
e'-e=-D_R+2\lVert u_I\rVert_1
-2a_I^{\mathsf T}A[I,R]y_R
+\sigma a_I^{\mathsf T}A[I]a_I.
}
\tag{10.1334}
~~~

A uniform (R) deletes each internal edge with probability
(\alpha_s=s(2m-s-1)/(m(m-1))), so the minimizing rule gives

~~~math
e'-e\ge2\lVert u_I\rVert_1-\alpha_se-(3s^2-s).
\tag{10.1335}
~~~

If (q_S=Q(A[S])), (\delta=q_S-e), and
(L_s=s(2m-s-1)), the common principal block gives

~~~math
\boxed{
h'-h=q_{S'}-q_S,\qquad
|h'-h|\le L_s,\qquad
h'-h\ge(e'-e)-\delta,\qquad
g'-g=(1-p_2)(e'-e)-(h'-h).
}
\tag{10.1336}
~~~

For (U=T\setminus I), the new cross profile is exactly

~~~math
\boxed{
u'_U=u_U-A[U,R]y_R+A[U,I]z_I,\qquad
u'_R=A[R,K]y_K+A[R,I]z_I.
}
\tag{10.1337}
~~~

Thus removing the old head can add (2s) to every retained outside
coordinate and creates an unrestricted profile on (R); no monotone
K-functional comparison follows.  If one retains the same full spin while
only changing selector membership, its project row is invariant and the old
and new completion increments satisfy

~~~math
\boxed{
Z'=Z-(e'-e),\qquad
p_2Z'-g'=p_2Z-g+(h'-h)-(e'-e).
}
\tag{10.1338}
~~~

The final correction is precisely the change in local deficit.  Since
(L_s=\Theta(sm)) can exceed (T_n) when (s) is near (H), and
(\delta) need not be small under the uniform local-state law, head energy
gain does not control selector excess.  More fundamentally, (10.1316) and
(10.1320) bound the far population before a head is selected.  Canonical-head
exchange is therefore **Retired as a rescue of the uniform K route**; a
future exchange theorem would have to prove the local recurrence after
charging or operate under a genuinely selected nonuniform law.

#### 10.108.6 Updated frontier

Wave 55 leaves the rigorous convergence interval unchanged but decisively
changes the strategic frontier:

- the full uniform arbitrary-cut annealed incidence (Z_t) is recurrence-
  local at entropy (H): (10.1317)--(10.1318) show that every
  (e^{-O(H)}) lower bound already contains saved (g\ge-O(T_n)) mass.
  This retires the uniform-completion detour as an independent route, while
  leaving the selected fixed-cut tail (10.795) intact;

- the orientation-paired K route is falsified, not merely weakened.
  Equation (10.1320) makes its proposed far population
  (e^{-\omega(H)}), and fixed-depth K abundance truncates to the local
  band.  Scalar profiles, canonical clipped heads, and exchange of those
  heads cannot evade the Euclidean or full-completion concentration walls;

- the leading concrete implementation is again selected-cut box plus
  spectral excess: establish a project-row box witness at
  (R_0=O(n^{9/4-c})), then at controlled cap find some nontrivial core
  with (P_\ell-\lambda_2(\ell)\ge e^{-O(H)}).  The overlap proof must use
  the complete weighted histogram.  Every linear-core one-threshold
  replacement is now impossible, and nonnegative mixtures cannot beat their
  best component;

- exact minimality now supplies real structure against the box-field spike.
  If (I\gg n^{9/4-c}), (10.1326) forces an (o(T_n))-deficit parent
  witness which reverses a fixed fraction of a localized field-excess block.
  The missing theorem is aggregate child-fibre/parent-witness compatibility:
  charge these migrated states back to low-cross favorable incidences or turn
  their multiplicity into a controlled-row box witness.  Pointwise common-
  witness compatibility is false finitely;

- at any fixed linear core, saved weighted excess requires
  (e^{nv_*(p,\alpha)+o(n)}) partners for a positive fraction of biased
  bases.  This capacity is abstractly feasible only by retaining the
  (\Theta(\sqrt n)) typical intersection levels.  A next attack must use
  signing-specific compatibility and row control to generate that full
  histogram, rather than another threshold or kernel mixture;

- these are decisive falsifications and a change of leading route, so
  `STEERING.md` is refreshed early at this boundary.  Wave 56 should compare
  aggregate migration compatibility, a direct exact-minimizer inequality for
  the full weighted histogram, and a selected-cut construction which bypasses
  the uniform product local-state law.  The next mandatory steering refresh
  is Wave 60 and must include a blank-slate abstraction audit.
