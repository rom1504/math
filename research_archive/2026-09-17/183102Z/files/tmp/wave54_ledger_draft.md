### 10.107 Wave 54: orientation pairs, flat box obstructions, and the linear-core feasibility wall

Wave 54 attacked the joint K-profile population, direct cancellation of
low-cross child spikes, and all-seed production of linear-core partners.
Orientation pairing strictly weakens the sufficient K-abundance lemma by
removing every far state whose opposite orientation already enters the local
recurrence band.  The box audit then constructs a competitive full-scale
hidden spike with \(d=\mathsf H_T=0\), showing that generic block algebra
cannot force cancellation.  Finally, an entropy calculation corrects the
Wave 53 one-threshold proposal: away from one tuned surface its required
partner count exceeds the entire Johnson ball.  An arbitrary-seed
positive-threshold star survives, but remains polynomial.

All algebraic identities and asymptotic entropy comparisons below are
**Verified**.  Stored-instance enumerations are **Numerical exact finite
diagnostics**.  The scalable signings are **Verified constructions satisfying
the stated generic bounds** but are not asserted to be exact minimizers.  All
four checkers were rerun independently.  Nothing in this wave proves
convergence.

#### 10.107.1 Orientation pairing isolates the genuinely selector-driven K-profile branch

Retain

~~~math
h=Q(A[S])-p^{3/2}q_n-t,\qquad
e=\sigma y^{\mathsf T}A[S]y,\qquad
g=(1-p_2)e-h.
~~~

Let \(\iota\) flip only the orientation \(\sigma\).  It is a
measure-preserving involution of the uniform selector/local-state law and
preserves \(u=A[T,S]y\), \(G\), \(J\), \(K_H(u)\), and \(b_H\).
Directly,

~~~math
\boxed{
g^\iota=-g-2h.
}
\tag{10.1296}
~~~

Equivalently, for the unoriented pair put
\(a=|y^{\mathsf T}A[S]y|\) and \(d=(1-p_2)a\).  Its two margins are
\(g_\pm=\pm d-h\), so both lie below \(-B\) exactly when
\(h-d>B\).  Thus the only pairs on which a genuinely nonlocal completion
tail is needed are those whose selector excess dominates even the larger
local orientation.

This gives a strictly weaker sufficient abundance theorem than (10.1285).
Fix \(H=\lceil n^{3/4-c}\rceil\), \(T_n=n^{3/2-c}\), a sufficiently
large fixed \(C_B\), and a slowly diverging \(\omega_n\).  For a state with
\(g<0\), put \(r=-g/p_2\) and define

~~~math
\mathcal R_n=
\{g<0,\ r\ge\omega_nT_n:
\ g^\iota\ge-C_BT_n
\ \hbox{or}\ K_H(u)\ge r+b_H\}.
~~~

Then the **Verified conditional convergence theorem** is

~~~math
\boxed{
\nu_m(\mathcal R_n)\ge e^{-C_0H}
\quad\Longrightarrow\quad
q_m\le p^{3/2}q_n+O(T_n).
}
\tag{10.1297}
~~~

Indeed one of the two branches has at least half the saved mass.  On the
first, apply \(\iota\); saved mass in \(g\ge-O(T_n)\) gives (10.1256) by
the same inverse Hanson--Wright argument, without a completion profile.
On the second, the dependence-safe K theorem (10.1276) gives
\(\exp\{-O(H)\}\) completion probability.  For all large \(n\), the
original far orientation already lies below the local band, so only

~~~math
h-(1-p_2)|y^{\mathsf T}A[S]y|>C_BT_n
~~~

must pay the K-profile.  This selector-driven both-orientations-far branch
is genuinely present in stored finite data.

There is also an exact canonical description of profile failure.  The dual
formula (10.1280) is equivalently

~~~math
K_H(u)=\min_{\theta\ge0}
\left\{\frac{H\theta}{2}+\sum_i
\begin{cases}
|u_i|^2/(2\theta),&|u_i|\le\theta,\\
|u_i|-\theta/2,&|u_i|\ge\theta.
\end{cases}\right\}.
~~~

At its positive minimizer, \(I=\{i:|u_i|\ge\theta\}\), \(s=|I|<H\),
and

~~~math
K_H(u)=\lVert u_I\rVert_1+
\sqrt{(H-s)\lVert u_{I^c}\rVert_2^2}.
~~~

Hence an energy-qualified but profile-bad state has a canonical exceptional
set of fewer than \(H\) outside vertices carrying the missing energy.  This
repackages (10.1281)--(10.1282); it removes an arbitrary choice of head but
does not prove saved mass.

Exhaustive **Numerical exact finite diagnostics** over the stored
\(A_6,A_8,A_9,A_{10}\) target-density cases check 53,878 negative states.
With the deliberately finite normalization \(H\le3\), \(t=b_H=0\),
41,640 have both orientations negative and 30,374 of those fail
\(K_H(u)\ge r\).  These counts show that both obstructions are real
pointwise, but have no asymptotic population force.

#### 10.107.2 Generic child spikes need force neither linear nor quadratic `A^2` cancellation

Write the signing in blocks

~~~math
A=\begin{pmatrix}P&B^{\mathsf T}\\B&C\end{pmatrix},
\qquad P=A[S],\quad B=A[T,S],\quad C=A[T],
~~~

and let \(y\) be a child ground.  Put \(u=By\), \(I=\lVert
Py\rVert_2^2\), \(X=\lVert u\rVert_2^2\), and
\(\mathsf H=A^2-(n-1)I_n\).  Direct block multiplication gives

~~~math
\boxed{
d=\mathsf H[T,S]y=BPy+Cu,\qquad
\mathsf H_T=BB^{\mathsf T}+C^2-(n-1)I_k,
}
\tag{10.1298}
~~~

and

~~~math
I+X=m(n-1)+y^{\mathsf T}\mathsf H[S]y,\qquad
\lVert\mathsf H[S]y\rVert_2^2+\lVert d\rVert_2^2
\ge\frac{\{I+X-m(n-1)\}^2}{m}.
\tag{10.1299}
~~~

These are **Verified frozen-curvature identities**.  They expose the gap in
a trace argument: all positive curvature caused by \(I\) may remain in the
inaccessible principal block \(\mathsf H[S]\).  Both principal blocks have
zero trace, but this does not transfer a negative direction to
\(\mathsf H_T\).

There are exact positive certificates when the accessible quantities are
large.  For every \(z\in[-1,1]^T\), product-bias rounding and conditional
expectation give

~~~math
\boxed{
\mathcal V(S,y)\le
\mu+2d^{\mathsf T}z+z^{\mathsf T}\mathsf H_Tz,
\qquad \mu=I+X+k(n-1).}
\tag{10.1300}
~~~

Taking \(z=-\theta\operatorname{sgn}(d)\) yields an \(\ell _1\) linear
gain, less its signed quadratic curvature.  If
\(\mathsf H_Te=\lambda e\), \(\lambda<0\), and \(\lVert e\rVert_2=1\),
the two choices \(z=\pm e/\lVert e\rVert_\infty\) give

~~~math
\mathcal V(S,y)\le\mu-
\frac{|\lambda|}{\lVert e\rVert_\infty^2}
-\frac{2|d^{\mathsf T}e|}{\lVert e\rVert_\infty}.
\tag{10.1301}
~~~

Thus a negative mode must be both large and delocalized.  A trace-zero
statement alone is insufficient.

The hoped-for spike contribution can also cancel exactly.  Orient
\(r_i=\sigma y_i(Py)_i\), let \(b_i=Be_i\),
\(c_i=\sigma y_ib_i\), \(a_i=c_i^{\mathsf T}d\), and
\(h_i=c_i^{\mathsf T}\mathsf H_Tc_i\).  Then

~~~math
\boxed{
\begin{aligned}
a_i
&=kr_i+\sum_{j\ne i}y_iy_jr_j\langle b_i,b_j\rangle
+\sigma y_i\langle b_i,Cu\rangle,\\
\mathcal V(S,y)
&\le\mu+\min_{0\le\theta\le1}
\{-2\theta a_i+\theta^2h_i\},\\
\sum_i r_i a_i
&=\frac12\{\lVert d\rVert_2^2+\lVert BPy\rVert_2^2
-\lVert Cu\rVert_2^2\}.
\end{aligned}}
\tag{10.1302}
~~~

All three identities are **Verified**.  The off-diagonal cross-column
correlations in the first line can erase \(kr_i\); low \(X\) does not
prevent this.

This failure occurs at the full asymptotic scale.  An explicit scalable
family uses a hub-adjoined competitive signing \(D\) of order \(s\),
pair-duplicates it to an order-\(m=2s\) child \(P\), takes
\(B=(W_0,-W_0)\) from \(k\) rows of a Sylvester Hadamard matrix, and takes
an order-\(k\) symmetric Paley conference signing \(C\).  Choosing
\(k=q+1\) for primes \(q\equiv1\pmod4\) and the least power of two
\(s\ge k\) gives a compact fixed-density sequence
\(2/3\le m/n<4/5\).  Direct calculation proves

~~~math
\boxed{
Q(A)=O(n^{3/2}),\quad \lVert A\rVert_{\rm op}=O(n^{3/4}),\quad
X=0,\quad d=0,\quad\mathsf H_T=0,\quad
I=\Theta(n^{5/2}),\quad
\mathcal V(S,y)=\Theta(n^{5/2}).}
\tag{10.1303}
~~~

This is a **Verified complete-signing construction**, but it is not asserted
to be an exact minimizer.  It falsifies every implication based only on the
generic cap, operator, child-ground, and block identities.  Exact finite
minimizers already show the pointwise wall: at \(A_9,m=8\) and
\(A_{10},m=8\), stored child grounds have \(X=d=0\),
\(\mathsf H_T=0\), and nonuniform fields with maximum \(7\).  Conversely,
the best stored \(A_{10},m=6\) box has \(\mathsf H_T=0\) and
\(d=(-8,8,8,-8)\), and the signed-\(d\) certificate is exact:
\(\mathcal V=74-64=10\).  These are **Verified exact finite diagnostics**,
not fixed-density asymptotic counterexamples.

#### 10.107.3 A one-threshold linear-core cluster is at or beyond the Johnson-ball capacity

Take \(m/n\to p\), \(\ell/n\to\alpha\), and \(s/n\to\beta\), with
\(0<\alpha<\beta<p<1\).  For a fixed \(m\)-selector, the number of
distinct partners with intersection at least \(s\) is

~~~math
\mathcal B_{n,m}(s)
=\sum_{j=s}^{m-1}\binom mj\binom{n-m}{m-j}
=\exp\{nV(p,\beta)+o(n)\},
\tag{10.1304}
~~~

where

~~~math
V(p,\beta)=
\begin{cases}
H(p),&\beta\le p^2,\\
pH((p-\beta)/p)+(1-p)H((p-\beta)/(1-p)),&\beta>p^2.
\end{cases}
~~~

The upper limit \(m-1\) removes the self-loop.  If
\(\epsilon_n=\exp\{-o(n)\}\), the required distinct-partner count on the
right side of the one-threshold lemma (10.1294) has exponent

~~~math
E(p,\alpha,\beta)
=pH(\alpha/p)
+(1-\alpha)H((p-\alpha)/(1-\alpha))
-\beta H(\alpha/\beta).
\tag{10.1305}
~~~

Exact differentiation gives

~~~math
\frac{\partial E}{\partial\alpha}
=\log\frac{(p-\alpha)^2}{(1-\alpha)(\beta-\alpha)}.
~~~

It follows that the **Verified sharp feasibility comparison** is

~~~math
\boxed{E(p,\alpha,\beta)\ge V(p,\beta).}
\tag{10.1306}
~~~

For \(\beta\le p^2\) the inequality is strict for every
\(\alpha>0\).  For \(\beta>p^2\), equality occurs only at

~~~math
\boxed{\alpha_*=\frac{\beta-p^2}{1+\beta-2p}.}
\tag{10.1307}
~~~

Thus (10.1294) is exponentially impossible away from this tuned surface,
even if the favorable family is the entire slice.  At equality it requires
the full Johnson-ball rate \(\exp\{nV+o(n)\}\), not merely some exponential
cluster.  This **Falsifies the general one-threshold implementation**, not
the exact weighted histogram identity (10.1292) or saved spectral excess.
In particular the Wave 53 parameters
\((p,\alpha,\beta)=(0.6,0.3,0.45)\) have
\(E=0.607492605894\ldots>V=0.602026382034\ldots\); the exact
\(\epsilon=0\) demand already exceeds the whole distinct Johnson ball at
\(n=160\).

There is nevertheless a new exact all-seed theorem.  Let \(S\) be any
selector, \(y\) any child ground, \(z\) any completion, and \(k=n-m\).
For every exchange \(S-i+j\), choose between \(z\) and \(z^j\) to maximize
the incoming port.  Adjacent caps and the assigned deficits obey

~~~math
|q_{S-i+j}-q_S|\le2(m-1),\qquad
\sum_{i\in S,\ j\notin S}D_{c_{ij}}(S-i+j)
\le4mk(m-1).
\tag{10.1308}
~~~

Hence at least half the ports have deficit at most \(8(m-1)\), and one of
the \(k+1\) centers \(z,z^j\) contains at least

~~~math
\boxed{\frac{mk}{2(k+1)}}
\tag{10.1309}
~~~

distinct favorable adjacent selectors, including no self-loop.  Its row is
at most \(\{\sqrt{R_2(z)}+2\sqrt n\}^2\).  Aggregating an incidence set
\(I_R\), with the at-most-\(k+1\) preimage multiplicity retained, gives at
least

~~~math
\boxed{\frac{mk}{2(k+1)^2}|I_R|}
\tag{10.1310}
~~~

ordered off-diagonal center/base/partner triples at the enlarged cap.
These statements are **Verified**; exhaustive \(A_{10},m=6\) checking
covers 8,320 child-ground/completion incidences.  But they give only
polynomial partners per base.  A greedy constant-weight-code model can
place exponentially many such stars in mutually separated components, so
total seed mass and local-star axioms alone cannot force a tuned full-ball
cluster.  That model is an abstract incidence obstruction, not an
exact-minimizer family.

#### 10.107.4 Edge-sign minimality gives a migration/entropy dichotomy, not field regularity

There is one exact consequence of discrete signing minimality, but it does
not see the subleading spike excess.  For an exact minimizer \(A\), let
\(\Delta_\omega=q_n-\sigma x^{\mathsf T}Ax\) and
\(s_e(\omega)=\sigma a_ex_ix_j\).  Fix any edge block \(P\), \(|P|=h\),
and write \(C_\omega(P)=\{e\in P:s_e(\omega)=-1\}\).  Applying the exact
simultaneous-flip certificate to every \(a\)-subset \(F\subseteq P\)
proves

~~~math
\boxed{
\Delta_\omega\le4a,\qquad
|F\cap C_\omega(P)|\ge\frac a2+\frac{\Delta_\omega}{8}
}
\tag{10.1311}
~~~

for at least one state \(\omega\), depending on \(F\).  Consequently,

~~~math
\binom ha\le
\sum_{\omega:\Delta_\omega\le4a}
\sum_{j\ge\lceil a/2+\Delta_\omega/8\rceil}
\binom{|C_\omega(P)|}{j}
\binom{h-|C_\omega(P)|}{a-j}.
\tag{10.1312}
~~~

Hoeffding's hypergeometric bound gives the **Verified migration/entropy
dichotomy**: for every \(0<\eta<1/2\), either some state with
\(\Delta_\omega\le4a\) is negative on more than
\((1/2-\eta)h\) edges of \(P\), or there are at least
\(\exp(2\eta^2a)\) oriented states with deficit at most \(4a\).

For a child ground, take \(P\) to be the positive internal star at vertex
\(i\).  Its size is exactly

~~~math
|P|=\frac{m-1+r_i}{2}.
\tag{10.1313}
~~~

At \(a=H=n^{3/4-c}\), (10.1311)--(10.1313) connect a fixed child spike to
parent near-ground geometry.  They still do not bound \(r_i\):
\(|P|=\Theta(n)\) even without the subleading spike, absolute state count
is not saved normalized project mass, and the witnessing restriction need
not be a child ground.  Stored \(A_8,A_9,A_{10}\) cases take the migration
branch strongly; the full positive star is certified by one parent state
of deficit \(0,8,4\), respectively.  This is **Verified exact finite
evidence**, not an asymptotic wall.  A useful continuation would need a
baseline-canceling perturbation that retains the excess \(r_i/2\), plus a
common-witness or favorable-restriction theorem.

#### 10.107.5 Updated frontier

Wave 54 leaves the rigorous convergence interval unchanged and proves no
asymptotic restriction estimate.  It changes the exact frontier as follows:

- the direct completion route now has the relaxed orientation-paired
  sufficient lemma (10.1297).  No K-profile is needed when the opposite
  orientation has \(g\ge-O(T_n)\); only the selector-driven region
  \(h-(1-p_2)|e|>C_BT_n\) must pay the exact clipped profile.  Profile
  failure has a canonical exceptional head of fewer than \(H\) outside
  vertices.  The missing theorem is saved mass for this union, or an
  exchange argument showing that canonical heads cannot carry all of the
  both-orientations-far population;

- the low-cross box route cannot obtain \(A^2\) cancellation from generic
  cap, operator, trace, or child-ground structure.  A scalable competitive
  signing has \(I=\Theta(n^{5/2})\) but \(X=d=0\) and
  \(\mathsf H_T=0\), and stored exact minimizers show the same pointwise
  hiding at finite codimension.  A positive result must now be an
  exact-minimizer localized-curvature exclusion, prove the mean
  regularity target (10.1290) by different means, or select another child
  incidence;

- the one-threshold cluster lemma (10.1294) should be retired at arbitrary
  linear parameters.  It is exponentially feasible only on the tuned
  surface (10.1307), where it asks for essentially the full Johnson-ball
  rate.  The large-core spectral route itself remains open through the exact
  weighted histogram (10.1292); it must control several intersection levels
  together rather than replace them by one threshold;

- every child-ground seed does produce a positive-threshold adjacent star
  and the exact aggregate (10.1310), but this supplies only polynomial
  partners per base.  Global seed count, local stars, and unweighted
  simultaneous edge-flip covers admit separated-component obstructions.
  Future overlap work must use signing-specific compatibility, controlled
  row geometry, or a direct weighted-histogram inequality;

- these corrections materially weaken the former leading structured
  implementation and strengthen the direct sufficient lemma, so
  \(STEERING.md\) is refreshed early at this boundary.  Wave 55 should
  compare canonical-head exchange on the both-orientations-far branch,
  exact-minimizer exclusion of hidden frozen curvature, and a genuinely
  multilevel weighted-intersection theorem.  The next mandatory steering
  refresh is Wave 59 and must include a blank-slate abstraction audit.
