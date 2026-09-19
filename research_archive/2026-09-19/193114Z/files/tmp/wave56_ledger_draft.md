### 10.109 Wave 56: a common active-face law and the mesoscopic core window

Wave 56 followed the selected-cut frontier after the uniform-completion
retirement.  Exact-minimizer perturbation can be upgraded from one migrated
witness per edge block to a single near-ground parent law which simultaneously
balances all bounded-size edge blocks.  This is genuine common-law structure,
although orientation cancellation, bare favorability, and project row remain
uncontrolled.  Independently, the common-core phase diagram identifies the
entropy-matched scale `ell~sqrt(nH)` and a sharp second boundary at
`ell~n^(5/6)`: unlike fixed linear cores, a mesoscopic one-threshold
implementation is count-feasible when `c>1/12`.  Finally, arbitrary selected
cut laws admit an exact collision/row extraction theorem, together with a
sharp matching obstruction to selector-wise existential witnesses.

The active-face LP checker, two independent mesoscopic factorial/log-domain
checkers, and the selected-channel checker were rerun independently.  All
analytic identities and finite-law theorems below are **Verified**;
asymptotic statements are uniform on compact fixed-density windows.  Nothing
in this wave proves the box witness, the fixed-cut tail, or convergence.

#### 10.109.1 Exact minimality gives one common near-ground response law

Let `A` be an exact minimizer, `q=Q(A)`, and on oriented parent states put

~~~math
E_A(\omega)=\tau x^{\mathsf T}Ax,\qquad
\Delta(\omega)=q-E_A(\omega),\qquad
s_e(\omega)=\tau a_ex_ix_j.
~~~

Fix an edge universe `P` incident to `v` vertices and an integer budget `r`.
For

~~~math
\mathcal P_r=\{p\in[0,1]^P:\ \sum_ep_e\le r\},
~~~

fractionally attenuate edge `e` by replacing `a_e` with
`(1-2p_e)a_e`, then independently round it to either sign with mean equal to
that coefficient.  For each oriented edge pattern the energy error is a sum
of centered variables bounded by four with total variance at most `16r`.
There are at most `2^v` patterns.  Bernstein and a union bound, with
`L=(v+1)log 2`, therefore give a simultaneous rounding error

~~~math
\boxed{
\eta_r=\sqrt{32r(v+1)\log2}+\frac83(v+1)\log2.
}
\tag{10.1339}
~~~

The rounded integral signing has cap at least `q` by exact minimality.  Hence
the fractional signing has cap at least `q-eta_r`, or, equivalently, for
every `p` in `mathcal P_r`,

~~~math
\min_\omega\left\{
\Delta(\omega)+4\sum_{e\in P}p_es_e(\omega)
\right\}\le\eta_r.
~~~

Finite minimax now gives the **Verified common active-face theorem**.  If
`h_r(z)` denotes the sum of the `r` largest positive coordinates of `z`,
then one law `mu` on parent states, independent of the subsequently chosen
edge block, satisfies

~~~math
\boxed{
\mathbb E_\mu\Delta
+4h_r\big((\mathbb E_\mu s_e)_{e\in P}\big)\le\eta_r.
}
\tag{10.1340}
~~~

This is stronger than selecting a new migrated witness for each block.  It
simultaneously controls every integral or fractional block of budget `r`.

Condition on `G={Delta<=eta_r/epsilon}`.  Markov gives
`mu(G)>=1-epsilon`; since a discarded state contributes at worst `-r` to
any fractional block, (10.1340) gives

~~~math
\boxed{
h_r\big((\mathbb E_\mu[s_e\mid G])_e\big)
\le\frac{\eta_r/4+\epsilon r}{1-\epsilon}.
}
\tag{10.1341}
~~~

When `v=O(n)`, `n<<r<=C n^(3/2)`, and `0<c<1/4`,
`eta_r=O(sqrt(rn)+n)=o(r)` and `eta_r=o(T_n)`.  Choosing
`epsilon=o(1)` but `eta_r/T_n=o(epsilon)` leaves one law supported on
`Delta=o(T_n)` with `h_r((E_mu s_e)_e)=o(r)`.

Apply this with `P` equal to all parent edges and
`r=C_1n^(3/2)>q_n`.  For every selector choose one oriented child ground.
Its nonnegative internal fields satisfy `sum_i r_i=Q(A[S])`; choosing `r_i`
child-positive incident edges at every vertex and taking their union, as in
(10.1324), produces an excess block `E_S` with

~~~math
q_m/2\le Q(A[S])/2\le |E_S|\le Q(A[S])\le q_n<r.
~~~

Thus `|E_S|=Theta(r)` uniformly and (10.1341) implies
`E_mu S_{E_S}=o(|E_S|)`.  Since the block sum lies between
`-|E_S|` and `|E_S|`,

~~~math
\boxed{
\mu\{S_{E_S}\le |E_S|/2\}\ge\frac13-o(1)
}
\tag{10.1342}
~~~

for every selector.  Fubini proves the **Verified common escape corollary**:

~~~math
\boxed{
\text{one parent state with }\Delta=o(T_n)
\text{ escapes }(1/3-o(1))\text{ of the chosen child fibres.}
}
\tag{10.1343}
~~~

Here escape means reversing more than one quarter of the selected
child-positive edges.  This closes the earlier separate-witness defect at
the level of block escape, but not at the level needed for recurrence.  A
one-state response with `Delta=0` and every `s_e=-1` already satisfies the
minimax conclusion: orientation cancellation permits tiny support.  Moreover
one selected block sum controls neither the complete local deficit
`delta_S(d)` nor the project row `x^T A^2x`.  The common escape incidence
therefore has constant selector--state information cost but is not known to
be bare-favorable and retains only the generic `O(n^(5/2))` row cap.  A
successor must turn a saved fraction of these common incidences into
bare-favorable low-row ones, or obtain a two-sided/sector-sensitive
perturbation certificate excluding the opposite-orientation response.

#### 10.109.2 The entropy-matched spectral core is mesoscopic

Let `N=binom(n,m)`, `p=m/n` stay in a compact subinterval of `(1/2,1)`,
put `a=(1-p)/p`, and suppose

~~~math
H=o(\ell),\qquad \ell=o(n).
~~~

The common-core density, first two nonconstant eigenvalues, and extraction
denominator have the **Verified uniform asymptotics**

~~~math
\boxed{
\begin{aligned}
-\log\rho_\ell
&=\ell\log(1/p)+\frac{1-p}{2p}\frac{\ell^2}{n}
+O\!\left(\frac\ell n+\frac{\ell^3}{n^2}+\log n\right),\\
\lambda_1(\ell)&=(1+o(1))a\frac\ell n,\\
\lambda_2(\ell)&=(1+o(1))a^2\frac{\ell^2}{n^2},\\
D_\ell&=(1+o(1))a\ell.
\end{aligned}}
\tag{10.1344}
~~~

The exact denominator identity behind the last line is

~~~math
D_\ell=(1-\lambda_1)
\left\{1+\frac{\ell(n-m)(n-2)}{(m-1)(n-\ell-1)}\right\}.
~~~

Thus `-log rho_ell=omega(H)` removes the direct collision term, while the
spectral extraction denominator costs only a polynomial.

For a fixed base selector, the kernel intersection has the exact law

~~~math
J_\ell-\ell\sim
\operatorname{Hypergeom}(n-\ell,m-\ell,m-\ell).
~~~

Relative to the ordinary Johnson mean, its mean shift and variance are

~~~math
\boxed{
\mathbb E J_\ell-\frac{m^2}{n}
=\frac{\ell(n-m)^2}{n(n-\ell)},\qquad
\operatorname{Var}J_\ell
=\frac{(m-\ell)^2(n-m)^2}
{(n-\ell)^2(n-\ell-1)}.
}
\tag{10.1345}
~~~

Hence the standardized shift is
`(1+o(1))a ell/sqrt(n)`, which diverges in the present regime.

Let `r_ell^*` be the least number of largest-intersection distinct partners
whose complete `K_ell` mass reaches the off-diagonal degree-two baseline
`lambda_2-h_ell`.  Closest-shell rearrangement, the exact hypergeometric law,
and uniform entropy expansion prove the **Verified sharp exponential
capacity**

~~~math
\boxed{
\log\frac{N}{r_\ell^*}
=\frac{a^2\ell^2}{2n}\{1+o(1)\}.
}
\tag{10.1346}
~~~

More precisely, if `barPhi(q_ell)=lambda_2`, then

~~~math
q_\ell=2\sqrt{\log(n/\ell)}
+O\!\left(\frac{\log\log(n/\ell)+1}{\sqrt{\log(n/\ell)}}\right),
~~~

and, with `tau=a ell/sqrt(n)` and the exact mode rate `I_*`,

~~~math
\log\frac N{r_\ell^*}
=I_*+\tau q_\ell
+O(q_\ell^2+(\ell/n)\tau q_\ell+\log n).
\tag{10.1347}
~~~

The quantile correction is superpolynomial but lower order than
`ell^2/n`.  For an actual favorable-family entropy loss `L=Theta(H)`, the
first full-histogram count-feasible scale is therefore

~~~math
\boxed{
\ell_{\rm count}\asymp\sqrt{nH}=n^{7/8-c/2}.
}
\tag{10.1348}
~~~

The statement is deliberately conditional on the actual loss: a lower bound
`|F|/N>=e^{-CH}` does not say the family is not denser.  Rather,
`ell>>sqrt(nH)` makes the worst guaranteed saved count sufficient at the
capacity level.  No signing correlation is proved.

#### 10.109.3 A second boundary reopens one threshold only mesoscopically

Let

~~~math
R_\ell=\max_{\ell\le s<m}
\left\{\sum_{j=s}^{m-1}\binom mj\binom{n-m}{m-j}\right\}
\frac{\binom s\ell}{d_\ell}
~~~

be the largest load certifiable by one overlap threshold even if the whole
slice is available.  Exact shell ratios and the hypergeometric local limit
theorem give

~~~math
\boxed{
R_\ell=(1+o(1))
\frac{p}{\sqrt{2\pi}(1-p)}\frac{\sqrt n}{\ell},
\qquad
\frac{R_\ell}{\lambda_2(\ell)}
=(1+o(1))\frac1{\sqrt{2\pi}a^3}
\frac{n^{5/2}}{\ell^3}.
}
\tag{10.1349}
~~~

Thus one-threshold compression is capacity-impossible above

~~~math
\ell=(2\pi)^{-1/6}\frac p{1-p}n^{5/6}
~~~

at the critical constant.  Combining this with (10.1348), the only
mesoscopic one-threshold feasibility window at loss `Theta(H)` is

~~~math
\boxed{
\sqrt{nH}\lesssim\ell\lesssim n^{5/6}.
}
\tag{10.1350}
~~~

At exponent level it is nonempty exactly for `c>=1/12`, with constants
decisive at equality.  This does not contradict the Wave 55 falsification,
which was uniform on compact **linear** core ranges.  For `c>1/12`, a
genuinely new mesoscopic one-threshold experiment survives.  For `c<1/12`,
every mesoscopic one-threshold choice is impossible and only the complete
histogram remains count-feasible above `sqrt(nH)`.

Even inside the window, the full histogram is stronger.  It reaches the
upper transition quantile `q_ell`, whereas the highest usable rectangle has
standardized height

~~~math
u_{\rm rect}=\sqrt{2\log(R_\ell/\lambda_2)}+o(\sqrt{\log n}).
~~~

It therefore tolerates fewer partners by

~~~math
\exp\{\tau(q_\ell-u_{\rm rect})+o(\tau\sqrt{\log n})\},
\tag{10.1351}
~~~

a superpolynomial but `e^{o(H)}` factor at the matched scale.  Capacity is
only a necessary condition: migration still supplies neither
`N e^{-O(H)}` distinct selectors nor bounded witness multiplicity.

#### 10.109.4 Exact selected-prior extraction isolates the collision obligation

For an oriented global cut `d`, retain its parent deficit `Delta(d)`, row
`R(d)=x^T A^2x`, and local deficit `delta_S(d)`.  Direct substitution into
(10.792) gives

~~~math
\boxed{
\widehat\ell(S,d)=\delta_S(d)-p_2\Delta(d)-B,
\qquad B=(p^{3/2}-p_2)q_n.
}
\tag{10.1352}
~~~

Let `F_t` be this bare favorable event and
`u_t(d)=U_m{S:(S,d) in F_t}`.  For **any** selector-independent global-cut
prior `nu`, put `Z_nu=E_nu u_t`, condition `nu` by the favorable weight, and
write the captured means as `Rbar,Deltabar`.  Markov and one joint averaging
give the **Verified selected-prior extraction**: for all `a,b>1` with
`theta=1-1/a-1/b>0`, some deterministic cut satisfies

~~~math
\boxed{
u_t(d)\ge\theta Z_\nu,qquad
R(d)\le a\overline R,qquad
\Delta(d)\le b\overline\Delta.
}
\tag{10.1353}
~~~

For a ground-face prior, omit the deficit clause; `a=2` preserves coverage
`Z_nu/2` and row `2Rbar`.  Conditioning `U_m tensor nu` on `F_t` gives the
exact cost identity

~~~math
\boxed{
-\log Z_\nu
=I(S;D)+D(P_S\Vert U_m)+D(P_D\Vert\nu).
}
\tag{10.1354}
~~~

Thus a selected prior with `-log Z_nu=O(H)` and captured
`Rbar=O(n^(9/4-c))` proves (10.795), without using the uniform local-state
completion law retired in Wave 55.

There is also a reference-free exact formulation.  For any joint law `P`
supported on `F_t`, put

~~~math
K(P)=I_P(S;D)+D(P_S\Vert U_m).
~~~

KL projection onto each favorable fibre gives

~~~math
\boxed{
K(P)\ge\mathbb E_D\log\frac1{u_t(D)},
}
\tag{10.1355}
~~~

and, for all nonnegative `lambda,eta`, the **Verified variational equality**

~~~math
\boxed{
\inf_{P:\,P(F_t)=1}
\{K(P)+\lambda\mathbb ER(D)+\eta\mathbb E\Delta(D)\}
=\min_{d:u_t(d)>0}
\left\{\log\frac1{u_t(d)}+\lambda R(d)+\eta\Delta(d)\right\}.
}
\tag{10.1356}
~~~

Equality uses a constant output and the uniform selector law on its favorable
fibre.  In particular, `K(P)<=K`, `ER<=R_0`, and `EDelta<=D_0` give one cut
with

~~~math
u_t(d)\ge e^{-3K},\qquad R(d)\le3R_0,qquad\Delta(d)\le3D_0.
~~~

The new obligation is exactly a low-information favorable collision, not
selector-wise witness existence.  In an abstract matching with one distinct
active row-`R_0` cut per selector, every selector has a perfect witness but

~~~math
\boxed{I(S;D)+D(P_S\Vert U_m)=\log\binom nm}
\tag{10.1357}
~~~

for every supported law, and every output fibre has mass `e^{-Theta(n)}`.
This is an incidence-level obstruction, not a signing construction.

Even the common near-ground minimax law from this wave gives `O(1)`
information only for its **escape** events.  Escape controls one selected
internal edge sum; it neither implies bare favorability in (10.1352) nor
controls captured row below the generic `O(n^(5/2))`.  These are the two
exact remaining bridges to the selected-prior package.

#### Updated frontier

1. **The migration multiplicity gap has narrowed to two exact bridges.**
   Exact minimality now produces one selector-independent law on an
   `o(T_n)` parent face which escapes a constant fraction of every chosen
   child excess block.  What remains is not common-witness existence: it is
   to prove that a saved fraction is bare-favorable and to reduce its captured
   project row from the generic `O(n^(5/2))` to `O(n^(9/4-c))`.  The
   opposite-orientation one-state response is the sharp obstruction to any
   marginal-only argument.

2. **The leading spectral scale is now mesoscopic.**  At actual entropy loss
   `Theta(H)`, the complete histogram becomes count-feasible at
   `ell~sqrt(nH)`.  For `c>1/12`, even a one-threshold implementation has a
   previously missed capacity window
   `sqrt(nH) lesssim ell lesssim n^(5/6)`; for `c<1/12` it is impossible and
   only the full histogram survives.  This reopens one-threshold work only at
   the mesoscopic scale and does not undo the linear-core falsification.

3. **The selected-prior obligation is exact.**  It suffices to construct a
   selector-independent cut prior with favorable mass `e^{-O(H)}` and
   favorable-captured row `O(n^(9/4-c))`.  Equivalently one needs a
   favorable joint law of information cost `O(H)` and that row cost.
   Selector-wise active witnesses alone can have information cost `Theta(n)`.

4. **The box remains open.**  The common active-face law is the strongest
   exact-minimality consequence so far, but escape does not yet imply box
   cancellation.  Wave 57 should test: (i) sector-sensitive or two-sided
   perturbations converting common escape to bare favorability; (ii) a
   mesoscopic threshold cluster built from migration/common-law structure;
   and (iii) direct captured-row control for the selected prior.

Because (10.1350) changes the viable implementation and leading scale,
`STEERING.md` receives an early refresh at this boundary.  The next mandatory
refresh is Wave 61 and must include a blank-slate abstraction audit.
