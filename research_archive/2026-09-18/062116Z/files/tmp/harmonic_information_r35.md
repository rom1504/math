# Wave 35 route C: endpoint erasure information, orientation loss, and migration variation

## Status

- **Verified:** the endpoint cost `C_j` in (10.980) is exactly the relative-
  entropy loss caused by deleting coordinate `j`:
  `D(mu_1||nu_beta)-D((mu_1)_(-j)||(nu_beta)_(-j))`.
- **Verified:** after lifting the endpoint to its canonical selector mixture,
  the corresponding component cost is exactly `C_j` plus the conditional
  selector information `I(S;D_j|D_(-j))`.  Thus selector information is the
  cancellation term, not an upper bound that can simply be discarded.
- **Verified:** for a vertex coordinate, components whose selector omits that
  vertex contribute zero.  The resulting inclusion-budget bound loses a
  factor `m+1=Theta(n)` and therefore does not reach the project scale.
- **Verified:** marginalizing the global orientation strictly loses an
  endpoint conditional KL.  An exact `A_6,m=3,beta=(log 2)/2` certificate has
  a nontrivial orientation likelihood ratio, so the orientation bit is not a
  redundant coordinate of `f_beta`; it cannot be deleted without retaining
  or separately bounding that conditional KL.
- **Verified:** the migration covariance in (10.980) is an exact information
  derivative.  Its adverse part is bounded either by the positive variation
  of a marginal KL path or by a weighted backward variation of the edge-
  context masses.  These give precise sufficient bounded-migration lemmas,
  but no such minimizer-specific lemma is proved here.
- **Open:** neither `sum_j C_j=O(n^(1/2-2c))` nor the required bounded
  migration follows from the identities alone.  The restoring comparison and
  adjacent-selector Hellinger estimate remain separate.

All finite calculations are reproduced by
`tmp/harmonic_information_r35_check.py`, using the repository virtual
environment.  Exact rational arithmetic is used for the orientation
certificate; quadrature-based migration decimals are explicitly numerical.

## 1. Endpoint and selector notation

Use the `n`-bit chart of the oriented-cut cube from Waves 33--34: coordinate
`0` is global orientation and coordinates `1,...,n-1` are projective vertex
flips.  Put

\[
 \nu=\nu_\beta,\qquad \mu=\mu_1=\nu f,\qquad
 \mu_t(d)=\frac{\nu(d)e^{t g(d)}}{Z_t},\qquad g=\log f.
\]

The canonical selector lift (10.938) is

\[
 f(d)=\sum_S q(S)h_S(d),\qquad
 \mu_S(d)=\nu(d)h_S(d),\qquad
 P(S,d)=q(S)\mu_S(d)=\mu(d)\pi_d(S).
 \tag{R35.1}
\]

Every `mu_S` is normalized.  In the matched harmonic construction its
restriction to `S` has the child Gibbs marginal, while its conditional law
outside `S`, given the restriction, is the parent Gibbs conditional.  Only
the likelihood identity in (R35.1) is needed below.

## 2. `C_j` is exactly an erasure relative entropy

Let `D_(-j)` denote all cube coordinates except `j`.  The endpoint binary
cost from (10.980) is

\[
 C_j=\mathbb E_{\mu_{-j}}
 D\!\left(\mu(D_j\mid D_{-j})\,\middle\|\,
          \nu(D_j\mid D_{-j})\right).
\]

The ordinary KL chain rule gives the exact identity

\[
 \boxed{
 C_j=D(\mu\Vert\nu)-D(\mu_{-j}\Vert\nu_{-j}).
 }
 \tag{R35.2}
\]

This also independently verifies that the endpoint expression
`sum_e M_e(1)K_e` in (10.980) has the `mu_1` context weight, not the base
context weight.  Summing over the `n` chart coordinates gives

\[
 \boxed{
 \sum_j C_j
 =nD(\mu\Vert\nu)-\sum_jD(\mu_{-j}\Vert\nu_{-j}).
 }
 \tag{R35.3}
\]

Thus `sum C_j` is an erasure-relative-entropy functional.  Formula (R35.3)
does not by itself upper-bound it at the desired scale; doing so would be an
endpoint-specific approximate-tensorization theorem.

## 3. Exact selector-information cancellation

Define the lifted conditional cost

\[
 J_j:=\mathbb E_{P_{S,D_{-j}}}
 D\!\left(P(D_j\mid S,D_{-j})\,\middle\|\,
          \nu(D_j\mid D_{-j})\right).
 \tag{R35.4}
\]

There are two ways to apply the KL chain rule to
`P` relative to `q tensor nu`:

\[
 \begin{aligned}
 D(P\Vert q\otimes\nu)
   &=D(\mu\Vert\nu)+I_P(S;D),\\
 D(P_{S,D_{-j}}\Vert q\otimes\nu_{-j})
   &=D(\mu_{-j}\Vert\nu_{-j})+I_P(S;D_{-j}).
 \end{aligned}
\]

Subtracting and using (R35.2) proves

\[
 \boxed{
 J_j=C_j+I_P(S;D_j\mid D_{-j}),
 \qquad
 C_j=J_j-I_P(S;D_j\mid D_{-j}).
 }
 \tag{R35.5}
\]

Equivalently, with

\[
 J_{S,j}=D(\mu_S\Vert\nu)
          -D((\mu_S)_{-j}\Vert\nu_{-j}),
\]

one has `J_j=sum_S q(S)J_(S,j)`.  This is the requested exact relation to
selector-posterior information.  In particular, bounding `J_j` and throwing
away the posterior information goes in the wrong quantitative direction:
the information is the cancellation which makes the mixture cheaper than
its components.

For a vertex coordinate `i>=1`, `h_S` is invariant under the flip whenever
`i notin S`.  Hence

\[
 \boxed{J_{S,i}=0\quad\text{if }i\notin S.}
 \tag{R35.6}
\]

The orientation coordinate has no analogous omitted-selector class.  For a
fixed `S`, at most
`m+1-1_{\{0\in S\}}` chart coordinates can have nonzero component cost.
Also `0<=J_(S,j)<=D(mu_S||nu)`.  Therefore the full consequence of the bare
inclusion budget is only

\[
 \boxed{
 \sum_jJ_j
 \le\mathbb E_{S\sim q}
 [(m+1-\mathbf1_{\{0\in S\}})D(\mu_S\Vert\nu)]
 \le(m+1)\{D(\mu\Vert\nu)+I_P(S;D)\}.
 }
 \tag{R35.7}
\]

The last equality used the mixture identity
`E_q D(mu_S||nu)=D(mu||nu)+I_P(S;D)`.  At fixed density the factor `m+1` is
linear in `n`.  Even if both endpoint and selector information already had
the target order, (R35.7) would be a factor `n` too large.

The cancellation is substantial in exact finite minimizers.  For
`A_6,m=3,beta=1/2`, exhaustive enumeration gives, numerically,

\[
 \sum_jJ_j=0.1905952906,\qquad
 \sum_j I(S;D_j\mid D_{-j})=0.1853778460,\qquad
 \sum_jC_j=0.005217444607.
 \tag{R35.8}
\]

Thus retaining only the inclusion-supported component costs loses a factor
about `36.5` already here.  This is a finite diagnostic, not an asymptotic
falsifier of a stronger signing-specific cancellation theorem.

There is a useful contextwise interpretation of the obstruction.  For a
vertex `i`, let `Z=1_{\{i\in S\}}`, condition on `D_(-i)=e`, and write
`alpha_e=P(Z=1|e)`.  Because the omitted-selector component has exactly the
base bit conditional,

\[
 \mu(D_i\mid e)=(1-\alpha_e)\nu(D_i\mid e)
                  +\alpha_e Q_e.
 \tag{R35.9}
\]

Consequently the endpoint conditional KL is the divergence of this mixture
from its base, while (R35.5) says that the component divergence also pays the
information revealing `S`.  If `alpha_e` is close to one, the omitted base
component supplies no uniform information comparison; this is exactly the
posterior-inclusion regime forced by vertex crossings in (10.981).  A valid
project-scale charge therefore needs a weighted theorem controlling those
near-total-inclusion contexts, not just the identity `sum_i 1_{i in S}=m`.

## 4. Orientation marginalization loses a real endpoint term

Let `X=D_(-0)` be the quotient which forgets global orientation, and denote
the corresponding marginals by `bar mu,bar nu`.  Applying (R35.2) to
coordinate zero gives

\[
 \boxed{
 D(\mu\Vert\nu)
 =D(\bar\mu\Vert\bar\nu)
  +\mathbb E_{\bar\mu}
    D(\mu(D_0\mid X)\Vert\nu(D_0\mid X))
 =D(\bar\mu\Vert\bar\nu)+C_0.
 }
 \tag{R35.10}
\]

Thus quotienting is valid data processing, but it decreases the very entropy
which must be upper-bounded.  It proves nothing about the missing `C_0`
unless that term is separately controlled.

The loss can be strict for the actual matched harmonic likelihood.  Take the
exact order-six minimizer `A_6`, `m=3`, and `beta=(log 2)/2`.  Let

\[
 x=(1,-1,-1,-1,-1,-1),\qquad d=-xx^{\mathsf T}
\]

with the diagonal set to zero.  Then
`<A_6,d>=10` and `<A_6,-d>=-10`.  Since all external energies are even, the
definition of `U_beta` can be evaluated exactly with powers of two.  The
checker obtains

\[
 U_\beta(d)=\frac{847888}{3590575},\qquad
 U_\beta(-d)=\frac{111448}{963325},\qquad
 \boxed{
 \frac{f_\beta(d)}{f_\beta(-d)}
 =\frac{U_\beta(d)}{U_\beta(-d)}
 =\frac{1165846}{571171}\ne1.
 }
 \tag{R35.11}
\]

Both states have positive base mass.  Their endpoint conditional orientation
odds therefore differ from their base odds by the nonunit factor in
(R35.11), proving `C_0>0` exactly.  The orientation response is also nonzero,
so this edge contributes strictly positive heat-bath energy at every
interior interpolation time.  Marginalization cannot remove it while
retaining either the original entropy identity or the original screened
energy.  The vanishing orientation cost in some one-deletion or symmetric
finite cases is accidental, not an algebraic property of `f_beta`.

## 5. Migration is an exact marginal-information derivative

Fix coordinate `j` and an unordered edge/context `e={x,y}`.  In addition to
the Wave 34 notation, write

\[
 \psi(t)=\log\mathbb E_\nu e^{tg},\qquad
 \ell_e(t)=\log\frac{M_e(t)}{M_e(0)}=A_e(t)-\psi(t),
\]

and define the marginal relative entropy

\[
 D^-_j(t):=D((\mu_t)_{-j}\Vert\nu_{-j})
            =\sum_eM_e(t)\ell_e(t).
\]

Since

\[
 M'_e(t)=M_e(t)\{A'_e(t)-\psi'(t)\},
 \qquad \psi'(t)=\mathbb E_{\mu_t}g,
\]

direct differentiation yields

\[
 \boxed{
 (D^-_j)'(t)=\operatorname{Cov}_{M_{t,j}}(A'_e(t),A_e(t)).
 }
 \tag{R35.12}
\]

Recall `k_e(t)=tA'_e(t)-A_e(t)`.  Hence, with all signs fixed,

\[
 \boxed{
 \operatorname{Cov}_{M_{t,j}}(A'_e,k_e)
 =t\operatorname{Var}_{M_{t,j}}(A'_e)-(D^-_j)'(t).
 }
 \tag{R35.13}
\]

This is the exact information-derivative form of the migration term.  It is
consistent with the conditional chain rule at every `t`.  Put

\[
 C_j(t)=D(\mu_t\Vert\nu)-D^-_j(t).
\]

The full exponential-family identity and conditional variance decomposition
are

\[
 \frac d{dt}D(\mu_t\Vert\nu)=t\operatorname{Var}_{\mu_t}(g),
 \qquad
 \operatorname{Var}_{\mu_t}(g)
 =\mathcal E_{t,j}(g)
  +\operatorname{Var}_{M_{t,j}}(A'_e).
\]

Together with (R35.13), these give

\[
 \boxed{
 C'_j(t)=t\mathcal E_{t,j}(g)
 +\operatorname{Cov}_{M_{t,j}}(A'_e,k_e).
 }
 \tag{R35.14}
\]

Integrating (R35.14), using `C_j(0)=0`, recovers (10.980) exactly.  This
independent derivation checks both the sign and the single factor `t`.

## 6. Two exact sufficient controls for adverse migration

Define

\[
 \mathcal A_j=\int_0^1
 [-\operatorname{Cov}_{M_{t,j}}(A'_e,k_e)]_+\,dt.
\]

Because the first term on the right of (R35.13) is nonnegative,

\[
 \boxed{
 \mathcal A_j
 \le\int_0^1[(D^-_j)'(t)]_+\,dt
 =:\operatorname{Var}^+_{[0,1]}D^-_j.
 }
 \tag{R35.15}
\]

Therefore a precise sufficient information-growth package is

\[
 \sum_j C_j(1)=O(n^{1/2-2c}),\qquad
 \sum_j\operatorname{Var}^+_{[0,1]}D^-_j
 =O(n^{1/2-2c}).
 \tag{R35.16}
\]

It implies both endpoint-cost and adverse-migration bounds in (10.983).
A stronger relative version,

\[
 \sum_j\operatorname{Var}^+D^-_j\le K_{\rm bt}\sum_jC_j(1),
 \tag{R35.17}
\]

for a fixed constant `K_bt` uniform in `n`, the relevant exact minimizer,
the ratio window, and the prescribed temperature, would make the endpoint
cost alone sufficient for migration.

Monotonicity of every `D^-_j(t)` is not by itself enough at the project
scale.  It only changes the left side of the second estimate in (R35.16) to
`sum_jD^-_j(1)`, while (R35.3) says

\[
 \sum_jD^-_j(1)=nD(\mu\Vert\nu)-\sum_jC_j.
\]

That is circular and carries a factor `n` unless a separate, much stronger
marginal endpoint estimate is proved.  What is needed is bounded positive
variation at the target scale, or a charge such as (R35.17), not mere
monotonicity.

There is a more local endpoint-cost criterion.  Since

\[
 \operatorname{Cov}_{M_{t,j}}(A'_e,k_e)
 =\sum_eM'_e(t)k_e(t),
 \qquad 0\le k_e(t)\le k_e(1)=K_e,
\]

one has

\[
 \boxed{
 \mathcal A_j
 \le\sum_eK_e\int_0^1[-M'_e(t)]_+\,dt.
 }
 \tag{R35.18}
\]

Consequently the weighted bounded-backtracking lemma

\[
 \boxed{
 \sum_{j,e}K_e\int_0^1[-M'_e(t)]_+dt
 \le K_{\rm bt}\sum_{j,e}M_e(1)K_e
 =K_{\rm bt}\sum_jC_j
 }
 \tag{R35.19}
\]

with `K_bt` fixed and uniform over the same project family would prove the
migration half of (10.983) from its endpoint-cost half.  A convenient
stronger pointwise version is
`int_0^1[-M'_e]_+dt<=K_bt M_e(1)` for every edge context.  The exact obstruction
is now explicit: a high-`K_e` context may carry substantial intermediate
mass and then lose it before `t=1`, where `C_j` weights it only by `M_e(1)`.

Finite enumeration does not settle (R35.17) or (R35.19).  For
`A_9,m=7,beta=2`, the checker finds numerically

\[
 \sum_jC_j\approx0.05086819103,\quad
 \int_0^1t\mathcal E_t(g)dt\approx0.05233613551,\quad
 \sum_j\mathcal A_j\approx0.002857229663.
 \tag{R35.20}
\]

Thus coefficient-one domination of weighted energy by endpoint cost is
already false, while a finite constant remains plausible on this data.  The
sampled marginal-KL derivatives for the `A_6/A_8/A_9` paths were nonnegative
at the tested Gauss--Legendre nodes, but this is numerical only and does not
prove monotonicity on the whole interval; as explained above, even exact
monotonicity would not close the exponent.

## Frontier

The endpoint part of the harmonic route now has an exact three-level
decomposition:

1. `C_j` is the endpoint erasure KL (R35.2).
2. Its selector lift is `C_j+I(S;D_j|D_(-j))` (R35.5); local inclusion merely
   bounds the larger lifted term with a fatal linear factor.
3. Adverse interpolation is positive marginal-information variation
   (R35.15), or weighted backward edge-context movement (R35.18).

The clean next sufficient lemma is therefore the pair (R35.16), or the
endpoint erasure bound together with the relative backward-migration charge
(R35.19).  Proving either still requires exact-minimizer structure.  The
orientation cost `C_0` must be retained, the restoring constants in (10.983)
remain independent inputs, and the adjacent-selector Hellinger estimate is
unchanged.  No convergence proof or scalable harmonic counterexample is
obtained in this wave.
