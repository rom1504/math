# Exact and approximate response states for heterogeneous binary mean field

**Date:** 2026-08-16  
**Status:** independent proof-level benchmark report. Only this requested
temporary report is created.

## 1. Main result and scope

Let a block \(A\) have \(n\) binary sites,

\[
H_A(x)=\sum_{i\in A}h_i x_i,\qquad
x_i\in\{0,1\},\quad h_i\in[-B,B],                         \tag{1.1}
\]

where \(B>0\). Write the fields in decreasing order
\(a_1\ge\cdots\ge a_n\), and define the fixed-occupancy optimum

\[
p_A(k)=\max_{\sum_i x_i=k}H_A(x)=\sum_{j=1}^k a_j,
\qquad 0\le k\le n.                                      \tag{1.2}
\]

The candidate survives, with one important qualification.

* In the local-field model, \(p_A\) is discrete concave. A single uniform
  terminal field \(\lambda\sum_i x_i\) recovers every coordinate; arbitrary
  lookup potentials \(g(k)\) are unnecessary. Thus \(p_A\), equivalently its
  slope multiset \(\{h_i\}\), is the coarsest exact state.
* Disjoint blocks compose by max-plus convolution of their profiles. In
  slope coordinates this is simply sorted multiset union.
* Quantizing slopes on a common mesh \(\eta\) produces an exactly mergeable
  histogram with

  \[
  \log_2|\mathcal S_{n,\eta}|
  =O\!\left((1+B/\eta)\log(n+1)\right),               \tag{1.3}
  \]

  and uniform additive response error at most \(\eta n/2\). In the intended
  regime \(0<\eta\le B\), this is
  \(O((B/\eta)\log n)\) bits. Repeated merging charges each site once, so the
  root error is proportional to total mass, not tree depth.
* If a known aggregate term destroys concavity, linear fields recover only
  the least concave majorant (the response roof). For a fixed uniform
  quadratic interaction the roof is still an exact composable quotient,
  because its cross term is bilinear. It can be strictly smaller than the
  raw conditional profile.

The histogram rate is a robust constructive upper bound, not a proved
minimax rate at error \(\Theta(\eta n)\). Exact grid histograms have a
matching binomial state-count lower bound, but a general macroscopic packing
proved below gives only

\[
\log_2 K_{c\eta n}
\ge c'\min\{n,\sqrt{B/\eta}\}.                       \tag{1.4}
\]

The gap is real in the present proof and is stated rather than hidden.

## 2. The non-tautological context family

A system is a finite multiset \(A=\{h_i\}\subset[-B,B]\), composed by
disjoint union. The **local aggregate contexts** are:

1. append an arbitrary future block \(C\) from the same class; and
2. apply one scalar chemical potential \(\lambda\in\mathbb R\) to the total
   occupancy \(K=\sum_i x_i\).

The response is

\[
\mathcal R_A(C,\lambda)
=\max_{x_A,x_C}
 \{H_A(x_A)+H_C(x_C)+\lambda(K_A+K_C)\}.             \tag{2.1}
\]

This family contains no freely specified table
\(g:\{0,\ldots,n\}\to\mathbb R\). Such a table could pin a fibre by
definition and would make profile recovery tautological. Here the terminal
query has only one real parameter. Because the local energy is additive,

\[
\mathcal R_A(C,\lambda)=R_A(\lambda)+R_C(\lambda),\qquad
R_A(\lambda):=\max_x\{H_A(x)+\lambda K_A\}.           \tag{2.2}
\]

The empty future is allowed. Hence contextual equivalence is exactly
equality of the one-variable functions \(R_A\).

Two conventions matter.

* After merging, only total occupancy is visible. If a future can still
  apply separate fields \(\lambda_1K_1+\lambda_2K_2\) to named old blocks,
  an anonymous union histogram is insufficient; the state must retain the
  tuple of labelled histograms or an equivalent joint profile.
* Equation (1.1) anchors the all-zero configuration at energy zero. For
  general local energies \(e_i(0),e_i(1)\), retain the additive baseline
  \(\sum_i e_i(0)\) and apply the theorem to increments
  \(h_i=e_i(1)-e_i(0)\). If responses are projective (modulo constants), the
  baseline is instead quotiented out.

For Ising spins, set \(x_i=(1+\sigma_i)/2\). Then
\(g_i\sigma_i=2g_i x_i-g_i\) and magnetization is \(M=2K-n\).
Thus the occupancy and magnetization versions differ only by an affine
coordinate change and an additive baseline; mass must be retained.

## 3. Exact theorem for local fields

### Theorem 3.1 (minimal state, metric, and merge law)

For the contexts in Section 2:

1. \(p_A\) is discrete concave, and its slopes are the decreasingly sorted
   local fields.
2. Uniform linear fields determine all values of \(p_A\):

   \[
   R_A(\lambda)
   =\sum_{i\in A}(h_i+\lambda)_+
   =\max_{0\le k\le n}\{p_A(k)+\lambda k\},             \tag{3.1}
   \]

   \[
   p_A(k)=\inf_{\lambda\in\mathbb R}
          \{R_A(\lambda)-\lambda k\}.                   \tag{3.2}
   \]

   At known mass and with \(h_i\in[-B,B]\), the restricted range
   \(\lambda\in[-B,B]\) already suffices.
3. For two blocks of the same mass,

   \[
   \sup_{\lambda\in\mathbb R}|R_A(\lambda)-R_{A'}(\lambda)|
   =\max_{0\le k\le n}|p_A(k)-p_{A'}(k)|.               \tag{3.3}
   \]

4. Disjoint union obeys

   \[
   p_{A\sqcup C}(t)
   =\max_{k+\ell=t}\{p_A(k)+p_C(\ell)\}.                \tag{3.4}
   \]

   The slopes on the right are the sorted multiset union of the two child
   slope multisets.
5. Up to injective recoding, \(p_A\), its slope multiset, and \(R_A\) are
   the coarsest exact deterministic state closed under these contexts.

#### Proof

For fixed \(k\), an optimizer must use the \(k\) largest fields: exchanging
an included field for a larger omitted one cannot decrease the energy and
strictly improves it when the inequality is strict. This proves (1.2), and

\[
p_A(k)-p_A(k-1)=a_k.                                \tag{3.5}
\]

The increments are nonincreasing, so \(p_A\) is discrete concave.

With a terminal field, each site is independently occupied precisely when
\(h_i+\lambda>0\), with either choice allowed at equality. This proves
(3.1). For \(1\le k<n\), every

\[
-a_k\le\lambda\le-a_{k+1}                           \tag{3.6}
\]

makes \(k\) a maximizer in (3.1). For \(k=0\), take
\(\lambda\le-a_1\); for \(k=n\), take \(\lambda\ge-a_n\).
Thus \(R_A(\lambda)-\lambda k\ge p_A(k)\) for all \(\lambda\), with equality
for one of these choices, proving (3.2). Ties may prevent \(k\) from being a
unique optimizer, but its point lies on the same supporting line and is
still recovered. All choices can be made in \([-B,B]\).

Let \(\delta=\max_k|p_A(k)-p_{A'}(k)|\). Taking maxima in (3.1) gives response
distance at most \(\delta\). Conversely, (3.2) applied to two response
functions at sup distance \(d\) gives
\(|p_A(k)-p_{A'}(k)|\le d\) for every \(k\). This proves (3.3).

A \(t\)-site subset of \(A\sqcup C\) uses \(k\) sites from \(A\) and
\(t-k\) from \(C\), and the two choices optimize independently once \(k\)
is fixed. This proves (3.4). Equivalently, the best \(t\) fields in the
union are the \(t\) largest elements of the union multiset, proving the
slope-union assertion. Multiset union is associative.

Finally, any exact state must answer the empty-future queries and hence
determine \(R_A\). Formula (3.2) then forces it to determine \(p_A\), and
(3.5) forces it to determine the field multiset. Conversely, that multiset
answers all contexts and composes by union. This proves minimality and
closure. \(\square\)

The point is stronger than “a conditional table is sufficient”: concavity
identifies a one-parameter family that exposes it, and converts tropical
convolution into ordinary histogram addition.

## 4. Quantized histogram theorem

Let

\[
M=1+\left\lceil\frac{2B}{\eta}\right\rceil,\qquad
\Delta=\frac{2B}{M-1}\le\eta,                       \tag{4.1}
\]

and use the \(M\)-point equally spaced grid from \(-B\) to \(B\). Let \(Q(h)\)
be a deterministic nearest grid point, so
\(|Q(h)-h|\le\Delta/2\le\eta/2\). Define

\[
S_\eta(A)=(c_1,\ldots,c_M),\qquad
c_j=|\{i:Q(h_i)=\gamma_j\}|.                         \tag{4.2}
\]

(For \(B=0\), the one-bin state is trivial.)

### Theorem 4.1 (depth-uniform response approximation)

For every block of mass \(n\):

1. \(S_\eta(A\sqcup C)=S_\eta(A)+S_\eta(C)\)
   coordinatewise.
2. At fixed mass the exact state count is

   \[
   |\mathcal S_{n,\eta}|=\binom{n+M-1}{M-1}
   \le(n+1)^M.                                        \tag{4.3}
   \]

3. If \(\widetilde p_A,\widetilde R_A\) are decoded from the grid multiset,

   \[
   |p_A(k)-\widetilde p_A(k)|
   \le\frac{\eta k}{2}\le\frac{\eta n}{2},            \tag{4.4}
   \]

   \[
   \sup_\lambda|R_A(\lambda)-\widetilde R_A(\lambda)|
   \le\frac{\eta n}{2}.                               \tag{4.5}
   \]

4. On any binary merge tree with total leaf mass \(N\), the root response
   error is at most \(\eta N/2\), independently of depth and bracketing.

The same error bound holds after adding any fixed known aggregate term that
is evaluated identically in the exact and decoded systems, including a
uniform quadratic mean-field term and terminal linear fields.

#### Proof

Part 1 holds because all leaves use the same sitewise map \(Q\), and merging
only adds counts. A mass-\(n\) histogram is a weak composition of \(n\) into
\(M\) parts, proving the equality in (4.3); the displayed upper bound follows
because every count lies in \(\{0,\ldots,n\}\).

For every \(k\)-site set \(T\),

\[
\left|\sum_{i\in T}h_i-\sum_{i\in T}Q(h_i)\right|
\le\frac{\eta k}{2}.                               \tag{4.6}
\]

Taking the maximum over all such sets in the two models proves (4.4).
Equation (4.5) follows from (3.3), or from the one-Lipschitz property

\[
|(h+\lambda)_+-(Q(h)+\lambda)_+|\le|h-Q(h)|.         \tag{4.7}
\]

More generally, every complete configuration on \(N\) sites changes its
local energy by at most

\[
\sum_i|h_i-Q(h_i)|x_i\le\frac{\eta N}{2}.            \tag{4.8}
\]

Any identical aggregate interaction cancels in this comparison, and taking
maxima preserves the bound. A grid value is never rounded again, so each
site contributes to (4.8) once. This proves the depth-independent assertion.
\(\square\)

This yields a useful joint limit. If

\[
\eta_N\to0,\qquad \frac{B\log N}{N}\ll\eta_N,         \tag{4.9}
\]

then both response error per site and state bits per site tend to zero. For
example, \(\eta_N=B/\sqrt N\) gives

\[
\text{absolute error }O(B\sqrt N),\qquad
\log_2|\mathcal S|=O(\sqrt N\log N).                 \tag{4.10}
\]

## 5. Lower bounds and the rate caveat

### 5.1 Exact and microscopic accuracy

Fix the arithmetic grid in Section 4 and a common mass \(n\). Distinct
histograms have distinct sorted slope lists, hence distinct profiles and
responses by Theorem 3.1. Every exact summary on this grid subclass
therefore needs exactly at least

\[
\binom{n+M-1}{M-1}                                  \tag{5.1}
\]

states, matching the histogram construction.

There is also a separated version. For two distinct grid histograms, choose
the first rank at which their sorted slopes differ. The preceding prefix
sums agree and the next difference is a nonzero integer multiple of
\(\Delta\). Hence (3.3) gives response distance at least \(\Delta\).
Any decoder with uniform error \(<\Delta/2\) still needs all the states in
(5.1).

This is not a lower bound at the much larger allowance
\(\Theta(\eta n)\).

### 5.2 Macroscopic packing

Here is an elementary packing at normalized tolerance \(\varepsilon\). Choose
integers \(q,s\) satisfying

\[
q\le\frac n8,\qquad q^2\le\frac{B}{64\varepsilon},\qquad
s=\left\lfloor\frac{n}{4q}\right\rfloor.             \tag{5.2}
\]

Let \(w=B/(2q)\), and take \(q\) intervals
\([c_j-w,c_j+w]\) with disjoint interiors tiling \([-B/2,B/2]\).
For each \(z\in\{0,1\}^q\), construct a mass-\(n\) field multiset \(A_z\):

* if \(z_j=0\), place \(2s\) fields at \(c_j\);
* if \(z_j=1\), place \(s\) fields at each of \(c_j-w,c_j+w\);
* place all unused fields at zero.

Every alternative has the same mass and total field sum. Flipping bit \(j\)
changes the response by

\[
D_j(\lambda)=s\big[
(c_j-w+\lambda)_++(c_j+w+\lambda)_+
-2(c_j+\lambda)_+\big].                            \tag{5.3}
\]

This is a triangular function supported on
\([-c_j-w,-c_j+w]\), with height \(sw\) at \(\lambda=-c_j\).
The supports have disjoint interiors. If \(z,z'\) differ at \(j\), evaluating
their response difference at \(-c_j\) eliminates every other changed bit and
gives

\[
d(A_z,A_{z'})\ge sw.                                \tag{5.4}
\]

Since \(q\le n/8\), one has \(s\ge n/(8q)\), and therefore

\[
sw\ge\frac{Bn}{16q^2}\ge4\varepsilon n.             \tag{5.5}
\]

The \(2^q\) systems are separated by more than \(2\varepsilon n\).
Two such systems cannot share a state decoded with error at most
\(\varepsilon n\), by the triangle inequality. Thus

\[
\log_2K_{\varepsilon n}\ge q
=\Omega\!\left(\min\{n,\sqrt{B/\varepsilon}\}\right) \tag{5.6}
\]

whenever the right side is above a fixed constant.

Taking \(\varepsilon\) to be a constant multiple of \(\eta\) proves (1.4).
This does not match the histogram upper bound. The histogram retains every
bin count exactly, while \(R_A\) is an integrated empirical-distribution
query and can be less sensitive to fine count changes. No
\(\Omega((B/\eta)\log n)\) macroscopic lower bound is proved here.

## 6. Known aggregate energy and the response roof

Let a known aggregate term \(\Phi_n(k)\) be added and put

\[
q_A(k)=p_A(k)+\Phi_n(k),\qquad
L_A(\lambda)=\max_k\{q_A(k)+\lambda k\}.             \tag{6.1}
\]

Define the roof on \([0,n]\) by

\[
\overline q_A(u)=
\max\left\{\sum_k\alpha_kq_A(k):
\alpha\in\Delta_{n+1},\ \sum_k\alpha_k k=u\right\}.  \tag{6.2}
\]

It is the least concave majorant of the points \((k,q_A(k))\), and

\[
L_A(\lambda)=\max_{u\in[0,n]}
\{\overline q_A(u)+\lambda u\},\qquad
\overline q_A(u)=\inf_\lambda\{L_A(\lambda)-\lambda u\}. \tag{6.3}
\]

These identities follow by maximizing a linear functional over convex
mixtures of the lifted points; a supporting line at each point of the
concave roof gives the inverse identity.

Therefore linear fields recover all raw coordinates exactly when \(q_A\) is
discrete concave. If \(q_A\) is nonconcave, they recover only its roof.
The underlying local-field histogram remains sufficient because it
reconstructs \(p_A\), after which the known \(\Phi_n\) can be added. Whether
the smaller roof is closed under future composition depends on the cross
interaction.

## 7. Fixed quadratic mean field

Consider the genuinely closed fixed-coefficient model

\[
H_A^J(x)=\sum_{i\in A}h_i x_i+J\binom{K_A}{2}.       \tag{7.1}
\]

Its conditional profile and merge law are

\[
q_A(k)=p_A(k)+J\binom{k}{2},                         \tag{7.2}
\]

\[
q_{A\sqcup C}(t)
=\max_{k+\ell=t}\{q_A(k)+q_C(\ell)+Jk\ell\}.         \tag{7.3}
\]

The cross term obeys

\[
Jk\ell+J(k+\ell)r
=J\ell r+Jk(\ell+r)
=J(k\ell+kr+\ell r),                                \tag{7.4}
\]

so the raw profile law is bracket-independent.

### Theorem 7.1 (coarsest roof quotient and associative law)

For concave roofs \(f\) on \([0,n]\), \(g\) on \([0,m]\), define

\[
G_{f,g}(t)=
\max_{\substack{u\in[0,n],\,v\in[0,m]\\u+v=t}}
\{f(u)+g(v)+Juv\},                                  \tag{7.5}
\]

\[
f\star_Jg=\operatorname{cav}G_{f,g}.                \tag{7.6}
\]

Then

\[
\overline q_{A\sqcup C}
=\overline q_A\star_J\overline q_C.                 \tag{7.7}
\]

The operation is associative on realizable roofs. For contexts generated by
repeated same-\(J\) block merges followed by a linear terminal field,
\((n,\overline q_A)\) is the coarsest exact state.

#### Proof

For any means \(u,v\), choose child occupancy distributions attaining the
two roofs and take their product. Independence gives expected cross energy
\(Juv\). The expected energy cannot exceed the largest pure-pair energy.
Conversely, pure occupancies are among the candidates. Thus, for every
\(\lambda\),

\[
\begin{aligned}
&\max_{u,v}\{\overline q_A(u)+\overline q_C(v)
             +Juv+\lambda(u+v)\}\\
&\quad=\max_{k,\ell}\{q_A(k)+q_C(\ell)+Jk\ell
                      +\lambda(k+\ell)\}.             \tag{7.8}
\end{aligned}
\]

The right side is the parent response. Grouping the left side by \(t=u+v\)
gives the linear response of \(G\); applying (6.3) gives its concave
majorant, proving (7.7).

For three roofs, either bracketing has terminal response

\[
\max_{u,v,w}\{
f(u)+g(v)+h(w)+J(uv+uw+vw)+\lambda(u+v+w)\}.         \tag{7.9}
\]

To see this directly, apply (7.8) to the inner pair with induced field
\(\lambda+Jw\), or symmetrically with field \(\lambda+Ju\).
The two responses agree for every \(\lambda\); (6.3) then makes the roofs
equal. This proves associativity.

The empty future forces every exact state to determine \(L_A\), hence the
roof. Conversely, (7.7) answers every generated future context. Therefore
the roof is coarsest. \(\square\)

### 7.2 Concavity-preserving and concavity-destroying regimes

From the sorted fields,

\[
q_A(k)-q_A(k-1)=a_k+J(k-1).                         \tag{7.10}
\]

Hence \(q_A\) is discrete concave exactly when

\[
a_k-a_{k+1}\ge J\qquad(1\le k<n).                   \tag{7.11}
\]

In particular, \(J\le0\) preserves concavity for every heterogeneous block
and every union. Then linear fields recover \(q_A\); subtracting the known
quadratic recovers \(p_A\), and differences recover the histogram. Thus the
histogram is again the coarsest exact state. Positive \(J\) can preserve
concavity for a special well-separated list, but not uniformly over
arbitrary bounded fields and merges.

For a strict-quotient example, let \(n=2\), \(J>0\), and
\(0<a<\min\{B,J/2\}\). The field multisets

\[
A=\{0,0\},\qquad A'=\{a,-a\}                         \tag{7.12}
\]

have raw profiles

\[
q_A=(0,0,J),\qquad q_{A'}=(0,a,J).                   \tag{7.13}
\]

Both roofs are the line from \((0,0)\) to \((2,J)\), although their
\(k=1\) optima differ. By Theorem 7.1, no sequence of same-\(J\) block
futures and linear terminal fields distinguishes them.

There is a complete collapse at strong positive coupling. If \(J\ge2B\),

\[
[q_A(k+1)-q_A(k)]-[q_A(k)-q_A(k-1)]
=a_{k+1}-a_k+J\ge0.                                 \tag{7.14}
\]

Thus \(q_A\) is convex, and its least concave majorant is the endpoint chord

\[
\overline q_A(u)=\frac{u}{n}
\left(\sum_i h_i+J\binom n2\right).                 \tag{7.15}
\]

At fixed mass the minimal state collapses to the single number
\(\sum_i h_i\). This is a sharp warning against asserting a histogram lower
bound after concavity has been destroyed.

The quantized-histogram upper bound remains valid for every sign and size of
\(J\): the configurationwise comparison (4.8) does not involve \(J\), so the
quadratic term cannot amplify field-rounding error.

### 7.3 Other quadratic normalizations

A fixed aggregate quadratic
\(\Phi(k)=\alpha k^2+\beta k\), with \(\Phi(0)=0\), has cross increment

\[
\Phi(k+\ell)-\Phi(k)-\Phi(\ell)=2\alpha k\ell,       \tag{7.16}
\]

so Theorem 7.1 applies with \(J=2\alpha\). The same holds in magnetization
coordinates because mass and magnetization add.

For Curie--Weiss normalization \(J K^2/N\), the meaning of \(N\) during
composition is essential.

* If a final target mass \(N\) is fixed and every partial block uses the same
  coefficient \(J/N\), the cross term is bilinear and the roof closes.
* If a mass-\(n\) child was summarized using \(J/n\) and is later
  reinterpreted inside mass \(n+m\) with coefficient \(J/(n+m)\), its
  internal curvature has changed. A roof made at the old curvature may have
  discarded newly exposed points. Retaining the underlying \(p\)/histogram
  and mass is safe; blindly merging native normalized roofs is not.

More generally, with known \(\Phi_n\), raw profiles always obey

\[
q_{n+m}(t)=\max_{k+\ell=t}
\{q_n(k)+q_m(\ell)+C_{n,m}(k,\ell)\},               \tag{7.17}
\]

\[
C_{n,m}(k,\ell)
=\Phi_{n+m}(k+\ell)-\Phi_n(k)-\Phi_m(\ell).          \tag{7.18}
\]

This raw law is associative by telescoping. Passing to roofs before merging
is justified by the product-mixture proof when the cross term is separately
affine in the two child aggregate variables. A generic size-changing
\(C_{n,m}\) need not be, so roof closure must be proved rather than assumed.

## 8. Further context and closure caveats

### 8.1 Which fields expose which coordinates?

The exact statements are:

* one-parameter linear fields expose every coordinate of the concave top-\(k\)
  profile \(p\), via (3.2);
* the same fields expose only \(\operatorname{cav}q\) when an aggregate
  energy makes \(q\) nonconcave;
* if the allowed terminal family is enlarged to tunable negative quadratic
  penalties

  \[
  -L(k-k_0)^2,\qquad L>0,                            \tag{8.1}
  \]

  every coordinate of any finite bounded \(q\) can again be exposed: choose
  \(L\) larger than all competing energy differences.

The last family is still simple and aggregate, not an arbitrary lookup
table, but it is a strictly richer experiment. Thus “aggregate field” is
not a complete specification; functional form and coefficient range matter.

### 8.2 Non-biaffine interactions can resurrect hidden data

Theorem 7.1 uses
\(\mathbb E[K_AK_C]=\mathbb E K_A\,\mathbb E K_C\) under independent child
mixtures. A nonlinear cross term such as \(K_A^2K_C\), or a coefficient that
changes after a merge, requires moments not retained by a scalar roof.
Previously hidden raw coordinates can then matter. The safe choices are to
retain the raw conditional profile, enlarge the feature vector until the
cross term is separately affine, or prove a model-specific closure identity.

### 8.3 Value error is not optimizer error

Theorem 4.1 controls optimum values. Near a tie, an arbitrarily small field
perturbation can change the maximizing occupancy or selected sites. No
optimizer-stability statement follows without a margin assumption. Likewise,
an unlabelled histogram does not preserve site identities or old block
labels; those are intentionally outside the declared context family.

## 9. Benchmark conclusion

The partial mean-field benchmark can be closed rigorously as follows.

1. For arbitrary bounded heterogeneous binary local fields, the coarsest
   exact state under uniform aggregate fields and unlabelled block union is
   the discrete-concave top-\(k\) profile, equivalently its slope multiset.
2. Its exact composition law is max-plus convolution, or multiset union in
   slope coordinates.
3. A fixed-grid histogram is an exact merge homomorphism with
   \(O((B/\eta)\log N)\) bits in the intended mesh regime and
   \(O(\eta N)\) response error. Error grows with total mass, not merge depth;
   choosing \(\eta_N=B/\sqrt N\) gives simultaneous sublinear error and
   sublinear state bits.
4. Exact grid recovery has a matching binomial lower bound. At macroscopic
   distortion the proved general lower bound is only
   \(\Omega(\min\{N,\sqrt{B/\eta}\})\) bits, so minimax optimality of the
   histogram rate remains open.
5. Fixed quadratic mean-field energy is legitimately closed through its
   bilinear cocycle. Concavity-preserving curvature leaves the histogram
   minimal; concavity-destroying curvature replaces it by a potentially
   strict concave-roof quotient. Strong positive curvature can collapse the
   state to mass and total field.
6. Persistent block labels, tunable additional curvature, non-biaffine cross
   terms, or changing normalization across scales alter the quotient and
   must be declared separately.

This gives a non-tautological exact theorem, a stable macroscopic
approximation theorem, and explicit lower-bound limits for the mean-field
campaign.
