# A compositional carrier--exposure law for benchmark state growth

**Status.** Independent proof-level report for the benchmark campaign.  No
repository file is edited.  All logarithms in state-count statements are
base two.  Covering numbers below are *external* unless explicitly stated;
an internal cover is obtained at twice the radius.

## 1. The point of the theorem

The four benchmarks currently use states which look superficially
different:

* a table on projective boundary assignments for separator Max-Cut/CSP;
* a two-entry (or finite-width) transfer state for one-dimensional Ising;
* a histogram of quantized site types for heterogeneous mean field;
* a vector of tropical aggregates for a lumpable weighted automaton.

There is one operational distinction which predicts their growth.  A state
must satisfy **both** of the following tests.

1. Its coordinates must carry every independently realizable response
   direction exposed by the allowed future queries.
2. Its coordinates must form a congruence for the future-composition
   algebra; a small one-time cover is not enough.

The theorem below makes the first test quantitative and the second exact.
It also distinguishes a semantic carrier from a smaller syntactic
max-affine presentation.  Its additive-atom clause gives a sharp polynomial
law and, as a new benchmark consequence, determines the response-metric
entropy of a finite-grid heterogeneous mean-field model.

## 2. Response systems and two kinds of presentation

Let \(\mathcal A\) be a class of partial objects, \(\mathcal Q\) a declared
class of future queries, and

\[
 F_A:\mathcal Q\longrightarrow \mathbb R
\]

the future-response function.  Write

\[
 d_{\rm rsp}(A,A')=\|F_A-F_{A'}\|_{\infty}.
 \tag{2.1}
\]

Everything below also applies after quotienting response functions by
constants, replacing the sup norm by \(\operatorname{osc}/2\).  A fixed
anchor query converts that version to the literal version.

Assume a monoid \(\mathsf C\) of continuations acts on an ambient class of
partial objects.  An **exact stable \(r\)-carrier** consists of

\[
 s:\mathcal A\to\mathbb R^r,
 \qquad D_q:\mathbb R^r\to\mathbb R\quad(q\in\mathcal Q),
 \tag{2.2}
\]

and updates \(U_C\) such that

\[
 F_A(q)=D_q(s(A)),\qquad
 |D_q(z)-D_q(z')|\le\|z-z'\|_\infty,
 \qquad s(CA)=U_C(s(A)).
 \tag{2.3}
\]

The last identity is the congruence requirement.  We say a test subfamily
has **range \(R\)** when its carrier image is contained in a translate of
\([-R,R]^r\); continuations may move it to a different resource slice or a
larger range.  No Lipschitz hypothesis
on \(U_C\) is needed for exact reuse, although one is useful for approximate
reuse.  A \(\zeta\)-accurate carrier may replace the first equality in
(2.3) by uniform error at most \(\zeta\); this is only a static approximation
unless its updates come with a separate depth-uniform error law.

A different notion is a **homogeneous binary max-affine grammar with \(g\)
shared parameters**:

\[
 F_A(q)=\max_{a\in E_q}\langle a,\theta(A)\rangle+b(q),
 \qquad E_q\subseteq\{0,1\}^g,
 \tag{2.4}
\]

where the response shapes have radius at most \(R\).  The same grammar is
called stable if all continuation derivatives remain of the form (2.4)
with the same \(g\) parameters and admit an exact parameter update.  A
grammar parameter is not a semantic response coordinate: optimizer changes
can make \(g\) parameters generate substantially more than \(g\) robust
response directions.

Finally, an **\((a,\alpha,k)\) exposed response cube** is a map

\[
 \iota:[-a,a]^k\to\mathcal A
 \quad\hbox{such that}\quad
 d_{\rm rsp}(\iota(u),\iota(v))
 \ge \alpha\|u-v\|_\infty .
 \tag{2.5}
\]

In applications (2.5) is not assumed abstractly: it is certified by actual
future queries which read the selected coordinates, or by balanced
two-sided exposure in the projective setting.

## 3. Main theorem: carrier, grammar, and exposure

### Theorem 3.1 (compositional carrier--exposure sandwich)

Let a response test family have range \(R\) inside an exact stable
\(r\)-carrier.
Then for every \(\epsilon>0\),

\[
 \log \operatorname{Cov}^{\rm ext}_{\epsilon}
       \{F_A:A\in\mathcal A\}
 \le
 r\log\!\left(1+{2R\over\epsilon}\right).
 \tag{3.1}
\]

For a \(\zeta\)-accurate carrier the same right side bounds
\(\log\operatorname{Cov}^{\rm ext}_{\zeta+\epsilon}\).

If the response system contains an \((a,\alpha,k)\) exposed cube and
\(0<\epsilon<\alpha a/2\), then

\[
 \log \operatorname{Cov}^{\rm ext}_{\epsilon}
       \{F_A:A\in\mathcal A\}
 \ge
 k\log\!\left(1+left\lfloor {\alpha a\over2\epsilon}
                         \right\rfloor\right).
 \tag{3.2}
\]

Consequently every exact stable carrier with the stated range satisfies

\[
 r\ge
 {k\log(1+\lfloor\alpha a/(2\epsilon)\rfloor)
  \over \log(1+2R/\epsilon)}.
 \tag{3.3}
\]

If instead the response class has a stable homogeneous binary max-affine
grammar with \(g\) shared parameters, then

\[
 \begin{split}
 \log \operatorname{Cov}^{\rm ext}_{\epsilon}
 &\le
 g\log\!\left(4\,{3^g+1\over2}\right)
 +g\log\!\left(1+{2R\over\epsilon}\right).
 \end{split}
 \tag{3.4}
\]

Thus an exposed packing with \(L\) bits forces

\[
 L\le O\!\left(g^2+g\log(1+R/\epsilon)\right).
 \tag{3.5}
\]

In particular, at fixed relative precision, \(L=2^{\Omega(w)}\)
forces carrier rank \(r=2^{\Omega(w)}\), whereas it forces only the weaker
but still exponential grammar bound \(g=2^{\Omega(w)}\) (with a possibly
halved constant in the exponent).

The carrier or grammar remains valid at arbitrary continuation depth because
of the exact update in (2.3), not because of the static covering argument.

#### Proof

Quantize every carrier coordinate on a mesh of width at most
\(2\epsilon\).  There are at most \(1+2R/\epsilon\) choices per coordinate
(ceilings change this by an inessential endpoint convention).  A nearest
grid point is within \(\epsilon\) in sup norm, and every readout is
one-Lipschitz, proving (3.1); a \(\zeta\)-accurate decoder adds \(\zeta\).
More formally one may use a maximal
\(\epsilon\)-separated subset of the interval, yielding exactly the same
display with a ceiling.

For (3.2), put a grid in each coordinate of \([-a,a]^k\) with spacing
\(4\epsilon/\alpha\).  It has at least

\[
 \left(1+\left\lfloor{\alpha a\over2\epsilon}\right\rfloor\right)^k
\]

points whose response distance is at least \(4\epsilon\).  One
external radius-\(\epsilon\) ball contains at most one of them.  This proves
(3.2), and (3.3) follows by comparison with (3.1).

For (3.4), every optimizer change in (2.4) lies on a central hyperplane with
normal in \(\{-1,0,1\}^g\setminus\{0\}\), modulo reversal.  Hence at most

\[
 N_g={3^g-1\over2}
\]

distinct hyperplanes occur, independently of the number of queries and
witnesses.  All relatively open faces of their arrangement number at most
\([4(N_g+1)]^g\).  On each face a fixed tie rule makes the full response one
linear map of \(\theta\); after quotienting constants its image has
dimension at most \(g\).  A radius-\(R\) subset of a \(g\)-dimensional
normed space has an external \(\epsilon\)-cover of size at most
\((1+2R/\epsilon)^g\).  Multiplying by the face count proves (3.4).
Exact closure under continuation is simply induction in (2.3), and is
logically separate from all covering estimates. \(\square\)

### Why Theorem 3.1 is not merely Myhill--Nerode

Myhill--Nerode identifies the exact contextual equivalence relation, but it
does not supply (3.2)--(3.5): robust scale-dependent lower bounds, the
precision cost per exposed direction, or a lower bound on a real-parameter
max-affine grammar despite unlimited parameter precision.  Nor does an
ordinary dimension count distinguish the semantic carrier rank \(r\) from
the syntactic grammar rank \(g\).  The future algebra enters twice: actual
queries certify (2.5), while the exact derivative update is what makes the
same coordinates reusable.

## 4. A sharp polynomial branch: finite additive atoms

The preceding theorem says when an exponential exposed family rules out a
small carrier.  The next clause identifies a broad algebra where the carrier
provably stays polynomial.

Let \(\phi_1,\ldots,\phi_d\in\ell_\infty(\mathcal Q)\), let

\[
 \mathcal H_{n,d}=\{c\in\mathbb N^d:\sum_{j=1}^dc_j=n\},
 \qquad F_c=\sum_{j=1}^dc_j\phi_j,
 \tag{4.1}
\]

and compose histograms by addition.  On

\[
 V=\{z\in\mathbb R^d:\sum_jz_j=0\}
\]

define \(Tz=\sum_jz_j\phi_j\),

\[
 \alpha=\inf_{z\in V,\ \|z\|_\infty=1}\|Tz\|_\infty,
 \qquad L=\max_j\|\phi_j\|_\infty,
 \tag{4.2}
\]

and the lattice margin

\[
\sigma=\inf_{0\ne z\in V\cap\mathbb Z^d}\|Tz\|_\infty.
\tag{4.3}
\]

There is also an exact arithmetic rank which must not be confused with
ordinary affine dimension.  Put

\[
 \Gamma_\Phi=
 \left\langle\phi_1-\phi_d,\ldots,
                    \phi_{d-1}-\phi_d\right\rangle_{\mathbb Z}
 \subset \ell_\infty(\mathcal Q),
 \qquad r_{\mathbb Z}=\operatorname{rank}_{\mathbb Z}\Gamma_\Phi.
 \tag{4.3a}
\]

This is a finitely generated torsion-free abelian group, so its rank is
well-defined and lies between zero and \(d-1\).
For projective responses, the same definition is made in
\(\ell_\infty(\mathcal Q)/\mathbb R\mathbf1\); otherwise a query-independent
baseline would be counted as an exposed generator.

### Theorem 4.1 (additive-generator response law)

For arbitrary \(\phi_1,\ldots,\phi_d\), contextual equivalence is exactly

\[
 c\sim c'\quad\Longleftrightarrow\quad T(c-c')=0,             \tag{4.3b}
\]

this equivalence is a congruence for histogram addition, and the number
\(N_n\) of exact mass-\(n\) contextual states satisfies

\[
 \boxed{N_n=\Theta_\Phi(n^{r_{\mathbb Z}}).}                 \tag{4.3c}
\]

Suppose in addition that \(d\ge2\) and \(\alpha>0\).  Then:

1. the histogram \(c\) is the coarsest exact contextual state and
   \(r_{\mathbb Z}=d-1\);
2. composition is the exact depth-independent update \(c\leftarrow c+c'\);
3. for every \(\epsilon>0\), with

   \[
   s_\epsilon=1+\left\lfloor{2\epsilon\over\alpha}\right\rfloor,
   \tag{4.4}
   \]

   one has

   \[
   \boxed{
   \left(1+\left\lfloor{n\over(d-1)s_\epsilon}\right\rfloor\right)^{d-1}
   \le \operatorname{Cov}^{\rm ext}_{\epsilon}\{F_c:c\in\mathcal H_{n,d}\}
   \le
   \left(2+\left\lceil{L(d-1)n\over\epsilon}\right\rceil\right)^{d-1}.}
   \tag{4.5}
   \]

4. If \(\sigma>0\), then for every \(0<\epsilon<\sigma/2\),

   \[
   \boxed{
   \operatorname{Cov}^{\rm ext}_{\epsilon}\{F_c:c\in\mathcal H_{n,d}\}
   =\binom{n+d-1}{d-1}.}
   \tag{4.6}
   \]

Thus every finite additive atom algebra has polynomial exact contextual
growth, with exponent equal to the arithmetic rank of the atom-response
differences.  For a nondegenerate query family the exponent is \(d-1\), and
throughout the
mesoscopic range its metric entropy has the expected
\((d-1)\log(n/\epsilon)\) law up to constants depending on the conditioning
of the query features.  This polynomiality is a theorem about the additive
future algebra, not a claim that the full microscopic landscape is small.

#### Proof

At fixed mass,

\[
 F_c=n\phi_d+\sum_{j<d}c_j(\phi_j-\phi_d).        \tag{4.3d}
\]

This proves (4.3b), and adding the same histogram preserves equality, so the
relation is a congruence.  Choose an isomorphism from \(\Gamma_\Phi\) to
\(\mathbb Z^{r_{\mathbb Z}}\).  The finitely many atom differences have
bounded integer coordinates in this basis.  Every sum in (4.3d) therefore
lies in a box of side \(O_\Phi(n)\), giving \(N_n=O_\Phi(n^{r_{\mathbb Z}})\).
Conversely, choose \(r_{\mathbb Z}\) of the displayed atom differences which
are linearly independent over \(\mathbb Q\).  Letting each of their counts
range independently from zero to
\(\lfloor n/r_{\mathbb Z}\rfloor\), with all remaining mass put in type
\(d\), gives \(\Omega_\Phi(n^{r_{\mathbb Z}})\) distinct sums.  (For rank
zero there is one.)  This proves (4.3c).

Since \(\alpha>0\), the restriction of \(T\) to \(V\) is injective.  At
fixed mass, equality \(F_c=F_{c'}\) therefore forces \(c=c'\).  Conversely
the histogram plainly answers every declared query.  Appending a block adds
its histogram and its response, proving exact closure and coarseness.  The
\(d-1\) differences are then linearly independent over \(\mathbb R\), hence
also over \(\mathbb Q\), so \(r_{\mathbb Z}=d-1\).

For the lower bound in (4.5), let

\[
 K=\left\lfloor {n\over(d-1)s_\epsilon}\right\rfloor
\]

and take all histograms

\[
 c(u)=\bigl(s_\epsilon u_1,\ldots,s_\epsilon u_{d-1},
 n-s_\epsilon\sum_{j<d}u_j\bigr),
 \qquad u_j\in\{0,\ldots,K\}.
 \tag{4.7}
\]

They are nonnegative.  Two distinct ones differ by at least
\(s_\epsilon\) in sup norm, so their response distance is at least
\(\alpha s_\epsilon>2\epsilon\).  They form the claimed packing.

For the upper bound, put \(h=\epsilon/[L(d-1)]\) (the case \(L=0\) is
trivial).  Round the first \(d-1\) coordinates of \(c\) to the nearest
multiple of \(h\), and set the last coordinate so the total remains \(n\).
The external centre need not be a nonnegative integer histogram.  The
\(\ell_1\) rounding error is at most \((d-1)h\), and hence its response error
is at most \(L(d-1)h=\epsilon\).  Each stored coordinate has at most the
number of values displayed in the right side of (4.5).

Finally, distinct mass-\(n\) integer histograms differ by a nonzero element
of \(V\cap\mathbb Z^d\), so their responses are \(\sigma\)-separated.  At
radius below \(\sigma/2\) every external ball contains at most one response.
There are exactly \(\binom{n+d-1}{d-1}\) histograms, proving (4.6).
\(\square\)

The condition \(\alpha>0\) may be weakened for exact recovery to injectivity
on the integer difference lattice.  It is deliberately retained in the
metric clauses because it is exactly the robust condition needed for
approximate response complexity.  Formula (4.3c) also exposes an important
precision pathology: a one-real-parameter query can assign rationally
independent atom responses and have \(r_{\mathbb Z}=d-1\).  Ordinary query
dimension alone therefore does not predict exact contextual growth.

## 5. New benchmark consequence: heterogeneous mean field

Take equally spaced fields

\[
 \gamma_j=-B+(j-1)\Delta,qquad
 \Delta={2B\over d-1},\qquad 1\le j\le d,
 \tag{5.1}
\]

and uniform chemical-potential queries \(\lambda\in[-B,B]\).  One site of
type \(j\) contributes

\[
 \phi_j(\lambda)=(\gamma_j+\lambda)_+.
 \tag{5.2}
\]

For \(z\in V\), let

\[
 G_z(\lambda)=\sum_jz_j(\gamma_j+\lambda)_+.
\]

At the knot \(\lambda=-\gamma_k\),

\[
 {G_z(-\gamma_k)\over\Delta}
 =S_k:=\sum_{j>k}z_j(j-k).                       \tag{5.3}
\]

The tail sums obey

\[
 \sum_{j\ge k+1}z_j=S_k-S_{k+1}.                \tag{5.4}
\]

Consequently

\[
 \|z\|_\infty
 \le4\max_k|S_k|
 \le {4\over\Delta}\|G_z\|_\infty,             \tag{5.5}
\]

so \(\alpha\ge\Delta/4\).  If \(z\) is a nonzero integer vector, some
\(S_k\) is a nonzero integer, hence \(\|G_z\|_\infty\ge\Delta\).  Moving
one site between adjacent bins attains equality, so

\[
 \boxed{\sigma=\Delta.}                          \tag{5.6}
\]

Since \(L\le2B\), Theorem 4.1 gives

\[
 \left(1+left\lfloor
 {n\over(d-1)(1+\lfloor8\epsilon/\Delta\rfloor)}
 \right\rfloor\right)^{d-1}
 \le \operatorname{Cov}^{\rm ext}_\epsilon
 \le
 \left(2+\left\lceil{2B(d-1)n\over\epsilon}\right\rceil\right)^{d-1},
 \tag{5.7}
\]

and, more sharply, for \(0<\epsilon<\Delta/2\),

\[
 \boxed{
 \operatorname{Cov}^{\rm ext}_\epsilon
 =\binom{n+d-1}{d-1}.}                           \tag{5.8}
\]

Thus the contextual framework independently recovers a histogram state,
proves it minimal for the declared uniform-field futures, and obtains a new
exact/approximate response-rate law rather than merely renaming the usual
mean-field aggregate.

For arbitrary fields in \([-B,B]\), round each field once to a common mesh
of spacing at most \(\eta\).  The scalar hinge is one-Lipschitz in its field,
so the response error is at most \(\eta n/2\) under nearest rounding.  Counts
add exactly under disjoint union.  On a merge tree, each site is charged
once; the root error is at most \(\eta N/2\), independent of tree depth.
The number of exact histogram states at mass \(n\) is

\[
 \binom{n+d-1}{d-1},\qquad d=1+\left\lceil{2B\over\eta}\right\rceil.
 \tag{5.9}
\]

The same carrier survives any known aggregate interaction evaluated from
the total occupancy, including a fixed uniform quadratic term, because the
configuration-wise local-field error is unchanged.

## 6. Benchmark score under the theorem

| benchmark | carrier predicted from futures | exact update | exposed/generator rank | growth verdict |
|---|---|---|---:|---|
| width-\(w\) Max-Cut/CSP | projective boundary profile | table gluing / max-plus elimination | \(2^{w-1}\) in pure Max-Cut | exponential in \(w\) |
| 1D Ising | projective two-spin profile | \(2\times2\) max-plus transfer | one | fixed-dimensional |
| \(d\)-type heterogeneous mean field | atom histogram | count addition | \(r_{\mathbb Z}=d-1\) on the uniform grid | \(\Theta(n^{d-1})\) exact states |
| \(r\)-class lumpable automaton | tropical aggregate vector | quotient max-plus transition | reachable pins, at most \(r\) | fixed-dimensional in word length |

### 6.1 Bounded-separator Max-Cut/CSP: exponential branch

For a boundary of width \(w\), let

\[
 q=|\{\pm1\}^w/\{\pm1\}|=2^{w-1}.
\]

Arbitrary future attachments read the projective boundary profile.  The
private Max-Cut compiler realizes every profile table (up to its calibrated
constant), so at scalar precision the response class contains a
\(q\)-coordinate cube.  The profile table is a stable \(q\)-carrier under
gluing.  Theorem 3.1 therefore gives matching

\[
 \log\operatorname{Cov}_\epsilon
 =\Theta\bigl(q\log(W/\epsilon)\bigr),            \tag{6.1}
\]

which makes exponential state growth a consequence of realizable query
exposure, not an assumption imported from treewidth dynamic programming.

Under unit boundary load, the response class is the full Lipschitz ball on
the projective Hamming cube.  At error \(\epsilon w\), the existing balanced
exposure theorem gives

\[
 \log\operatorname{Cov}_{\epsilon w}
 \ge2^{(1-H_2(2\epsilon)+o(1))w}.                \tag{6.2}
\]

Since the natural carrier range is \(O(w)\), (3.3) forces exponentially
many stable bounded-range carrier coordinates.  Combining (6.2) with (3.4)
also forces exponentially many binary shared parameters, though with the
expected square-root loss in the response-bit exponent.  This separates
semantic table size from the number of edge-weight parameters in a compiler.

### 6.2 One-dimensional / finite-width Ising: fixed branch

For a chain with one exposed spin, the conditional profile has two entries;
modulo constants it has one coordinate.  Arbitrary endpoint fields expose
that coordinate, and a \(2\times2\) max-plus transfer kernel updates it.
Thus \(r=k=1\) and

\[
 \log\operatorname{Cov}_\epsilon
 =\Theta(\log(1+R/\epsilon)).                    \tag{6.3}
\]

For fixed strip width \(w_0\), the same argument has at most
\(2^{w_0}-1\) projective coordinates, a constant independent of chain
length.  This is exactly why repeated composition does not create growing
microscopic information at fixed width.

### 6.3 Heterogeneous mean field: polynomial branch

Sections 4--5 give the precise result.  A fixed quantized atom alphabet has
\(d-1\) independent count generators and exactly
\(\binom{n+d-1}{d-1}\) states at sub-grid response error.  Composition is
histogram addition.  The polynomial exponent is forced by uniform-field
queries through the hinge finite-difference identity (5.3), rather than by
assuming magnetization or a histogram in advance.

### 6.4 Lumpable weighted automata: fixed exposed quotient

Under exact tropical lumpability into \(r\) aggregate classes, the aggregate
vector is a stable \(r\)-carrier.  If reachable suffixes robustly pin all
aggregate coordinates, they expose an \(r\)-cube.  Theorem 3.1 gives the
matching law

\[
 \log\operatorname{Cov}_\epsilon
 =\Theta(r\log(1+B/\epsilon)).                   \tag{6.4}
\]

on the reachable aggregate box.  If only \(k<r\) aggregates are reachable,
the lower bound is \(k\), which is the correct limitation: formal states
not exposed by a legal suffix need not be stored for this experiment.

## 7. What the theorem predicts, and what it does not

The common law is now concise:

> **State growth is the metric entropy of the realizable exposed carrier,
> subject to that carrier being a congruence for the future algebra.**

This predicts two different mechanisms for small states.

* Fixed boundary/aggregate rank gives a bounded-dimensional carrier
  (Ising, lumped automata).
* A finite additive atom algebra gives a histogram simplex, hence polynomial
  growth in mass, with exact exponent equal to its arithmetic response
  generator rank (mean field).

It also predicts the separator obstruction: arbitrary lookup futures and a
universal compiler expose one independent response coordinate for each
boundary assignment, so the carrier rank itself is exponential in width.

The following limitations are essential and sharp.

1. **Static entropy is not reuse.**  The transition-toll example in the
   repository has a vanishing one-step defect and fixed long-depth drift.
   Omitting (2.3) would make Theorem 3.1 a one-time sketching theorem, not a
   compositional theorem.
2. **Query affine dimension is not exposure dimension.**  A one-dimensional
   curve of weighted-automaton suffixes can robustly expose every one of
   \(p\) coordinates.  The lower certificate must use realizable query
   action, not the dimension of the query parameter space.
3. **Semantic rank is not grammar rank.**  Optimizer fans let \(g\) binary
   parameters generate \(\exp(O(g^2))\) robust response cells.  Hence the
   carrier lower bound cannot be transferred to syntax without the weaker
   grammar estimate (3.4).  That estimate is genuinely homogeneous:
   arbitrary witness-specific offsets could create unboundedly many
   parallel comparison hyperplanes and require an additional offset/resource
   count.
4. **Range and regularity cannot be omitted.**  One unrestricted real number
   can encode an arbitrary finite state set discontinuously.  The bounded
   response-scale range and Lipschitz readouts prevent this fake
   compression; homogeneous binary grammar gives a different precision-free
   safeguard.
5. **Finite atoms and conditioning are substantive.**  If \(d\) grows with
   system size, the histogram law may itself become exponential.  If
   \(\alpha\) tends to zero, exact states can remain distinct while their
   macroscopic response entropy collapses.  If \(T\) has a kernel, the true
   state is the corresponding quotient of the histogram, not all counts.
   Exact growth is controlled by the arithmetic rank \(r_{\mathbb Z}\),
   while robust approximate growth is controlled by the real conditioning
   \(\alpha\); these can behave very differently under unlimited precision.
6. **Exact and approximate mean-field claims differ.**  The finite-grid law
   is exact.  Continuous fields require quantization and incur
   \(\eta n/2\) error.  A still smaller static metric net of histograms need
   not be closed under arbitrary-depth merging; the exact count carrier is.
7. **Declared futures matter.**  Anonymous mean-field union permits one
   histogram.  Futures which retain separate labels for old blocks require
   their tuple.  Similarly, inaccessible automaton suffixes do not
   contribute to exposure.

## 8. Director judgment and next theorem

This theorem is stronger than a taxonomy: it produces the previously
missing approximate response-rate law (5.7) for the mean-field histogram,
recovers the correct polynomial exponent at microscopic resolution, and
puts that result under the same operational lower certificate which proves
separator exponentiality and automaton minimality.  The framework is
therefore generative at least for finite additive-atom systems.

The remaining conceptual gap is dynamic rather than static.  The strongest
next theorem is:

> **Approximate congruence dichotomy.**  Characterize when an
> \(\epsilon\)-sufficient exposed carrier admits a derivative-compatible
> update whose total error is \(O(\epsilon\,\text{mass})\), rather than
> \(O(\epsilon\,\text{depth})\).  Prove that, on a finite recurrent future
> graph, every such mechanism factors (up to controlled error) through an
> additive atom homomorphism, a tropical gauge/coboundary, or a syndetic
> reset; or exhibit a fourth mechanism.

That statement would explain not only how many response coordinates are
needed, but exactly when an approximate response quotient can be reused.
