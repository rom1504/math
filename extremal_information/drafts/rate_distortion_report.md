# Extremal rate--distortion and query sufficiency for finite landscapes

## Frozen setup (fixed before consulting repository no-go artifacts)

Let \(X_n=\{-1,+1\}^n\), let \(\mathcal H_n\subset\mathbb R^{X_n}\) be a finite class of labelled energy landscapes, and for an external field \(h\in\mathbb R^n\) define the perturbed optimum-value query
\[
V_H(h):=\max_{x\in X_n}\{H(x)+\langle h,x\rangle\}.
\]
A possibly randomized \(R\)-bit sketch is a channel \(P_{Z\mid H}\) with \(|\operatorname{supp}Z|\le 2^R\), followed by a decoder \(z\mapsto \widehat V_z\) on a declared query family \(\mathcal U_n\).  Two distortions will be used and not changed after importing the repository evidence:
\[
d_\infty(H,z)=\sup_{h\in\mathcal U_n}|\widehat V_z(h)-V_H(h)|,
\qquad
d_2^2(H,z)=\mathbb E_{h\sim\nu_n}|\widehat V_z(h)-V_H(h)|^2.
\]
The induced query pseudometric is
\[
\rho_{\mathcal U}(H,H')=\sup_{h\in\mathcal U_n}|V_H(h)-V_{H'}(h)|.
\]
The rate--distortion function under a prior \(\Pi\) is
\[
R_\Pi(D)=\inf_{P_{Z\mid H},\widehat V:\ \mathbb E d(H,Z)\le D} I_\Pi(H;Z).
\]
“Query-sufficient” always means sufficient for this full perturbed query experiment (uniformly, or under the explicitly stated \(\nu_n\)); it never means merely retaining \(V_H(0)=\max H\).  Landscape-description information and optimizer-state information are counted separately: quantized dense pair couplings have \(\binom n2=\Theta(n^2)\) independent coordinates, while a labelled spin state has only \(n=\Theta(n)\) bits.

The report below will instantiate this setup, prove lower bounds, and then compare them with the repository's exact no-go ingredients.

## Result in one sentence

There is a \(2^{\binom n2}\)-element family of normalized Ising landscapes for which every unperturbed maximum is **the same number, zero**, yet an average-square accurate sketch of the optimum under a fixed family of external-field interventions needs
\[
\binom n2\,[1-h_2(D)]
\]
bits; under additive uniform error below one coupling quantum it needs exactly
\(\binom n2\) bits.  In contrast, the analogous state-location problem costs only
\(n[1-h_2(D)]\) bits.  The gap is landscape/counterfactual information versus
state/extremizer information, not a disguised demand to retain \(\max H\).

This is a finite, nonasymptotic information theorem.  It neither targets nor
implies convergence of \(M_n\).

## 1. Four notions that should not be conflated

For a finite class \(\mathcal H\), write
\[
\mathsf a(H)=(V_H(u))_{u\in\mathcal U}
\]
for its complete query-answer vector.

1. **A scalar optimum** is \(V_H(0)\).  It can be constant on a very rich
   class.
2. **A state certificate** is one \(x^*(H)\in\arg\max H\).  On
   \(\{\pm1\}^n\), its labelled description has at most \(n\) bits.
3. **A landscape description** identifies the interaction data.  A dense
   one-bit pair-interaction matrix has
   \(N=\binom n2=\Theta(n^2)\) bits (and \(bN\) bits at \(b\)-bit
   quantization).  This statement concerns the dense Ising parameterization;
   an unrestricted lookup table would of course be exponentially larger.
4. **A query-sufficient transcript** supports counterfactual optimum queries.
   Its necessary rate is governed by the geometry of
   \(\{\mathsf a(H):H\in\mathcal H\}\), and can be as large as the landscape
   description even when a scalar optimum and one optimizer are retained.

The operational hierarchy is therefore
\[
\text{scalar optimum}\quad\ll\quad
\text{one state}\quad\ll\quad
\text{uniform counterfactual response/landscape}
\]
on the hard family below.  The last two inequalities are information
inequalities, not computational-complexity assertions.

## 2. Literature coordinates and the exact model used here

The framework combines five classical coordinates, but the theorem below is
self-contained.

* [Shannon's 1959 fidelity paper](https://ieeexplore.ieee.org/document/5311476)
  defines rate--distortion by minimizing mutual information subject to an
  expected distortion constraint and computes the equiprobable binary/Hamming
  curve \(R(D)=1-h_2(D)\).  The proof below uses its one-shot converse, not an
  asymptotic coding theorem.
* [Blackwell's comparison of experiments](https://doi.org/10.1214/aoms/1177729032)
  orders experiments by all decision problems and characterizes exact
  simulation by a garbling kernel.  [Le Cam's approximate sufficiency
  paper](https://doi.org/10.1214/aoms/1177700372) supplies the approximate
  experiment-comparison viewpoint.  In the finite convention used here,
  \[
  \delta(\mathsf E,\mathsf F)
   :=\inf_K\sup_{\theta}
     \|K P^{\mathsf E}_\theta-P^{\mathsf F}_\theta\|_{\rm TV}.
  \]
  We do not rename value error as total-variation deficiency: we first decode
  the finite parameter and only then invoke this standard deficiency.
* [Kolmogorov and Tikhomirov's \(\varepsilon\)-entropy and
  \(\varepsilon\)-capacity](https://www.mathnet.ru/php/archive.phtml?jrnid=rm&option_lang=eng&paperid=7289&wshow=paper)
  are respectively logarithmic covering and packing numbers.  The relevant
  metric here is the query pseudometric \(\rho_{\mathcal U}\), rather than a
  norm on the raw coupling matrix.
* For quadratic-form and cut sketches, the distinction between one fixed query
  and all queries is essential.  Andoni--Krauthgamer--Woodruff prove an
  \(\Omega(n/\varepsilon^2)\)-bit lower bound for simultaneous cut
  approximation in [*The Sketching Complexity of Graph
  Cuts*](https://arxiv.org/abs/1403.7058) (in its stated
  \(\varepsilon>n^{-1/2}\) regime).  Andoni--Chen--Krauthgamer--Qin--
  Woodruff--Zhang give, among other regimes, \(\Omega(n^2)\) for general
  quadratic-form sketches and \(\Omega(n/\varepsilon^2)\) for for-all
  Laplacian cut queries in [*On Sketching Quadratic
  Forms*](https://doi.org/10.1145/2840728.2840753).  For one global objective
  in a different model, Kapralov--Krachun prove that a single-pass streaming
  algorithm beating the factor-two barrier for Max-Cut needs \(\Omega(n)\)
  space in [*An Optimal Space Lower Bound for Approximating
  MAX-CUT*](https://doi.org/10.1145/3313276.3316364).
  These results are precedents, not black boxes in the proof: their
  multiplicative-error, streaming, and for-each/for-all quantifiers differ
  from the additive intervention oracle below.
* [Kullback--Leibler's original information/sufficiency
  paper](https://doi.org/10.1214/aoms/1177729694) and
  [Sanov's primary large-deviation report](https://repository.lib.ncsu.edu/items/8f909775-ba1b-4874-acc2-362a8221edb0)
  identify relative entropy as the information price and rare-event rate.
  [Bretagnolle--Huber's primary minimax paper](https://doi.org/10.1007/BF00535278)
  is a classical testing/event-probability counterpart.
  Section 6 gives an exact finite event bound; no asymptotic large-deviation
  approximation is needed.

Throughout, \(h_2(p)=-p\log_2p-(1-p)\log_2(1-p)\), with the usual endpoint
convention; \(H_2\) and \(D_2\) denote entropy and relative entropy in bits.

## 3. General finite query-sufficiency bounds

### 3.1 Deterministic metric-entropy converse

Let \(\mathcal P\subset\mathcal H\) satisfy
\(\rho_{\mathcal U}(H,H')>2\varepsilon\) for all distinct
\(H,H'\in\mathcal P\).  If an \(R\)-bit deterministic sketch obeys
\(d_\infty(H,z(H))\le\varepsilon\) for every \(H\), then
\[
\boxed{R\ge\log_2|\mathcal P|.}                         \tag{3.1}
\]
Indeed, two members assigned the same message would both lie within
\(\varepsilon\) of the same decoded answer vector, contradicting the
triangle inequality.  Equivalently,
\[
R\ge \log_2 \mathsf{Pack}
   (\mathcal H,\rho_{\mathcal U},2\varepsilon).
\]
This is the worst-case, zero-error face of query rate--distortion.

### 3.2 Bayesian and Blackwell--Le Cam faces

Let \(\Theta\sim\Pi\), let \(Z\) be any transcript, and suppose an action
\(\widehat\Theta(Z)\) has expected distortion at most \(D\).  Then any lower
bound on the ordinary source rate--distortion function
\[
\inf I(\Theta;Z)
\]
is automatically a sketch lower bound, because a fixed-length \(R\)-bit
message satisfies \(I(\Theta;Z)\le H(Z)\le R\).  Public randomness independent
of \(\Theta\) does not change this inequality after conditioning on it.

Let \(\mathsf I\) be the identity experiment that outputs \(\Theta\), and let
\(\mathsf S\) be the sketch experiment.  If
\(\delta(\mathsf S,\mathsf I)\le\eta\), some kernel decodes \(\Theta\) with
worst-case error at most \(\eta\).  Under a uniform prior on \(L\) parameters,
Fano's inequality gives
\[
R\ge I(\Theta;Z)
 \ge \log_2L-h_2(\eta)-\eta\log_2(L-1).              \tag{3.2}
\]
Thus exact Blackwell sufficiency for an injective query-answer map is a
lossless source code.  Approximate value fidelity and Le Cam deficiency are
not identical notions, but a stable decoder transfers the former into the
latter.

## 4. Main theorem: pinned optimum queries have quadratic rate

Let \(n\ge2\), let \(E_n=\{\{i,j\}:1\le i<j\le n\}\), and let
\(N=|E_n|\).  Fix an energy quantum \(a>0\).  For
\(A\in\{\pm1\}^{E_n}\), put
\[
q_A(x)=a\sum_{i<j}A_{ij}x_ix_j,\qquad
c_A=\max_x q_A(x),\qquad
H_A(x)=q_A(x)-c_A.                                    \tag{4.1}
\]
Every member of this \(2^N\)-element class has
\[
\boxed{\max_xH_A(x)=0.}                               \tag{4.2}
\]
For each \(u\in X_n\), query the external field
\(h^u=Mu\), where
\[
M>a(n-1).                                             \tag{4.3}
\]
Let \(U\) denote a uniform spin and
\(\chi_{ij}(u)=u_i u_j\).

> **Theorem 4.1 (finite Ising query rate--distortion).**  Let \(A\) be uniform
> on \(\{\pm1\}^{E_n}\).  From an \(R\)-bit transcript \(Z\), a decoder
> supplies a function \(u\mapsto\widehat V_Z(h^u)\).  Define
> \[
> d_Q(A,Z):={1\over a^2N}\,
>   \mathbb E_U\bigl[
>    \widehat V_Z(h^U)-V_{H_A}(h^U)
>   \bigr]^2.                                         \tag{4.4}
> \]
> If \(\mathbb E d_Q(A,Z)\le D\le1/2\), then
> \[
> \boxed{
> I(A;Z)\ge N[1-h_2(D)],\qquad
> R\ge N[1-h_2(D)].}                                  \tag{4.5}
> \]
> In particular, if the same decoded answer function has uniform additive
> error strictly below \(a\) for every \(A,u\), then \(A\) is recovered
> exactly and \(R\ge N=\binom n2\).

The normalization in (4.4) is deliberately weak: constant \(D\) permits
root-mean-square value error as large as \(a\sqrt{ND}\).  The theorem still
forces a positive quadratic rate.

### Proof

Fix \(A,u\), and let \(x\) differ from \(u\) in \(k\ge1\) coordinates.  The
field loses \(2Mk\).  Only the \(k(n-k)\) edges crossing the flipped set can
change their interaction contribution, and their total possible gain is at
most \(2ak(n-k)\le2ak(n-1)\).  By (4.3), \(u\) is the unique maximizer, hence
\[
V_{H_A}(h^u)=Mn+q_A(u)-c_A.                            \tag{4.6}
\]

Set
\[
f_A(u)=V_{H_A}(h^u)-Mn=q_A(u)-c_A,
\qquad
\widehat f_z(u)=\widehat V_z(h^u)-Mn.
\]
The degree-two Walsh coefficient of \(f_A\) is
\[
\mathbb E_U f_A(U)\chi_{ij}(U)=aA_{ij};               \tag{4.7}
\]
the unknown normalization \(c_A\) is a degree-zero coefficient and
disappears.  Define
\[
\widehat b_{ij}(z)=\mathbb E_U\widehat f_z(U)\chi_{ij}(U),
\qquad
\widehat A_{ij}(z)=\operatorname{sgn}\widehat b_{ij}(z),
\]
with an arbitrary convention at zero.  Bessel's inequality for the Walsh
characters gives
\[
\sum_{i<j}\bigl(\widehat b_{ij}(z)-aA_{ij}\bigr)^2
\le
\mathbb E_U\bigl(\widehat f_z(U)-f_A(U)\bigr)^2.       \tag{4.8}
\]
If the decoded sign on an edge is wrong, the corresponding squared
coefficient error is at least \(a^2\).  Therefore, pointwise in \(A,z\),
\[
{d_{\rm Ham}(A,\widehat A(z))\over N}\le d_Q(A,z).    \tag{4.9}
\]

Let \(p_{ij}=\Pr\{A_{ij}\ne\widehat A_{ij}(Z)\}\).  The mean of the \(p_{ij}\)
is at most \(D\).  Since the \(N\) source bits are independent and uniform,
\[
\begin{aligned}
I(A;Z)
&=N-H(A\mid Z)\\
&\ge N-H(A\mid\widehat A)\\
&\ge N-\sum_{i<j}h_2(p_{ij})\\
&\ge N[1-h_2(D)],
\end{aligned}                                         \tag{4.10}
\]
where the last step is concavity of binary entropy.  This proves (4.5).

For the uniform claim, an error bound
\(\|\widehat f_z-f_A\|_\infty< a\) makes every coefficient error in
(4.7) strictly smaller than \(a\), so every sign is exact.  Equivalently, if
\(A\ne B\), a differing degree-two coefficient gives
\[
\rho_{\mathcal U}(H_A,H_B)
\ge \left|\mathbb E_U(f_A-f_B)\chi_{ij}\right|=2a.    \tag{4.11}
\]
Thus the whole class is \(2a\)-separated, equivalently an `r`-packing for
every `r<2a`, and (3.1) also gives \(R\ge N\).
Storing \(A\) itself uses \(N\) bits and permits exact (computationally
unbounded) decoding, so the lossless rate is exactly \(N\).
\(\square\)

### Experiment-comparison corollary

Let \(\mathsf Q\) be the deterministic experiment that outputs the complete
vector \((V_{H_A}(h^u))_{u\in X_n}\), and let \(\mathsf I\) output \(A\).
Forward evaluation gives a kernel \(\mathsf I\to\mathsf Q\), while (4.7)
gives a deterministic kernel \(\mathsf Q\to\mathsf I\).  Thus
\(\mathsf Q\) and \(\mathsf I\) are Blackwell equivalent on this family.
If a sketch experiment \(\mathsf S\) has
\(\delta(\mathsf S,\mathsf Q)\le\eta\), then (3.2) specializes to
\[
R\ge N-h_2(\eta)-\eta\log_2(2^N-1)
 \ge (1-\eta)N-h_2(\eta).                             \tag{4.12}
\]
This is the Le Cam-deficiency face of the same lower bound.  It is stronger
than asking for one decision problem only because the declared query
experiment is injective on the hard family.

### Bounded coupling interventions give the same decoder

The growing field magnitude in (4.3) is not essential if pair-coupling
queries are allowed.  For \(L>a\), add the rank-one coupling
\[
K^u_{ij}=L u_i u_j.
\]
Writing \(z_i=u_ix_i\), the perturbation loses \(2Lk(n-k)\) when \(z\) is
not constant, whereas \(q_A\) can gain at most \(2ak(n-k)\).  Hence the
maximizers are \(x=\pm u\), and
\[
\max_x\left\{H_A(x)+\sum_{i<j}K^u_{ij}x_ix_j\right\}
=LN+q_A(u)-c_A.                                       \tag{4.13}
\]
The Walsh decoder and every bound above apply unchanged.  Here each added
coupling has magnitude \(L=O(a)\), independent of \(n\); the query family has
only \(2^{n-1}\) distinct members because \(K^u=K^{-u}\).

## 5. State information is linear, not quadratic

For \(s\in\{\pm1\}^n\), define the normalized linear landscape
\[
G_s(x)=a\langle s,x\rangle-an.
\]
Again \(\max_xG_s(x)=0\), now with unique optimizer \(s\).  At the \(n\)
external-field queries \(h^i=ae_i\), direct coordinatewise maximization gives
\[
V_{G_s}(h^i)=a s_i.                                    \tag{5.1}
\]
If
\[
d_{\rm state}(s,z)={1\over a^2n}\sum_i
 \bigl(\widehat V_z(h^i)-V_{G_s}(h^i)\bigr)^2,
\]
then thresholding the answers produces
\[
{d_{\rm Ham}(s,\widehat s(z))\over n}\le d_{\rm state}(s,z).
\]
The same entropy proof yields
\[
\boxed{R\ge n[1-h_2(D)]}                              \tag{5.2}
\]
under uniform \(s\) and expected distortion \(D\le1/2\).  This is the
\(\Theta(n)\) state/extremizer scale.

The distinction remains sharp inside the quadratic class.  Choose one
labelled optimizer \(x^*(A)\) for every \(A\).  By pigeonhole, some common
optimizer fiber contains at least \(2^{N-n}\) different matrices.  On that
fiber both \(\max H_A=0\) and \(x^*(A)\) are identical, yet (4.11) separates
every pair of full intervention-response vectors.  Retaining a maximum and
one maximizer can therefore leave \(N-n=\Theta(n^2)\) counterfactual bits
unresolved.

## 6. Rare events and the repository's sign-near entropy boundary

### 6.1 A general rare-success information lemma

Let \(\Theta\sim\Pi\), and suppose each transcript \(z\) specifies a success
set \(G_z\) with \(\Pi(G_z)\le p\).  If
\(\Pr\{\Theta\in G_Z\}\ge1-\delta\), \(0\le\delta\le1/2\), then
\[
\boxed{
I(\Theta;Z)\ge(1-\delta)\log_2{1\over p}-h_2(\delta).} \tag{6.1}
\]
To prove it, compare the joint law \(P_{\Theta Z}\) with
\(\Pi\otimes P_Z\) and apply data processing to the indicator
\(1_{\{\Theta\in G_Z\}}\).  Under the product law the event has probability at
most \(p\); binary relative entropy is at least
\((1-\delta)\log_2(1/p)-h_2(\delta)\).  Thus an algorithm that reliably lands
in an exponentially rare good set must acquire a proportional number of
bits, independent of how that good set was found.

### 6.2 Exact posterior-polarization converse

The strongest information ingredient in
[the repository's sign-near weighted-recovery note](../../artifacts/sign_near_weighted_recovery.md)
has a direct query-sketch interpretation.  Let \(A\) be uniform on its \(N\)
edge signs, and for a transcript \(Z=z\) put
\[
w_e(z)=\mathbb E[A_e\mid Z=z],\qquad
V_z=\sum_e(1-w_e(z)^2).                                \tag{6.2}
\]
For each posterior \(\mu_z\), entropy subadditivity, monotonicity of \(h_2\),
and Jensen give
\[
H_2(\mu_z)
\le \sum_e h_2\!\left({1-|w_e(z)|\over2}\right)
\le N h_2\!\left({V_z\over2N}\right).                \tag{6.3}
\]
Averaging once more and using concavity yields the nonasymptotic converse
\[
\boxed{
I(A;Z)\ge
N\left[1-h_2\!\left({\mathbb E V_Z\over2N}\right)\right].} \tag{6.4}
\]
Consequently, a transcript whose posterior barycenter is globally sign-near,
\(\mathbb E V_Z=o(N)\), has rate \(N-o(N)=\Theta(n^2)\).  This makes precise
the statement in [ACTIVE_STATE.md](../../ACTIVE_STATE.md) that sign-near
weighted recovery exposes almost the entire sign skeleton.  It also imports
the note's sharp bound rather than the coarser Pinsker estimate.

There is an exact rare-event corollary.  If \(U\) is uniform on edge signings,
\(E\) is nonempty, \(\mu=U(\cdot\mid E)\), \(w_e=\mathbb E_\mu A_e\), and
\(V=\sum_e(1-w_e^2)\), then
\[
-\log_2U(E)=D_2(\mu\|U)
\ge N\left[1-h_2\!\left({V\over2N}\right)\right].    \tag{6.5}
\]
Thus \(V=o(N)\) forces
\[
U(E)\le2^{-N+o(N)}.                                   \tag{6.6}
\]
Near-polarization is not a low-cost conditioning trick: at the edge scale it
is nearly a singleton rare event.

The repository note also proves that \(V=o(N)\) permits biased rounding with
\(o(n^{3/2})\) Boolean-objective error and vanishing normalized bilinear/action
error.  Combining that rounding fact with (6.4) gives a clean separation:
the **rounding computation** can be cheap once a sign-near barycenter exists,
but a generic transcript that identifies such a barycenter is
information-heavy.  This observation is a boundary on recovery schemes, not
a proposed route to the asymptotics of \(M_n\).

## 7. Finite falsifier

The proof is algebraic, but its pinning and Fourier identities have a complete
finite falsifier.  The following dependency-free program enumerates every
landscape through \(n=5\), every pinning field, and every Walsh-recovered edge.

```python
from itertools import product, combinations

for n in range(2, 6):
    edges = list(combinations(range(n), 2))
    N, a, M = len(edges), 1, n       # M > a(n-1)
    X = list(product((-1, 1), repeat=n))
    seen = set()
    for A in product((-1, 1), repeat=N):
        def q(x):
            return sum(A[k]*x[i]*x[j]
                       for k, (i, j) in enumerate(edges))
        c = max(q(x) for x in X)
        V = {u: max(q(x)-c+M*sum(ui*xi for ui, xi in zip(u, x))
                    for x in X)
             for u in X}
        assert all(V[u] == M*n + q(u) - c for u in X)
        recovered = []
        for i, j in edges:
            b = sum((V[u]-M*n)*u[i]*u[j] for u in X)/(2**n)
            recovered.append(1 if b > 0 else -1)
        assert tuple(recovered) == A
        seen.add(tuple(V[u] for u in X))
    assert len(seen) == 2**N
    print(n, 2**N, "landscapes checked")
```

The executed output was

```text
2 2 landscapes checked
3 8 landscapes checked
4 64 landscapes checked
5 1024 landscapes checked
```

Any failed assertion falsifies a specific step of Theorem 4.1.  The field
threshold is substantive, not decorative: for \(n=3,a=1,M=1\), taking all
three couplings \(-1\) and \(u=(-1,-1,-1)\) already violates (4.6).  At
\(M=2=a(n-1)\) the value identity holds (possibly with ties); the strict
hypothesis was chosen to guarantee a unique pinned maximizer.

## 8. Application outside \(M_n\): vertex-prize Max-Cut oracles

Let \(B\in\{0,1\}^{E_n}\) encode a dense unweighted graph and
\[
C_B(x)=\sum_{i<j}B_{ij}{1-x_ix_j\over2}
\]
be its cut value.  Consider the intervention oracle
\[
W_B(u)=\max_x\{C_B(x)+M\langle u,x\rangle\},
\qquad M>{n-1\over2}.                                  \tag{8.1}
\]
The same flip-set argument pins \(x=u\), so
\[
W_B(u)=Mn+C_B(u),qquad
\mathbb E_U[W_B(U)-Mn]\,U_iU_j=-{B_{ij}\over2}.       \tag{8.2}
\]
Therefore a single sketch that answers **all** vertex-prize Max-Cut optimum
queries with additive error \(<1/4\) identifies every edge by thresholding
the Walsh coefficient.  It must use
\[
\boxed{\binom n2\text{ bits in the worst case}.}       \tag{8.3}
\]
This is not the usual problem of approximating one Max-Cut value, nor the
usual relative-error all-cut sketch: the query asks for an optimum after a
strong but succinctly described vertex intervention.  It is a concrete
inverse-optimization statement: uniform counterfactual optimum access is a
lossless graph description.  This one-sided nonnegative cut landscape is
intentionally distinct from the absolute signed quadratic objective defining
\(M_n\).

## 9. Scope and sharp boundaries

* The quadratic lower bound is information-theoretic and finite.  It says
  nothing about encoding/decoding time.
* The query family is exponentially large but each query has an \(n\)-bit
  sign description.  Uniform reconstruction, or the mean-square fidelity in
  (4.4) for one common decoded answer function, is essential.  A guarantee
  for one preselected query is a different model.
* External-field pinning uses \(M=\Theta(an)\).  If only bounded external
  fields are allowed, this proof does not apply.  Rank-one coupling
  interventions (4.13) recover the result with bounded per-edge magnitude.
* The additive threshold tracks the coupling quantum.  Coarser additive or
  purely multiplicative error can merge landscapes; one should then compute
  the packing or rate--distortion function of the merged response class
  rather than quote (4.5) unchanged.
* The posterior result (6.4) concerns information needed to identify a
  near-polarized edge phase.  It does not say that every useful summary of a
  landscape has quadratic entropy; Section 5 gives an explicit linear-rate
  task.
* Nothing here compares different orders, constructs extremal signings, or
  proves a limit.  In particular, the report deliberately makes no claim
  about convergence or nonconvergence of \(M_n/n^{3/2}\).
