# A high-original-degree probe of the actual first marked feedback

2026-09-06. Independently audited at the fixed-first-history, fixed-L
scope in `transfer_adversary_high_degree_actual_feedback_audit_2026_09_06.md`.
The finite Gaussian rank
inequality is separately proved in
`transfer_director_inhomogeneous_gaussian_return_rigidity_2026_09_06.md`.
This note supplies an actual mixed comparison and its quantitative
consequence. Its scope is the ALREADY AUDITED first marked history, not
an arbitrary later computation and not yet a better unrestricted decimal.

## 1. Frozen actual fields and the fixed high-degree choice

Let `B=A/sqrt(n-1)` be an actual symmetric hollow signing, with fixed
`||B||op<=L`, and `Q=B^2`. On independent uniform signs S retain literally

```
G=BS, D_i=S_i h_2(G_i), Y=BD,
W=(S,G,Y,QS,QD).
```

Let `f(g,y)` be fixed, odd, Gaussian-a.e.-continuous and ternary, with
`0<E f(N_0,N_1)^2<1`. Set `H=1-f^2`, `mu=E H>0`, and

```
f=b_0 g+b_1 y+r,
V=b_0 QS+b_1 QD, Z=Br(G,Y), BF=V+Z, F=f(G,Y),
C=H(G,Y) sign(V+Z).
```

The sign at zero can be fixed arbitrarily to be Boolean. Its choice is
irrelevant in averaged limits by the audited small-ball bound. In
particular `u_+=F+C` and `u_-=-F+C` are actual Boolean vectors.
With, for example, sign(0)=+1, C need NOT be exactly odd or centered at
finite n. The probe U_k and E_k are exactly odd and centered. Thus each
mixed covariance with the probe equals the uncentered product expectation
used below. The hard-tie modification costs o(1) in averaged L2 relative
to the odd sign(0)=0 comparison; no exact centering of C is required.

Expand r in normalized local Hermites. Give `h_a(g)h_b(y)` ORIGINAL
degree `a+3b`. Since a bounded nonconstant function cannot be a finite
Gaussian polynomial, there is a FIXED odd `P>3` for which

```
r_P=sum_{a+3b=P} r_ab h_a(g)h_b(y),
tau_P^2=sum_{a+3b=P} r_ab^2>0.
```

Only terms with local degree `a+b>=3` occur here. Choose this P once,
before taking n to infinity. No uniform lower bound over an increasing
primitive cutoff is asserted. Let `Z_P` be the transported exact
injective original-P forest main of `B r_P(G,Y)`. At fixed P the raw
replacement has vanishing averaged L2 error. Define

```
R_P=sum_{a+3b=P} r_ab^2 Q^{circ(a+b)},
T_P=B R_P B, q_i=(T_P)_ii,
T=B[sum_{a+b>=3 odd} r_ab^2 Q^{circ(a+b)}]B, sigma_i^2=T_ii.
```

The established nonlinear covariance theorem, applied to each fixed
original-degree block before summation, gives normalized nuclear errors
`o(n)`. Its exact Schur main gives

```
0 <= T_P <= L^2 tau_P^2 Q,
q_i <= L^2 tau_P^2, sigma_i^2 <= L^2 ||r||_2^2 <= L^2,
#{i:q_i<=epsilon tau_P^2} <= epsilon n.                 (1)
```

The last statement is the exact odd-Schur subset bound in
`continued_feedback_threshold_without_variance_floor_2026_09_06.md`.
Also `n^-1 E||V||^2 <= L^2+o(1)`. No law of QS or QD is Gaussianized.

## 2. Actual mixed-covariance lemma

Fix odd `k>=3` and a fixed variance cutoff `epsilon>0`. Put
`J={i:q_i>epsilon tau_P^2}`. With an auxiliary Gaussian independent of
the ACTUAL W, define the bounded odd increasing scalar function

```
g_i(t)=E_{S,N}[H(G_i,Y_i)
       sign(V_i+sqrt(q_i)t+sqrt(sigma_i^2-q_i)N)],
alpha_ik=E_G[g_i(G)h_k(G)],
d_i=1_J(i) alpha_ik/q_i^(k/2), D_k=diag(d_i),
U_ki=q_i^(k/2) h_k(Z_Pi/sqrt(q_i)), E_k=B D_k U_k,
C^J=1_J C.
```

Here h_k is the ORTHONORMAL probabilists' Hermite polynomial. Equivalently
U is the standard variance-Hermite `He_k(Z_Pi;q_i)/sqrt(k!)`, where
`He_k(z;q)=q^(k/2) He_k(z/sqrt(q))` on the right uses the unnormalized
probabilists' polynomial.
The exact degree-kP squarefree main can equivalently be used for U_k;
the difference is o(1) in averaged L2 at these fixed parameters.
The coefficients d_i are deterministic and uniformly bounded at fixed
epsilon, P and k. The lemma is

```
|| Cov(U_k)-T_P^{circ k} ||_*/n -> 0,
|| Cov(C^J,U_k)-D_k T_P^{circ k} ||_*/n -> 0.             (2)
```

Thus, with every occurrence of B and C actual,

```
n^-1 E[(BC^J)^T E_k]
 = n^-1 tr(Q D_k T_P^{circ k} D_k)+o(1),
n^-1 E||E_k||^2 <= L^2+o(1).                            (3)
```

The ideal row variances of E_k are bounded by the fixed number
`M_E=L^(2k+2) epsilon^(-k)`; their difference from actual variances is
o(1) in averaged absolute value. This follows by Schur multiplication,
`||T_P||op<=L^4 tau_P^2`, and `|alpha_ik|<=1`. The averaged bound in
(3) is much better than this row cap: its source has ideal diagonal
variance `1_J alpha_ik^2<=1`.

### 2.1 Fixed-stage diagram proof of (2)

Use the exact local surgery and coherent equality bookkeeping of
`continued_feedback_first_marked_history_energy_projection_2026_09_06.md`,
Sections 2--5. At a fixed polynomial stage, expand the LEFT response in
variance-Hermites of its finitely many transported residual channels,
with polynomial coefficients in the primitive list W. Its internal
coherent Boolean contractions remain exact. Remove only noise-touching
local collisions and the noise-noise contractions canceled by normal
ordering. Those removals are o(1) in averaged L2.

The RIGHT kernel is k independent formal copies of the SAME original-P
noise, projected onto distinct seed labels. Every such noise has all
proper fixed-root cuts `O(n^-1/2 polylog n)` and all global cuts bounded
by a fixed polylogarithm. Every coherent primitive has at most THREE
original slots, strictly fewer than P. This strict inequality is the
reason for choosing P>3; orthogonality alone would not suffice.

In a surviving two-root pairing, an isolated whole pair between right
noise and a left noise is its exact channel covariance. A noise meeting
several primitive vertices has a proper-cut gain; a retained Boolean
hyperedge has the corresponding fixed-slot influence gain. Expand
coherent coefficients into their primitive diagrams BEFORE making this
claim. The hypergraph first-merge inequality from
`continued_feedback_marked_local_noise_separation_2026_09_06.md`,
Section 4, keeps internal coherent collisions and all open labels.

The nonleading diagrams have the following exhaustive disposition.

* If there are two vertex-disjoint partial merges, perform them first.
  Each costs `n^-1/2 polylog n`; all later merges are Hilbert-contractive.
  The resulting entry is `O(n^-1 polylog n)`, hence its two-root matrix
  has Frobenius norm `O(polylog n)` and nuclear norm o(n).
* If a single nontrivial star is accompanied by an isolated whole-noise
  edge, retain that edge as a covariance matrix of bounded operator
  norm. The star is entrywise `O(n^-1/2 polylog n)`, while the retained
  covariance has Frobenius norm `O(sqrt(n))`. Their entrywise product
  again has Frobenius norm `O(polylog n)`.
* A star consuming ALL k right branches cannot have an old coherent
  primitive as its center: that primitive has at most three slots.
  Its only remaining center is one left transported noise. Pull this
  center back through its B to its original old G/Y source. The k right
  noises must then meet primitive tensors of degree at most three,
  either in this old source or in the left coherent coefficient.
  Hall matching assigns them k distinct primitive neighbors, because
  k' selected right noises have P k' slots and a primitive has at most
  three. Every assigned merge is proper or has a retained-label
  influence gain. The source contraction is `O(n^(-k/2) polylog n)`
  uniformly in the three roots. Summing the pulled-back B row costs at
  most sqrt(n). Since k>=3, its final entry is `O(n^-1 polylog n)`.

The isolated-edge factorization in the second bullet is applied only to
an UNRESTRICTED equality diagram after finite inclusion-exclusion, not
to the original globally distinct source sum. Its correction accounting
is as follows. Identifying a whole-pair label with a different component
forces that label to be RETAINED when the whole pair is merged. For
squarefree noise kernels K,K', the resulting slice tensor has Hilbert
norm at most the maximum fixed-slot influence square root of K times
the Hilbert norm of K'. This supplies a new `n^-1/2 polylog n` gain.
If a whole pair is joined to the star, use this gain and the star's
original proper merge on disjoint primitive vertices. If two whole
components are joined, both retained-pair merges supply gains. An
identification of two distinct slots within any primitive vanishes.
In particular two different labels of a star center cannot be merged
to eliminate its proper cut: that correction is zero. If an unaffected
whole component remains, its covariance still factors and supplies the
Frobenius sqrt(n) bound. Repeating this finite classification handles
several identifications; no distinctness restriction is silently dropped
from a factored main term.

These cases also cover a collection of stars: two nontrivial components
give two disjoint gains. A connected bipartite graph that is not an
edge or a star contains a P4, whose outer edges give the two disjoint
merges. A non-pair Boolean hyperedge is handled by a retained-label
merge; any equality between distinct declared free right labels is
excluded, and its inclusion-exclusion corrections obey the same bounds.
In particular a polynomial in W is not replaced by an independent
Gaussian coefficient.

The only remaining diagrams match each of the k right noises wholly
to a left original-P noise, and leave the coherent coefficient at its
ACTUAL mean. Their sum is precisely `D_k T_P^{circ k}`, with the stated
normalization of He_k. For Cov(U_k), the two sides have the same number
and degree of branches; a lone all-consuming star is impossible, and
the same proof gives `T_P^{circ k}`.

First use exact main channel covariances and their exact row variances.
Replacing them by T_P uses their normalized nuclear comparison and the
bounded Gram Schur-multiplier inequality; it does NOT require operator
convergence of the raw covariance. Row-variance replacement is controlled
in averaged absolute value. Bounded d_i at the fixed cutoff permits it.

### 2.2 Ordered bounded and hard-threshold passage

For this primitive list, QS is subgaussian and QD has the already audited
uniform exponential moment; see the marked local-noise proof, Section 7.
Thus the existing finite-catalog polynomial approximation of bounded
functions of the LITERAL W applies. Approximate the bounded old f,H and
the response sign in the audited order: fixed Hermite/response stage,
then n to infinity, then approximation error to zero. The exact odd-Schur
inverse-variance small-ball bound removes the hard threshold. On J there
is also the simpler fixed Gaussian variance floor.

For random vectors X,Y, `||E XY^T||_* <= (E||X||^2 E||Y||^2)^(1/2)`.
Consequently every o(1) averaged L2 replacement made here is o(n) in
the required nuclear norm. The fixed cutoff keeps all deterministic
coefficients bounded. No unbounded diagonal multiplier is applied to
a raw covariance error. This proves the ordered version of (2).

## 3. A uniformly positive FULL return correlation

Let `gbar(t)=n^-1 sum_i sqrt(q_i) g_i(t)`. Its amplitude is at most
`M0=L tau_P`, and it is odd and nondecreasing. Its first coefficient is

```
<gbar,h_1>=n^-1 sum_i 2 q_i E[H phi(V_i/sigma_i)/sigma_i].
```

Set `q0=tau_P^2/2` and `R0^2=8 L^2/mu`. By (1), at least half the
rows have q_i>q0. The old local marginal gives `E H_i=mu+o(1)` uniformly;
the averaged V second moment bounds the H-mass with |V|>R0 by mu/8+o(1).
Since `sqrt(q0)<=sigma_i<=L` on the retained rows, for all large n

```
<gbar,h_1> >= c0 := mu q0 phi(R0/sqrt(q0))/(2L)>0.       (4)
```

Apply the director's finite Gaussian rigidity theorem with c0,M0. Put
`delta=E(c0|N|-M0)_+^2`, and fix
`K0>=max(2,36 M0^4/delta^2)`. There is some odd `3<=k<=K0` such that

```
|<gbar,h_k>| >= a0:=sqrt(delta/(2K0)).
```

The choice may depend on n, but belongs to this FIXED finite set. All
comparison errors in (2) are therefore uniform over the needed choices.
Removing rows q_i<=epsilon tau_P^2 changes any coefficient of gbar by
at most `tau_P epsilon^(3/2)`. Choose epsilon once so that

```
tau_P epsilon^(3/2)<=a0/2,
L^2 sqrt(epsilon)<=a0^2/(8 L^2 tau_P^2).                (5)
```

The signed-diagonal rank identity, with `Q>=T_P/(L^2 tau_P^2)`, now gives

```
n^-1 E[(BC^J)^T E_k] >= a0^2/(4 L^2 tau_P^2)+o(1).
```

This initially concerns C^J, not the full C. The distinction matters.
But (3), `||B||op<=L`, `|C|<=1`, and |J^c|<=epsilon n imply

```
n^-1 |E[(B(C-C^J))^T E_k]| <= L^2 sqrt(epsilon)+o(1).
```

Hence the ACTUAL unmodified return obeys the closed quantitative test

```
n^-1 E[(BC)^T E_k] >= a_*+o(1),
a_*:=a0^2/(8 L^2 tau_P^2)>0.                           (6)
```

The first coefficient, P, K0 and epsilon all precede the signing limit.
This is stronger than saying that every finite cutoff leaves some tail:
it supplies one finite tested degree with a fixed positive correlation.

## 4. Exact local queries needed for an energy gain

The next step must compare E_k with the following literal query list:

```
W, every fixed-degree old nonlinear channel of Z,
BC-beta_i E_ki,
beta_i=E[(BC)_i E_ki]/E[E_ki^2].
```

The beta coefficients require a positive probe-variance cutoff and a
root second-moment cap, as in the audited ordered next-gain theorem.
Primitive orthogonality by itself is NOT sufficient for this list.

Here is the fixed-stage full-contraction check for that list.
The probe's exact main has original degree kP and small proper cuts by
the two-factor flat-transport theorem. Higher degrees of an old channel
Z have small proper cuts on the RIGHT. Equal-degree probe/old-channel
covariances vanish after pulling the old channel to its G/Y source and
Hall-matching k probe branches to distinct degree-at-most-three source
primitives. Lower degrees are handled by proper cuts of the probe.

For the returned response decompose at a fixed stage

```
BC=Bc0(W)+B[(A(W)-a)Z+R_{>=2}]+BD_aB r(G,Y).
```

Higher degrees of Bc0 have small proper cuts by the exact boundary-graph
theorem, since those degrees exceed kP>3. Each middle term has at least
two rooted factors and thus small proper cuts after B transport. For
the last term, full contraction of the probe into a larger original
degree is pulled back on BOTH output-root axes. Open-mark Hall matching
gives k disjoint proper/influence gains against its old G/Y source.
Its two-root open tensor has Frobenius norm `O(n^(1-k/2) polylog n)`;
bounded outer operators followed by the diagonal-root projection keep
it o(sqrt(n)). Equal-degree covariance is retained by beta, not discarded.
All positive coherent equality diagrams have the boundary-graph global
cut bounds, including diagrams with internal Boolean collisions.

These are exactly the proper/full derivative contraction conditions of
`continued_audit_boolean_stable_noise_stein_2026_09_06.md`. Applying its
characteristic-function argument retains the actual joint old-query law
and makes the probe Gaussian independent of the REGRESSED query list.
The bounded passage uses fixed-degree approximations, root caps, a fixed
probe-variance floor and then ordered removal, exactly as in
`continued_audit_ordered_bounded_next_gain_2026_09_06.md`. It does not
assert a Gaussian law for BC-beta E, QS or QD.

Equation (6) alone is not yet a spin-flip gain: correlation and Gaussian
marginals would not suffice without this exact literal-query comparison.
The following section records the exact implication of that comparison.

## 5. Quantitative actual cap consequence of the local comparison

Set `Lambda(B)=max_{x in {-1,1}^n}|x^T Bx|/(2n)` and
`j_n=n^-1 E sum_i H_i |(BF)_i|`. The mixed and local comparisons in
Sections 2 and 4 imply a constant `Delta(f,L)>0` such that

```
liminf_n [Lambda(B_n)-j_n] >= Delta(f,L).                 (7)
```

Here is an explicit, deliberately conservative choice. Use a_* from (6),
choose the k from the finite set there, and replace the probe row cap by
the maximum `M_E=max_{3<=k<=K0 odd} L^(2k+2) epsilon^(-k)`. Put

```
A0=a_*/(2 sqrt(2 M_E)), s0=A0/2,
theta=A0^2/(4 L^2), R1=sqrt(2)*L/sqrt(theta),
z0=s0 phi(R1/s0)-R1 barPhi(R1/s0)>0,
w0=theta*z0/4, Delta(f,L)=w0^2/(2L).                    (8)
```

The positivity of z0 follows either from strict Mills' inequality or
from `z0=E[-R1-s0 N]_+`.

To prove (7), take either Boolean endpoint `u=+F+C` or `u=-F+C`.
The local comparison gives vanishing average covariance of the probe
with BF, so (6) applies to Bu as well. Write

```
v_i=E E_ki^2, b_i=E[(Bu)_i E_ki],
s_i=|b_i|/sqrt(v_i), R_i=(Bu)_i-(b_i/v_i) E_ki.
```

Cut first to `v_i>nu=(a_*/(4L))^2` and `v_i<=2 M_E`.
The omitted low-variance contribution to the average |b_i| is at most
`L sqrt(nu)=a_*/4`, by Cauchy--Schwarz and
`n^-1 E||Bu||^2<=L^2`. The high-variance contribution is o(1): the
ideal variances are at most M_E and their averaged absolute error is
o(1), so the total actual variance on those rows is o(n).
Consequently the retained roots have average s_i at least A0 for all
sufficiently large n, while their average s_i^2 is at most L^2 by
Bessel's inequality. At least a fraction theta of all roots therefore
has `s_i>=s0`.

On each retained root, Section 4 makes E_ki independent Gaussian of
variance v_i relative to the LITERAL joint law of `(u_i,R_i)`.
Regression gives `n^-1 sum_i E R_i^2<=L^2`. Hence the event consisting
of `s_i>=s0` and `|R_i|<=R1`, with a uniform random root, has probability
at least theta/2. For either FIXED energy orientation e in {-1,1}, the
conditional negative local field has expectation

```
E_N[-e u_i(R_i+(b_i/sqrt(v_i))N)]_+ >= z0
```

on this event. The sign e u_i only changes the Gaussian orientation.
The comparison applies at fixed variance and root-moment cutoffs; remove
the root caps afterwards. Their contribution to the relevant square-root
coefficients is O(R^(-1/2)), using the averaged second moments, exactly
as in the audited ordered-gain proof. Nonnegative omitted negative-field
contributions can be discarded. Uniform integrability for the remaining
linear-growth tests follows from the same averaged second moments.
Thus for BOTH endpoints and BOTH energy orientations,

```
liminf_n n^-1 E sum_i [-e u_i(Bu)_i]_+ >= w0.            (9)
```

The final step is an exact finite signing inequality, with no comparison.
For any Boolean u define

```
w(S)=n^-1 sum_i [-e u_i(Bu)_i]_+,
d_i(S)=-u_i 1{e u_i(Bu)_i<0}.
```

For `0<=alpha<=1`, the vector `u+alpha d` belongs to the cube, and

```
e[(u+alpha d)^T B(u+alpha d)-u^T Bu]/(2n)
 >= alpha w(S)-L alpha^2/2.
```

The random vector need not be Boolean: hollowness makes its quadratic
energy a multilinear cube functional, whose maximum absolute value is
attained at a Boolean vertex. Choose the deterministic
`alpha=E w/L<=1`. It improves expected oriented energy by at least
`(E w)^2/(2L)`.

Finally the two endpoints have expected energies
`e_n(F)+e_n(C) +/- j_n`. One of the four FIXED choices of endpoint and
energy orientation therefore has expected oriented energy at least j_n.
Apply (9) and the finite ascent to that choice. This proves (7).
There is no seed-dependent choice of the initial energy orientation in
this argument, and no replacement of a signed energy by its expectation
inside an absolute value.

## 6. Scope barriers kept explicit

The bounded approximation in Section 2.2 uses the exponential moments of
this particular first marked coherent frame. For a general rich finite
tree frame, a high-degree literal return need not have such an exponential
moment. A primitive-degree inequality alone does not prove uniform
polynomial approximation of bounded functions of that return. Therefore
no arbitrary-rich-frame extension is being silently imported here.

Despite the audited Section 4 comparison, a positive gain depending on fixed L
does not automatically improve the unrestricted lower bound: the
spectral-core deletion loses O(1/L). The constants in (4)--(6) may be far
smaller. No numerical improvement, original scale transfer, convergence,
or nonconvergence follows from that unverified comparison of rates.

## 7. Reproducible finite checks

`computations/transfer_seed_high_degree_feedback_checks_2026_09_06.py`
checks the signed rank inequality, Gaussian mixed-Hermite normalization
with a two-atom NON-Gaussian coherent shift, the edge/star/P4 graph
dichotomy, and the exact feasible spin-ascent algebra on complete finite
Boolean cubes. The run passed 196 signed rank instances, 3240 mixed
Gaussian instances, 2452 bipartite graphs and 72 actual signing/source
instances (70 had a strictly positive finite ascent).

Those checks do not certify the asymptotic Boolean mixed comparison;
its tensor and ordered-limit arguments must be audited mathematically.
