# Independent audit of the first marked-history energy closure

Date: 2026-09-06. The full theorem in
`continued_feedback_first_marked_history_energy_projection_2026_09_06.md`
has been read and independently verified, including its bounded closure.
The director independently reconstructed this audit's two final crosses.
This note concerns one marked history, not an unrestricted-depth
recursion or a convergence proof.

## 1. Setup and normalization

Let B be a normalized hollow symmetric signing with bounded operator
norm, Q=B^2, and independent Boolean seed S. Put

```math
G=BS,\qquad D_i=S_i h_2(G_i),\qquad Y=BD,\qquad
W_i=(S_i,G_i,Y_i,(QS)_i,(QD)_i).
```

The literal coherent fields QS and QD are retained; neither is replaced
by a Gaussian or by a local conditional expectation. Write a fixed odd
old-frame response as

```math
f(g,y)=b_0g+b_1y+r(g,y),\qquad
r=\sum_{p+q\ge3\ {m odd}}r_{pq}h_p(g)h_q(y).
```

For fixed finite polynomial work let

```math
Z=Br(G,Y),\quad
T=B\left[\sum_{p,q}r_{pq}^2Q^{\circ(p+q)}\right]B,
\quad \sigma_i^2=T_{ii},\quad V=b_0QS+b_1QD.
```

Take bounded even H and an odd smooth bounded response psi. The intended
energy identity is, with e(X)=E[X^T B X]/(2n),

```math
C_i=H(G_i,Y_i)\psi(V_i+Z_i),
\qquad c_{0,i}=H(G_i,Y_i)\mathbb E_N\psi(V_i+\sigma_iN),
```

```math
a_i=\mathbb E_{S,N}
 H(G_i,Y_i)\psi'(V_i+\sigma_iN),
```

```math
e(C)=e(c_0)+\frac1n\mathbb E[c_0^TBD_aZ]
 +\frac1{2n}\operatorname{Tr}(BD_aTD_a)+o(1).       (1)
```

The cross coefficient is 1/n, while each self-energy has denominator
2n. The deterministic diagonal D_a cannot be silently replaced by a
constant. Equation (1) is an energy projection, not a covariance-matrix
or random-vector approximation.

## 2. Structural estimates checked independently

The complete file
`continued_feedback_all_global_walsh_cuts_2026_09_06.md` has been read.
Its exact higher Boolean derivative cover formula is correct: the
alternating subset sum retains covers of the differentiated label set;
the cover makes the repeated-input embedding injective; accumulated
seed characters are column signs, not uncontrolled row factors. Fourier
orthogonality identifies every global Walsh cut, up to fixed factorials.
Consequently fixed centered polynomial circuits have polylogarithmically
bounded global cuts. All exact internal coherent Walsh contractions stay
inside those circuits.

The director's two-factor transport lemma was reconstructed in Section 7
of `continued_audit_transported_tensor_powers_2026_09_06.md`: if at least
two rooted input tensors have bounded global cuts, every proper local
cut after transport by a flat row of B is at most

```math
\max_j|B_{ij}|\ \prod_b C_b.
```

Opposite root orientations give an exact composition with a middle
diagonal B-row; no partial trace is taken. Exact distinctness is a finite
composition of rectangular pinching complements, so it enlarges these
cut bounds by only a degree-dependent constant.

The old residual source forest has bounded absolute covariance row and
column sums. In particular its covariance followed by B has entries
O(n^-1/2). The exceptional complete degree-three comparison is also
small: the only new residual of original degree three is transported
h_3(G), and its full covariance with old Y has operator norm
O(n^-1/2). These facts were checked before using the argument below.

Write epsilon=n^-1/2, suppressing fixed powers of log(n) when global
coherent cuts occur. The local noise/coherent surgery is performed FIRST:
remove contractions touching a new-noise branch, retaining its whole
Gaussian pairings and charging the rest in averaged L2. All internal
coherent Walsh contractions remain exact. The previously audited
first-merge argument uses old primitive degrees at most three; its sole
equal-degree exception is the explicit h_3(G)/marked-input covariance.
Global total-degree counts alone do not justify this surgery.

## 3. Local Gaussian testing after a multi-factor transport

Consider an exact squarefree main tensor X made from at least two
rooted factors, with at least two noise factors or with one noise factor
and a positive-degree centered coherent coefficient. Its transport BX
has proper cuts epsilon times a fixed logarithmic power. Its original
degree is greater than three in the applications below.

The local stable-noise lemma therefore compares (BX,Z,W) with (L,N,W),
where (L,N) is jointly centered Gaussian and independent of the ACTUAL
Boolean W. Their full mutual covariance is retained. Mixed original
noise degrees are not declared orthogonal; only their actual full
covariances determine the Gaussian comparison. This statement follows
either from the audited contraction/first-merge moment proof or the
self-contained stable Boolean Stein module.

Transported row variances need not be uniformly bounded. Their average
is bounded by the bounded-op transport and the bounded average source
variance. A fixed variance cutoff followed by its removal handles this
point in averaged tests; alternatively all fixed polynomial comparison
errors are epsilon times logarithmic powers. No uniform row-variance
claim follows from global-cut bounds alone.

## 4. Higher conditional-noise terms have no energy contribution

At fixed polynomial approximation, expand using unnormalized Gaussian
Wick polynomials at variance sigma_i^2:

```math
C=c_0+A(W)Z+R_{\ge2},\qquad
A_i(W)=H(G_i,Y_i)\mathbb E_N\psi'(V_i+\sigma_iN).
```

Thus E A_i=a_i. Each term of R_{>=2} has at least two noise branches
after the local Wick surgery. For every such term X, Section 3 makes
BX jointly Gaussian with Z and independent of W. A centered Gaussian
linear variable has zero expectation against every Wick polynomial in
Z of degree at least two. It follows that

```math
e(R_{\ge2})=o(1),\qquad
\mathbb E[R_{\ge2}^TBc_0]/n=o(1).
```

The remaining test against A(W)Z reduces, by the same comparison, to
the mean a_i times Cov((BX)_i,Z_i). Pulling Z back to its old source
therefore leaves a deterministic bounded-op matrix M=BD_aB paired with
Cov(X_i,r(G_j,Y_j)). This covariance is O(epsilon^2 polylog(n)) entrywise.

Here is a precise two-gain proof. After the internal collision surgery,
the covariance diagram is bipartite: every left new-noise branch has
at least three marked slots, and every right old primitive G or Y has
at most three. Any k left noise branches therefore have at least k
right neighbors, by slot counting. Hall's condition gives distinct
matched old primitives for all the left noise branches. A chosen match
gives a small proper-cut contraction unless it is a full degree-three
pair; that isolated exceptional pair has the already established small
h_3(G)/Y covariance. For two chosen noise branches the matched merges
are vertex-disjoint, so the two small factors multiply. All remaining
merges are Hilbert-contractive, even when matched vertices have edges
to other vertices. Consequently at least two noise branches give the
claimed epsilon^2 bound. Since each row of M has absolute sum O(sqrt(n)),
the normalized weighted covariance tends to zero.

Hall is being applied only AFTER surgery. It cannot be inferred from
the degree of an unexpanded coherent polynomial with internal collisions.

## 5. The centered random linear coefficient

Put Atilde=A-a. It is even and has only positive even Walsh degrees.
The exact product main of Atilde Z has original degree at least five.
The two-factor transport lemma and Section 3 imply

```math
e(Atilde Z)=o(1),\qquad
\mathbb E[(Atilde Z)^TBc_0]/n=o(1).
```

For the self-energy, the Gaussian comparison gives a covariance times
E Atilde_i=0, not a claim that the transported variable and Z have zero
covariance.

It remains to bound its cross with aZ. Pull back the right Z and put
M=BD_aB again. Expand
Cov(Atilde_i Z_i,r(G_j,Y_j)) in old G/Y primitive branches on the right.

If some old Y splits between the left coherent Atilde blob and the
left noise, merge that Y first with Atilde. Its proper cut gives
epsilon. Since the noise touches this Y only partly, it has another
old neighbor; merge it with that different old primitive using its
proper cut. These are disjoint vertex pairs and give a second epsilon.
The remaining network is Hilbert-contractive. Such entries are
O(epsilon^2 polylog(n)).

If no old Y splits, every right old primitive belongs wholly to either
the noise group or the coherent group. The noise group has odd size:
at least three branches give the bounded-absolute-row source covariance
followed by B, while a single branch is precisely the small h_3(G)/Y
exception. Its factor S has max-entry O(epsilon). The centered even
coherent coefficient leaves a nonempty even group of old branches.
Its covariance factor Gamma has Frobenius norm O(sqrt(n) polylog(n))
by the global root-map estimates. Thus

```math
\|S\circ\Gamma\|_F=O(\operatorname{polylog}(n)),
\qquad
|\langle M,S\circ\Gamma\rangle|/n=o(1).
```

Cross-group distinctness exclusions pull back to a local collision
between the left noise and coherent coefficient; those are among the
errors already charged before factorization. One must not restore an
arbitrary coherent pre-Walsh tensor here. This proves

```math
\mathbb E[(Atilde Z)^TBD_aZ]/n=o(1).
```

Combining Sections 4 and 5 leaves exactly (1), using the previously
audited normalized-nuclear comparison of Cov(Z) with T against the
bounded matrix BD_a.

## 6. Bounded closure and scope checks

The fixed-polynomial proof uses constants depending on the polynomial
degrees; approximation precedes the n limit and degree limits follow it.
Canonical W has uniformly bounded exponential moments near zero, by
the special decoupled bilinear representation of sums of S_i h_2(G_i).
This is stronger than generic degree-three hypercontractivity and makes
its subsequential laws moment-determinate.

A safe polynomial approximation argument uses finite catalogs, not an
arbitrary-frequency Taylor series. The family of coherent row laws is
weakly precompact and uniformly integrable at every fixed polynomial
degree. Under each limiting law polynomials are L2-dense: an orthogonal
signed measure has an analytic Fourier transform near zero with every
derivative zero, hence vanishes. The error for any fixed polynomial is
continuous across the family by uniform integrability. A finite cover
therefore supplies finitely many polynomials at each target accuracy;
deterministic row choices among this finite catalog preserve the
structural cut bounds. Gaussian variance parameters lie in a compact
interval and can be included in this finite-cover argument, including
variance zero through unnormalized Wick polynomials.

For fixed smooth bounded responses, local noise separation first controls
the actual and comparison L2 approximation errors; bounded-op energy
continuity then removes them. Extending r from fixed finite Hermite
polynomials is a separate ordered L2 approximation in the old Gaussian
frame. Arbitrary bounded measurable functions without the requisite
actual-law approximation are not silently covered.

This audit does not Gaussianize QS or QD, omit e(c0), evaluate the raw
bare-star cross, assert cutoff-free unbounded hard-threshold derivatives,
or prove any universal improvement to the established lower bound.

The final theorem's precise smooth hypotheses pass: f is bounded, odd
under simultaneous sign reversal, and Gaussian-a.e.-continuous; H is
bounded, even, and Gaussian-a.e.-continuous; psi is odd C2 with psi and
its first two derivatives bounded. The polynomial proof normal-orders at
the actual squarefree-main variance before replacing it by T_ii. In the
bounded closure the first-noise coefficient is controlled in its weighted
Gaussian L2 norm; fixed bounded deterministic coefficient approximants
allow the old nuclear covariance comparison to transfer that norm to the
raw channel. Uniform ideal-channel Hermite-tail coupling then controls
the final f approximation, without multiplying an uncontrolled raw row
variance by a merely averaged coefficient error.

The optional hard-threshold extension requires positive residual norm
tau^2=||r||_2^2>0. The smooth theorem has no such restriction. If tau=0,
the inverse-variance argument alone cannot justify applying a discontinuous
sign to an asymptotically vanishing raw residual. This qualifier was sent
to the theorem author for the optional threshold paragraph.

## 7. Independent finite stress test of the individual remainders

`computations/continued_audit_marked_energy_terms_2026_09_06.py` uses
actual positive-twin Steiner hollow signings, f=sin(Y), H=1, and the sine
feedback response. The conditional Gaussian sine formulas are exact, so
there is no quadrature error. It separately measures every remainder in
Sections 4 and 5, using an independent pilot for the deterministic mean
coefficient. The exact samplewise algebra summing those remainders to
the full projection difference was checked to error below 2e-14.

Reproducible command:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_audit_marked_energy_terms_2026_09_06.py --samples 6000
```

At n=56,240,992 the final centered-coefficient/deterministic-noise cross
was respectively -0.455805, -0.129170, -0.021268. The higher-noise
self-energy was 0.035293, 0.008244, 0.000741. The complete projection
error was -0.086339, -0.007488, -0.001151; its conditional pilot standard
errors were 0.003482, 0.001959, 0.000965. The asymptotic noise trace itself
remained nonzero, 0.026060, 0.027992, 0.028615.

Thus this test exposes substantial finite-size remainders while finding
no persistent falsifier. The standard errors omit pilot uncertainty;
these observations do not prove a convergence rate or a finite-n bound.

## 8. A stronger local covariance consequence for the next return

At a fixed polynomial stage define the exact returned remainder
eta=B(C-c0-D_a Z), using the already specified squarefree-main
approximations. Sections 4 and 5 bound its source crosses against
BD_dB for EVERY bounded deterministic diagonal d, not only d=a.
Their proofs also apply separately to each fixed original-degree main
of the left remainder and each original-degree old residual source on
the right. Choosing d_i to be the sign of the corresponding same-root
covariance therefore gives

```math
\frac1n\sum_i|\operatorname{Cov}(\eta_{p,i},Z_{q,i})|
\longrightarrow0
```

for every fixed pair of original degrees p,q. Different degrees are
already orthogonal. Combining this with the joint local Gaussian
comparison shows that the finite eta-component vector is Gaussian
independent of the OLD noise Z, after averaging over roots. This is
stronger than merely retaining a possibly nonzero eta/Z covariance.

It does not prove independence from the NEW literal return
L=B(c0+D_a Z). The condition in
`continued_feedback_second_return_boundary_2026_09_06.md` correctly
isolates the unproved full contraction of a lower-degree eta component
into a higher-degree L component. Proper cuts alone do not control that
contraction. Regressing equal-degree covariance and keeping appropriate
variance cutoffs are legitimate conditional steps; the needed higher-
degree full-contraction condition is not established by the present audit.
