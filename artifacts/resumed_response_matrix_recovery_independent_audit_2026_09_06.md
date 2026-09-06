# Independent reconstruction of exact finite-frame recovery

Date: 2026-09-06. Status: the director's constructive converse is
independently reconstructed below. The matrix theorem and the anchored
relative-height extension are new in this campaign; the unanchored
scalar critical existence criterion was already proved in
`fresh_tree_fixed_point_critical_independent_2026_09_05.md` and is not
claimed as new here.

## 1. Exact necessary and sufficient condition

Let U be the canonical creation isometry from jointly even Gaussian L2
to first chaos. Let g: R^d -> R^d be jointly even and measurable, with

    E[g(N)g(N)^T]=I_d,   N standard Gaussian.

For 0<=P<=I put K(P)=E[g(X)g(Y)^T], where (X,Y) have standard marginal
laws and cross covariance P. This is the UNCENTERED Gram kernel.
Then the following conditions are equivalent:

1. There exists a standard Gaussian first-chaos frame G in the actual
   canonical space such that G=U g(G), componentwise.
2. The matrix iteration P_0=0, P_(h+1)=K(P_h) tends to I.
3. There is no proper PSD subfixed matrix: no 0<=P<=I, P!=I with K(P)<=P.

When the frame exists, it is unique as a random vector on this fixed
canonical space, not merely unique in distribution.

The necessary direction and Loewner monotonicity are proved in
`resumed_response_matrix_height_feedback_obstruction_2026_09_06.md`.
Here is a direct construction proving the converse, without importing
an RDE endogeny theorem.

## 2. Construction and the compatibility identity

Let F_h be the canonical Gaussian sigma-field through tree height h,
Pi_h its first-chaos projection, and F_0 trivial. For an even f,

    Pi_(h+1) U f = U E[f|F_h].

Set Y_0=0. Recursively define the even F_h-measurable vector

    t_h = E_Z g(Y_h + (I-P_h)^(1/2) Z),
    Y_(h+1) = U t_h,

where Z is an external standard Gaussian, integrated out rather than
adjoined to the canonical probability space. We claim inductively:

    Y_h belongs to first chaos through height h,
    E[Y_h Y_h^T]=P_h,
    Pi_h Y_(h+1)=Y_h.                              (1)

The parity issue is material but harmless: replacing the entire
underlying Gaussian family by its negative negates Y_h; evenness of g
and symmetry of external Z imply that t_h is unchanged. Thus U t_h
is always defined.

The covariance formula follows directly from the shared-Y_h Gaussian
coupling:

    E[t_h t_h^T]=K(P_h).

For compatibility, condition t_h on F_(h-1). By the preceding inductive
identity, Y_h=Y_(h-1)+Delta_h, where Delta_h is a first-chaos Gaussian
independent of F_(h-1), with covariance P_h-P_(h-1). Integrating both
Delta_h and external Z leaves covariance I-P_(h-1), hence

    E[t_h|F_(h-1)]=t_(h-1).

The creation/projection identity now gives
Pi_h Y_(h+1)=U t_(h-1)=Y_h. The base case is
Y_1=(Eg)G_edge and Pi_0 Y_1=0. This proves (1), including singular
P_h; no invertibility is used.

If P_h tends to I, the compatible first-chaos projections form an
L2-bounded Gaussian martingale and converge in L2 to a standard frame G.
Its projection on F_h is Y_h and its independent residual covariance
is I-P_h. Therefore

    E[g(G)|F_h]=t_h.

Martingale convergence and the isometry give

    U g(G)=lim_h U t_h=lim_h Y_(h+1)=G.

For uniqueness, any candidate frame must have projection Y_0=0; its
conditional Gaussian law forces each succeeding projection to be the
same recursively defined Y_h. Height exhaustion then forces the frame
itself to be G.

## 3. Why no proper subfixed matrix is exactly sufficient

K is continuous on the compact PSD interval [0,I], even for merely L2
g. One proof approximates each g_i by a finite Gaussian polynomial.
Gaussian conditional expectation is L2-contracting, so the error in
each Gram entry is uniformly bounded in P by a constant times the L2
approximation error. Polynomial kernels are continuous in P.

The iteration P_h is Loewner increasing and bounded by I. Its limit
P_infinity is therefore a fixed point by continuity. If it is proper,
it is itself a forbidden subfixed matrix. Conversely, any proper
subfixed matrix bounds all iterates by monotonicity. This proves the
equivalence of conditions 2 and 3, and closes the exact recovery theorem.

The construction still makes sense when P_infinity<I; it then gives a
deficient-covariance Gaussian vector satisfying a smoothed equation.
It does NOT turn that vector into a standard endogenous solution by
adding an external Gaussian. Such an addition would change the required
same-space fixed-point equation.

## 4. Anchored recovery, including actual cyclic anchor frames

Let A be an already actual finite standard creation-closed frame:
A=U a(A), with jointly even inverse vector a and E[aa^T]=I. The frame
need not be finite canonical or causal; it can be a previously constructed
cyclic frame with countable canonical support.

Let h(A,z) be jointly even in all inputs, scalar for now, and suppose
under independent standard (A,Z) that

    E h(A,Z)^2=1,   E[a(A)h(A,Z)]=0.

Define K_A(q)=E[h(A,X)h(A,Y)], with A shared and Corr(X,Y)=q. There exists
a standard first-chaos Z independent of A satisfying Z=U h(A,Z) if and
only if K_A^h(0) tends to 1, equivalently K_A(q)>q for all q<1. The
extension is unique with A fixed.

One proof augments the full inverse vector by h and applies Sections
1--3: the limiting anchor covariance block must be I because the
anchor subsystem is already actual, and then PSD forces its cross
block to vanish. Any deficient innovation limit would be a scalar
fixed point of K_A below one.

There is also an exact anchored-depth construction giving useful error
bounds. Let V_0=span(A) inside canonical first chaos, and recursively

    V_(h+1)=U L2_even(sigma(V_h)).

These are closed first-chaos subspaces and are nested: creation closure
puts V_0 inside V_1, and induction gives all later inclusions. Their
union is dense in canonical first chaos, since U1 enters V_1 and each
subsequent step contains the next canonical tree height. Moreover

    Proj_(V_(h+1)) U f = U E[f|sigma(V_h)].

Start Z_0=0, q_0=0 and define

    Z_(h+1)=U E_N h(A,Z_h+sqrt(1-q_h)N).

Orthogonality to a(A) gives Proj_(V_0)Z_1=0. The compatibility induction
of Section 2, now with A already measurable in every level, proves

    Proj_(V_h) Z_(h+1)=Z_h,
    E Z_h^2=q_h,   q_(h+1)=K_A(q_h),
    E[A Z_h]=0.

Thus this is an actual canonical-space construction, not a distributional
recursion on hypothetical independent child copies. In the successful
case its exact truncation error is E|Z-Z_h|^2=1-q_h.

The same argument works for a finite vector innovation, with scalar q
replaced by a PSD matrix Q and K_A(Q); no proper PSD subfixed innovation
matrix is again necessary and sufficient.

## 5. Critical scalar classification and convergence rate

Expand only in the scalar innovation:

    h(A,z)=sum_(ell>=0) b_ell(A) h_ell(z),
    w_ell=E b_ell(A)^2,  K_A(q)=sum_ell w_ell q^ell,
    sum_ell w_ell=1,   D=sum_ell ell w_ell.

Here odd ell may occur, unlike the purely unanchored even scalar case.
Consequently the exceptional deterministic-one-child kernel K_A(q)=q
must be stated explicitly. Exact endogenous scalar extension exists iff

    D<1, or D=1 and some w_ell>0 for ell>=2.       (2)

If D<1, convexity gives K_A(q)>=1+D(q-1)>q. If D=1 with a genuine
nonlinear term, strict convexity above the tangent at 1 gives the same
strict inequality. If D=1 without such a term, normalization forces
w_1=1 and K_A(q)=q, which never leaves zero. If D>1, including infinite
D, the left derivative at one forces K_A(q)<q somewhere below one.

For D<1, 1-q_h<=D^h. At D=1, assume the finite second factorial moment

    0<K_A''(1)=sum_ell ell(ell-1)w_ell<infinity.

Then with delta_h=1-q_h,

    delta_(h+1)=delta_h-(K_A''(1)/2)delta_h^2+o(delta_h^2),
    1/delta_(h+1)-1/delta_h -> K_A''(1)/2,

and Cesaro averaging gives the sharp anchored-depth rate

    1-q_h ~ 2/[K_A''(1) h].                       (3)

Each delta_h is positive: K_A(q)<1 for q<1 because some positive-degree
coefficient is present. Finite polynomials satisfy the second-moment
hypothesis, but polynomiality is not needed for (3).

Reconciliation with the archive: the 2026-09-05 note already proved
unanchored scalar existence and uniqueness for D<=1, using a forced
tree-coefficient/Galton--Watson recursion, and the subcritical D^h tail.
The unanchored even case cannot have the exceptional w_1=1 kernel.
What is added here is exact multivariate recovery from PSD kernel
iteration, extension relative to arbitrary actual finite cyclic anchors,
and the explicit anchored critical rate. These are structural feasibility
results, not a new numerical lower certificate or a convergence theorem
for M_n/n^(3/2).
