# Original convergence: lower sampling and shared optimizer laws

Date: 2026-09-06. Status: proof-obligation audit and two exact shared-law
obstructions. No convergence or nonconvergence theorem is claimed.

Write H_A(x)=sum_(i<j)a_ij x_i x_j, Q(A)=max|H_A|,
M_n=min_A Q(A), and c_n=M_n/n^(3/2). The current certified interval is

    .4314603928237005 <= liminf c_n <= limsup c_n <= .5.

## 1. Exact scope of the proposed lower-sampling envelope

The proposed statement

    c_n >= c_k-e(k),  n>=k,   e(k)->0

would prove convergence. However existence of some unspecified e(k)
is exactly equivalent to convergence of the bounded sequence c_n,
not an independently weaker scalar obligation. The smallest error is

    e_*(k)=max(0,c_k-inf_(n>=k)c_n).

If c_n converges, e_*(k)->0. If limsup c_n exceeds liminf c_n, choosing
k near the limsup and a later n near the liminf keeps e_*(k) bounded
away from zero. A proved explicit error such as C/sqrt(k) would be
substantially stronger and would require actual signing structure.

The pointwise scale-preserving principal restriction statement in
`nested_restriction_paving.md` is already false without error at small
orders, including exact conference examples. The minimizing-fibre
restriction examples in `cross_order_outward_director_review.md` are
also only finite warnings. None produces a growing-order gap between
the ACTUAL minima, so none falsifies the displayed asymptotic envelope.

## 2. Improved range bounds do not eliminate the low-profile obstacles

The paired full-response construction bounds half the total range, not
only Q. Let W=(max H-min H)/2 and let c_*=.4314603928237005. At any
asymptotic near-minimizer,

    W>=(c_*-o(1))n^(3/2),  Q<=(.5+o(1))n^(3/2).

Thus both oriented extrema are at least
(2c_*-.5-o(1))n^(3/2). Exact range superadditivity under a vertex
partition gives every sublinear principal block the range budget

    range(A[S]) <= (1-2c_*+o(1))n^(3/2).

In particular a coherent clique of size gamma n^(3/4) is limited to
gamma<=sqrt(2-4c_*), approximately .5236, not to gamma=o(1). Every
fixed lower constant below .5 leaves a nonzero allowance of this type.

More decisively, the exact-sign near-minimizer construction in
`cross_order_exposed_shell_projection_barriers.md` survives every
improvement of the lower constant. Start with an exact order-(n-s)
minimizer, s=floor(sqrt(n)/log n), and attach s vertices with all-positive
cross edges. Its cap is at most M_n+sn+O(s^2)=M_n+o(n^(3/2)), while its
operator norm is at least sqrt(s(n-s))~n^(3/4)/sqrt(log n).

This example does not claim exact optimality. It does prove that
asymptotic near-minimality, even with the new universal lower bound,
does not force a bounded normalized operator norm or the needed raw
operator concentration. Fixed local-profile indistinguishability also
does not compare actual minima: a planted high-cap comparison matrix
is not thereby another minimizer. These distinctions prevent promoting
the archive's proof-class barriers into a scalar-envelope falsifier.

## 3. What exact edge-flip optimality really says

Use oriented features v_e=σx_i x_j, σ in {+1,-1}, and write

    g(v)=M_n-sum_e a_e v_e>=0

for an exact minimizing signing A. Flipping a set F of coefficients
changes its score at v by -2 sum_(e in F)a_e v_e. Since the resulting
signing cannot have cap below M_n, one has the exact statement

    for every F there exists v such that
       g(v)+2 sum_(e in F)a_e v_e <= 0.             (1)

For a single edge, its witness necessarily has g<=2 and a_e v_e=-1.
For |F|=r, a witness has g<=2r. Different flip sets may require
different states. Thus (1) is not a common-law stationarity condition.
Convexifying the flip choices and exchanging min and max would require
an additional rounding estimate, which is precisely where a critical
scale error can enter.

There is an unavoidable budget for any attempted common law. For EVERY
law mu on oriented states with E_mu g<=eta,

    sum_e (a_e E_mu v_e)_+ >= M_n-eta.             (2)

This follows by summing a_e E_mu v_e=M_n-E_mu g and dropping its
negative coordinates. More generally, if nonnegative weights lambda_T
give a fractional edge cover, sum_(T containing e)lambda_T>=1, then

    sum_T lambda_T sum_(e in T)(a_e E_mu v_e)_+
       >= M_n-eta.                                (3)

Hence one near-active law cannot have a vanishing aggregate positive-
alignment budget on an edge-covering family. At eta=o(n^(3/2)), the
budget is already at least (c_*-o(1))n^(3/2). This does not rule out
more structured joint laws, but it excludes the naive simultaneous
antialignment upgrade of the separate edge witnesses.

## 4. Exact minimizers need not possess a radial ground-state law

An apparently milder shared-law proposal would ask for a law on exact
oriented ground states such that E v_e=rho a_e for all edges. Activity
would force rho=M_n/binom(n,2). This is false even for an actual exact
minimizer, not merely for a low-cap comparison matrix.

Take n=4, with all a_ij=+1 except a_23=-1 (indices 0 through 3).
Direct enumeration gives Q(A)=4. Every four-vertex signing has
E H_A^2=6 and only even energy values, so Q(A)>=4; this A is therefore
an exact minimizer without relying on an archived classification.

For every spin x,

    x_0x_1+x_0x_3+x_1x_2-x_2x_3
      =x_1(x_0+x_2)+x_3(x_0-x_2) in {-2,2}.

The same absolute bound holds for oriented features. If E v=rho a,
taking expectation gives 4rho<=2. Thus rho<=1/2, whereas activity
would require rho=4/6. In fact rho=1/2 is attained by the equally
weighted four oriented states

    (σ,x)=(+,(1,1,-1,1)), (+,(1,1,1,-1)),
          (+,(1,1,1,1)), (-,(1,1,-1,-1)).

Their mean oriented energy is 3, forcing mean gap 1. Consequently
even exact minimization does not supply a radial subgradient supported
on exact ground states. This is a finite falsifier of the exact claim;
it does not refute an asymptotic radial approximation using a growing
critical gap window.

The exact statements above can be replayed with
`computations/resumed_bound_audit_radial_ground_law_exact.py`.

## 5. A scalable edit-distance consequence of a genuine radial law

The radial-law obstruction has a useful positive counterpart for a
specific structured class. Let H be a symmetric full sign Hadamard,
H^2=nI, possessing an orthogonal Boolean eigenbasis x^(1),...,x^(n).
Write its eigenvalues as sigma_j sqrt(n), sigma_j in {+1,-1}. Uniformly
sample j and use the oriented spin (sigma_j,x^(j)). Spectral expansion
gives the exact oriented covariance

    E[sigma_j x^(j)(x^(j))^T]=H/sqrt(n).

If A is any hollow signing differing from H on d UNORDERED off-diagonal
edges, this law immediately gives

    Q(A)>=E[sigma_j H_A(x^(j))]
         =(binom(n,2)-2d)/sqrt(n),

or

    Q(A)/n^(3/2)>=1/2-1/(2n)-2d/n^2.              (4)

For example H=(J_4-2I_4)^(tensor a) has the required Boolean eigenbasis,
obtained by tensoring the four Walsh eigenvectors of J_4-2I_4. Switching
and relabeling preserve the property. Thus a hypothetical sequence of
actual minimizers with c_n<=1/2-delta must remain at Hamming distance at
least (delta/2-o(1))n^2 from every member of this class. An o(n^2)-edit
recovery into these resonant Hadamards cannot preserve a genuinely
sub-1/2 minimizing subsequence.

The archive already established that a complete Boolean eigenbasis
forces the 1/2 scale; (4) records its exact robust edit-distance
consequence. It is not a theorem about arbitrary conference matrices,
which need not have such an eigenbasis. The displayed radial law also
need not live inside an O(sqrt(n)) exact-cap window: for the hollowed
tensor examples the two eigenvalue orientations differ in cap energy
by O(n). Consequently (4) does not supply the critical-window insertion
law ruled open in Sections 3--4.

## 6. Remaining actual obligation

The old fixed-block insertion criterion remains a sufficient variant
of the archived cap-discrepancy route, not new progress by itself.
Likewise the lower-sampling envelope cannot be obtained by renaming
convergence or by assuming that regularized phase norms coincide with
original minimizing caps. The present nonlocal phase modules give
LOWER tests on larger amplifications; they do not force an original
minimizer to be upper-stable under those amplifications.

A useful new optimizer theorem must control conditional geometry of
the jointly labelled near-active states, while respecting the budget
(3) and allowing the radial failure in Section 4. Neither separate
edge witnesses, unconditional permutation averaging, low cap, nor
asymptotic near-minimality supplies that missing information.
