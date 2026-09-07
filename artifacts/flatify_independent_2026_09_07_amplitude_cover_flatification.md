# Diagonal amplitude-cover flatification

2026-09-07. A sufficient actual-sign operation, not a convergence theorem.
The underlying bounded-coefficient Bernstein rounding lemma is already in
`second_phase_independent_abstraction.md`; the point here is an explicit
condition paying for coefficients outside the sign cube without clipping.

Let B be a real symmetric hollow n by n matrix and write
Q(B)=max_x |sum_{i<j} B_ij x_i x_j|. Suppose nonnegative numbers ell_i obey

    |B_ij|^2 <= (1+ell_i)(1+ell_j),  i != j.

Put L=sum_i ell_i, s_i=(1+ell_i)^(-1/2), C=SBS and
V=sum_{i<j}(1-C_ij^2). Then an actual hollow full signing A exists with

    Q(A) <= Q(B) + sqrt(2 V h) + 4h/3,
    h=(n+2) log 2.

Indeed C is in the coefficient cube. Multilinearity implies Q(SBS)<=Q(B):
every Sx is in [-1,1]^n, and a multilinear quadratic has its extrema on
vertices. Independently round each C_ij to a sign of mean C_ij. Each fixed
spin's centered energy has variance V and summands of absolute value at
most two. Bernstein at sqrt(2Vh)+4h/3, followed by the two-sided union bound
over all 2^n spins, has failure probability at most 1/2.

The variance deficit has the useful deterministic estimate

    V <= binom(n,2) - sum_{i<j} B_ij^2 + (n-1)L + L^2/2.

To see this, set p_ij=(1+ell_i)(1+ell_j). Since B_ij^2<=p_ij,

    B_ij^2 - B_ij^2/p_ij <= p_ij-1
                              =ell_i+ell_j+ell_i ell_j.

Sum and use sum_{i<j}ell_i ell_j <= L^2/2. No uniform bound on the
uncontracted amplitudes B_ij is required. The right side is automatically
nonnegative whenever the cover hypothesis holds.

Consequently, if L=o(n) and sum B_ij^2 >= binom(n,2)-o(n^2), then V=o(n^2)
and Q(A)<=Q(B)+o(n^(3/2)). Thus a real-coefficient seed-transfer operation
meeting these two concrete hypotheses would transfer the ORIGINAL cap.
Neither small Frobenius distance to the sign cube alone nor a low cap alone
establishes the amplitude-cover condition. A matching of constant-size
amplitude excesses already forces L=Omega(n), although it may be harmless
by a separate O(n) direct correction. This is a sufficient, not necessary,
criterion.

The cover can be checked by a convex program in r_i=log(1+ell_i): minimize
sum_i(exp(r_i)-1), subject to r_i>=0 and
r_i+r_j>=2 log^+|B_ij|. This computational formulation does not provide the
missing seed-preserving B.

## Scope relative to prior surgery

For B=A-Delta with A a signing, the previously proved target-contraction
surgery obtains a cover from ell_i=|Delta|_ii, since
|Delta_ij|<=sqrt(ell_i ell_j) and
1+sqrt(ell_i ell_j)<=sqrt((1+ell_i)(1+ell_j)). That route also estimates
the variance directly using the nuclear budget. The present criterion
allows arbitrary real B and covers supplied by entry geometry, but does
not supersede or independently reprove a new nuclear-budget operation.

No actual all-order low-cap B meeting this condition has been constructed.
