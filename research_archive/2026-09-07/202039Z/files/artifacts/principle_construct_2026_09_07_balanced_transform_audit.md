# Balanced-transform and independent-row repair audit

2026-09-07. **PASS** for the balanced-transform theorem in
`principle_invent_2026_09_07_balanced_transform_compiler.md`, read in full.

For a column of sum s, majority sampling has N=(q+|s|)/2 positions and
r=|s|/2 selected flips. The exact mean is q/(q+|s|) times the centered
column. The centered fluctuation covariance is four times the uniform
fixed-size-subset covariance, with nonzero eigenvalue
`4p(1-p)N/(N-1)`, p=r/N. This is bounded by 8S/q; its trace is at most
2|s|. The centered column norm is at most 2sqrt(2S). Independent-column
rectangular Bernstein has variance parameter at most max(8mS/q,2S),
giving precisely the stated operator error. The zero-excess case is
deterministic and requires no division by an empty majority set.

For selected rows of a dephased full Hadamard with its constant column
omitted, HH^T=MI-J. Thus (PH)(PH)^T=MP. The repaired TT^T estimate follows
from the operator perturbation and the O(sqrt(M)) factor norms.

The signed-swap factorization transfers the repair error uniformly over
the edge mask. Its o(M) operator error costs o(n^(3/2)) in quadratic
energy when n and M^2 are comparable. The original projected Boolean cap
can be invoked only on fibre-balanced words; the artifact correctly does
not apply it to arbitrary centered nonbalanced words. Hollowing a diagonal
rank-one fibre block produces only an O(1) operator diagonal correction.

The conditional deletion inequality follows exactly from convexity when
the deleted signs remain independent fair signs given the retained data.
Conditioning on phase-independent selectors or excesses preserves this.
Conditioning on the whole cap would not automatically preserve it and
must not be silently inserted. With this stated scope, the balanced-face
certificate transfer is valid; the mixed-mean response remains open.

## Related independent-row repair, reconstructed for synthesis

For the k=2 directed-sign compiler, start with independent fair rows eta
and independently repair each row's majority to exact balance (using an
extra diagonal port if needed). Conditional on the original eta, write
the repaired row mean as

    mu_ij=(1-p_i)eta_ij-a_i,
    p_i=|s_i|/(n+|s_i|),  a_i=sign(s_i)p_i.

The repaired rows remain independent, with centered covariance at most
C p_i I and centered norm at most C sqrt(r_i), r_i=|s_i|/2.

For D'_ij=S_ij eta'_ij eta'_ji, expose repaired rows in order. The exact
Doob increment at row i is `e_i w_i^T+w_i e_i^T`, where

    w_ij=S_ij(eta'_ij-mu_ij) eta'_ji,   j<i,
    w_ij=S_ij(eta'_ij-mu_ij) mu_ji,     j>i,
    w_ii=0.

The partner multipliers are predictable and have magnitude at most one.
Hence conditional covariance of w_i is at most Cp_i I, its conditional
squared norm expectation is at most Cnp_i=O(r_i), and its norm is at
most Csqrt(r_i). The sum of predictable star-square matrices is bounded
by `C(max_i r_i+sum_i p_i) I`. Thus the proposed matrix-Freedman argument
has exactly the advertised variance and increment parameters; it omits
no hidden quadratic dependence between repaired rows.

This operation produces independent uniform BALANCED rows, by coordinate
permutation symmetry, not orthogonal rows. Any consequence for genuinely
orthogonal incidence arrays would need a separate argument.

For clarity, the intended physical frame has only two rows, 1 and eta_i.
Exact balance DOES make those two rows orthogonal, so the repair covers
that intended independent two-row-frame law. Only mutual orthogonality
of the different macro incidence rows eta_i is excluded by the caution.
