# Coding and Boolean polynomial literature: exact usable scope

2026-09-06. Primary-source reconnaissance. No imported existence-of-limit
theorem has been found. This note separates actual equalities from
asymptotic-order estimates and records a concrete recurrence mapping.

## 1. Recent Sidon asymptotics are not a normalized limit

[Defant--Galicer--Mansilla--Mastylo--Muro, arXiv:2607.05594](https://arxiv.org/pdf/2607.05594)
Theorem 3.3 (printed p.35) estimates the Sidon, unconditional-basis and
Gordon--Lewis constants of fixed-degree Hamming-scheme spaces. Crucially,
the authors explicitly define their decorated asymptotic symbol to mean
upper and lower inequalities with factors c_1(q)^d and c_2(q)^d. It does
NOT mean ratio tending to one. Their precise Gaussian/Hermite limits
concern projection kernels and projection constants, a different invariant.
For q=2,d=2, the Sidon conclusion gives order sqrt(n), not its limiting
coefficient or existence of such a coefficient. These exact passages were
read directly in the primary PDF.

The predecessor
[arXiv:2302.00233v2](https://arxiv.org/pdf/2302.00233),
Theorem 5.1 (printed p.24), bounds the degree-d Boolean Sidon constant by
a degree-dependent constant times the degree-(d-1) projection constant.
Its introductory Theorem 4.1 discussion gives an exact Hermite-integral
limit for the PROJECTION constant. Neither statement supplies an
asymptotically equality-preserving transfer for all-moduli-one coefficients.

The exact relation to our problem is elementary. If Sid_n is the Sidon
constant of the weight-two Walsh system and N=binom(n,2), then

```
Sid_n=sup_(a!=0) [sum_(i<j)|a_ij|]/Q(a),
N/Sid_n=min_{sum|a_ij|=N} Q(a) <= M_n.                (1)
```

Thus even a limit for Sid_n/sqrt(n) would give only one side unless the
coefficient-modulus relaxation in (1) were proved asymptotically exact.
At n=4 it is already not exact: a four-cycle CHSH signing has Q=2 and
l1 mass4; rescaling it to mass6 gives Q=3, whereas M_4=4. This finite
example is not an asymptotic obstruction, but forbids equating the
objects by definition. No flatness theorem is silently assumed.

## 2. Gale--Berlekamp remains a different optimization

[Pellegrino--Raposo, arXiv:2111.00445](https://arxiv.org/pdf/2111.00445),
Section 6 (printed p.15), records sqrt(2/pi)+o(1)<=G_n/n^(3/2)<=1+o(1)
for the square Gale--Berlekamp game. Here G_n minimizes the bilinear norm
max_(x,y)|x^T A y| over full rectangular sign matrices, with two independent
spin vectors. This is not the hollow symmetric same-spin absolute cap.
The displayed asymptotic bounds do not state existence of a normalized
limit, either for G_n or for M_n. Hadamard-based constructive bounds
therefore do not furnish the missing original seed-transfer theorem.

## 3. A spherical secondary construction is exactly a cavity problem

The weight-two Reed--Muller/cut-code mapping is already audited in
`critical_scale_code_audit_2026_09_05.md`. The new equatorial-switching
lemma in `decisive_audit_equatorial_switching_2026_09_06.md` shows that
restricting to balanced weight-two truth tables changes M_n by only
O(sqrt(n log n)), at EVERY order.

[Gini--Meaux, WAPB secondary constructions](https://orbilu.uni.lu/bitstream/10993/52556/1/constrWAPB.pdf)
uses the standard Boolean selector construction
f(x,z)=(1+z)f_0(x)+z f_1(x), with arithmetic modulo two. Specialize
DIRECTLY to weight two: f_0 supplies the old pair table and f_1 on
weight one supplies the appended signed row v. Therefore its exact
restricted Walsh maximum is

```
max_s [|H_A(s)|+|v^T s|].                             (2)
```

Indeed the new affine query bit independently chooses the sign of the
linear term. This independent reconstruction shows that the standard
one-variable spherical recursion retains the full cavity/slack problem.
The triangle bound gives only Q(A)+n. Balancedness of f_1 restricts the
row's total sum, not the correlations v^T s on near-extreme old states.
No new upper recurrence follows merely from the secondary construction's
existence. Its paper's abstract and definitions were read; its full
quantitative nonlinearities theorem has NOT been imported here.

## 4. Access and preservation

The cited arXiv PDFs and the WAPB author's introductory PDF pages were
read through the web tool. Follow-up full-PDF retrieval from ORBilu
timed out; a local ordinary TLS fetch reported an expired certificate.
Certificate verification was not disabled. No downloaded upstream file,
credential, ignored input, or unrecorded external code is a dependency.
The elementary formulas (1)--(2) were reconstructed directly, not inferred
from unavailable pages or from secondary summaries.
