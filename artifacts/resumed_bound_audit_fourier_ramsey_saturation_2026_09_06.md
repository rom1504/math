# Independent audit: Fourier Ramsey saturation and its actual-signing scope

Date: 2026-09-06. The exact finite-field theorem, robust finite-field
version, spectrally flat additive-Cayley consequence, and exact finite
abelian-group extension all pass independent reconstruction. No general
signing theorem or robust arbitrary-group extension is asserted.

Main audited artifact:
`resumed_convergence_fourier_involution_ramsey_saturation_2026_09_06.md`.
The analytic dependency
`resumed_convergence_paley_squarefree_audit_2026_09_06.md` was read fully
and its squarefree and finite-moment steps independently reconstructed.

## 1. Primary Ramsey inputs and a coefficient convention

The primary research article
[Frankl--Graham--Rodl, Iterated Combinatorial Density Theorems (1990)](https://fanchung.ucsd.edu/ron/papers/90_04_iterated.pdf)
states both needed classical results: its page 99 gives the finite
Ramsey theorem for coloring fixed-dimensional vector subspaces; page
105 gives Deuber's finite theorem for (m,p,c)-sets. The latter permits
any positive target c, including c=1. Its definition requires the
whole generated set, not just its generators, to lie in the positive
integers. These hypotheses were checked directly, not inferred from
an informal description of the theorem.

To remove any strict-versus-weak coefficient-bound convention, request
the Deuber target parameter K+1. It then includes every coefficient
lambda with |lambda|<=K. For generators z_1,...,z_r, positivity at
lambda_j=-K gives

    z_i>K sum_(j<i) z_j.                               (1)

Fix one positive host (M,P,C)-set supplied by the finite theorem and
let N be its maximum. Every coloring of [1,N] restricts to this host,
so it contains the desired target set. Each z_i itself is in that set
because c=1, hence 1<=z_i<=N. This is the finite uniform N used below.

## 2. Analytic finite-pattern lemma

For every error epsilon>0 one can fix an odd polynomial degree D,
a number r of coordinates, and a finite moment cutoff K, in that
order, with constants independent of the odd characteristic.

For iid Y_j uniform on F_p, X_j=sqrt(2)cos(2pi Y_j/p) have mean zero,
variance one and absolute value at most sqrt(2), uniformly in p.
The same assertions hold on every odd cyclic group Z_m with m>=3.
Let S=r^(-1/2)sum X_j and define the ordered squarefree statistic W_l.
The exact identities

    E W_l^2=E[He_l(S)W_l]=l!(r)_l/r^l

follow because every distinct index in W_l must be covered in the
other factor; all lower-degree terms vanish. Gaussian moment
convergence for bounded independent coordinates is uniform here.
Thus a finite odd Hermite approximation to sign(S) may be replaced
in L2 by a sum of odd-degree W_l, with arbitrarily small error after
r is fixed sufficiently large.

Substituting X_j(t)=sqrt(2)Re chi_(a_j)(t) leaves a polynomial P whose
frequencies are nonempty signed squarefree sums of odd support size
at most D. It suffices to arrange that these frequencies are nonzero
and all carry one multiplier sign. No whole trace-label independence
is necessary. It is enough that all integer relations of total
coefficient size at most K agree with those of the iid model.

The discontinuity of sign creates no hidden transfer step: if h is
odd, (sign(u)-h(u))^2 is continuous at zero, even with sign(0)=1.
Approximate that continuous function uniformly on the fixed interval
[-sqrt(2r),sqrt(2r)] by a polynomial and enlarge K. The square of the
Hermite-to-squarefree difference is already a polynomial of degree
at most 2D. These two finite-moment transfers give ||f-P||_2<epsilon
for a genuinely Boolean f. All choices precede the group-size limit.

## 3. Exact saturation: uniform characteristic and group thresholds

In large characteristic p>K r N, use the positive integer Deuber
generators from Section 1 as prime-field frequencies. A nonempty
signed squarefree sum can be negated so its last coefficient is +1;
it is then in the monochromatic target set. Evenness of the
multiplier handles the possible global negation.

For any nonzero short integer relation, its largest-index term and
(1) show its integer value is nonzero, while its absolute value is
at most K r N<p. Hence it also remains nonzero in the field. This
works even in prime fields, where full coordinate independence
would be impossible.

For each of the finitely many remaining odd primes, color each
projective line by the multiplier of its unique vector whose first
nonzero ambient coordinate is one. Vector-space Ramsey supplies a
homogeneous r-subspace once the extension degree is sufficiently
large. An echelon basis has increasing pivots equal to one. The
first nonzero coordinate of any signed squarefree sum is therefore
exactly +1 or -1, so line normalization changes it only by a sign.
This is why even multipliers suffice; invariance under all nonzero
scalars is not being assumed. Linear independence gives every
required short-relation condition automatically.

Take the maximum of the finitely many small-prime dimension
thresholds and the large-characteristic threshold. This proves one
threshold in q for every fixed epsilon, uniformly over all even
off-zero sign multipliers. The resulting P is mean zero and belongs
to one exact eigenspace. Orthogonal projection gives

    max_(Boolean f) |<f,Uf>|>=1-2epsilon^2,
    |E f|<=epsilon.

The exact arbitrary finite abelian-group extension is also sound.
An even group has a nontrivial order-two Boolean character. For odd
groups, an unbounded exponent supplies a cyclic dual subgroup whose
order replaces p in the Deuber short-relation proof. Bounded exponent
and unbounded order force an elementary abelian subgroup of growing
rank for one of finitely many primes. More explicitly, bounding the
exponent and all those ranks bounds the group order by the cyclic
factor decomposition. Thus this case split also yields a uniform
finite threshold, rather than just a subsequence argument.

## 4. Robust finite-field step and actual Cayley signings

Fix epsilon, eta and all analytic/Ramsey parameters FIRST. In large
characteristic the host is [1,N]; in each small characteristic it is
the nonzero part of a fixed finite-dimensional subspace. Their sizes
are bounded by one finite constant depending only on the accuracy.

If a bad frequency set has density tending to zero, multiplication
of the entire host by a uniformly random nonzero field element has
expected bad intersections at most host_size times bad_density
(with the harmless q/(q-1) adjustment). Hence some scalar dilation
avoids the entire host. All short relations and linear independence
are preserved. Coloring this dilated host by multiplier signs now
gives P with every supporting multiplier satisfying s*m>=1-eta.

If ||U||op<=C and ||f-P||2<epsilon, then

    s<f,Uf>>=(1-eta)(1-epsilon)^2-C epsilon(2+epsilon).

First send q to infinity, then eta and epsilon to zero. This proves
the robust lower Rayleigh limit. It does not require any independence
of the bad set, which is an arbitrary deterministic spectral set.

For a symmetric hollow additive-Cayley signing A on F_q, put
U=A/sqrt(q). Its Fourier multipliers are real and even and

    q^(-1)Tr(U^2)=1-1/q.

If ||U||op<=1+o(1), this second-moment identity forces the density of
frequencies with ||m|-1|>eta to vanish for every fixed eta>0. The
robust theorem gives max_Boolean |<f,Uf>|>=1-o(1), while the operator
norm gives the matching upper bound. In the original half-energy
normalization this is exactly

    Q(A)/q^(3/2)->1/2.

This rules out spectrally near-Frobenius-optimal additive-Cayley
families over finite fields as sub-1/2 upper constructions. It does
not rule out Cayley families with fixed spectral excess, arbitrary
non-Cayley signings, or original nonconvergence. The robust argument
uses scalar transitivity on field frequencies and is not silently
extended to arbitrary finite abelian groups.
