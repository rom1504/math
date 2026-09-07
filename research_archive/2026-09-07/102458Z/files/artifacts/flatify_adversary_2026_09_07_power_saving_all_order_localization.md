# Power-saving all-order subspace localization: prime gap and zero embedding

2026-09-07. This audits the all-order extension of the actual-sign
orthogonal-target rounding operation. It needs no prime-in-progressions
short-interval theorem and no assumed Hadamard conjecture.

## Published short-interval input

Baker, Harman and Pintz, *The Difference Between Consecutive Primes, II*,
Proc. London Math. Soc. (3) 83 (2001), 532--562, Theorem 1 on the first
page, states that every sufficiently large x has a prime between
x-x^0.525 and x. The actual paper was read directly through the
[author-uploaded published PDF](https://www.researchgate.net/profile/Roger-Baker-2/publication/228991240_The_Difference_Between_Consecutive_Primes_II/links/53e14c630cf2235f35296d79/The-Difference-Between-Consecutive-Primes-II.pdf).
The [publisher record](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/plms/83.3.532)
identifies the 2001 volume and pages; its later online publication metadata
should not be mistaken for the theorem's original date.

Put t=n/2-1 and x=t+2t^0.525. For large t, x-x^0.525>=t, so that theorem
gives a prime p>=t with p-t=O(t^0.525). Only an unrestricted prime is used.

## Hadamard order 2(p+1) for EVERY odd prime

Here is a direct verification of the needed Paley construction, rather
than importing a congruence restriction incorrectly. Let Q_ij=chi(i-j)
on F_p, with chi(0)=0. Then Q1=0, QQ^T=pI-J, and
Q^T=chi(-1)Q. These follow from the elementary quadratic-character
correlation sum; for distinct shifts its value is -1, obtained by counting
solutions of v^2=t(t-1) via (2t-1-2v)(2t-1+2v)=1.

If p=3 modulo 4, [[1,1^T],[1,Q-I]] is a sign Hadamard of order p+1;
tensor with H2 to obtain order 2(p+1). If p=1 modulo 4, the symmetric
conference matrix S=[[0,1^T],[1,Q]] obeys S^2=pI. The sign block matrix

    [[S+I,S-I],[S-I,-S-I]]

has square 2(p+1)I. Thus it has the same required order 2(p+1). These
are the two constructions originating in [Paley's 1933 paper](https://onlinelibrary.wiley.com/doi/10.1002/sapm1933121311);
the algebra above supplies the relevant proof even though the publisher's
accessible page exposed only its bibliographic record.

Consequently every sufficiently large n has an available Hadamard order

    h>=n, h-n=O(n^(21/40)).

## Zero embedding avoids a square-root padding loss

Embed the original arbitrary subspaces E_L,E_R of R^n in R^h by ZERO
coordinates, not by padded Boolean centers. Their dimensions and original
orthogonal projections are unchanged. Apply the finite-rank orthogonal
localization and full-sign rounding in dimension h, then restrict its
actual sign output to the original n-by-n rectangle.

For physical x,y in {+/-1}^n, their zero extensions have norm sqrt(n)
and belong to [-1,1]^h. A uniform bilinear difference bound on the Boolean
cube extends to this box by separate convexity, so the full rounding
error applies with no extra term for deleted coordinates. The target
residual bound is exactly

    sqrt(h)*sqrt(n-||P_Lx||^2)*sqrt(n-||P_Ry||^2)
             +2sqrt(h*n*d),

where d=dim(E_L)+dim(E_R). The zero cross-projection construction in
`flatify_adversary_2026_09_07_rectangular_orthogonal_rounding_audit.md`
justifies the localized terms without requiring projections to preserve
the cube.

Changing the leading scale from sqrt(h) to sqrt(n) costs at most

    n*(sqrt(h)-sqrt(n))=O(sqrt(n)*(h-n))=O(n^(41/40)).

Thus the resulting full sign bridge at EVERY original order satisfies
the desired residual bound with error

    O(n^(5/4)*sqrt(d)+n*sqrt(d)+n+n^(41/40)).

For d>=1 this is simply O(n^(5/4)*sqrt(d)). In particular
d=O(n^(1/2-eta)) gives a fixed normalized power saving O(n^(-eta/2)).
This is stronger than a sign-padding/deletion estimate
O(n*sqrt(h-n))=O(n^(101/80)); the latter is unnecessary for this
arbitrary-subspace operation because zero embedding preserves the exact
projection geometry.

The theorem constructs an actual full sign bridge with an explicit
all-order rate. It does not furnish subspaces that control all actual
child energies, and therefore does not itself establish recurrence.
