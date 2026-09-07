# Paley nonclosure: a precise input, not a pressure theorem

Status: partial route, set aside after the exact Gaussian-phase theorem.

For quadratic-character Paley columns indexed by distinct field elements,
the correlation of a product of s>=2 distinct columns with any one column
is O_s(q^-1/2), apart from O_s(1/q) border/diagonal corrections. If the test
column repeats one factor, remove its square first; the remaining polynomial
has distinct roots and is not a square. The only exceptional linear case
is s=1 with the same test column. Thus a fixed nonlinear multilinear function
of r distinct Paley columns has o(1) energy in ANY fixed-size set of spectral
columns. Its linear terms are the only finitely concentrated part.

The exact character-sum hypothesis used here is the non-square-polynomial
bound with d distinct roots, magnitude at most (d-1)sqrt(q). A primary
research-paper statement and application is Appendix A, Theorem A.1 in
[Ample simplicial complexes](https://link.springer.com/article/10.1007/s40879-021-00521-5).
The original Weil paper was located at DOI10.1073/pnas.34.5.204, but its
PDF could not be opened in the available web reader; no uninspected detail
of that PDF is used.

To quantify the conclusion, write a nonlinear part f_nl as a fixed sum of
Walsh monomials of the chosen column signs. Each normalized Fourier
coefficient is O_(r,f)(q^-1/2), hence projection onto L arbitrary columns has
normalized energy O_(r,f)(L/q). The empirical input patterns tend to uniform
independent signs by the same fixed-degree character-sum estimate, so the
total nonlinear energy tends to its Walsh Parseval energy. This sharply
contrasts with the diagonal finite-spike routing possible in Walsh bases.

WHAT IS MISSING: no theorem turns this fixed-complexity nonclosure statement
into a uniform upper bound on the entire spin-summed partition. One needs an
inverse pressure theorem: every source family with pressure within epsilon
of E must admit a finite-complexity latent description whose conditional-mean
energy is concentrated in finitely many actual columns, with quantitative
loss otherwise. Fixed-moment CLT, fixed-degree Weil bounds, and typical-spin
Gaussianity do not supply such an exponential-scale inverse theorem.

There is also a parameter obstruction. At p31/32,t4 the exact new phase
theorem gives E=g4(1), attained by an uninformative channel. The universal
Gaussian/BSC routing lower theorem applies to Paley bases too. Therefore
nonlinear-mean delocalization CANNOT strictly improve below E there. A
potential strict effect must occur in a genuinely nonlinear latent phase or
in a different construction/pressure architecture.

An elementary diagnostic for a nonlinear phase is p<1/3. If a posterior
mean M is a finite linear combination of independent Rademachers and
|M|<=1/sqrt(p), then E M^4<=3(E M²)²<=3. For normalized ternary X,
E X^4=1/p. Since |x^4-y^4|<=4p^(-3/2)|x-y| on this interval, any such
posterior representation has distortion D=E(X-M)² satisfying
`D>=p(1-3p)²/16`. Thus this restricted affine-label family cannot approach
zero distortion when p<1/3, while full revelation is always a finite-information
channel. At sufficiently large t its g_t(D) penalty diverges to minus infinity,
whereas E stays at least -H(nu_p). This proves a gap for that RESTRICTED
affine-label model, not an actual Paley pressure gap.
