# Preserved director reconnaissance: unresolved comparisons

2026-09-06, research notes retained without polishing. **Not a proof of
convergence and not a new canonical route.** Exact identities are explicitly
separated from discarded comparison hypotheses.

## Threshold replicas: an exact identity and an unproved replacement

For one standard Gaussian threshold bit `X`, put
`u(g)=(X(g)-m)/sqrt(v)`. Draw two Gaussians independently conditional on
the same threshold bit. Their joint law relative to two independent
standard Gaussians has density `1+u(g)u(h)` (verify separately on the
two equal cells and on unequal cells). Let this law be `nu_delta`.
If `L_R` is the Gaussian likelihood at covariance `R`, exactly

```math
 1+\chi^2(p_R,\mu_m)
 =\mathbb E_{\nu_\delta^{\otimes N}}L_R(G)L_R(H).
```

Replacing this pair law by a Gaussian pair of covariance `eta=a^2/v`
would give `det(I-eta^2 s^2 B^2)^(-1/2)` at `R=I+sB`.
The Gaussian integral is elementary: in the symmetric/antisymmetric
replica coordinates its two determinant factors reduce to
`det(I-eta sB)` and `det(I+eta sB)`. **No valid comparison in the required
direction was proved.** This is a reference calculation, not a formula
for threshold entropy.

Entrywise flatness of `B` does not imply entrywise flatness of `B^2`:
duplicated Hadamard rows produce macroscopic twin Gram entries while
the original entries vanish with order. Thus connected Hermite diagrams
with repeated local branches cannot be discarded just by saying that
all covariance entries are small. A fixed-order Taylor calculation also
does not control the fixed-parameter order limit. These problems motivated
the subsequently proved OU/projection/resampling inequality, which avoids
both replacements rather than assuming either.

Rough exploratory guide, not a certificate: at `delta=2^-24`, Gaussian
threshold squared location is about 28. An unattenuated correlation near
`1/sqrt(28)` could produce a useful energy credit only if the discrete
KL multiplier were small. This estimate suggested the question; the
actual proof instead uses exact `s=1/32,r=1/25` and certified inequalities.

## Banach-space literature mapping: no convergence import

Primary source read: Defant, Galicer, Mansilla, Mastylo, Muro,
[arXiv:2302.00233v2](https://arxiv.org/pdf/2302.00233), introduction,
definitions in Section 2 and the stated projection/Sidon results.
For the complete quadratic Walsh support, the projection constant is
the expected absolute all-positive quadratic sum. Its normalized limit
follows from the ordinary CLT and is `sqrt(2/(pi e))` after division by
the number of variables. This is an averaged all-positive kernel, not
the minimum signing supremum. The Sidon constant instead optimizes
arbitrary real Fourier coefficients and is comparable to a degree-one
projection constant only up to fixed factors.

Writing `S_n` for that real quadratic Sidon constant gives only
`M_n >= binom(n,2)/S_n`. Neither asymptotic equality with flat coefficients
nor convergence of `S_n/sqrt(n)` follows from the inspected theorems.
Replacing the original problem by either assertion would insert an
unproved obligation, not import a known solution. No such route is promoted.

## Other literature candidates not used as theorem inputs

- Attard, Jepps, Marcelja (1997), Phys. Rev. E 56, 4052: search located
  a correlation expansion for entropy; no dimension-uniform threshold
  remainder was checked, so it is not a proof dependency.
- A direct publisher PDF retrieval of Finner's 1992 generalized Holder
  paper failed. The weighted information note instead contains a complete
  finite-product Holder proof; it does not depend on an unread theorem.
- Generic Gaussian information inequalities and fixed common-input
  threshold models do not by themselves provide the rare-threshold
  factor required here. This is the bounded search's outcome, not a
  claim that the literature contains no related result.
