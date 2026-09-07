# Independent reconstruction of the pinned G/P Gaussian-pair sector

Scope update: all kernel audits remain valid, but the fixed-center
energy-window conclusion is dominated by the elementary deterministic
bound in `flatify_adversary_2026_09_07_dephased_hadamard_energy_window_audit.md`.

2026-09-07. Final source
`flatify_independent_2026_09_07_mixed_paired_bridge_sector.md` has now been
read completely: PASS. The covariance, exact-law coupling, pin-derived
physical entropy, and same-ensemble type union all pass reconstruction.
The separate independent checker `flatify_adversary_2026_09_07_mixed_pinned_certificate.py`
also verifies all rational determinants and directed margins.

Use alpha=4/5, s=3/5, t=10/3, lambda=35/18, and O=H2/sqrt(2). The
exponent is t a^T K b-lambda(||a||^2+||b||^2). A regular G pair has
covariance s^2 I2; a regular P pair has covariance diag(2s^2,0), with
its active coordinate fixed in bit zero. Set

    a0=1+2lambda*s^2=12/5, aP=1+4lambda*s^2=19/5.

## Regular Gaussian factors

On one edge the exact factors are

    GG: D^-1,       D=a0^2-t^2*s^4=108/25;
    GP: E^(-1/2),   E=a0*(a0*aP-2t^2*s^4)=1872/125;
    PP: F^(-1/2),   F=aP^2-2t^2*s^4=289/25.

For GP the scalar nonzero P coordinate couples to a unit vector on the
G side, giving cross-variance 2s^4. For PP the active-active entry of O
is 1/sqrt(2), giving the SAME cross-variance 2s^4. These are matrix/pair
calculations, not scalar magnitude substitutions.

## Pins and arbitrary center means

Pin each center's unique Walsh coefficient to one designated column pair.
The pinned target mean is mu=alpha*sqrt(k) in bit zero, while the other
coordinate has Gaussian variance s^2 for G and zero for P. Point every
center pin into a GOOD opposite fibre and forbid reciprocal pins. The
opposite pair on such an edge is therefore an ordinary isotropic G pair.

Integrating it first gives the exact mean factor

    exp[mu^2*(t^2*s^2/(2a0)-lambda)]
      =exp[-(10/9)mu^2].

Orthogonality of O kills the mean-partner cross term, so this is independent
of the pinning row's own G/P type. There are 2m nonreciprocal pinned
edges. Since their squared means sum to 2alpha^2 n, the physical norm
compensation 2lambda*n plus the mean penalty is exactly (37/15)n.
This equals t*(37/50)n.

Relative to the ordinary edge determinant, replacing the center variable
by its fixed mean costs sqrt(9/5) on G->G and sqrt(13/5) on P->G.
For example the latter ratio is sqrt(E)/a0=sqrt(13/5); the former is
D/[a0*sqrt(9/5)]=sqrt(9/5). Thus all pins cost at most m log(13/5)
in the log moment. This is sublinear in n and has no hidden large-mean
factor. Replacing t by -t leaves these factors unchanged and pays the
other bridge polarity.

## Entropy and heterogeneous fractions

Write eL,eR for the fractions of P fibres on the two sides. Relative to
the all-G log moment, the leading increase per n is

    A*(eL+eR)/2+B*eL*eR,
    A=log D-(1/2)log E,
    B=-(1/2)log D+(1/2)log E-(1/4)log F.

Separate typed noise windows with fraction delta=1/10 among G single
variables and among P paired variables give physical pair entropy

    [2-(eL+eR)/2]*h(delta)*n+o(n).

This deficit is NOT justified from a bare global Hamming sphere: groups
of sizes one and two can have different optimizing flip frequencies.
The explicit typed windows, or an explicit implication from sufficiently
small pin errors, are essential. The dominant paired-exception family
constructed earlier does satisfy these typed windows.

Here A=.11003094..., B=.00970904..., so B>=0 and h(delta)-A>B. Therefore
the entropy-versus-moment margin for every eL,eR in [0,1] is at least the
all-G margin I0-2h(delta)>.08, with I0=(1/2)log(108/25). Indeed
eL*eR<=(eL+eR)/2. A directed certificate should verify these elementary
logarithmic inequalities, not rely on the displayed decimals. Both the
source's checker and our independent checker now do so.

The SAME source also certifies threshold 18/25=.72: reducing 37/50 by
1/50 spends only t/50=1/15 of the margin, leaving more than 1/75.
This is below the c=47/100 parent bridge budget
(94/100)*(sqrt(2)-16/25). Thus the mixed sector can be paid at that same
seed constant without a new source optimization.

## Actual construction-law and type-union requirements

Form active/inactive Walsh column pairs in every row in advance, even
when its declared type is G. Fix the center pair. Randomize only the
remaining pair labels and independent column signs. The relevant regular
profile is the empirical law of ORDERED absolute pairs, not the separate
scalar histogram and not the full arbitrary column-permutation orbit.
Match this pair law to independent Gaussian G/P pairs; then apply an
independent uniform signed-pair permutation simultaneously to both words.
This gives exactly the physical construction law and the conditioned
Gaussian law. P's zero second coordinate is preserved. The center and
partner pin errors must be paid separately in squared Euclidean norm.
A joint row event controlling the empirical regular-pair W2 error and
the O(1)-variance pinned partner has probability at least 1/2 for large k.
Its conditioning cost over all rows is at most 2m log2. Discard conditioning
before Gaussian integration. Fixed positive margins absorb all o(n)
transport errors and these O(m) finite factors.

For a genuine simultaneous type union, reserve two G fibres a,b on the
left and c,d on the right in advance. Set right pin destinations psi(j)=a
except psi(c)=b; set left destinations phi(i)=c except phi(b)=d. There
are no reciprocal pins, and every destination is reserved G. This layout
does not depend on any other type labels. The SAME paired-frame ensemble
then supports all 2^(2m-4) G/P assignments on the other rows, costing only
O(m) additional entropy. The unequal-fraction formula above covers all
such assignments. This does not union over translation constraints or
over pairings selected after observing candidate spins.

## Scope

Subject to these explicit finite coupling and typed-window conditions,
the operation pays a real family of non-Gaussian paired exceptions rather
than deleting them from the count. It still leaves other profile laws,
the reserved fibres' possible exceptional states, unequal noise windows,
and the broader whole-parent obligation uncontrolled. It keeps child
internal edges fixed and therefore says nothing negative about broader
flatification operations that change those edges.

## Final-source details checked

The final source's pin error gives average (alpha_i-alpha)^2<=tau^2.
Binary entropy continuity and Jensen therefore derive the typed entropy
bound directly, with error 2h(tau/2)n and overlap enumeration 2m log(k+1).
No separate typed-window assumption is needed once that pin constraint
is present. The Euclidean matching error eta=tau+2epsilon is conservative:
one term pays empirical pair matching and the other the reference partner.
The exponent error is t+2lambda=65/9 times eta(2+eta)n, as stated.

For typicality of ordered PAIR profiles one uses a two-dimensional
marginal replacement and a four-dimensional two-pair replacement to
control empirical test variance. The source's one/two-coordinate wording
should be read with this extension; fixed dimension four follows from
the same bounded-coefficient argument. No theorem depends on a growing-
dimensional replacement or an exponential typicality assertion.

## Stronger root upgrade: no approximation of P residuals

Root's subsequent stronger operation passes independent reconstruction.
Keep t=10/3 and lambda_G=35/18, but take lambda_P=5/4. A P row now
requires only its exact translation constraint and small pinned-mean
error; there is NO empirical Gaussian condition on its other coefficients.
In the comparison vector replace only its pinned coefficient by
alpha*sqrt(k), leaving every residual coefficient and its random signs
and pair labels exactly as in the actual physical word.

Condition on ALL these P words and their construction randomness first.
The G comparison Gaussian coordinates remain independent. Every regular
GP edge integrates to

    (5/12)*exp[-(5/12)u^2] <= 5/12,

where u is its arbitrary actual active P coefficient. Every PP edge has
tilted exponent

    (t/sqrt(2))*uv-lambda_P*(u^2+v^2) <= 0,

since t/(2sqrt(2))<lambda_P. Thus its integrand is at most one pointwise.
These estimates do not require independent P coordinates, a P moment
bound beyond the exact physical norm, or any P-profile typicality.

On a P-pinned edge, Gaussian integration into the reserved G neighbour
still has mean coefficient t^2 s^2/(2a0)-lambda_P=-5/12. The regular
constant is 5/12, so there is no additional P-pin determinant ratio.
The G-pin ratio remains sqrt(9/5). Relative to a G row, changing the
norm compensation and its pinned mean together changes the log moment by

    -(lambda_G-lambda_P)*(1-alpha^2)*k = -k/4.

The independent check is important: the center mean is not discarded,
and the whole changed norm penalty cannot be claimed as a saving.

Distribute this saving and the pin-derived physical entropy deficit over
the affected edges. Relative to the all-G entropy-plus-moment expression,
one GP edge contributes

    log(9/5)-1/2-h(1/10) < 0,

and one PP edge contributes

    log(108/25)-1-2h(1/10) < 0.

Hence every type mixture has at least the all-G margin, including at
bridge threshold .72. The same reserved-G pin map and same-ensemble union
apply. Only G rows need the ordered-pair Gaussian coupling and pinned
partner error. P rows need their exact support constraint and mean-pin
error. Both bridge polarities are covered because the integrated terms
depend on t^2 and the PP estimate holds for either sign of t.

This strictly strengthens the earlier Gaussian-P theorem: it pays all
spectral exceptions within the exposed P support, not only typical P
spectra. It still does not expose arbitrary constraints chosen after
seeing a candidate, nor control arbitrary unstructured G exceptions.
