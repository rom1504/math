# Optimized arbitrary-center Gaussian-noise shell rate

Scope update: the fixed-center energy-window application is superseded by
`flatify_adversary_2026_09_07_dephased_hadamard_energy_window_audit.md`,
which controls all profiles and unequal noise levels by an elementary
deterministic bridge. The certificates below remain valid sector results.

2026-09-07. Exact variational reduction; the numerical all-shell maximum
below is exploratory, not a directed continuum certificate.

The arbitrary-center actual-sign theorem of the independent track extends
structurally beyond its rational alpha=4/5, b=37/50 point. Fix 0<alpha<1,
G=[[1,alpha],[alpha,1]], and constrain a tilted four-vector covariance to

    Sigma=[[G,J],[J,G]], J=[[0,z],[z,b]].

The transpose-symmetric form loses nothing when maximizing log determinant:
average any feasible cross-block with its transpose and use concavity.
Put d=1-alpha^2. Then

    det Sigma=(d-z^2)^2-(b-2alpha*z)^2.

The feasible domain is d-z^2>|b-2alpha*z|. The unique maximizing z obeys

    z^3-(1+alpha^2)z+alpha*b=0.

Uniqueness follows from strict concavity of logdet on the affine covariance
slice; the relevant root lies on the decreasing branch of this cubic.
At this stationary point the off-diagonal entries of the cross precision
block vanish, because the corresponding derivative of logdet vanishes.
Thus P=Sigma^-1 has exactly the admissible form
[[G^-1+L,-T],[-T,G^-1+L]] with T diagonal.

Both center marginal variances are one and their cross covariance is
zero. Hence the center Schur complement of P is I. Integrating the two
replica noises is therefore pointwise constant in arbitrary actual center
arrays, exactly as in the rational certificate. Moreover tracing P Sigma
gives Tr(LG)=Tr(TJ)=t1*b. The resulting conditional rate is exactly

    I_alpha(b)=-1/2 log(((d-z^2)^2-(b-2alpha*z)^2)/d^2).

For alpha=0 the formula continuously reduces to
I_0(b)=-1/2 log(1-b^2). The optimal source t1 is positive when b>0:
I increases with z, while db/dz=(1+alpha^2-3z^2)/alpha>0 on the feasible
interval z^2<min(alpha^2,1-alpha^2). For alpha^2<=1/2 the numerator is
at least 1-2alpha^2, and for alpha^2>=1/2 it is at least 4alpha^2-2,
with strictness supplied by strict feasibility when alpha^2=1/2.
This is an actual-sign sector probability bound via the same transport
and conditioning proof, not a substitution of a Haar bridge.

For Hamming fraction delta=(1-alpha)/2, solving
I_alpha(b)=2 h(delta) gives the threshold for a full cloud-pair union.
Typical child energies contribute 2 alpha^2 c, so a restricted parent
improvement at seed constant c requires

    c > b/(2 sqrt(2)-2 alpha^2).

The preserved floating script
`computations/flatify_adversary_2026_09_07_conditional_shell_rate.py`
and result JSON evaluate this equation on a 1000-point grid and locally
optimize its largest peak. They find .4648075500382 at delta=.1192912491;
at delta=.1 they give b=.7137683295208 and c=.4609634629320.
These values are not directed certificates of a uniform maximum.
The first implementation tested determinant positivity alone and failed
outside the positive-definite component; the corrected script explicitly
checks d-z^2>|b-2alpha*z|. A positive determinant alone is insufficient.

Even a certified uniform scalar threshold would cover only the joint
Gaussian-noise profile sectors. Typicality leaves an exponentially large
exceptional family, and no estimate here pays a covering by centers or
those exceptions. In particular no original-value recurrence follows.

## Directed uniform certificate at c=47/100

`computations/flatify_adversary_2026_09_07_all_noise_certificate.py` and its
result JSON now certify the entire SAME-NOISE-LEVEL parameter interval.
This is no longer an extrapolation of the preceding floating maximum.
The middle interval alpha in [59/100,99/100] is covered by 522 closed
rational cells, with endpoints checked to tile that interval exactly.
On each cell one rational z is supplied. Define

    b(alpha)=z*(1+alpha^2-z^2)/alpha.

This makes z stationary EXACTLY for every alpha in the cell, without
requiring any approximate polynomial-root assertion. The useful identities
are

    det Sigma=(d-z^2)^2*(1-z^2/alpha^2),
    I=-log(1-z^2/d)-(1/2)log(1-z^2/alpha^2).

Directed interval evaluation certifies d>z^2 and alpha^2>z^2, and both

    I > 2h((1-alpha)/2),
    b(alpha) < (94/100)*(sqrt(2)-alpha^2)

uniformly in every cell. The script uses a rational LOWER bound on sqrt(2),
checks its square exactly, and accepts only positive rational lower
endpoints for both gaps. Its floating root finder merely proposes z;
all acceptance conditions are directed inequalities. Every rational z,
cell, and lower gap is stored. There were 1043 total subdivision nodes.

For alpha<=59/100, the operator bound |x^T C y|<=n^(3/2) suffices:
the right-hand bridge budget exceeds one by more than .00214. For
alpha>=99/100, write normalized physical vectors as x=alpha*x0+s*u and
y=alpha*y0+s*v with s=sqrt(1-alpha^2), u,v unit and orthogonal to their
respective centers. If the center bridge has normalized magnitude c0,
orthogonality of C/sqrt(n) gives

    |h1| <= c0+2s+s^2.

The latter is at most c0+2sqrt(1-(99/100)^2)+1-(99/100)^2,
strictly below (94/100)*(sqrt(2)-1) by more than .087.
Thus both endpoint intervals are paid by deterministic geometry, not
Gaussian-profile assumptions. Negative overlaps reduce by global spin
flips; the same-level condition means equal absolute overlaps on the
two sides.

For the middle interval, the finite family of compact, positive-definite
precision ranges bounds all transport error coefficients uniformly. The
reference noise variances stay bounded away from zero and the arbitrary-
center empirical transport argument is uniform there. Union over the
at most n+1 exact noise levels costs only O(log n). The positive certified
rate and cap gaps absorb this, the 2m log2 conditioning cost, and uniformly
vanishing profile errors. The one shared center constraint has failure
probability o(1), as before.

Consequently, at every compatible growing Hadamard order there exist
actual full sign bridges which, simultaneously for every same noise level,
control all pairs in its stated arbitrary-center Gaussian-noise sector
(and all pairs at the two endpoint intervals). In the corresponding
typical child-energy windows, for every seed constant c>=47/100,

    |H_A(x)+H_D(y)|+|x^T C y| < (2sqrt(2)c-epsilon) n^(3/2)

for some fixed epsilon>0 once the window errors vanish. This follows from
2c alpha^2 plus the certified bridge budget; increasing c only enlarges
the allowed budget. The centers may be arbitrary actual spins; neither
their coefficient profiles nor their bridge polarity is prescribed.

This pays the continuous common noise-level parameter, NOT unequal noise
levels, exceptional profiles, arbitrary energy shells, or a covering by
centers. It proves neither an all-order child-preserving interpolation
nor an original-value recurrence. No further constant tuning is needed
for this structural conclusion.
