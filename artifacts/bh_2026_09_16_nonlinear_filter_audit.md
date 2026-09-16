# Nonlinear polynomial filters: what polynomial BH can and cannot detect

2026-09-16. Status: DIRECTOR DERIVATION; weighted, coefficient-filter,
and rare-tail arguments have independent power-agent audit PASS. The
independent paper verifier also checked the all-signing small-degree
extension and the full rare-tail argument in Section7. This is a scoped obstruction, not a theorem
about the moment law of actual exact minimizers.

## 1. Setting and elementary inputs

Let A_n be a hollow symmetric full signing with
`||A_n||op<=L sqrt(n)` for a FIXED L. Write
`H_A(x)=sum_(i<j) a_ij x_i x_j`, `Y_n=H_A/n^(3/2)`.
Such actual signings exist at every order: take a principal submatrix of
the next Sylvester matrix and remove its diagonal. A standard elementary
decoupling/Gaussian-mgf argument, reconstructed in the companion power
artifact, gives, uniformly j>=1,

    ||Y_n||_(2j)<=C_L(sqrt(j/n)+j/n).                 (1)

For j<=n this is at most `2C_L sqrt(j/n)`.

Fix c>0. Let P be ANY real or complex univariate polynomial of degree d
with `sup_(t in [-c,c])|P(t)|<=1`. No positivity of its coefficients,
monotonicity, or Chebyshev form is assumed. Writing P(t)=sum a_j t^j,
Bernstein--Walsh and Cauchy's coefficient formula imply

    |a_0|<=1,    |a_j|<=(e d/(c j))^j,   1<=j<=d.    (2)

Indeed the Green-function bound on |z|=R is
`|P(z)|<=(R/c+sqrt(1+(R/c)^2))^d<=exp(dR/c)`.
Optimize Cauchy's bound at R=cj/d. The primary-paper audit also gives
an elementary proof of this complex-plane estimate by the Joukowski map.

## 2. Weighted square-function theorem

For a Walsh polynomial f on n spins, put

    S(f)^2=sum_(r>=1) r^(-10)
                (sum_(|T|=r)|fhat(T)|^(2r/(r+1)))^((r+1)/r).

This is a seminorm. Layerwise norm comparison and Parseval give

    S(f)^2<=e n sum_(T nonempty)|fhat(T)|^2
           <=e n ||f||_2^2.                          (3)

More generally any fixed nonnegative layer weights bounded by one give
the same upper bound. Using (1)--(3), for d<=n,

    S(P(Y_n))
    <=sqrt(e n) sum_(j=1)^d (z/sqrt(j))^j,
    z=C_(c,L) d/sqrt(n).

The constant term drops out. Since `j^(j/2)>=sqrt(j!)`, Cauchy--Schwarz
with geometric weights proves

    sum_(j>=1) z^j/sqrt(j!)<=sqrt(2) z exp(z^2).

For example, divide out z, use `(l+1)!>=l!`, and apply Cauchy--Schwarz
to `2^(-l/2)` and `(sqrt(2)z)^l/sqrt(l!)`. Consequently

    S(P(Y_n)) <= C_(c,L) d exp(C'_(c,L) d^2/n).       (4)

The SAME conclusion holds for the improved r^(-5) square weights in the
independent audit: the layer factor becomes r^(-6)<=1 in (3).

In particular, for growing d with `d^2/n=o(log d)`, the right side of
(4) is `d^(1+o(1))`. Neither the paper's bound `K(2d)^22` nor the
reconstructed stronger bound `K'(2d)^3` can be contradicted by these
filters. This includes all growing `d=o(sqrt(n log n))`.

## 3. Full coefficient-norm theorem

Suppose `2<=d<=n/4`, and set `p=4d/(2d+1)`. The Walsh support of Y_n^j
has at most `D_(n,j)=sum_(r<=2j) binom(n,r)` elements, hence

    ||hat(Y_n^j)||_p <=D_(n,j)^(1/(4d))||Y_n||_(2j)^j
                    <=(e n)^(j/(2d))||Y_n||_(2j)^j.

This uses only the Boolean symmetric-difference product; all cancellations
are retained inside each final Fourier coefficient. Triangle inequality,
(1),(2), and the previous exponential sum now give

    ||hat(P(Y_n))||_p <=1+sqrt(2) z_d exp(z_d^2),
    z_d=C_(c,L) d (e n)^(1/(2d))/sqrt(n).             (5)

Uniformly for every sequence `2<=d=o(sqrt(n))`, z_d tends to zero.
To see this, split at d=log n. Below the split, d>=2 gives
`z_d<=C log(n) n^(-1/4)`; above it, `(en)^(1/(2d))` is bounded and
d/sqrt(n) tends to zero. Thus

    sup_(P as above) ||hat(P(Y_n))||_(4d/(2d+1))
                                                    <=1+o(1). (6)

The lower bound1 is attained by the constant polynomial if degree means
at most d. If degree is EXACTLY d, arbitrary small nonzero leading
perturbations approach it. For a full signing with 2d<=n, a nonzero
leading coefficient produces a nonzero top Walsh level by the odd-hafnian
identity; degree cannot silently collapse to two.

For growing d=o(sqrt(n log n)), (5) more generally is `d^(o(1))`:
when d>=log n, the factor `(en)^(1/(2d))` is bounded, and
`log(1+sqrt(2) z_d exp(z_d^2))<=C(1+z_d^2)` for all z_d>=0.
The latter is `O(1+d^2/n)=o(log d)`. Below log n use (6).
In particular it remains negligible compared with either degree power
27 or11/2. The sharper statement (6) is the primary claimed conclusion.

## 4. Scope: a failed certificate is not a failed optimization theorem

If one assumes `Q(A)<=c n^(3/2)`, then `||P(Y_n)||infty<=1`, so the
paper supplies necessary coefficient inequalities. Equations (4)--(6)
show that these inequalities alone cannot detect the violation of that
cap assumption on the bounded-operator control family at the stated
filter degrees. Their coefficients are too small, even though a rare
actual extremum can lie outside [-c,c]. Signs in P and Walsh reduction
do not rescue the filter.

This does NOT prove actual minimizing signings have bounded operator norm,
does NOT exclude an optimizer-specific high-order coefficient theorem,
and does NOT exclude filters of linear or larger degree. Such a theorem
would be additional mathematical input, not supplied by polynomial BH.
No lower or upper constant for M_n is improved here.

## 5. Explicit full-sign false positives

The preceding conclusion can be tested against a cap hypothesis that is
actually false, without assuming anything about minimizers. Let

    R=J_4-2I_4,  H=R^(tensor k),  n=4^k,
    delta=(-1)^k,  A=H-delta I_n.

Then A is symmetric, hollow and full-sign, `H^2=nI`, and
`||A||op=sqrt(n)+1`. In fact its cap is EXACTLY

    Q(A)=n(sqrt(n)+1)/2.

The spectral norm gives the upper bound. For odd k, the all-ones vector
is a Boolean eigenvector of A with eigenvalue sqrt(n)+1. For even k,
take a balanced Boolean eigenvector of R with eigenvalue -2 in one
factor and all-ones vectors (eigenvalue +2) in the remaining factors.
Their tensor is a Boolean eigenvector of A with eigenvalue -sqrt(n)-1.
Both cases attain the displayed cap.

Thus `Q(A)/n^(3/2)->1/2`. For every fixed c<1/2 the hypothesis
`Q(A)<=c n^(3/2)` is false on this family for all large k. Nevertheless
ALL normalized polynomial filters of growing degree
`d=o(sqrt(n log n))` satisfy the necessary coefficient inequalities
supplied by the paper's degree growth (and by the refined exponents in
the independent audit). The obstruction concerns precisely those
certificates; it is not an obstruction to stronger, model-specific
coefficient inequalities or to filters at larger degree.

The order limit in this statement is along n=4^k. No all-order
realization or assertion that these signings are near-minimizing is
needed. The finite replay uses exact integer Walsh transforms of
Chebyshev filters, with floating-point coefficient norms clearly
separated from exact identities.

## 6. Stronger all-signing small-degree theorem

The low-degree assertion does NOT in fact need an operator assumption.
For every full signing, degree-two hypercontractivity and Parseval give

    ||Y||_(2j) <= (2j-1)||Y||_2
               <= sqrt(2) j/sqrt(n).

Combined with (2), put `z=sqrt(2)e d/(c sqrt(n))`. Whenever z<1,

    W_s(P(Y)-a_0-a_1Y) <= sqrt(en) z^2/(1-z),

where `W_s^2=sum_(r>=1) r^(-2s) S_r^2`, for every fixed s>=0.
In particular for d=o(sqrt(n)) this error is `O_c(d^2/sqrt(n))`,
uniformly over ALL full signings and ALL interval-normalized P.
The reverse triangle inequality for this seminorm proves

    W_s(P(Y))
      =2^(-s)|P'(0)| binom(n,2)^(3/4)/n^(3/2)
                         +O_c(d^2/sqrt(n)).

Thus low-degree filters have asymptotically only a linear coefficient
response even on signings with very large cap. Bernstein's interior
inequality `|P'(0)|<=deg(P)/c` converts this into a sharp normalized
efficiency theorem; see the independently written
`bh_2026_09_16_linear_weighted_filter_analysis.md` for the quantifiers,
the exact degree check, and the distinction between an all-degree
linear bound and an asymptotic slope with an additive intercept.

The full coefficient estimate improves in the same way. For
`p=4d/(2d+1)` let `z_d=sqrt(2)e d(en)^(1/(2d))/(c sqrt(n))`.
Whenever z_d<1,

    ||hat(P(Y))||_p <=1+z_d/(1-z_d).

As in Section3, z_d tends to zero uniformly along every sequence
`2<=d=o(sqrt(n))`. Equation (6) therefore holds for ALL full signings,
not just the bounded-operator family.

## 7. Rare-tail extension to EVERY sublinear filter degree

The preceding Taylor estimates are deliberately local. For larger d,
one must keep cancellations in P rather than separately pay its Taylor
terms. On bounded-operator full signings this gives a much stronger
conclusion.

**Theorem.** Fix L,c>0. There is eta=eta(L,c)>0 such that the following
holds for every full signing with `||A||op<=L sqrt(n)`, every polynomial
P of its ACTUAL univariate degree d with `||P||_[-c,c]<=1`, and every
sequence `d->infinity`, `d<=eta n`:

    W_s(P(Y))/d^a ->0               for every fixed a>1, s>=0;
    ||hat(P(Y))||_(4d/(2d+1))/d^b ->0  for every fixed b>1/2.

Both conclusions are uniform over A and P. In particular they hold for
every growing sublinear degree d=o(n). The actual reduced Walsh degree
is 2d: choose eta<=1/4 and use the nonzero odd-hafnian top coefficient.

**Proof of the tail step.** The spectral bound gives `|Y|<=R=L/2`.
The elementary decoupling/mgf proof in the power artifact yields a
constant gamma=gamma(c,L)>0 such that

    Pr(|Y|>c) <=2 exp(-gamma n).

The normalized interval bound and the Joukowski estimate give, on the
whole possible energy interval,

    |P(t)|<=exp(K d),  |t|<=R,
    K=log(R/c+sqrt(1+(R/c)^2))>0.

Consequently, keeping P intact rather than summing Taylor magnitudes,

    ||P(Y)||_2^2 <=1+2 exp(2Kd-gamma n).

Choose `eta=min(1/4,gamma/(4K))`. For d<=eta n the error is at most
`2 exp(-gamma n/2)`. This also covers c>=R, where the tail is actually
empty. The constants are fixed before n and d are sent to infinity.

**Weighted conclusion.** Choose
`1/(2a)<theta<1/2`. If `d<=n^theta`, (4) gives
`W_s(P(Y))<=C d exp(C'd^2/n)=O(d)`, uniformly in this range.
Dividing by d^a tends to zero when d grows. If `d>=n^theta`, use
`W_s<=sqrt(en)||P(Y)||_2`; its ratio to d^a is
`O(n^(1/2-a theta))=o(1)`.

**Full coefficient conclusion.** Choose
`1/(2b+1)<theta<1/2`. For d<=n^theta, Section6 gives coefficient
norm `1+o(1)` uniformly once d>=2. For d>=n^theta, the support has at
most `sum_(r<=2d)binom(n,r)<=(en/(2d))^(2d)` elements because
`2d<=n/2`. Norm comparison therefore gives

    ||hat(P(Y))||_(4d/(2d+1))
       <=sqrt(en/(2d)) ||P(Y)||_2.

Its ratio to d^b is `O(n^(1/2-(b+1/2)theta))=o(1)`.
This proves both statements.

**Exact consequence and scope.** For any fixed positive constants in
the paper's weighted degree bound22 and full degree bound27, these
necessary BH inequalities are eventually satisfied on the explicit
false-positive family of Section5 for EVERY growing sublinear-degree
filter. The same is true for the refined exponents3 and11/2. Thus the
gap between square-root and linear filter degree is not a loophole.
Both older subexponential and newer polynomial constants still leave
the rare-event obstruction intact. The theorem even covers a small
fixed linear-degree interval, but does not cover arbitrary linear
degrees. Nor does it address stronger model-specific inequalities,
all possible uses of the new proof mechanism, or actual minimizer laws.

The conclusion charges P its own degree. It does not falsely charge a
linear filter the much larger cutoff degree of a family containing it.
Bounded-degree tests remain their separate finite-degree inequalities.

## 8. Finite replay and evidence status

Run:

```sh
.venv/bin/python computations/bh_2026_09_16_filter_checks.py --output computations/results/bh_2026_09_16_filter_checks.json
```

The replay verifies24 Chebyshev filters on four actual sign matrices:
regular Hadamard orders4 and16, the stored order9 cap12 witness, and
the order16 all-positive matrix. Exact integer checks include energies,
every Walsh transform, Parseval, parity, and nonvanishing top degree
when2d<=n. The rational normalization is c=2/5. Matrix inputs, hashes,
energy histograms, and all nonzero coefficient-magnitude histograms
are retained. Coefficient norms are60-digit numerical diagnostics, not
interval certificates. No asymptotic theorem is inferred from these
small orders; the all-signing theorem and rare-tail proof supply that.
