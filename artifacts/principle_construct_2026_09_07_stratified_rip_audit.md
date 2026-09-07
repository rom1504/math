# Independent audit: stratified RIP and near-linear sparse activity

2026-09-07. **PASS**, mathematical derivation coordinated with synthesis,
which independently read the complete original Rudelson--Vershynin lemma
and proof. The actual law is MODIFIED by per-fibre selector conditioning.

Primary input: Rudelson--Vershynin, *On sparse reconstruction from Fourier
and Gaussian measurements*, Lemma3.6,
<https://public.websites.umich.edu/~rudelson/papers/convex-relaxation.pdf>.
An independently checked primary exposition is Tropp et al., *Beyond
Nyquist: Efficient Sampling of Sparse Bandlimited Signals*, Lemmas17--18,
<https://users.cms.caltech.edu/~jtropp/papers/TLDRB10-Beyond-Nyquist.pdf>.
The latter explicitly notes that the independent row distributions may
differ. The deterministic Rademacher estimate has factor at most
C sqrt(s) log^2(m) when coordinates are bounded by1 and row count<=m.

## Exact adaptation, with no identical-distribution premise

Fix a full unnormalized Hadamard H of order m and ANY partition into
M=m/L groups of size L. For a one-hole selector, let X_g be the full
Hadamard row omitted from group g. These vectors are independent, each
has infinity norm1, and aggregate isotropy is exact:

    sum_g E X_g X_g^T = (1/L) H^T H = M I.

Let ||B||_s denote the maximum absolute quadratic form on s-sparse unit
vectors, and put D=E||sum_g X_g X_g^T-MI||_s. Symmetrization is valid
for independent nonidentical random vectors. The deterministic RV lemma
and Jensen give

    D <=2a sqrt(D+M),     a<=C sqrt(s) log^2(m).

Solving yields D<=4a^2+2a sqrt(M). Thus D/M tends to zero uniformly over
H and its partition when s log^4(m)/m tends to zero (L fixed).
Markov gives an error epsilon_m tending to zero and a uniform probability
at least1/2 of ||missingGram-MI||_s<=epsilon_m M. Intersect with the
uniform column-excess good event, whose failure probability tends to zero;
adjust constants so the intersection still has probability at least1/2.

The retained Gram is mI minus the missing Gram, hence its restricted
upper bound is q+epsilon_m M, q=(1-1/L)m. No independent retained-row
assumption is needed.

## Conditioning and exact balance are paid

For each top-level physical fibre, condition its selector GIVEN its full
frame and partition on this RIP/excess event. The conditional selector
density relative to its old law is at most2 uniformly. Consequently every
nonnegative one-row partition expectation is at most twice its previous
one, even though conditioning may couple the child variables within that
row. No conditional child independence is asserted or needed. Across
physical fibres the entire sampling remains independent.

Restricted Gram norms and maximum non-DC column excess are invariant
under available output column signed permutations. Thus their conditional
group law, required by graph Finner, remains intact. The added log2 per
row is o(m) and cannot change a strict leading row certificate.

Use the already proved random exact balancing with the selected raw matrix
fixed. Its operator error relative to P H_retained is o(sqrt(m)), uniformly
under the imposed excess bound, with sufficiently strong polynomial tail
to hold over all fibres. Since P is an orthogonal contraction, for every
s-column submatrix the repaired upper singular value squared is q+o(m).

Applying the deterministic directed-port reciprocal-swap argument now
gives the RELATIVE coefficient sqrt(p)/2+o(1) for all words with at most
s+1 active fibres, uniformly over every outer sign seed and symmetric
mask. Thus any chosen s=o(m/log^4 m) is allowed. In particular every
fixed power-sublinear active-fibre range is covered. The older coherence
argument remains a fully elementary fallback without RV.

This is a real stronger actual ensemble operation, not a claim that the
unconditioned law automatically has simultaneous RIP over all fibres.
It does not control a positive fraction of active fibres or establish
selected-minimizer cap transfer.

The final assembled canonical Section6 was read completely and passes
also for general fixed k<L. Order each group uniformly and discard its
first L-k slots. Each slot separately has independent rows across groups
and the same aggregate isotropy. Its finitely many dependent companions
can be handled by a union and the restricted-norm triangle inequality.
No independence between slots is needed. The canonical proof uses a
near-one good probability, reducing the whole-array density cost further
to exp(o(m)). Its explicit warning that conditioning is top-level only
and need not preserve intra-row child independence is correct.

## Quantitative bounded-density variant (director derivation, audited PASS)

Let theta=s log^4(m)/m tend to zero and use the uniform expectation bound
B_m=4a^2+2a sqrt(M). Conditioning the restricted deviation to be at most
4B_m has failure probability at most1/4. Intersecting the excess event
therefore leaves success probability at least1/2 for all large orders.
Use density at most2 per fibre rather than the near-one normalization.
The entire array density costs exp(O(m)), still subleading to each fixed
accuracy leading certificate.

The restricted Gram error divided by m is O(sqrt(theta)+theta). The
balanced repair operator error divided by sqrt(m) is
O(m^(-1/4) log^(5/4)m), including a smaller mean-contraction term.
Squaring the repaired singular-value bound therefore gives relative
phase coefficient error

    O(sqrt(theta)+m^(-1/4)log^(5/4)m).

For s=m^(1-eta), fixed eta>0, this is
O(m^(-eta/2)log^2m+m^(-1/4)log^(5/4)m). These are paid power-saving
errors for this sparse-active phase, not an original minimax recurrence.

Final assembled canonical Section6.5 read PASS. For d=L-k discarded
slots, triangle bounds the expected TOTAL restricted deviation by d B_m.
Condition total deviation<=4d B_m directly. Markov failure is at most1/4,
without independence between slots; intersecting the excess event gives
the claimed success probability at least1/2. This makes the fixed-slot
constant enlargement completely explicit.
