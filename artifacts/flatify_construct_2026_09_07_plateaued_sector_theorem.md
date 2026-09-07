# A uniform strict-subhalf theorem for exact flat-spectrum row sectors

Date: 2026-09-07. Constructive-track proof; independently reconstructed by
the adversarial track. This is a proved sector bound for the actual rank-two
construction, not a cap bound on its entire Boolean cube.

## 1. General Hadamard theorem

Let k=2m=2^d, N=mk, and use arbitrary Hadamard frames F_i of order k in
the rank-two cross construction. Column signs/permutations and outer edge
signs are arbitrary. Complete each fibre by actual signs of cap O(k^(3/2));
the total internal contribution is O(N^(5/4)). For a row spin x_i define
w_i=F_i^T x_i/sqrt(k), so ||w_i||^2=k.

Call a row EXACT-FLAT if all its nonzero transformed coefficients have
one common absolute magnitude. If the unnormalized magnitude is M_i and
there are r_i nonzero coefficients, then

    r_i M_i^2=k^2.

Since M_i is an integer and k is a power of two, M_i is a power of two.
Hence normalized amplitudes A_i=M_i/sqrt(k) have pairwise ratios in2^Z.
This arithmetic conclusion uses only the sign Hadamard identity, not a
Walsh-specific classification of plateaued functions.

THEOREM. Uniformly over every spin whose rows are all exact-flat,

    |H_C(x)| <= (sqrt(2)/3) N^(3/2)+O(N^(5/4)).             (1)

The coefficient sqrt(2)/3=.4714045207910317 is strictly below one half.
Different rows may have completely different support sizes and magnitudes;
they may vary with k. The bound is deterministic, so it holds simultaneously
for all such rows for every realization of the construction.

PROOF. Write the two signed feature coordinates at one endpoint of an
edge as a and at the other as b. The cross energy normalized by sqrt(N)
is exactly a^T R b, R=H2/sqrt(2). Every nonzero coordinate of a has
magnitude A and every nonzero coordinate of b has magnitude B, with A/B
dyadic. If either endpoint is zero, the assertion is trivial. Let u,v in
{1,2} count active coordinates. The maximum absolute UNNORMALIZED H2
bilinear numerator is AB for (u,v)=(1,1), and2AB for the other three cases.

For (1,1) and (2,2), dividing by the four-square budget gives at most1/2
before the factor1/sqrt(2). For (1,2), the needed inequality is

    2AB <= (2/3)(A^2+2B^2),

equivalent to (A-B)(A-2B)>=0. It holds because no power of two lies
strictly between1 and2. The (2,1) case is symmetric. Therefore

    |a^T R b| <= (sqrt(2)/3)(||a||^2+||b||^2).             (2)

Every cross feature coordinate appears in exactly one unordered edge.
The sum of their squared values is at most sum_i||w_i||^2=N; the two
unused loop coordinates per fibre only reduce that budget. Sum (2), then
apply the actual internal completion bound. This proves (1).

If every row is fully flat (no zero coefficients), its normalized magnitude
is1 and every edge has active counts(2,2). The stronger coefficient is
1/(2sqrt(2))=.3535533905932738. This includes the all-bent Walsh sector.

## 2. Quantitative sector with some unrestricted rows

Suppose at most a fraction eta of the fibres are not exact-flat. Put
c0=sqrt(2)/3. The normalized cross cap is at most

             c0+eta/(12sqrt(2)),             0<=eta<=8/17;
    q(eta) = sqrt(eta(1-eta)),              8/17<=eta<=1/2;
             1/2,                          1/2<=eta<=1.    (3)

To verify this, measure feature-square mass in units of N. Good rows have
mass1-eta and bad rows masseta. Let G be the good-coordinate mass assigned
to good-bad edges, and B the bad-coordinate mass assigned to those edges.
The remaining good-good contribution is at most c0(1-eta-G), the cross
contribution is at most sqrt(GB) by Cauchy--Schwarz, and bad-bad is at most
(eta-B)/2. Thus maximize

    c0(1-eta-G)+sqrt(GB)+(eta-B)/2,
    0<=G<=1-eta, 0<=B<=eta.

For eta<=8/17 the maximizer is B=eta,G=9eta/8. For8/17<=eta<=1/2 it is
B=eta,G=1-eta. For eta>=1/2 choose B=G=1-eta and the value is1/2.
This gives (3), again with only O(N^(5/4)) added for actual completion.

Consequently, a spin exceeding normalized coefficient .493608093588748653
must have at least .3768071254040348+o(1) non-exact-flat fibres. This is a
restriction on possible large-energy witnesses, not an exclusion of them.

## 3. A separate metric extension, with NO approximate-spectrum counting

Let v be any comparison feature vector with norm sqrt(N), whose rows have
the exact dyadic flat structure used in (2) and norm sqrt(k) each. It need
not be the transform of a Boolean row. If ||w-v||^2<=epsilon N, then

    |H_cross(w)|/N^(3/2) <= c0+sqrt(epsilon).              (4)

Indeed the normalized cross operator T has norm1 and
H_cross(w)/sqrt(N)=w^T T w/2. Its quadratic difference is at most
||w-v||(||w||+||v||)/2<=sqrt(epsilon)N. This metric statement is proved
directly; no Potapov counting theorem is applied to approximate spectra.

## 4. Actual primary counting bounds and how they enter

For Walsh frames, exact-flat rows are plateaued. Potapov,
https://arxiv.org/html/2303.16547 (v3,18November2024), Theorems1(a),2, gives

    log2 #bent(d) <= (11/32+o(1)) k;
    log2 #s-plateaued(d)
      <= [alpha/8+(h2(2^-s)+2^-s)/4+o(1)] k,             (5)

where alpha=1+(3/8)log2 6, h2 is binary entropy, and s>0 is FIXED as
d tends to infinity. The same source gives degree <=(d-s)/2+1, hence a
leading1/2 bound for fixed positive s. Theorem1(b)'s sharper near-bent
count applies only to restrictions of bent functions, not all near-bent rows.

Combining (5) with the elementary support count proved in the rank-two
record gives admissible entropy coefficients (in bits per coordinate):

    e0=11/32;
    e_s=min{1/2, alpha/8+(h2(2^-s)+2^-s)/4, h2(2^-s)}.

In particular e1=e2=1/2 and e3=.4133112280212033. For a FIXED finite set
of levels and row fractions p_s, the log row-tuple count per N is at most
(log2)sum_s p_s e_s+o(1). Choosing which fibres have each type costs o(N).
Growing s must NOT be inserted into the fixed-s asymptotic in (5); the
uniform elementary support bound handles the large-s tail instead.

These are legitimate row-count inputs to the director's finite typed
certificate. However (1) is stronger for excluding the entire exact-flat
sector at thresholds above c0: its probability of producing such a witness
is ZERO, so no entropy payment is necessary there. The primary counts are
retained for any attempted finer threshold or mixed unrestricted sector.

## 5. Exact remaining obligation

The all-exact-flat heterogeneous sector is fully controlled with a uniform
strict-subhalf coefficient. This is a genuine proof step, not a named
condition or numerical proxy. The remaining parent cube includes rows
with two or more nonzero magnitudes, including typical Gaussian-like
Walsh profiles. To improve the current global coefficient using this route,
one still must control witnesses with a positive fraction of those rows,
or quantify their distance from the exact-flat comparison class.

No claim is made that every near-ground state is close to plateaued, that
Potapov's exact counts apply to approximate spectra, or that the new sector
theorem transfers the input seed cap. Independent column signs still erase
the outer seed in the randomized construction.

A concrete limit of the local argument: a=(1,4), b=(4,-2) give
|a^T R b|/(||a||^2+||b||^2)=26/(37sqrt(2))=.496886... . Thus simply
requiring all four displayed coefficients to be dyadic does not extend
the theorem. Such pairs occur in actual order64 Walsh transforms: flip16
entries of the all-positive row, split6/10 between the two level sets of
a nontrivial character to obtain coefficients(8,32), or split12/4 to
obtain(32,-16), and normalize by8. This is only a local compatibility
example; it is not a global bad parent construction.
