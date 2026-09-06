# Audit of rare-selector obstructions and anti-sparsity conditioning

Date: 2026-09-06. Full independent read and reconstruction of
`continued_feedback_weave_rare_selector_obstruction_2026_09_06.md` and
`continued_feedback_selector_antisparsity_2026_09_06.md`. Their proved
statements pass; all numerical optimizations remain uncertified probes.

## 1. Rare exact Walsh selectors

At exact dyadic retention p=a/b, a selector pulled back from an a-point
subset of a fixed b-point quotient has at most b nonzero Walsh spectral
coordinates with aligned spins. Fixing every remaining nonzero labeled
coordinate in a permutation leaves (d-r)! permutations, d=m-1 and r<=b.
Each fixed nonzero diagonal kernel entry is at least 1/2, independently
of its magnitude and of t. Thus

```math
\operatorname{per}K/d!\ge2^{-r}(d-r)!/d!\ge2^{-b}d^{-b}.
```

The single selector has probability 1/binomial(m,pm). The resulting
one-row pressure lower bound -h(p) is UNIFORM in t. It therefore rules
out the unconditioned p=15/16 Gaussian-optimal tilt 7.74596669: the
necessary upper threshold is h(p)/(1-sqrt(p))=7.36255134.

The nearby-selector count also checks. Removing alpha m elements of
the base selector and adding alpha m from its complement gives rarity
rate h(p)-p h(alpha/p)-(1-p)h(alpha/(1-p)). Every nonexceptional Walsh
character is balanced in each population. Their two finite-population
variance contributions give

```math
v_p(\alpha)=\frac{2\alpha-\alpha^2/[p(1-p)]}{p}.
```

Distinct character cosets of the fixed annihilator have balanced four-
class populations and independent normal limits; only O(m) exceptional
pairs remain. Fixing all finitely many exceptional spectral coordinates
in the permanent costs only a polynomial factor. With spin flips of
fraction theta in each retained population, the row-count rate is
p h(theta), and the variance becomes
1-(1-2theta)^2(1-v_p(alpha)). This reconstructs the note's full lower
envelope, without treating it as a pressure equality.

The fixed-input-fibre alternative also has the stated count and mixture:
only hit fibres carry distinct spin choices, while zero internal
frequency has occupancy second moment bp(1-p)+b^2p^2 and every other
internal frequency has second moment bp(1-p). The Sinkhorn evaluations
in the note are expressly approximate and are not certified bounds.

## 2. Exact and approximate anti-sparsity

An s-dimensional real subspace has an injective restriction to some s
coordinate evaluations. It therefore contains at most 3^s ternary
vectors. Union over the binomial(m,s) spectral supports and divide by
binomial(m,k) selectors. This proves the exact span count, including
all 2^k supported spin choices without paying a separate 2^k factor.

For the approximate version, the projected vector lies in the radius
sqrt(k) ball of the s-dimensional subspace. A delta sqrt(m) net has
size at most (1+2sqrt(p)/delta)^s. A nearest ternary rounding disagreement
costs at least 1/4 in squared distance from the unrounded center, including
ties. Therefore at most 4(rho+delta)^2 m coordinates differ. The ternary
Hamming-ball count is sum_(j<=rm) binomial(m,j)2^j, with exponent
h(r)+r log2 for r<2/3. All factors in the displayed approximate bound
and its p=15/16 numerical entropy gap consequently check.

The normalization from O=H^T/sqrt(m) to amplitudes a=|H^T xi|/sqrt(k)
is important: the tail energy bound becomes epsilon m with
epsilon=rho^2/p. It is not rho^2 m without that factor.

## 3. Conditioning and the surviving obstructions

Independently drawing fibre selectors from the uniform law conditioned
on a good event preserves the exact fixed-selector soft weave proof.
Its one-row expectation is E[Z_T 1_good]/P(good), exactly as stated.
For P(good)=1-exp(-Omega(m)), the normalization has no exponential cost.
Removing the rare finite-spike selectors is therefore a legitimate
existence-construction choice, not a contradiction with their lower bound.

The TYPICAL Gaussian-profile lower obstruction survives any conditioning
whose probability tends to one. More generally, for a selector class
of probability exp(-I m+o(m)), removal probability at most exp(-c m)
with I<c deletes a vanishing relative fraction of that class. Its
conditional typical-profile lower proof, when its row-choice count is
uniform across selectors in the class, survives too. This does not
assert that an arbitrary exponentially weighted partition contribution
survives merely from a selector mass comparison.

Thus the displayed anti-sparsity event may remove the exact spike class
while leaving the Gaussian p<=0.9225232669 obstruction untouched. The
remaining expectation still needs a genuine upper bound.

## 4. The permanent limitation is genuine

PSD Cauchy--Schwarz gives K(a,b)<=sqrt(K(a,a)K(b,b)); every permutation
product is at most the product of diagonal entries. The function
phi_t(x)=-log[(1+exp(-4tx))/2]/2 is increasing and concave with phi_t(0)=0.
After one deletion, remove the s-1 largest remaining magnitudes as well.
All outside squares are at most m/(s-1), and the concave chord lower
bound yields the stated decay epsilon q phi_t(1/q) per coordinate.
This argument is asymptotic with s>=2.

The flat profile obeys the anti-sparsity constraint whenever
epsilon<1-q, yet its permanent square root has decay at most log2/2
per coordinate. Therefore relaxing to all profiles satisfying this
constraint and then paying 2^k cannot give the desired upper bound for
p>1/2. This is a limitation of that relaxation, not a claim that an exact
flat spectrum occurs on every good selector.

## 5. A narrow extension to matrix-dependent tilts

For UNCONDITIONED uniform selectors at EXACT dyadic p in the Walsh
setting, combine the uniform spike bound with the typical-profile bound.
This excludes even arbitrary sequences of tilts t_m when
0<p<=p_*=0.9225232669048273 from satisfying the strict negative-exponent
certificate.

Indeed t_m tending to infinity makes the tilt penalty
t_m(1-sqrt(p)) dominate the uniform -h(p) spike lower bound. If t_m
tends to zero, K_t(a,b)>=exp(-t(a^2+b^2)) and sum a_j^2=m imply
L_t>=exp(-tm), so the row-count pressure tends to at least p log2.
For a positive finite subsequential tilt t, monotonicity of Z_T(t) and
the fixed-tilt lower bound at t+delta give, after delta decreases to zero,
the same lower exponent g_p(t) as before. Its infimum is nonnegative
for p<=p_*. A desired additional positive gap only increases the tilt
penalty. These three subsequence cases exhaust all t_m.

This extension uses the rare exact selector and thus does NOT carry over
automatically to conditioned selectors or nonexact dyadic proportions.

## 6. Walsh divisibility and translated RM information sets

The complete director file
`continued_director_walsh_divisibility_selector_count_2026_09_06.md`
has now been independently reconstructed. Its degree bound is precisely
d-r, NOT d-r+1. The ternary input xi reduces modulo two directly to
the support indicator. The algebraic coefficient at a marked subset J
is the integer subcube sum modulo two, and Fourier inversion gives

```math
\sum_{x:\ x_{J^c}=0}\xi(x)
=2^{|J|-d}\sum_{y:\ y_J=0}(H\xi)(y).
```

If every Walsh coefficient is divisible by 2^r, the right side is even
as soon as |J|>d-r. There is no division by two from converting a sign
function to a Boolean function; importing that unrelated convention
would introduce the off-by-one error the director correctly avoids.

The weight-at-most-s Hamming ball is an information set for RM(s,d)
by triangular Boolean Mobius inversion, and all its translations remain
information sets. A codeword of weight at most delta m has some translated
information set with at most floor(delta K) ones, where K is the code
dimension. Counting these patterns and translations gives exactly
m sum_(j<=delta K) binomial(K,j). Complementing a degree-s selector
indicator stays in RM(s,d), so the displayed selector probability follows.
For s=d/2+O(1), K/m tends to 1/2, giving the rate -h(p)/2.

A fixed-weight swap ball has size
sum_(a<=alpha m) binomial(k,a)binomial(m-k,a). Its entropy increases up
to alpha=p(1-p), since its derivative is
log[(p-alpha)(1-p-alpha)/alpha^2]. Thus the director's specified range
and neighborhood exponent are correct. At p=15/16 the largest strict
radius allowed by that displayed exponent is
alpha<0.0160686461175...; for alpha=0.015 the added rate is
0.1113493487<h(p)/2=0.1168958294.

Finally, flipping one nonzero input spin changes every Walsh coefficient
by plus or minus two. A spectrum divisible by four becomes identically
two modulo four, while its normalized mean squared spectral displacement
is only 4/k. The unchanged support is still rejected by the existence-
based filter, so this is a correct limitation on spectral divisibility
stability, not a counterexample to the selector count.
