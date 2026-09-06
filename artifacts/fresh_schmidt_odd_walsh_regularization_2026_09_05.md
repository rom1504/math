# Schmidt's theorem settles the odd-Walsh regularized seed

Date: 2026-09-05. This is an exact primary-theorem import and scope correction,
not a proof of convergence of the original minimum.

## 1. Imported theorem and normalization

Kai-Uwe Schmidt, *Asymptotically optimal Boolean functions*, Journal of
Combinatorial Theory, Series A 164 (2019), 50–59,
[DOI](https://doi.org/10.1016/j.jcta.2018.12.005),
[author manuscript](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Diskrete_Mathematik/Publications-schmidt/pw.pdf).
Theorem 1 and the Fourier normalization immediately preceding Proposition 3
give, for the Sylvester matrix `W_N`, `N=2^d`,

`mu_d := min_{f in {±1}^N} ||W_N f||_infinity / sqrt(N) -> 1`.

This limit includes **all odd dimensions**, not merely even dimensions or a
subsequence. The primary proof first constructs a sequence of odd dimensions,
then uses `1<=mu_(d+2)<=mu_d` to obtain the full odd-dimensional limit.
Theorem 1, Proposition 3, and their proofs were read in the author manuscript.

## 2. Exact bilinear and same-spin consequences

Take a minimizing `f`. Since `W_N^T W_N=N I`,

`||W_N f||_1 >= ||W_N f||_2^2 / ||W_N f||_infinity
               = N^(3/2)/mu_d`.

Consequently, writing `beta(K)=max_{x,y signs}|x^T K y|`,

`1/mu_d <= beta(W_N)/N^(3/2) <= 1`.

Thus the normalized **maximal** bilinear excess tends to one. This conclusion
uses both Parseval and the peak bound; it is not a claim about a prescribed
Boolean input, greedy endpoint, or attraction basin.

For the same-spin objective put `g=sign(W_N f)` and
`z=(f,g,g,-f)`. Direct block multiplication, with
`W_4=H_2 tensor H_2`, gives

`z^T (W_4 tensor W_N) z = 8 g^T W_N f = 8 ||W_N f||_1`.

It follows that

`max_z |z^T W_(4N) z| / (4N)^(3/2) >= 1/mu_d -> 1`.

The spectral upper bound is one. Shifting `d` by two preserves odd parity,
so the normalized maximal same-spin Rayleigh quotient tends to one in **every
odd dimension** as well.

## 3. Exact landing in the campaign regularization

Use `q(K)=0.5 max_x |x^T K x|` and the actual regular outer
`R4=J4-2I4`. The independently checked signed-congruence identity from
`fresh_limit_algebra_phase2_2026_09_05.md` identifies
`R4^tensor a tensor H2` with `W_(2^(2a+1))`, up to signed permutation and
global sign. Hence

`lim_a q(R4^tensor a tensor H2)/4^(3a/2) = sqrt(2)`.

Therefore the Walsh-only regularization and the larger regularization allowing
the order-144 generator both satisfy

`R4_regularized(H2) = R(H2) = sqrt(2)`.

The upper bound for the larger family is the elementary spectral bound.
This rules out every positive uniform phase-defect upper certificate for the
H2 seed, including a defect restricted only to the Walsh branch. It does not
identify `R(B)` for a general seed, and it does not prove that every Hadamard
matrix becomes asymptotically regularizable under fixed R4 amplification.

## 4. What the primary proof does and does not generalize

Schmidt's construction uses a subgroup of index `v=7^e` in the multiplicative
group of a binary finite field. On almost all cosets the sign function is
constant. Special index-two Gauss-sum evaluations and Davenport–Hasse lifting
permit all relevant normalized Gauss sums to approach one along odd extension
degrees. The resulting Fourier transform acts approximately by inversion on
the coset values. Spencer's discrepancy theorem fills the exceptional coset
with Fourier error `O(sqrt(log(v)/v))`. Proposition 3 gives the explicit bound
`mu_d <= 1+12 sqrt(log(2v)/v)` in a suitable odd dimension.

The algebraic character identity is specific to the Walsh/finite-field
transform. Replacing it by an arbitrary fixed Hadamard tensor factor is a
new realization problem; the imported theorem does not assert that extension.

## 5. Archive reconciliation

`walsh_bent_stability_literature_audit.md`, Section 5, already cites Schmidt
but discusses only stability and greedy-basin information. That scoped
observation remains correct. It does **not** justify treating the global
odd-Walsh maximum or the regularized H2 seed as unresolved: Sections 2–3 above
are an overlooked direct consequence of the cited theorem.

Likewise, the permanent-stagnation examples in the phase-2 algebra report
remain correct for those recursive starting states. They do not contradict
the existence of different asymptotically saturating states.
