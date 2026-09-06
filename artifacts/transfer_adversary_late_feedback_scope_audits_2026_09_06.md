# Late feedback scope audits

2026-09-06. Independent checks of two companion results. Both pass at
their explicit scopes; neither decides the original convergence question.

## 1. Exact twin-Hadamard rich holes without a high-degree tail

The director's construction in
`transfer_director_rich_holes_without_high_degree_2026_09_06.md` was
checked both in its initial perturbative form and in its stronger exact
form. Let k>=2 be a power of two, H_k the symmetric Sylvester matrix,
K=J_2 tensor H_k, A=K-diag(K), B=A/sqrt(2k-1), Q=B^2.
Writing d_a=H_aa and t_a=S_(1,a)+S_(2,a), direct expansion gives

```
(A^2 S)_(s,a)=(2k-2)t_a+S_(s,a)-R_a,
R_a=sum_(b != a) H_ab(d_a+d_b)t_b.
```

The diagonal d_a=(-1)^(Hamming weight of a) is balanced. There are
k/2-1 indices b!=a with d_b=d_a, and each summand then has magnitude
at most four. Hence |R_a|<=2k-4. If t_a=0, the numerator has absolute
value at most 2k-3. If t_a=2z, then S_(s,a)=z and its z-oriented value
is at least 2k+1. Both inequalities are STRICT relative to the threshold
2k-1. Therefore the fixed response satisfies the exact cube identity

```
sign((QS)_i) 1{|(QS)_i|>1} = (S_i+S_mate(i))/2.
```

It has exactly one-half expected mark mass and one-half expected holes,
and NO original Fourier degree above one. The matrix is an actual hollow
signing with asymptotic normalized operator norm sqrt(2).

This genuinely disproves an implication from arbitrary rich-frame holes
to positive high-original-degree mass. It does not challenge that
implication for a bounded nonconstant function of Gaussian old variables.
Nor is the parent near-minimizing: tr(H_k)=0 and evaluating at (z,z)
gives `H_A((z,z))=4 H_(hollow H_k)(z)`. Thus its normalized cap is at
least sqrt(2) times the order-k minimum's normalized cap, separated from
the known sub-half upper by the banked lower bound.

## 2. Uniform fixed-complexity first-marked responses

Read all of
`transfer_fresh_uniform_first_marked_family_and_core_barrier_2026_09_06.md`.
The finite-dimensional residual-distance argument, explicit coefficient
bounds, noise-sensitivity tail and finite degree catalog pass.

For the important ACTUAL uniformity step, averaged marginal weak
convergence to a continuous Gaussian CDF gives uniform CDF convergence
by a finite threshold grid. If the CDF error is zeta_n, every averaged
atom has mass at most 2*zeta_n. Moving threshold strips for a fixed
finite grid thus have precisely the required vanishing actual mass.
After taking one fixed label pattern and convergent thresholds, the
actual ternary vectors converge in averaged L2. No assertion about an
arbitrary Gaussian-equivalence class of responses is made.

For actual ternary F,G and their hole masks, the exact estimate is

```
|j(B,F)-j(B,G)| <= 2L sqrt(average E|F-G|^2).
```

One term uses |H_F-H_G|<=|F-G| and the old transported L2 bound; the
other uses bounded-op transport of F-G. This transfers the fixed-limit
response theorem to the moving compact family without uniform diagram
approximations of growing complexity.

## 3. The conservative certificate cannot pay its own core loss

The same companion artifact's constants also pass independently.
With its endpoint parameters,

```
d0 <= a^2/(64 M sqrt(R) sqrt(2pi)),
R=128 M L^4/a^2,
Delta=d0^2/(2L)
 <= a^6/(2^21 pi M^3 L^5)
 <= 1/(2^24 pi L^5),
```

using M>=2 and the source's a<=1. If p=1-epsilon and
`L>=D/(epsilon sqrt(p))`, with D=4 K_G C>=1, then

```
p^(3/2) Delta <= epsilon^5 p^4/(2^24 pi D^5)
              < j epsilon <= j(1-p^(3/2))
```

for every epsilon in (0,1) and j>=1/4. This is an analytic all-epsilon
barrier for THIS explicitly certified positive number, not an upper
bound on the actual possible spin-ascent improvement.

The finite core normalization in the companion is correctly
`m sqrt(m-1)/N^(3/2)`. Its use of a strict operator-cap slack, rather
than replacing the exact denominator by its limit prematurely, is
also necessary and correctly stated.
