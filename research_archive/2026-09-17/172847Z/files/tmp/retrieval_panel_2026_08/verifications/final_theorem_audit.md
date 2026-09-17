# Independent final theorem audit

## Verdict

Both artifacts are mathematically verified as stated.  I found no
normalization, sign, quantifier, finite-field coverage, continuity-constant,
or closure error.  There are two harmless exposition omissions noted below;
neither requires a theorem correction.

## 1. Square-field Paley terminal-drift obstruction

1. **Conference normalization.**  For odd prime-power `r`, `q=r^2` is
   `1 mod 4`, so the displayed Paley matrix is symmetric and satisfies
   `C^2=r^2 I` at order `n=r^2+1`.  Therefore
   `(1/2)|x^TCx| <= (r/2)||x||_2^2=rn/2`.

2. **Boolean eigenvectors.**  Write `F=K+Kt`.  For
   `x_infinity=1`, `x_(a+bt)=f(b)`, and `sum_K f=1`, the infinity row is
   `r`.  At a finite coordinate in fibre `b_0`, the within-fibre character
   sum is `r-1`, while every other affine `K`-line has sum `-1`.  Hence
   `(Cx)_(a+bt)=1+(r-1)f(b)-sum_(d!=b)f(d)=rf(b)`.

   Thus `Cx=rx`, attaining the spectral upper bound and proving
   `Q(C)=rn/2` in the one-copy convention.

3. **Every finite edge is covered.**  For `z!=w`, more than `r-1`
   squares `u` are available, so one can choose square `u` with
   `u(z-w) notin K`.  Multiplication by `u` is a Paley automorphism and puts
   the two endpoints in distinct fibres `b,d`.  A sign function with
   `(r+1)/2` positive values can impose either requested product
   `f(b)f(d)=-chi(z-w)` for every odd `r>=3`.  Pullback gives
   `C_ij x_i x_j=-1`.

4. **Infinity edges are covered.**  Prescribing the relevant fibre value
   `f(b)=-1` is compatible with `sum f=1` because `r>=3`; hence the same
   opposition identity holds for `{infinity,z}`.

5. **Edge-flip sign and exact cap.**  Flipping the unordered edge changes
   the one-copy energy by `H_(C^e)(x)-H_C(x)=-2C_ijx_ix_j=+2`
   at the opposing `r`-eigenvector.  Conversely a one-edge flip changes
   every one-copy energy by at most `2`.  Therefore, for every edge,
   `Q(C^e)=Q(C)+2`, not `Q(C)-2` or merely a lower bound.

6. **Coset direction.**  Since `R=(binom(n,2)-Q)/2`, every flip decreases
   `R` by one.  Thus `b(U_C)=0` and
   `z(U_C)=Q(C)/n^(3/2)=r/(2 sqrt(r^2+1))->1/2`.

7. **Drift implication.**  The Paley roots eventually lie in every `I`
   containing `[0.33,0.51]`.  Uniformity over every eligible root and
   `epsilon_n->0`, followed by continuity, gives `beta(1/2)=0`; the unique
   zero is therefore `c=1/2`.  The established frontier puts every deepest
   root in `I` eventually, and deepest roots have `b=0`.  Compactness plus
   uniqueness then gives `M_n/n^(3/2)->1/2` along all orders.  The theorem is
   a scalable obstruction to calling `L_drift` a weaker convergence lemma;
   it is not a falsification of `L_drift`.

## 2. Finite certificate script

Running `.venv/bin/python computations/audit_paley_edge_traps.py` succeeded
and reproduced
`computations/results/paley_conference_edge_trap_audit.json` exactly.
It verifies all `45` edges at `(n,Q)=(10,15)` and all `325` edges at
`(26,65)`, certifying flipped caps `17` and `67`.  The script's oriented
condition `sign*C_ij*x_i*x_j=-1` is correct for both top and bottom
extremizers.  Its claim of exact flipped caps also legitimately uses the
one-edge Lipschitz upper bound; it is not inferred from coverage alone.

## 3. Quantitative action continuity

1. **Finite-matrix factor.**  On the uniform `n`-point probability space,
   `T_A=A/sqrt(n)` gives `<f,T_Af>=f^TAf/n^(3/2)`.
   Hollowness makes this expression affine in each coordinate, so its
   absolute maximum on `[-1,1]^n` occurs at a cube vertex.  Hence
   `Phi(T_A)=2Q(A)/n^(3/2)` exactly.

2. **Closed-profile step.**  Every actual profile law has `|X|<=1` and
   `E Y^2<=C^2`.  Portmanteau/lower semicontinuity preserves these facts in
   the weak closure.  The common second-moment bound also gives uniform
   integrability of `XY`, so taking the profile closure neither changes the
   energy supremum nor invalidates passage to weak limits.

3. **Constants.**  Strassen coupling at LP distance `delta` gives a good
   event of probability at least `1-delta`.  For
   `g_R(x,y)=x theta_R(y)`, the good-event cost is `(R+1)delta`, the bad-event
   cost is `2R delta`, and the two tails cost `2C^2/R`.  Thus
   `|int xy dmu-int xy dnu| <= 2C^2/R+(3R+1)delta`.
   Setting `R=C/sqrt(delta)` gives exactly `5C sqrt(delta)+delta`.
   Hausdorff matching in both directions justifies taking absolute suprema.

4. **Action metric.**  The `k=1` term has weight `1/2`, so
   `delta<=2d_M`; monotonicity gives the stated
   `5C sqrt(2d_M)+2d_M` modulus.  This proves continuity on every common
   `2->2`-bounded action class.  It does not provide spectral regularization
   or all-large-order exact signed realizers, so the artifact's boundary is
   correctly stated.

## Harmless wording points

- In the proof, “take laws `mu,nu`” should be read as: for each law in either
  closed profile, choose a Hausdorff-matched law in the other profile (using
  `delta+o(1)` if needed).
- The displayed theorem states `0<delta`; at `delta=0`, equality follows
  directly from identical closed profiles plus the same uniform-integrability
  argument.  Thus equation (2) also covers `d_M=0`.
