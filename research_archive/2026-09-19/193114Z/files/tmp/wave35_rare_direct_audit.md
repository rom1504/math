# Independent Wave 35 audit: rare-center Euler and direct-cut hinge memos

## Overall result

The block increments, one-bit row sums, Gibbs signs, hinge sandwich, and
strict soft-to-hard margin were independently re-derived and are correct.
Two strategic conclusions require substantive weakening, and one checker
description overstates its coverage.

## Rare-center memo

### Accepted identities

- From `F=-log q`, global minimality of `F+gamma R` gives
  `log E_pi exp(-lambda delta_B)<=gamma Delta_B R`; the sign in (A.3) is
  correct.
- Expanding `z^B=z-2P_Bz` gives (A.2), and summing the one-bit increments
  gives `4[n(n-1)-R]`.
- A chosen nearest mismatch block has size `a_z(S)` and each of its bits
  reduces projective distance by one.  Hence (A.6) and the `m-2a` inequality
  in (A.7) are correct.
- For `G_B`, the inequalities are
  `U(G_B)/q(z)<=q(z^B)/q(z)<=exp(gamma Delta_B R)`.  Summing gives (A.8)
  with the displayed sign.
- The Gibbs identity is `F=D(pi||U)+lambda E_pi a`; all signs in (A.10) are
  correct.
- The row descent gives
  `R<=n(n-1)+lambda n/(4gamma)` for `gamma>0`, hence (A.15).  The scale
  ratio in (A.16) is `lambda n/(TL0)=Theta(n/d)`, with
  `d=Theta(T n^(1/2-2c0)/(log n)^2)=o(n)`.
- The strict conditions in (A.17) are correct: `A<2Gamma` forces the row
  cap and `A<Lambda` supplies the soft-hard margin.

### Required correction

The verdict calls (A.12) “equivalent” to
`D(pi_z||U_m)=O(TL0)`.  This is false.  The exact identity is

```math
F_lambda(z)=D(pi_z||U_m)+lambda E_pi a_z.
```

KL control alone does not control `q=e^{-F}`.  The correct alternatives are:

- (A.12) is a sufficient mismatch-compression certificate; or
- prove both `D(pi_z||U_m)=O(TL0)` and
  `lambda E_pi a_z=O(TL0)` at a low-row center.

For the unpenalized soft maximizer, (A.11) can make the energy term small,
but that optimizer need not satisfy the row cap.  Penalized low-row
optimization does not automatically retain (A.11).

### Wording corrections

- The mismatch count automatically **certifies** only `q>=exp(-O(n))`; it
  does not prove that the best center is `exp(-Theta(n))` in every signing.
- “Strongest elementary use of all block inequalities” is not a proved
  optimality statement and should be softened to “a direct joint use”.
- State `gamma>0` explicitly in (A.14).

## Direct-cut memo

### Accepted identities

- Applying the finite Gibbs principle for each cut and commuting two finite
  maxima gives (M35.1); this is existential, not an adversarial-law target.
- (M35.2) has the correct KL sign.
- Splitting at `h<=t+u` gives
  `K<=U{h<=t+u}+exp(-lambda u)`.  Thus `A<Lambda` in (M35.4)--(M35.5) is the
  necessary strict margin.
- With `t+u=O(n^(3/2-c'))`, `TL0=O(n^(3/4-c'))`, and
  `C=2n(n-1)`, the resulting cut meets the stronger row version of (10.795).
  The constants and witnesses must be quantified uniformly for every
  relevant `(A,n,m)`, while `w,d,t,u` may depend on that target pair.
- A physical block flip changes the child and parent crossing sums by
  `-4sigma W_(H,S)` and `-4sigma W_H`; substituting into `h` gives (M35.8).
- Expanding `x^H=x-2P_Hx` gives (M35.9), independent of orientation.
  The summed one-bit relations in (M35.10) and orientation increment `2X_d`
  have the correct factors.
- Optimality of `-log K+gamma R` gives (M35.11) with the displayed sign.
- The one-bit loss bound (M35.12) and generic enforcement coefficient
  `gamma=Theta(lambda)` in (M35.13) are correct.

### Required strategic correction

Section 4 asserts
`gamma R_2=Omega(n^(5/4))` and concludes that simultaneous generic row
enforcement and a target-scale objective bound are impossible.  No available
fact gives `R_2(d)=Omega(n^2)` for the optimizing cut.  Since
`R_2(d)=||Ax||^2`, it may be only `O(n)` or even vanish for a singular
signing.

The rigorous scoped conclusion is:

- `gamma=Theta(lambda)=Theta(n^(-3/4))` charges a cut with
  `R_2=Theta(n^2)` by `Theta(n^(5/4))`, far above `TL0`;
- a target-scale bound on the whole penalized objective would instead force
  the much stronger row estimate
  `R_2=O(TL0/lambda)=O(T n^(3/2-c0))`;
- generic bit descent enforces only `R_2<=2n(n-1)` and does not produce this
  stronger estimate.

This is a genuine scale mismatch and mechanism wall, but not an
impossibility theorem.  An anomalously low-row soft optimizer is not ruled
out.

### Checker-scope correction

The direct-cut script exhaustively constructs all cuts and selector losses,
but it does not exhaustively verify every block identity:

- `audit_block_increments` samples 300 RNG-selected `(cut,H)` pairs;
- `audit_free_energy` checks Euler inequalities only for the `n` one-vertex
  flips and orientation, not all blocks.

The identities remain algebraically verified, but §5 should call these
randomized block checks and exhaustive-state/one-bit checks, or the checker
should be extended to enumerate all blocks.

The reported `A_8` and `A_9` finite optimizer values were reproduced.
