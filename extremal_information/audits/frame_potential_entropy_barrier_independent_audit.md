# Independent audit: frame-potential entropy barrier

**Frozen source:** `extremal_information/drafts/frame_potential_entropy_barrier.md`

**SHA-256:**
`fc3733f33ca316d4e9a2736889fc730f9db9f368a7e5c3e9bc7f3dad61256095`

**Verdict:** **PASS.**  The probability, entropy, projection-splitting, and
scope claims checked below are correct as stated.  I found no required
repair.

## Claims checked

1. **Basic frame identities.**  For independent sign rows `R_i`,

   ```math
   \|BB^T\|_F^2=r^3+2\sum_{k=1}^{r-1}
   \sum_{j\le k}\langle R_{k+1},R_j\rangle^2,
   ```

   `E||BB^T||_F^2=2r^3-r^2`, and the deterministic floor is `r^3`.

2. **Lemma FE.2.**  With `G=sum_{j<=k}R_jR_j^T` and spectral cutoff
   `H=G 1_[0,Lr](G)`, `L=8/(a eta)`, the discarded trace satisfies

   ```math
   \operatorname{Tr}(G-H)
   \le {\operatorname{Tr}G^2\over Lr}
   \le {\eta\over4}kr.
   ```

   Hence the asserted low-row event forces a downward deviation of at least
   `3 eta k r/4` for `X^T H X`.  The bounds
   `||H||_op<=Lr` and `||H||_F^2<=L k r^2` give an
   `exp(-c_(a,eta) r)` Hanson--Wright charge uniformly over every admissible
   deterministic prefix.  Both branches of the Hanson--Wright minimum have
   the claimed linear-in-`r` lower bound.

3. **Low-row counting in FE.1.**  On
   `F(B)<=(2-delta)r^3`, all prefix potentials are at most `2r^3`, and
   `sum Z_{k+1}<=(1-delta)r^3/2`.  If fewer than `rho r` late rows were
   low, the lower bound in (FE.19) would have leading coefficient

   ```math
   (1-\delta/8)(1-\delta/4)>1-\delta.
   ```

   Thus at least `ceil(rho r)` charged indices are necessary.

4. **Adaptive exposure and union bound.**  For each fixed charged index
   set `I`, rows may be exposed sequentially.  Conditional on every realized
   admissible prefix, FE.2 applies uniformly to the next independent row; if
   a prefix is inadmissible, the target event has already failed.  Iteration
   therefore gives `p_r^|I|` without an independence assumption between the
   low-row events.  The at-most-`2^r` union bound then yields
   `exp(-c_delta r^2)`.  There is no hidden conditioning gap here.

5. **Entropy corollary to FE.1.**  The deterministic floor and the threshold
   `(2-delta/2)r^3` force the low-frame event to have `q_r`-mass at least
   `delta/(2-delta)` in the zero-error case.  Binary data processing against
   its `U_r`-mass `exp(-Omega(r^2))` proves the stated quadratic relative-
   entropy lower bound.  Replacing `delta` by a fixed smaller value absorbs
   the `o(r^3)` error.

6. **Conference projection and FE.3.**  Writing
   `C=A/sqrt(r-1)`, left and right multiplication by `C` are commuting
   self-adjoint involutions.  Thus

   ```math
   P_\epsilon={1\over2}(I+\epsilon L_C R_C)
   ```

   is an orthogonal projection, `Tr P_epsilon=r^2/2`, and

   ```math
   \|AB+\epsilon BA\|_F^2
   =4(r-1)\,\operatorname{vec}(B)^T
     P_\epsilon\operatorname{vec}(B).
   ```

   The rank follows from `Tr C=0`.  The threshold in (FE.32) is, for all
   sufficiently large `r`, at least `delta r^2/16` below the projection's
   mean.  Hanson--Wright therefore gives `exp(-Omega_delta(r^2))`.
   Splitting `J_epsilon<=(4-delta)r^3` between the frame and intertwiner
   deficits is valid, and the entropy consequence follows from the floor
   `J_epsilon>=r^3`; the displayed mass `delta/(6-delta)` is correct.

7. **Pressure scope.**  The note does not infer a quartic or frame deficit
   from bare low pressure.  It invokes only the separately proved
   small-`beta`, FMW-power-regular implication
   `J_epsilon/r^3<=1+O(beta^2)`, which for sufficiently small fixed `beta`
   lies inside a fixed FE.1 deficit event.  The warning that higher cumulants
   may compensate at one fixed temperature, and that power-irregular
   target-reaching bridges remain uncontrolled, is logically necessary and
   correctly retained.

## Corrections

None required.
