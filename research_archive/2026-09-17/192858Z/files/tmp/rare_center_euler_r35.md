# Wave 35 A: exact Euler calculus for a row-penalized rare center

Status: the identities and implications labeled **Verified** below were
derived directly and checked by exhaustive enumeration in
`tmp/rare_center_euler_r35_check.py`.  The finite decimals are **Numerical**.
The conclusion is a scoped negative one: the bare Euler-summing mechanism
does not supply the strict soft margin.  This does **not** falsify (10.967),
nor does it falsify the possibility that exact-signing minimality implies an
additional selector-agreement or mismatch-compression theorem.

## 1. Setup and exact block Euler inequalities

Fix the uniform selector slice `Omega_m`, favorable projective child fibers
`F_S`, and

```math
a_z(S)=\min_{[y]\in F_S}d_{\rm pr}(z_S,y),\qquad
q_\lambda(z)=\mathbb E_{S\sim U_m}e^{-\lambda a_z(S)},\qquad
F_\lambda(z)=-\log q_\lambda(z).
```

The center is a projective spin `z in {+-1}^n/{+-1}` and

```math
R(z)=R_2(z)=z^{\mathsf T}A^2z=\lVert Az\rVert_2^2.
```

For `B subseteq [n]`, write `z^B` for the spin obtained by flipping `B`, and
put

```math
\delta_B(S)=a_{z^B}(S)-a_z(S),\qquad
\pi_z(S)=\frac{U_m(S)e^{-\lambda a_z(S)}}{q_\lambda(z)}.
```

The projective triangle inequality gives
`|delta_B(S)| <= |B cap S|`.  Direct division of the two partition functions
gives the first exact identity

```math
\frac{q_\lambda(z^B)}{q_\lambda(z)}
=\mathbb E_{\pi_z}e^{-\lambda\delta_B(S)}. \tag{A.1}
```

Also, expanding `z^B=z-2P_Bz` gives

```math
\Delta_BR:=R(z^B)-R(z)
=-4\sum_{i\in B}z_i(A^2z)_i
+4\sum_{i,j\in B}z_i(A^2)_{ij}z_j. \tag{A.2}
```

Let `z=z_gamma` be a global minimizer of the smooth absolute-row objective

```math
\Psi_{\lambda,\gamma}(z)=F_\lambda(z)+\gamma R(z),\qquad \gamma\ge0.
```

Then (A.1) and exact cube optimality give, for **every** block `B`,

```math
\boxed{
\log\mathbb E_{\pi_z}e^{-\lambda\delta_B(S)}
\le \gamma\Delta_BR.} \tag{A.3}
```

This is the complete block Euler system; there is no differentiability or
asymptotic qualification.

For a one-bit block, `(A^2)_(ii)=n-1`, so

```math
\boxed{
\Delta_iR=4[(n-1)-z_i(A^2z)_i],\qquad
\sum_{i=1}^n\Delta_iR=4[n(n-1)-R(z)].} \tag{A.4}
```

Here `delta_i(S) in {-1,0,1}`.  If

```math
\beta_i=\pi_z\{\delta_i=-1\},\qquad
\chi_i=\pi_z\{\delta_i=+1\},
```

then the one-bit left side is exactly

```math
\log\{e^\lambda\beta_i+1-\beta_i-\chi_i+e^{-\lambda}\chi_i\}. \tag{A.5}
```

Choose deterministically, for each `S`, a closest oriented representative of
a fiber label, and let `B(S) subseteq S` be its mismatch set.  Then
`|B(S)|=a_z(S)`, and every `i in B(S)` has `delta_i(S)=-1`.  Consequently

```math
\mathbb E_{\pi_z}a_z(S)\le\sum_i\beta_i,
\qquad
\beta_i\le \exp\{-\lambda+\gamma\Delta_iR\}. \tag{A.6}
```

Jensen in (A.3), followed by
`sum_i delta_i(S) <= m-2a_z(S)`, also gives the exact summed inequality

```math
\boxed{
\lambda\{2\mathbb E_{\pi_z}a_z(S)-m\}
\le\sum_i\log\mathbb E_{\pi_z}e^{-\lambda\delta_i(S)}
\le4\gamma[n(n-1)-R(z)].} \tag{A.7}
```

The `m/2` term makes (A.7) much too weak at the desired scale.

## 2. The exact entropy wall

One direct use of **all** block inequalities still exposes an
exponential mismatch-pattern cost.  Group selectors according to the
deterministic nearest mismatch block:

```math
G_B=\{S:B(S)=B\}.
```

For `S in G_B`, flipping `B` lands exactly in its chosen favorable fiber, so
`a_(z^B)(S)=0`.  Thus

```math
\frac{q_\lambda(z^B)}{q_\lambda(z)}
\ge\frac{U_m(G_B)}{q_\lambda(z)}.
```

Combining this with (A.3) and summing the partition `{G_B}` proves

```math
\boxed{
q_\lambda(z)
\ge
\left(\sum_{B:G_B\ne\varnothing}e^{\gamma\Delta_BR}\right)^{-1}.} \tag{A.8}
```

At `gamma=0`, this is

```math
q_\lambda(z)\ge \frac1{L_z},\qquad
L_z=|\{B(S):S\in\Omega_m\}|
\le\sum_{j\le m/2}\binom nj. \tag{A.9}
```

For a fixed-density slice, the logarithm of the last quantity is
`Theta(n)`.  On the Wave 35 scales,

```math
TL_0=O(n^{3/4-c_0+\eta})=o(n),\qquad \eta<c_0,
```

so (A.9) misses the needed exponent.  Positive `gamma` does not generically
help (A.8), since mismatch flips that increase row square receive weights
larger than one.

There is an equivalent information-theoretic diagnosis.  The tilted selector
law satisfies the exact Gibbs identity

```math
\boxed{
F_\lambda(z)
=D(\pi_z\Vert U_m)+\lambda\mathbb E_{\pi_z}a_z(S).} \tag{A.10}
```

For an unpenalized global maximizer of `q_lambda`, (A.3) has right side zero.
Equations (A.5)--(A.6) then imply

```math
\mathbb E_{\pi_z}a_z(S)\le ne^{-\lambda}. \tag{A.11}
```

Thus the cube Euler conditions can make the tilted **energy** essentially
zero while leaving all of `F_lambda` in the selector KL term.  A center which
matches an exponentially rare selector subfamily exactly has precisely this
behavior.  The independent-label wall (10.930) shows that this phenomenon is
not removable for an arbitrary fiber system.  Such fibers need not arise
from one exact signing, so this is a mechanism wall, not an actual-minimizer
counterexample.

Equation (A.8) states the missing input sharply: an estimate

```math
\sum_{B:G_B\ne\varnothing}e^{\gamma\Delta_BR}
\le e^{O(TL_0)} \tag{A.12}
```

would prove the soft margin.  But (A.12) is a minimizer-specific
cross-selector mismatch-compression/agreement theorem; it does not follow by
summing the Euler inequalities themselves.

## 3. Row enforcement and the incompatible smooth-penalty scales

One-bit Lipschitzness gives

```math
|F_\lambda(z^i)-F_\lambda(z)|\le\lambda. \tag{A.13}
```

If a local minimizer of `F_lambda+gamma R` had

```math
R(z)>n(n-1)+\frac{\lambda n}{4\gamma},
```

then (A.4) would supply an `i` with
`Delta_i R < -lambda/gamma`, and (A.13) would give a strict descent.  Hence

```math
\boxed{
R(z)\le n(n-1)+\frac{\lambda n}{4\gamma}.} \tag{A.14}
```

In particular, the generic bit-descent conversion to
`C_2={R<=2n(n-1)}` requires

```math
\gamma\ge\frac{\lambda}{4(n-1)}. \tag{A.15}
```

By contrast, an absolute row penalty which is affordable at project exponent
for a center with `R=Theta(n^2)` has natural scale

```math
\gamma=O(TL_0/n^2).
```

Using

```math
\lambda=\Theta(rL_0/D),\qquad
d=D/s,\qquad s\sim r/T,
```

the ratio between (A.15) and this target-compatible scale is

```math
\frac{\lambda/n}{TL_0/n^2}
=\Theta\left(\frac{\lambda n}{TL_0}\right)
=\Theta\left(\frac nd\right)\longrightarrow\infty. \tag{A.16}
```

Indeed
`d=Theta(T n^(1/2-2c_0)/(log n)^2)=o(n)`.  Therefore the standard smooth
penalty has no parameter which simultaneously (i) enforces `C_2` by generic
cube descent and (ii) keeps its row charge at the desired exponent.  This is
a failure of this conversion, not a theorem that every possible Pareto or
barrier argument must fail.

There are two exact but non-improving conversions worth recording.

1. Since `E_z R(z)=tr(A^2)=n(n-1)`, a center with
   `R<=n(n-1)` exists.  Also `F_lambda(z)<=lambda floor(m/2)` for every
   center.  Thus minimizing

   ```math
   F_\lambda(z)+K\mathbf1_{\{R(z)>2n(n-1)\}},
   \qquad K>\lambda\lfloor m/2\rfloor,
   ```

   produces a center in `C_2`.  But its Euler system is useful only for flips
   that remain in `C_2`; outgoing flips receive the enormous allowance `K`.
   This is just the hard constrained problem in penalty form.

2. A genuinely sufficient but stronger smooth lemma is available.  Put
   `gamma=Gamma TL_0/[n(n-1)]`.  If constants satisfy
   `A<min{Lambda,2Gamma}` and one proves

   ```math
   \min_z\{F_\lambda(z)+\gamma R(z)\}\le A TL_0, \tag{A.17}
   ```

   then its minimizer cannot have `R>2n(n-1)`, because the row term alone
   would exceed `2Gamma TL_0`.  It also has
   `q_lambda(z)>=e^{-A TL_0}`, so (10.968) yields (10.967).  However, (A.17)
   charges row square even inside `C_2`, is strictly stronger than (10.967),
   and the Euler identities above do not prove it.

## 4. Exhaustive `A_6/A_8/A_9` diagnostic

The checker enumerates all projective centers, all `(n-1)`-selectors, all
exact projective child-ground fibers, every one-bit flip, and every flip
block, at `lambda=5`.  It verifies (A.1)--(A.8) numerically to `2e-10` after
all distances and row squares have been computed exactly.

For the unpenalized soft-degree maximizer:

| signing | row range (all centers) | `C_2` cap | maximizing `(R,q)` | distance histogram |
|:--|:--|:--|:--|:--|
| `A_6`, `m=5` | `30` | `60` | `(30, 0.8333409000)` | `0^5,2^1` |
| `A_8`, `m=7` | `8..72` | `112` | `(64, 0.3750058279)` | `0^3,2^1,3^4` |
| `A_9`, `m=8` | `16..128` | `144` | `(128, 0.3355843943)` | `0^3,1^3,2^1,3^1,4^1` |

Thus the best finite centers are already in `C_2`.  Nevertheless, at the
generic enforcement coefficient

```math
\gamma=1.01\,\lambda/[4(n-1)],
```

the `A_9` objective chooses the unique row-`24` center with

```math
q=3.0334827152\cdot10^{-5},\qquad
F=10.4032140946,
```

and distance histogram `2^6,3^2,4^1`; it has no exact child-fiber hit.  The
unpenalized row-`128` center has `q=0.3355843943`.  Hence the generic smooth
coefficient can severely over-regularize soft degree even when the desired
row constraint is already vacuous.  This is finite evidence only, and these
exact-ground fibers are narrower than the asymptotic favorable fibers.

## 5. Verdict and viable successor

**Verified verdict.**  Exact one-bit and block optimality reduce the proposal
to (A.8).  Their automatic consequence certifies only an `e^{-O(n)}` center;
they control the tilted distance energy but leave a selector KL/mismatch-code
entropy of order `n`.  Since `TL_0=o(n)`, merely summing these inequalities
cannot establish the strict soft margin.  Smooth absolute-row enforcement by
the generic descent argument operates at an incompatible coefficient.

One useful surviving sufficient statement is (A.12).  An alternative package
is `D(pi_z||U_m)=O(TL_0)` **and**
`lambda E_(pi_z)a_z(S)=O(TL_0)` at a suitably low-row optimized center; KL
control alone is insufficient by (A.10).  A proof would have to use the fact
that all favorable fibers come
from principal restrictions of the **same exact minimizer**, for example via
a two-selector/list-agreement theorem.  That is a distinct structural input,
not another Euler summation.  The direct exceptional-center target (10.967)
therefore remains open and unfalsified.
