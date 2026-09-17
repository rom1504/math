# Wave 27 Route 4: the signed deletion deck and deficit-square obstruction

## Status

The signed deletion-deck identity, parent/child change of measure, exact
extension-deficit square identity, child-profile upper bound, and the
Krawtchouk consequence of all star-competitor norm floors below are
**Verified**.  `tmp/cavity_moment_r27.py` audits every normalization on
`A_5,A_6,A_8,A_9` at `beta in {1/4,1/2,1}`.

The deck condition gives a new exact sufficient test for a selectable cavity
moment.  It is **not forced** by finite exact minimality: all four audited
minimizers fail it at `beta=1`, and `A_9` fails by a large signed margin.  The
beta-derivative identity points in the opposite direction and gives a sharp
falsification template, but no unbounded exact-minimizer family realizing the
template is known.  Thus there is no convergence proof and no asymptotic
counterexample.

## 1. Exact signed deletion deck

Let `B` have order `r`, put `h=r-1` and `rho=tanh(2 beta)`, and retain the
signed even-Eulerian polynomial of (10.814).  For an even Eulerian edge set
`F`, let

```math
n_j(F)=|\{i:\deg_F(i)=j\}|.
```

Write `P_i=P_{B[-i]}(rho)`.  The star interpolation from (10.815) has the
two simultaneous expansions

```math
\mathscr P_i(t)
=\sum_F b_F\rho^{|F|}t^{\deg_F(i)}
=P_i\left[1+\frac{\rho^2}{2}
  \{v_i(\beta;B)-h\}t^2+O(t^4)\right].
```

The first expansion keeps the original signs `b_F`; no absolute cycle
coefficient occurs.  Equating the quadratic coefficients gives

```math
\boxed{
P_i\{v_i(\beta;B)-h\}
=\frac2{\rho^2}
\sum_{F:\deg_F(i)=2}b_F\rho^{|F|}.}
\tag{R27.1}
```

Consequently, with the signed deletion decks

```math
N_0(\rho)=\sum_Fn_0(F)b_F\rho^{|F|}=\sum_iP_i>0,
\qquad
N_2(\rho)=\sum_Fn_2(F)b_F\rho^{|F|},
```

one has the exact weighted average

```math
\boxed{
\frac{\sum_iP_iv_i}{\sum_iP_i}
=h+\frac2{\rho^2}\frac{N_2(\rho)}{N_0(\rho)}.}
\tag{R27.2}
```

Thus, at a frozen root scale `alpha`, some vertex has
`v_i >= (9/16)alpha^2 r` whenever

```math
\boxed{
\frac{N_2(\rho)}{N_0(\rho)}
\ge\frac{\rho^2}{2}
\left\{\frac9{16}\alpha^2r-(r-1)\right\}.}
\tag{R27.3}
```

This is a path-selectable criterion if (R27.3) holds at every restriction
visited by the path.  Exact minimality of the root alone does not make later
restrictions exact minimizers, so a theorem only for isolated exact roots
would not settle the path compatibility.

There is an equivalent primal formula which checks all orientation and
projective factors.  Under the parent Gibbs law let

```math
h_i(\sigma,x)=\sigma x_i(Bx)_i,
\qquad A_i=e^{-\beta\kappa_{\beta,i}(B)}.
```

Deleting `i` and then extending its spin uniformly gives exactly

```math
\boxed{
A_i=\mathbb E_{\nu_{\beta,B}}e^{-2\beta h_i},
\qquad
A_iv_i=\mathbb E_{\nu_{\beta,B}}
   [h_i^2e^{-2\beta h_i}].}
\tag{R27.4}
```

In particular, (R27.2) is genuinely a positive mixture of the exact child
Gibbs moments even though its cycle numerator is signed.  The familiar
identity `sum_i h_i=e_B` does not lower-bound (R27.4): the exponential tilt
suppresses the large positive fields that carry the parent energy.

## 2. What all star-competitor norm floors imply

For a fixed vertex, let

```math
f_i(S)=\mathcal P_{B^S}(\rho),
\qquad S\subseteq\delta(i),
```

where `B^S` flips precisely the star edges in `S`.  The Boolean-noise
identity says that `P_i` is the uniform mean of `f_i`, while its summed
Fourier coefficient at star degree two is the numerator in (R27.1).  If
`k=|S|`, the degree-two Krawtchouk kernel is

```math
K_2(k)=\frac{(h-2k)^2-h}{2},
\qquad \min_kK_2(k)=-\lfloor h/2\rfloor.
```

Suppose `B` is an exact minimizer with `Q(B)=q`.  Every actual complete
signing `B^S` has norm at least `q`, hence

```math
f_i(S)\ge
L_q:=\frac{2\cosh(\beta q)}
 {2^r\cosh(2\beta)^{\binom r2}}.
```

Applying the minimum of `K_2` to the nonnegative function `f_i-L_q` gives

```math
\boxed{
v_i\ge h-
\frac{2\lfloor h/2\rfloor}{\rho^2}
\left(1-\frac{L_q}{P_i}\right).}
\tag{R27.5}
```

This is the strongest bound obtainable from the common pointwise floor and
the known mean by this exact degree-two extraction.  It is asymptotically
nonpositive at fixed `beta` because `L_q/P_i` can be exponentially small.

One can retain every competitor norm without losing cancellation.  Define

```math
L_i(S)=\frac{2\cosh\{\beta Q(B^S)\}}
 {2^r\cosh(2\beta)^{\binom r2}}.
```

The same argument applied pointwise to `f_i-L_i` yields

```math
\boxed{
v_i\ge h+\frac2{\rho^2P_i}
\left[-\lfloor h/2\rfloor P_i
+2^{-h}\sum_S L_i(S)
  \{K_2(|S|)+\lfloor h/2\rfloor\}\right].}
\tag{R27.6}
```

Thus even the full list of complete-signing star norms has a precise,
cancellation-preserving interface.  On the finite audit it remains far below
the required `9alpha^2r/16`; at `beta=1` it gives respectively
`-0.1601, 0.7687, 1.1266`, and at most `0.0167` on
`A_5,A_6,A_8,A_9`, against thresholds
`1.44,1.5625,3.515625,4`.  In particular, norm floors alone do not control
the signed deck at the necessary linear scale.

## 3. Exact beta-derivative obstruction

Let `C_i=B[-i]`, `d_i=Q(B)-Q(C_i)`, and under its exact Gibbs law put

```math
\delta_i=Q(C_i)-e_{C_i}(\omega),
\qquad
D_i(\gamma)=\sum_\omega e^{-\gamma\delta_i(\omega)},
\qquad \ell_i=\log D_i.
```

The two extensions of a child state have parent deficits

```math
\Delta_{i,\pm}=d_i+\delta_i\pm2L_i.
```

Both are nonnegative.  More precisely their product retains the missing
signed field exactly:

```math
\Delta_{i,+}\Delta_{i,-}=(d_i+\delta_i)^2-4L_i^2.
```

Averaging under the child Gibbs law and using
`ell_i'=-E delta_i`, `ell_i''=Var(delta_i)` proves

```math
\boxed{
4v_i
=(d_i-\ell_i')^2+\ell_i''
-\mathbb E_{\mu_{\beta,C_i}}
 [\Delta_{i,+}\Delta_{i,-}]
\le(d_i-\ell_i')^2+\ell_i''.}
\tag{R27.7}
```

This is the natural beta-derivative correction, and its direction is the
opposite of the desired moment floor.  It gives a sharp new falsification
template: if an unbounded family has, for every eligible vertex,

```math
d_i=o(\sqrt r),
\qquad -\ell_i'(\beta)=o(\sqrt r),
\qquad \ell_i''(\beta)=o(r),
```

then `max_i v_i=o(r)` and the cavity-moment criterion fails whenever the
root scales stay bounded below.  At `beta=1`, `A_9` has
`4v_i/[(d_i-ell_i')^2+ell_i'']` between `0.9305` and `0.9554`, so the upper
envelope is close to exact rather than merely formal.  No unbounded exact
minimizer family with these derivative properties is known.

The same extension cap gives a child-only strengthening of the deficit-ratio
falsifier.  Since `2|L_i|<=d_i+delta_i`,

```math
\boxed{
e^{\beta\kappa_{\beta,i}(B)}
\le
\frac{e^{\beta d_i}D_i(0)+e^{-\beta d_i}D_i(2\beta)}
     {2D_i(\beta)}.}
\tag{R27.8}
```

Equivalently, the parent deficit sum in (10.820) obeys

```math
D_\beta(B)\le D_i(0)+e^{-2\beta d_i}D_i(2\beta).
```

Therefore a realizable exact-minimizer family for which the logarithm of the
ratio on the right of (R27.8) is `o(sqrt(r))`, uniformly over the relevant
deletion cutset, would directly falsify fixed-temperature constant shortfall.
Unlike (10.820), this test needs only each child's three-temperature deficit
profile, not the parent/child deficit ratio.  The finite matrices satisfy the
inequality, but it is not asymptotically sharp enough there to supply such a
family.

## 4. Finite audit

For the root-normalized one-step threshold, the key deck values are:

| `beta` | signing | `N_2/N_0` | required RHS in (R27.3) | weighted `v` | outcome |
|---:|:---|---:|---:|---:|:---|
| `1/4` | `A_5` | `-0.150318` | `-0.273347` | `2.592217` | pass |
| `1/4` | `A_6` | `-0.270448` | `-0.367043` | `2.467151` | pass |
| `1/4` | `A_8` | `-0.334719` | `-0.372048` | `3.865230` | pass |
| `1/4` | `A_9` | `-0.391491` | `-0.427105` | `4.333538` | pass |
| `1/2` | `A_5` | `-0.851707` | `-0.742433` | `1.063209` | fail |
| `1/2` | `A_6` | `-1.093545` | `-0.996919` | `1.229323` | fail |
| `1/2` | `A_8` | `-1.393191` | `-1.010513` | `2.196105` | fail |
| `1/2` | `A_9` | `-1.868434` | `-1.160051` | `1.557408` | fail |
| `1` | `A_5` | `-1.790635` | `-1.189567` | `0.146476` | fail |
| `1` | `A_6` | `-1.856704` | `-1.597319` | `1.004292` | fail |
| `1` | `A_8` | `-2.700049` | `-1.619101` | `1.189376` | fail |
| `1` | `A_9` | `-3.624087` | `-1.858698` | `0.200807` | fail |

The checker independently enumerates all star perturbations, all oriented
projective parent and child states, and every actual competitor norm.  It
verifies (R27.1)--(R27.8), the signed weighted average, both Krawtchouk
bounds, and the displayed numerical data.

## Route conclusion

The cavity moment has an exact signed-deck selector: proving (R27.3) at each
visited restriction would give precisely the `9/16` premise without taking
absolute cycle coefficients.  The finite audit shows that exact minimality
does not force the signed ratio even at one exact root, while the strongest
direct degree-two use of all star-competitor norm floors remains far too
small.

The beta derivatives do not repair the sign.  They produce the nearly sharp
upper envelope (R27.7) and the child-profile falsifier (R27.8).  The next
genuine positive step would have to control the signed ratio `N_2/N_0` using
more than competitor maxima, or prove that some path avoids the low-deck
roots.  The clean negative alternative is an unbounded exact-minimizer
family satisfying the child derivative/profile criteria above; no such
realizable family is presently available.
