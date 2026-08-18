# Second independent audit: rowwise antipodal halfcube pressure

**Frozen source:**
`extremal_information/drafts/conference_row_halfcube_pressure.md`

**SHA-256:**
`d398b66308c63cdc2c4b00850c91ca769aceb00b66b29b5aec6d9a7a5b172c78`

**Verdict: PASS.**  The exact halfcube count, tie handling, uniform
pushforward, switching identity, common operator event, two-sided
convex-Lipschitz concentration, uniform center, orbit union, `L^1`
conclusion, and stated scope are correct.  No source repair is required.

The lower-tail part overlaps the more general adaptive-gauge theorem CB.5.
The additional content of RH.1 is the exact nonlinear halfcube
specialization, the upper tail, and hence two-sided convergence in
probability and in `L^1`.

## 1. Exact count, ties, and pushforward

For a nontied row `R`, exactly one of `R,-R` satisfies

```math
u_i\langle R,v\rangle>0.
```

For a tied row, oddness gives `chi(-R)=-chi(R)`, so exactly one member of
the pair satisfies `u_i chi(R)=1`.  Sign rows have no fixed point under
antipodal negation.  Thus, including the even-`r` tie mass, `H_i` contains
one member from each of the `2^(r-1)` antipodal pairs.  Independence of the
row constraints gives

```math
|\mathcal H_r(u,v)|=(2^{r-1})^r=2^{r^2-r}.
```

For each input row `W_i`, there is a unique `s_i(W)` such that
`s_i(W)W_i` is the selected representative.  Conversely, if
`B in H_r(u,v)`, its preimages are exactly

```math
W=D_tB,\qquad t\in\{+-1\}^r.
```

The canonical selector returns `s(W)=t`, so every output has exactly
`2^r` preimages.  Therefore `D_(s(W))W` is exactly uniform on the halfcube;
there is no weighting from ties or row margins.

## 2. Switching identity and common regular event

In the pressure with bridge `D_sW`, substitute `x=D_sx'`.  Then

```math
H_A(D_sx')=H_{D_sAD_s}(x'),
\qquad
(D_sx')^TD_sWy=(x')^TWy,
```

which proves RH.12 with no sign or orientation change.  Equivalently,
conjugation by `diag(D_s,I)` identifies the two parent matrices and hence
their operator norms.

Because a symmetric conference signing satisfies
`||A_r||_(op)=sqrt(r-1)`, block triangle inequality gives, simultaneously
for every switch `s`,

```math
\begin{aligned}
\left\|{\beta\over\sqrt{2r}}
\begin{pmatrix}D_sA_rD_s&W\\W^T&\epsilon A_r\end{pmatrix}
\right\|_{op}
&\le {\beta\over\sqrt{2r}}
  (\sqrt{r-1}+\|W\|_{op})\\
&\le {\beta\over\sqrt2}
 (\sqrt{1-1/r}+2+\delta)<\kappa
\end{aligned}
```

on `||W||_(op)<=(2+delta)sqrt(r)`, for all large `r`.  The same single
norm event works for all `2^r` switches.  Its complement has probability
`e^(-c_0r)` by the rectangular Rademacher norm tail.

## 3. The Talagrand input is genuinely two-sided

Let `K_1` be the convex operator-regular set for the unswitched children,
and let `g_1` be the supporting-plane extension from the archived regular
theorem.  It is convex and has a Frobenius Lipschitz constant independent of
`r`.  Talagrand's convex-Lipschitz theorem on the Rademacher cube gives a
two-sided Gaussian tail about its median, and equivalently about its mean:

```math
\Pr\{|g_1-Eg_1|>z\}\le C e^{-c z^2}
```

with constants depending only on the fixed high-temperature parameters.
The fact that the archived negative-moment application displayed only the
lower exponential moment does not make the underlying concentration
one-sided.

For a fixed switch define

```math
K_s=\left\{W:
 \left\|{\beta\over\sqrt{2r}}
 \begin{pmatrix}D_sAD_s&W\\W^T&\epsilon A\end{pmatrix}
 \right\|_{op}\le\kappa\right\}.
```

Conjugation shows

```math
W\in K_s\quad\Longleftrightarrow\quad D_sW\in K_1.
```

Consequently one may choose the extension

```math
g_s(W)=g_1(D_sW).
```

It is convex, has the same Lipschitz constant, and agrees with the switched
pressure on `K_s`.  Since `D_sW` is uniform whenever `W` is uniform,

```math
E g_s(W)=E g_1(W)
```

exactly.  Thus the center `h_beta r+o(r)` and its `o(r)` error are uniform
over all `s`; this is stronger than obtaining a separate asymptotic for
each member of a growing family.

On the common norm event of Section 2, `W in K_s` for every `s`, so the
two-sided Talagrand tail at deviation `eta r/2`, together with the uniform
center, proves RH.16.  No unproved two-sided pressure theorem is being
imported beyond the standard two-sided convex-Lipschitz inequality.

## 4. Union constants and the adaptive selector

For one switch and fixed orientation, RH.16 costs at most
`2e^(-c_(beta,eta)r^2)`.  Union over exactly `2^r` switch vectors costs

```math
2^{r+1}e^{-c_{\beta,\eta}r^2}.
```

Adding the once-only common norm exception gives RH.17.  For all large
`r`, the union term is at most `e^(-c_1r^2)` after decreasing the positive
constant.  Because the event is simultaneous in `s`, it remains true for
the fully bridge-dependent choice `s=s(W)`.  Combining this with the exact
uniform pushforward from Section 1 proves RH.19 and convergence in
probability under the halfcube law.

## 5. The `L^1` step and the child target

For a sign parent on `2r` vertices,

```math
0\le f_{\epsilon,r}(B)\le C_\beta r^{3/2}
```

uniformly: `cosh>=1` gives the lower bound, while the number of signed
quadratic terms is `O(r^2)` and the inverse-temperature scale is
`O(r^(-1/2))`.  After division by `r`, the contribution of the exceptional
event in RH.19 is bounded by

```math
O_\beta(\sqrt r)
 (e^{-c_0r}+e^{-c_1r^2})=o(1).
```

On its complement the normalized absolute error is at most `eta`.
Letting `eta` tend to zero proves the claimed `L^1` convergence.  Finally,
for any deterministic error `e_r=o(r)`, the fixed gap

```math
h_\beta-\tau_\beta=\gamma(\beta)>0
```

puts `tau_beta r+e_r` below `(h_beta-gamma(beta)/2)r` eventually.  RH.20
follows from the lower half of RH.19.

## 6. Archive comparison and scope

1. The small-Hamming retraction theorem does not apply: a majority repair
   can be of order `r^(3/2)`.  RH.1 uses an exact orbit representation, not
   an edit estimate.
2. The switching-gauge quotient theorem supplies the algebraic conjugacy
   behind RH.12, but its optimized response statement does not imply a
   quenched pressure limit on this fixed nonlinear fibre.
3. The later general adaptive-gauge result CB.5 already proves the
   orientation-specific **lower-tail** bound for every row/column switching
   selector.  RH.1's lower-tail conclusion is therefore subsumed by CB.5.
   Its two-sided Talagrand use adds the upper tail and `L^1` convergence,
   while its exact selector construction identifies the declared halfcube
   law.
4. The orbit union loses only `O(r)` in an `Omega(r^2)` exponent, but the
   common operator-norm exception itself has probability only
   `e^(-Omega(r))`.  The theorem therefore does not exclude an
   `e^(-Theta(r))` low-pressure subfamily inside the halfcube's
   operator-irregular remainder.
5. A row-sign-invariant magnitude fibre is not a selector of antipodal
   representatives and is outside RH.1.  The final proposed surviving
   direction is consequently not ruled out by the proof.

The source's conclusion should thus be read literally: the whole exact
rowwise halfcube has the typical pressure rate and cannot itself be a
same-temperature favorable basin.  It does not establish a full
superexponential lower-pressure tail for that halfcube or for arbitrary
rowwise nonlinear constraints.
