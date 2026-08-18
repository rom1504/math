# Adversarial audit: projected-coupling no-gain criterion

**Frozen source:**
`extremal_information/drafts/conference_projected_coupling_no_gain.md`

**SHA-256:**
`c2cf6feaabf3e3cdc6b1b4687ff71da51c5510a2956793ea84f1b8a3ed3b791c`

**Verdict:** **PASS, with one proof clarification and a useful
strengthening.**  Conditions PC.4--PC.7 imply PC.8 with all constants as
written.  In fact PC.8 has a direct one-sided proof that needs only
`Pr(G_r^c)=o(1)` and does not need condition 4 at all.  The source instead
passes through the stronger two-sided `L^1` assertion PC.11; for that
intermediate assertion, condition 4 must be read as a uniform crude bound
for all projected bridges occurring in the comparison, not only for the
full endpoint `B_r`.  This reading is automatic for exact-sign `B_r`, as
claimed.  Clarifying the referent of condition 4 would remove the only
ambiguity.

I found no circular pressure assumption.  PC.5 and PC.6 can be difficult to
verify, but they are geometric coupling hypotheses and do not encode the
desired lower pressure conclusion.

## 1. Operator normalization

Put `t=beta/sqrt(2r)`.  For a real bridge `C`, the block interaction is

```math
X(C)=t\begin{pmatrix}A_r&C\\C^T&\epsilon A_r\end{pmatrix}.
```

Since `||A_r||op=sqrt(r-1)`, PC.5 gives on `G_r`

```math
\|X(B_r^\circ)\|_{op},\ \|X(W_r^\circ)\|_{op}
\le {\beta\over\sqrt2}
 \left(\sqrt{1-1/r}+2+\delta\right).
```

The right side is eventually below a fixed `kappa<1/2` by PC.4.  Also
`||W_r(I-P_r)||op<=||W_r||op`, so the omitted projected-iid bound truly is
automatic.  The operator ball is convex, hence each comparison segment used
by pressure stability remains inside it.

The numeric hypothesis is slightly redundant in notation: PC.4 already
says the displayed coefficient is below `1/2`, after which `kappa` is chosen
between them.  It is nevertheless correct and is equivalent to the strict
campaign range after allowing a sufficiently small `delta>0`.

## 2. Nuclear constants

For a bridge increment `E`, the symmetric off-diagonal dilation has nuclear
norm `2||E||_*`.  The archived bound

```math
|d\log\overline Z_X[H]|\le {K_\kappa\over2}\|H\|_*
```

therefore gives

```math
K_\kappa t\|E\|_*
\le {K_\kappa\beta\over\sqrt2}\|E\|_F
```

when `rank(E)<=r`.  This is PC.9 with
`E=(B_r-W_r)(I-P_r)`.  There is no missing `sqrt(r)` or factor two.

For the removed part no Frobenius conversion is made, so the corresponding
bound is

```math
K_\kappa t\|CP_r\|_*
={K_\kappa\beta\over\sqrt{2r}}\|CP_r\|_*.
```

This checks both PC.10 and PC.12.  Dividing pressure by `r` shows why the
natural threshold in PC.7 is exactly `o(r^(3/2))`.

For PC.12 only the base `X(B_r^circ)` must be regular.  Indeed

```math
g(s)=\log\overline Z(X(B_r^\circ)+sY)
```

is globally convex, and the covariance theorem bounds `g'(0)` at the base.
The supporting-line inequality at `s=1` remains valid even when the full
endpoint is arbitrarily operator-irregular.  Thus the theorem does not
silently assume regularity of `B_r`.

## 3. Direct derivation of PC.8

The cleanest audit is to avoid the stronger intermediate PC.11.  On `G_r`,
apply PC.9, PC.10, and the one-sided PC.12 consecutively to obtain

```math
f_r(B_r)\ge f_r(W_r)-L_r,
```

where

```math
L_r={K_\kappa\beta\over\sqrt2}
       \|(B_r-W_r)(I-P_r)\|_F
 +{K_\kappa\beta\over\sqrt{2r}}
       (\|W_rP_r\|_*+\|B_rP_r\|_*).
```

Therefore

```math
1_{G_r}\left(h_\beta-{f_r(B_r)\over r}\right)_+
\le
\left|h_\beta-{f_r(W_r)\over r}\right|
+{L_r\over r}.
```

The iid conference theorem makes the first expectation tend to zero.
PC.6 and PC.7 give

```math
{\mathbb E L_r\over r}=o(1).
```

On `G_r^c`, every cosh pressure is nonnegative, so

```math
\left(h_\beta-{f_r(B_r)\over r}\right)_+\le h_\beta.
```

Consequently the exceptional contribution is at most
`h_beta Pr(G_r^c)=o(1)`.  This proves PC.8 directly.  It shows that the
theorem's conclusion remains valid if PC.5's probability assumption is
weakened from `o(r^(-1/2))` to `o(1)` and condition 4 is deleted.

This direct route also defeats a possible rare-event counterexample: even
an enormous real endpoint on `G_r^c` cannot enlarge a *lower* shortfall,
because pressure is nonnegative.

## 4. The stronger PC.11 and exceptional-event uniform integrability

As written, the proof claims the two-sided statement

```math
f_r(B_r^\circ)/r\longrightarrow h_\beta
\quad\text{in }L^1.
```

To justify it through PC.9, the bad-event contribution must control
`f(B_r^circ)` and `f(W_r^circ)`, and the PC.10 comparison also involves
`f(W_r)`.  Thus condition 4 should explicitly mean

```math
\sup_{G_r^c}
\max_{C\in\{B_r,B_r^\circ,W_r,W_r^\circ\}}f_r(C)
=O(r^{3/2}),
```

or, more minimally, that the corresponding bad-event expected pressure
sum is `o(r)`.

Under that reading, PC.5's
`Pr(G_r^c)=o(r^(-1/2))` yields

```math
O(r^{3/2})\Pr(G_r^c)=o(r),
```

which is exactly the required uniform integrability.

For an exact-sign `B_r`, the claimed automatic bound includes its
projection.  Orthogonal projection is a Frobenius contraction, so

```math
\|B_r^\circ\|_F\le\|B_r\|_F=r,
\qquad
\|W_r^\circ\|_F\le r.
```

For any such bridge `C`,

```math
|x^TCy|\le r\|C\|_{op}\le r\|C\|_F\le r^2,
```

while the two conference child energies also total `O(r^2)`.  Multiplication
by `beta/sqrt(2r)` gives the uniform `O_beta(r^(3/2))` pressure bound.  Hence
PC.11 is valid for the intended exact-sign application.

If condition 4 were read as bounding only `f(B_r)`, its use in the proof of
PC.11 would be incomplete for arbitrary real bridges.  This is a proof-text
ambiguity, not a counterexample to PC.8, because the direct argument in the
previous section proves the declared theorem without PC.11.

## 5. Rank interpretation

If `rank(P_r)=k_r`, then `rank(CP_r)<=k_r`, and Cauchy--Schwarz on singular
values gives

```math
\|CP_r\|_*\le\sqrt{k_r}\|CP_r\|_F.
```

For an iid Rademacher row `w`,

```math
\mathbb E\|wP_r\|_2^2=\operatorname {Tr}P_r=k_r.
```

Summing over rows and applying Jensen gives

```math
\mathbb E\|W_rP_r\|_*
\le\sqrt{k_r}\,\mathbb E\|W_rP_r\|_F
\le k_r\sqrt r,
```

which verifies PC.14.

If each sign row is conditioned on an event of probability at least `p_0`,
then, without requiring coordinate independence after conditioning,

```math
\mathbb E[\|RP_r\|_2^2\mid E]
\le p_0^{-1}\mathbb E\|wP_r\|_2^2
={k_r\over p_0}.
```

The same row-sum/Jensen argument yields
`E||B_rP_r||_*<=p_0^(-1/2)k_r sqrt(r)`.  Thus `k_r=o(r)` indeed implies
PC.7.  Cross-row independence is not needed for this particular estimate;
the marginal conditional-energy bound suffices.

The word “rank” is therefore literal: the removed component has matrix rank
at most `k_r`.  Rank alone does not imply harmlessness for arbitrary row
laws; the fixed-energy/constant-mass qualification is essential, and the
source states it in the paragraph deriving automatic PC.7.

## 6. Counterexample and circularity checks

1. **Huge irregular endpoint.**  Choosing `B_rP_r` with enormous operator
   norm does not violate the argument if its expected nuclear norm is
   `o(r^(3/2))`: convexity bounds its possible downward effect at the regular
   base.  It may raise pressure dramatically, consistent with the one-sided
   theorem.
2. **Rare bad coupling event.**  For a two-sided `L^1` statement, large real
   projected bridges can defeat uniform integrability unless condition 4 is
   read broadly.  For PC.8 itself, nonnegativity of pressure makes every bad
   event harmless once its probability tends to zero.
3. **Low rank with large amplitude.**  A rank-one component can have leading
   nuclear mass and is not covered merely because its rank is one.  PC.7
   explicitly excludes this, so the rank corollary is not making that false
   inference.
4. **Coupling circularity.**  PC.6 asks for a geometric coupling with small
   projected Frobenius edit cost, not for pressure closeness.  PC.5 asks for
   a sharp operator-regular projected bulk.  Neither contains `f_r` or the
   target `h_beta`; they may be strong, but they are falsifiable and are not
   equivalent by definition to PC.8.

No candidate satisfying the four stated conditions violates PC.8.

## 7. Recommended source clarification

No theorem repair is necessary.  For maximal precision, either:

- prove PC.8 by the direct one-sided chain above, weaken
  `Pr(G_r^c)=o(r^(-1/2))` to `o(1)`, and remove condition 4; or
- retain PC.11 and say explicitly that condition 4 bounds the pressures of
  `B_r`, `B_r^circ`, `W_r`, and `W_r^circ` on `G_r^c`.

The first version is both stronger and better aligned with the theorem's
one-sided conclusion.
