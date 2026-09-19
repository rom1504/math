# Wave 49B independent cross-audit

## Verdict

The annealed-incidence theorem is correct.  Its column extraction can be
sharpened from `u_t(d)>=exp(-2L_t)` to `u_t(d)>=Z_t/2`, with the same factor-two
row cutoff.  There is no projective-chart normalization loss.

The cylinder wall `Y_S<=t+K(K-1)` is also correct provided `K` is defined as
the **total Boolean-chart codimension**: one orientation bit plus the number
of fixed independent projective spin bits.  Equivalently, `K` is the size of
the anchored physical vertex set after fixing the gauge spin.  If `r` instead
means the number of additional fixed spin bits excluding the gauge anchor,
the same bound reads `Y_S<=t+r(r+1)` and the cylinder mass is
`2^{-(r+1)}`.  This is a notation/quantifier clarification, not a missing
factor in the intended `K(K-1)` statement.

## 1. KL and row audit

For `F={(S,d):widehat ell(S,d)<=t}`, `Z=(U_m x Pi_n)(F)>0`, and the
conditioned law `P`, the density is exactly `1_F/Z`.  Hence

```math
D(P\Vert U_m\otimes\Pi_n)=-\log Z=:L.
```

The chain rule gives

```math
L=D(P_D\Vert\Pi_n)+D(P_S\Vert U_m)+I_P(S;D),
```

so in particular `D(P_D||Pi_n)<=L`.

For `C=A^2-(n-1)I`, `tr C=0` and

```math
R_2(d)=n(n-1)+x^{\mathsf T}Cx.
```

Although `Pi_n` is expressed using a projective chart and an orientation,
the `x`-quadratic is even.  Data processing to the projective `x` marginal,
followed by the symmetric lift to iid Rademacher `x`, preserves or decreases
KL.  Entropy duality and the centered Hanson--Wright mgf therefore give

```math
E_PR_2(D)\le n(n-1)+O(\lVert C\rVert_F\sqrt L+
\lVert C\rVert_{op}L).
```

The exact trace calculation is

```math
\lVert C\rVert_F^2=\operatorname{tr}A^4-n(n-1)^2.
```

Using `||A||_op^2<=2q_n`, `tr A^2=n(n-1)`, and
`q_n=O(n^(3/2))` gives

```math
\lVert C\rVert_F=O(n^{7/4}),\qquad
\lVert C\rVert_{op}=O(n^{3/2}),
```

and thus

```math
\overline R:=E_PR_2(D)
\le n(n-1)+O(n^{7/4}\sqrt L+n^{3/2}L).
```

There is no lost factor in these exponents.  At
`L=O(n^(3/4-c))`, the two errors are
`O(n^(17/8-c/2))` and `O(n^(9/4-c))`; the former and the `n^2` baseline are
bounded by the latter precisely when `c<=1/4`.

## 2. Sharper simultaneous extraction

Write

```math
\alpha(d)=U_m\{S:(S,d)\in F\},\qquad
\mu(d)=P_D(d)=\frac{\Pi_n(d)\alpha(d)}Z.
```

Markov under `mu` gives

```math
\mu\{d:R_2(d)\le2\overline R\}\ge\frac12.
```

For this row-good set `G`, therefore

```math
\sum_{d\in G}\Pi_n(d)\alpha(d)=Z\mu(G)\ge Z/2.
```

Since `Pi_n(G)<=1`, some `d in G` satisfies

```math
\boxed{R_2(d)\le2\overline R,\qquad \alpha(d)\ge Z/2.}
```

This is valid as written for the finite uniform oriented-projective law; no
chart multiplicity enters.  More generally, a row factor `a>1` gives mass
at least `Z(1-1/a)`.  Thus the separate logarithmic-cost averaging in
(R49.5)--(R49.6) is correct but unnecessary and loses a factor two in the
exponent.

## 3. Cylinder conditional-mean wall

Use the chart `x_0=1`.  Suppose a favorable cylinder fixes the orientation
and `K-1` independent tail-spin bits.  Let `J` be those vertices together
with the gauge anchor `0`, so `|J|=K`; the product mass is exactly `2^{-K}`.
On averaging the favorable inequality

```math
Y_S-\sigma x^{\mathsf T}W_Sx\le t
```

over all free chart bits, every quadratic monomial vanishes except those
with both endpoints in `J`.  Since `W_S` has zero diagonal and every
off-diagonal entry has absolute value at most one,

```math
\sigma E(x^{\mathsf T}W_Sx)
\le\sum_{i\ne j\in J}|(W_S)_{ij}|\le K(K-1).
```

Hence exactly

```math
\boxed{Y_S\le t+K(K-1).}
```

If the orientation is not fixed and `Y_S>t`, the two orientations give a
contradiction, so the only orientation-free case already has `Y_S<=t`.
Finally `Q(A[S])>=q_m` turns the wall into the desired restriction recurrence
when `K=O(n^(3/4-c))`.  The wall and its circularity claim are therefore
correct with the codimension convention made explicit.
