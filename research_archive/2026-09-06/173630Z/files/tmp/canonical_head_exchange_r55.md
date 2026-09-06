# Wave 55: joint-completion concentration retires the far-tail implementation

## Scope and conclusion

- **Verified:** the full annealed arbitrary-cut incidence has a universal
  joint Hanson--Wright decomposition.  At depth `D T_n` its contribution is
  at most `2 exp(-c D H)`, uniformly in the selector.
- **Verified consequence:** any lower bound `Z_t >= exp(-C_0 H)` forces saved
  local-margin mass after choosing one sufficiently large fixed `D`; that
  local mass already proves the restriction recurrence.  Thus the uniform
  annealed completion-tail route is recurrence-circular as an independent
  implementation.  This does not say that `Z_t` itself is small.
- **Verified sharper K-specific consequence:** the K-qualified contribution
  at depth `D T_n` is at most `2 exp(-c D^2 H)`.  In particular, the proposed
  diverging-depth abundance (10.1277), and hence its scalar subevent
  (10.1285), is impossible: at `D=omega_n -> infinity` it is
  `exp(-omega(H))` for every exact minimizer.
- **Verified:** a fixed large depth does not rescue the K route.  Any
  `exp(-O(H))` K-qualified population can be truncated above at another fixed
  multiple of `T_n`, leaving saved mass in an `O(T_n)` local band, which again
  proves the recurrence directly.
- **Verified exact algebra:** canonical-head exchange has clean formulas, but
  selector-cap change is controlled only by an `O(sm)` Lipschitz term and the
  tracked state's local deficit.  More importantly, exchanging the head
  cannot evade the joint-completion theorem, which already bounds the entire
  far contribution before a head is selected.

Throughout, `A` is an exact order-`n` minimizer, `q=Q(A)=q_n`,
`p=m/n` stays in a compact subinterval of `(1/2,1)`, `0<c<1/4`, and

```math
H=\lceil n^{3/4-c}\rceil,
\qquad T_n=n^{3/2-c},
\qquad p_2=\frac{m(m-1)}{n(n-1)}.
```

We use the already verified exact-minimizer estimates

```math
q_n=O(n^{3/2}),
\qquad \lVert A\rVert_{\rm op}^2\le 2q_n.
```

## 1. A joint Hanson--Wright theorem for the complete incidence

Fix `S`, put `T=S^c`, and write

```math
A=\begin{pmatrix}P&B^{\mathsf T}\\ B&C\end{pmatrix}.
```

For a local projective spin `y`, an independent outside spin `w`, and an
orientation `sigma`, the centered completion increment in (10.1244) is

```math
Z=\sigma\{2y^{\mathsf T}B^{\mathsf T}w+w^{\mathsf T}Cw\}.
```

Let `x=(y,w)` and

```math
C_S=A-P_SAP_S
=\begin{pmatrix}0&B^{\mathsf T}\\B&C\end{pmatrix}.
```

The choice of a representative of projective `y` causes no distributional
problem: changing that representative also changes `w` to `-w`.  Hence
uniform projective `y` and uniform `w` map bijectively to uniform projective
`x`.  Therefore

```math
\boxed{Z=\sigma x^{\mathsf T}C_Sx}
\tag{R55.1}
```

under a uniform full oriented projective spin.  Moreover,

```math
\boxed{
\lVert C_S\rVert_F^2=2mk+k(k-1)<n^2,
\qquad
\lVert C_S\rVert_{\rm op}
\le 2\lVert A\rVert_{\rm op}
\le2\sqrt{2q_n}.}
\tag{R55.2}
```

The first identity just counts retained ordered off-diagonal entries.  The
operator bound follows from `C_S=A-P_SAP_S` and compression.  Since `C_S`
has zero diagonal, Hanson--Wright, followed if desired by the independent
orientation symmetry, gives universal constants `c_0,c_1>0` such that for
every fixed `D>=1` and all sufficiently large `n`,

```math
\begin{aligned}
\Pr_{y,w,\sigma}\{Z\le-DT_n\}
&\le2\exp\left[-c_0\min\left{
\frac{D^2T_n^2}{n^2},
\frac{DT_n}{\sqrt{q_n}}
\right}\right]\\
&\le 2e^{-c_1DH}.
\end{aligned}
\tag{R55.3}
```

Indeed the two exponents have orders `D^2 n^(1-2c)` and
`D n^(3/4-c)`; their ratio tends to infinity because `c<1/4`.

Now retain

```math
h=Q(A[S])-p^{3/2}q_n-t,
\qquad e=\sigma y^{\mathsf T}Py,
\qquad g=(1-p_2)e-h.
```

The exact normal form (10.1244) says that the arbitrary-cut incidence is
the event `Z<=g/p_2`.  Splitting according as `g>=-p_2DT_n` and applying
(R55.3) on the complement proves

```math
\boxed{
Z_t\le
\nu_m\{g\ge-p_2DT_n\}+2e^{-c_1DH}.}
\tag{R55.4}
```

This is uniform in `S`, so averaging selectors costs nothing.

Suppose the proposed annealed target `Z_t>=e^{-C_0H}` holds.  Choose one
fixed `D` with `c_1D>C_0+1`.  For all large `n`, (R55.4) gives

```math
\nu_m\{g\ge-p_2DT_n\}\ge\tfrac12e^{-C_0H}.
```

The inverse-Hanson--Wright local-band extraction (R54.5)/(10.1256) now gives

```math
q_m\le p^{3/2}q_n+O(T_n)
```

directly, without using any completion or row extraction.  Thus a saved
uniform arbitrary-cut incidence remains a correct sufficient condition, but
it is not an independent route around the local recurrence.

## 2. The sharper K-qualified upper tail

For fixed `S`, let

```math
u=By,\qquad G=\lVert u\rVert_2^2=y^{\mathsf T}Ry,
\qquad R=B^{\mathsf T}B.
```

Exactly,

```math
\boxed{
\mathbb E_yG=\operatorname{tr}R=mk,
\quad \lVert R\rVert_{\rm op}\le2q_n,
\quad \lVert R\rVert_F^2
\le\lVert R\rVert_{\rm op}\operatorname{tr}R
\le2q_nmk.}
\tag{R55.5}
```

The exact dual formula gives `K_H(u)<=sqrt(HG)`.  Since `b_H>=0`, a
K-qualified state at depth `r>=DT_n` must satisfy

```math
G\ge\frac{(r+b_H)^2}{H}\ge\frac{D^2T_n^2}{H}.
```

The latter threshold is `Theta(D^2 n^(9/4-c))`, while `mk<=n^2`.
Hanson--Wright for `y^T R y`, using (R55.5), therefore yields another
universal `c_2>0` such that

```math
\boxed{
\nu_m\{g<0, r\ge DT_n, K_H(u)\ge r+b_H\}
\le2e^{-c_2D^2H}.}
\tag{R55.6}
```

The intermediate two exponents are
`D^4 n^(1-2c)` and `D^2 n^(3/4-c)`; the second is the smaller one for all
large `n`.  Correlation between `r` and `G` is irrelevant because the event
is contained pointwise in this fixed upper tail.

Taking `D=omega_n -> infinity` in (R55.6) gives

```math
\nu_m\{r\ge\omega_nT_n, K_H(u)\ge r+b_H\}
\le e^{-\omega(H)}.
\tag{R55.7}
```

Thus (10.1277) is **falsified**, and (10.1285), being a subevent of it, is
also falsified.  In the orientation-relaxed union (10.1297), the K branch is
`e^{-omega(H)}`; any `e^{-C_0H}` union mass must therefore come from
`g^iota>=-C_BT_n`, whose orientation image is already saved local-band mass.

A fixed threshold `C T_n` does not avoid this conclusion.  If

```math
E_C=\{g<0, r\ge CT_n, K_H(u)\ge r+b_H\}
```

has mass at least `e^{-C_0H}`, choose a fixed `D>C` with
`c_2D^2>C_0+1`.  Equation (R55.6) leaves at least half that mass on
`C T_n<=r<D T_n`.  There `g=-p_2r>=-DT_n`, so the same local-band
extraction proves the recurrence directly.  This is stronger than merely
saying that diverging depth was too ambitious.

## 3. Exact canonical-head exchange algebra

Although the concentration theorem already controls every far completion,
the proposed exchange can be written exactly.  Let `I subset T` be the
canonical clipped head, `s=|I|<H`.  For `s>0`, put

```math
a_i=\operatorname{sgn}(u_i)\quad(\operatorname{sgn}0:=1),
\qquad z_I=\sigma a_I.
```

For each `s`-set `R subset S`, let `K=S\setminus R` and define the oriented
deleted energy

```math
D_R=\sigma\{y^{\mathsf T}A[S]y-y_K^{\mathsf T}A[K]y_K\}.
```

Use the specified removal rule: choose a minimizer of `D_R`, breaking ties
lexicographically.  Put `S'=K union I`, `T'=(T\setminus I) union R`, and
`y'=(y_K,z_I)`.  If `e'=sigma y'^T A[S']y'`, then direct expansion gives

```math
\boxed{
e'-e=-D_R+2\lVert u_I\rVert_1
-2a_I^{\mathsf T}A[I,R]y_R
+\sigma a_I^{\mathsf T}A[I]a_I.}
\tag{R55.8}
```

For a uniform `s`-set `R`, every internal edge is deleted with probability

```math
\alpha_s=\frac{s(2m-s-1)}{m(m-1)}.
```

Hence the removal rule satisfies `D_R<=alpha_s e`, and (R55.8) implies

```math
e'-e\ge2\lVert u_I\rVert_1-\alpha_se-(3s^2-s).
\tag{R55.9}
```

Let `q_S=Q(A[S])`, `h'=q_{S'}-p^(3/2)q_n-t`, and
`delta=q_S-e`.  With

```math
L_s=s(2m-s-1),
```

monotonicity through the common principal block `A[K]` gives the exact
deterministic controls

```math
\boxed{
h'-h=q_{S'}-q_S,
\quad |h'-h|\le L_s,
\quad h'-h\ge(e'-e)-\delta.}
\tag{R55.10}
```

The margin change is consequently

```math
\boxed{g'-g=(1-p_2)(e'-e)-(h'-h).}
\tag{R55.11}
```

The cross vector changes, for `U=T\setminus I`, by

```math
\boxed{
u'_U=u_U-A[U,R]y_R+A[U,I]z_I,
\qquad
u'_R=A[R,K]y_K+A[R,I]z_I.}
\tag{R55.12}
```

After zero-extending both profiles to `[n]`, the exact squared change is

```math
\lVert\widetilde u'-\widetilde u\rVert_2^2
=\lVert u_I\rVert_2^2+\lVert u'_R\rVert_2^2
+\lVert-A[U,R]y_R+A[U,I]z_I\rVert_2^2.
\tag{R55.13}
```

Thus the old canonical head is removed, but up to `2s` is added to every
retained outside coordinate and an unrestricted new `R`-profile appears;
there is no monotone K-functional comparison.

Finally choose arbitrary spins on `U` and retain the same full spin `x`
when only selector membership changes.  Its project row is exactly invariant:

```math
\boxed{R_2(x)=\lVert Ax\rVert_2^2\quad\hbox{before and after exchange}.}
\tag{R55.14}
```

If `Z` and `Z'` are the old and new completion increments, then

```math
Z'=Z-(e'-e),
\qquad
p_2Z'-g'=p_2Z-g+(h'-h)-(e'-e).
\tag{R55.15}
```

The last correction is exactly the change in local deficit
`(q_{S'}-e')-(q_S-e)`.  Equations (R55.10)--(R55.13) expose the scoped
exchange obstruction: for the uniform local states used by the incidence,
`delta` need not be small, and `L_s=Theta(sm)` can exceed `T_n` by a power
when `s` is near `H`.  A state-energy gain from the canonical head therefore
does not by itself charge selector excess or preserve the K-profile.  In any
case, (R55.3)--(R55.7) show that an exchange proof cannot rescue the
K-qualified far population as an independent uniform-completion mechanism;
it would have to establish the local recurrence after charging, or introduce
a genuinely different non-uniform mechanism.

## 4. Reproducibility

`tmp/canonical_head_exchange_r55_check.py` verifies (R55.1)--(R55.2),
(R55.5), and (R55.8)--(R55.15) on deterministic samples from the stored exact
`A_6,A_8,A_9,A_10` minimizers.  It also checks the exact K upper bound and
the set inclusion underlying (R55.6).  This finite checker audits algebra;
the asymptotic tail estimates themselves follow from Hanson--Wright and the
displayed norm bounds.
