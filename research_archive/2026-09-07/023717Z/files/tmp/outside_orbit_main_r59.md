# Wave 59 main audit: entropy-visible outside switching

## Status

The identities, KL lower bound, and marginal-visibility estimates below are
verified.  They give an exact global interface between a bare favourable
fibre and row-descending block switches.  The result is mixed: the parent
deficit is automatically affordable for an `O(H)` outside block, and a saved
fibre cannot hide a row much larger than the project cap in its average
complement.  But exact restriction preservation has an unavoidable linear
information cost per flipped coordinate, while the universal orbit contracts
row by only that same linear fraction.  Hence this mechanism alone does not
regularize a generic high row to the project cap.

Independent enumeration: `tmp/outside_orbit_main_r59.py`.

## 1. A restriction-preserving joint law

Fix an oriented cut `d=(tau,x)` and a nonempty bare fibre `F=F_t^d`.  Draw
`S` uniformly from `F`, put `T=S^c`, `|T|=r=n-m`, and draw a `k`-set
`U subset T` uniformly.  Output `d^U`, the physical-spin flip on `U`.

Since `U cap S` is empty, the induced local state is unchanged.  Therefore

```math
\widehat\ell(S,d^U)-\widehat\ell(S,d)
=-p_2\{\Delta(d^U)-\Delta(d)\}
\le p_2\Delta(d).
\tag{M59.1}
```

The last inequality uses only `Delta(d^U)>=0`.  Thus a ground column retains
the whole sampled incidence with no threshold loss; a near-ground column
with `Delta(d)=O(T_n)` retains it after `O(T_n)` enlargement.

If `h=-log U_m(F)`, the reference-free information of this joint law is at
most

```math
K(S;d^U)\le h+log\frac{\binom nk}{\binom{n-m}k}
=h+k\log\frac1{1-p}+O(k^2/n).
\tag{M59.2}
```

There is also an exact converse for every joint law supported on
`U cap S=empty`, with no fixed-size assumption:

```math
\boxed{
I(S;U)+D(P_S\Vert U_m)
=\mathbb E_U D(P_{S|U}\Vert U_m)
\ge\mathbb E_U\log
\frac{\binom nm}{\binom{n-|U|}m}.}
\tag{M59.3}
```

This is KL projection onto the family of selectors avoiding `U`.  On a
compact fixed-density window it implies `E|U|=O(H)` whenever the information
cost is `O(H)`.  Exact preservation of a local state therefore cannot
randomize a linear number of parent coordinates at the saved entropy scale.

## 2. Exact conditional row and deficit formulas

Gauge `A` by `x`.  Put

```math
C=diag(x)A^2diag(x),\qquad g=C1,\qquad R=1^TC1,
```

and let `C_T=1_T^TC1_T`, `g(T)=sum_(i in T)g_i`.  For the signed edge matrix
`B=tau diag(x)Adiag(x)`, put `rho=B1`, `E=1^TB1=q-Delta`,
`E_T=1_T^TB1_T`, and `rho(T)=sum_(i in T)rho_i`.  Direct fixed-size sign
averaging gives, with

```math
a_k=4k/r,\qquad d_k=4k(k-1)/(r(r-1)),
```

the exact identities

```math
\boxed{
\begin{aligned}
\mathbb E[R(d^U)\mid S]
&=R+a_k\{r(n-1)-g(T)\}
-d_k\{r(n-1)-C_T\},\\
\mathbb E[\Delta(d^U)\mid S]
&=\Delta+a_k\rho(T)-d_kE_T.
\end{aligned}}
\tag{M59.4}
```

For `k=1` these reduce to the known one-spin Euler formulas.  The point here
is their conditioning on an arbitrary favourable fibre.

## 3. Fibre entropy forces average visibility of high row

Let `mu_i=P(i in S | F)` and `p=m/n`.  Subadditivity of entropy and the
global strong concavity of binary entropy give

```math
\boxed{
\sum_i(\mu_i-p)^2
\le\frac12\{h+n h_b(p)-\log\binom nm\}
=O(h+\log n).}
\tag{M59.5}
```

Since `sum rho_i=E`, `||rho||^2=R`, `sum g_i=R`, and exact minimality gives

```math
||g||^2=x^TA^4x\le||A||_op^2R\le2qR,
```

Cauchy--Schwarz yields

```math
\boxed{
\begin{aligned}
\left|\mathbb E_F\rho(S^c)-(1-p)E\right|
&\le O(\sqrt{(h+\log n)R}),\\
\left|\mathbb E_F g(S^c)-(1-p)R\right|
&\le O(\sqrt{(h+\log n)qR}).
\end{aligned}}
\tag{M59.6}
```

Thus an `e^{-O(H)}` fibre cannot hide a row `R >> qH=R_*` from its average
complement.  The project scale is exactly where the entropy bias becomes as
large as the row itself.

Also `C_T<=2qr` and `|E_T|<=q`.  Averaging (M59.4), taking `h=O(H)` and
`k=Theta(H)`, gives

```math
\boxed{
\begin{aligned}
\mathbb E R(d^U)
&=R-\frac{4k}{n}\{R-n(n-1)\}
+O\left(\frac{k}{n}\sqrt{HqR}+\frac{k^2q}{n}\right),\\
\mathbb E\Delta(d^U)
&\le\Delta+O\left(\frac{kq}{n}
+\frac{k}{n}\sqrt{HR}+\frac{k^2q}{n^2}\right).
\end{aligned}}
\tag{M59.7}
```

The row equality uses `O` for a signed error.  Uniformly for
`R=O(nq)=O(n^(5/2))`, the second line is `Delta+o(T_n)`.  Parent deficit is
therefore not the obstruction.

At the row target, both errors in the first line have scale
`n^(2-2c)`, the same scale as the universal decrement `(H/n)R_*`; no sign is
left.  Far above the target the main decrement is visible, but it is only an
`H/n` fraction of the current row.  By (M59.3), applying a constant-fraction
universal randomization would cost `Theta(n)` information.  The exact open
successor is consequently a **high-leverage hub theorem**: use at most
`O(H)` fibre-visible coordinates whose block switch removes much more row
than the uniform `H/n` contraction, while retaining the signed relative-row
and deficit controls.  Neither entropy nor the PSD identities alone proves
such leverage.
