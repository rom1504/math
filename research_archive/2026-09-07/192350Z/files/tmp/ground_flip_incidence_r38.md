# Wave 38 main attack: structured-flip ground compatibility and its slack wall

## Outcome

The complement-flip cover has a stronger exact subcover obtained by choosing
an actual ground of each flipped signing.  In the gauge of such a state, every
switched cut has nonnegative weight.  Translated back to the parent, this
forces the selector to contain every vertex of positive parent signed degree,
and gives an all-subset cut corridor.  For a threshold certificate which is
not an exact flipped ground, the same statement survives with a sharp slack
equal to its flipped-ground deficit.

This structure does **not** by itself solve row compression.  Exact `A_8`
enumeration shows that replacing threshold incidences by actual flipped
grounds deletes every witness at the mean row cap for all 28 selectors, even
though 24 of those selectors have a threshold witness at or below that cap.
Thus ground compatibility is a strictly stronger and row-hostile
restriction at the mean row cap.  Any use of the clean zero-slack corridor
must separately prove low-row exact grounds at the larger project cap;
otherwise it must retain the flipped-cap surplus or another source of
near-ground slack.

All algebraic statements and finite incidence/miss counts below are exact and
checked by `tmp/ground_flip_incidence_r38_check.py`.  Fractional-cover values
reported by the checker use floating-point linear programming and are
numerical only.

## 1. Two simultaneous cut-stability systems

Fix an exact minimizer `A`, with `Q(A)=q`, an oriented parent state
`d=(sigma,x)`, and put

```math
s_{ij}=\sigma a_{ij}x_ix_j,
\qquad
r_i(d)=\sum_{j\ne i}s_{ij},
\qquad
E_d=\sum_i r_i(d)=q-\Delta_d.
```

For an `m`-selector `S`, let `B_S=A^(F_S)` be the signing obtained by
flipping every edge outside `E(S)`.  In the same gauge its signed edges are

```math
t^S_{ij}=
\begin{cases}
s_{ij},&i,j\in S,\\
-s_{ij},&\text{otherwise}.
\end{cases}
```

Its energy at `d` is the complement-flip payoff

```math
L_S(d)=\langle B_S,d\rangle=2c_S(d)-E_d.
```

Write

```math
M_S=Q(B_S),
\qquad
\Gamma_S(d)=M_S-L_S(d)\ge0.
```

For a vertex subset `U`, let `C_s(U)` and `C_t^S(U)` denote the unordered
signed weights across its cut, and put

```math
C_s^S(U)=
\sum_{\substack{i\in U\cap S\\j\in S\setminus U}}s_{ij}.
```

Flipping the spins in `U` changes the two oriented energies by four times
their respective cut weights.  The parent and flipped cap inequalities give
exactly

```math
\boxed{
C_s(U)\ge-\frac{\Delta_d}{4},
\qquad
C_t^S(U)=2C_s^S(U)-C_s(U)
\ge-\frac{\Gamma_S(d)}4
\quad\text{for every }U\subseteq[n].
}
\tag{R38.G1}
```

The second identity is just the fact that the edges internal to `S` retain
their sign while every other crossing edge reverses it.  It contains no
averaging or asymptotics.

If `U` is contained in the omitted set `S^c`, its internal-`S` term vanishes,
so (R38.G1) becomes the two-sided corridor

```math
\boxed{
-\frac{\Delta_d}{4}
\le C_s(U)\le
\frac{\Gamma_S(d)}4,
\qquad U\subseteq S^c.
}
\tag{R38.G2}
```

For a singleton omitted vertex this says

```math
r_i(d)\le\Gamma_S(d)/4,
\qquad i\notin S.
\tag{R38.G3}
```

## 2. Exact grounds and surplus-stratified threshold certificates

Let

```math
\mathcal G_d=\{S:d\text{ is a positive ground of }B_S\}.
```

For `S in G_d`, `Gamma_S(d)=0`.  Therefore, with

```math
P(d)=\{i:r_i(d)>0\},
```

equation (R38.G3) forces

```math
\boxed{
S\in\mathcal G_d\quad\Longrightarrow\quad P(d)\subseteq S,
\qquad
|\mathcal G_d|
\le\binom{n-|P(d)|}{m-|P(d)|}.
}
\tag{R38.G4}
```

The binomial is understood as zero when `|P(d)|>m`.  Choosing an actual
ground of every `B_S` gives a cover by these structured incidences because
`M_S>=q` by exact signing minimality.  Hence (R38.G4) is a genuine
selector-incidence constraint on one valid complement-flip subcover.

The useful version for the larger threshold incidence keeps the surplus.
Put

```math
\delta_S=M_S-q\ge0.
```

If `L_S(d)>=q`, then `Gamma_S(d)<=delta_S`.  Consequently, on the selector
stratum `delta_S<=gamma`,

```math
\boxed{
\{i:r_i(d)>\gamma/4\}\subseteq S.
}
\tag{R38.G5}
```

More strongly, every `U subseteq S^c` obeys (R38.G2) with upper endpoint
`gamma/4`.  This yields an exact surplus-or-support dichotomy: a threshold
witness may omit a large positive-degree vertex only by paying at least four
times that degree in flipped-cap surplus.

The direction of (R38.G4)--(R38.G5) is important.  It upper-bounds how many
selectors a highly constrained ground state can serve; it does not construct
a low-row state or upper-bound a fractional cover.  Moreover, at the
project-scale low-cap threshold the allowed surplus can exceed every
singleton degree, so the all-subset corridor, rather than only (R38.G5),
would have to do the work.

## 3. Exact `A_8` wall to imposing zero slack

The checker enumerates every oriented projective state and every selector for
`A_6,m=5`, `A_8,m=6`, and `A_9,m=7`.  It separately records

```math
I_d=\{S:L_S(d)\ge q\}
```

and `G_d`, verifies (R38.G1)--(R38.G4) on every exact flipped ground, and
computes the minimum parent row square in each selector fiber.

For `A_8,m=6`, every complement-flipped signing has

```math
M_S=24,
\qquad q_8=20.
```

Every one of the 28 selectors has minimum exact-ground row square `64`, so
all 28 lose every ground witness at the mean cap
`R_2<=8*7=56`.  In contrast, minimum threshold-incidence row squares have
the exact histogram

```text
40^8, 56^16, 64^4.
```

Thus only four selectors lose all threshold witnesses at the same cap.  The
unrestricted fractional cover values are numerically `20` for ground
incidence and `10` for threshold incidence.  The numerical LP values are not
needed for the exact wall: the 28-versus-4 miss counts already prove that
zero-slack ground restriction can discard the low-row certificates.

For context, `A_6,m=5` has `M_S=q_6=10` on all selectors and ground incidence
coincides with threshold incidence at minimum row `30`.  On `A_9,m=7`, all
flipped caps equal `28` versus `q_9=24`; the selectorwise minimum row
histogram happens to agree for the two incidences, but this does not repair
the exact `A_8` failure.

These are finite mechanism walls, not asymptotic counterexamples: the project
row cap is larger than the mean scale.  They do show that a proof cannot
silently choose exact flipped grounds and inherit the needed row property
merely in order to gain the clean positive-cut system.  A viable low-cap
continuation must either prove a new project-cap row theorem for exact flipped
grounds, use the slack corridor (R38.G1) and trade surplus against row
regularization quantitatively, or find a different structured family of
threshold witnesses.
