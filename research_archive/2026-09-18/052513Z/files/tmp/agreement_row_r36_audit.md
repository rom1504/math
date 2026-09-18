# Independent audit of `agreement_row_r36.md`

## Verdict

**PASS with one notation-only correction.**  In (A36.5), replace the comma
between the inclusion probability and `c_{ijk}(z)` by multiplication (for
example `\,c_{ijk}(z)`).  All mathematical identities, inequalities,
exponent comparisons, scoped examples, and the sufficient package are
correct as stated.

The checker

```text
.venv/bin/python tmp/agreement_row_r36_check.py
```

passes independently.

## 1. Uniform and anchored identities

After gauging by `z`, write `W=D_zAD_z`.  For every fixed `m`-set,

```math
L_S(z)=m(m-1)+2\sum_{\{i,j,k\}\subset S}c_{ijk}.
```

The diagonal inner-index terms give `m(m-1)`.  Each unordered triple gives
twice the displayed three-term coefficient.  Since

```math
R_2(z)=n(n-1)+2\sum_{i<j<k}c_{ijk},
```

uniform inclusion probability `p_3` immediately gives (A36.1), using
`m(m-1)=p_2n(n-1)`.

Conditionally on `v in S`, a triple containing `v` has probability
`alpha_2`, while a triple avoiding `v` has probability `alpha_3`.  Splitting
the triple sum at `v` gives exactly (A36.2), including the factor two and the
base term.

For (A36.3),

```math
\sum_{j<k\ne v}w_{vj}w_{vk}
=\frac{r_v^2-(n-1)}2,
```

and

```math
\sum_{j<k\ne v}w_{jk}(w_{vj}+w_{vk})
=(Wr)_v-(n-1).
```

Their sum is the asserted formula.  Moreover
`Wr=D_zA^2z`, so its second displayed form is correct.  The bound on `C_v`
follows because each triple coefficient has absolute value at most three.
At fixed positive selector density, `alpha_3` is bounded below and every
remaining term is `O(n^2)`, proving (A36.4).  The target exponent is strictly
above two, so inversion of (A36.2) is legitimate at that scale.

For an arbitrary slice law, the correct typeset form of (A36.5) is

```math
\mathbb E_PL_S(z)
=m(m-1)+2\sum_{i<j<k}P\{i,j,k\in S\}\,c_{ijk}(z).
```

The memo currently has a comma in place of the final multiplication; this
does not affect the argument.

## 2. Hamming transfer and exponents

Principal compression and `||z_S-y^S||_2=2sqrt(e_S)` give

```math
\big|\,||A[S]z_S||_2-||A[S]y^S||_2\,\big|
\le2||A||_{op}\sqrt{e_S},
```

which is (A36.7).  Applying scalar `L^2(P)` Minkowski and then
`E_Pe_S<=C_(v_*)` gives (A36.8); no Jensen direction is reversed.

With `T<=n^eta` and `c'=c_0-eta`, the two terms in (A36.9) have powers

```text
3/2
and at most 2-2c_0+eta,
```

whereas the target power is

```text
9/4-c'=9/4-c_0+eta.
```

The margins are respectively at least `1/2` and exactly `1/4+c_0` (up to
the harmless logarithmic gain).  Thus the asserted little-`o` comparison is
uniform.  The full-word perturbation bound (A36.10) is the same argument and
has the correct constant.

## 3. Ground stability and the sufficient package

Flipping `F` changes the gauged oriented energy by `-4a_S(F)`.  Comparing
with the child absolute cap gives `a_S(F)>=-delta_S/4`; singleton cuts give
the row-field form in (A36.11).  At zero deficit, every row field lies in
`[0,m-1]`, sums to `Q_S`, and therefore

```math
\sum_i(r_i^S)^2\le(m-1)\sum_i r_i^S=(m-1)Q_S.
```

Extending a child ground by uniform independent outside spins proves
`Q_S<=q_n`, so (A36.12)--(A36.13) and their exponent gap are correct.

The package (A36.14) is sufficient exactly as claimed:

1. its first line, conflict `O(d)`, and (A36.8)--(A36.9) give
   `E_P L_S(z)=O(n^(9/4-c'))`;
2. its second line transfers this to the uniform anchored slice;
3. (A36.2)--(A36.4), with the `O(n^2)` correction absorbed, give
   `R_2(z)=O(n^(9/4-c'))`.

The combined anchored row-transfer lemma is a valid alternative sufficient
statement.  It is explicitly only the weaker arbitrary-cut row clause, not
the stronger `2n(n-1)` exceptional-center cap.

## 4. Scoped examples

The conference-core example is correct.  The `h` hub rows alone give
`Omega(hn^2)`, while selectors avoiding the hubs have local row square at
most `(ell-1)m`.  Conditional hub avoidance has

```math
\log\beta^{-1}
=h\log(1/(1-rho))+O(h^2/n),
```

so its entropy cost is `Theta(n^(1/4))`.  The full anchored mean is
`Omega(hn^2)` by (A36.2), since the anchor correction is only `O(n^2)`.
The cap `Q(A)=O(n^(3/2))` follows from the conference spectral norm plus the
`O(hn)` hub contribution.  The memo correctly disclaims both favorability
and exact minimality.

For the positive-hub join, the triangle inequality gives

```math
Q(A)\le Q(B)+2hell+h(h-1),
```

and the all-one word attains equality, proving (A36.18).  For
`S=H union R`, the positive hub part contributes its full edge cap and

```math
Q(A[S])-1^TA[S]1
\le Q(B[R])-1^TB[R]1
\le2Q(B[R])\le2Q(B),
```

so (A36.20) is correct.  At fixed selector density the favorable coefficient
is a positive constant; choosing fixed `K` large enough makes its own-parent
allowance dominate `2Q(B)`.  General norm inequality (10.67) gives
`||A||_op=O_K(n^(3/4))`, and (A36.10) then preserves the
`Omega_K(n^(5/2))` row square throughout every `o(n)` Hamming ball.

The example is correctly scoped: replacing the true project value `q_n` by
the larger cap `Q(A)` is not valid for an actual-minimizer argument, so it is
only a mechanism wall.

