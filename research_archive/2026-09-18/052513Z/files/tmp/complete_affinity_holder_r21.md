# Wave 21: complete-energy Hölder identities for coordinate affinity

This note starts from (10.687) and retains the complete-signing identity
`Delta=q-2 sum_e s_e` at every state.  It produces exact aggregate and
product formulas, but every direct Hölder consequence has the wrong sign for
(10.624).  A separate pair-field argument gives a genuine universal lower
bound on the best cavity reward, of order only `r^(-1/2)` at fixed
temperature.  No complete-signing asymptotic counterexample to (10.617) is
claimed.

## 1. Exact statewise product law

Let `B` be an order-`r` complete signing, `q=Q(B)`, and use the lifted
oriented state space.  Put

```math
\Delta_\omega=q-e_B(\omega),
\qquad
w_\omega=e^{-\beta\Delta_\omega},
\qquad
D_\gamma(B)=\sum_\omega e^{-\gamma\Delta_\omega}.
```

For coordinate `i`, retain

```math
h_i(\omega)=\sigma x_i(Bx)_i,
\qquad
A_i=\operatorname{BC}_i=e^{-\beta\kappa_{\beta,i}}.
```

Equation (10.687) is

```math
A_i
=\frac1{D_\beta(B)}
\sum_\omega w_\omega e^{-2\beta h_i(\omega)}.
```

Completeness and the single-signing energy law enter through the exact
identity

```math
\boxed{
\sum_{i=1}^r h_i(\omega)
=e_B(\omega)=q-\Delta_\omega.
}
```

Thus, state by state,

```math
\boxed{
\prod_{i=1}^r e^{-2\beta h_i(\omega)}
=e^{-2\beta(q-\Delta_\omega)},
}
```

and for the directed geometric masses from (10.687),

```math
\boxed{
\prod_{i=1}^r
e^{-\beta[\Delta_\omega+2h_i(\omega)]}
=e^{-2\beta q-\beta(r-2)\Delta_\omega}.
}
```

These are the strongest products available before summing over states.

## 2. Hölder gives a lower affinity product

Apply generalized Hölder to the parent Gibbs law and the positive functions
`exp(-2 beta h_i)`.  The preceding product law gives

```math
\boxed{
\left(\prod_{i=1}^r A_i\right)^{1/r}
\ge
R_{\beta,r}(B)
:=
e^{-2\beta q/r}
\frac{D_{\beta(1-2/r)}(B)}{D_\beta(B)}.
}
```

Pointwise AM--GM, followed by the same Gibbs average, also gives

```math
\boxed{
\frac1r\sum_{i=1}^r A_i
\ge R_{\beta,r}(B).
}
```

Since `D_{beta(1-2/r)}>=D_beta`, a simpler consequence is

```math
\left(\prod_iA_i\right)^{1/r}\ge e^{-2\beta q/r}.
```

In reward form the exact inequality is

```math
\boxed{
\frac1r\sum_i\kappa_{\beta,i}(B)
\le
\frac{2q}{r}
-\frac1\beta
\log\frac{D_{\beta(1-2/r)}(B)}{D_\beta(B)}
\le\frac{2q}{r}.
}
```

Thus completeness forces at least one **large affinity / small reward**.  It
does not force the small affinity required by (10.624).  One-edge witnesses
only add low-deficit terms to the sums and cannot reverse this inequality.

There is a second exact product form.  Let

```math
d_i=q-Q(B[-i]).
```

Using (10.623) coordinate by coordinate gives

```math
\boxed{
\prod_{i=1}^r A_i
=2^r e^{-\beta\sum_i d_i}
\frac{\prod_iD_\beta(B[-i])}{D_\beta(B)^r}.
}
```

Equivalently, Hölder applied after lifting every child deficit to the parent
cube reproduces the same lower bound.  Exact `Q`-minimality says only
`d_i<=q_r-q_{r-1}` at a minimizing root; it gives no upper bound on the
child-deficit product capable of reversing Hölder.

## 3. Exact unordered-pair representation

There is a complementary identity which does yield a lower bound on some
reward.  Group the parent states into unordered pairs
`{omega,tau_i omega}`.  In each pair choose the endpoint with
`h_i=u>=0`, and give the pair its normalized total Gibbs mass `pi_i`.  The
other endpoint has field `-u` and weight ratio `e^(-4 beta u)`.  Therefore

```math
\boxed{
A_i
=\mathbb E_{\pi_i}\operatorname{sech}(2\beta u),
\qquad
m_i:=\mathbb E_{\nu_B}h_i
=\mathbb E_{\pi_i}[u\tanh(2\beta u)].
}
```

Use

```math
1-\operatorname{sech}(2x)
=\tanh(2x)\tanh x
```

and concavity of `tanh` on `[0,infinity)`.  Since `0<=u<=r-1`,

```math
\tanh(\beta u)
\ge\frac{u}{r-1}\tanh(\beta(r-1)).
```

It follows exactly that

```math
\boxed{
1-A_i
\ge
\frac{\tanh(\beta(r-1))}{r-1}\,m_i.
}
```

The complete energy law now gives

```math
\sum_i m_i=\mathbb E_{\nu_B}e_B.
```

Hence some coordinate satisfies the genuine lower bound

```math
\boxed{
\max_i\kappa_{\beta,i}(B)
\ge
-\frac1\beta\log\left(
1-
\frac{\tanh(\beta(r-1))}{r(r-1)}
\mathbb E_{\nu_B}e_B
\right).
}
```

Finally, the Gibbs variational identity gives
`E_nu e_B>=P_beta(B)>=q-r log(2)/beta`.  Therefore

```math
\boxed{
\max_i\kappa_{\beta,i}(B)
\ge
-\frac1\beta\log\left[
1-
\frac{\tanh(\beta(r-1))}{r(r-1)}
\left(q-\frac{r\log2}{\beta}\right)_+
\right].
}
```

For fixed `beta` and `q=alpha r^(3/2)+o(r^(3/2))`, this is only

```math
\frac{\alpha+o(1)}{\beta\sqrt r},
```

not `c sqrt(r)`.  The loss is sharp for the information used: a pair law may
put probability `Theta(r^(-1/2))` at `u=Theta(r)` and all remaining mass at
`u=0`.  It then has mean field `Theta(sqrt(r))` but affinity
`1-Theta(r^(-1/2))`.  Ruling out precisely this positive-heavy
concentration requires new minimizer structure; total energy and Hölder do
not do so.

## 4. Why edge witnesses do not improve the aggregate

For every low witness from (10.685), nonnegativity of its coordinate-flipped
deficit gives

```math
h_i(\omega)\ge-\Delta_\omega/4\ge-1.
```

Also

```math
\sum_i h_i(\omega)=q-\Delta_\omega.
```

Thus each witness has a coordinate with field at least `(q-4)/r`, and its
geometric contribution at that coordinate is small.  The coordinate may
depend on the witness.  Even a constant-size edge cover can put its large
field on different coordinates, while other low states dominate the Gibbs
denominator.  Summing or multiplying these witness contributions invokes
the lower-product Hölder inequality above.  No common-coordinate upper bound
results.

Exact global minimality beyond its one-edge witnesses could conceivably rule
out this field migration or the heavy-field mixture.  No such implication is
proved here.  In particular, this note neither constructs a complete-signing
family with bounded reward nor falsifies (10.617).

## 5. A8/A9 audit

The verifier evaluates every lifted state, every coordinate partner, every
child deficit sum, and every complete edge witness.  Representative values
are:

| matrix | `beta` | Holder `R` | geometric mean `BC` | actual `max kappa` | pair-field lower bound |
|---|---:|---:|---:|---:|---:|
| `A_8` | 0.5 | 0.103070 | 0.323374 | 2.257890 | 0.806682 |
| `A_9` | 0.5 | 0.074182 | 0.341042 | 2.665415 | 0.795337 |
| `A_8` | 1.0 | 0.007289 | 0.104231 | 2.261149 | 0.436763 |
| `A_9` | 1.0 | 0.004847 | 0.248269 | 1.760116 | 0.405331 |

Both exact minimizers satisfy all one-edge certificates.  The Hölder lower
bound can be very slack and always has the wrong direction.  The pair-field
bound is valid but decays at the asymptotic competitive scale.

## 6. Scoped conclusion

- **Verified:** the statewise products, exact child product, Hölder and
  arithmetic lower-affinity bounds, unordered-pair identities, and weak
  max-reward lower bound.
- **Closed mechanism:** aggregate/product/Hölder use of
  `Delta=q-2 sum_e s_e`, even supplemented by raw one-edge witnesses, cannot
  produce the desired direction.  It controls average reward from above.
- **Still open:** a genuinely minimizer-specific theorem excluding migration
  of the heavy coordinate across near states or excluding a rare
  `u=Theta(r)` positive-field component.  Such a theorem would have to use
  more than complete energy accounting and edge-flip incidence.
- **Quantifiers:** no asymptotic minimizer counterexample is supplied, so
  neither the literal every-minimizer version nor the target-specific
  convergence version of (10.617) is falsified.

All claims are audited in `complete_affinity_holder_r21.py`.
