# Wave 18 Route 2: finite-window pressure descent

## Status

- **Verified:** the exact fixed-temperature window telescope, its sharp
  entropy-corrected form, the hierarchy between pairwise endpoints and a
  single nested deletion chain, and the exact changing-temperature correction.
- **Verified no-go:** at a fixed temperature, cavity rewards form a gradient
  on the subset lattice.  Their sum depends only on the terminal restriction,
  so reordering or shortest-path ``banking'' cannot create a second scalar.
- **Verified no-go:** both the fixed-orientation and orientation-inclusive
  Finner inequalities give a lower bound on average restricted pressure at a
  hotter temperature.  This is the wrong direction for selecting a small
  restriction.  The reverse contraction already fails on the order-nine
  minimizer.
- **Open target:** prove a uniform, all-fixed-proportion upper bound on the
  minimum restricted pressure (or norm).  A cumulative power-saving error is
  enough; the pointwise constant-shortfall condition (10.617) is much stronger
  than necessary.

## 1. The exact weakest window telescope

Let `A` be an order-`n` global minimizer and put

```math
\alpha=\frac{q_n}{n^{3/2}}.
```

Fix one inverse temperature `beta` for the entire root window.  Let

```math
B_n=A\supset B_{n-1}\supset\cdots\supset B_m
```

be any induced deletion chain, and let `i_r` be the vertex removed from
`B_r`.  Define the root-scale edge shortfall and its cumulative value by

```math
\delta_r
=\alpha\bigl(r^{3/2}-(r-1)^{3/2}\bigr)
-\kappa_{\beta,i_r}(B_r),
\qquad
\Delta_{n,m}=\sum_{r=m+1}^n\delta_r.
\tag{WP1}
```

The cavity identity telescopes exactly:

```math
\sum_{r=m+1}^n\kappa_{\beta,i_r}(B_r)
=P_\beta(A)-P_\beta(B_m),
```

and therefore

```math
\Delta_{n,m}
=\alpha(n^{3/2}-m^{3/2})
-\bigl(P_\beta(A)-P_\beta(B_m)\bigr).
\tag{WP2}
```

Since `P_beta(A)<=Q(A)=q_n` and
`Q(B_m)<=P_beta(B_m)+m log(2)/beta`, one obtains

```math
\boxed{
\frac{q_m}{m^{3/2}}
\le
\frac{q_n}{n^{3/2}}
+\frac{\Delta_{n,m}}{m^{3/2}}
+\frac{\log2}{\beta\sqrt m}.
}
\tag{WP3}
```

Consequently, for one pair `(n,m)` the weakest useful pressure hypothesis is
only the existence of a minimizer `A`, an order-`m` restriction `B_m`, and a
fixed-temperature chain to it for which

```math
\Delta_{n,m}=o(m^{3/2}),
\qquad
\beta\sqrt m\longrightarrow\infty.
\tag{WP4}
```

No stepwise inequality is required.  Positive reward surplus can be banked
against later shortfalls; only the cumulative balance at the landing order
matters.  Condition (10.617) implies
`Delta_{n,m}<=K(n-m)=O(n)`, whereas (WP4) permits any cumulative
`o(n^{3/2})` error in a fixed-proportion window.

There are three logically distinct uniform versions.

1. **Pairwise endpoint version (weakest):** for every
   `m in [rho n,n]`, one may choose a different minimizer and a different
   order-`m` restriction satisfying (WP4).
2. **Fixed-root cardinality profile:** one minimizer works for every such
   `m`, but the good restrictions need not be nested.
3. **One-chain prefix version:** one deletion order obeys

   ```math
   \sup_{\rho n\le m\le n}\Delta_{n,m}=o(n^{3/2}).
   \tag{WP5}
   ```

For inequalities between the numbers `q_n`, the pairwise version is enough:
each directed comparison may use a fresh minimizer.  A supermartingale or a
single recursive restriction construction needs the stronger nested version.
In every case the condition must hold at every requested fixed-proportion
landing order.  The telescope at `m=0` alone says nothing about where the
reward is located.

For the adaptive convergence tail, plain `o(1)` comparison costs need not be
summable.  If `eta(N)` bounds the right-hand error in (WP3) for every edge in
a fixed-ratio window, a sufficient rate is

```math
\sum_{j\ge0}\eta\!\left(\left\lceil\rho^{-j}N\right\rceil\right)
\longrightarrow0.
\tag{WP6}
```

Any power saving `eta(N)=O(N^{-c})`, `c>0`, suffices.  Equivalently, a
cumulative pressure shortfall `Delta=O(n^{3/2-c})` is already enough and is
strictly weaker than the `O(n)` consequence of (10.617).

## 2. Exact entropy correction: pressure banking is endpoint selection

For an order-`r` matrix `B`, define

```math
D_\beta(B)=\sum_{\omega\in\Omega_r}
e^{-\beta(Q(B)-e_B(\omega))},
\qquad
L_\beta(B)=Q(B)-P_\beta(B)
=\frac{r\log2-\log D_\beta(B)}{\beta}.
\tag{WP7}
```

Thus `0<=L_beta(B)<=r log(2)/beta`.  Because the root is a minimizer,
(WP2) has the exact sharpening

```math
\boxed{
Q(B_m)-\alpha m^{3/2}
=\Delta_{n,m}+L_\beta(B_m)-L_\beta(A).
}
\tag{WP8}
```

The right side is independent of temperature after all terms are combined.
At one deletion, if `d_i=Q(B)-Q(B[-i])`, then

```math
\kappa_{\beta,i}(B)
=d_i+L_\beta(B[-i])-L_\beta(B).
\tag{WP9}
```

Thus every apparent finite-temperature gain beyond the norm decrement is
exactly an increase of the soft-max entropy gap.  Finite temperature provides
a smoother objective on which one might prove an endpoint-selection theorem,
but it does not provide an independent terminal credit.

There is also an exact subset-lattice formulation.  Put

```math
F_{\alpha,\beta}(S)=P_\beta(A[S])-\alpha|S|^{3/2}.
```

For an edge `S -> S\setminus{i}` its increment is precisely `delta` from
(WP1).  Hence, for every two deletion orders with the same endpoint `T`,

```math
\sum\delta=F_{\alpha,\beta}(T)-F_{\alpha,\beta}(V).
\tag{WP10}
```

The shortest-path problem at level `m` therefore collapses to
`min_{|T|=m}F(T)`; order-dependent reward banking cannot improve a fixed
endpoint.  Controlling all levels along one chain remains a genuine nesting
compatibility problem because minimizers of `F` at different cardinalities
need not be nested.

## 3. Temperature regimes and the exact price of changing temperature

Suppose `m>=rho n` and `beta_n` is frozen throughout the root window.
The entropy term in (WP3) has the following exact scales.

- If `beta_n=beta_0>0`, it is at most
  `log(2)/(beta_0 sqrt(rho n))=O(n^{-1/2})`.
- If `beta_n=b/sqrt(n)` with fixed `b`, it is at most
  `log(2)/(b sqrt(rho))`, which does **not** tend to zero.  If a cumulative
  lemma were available for every fixed `b` with an `o_n(1)` error, one could
  take `n -> infinity` first and then `b -> infinity`, but a single fixed `b`
  does not prove (WP4).
- If `beta_n=b_n/sqrt(n)`, then the entropy term is
  `O_rho(1/b_n)` and is `o(1)` exactly when `b_n -> infinity`.  For the
  adaptive tail (WP6), one additionally needs the geometric tail of
  `1/b_n` to be summable; for example `b_n=n^c`, `c>0`, works.

Freezing the root temperature avoids an extra drift.  If instead a chain uses
temperatures `beta_r`, set

```math
\tau_r
=P_{\beta_r}(B_{r-1})-P_{\beta_{r-1}}(B_{r-1}).
```

Then the exact telescope is

```math
P_{\beta_n}(A)-P_{\beta_m}(B_m)
=\sum_{r=m+1}^n\bigl(\kappa_{\beta_r,i_r}(B_r)+\tau_r\bigr).
\tag{WP11}
```

The pressure derivative is

```math
\partial_\beta P_\beta(B)
=\frac{\operatorname{Ent}(\mu_{\beta,B}\mid\mathrm{unif})}{\beta^2},
\qquad
0\le\partial_\beta P_\beta(B)
\le\frac{|B|\log2}{\beta^2}.
\tag{WP12}
```

Consequently, if `beta_{r-1}>=beta_r`,

```math
-(r-1)\log2
\left(\frac1{\beta_r}-\frac1{\beta_{r-1}}\right)
\le\tau_r\le0.
\tag{WP13}
```

For the natural but changing choice `beta_r=b/sqrt(r)`, the total omitted
temperature correction can be as large as

```math
\frac{\log2}{b}
\sum_{r=m+1}^n(r-1)(\sqrt r-\sqrt{r-1})
=O\!\left(\frac{n^{3/2}-m^{3/2}}b\right).
\tag{WP14}
```

It is a leading fixed-window error for fixed `b`; it becomes subleading if the
corresponding `b_n -> infinity`.  A telescope that changes temperature but
omits `tau_r` is invalid.

## 4. Random fixed-size restrictions and the Finner direction

Let

```math
Z_A^\sigma(\beta)=\mathbb E_x
e^{\beta\sigma x^{\mathsf T}Ax}.
```

For all order-`m` subsets `T`, put
`d=binom(n-1,m-1)`.  Each coordinate belongs to exactly `d` subsets and
each edge belongs to `binom(n-2,m-2)` subsets.  Finner with weight `1/d`
therefore gives the fixed-orientation inequality

```math
\boxed{
Z_A^\sigma(\beta)
\le
\prod_{|T|=m}
Z_{A[T]}^\sigma\!\left(\beta\frac{n-1}{m-1}\right)^{1/d}.
}
\tag{WP15}
```

The shared orientation variable cannot be added with these weights: it lies
in every factor and its total cover weight is
`binom(n,m)/d=n/m>1`.  It can be included only after reducing every weight to
`1/binom(n,m)`.  Viewing the oriented projective partition function as the
normalized expectation over independent `(sigma,x)` then yields the valid
orientation-inclusive statement

```math
\boxed{
\mathcal Z_\beta(A)
\le
\prod_{|T|=m}
\mathcal Z_{\beta n(n-1)/(m(m-1))}(A[T])^{1/\binom nm}.
}
\tag{WP16}
```

Equivalently,

```math
\mathbb E_T
P_{\beta n(n-1)/(m(m-1))}(A[T])
\ge
\frac{m(m-1)}{n(n-1)}P_\beta(A).
\tag{WP17}
```

This is a **lower** bound on average restricted pressure, at a hotter
temperature.  Selected descent requires an upper bound that ensures at least
one small endpoint.  Thus Finner/Shearer has the wrong direction.  At zero
temperature (WP17) reduces to the elementary lower bound obtained by
restricting a root ground state.

A hoped-for reverse normalized contraction is genuinely false, not merely
unproved.  For the exact `A_9` minimizer all nine order-eight restrictions
have norm `24`, whereas

```math
\left(\frac89\right)^{3/2}Q(A_9)
=\frac{128}{9}\sqrt2
=20.11325955\ldots<24.
```

Hence the reverse pressure contraction fails for all sufficiently large
`beta` by continuity.  Random fixed-size averaging alone cannot bypass the
flat child wall.

## 5. Exact `A_9` flat-chain audit

Take the certified chain

```math
A_9\longrightarrow B_8=A_9[-7]
\longrightarrow B_7=A_9[-\{7,0\}].
```

It has

```math
Q(A_9)=24,
\qquad Q(B_8)=24,
\qquad Q(B_7)=22,
\qquad \alpha=\frac89.
```

Thus the total norm decrements are `0` and `2`, while the exact endpoint
excesses over the root scale are

```math
Q(B_8)-\alpha8^{3/2}
=24-\frac{128}{9}\sqrt2
=3.886740446\ldots,
\tag{WP18}
```

and

```math
Q(B_7)-\alpha7^{3/2}
=22-\frac{56}{9}\sqrt7
=5.537547398\ldots.
\tag{WP19}
```

At every fixed temperature, (WP9) makes the cancellation transparent:

```math
P_\beta(A_9)-P_\beta(B_8)
=L_\beta(B_8)-L_\beta(A_9),
\tag{WP20}
```

and

```math
P_\beta(A_9)-P_\beta(B_7)
=2+L_\beta(B_7)-L_\beta(A_9).
\tag{WP21}
```

So the pressure reward on the zero-decrement step is *entirely* paid back by
the terminal entropy gap.  For example, with the root-frozen temperature
`beta=(1)/sqrt(9)=1/3`, exact energy enumeration gives

| endpoint | cumulative pressure reward | pressure shortfall `Delta` | `L(endpoint)-L(A_9)` | exact excess |
|---|---:|---:|---:|---:|
| `B_8` | 2.024216818 | 1.862523628 | 2.024216818 | 3.886740446 |
| `B_7` | 5.360214382 | 2.177333016 | 3.360214382 | 5.537547398 |

The final column is exactly the sum prescribed by (WP8).  The checker
`tmp/check_window_pressure_r18.py` verifies these identities for six
temperatures using exact integer energy spectra and floating-point
log-sum-exp only.

## 6. Ledger-ready conclusion

The stepwise constant-shortfall target (10.617) can be weakened rigorously to
an all-landing-order cumulative condition: for every fixed-ratio pair
`m in [rho n,n]`, select an order-`m` restriction for which the total
root-scale cavity shortfall is `o(n^{3/2})`, uniformly, while freezing a
temperature satisfying `beta_n sqrt(n) -> infinity`.  A power-saving
cumulative error gives the summable adaptive tail and would prove
convergence.  This permits arbitrary bad individual steps and genuine reward
banking across a window.

However, fixed-temperature cavity rewards are a gradient, so shortest-path
banking is exactly endpoint selection.  The entropy-corrected pressure
telescope is identically the zero-temperature endpoint norm, and `A_9` shows
the cancellation on a flat chain at every temperature.  Random restriction
plus Finner supplies only a lower average-pressure bound at a hotter
temperature, the opposite of what endpoint selection needs.  The surviving
new target is therefore a matrix-specific upper bound on
`min_{|T|=m}P_beta(A[T])` (or directly on `min Q(A[T])`), uniform at every
fixed-proportion cardinality, with any geometrically summable normalized
error.
