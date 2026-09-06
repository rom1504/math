# Wave 39 root audit: row-priced complement witnesses

## Status

This note gives an exact pointwise bicriteria reformulation and a scoped
finite warning.  It does **not** control cross-selector congestion and hence
does not prove (10.1047).

For a selector `S`, retain

```math
B_S=-A+2P_SAP_S,qquad L_S(d)=\langle B_S,d\rangle,qquad R(d)=R_2(d).
```

For `alpha>=0`, choose an oriented cut maximizing

```math
\Phi_{S,\alpha}(d)=\alpha L_S(d)-R(d).
```

Then every maximizer satisfies the exact row inequality

```math
\boxed{R(d)\le n(n-1)+\alpha L_S(d)
             \le n(n-1)+3\alpha q_n.}
\tag{R39.P1}
```

Define the critical row price

```math
\alpha_S^*=\inf\{\alpha\ge0:
 \text{some maximizer of }\Phi_{S,\alpha}\text{ has }L_S(d)\ge q_n\}.
\tag{R39.P2}
```

The infimum is attained, and maximizing energies are nondecreasing with the
price.  Consequently

```math
\max_{|S|=m}\alpha_S^*=O(n^{3/4-c})
\tag{R39.P3}
```

would give a row-good complement certificate for every selector at the
project cap.  This still gives only pointwise coverage: a proof of (10.1047)
also needs the resulting witnesses to collide with fractional weight
`exp{O(n^(3/4-c))}`.

## Proof of the row inequality

Write `d=(sigma,x)` and set

```math
T_i=x_i(A^2x)_i-(n-1),
\qquad v_i=\sigma x_i(B_Sx)_i.
```

Flipping physical vertex `i` changes the two objectives by

```math
R(d^i)-R(d)=-4T_i,
\qquad L_S(d^i)-L_S(d)=-4v_i.
```

Global optimality of `d` gives `T_i<=alpha v_i` for every `i`.  Summing and
using

```math
\sum_iT_i=R(d)-n(n-1),
\qquad \sum_iv_i=L_S(d)
```

proves the first inequality in (R39.P1).  Orientation reversal negates `L`
without changing `R`, so a positive-price maximizer has `L>=0`; (10.1071)
then gives `L<=3q_n`.

For `alpha_2>alpha_1`, add the two optimizer inequalities for maximizers
`d_1,d_2` to obtain

```math
(\alpha_2-\alpha_1)(L_S(d_2)-L_S(d_1))\ge0.
```

Thus optimizer energy is monotone (including across ties in the displayed
paired sense), and finiteness of the cut cube proves attainment in (R39.P2).

The unconditional bound is far weaker.  A `B_S`-ground has `L=M_S>=q_n`
and `R<=n||A||_op^2<=2nq_n`.  Every oriented energy lies on an integer
lattice of fixed spacing, so comparison with every lower-energy profile
only gives

```math
\alpha_S^*=O(nq_n)=O(n^{5/2}).
```

The tempting `O(n)` estimate obtained by comparing the ground with one
typical low-energy cut is invalid: a different profile can have energy only
one lattice step below `M_S` and much smaller row.  Thus (R39.P3) is a real
power-scale slope theorem, not a consequence of ordinary averaging.

## Exact finite audit and congestion warning

`tmp/complement_row_price_r39.py` enumerates all profiles using rational
critical prices.  At the common price `max_S alpha_S^*` it also computes the
fractional cover of price-optimal threshold witnesses (the LP values are
floating-point diagnostics):

| case | critical-price histogram | common-price active oriented cuts | common-price fractional cover |
|---|---|---:|---:|
| `A_6,m=5` | `0^6` | 32 | 2 |
| `A_8,m=6` | `2^8,4^16,(14/3)^4` | 32 | 16 |
| `A_9,m=7` | `2^10,3^18,(7/2)^1,5^2,6^2,7^2,8^1` | 50 | 36 |

For comparison, the unrestricted numerical `A_9` cover is `13.75`.
Thus a common row price can select markedly less congested witnesses even at
finite order.  The data do not falsify (R39.P3) or (10.1047), but they show
that pointwise row-price control is not a substitute for the fractional
collision theorem.

Run:

```bash
.venv/bin/python tmp/complement_row_price_r39.py
```

The current run ends with `PASS complement_row_price_r39`.
