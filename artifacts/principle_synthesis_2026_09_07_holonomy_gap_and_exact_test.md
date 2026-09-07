# The quantitative holonomy obligation and an exact selectable-child test

2026-09-07. Status: exact criterion and exact finite contraction calculation.
No asymptotic seed-value transfer is claimed. This follows the actual
fixed-seed construction in `principle_synthesis_2026_09_07_fixed_seed_mixed_orbit_weave.md`.

## 1. A concrete port-tensor gap, not a covariance condition

Fix the physical retention `k`, macro order `m`, depth `r`, and `t>0`.
Take the actual recursive row basis with NO final column signs, a uniform
physical spin, and an ordinary output permutation. Normalize its spectrum
by `sqrt(k)` and delete the coordinate assigned to the macro diagonal.
Let `nu_m` be the resulting exchangeable law on `R^(m-1)`. Define the
Gaussian port tensor

```math
T_m=\mathbb E_{v\sim\nu_m}\bigotimes_{j=1}^{m-1}\phi_t(v_j),
\quad
\rho_m=\|T_m\|,
\quad
\rho_m^2=\mathbb E_{v,w\sim\nu_m}
                   e^{-t\|v-w\|^2}.                       \tag{1}
```

The two vectors in the last expectation are independent. For a fixed
macro signing `S`, put

```math
\mathcal C_m(S)=\mathbb E_{v_1,\ldots,v_m\sim\nu_m}
 \prod_{i<j}e^{-t(v_i(j)-S_{ij}v_j(i))^2},
\quad
\overline{\mathcal C}_m(S)
 =\frac{\mathcal C_m(S)+\mathcal C_m(-S)}2.                  \tag{2}
```

This is a finite, nonnegative, explicitly specified Gaussian-kernel
contraction of the common row tensor. Its definition averages the actual
row law before the spin maximum; no maximizing parent configuration enters
it. Reflection-twisted graph Cauchy--Schwarz gives
`overline C_m(S)<=rho_m^m`. Define the nonnegative gap

```math
\Gamma_m(S)=\frac1{m^2}
 \log\frac{\rho_m^m}{\overline{\mathcal C}_m(S)},
\qquad
a_m=\frac{k}{m}\log2+\frac1m\log\rho_m.                   \tag{3}
```

The spin sum and two energy orientations are exactly
`2^(mk+1) overline C_m(S)`. Thus, for any `gamma>=0`, the probability
that the fixed-seed weave has a spin with
`D_sigma^off<=2 gamma m^2 k` is at most

```math
2\exp\{m^2[t\gamma+a_m-\Gamma_m(S)]\}.                   \tag{4}
```

If `c_m=Q(S)/m^(3/2)` and `p_m=k/m`, the direct selectable seed target
is therefore supplied by the explicit inequality

```math
\Gamma_m(S)\ge a_m+t(1-2\sqrt{p_m}\,c_m)+\epsilon_m,
\qquad \epsilon_m m^2\longrightarrow\infty.               \tag{5}
```

It implies actual hollow output cap at most
`(c_m+o(1))(mk)^(3/2)`. For selectable exact-child scope, (5) need only
hold for ONE signing with `Q(S)=M_m`, not for all of them. Equivalently
the maximum of `Gamma_m(S)` over that exact optimizer set must meet the
displayed threshold. Parameters and depth may depend on a fixed desired
accuracy, but a growing-depth implementation would require its own uniform
error proof.

The two-replica row norm `rho_m` can be smaller than the sum of individual
orbital norms used in the old certificate. The gap in (3) separates this
possible row-level improvement from the genuinely joint contraction gap.
The coarse old exponent can replace `a_m` by its known upper bound, but
that gives a stronger sufficient condition, not an equivalence.

There is also a genuine architectural restriction: reciprocal-column
matching gives the actual fixed-retention family floor
`sqrt(p)/2-o(1)` (see the archived reciprocal-weave matching theorem).
Consequently reaching a child value `c` requires `p<=4c^2+o(1)` within
this family. At `p=.96` the floor is about `.489898`; a possible value
near the original lower bound would require substantially smaller
retention. No improvement of the existing fixed-p contraction alone
removes this floor.

## 2. An actual Hadamard row sector with a finite exact contraction

Take macro order `m=8`, depth zero, the Sylvester Hadamard `H_8`, and
physical retention `k=4`. Under the fresh signed input permutation, the
physical word is uniform on the `binom(8,4)2^4=1120` words having four
entries in `{+-1}` and four zero entries. Exactly 112 of them have a
flat normalized spectrum

```math
v=H_8^Tf/2\in\{\pm1\}^8.                                \tag{6}
```

Among these 112 words, 28 spectra have two positive coordinates, 56
have four, and 28 have six. Apply the actual ordinary output permutation
and condition on this row sector. Its exchangeable distinct-coordinate
moments, up through degree seven, are exactly

```math
\mu_0=1,\qquad \mu_4=-1/35,
\qquad \mu_1=\mu_2=\mu_3=\mu_5=\mu_6=\mu_7=0.             \tag{7}
```

The values follow either from the enumerated 112 integer words or from
the coefficient of `z^d` in
`(1+z)^a(1-z)^(8-a)` averaged over the three displayed positive counts.
Deleting the assigned diagonal coordinate preserves these moments on
the seven remaining distinct coordinates.

For `a,b` signs, write

```math
e^{-t(a-Sb)^2}=c_t(1+u_t Sab),
\quad c_t=\frac{1+e^{-4t}}2,
\quad u_t=\tanh(2t).                                     \tag{8}
```

Expanding over edge subsets of `K_8`, (7) kills every term except those
whose vertices have degree zero or four. If `A_s(S)` is the sum of
`prod_{e in F} S_e` over all simple four-regular subgraphs on `s` active
vertices, the EXACT contraction is therefore

```math
\frac{\mathcal C_8^{\rm flat}(S)}{c_t^{28}}
 =1+\sum_{s=5}^{8}A_s(S)z_t^s,
\qquad z_t=-\frac{u_t^2}{35}\in(-1/35,0).                 \tag{9}
```

The denominator is exactly the contraction after averaging independent
macro signs. Also `C_8^flat(-S)=C_8^flat(S)`, because each surviving
subgraph has `2s` edges. This computation concerns a real row sector of
the fixed-seed compiler, not an invented covariance model.

Only 23,551 subgraphs are needed: 56 with five active vertices, 420 with
six, 3,720 with seven, and 19,355 with eight. There is no full parent-spin
enumeration and no summation over all incidence words.

## 3. Actual exact minimizers have opposite strict responses

Use the four rooted order-eight optimizer representatives indexed by
seven-vertex graph-atlas entries 580, 722, 870 and 1005. The root has
positive incident signs and the atlas graph supplies negative internal
edges. Their cap is exactly `M_8=10`; the archived order-eight optimum
certificate supplies the lower bound, and all 128 projective spin tests
independently verify the four witnesses. The complete integer coefficients
in (9) are:

| Macro seed | Q(S) | A_5 | A_6 | A_7 | A_8 |
| --- | ---: | ---: | ---: | ---: | ---: |
| all positive | 28 | 56 | 420 | 3720 | 19355 |
| atlas580 or atlas1005 | 10 | -8 | 4 | 72 | -5 |
| atlas722 or atlas870 | 10 | 24 | 68 | -24 | -5 |

These signs are not merely small-temperature expansions. Factoring
`R_S(z)-1=z^5 q_S(z)`, exact Bernstein coefficients of the cubic `q_S`
on `[-1/35,0]` show

```math
R_{580}(z)>1,\qquad R_{722}(z)<1
\quad\text{for every }z\in[-1/35,0).                     \tag{10}
```

For example the four Bernstein coefficients for atlas580 are
`-2763/343, -29608/3675, -844/105, -8`, all negative; those for atlas722
are `188973/8575, 83416/3675, 2452/105, 24`, all positive. Since `z^5<0`,
the signs in (10) follow. The all-positive seed also has `R(z)<1` on
this interval, despite having a much larger original cap.

Thus selecting an actual exact child CAN improve this exact reflection
contraction relative to independent macro signs, while another exact
child worsens it, for every positive kernel parameter. The input cap
alone does not determine even the sign of this finite response.
This does not defeat the selectable quantifier: atlas722 is a favorable
selection for this sector. It also does not prove a favorable choice
simultaneously for every other row sector.

For this sector the port norm itself is explicit:

```math
\rho_8^2=c_t^7(1+u_t^4/35),
\quad
\Gamma_8^{\rm flat}(S)
 =\frac1{64}\left[4\log(1+u_t^4/35)-\log R_S(z_t)\right].  \tag{11}
```

The first part of (11) is seed-independent. The exact polynomial shows
the additional seed-dependent correction separately.

## 4. Verification and the unclosed scale

The full integer calculation is reproducible with
`computations/principle_synthesis_2026_09_07_fixed_seed_holonomy_polynomial.py`.
Its complete matrices, moment counts, subgraph counts, coefficients and
rational Bernstein checks are in the corresponding results JSON. All
checks PASS. This is an exact finite calculation, not an interval or
floating-point optimization.

This sector has probability `1/10` for one row. Controlling its
contraction does not control all row sectors in the parent spin union.
No original cap bound or asymptotic recurrence is inferred from (10).

The director's fourth-moment scale audit also remains important: for
fixed positive retention, the full-spin fourth moment in four distinct
spectral coordinates is only `O_r(m^-2)`, so a fixed `K_5` term is
`O_r(m^-10)` and the sum of all such terms is still `O_r(m^-5)`.
A leading gap in (5) requires more than preserving or optimizing a
bounded collection of low-degree cycles.

The independent seed-pressure theorem now extends the earlier finite-type
audit to ALL uniformly bounded profiles, including asymmetric ones, and
then to uniformly square-integrable profile families by quantization and
clipping. See `principle_construct_2026_09_07_uniform_integrable_seed_pressure.md`.
Thus asymmetry by itself does not escape at fixed `t` for low-cap seeds.
The actual all-spin row family fails uniform square-integrability because
restricted column words have `sqrt(k)` spikes carrying a positive energy
fraction. The unresolved possibility lies in those coherent tails,
growing temperature, cross-row dependence, or a different small-retention/
constant-multiplier architecture. None has been replaced by a theorem here.
