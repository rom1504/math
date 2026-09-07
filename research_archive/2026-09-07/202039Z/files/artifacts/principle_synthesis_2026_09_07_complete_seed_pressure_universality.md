# Complete fixed-temperature seed-pressure universality

2026-09-07. **Proved jointly with the construction auditor; full independent
reconstruction PASS, director audit requested.** The new step is the removal
of ALL square-integrability-of-tails restrictions. The proof combines the
auditor's descending signed-port matching with the category-entropy and
tagged-zero conditioning argument developed here.

This is a theorem about an actual positive-kernel permutation partition,
not about actual parent maxima. It closes leading seed sensitivity in the
complete fixed-temperature compiler of
`principle_synthesis_2026_09_07_fixed_seed_mixed_orbit_weave.md`, despite that
compiler's exact finite cycle dependence. It does not rule out a different
seed-dependent basis law, correlated port placement, a singular temperature,
or another construction under the user's allowance of global edge changes.

## 1. The theorem

At vertex `i` prescribe a real multiset `u_i` of length `d=m-1` and assume

```math
\sum_{a\in u_i}a^2\le C m\quad\text{for every }i.          \tag{1}
```

Independently uniformly permute each row among its directed ports. For a
hollow macro signing `S`, define

```math
Z_S(u;t)=\mathbb E\prod_{i<j}
       e^{-t(u_{ij}-S_{ij}u_{ji})^2},\qquad t>0.            \tag{2}
```

Fix finite constants `C,K,t`. Uniformly over EVERY profile array obeying
(1), including arbitrarily many atoms, signed asymmetry, and coherent
spikes as large as `sqrt(Cm)`, and every pair of seeds obeying

```math
Q(S),Q(T)\le K m^{3/2},
```

one has

```math
\boxed{\ |\log Z_S(u;t)-\log Z_T(u;t)|=o(m^2).\ }          \tag{3}
```

The little-oh is uniform in the profile array and the seeds. No uniform
integrability premise is needed. The per-ROW bound in (1), rather than
only a total bound, is used in the heavy matching lemma.

The bounded-profile version of (3) is the earlier finite-type theorem
plus the auditor's uniform quantization/pruning extension, preserved in
`principle_construct_2026_09_07_uniform_integrable_seed_pressure.md`.
The rest of this note proves the missing passage from (1) to bounded
profiles. It does not repeat or rebrand the earlier finite-type argument.

**Dependency hypothesis audit.** I read the original finite-type proof
`principle_invent_2026_09_07_finite_type_seed_universality.md` completely.
Its alphabet is one fixed global involutive alphabet, but each row may
have an arbitrary support, with zero masses handled by restriction to
that LOCAL support. The common bounded grid used by the auditor includes
zero and both signs. Rowwise pruning merges rare atoms into a row's own
most frequent atom; it therefore preserves membership in the common grid
while permitting different local supports. Every retained positive mass
is bounded below by the chosen fixed pruning threshold. A row consisting
only of zeros has a singleton support and satisfies the proof trivially.
The Gaussian kernel is strictly positive on the common finite grid and
has the required simultaneous-sign symmetry. Its smooth dual comparison
is on the full simplex, including the boundaries containing those local
supports. Thus row-specific zero atoms and asymmetric pruned supports do
not add a missing hypothesis to the bounded-profile input used below.

## 2. A deterministic perturbation bound

Let `v_i` be another list on the same labeled coordinate copies, obtained
by deleting some values or otherwise changing them. Suppose both total
energies are at most `Cm^2`, and
`sum_ports (u-v)^2<=eta m^2`. Use the SAME row permutations for the two
arrays. The edge-difference map has norm at most `sqrt(2)`, so

```math
\left|\log Z_S(u;t)-\log Z_S(v;t)\right|
 \le4t\sqrt{C\eta}\,m^2.                                 \tag{4}
```

This is pointwise before averaging and holds for every `S`. Treating
copies as labeled causes no change to the uniform multiset law. A fixed
relabeling/deletion of copies pushes a uniform permutation to the uniform
permutation law of the modified multiset, even when some equal copies are
deleted and others are not.

## 3. Descending heavy-bin matching

Assume first that all nonzero magnitudes are either at most `L` or at
least `V`, where `V>2L` and `V>=2`. Fix `0<delta<=1`. Put

```math
a_b=V(1+\delta)^b,
```

and bin the heavy magnitudes in `[a_b,(1+delta)a_b)`. There are at most
`B=O_delta(log m)` nonempty bins, since every coordinate is at most
`sqrt(Cm)`. Process the bins in DECREASING magnitude order.

For each bin, greedily match pairs of still available signed stubs at
distinct vertices using a previously unused macro edge. A pair of stub
signs `epsilon_i,epsilon_j` may use edge `ij` exactly when
`S_ij=epsilon_i epsilon_j`. Continue until the matching is maximal.
This is a finite constructive procedure; it need not find a maximum
matching. Every physical edge is used at most once over all bins.

At a bin with lower endpoint `a`, the degree already used at any vertex,
including the current bin, is at most

```math
D_a=Cm/a^2,                                               \tag{5}
```

because every used stub has magnitude at least `a` and (1) holds. Let
`U` be the vertices still having an unmatched stub in this bin. Choose
one available sign `epsilon_i` at each `i in U`. Maximality implies that
every UNUSED edge in `U` has `S_ij epsilon_i epsilon_j=-1`. There are at
most `|U|D_a/2` used edges in `U`. Hence the signed energy on this partial
spin is at most `|U|D_a-binom(|U|,2)`. Its absolute value cannot exceed
`Q(S)`, by principal cap monotonicity. Therefore

```math
\binom{|U|}{2}-|U|D_a\le Q(S),
\qquad |U|\le2D_a+1+\sqrt{2Q(S)}.                        \tag{6}
```

This argument allows either sign at a residual vertex and does not
require separate positive or negative degree assumptions on the seed.
Descending order is essential for the used-degree bound (5).

Delete all unmatched heavy stubs, replacing them by zero. The total
deleted energy `E_bad` obeys

```math
E_{\rm bad}
 \le Cm\sum_b|U_b|
 \le\frac{2C^2m^2}{V^2[1-(1+\delta)^{-2}]}
       +CmB\,[1+\sqrt{2K}\,m^{3/4}].                    \tag{7}
```

Thus, for fixed `C,K,delta,V`,

```math
\eta_{\rm bad}:=E_{\rm bad}/m^2
 \le \frac{A_{C,\delta}}{V^2}
       +O_{C,K,\delta,V}(m^{-1/4}\log m).                \tag{8}
```

The remaining heavy stubs have an actual compatible placement in pairs
on distinct macro edges. Within each pair the magnitudes lie in the
same relative bin and the signs agree through `S`. Consequently their
total squared edge defect is at most

```math
\frac12\delta^2 C m^2.                                  \tag{9}
```

This is the auditor's sparse coherent-port compiler in the form needed
here. It handles arbitrary signed, heterogeneous, and growing heavy
types, not just an exactly balanced incidence profile.

## 4. Heavy category placement has small entropy cost

Let `u^S` be the profile after the deletion in Section 3. At row `i`,
let `h_{ib,+},h_{ib,-}` count its remaining heavy stubs in the indicated
bin and sign. The matching prescribes the PORT SET of each such category;
it does not prescribe an ordering of individual values within the bin.
Under the original independent uniform row permutations the probability
of these placements is

```math
\prod_i\frac{(d-h_i)!\prod_{b,\sigma}h_{ib,\sigma}!}{d!},
\qquad h_i=\sum_{b,\sigma}h_{ib,\sigma}.                  \tag{10}
```

This category rather than copy-by-copy event is important. A bound of
the form `H log m` would be inadequate.

Here is a uniform entropy estimate for (10). Use the reference masses

```math
q_{b,+}=q_{b,-}=\frac{\delta}{8a_b^2},\qquad
q_0=1-\sum_{b,\sigma}q_{b,\sigma}.                       \tag{11}
```

Since `1-(1+delta)^(-2)>=delta/2`, their heavy sum is at most
`1/(2V^2)`, so `-log q_0<=1/V^2`. If `p` is any row's empirical
category law, its entropy is at most its cross-entropy against (11).
Moreover `sum p_{b,sigma}a_b^2<=Cm/d<=2C`. The function
`[log(8/delta)+2log a]/a^2` decreases for `a>=1`. It follows that

```math
H(p)\le h_{C,\delta}(V)
 :=\frac{1+2C[\log(8/\delta)+2\log V]}{V^2}.              \tag{12}
```

The multinomial upper bound `multinomial(d;dp)<=exp(dH(p))` is exact.
Therefore the probability in (10) is at least

```math
\exp[-h_{C,\delta}(V)m^2].                              \tag{13}
```

There is no minimum atom-mass assumption, no bounded-number-of-bins
assumption, and no Stirling error in this step.

## 5. Fixing the zero ports costs only bounded-light swaps

Let `ell` be the COMMON bounded profile obtained by replacing ALL heavy
values by zero and leaving the light values unchanged. It is independent
of `S`. Conditional on the event (10), collapse the placed heavy values
to zero. The remaining light values are still independently uniformly
permuted in their free ports.

Couple this conditional bounded profile with an unconditioned uniform
permutation of `ell` as follows. Artificially tag `h_i` zero copies in
row `i` and assign them bijectively to the prescribed heavy ports.
Sequentially swap each tag into its assigned port. A uniform permutation
pushes to the uniform law conditional on the tagged placements: at each
stage the unconditioned tag is uniform among the remaining positions,
and every resulting conditioned arrangement has the same number of
preimages. The previously fixed tags are not disturbed. After forgetting
tags, this is exactly the required conditional light law.

There are at most `H=sum_i h_i` swaps, and

```math
H\le Cm^2/V^2.                                          \tag{14}
```

All values in these two bounded profiles have magnitude at most `L`.
Changing one directed port changes its log kernel by at most `8tL^2`;
a swap changes at most two ports. Thus the coupled log weights differ
by at most `16tL^2H`, and

```math
\mathbb E[K_S(\ell)\mid\text{prescribed zero ports}]
 \ge e^{-16tCL^2m^2/V^2}\,Z_S(\ell;t).                   \tag{15}
```

This comparison concerns normalized permutation expectations. It does
not replace a rare conditioning event by its unnormalized weight.

On (10), every placed heavy edge has two heavy endpoints; all other
edges are light-light. Combining (9), (13), and (15) yields

```math
Z_S(u^S;t)\ge Z_S(\ell;t)
 \exp\left\{-m^2\left[
 h_{C,\delta}(V)+\tfrac12t\delta^2C
                 +16tCL^2/V^2\right]\right\}.             \tag{16}
```

Finally, by (4),(8), replacing `u` by `u^S` costs at most
`4t sqrt(C eta_bad)m^2` in log partition.

## 6. The separated-tail comparison

Removing the heavy values can only INCREASE every edge kernel when
`V>2L`. For a heavy-light edge,

```math
(a-Sb)^2-b^2=a^2-2Sab\ge a^2-2|a|L\ge0;
```

for a heavy-heavy edge the new squared difference is zero. The labeled
permutation pushforward is exact, so `Z_S(u;t)<=Z_S(ell;t)`.
Together with (16) this proves the quantitative, uniform estimate

```math
0\le\frac{\log Z_S(\ell;t)-\log Z_S(u;t)}{m^2}
\le4t\sqrt{C\eta_{\rm bad}}
  +h_{C,\delta}(V)+\tfrac12t\delta^2C+16tCL^2/V^2.          \tag{17}
```

The right side can be made arbitrarily small by first taking `delta`
small, then the separation ratio `V/L` large, then `V` large, and finally
`m` large. The descending matching depends on `S`, but the bounded
comparison profile `ell` does not.

## 7. Every bounded-energy array has a usable low-energy annulus

Fix an accuracy. Choose a small `delta`, a large separation ratio `R>2`,
a large integer `J`, and a large initial radius `L_0`. Consider the
disjoint amplitude annuli

```math
(L_j,V_j),\qquad L_j=L_0R^j,
\quad V_j=RL_j,\quad j=0,\ldots,J-1.                     \tag{18}
```

Their total energy is at most `Cm^2`; hence one annulus carries at most
`Cm^2/J`. Choose it using only the profile array, not the seed. Delete
its values. Equation (4) changes each log partition by at most
`4tC m^2/sqrt(J)`. The remaining array is separated into light values
at most `L_j` and heavy values at least `V_j`.

Apply (17) for `S` and `T`, with the SAME bounded light profile. The
earlier bounded-profile theorem compares their bounded-light pressures
by `o(m^2)`. It applies uniformly because `L_j<=L_0R^J` is a fixed
constant before `m` grows and there are only finitely many choices.

To make every error arbitrarily small, choose `delta` small and `R`
large, choose `J` large, then choose `L_0` large for (8),(12), and only
then let `m` grow. This proves (3). The annulus index may depend on the
array; all comparison errors are uniform over its finite range.

## 8. Exact consequence for the cycle-sensitive weave

The actual normalized Hadamard row spectra obey
`sum_j v_i(j)^2=m`; deletion of the diagonal coordinate only lowers this.
Thus (1) holds with `C=1` for EVERY physical spin query, including the
coherent column-word spikes. Conditional on the bases and physical spins,
the final ordinary port permutations are independent. The comparison
(3) is uniform over the resulting profile arrays. Consequently it may
be integrated over ANY COMMON profile-array law independent of the seed,
even if those profiles themselves are correlated across rows:

More explicitly, first condition on the coordinate assigned to the physical
macro diagonal at each row, then use independent permutations on the
remaining ports. A diagonal-dependent but seed-independent nonnegative
prefactor may be included in the common profile measure. No seed-dependent
diagonal factor is silently included in this assertion.

```math
\left|\log\mathbb E_{\rm profiles}Z_S
       -\log\mathbb E_{\rm profiles}Z_T\right|=o(m^2).     \tag{19}
```

In particular the complete full-spin, fixed-`t` raw-kernel counting
pressure of the new mixed-orbit weave is leading-order identical for all
bounded-cap macro seeds. This strengthens the earlier finite-type and
uniformly-integrable statements to the entire actual query family.

Its exact `K_5` statistic and the exact finite optimizer contraction
differences remain valid. They distinguish LAWS and finite responses,
but cannot by themselves produce a leading `m^2` gap separating favorable
low-cap seed values under this common fixed-temperature permutation law.

This is NOT an assertion that all such seeds have the same original cap,
that actual parent maxima have been compared, or that every possible
construction must be seed-blind. Seed-dependent profile laws, nonproduct
conditional port placements, or a temperature growing with order are
outside (19). A small-retention regime with a different pressure scale
must also be audited separately. Global old-edge changes and selection
among optimal children remain allowed in the original open problem.

## 9. Zero-temperature minimum port-defect universality

Define the actual minimum over row permutations, for a PRESCRIBED profile
array, by

```math
G_S(u)=\min_{\pi_1,\ldots,\pi_m}
             \sum_{i<j}(u_{ij}-S_{ij}u_{ji})^2.             \tag{20}
```

Under the same per-row energy and bounded-cap seed hypotheses,

```math
\boxed{\ |G_S(u)-G_T(u)|=o(m^2)\ }                       \tag{21}
```

uniformly over ALL profile arrays. This is stronger than taking the
fixed-temperature limit in (3) without controlling its entropy error.

First suppose the profiles are uniformly bounded. Quantize to a common
fixed grid with `L_grid` symbols. The number of distinct row-profile
arrays is at most `L_grid^(m(m-1))`, and every array has equal probability
under the uniform multiset-permutation law. Therefore

```math
G_S\le-\frac1t\log Z_S
 \le G_S+\frac{m(m-1)\log L_{\rm grid}}t.                 \tag{22}
```

For fixed grid and fixed `t`, (3), or just its earlier bounded-profile
version, shows that the soft quantities for `S,T` differ by `o(m^2)`.
Equation (22) gives a limsup for the ground-value difference of at most
`log L_grid/t`. Let `t` tend to infinity after the order limit.
Quantization changes every squared-defect sum by at most
`4V delta_grid m(m-1)` for profiles bounded by `V`; then let the grid
spacing tend to zero. This proves the uniform bounded-profile case of
(21).

For separated heavy and light profiles, the proof above is already
constructive at zero temperature. Start from a minimizing light
permutation. Move tagged zeros to the prescribed heavy ports using at
most `H` swaps, costing at most `16L^2H`. Insert the matched heavy values,
costing at most `delta^2Cm^2/2`. Restore the unmatched values using the
deterministic perturbation bound. The category probability is no longer
needed. Heavy removal is still pointwise favorable. Thus

```math
0\le\frac{G_S(u)-G_S(\ell)}{m^2}
\le4\sqrt{C\eta_{\rm bad}}
       +\tfrac12\delta^2C+16CL^2/V^2.                    \tag{23}
```

The same low-energy annulus, chosen by the profiles alone, reduces the
general case to (23) and the bounded-profile ground theorem. The same
accuracy-first parameter order proves (21).

The row permutations in (20) are allowed to depend on the prescribed
query profile. Thus this result does NOT interchange a single parent
construction choice with its maximum over ALL physical spins. In
particular it is not a comparison of the optimal parent caps.

## 10. Uniform Laplace principle for every diverging temperature

The elementary full permutation-count bound would require
`t_m/log m -> infinity` to compare soft and hard minima at scale `m^2`.
The heavy compiler removes that unnecessary condition. In fact, for
EVERY `t_m -> infinity`,

```math
\sup_{u,S}
 \frac{-t_m^{-1}\log Z_S(u;t_m)-G_S(u)}{m^2}
 \longrightarrow0,                                    \tag{24}
```

where the supremum obeys the same fixed `C,K` hypotheses.

Here are the uniformity details. Divide (17) by `t` before interpreting
its soft quantities as energies. The separated-tail error becomes

```math
4\sqrt{C\eta_{\rm bad}}+
  h_{C,\delta}(V)/t+\tfrac12\delta^2C+16CL^2/V^2.          \tag{25}
```

The corresponding hard comparison is (23). Middle-annulus deletion
has the same `4C/sqrt J` error in both energy units. For the bounded
light profile, fix a finite quantization grid. Its soft-hard gap is
at most `m(m-1)log L_grid/t`, plus the deterministic quantization error.
Choose the annulus, separation, matching and grid accuracies first;
then let `m` and `t_m` grow. The only entropy terms in (25),(22) are
fixed constants divided by `t_m`. This proves (24) uniformly, without
any condition on the rate of divergence.

Combining (21),(24) gives

```math
\frac{|\log Z_S(u;t_m)-\log Z_T(u;t_m)|}{t_m m^2}
 \longrightarrow0.                                    \tag{26}
```

The denominator in (26) is IMPORTANT. It does not prove
`o(m^2)` comparison of unscaled log pressures at a growing temperature.
In particular it does not eliminate a possible critical-temperature
effect in a smaller subleading pressure scale.

## 11. A positive finite-accuracy entropy consequence

The same proof shows that for every fixed `epsilon>0` there is a finite
constant `A_epsilon`, independent of `m`, profiles and allowed seeds,
such that for all sufficiently large `m`,

```math
\Pr_{\rm row\ permutations}
 \{D_S(u)\le G_S(u)+\epsilon m^2\}
 \ge e^{-A_\epsilon m^2}.                              \tag{27}
```

Thus approximating the best port defect to fixed accuracy does not pay
the apparent `m^2 log m` cost of specifying every row permutation.

One can derive (27) directly from the proof's fixed heavy categories and
fixed bounded grid. Alternatively choose a fixed sufficiently large `t`
so that the uniform soft-hard estimate is at most `epsilon m^2/2`.
Writing `p` for the probability in (27), split its partition into that
event and its complement to obtain

```math
e^{-t\epsilon m^2/2}
 \le e^{tG_S}Z_S\le p+e^{-t\epsilon m^2}.
```

For large `m`, this gives
`p>=one-half exp(-t epsilon m^2/2)`, which is (27).
This is a positive approximation-and-entropy theorem for the specified
port problem. It is not a finite response state for the original
simultaneous Boolean maximum.
