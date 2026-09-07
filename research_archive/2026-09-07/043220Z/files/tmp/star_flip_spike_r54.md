# Wave 54 main audit: star-flip minimality gives a migration/entropy dichotomy

## Scope and outcome

This note asks whether exact discrete edge-sign minimality can bound a large
local field of a child ground, as required by (10.1290).  It proves an exact
covering identity and a quantitative dichotomy for every edge block.  Applied
to the positive star of a spiky child ground, the result gives either many
near-parent-ground witnesses or one near-parent-ground witness that
anti-aligns with roughly half the star.

This is a genuine bridge from a fixed child spike to parent near-ground
geometry, but it does **not** prove field regularity or `A^2` cancellation.
The witnesses may migrate with the perturbation, their restrictions need not
be child grounds, and absolute state count is not saved project incidence.
Exact `A8`, `A9`, and `A10` data take the anti-aligned branch very strongly.

## 1. Exact block-cover theorem

Let `A` be an exact order-`n` minimizer, `q=Q(A)`, and use oriented full
states `omega=(sigma,x)`.  Put

```math
E_A(\omega)=\sigma x^{\mathsf T}Ax,
\qquad
\Delta_\omega=q-E_A(\omega),
\qquad
s_e(\omega)=\sigma a_ex_ix_j.
```

Fix an arbitrary edge block `P`, `|P|=h`, and set

```math
C_\omega(P)=\{e\in P:s_e(\omega)=-1\},
\qquad t_\omega=|C_\omega(P)|.
```

For every `F subseteq P`, flipping exactly the signing edges in `F` and
using exact global minimality gives

```math
0\le Q(A^F)-q
=\max_\omega\left\{-\Delta_\omega-4\sum_{e\in F}s_e(\omega)\right\}.
```

Consequently, for every `k`-set `F subseteq P`, some `omega` satisfies

```math
\boxed{
\Delta_\omega\le4k,
\qquad
|F\cap C_\omega(P)|\ge \frac{k}{2}+\frac{\Delta_\omega}{8}.}
\tag{R54.1}
```

This gives the following exact covering count.  Define

```math
N_{h,k}(t,\Delta)
=\sum_{j\ge\lceil k/2+\Delta/8\rceil}
\binom tj\binom{h-t}{k-j}.
```

Then

```math
\boxed{
\binom hk
\le
\sum_{\omega:\Delta_\omega\le4k}
N_{h,k}(t_\omega,\Delta_\omega).}
\tag{R54.2}
```

The proof is only the union bound over the exact certificates (R54.1), so
there is no hidden interchange of `min`, `max`, or expectation.

For `0<eta<1/2`, Hoeffding's upper-tail estimate for sampling without
replacement implies the rigorous dichotomy

```math
\boxed{
\begin{array}{l}
\text{either some }\omega\text{ with }\Delta_\omega\le4k\text{ has}
t_\omega>(1/2-\eta)h,\\[1mm]
\text{or there are at least }\exp(2\eta^2k)
\text{ oriented states with }\Delta_\omega\le4k.
\end{array}}
\tag{R54.3}
```

Indeed, under `t_omega/h <= 1/2-eta`, its normalized contribution to
(R54.2) is at most `exp(-2 eta^2 k)`; the additional positive deficit in
the threshold only decreases it.

## 2. Application to a child-field spike

Let `S` be an `m`-selector, orient a child ground `y` by `sigma_y`, and put

```math
r_i=\sigma_yy_i(A[S]y)_i\ge0.
```

For a fixed `i in S`, take the positive internal star

```math
P_i=\{\{i,j\}:j\in S\setminus\{i\},\ 
\sigma_ya_{ij}y_iy_j=1\}.
```

Its size records the spike exactly:

```math
\boxed{|P_i|=\frac{m-1+r_i}{2}.}
\tag{R54.4}
```

Thus (R54.1)--(R54.3) apply to every spiky child coordinate.  With
`k=H=n^(3/4-c)` and fixed selector density, they yield either
`exp(Omega(H))` parent states within deficit `4H` or a parent near-ground
state whose signed star disagrees with almost half of `P_i`.

Neither conclusion presently pays (10.1290):

1. `exp(Omega(H))` is only an absolute count inside a state space of size
   `exp(Theta(n))`; it can still have exponentially unsaved normalized mass.
2. The anti-aligned state may depend on the flipped subset, and its
   restriction to `S` need not be a child ground or even a favorable label.
3. Equations (R54.1)--(R54.4) contain neither
   `d=(A^2)[T,S]y` nor a negative direction of `H_T`; hence they give no
   completion cancellation.
4. The size `|P_i|` is already `Theta(n)` at fixed density even when `r_i`
   is small.  The theorem detects the positive star, but does not extract
   the subleading excess `r_i` at scale `n^(3/4-c)`.

The last point is the sharp obstacle to using this exact cover as a spike
bound.  A continuation would need a *biased* perturbation comparison that
cancels the baseline `(m-1)/2` and retains the excess `r_i/2`, together with
a common-witness or favorable-restriction theorem.  Plain simultaneous
edge-flip minimality has neither feature.

## 3. Exact finite diagnostic

`tmp/star_flip_spike_r54_check.py` exhausts all oriented parent states and
all subsets of the selected positive star.  It independently verifies the
simultaneous-flip certificates.  On the spikiest audited child grounds:

| instance | `m` | `r_max` | `|P_i|` | full-star best `(Delta,sum_P s)` |
|---:|---:|---:|---:|---:|
| `A8` | 5 | 4 | 4 | `(0,-2)` |
| `A9` | 6 | 5 | 5 | `(8,-5)` |
| `A10` | 6 | 5 | 5 | `(4,-5)` |

For `A9` and `A10`, a single near-ground state is negative on every edge
of `P_i`; for `A8` the optimum is already an exact parent ground and is
negative on three of four edges.  These are exact finite migration warnings,
not asymptotic counterexamples.

## Research judgment

The higher-order star-flip idea produces the useful exact theorem
(R54.1)--(R54.3), but closes as a direct route to field regularity.  It
should be retained only if coupled to a new baseline-canceling perturbation
or a theorem transferring the parent witnesses into the same favorable
project family.  Reusing the unweighted global edge-cube cover alone would
repeat the old migration obstruction.
