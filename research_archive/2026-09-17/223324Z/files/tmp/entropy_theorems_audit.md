# Adversarial audit: landmark/exposure and binary max-affine entropy

## Verdict

The mathematical core now passes. The balanced projective identity has the
advertised constants, and facet deletion really gives shape radius `1/2` and
pairwise distance `1`. The amended radius truncation also fixes the earlier
major defect in the max-affine entropy comparison: without it, positive
homogeneity makes the response cone unbounded and its fixed-scale cover
infinite.

Remaining corrections are mostly convention, proof-completeness, and scope:
the weighted-automaton sharp rate needs an explicit epsilon range and a
corner packing; strict/non-strict packing terminology should be uniform; and
the `exp(Theta(m^2))` common-fan assertion needs a classical threshold-function
citation. The advertised “exact remaining problem” is not an exhaustive
dichotomy.

## 1. Landmark and balanced exposure

### Constants and covers

* **(L.2) passes as an external cover.** Upward rounding gives
  `g-f in [-r,r+eta)`, hence `d_sh([g],[f])<r+eta/2`. There are at most
  `ceil(D/eta)+1` values per landmark. For a proper subclass the constructed
  centre need not belong to the subclass, so `Cov^ext` is essential.
* **(L.6) also passes only externally.** Nearest-grid error `eta/2` plus
  nearest-landmark error `Lr` gives `Lr+eta/2`, with
  `ceil(2M/eta)+1` symbols. The current `ext` qualifier is correct.
* Selecting a class member from each nonempty external ball gives internal
  covers at radii `2r+eta` and `2Lr+eta`. For the full Lipschitz ball, the
  (L.2) envelope centre is itself one-Lipschitz, so no doubling is needed.
* **(L.3) is valid even for external covers.** The exposed family is
  non-strictly `2gamma`-separated, while an external epsilon-ball has
  diameter at most `2epsilon<2gamma`. Writing `Cov^ext` would record the
  strongest statement. The strict hypothesis `epsilon<gamma` is necessary.

### Exact identity

Under the weak-margin convention and the draft's non-strict query packing
convention,

```math
E_gamma(Lip_1(X)/R1)=2 floor(P_X(2gamma)/2)
```

is correct. Necessity gives, for every exposed pair `i,j`,

```math
d(x_i,x_j)>=max(a_i-a_j+2gamma,a_j-a_i+2gamma)>=2gamma.
```

Conversely, `+gamma/-gamma` labels on a constant-weight subset of a
non-strict `2gamma`-packing are one-Lipschitz and extend by McShane. Define
`E_gamma=0` explicitly when no positive even size is feasible. At line 138,
“`2gamma`-packing” should be “non-strictly `2gamma`-separated” (or a
`rho`-packing for every `rho<2gamma`) if manuscript `Pack_r` means strict
separation.

(L.5) is correct for `L>0`, `gamma>0`; state those hypotheses. The
fat-shattering/metric-packing connection for Lipschitz functions is
classical (Gottlieb--Kontorovich--Krauthgamer, *IEEE TIT* 60 (2014),
Theorem 1 and its zero-threshold lemma); McShane supplies the matching lower
bound here. “Balanced exposure dimension” should be presented as a custom
projectively invariant, constant-weight analogue, not an established
standard notion.

The opening example is now correctly labelled nonprojective. For its
suggested projective replacement, state the geometry: take `X={0,1}` with
`d(0,1)=B` and `f_v=(0,v)`, `v in [-B,B]`. Then
`d_sh([f_v],[f_u])=|v-u|/2`, the entropy is
`Theta(log(1+B/epsilon))` for `0<epsilon<B`, but `E_gamma=2` for every
`0<gamma<=B/2`.

## 2. Weighted automata

The structural claim

```math
k_pin <= fat_gamma({g_v}) <= min(p,P_H(gamma)),  0<gamma<B,
```

passes. The function `g_v` is 2-Lipschitz in the projective suffix metric,
giving the packing scale `gamma`. Thresholding transforms a query-threshold
pair into a point of `R^p`; the subgraph class is the complement of an upper
orthant, whose VC dimension is exactly `p`. Margin `2B` pins the chosen
coordinate, including endpoint ties, and gives `g_v(h^(j))=v_(i_j)`.

(WA.E)'s volume bounds are correct, but alone do not prove the following
sharp-rate sentence uniformly as `epsilon` approaches `B`. Add the range and
corner packing:

```math
0<epsilon<B,
max{(B/epsilon)^k_pin,2^k_pin}
 <= Cov^ext_epsilon <= Cov_epsilon <= (1+2B/epsilon)^p.
```

The `2^k_pin` corners are pairwise distance `2B`. Thus, when `k_pin=p`,

```math
log Cov_epsilon=Theta(p log(1+B/epsilon)),  0<epsilon<B.
```

Without the range restriction this is false: at `epsilon>=B` the cube has a
one-centre cover.

## 3. Max-Cut/CSP corollary

(MC.1) is correct. A projective Hamming code of distance `2epsilon w` has
exponent `1-H_2(2epsilon)`, and forcing even cardinality is negligible. To
lower-bound a cover at radius exactly `epsilon w`, take
`gamma_w>epsilon w` with `gamma_w/w -> epsilon`; the current “from above”
sentence handles this strictness correctly.

For `k=E_(epsilon w)(C_(w,m))`, (L.3) at radius `epsilon w/2` and Theorem
16.8 give

```math
k-log_2(k+1) <= C_epsilon [m^2+m log(w+m)].
```

Hence `k=O_epsilon(m^2+m log(w+m))`. This verifies (MC.3); `<=` is clearer
than `=O(...)`. The constants `(3^m-1)/2`,
`[4((3^m-1)/2+1)]^m`, and `1+2R/delta` agree with Theorem 16.8. For one fixed
topology the `log T` term is absent; it enters only after union over `T`
topologies.

## 4. Facet deletion and generic perturbation

The construction and all metric constants pass. With `V_F=V cap F`,
`A_0=V`, `A_F=V setminus V_F`, and
`gamma_F=h_V(u_F)-h_(V setminus V_F)(u_F)`, one gets

```math
f_(theta_F)=c_F 1-e_F,
||[f_(theta_F)]||_sh=1/2,
d_sh([f_(theta_F)],[f_(theta_G)])=1.
```

Distinct facets are incomparable, so `F setminus G` contains a vertex; this
justifies every off-diagonal coordinate. The same argument works for an
antichain of nonempty proper faces.

The tie perturbation passes. The complement of finitely many comparison
hyperplanes is dense, support functions are continuous, and positive
rescaling preserves tie-freeness because the hyperplanes are central. The
perturbed shapes can therefore be restored to norm exactly `1/2` and kept at
distance `>1-2eta`. This proves response-metric robustness and tie-freeness,
not a dimension-uniform parameter-space distance from the tie arrangement;
the draft now makes that distinction.

With strict `Pack_r`, the exact targets are not a “1-packing”: they are
non-strictly 1-separated and are a `rho`-packing for every `rho<1`. The
external-cover lower bound for `delta<1/2` is correct. The endpoint is truly
excluded: the zero shape lies at distance `1/2` from every exact target.

The new radius truncation in `C_m(R,delta)` is essential and correct. It
would be cleaner for (7) to display the truncated set explicitly and for
(8) to use the sharper `C_m(1/2,delta)` instead of `C_m(1,delta)`, but the
current inequality is true.

## 5. Sources and overclaims

The cited 0/1-polytope results support the exact asymptotic scale:

* Gatzouras--Giannopoulos--Markoulakis prove a full-dimensional example with
  at least `(c m/(log m)^2)^(m/2)` facets, yielding
  `(m/2)(log m-2 log log m-O(1))` response bits.
* Fleiner--Kaibel--Rote, Corollary 8, give
  `f_(m-1)<=C(m-2)!` and
  `f_i<=C(m-2)![2(i+1)]^(m(m-1)/(m+1))` for `i<m-1`. This is uniformly
  `exp(O(m log m))`; summing over `i` controls the whole face lattice. It
  supports only the stated single-polytope, one-shape-per-face mechanism,
  not arbitrary common refinements.

Two claims still need adjustment:

1. Robust-draft line 197 says a common normal fan *can have*
   `exp(Theta(m^2))` cells. This is true but non-elementary. Ternary normals
   contain `(1,x)`, `x in {+-1}^{m-1}`; those chambers encode Boolean linear
   threshold functions, whose logarithmic count is asymptotic to `(m-1)^2`
   by Zuev's classical threshold-function theorem. Add that citation and
   reduction, or weaken the sentence to the elementary `exp(O(m^2))` upper
   bound.
2. Lines 202--208 are not an “exact” dichotomy. The optimum could be strictly
   intermediate between `m log m` and `m^2`; neither endpoint theorem need
   hold. Call it an open gap and present the two alternatives only as
   possible endpoint-resolving routes.

Finally, landmark quantization, Lipschitz fat-shattering, orthant VC
dimension, volume covering, finite hyperplane arrangements, and threshold
chamber enumeration are classical ingredients. Absent a dedicated novelty
search, the drafts should claim the benchmark-specific synthesis and
consequences, not external novelty for those ingredients.
