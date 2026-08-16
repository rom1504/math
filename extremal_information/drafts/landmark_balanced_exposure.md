# A landmark--exposure dimension law

## Verdict

There is a rigorous cross-benchmark law, but it is a **two-sided
certificate**, not a one-number formula. At scale `epsilon`, query landmarks
give an upper bound and independently exposed queries give a lower bound. For
a response language realizing the whole Lipschitz ball (unit-load Max-Cut),
the exposure dimension is exactly a packing number of the query interface.
For max-plus residuals it is bounded by the number of raw coordinates and is
equal to the number of robustly pinned coordinates when all coordinates are
exposed.

A proposed law of the form

```math
\log \operatorname{Cov}_\epsilon(\mathcal F)
\asymp \operatorname{fat}_{c\epsilon}(\mathcal F)
```

is false already in the nonprojective companion: one query and the one-state
max-plus family `f_v(*)=v`, `v in [-B,B]`, have fat dimension one at every
scale below `B`, but covering entropy `Theta(log(1+B/epsilon))`.  In the
projective language one may instead use the two-query interval
`f_v=(0,v)`.  The missing factor is scalar precision.  The corrected law
below retains it explicitly.

## Exact theorem (landmarks and balanced exposure)

Let `(X,d)` be a finite metric space of diameter `D`. Write

```math
\mathcal L(X)=\operatorname{Lip}_1(X,d)/\mathbb R\mathbf 1,
\qquad
d_{\rm sh}([f],[g])={1\over2}\operatorname{osc}(f-g).
```

Let `N_X(r)` be the least size of an `r`-net and `P_X(s)` the largest
cardinality of a set whose distinct points have distance at least `s`.
For a subclass `F subseteq mathcal L(X)`, define its **balanced exposure
dimension** `E_gamma(F)` to be the largest even `k` for which there are
queries `x_1,...,x_k` and thresholds `a_1,...,a_k` such that, for every
`U subseteq [k]` with `|U|=k/2`, some representative `f_U` of a member of
`F` satisfies

```math
f_U(x_i)\ge a_i+\gamma\quad(i\in U),
\qquad
f_U(x_i)\le a_i-\gamma\quad(i\notin U).          \tag{L.1}
```

This definition is invariant under addition of constants and is the
projective version of fat shattering appropriate to separator profiles.
Then, for `r,eta>0`, with external cover centers,

```math
\boxed{
\log_2\operatorname{Cov}^{\rm ext}_{r+\eta/2}(\mathcal F,d_{\rm sh})
\le N_X(r)\log_2\left(\left\lceil {D\over\eta}\right\rceil+1\right).}
                                                               \tag{L.2}
```

For every subclass `F` and every `0<epsilon<gamma`,

```math
\boxed{
\log_2\operatorname{Cov}_{\epsilon}(\mathcal F,d_{\rm sh})
\ge
\log_2 {E_\gamma(\mathcal F)\choose E_\gamma(\mathcal F)/2}.}
                                                               \tag{L.3}
```

For the complete Lipschitz language the exposure dimension is computed
exactly from the *query* space:

```math
\boxed{
E_\gamma(\mathcal L(X))
=2\left\lfloor {P_X(2\gamma)\over2}\right\rfloor.}           \tag{L.4}
```

For `L>0` and `gamma>0`, the nonprojective companion is

```math
\operatorname{fat}_\gamma(\operatorname{Lip}_L(X))
=P_X(2\gamma/L).                                             \tag{L.5}
```

If a literal response subclass lies in `[-M,M]^X` and is `L`-Lipschitz,
nearest-landmark quantization also gives

```math
\log_2\operatorname{Cov}^{\rm ext}_{Lr+\eta/2}(\mathcal F,\|\cdot\|_\infty)
\le N_X(r)\log_2\left(\left\lceil {2M\over\eta}\right\rceil+1\right).
                                                               \tag{L.6}
```

### Proof

Normalize each Lipschitz `f` by `min f=0`; then `0<=f<=D`. On an `r`-net
`S`, round `f(s)` upward to a multiple `q_s` of `eta` and form

```math
u(x)=\min_{s\in S}\{q_s+d(x,s)\},
\quad
ell(x)=\max_{s\in S}\{q_s-d(x,s)\},
\quad g=(u+ell)/2.
```

Both envelopes are one-Lipschitz. The Lipschitz inequalities and a landmark
within distance `r` give

```math
f\le u<f+\eta+2r,
\qquad f-2r\le ell<f+\eta.
```

Thus `osc(g-f)<2r+eta`, proving (L.2); there are at most
`ceil(D/eta)+1` choices per landmark. Nearest-grid quantization at radius
`eta/2`, followed by nearest-landmark decoding, proves (L.6).
The decoded functions need not lie in the subclass `F`, which is why these
are external covers.  Selecting a member of `F` from every nonempty external
ball gives internal covers at twice the displayed radii.

If (L.1) holds and `U,V` are different balanced subsets, there are
`i in U\V` and `j in V\U`. Hence `f_U-f_V` is at least `2gamma` at `x_i`
and at most `-2gamma` at `x_j`, so their shape distance is at least
`2gamma`. This proves (L.3).

Conversely, if (L.1) holds, choosing balanced patterns which respectively
contain `i` and `j` shows

```math
d(x_i,x_j)\ge
\max\{a_i-a_j+2\gamma,a_j-a_i+2\gamma\}\ge2\gamma.
```

Thus exposed queries are non-strictly `2gamma`-separated.  If `C` is any
such set,
assign `+gamma` on a balanced subset of `C` and `-gamma` on its complement.
These data are one-Lipschitz and extend to `X` by the McShane formula. This
proves (L.4). The same two-pattern necessity and the same extension, without
the balancing restriction, prove (L.5).

## Weighted-automaton specialization

For suffix vectors `H subset R^p`, put

```math
F_v(h)=\max_i(v_i+h_i),
\qquad
g_v(h)=F_v(h)-F_0(h),
\qquad v\in[-B,B]^p.
```

Subtracting the fixed baseline does not change distances between residuals.
Moreover `|g_v|<=B`, and `g_v` is two-Lipschitz for the projective suffix
metric. If `fat_gamma` denotes ordinary fat dimension, then

```math
\boxed{
k_{\rm pin}\le \operatorname{fat}_\gamma\{g_v:v\in[-B,B]^p\}
\le \min\{p,P_H(\gamma)\},
\qquad 0<\gamma<B,}                               \tag{WA.L}
```

where `k_pin` is the number of distinct coordinates having suffixes with
margin `2B` as in (16.62).

The metric upper bound follows from (L.5). The bound by `p` is structural,
not parameter counting: for a query-threshold pair `(h,t)`,

```math
g_v(h)\ge t
\quad\Longleftrightarrow\quad
\exists i:\ v_i\ge t+\max_jh_j-h_i.
```

Thus subgraph labelings are complements of upper orthants in `R^p`. Upper
orthants have VC dimension `p`: if `p+1` points were shattered, excluding
each point while retaining all the others would make each point a unique
coordinate minimum, requiring `p+1` distinct coordinates. Robust pinning
gives `g_v(h^{(j)})=v_{i_j}`, proving the lower bound.

More strongly, the pinned subfamily is an isometric `[-B,B]^{k_pin}`
sup-cube. Consequently

```math
\max\{2,B/\epsilon\}^{k_{\rm pin}}
\le \operatorname{Cov}_\epsilon
\le (1+2B/\epsilon)^p,                            \tag{WA.E}
```

for `0<epsilon<B`, up to universal changes of radius.  The `B/epsilon`
term is the usual volume bound, while the `2^k` term is the corner packing
of the exposed cube; the latter is needed as `epsilon` approaches `B`.
When all `p` coordinates are exposed this is the sharp
`Theta(p log(1+B/epsilon))` law. Query landmarks give the complementary upper
bound (16.61); max-plus covariance improves the generic factor `2r` in
(L.6) to `r`. The affine-line example does not contradict the law: its
projective packing/covering number at the exposure scale grows with `p`, even
though its affine dimension is one.

## Max-Cut/CSP specialization

For unit boundary load, the normalized pure-Max-Cut compiler proves

```math
\{[h_G]\}=\mathcal L(X_w),
\qquad X_w=\{\pm1\}^w/\{s\sim-s\}.
```

Hence (L.4) computes the response exposure dimension exactly. At
`gamma=epsilon w`, `0<epsilon<1/4`, a greedy projective Hamming code gives

```math
E_{\epsilon w}(\mathcal L(X_w))
\ge 2^{(1-H_2(2\epsilon)+o(1))w}.                 \tag{MC.1}
```

Together with (L.2)--(L.3), this recovers the double-exponential cover size
(exponentially many response bits) from a scale-dependent packing of the
separator query interface.  For a strict lower bound on a cover at radius
`epsilon w`, take `gamma=(epsilon+o(1))w` from above; this leaves the entropy
exponent unchanged.

Component size sharply limits this dimension. A fixed `m`-edge Max-Cut or
Boolean-CSP topology is a binary-witness max-affine grammar in `m` shared
weight parameters. There are at most `(3^m-1)/2` distinct optimizer
comparison hyperplanes. Their arrangement has at most

```math
[4((3^m-1)/2+1)]^m
```

faces; on every face the response image modulo constants is linear of rank
at most `m`. Covering each radius-`R` image and taking the union proves

```math
\log_2\operatorname{Cov}_\delta
\le \log_2 T
+m\log_2(4((3^m-1)/2+1))
+m\log_2(1+2R/\delta),                            \tag{MC.2}
```

where `T` is the number of allowed topologies. This is the arrangement proof
of Theorem 16.8, included to make the dimension calculation explicit. For
unit-load `m`-edge components it gives

```math
E_{\epsilon w}(\mathcal C_{w,m})
-\log_2(E_{\epsilon w}(\mathcal C_{w,m})+1)
\le O_\epsilon(m^2+m\log(w+m)).                  \tag{MC.3}
```

Indeed (L.3) at cover radius `epsilon w/2` and
`binom(k,k/2)>=2^k/(k+1)` give (MC.3). The same statement holds for a
bounded weighted Boolean-CSP presentation after inserting its topology
count and response-oscillation bound. Thus the exposure dimension is full
query-packing size for the exponential private compiler, but only polynomial
in the shared presentation size for bounded components.

## What this is, and what it is not

* It is not tautological metric entropy. `N_X(r)` and `P_X(2gamma)` are
  computed in the **query interface before covering responses**; (WA.L) and
  (MC.2) compute realized exposure from automaton coordinates and shared
  optimizer geometry. No covering number is hidden in the definition of
  `E_gamma`.
* It is not classical rate--distortion. The statements are deterministic,
  worst-case, and uniform in the declared queries. A prior and Fano can be
  added afterward, but no expected-loss or mutual-information quantity is
  used here.
* It is more than attaching the phrase “fat shattering.” The substantive
  facts are the exact identity (L.4), the orthant computation (WA.L), the
  isometric exposed cube (WA.E), and the finite-normal arrangement ceiling
  (MC.2). The one-query counterexample shows why a bare fat dimension at one
  scale is insufficient; scalar precision or a genuinely multiscale profile
  is unavoidable.

## Composition verdict: static entropy only

The theorem controls a response class at one declared scale. It does **not**
make a cover into a derivative congruence. A once-decoded full residual can
be passed through a fixed nonexpansive continuation without amplifying its
one error, but rounding back to a finite landmark code after each update can
accumulate error. Exact tropical lumpability (16.64)--(16.67), or an exact
metric-shell semilattice, is separate algebraic information and is what gives
arbitrary-depth closure.

The exact path-kernel family makes the limitation decisive. For `q>=3`, on
`{0,...,q-1}`, let

```math
K_\delta(i,j)=0\ (i=j),
\qquad K_\delta(i,j)=a|i-j|-\delta\ (i\ne j),
\qquad a>\delta.
```

It is within `delta` of the path metric and has one-step idempotence defect
`delta`, but for `1<=T<=q-1`,

```math
K_\delta^{\star T}(i,j)
=a|i-j|-\delta\min\{T,|i-j|\}\quad(i\ne j).
```

Thus `||K_delta^{star T}-K_delta||_infinity=(T-1)delta`. Taking
`delta=c/q` and `T=q-1` makes the local/static error vanish while the
long-depth drift tends to `c`. Consequently no query/fat/exposure dimension
law alone supplies nonaccumulating composition; one must additionally prove
that the chosen state descends under the future semigroup.
