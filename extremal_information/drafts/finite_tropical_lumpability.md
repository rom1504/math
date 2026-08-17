# Finite tropical lumpability, dynamic response memory, and exact benchmarks

**Status.** Proof source for the finite-lumpability campaign.  Statements in
Sections 2--7 have complete proofs below; the exact arithmetic checks are in
[`../experiments/verify_pullback_holonomy.py`](../experiments/verify_pullback_holonomy.py)
and
[`../experiments/verify_switching_benchmarks.py`](../experiments/verify_switching_benchmarks.py),
with strict-strip and weighted-cycle certificates in
[`../experiments/verify_strict_strip_response_certificates.py`](../experiments/verify_strict_strip_response_certificates.py)
and
[`../experiments/verify_weighted_contextual_defect_cycles.py`](../experiments/verify_weighted_contextual_defect_cycles.py).
The purpose is to separate four resources which a single word such as
"state" can hide:

1. an exact symbolic control language;
2. approximation of the numerical response;
3. reusable dynamic memory;
4. the information size of the response carrier.

The abstract refinement theorem is a finite-state
Myhill--Nerode/bisimulation fact.  Likewise, contractive approximate
abstractions and weighted-automaton minimization are classical.  The useful
new content here is their response-theoretic separation, the finite
pullback-normal certificate for selector systems, the sharp strip response
metric, and the counterexamples showing that a small static response cover
need not give any small reusable dynamic quotient.

## 1. Conventions

Let `E` be a finite input alphabet and let `F_e:X->X`, `e in E`, be maps.
The finite observation partition is `P_0`; it may record an optimizer cell,
a tie-consistent selector, or a declared finite output.  For a word
`w=e_1...e_t`, write

```math
F_w=F_(e_t) circ ... circ F_(e_1).
```

A deterministic quotient is **sound and path-realizing** when every quotient
transition is induced by a whole source class.  Thus every labelled path from
a nonempty class is followed by every seed in that class.  This is stronger
than drawing an edge whenever two local cells have a nonempty compatible
intersection.

For max-plus maps use the column convention

```math
(F_Kx)_i=max_j(K_(ij)+x_j).
```

On `R^r/R1`, use `||[x]||_H=osc(x)/2`.

## 2. The intrinsic exact quotient

Define a sequence of equivalence relations by

```math
x equiv_t y
iff
F_wx and F_wy lie in the same P_0 atom for every |w|<=t.       \tag{2.1}
```

Equivalently,

```math
P_(t+1)=P_t wedge bigwedge_(e in E) F_e^(-1)P_t.                \tag{2.2}
```

### Theorem 2.1 (finite contextual lumpability)

There is a finite deterministic forward congruence refining `P_0` if and
only if the refinement sequence (2.2) stabilizes after finitely many steps.
At its first stable stage, `P_t` is the coarsest such congruence.  If a
congruence with `k` classes exists, stabilization occurs after at most
`k-|P_0|` strict refinements.

#### Proof

If `P_(t+1)=P_t`, (2.2) says exactly that `x equiv_t y` implies
`F_ex equiv_t F_ey` for every `e`; hence `P_t` is a forward congruence.
Conversely, a forward congruence `Q` refining `P_0` distinguishes any pair
distinguished after a word, so `Q` refines every `P_t`.  Therefore
`|P_t|<=|Q|`.  Every strict refinement of a finite partition increases its
number of classes, so the sequence stabilizes.  The same induction shows
that every forward congruence refining `P_0` refines the stable partition,
which proves coarseness. `square`

This theorem is a characterization, not yet a compression theorem: computing
`P_t` can enumerate the full orbit language.  The next result gives a finite
dual certificate whose size is controlled before orbit cells are generated.

When the observation itself is robustly query-exposed, refinement growth is
already a quantitative obstruction. Let it take values in a finite subset of
a metric space with pairwise distances at least `Delta`, and let `N_t=|P_t|`.

### Theorem 2.2 (observable refinement-growth sandwich)

For every `0<=epsilon<Delta/2`, the finite-horizon predictor complexity of
Section 5 obeys

```math
N_T<=C_T(epsilon)<=sum_(t=0)^T N_t.                              \tag{2.3}
```

The upper predictor is exact. Consequently, if `N_T` is unbounded, no finite
`epsilon`-predictor works at every depth; if
`liminf T^(-1)log N_T=h>0`, at least `(h/log 2-o(1))T` bits are necessary.
If the refinement stabilizes at `N` classes, one exact `N`-state congruence
works forever.

#### Proof

Different `P_T` atoms have observation trees which differ after some word of
length at most `T`; their values there are `Delta`-separated. A common
predictor state would put them within `2epsilon`, proving the lower bound.

For the upper bound, use layered states `(t,[x]_(P_t))`, `0<=t<=T`. Initialize
at `(T,[x]_(P_T))`; on input `e`, send

```math
(t,[x]_(P_t))->(t-1,[F_ex]_(P_(t-1))).                           \tag{2.4}
```

This is well-defined by (2.1). Decode the current observation and make the
level-zero states absorbing after the declared horizon. The number of states
is the right side of (2.3). Stabilization invokes Theorem 2.1 and removes the
time layer. `square`

This is not automatically a lower bound for numerical response. A branch
label near a tie may not be separated by any declared value query. The
observation must itself belong to the response interface, or be recovered
from it with margin `Delta`. The theorem prices symbolic path existence and
response approximation separately rather than identifying them.

## 3. A finite pullback-normal certificate

Let `V=R^r/R1`, of dimension `d=r-1`, and suppose each `F_e:V->V` is
piecewise affine.  Fix a deterministic tie rule.  On a selector region its
formula is

```math
A_(e,sigma)z=P_sigma z+b_(e,sigma).                         \tag{3.1}
```

Let `H={h_1,...,h_m}` be gauge-invariant affine forms whose sign partition
refines the declared observation partition. Write `C_s` for the
possibly empty relatively open sign face on which
`sgn h_j=s_j in {-1,0,1}`.

### Theorem 3.1 (finite invariant-arrangement certificate)

Assume every declared generator is enabled on the whole sign atom on which it
is used, and:

1. the signs of `H` determine, for each input `e`, a selector formula
   (3.1) which agrees with `F_e` on the whole face, including ties; and
2. for every `h in H` and every selector formula used on a nonempty face,
   `h circ A_(e,sigma)` is zero, a nonzero constant, or a nonzero scalar
   multiple of a member of `H`.

Then the nonempty sign faces form a finite forward congruence:

```math
F_e(C_s) subseteq C_(tau_e(s)).                              \tag{3.2}
```

It is sound and path-realizing.  In particular, every graph-cycle repetition
starting at a nonempty face is a genuine finite trajectory.  The number of
faces is at most

```math
N_H<=sum_(j=0)^d 2^j {m choose j}<=3^m.                     \tag{3.3}
```

The certificate is checked on the finitely many generator branches, not on
all words.  The reverse witness graph of Theorem 16.19 consequently has at
most `N_H r^2` vertices for one carrier and `N_H r^4` vertices for the paired
cross carrier.

#### Proof

On `C_s`, assumption 1 fixes a valid affine formula.  For each `h in H`,
assumption 2 makes the sign of `h(F_ez)` a function only of `s`: it is fixed
for a nonzero constant, zero for the zero form, and otherwise the sign of one
coordinate of `s`, possibly reversed.  This proves (3.2).  Induction proves
that every seed in a source face follows every labelled quotient path from
that face, including repeated cycles.  The standard maximal count for all
relatively open faces of `m` affine hyperplanes in dimension `d` is the first
bound in (3.3); `3^m` is immediate from sign vectors. `square`

The signs need not determine independently chosen tangent selectors on a tie
face.  The theorem is exact for map values because all selected maximizing
coordinates give the same value there.  A claim about joint tangent or
secant paths must instead apply the certificate to the appropriate paired
lift `(x,y)->(F_ex,F_ey)`, or prove a common perturbation cone.  Independent
tie branching is only an upper relaxation and can create false cycles.
Likewise, a raw observable transfers through a block quotient only when it
factors through that quotient.  Global surjectivity of `Lambda_c` does not
control variation inside its fibres: a raw coordinate difference can vary
while all block maxima remain fixed.  Microscopic witness claims therefore
need either quotient-factorable terminal observations or a separate uniform
fibre-oscillation bound.  For control-dependent partial maps, enabling domains
must be included among the forms defining the sign atoms; otherwise use the
edge-domain certificate below.

### Proposition 3.2 (finite dual closure from zero affine holonomy)

For selector-affine branches, start with finitely many coordinate-comparison
forms

```math
h_(ijc)(z)=z_i-z_j-c.
```

Create the finite ordered-pair transport multigraph whose branch edge sends
`(i,j)` to `(sigma(i),sigma(j))` with weight `b_i-b_j`; diagonal pairs are
terminal constants.  If every directed cycle reachable from a starting pair
has total weight zero, the set of pullback hyperplanes reachable from the
starting comparisons is finite.  Its oriented affine closure (retaining the
sign of every proportionality multiplier) is an `H` satisfying
assumption 2 of Theorem 3.1.  More precisely, if `P` is the finite collection
of simple transport paths, there are at most `|H_0||P|` pulled-back
hyperplanes.

If a reachable cycle has nonzero weight and can be iterated syntactically,
its repeated pullbacks give infinitely many distinct parallel hyperplanes.
This obstructs this particular finite-normal certificate, but not every
possible finite quotient: the branch cycle may be dynamically unrealizable
or later cells may saturate.

#### Proof

One pullback is

```math
h_(ijc)(P_sigma z+b)
=z_(sigma(i))-z_(sigma(j))-[c-(b_i-b_j)].                    \tag{3.4}
```

Along a path, the threshold is the original threshold minus the path weight.
Delete cycles from a path.  Zero cycle weights do not change its pulled-back
form, leaving one of finitely many simple paths.  Conversely, repeating a
nonzero cycle which returns to the same ordered pair changes the threshold
by a nonzero multiple on every repetition, so the resulting parallel
hyperplanes are distinct. `square`

Both restrictions in Proposition 3.2 are essential.  Normal directions alone
do not determine affine cells: for `F(x)=x+1` and `H={x}`, the negative cell
contains points whose images have different signs, although the normal is
fixed.  And zero *additive* holonomy is insufficient when pullbacks rescale
normals: `F(x)=2x`, `h=x-1` has zero translation but produces the distinct
hyperplanes `x=2^(-k)`.  Proposition 3.2 works because selector pullback sends
an oriented coordinate-difference normal to another such normal with unit
coefficient.  Theorem 3.1 itself is stated for full affine forms, not merely
their linear directions.

A common invariant normal fan is therefore sufficient but not necessary as
an *initial* description: finitely many pullback refinements may stabilize.
The invariant arrangement at stabilization is the real object.

### Proposition 3.3 (surjective edge-domain alternative)

Whole-cell determinism is not necessary for the finite-orbit realization
used by Theorem 16.19.  Let an edge `e:q->q'` have a domain
`E_e subseteq C_q` on which its displayed affine formula is exact.  Suppose
the outgoing domains cover every actual one-step behavior and

```math
F_e(E_e)=C_(q')                                                     \tag{3.5}
```

for every declared edge.  Then the graph is sound and every finite graph path
is realizable.  Indeed, choose any terminal point and lift it backward
edge-by-edge using (3.5).  The seed may depend on the path length, exactly as
allowed in the finite-orbit pumping conclusion of Theorem 16.19.

## 4. Two dual max-plus block quotients

For a partition `Pi={I_1,...,I_r}` and a gauge `c`, put

```math
J_c(u)_i=u_a-c_i  (i in I_a),
Lambda_c(x)_a=max_(i in I_a)(x_i+c_i),
Lambda_c J_c=id.                                                   \tag{4.1}
```

The following distinction prevents a common false lumping argument.

### Proposition 4.1 (embedded versus aggregate congruence)

If

```math
max_(j in I_b)(K_(ij)-c_j+c_i)=S_K(a,b)  (i in I_a),               \tag{4.2}
```

then `F_KJ_c=J_cF_(S_K)`: the block-constant embedded suffix space is
invariant.  If instead

```math
max_(i in I_a)(K_(ij)+c_i-c_j)=S~_K(a,b)  (j in I_b),              \tag{4.3}
```

then `Lambda_cF_K=F_(S~_K)Lambda_c` on every raw state.  Both are simple
regroupings of maxima.  Neither condition implies the other.  Under the
row-vector convention, (4.3) is the transpose of (16.64).

For a fixed gauge and initial coloring, repeated splitting by the block-max
signatures in (4.2), or by the transposed signatures for (4.3), produces the
coarsest stable refinement.  This is a finite partition-refinement algorithm;
optimizing over gauges is a separate feasibility problem and can give
incomparable partitions.

If the coefficients in the appropriate independence condition differ from
their block value by at most `epsilon`, the associated one-step coordinate
sup error, hence its Hilbert error, is at most `epsilon`.

The aggregate law gives a finite approximate state without constructing any
optimizer-cell language.

### Theorem 4.2 (approximate block lumpability with switching reset)

Let `K_e` be a finite alphabet of finite max-plus matrices, and suppose a fixed
partition, gauge, and `r by r` matrices `S_e` satisfy

```math
|max_(i in I_a)(K_e(i,j)+c_i-c_j)-S_e(a,b)|<=epsilon
quad(j in I_b) .                                                   \tag{4.4}
```

Let `Y subset R^r/R1` contain all actual aggregate states under consideration
and be forward invariant under every `F_(S_e)`. Assume every quotient map is
defined and nonexpansive on `Y`, and every contiguous legal length-`L`
quotient composition is `rho`-contractive on `Y` in Hilbert norm, `rho<1`.
If `Y`
lies in a Hilbert ball of radius `R`, then for every `h>0` there is a finite
deterministic response state with at most

```math
(1+2R/h)^(r-1)                                                     \tag{4.5}
```

states whose aggregate shadow error, at every depth and under every legal switch
word, is at most

```math
h+{L(epsilon+h) over1-rho}.                                      \tag{4.6}
```

If a raw response is within `kappa` of an `L_O`-Lipschitz function of
`Lambda_c(x)`, its response error is at most

```math
kappa+L_O[h+{L(epsilon+h) over1-rho}].                            \tag{4.7}
```

Thus the prequantized quotient has dimension `r-1`; after quantization the
reusable simulator has the finite state count (4.5), with no active-cell
enumeration in the defect test. The block-defect certificate beyond the given
raw matrices stores `O(|E|pr+|E|r^2+p)` coefficients for the defects, quotient
matrices, gauge, and partition. The realized transducer additionally stores
the net centers, transitions, and decoder data; certifying every legal
length-`L` contraction can itself be expensive unless an intrinsic reset or
scrambling certificate is available. At `epsilon=0`, before net quantization,
(4.3) is an exact semiconjugacy for arbitrary raw states.

#### Proof

Regrouping maxima gives

```math
Lambda_c(F_(K_e)x)_a
=max_b max_(j in I_b){x_j+c_j+L^e_(a j)},                         \tag{4.8}
```

where `L^e_(a j)` is the coefficient on the left of (4.4).  Replacing every
such coefficient by `S_e(a,b)` changes each output coordinate by at most
`epsilon`, proving the one-step Hilbert defect. Take an internal `h`-net `C`
of `Y` with deterministic nearest map `Q`, initialize
`c_0=Q(Lambda_cx_0)`, and set
`c_t=Q(F_(S_(e_t))c_(t-1))`. The fresh
residual is at most `epsilon+h`.  Group its transported secants into
length-`L` blocks and sum the geometric series, as in Theorem 16.18; the
initial net error contributes at most `h`.  The standard volumetric net bound
in the `(r-1)`-dimensional Hilbert norm gives (4.5), and the response decoder
gives (4.7). `square`

The finite input path is exact, but (4.6) is metric shadowing rather than an
exact claim about tie-selector itineraries.  Without contraction, the same
construction incurs `h+T(epsilon+h)`.  The repeated-rounding strip example
in Section 7 shows that this linear alternative can be sharp.

The contraction is deliberately restricted to `Y`. Globally on the whole
projective space, a nonconstant finite max-plus map has Hilbert Lipschitz
coefficient one; a strict global coefficient is a projective reset. Legal-word
recognition, if not the full shift, also requires its own finite control state
in addition to the net centers counted in (4.5).

## 5. Static response entropy times dynamic forgetting

Let `(Y,d)` be a forward-invariant metric quotient, with maps `G_e:Y->Y`.
Suppose every `G_e` is nonexpansive and every contiguous legal length-`L`
composition occurring in a trajectory is `rho`-Lipschitz, `rho<1`. Let
`pi:X->Y` have semiconjugacy defect

```math
d(pi F_e x,G_e pi x)<=zeta.                                      \tag{5.1}
```

Assume the response `R:X->l_infinity(Theta)` admits an `L_R`-Lipschitz
decoder `Psi:Y->l_infinity(Theta)` with

```math
||R(x)-Psi(pi x)||_infinity<=kappa.                               \tag{5.2}
```

### Theorem 5.1 (depth-uniform response net under block contraction)

Let `C` be an internal `eta`-net of `Y`, and define finite transitions by a
deterministic nearest-net map

```math
c -> Q(G_e c).
```

Starting with `d(pi x_0,c_0)<=eta`, for every input word and every depth,

```math
d(pi x_t,c_t)
<=eta+{L(zeta+eta) over 1-rho},                                  \tag{5.3}
```

and therefore

```math
||R(x_t)-Psi(c_t)||_infinity
<=kappa+L_R[eta+{L(zeta+eta) over 1-rho}].                        \tag{5.4}
```

The reusable state count is at most `Cov_eta(Y)`, independently of depth.
Without block contraction, the same proof gives only
`eta+t(zeta+eta)`.

#### Proof

At each step, (5.1), nonexpansiveness, and quantization contribute a fresh
residual at most `zeta+eta`.  Group the transported residuals into blocks of
length `L`.  Each older block is attenuated by another factor `rho`, so their
sum is at most `L(zeta+eta)sum_(j>=0)rho^j`; the initial residual contributes
at most `eta`.  Apply (5.2) and Lipschitz continuity of `Psi`. `square`

The theorem is an approximate-bisimulation theorem written in response
coordinates.  Its nonclassical content is not the contraction argument; it
is the accounting identity it suggests:

```math
reusable response bits
<=log_2 Cov_eta(Y),
```

at distortion (5.4).  Static response entropy becomes dynamically valid
only after a congruence or a forgetting estimate controls the reuse error.
The abstract architecture is close to incremental-stability symbolic models;
see [Girard--Pola--Tabuada](https://arxiv.org/abs/0807.5022).

The preceding estimate has an intrinsic rate--distortion form. Let `R:X->Z`
be an `L_R`-Lipschitz response into a metric space, with `L_R>0`, and define
the horizon-`T` contextual
pseudometric

```math
d_T(x,y)=max_(|w|<=T)d_Z(R(F_wx),R(F_wy)).                       \tag{5.5}
```

An `epsilon`-predictor is a finite set `S`, an encoder `q:X->S`, deterministic
updates `delta_e:S->S`, and a decoder `g:S->Z` such that

```math
d_Z(R(F_wx),g(delta_wq(x)))<=epsilon                             \tag{5.6}
```

for all declared words of length at most `T`.  This is a word-consistent
simulator from an initially encoded state; it does not require the exact
congruence identity `qF_e=delta_eq`.  Write `C_T(epsilon)` for its smallest
number of states.  For maps on a metric space, define the suffix
memory gain

```math
G_T=sup_(t<=T,e_1...e_t)
 [Lip(F_(e_t)...F_(e_1))
  +sum_(s=1)^t Lip(F_(e_t)...F_(e_(s+1)))],                      \tag{5.7}
```

where the empty suffix has Lipschitz constant one.

### Theorem 5.2 (response packing--memory sandwich)

Then

```math
Pack_(>2epsilon)(X,d_T)
<=C_T(epsilon)
<=Cov^int_(epsilon/(L_R G_T))(X).                               \tag{5.8}
```

The upper bound uses internal covers; the lower bound is valid for every
external predictor.  In particular:

```math
G_T<=T+1                                                        \tag{5.9}
```

for nonexpansive maps. If, in addition, every contiguous legal length-`L`
block occurring in a trajectory is `rho`-Lipschitz, `rho<1`, then

```math
G_infinity<={L over1-rho}.                                      \tag{5.10}
```

Thus the correct general interaction is not literally a product of static
state count and mixing time.  Dynamic memory changes the **resolution** at
which static metric entropy must be paid.  If
`log Cov_eta(X) asymp d log(1/eta)`, it contributes an additive
`d log G_T` bits; an expanding inverse such as Theorem 6.1 can instead make
the contextual packing itself exponential in depth.

#### Proof

If `q(x)=q(y)`, the predictor supplies the same decoded answer after every
word, so (5.6) and the triangle inequality give `d_T(x,y)<=2epsilon`.
This proves the packing lower bound.  For the upper bound, choose an
`eta`-net `C`, encode by a nearest center, and after each exact map quantize
again to a nearest center.  The initial and every fresh quantization residual
have norm at most `eta`.  At time `t` they are transported by exactly the
full product and the suffix products in (5.7).  Hence state error is at most
`eta G_T`, and response error is at most `L_R eta G_T`.  Set
`eta=epsilon/(L_R G_T)`.  Summing unit suffix bounds proves (5.9); grouping
suffixes into length-`L` blocks proves (5.10). `square`

Equation (5.8) is useful only when the internal metric itself is a strict
response-sufficient quotient.  Replacing `X` by the full raw landscape would
be tautological.  Exact synchronization, an invariant carrier, or a proved
response isometry such as Theorem 7.1 is what licenses the smaller space.

There is a complementary law for the external response tree.  It does not by
itself give transition updates, but it measures exactly how contraction and
input branching turn local static response entropy into future-query entropy.
For a response `h:X->Z`, set

```math
S_h(r,epsilon)=sup_(A subseteq X, diam A<=r)
 Cov^ext_epsilon(h(A)).                                           \tag{5.11}
```

Let `q=|E|`, `diam X<=D`, and let every `F_e` be `rho`-Lipschitz for
`0<rho<1`.  Define the full response-tree map

```math
Phi_T(x)=(h(F_wx))_(|w|<=T)                                     \tag{5.12}
```

with the sup metric over words.

### Theorem 5.3 (contraction-weighted context-tree transform)

For every horizon,

```math
log Cov^ext_epsilon(Phi_T(X))
<=sum_(k=0)^T q^k log S_h(D rho^k,epsilon).                       \tag{5.13}
```

If `epsilon>0` and `h:X->R^p` is `L`-Lipschitz in the sup norm, with finite
`L,D` and `q>=1`, put

```math
n_k=max{1,ceil(LD rho^k/(2epsilon))}.
```

Then

```math
Cov^ext_epsilon(Phi_T(X))
<=prod_(k=0)^T n_k^(p q^k).                                     \tag{5.14}
```

Every factor and exponent in (5.14) is sharp over a horizon-dependent family
of finite affine systems. For each `T`, take
`X=prod_(|u|<=T)[0,D]^p`, use the sup metric, let `h(x)=Lx_empty`, and define

```math
(F_a x)_u=cases(rho x_(au),&|u|<T;
                0,&|u|=T)                                      \tag{5.15}
```

when `au` is in the finite tree, filling the boundary by zero.  Then

```math
Phi_T(X)=prod_(k=0)^T[0,LD rho^k]^(p q^k),                       \tag{5.16}
```

whose exact external sup-covering number is the right-hand side of (5.14).
This proves distribution-free sharpness; it is not an asymptotic statement
about one fixed finite-dimensional system.

Writing `A=LD/(2epsilon)>1` and
`m=ceil(log A/log(1/rho))`, the transform stops after scale `m`.  For `q>=2`
and `T>=m`, the sharp example has

```math
log Cov^ext_epsilon(Phi_T(X))
=Theta_(p,q,rho)(A^(log q/log(1/rho))).                           \tag{5.17}
```

For `q=1`, it is `Theta_(p,rho)(m^2)`.  Thus contraction makes the
infinite-depth response tree finitely compressible at each fixed accuracy,
but branching can turn logarithmic mixing depth into a power of
`1/epsilon`.

#### Proof

For `|w|=k`, the set `F_wX` has diameter at most `D rho^k`; cover its response
image using (5.11).  Taking the Cartesian product of these external covers
over the `q^k` words at every level proves (5.13).  A set of diameter `r` has
response image inside a `p`-box of side at most `Lr`, proving (5.14).

In (5.15), the different context-tree coordinates are independent and the
word `a_1...a_k` reads
`h(F_(a_k)...F_(a_1)x)=L rho^k x_(a_1...a_k)`.  A rectangular box in the sup
metric has covering number equal to the product of the one-dimensional
covering numbers: the regular grid proves the upper bound, and a Cartesian
product of points separated by more than `2epsilon` proves the lower bound.
Indeed, for an interval of length `ell` and
`n=ceil(ell/(2epsilon))>=2`, the `n` equally spaced points including both
endpoints have spacing `ell/(n-1)>2epsilon`, while `n` equal subintervals
give the matching cover.  This proves (5.16)
and sharpness.  Only levels `k<m` have `n_k>1`.  For `q>=2`, the last active
level gives the lower bound `(log 2)pq^(m-1)`, while summing
`log ceil(A rho^k)` gives `O_(p,q,rho)(q^m)`.  Since
`q^(m-1)<A^(log q/log(1/rho))<=q^m`, (5.17) follows.  For `q=1`, summing the
linear sequence `log A-k log(1/rho)` gives quadratic order in `m`. `square`

Every predictor supplies an external response-tree center, so

```math
Cov^ext_epsilon(Phi_T(X))<=C_T(epsilon).                          \tag{5.18}
```

The reverse can fail because arbitrary external table centers need not be
shift-consistent.  Theorem 5.2 adds word-consistent deterministic updates
through the declared horizon; it still does not impose exact semiconjugacy.
The Cantor system below has a tiny
one-step response image but an expanding inverse which makes its contextual
packing exponential; it therefore separates the two laws sharply.

The finite-horizon notion nevertheless has an exact compactness law.  Let
`C_infinity(epsilon)` be the least number of states in one predictor satisfying
(5.6) for every finite word.

### Theorem 5.4 (finite predictive compactness)

If the response space `Z` is compact, then

```math
\boxed{C_infinity(epsilon)=sup_(T>=0) C_T(epsilon).}              \tag{5.19}
```

In particular, uniformly bounded finite-horizon response memory cannot evade
one finite predictor by changing its transition graph at every horizon.  If
no finite infinite-depth predictor exists, `C_T(epsilon)` must diverge.

#### Proof

One infinite predictor restricts to every horizon. Conversely, choose a
strictly increasing sequence of horizons with `C_T(epsilon)<=S`; after
relabelling, pad predictors with unused self-loop states and one fixed decoder
value to use the common set `[S]`. There are only finitely many families of
deterministic transition maps `delta_e:[S]->[S]`, so along a subsequence of
horizons they are one fixed family.  Compactness of `Z^S` gives a further
subsequence on which every decoder value `g_T(s)` converges to `g(s)`.

On this final subsequence, for each `x`, at least one state `s_x` occurs
infinitely often among its
initial encodings `q_T(x)`.  Define `q_infinity(x)=s_x`.  Fix any finite word
`w` and pass along that infinite occurrence subsequence.  Eventually
`T>=|w|`, so

```math
d_Z(R(F_wx),g_T(delta_ws_x))<=epsilon.
```

Taking the decoder limit proves the same inequality with `g`.  Since `x` and
`w` were arbitrary, this is one `S`-state infinite-depth predictor. `square`

This is compactness of predictive behavior, not realizability of a quotient
partition.  The selected state `q_infinity(x)` need not satisfy
`q_infinity(F_ex)=delta_e q_infinity(x)`.  Requiring that stronger identity
is the exact lumpability problem of Theorems 2.1 and 3.1.

## 6. Static compression can fail dynamically

### Theorem 6.1 (compact encoder--decoder memory explosion)

On `X=[0,1]`, let `h(x)=x`,

```math
E_0(x)=x/3,
E_1(x)=(x+2)/3,
```

and let the continuous piecewise-affine decoder be

```math
R(x)=cases(3x,&0<=x<=1/3;
           2-3x,&1/3<=x<=2/3;
           3x-2,&2/3<=x<=1).
                                                                    \tag{6.1}
```

Then `R E_b=id`.  After `t` encoders, the `2^t` states reached from zero are
pairwise `1/3`-separated by one of the future queries
`h,h circ R,...,h circ R^(t-1)`.  Consequently every reusable deterministic
summary which answers all those future queries to error `epsilon<1/6` has at
least `2^t` states.  For total encode-plus-query depth `T=2t-1`, it needs at
least `(T+1)/2` bits.

In contrast, the static output has an `epsilon`-cover of size
`ceil(1/(2epsilon))`, and even the numerical one-step response vector

```math
(h,h circ E_0,h circ E_1,h circ R)
```

has an `epsilon`-cover of size at most `ceil(3/(2epsilon))+1`.  Thus no bound
depending only on a finite static response cover can control reusable dynamic
memory, even for five affine pieces on a compact invariant interval.

#### Proof

Writing the newest ternary digit first, an encoded state is

```math
x_b=sum_(r=1)^t 2b_(t-r+1)3^(-r).
```

For two different words, let `r` be the first reversed-digit position where
they differ.  Applying `R^(r-1)` strips their common prefix.  One remaining
state lies in `[0,1/3]`, the other in `[2/3,1]`, so their observed values
differ by at least `1/3`.  Two states in one summary fibre would have decoded
responses at distance at most `2epsilon`, proving the packing.  The one-step
response vector is `3`-Lipschitz, proving its scalar-grid cover. `square`

### Proposition 6.2 (Helly-sharp tie incompatibility)

Modulo common constants, let cyclic coordinates satisfy `x_(m+1)=x_1`,
choose `sum_i c_i=1`, and set

```math
phi_i(x)=max{x_(i+1),x_i+c_i}.
```

Every proper subfamily of the `m` tie faces intersects, but their total
intersection is empty.  More sharply, if

```math
r_i=x_(i+1)-x_i-c_i,
```

then

```math
min_x max_i |r_i|=1/m.                                           \tag{6.2}
```

Indeed, deleting one edge leaves a forest on which the tie equations can be
integrated.  On the full cycle, `sum_i r_i=-1`, proving the lower bound;
setting every `r_i=-1/m` attains it.  Hence compatibility of every
`m-1` local faces does not imply a common exact trajectory, and failure can
be postponed to arbitrary dimension.  The case `m=3` is the three tie lines
of `max(u,0)`, `max(v,u)`, and `max(0,v+1)`, with robust gap `1/3`.

### Proposition 6.3 (one-step observations are not a congruence)

On `R^3`, observe `h(x)=x_1` and let `A=P_(12)`, `B=P_(23)`.  The states
`0` and `e_3` agree under the empty word, `A`, and `B`, but

```math
h(ABe_3)=1 !=0=h(AB0).                                           \tag{6.3}
```

Thus even every individual-mode response can miss a switching-created
feature.  Dimension three is minimal for this linear two-switch pattern: in
dimension two, `ker h` is one-dimensional, so if `v,Av,Bv in ker h`, then
`Bv=lambda v` and `ABv=lambda Av in ker h`.

An irrational circle rotation supplies a complementary nonexpanding example:
its pullback boundaries are `-k alpha mod 1`, all distinct.  The contextual
refinement never stabilizes.  Its Sturmian branch language is non-sofic, so
no finite graph can be both sound and finite-path-realizing.

## 7. Width-three Ising as a switching benchmark

Let `X={-1,1}^w`.  A strip prefix has boundary profile

```math
f_P(x)=max{E_P(sigma): sigma_boundary=x}.
```

A new column updates it by

```math
(T_cf)(y)=A_c(y)+max_x{f(x)+sum_i J_i^c x_i y_i}.                 \tag{7.1}
```

This state was derived from contextual response before comparison with the
classical transfer matrix.

### Theorem 7.1 (sharp strip response state and metric)

Under arbitrary real strip continuations, the exact contextual quotient is
the complete `2^w`-entry table `f`.  More precisely,

```math
D_ctx(f,g)
:=sup_C |Opt_f(C)-Opt_g(C)|
=||f-g||_infinity.                                                \tag{7.2}
```

Modulo additive constants, the metric is
`osc(f-g)/2`.  Therefore rounding a current profile once to sup error
`epsilon` preserves every future optimum, at every depth, to exactly
`epsilon`.

#### Proof

The upper bound is max-plus nonexpansiveness.  To expose coordinate `a`, add
one column with `J_i=M>0`, fields `h_i=La_i`, no vertical edges, and `L>M`.
For fixed old `x`, the new spin is uniquely `y=a`.  If

```math
2M>max_(x!=a) max{f(x)-f(a),g(x)-g(a)},                            \tag{7.3}
```

the old boundary is uniquely `x=a` for both profiles.  The difference of
continued optima is therefore `f(a)-g(a)`.  Choose `a` attaining the sup
norm. `square`

This proves an exact minimality statement, not merely the familiar transfer
update.  Width two realizes every projective table by Walsh expansion.  A
strict width-three strip also has a full-dimensional reachable response
neighborhood: the explicit integer certificate recorded in the benchmark
verification has determinant `-1024`, inverse sup-operator norm `2`, and
optimizer margin at least `2`.  Its fixed affine chart realizes the closed
sup-cube of radius `1/2` (with ties permitted on the boundary), and the
radius-`1/8` cube has unique predecessors with margin at least `3/2`.  It
yields eight absolute and seven projective continuous response degrees.

### Proposition 7.2 (finite restricted strip quotient)

For a finite rational column alphabet on a fixed width, normalize a profile
as `(b,p)`, where `b=max f` and `max p=0`.  If all coefficients and the seed
lie in `eta Z`, every normalized successor lies in

```math
S_B={p in (eta Z)^(2^w): max p=0, min p>=-B},                    \tag{7.4}
```

where

```math
B=max_c max_(y,z)[A_c(y)-A_c(z)
                  +2 sum_(i:y_i!=z_i)|J_i^c|].                  \tag{7.5}
```

Thus, with `K=floor(B/eta)`, there are at most

```math
(K+1)^(2^w)-K^(2^w)                                             \tag{7.6}
```

normalized control states.  Partition refinement by each column's pair
`(baseline toll, successor block)` gives the coarsest exact weighted quotient.

The bound follows because the difference of two maxima lies between the
minimum and maximum of the pointwise kernel difference.  Integrality and
normalization then give (7.4)--(7.6).  This is a genuine finite
path-realizing weighted state, with an exact unbounded scalar register kept
separate.

### Proposition 7.3 (fresh quantization can drift linearly)

Embed width one in width two and repeat a column with horizontal coupling
`J=-K`, field `h=s`, and `0<s<K`.  From a flat message, the exact optimum
after `n` columns is

```math
V_n=nK+(n mod 2)s.                                                \tag{7.7}
```

A scheme which normalizes and rounds the shape to a grid of mesh `Delta`,
where `s<Delta/4`, rounds every raw shape `(-2s,0)` back to flat and predicts
`n(K+s)`.  Its error is `ns` for even `n` and `(n-1)s` for odd `n`.
Thus the one-shot isometry (7.2) does not justify repeated projection onto a
static response net.  A congruence, telescoping gauge, or recurrent
contraction is genuinely additional information.

More generally, put the exact and approximate controls in one finite product
graph and label each edge by the difference of their scalar baseline tolls.
If a reachable repeatable cycle has nonzero total label `eta`, repeating it
produces absolute error `k eta`.  Uniform absolute accuracy therefore forces
zero discrepancy on every reachable cycle, equivalently a scalar coboundary
on each reachable strongly connected component.  In (7.7), both projective
controls return exactly to flat after two steps, but the reward-cycle
discrepancy is `2s`.  Projective lumpability can therefore coexist with
extensive error in the absolute optimum; the additive register and its
cocycle are a separate resource.

## 8. Weighted automata: exact refinement and a sharp defect cycle

The weighted-automaton benchmark uses the row action

```math
(vT_e)_j=max_i(v_i+T_e(i,j)).                                    \tag{8.1}
```

For blocks `I_a` and fixed gauges `c_i`, define

```math
R_e(i,b)=max_(j in I_b){T_e(i,j)+c_j}-c_i.                       \tag{8.2}
```

Starting from equal terminal values `beta_i-c_i`, repeatedly split each
source block by all signatures `(R_e(i,B))_(e,B)`. The stable result is the
coarsest strong block-max quotient for this fixed gauge. Any other stable
refinement refines every stage, because a current target block is a union of
its target blocks and maxima regroup over that union. This derivation used
contextual signatures before comparison with weighted-automaton minimization.

For an approximate quotient `S_e`, put

```math
epsilon_e=max_(a,b,i in I_a)|R_e(i,b)-S_e(a,b)|.                 \tag{8.3}
```

The aggregate `pi(v)_a=max_(i in I_a)(v_i+c_i)` satisfies

```math
||pi(vT_e)-pi(v)S_e||_infinity<=epsilon_e,                       \tag{8.4}
```

and hence a word `e_1...e_n` has error at most
`sum_t epsilon_(e_t)`. If a microscopic maximizing cycle and its quotient
block cycle are both repeatable, their per-cycle discrepancy `D` has
telescoping gauges and produces error exactly `kD` after `k` repetitions.
On a finite joint control graph, bounded scalar error is therefore equivalent
to zero discrepancy on every repeatable cycle, or to a coboundary on every
reachable strongly connected component.

The exact four-state verifier discovers a two-block quotient and checks 4,802
aggregation identities. Perturbing one maximizing microscopic self-loop by
`delta` makes exact refinement split all four states. The old quotient still
has one-step defect exactly `delta`, while the perturbed loop gives error
`n delta` at depth `n`. Exact quotient cardinality is discontinuous, but
finite-depth response distortion is continuous and the general sum bound is
sharp.

Strong weighted partition refinement is classical; the benchmark-level
addition is the exact quantitative relation between its failed equality, the
finite-horizon response modulus, and the repeatable discrepancy cocycle.

## 9. What the campaign establishes

The finite dynamic law is now precise:

```math
usable future compression
= finite response cover
+ path-realizing congruence or quantitative forgetting.          \tag{8.1}
```

Theorem 3.1 supplies a checkable exact congruence; Theorem 5.1 supplies a
depth-independent approximate substitute; Theorem 6.1 proves that static
entropy alone cannot replace either.  Theorem 16.19 then turns an exact
presentation into a finite pumpable-cycle decision, while the strip
benchmark proves that exact contextual minimality and approximate dynamic
reuse are different questions even in the simplest switching model.

This is a dynamic/compositional component, not just static response
geometry.  It still leaves a sharp converse open: characterize when failure
of pullback stabilization forces either observable cycle drift or horizon-
growing contextual entropy, rather than merely a larger but bounded exact
quotient.

The boundary with existing theory is important. Weighted block refinement is
classical; compare
[Lombardy--Sakarovitch](https://arxiv.org/abs/2112.09387). Finite approximate
symbolic models under incremental stability and finite exact bisimulations
from contractive polyhedral Lyapunov slices are also established mechanisms;
compare
[Girard--Pola--Tabuada](https://arxiv.org/abs/0807.5022) and
[Ding--Lazar--Belta](https://arxiv.org/abs/1203.6408). The project-level
theorems are the response-entropy/memory bounds, their exact extremizers, the
selector pullback-holonomy certificate, and the benchmark isometries and
falsifiers which place those classical mechanisms in one information law.

A completely general effective converse is not credible: rational
piecewise-affine reachability is undecidable already in dimension two.
[Varonka--Watanabe](https://arxiv.org/abs/2502.19923) recover decidability for
two-dimensional Bellman operators, illustrating why the next converse should
stay inside selector/Bellman subclasses with explicit enabling geometry.
