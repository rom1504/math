# A uniform matching-scale bound for sparsely active balanced bulk

2026-09-07. **Positive actual-construction phase bound; independent
construction-agent reconstruction PASS, including the primary-backed
conditioned-RIP extension in Section 6.** Arbitrary bias scales and coherent spectra are permitted.
The conclusion concerns a restricted set of physical spin words, not the
whole parent maximum or original convergence.

## 1. Deterministic theorem, including every macro seed

There are `m` fibres of size `q`, with `N=mq`. At fibre `i`, let `T_i`
have balanced columns in `R^q`. Assume that every active column has
squared norm at most `q`, and distinct active columns satisfy

```
|<T_i(:,j),T_i(:,k)>|<=mu.                              (1)
```

Self columns and symmetrically masked columns may be zero. For ANY
fixed symmetric outer signs `S_ij`, form the hollow bulk

```
W_ij=S_ij T_i(:,j) T_j(:,i)^T  (i!=j),    W_ii=0.
```

Let `P_i=I-J/q`. If a physical spin word `x` has exactly `ell`
nonconstant fibres, then `H_W(x)=0` when `ell<=1`, and for `ell>=2`,

```
|H_W(x)| <= [q+(ell-2)mu]/2 * sum_i ||P_i x_i||^2.       (2)
```

This is simultaneous over EVERY such word, EVERY outer seed, and EVERY
symmetric mask satisfying the assumptions. It requires no local
stability, random signs, biased-slice entropy bound, or restriction on
the ratios of the positive minority densities.

Proof. Let `I` be the set of nonconstant fibres and put `y_i=P_i x_i`.
Balanced columns imply `W=P_all W P_all`, and every response of a
constant fibre is zero. Define the directed-port vector
`h_i(j)=T_i(:,j)^T y_i` for distinct `i,j` in `I`, with masked entries
zero. The signed reciprocal swap

```
(Rh)_i(j)=S_ij h_j(i)
```

has operator norm at most one. Thus

```
|H_W(x)|=|h^T R h|/2<=sum_i ||h_i||^2/2.
```

Each `h_i` uses at most `ell-1` columns. The Gram matrix of those
columns has diagonal at most `q` and off-diagonal magnitude at most
`mu`, so Gershgorin gives operator norm at most `q+(ell-2)mu`.
Consequently `||h_i||^2<=[q+(ell-2)mu]||y_i||^2`; summing proves (2).
With at most one nonconstant fibre no surviving macro edge exists.

The coefficient at `ell=2` is sharp. If an unmasked macro edge `ij`
has two full balanced sign columns, choose `x_i=T_i(:,j)` and
`x_j=S_ij T_j(:,i)`, and make every other fibre constant. Then both
active words are balanced,

```
H_W(x)=q^2,        sum_i||P_i x_i||^2=2q.
```

Their variance quotient is EXACTLY `q/(2sqrt(N))=sqrt(q/m)/2`.
This is a sharp relative-phase witness, not a leading original-cap
witness: its energy is only order `N` when `q` is comparable to `m`.
In particular it must not be identified with a homogeneous
near-constant-band lower bound.

## 2. The actual stratified selector has the required coherence

Use the fixed-`L` selector of
`principle_invent_2026_09_07_stratified_marked_selector.md`:
ambient order `m=LM`, physical size `q=kM`, `p=k/L<1`, and exactly
`k` coordinates chosen independently and uniformly in each of `M`
groups of size `L`. Condition on the full Hadamard transform and its
partition. The partition may depend on that transform as prescribed
in the construction; no independence between them is needed here.

Two distinct full Hadamard columns have product sum zero. Their
selected inner product is a sum of independent group contributions,
each bounded in absolute value by `L`, and its expectation is `p`
times the full inner product, hence zero. Hoeffding therefore gives

```
Pr(|selected column inner product|>u)
                 <=2 exp[-u^2/(2mL)].                 (3)
```

The same estimate holds for the selected sum of any nonconstant full
column. There are at most `m^3` column pairs across all `m` fibres.
For fixed `L`, a union bound at `u=C_L sqrt(m log m)` therefore gives,
with probability `1-o(1)`, simultaneously:

```
all selected nonconstant column sums are O_L(sqrt(m log m));
all distinct selected column inner products are O_L(sqrt(m log m)).
```

Assume `q` is even, as in exact balanced repair. Balancing each
nonconstant sign column flips exactly half the absolute column sum,
so at most `O_L(sqrt(m log m))` entries are changed in each column.
Whatever allowable locations the repair chooses, the inner product of
two columns changes by at most twice the total number of their flips.
Thus the EXACT repaired sign columns have squared norms `q` and

```
mu=O_L(sqrt(m log m)).                                (4)
```

Deleting the marked constant self-port and any subsequent symmetric
macro mask preserves this bound on the remaining columns. The random
repair law's finer dependence or operator estimate is not needed for
(4). These conclusions hold conditional on every full transform array
and its prescribed partitions, uniformly; the union bound can therefore
be averaged over the recursive transform law.

## 3. Actual strict-subhalf heterogeneous phase

Combining (2)--(4), the actual repaired stratified bulk satisfies,
with probability `1-o(1)`, simultaneously over every outer seed and
all words with

```
ell=o(sqrt(m/log m))
```

the RELATIVE bound

```
|H_W(x)| <= [sqrt(p)/2+o(1)] sqrt(N)
                                  *sum_i ||P_i x_i||^2. (5)
```

The error is uniform over that chosen active-fibre range and multiplies
the actual variance; it is not only an additive `o(N^(3/2))` statement.
Every nonconstant fibre may have an arbitrary allowed bias. In
particular highly coherent spectra, balanced fibres, nearconstant
fibres, and arbitrarily many minority-density scales can coexist.

For a fixed target `b>sqrt(p)/2`, (2) more explicitly controls all
`ell` for which

```
(ell-2)mu <= (2b-sqrt(p))sqrt(N).                      (6)
```

When `p<1` one may choose `b<1/2`, so this is a strict-subhalf phase
of the ACTUAL bulk. Unlike a pointwise seed statement obtained by
absorbing signs into independent output gauges, (2) is deterministic
and genuinely holds for every seed at once.

## 4. What this adds, and what remains

The bounded-ratio packet theorem controls many active fibres but keeps
their nonzero densities in a common finite-ratio band. The present
coherence theorem controls a sparse set of active fibres with NO bias
or spectral-tail restriction. The companion local-stability theorem
in `principle_synthesis_2026_09_07_balanced_bulk_sparse_stability.md`
handles many heterogeneous fibres when enough minority mass has
controlled incoming tails.

These are complementary actual constructions/estimates, not a complete
cover: profiles with many active fibres, widely separated density
scales, and coherent incoming mass remain possible. Separate phase
bounds cannot be added without paying their mixed cross interaction.
No selected-child transfer, complete variance-ratio theorem, or
convergence claim follows yet.

## 5. Finite exact checks

`computations/principle_synthesis_2026_09_07_sparse_active_check.py`
builds actual one-hole selectors at `L=4`, ambient orders `m=8,16,32`,
and exact individually balanced sign columns. It checks 1,200 instances
of (2), with the exact computed column coherence and rational variance,
over `ell=1,2,3,4`. The two-fibre witness attains energy `q^2` and
variance `2q` exactly at all three orders. All checks PASS; the
corresponding `_check_results.json` preserves the output. The small
orders are not a numerical test of asymptotic coherence decay.

## 6. Primary-backed extension to `o(m/log^4 m)` active fibres

There is a stronger actual phase theorem. Fix `L` and `k<L`, and
integers `s=s_m>=2` with `s log^4(m)/m ->0`. A conditioned version
of the top-level stratified selector and the same operator-controlled
balanced repair satisfy (5) uniformly for ALL words with at most
`s+1` nonconstant fibres and EVERY seed/mask. The conditioning costs
`o(m)` in total log density across all fibres, preserving the existing
fixed-accuracy leading certificates. This is an explicitly changed
selector law, not an unproved all-fibre high-probability assertion for
the original unconditioned law.

### 6.1 Exact primary input

We use Lemma 3.6 of Rudelson--Vershynin, *On sparse reconstruction
from Fourier and Gaussian measurements*, in the
[author-hosted primary paper](https://public.websites.umich.edu/~rudelson/papers/convex-relaxation.pdf).
The lemma and its complete proof on printed pages 7--11 were read
directly. For deterministic `z_1,...,z_M in R^m`, `M<=m`, with
`||z_g||infty<=1`, it bounds the independent Rademacher average by

```
E_epsilon ||sum_g epsilon_g z_g z_g^T||_(s)
 <=a ||sum_g z_g z_g^T||_(s)^(1/2),
a=C sqrt(s) log(2s) sqrt(log(2m) log(2M)).               (7)
```

Here `||B||_(s)=max_(|J|<=s)||B_(J,J)||op`. Enlarging logarithms
avoids small-order conventions. The lemma imposes NO distribution,
identical-law, orthogonality, or isotropy assumption on the rows.
We do not import its Bernoulli-selector corollary for our different law.

Exact opening hypothesis, Lemma 3.6 on printed page 7:
“Let x_1,...,x_k, k<=n, be vectors in C^n with uniformly bounded entries”.
Its displayed equation (10) is the imported bound; the following
hypothesis is `||x_i||infty<=K`. Our use has `k=M`, `n=m`, `r=s`, `K=1`.

### 6.2 Independent nonidentical missing rows

First let `k=L-1`. Condition on ANY full sign Hadamard `H` and
ANY partition into `M=m/L` physical groups of size `L`. The missing
row `X_g` is chosen uniformly in each group, independently between
groups. Its distribution can depend on `g`, but EXACTLY

```
||X_g||infty=1,
sum_g E X_g X_g^T=(1/L)H^T H=M I_m.                   (8)
```

Set `D=E||sum_g X_g X_g^T-M I||_(s)`. Independent copies and Jensen
give `D<=E||sum_g(X_g X_g^T-X'_g X'_g^T)||_(s)`.
The independent differences are symmetric, so inserting independent
Rademacher signs preserves their law. The triangle inequality then
bounds this by `2E||sum_g epsilon_g X_g X_g^T||_(s)`.
No identical-distribution step is used. Apply (7) conditional on the
rows, and then Jensen and the triangle inequality:

```
D<=2a E sqrt(||sum_g X_g X_g^T||_(s))
 <=2a sqrt(D+M),
D<=4a^2+2a sqrt(M).                                   (9)
```

Since `a^2=O(s log^4 m)=o(m)`, the last bound is `o(m)`, uniformly
over the transform and partition. Put
`u_m=(4a^2+2a sqrt(M))/m ->0`. Markov gives, with probability at
least `1-sqrt(u_m)`,

```
||sum_g X_g X_g^T-M I||_(s)<=sqrt(u_m)m.
```

The retained matrix has Gram `H_T^T H_T=mI-sum_g X_g X_g^T`, hence

```
sup_(|J|<=s)||H_T(:,J)||op^2<=q+sqrt(u_m)m.             (10)
```

For general fixed `k<L`, generate a uniform ordering inside each group
and discard its first `L-k` slots. Each slot separately consists of
independent uniform missing rows across groups, so (8)--(10) apply.
The slots are dependent, but summing their restricted Gram deviations
and taking a union over their FIXED number needs no slot independence.
Thus (10), with an `o(m)` error and probability `1-o(1)`, holds for
every fixed `k<L`.

### 6.3 Paid top-selector conditioning and exact repair

Choose `delta_m->0` so that the event `G_m` comprising (10) with
upper error `delta_m m` and the column-excess bounds in Section 2 has
probability at least `1-delta_m`, uniformly conditional on every full
transform and partition. At EACH TOP-LEVEL physical fibre, keep the
original transform law and sample the selector conditional on `G_m`.
Different fibres are still independent. Relative to the original joint
top-fibre law, the density is at most `1/(1-delta_m)`, and across the
whole array it is at most

```
(1-delta_m)^(-m)=exp(o(m)).                            (11)
```

Every nonnegative row or whole-array expectation therefore increases
by at most the corresponding explicit factor. Existing strict leading
`m^2`-scale certificates keep their margins. In the rare-packet theorem,
bias and recursive depth are fixed before ambient order, so its strict
row margin also survives. We do NOT infer a certificate uniform under
a simultaneously vanishing bias from this order of limits.

The event is invariant under spectral-column signed permutations,
since restricted Gram norm and absolute column excess are invariant.
Child matrices need NOT remain independent after conditioning: the
density comparison is for the WHOLE top-fibre law and uses the old
unconditioned certificate. The conditioning is not inserted afresh at
every recursive child node.

Apply the exact balanced-column repair from
`principle_invent_2026_09_07_balanced_transform_compiler.md`. Its matrix
Bernstein bound is uniform conditional on each raw matrix and its
column-excess bound. Choose the probability constant large enough to
hold simultaneously in all `m` fibres. Then, with probability `1-o(1)`,

```
max_i ||T_i-P_i H_(i,T)||op=epsilon_m sqrt(m),
epsilon_m->0.                                        (12)
```

The marked constant self-column may be omitted: it projects to zero.
Because `P_i` is a contraction, (10),(12) imply uniformly in `|J|<=s`,

```
||T_i(:,J)||op<=sqrt(q+delta_m m)+epsilon_m sqrt(m),
||T_i(:,J)||op^2<=q+o(m).                              (13)
```

The OPERATOR-controlled repair is essential here. Unlike Section 2,
we do not assert (13) for arbitrary locations of all allowed flips.

### 6.4 Actual relative phase bound and remaining scope

With at most `s+1` active fibres, every response uses at most `s`
columns. Replace Gershgorin by (13) in the proof of (2), obtaining

```
|H_W(x)|<=[q+o(m)]/2 *sum_i||P_i x_i||^2
 =[sqrt(p)/2+o(1)]sqrt(N)*sum_i||P_i x_i||^2.           (14)
```

This bound is simultaneous in every seed, mask, bias scale, and
coherent profile in the stated active-fibre range. Its error is
RELATIVE even at microscopic total variance. The two-active witness
shows its leading coefficient is sharp in this phase. The conditioned
law retains the previously paid balanced-face and fixed-band bounds,
but (14) does not handle macroscopically many active fibres or pay
their mixed cross interaction. Original convergence is still open.

### 6.5 Quantitative power saving (director derivation, independent audit PASS)

For fixed L,k put theta_m=s log^4(m)/m and
e_m=m^(-1/4)log^(5/4)m. The exact operator-controlled repair has
normalized operator error O(e_m). Instead of near-one conditioning,
condition each selector on restricted deviation at most four times the
UNIFORM deterministic upper bound B_m=4a^2+2a sqrt(M), and on the
column-excess event. Markov and the excess estimate give probability
at least1/2 for all sufficiently large m. The density cost is at most
2 per fibre, hence exp(O(m)) globally, still subleading for all previous
fixed-accuracy leading certificates. For multiple discarded slots use
their fixed finite union, enlarging the constant in B_m as needed.

The restricted Gram error divided by m is
O(sqrt(theta_m)+theta_m). Squaring the repaired singular-value bound
adds O(e_m+e_m^2). Consequently, when theta_m tends to zero, uniformly
all the same words/seeds/masks,

    |H_W(x)| <= [sqrt(p)/2+O(sqrt(theta_m)+e_m)] sqrt(N)
                                      *sum_i||P_i x_i||^2.

In particular s=floor(m^(1-eta)), fixed0<eta<1, gives coefficient error
O(m^(-eta/2)log^2(m)+m^(-1/4)log^(5/4)(m)). This is a quantitative
RELATIVE sparse-phase theorem for a modified actual selector law.
It is not a complete-parent defect estimate or an original M_n recurrence.
