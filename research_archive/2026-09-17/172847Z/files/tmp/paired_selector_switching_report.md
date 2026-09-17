# Paired-selector switching dynamics: an exact witness-cycle theorem and its boundary

## Executive result

There is a clean finite theorem for max-plus systems with switching cells and ties, provided the switching dynamics admits an **exact regular paired-selector presentation**. The required state is not a full trajectory or a suffix-set lift: it is one finite control state together with one ordered coordinate witness.

For a presentation with control set `Q` and `d` affine-selector coordinates, a weighted graph on at most

\[
|Q|d^2
\]

vertices decides uniform boundedness of every declared directed response. A positive relevant cycle is exactly a pumpable linear-drift certificate. In the symmetric/Hilbert case, boundedness is equivalent to zero weight on every relevant cycle, equivalently to a coboundary condition on each relevant strongly connected component.

For two `r`-coordinate max-plus systems with different selectors, the correct closed joint state is the `r^2`-coordinate cross-difference array `w_(u,p)=y_u-x_p`; its ordered-witness graph therefore has at most `|Q|r^4` vertices. The smaller `|Q|r^2` lift applies when an `r`-coordinate error recursion already closes (for example under common selectors), not to a generic different-selector Hilbert comparison.

This theorem handles selector changes, different left/right selectors, and ties through **jointly realizable** resolved edges. It does **not** justify replacing the exact presentation by the usual local face-adjacency graph: a two-dimensional all-finite clamp gives a concrete spurious-cycle counterexample. Moreover, no unrestricted algorithmic theorem of this kind is possible for arbitrary sparse tropical rational systems, by a reduction from Krob's undecidable equality problem.

The useful next target is therefore not a larger generic state. It is a structural theorem guaranteeing an exact finite Markov presentation for a natural subclass (for example an invariant finite normal fan, a genuinely lumpable tropical quotient, or a bounded finite projective orbit).

## 1. Exact regular paired-selector presentations

Let `Q=(V,E)` be a finite directed multigraph with a set `I` of allowed initial vertices. Fix a coordinate dimension `d >= 2`. Each edge

\[
e:q\longrightarrow q'
\]

carries an affine coordinate-selector map

\[
A_e z=P_{\sigma_e}z+b_e,
\qquad (P_{\sigma_e}z)_i=z_{\sigma_e(i)},
\]

where `sigma_e:[d]->[d]` need not be injective and `b_e` is a real vector. Parallel edges are allowed; in particular, distinct resolutions of a max-plus tie may be represented separately.

Each terminal control state `q` has a finite set of ordered observation pairs

\[
O_q\subseteq[d]^2.
\]

For a forward path `p=e_1...e_t` from an allowed start to `q`, put

\[
A_p=A_{e_t}\circ\cdots\circ A_{e_1}
\]

and define its directed output from `z_0` by

\[
D_p(z_0)=\max_{(i,j)\in O_q}
\big((A_pz_0)_i-(A_pz_0)_j\big).
\]

Assume the allowed initial set has bounded oscillation:

\[
\sup_{z_0}\max_{i,j}|(z_0)_i-(z_0)_j|\le R_0<\infty.
\]

The presentation is **sound** if every actual resolved orbit of the system determines a graph path and obeys its stated affine selector updates. It is **path-realizing** if every relevant finite graph path, including every finite repetition of a relevant cycle with fixed entry and exit paths, is realized by some allowed initial state. An **exact regular paired-selector presentation** is both sound and path-realizing. (If one wants a single infinite pumping orbit rather than a sequence of finite witnesses, add the stronger nested-cylinder realization condition; it is automatic in the invariant-cell criterion of Section 7.)

The path-realizing clause is substantive. It is the finite-memory/Markov property needed to concatenate local selector choices. Merely knowing that every individual edge is locally feasible is not enough.

## 2. The reverse ordered-witness graph

Define a finite weighted graph `G^<->` with vertices

\[
(q,i,j),\qquad q\in V,\quad i,j\in[d].
\]

For every presentation edge `e:q->q'` and every ordered pair `(i,j)`, insert the reverse edge

\[
(q',i,j)\longrightarrow(q,\sigma_e(i),\sigma_e(j))
\]

with weight

\[
\omega_e(i,j)=b_e(i)-b_e(j).
\]

A lifted vertex or edge is **relevant** when it lies on a reverse path from a declared terminal observation `(q,i,j)`, `(i,j) in O_q`, to an allowed initial control state. Thus irrelevant selector cells and coordinate witnesses can be discarded by ordinary reachability and co-reachability.

The key identity is one-dimensional even though the dynamics is not:

\[
(A_ez)_i-(A_ez)_j
=z_{\sigma_e(i)}-z_{\sigma_e(j)}+b_e(i)-b_e(j).
\tag{1}
\]

Iterating (1) backward along a path shows that every terminal coordinate difference equals an initial coordinate difference plus the weight of its unique lifted reverse witness path.

## 3. Witness-cycle theorem

### Theorem (finite paired-selector drift dichotomy)

For a sound regular paired-selector presentation:

1. If the relevant part of `G^<->` has no positive-weight directed cycle, then all declared directed outputs are uniformly bounded above, independently of path length.
2. More quantitatively, the bound is `R_0+K`, where `K` is the maximum weight of a simple relevant reverse path (with the harmless finite entry/exit convention absorbed into `K`). A Bellman-Ford potential computes such a bound.

If the presentation is also path-realizing, the converse holds:

3. If a relevant lifted cycle has weight `c>0`, there are fixed forward entry, cycle, and exit words `u,v,w` such that, for every `k`,

\[
\sup_{z_0\ {m allowed}}D_{uv^kw}(z_0)\ge kc-O(1).
\]

Thus positive cycle weight is precisely a pumpable linear-drift certificate for the worst-case response. Under nested-cylinder realization, the same conclusion holds along one infinite orbit.

Consequently, for an exact regular presentation, uniform upper boundedness is equivalent to absence of positive relevant cycles.

If the observation family and relevant graph are closed under pair reversal `(i,j)<->(j,i)`, then the corresponding two-sided/Hilbert response is uniformly bounded if and only if every relevant cycle has weight zero. On every relevant strongly connected component this is equivalent to the existence of a vertex potential `h` with

\[
\omega(u\to v)=h(v)-h(u)
\]

for every lifted edge in the component.

### Proof

For a forward path and terminal witness `(i,j)`, repeatedly apply (1). This gives

\[
(A_pz_0)_i-(A_pz_0)_j
= (z_0)_{i_0}-(z_0)_{j_0}+W(\widehat p),
\tag{2}
\]

where `widehat p` is the associated reverse witness path and `W` is the sum of its edge weights.

If every relevant cycle has nonpositive weight, delete closed subwalks from any relevant reverse path. Removing a nonpositive cycle cannot decrease its weight, so the remaining simple path has weight at least that of the original path. There are finitely many simple relevant paths, hence their weights have a finite maximum `K`. Equation (2) and the initial oscillation bound give `D_p(z_0)<=R_0+K`.

Conversely, let a relevant reverse cycle have weight `c>0`. Relevance supplies fixed reverse entry and exit paths. Reading the entire lifted path backward gives a forward control path containing the corresponding cycle word. Path realization allows the same cycle word to be repeated arbitrarily often without losing feasibility. Equation (2) then gives a fixed endpoint contribution plus `kc`, proving linear drift.

Under pair reversal every edge weight is negated. Hence any negative relevant cycle yields a positive reversed cycle. Two-sided boundedness is therefore equivalent to zero weight on every relevant cycle. Finally, the standard path-independence argument on a strongly connected graph says that a weight function has zero integral around every directed cycle exactly when it is a coboundary: fix a base vertex and define `h(v)` as the weight of any path from the base to `v`; zero cycle weights make this independent of the chosen path. QED.

### Complexity and information content

The lifted graph has at most `|V|d^2` vertices and `|E|d^2` edges. Its state is a finite control cell plus one ordered coordinate witness. **Conditional on `Q` itself being a genuine sub-landscape quotient**, this is smaller than storing trajectories, reachable subsets, or the suffix-set lift exponential in the number of coordinate maps. The theorem does not by itself prove that `Q` is small; that is exactly the lumpability obligation isolated below. It does identify what the additional response memory must retain: the coordinate pair currently certifying the directed response, not all possible future witnesses simultaneously.

## 4. Application to paired max-plus systems with switching and ties

On a polyhedral selector face, a max-plus map is affine with a coordinate selector. Let two `r`-coordinate systems have the resolved joint update

\[
x'=P_\sigma x+a,
\qquad
y'=P_\tau y+c.
\tag{4}
\]

When `sigma` and `tau` differ, the ordinary error `y-x` does not close. The exact smaller-than-trajectory joint state is instead the cross-difference array

\[
w_{u,p}=y_u-x_p,
\qquad (u,p)\in[r]^2.
\]

Equation (4) gives

\[
w'_{u,p}=w_{\tau(u),\sigma(p)}+c_u-a_p.
\tag{5}
\]

Thus `w` is itself an affine coordinate-selector system of dimension `d=r^2`. Its ordered coordinate witnesses are

\[
w_{u,p}-w_{v,q}
=(y_u-y_v)-(x_p-x_q),
\]

and propagate jointly by

\[
(u,v;p,q)
\longmapsto
(\tau(u),\tau(v);\sigma(p),\sigma(q))
\]

with edge weight

\[
(c_u-c_v)-(a_p-a_q).
\tag{6}
\]

In particular,

\[
2d_H(x,y)
=\max_{i,j}\big[w_{i,i}-w_{j,j}\big],
\]

so the terminal observation set consists of `(i,j;i,j)`. The generic paired witness graph has at most `|Q|r^4` vertices. This four-index lift is essential when selectors differ: it retains cancellation between the two channels before taking the absolute/projective maximum. When `sigma=tau` along every transition, the diagonal cross differences close as the ordinary `r`-coordinate error, reducing the lift to `|Q|r^2`.

The initial-oscillation hypothesis must still be checked. If both systems start from the same but otherwise arbitrary projective vector, the off-diagonal entries of `w` can have unbounded spread even though the initial Hilbert error is zero; differing selectors can expose those entries. A bounded-image prefix, a reset, or a declared bounded initial projective domain is therefore required for a uniform-input conclusion. This is an endpoint obligation, not a cycle obstruction.

Switching between paired faces becomes movement of the finite control state. Ties may become parallel resolved edges only after a common tie-direction cone or equivalent exact active-set state has been retained. Arbitrarily and independently choosing one maximizer in each tied output is generally not sound: the true tangent update takes a max over the complete active set, and different formal selector choices need not be simultaneously realizable by one perturbation. The exact-language hypothesis includes this compatibility requirement.

Thus the theorem genuinely covers:

- different selectors on the two channels;
- selector switching across time;
- simultaneous cancellation of the two channel increments in (6), before absolute values;
- ties, provided each resolved tied edge is jointly realizable and the global path language is exact.

For a merely sound over-approximation, absence of positive cycles is still a valid sufficient condition for bounded directed error, because every actual orbit is represented. But a positive graph cycle is then only a candidate obstruction: without path realization it can be spurious.

This separates two notions that were blurred in the earlier finite-fibre discussion:

1. **fresh arbitrary residuals**, for which the robust tropical lower bounds and suffix-set obstruction remain decisive;
2. **coherent finite switching kernels**, for which one coordinate-witness cocycle can telescope and the finite cycle criterion is exact.

The theorem therefore escapes the scalar/independently-paid-channel no-go results only through a real new hypothesis: temporal coherence of the joint selector state.

## 5. Why the naive local face graph is not exact

There is a minimal all-finite counterexample. Fix `0<delta<1/3` and use the max-plus matrix

\[
S_\delta=
\begin{pmatrix}
0&0\\
-1-\delta&-\delta
\end{pmatrix},
\qquad
(F_{S_\delta}u)_b=\max_a(u_a+(S_\delta)_{ab}).
\]

In projective coordinate `z=u_2-u_1`, its induced map is

\[
F_\delta(z)
=\max(0,z-\delta)-\max(0,z-1-\delta)
=\operatorname{clip}(z-\delta,0,1).
\tag{3}
\]

The maximal slope-one selector cell is

\[
C=(\delta,1+\delta).
\]

It has a local self-adjacency in the usual face graph: for every `z in (2delta,1+delta)`, both `z` and `F_delta(z)=z-delta` lie in `C`, so `F_delta(C) intersect C` is nonempty. Nevertheless every orbit in `C` loses `delta` at each step and leaves the cell after finitely many iterations (uniformly in the starting point of this bounded cell). There is no orbit realizing arbitrary repetitions of the apparent loop.

Therefore a graph built only from nonempty one-step face intersections can contain a false pumpable cycle. The missing datum is not another local label; it is compatibility of arbitrarily long path cylinders. Ties can create further false concatenations, not repair this one.

The exact theorem above must consequently be stated either:

- for a genuinely Markov/path-realizing partition;
- for a symbolic presentation whose path language is independently proved exact; or
- one-sidedly, with a sound graph used only to prove boundedness.

There is an equally small tie warning. The all-finite matrix

\[
S_{\rm tie}=
\begin{pmatrix}
0&1\\
-1&0
\end{pmatrix}
\]

induces the constant projective map `z->1`. At `z=1`, both inputs tie in both output coordinates. Independent symbolic choices include the identity and swap selectors, neither of which is a projective reset, even though the exact map is a complete reset. Its directional update is

\[
e\longmapsto
\big(\max(e_1,e_2),\max(e_1,e_2)\big),
\]

which is projectively zero. Thus a tie-aware exact presentation must retain the common perturbation cone or the full active-set tangent map; treating tied selectors as independent letters can manufacture false surviving directions and false cycles.

## 6. An unrestricted exact algorithm is impossible

The need for a structural subclass is not just a technical inconvenience. Krob proved that equality of rational series over the tropical semiring is undecidable (for alphabets with at least two letters): Daniel Krob, [*The equality problem for rational series with multiplicities in the tropical semiring is undecidable*](https://doi.org/10.1142/S0218196794000063), *International Journal of Algebra and Computation* 4 (1994), 405-425; an accessible record is also available at [HAL](https://hal.science/hal-00017323).

That theorem yields the following sharper boundary.

### Proposition (undecidability of bounded-versus-pumpable tropical response)

There is no algorithm for arbitrary integer-weight max-plus rational series `F,G` that, under the promise that one of the following holds, always decides which one holds:

1. `sup_w |F(w)-G(w)|<infinity`;
2. there are fixed words `u,v,w`, a constant `c>0`, and arbitrarily large `k` such that

\[
|F(uv^kw)-G(uv^kw)|\ge kc-O(1).
\]

The claim remains a directed/projective-response statement after anchoring each scalar output as the two-coordinate vector `(F(w),0)`.

### Reduction

Start with two max-plus rational series

\[
f(w)=\alpha_f A_w\beta_f,
\qquad
g(w)=\alpha_g B_w\beta_g.
\]

First totalize them without changing their equality problem. Choose `M` bounding the absolute values of all finite presentation weights in both systems and choose `C>M`. Add the same constant `C` to every transition and compatible common constants to initial and final weights. Every finite accepting value then becomes positive, while a nonaccepted word still has value `-infinity`. Take the max with the constant-zero series. The resulting series are finite everywhere, and equality of the transformed pair is equivalent to equality of the original pair because both received the same length-dependent shift before the zero fallback.

For each series add a delimiter letter `#` represented by the rank-one tropical matrix

\[
R_f=\beta_f\otimes\alpha_f,
\qquad (R_f)_{ij}=\beta_f(i)+\alpha_f(j),
\]

and similarly for `g`. For every row vector `x`,

\[
xR_f=(x\beta_f)+\alpha_f.
\]

Consequently, on a segmented word

\[
u_1\#u_2\#\cdots\#u_k,
\]

the extended series evaluate to

\[
F=\sum_{t=1}^k f(u_t),
\qquad
G=\sum_{t=1}^k g(u_t).
\]

If `f=g`, the response difference is identically zero. If `f(u)-g(u)=c` is nonzero for one word `u`, then the word

\[
(u\#)^{k-1}u
\]

has difference `k c`, hence absolute difference `k|c|`. This is of the promised fixed entry/cycle/exit form (for example, take the repeat block `u#` and absorb the final delimiter convention into the fixed suffix). Thus an algorithm for the promised bounded-versus-pumpable dichotomy would decide tropical rational-series equality, contradicting Krob's theorem.

This reduction uses sparse tropical matrices with `-infinity` entries. It does **not** establish undecidability for the all-finite, bounded-projective-image, or fixed-low-dimensional subclasses relevant to some concrete max-plus models. That distinction is important: those are precisely the subclasses in which an exact finite Markov presentation might still be forced by geometry.

## 7. A checkable sufficient realization criterion

The exact-presentation hypothesis is automatic under a strong but concrete finite lumpability condition.

### Proposition (invariant finite-cell realization)

Let `X` be an invariant projective state space partitioned into finitely many nonempty cells `(C_q)_(q in V)`. Let each continuation letter `a` act by a map `F_a:X->X`. Suppose that for every cell `q` and letter `a` there are a uniquely declared next cell `delta(q,a)`, a selector `sigma_(q,a)`, and a vector `b_(q,a)` such that

\[
F_a(C_q)\subseteq C_{\delta(q,a)}
\]

and, throughout `C_q`,

\[
F_a z=P_{\sigma_{q,a}}z+b_{q,a}.
\]

Then the deterministic cell automaton is an exact regular paired-selector presentation. In particular, its ordered-witness graph gives an exact boundedness-versus-linear-drift criterion for every finite family of coordinate-difference observations.

This remains valid when cells are lower-dimensional tie faces: different active selectors may be retained as parallel labels if their affine formulas agree on the entire face. What is forbidden is a purported transition that is feasible only on a proper subregion whose future feasibility was forgotten.

#### Proof

Choose any initial point in a permitted nonempty cell. Cell inclusion inductively forces its trajectory under every word to follow the automaton path, and the affine identities give the displayed selector recursion. Conversely, every word-defined automaton path from a permitted cell is realized by applying that word to any point of the cell; whole-cell inclusion ensures that no later transition loses feasibility. Thus soundness and path realization both hold, and the witness-cycle theorem applies. QED.

This proposition gives a precise, checkable version of the proposed “invariant normal fan” mechanism. It is stronger than necessary, but unlike local face adjacency it actually proves the concatenation property. It also exposes the approximation problem sharply: replacing whole-cell inclusion by nonempty intersection destroys path realization, as the clamp example shows.

## 8. Research judgment and next theorem

The exact finite theorem is now clear, as is its limitation:

- If an exact regular paired-selector presentation is supplied, the ordered-witness graph gives a complete and computationally small boundedness/drift criterion.
- A local selector-face graph is not automatically such a presentation.
- No procedure can construct an exact finite criterion for all tropical rational systems.

The strongest next theorem is therefore a **finite tropical lumpability theorem**:

> Give checkable hypotheses on a family of max-plus maps and its declared response queries under which a finite polyhedral quotient is invariant and path-realizing; then the ordered-witness cycle criterion decides bounded response error.

Promising concrete hypotheses, in decreasing order of conceptual value, are:

1. a common invariant finite normal fan whose image of every cone is contained in a single cone (with boundary ties resolved compatibly);
2. a finite projective orbit or finite collection of invariant metric shells;
3. an explicit synchronizing/reset word making all post-reset selector paths Markov;
4. a bounded integer-weight, fixed-dimension all-finite subclass admitting an effective projective partition.

The first is the cleanest unifying target. It says exactly when a switching-cell description is a true quotient rather than an over-approximation. Once proved, the witness-cycle theorem supplies the interacting-composition result automatically and without separately paying scalar channels.

## Status classification

- **Proved:** ordered-witness cycle theorem under a stated exact regular presentation.
- **Proved:** invariant finite-cell lumpability implies such an exact presentation.
- **Proved:** all-finite two-dimensional clamp falsifying naive local face-cycle pumping.
- **Proved from published undecidability:** no unrestricted bounded-versus-pumpable decision algorithm for arbitrary sparse tropical rational series.
- **Open:** natural intrinsic hypotheses guaranteeing an exact finite path-realizing quotient for the benchmark systems of interest.
