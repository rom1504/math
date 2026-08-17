# Generic low-dimensional deep holes: three rigorous localization barriers

**Status:** theorem draft; no canonical edits.  The results below rule out
ordinary direct-sum/hierarchical countermodels and near-full-dimensional
linear cap carriers.  They also prove a local cap-cover entropy inequality:
low dimension forces projective spread of order the covering-radius deficit
inside an `o(t)` shell.  They do **not** decide whether an arbitrary
`[N,Theta(sqrt N)]` antipodal code can have every `o(N^(3/4))` deep-hole
shell projectively collapsed.

## 1. Setup

Let `C<=F_2^N` be a binary linear code containing the all-one word.  Put

```math
\Delta(C)=\min_{a\in\{\mathord\pm1\}^N}
 \max_{c\in C}\sum_{e=1}^N a_e(-1)^{c_e}.
```

Antipodality makes the inner maximum an absolute maximum.  The exact
covering-radius dictionary is

```math
\rho(C)={N-\Delta(C)\over2},
\qquad
t(C):={N\over2}-\rho(C)={\Delta(C)\over2}.          \tag{GL.1}
```

For a deep hole `y`, write

```math
\mathcal L_s(y)=\{c\in C:d(y,c)\le\rho(C)+s\},
\qquad
d_{\rm pr}(c,c')=\min\{d(c,c'),N-d(c,c')\}.
```

The target scaling is

```math
k=\dim C=\Theta(\sqrt N),
\qquad t(C)=\Theta(N^{3/4}).                         \tag{GL.2}
```

## 2. Direct products activate before the target shell scale

### Theorem GL.1 (radius-one block activation)

For `1<=j<=b`, let `C_j<=F_2^{N_j}` contain its local all-one word, let
`y_j` be a deep hole, and put `rho_j=rho(C_j)`.  Form

```math
C=\bigoplus_{j=1}^b C_j,
\qquad y=(y_1,\ldots,y_b),
\qquad N=\sum_jN_j.
```

Then `rho(C)=sum_j rho_j`, `y` is a deep hole, and

```math
\prod_{j=1}^b\mathcal L_{s_j}(y_j)
 \subseteq\mathcal L_{\sum_js_j}(y).                \tag{GL.3}
```

For each `j`, choose `c_j^0,c_j^1 in mathcal L_1(y_j)` and put

```math
d_j=d(c_j^0,c_j^1).
```

If

```math
\max_jN_j=o(N),
\qquad
\sum_jd_j\ge N/3,                                   \tag{GL.4}
```

then

```math
\operatorname{diam}_{\rm pr}\mathcal L_b(y)
 \ge N/3-o(N).                                      \tag{GL.5}
```

More generally, (GL.5) holds with `b` replaced by `sum_j s_j` whenever
the chosen pairs lie in the corresponding `s_j` shells.

#### Proof

Distances and covering radii add under direct sums, proving (GL.3).  Hold
the zeroth word fixed in every block.  For a subset `J` of blocks, use the
first word precisely in `J`; its distance from the all-zeroth product is

```math
w_J=\sum_{j\in J}d_j.                                \tag{GL.6}
```

Greedily add blocks until this sum first reaches `N/3`.  Assumption (GL.4)
ensures that this happens, and the overshoot is at most
`max_j d_j<=max_jN_j=o(N)`.  Hence

```math
N/3\le w_J\le N/3+o(N)<2N/3,
```

so `min(w_J,N-w_J)>=N/3-o(N)`.  Both product words belong to the shell in
(GL.3), proving (GL.5). `square`

### Consequence for the proposed scaling

Every nonzero factor uses at least one code dimension, so `b<=k`.  Under
(GL.2),

```math
b\le k=Theta(\sqrt N)=o(t(C)).                       \tag{GL.7}
```

Thus any fine direct-sum construction satisfying (GL.4) already has a
fixed-scale projective pair in an `o(t)` shell.  In particular, any family
whose exhibited radius-one alternatives have aggregate activation mass at
least `N/3` cannot produce the desired countermodel.  A nested construction
does not evade the statement at a node which is itself a genuine direct sum
and satisfies (GL.4): GL.1 applies there verbatim.  Hence a direct-sum
hierarchy avoiding the theorem must arrange that every balanced direct-sum
node has total radius-one activation mass less than one third of that node's
length.  This is a precise recursive localization obligation, not a proof
that all hierarchies fail.

This is stronger than applying sphere covering separately to the factors.
It uses exact deep-hole stability at radius one and shows why adding many
small probabilistic factors is actively harmful: all factors switch before
the `N^(3/4)` scale.

## 3. A split discrepancy bound for localized subcodes

We next quantify a different common proposal: place most code dimensions in
a projectively small carrier and leave only a few global quotient bits.

### Lemma GL.2 (split discrepancy)

Let `C<=F_2^N` have dimension `k`.  Suppose it has an `ell`-dimensional
subcode `U` supported on a coordinate set `S` of size `D`.  Then

```math
\boxed{
\Delta(C)\le
 \sqrt{2D(k+2)\log2}
 +\sqrt{2(N-D)(k-\ell+2)\log2}.}                    \tag{GL.8}
```

Consequently

```math
t(C)\le {1\over2}\left[
 \sqrt{2D(k+2)\log2}
 +\sqrt{2(N-D)(k-\ell+2)\log2}\right].              \tag{GL.9}
```

#### Proof

For any binary code `B` of length `m` and dimension `r`, independent
Rademacher signs `epsilon_i` satisfy, for a fixed sign word `z`,

```math
\Pr\{|\epsilon\mathbin\cdot z|\ge T\}
 \le2\exp(-T^2/(2m)).                                \tag{GL.10}
```

A union bound over at most `2^r` restricted sign words shows that one
choice of `epsilon` has every absolute response below

```math
\sqrt{2m(r+2)\log2}.                                 \tag{GL.11}
```

Apply this independently on `S` and `S^c`.  The restriction of `C` to `S`
has dimension at most `k`.  Because `U` lies in the kernel of puncturing to
`S^c`, the restriction there has dimension at most `k-ell`.  For every
codeword, the absolute full response is at most the sum of the two restricted
absolute responses.  Equations (GL.8)--(GL.9) follow from (GL.1). `square`

The harmless `+2` makes the union-bound probability strictly below one and
also covers zero-length endpoint cases.  No independence between the two
restricted code images is asserted or needed.

### Lemma GL.3 (a linear projective cap is actually supported locally)

Let `W<=F_2^N` contain `mathbf1`, and suppose

```math
d_{\rm pr}(w,0)\le D<N/3\qquad(w\in W).              \tag{GL.12}
```

Then

```math
W=U\oplus\langle\mathbf1\rangle,                    \tag{GL.13}
```

where every word of `U` has weight at most `D` and

```math
|\operatorname{supp}U|\le2D.                        \tag{GL.14}
```

#### Proof

Let `U={w in W:wt(w)<=D}`.  Every projective class has exactly one member in
`U`.  If `u,v in U`, then `wt(u+v)<=2D<N-D`.  By (GL.12), `u+v` has weight
at most `D` or at least `N-D`; the latter is impossible.  Thus `U` is a
subspace and (GL.13) follows.  On every coordinate used by `U`, exactly half
of its words are one.  Averaging their weights gives

```math
|\operatorname{supp}U|/2
=\mathbb E_{u\in U}\operatorname{wt}(u)\le D,
```

which proves (GL.14). `square`

### Corollary GL.4 (a collapsed near-full-dimensional projective affine carrier is impossible)

Assume (GL.2), and let

```math
\pi:C\longrightarrow C/\langle\mathbf1\rangle
```

be the projective quotient.  Suppose the projective image of some thin shell
contains an affine subspace

```math
\bar c+\bar W\subseteq C/\langle\mathbf1\rangle
```

of dimension `r`, and suppose the chosen shell representatives of these
projective classes have pairwise projective distance at most `D=o(N)`.  Then

```math
k-r=\Omega(k).                                       \tag{GL.15}
```

In particular, no projectively collapsed shell can contain a projective
affine subspace of dimension `k-o(k)`.

#### Proof

Let `W<=C` be the inverse image of the direction space `bar W`.  It contains
`mathbf1` and has dimension `r+1`.  Differences of shell representatives
show that every projective class in `W` has weight at most `D`.  For all
large `N`, Lemma GL.3 therefore gives a subcode `U<=W<=C` of dimension
`ell=r` supported on at most `2D=o(N)` coordinates.  Apply GL.2 with this
support.  The first square root in (GL.8) is

```math
o(\sqrt{Nk})=o(N^{3/4}).                             \tag{GL.16}
```

Since `Delta(C)=2t(C)=Theta(N^(3/4))`, the second square root must have that
order.  Therefore

```math
k-\ell=\Omega(\sqrt N)=\Omega(k),
```

which is (GL.15). `square`

A quantitative version is immediate.  If `t(C)>=cN^(3/4)` and
`k<=Ksqrt N`, then, for `D=o(N)`,

```math
k-r+2\ge\left({2c^2\over\log2}-o(1)\right)\sqrt N.
                                                               \tag{GL.17}
```

## 4. Local cap-cover entropy and the exact two-cap ceiling

The preceding affine-carrier argument does not control a nonlinear shell.
There is nevertheless an exact counting inequality.  It proves diffusion at
the mesoscopic `t(C)` scale and shows why the proof stops exactly when an
oppositely oriented projective cap can enter.

### Lemma GL.5 (opposite-lift entry gap)

Let `C` be antipodal, let `y` be a deep hole, and abbreviate
`rho=N/2-t`.  If `c,c' in mathcal L_s(y)` and their actual Hamming distance
uses the large projective lift,

```math
d(c,c')=N-d_{\rm pr}(c,c'),
```

then

```math
d_{\rm pr}(c,c')\ge2(t-s).                         \tag{GL.18}
```

Consequently, if
`diam_pr mathcal L_s(y)<2(t-s)`, the whole shell lies in one ordinary
Hamming cap: for any fixed `c_0 in mathcal L_s(y)`,

```math
d(c,c_0)=d_{\rm pr}(c,c_0)
<2(t-s)\qquad(c\in\mathcal L_s(y)).                \tag{GL.19}
```

#### Proof

Put `E=supp(y+c)` and `E'=supp(y+c')`.  Their sizes are at most
`rho+s`, while

```math
N-d_{\rm pr}(c,c')=d(c,c')=|E\mathbin\triangle E'|
\le |E|+|E'|\le2(rho+s)=N-2(t-s).
```

This proves (GL.18), and the consequence follows by applying it to every
pair involving `c_0`. `square`

### Theorem GL.6 (local cap-cover entropy)

Let `C<=F_2^N` be antipodal of dimension `k`, let `y` be a deep hole, and
write `rho=N/2-t`.  Fix an integer `1<=r<t`, put

```math
m=N-rho=N/2+t,
\qquad h=\lceil r/2\rceil,
```

and suppose

```math
D:=\operatorname{diam}_{\rm pr}\mathcal L_r(y)<2(t-r).          \tag{GL.20}
```

Then the shell obeys the exact information lower bound

```math
\boxed{
|\mathcal L_r(y)|
\binom rh\left({D\over m}\right)^h\ge1.}           \tag{GL.21}
```

In particular,

```math
k\ge h\log_2(m/D)-\log_2\binom rh.                 \tag{GL.22}
```

(The hypotheses force `D>0`, so the logarithm is well defined.)

#### Proof

Choose a nearest word `c_0`, and set `E_0=supp(y+c_0)`, so
`|E_0|=rho` and `|E_0^c|=m`.  For every `r`-set `F subseteq E_0^c`, the
covering-radius property supplies `c_F in C` with

```math
d(y+F,c_F)\le rho.
```

Such a word cannot equal `c_0`, because
`d(y+F,c_0)=rho+r`.
Nor can it equal `c_0+mathbf1`, because

```math
d(y+F,c_0+\mathbf1)=N-(rho+r)=rho+2t-r>rho
```

when `r<t`.  Hence `c_F` represents a genuinely different projective class
from `c_0`, and therefore `D>0`.
The triangle inequality gives `c_F in mathcal L_r(y)`.  If
`E_F=supp(y+c_F)`, deepestness of `y` gives `|E_F|>=rho`, and hence

```math
rho\ge |E_F\mathbin\triangle F|
=|E_F|+r-2|E_F\cap F|
\quad\Longrightarrow\quad |E_F\cap F|\ge h.       \tag{GL.23}
```

By GL.5 and (GL.20), every shell word uses the small lift from `c_0`.
Thus the sets

```math
G_c=E_0^c\cap\operatorname{supp}(c+c_0),
\qquad c\in\mathcal L_r(y),
```

all have size at most `D`.  Equation (GL.23) says that every `r`-subset
of the `m`-point set `E_0^c` meets at least one `G_c` in at least `h`
points.

Choose that `r`-subset uniformly and then order it uniformly.  For fixed
`G_c`, a union bound over the `binom(r,h)` choices of `h` sample positions
gives

```math
\Pr\{|F\cap G_c|\ge h\}
\le\binom rh{(|G_c|)_h\over(m)_h}
\le\binom rh(D/m)^h.                               \tag{GL.24}
```

The events over all shell words cover every `F`, so another union bound
proves (GL.21).  Finally `|mathcal L_r(y)|<=|C|=2^k` gives (GL.22).
`square`

### Corollary GL.7 (generic mesoscopic diffusion)

For any sequence of antipodal `[N,k]` codes and deep holes satisfying

```math
k=o(t),\qquad t=o(N),\qquad t\longrightarrow\infty,              \tag{GL.25}
```

there is a shell width `r=o(t)` such that

```math
\operatorname{diam}_{\rm pr}\mathcal L_r(y)
\ge(2-o(1))t.                                      \tag{GL.26}
```

Indeed, (GL.26) holds for every integer sequence `r` with `k=o(r)` and
`r=o(t)`.  In the target scaling, for example,
`r=ceil(sqrt(kt))=Theta(N^(5/8))` works.

#### Proof

Suppose instead that the diameter `D` is below `2(t-r)`.  From (GL.21),
`|mathcal L_r(y)|<=2^k`, and `binom(r,h)<=2^r`,

```math
{D\over m}\ge2^{-(k+r)/h}=1/4-o(1),                \tag{GL.27}
```

because `h=(1/2+o(1))r` and `k=o(r)`.  Thus
`D>=(1/4-o(1))m=Theta(N)`, contradicting
`D<2(t-r)=o(N)`.  Therefore `D>=2(t-r)=(2-o(1))t`.
`square`

The theorem is a genuine low-dimensional consequence, but its scale is
sharp for this proof architecture.  At projective distance `2(t-r)`, an
opposite Hamming lift can enter and automatically answer many local flip
queries.  The PP.2 two-cap construction exhibits exactly this mechanism
(at higher dimension).  Therefore cap-cover entropy alone cannot upgrade
(GL.26) from `Theta(t)` to `Theta(N)` without new control of the opposite
cap.

## 5. What these theorems decide

Three natural generic constructions are now rigorously constrained.

1. **Many independent or hierarchical blocks.**  Radius-one alternatives in
   all blocks activate at total excess at most the number of blocks.  Since
   that number is at most `k=o(t)`, any macroscopic activation mass yields
   diffusion before the target scale (GL.1).
2. **One large linear cap carrier plus a small global quotient.**  A linear
   projective cap is supported on only `o(N)` coordinates (GL.3), and the
   split discrepancy inequality forces a linear fraction of all `k`
   dimensions to remain outside it (GL.4).
3. **A single nonlinear Hamming cap.**  Exact `r`-flip stability forces the
   cap family to cover all `r`-sets at half incidence.  Its entropy exceeds
   the available `k` bits once `k=o(r)=o(t)`, forcing spread to at least
   `(2-o(1))t`.  This is exactly the threshold beyond which GL.5 no longer
   certifies a common Hamming lift; it does not prove that the witness pair
   actually uses opposite lifts (GL.5--GL.7).

These conclusions are not consequences of PP.3 alone.  PP.3 treats a fully
separable repetition-plus-code block.  GL.1 handles arbitrary direct-product
factors through their actual radius-one shell geometry, while GL.4 handles an
overlapping ambient code whenever the proposed shell carrier contains a
large affine subcode.  GL.6 is the new low-dimension input absent from PP.1's
coordinate cover: it turns all simultaneous `r`-flip obligations into a
quantitative covering-design entropy bound.

## 6. Exact unresolved loophole

An unrestricted countermodel would now have to satisfy all of the following.

* It is genuinely overlapping: no balanced direct-sum level has macroscopic
  radius-one activation mass.
* Its `o(t)` shell has projective diameter between `Theta(t)` and `o(N)`;
  GL.7 excludes a smaller scale.
* It lies beyond the common-lift regime controlled by GL.5.  It may use two
  opposite Hamming caps, or one mesoscopic cap of diameter at least
  `(2-o(1))t`, and it is **affine-subspace-evasive**: it contains no affine
  subspace of dimension `k-o(k)` in its projective image.
* It contains no affine family capturing more than `(1-Omega(1))k`
  projective directions, even though the nonlinear shell may still span the
  full quotient.
* Although exact deep-hole stability supplies local leaders covering every
  `r=o(t)` flip query, those leaders never amplify their forced
  `Theta(t)` projective spread to a fixed fraction of `N`.

The theorems above do not rule out such a nonlinear clustered shell.
Conversely, no
probabilistic, direct-sum, Plotkin-style, or low-support-subcode construction
found in this audit realizes these five requirements.  The exact remaining
generic lemma is therefore narrower than the original question:

> **Nonlinear-cluster entry.**  In an antipodal `[N,k]` code with
> `t(C)>=c sqrt(Nk)`, does every deep hole have some `s=o(t(C))` for which
> `mathcal L_s(y)` either has fixed-scale projective diameter or contains an
> affine subspace of dimension `k-o(k)` in its projective image?

An affirmative answer, together with GL.4, forces diffusion.  A negative
answer must exhibit the affine-subspace-evasive clustered shell described
above.  This is a genuinely additive-combinatorial obligation; sphere
covering and direct-product constructions alone do not decide it.

## 7. Director classification

* **Full generic diffusion theorem:** open.
* **Low-dimensional collapsed countermodel:** not found.
* **Direct-sum/hierarchical route:** rigorously obstructed by GL.1 unless it
  recursively assumes localized factors.
* **Large linear carrier route:** rigorously obstructed by GL.2--GL.4.
* **Sub-threshold nonlinear route:** quantitatively obstructed by
  GL.5--GL.7; thin shells must already spread on the `Theta(t)` scale.
* **New precise bottleneck:** a nonlinear, affine-subspace-evasive shell
  beyond the common-lift threshold, and amplification from `Theta(t)` to
  `Theta(N)` projective spread.

This is a decisive construction barrier, not evidence that the augmented cut
code satisfies Cut-DH(3).  The remaining loophole is broad enough that no
cut-specific theorem should be inferred from these generic results.
