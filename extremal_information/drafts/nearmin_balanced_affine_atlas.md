# A balanced affine atlas inside every exact-minimizer shell

Date: 2026-08-17.

Status: proof draft for independent audit.  This note combines the
fractional-balance law FB.1--FB.2 with the multiscale affine-shell argument
MP.3.  The new point is that an arbitrary near-top atom can be thickened
*in place*: no local-ascent map, and hence no uncontrolled change of its
edge barycentre, is needed.

The output is a sublinear collection of large, one-sided affine response
charts whose uniform mixture is still edge-balanced.  It is a strict,
existential one-block response compression for a signing-dependent designed
query bank--not an efficient encoder or a common exogenous query quotient.
The last section gives the exact Walsh obstruction to promoting the atlas to
a language which mixes ports from different charts.

Throughout, `A` is a hollow symmetric signing of order `n`,

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q=Q(A)=\max_x|H_A(x)|,
\qquad E={n\choose2}.
```

For `sigma in {+-1}` and `x in {+-1}^n`, write

```math
z(sigma,x)=(sigma x_ix_j)_{i<j},
\qquad d(sigma,x)=Q-\langle a,z(sigma,x)\rangle.
```

## 1. Direct thickening of a nonstationary shell atom

Switch and orient at `(sigma,x)`:

```math
D=\sigma\,\operatorname {diag}(x)A\operatorname {diag}(x),
\qquad P=H_D(1)=Q-d,
```

and let `ell_i=sum_jd_ij` be its local fields.  They need not be
nonnegative.

### Lemma AA.1 (negative local-field mass)

Put

```math
L_-:=\sum_i(-\ell_i)_+.
```

Then

```math
\boxed{L_-\le2\sqrt{Qd}.}                                      \tag{AA.1}
```

#### Proof

Let `J={i:ell_i<0}` and include every vertex of `J` independently with
probability `p` in a random set `R`.  The exact flip identity gives

```math
\mathbb E[H_D(1^R)-P]
=2pL_-+4p^2\sum_{\{i,j\}\subseteq J}d_{ij}
\ge2pL_--4p^2Q.                                                \tag{AA.2}
```

The last inequality is the subset completion bound MP.2.  Every realization
has energy at most `Q`, so the left side is at most `d`.  Taking `p=1`
first gives `L_-<=(d+4Q)/2<=3Q`, since `d<=2Q`.  Hence
`p=L_-/(4Q)<=3/4` is admissible.  Substitution into (AA.2) gives
`d>=L_-^2/(4Q)`, which is (AA.1). `square`

### Theorem AA.2 (affine chart through every near-top atom)

For every integer `2<=q<=n`, there is a set `I subseteq[n]` of size

```math
k=\lfloor n/q\rfloor
```

such that, for every `S subseteq I`,

```math
\boxed{
 \sigma H_A(x^S)
 \ge Q-d-{8Q+4\sqrt{Qd}\over q}.}                            \tag{AA.3}
```

Thus the affine cube through the *original* atom has one common quadratic
orientation.  It is closed under every odd coordinatewise product.

#### Proof

Write `L_+=sum_i(ell_i)_+`.  Since `sum_iell_i=2P`, Lemma AA.1 gives

```math
L_+=2P+L_-
\le2Q+2\sqrt{Qd}.                                             \tag{AA.4}
```

Partition `[n]` into `q` cells whose sizes differ by at most one.  For a
cell `J_b`, put

```math
C_b=2\sum_{i\in J_b}(\ell_i)_++4Q_-(D[J_b]).
```

The partition cap budget MP.1 and (AA.4) imply

```math
\sum_bC_b\le2L_++4Q\le8Q+4\sqrt{Qd}.                         \tag{AA.5}
```

Choose a cell of cost at most the average and, if necessary, choose a
`k`-element subset `I` of it.  Principal restriction cannot increase the
one-sided negative cap: a spin on `I` can be completed randomly on the
rest of the cell, preserving its internal energy in expectation.  Hence
the cost does not increase on passing to `I`.

For `S subseteq I`, the exact flip identity and MP.2 give

```math
\begin{split}
P-H_D(1^S)
 &=2\sum_{i\in S}\ell_i
     -4\sum_{\{i,j\}\subseteq S}d_{ij}\\
 &\le2\sum_{i\in I}(\ell_i)_++4Q_-(D[I])\\
 &\le {8Q+4\sqrt{Qd}\over q}.
\end{split}                                                   \tag{AA.6}
```

Switch back to obtain (AA.3).  Odd products XOR the masks. `square`

The `sqrt(Qd)` term is exactly what avoids the archived local-ascent
obstruction: a near-top atom may have many negative local fields, but their
*total mass* is controlled, so the chart need not move its centre.

## 2. Thickening a balanced shell law

For the chart in AA.2, let `nu_(sigma,x,I)` be the uniform law on its
`2^k` signed cuts `z(sigma,x^S)`.  Uniform mask averaging erases precisely
the edges incident with `I`.  Therefore

```math
{1\over E}\sum_e
 \left|\mathbb E_\nu[a_ez_e]-a_ez(sigma,x)_e\right|
=\theta_(n,k):={kn-k(k+1)/2\over E}
\le {2k\over n-1}.                                          \tag{AA.7}
```

### Theorem AA.3 (balanced affine-atlas implication)

Suppose `mu` is a probability law on signed cuts satisfying

```math
d(z)\le d_0\quad(mu\hbox{-a.s.}),
\qquad
{1\over E}\sum_e|\mathbb E_\mu[a_ez_e]|\le\delta.          \tag{AA.8}
```

For every integers `K>=1` and `2<=q<=n`, there are at most `K` charts

```math
\mathcal C_r=\{z(\sigma_r,x_r^S):S\subseteq I_r\},
\qquad |I_r|=\lfloor n/q\rfloor,                            \tag{AA.9}
```

and a law `bar nu` which chooses a chart uniformly and then chooses a
uniform atom in that chart, such that:

1. every atlas atom has the chart's common one-sided orientation and

```math
d(z)\le
D_(q,d_0):=d_0+{8Q+4\sqrt{Qd_0}\over q};                    \tag{AA.10}
```

2. the atlas law retains fractional edge balance,

```math
\boxed{
 {1\over E}\sum_e|\mathbb E_(\bar\nu)[a_ez_e]|
 \le\delta+K^{-1/2}+{2\lfloor n/q\rfloor\over n-1}.}       \tag{AA.11}
```

#### Proof

Sample `K` independent atoms from `mu`.  Coordinatewise Jensen and the
variance bound show that the expected normalized `l_1` distance between
their empirical signed-edge mean and the mean of `mu` is at most
`K^{-1/2}`.  Fix a sample attaining that bound.  Apply AA.2, with the same
`q`, to each sampled atom and average the resulting uniform chart laws.
Equation (AA.10) follows from AA.3.  The triangle inequality and (AA.7)
give (AA.11). `square`

This is a genuine compressed representation of the two conclusions in
(AA.8): it stores `K` triples `(sigma_r,x_r,I_r)`, rather than listing any
of the `K2^k` generated shell atoms.  Its existence proof may inspect the
full minimax shell law and the relevant principal caps; it supplies no
efficient procedure for finding the atlas.  The information retained by the
certificate, rather than its computational discovery cost, is what is being
bounded.

## 3. Declared one-block response composition

After dropping at most one coordinate, take `I_r' subseteq I_r` even and
form the star frame

```math
W_r=(x_r,(x_r^{\{i\}})_{i\in I_r'}).
```

For the oriented one-sided trust response

```math
\mathcal R_(A,\sigma)(g)
=\max_y\{\sigma H_A(y)+g\mathbin\cdot y\},                  \tag{AA.12}
```

and the common absolute trust response

```math
\mathcal B_A(g)=\max_{y,\tau\in\{+-1\}}
 \{\tau H_A(y)+g\mathbin\cdot y\},
```

the majority-selector calculation of MP.3 now gives, simultaneously for
every chart `r`, every `epsilon in {+-1}^{|I_r'|+1}`, and every `m>=0`,

```math
\boxed{
0\le Q+m\|W_r\epsilon\|_1
 -\mathcal R_(A,\sigma_r)(mW_r\epsilon)
\le D_(q,d_0).}                                             \tag{AA.13}
```

The identical bound holds with `mathcal R_(A,sigma_r)` replaced by the
single response `mathcal B_A`, for all charts at once.

Indeed the sign selector of `W_repsilon` belongs projectively to
`{x_r^S:S subseteq I_r'}`, pays the full field norm, and has the one-sided
energy from AA.10; cap plus Holder is the upper bound.

Thus the atlas has an honest, but signing-dependent, one-block composition
payoff: it gives an `O(D_(q,d_0))` response roof for the declared union of
`K` growing star-interface languages designed from `A`.  It is not merely a
shell cover.  It is also not an efficient encoder or a quotient for one
common exogenous query class shared by different signings.  The state has
an `O(Kn)`-bit direct presentation, whereas the generated response bank has
`K2^{Theta(n/q)}` endpoints.  For `K=o(n)` this is subquadratic and hence
strictly smaller than an arbitrary edge landscape or a full Boolean
response table.

This existential response certificate includes the scalar baseline `Q`; it
does not provide an algorithm for computing `Q`, the centres, or the
supports.  Nor is `mathcal B_A` by itself an ordinary
fixed-orientation signing parent: the auxiliary sign `tau` is part of the
declared query.  The one-sided statement in (AA.13) is physical after the
chart orientation has been fixed, while the common absolute statement is
the chart-union response object used by the repository's roof formalism.

## 4. Exact-minimizer consequence

More generally, suppose `A_n in N_n(epsilon_n)`, `epsilon_n ->0`, and choose
`kappa_n ->0` so that

```math
{\epsilon_n+\eta_n(\kappa_n)\over\kappa_n}\longrightarrow0.
                                                                    \tag{AA.14a}
```

(For example one may take a slowly enlarged maximum of
`sqrt(epsilon_n)` and `n^(-1/3)`.)  If `K_n,q_n -> infinity` while
`K_n,q_n=o(n)`, then FB.1--FB.2 and AA.3 give a `K_n`-chart atlas with

```math
{D_(q_n,2\kappa_nn^{3/2})\over n^{3/2}}
=O(\kappa_n+q_n^{-1})=o(1),
```

normalized edge bias

```math
O\left(
 {\epsilon_n+\eta_n(\kappa_n)\over\kappa_n}
 +K_n^{-1/2}+q_n^{-1}
\right)=o(1),                                                \tag{AA.14b}
```

and chart dimension `floor(n/q_n)`.  Thus the atlas conclusion holds
uniformly throughout every genuinely vanishing near-minimizer class, not
only at exact minimizers.

Use FB.1 for an exact minimizer with `kappa=n^(-1/6)`, and FB.2 with

```math
K=\lceil n^{1/3}\rceil,
\qquad q=\lceil n^{1/6}\rceil.                              \tag{AA.14}
```

The standard random-sign union bound gives `M_n=O(n^(3/2))`, while the
known rigorous lower bound gives `M_n>=c n^(3/2)`.  The latter makes the
absolute-bias clause FB.6 applicable for all sufficiently large `n`.
Thus FB.1 gives `d_0=2n^(4/3)` and `delta=O(n^(-1/6))`.
Consequently AA.3 gives:

### Corollary AA.4 (sublinear balanced atlas for every exact minimizer)

Every sufficiently large exact minimizer has an atlas of at most
`O(n^(1/3))` charts, each of affine dimension `Theta(n^(5/6))`, such that

```math
\begin{array}{ll}
\text{common-oriented shell deficit in every chart}
   &=O(n^{4/3})=o(n^{3/2}),\\[2mm]
\text{normalized signed-edge bias of the atlas law}
   &=O(n^{-1/6}),\\[2mm]
\text{one-sided response-roof error on every declared chart endpoint}
   &=O(n^{4/3}).
\end{array}                                                  \tag{AA.15}
```

The existential, signing-dependent atlas has an `O(n^(4/3))`-bit direct
presentation and generates

```math
\exp(\Theta(n^{5/6}))
```

one-sided selector endpoints per chart.  This is a Level-5 theorem about
actual exact minimizers at arbitrary orders.  As a labelled generative
certificate it uses substantially fewer bits than an arbitrary edge
signing or full Boolean response table; it is not a common contextual
quotient and does not assert an efficient construction.

The exponents in (AA.14) merely balance the three errors.  More generally,
any `K,q -> infinity` with `K,q=o(n)` gives a sublinear number of charts,
vanishing empirical/erasure bias, and vanishing normalized shell/response
error; the affine dimension is `floor(n/q)`.

## 5. Exact ceilings and obstruction tests

The result survives the four requested stress tests, but they sharply limit
its interpretation.

1. **Local ascent.**  No ascent is used.  AA.1 is the quantitative repair:
   a centre of deficit `d` has negative local-field mass at most
   `2sqrt(Qd)`, and its signed-edge word is left unchanged before
   thickening.

2. **Singleton exact-active shells.**  The charts start in the thick FB.1
   shell, not in `S_0`.  Thus a matrix with one projective exact maximizer is
   compatible with AA.4; the theorem does not claim exact-active
   multiplicity.

3. **Sparse-flip/geodesic planted faces.**  A thinner shell can freeze a
   planted edge set.  AA.11 pays exactly the fraction of edges erased by its
   coordinate supports, while FB.1 supplies balance only at its prescribed
   thicker scale.  Hence no frozen-edge conclusion is silently assumed
   away.

4. **Walsh mixed orientations.**  The SC.1 Walsh signing has an optimally
   (and asymptotically) balanced full exact-shell law containing three active
   centres whose cross-chart majority has zero quadratic energy.  The three
   centres alone are not claimed to form the balanced law.  Each AA chart
   remains one-sided, but no
   conclusion is asserted for selectors mixing ports from different
   charts.  Indeed SC.1 proves that the implication

```text
balanced affine atlas + one-sided response inside each chart
    ==> one-sided response for arbitrary cross-chart selectors
```

   is false, already for three singleton subcharts.  The orientation word,
   or an equivalent cross-chart Fourier-coherence statistic, is missing.

There is a second physical ceiling.  The raw exact-sign port block of one
chart has an aligned endpoint of size `Theta(n^2/q)`, and independently
paying the omitted channel in a fixed-ratio compiler incurs the MC.4
`Theta(n^(3/2))` residual unless linear rank is retained.  Therefore AA.4
does **not** give a fixed-ratio recurrence, an all-context state, or a
cross-order congruence.  Its exact implication graph is

```text
FB.1--FB.2
   + direct nonstationary thickening AA.1--AA.2
       ==> sublinear balanced affine atlas
       ==> compressed one-block response on a union of chart languages
       X=> mixed-chart selector response
       X=> fixed-ratio physical composition.
```

The positive theorem is consequently strict but scoped: near-minimality
does force a large, low-information *atlas* of coherent local response
languages.  The next missing datum is not more within-chart shell entropy;
it is a reusable cross-chart transition/coherence law.  Walsh SC.1 proves
that such a law cannot be recovered from edge balance and chart-wise
one-sidedness alone.
