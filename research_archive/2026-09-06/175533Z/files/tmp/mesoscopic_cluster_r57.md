# Wave 57 memo: the mesoscopic route needs Johnson-correlated witness reuse

Status: **Verified combinatorial implication and abstract counterexample.**
The exact identities below were independently checked on the finite slice
`(n,m,ell)=(12,7,3)` by
`tmp/mesoscopic_cluster_r57_check.py`.  The counterexample is an abstract
favorable-incidence system, not a claim about an exact minimizing signing.

The conclusion is fairly sharp.  Distinct-selector production, balanced
fibres of the count-feasible size, and bounded pair-witness multiplicity do
not imply mesoscopic spectral excess.  What is missing is a correlation
statement: a high-overlap Johnson neighbour must retain the *same* witness
with a precisely quantified probability.

## 1. Incidence notation

Let

```math
\Omega=\binom{[n]}m,\qquad N=|\Omega|,
```

and let `z` range over a finite class of controlled-row centers.  Write

```math
F_z=\{S\in\Omega:(S,z)\text{ is favorable}\},\qquad r_z=|F_z|>0.
```

For `2<=ell<m`, put

```math
d_\ell=\binom m\ell\binom{n-\ell}{m-\ell},\qquad
K_\ell(S,T)=\frac{\binom{|S\cap T|}{\ell}}{d_\ell},
```

```math
h_\ell=K_\ell(S,S)=\binom{n-\ell}{m-\ell}^{-1},
\qquad
\lambda _2(\ell)=\frac{b_\ell}{d_\ell},
```

where

```math
b_\ell=\binom{m-2}{\ell-2}\binom{n-\ell-2}{m-\ell}.
```

The center retention and the project-row aggregate are

```math
p_\ell(z)=\frac1{r_z}\sum_{S,T\in F_z}K_\ell(S,T),
\qquad
P_\ell=\frac{\sum_z r_z^2p_\ell(z)}{\sum_zr_z^2}.
```

These are exactly the normalizations of (10.1266), (10.1292).

For `ell<=s<m`, let

```math
\Gamma_s(S)=\{T\in\Omega:T\ne S,\ |S\cap T|\ge s\},
\qquad D_s=|\Gamma_s(S)|
=\sum_{j=s}^{m-1}\binom mj\binom{n-m}{m-j}.
```

Define the ordered high-overlap pair count

```math
E_s(z)=|\{(S,T)\in F_z^2:S\ne T,\ |S\cap T|\ge s\}|.
```

## 2. Exact threshold-retention lemma

Choose a center with probability `r_z^2/sum_x r_x^2`, then choose `S`
uniformly from `F_z`, and finally choose `T` uniformly from
`Gamma_s(S)`.  The probability that the same center also witnesses `T` is

```math
\boxed{
\Pi_s=
\frac{\sum_z r_zE_s(z)}{D_s\sum_zr_z^2}.}
\tag{R57.1}
```

This size bias is forced by the spectral aggregate; it is not an arbitrary
choice of sampling law.

**Threshold-retention lemma.**  For every favorable incidence graph,

```math
\boxed{
P_\ell\ge h_\ell+R_{\ell,s}\Pi_s,
\qquad
R_{\ell,s}=\frac{D_s\binom s\ell}{d_\ell}.}
\tag{R57.2}
```

Consequently, for every `epsilon>0`,

```math
\boxed{
\Pi_s\ge
\frac{\lambda _2(\ell)-h_\ell+\epsilon}{R_{\ell,s}}
\quad\Longrightarrow\quad
P_\ell-\lambda _2(\ell)\ge\epsilon.}
\tag{R57.3}
```

**Proof.**  The diagonal contributes exactly `h_ell`.  Every ordered pair
counted by `E_s(z)` contributes at least `binom(s,ell)/d_ell`.  Hence

```math
\begin{aligned}
P_\ell
&=h_\ell+
\frac{\sum_zr_z
 \sum_{\substack{S,T\in F_z\\S\ne T}}
 \binom{|S\cap T|}{\ell}}
 {d_\ell\sum_zr_z^2}\\
&\ge h_\ell+
\frac{\binom s\ell\sum_zr_zE_s(z)}
 {d_\ell\sum_zr_z^2}
=h_\ell+R_{\ell,s}\Pi_s.
\end{aligned}
```

This proves both claims.  It is (10.1294) rewritten as an exact property of
the selector--witness incidence graph, but (R57.1) identifies the required
notion of witness reuse and removes all superfluous pointwise assumptions.

There is an equivalent codegree formulation.  Put

```math
c_r(S,T)=\sum_{z:S,T\in F_z}r_z,
\qquad d_r(S)=\sum_{z:S\in F_z}r_z.
```

Then

```math
\sum_Sd_r(S)=\sum_zr_z^2,
```

and exact double counting gives

```math
\boxed{
\Pi_s=
\frac{\sum_S\sum_{T\in\Gamma_s(S)}c_r(S,T)}
 {D_s\sum_Sd_r(S)}.}
\tag{R57.4}
```

The same codegrees give the complete weighted histogram without loss:

```math
\boxed{
P_\ell=h_\ell+
\frac{
 \sum_S\sum_{T\ne S}
 \binom{|S\cap T|}{\ell}c_r(S,T)}
 {d_\ell\sum_Sd_r(S)}.}
\tag{R57.4a}
```

Thus (R57.3) is the weakest *global average codegree condition* that uses
only one overlap threshold.  A convenient stronger, local hypothesis is

```math
\sum_{T\in\Gamma_s(S)}c_r(S,T)
\ge\theta D_sd_r(S)\quad\hbox{for every }S,
\tag{R57.5}
```

with

```math
\theta\ge
\frac{\lambda _2-h_\ell+\epsilon}{R_{\ell,s}}.
```

No lower bound on the number of witnesses or on the number of distinct
selectors in a fibre substitutes for this shared-witness codegree.

## 3. The exact requirement at `c=1/8`

Fix `p=m/n` in a compact subset of `(1/2,1)`, let

```math
a=\frac{1-p}{p},\qquad H=n^{5/8},\qquad
\ell=\kappa n^{13/16}(1+o(1)).
```

This is the entropy-matched scale `sqrt(nH)` and lies below `n^(5/6)`.
Let `s_*` maximize `R_{ell,s}`.  The Wave 56 local-limit calculation gives

```math
R_{\ell,s_*}
=(1+o(1))\frac1{\sqrt{2\pi}\,a\kappa}n^{-5/16},
```

whereas

```math
\lambda _2(\ell)
=(1+o(1))a^2\kappa^2n^{-3/8},
\qquad h_\ell=e^{-\Theta(n)}.
```

Therefore the zero-excess retention baseline is

```math
\boxed{
\theta_0:=\frac{\lambda _2-h_\ell}{R_{\ell,s_*}}
=(1+o(1))\sqrt{2\pi}\,a^3\kappa^3n^{-1/16}.}
\tag{R57.6}
```

To save `epsilon=exp(-C H)`, the exact reduced lemma is

```math
\boxed{
\Pi_{s_*}\ge\theta_0+
\frac{e^{-CH}}{R_{\ell,s_*}}
=\theta_0+e^{-(C+o(1))H}.}
\tag{R57.7}
```

Equivalently, under the size-biased incidence law, an average favorable
selector must have at least

```math
D_{s_*}\left\{
(1+o(1))\sqrt{2\pi}a^3\kappa^3n^{-1/16}
+e^{-(C+o(1))H}\right\}
```

distinct `s_*`-overlap partners carrying the *same* center.  A fixed
multiplicative improvement over `theta_0` is much stronger than needed: it
would give polynomial, rather than saved-exponential, spectral excess.

This quantifies the gap left by migration.  Polynomial adjacent stars do
not approach (R57.7), because `D_{s_*}` is exponentially large.  Dependent
random choice applied only to the selector--center density also does not
know that the pair found lies in `Gamma_{s_*}`.

## 4. Balanced partitions falsify degree and bounded-multiplicity packages

The following construction is exact and requires no asymptotic probability
estimates.

**Balanced-partition obstruction.**  Fix positive block sizes
`r_1,...,r_B` with `sum_i r_i=N`, and choose a uniformly random partition

```math
\Omega=F_1\mathbin{\dot\cup}\cdots\mathbin{\dot\cup}F_B,
\qquad |F_i|=r_i.
```

Regard each block as the favorable family of one center.  Set

```math
q_{\rm part}=
\frac{\sum_i r_i^2(r_i-1)}{(N-1)\sum_i r_i^2}.
```

Then, exactly,

```math
\boxed{
\mathbb E P_\ell=h_\ell+(1-h_\ell)q_{\rm part},
\qquad
\mathbb E\Pi_s=q_{\rm part}\quad\hbox{for every }s.}
\tag{R57.8}
```

In particular,

```math
q_{\rm part}\le\frac{r_{\max}-1}{N-1}.
\tag{R57.9}
```

**Proof.**  For each ordered distinct pair `(S,T)`, the probability that
both lie in block `i` is `r_i(r_i-1)/(N(N-1))`; its contribution to the
`r_i`-weighted numerator is therefore multiplied by `r_i`.  The stochastic
row identity

```math
\sum_{T\ne S}K_\ell(S,T)=1-h_\ell
```

gives the first expectation.  The overlap graph `Gamma_s` is regular of
degree `D_s`, so the same calculation with its `ND_s` ordered edges gives
the second.  Inequality (R57.9) is immediate.

Now choose

```math
B=\lceil e^H\rceil
```

and balanced sizes `r_i in {floor(N/B),ceil(N/B)}`.  Since `H=o(n)` and
`N=e^{Theta(n)}`,

```math
q_{\rm part}=e^{-(1+o(1))H},
\qquad
\mathbb E P_\ell=e^{-(1+o(1))H}\ll\lambda _2(\ell).
```

Hence at least one deterministic partition satisfies

```math
\boxed{P_\ell-\lambda _2(\ell)<-\tfrac12\lambda _2(\ell)}
\tag{R57.10}
```

for all sufficiently large `n`.  Its saved positive excess is zero.

This incidence graph has all of the following strong regularity:

- every selector has exactly one witness;
- every center has `(1+o(1))N e^{-H}` distinct selectors;
- every selector pair has codegree at most one;
- the center fibres partition the whole slice, with no repeated incidence.

Taking `d` labelled copies gives exact left degree `d`, right degrees of the
same size, pair codegree at most `d`, and the same `P_ell`.  Thus neither a
perfect cover, target-scale distinct fibres, nor bounded witness
multiplicity yields any spectral excess.  The block label is statistically
independent of Johnson geometry.  Moreover, for the deterministic partition
chosen in (R57.10), (R57.2) itself gives

```math
\Pi_{s_*}\le\frac{P_\ell-h_\ell}{R_{\ell,s_*}}
\le e^{-(1+o(1))H},
```

exponentially below the necessary `n^(-1/16)` in (R57.6).

This also explains why ordinary dependent random choice cannot close the
gap from degree data: the partition construction already saturates the
generic density-scale common-neighbour behavior while avoiding the needed
geometric correlation.  A successful signing theorem must force
`Gamma_{s_*}`-correlated reuse, not merely many incidences.

## 5. Assessment of the Johnson small-set-expansion literature lead

Khot--Minzer--Moshkovitz--Safra, *Small-Set Expansion in the Johnson Graph*,
Theory of Computing 21 (2025), Article 2,
<https://theoryofcomputing.org/articles/v021a002/>, is conceptually aligned
with (R57.4): nonexpanding Johnson families must correlate with intersections
or unions of basic coordinate sets.  Their Theorem 1.3, however, does not
currently supply a quantitative input here.

There are four precise mismatches.

1. Their theorem is stated for `J(n,k,alpha k)` with sequential quantifiers
   `k` sufficiently large and then `n>=n_0(k)`.  No uniformity is supplied
   for the fixed-density diagonal `k=pn` used here.

2. Their conclusion resolves a fixed expansion defect `eta>0`.  We need a
   vanishing retention of order `lambda_2=n^(-3/8)` and must distinguish it
   to additive accuracy `e^{-Theta(H)}`.  The paper gives no parameter
   dependence usable at that resolution.

3. Their pseudorandomness parameters are fixed after the fixed defect is
   chosen.  Our favorable densities and the natural restriction widths vary
   as `e^{-Theta(H)}` and `Theta(H)` or `ell`; the required regime is not a
   consequence of their constant-parameter statement.

4. Their graph has one exact intersection shell.  `K_ell` is a weighted
   down--up kernel (equivalently, a mixture of intersection shells).  A
   shell-by-shell import would require estimates uniform both in the
   fixed-density vertex size and across the moderate-deviation window; these
   are not stated.

Thus the paper supports the *research heuristic* that failure of expansion
should create coordinate-core structure, but invoking it here would leave
both the fixed-density extension and all decisive quantitative dependence
unproved.  The exact conditional-load identity already available in
(10.1266),

```math
\langle f,K_\ell f\rangle
=\mathbb E_R\Pr(F\mid R)^2,
```

is the applicable finite statement of that heuristic.  What is still needed
is a minimizer-specific reason that favorable incidences have the
size-biased codegree in (R57.7).

## 6. Research judgment

The one-threshold mesoscopic route is not falsified, but its missing theorem
is now exact:

> At some controlled row cap and `ell~sqrt(nH)`, prove that a
> size-biased favorable incidence retains its center under an optimal
> high-overlap Johnson move with probability at least the polynomial
> baseline (R57.6), plus the saved margin in (R57.7).

This is substantially stronger than producing `N e^{-O(H)}` selector
labels across witnesses or bounding preimage multiplicity.  Any proposed
migration-to-cluster argument should be tested directly against `Pi_s`.
If it only controls left/right degrees, ordinary pair codegrees, or local
adjacent stars, the balanced-partition obstruction applies.

For the full weighted histogram, the analogous exact target is obtained by
replacing the threshold edge indicator in (R57.4) with
`binom(|S cap T|,ell)/d_ell`; it can exploit several shells, but it still
requires Johnson-correlated same-witness reuse.  Nothing in the present
active-face law supplies that correlation.
