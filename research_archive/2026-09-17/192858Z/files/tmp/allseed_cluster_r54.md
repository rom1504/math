# Wave 54: all-seed exchange and the one-threshold feasibility wall

## Status and scope

This attack does **not** prove large-core excess, the restriction estimate, or
convergence.  It has two rigorous outputs.

1. The one-threshold sufficient condition (10.1294) has a sharp entropy
   feasibility wall.  After the diagonal self-loop has already been removed,
   its required number of partners is never exponentially smaller than the
   entire relevant Johnson ball.  Away from one tuned equality surface, it is
   exponentially larger and therefore impossible.  In particular the
   numerical linear-scale example following (10.1295) has an asymptotically
   impossible hypothesis.
2. The positive-threshold port star extends from globally maximizing selectors
   to **every** child-ground selector, at threshold `8(m-1)`.  Aggregating these
   stars gives an exact all-seed ordered-pair inequality, but only polynomially
   many partners per base.  Pure double counting does not upgrade it to the
   exponential multiplicity demanded by a feasible version of (10.1294).

The feasibility wall falsifies only the one-threshold implementation
(10.1294) away from its equality surface.  It does **not** falsify the exact
weighted histogram identity (10.1292), saved spectral excess, or a signing
theorem controlling several intersection levels together.  The code model in
Section 4 is explicitly an abstract incidence obstruction, not a family known
to arise from an exact minimizing signing.

The computations are reproduced by `tmp/allseed_cluster_r54_check.py`; saved
output is in `tmp/allseed_cluster_r54_check.out`.

## 1. Entropy feasibility of the one-threshold cluster lemma

Assume

```math
\frac mn\to p\in(1/2,1),\qquad
\frac\ell n\to\alpha,\qquad
\frac sn\to\beta,\qquad 0<\alpha<\beta<p.
```

For a fixed `m`-selector `S`, the total number of **distinct** `m`-selectors
having intersection at least `s` with `S` is

```math
\mathcal B_{n,m}(s)
=\sum_{j=s}^{m-1}\binom mj\binom{n-m}{m-j}.
\tag{R54C.1}
```

The upper limit `m-1` is important: the omitted `j=m` term is exactly the
self-loop `T=S`.  Standard largest-term estimates give

```math
\frac1n\log\mathcal B_{n,m}(s)\longrightarrow V(p,\beta),
\tag{R54C.2}
```

where, with `H(x)=-x log x-(1-x)log(1-x)`, 

```math
V(p,\beta)=
\begin{cases}
H(p),&\beta\le p^2,\\[1mm]
pH\!\left(\dfrac{p-\beta}{p}\right)
+(1-p)H\!\left(\dfrac{p-\beta}{1-p}\right),&\beta>p^2.
\end{cases}
\tag{R54C.3}
```

Indeed the exchange-distance summand with `d=m-j` has exponent

```math
pH(d/(pn))+(1-p)H(d/((1-p)n)),
```

whose maximum occurs at the typical distance `d=p(1-p)n`; if the ball stops
before that distance, its boundary term dominates.

Now let `epsilon_n=exp{-o(n)}`; this includes the saved target
`exp{-O(n^(3/4-c))}`.  Because `lambda_2(ell)` tends to a positive constant,
the right side of (10.1294), including its exact diagonal subtraction, has
rate

```math
\frac1n\log
\frac{b_\ell-\binom m\ell+d_\ell\epsilon_n}{\binom s\ell}
\longrightarrow E(p,\alpha,\beta),
\tag{R54C.4}
```

with

```math
E(p,\alpha,\beta)
=pH(\alpha/p)
+(1-\alpha)H\!\left(\frac{p-\alpha}{1-\alpha}\right)
-\beta H(\alpha/\beta).
\tag{R54C.5}
```

Here `b_ell/d_ell=lambda_2(ell)`, while
`binom(m,ell)/d_ell=h_ell` is exponential; hence the subtraction in
(R54C.4) is precisely the removed self-loop and does not alter the rate.

There is a sharp **Verified entropy comparison**:

```math
\boxed{E(p,\alpha,\beta)\ge V(p,\beta).}
\tag{R54C.6}
```

The proof is elementary and also gives all equality cases.  For fixed
`p,beta`, differentiation in `alpha` gives

```math
\frac{\partial E}{\partial\alpha}
=\log\frac{(p-\alpha)^2}{(1-\alpha)(\beta-\alpha)}.
\tag{R54C.7}
```

The sign of (R54C.7) is the sign of

```math
\Phi(\alpha)
=(p-\alpha)^2-(1-\alpha)(\beta-\alpha)
=p^2-\beta+(1+\beta-2p)\alpha.
\tag{R54C.8}
```

If `beta<=p^2`, then `Phi(alpha)>0` for `0<alpha<beta`: when its slope is
nonnegative this follows at `alpha=0`, and when its slope is negative it
follows at `alpha=beta` from `Phi(beta)=(p-beta)^2`.  Therefore
`E(alpha)>E(0)=H(p)=V` for every linear `alpha>0`.

If `beta>p^2`, then also `beta>2p-1`, so `1+beta-2p>0`; (R54C.8) has the
unique root

```math
\alpha_* =\frac{\beta-p^2}{1+\beta-2p}\in(0,\beta).
\tag{R54C.9}
```

Thus `alpha_*` is the unique global minimum.  Direct substitution, using

```math
1-\alpha_* =\frac{(1-p)^2}{1+\beta-2p},\quad
p-\alpha_* =\frac{(1-p)(p-\beta)}{1+\beta-2p},\quad
\beta-\alpha_* =\frac{(p-\beta)^2}{1+\beta-2p},
```

gives `E(p,alpha_*,beta)=V(p,beta)`.  This proves (R54C.6), with equality
for fixed positive `alpha` **only** on (R54C.9).

Since every `N_z(s)` in (10.1294) is bounded by (R54C.1), (R54C.6) has the
following exact asymptotic consequence.

- If `E>V`, the hypothesis of (10.1294) is impossible for all sufficiently
  large `n`, even for the entire slice as a favorable family.
- On the tuned surface `alpha=alpha_*`, the hypothesis is not ruled out by
  exponents, but it requires the full Johnson-ball partner rate
  `exp{nV+o(n)}`.  An exponentially thinned local cluster is insufficient.

Thus (10.1294) is much stronger than merely asking for exponentially many
partners.  The exact weighted condition (10.1292) remains viable because it
can distribute binomial weight over the complete intersection histogram.

## 2. Correction to the displayed Wave 53 parameter example

For the parameters displayed after (10.1295),

```math
p=3/5,\qquad \alpha=3/10,\qquad \beta=9/20,
```

the exact rates are

```math
E=0.6074926058936316\ldots,\qquad
V=0.6020263820344778\ldots,\qquad
E-V=0.0054662238591538\ldots.
\tag{R54C.10}
```

The equality core for this `(p,beta)` would instead be `alpha_*=0.36`.
Consequently the displayed `alpha=0.3` one-threshold hypothesis is
asymptotically impossible.  This does not affect the exact validity of
(10.1294) as an implication; it corrects the assessment that these parameters
give a potentially usable cluster input.

There is also a finite exact check.  On multiples of twenty, set
`m=3n/5`, `ell=3n/10`, and `s=9n/20`.  At `n=160`, already

```math
\frac{b_\ell-\binom m\ell}{\binom s\ell}
>\sum_{j=s}^{m-1}\binom mj\binom{n-m}{m-j}.
\tag{R54C.11}
```

The left side is the `epsilon=0` demand, so every positive saved `epsilon`
only strengthens the impossibility.  The checker finds `n=160` as the first
such multiple among `20,40,...,1000`; the positive rate gap proves the same
for all sufficiently large multiples.

## 3. Every child ground has a positive-threshold exchange star

The maximal-selector assumption in (10.1271)--(10.1273) can be removed if one
pays an `O(n)` deficit threshold.  Write

```math
q_S=Q(A[S]),\qquad |S|=m,\qquad T=S^c,\qquad |T|=k.
```

Fix any child ground `y` of `A[S]`, orient it by `sigma` so that
`sigma y^T A[S]y=q_S`, and complete it arbitrarily to `z`.  For `i in S`
put

```math
r_i=\sigma y_i\sum_{u\in S\setminus\{i\}}a_{iu}y_u.
```

One-spin optimality and summation give

```math
0\le r_i\le q_S/2,\qquad \sum_{i\in S}r_i=q_S.
\tag{R54C.12}
```

Adjacent child caps have the **Verified universal Lipschitz bound**

```math
\boxed{|q_{S-i+j}-q_S|\le2(m-1).}
\tag{R54C.13}
```

To prove it, use the common `(m-1)`-set `U=S-i`.  Adding one vertex can only
add `2(m-1)` in absolute energy, while choosing the added spin to agree with
the sign of a ground energy on `U` shows `q_U<=q_S,q_{U+j}`.

For a port `(i,j)`, put

```math
h_{ji}=\sum_{u\in S\setminus\{i\}}a_{ju}y_u.
```

Choose between `z` and `z^j` so that the incoming spin maximizes the absolute
energy on `S-i+j`; call the chosen center `c_{ij}`.  Since
`q_S-2r_i>=0`, direct expansion gives the **Verified exact deficit formula**

```math
D_{c_{ij}}(S-i+j)
=q_{S-i+j}-q_S+2(r_i-|h_{ji}|).
\tag{R54C.14}
```

In particular,

```math
D_{c_{ij}}(S-i+j)
\le[q_{S-i+j}-q_S]_++2r_i.
```

Summing (R54C.13)--(R54C.14), using (R54C.12) and the trivial
`q_S<=m(m-1)`, proves the **Verified arbitrary-seed port budget**

```math
\boxed{
\sum_{i\in S,j\in T}D_{c_{ij}}(S-i+j)
\le4mk(m-1).}
\tag{R54C.15}
```

At least `mk/2` ports therefore have deficit at most `8(m-1)`.  They are
assigned among only the `k+1` centers `z,z^j`, and distinct ports give
distinct exchanged selectors.  Hence some such center has at least

```math
\boxed{\frac{mk}{2(k+1)}}
\tag{R54C.16}
```

distinct `8(m-1)`-favorable one-exchange selectors.  The base selector `S`
remains an exact ground at every assigned center, and every pair counted in
(R54C.16) is off-diagonal.  Also

```math
R_2(c_{ij})\le\{\sqrt{R_2(z)}+2\sqrt n\}^2.
\tag{R54C.17}
```

The threshold `8(m-1)=O(n)` is below the active restriction tolerance scale,
but importing positive-deficit families into a particular spectral cap still
requires the usual tolerance bookkeeping; (R54C.16) is not asserted as a
zero-deficit theorem.

On the stored exact `A10` minimizer at `m=6`, exhaustive enumeration checks
all 8,320 child-ground/completion incidences.  The minimum assigned star at
threshold forty is six, the largest adjacent cap jump is eight (against the
universal bound ten), and every formula and budget above holds exactly.  This
is a finite mechanism check, not asymptotic evidence.

## 4. Exact aggregation and why it still does not make a global cluster

Let `I_R` be any collection of exact-ground incidences `(S,z)` with
`R_2(z)<=R`, and apply (R54C.16) to every incidence.  For a fixed base `S`
and output center `c`, at most `k+1` inputs can map to `(S,c)`, because the
input is one of `c,c^j` with `j notin S`.  Therefore the output contains at
least

```math
\frac{|I_R|}{k+1}
\tag{R54C.18}
```

distinct `(base selector, center)` flags.  Multiplying by the star size gives
the **Verified all-seed ordered-pair bound**

```math
\boxed{
\#\{(c,S,U):S\ne U,\ |S\cap U|=m-1,\ S\text{ is a flagged ground},
\ D_c(U)\le8(m-1)\}
\ge\frac{mk}{2(k+1)^2}|I_R|,}
\tag{R54C.19}
```

up to harmless integer rounding, with every output center under the enlarged
cap in (R54C.17).  No self-loop is counted.

This is the strongest conclusion available from star aggregation alone.  It
still gives only `Theta(n)` partners for each flagged base.  Repeated stars
may occupy different centers, and (R54C.18) supplies no signing-specific
reason for their bases to be mutually close.

There is a scalable **abstract counting obstruction**.  Let
`u=m-s=(p-beta)n+O(1)` and take `beta` sufficiently close to `p` that the
Johnson-ball rate `V(p,beta)<p log 2`.  A greedy constant-weight code with
minimum exchange distance greater than `u+2` has size at least

```math
\frac{\binom nm}{\sum_{d\le u+2}\binom md\binom kd}
=\exp\{n(H(p)-V(p,\beta))+o(n)\}.
\tag{R54C.20}
```

This is exponentially larger than the average base-flag rate obtainable by
putting one chosen ground for every selector through all its completions,

```math
\frac{\binom nm\,2^{n-m}}{(k+1)2^{n-1}}
=\exp\{n(H(p)-p\log2)+o(n)\}.
\tag{R54C.21}
```

If `H(p)-p log 2>0`, choose a subset of the code at the latter exponential
size and attach the full one-exchange star to every codeword.  If that rate is
nonpositive, (R54C.21) does not force even exponentially many base flags at
one center.  Stars around distinct codewords have no cross-pair with
intersection at least `s`, by the triangle inequality for exchange distance.
Each member therefore has only `O(n^2)` such partners, despite (in the
positive-rate case) the exponentially large supply of bases and a star at
every base.

This construction is exact as a Johnson-incidence model and falsifies a proof
using only total all-seed mass, the number of centers, and local-star axioms.
It is **not** claimed to be the ground-incidence system of an exact signing;
a future proof may still exploit compatibility of restrictions, row geometry,
or exact-minimizer structure to rule it out.

## 5. Research judgment and exact remaining step

All-seed aggregation improves the port theorem by removing the
maximal-selector restriction and yields the exact global count (R54C.19), but
it does not force the large-core partners needed by (10.1292).  Moreover,
(10.1294) should no longer be pursued at arbitrary linear parameters: it is
asymptotically impossible unless `alpha=alpha_*` in (R54C.9).

A live all-seed continuation must prove one of the following genuinely
signing-specific statements at a controlled project-row cap:

1. on the tuned equality surface, a saved biased mass of centers has nearly
   the full Johnson-ball exponential partner rate; or
2. directly in the exact self-loop-free identity (10.1292), the complete
   binomial-weighted intersection histogram exceeds its `lambda_2` threshold
   by `exp{-O(n^(3/4-c))}` without reducing to a single threshold.

Neither follows from the number of seeds, all completions, or polynomial
exchange stars.  Controlled-cap compatibility and a genuinely global
intersection theorem remain the exact missing inputs.
