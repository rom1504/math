# Thermal activation and the replica budget of a rare phase

**Status.**  Rigorous, scoped adversarial-statistical-mechanics theorem.
This note gives one exact rare-phase carrier and one sharp no-go theorem for
bounded-temperature observation.  It is deliberately separate from the
spectral-spike benchmark: the missing object here is a set of exponentially
small Gibbs mass, not a zero-density spectral atom.  It does not address the
dense signing convergence problem.

## 1. A deterministic marked-phase landscape

Let `Omega_n` be a finite set, let `P_n` be a nonempty proper subset, and put

```math
q_n={|P_n|\over|\Omega_n|},
\qquad
H_n^{(0)}(x)=0,
\qquad
H_n^{(1)}(x)=E_n\mathbf 1_{P_n}(x).                    \tag{TA.1}
```

The normalized partition function and Gibbs measure are

```math
Z_{n,s}(\beta)={1\over|\Omega_n|}
 \sum_x e^{\beta H_n^{(s)}(x)},
\qquad
\mu_{n,s,\beta}(x)
 ={e^{\beta H_n^{(s)}(x)}\over
   |\Omega_n|Z_{n,s}(\beta)}.                           \tag{TA.2}
```

For the marked landscape, write

```math
p_n(\beta)=\mu_{n,1,\beta}(P_n)
 ={q_ne^{\beta E_n}\over1-q_n+q_ne^{\beta E_n}}.       \tag{TA.3}
```

The pair

```math
\Theta_n=(q_n,E_n)                                      \tag{TA.4}
```

is an exact presented state for the declared queries consisting of the
uniform pressure and the Gibbs mass of the marked phase:

```math
Z_{n,1}(\beta)=1-q_n+q_ne^{\beta E_n},
\qquad p_n(\beta)={q_ne^{\beta E_n}\over Z_{n,1}(\beta)}.
                                                               \tag{TA.5}
```

It has an exact conjunctive product.  Given two marked landscapes, define

```math
(\Omega,P,E)\mathbin\wedge(\Omega',P',E')
=(\Omega\times\Omega',P\times P',E+E').                \tag{TA.6}
```

This is the natural operation for a hard reward obtained only when both
independent constraint systems are satisfied.  Then

```math
\boxed{\Theta\mathbin\wedge\Theta'=(qq',E+E').}         \tag{TA.7}
```

Thus the extensive rarity costs `J=-log q` and rewards add.  For factors of
sizes `n,m`, the product's intensive rate is the weighted average
`(nI+mI')/(n+m)`, not `I+I'`.  The state is strictly smaller than the marked
set when its geometry is not queried.  This exact carrier is not claimed for
ordinary additive Hamiltonian composition: that operation creates the two
intermediate phases and requires their phase table (or energy enumerator).

The same carrier exists over a nonflat background.  Let `K:Omega->R`, let
`mu_(K,beta)` be its Gibbs law, fix `E>=0`, and boost a marked event by

```math
K^+(x)=K(x)+E\mathbf1_P(x),
\qquad q_K(\beta)=\mu_{K,\beta}(P).                     \tag{TA.7a}
```

### Lemma TA.0 (exact rare-event tilt and background composition)

For every finite background and every `beta>=0`,

```math
{Z_{K^+}(\beta)\over Z_K(\beta)}
=1-q_K(\beta)+q_K(\beta)e^{\beta E},                   \tag{TA.7b}
```

```math
p_K(\beta):=\mu_{K^+,\beta}(P)
={q_K(\beta)e^{\beta E}\over
  1-q_K(\beta)+q_K(\beta)e^{\beta E}},                 \tag{TA.7c}
```

and

```math
\|\mu_{K^+,\beta}-\mu_{K,\beta}\|_{TV}
=p_K(\beta)-q_K(\beta).                                \tag{TA.7d}
```

Consequently the relative marked-phase state

```math
\Theta(K,P,E)=(q_K(\mathord\cdot),E)                    \tag{TA.7e}
```

answers the complete relative-pressure and marked-mass curves.  If `F>=0`,
`K\oplus L(x,y)=K(x)+L(y)` and the marked composition is
`P times Q` with reward `E+F`, then

```math
\boxed{
\Theta(K\oplus L,P\times Q,E+F)
=(q_K(\mathord\cdot)q_L(\mathord\cdot),E+F).}           \tag{TA.7f}
```

If absolute rather than relative pressures are declared, adjoining the base
curve `Z_K` gives the equally exact product update
`Z_(K\oplus L)=Z_KZ_L`.

For arbitrary `K`, the component `beta -> q_K(beta)` is an entire function
and need not have a finite description.  The finite two-coordinate collapse
occurs in the flat membership models used in Theorem TA.1 and Section 3; no
general compression claim is hidden in (TA.7e).

Moreover, for any adaptive experiment using a deterministic (or almost
surely bounded) budget of `K_samp` Gibbs samples at
temperatures in a set `mathcal B`, its two transcript laws for `K` and `K^+`
satisfy

```math
\|\mathcal L_{K^+}-\mathcal L_K\|_{TV}
\le K_{samp}\sup_{\beta\in\mathcal B}
       (p_K(\beta)-q_K(\beta)).                         \tag{TA.7g}
```

Conversely, at any fixed `beta`, phase counting distinguishes the two laws
whenever `K_samp p_K(beta)->infinity` and
`q_K(beta)/p_K(beta)->0`.

#### Proof

Split each partition sum over `P` and `P^c`; this gives (TA.7b)--(TA.7c).
The two laws have the same conditional distribution inside each of those
two cells, so their total-variation distance is the difference of their
cell masses, proving (TA.7d).  Gibbs laws of additive backgrounds factor,
so the base mass of `P times Q` is `q_K q_L`, proving (TA.7f).  Sequential
maximal coupling until the first differing sample proves (TA.7g).  For the
converse, count visits to `P` and apply Chernoff under `K^+` and Markov under
`K`, exactly as in (TA.19)--(TA.20) below. `square`

The conjunctive boosted composite in (TA.7f) is deliberately not the sum
`K^++L^+`: the latter has two intermediate marked phases.  This distinction
is one of the theorem's internal scope checks.

## 2. Sharp thermal-observation theorem

Assume

```math
-{1\over n}\log q_n\longrightarrow I\in(0,\infty),
\qquad {E_n\over n}\longrightarrow\delta\in(0,\infty). \tag{TA.8}
```

An observer is told the two candidate presentations but is given only Gibbs
samples from the unknown candidate `s in {0,1}`.  At round `t` it may choose
an inverse temperature `beta_t in [0,B]` as an arbitrary randomized function
of all previous samples, and then receives an independent sample from
`mu_(n,s,beta_t)`.  Thus this is stronger than observing any predeclared
finite collection of overlaps or replica statistics: those are
post-processings of the full configurations.  Energy values evaluated under
the two known candidate Hamiltonians are also post-processings; an oracle
that separately reveals which unknown Hamiltonian supplied an energy label
would be an additional observation channel and is not included.

### Theorem TA.1 (rare-phase activation equals the replica exponent)

Fix `B>0` with

```math
a_B:=I-B\delta>0.                                      \tag{TA.9}
```

Then the following hold.

1. **The whole bounded pressure curve is exponentially blind.**  Uniformly
   for `0<=beta<=B`,

   ```math
   0\le {1\over n}\log Z_{n,1}(\beta)
   -{1\over n}\log Z_{n,0}(\beta)
   \le {1\over n}\exp\{-na_B+o(n)\}.                  \tag{TA.10}
   ```

   At `beta=B` the right exponential rate is attained:

   ```math
   {1\over n}\log Z_{n,1}(B)
   = {1\over n}\exp\{-na_B+o(n)\}.                   \tag{TA.11}
   ```

2. **Every subcritical adaptive replica experiment is blind.**  Let `K_n`
   be the total number of samples.  If

   ```math
   \limsup_n{1\over n}\log K_n<a_B,                   \tag{TA.12}
   ```

   then the total-variation distance between the complete adaptive
   transcripts under `s=0` and `s=1` tends to zero.  Hence every statistic
   of those samples—including every sampled multi-overlap array—has the same
   asymptotic law under the two landscapes.  With a uniform prior on `s`, the
   mutual information between `s` and the transcript tends to zero.

3. **The exponent is sharp.**  If

   ```math
   \liminf_n{1\over n}\log K_n>a_B,                   \tag{TA.13}
   ```

   then `K_n` samples at the single temperature `B`, followed only by
   counting how many lie in `P_n`, distinguish `s=0` from `s=1` with error
   tending to zero.

4. **The zero-temperature response remains macroscopically separated.**

   ```math
   {1\over n}\max_xH_n^{(1)}(x)
   -{1\over n}\max_xH_n^{(0)}(x)\longrightarrow\delta. \tag{TA.14}
   ```

   Here is a precise robust-pressure version of the associated minimax
   statement.  Put

   ```math
   F_{n,s}(\beta)={1\over n}\log Z_{n,s}(\beta),
   \quad
   \Delta_n^P=\|F_{n,1}-F_{n,0}\|_{L^\infty[0,B]},
   \quad
   d_n(B)=\|\mu_{n,1,B}-\mu_{n,0,B}\|_{TV}.             \tag{TA.14a}
   ```

   Suppose a pressure oracle may adversarially return any curve within
   `eta_n` in `L^infinity[0,B]` of the true curve.  If

   ```math
   2\eta_n\ge\Delta_n^P,
   \qquad K_nd_n(B)\longrightarrow0,                    \tag{TA.14b}
   ```

   then the combined pressure-plus-replica experiment has minimax absolute
   error at least `delta/2-o(1)` for the normalized maximum.  Exact finite
   pressure values are not covered: they distinguish the candidates.

The quantity `a_B=I-B delta` is therefore simultaneously:

* the residual free-energy cost of activating the rare phase;
* the robust pressure-resolution exponent in (TA.14b); and
* the exponential number of Gibbs replicas needed to see the phase.

#### Proof

For `s=0`, `Z_(n,0)=1` and every Gibbs measure is uniform.  Equations
(TA.3)--(TA.5) are direct counting.  From (TA.8), uniformly on `[0,B]`,

```math
q_n(e^{\beta E_n}-1)\le\exp\{-na_B+o(n)\}.             \tag{TA.15}
```

Using `log(1+u)<=u` proves (TA.10).  At `B`, the expression in (TA.15) is
`exp(-na_B+o(n))`; since it tends to zero, `log(1+u)=u(1+o(1))`, proving
(TA.11).

Both Gibbs laws are uniform conditional on `P_n` and conditional on its
complement.  Their total-variation distance is therefore exactly

```math
d_n(\beta)
=\|\mu_{n,1,\beta}-\mu_{n,0,\beta}\|_{TV}
=p_n(\beta)-q_n.                                        \tag{TA.16}
```

It is increasing in `beta`, and (TA.3), (TA.8), and (TA.9) give

```math
d_n(B)=\exp\{-na_B+o(n)\}.                              \tag{TA.17}
```

Couple the two adaptive experiments sequentially.  As long as their
histories agree, the observer chooses the same next temperature, and a
maximal coupling makes the next samples disagree with probability at most
`d_n(B)`.  A union bound gives

```math
\|\mathcal L_1(\mathrm{transcript})
 -\mathcal L_0(\mathrm{transcript})\|_{TV}
\le K_nd_n(B).                                          \tag{TA.18}
```

This tends to zero under (TA.12).  Data processing gives the assertion for
all replica statistics.  For a binary prior, vanishing total variation of
the two conditional transcript laws also forces their Jensen--Shannon
divergence, hence the mutual information, to vanish.

For the converse, sample at `B` and let `N_n` count visits to `P_n`.  Under
`s=0` it is binomial with mean `K_nq_n`; under `s=1` it is binomial with mean
`K_np_n(B)`.  Assumption (TA.13) and (TA.3) give

```math
K_np_n(B)\longrightarrow\infty,
\qquad {q_n\over p_n(B)}=\exp\{-nB\delta+o(n)\}\longrightarrow0.
                                                               \tag{TA.19}
```

Chernoff's bound under `s=1` and Markov's inequality under `s=0` show

```math
\Pr_1\{N_n<K_np_n(B)/2\}\to0,
\qquad
\Pr_0\{N_n\ge K_np_n(B)/2\}
\le {2q_n\over p_n(B)}\to0.                             \tag{TA.20}
```

This proves (3).  Equation (TA.14) is immediate from (TA.1).  For
(TA.14b), the two `eta_n`-balls around the pressure curves intersect: the
oracle may return their midpoint under either candidate.  The remaining
transcript laws have total variation at most `K_nd_n(B)=o(1)`.  Two
parameters separated by `delta+o(1)` whose observation laws have total
variation `o(1)` have minimax absolute estimation risk at least
`delta/2-o(1)` (integrate the pointwise triangle inequality against the
common part of the two laws). `square`

### Boundary of the theorem

The subcritical condition is sharp in another elementary sense.  If
`B delta>I`, then `p_n(B)->1` while `q_n->0`, so one temperature-`B` sample
distinguishes the candidates with probability tending to one.  At
`B delta=I`, subexponential factors decide the answer and there is no
universal limit under (TA.8) alone.  More generally the exact statements are
`K_nd_n(B)->0` for blindness and
`K_np_n(B)->infinity`, `q_n/p_n(B)->0` for phase-counting recovery; slogans
of the form `exp(n(a_B plus/minus o(1)))` do not resolve the critical window.
Thus the theorem is not the false claim that finite temperature always loses
the ground state; it identifies the activation threshold and the observation
cost below it.

## 3. Concrete code and CSP benchmarks

Let `Omega_n=F_2^n` and let `C_n` be any binary linear code with

```math
{\dim C_n\over n}\longrightarrow R\in[0,1).
```

Take `P_n=C_n` and reward code membership by `E_n=n delta`.  Then

```math
I=(1-R)\log2.                                            \tag{TA.21}
```

For every

```math
B<{(1-R)\log2\over\delta},                              \tag{TA.22}
```

the sharp exponential replica threshold is

```math
\exp\{n((1-R)\log2-B\delta+o(1))\}                      \tag{TA.23}
```

in the following precise sense: every strict smaller exponent is blind and
every strict larger exponent permits consistent phase counting.  The
critical subexponential factor is not specified.  This holds even though the
normalized ground-state reward is `delta`.  Direct products of codes are
exactly the conjunctive composition (TA.6), so the two-coordinate state
`(codimension, reward)` adds and answers all declared thermal queries without
retaining a generator matrix, codeword list, or coset geometry.  This is a
membership-bonus model, not a covering-radius, coset-distance, or decoding-
energy theorem.

The same statement applies verbatim to a deterministic CSP family whose
satisfying set has density `exp(-nI+o(n))`, when the Hamiltonian pays a
macroscopic bonus only for satisfying every constraint.  Disjoint
conjunction multiplies satisfying densities and adds declared bonuses.  No
random-disorder assumption is used in either application.

## 4. What is classical and what this theorem adds

The identities (TA.3), the Gibbs likelihood ratio, maximal coupling,
Chernoff bounds, and the two-point minimax inequality are classical.  The
carrier `(q,E)` for a two-level energy histogram is also elementary.  The
entropy--reward balance and finite phase-count multiplication already appear
in repository Proposition 27.2 and Theorem 27.3, while the entropy-tilted
bridge audit obtains the analogous `log(1/p)` cost by changing measure over
bridges.  The genuine project-specific increment here is the exact adaptive
transcript theorem and its sharp Gibbs-replica threshold.  It synthesizes
four scales in one operational response statement:

```math
\boxed{
\text{rarity }I
-\text{ thermal reward }B\delta
=\text{ pressure-resolution exponent}
=\text{ replica-budget exponent}.}                       \tag{TA.24}
```

This is not the posterior-width theorem in different notation.  That theorem
lower-bounds bits needed to encode one of many separated response vectors.
Here the hidden family has only two members; the obstruction is that the
declared Gibbs query assigns exponentially small mass to the response-
distinguishing set.  Nor is this merely Proposition 27.2: that proposition
shows equality of limiting bounded-temperature pressures, whereas TA.1
shows that even the full adaptive Gibbs experiment and every replica
observable remain blind up to a sharp exponential sample threshold.

The spectral-spike benchmark has a different carrier and mechanism.  A
rank-one perturbation has vanishing empirical spectral mass but creates an
outlier through a resolvent equation.  Here an exponentially rare
configuration phase is activated when its Boltzmann reward pays its entropy
cost, and the exact state composes under hard conjunction.

## 5. Internal falsifiers and research relevance

The scope can be killed in three explicit ways.

1. **Raise the temperature budget.**  Once `B delta>I`, one sample exposes
   the marked phase.
2. **Raise the replica budget.**  Below the thermal transition, the exact
   phase-count criterion is `K_np_n(B)->infinity`.  A claim of blindness
   beyond that scale is false; the critical subexponential window is not
   fixed by (TA.8).
3. **Allow geometric future contexts.**  If `P,P' subset Omega` are disjoint
   and equally large, they have the same `(q,E)`.  The future reward
   `V=E' 1_P` gives

   ```math
   \max(H_P+V)=E+E',
   \qquad
   \max(H_{P'}+V)=\max\{E,E'\},                         \tag{TA.25}
   ```

   so the carrier fails by `min(E,E')`.  It is sufficient only for the
   two-phase thermal queries declared in (TA.5) and conjunctive composition,
   not arbitrary labelled couplings.  Ordinary additive composition likewise
   creates intermediate phases and forces the state to expand to a phase
   enumerator.

The theorem gives a precise rare-event complement to contextual response
geometry: a small realizable response image can still be operationally
invisible when every allowed bounded-temperature query puts exponentially
small mass on its exposed face.  It suggests that any statistical-mechanics
state intended to preserve zero-temperature response must specify not only
which extremal phases exist, but also their activation costs relative to the
largest allowed temperature and the available query/replica budget.

**Portfolio judgment:** keep this as a scoped Level-3 theorem/no-go.  It is a
strict, composable state for marked code/CSP phases and a sharp operational
falsifier for bounded-temperature summaries.  It does not justify a new full
theory branch, because arbitrary labelled futures immediately require the
geometry of `P_n` and ordinary additive interactions create a growing phase
table.
