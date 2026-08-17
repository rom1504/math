# Deterministic de Bruijn response versus anticipatory proof memory

Status: proof draft, independently derived by the weighted-carrier agent and
sent for adversarial audit.

This example answers the live converse question negatively in its strongest
form.  Weighted anticipatory supports admit an exact finite mean-payoff-game
characterization, but their smallest size is not determined by the semantic
response image: the latter can be one point while support-certificate memory
is exponential.

## 1. Exact support-game characterization

Use the notation of Theorem 17.1s.  For a finite coarse edge
`c=(a,e,b)` and nonempty supports `K subset pi^(-1)(a)`,
`L subset pi^(-1)(b)`, define

```math
d_c(K,L)=\max_{j\in L}
 \left[S_e(a,b)-\max_{i\in K}\widehat T_e(i,j)\right]_+,      \tag{DB.1}
```

with value `+infinity` if an endpoint has no finite predecessor.  This is
the least shortfall making the backward-surjective inclusion (17.7ay) true.

### Theorem DB.1 (finite support carriers are mean-payoff strategies)

Fix nonnegative budgets `beta_c`.  The following are equivalent:

1. a finite Theorem-17.1s support-plus-potential carrier exists;
2. finite nonempty support families and a successor selector can be chosen so
   that every directed cycle of the selected support graph has

   ```math
   \sum_c(d_c(K,\sigma(K,c))-\beta_c)\le0;                    \tag{DB.2}
   ```

3. in the finite game where the environment chooses the next coarse edge and
   the controller chooses the successor support, the controller can keep the
   accumulated adjusted shortfall uniformly bounded **above**, from every
   required selected start support;
4. those starts lie in the controller's threshold-zero mean-payoff winning
   region.

Finite memory gives no advantage over a positional support strategy.  With
a budget of at most `N` supports, the optimal uniform toll is therefore

```math
\beta_N=\min_{|W|\le N,\sigma}
 \max_{\Gamma\ {\rm selected\ cycle}}
 {1\over|\Gamma|}\sum_{\gamma\in\Gamma}d_\gamma.             \tag{DB.3}
```

#### Proof

For a fixed selector, the potential condition is

```math
d_c(K,L)-\beta_c\le\psi(L)-\psi(K).                           \tag{DB.4}
```

It implies (DB.2) by telescoping.  Conversely, if the selected graph has no
positive adjusted cycle, define `psi(L)` as the maximum adjusted weight of a
selected path ending at `L`, after fixing one zero-potential source in each
reachable component.  Deleting nonpositive cycles makes this maximum finite,
and appending one edge gives (DB.4).  Uniform upper path bounds and absence of
a positive selected cycle are equivalent.  Positional determinacy of finite
mean-payoff games proves that finite controller memory gives no advantage.
For edge-dependent budgets, replace the uniform toll in (DB.3) by adjusted
edge costs; the scalar cycle-mean formula is the uniform-budget case.
`square`

This is an exact characterization of the **17.1s certificate architecture**,
not of arbitrary wordwise spectral equivalence.  A failed support strategy
can be pumped against that strategy while a different future-dependent raw
witness may still produce the correct scalar response.

## 2. Deterministic de Bruijn shift

Let `E` have `q>=2` letters, let `m>=1`, `C>0`, let `I=E^m`, and set

```math
F_e(s_1...s_m)=s_2...s_me,

T_e(s,t)=\begin{cases}0,&t=F_e(s),\\-C,&t\ne F_e(s),\end{cases} \tag{DB.5}
```

with scalar coarse system `S_e=0`.

### Theorem DB.2 (trivial response, exponential exact support memory)

Every word has spectral response zero, but every exact anticipatory-support
carrier has at least `q^m` states.  Singleton supports attain this bound.

#### Proof

For a nonempty word `w`, take an `m`-window of the bi-infinite periodic word
`w^infinity`.  It is a fixed point of `F_w`, so `T_w` has a zero diagonal
entry.  No entry is positive, hence `rho(T_w)=0`.

First note that any gauge compatible with exact upper domination is constant.
Indeed, every deterministic edge gives `h(s)<=h(F_e(s))`, and the directed
de Bruijn graph is strongly connected.  Thus the tight and nontight gaps
remain exactly `0` and `C`.

Choose a carrier support of maximal potential.  At `beta=0`, its successor
has no larger potential, while `d<=psi(K')-psi(K)` and `d>=0`; hence `d=0`
and the successor also has maximal potential.  Inductively, every update from
this start is tight and obeys

```math
K'\subseteq F_e(K).                                            \tag{DB.6}
```

For any nonempty `K` and every `u in E^m`, the deterministic image `F_u(K)`
is the singleton `{u}`.  Following the `q^m` length-`m` words from the chosen
maximal-potential state therefore reaches the `q^m` distinct singleton
supports.  Hence the carrier has at least `q^m` states.  The raw singleton
supports, updated by `F_e`, form an exact carrier. `square`

This is not the repository's free-tail de Bruijn relation, whose appended
symbol is unconstrained and whose full support survives every letter.  Here
the input deterministically supplies the appended symbol.

### Theorem DB.3 (sharp approximate support-certificate tradeoff)

If a uniform-toll carrier has `N` states, then for every `ell<=m`,

```math
\ell\beta<C\quad\Longrightarrow\quad N\ge q^\ell.             \tag{DB.7}
```

Conversely, for every `0<=L<m` there is a carrier with

```math
N_L=1+q+\cdots+q^L={q^{L+1}-1\over q-1},
\qquad \beta_L={C\over L+1}.                                  \tag{DB.8}
```

Thus, for `1<=N<q^m`, before the discontinuous exact threshold,

```math
\beta_N=\Theta\left({C\over1+\log_qN}\right).                 \tag{DB.9}
```

#### Proof

Start at a carrier support maximizing its potential.  Along a word `u` of
length `ell`, (17.7az) telescopes to

```math
\sum_{s=1}^\ell d_s
\le\ell\beta+\psi(K_u)-\psi(K_0)\le\ell\beta<C.               \tag{DB.10}
```

Every shortfall is nonnegative, so every step has shortfall below `C` and
must use only the deterministic zero edge.  Therefore
`K_u subseteq F_u(K_0)`.  For distinct length-`ell` words these images lie
in disjoint suffix cylinders, proving (DB.7).

For the upper bound, index supports by words `u in E^(<=L)` and put

```math
K_u=\{s\in E^m:s\text{ ends in }u\}.
```

At depths below `L`, update `u` to `ue` at zero shortfall.  At depth `L`,
reset to `K_emptyset=I`; because every nontight entry in (DB.5) has weight
`-C`, this reset has shortfall `C`.  With

```math
\psi(K_u)=-{|u|C\over L+1},
```

the right side of the potential inequality is zero on an ordinary update
and exactly `C` on a reset.  This proves (DB.8).  Combining (DB.7) and
(DB.8) gives (DB.9). `square`

## 3. Consequence for the proposed unified law

The example has

```math
\text{semantic scalar response entropy}=0,
\qquad
\text{exact anticipatory proof memory}=m\log_2q.               \tag{DB.11}
```

Together with the existing free-tail de Bruijn and width-two Ising examples,
it realizes three distinct patterns:

| model | scalar response | anticipatory support | forward/path state |
|---|---:|---:|---:|
| free-tail de Bruijn | 1 | 1 | `2^m` |
| deterministic de Bruijn | 1 | `q^m` | 1 |
| width-two Ising | scalar | 2 | 4 |

Therefore a formula equating reusable semantic complexity with response
entropy times support-carrier size is false: carrier size can measure the
complexity of one proof architecture rather than an observable response.
Forward source-total and anticipatory target-surjective certificates are
incomparable: the free-tail and width-two examples favor the latter, while
the deterministic shift has the exact one-state rowwise lift
`s -> F_e(s)` but an exponential anticipatory carrier.
The missing third resource is only intrinsic when a declared rooted probe or
future interaction exposes the hidden witness phase.

The next theorem should characterize the smallest probe family that converts
support memory into an actual response packing.  That is more discriminating
than another twins-property variant: the deterministic example is already a
deterministic weighted automaton, so classical twinning is vacuous.
