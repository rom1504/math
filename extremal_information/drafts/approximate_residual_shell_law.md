# Approximate residual shells: terminal recovery, compatibility flux, and a sharp rate barrier

**Status.**  Proof draft.  The finite identities and the sharp two-letter
example are checked by
[`../experiments/verify_approximate_residual_shell_law.py`](../experiments/verify_approximate_residual_shell_law.py).

This note addresses the approximate question left after Theorem 17.1u.  Its
main conclusion is deliberately asymmetric:

* a bounded-delay approximate row residual gives a **depth-uniform rooted
  terminal-response theorem**;
* it does **not** give a depth-uniform additive or spectral-response theorem;
* even exact projective resets can hide a scalar compatibility cocycle whose
  error is paid once per switch.

Thus projective forgetting and scalar-toll synchronization are different
resources.  A greatest support core is a third, optional witness resource;
its failure is not used as evidence of scalar drift.

## 1. Conventions

All matrices are finite real max-plus matrices on a finite set `I`.  Products
act on row vectors from the right.  For a matrix `P` and a row profile `p`,
write

```math
\operatorname{rad}(P;p)
 :=\inf_{a\in\mathbb R^I}\max_{i,j}|P(i,j)-a_i-p_j|.       \tag{ARS.1}
```

This is a two-sided entrywise version of approximate max-plus row rank one.
For a row `u` and a terminal field `z`, let

```math
\mathcal R_z(u):=\max_j(u_j+z_j)-\max_j u_j.               \tag{ARS.2}
```

The subtraction removes the accumulated scalar amplitude.  It is exactly
the rooted terminal response of the normalized row.

For projective profiles use

```math
d_\infty([p],[q]):=\inf_{b\in\mathbb R}\|p-q-b\mathbf1\|_\infty.
                                                                    \tag{ARS.3}
```

## 2. A last-window theorem with no accumulated terminal error

### Theorem ARS.1 (bounded-delay residual-shell recovery)

Fix `D>=1`.  Suppose that for each legal word `v` of length `D` there is a
profile `p_v` with

```math
\operatorname{rad}(T_v;p_v)\le\epsilon.                    \tag{ARS.4}
```

Let `C` be a finite profile code and choose `r_c`, together with a map
`chi:v\mapsto c`, such that

```math
d_\infty([p_v],[r_{\chi(v)}])\le\eta.                       \tag{ARS.5}
```

If a word `w` of length at least `D` ends in `v`, then for every initial raw
state `i` there is a scalar `b_{w,i}` satisfying

```math
\max_j|T_w(i,j)-b_{w,i}-r_{\chi(v)}(j)|
 \le \epsilon+\eta.                                         \tag{ARS.6}
```

Consequently, for every terminal field `z`,

```math
|\mathcal R_z(T_w(i,\cdot))-\mathcal R_z(r_{\chi(v)})|
 \le 2(\epsilon+\eta),                                     \tag{ARS.7}
```

uniformly in the depth, the prefix before `v`, the initial state, and the
size of `z`.

The runtime state may always store the legal length-`D` suffix, so its
recurrent cardinality is at most the number of legal `D`-words.  It may store
only `c` when `chi` is a right congruence for the suffix shift:

```math
\chi(\operatorname{suf}_D(vc))
   =\delta(\chi(v),c).                                      \tag{ARS.8}
```

Shorter prefixes are a finite transient.  No path selector or active-cell
language is present in this state.

#### Proof

Write the factorization supplied by (ARS.4) as

```math
T_v(k,j)=a_k+p_v(j)+E_{kj},\qquad |E_{kj}|\le\epsilon.       \tag{ARS.9}
```

If `w=uv`, then

```math
T_w(i,j)
=\max_k\{T_u(i,k)+a_k+p_v(j)+E_{kj}\}.
```

Putting `b_{w,i}=\max_k(T_u(i,k)+a_k)` leaves an error in
`[-epsilon,epsilon]`.  Re-gauging (ARS.5) proves (ARS.6).
Each of the two maxima in (ARS.2) changes by at most
`epsilon+eta`, which proves (ARS.7). `square`

This is not Theorem 17.2 in disguise.  The maps may have Hilbert Lipschitz
coefficient one, there is no multiplicative contraction of nearby states,
and no rounding error is injected at every step.  The conclusion is possible
because the declared query removes the scalar amplitude and only the final
row profile is observed.  Conversely, it says nothing about a sum of scalar
tolls.

## 3. Exact rank-one switching exposes the missing scalar object

The obstruction is already present before any approximate dynamics.

### Proposition ARS.2 (directed compatibility law)

Suppose every generator on a finite directed legal-word presentation is
exactly max-plus row rank one,

```math
T_e(i,j)=a_e(i)+p_e(j).                                      \tag{ARS.10}
```

Define the directed compatibility table

```math
\varphi(e,f):=\max_j\{p_e(j)+a_f(j)\}.                       \tag{ARS.11}
```

More generally an exact rank-one product has the presented state
`(a_left,p_right,s)`, and composition closes by

```math
(a,p,s)\star(b,q,t)
   =(a,q,s+t+\max_j(p_j+b_j)).                                \tag{ARS.11a}
```

Thus the state retains two endpoint types and one scalar, not the full
matrix product.

For every nonempty cyclic word `w=e_1\cdots e_t`, with `e_{t+1}=e_1`,

```math
\rho(T_w)=\sum_{s=1}^t\varphi(e_s,e_{s+1}).                 \tag{ARS.12}
```

If `\widehat\varphi` is any proposed finite quotient table, its scalar
response has a depth-independent error on all cyclic words exactly when
the defect

```math
d(e,f)=\varphi(e,f)-\widehat\varphi([e],[f])                 \tag{ARS.13}
```

has zero sum on every directed cycle of the letter/quotient graph.  On each
strongly connected component this is equivalent to a potential

```math
d(e,f)=\psi(f)-\psi(e).                                      \tag{ARS.14}
```

Otherwise a violating cycle is a genuine pumpable scalar-response witness.

#### Proof

Multiplying rank-one matrices gives

```math
(a_e\otimes p_e)(a_f\otimes p_f)
=\varphi(e,f)+a_e\otimes p_f.                               \tag{ARS.15}
```

Iteration leaves the first left profile, the last right profile, and the
sum of interior compatibilities.  The max-plus eigenvalue of the resulting
rank-one matrix adds the closing compatibility `varphi(e_t,e_1)`, proving
(ARS.12).  Repeating a directed cycle proves necessity in the second claim.
Deleting cycles from an arbitrary path leaves a bounded simple path, proving
sufficiency.  The potential characterization is the finite cycle-coboundary
criterion. `square`

The cycle criterion is the specialization of Theorem 17.1l to the intrinsic
rank-one compatibility graph; it is recorded here to identify the relevant
quantity, not claimed as a new version of finite reward cohomology.  The
new point for the approximate-residual question is that a tiny projective
profile merge can change `varphi` by a tiny amount **at every composition**.
Projective reset does not contract that additive cocycle.

### Proposition ARS.2a (approximate rank-one block carrier)

Let `V` be a finite alphabet of legal blocks, all of the same original
length `D`, with concatenations declared by a finite directed graph, and
suppose

```math
\|T_v-a_v\otimes p_v\|_\infty\le\epsilon_v
\qquad(v\in V).                                               \tag{ARS.15a}
```

Put `phi(v,u)=max_j(p_v(j)+a_u(j))`.  For every legal cyclic block
word `W=v_1\cdots v_k`,

```math
\left|\rho(T_W)-\sum_{s=1}^k\varphi(v_s,v_{s+1})\right|
 \le\sum_{s=1}^k\epsilon_{v_s},
\qquad v_{k+1}=v_1.                                          \tag{ARS.15b}
```

More generally, let a proposed finite block quotient carry a compatibility
table `phihat`, and let

```math
\Delta=\max_C{1\over|C|}
 \left|\sum_{(v,u)\in C}
 [\varphi(v,u)-\widehat\varphi([v],[u])]\right|              \tag{ARS.15c}
```

over directed simple cycles of the legal block graph.  Its asymptotic
spectral distortion is at most

```math
{\max_v\epsilon_v+\Delta\over D}                             \tag{ARS.15d}
```

per original letter, up to a bounded transient when the word length is not
a multiple of `D`.  In the exact case `epsilon_v=0`, `Delta=0` is necessary
and sufficient for a depth-independent absolute error.

Indeed max-plus multiplication is one-Lipschitz in each factor in the
entrywise sup norm.  Replacing the `k` factors successively changes the
product, and hence its spectral radius, by at most `sum epsilon_v`.  Apply
ARS.2 to the rank-one reference product.  Cycle deletion gives (ARS.15c--d),
and repeated cycles give the exact-case converse.

This is the strongest generic scalar conclusion available from blockwise
approximate row rank one: a **rate** theorem.  A right-profile cover can make
the terminal state small, but the compatibility table also pairs each right
profile with the next block's left profile.  Unless that pairing descends to
the quotient modulo a cycle coboundary, the cover is not a reusable scalar
state.

### Theorem ARS.3 (sharp reset-versus-rate counterexample)

For every `delta>0` there are two `2 by 2` matrices `A,B` such that:

1. every generator and every product is exactly row rank one, hence every
   projective map is a reset with contraction coefficient zero;
2. their two right profiles have optimal one-point projective covering radius
   `delta/4` (and pairwise distance `delta/2`);
3. collapsing them to one dynamic state, while retaining letter-dependent
   scalar tolls, has optimal asymptotic
   response distortion exactly `delta/4` per letter, and no depth-uniform
   absolute error is possible.

Take

```math
A=\begin{pmatrix}0&0\\0&0\end{pmatrix}
  =\binom00\otimes(0,0),
\qquad
B=\begin{pmatrix}\delta&2\delta\\0&\delta\end{pmatrix}
  =\binom{\delta}0\otimes(0,\delta).                         \tag{ARS.16}
```

Since `d_infinity([(0,0)],[(0,delta)])=delta/2`, their midpoint
profile `(0,delta/2)` gives, and the triangle inequality forces, optimal
one-center radius `delta/4`.

Their compatibility table, in the order `A,B`, is

```math
\varphi=\delta\begin{pmatrix}0&1\\1&1\end{pmatrix}.        \tag{ARS.17}
```

Thus for a cyclic binary word of length `t`,

```math
\rho(T_w)=\delta\,[t-N_{AA}^{\rm cyc}(w)].                  \tag{ARS.18}
```

Powers of `A` force a bounded-error one-state toll to satisfy `g_A=0`,
powers of `B` force `g_B=delta`, while powers of `AB` force
`g_A+g_B=2delta`; the three requirements are inconsistent.

More quantitatively, if distortion is measured by

```math
D_1:=\inf_{g_A,g_B}\sup_{w\ne\varnothing}
 {1\over|w|}|\rho(T_w)-N_A(w)g_A-N_B(w)g_B|,                 \tag{ARS.19}
```

then

```math
D_1=\delta/4.                                                \tag{ARS.20}
```

The lower bound follows from the three cyclic words `A`, `B`, `AB`.  For the
upper bound take

```math
g_A=\delta/4,\qquad g_B=5\delta/4.                            \tag{ARS.21}
```

If both letters occur and `k` is the number of cyclic `A`-runs, then
`N_{AA}=N_A-k`, and the total defect equals

```math
\delta(k-|w|/4).
```

Since `0<=k<=|w|/2`, its magnitude is at most
`delta|w|/4`; the constant words have the same bound.

This is a width-one Ising transfer benchmark: every real `2 by 2` table is
a binary source field plus target field plus one Ising bond, and (ARS.16)
even has zero bond because it is rank one.  Therefore the obstruction is not
caused by a wide strip, incomplete supports, ties, or failure of mixing.  It
is a two-step compatibility flux omitted by the one-profile cover.

## 4. A fixed one-profile shell can contain arbitrary weighted dynamics

The preceding two-letter example is the smallest witness.  There is also a
general scale-separation statement.

### Proposition ARS.4 (small-shell universality)

Let `U_e` be any finite all-finite max-plus alphabet with every entry in
`[-1,0]`, and put `T_e=alpha U_e`, `alpha>0`.  Then

```math
\rho(T_w)=\alpha\rho(U_w)                                    \tag{ARS.22}
```

for every word.  Nevertheless every nonempty product has each one of its
rows in an interval of length at most `alpha`.  In particular,

```math
\operatorname{rad}(T_w;0)\le\alpha/2.                       \tag{ARS.23}
```

after one global scalar gauge, so the whole product semigroup has a
one-profile residual cover at radius `alpha`.

#### Proof

Positive scalar multiplication commutes with max and addition, proving
(ARS.22).  For a fixed initial state `i`, start with an optimal path from
`i` to `j`.  Changing only its last edge realizes any other terminal state
at loss at most `alpha`; reversing the comparison shows that the range of
row `i` is at most `alpha`.  Center each row separately through the free
left profile in (ARS.1).  This proves (ARS.23). `square`

Hence no theorem can bound scalar-response memory solely by delay, a fixed
projective residual-cover cardinality, and cover radius.  Below the natural
rate scale `alpha`, the shell contains a scaled copy of every bounded
all-finite weighted automaton.  This is stronger than saying that a local
error can accumulate: it says that the omitted shell can carry an arbitrary
finite response algebra.  The correct alternatives are to

* tolerate a rate error comparable to the shell radius;
* retain the directed compatibility/cycle data;
* prove an exact or summably accurate scalar cocycle;
* or impose a genuinely stronger synchronization hypothesis on the model.

## 5. Optional support leakage is a separate charge

Suppose an approximate residual presentation also declares threshold
relations and a nonempty greatest support core.  Theorem 17.1s applies to
that support system: a per-edge support shortfall `ell_e`, modulo a finite
support potential, contributes `sum ell_e` to the lower spectral bound.
This charge is independent of the terminal profile error in (ARS.7) and of
the scalar compatibility defect in (ARS.13).

There is no converse from core failure.  The mandatory Theorem-17.1u pair

```math
p=(0,-1),\qquad
T_a=\begin{pmatrix}0&-1\\-2&-3\end{pmatrix},\qquad
T_b=\begin{pmatrix}-2&-3\\1&0\end{pmatrix}                  \tag{ARS.24}
```

satisfies `pT_a=pT_b=p`, and hence every word has spectral radius zero, but
the zero-threshold relations have empty common core.  Thus an empty support
core remains a failure of one witness presentation, not semantic drift.

## 6. Benchmark score and relation to the existing theory

| system | terminal state | scalar state | conclusion |
|---|---|---|---|
| bounded-delay near-rank-one max-plus product | last `D`-suffix residual, or a right-congruent cover | compatibility data not supplied by the profile cover | rooted terminal error is depth-uniform; scalar error can have positive rate |
| width-one switching Ising (ARS.16) | one profile gives error at most `delta/2` in (ARS.7) | two endpoint residual types and the `2 by 2` compatibility table are exact | contraction zero does not repair the one-state rate `delta/4` |
| bounded all-finite weighted automaton | one profile at radius `alpha/2` after scaling | arbitrary original weighted response algebra below scale `alpha` | no universal memory bound from residual radius alone |
| Theorem-17.1u empty-core pair | one exact residual | scalar toll zero | support-core failure has no scalar implication |

The result sharpens, rather than repeats, the earlier theorems:

* **versus Theorem 17.2:** ARS.1 needs no Lipschitz contraction and incurs no
  repeated quantization, but it answers only scalar-normalized terminal
  queries.  ARS.3 proves why the same argument cannot answer accumulated
  spectral response;
* **versus Theorem 17.1l:** the cycle LP remains the right tool after a
  finite reward state is known.  ARS.4 shows that an approximate row-residual
  cover does not itself provide such a compressed reward state;
* **versus Theorem 17.1u:** exact residual profiles synchronize scalar tolls.
  Approximate profile merging can erase the directed compatibility table,
  even when every microscopic map is already an exact reset.

The clean approximate law is therefore

```math
\boxed{
\text{terminal residual error is paid once, whereas unresolved scalar
compatibility is paid once per exposed cycle step.}}
```

Any future positive depth-uniform scalar theorem must make the latter a
coboundary (or make its errors summable); projective forgetting alone cannot
do so.
