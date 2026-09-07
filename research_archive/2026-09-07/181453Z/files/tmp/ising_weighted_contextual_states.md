# Contextual states for two compositional optimization models

This note was derived before inspecting the repository.  It uses maximization;
replace `max,+` by `min,+` to obtain the minimization version.

## 1. General criterion

For a partial object `P` and a compatible future context `C`, let

\[
  \operatorname{Ans}(P,C)
\]

be the exact optimum of their composition.  Define

\[
  P\equiv P'\quad\Longleftrightarrow\quad
  \operatorname{Ans}(P,C)=\operatorname{Ans}(P',C)
  \text{ for every allowed }C.
\]

The map `P -> (C -> Ans(P,C))` is an exact state, and it is coarsest:
any exact summary assigning the same state to `P,P'` must make every future
answer equal, hence can identify only contextually equivalent objects.  The
substance in each model is to give this residual response function a finite
description and an exact composition law.

## 2. Ising chain and bounded-width fragment

### 2.1 Candidate state

Let the exposed ordered boundary have `w` spins and put
`S={-1,+1}^w`, `n=|S|=2^w`.  For a prefix fragment `P`, with internal
spins `x`, define its boundary value table

\[
 f_P(s)=\max_x E_P(x,s),\qquad s\in S,
\]

using `-infinity` for an infeasible boundary condition.  A future fragment
`C` similarly induces

\[
 g_C(s)=\max_y E_C(s,y).
\]

If all cross-cut dependence is through the exposed spins, separation gives

\[
 \operatorname{Ans}(P,C)=\max_{s\in S}\{f_P(s)+g_C(s)\}.       \tag{1}
\]

Thus `f_P` is an exact one-sided state.

### 2.2 Coarseness and minimality

Assume future contexts may apply arbitrary finite fields to boundary spins
(equivalently, may use a pinned future spin to induce such a field).  Then

\[
       P\equiv P' \quad\Longleftrightarrow\quad f_P=f_{P'}.    \tag{2}
\]

The reverse implication follows from (1).  For the forward implication,
suppose the tables differ at `s0`.  Use the future boundary field

\[
       g_M(s)=M\sum_i s_i(s_0)_i.
\]

For finite tables, choosing `M` larger than half the range of both tables
makes `s0` the unique maximizer in (1) for both prefixes.  The two answers
then differ by `f_P(s0)-f_P'(s0)`.  An infeasible entry is distinguished by
the same construction as `M` tends large: the feasible fragment can use
`s0`, whereas the other must lose at least `2M`.  Therefore every coordinate
is observable by some future; no boundary condition may be dropped merely
because it is currently suboptimal.

This is also a useful qualification: with a restricted set of futures, the
coarsest state is the quotient of tables under equality of (1) only for that
restricted set.  Equality of the raw tables is forced by selectable boundary
conditions, not by separator size alone.

If absolute optimum values matter, adding a constant to every entry changes
every future answer by that constant and is not an exact equivalence.  For a
projective/relative state one may quotient by `f ~ f+c*1`, but an exact
algorithm must carry the discarded scalar separately.

### 2.3 Exact composition

Give a two-ended fragment `A` the kernel

\[
 K_A(s,t)=\max_x E_A(s,x,t).
\]

Assign every local term to exactly one fragment.  Gluing `A` and `B` along
their common boundary gives

\[
 K_{A\circ B}(s,u)
   =\max_t\{K_A(s,t)+K_B(t,u)\}.                              \tag{3}
\]

For a prefix,

\[
       f_{P\circ A}(t)=\max_s\{f_P(s)+K_A(s,t)\}.             \tag{4}
\]

Equations (3)-(4) are associative because maximization distributes over
addition.  If arbitrary left and right boundary fields are allowed, the same
forcing proof selects any pair `(s,t)`, so a two-ended kernel is itself the
coarsest exact state for arbitrary two-sided composition.

### 2.4 Width one: a scalar relative state and an exact recurrence

For a chain, write `a=f_-=f(-1)`, `f_+=f(+1)` and

\[
       d=f_+-f_-.
\]

The exact state is `(a,d)`; the coarsest projective state is just `d`.
Appending a spin `y` with interaction `Jxy` and field `hy`
gives

\[
 f'(y)=hy+\max_x\{f(x)+Jxy\}.
\]

Consequently

\[
 a'=a-h+\max(d-J,J),                                         \tag{5a}
\]

and

\[
 d'=2h+\max(d+J,-J)-\max(d-J,J)
    =2h+\operatorname{sgn}(J)
       \operatorname{clip}(d,-2|J|,2|J|).                    \tag{5b}
\]

The clipped term is zero for `J=0`.  Formula (5b) is quantitative: a fixed
edge transmits only the old gap clipped to magnitude `2|J|`, even though no
two old gaps are equivalent against *all* possible future couplings.

If `|h|<=H` and `|J|<=K`, every post-transition gap lies in
`[-2(H+K),2(H+K)]`.  Hence an `epsilon`-net for the relative chain state has
`O((H+K)/epsilon)` cells.  Max-plus updates are nonexpansive in sup norm, so
an error introduced once does not grow under exact later updates; rounding
after every one of `L` updates can in the worst case accumulate `L*epsilon`.

There is also a sharp two-ended chain statement.  Every real `2 by 2`
kernel is already an Ising edge:

\[
       K(x,y)=c+u x+v y+Jxy.                                 \tag{5c}
\]

Indeed, `c,u,v,J` are the four Walsh coefficients of `K`.  Thus arbitrary
endpoint kernels are not merely an ambient relaxation for a chain: they are
realizable with one pairwise coupling and two endpoint fields.  Under
arbitrary endpoint contexts, a reusable chain segment therefore has four
exact scalar coordinates, or three normalized/projective coordinates.  For
normalized range `B`, its projective covering number is
`Theta((B/epsilon)^3)`.

### 2.5 Quantitative approximate state complexity at width `w`

For tables define the exact contextual metric

\[
 d_{ctx}(f,f')=\sup_g
 \left|\max_s(f(s)+g(s))-\max_s(f'(s)+g(s))\right|.
\]

If boundary conditions are selectable, the elementary max inequality gives
the upper bound and a selecting field gives the reverse bound, so

\[
       d_{ctx}(f,f')=\|f-f'\|_\infty.                         \tag{6}
\]

After quotienting a common offset, the induced metric is

\[
 \inf_c\|f-f'-c\mathbf 1\|_\infty
   ={1\over2}\left(\max_s\Delta(s)-\min_s\Delta(s)\right),
 \quad \Delta=f-f'.                                         \tag{7}

Suppose normalized tables have oscillation at most `B`.  Anchoring one of
the `n=2^w` coordinates and quantizing the other `n-1` gives an upper cover
of size

\[
       (1+O(B/\epsilon))^{2^w-1}.                             \tag{8}

\]

The exponent is necessary for general width-`w` Ising fragments, rather
than merely an artifact of tabulation.  To realize an arbitrary table
`F:S->R`, subtract `c=min_s F(s)` and put `lambda_a=F(a)-c>=0`.
For every assignment `a`, introduce an independent binary internal spin
`y_a in {0,1}` and the pairwise Ising/QUBO energy

\[
 y_a\lambda_a\left(\sum_i a_i s_i-(w-1)\right).
\]

Maximizing over `y_a` contributes `lambda_a` exactly when `s=a`, and zero
otherwise.  Summing the independent gadgets and restoring `c` realizes
`F(s)`.  Replacing `y_a` by an Ising spin is an affine change of variables,
so all terms remain unary or pairwise.  The bags `S union {y_a}`, processed
one at a time, have size `w+1`; the construction has treewidth `w`.

Restrict `F(s_ref)=0` and every other coordinate to a grid with spacing
`4*epsilon` in an interval of length `B`.  By (7), distinct tables are more
than `epsilon` apart.  Therefore every projective `epsilon`-summary needs at
least

\[
 \left(1+\left\lfloor {B\over4\epsilon}\right\rfloor\right)^{2^w-1}
                                                                    \tag{9}
\]

states for this gadget-complete family.  Combining (8)-(9), the precision
cost is

\[
       \Theta\bigl(2^w\log(B/\epsilon)\bigr)\text{ bits},
\]

while the number of distinguishable normalized states is doubly exponential
in `w`.  For a strict nearest-neighbor lattice strip, (8) remains an upper
bound and every differing reachable coordinate remains distinguishable; the
gadget lower bound applies when the width model permits general pairwise
interactions within a width-`w` tree decomposition.

With coefficient bounds, `B` need not grow with prefix length.  Changing a
boundary condition affects only terms incident to the current boundary, so

\[
 \operatorname{osc}(f_P)
 \le 2\sum_{\text{terms incident to boundary}} |\text{coefficient}|,
\]

which is `O(w(H+Delta*K))` at bounded local degree `Delta`.

## 3. Weighted-language optimization

### 3.1 Candidate residual state

Consider a max-plus weighted automaton with `m` control states, initial row
`alpha`, final column `beta`, and a transition matrix `T_a` for each letter.
For a word `x`, set

\[
 T_x=T_{x_1}\otimes\cdots\otimes T_{x_k},\qquad
 L(x)=\alpha\otimes T_x\otimes\beta,
\]

where `(A tensor B)_{ij}=max_k(A_{ik}+B_{kj})`.

The forward vector of a prefix and backward vector of a suffix are

\[
 p_x=\alpha\otimes T_x,\qquad h_y=T_y\otimes\beta,
\]

and

\[
       L(xy)=\max_q\{p_x(q)+h_y(q)\}.                         \tag{10}

\]

The coarsest exact prefix state is the weighted residual

\[
       \rho_x:y\longmapsto L(xy).                             \tag{11}

\]

Equivalently, it is the quotient of forward vectors under

\[
 p\sim p'\quad\Longleftrightarrow\quad
 p\otimes h_y=p'\otimes h_y\quad\text{for every suffix }y.   \tag{12}

\]

The raw `m`-vector can be strictly finer for a fixed automaton: a coordinate
that no accepting suffix can expose, or two coordinates with identical
future behavior, need not survive (12).  If future contexts may choose an
arbitrary terminal vector, coordinate-selecting vectors prove `p~p'` iff
`p=p'`.

### 3.2 Minimality and exact composition

If two prefixes have different residuals, by definition some suffix `y`
gives different final scores; hence no exact online summary can identify
them.  Conversely (11) directly answers every suffix query.  This proves
coarseness without assuming a particular automaton representation.

Appending a block `u` gives

\[
       p_{xu}=p_x\otimes T_u,\qquad T_{uv}=T_u\otimes T_v.     \tag{13}

\]

Directly on residuals, define `(D_u rho)(y)=rho(uy)`.  Then
`rho_{xu}=D_u rho_x` and `D_v(D_u rho)=D_{uv}rho`.  Thus both the matrix and
the representation-independent residual have exact associative composition.

Absolute and projective states must again be distinguished.  Replacing `p`
by `p+c*1` adds `c` to every suffix answer, so exact weighted-language values
distinguish it.  If only relative continuation behavior matters, `[p]` is a
projective state and `c` may be accumulated separately.  Even a one-control-
state automaton with a loop of weight one has infinitely many exact residuals
`rho_{a^k}(y)=k+|y|`, although all normalized residuals coincide.

### 3.3 Quantitative approximate state complexity

Max-plus multiplication is nonexpansive:

\[
       \|p\otimes T-p'\otimes T\|_\infty
       \le \|p-p'\|_\infty.                                  \tag{14}

\]

Therefore an `epsilon` approximation of a forward vector changes every
later suffix answer by at most `epsilon` if later operations are exact.  With
arbitrary terminal vectors, selector vectors make (14) tight and the exact
contextual metric is sup norm; modulo offsets it is again (7), now over `m`
coordinates.

If finite normalized coordinates lie in an interval of length `B`, the
projective covering number is

\[
       (1+O(B/\epsilon))^{m-1}.                               \tag{15}

\]

Arbitrary source-to-state weights realize every vector in such a box, and
arbitrary terminal weights expose every coordinate, so the same grid packing
as in (9) gives

\[
       \left(1+\left\lfloor {B\over4\epsilon}\right\rfloor
       \right)^{m-1}                                         \tag{16}

\]

as a lower bound.  Thus `Theta(m log(B/epsilon))` bits are necessary and
sufficient in this fully observable model.  For a fixed automaton, replace
sup norm by

\[
 d_A(p,p')=\sup_y|p\otimes h_y-p'\otimes h_y|;
\]

its zero classes are precisely (12), and the approximate state count is the
covering number of reachable forward vectors in this residual metric.  For
a horizon `H` and transition weights bounded by `W`, one may take
`B=O(HW)` (plus the initial spread), giving the explicit upper bound
`2^m(1+O(HW/epsilon))^(m-1)` when reachability masks are counted.

## 4. Post-derivation comparison

No repository material was used in Sections 1-3.  This section records the
comparison made afterward.

### 4.1 Classical transfer matrices and Myhill--Nerode

Equations (3)-(4) are precisely the zero-temperature transfer-matrix or
max-sum dynamic-programming construction.  At positive temperature, summing
Boltzmann weights gives ordinary nonnegative transfer-matrix multiplication;
in log coordinates the `max` in (3) is replaced by `log-sum-exp`.  The
ground-state limit is max-plus multiplication.  The width-one gap `d` is the
projective max-product/min-sum message (a zero-temperature cavity field), and
(5b) is its Ising-specific saturation law.  The classical machinery proves
sufficiency and computation; the contextual argument identifies the precise
query class under which the table is also minimal.

For the second model, the bi-infinite matrix

\[
       H(x,y)=L(xy)
\]

is the weighted Hankel matrix, and `rho_x` is exactly its row indexed by
`x`.  Equality of these rows is the right Nerode relation.  In the Boolean
case, the number of distinct residual rows is the state count of the minimal
DFA (Nerode, *Linear Automaton Transformations*, 1958,
DOI `10.1090/S0002-9939-1958-0135681-9`).

Two weighted qualifications matter.

1. A finite `m`-dimensional weighted linear representation can have
   infinitely many distinct exact residual rows.  The one-state unit-weight
   loop above is the smallest example.  Thus finite residual *cardinality*
   is not the right characterization of nondeterministic/linear weighted
   automaton dimension.  Over fields, that role is played by Hankel rank;
   over max-plus it becomes a semimodule/factorization problem rather than
   ordinary rank.
2. Deterministic weighted-automaton minimization normally redistributes
   path weights (weight pushing) before merging states with the same
   continuation behavior.  This is the algorithmic counterpart of splitting
   an absolute offset from the projective residual.  Mohri's
   [weighted-automata survey](https://cs.nyu.edu/~mohri/pub/hwa.pdf) gives
   the tropical determinization, weight-pushing, and minimization versions.

Therefore (11)-(12) are accurately called a weighted Nerode residual, but it
would be inaccurate to claim the unqualified Boolean theorem “finite
weighted automaton iff finitely many exact residuals.”

### 4.2 Repository collisions

The main generic conclusions already occur in the repository.

- `extremal_information/drafts/phase2_feature_growth.md`, Theorems FG.1--FG.4,
  already proves that the boundary response kernel is the coarsest exact
  endpoint statistic, that endpoint-response distance is entrywise sup
  distance, that gluing is max-plus multiplication, and that universal
  one- and two-sided separators require respectively `q^w` and `q^(2w)`
  scalar entries (with matching cube metric entropy).  Its Sections 5.1--5.2
  already name finite-state spin chains and trellis/weighted-automaton
  decoding as applications.
- `extremal_information/drafts/phase3_contextual_response_law.md`, Proposition
  CRL.1, already gives the abstract future-response pseudometric, its
  translation contraction, the zero-distance syntactic congruence, and its
  coarsest-state property.  It also states the same repeated-rounding error
  warning used after (14).
- `extremal_information/drafts/phase2_known_model_validation.md` adds the
  classical max-plus cyclicity theorem for homogeneous repetition.  That is
  a stronger asymptotic saturation statement for a fixed repeated kernel and
  was not independently rederived here.

So the boundary table, max-plus composition, sup-norm isometry, generic cube
packing, and abstract contextual quotient are direct collisions, not new
results for this repository.

### 4.3 Model-specific additions not found in the collision search

The parts that appear to add something are narrower.

1. Equations (5a)-(5b) reduce the exact one-sided Ising-chain state to an
   offset and one gap, and give the exact signed clipping recurrence and the
   `O((H+K)/epsilon)` normalized state count.
2. Equation (5c) shows that the repository's arbitrary `2 by 2` universal
   kernel lower bound is already realizable by one genuine Ising edge, not
   merely by an unrestricted landscape.
3. The lookup gadget following (8) realizes every one-sided boundary table
   using only unary and pairwise Ising terms at treewidth `w`.  This closes
   the local-realizability caveat for general treewidth-`w` pairwise Ising
   fragments and yields the lower bound (9).  It does **not** by itself close
   that caveat for a strict bounded-degree nearest-neighbor rectangular
   lattice.
4. Equations (11)-(12) distinguish the raw trellis vector/kernel, which is
   minimal under arbitrary terminal metrics, from the potentially smaller
   residual quotient for the suffix set of one fixed weighted automaton.
   This resolves an otherwise easy context-class ambiguity.
5. The normalized/projective estimates (8)-(9) and (15)-(16) sharpen the
   ambient entry count to the exact exponent after a common additive offset
   is separated.

## 5. Finite audit

`tmp/verify_ising_weighted_contextual_states.py` checks the derivations with
exact integer arithmetic.  Its current run passes:

- `29,575` full `(offset,gap)` chain updates, including both signs of `J`;
- `168` brute-force boundary evaluations of arbitrary-table pairwise lookup
  gadgets through width three;
- `750` random selector-field instances of the sup-norm response isometry;
- `500` each of max-plus associativity, nonexpansiveness, and prefix/suffix
  factorization; and
- all `127` suffixes of length at most six in a fixed two-state automaton
  where two distinct forward vectors have identical residuals.
