# Mesoscopic shell pairs admit low-cap exact-sign selectors

**Status.** Task-local theorem draft.  The archive comparison in Section 6
is part of the statement's scope: the first result is only a labelled
evaluation theorem, while the second gives a genuine all-spins-free parent
selector but not a contextual packing of distinct child signings.

The input is deliberately much smaller than a Boolean landscape.  It is one
pair of oriented cut words and, for the physical result, their common child
signing.  No claim below transfers the resulting selector to a cross-order
recurrence.

## 1. Notation

Put

```math
E={n\choose2},\qquad
c(x)=(x_ix_j)_{i<j},\qquad
H_A(x)=\langle A,c(x)\rangle .                 \tag{MQ.1}
```

An augmented cut is `z=sigma c(u)`, with `sigma in {+-1}`.  Its positive
energy is `\langle A,z\rangle`.  Write

```math
d_P(z,z')=\min\{d_H(z,z'),E-d_H(z,z')\}.        \tag{MQ.2}
```

For an exact sign matrix `B in {+-1}^{n times s}` define its free-shore
response

```math
F_B(x)=\max_{y\in\{+-1\}^s}|x^TBy|=\|B^Tx\|_1. \tag{MQ.3}
```

## 2. Every mesoscopic pair has a low-cap labelled sign query

### Lemma MQ.1 (pair-conditioned exact-sign query)

Let `z,z' in {+-1}^E`, put `h=d_H(z,z')`, and let `0<=t<=h`.  There is an
exact edge signing `b in {+-1}^E` such that

```math
Q(b)\le t+\sqrt{2E(n+3)\log2},                  \tag{MQ.4}
```

and

```math
\langle b,z\rangle-\langle b,z'\rangle
\ge2t-2\sqrt{2h\log4}.                          \tag{MQ.5}
```

Consequently, if `M=Theta(n^(3/2))` and `h>=delta M`, taking
`t=delta M/2` gives

```math
Q(b)=O(n^{3/2}),\qquad
\langle b,z\rangle-\langle b,z'\rangle
\ge delta M-O(n).                               \tag{MQ.6}
```

For a finite family whose oriented words have pairwise Hamming distance at
least `delta M`, applying the lemma separately to every unordered pair gives
one common bank of exact-sign queries with the uniform bounds (MQ.6).

#### Proof

Let `D={e:z_e ne z'_e}` and `mu=t/h` (the case `h=0` is vacuous).  Choose
the coordinates of `b` independently.  On `D`, impose

```math
\mathbb E b_e=\mu z_e,
```

and off `D` make `b_e` unbiased.  For every augmented cut `w`,

```math
|\mathbb E\langle b,w\rangle|
=\mu\left|\sum_{e\in D}z_ew_e\right|\le t.     \tag{MQ.7}
```

Hoeffding and a union bound over at most `2^n` augmented cuts show that,
with failure probability at most `1/4`, all their centred sums have modulus
at most the square root in (MQ.4).  On the other hand,

```math
{1\over2}\{\langle b,z\rangle-\langle b,z'\rangle\}
=\sum_{e\in D}b_ez_e
```

has mean `t`.  Its lower deviation exceeds `sqrt(2h log 4)` with
probability at most `1/4`.  The two good events therefore intersect, proving
(MQ.4)--(MQ.5). `square`

The query in MQ.1 is genuinely exact-sign and has the right cap scale.  It
is nevertheless a **labelled/pinned** query: its two outputs are evaluations
at the prescribed words.  Adding `b` coefficientwise to `A` leaves the
exact-sign class, and maximizing the sum forgets which one of the two words
was declared as the state.  The next theorem identifies what can still be
realized physically.

This first failure has an exact same-order form.  Suppose both relation
classes `{e:A_e=b_e}` and `{e:A_e=-b_e}` are nonempty.  If an exact signing
`C`, a scalar `lambda`, and a constant `r` satisfied

```math
H_C(x)=H_A(x)+\lambda H_b(x)+r\quad\hbox{for every Boolean }x, \tag{MQ.6a}
```

uniqueness of the Boolean Fourier expansion would give `r=0` and
`C_e=A_e+lambda b_e`.  Exact signs would then require both
`|1+lambda|=1` and `|1-lambda|=1`, hence `lambda=0`.  If only one relation
class occurs, the alleged query is merely a scalar copy of `A`.  Thus a
nondegenerate MQ.1 query has no same-order exact-sign additive overlay;
auxiliary spins or a different composition law are genuinely necessary.

## 3. A diffuse public shore turns a mesoscopic pair into a physical selector

The elementary inequality below is the useful source of diffuseness.  If
`c_j>=0`, `L=sum_j c_j`, and `V=sum_j c_j^2`, then

```math
\sum_j\min\{k,c_j\}
\ge {kL^2\over kL+V}.                            \tag{MQ.8}
```

Indeed `min{k,c}>=kc/(k+c)`, and Cauchy--Schwarz gives

```math
\sum_j{c_j\over k+c_j}
\ge{L^2\over kL+V}.
```

### Theorem MQ.2 (linear-shore exact-sign pair selector)

Fix constants `0<kappa<=K<infinity`.  For all sufficiently large `n`, let
`u,v in {+-1}^n` have projective vertex distance

```math
k=\min\{d_H(u,v),n-d_H(u,v)\},\qquad
kappa\sqrt n\le k\le K\sqrt n.                 \tag{MQ.9}
```

There are an exact sign cross block `B in {+-1}^{n times n}` and a Boolean
shore state `y_0` such that

```math
F_B(u)=\|B\|_{infinity to1}
\le (K+C_0)n^{3/2},                              \tag{MQ.10}
```

while

```math
F_B(u)-F_B(v)\ge c_(kappa,K)n^{3/2},             \tag{MQ.11}
```

where `C_0` and `c_(kappa,K)>0` are absolute apart from the displayed
parameters.  Moreover, one may prescribe any exact-sign shore signing `D`
and any positive ground state `y_0` of `D`; the columns of `B` can be gauged
so that the same `y_0` attains `u^TBy_0=F_B(u)`.

One explicit constant accounting is the following.  If a public
Rademacher matrix in the proof is chosen with

```math
\|C\|_(2 to2)\le L_0(\sqrt m+\sqrt n),           \tag{MQ.12}
```

then one may take

```math
C_0=2L_0,
\qquad
c_(kappa,K)={kappa\over4KL_0+8L_0^2}.            \tag{MQ.13}
```

#### Proof

Replace `v` by `-v` if necessary and let

```math
S=\{i:u_i ne v_i\},\qquad |S|=k,qquad m=n-k.
```

Choose a public exact sign matrix `C in {+-1}^{m times n}` satisfying
(MQ.12).  Such a matrix exists by the standard rectangular Rademacher
operator-norm bound.  Choose `x_0` maximizing `\|C^Tx\|_1` and put

```math
c_j=|(C^Tx_0)_j|,
\qquad L=\sum_jc_j,
\qquad V=\sum_jc_j^2.                            \tag{MQ.14}
```

Row-gauge `C` so that the restriction of `u` to `S^c` maps to `x_0`, and
column-gauge it so that

```math
C^Tu_(S^c)=(c_jy_(0,j))_j.                       \tag{MQ.15}
```

On every row `i in S`, set

```math
B_(ij)=u_i y_(0,j),                              \tag{MQ.16}
```

and use the gauged `C` on `S^c`.

For arbitrary `x`, put `a=sum_(i in S)u_ix_i`.  The triangle inequality and
the definition of `x_0` give

```math
\|B^Tx\|_1
\le n|a|+\|C^Tx_(S^c)\|_1
\le nk+L.                                       \tag{MQ.17}
```

At `x=u`, (MQ.15)--(MQ.16) attain equality and the maximizing shore spin is
`y_0`.  At `x=v`,

```math
F_B(v)=\sum_j|c_j-k|,
```

so

```math
F_B(u)-F_B(v)=2\sum_j\min\{k,c_j\}
\ge {2kL^2\over kL+V}.                          \tag{MQ.18}
```

It remains only to check scale.  Averaging `\|C^Tx\|_1` over uniform
Boolean `x` and using the `p=1` Khintchine inequality gives

```math
L\ge n\sqrt{m/2}\ge {1\over2}n^{3/2}.           \tag{MQ.19}
```

Also

```math
L\le\sqrt n\|C\|\sqrt m\le2L_0n^{3/2},
\qquad
V\le\|C\|^2m\le4L_0^2n^2.                    \tag{MQ.20}
```

Equations (MQ.9), (MQ.18)--(MQ.20) give (MQ.11)--(MQ.13), while
(MQ.17) gives (MQ.10). `square`

The ground-pair gauge in this proof is public auxiliary information; it
does not optimize the child `A` or use `M_n` at the target order.  It does,
however, make MQ.2 existential rather than an efficient compiler theorem.
On structured orders a spectrally flat public matrix with a known ground
pair can replace this step.  A uniformly explicit arbitrary-order version
is open.

## 4. Exact completion and optimizer exclusion

### Corollary MQ.3 (same-spin parent selection)

Let `A` be a hollow exact signing with `Q(A)=O(n^(3/2))`, and let
`z=sigma c(u)` obey

```math
\langle A,z\rangle\ge Q(A)-d.                   \tag{MQ.21}
```

Let `v` satisfy (MQ.9), construct `B,y_0` by MQ.2, and choose a hollow exact
shore signing `D` of order `n` with

```math
Q(D)=O(n^{3/2}),\qquad H_D(y_0)=Q(D).
```

(A random signing supplies such a `D`; gauge its ground to the prescribed
`y_0`.)  Form the complete exact-sign parent

```math
P=
\begin{pmatrix}
A&\sigma B\\
\sigma B^T&\sigma D
\end{pmatrix}.                                  \tag{MQ.22}
```

Then `P` has order `2n`, cap `O(n^(3/2))`, and the prescribed old spin `u`
is within `d` of the full parent cap:

```math
Q(P)-\max_y|H_A(u)+\sigma u^TBy+\sigma H_D(y)|
\le d.                                           \tag{MQ.23}
```

If `d=0`, `(u,y_0)` (with the sign dictated by `sigma`) is an exact global
parent maximizer and

```math
Q(P)=Q(A)+F_B(u)+Q(D).                           \tag{MQ.24}
```

At the same time the rival old spin is physically excluded:

```math
\max_y|H_A(u)+\sigma u^TBy+\sigma H_D(y)|
-\max_y|H_A(v)+\sigma v^TBy+\sigma H_D(y)|
\ge c_(kappa,K)n^{3/2}-d.                       \tag{MQ.25}
```

#### Proof

Every parent energy has absolute value at most

```math
Q(A)+F_B(u)+Q(D),                                \tag{MQ.26}
```

because `u` maximizes the entire bipartite block.  At `(u,y_0)` all three
components have sign `sigma`, and (MQ.21) therefore attains the right side
of (MQ.26) up to `d`.  This proves (MQ.23)--(MQ.24).  For fixed old spin
`v`, its rooted parent response is at most

```math
Q(A)+F_B(v)+Q(D).
```

Subtract from the target value and use MQ.2. `square`

Now suppose `M_n` is between two positive constant multiples of `n^(3/2)`
and two augmented cut words obey

```math
aM_n\le d_P(z,z')\le bM_n.                       \tag{MQ.27}
```

For large `n`, their projective vertex distance satisfies (MQ.9), with
constants depending only on `a,b` and the two cap constants.  Indeed, in
this `o(E)` regime the exact identity is

```math
d_P(z,z')=k(n-k).                                \tag{MQ.28}
```

Thus every `Theta(M_n)` pair from the new fractional-reservoir packing has
an exact-sign bounded-cap parent selector.  Applying MQ.3 separately to all
ordered pairs gives a common externally selected bank of pairwise exclusion
contexts.

## 5. Why this still does not prove a contextual response packing

The distinction is sharp.

1. **MQ.1 is labelled.**  It separates the prescribed evaluations
   `z,z'`, not two unconstrained optima.
2. **MQ.3 is physical but selector-valued.**  It makes `u` a genuine (or
   near) parent optimizer and puts `v` behind it.  Both roots belong to the
   same child matrix `A`.  The scalar answer `Q(P)` need not reveal which
   rival was excluded, and all pair queries can have the same scalar cap.
   Equivalently, MQ.3 does give a `K`-point packing of the **rooted-slice**
   systems `(A,z^u)` under the pair-context bank, but freeing the old spin
   erases the system label and leaves one unrooted parent.
3. **A contextual packing needs systems, not just roots.**  To infer
   `log K` response bits one still needs exact-sign child states `A^u` for
   which one common query bank produces pairwise separated scalar parent
   responses.  Pairwise shell geometry supplies no such map
   `u mapsto A^u`.
4. **The loss is not polarization.**  MQ.3 keeps the cross block and shore
   joint and aligns their exact ground pair.  The remaining gap is rooted
   versus unrooted information: after maximization a parent reports its
   winner's value, not the identities of all suppressed competitors.

There is a formal non-identifiability behind item 3.  Any unrooted parent
constructor `P(A,C)` whose state input is only `A` assigns exactly the same
response vector to `(A,z)` and `(A,z')`; the root label is absent from the
physical system.  Encoding the root into a new child matrix is therefore a
logically additional operation.  Universal locking is ruled out at low cap
by the archived universal-pin theorem, while the archived sparse-flip
compiler pays the separation fraction and gives only `O(n)` signal at
edge distance `Theta(M_n)`.  MQ.3 evades both facts as a selector, but not as
a state encoder.

## 6. Archive comparison and frontier impact

The closest archived results are:

- Theorem 36.3, which gives an exponential **pinned-field response packing**
  but pays a quadratic physical calibration.
- Theorem 36.11/AO.2, whose `Theta(sqrt n)` rank-one shore gives physical
  scalar-response separation `gamma n^(3/2)`.  At
  `gamma=Theta(M_n/E)` this is only `Theta(n)`.
- the universal-pin cap barrier, which applies to a context that locks one
  old state for every possible child;
- the rowwise microcanonical compiler, which realizes a prescribed field
  but requires an additional uniform affine endpoint law to preserve a
  scalar response gap.

MQ.1 is not equivalent to any of these: it records the simple but previously
unstated fact that exact-sign **quadratic** labelled queries already see the
energy-scale packing.  MQ.2--MQ.3 use a linear-width diffuse shore rather
than a repeated `sqrt(n)`-width rank-one shore.  They physically select a
mesoscopic shell word at the correct scale without universal locking and
without reconstructing the child optimum.

The frontier impact is therefore precise but limited:

```text
PROVED: Theta(M_n) projective shell geometry is visible to exact-sign,
        O(n^(3/2))-cap labelled queries and to pairwise all-spins-free
        optimizer selectors.

NOT PROVED: a packing of scalar parent responses, a reusable child-state
            encoding, a cross-order congruence, or any recurrence for M_n.
```

The next non-equivalent lemma would have to turn the optimizer-exclusion
table into scalar response separation while using strictly less than a
root-dependent near-minimizer child for every word.  Merely building more
pair queries cannot do this.

## 7. Exact finite verification

[`../experiments/verify_mesoscopic_pair_query_selector.py`](../experiments/verify_mesoscopic_pair_query_selector.py)
exhaustively checks (MQ.8) and the finite identities (MQ.14)--(MQ.25) on
twelve seeded instances of orders `4<=n<=7`, including the complete
`2n`-spin parent.  The frozen output is
[`../experiments/mesoscopic_pair_query_selector_results.json`](../experiments/mesoscopic_pair_query_selector_results.json).
These computations check algebra and signs only; the asymptotic existence
of the spectrally flat public block is proved probabilistically above.
