# Two-fiber conference systems: cyclic certificates and the centered-module barrier

Date: 2026-07-31. This note classifies every statement as an exact theorem,
an exact finite certificate, or an open construction problem.

## 1. Boolean invariant-plane formulation

Put

```math
N=4k^2+2,qquad s=2k^2+1,qquad Q=4k^2+1=N-1.
```

Let `S` be a symmetric conference signing of order `N`, so
`S^2=Q I`. The target equal two-fiber quotient exists after switching if and
only if there are Boolean vectors `p,q` such that

```math
p^{\mathsf T}q=0,qquad
Sp=p+2kq,qquad Sq=2kp-q.                         \tag{1}
```

In fact it is enough to find a Boolean `p` for which

```math
q={Sp-p\over 2k}                                   \tag{2}
```

is Boolean. Comparing `||Sp||^2=QN` with the squared norm of
`p+2kq` forces `p^Tq=0`, and applying `S` to (2) gives the second equation
in (1).

To see the quotient exactly, switch by `diag(p)` and put `r=p\circ q`.
Then

```math
(\operatorname{diag}p)S(\operatorname{diag}p)\mathbf1
  =\mathbf1+2kr,
\qquad
(\operatorname{diag}p)S(\operatorname{diag}p)r
  =2k\mathbf1-r.                                    \tag{3}
```

The two signs of `r` occur `s` times each. Solving the sum and difference
of (3) on the two fibers gives internal row sums `2k,-2k` and both cross
row sums equal to one. Thus (1) is exactly the desired equitable partition,
not merely a necessary spectral shadow of it.

There is no immediate arithmetic obstruction. The quotient

```math
T_k=\begin{pmatrix}2k&1\\1&-2k\end{pmatrix}
```

has characteristic polynomial `t^2-Q`, while the full conference matrix
has characteristic polynomial `(t^2-Q)^s`. The orthogonal complement simply
has characteristic polynomial `(t^2-Q)^(s-1)`. Row-sum parities and the
conference identity modulo two are also compatible. The usual necessary
sum-of-two-squares condition is vacuous here because `Q=(2k)^2+1`.

## 2. A cyclic subfamily and one exact identity

A much more explicit sufficient construction uses two length-`s` cyclic
sequences. Let `a_0=0`, let every other `a_j` be a sign, and let every `c_j`
be a sign. Assume

```math
a_{-j}=a_j,qquad \sum_ja_j=2k,qquad \sum_jc_j=1.   \tag{4}
```

Let `A=circ(a)` and `C=circ(c)`. Because circulant matrices commute, the
single complementary-autocorrelation identity

```math
AA^{\mathsf T}+CC^{\mathsf T}=Q I_s                \tag{5}
```

implies that

```math
S=\begin{pmatrix}A&C\\C^{\mathsf T}&-A\end{pmatrix} \tag{6}
```

is a symmetric conference signing with exactly the target quotient. This
state uses two cyclic sign sequences and (5); verification takes polynomial
time and does not encode the full Boolean bridge response.

There is an exact almost-difference-set form. Define

```math
P=\{j:a_j=-1\},\qquad R=\{j:c_j=-1\},
```

and let `N_X(h)=|X cap (X+h)|`. Conditions (4)--(5) are equivalent to

```math
|P|=k^2-k,quad |R|=k^2,quad P=-P,
```

and

```math
\boxed{N_P(h)+N_R(h)+\mathbf1_P(h)=k^2-k
       \quad(h\ne0).}                               \tag{7}
```

Thus the low-multiplicity shifts of the almost supplementary difference
pair must be the first block `P` itself. Call (7) the **self-indexed ASDS
condition**.

Armario and Flannery, *Almost supplementary difference sets and quaternary
sequences with optimal autocorrelation* (2020), prove the general incidence-
matrix correspondence for ASDS (their Theorem 2) and an equivalence between
certain amicable ASDS and optimal quaternary sequences (their Theorem 3):
[arXiv:1911.08828](https://arxiv.org/abs/1911.08828). Their existence remarks
give broad ASDS families when `2s-1` is a prime power. The cited results do
**not** impose that the low-multiplicity shift set equals `P`. Therefore they
do not, as stated, prove (7) or an infinite family of (6). Treating a generic
ASDS/OQS as a solution would drop the decisive self-indexing condition.

### 2.1 Prime Paley conferences give a deterministic two-cycle construction

There is nevertheless a clean algebraic source of target two-fiber systems.
Let `Q=4k^2+1` be prime and view the Paley conference matrix on the projective
line `P^1(F_Q)`. A nonsplit torus in `PSL(2,Q)` has order

```math
s={Q+1\over2}=2k^2+1
```

and acts without fixed points on the projective line. It therefore has two
orbits of length `s`. A Mobius generator preserves the Paley two-graph; on a
conference representative this is a signed permutation symmetry. Since `s`
is odd, its sign cocycle on each cyclic orbit can be gauged away. Ordering the
two cycles consecutively makes all four blocks circulant.

The quotient entries are integers and its square is `Q I_2`. Because a prime
`Q=4k^2+1` has, up to signs and order, the unique sum-of-two-squares
representation `Q=(2k)^2+1`, the two-fiber row sums are the target
`2k,1,-2k` after harmless fiber switches. Thus the general two-fiber block
system exists whenever `4k^2+1` is prime.

The verifier makes the construction completely explicit using

```math
g_t(x)=-{1\over x+t}
```

on the projective line, selecting the first `t` for which `g_t` has order
`s`, factoring its exact Paley switching cocycle, and gauging it around the
two cycles. For `Q=5,17,37,101` it further verifies that the chosen gauge has
the stronger form (6), including `B=-A`. This gives the new exact positive
control `k=5`, `s=51`, `N=102` (with `t=6`).

This is an unbounded-looking formula but not a proved infinite sequence of
orders: infinitude of primes of the form `4k^2+1` is an instance of the open
quadratic-prime problem. More importantly for the present project, these are
ordinary Paley conference matrices in a new gauge. The already proved prime-
Paley cosine saturation theorem gives their normalized Boolean cap a
subsequence approaching `1/2`. Hence the construction proves existence of
the quotient algebra, not useful landing near a hypothetically smaller true
optimum.

This exact interpretation is consistent with the adjacent two-circulant
literature: Barrera Acevedo--O Cathain--Dietrich prove that prime Paley type-II
Hadamard matrices admit two-circulant-core presentations
([Proposition 4.2](https://doi.org/10.1007/s10801-021-01033-x)), and
Iverson--Jasper--Mixon discuss two-circulant Paley/ETF representations
([arXiv:2410.17379](https://arxiv.org/abs/2410.17379)). The projective-line
calculation above checks the mapping at the conference normalization actually
used here rather than importing the Hadamard or ETF statement by analogy.

## 3. Exact finite certificates

The program
[`computations/two_fiber_cyclic_conference.py`](../computations/two_fiber_cyclic_conference.py)
independently verifies all sequence, ASDS, block, and full conference
identities with integer arithmetic. Its result file is
[`computations/results/two_fiber_cyclic_conference.json`](../computations/results/two_fiber_cyclic_conference.json).

Besides the order-six base, the saved non-Paley search branch certifies:

```text
k=2, s=9, N=18
a=(0,1,1,1,-1,-1,1,1,1)
c=(1,1,1,-1,1,-1,-1,1,-1)

k=3, s=19, N=38
a=(0,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1)
c=(1,1,1,1,-1,1,-1,1,-1,-1,-1,1,1,1,-1,-1,1,-1,-1).
```

These are exact finite existence results, not evidence by floating-point
residuals. The `k=2` solution has `P={4,5}` and
`R={3,5,6,8}` in `Z_9`; the `k=3` solution has the sets recorded in the JSON
certificate. An exhaustive signature match over all admissible symmetric
`a` and all row-sum-one `c` found the displayed solutions. This establishes
that the cyclic restriction is nonempty through order 38.

The verifier now classifies their conference graphs exactly:

- `k=2` is switching/permutation equivalent to the Paley conference of order
  18. The JSON records and checks an explicit 17-vertex conference-graph
  isomorphism.
- The saved `k=3` certificate is **not** Paley. Its normalized conference graph
  has 615 four-cliques and 65 five-cliques, while the Paley graph on 37
  vertices has 555 four-cliques and no five-cliques. Clique counts are graph
  invariants, so this is a short exact nonisomorphism certificate.
- Separately, the nonsplit-torus construction gives a Paley cyclic certificate
  at `k=3`. Thus the self-indexed ASDS solution space already contains both
  Paley and non-Paley equivalence classes at order 38.

An exhaustive enumeration sharpens the last statement. There are exactly 627
oriented `k=3` pairs `(a,c)`. Under the natural operations of simultaneously
multiplying both index sets by a unit of `Z_19`, cyclically shifting `c`, and
reversing `c`, they form exactly three orbits, of sizes 171, 342, and 114.
The corresponding normalized conference graphs have respective
`(K_4,K_5)` counts

```text
(555,0), (615,65), (610,48).
```

Consequently all three natural orbits are graph-inequivalent: the first is
the explicit Paley nonsplit-torus orbit, and the other two are distinct
non-Paley classes. This is an exact finite classification under the stated
natural equivalence only; it does not supply an asymptotic family or a cap
bound.

### 3.1 The strong cyclic subclass is empty at `k=4`

The next discriminator now has a certified answer. For `k=4`, `s=33`, there
are `binom(16,6)=8008` possible symmetric internal negative supports `P`.
For 7,998 of them the audit constructs an explicit integer vector `z` with

```math
z^{\mathsf T}(65I-AA^{\mathsf T})z<0.
```

This is an exact obstruction to a bridge `C`, since
`CC^{\mathsf T}=65I-AA^{\mathsf T}` would be positive semidefinite. The ten
surviving supports form one orbit under multiplication by units of `Z_33`.
For a representative with negative half-support

```text
{1,3,4,5,12,15},
```

the required bridge negative support `R` has a fixed 16-component circular
distance profile. Summing the profile at multiples of three forces its three
residue-class sizes to be a permutation of `(4,6,6)`. Translation puts the
size-four class at residue zero and then puts an omitted point at coordinate
zero. An exhaustive integer backtrack over the resulting search space visits
43,268,109 nodes and finds no `R`. Unit multiplication transfers this result
to all ten survivors.

Thus there is **no strong cyclic solution at `k=4`**, although general
conference matrices of order 66 exist. This is a finite falsifier for the
strong cyclic/self-indexed ansatz at one parameter, not an obstruction to
the unrestricted two-fiber block system or to an infinite subsequence that
skips `k=4`. The exact audit and stable result are
[`audit_k4_strong_cyclic_nonexistence.cpp`](../computations/audit_k4_strong_cyclic_nonexistence.cpp)
and
[`k4_strong_cyclic_nonexistence.json`](../computations/results/k4_strong_cyclic_nonexistence.json).

Exhaustive enumeration of all `2^17` Boolean states modulo global negation
gives cap 33 for the saved order-18 signing (204 projective maximizers). This
is an exact cap of that signing and hence a rigorous `M_18 <= 33` construction
bound; it is not a proof that `M_18=33`. The order-38 non-Paley certificate
supplies finite structural diversity but, without a scalable cap theorem,
does not itself count as a landing result.

## 4. What the block identities alone prove about Boolean cap

There is a sharp certificate barrier. Write arbitrary Boolean fiber spins as

```math
x=\alpha\mathbf1+u,qquad y=\beta\mathbf1+v,qquad
u,v\perp\mathbf1,qquad
\alpha={\mathbf1^{\mathsf T}x\over s},\quad
\beta={\mathbf1^{\mathsf T}y\over s}.               \tag{8}
```

Regularity makes the fiber-constant plane and its orthogonal complement
invariant. On the first plane the quadratic energy is exactly

```math
H_0=s\{k(\alpha^2-\beta^2)+\alpha\beta\}.           \tag{9}
```

On the complement, the conference identity gives operator norm `sqrt(Q)`.
Consequently every such two-fiber conference satisfies

```math
|H_S(x,y)|\le
s|k(\alpha^2-\beta^2)+\alpha\beta|
+{s\sqrt Q\over2}(2-\alpha^2-\beta^2).             \tag{10}
```

This envelope never exceeds the ordinary spectral bound `s sqrt(Q)`, since

```math
|k(\alpha^2-\beta^2)+\alpha\beta|
\le {\sqrt Q\over2}(\alpha^2+\beta^2);             \tag{11}
```

(11) is just the operator norm of `T_k/2`.

More importantly, maximizing (10) over Boolean means cannot certify a
leading improvement. Since `s` is odd, choose spins with
`alpha=beta=1/s`. Their right-hand side is

```math
s\sqrt Q-{\sqrt Q-1\over s}.                        \tag{12}
```

Thus even the best use of (8)--(11) can shave at most `O(s^(-1/2))` in
absolute energy from a quantity of order `s^(3/2)`. The irrational quotient
eigenline can give strict finite-order slack, but not a sub-`1/2` leading
constant.

This is a **certificate barrier**, not a falsification of the cyclic family.
The actual Boolean cap can be much smaller than the right-hand side of (10).
Proving that requires a nontrivial-module Boolean inequality using the
specific autocorrelation/difference structure, not merely conference
orthogonality and the two-dimensional quotient.

## 5. Exact open problem and falsifier

The clean scalable construction question is now:

> Construct **non-Paley** self-indexed ASDS pairs (7) for an infinite sequence
> of `k`, or prove a scalable obstruction to them, and obtain a Boolean cap
> bound using more than the two-fiber spectral decomposition.

Existence alone does not prove a useful landing theorem: all members are
conference matrices and (10) returns the `1/2` spectral constant. The prime
Paley branch is now explained and is asymptotically saturating, so further
Paley examples are not progress. The saved non-Paley `k=3` solution and the
certified `k=4` nonexistence show that the strong cyclic subclass is
arithmetically selective rather than an all-parameter construction. The next
existence target must therefore supply an explicit infinite parameter set;
another isolated solver search would not establish scalability.

Conversely, nonexistence in the cyclic subclass would falsify only this
polynomial-state candidate, not the general two-fiber block system (1). The
self-indexed ASDS condition remains a concrete, independently checkable
target strictly smaller than unrestricted bridge optimization, but a useful
research continuation must connect its non-Paley solutions to Boolean cap or
cross-order landing rather than only enumerate more conference gauges.
