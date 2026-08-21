# Exposed-shell sparse repair: distance and information no-go

**Status.**  Rigorous proof-class obstruction and exact finite audit.  The
natural extension of the sparse-repair condition (DR.25) from exposed
shells to arbitrary exact quadratic-cap shells is false on a scalable
family.  On a genuinely exposed shell, almost every root is intrinsically
farther than the proposed repair radius.  Finally, an unconstrained
low-label injection is exactly equivalent to the desired bad tail, rather
than a weaker lemma.  None of these statements falsifies DR.25 in its
original exposed-shell scope or falsifies ERSR: the scalable trap shell is
not known to be selected at any fixed temperature with positive mass.

Throughout,

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|,
\qquad D_N={N\choose2}.
```

## 1. Exposed shells are overwhelmingly far from the lower layer

Fix `beta>0`, let `q` be an exact cap selected by the exposed-layer
variational problem, and put

```math
\mathcal E=\{A:Q(A)=q\},
\qquad
\mathcal L=\{A:Q(A)\le q-2\}.
```

Proposition DR.2 gives, with
`rho=|L|/(|E|+|L|)`,

```math
\rho\le e^{-2\beta\sqrt N}.                           \tag{SR.1}
```

Let `V_D(k)=sum_(i<=k) binom(D,i)`, and write
`d_H(A,L)=min_(B in L)d_H(A,B)`.

### Proposition SR.1 (exposed-shell distance bound)

For every integer `k>=0`,

```math
\boxed{
\Pr_{A\sim U(\mathcal E)}
 \{d_H(A,\mathcal L)\le k\}
\le V_{D_N}(k){\rho\over1-\rho}.}                    \tag{SR.2}
```

Consequently, if

```math
k=\left\lfloor{\gamma\sqrt N\over\log N}\right\rfloor,
\qquad 0\le\gamma<{4\beta\over3},
```

then

```math
\boxed{
\Pr\{d_H(A,\mathcal L)\le k\}
\le
\exp\left[-\left(2\beta-{3\gamma\over2}+o(1)\right)
\sqrt N\right].}                                    \tag{SR.3}
```

In particular, the bound is `exp(-2 beta sqrt(N)+o(sqrt(N)))` for
`k=o(sqrt(N)/log N)`.

**Proof.**  The union of radius-`k` Hamming balls about `L` has size at
most `|L|V_(D_N)(k)`.  Intersect with `E` and divide by
`|E|=(1-rho)(|E|+|L|)` to obtain (SR.2).  For the displayed `k`, the
standard binomial-ball estimate gives

```math
\log V_{D_N}(k)
\le k\log(eD_N/k)
=\left({3\gamma\over2}+o(1)\right)\sqrt N.
```

Combine this with (SR.1). `square`

Thus sparse-repairable parents are already an exponentially negligible
part of the exact exposed shell.  A theorem asserting that every bad
incidence has such a repair is a strong support-localization theorem, not a
generic consequence of exposure.

## 2. Low-label injection without structure is the tail bound itself

Let `S_m` be the family of `m`-subsets, and let
`B subset E times S_m` be any permutation-invariant bad-incidence
relation.  Write

```math
p={|\mathcal B|\over|\mathcal E||\mathcal S_m|}.
```

### Proposition SR.2 (capacity equivalence)

For a positive integer `R`, there exists an injection

```math
(A,S)\longmapsto(R(A,S),S,\kappa(A,S))
\in\mathcal L\times\mathcal S_m\times[R]             \tag{SR.4}
```

if and only if

```math
\boxed{
p\le R{|\mathcal L|\over|\mathcal E|}
  =R{\rho\over1-\rho}.}                              \tag{SR.5}
```

**Proof.**  Permutation invariance makes the number of bad parents the same
for every fixed `S`, namely `p|E|`.  The target fibre at that same `S` has
capacity `R|L|`, proving necessity.  If the inequality holds, enumerate the
two finite fibres and inject the first into the second independently for
each `S`. `square`

Therefore the low-label condition (DR.22), if no locality or algebraic
constraint is imposed on `R`, is exactly the desired rare-tail estimate.
It is not a strict reduction.  Local sparse repair supplied concrete
structure, but Proposition SR.1 shows how exceptional that structure must
be, and the next section shows that it can fail at the full project scale.

## 3. A scalable conference trap with a bad comparable restriction

Let `r` be an odd prime power, `F=GF(r^2)`, `K=GF(r)`, and let `C` be the
symmetric square-field Paley conference signing on `F union {infinity}`.
Its order and cap are

```math
N=r^2+1,
\qquad C^2=r^2I,
\qquad Q(C)=q={rN\over2}.                             \tag{SR.6}
```

We first isolate the general mechanism.

### Proposition SR.3 (edge-balanced ground-orbit trap)

Suppose an order-`N` signing `A` of cap `q` admits a probability law on
oriented absolute ground cuts `v`.  Concretely, if `x` is a Boolean ground
spin and `sigma=sign(H_A(x))`, then `v_e=sigma x_i x_j` on
`e={i,j}`, so that `sum_e a_ev_e=q`.  Suppose this law satisfies

```math
\mathbb E[a_ev_e]={q\over D_N}
\quad\hbox{for every edge }e.                        \tag{SR.7}
```

If `A^F` is obtained by flipping `k` edges, then

```math
\boxed{
Q(A^F)\ge q\left(1-{2k\over D_N}\right).}           \tag{SR.8}
```

In particular,

```math
Q(A^F)\le q-2\quad\Longrightarrow\quad k\ge {D_N\over q}. \tag{SR.9}
```

**Proof.**  At an oriented ground cut, flipping `F` changes its energy from
`q` to `q-2 sum_(e in F)a_ev_e`.  Average and use (SR.7); one member of the
ground law attains at least the average.  If `k<D_N/q`, the result is
strictly greater than `q-2`.  Every cap is congruent to `D_N` modulo two,
so it is then at least `q`. `square`

For the conference signing, (SR.7) holds with value exactly `1/r=q/D_N`.
A ground law is obtained as follows.  Fix `t in F setminus K`, choose a
square `u in F^*` uniformly, and independently choose uniformly
`f:K->{+-1}` with `sum f=1`, put

```math
x_{infinity}=1,
\qquad x_z=f(b)\quad\hbox{when }uz=a+bt.
```

Multiplication by a square preserves the Paley signing, and the standard
subfield character sum gives `Cx=rx` after pulling this vector back through
that multiplication.  For a finite edge of
negative Paley sign, its endpoints always lie in different fibres and
`E[f(b)f(d)]=-1/r`, giving oriented load `1/r`.  For a positive finite
edge, a fraction `2/(r+1)` of square multipliers place its endpoints in one
fibre; hence its average load is

```math
{2\over r+1}-{r-1\over r+1}{1\over r}={1\over r}.
```

For a border edge, `E f(b)=1/r`.  Thus every edge has the required load.
Since `D_N/q=r`, Proposition SR.3 proves

```math
\boxed{
d_H\bigl(C,\{B:Q(B)\le q-2\}\bigr)\ge r,}          \tag{SR.10}
```

Since `r=sqrt(N-1)`, this is an `Omega(sqrt(N))` lower bound.

This trap also has exact bad restrictions at a comparable scale.  Choose
an even `h` with `h/r->1/2`, take the union `S` of `h` parallel affine
`K`-fibres, and put `m=rh`.  If the spin is constant with sign
`epsilon_b` on fibre `b`, the same subfield character sum gives

```math
H_{C[S]}(x)={r\over2}\left(hr-
                    \left(\sum_b\epsilon_b\right)^2\right).
```

Choose the fibre signs to be balanced.  The displayed energy is then

```math
H_{C[S]}(x)={rm\over2}.
```

Conversely, `||C[S]||_op<=||C||_op=r`, so this is the exact cap:

```math
Q(C[S])={rm\over2}.                                  \tag{SR.11}
```

Therefore

```math
\boxed{
{Q(C[S])\over m^{3/2}}\longrightarrow{1\over\sqrt2},
\qquad
{Q(C)\over N^{3/2}}\longrightarrow{1\over2}.}       \tag{SR.12}
```

The incidence is bad by a fixed margin, yet lowering the parent by one cap
step requires at least `r` flips, hence `Omega(sqrt(N))` flips.  This
disproves the natural
uniform extension of (DR.25) to arbitrary exact quadratic-cap shells.

It does **not** disprove the exposed-shell target.  The conference cap is
not known to be selected by the fixed-temperature variational problem.
The known bad fibre unions number only `exp(Theta(sqrt(N)))` among
`exp(Theta(N))` comparable subsets, and the conference equivalence orbit
has unknown density in its cap shell.  No positive exposed-shell bad mass
follows.

## 4. Exact order-eight distance audit

The finite exposed shell behaves differently.  Exhaustive breadth-first
search in the order-eight switching quotient gives

```math
\begin{array}{c|rrr}
\text{distance from cap 12 to cap at most 10}&1&2&3\\ \hline
\text{number of root-gauge signings}&97440&168840&70560.
\end{array}                                          \tag{SR.13}
```

Thus every cap-12 signing, including every bad restriction incidence for
child orders three through five, has an arbitrary-edge repair of radius at
most three.  By contrast, cap ten is the bottom layer and all 4,200 shell
points are unreachable from cap at most eight.  Restricting repair edges to
the child set also leaves most bad incidences unreachable.  These facts are
exact finite evidence, not an asymptotic claim.

The computation is reproduced by

```bash
.venv/bin/python computations/audit_exposed_shell_repair_distances.py \
  --output computations/results/exposed_shell_repair_distance_audit.json
```

## 5. Consequence for the cross-order campaign

The sparse-repair implementation does not improve the actual defect
`E_(m,n)(beta)=O_beta(N)`:

1. unconstrained low-label injection is equivalent to the desired tail;
2. almost all exposed-shell roots are outside the proposed sparse radius;
3. natural exact cap shells contain bad roots requiring at least an
   `Omega(sqrt(N))` repair radius;
4. exposure or positive shell mass for those roots is unproved.

Accordingly this is a **STRIKE**, not Level 6 and not a RESET.  The
permanent SML remains unchanged.  Any surviving rare-event route must
control the bad-incidence mass directly on the overwhelming Hamming-far
part of the actual exposed shell; calling that control a remote injection
or sparse repair does not make it a weaker theorem.
