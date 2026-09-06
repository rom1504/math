# A seed-only conditional tail need not have the weave's extensive speed

Date: 2026-09-06. A finite actual-signing counterexample to a uniform
conditional-tail proposal. It concerns fixed bases after exposure, with
only the independent seed signs left random. It does not contradict the
proved annealed bound that also averages randomized recursive bases.

## 1. Statement and precise scope

For every Hadamard order `t`, put `m=64t`, `k=62t`, and `N=mk`.
There are deterministic Hadamard bases in the rank-one weave, with this
exact retention `p=k/m=31/32`, such that independent uniform seed
off-diagonal signs and fixed seed diagonal `+1` have an event of probability

```math
2^{-2016t}=2^{-(63/2)m}                                    (1)
```

on which the retained hollow signing `A` has

```math
Q(A)=\frac{(m+1)N}{2},\qquad
\frac{Q(A)}{N^{3/2}}\longrightarrow\frac1{2\sqrt{31/32}}>\frac12.
```

It has a genuine global absolute maximizer whose aligned local fields
are exactly constant. Thus a high-cap, low-field-dispersion event cannot
have a uniform upper probability `exp(-c m^2)` after arbitrary bases
are exposed. The surviving seed signs really are independent before the
event is imposed; the failure is not a misuse of conditional independence.

IMPORTANT SCOPE: these examples have normalized mean field at the extreme
spectral edge, up to the hollow diagonal correction. They do NOT disprove
a conditional tail additionally restricted to mean fields a fixed distance
inside the extreme spectrum, such as the middle-space net candidates.

## 2. A fixed base with a retained Boolean extreme eigenvector

Let `J4-2I4` be indexed by pairs `(i,a)` in `{0,1}^2`. It is the
order-two-fibre weave from

```math
H_0=\begin{pmatrix}1&1\\-1&1\end{pmatrix},\quad
H_1=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\quad
S_2=\begin{pmatrix}-1&1\\1&-1\end{pmatrix}.
```

For every function `f` on `{0,1}`, the vector `(-1)^i f(a)` is a
`-2` eigenvector of `J4-2I4`: its coordinates sum to zero. Therefore,
in the first five factors of `(J4-2I4)^(tensor 6)`, every vector

```math
(-1)^{i_1+\cdots+i_5}F(a_1,\ldots,a_5)
```

is a `-32` eigenvector. In the last factor use the full constant
eigenvector, of eigenvalue `+2`. Choose `F` to be one except at
`(a_1,...,a_5)=(1,...,1)`, where it is zero. Reordering the six pairs
into a fibre coordinate `i=(i_1,...,i_6)` and a row coordinate
`a=(a_1,...,a_6)` gives a full order-4096 weave `W0`, with 64 fibres
of size64 and diagonal `+1`. Its tensor seed is `S0=S2^(tensor 6)`.

There are exactly two zero coordinates in each fibre, namely the rows
whose first five row bits are all one. Retain the other62 rows. The
restricted vector `x0` is Boolean, and its zero extension satisfies
`W0 x0_extended=-64 x0_extended`. Hence its retained matrix `K0`
satisfies `K0 x0=-64 x0`. The zero extension is an actual eigenvector,
so this does not approximate an eigenvector by silently dropping boundary
terms.

Let `H_i^0` denote the 64-by-64 tensor Hadamard in fibre `i` of this
base weave. These are explicit tensor products of `H_0,H_1` above.

## 3. Amplification retaining only O(m) effective random seed edges

Fix any Hadamard `H_t`. Index large fibres by `(i,u)`, where
`i in [64]` and `u in [t]`. In every such fibre use
`H_(i,u)=H_i^0 tensor H_t`, retaining the62 allowed base rows times
all `t` rows. This completely fixes the bases and selectors.

The only remaining randomness is the symmetric seed `S`, with iid
uniform off-diagonal signs and diagonal `+1`. Let `E_t` be the event

```math
S_{(i,u),(j,u)}=S^0_{ij}
\quad\hbox{for all }u\in[t],\ i<j.
```

It fixes exactly `t binom(64,2)=2016t` independent signs, proving (1).
All signs between different `u` groups remain unrestricted and random.

Use the explicit retained Boolean vector

```math
x_{(i,u),(a,v)}=x^0_{i,a}H_t(v,u).
```

Its row spectrum at fibre `(i,u)`, coordinate `(j,w)`, is exactly

```math
h_{(i,u)}(j,w)=t\,\boldsymbol1_{w=u}\,h_i^0(j),
\qquad h_i^0=(H_i^0[T,:])^Tx_i^0.                         (2)
```

This is column orthogonality of `H_t`. Formula (2) makes every cross-group
seed sign disappear from the field of this spin. On `E_t`, the surviving
sum in each group is the base eigenvector equation, so

```math
Kx=-64t\,x=-m x.                                         (3)
```

All diagonal entries of `K` are `+1`, because the seed diagonal is
`+1`. Consequently its hollowing `A=K-I` satisfies `Ax=-(m+1)x`.
Since the full weave has spectral norm `m`, principal compression gives
`||K||op<=m` and `||A||op<=m+1`. The Boolean vector in (3) attains that
bound, proving the exact global absolute cap in Section1. After aligning
the negative objective sign, every local field is `(m+1)/sqrt N`.

## 4. Conditional conclusion and what remains open

For every fixed `c<1/(2 sqrt(31/32))`, these fixed bases satisfy

```math
\Pr_S\{Q(A)>cN^{3/2}\ \hbox{and an absolute maximizer has }d_1=0\}
\ge2^{-(63/2)m}
```

for sufficiently large `t`. Therefore no basis-uniform `m^2`-speed
estimate for that event can replace the recursive-base average.

The exhibited bases may be rare or absent under additional restrictions
on a specified recursive ensemble; their relevance is to a purported
uniform conditional estimate for arbitrary revealed Hadamards. A proof
using typical exposed bases, an explicit interior mean-field restriction,
or a partially exposed recursive ensemble remains possible. Such a proof
must state and establish its conditional law separately.

Reproduction:
`computations/transfer_adversary_seed_only_conditional_tail_2026_09_06.py`
checks the exact field identities without building the large full weave.
