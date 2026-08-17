# Microcanonical hypographs as a rare-event compactness state

Status: literature-grounded scoped theory card, independently audited after
the repairs below. Promote only for leading exponential
density-of-states queries.  It does not preserve subexponential decorations
or full extremal point processes.

## 1. Object and topology

Let `K` be a compact metric descriptor space.  For a finite deterministic
landscape `X_n`, a speed `a_n -> infinity`, and a descriptor
`Q_n:X_n->K`, define on the finite image

```math
s_n(q)={1\over a_n}\log|\{x:Q_n(x)=q\}|,
\qquad\log0=-\infty.                                    \tag{MH.1}
```

More generally allow bounded-above upper-semicontinuous profiles.  Store the
closed hypograph.  Sequential hypograph convergence is

```math
q_n\to q\Longrightarrow\limsup_ns_n(q_n)\le s(q),
```

together with, for every `q`, one recovery sequence `q_n->q` satisfying

```math
\liminf_ns_n(q_n)\ge s(q).                              \tag{MH.2}
```

Sign reversal is epi/Gamma convergence. After imposing one common upper
normalization `s_n<=M`, use the compactified ordinate `[-infinity,M]`.
Closed downward hypographs in `K times [-infinity,M]` then form a compact
hyperspace. The standard references are
[Beer--Rockafellar--Wets (1992)](https://doi.org/10.1090/S0002-9939-1992-1119262-6)
and [Rockafellar--Wets, *Variational Analysis*, Section 7](https://sites.math.washington.edu/~rtr/papers/rtr169-VarAnalysis-RockWets.pdf).

## 2. Exact composition and recovery

### Theorem MH.1 (compact hypograph sup-convolution)

Let `K_1,K_2,K` be compact metric spaces, let
`m:K_1 times K_2->K` be continuous, and suppose bounded-above usc functions
obey

```math
f_n\xrightarrow{h}f,
\qquad g_n\xrightarrow{h}g.                             \tag{MH.3}
```

Define

```math
(f\star_mg)(z)=
\max_{m(x,y)=z}\{f(x)+g(y)\},                           \tag{MH.4}
```

with value `-infinity` on an empty fibre.  Then

```math
f_n\star_mg_n\xrightarrow{h}f\star_mg.                 \tag{MH.5}
```

For every continuous tilt `V:K->R`,

```math
\max_z\{(f_n\star_mg_n)(z)+V(z)\}
\longrightarrow
\max_z\{(f\star_mg)(z)+V(z)\}.                        \tag{MH.6}
```

Every cluster point of exact maximizers is a limit maximizer. If the limiting
maximum is finite, every limit maximizer has an asymptotically maximizing
recovery sequence.

#### Proof

For the hypograph upper bound, take a subsequence realizing the limsup and
maximizing decompositions `(x_n,y_n)`. Compactness supplies a convergent
subsequence. Continuity of `m` keeps its limit in the target fibre, and the
two hypograph upper bounds give (MH.5). For recovery, choose a maximizing
decomposition `(x,y)` of the limit value and separate recovery sequences for
`f_n` and `g_n`; their images `m(x_n,y_n)` converge to `m(x,y)`. This is the
hypograph lower bound. Apply the same two arguments after adding continuous
`V` to obtain (MH.6) and the optimizer statements. `square`

This is the compact sign-dual case of the epi-sum/infimal-convolution theorem
in Rockafellar--Wets, Proposition 7.56(a).  Noncompact versions need total
epi-convergence and a recession transversality condition; compactness removes
the escape-to-infinity issue.

### Corollary MH.2 (finite count-convolution)

Let `A_n(x),B_n(y)` be nonnegative integer multiplicities on finite grids,
assume their normalized log profiles hypographically converge as in (MH.3),
and

```math
C_n(z)=\sum_{m(x,y)=z}A_n(x)B_n(y).                     \tag{MH.7}
```

Let `D_n(z)` count decompositions with positive summand. If

```math
\sup_z\log\max(1,D_n(z))=o(a_n),                         \tag{MH.7a}
```

then the normalized log-count profile of `C_n` has the same hypograph limit
as the supremal convolution in (MH.4). On the common effective domain,

```math
0\le {1\over a_n}\log C_n(z)
-\max_{m(x,y)=z}\left{{\log A_n(x)\over a_n}
                       +{\log B_n(y)\over a_n}\right}
\le {\log\max(1,D_n(z))\over a_n}=o(1).                \tag{MH.8}
```

A zero count gives `-infinity` to both profiles and is excluded from the
ordinary subtraction in (MH.8).

If, in addition, the number of occupied descriptor fibres is
`exp(o(a_n))` (or an equivalent uniform coarse-bin Laplace principle holds),
then a recovered fibre event lying `Delta` below the maximum has uniform
probability `exp(-a_n Delta+o(a_n))`; it is not erased merely because it is
exponentially rare. Without that descriptor-complexity condition, exact
fibre hypographs need not control total mass.

Conversely, modulo an additive normalization, every bounded nonnegative usc
profile on a compact space has an **unstructured abstract-landscape**
finite-grid recovery sequence: approximate its hypograph by finite downward
hypographs and assign multiplicity `ceil(exp(a_n r))` at sampled top heights,
choosing the mesh with logarithmic cardinality `o(a_n)`.

## 3. Benchmark and falsifier

For the deterministic mean-field Blume--Emery--Griffiths model with
`x_i in {-1,0,1}`, the empirical occupation vector `L` has raw-count entropy

```math
s_{count}(L)=-\sum_jL_j\log L_j.                         \tag{MH.9}
```

The cited uniform-prior convention subtracts the harmless constant `log 3`,
giving `s_prob=s_count-log3`.

and energy

```math
u(L)=L_++L_- -K(L_+-L_-)^2.                             \tag{MH.10}
```

The microcanonical profile is the constrained supremum of (MH.9) over the
level set of (MH.10).  At parameters exhibiting ensemble nonequivalence it
has a nonconcave interval invisible to the canonical temperature family.
See [Ellis--Touchette--Turkington (2004)](https://doi.org/10.1016/j.physa.2003.11.028)
and the general large-deviation treatment in
[Ellis--Haven--Turkington (2000)](https://doi.org/10.1023/A:1026446225804).
This benchmark is naturally a density-of-states compactness problem rather
than a finite contextual quotient.

Bounded-temperature pressure can miss even the maximum. Fix `B>0` and
`0<delta<1/B`. Let landscape `A_n` have `ceil(e^n)` states of energy density
zero and one state of energy density `delta`; let `B_n` have only the bulk.
For every `|beta|<=B`,

```math
{1\over n}\log(\lceil e^n\rceil+e^{\beta n\delta})
-{1\over n}\log\lceil e^n\rceil\longrightarrow0,       \tag{MH.11}
```

uniformly, although the normalized maxima differ by `delta`. Their
hypographs differ by the isolated branch at `(delta,0)`.

The limit object has a sharp ceiling. One maximal state and `e^(sqrt(n))`
maximal states have the same speed-`n` profile. It cannot determine extremal
spacings, Poisson/Cox decorations, or genealogy-marked branches. Those need
a finer extremal process or a multi-speed hierarchy.

Algebraically this is a speed-sensitive contextual response roof under all
continuous descriptor tilts, not a wholly orthogonal theory. Its distinct
value is that the full usc hypograph retains nonconcave finite-rate branches
which a linear-temperature/Legendre summary convexifies away.
