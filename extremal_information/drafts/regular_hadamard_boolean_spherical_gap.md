# A fixed Boolean--spherical gap at bounded total port mass

**Status.** Rigorous scalable counterexample with an exact finite/tensor
verifier.  This settles the proposed uniform Boolean-recovery question for
the regular-Hadamard SA.3 trust carrier in the negative.

The obstruction already occurs with two ports, a single common Hadamard
factor, and total repeated-port mass `mp/r=1`.  It is not a generic
Grothendieck loss: both optima and their difference are computed exactly.

## 1. Setup

Let `H` be a symmetric entrywise sign matrix satisfying

```math
H^2=r^2I,
\qquad \operatorname{tr}H=0,                       \tag{BG.1}
```

and suppose it has two orthogonal Boolean top eigenvectors

```math
a,b\in\{+-1\}^n,
\qquad Ha=ra,
\qquad Hb=rb,
\qquad a^Tb=0.                                      \tag{BG.2}
```

Put `A=H-diag(H)`, an exact hollow signing.  Append two disjoint shores,
each containing an integer `m`
Boolean spins, and connect them to the old block by the repeated columns
`a` and `b`.  Before adding interactions among the new spins, endpoint
optimization gives the exact Boolean response

```math
\mathcal B_m(a,b)
=\max_{u\in\{+-1\}^n}
 \left\{|H_A(u)|+m|a^Tu|+m|b^Tu|\right\},            \tag{BG.3}
```

where, because `tr H=0`,

```math
H_A(u)={1\over2}u^THu.                              \tag{BG.4}
```

The SA.3 spherical trust relaxation of the same Boolean representation is

```math
\mathcal S_m(a,b)
=\max_{\|u\|_2^2=n}
 \left\{\left|{1\over2}u^THu\right|
       +m|a^Tu|+m|b^Tu|\right\}.                    \tag{BG.5}
```

As in SA.3, (BG.5) relaxes the trace-zero Boolean representation (BG.4),
not the hollow quadratic evaluated literally off the cube.

## 2. Exact gap theorem

### Theorem BG.1 (orthogonal top ports have a fixed trust gap)

Under (BG.1)--(BG.2), for every integer `m>=0`,

```math
\boxed{
\mathcal B_m(a,b)={rn\over2}+mn,}                   \tag{BG.6}
```

whereas

```math
\boxed{
\mathcal S_m(a,b)={rn\over2}+\sqrt2\,mn.}           \tag{BG.7}
```

Consequently the exact Boolean--spherical integrality gap is

```math
\boxed{
\mathcal S_m(a,b)-\mathcal B_m(a,b)
=(\sqrt2-1)mn.}                                    \tag{BG.8}
```

#### Proof

Orthogonality of two Boolean vectors means that they agree on exactly half
the coordinates and disagree on exactly half.  Hence, for every endpoint
sign pair,

```math
\|\epsilon_1a+\epsilon_2b\|_1=n.                   \tag{BG.9}
```

Duality of `l_1` and the Boolean cube gives

```math
\max_{u\in\{+-1\}^n}(|a^Tu|+|b^Tu|)
=\max_{\epsilon_1,\epsilon_2}
  \|\epsilon_1a+\epsilon_2b\|_1=n.                 \tag{BG.10}
```

Also `|u^THu|/2<=rn/2` for Boolean `u`.  Thus the right side of (BG.6) is
an upper bound.  It is attained at `u=a`, because `a` is a top eigenvector,
`a^Ta=n`, and `b^Ta=0`.  This proves (BG.6).

On the Euclidean sphere,

```math
|a^Tu|+|b^Tu|
\le\max_{\epsilon_1,\epsilon_2}
 \|\epsilon_1a+\epsilon_2b\|_2\|u\|_2
=\sqrt2\,n.                                         \tag{BG.11}
```

Together with the spectral quadratic bound this proves that the right side
of (BG.7) is an upper bound.  It is attained at

```math
u={a+b\over\sqrt2}.                                  \tag{BG.12}
```

Indeed, (BG.12) has squared norm `n`, remains in the `+r` eigenspace, and
has `a^Tu=b^Tu=n/sqrt(2)`.  Thus (BG.7), then (BG.8), follow. `square`

The failure is geometric and exact.  The spherical carrier uses the
Euclidean midpoint of two Boolean poles; the cube pays the `l_1` rather than
the `l_2` support of their sum.

## 3. Bounded total repeated-port mass

There are `p=2` ports, so the collective mass parameter from GE.19 is

```math
c={mp\over r}={2m\over r}.                            \tag{BG.13}
```

Choose

```math
m={r\over2}.                                         \tag{BG.14}
```

Then `c=1`, while (BG.8) becomes

```math
\mathcal S_m-\mathcal B_m
={\sqrt2-1\over2}rn.                                 \tag{BG.15}
```

Thus bounded total port mass, which is sufficient for metric continuity of
the **spherical** response, does not make the Boolean integrality gap
`o(rn)`.  Even the smallest nontrivial common-factor Gram state

```math
G=R=I_2                                               \tag{BG.16}
```

has a fixed normalized gap.

More generally, fixing any `c>0` and taking `m=cr/2` gives gap

```math
{c(\sqrt2-1)\over2}rn.                               \tag{BG.17}
```

up to the harmless integral rounding of `m`.

## 4. A scalable regular-Hadamard family

Let `H_16` be the order-16 regular Walsh signing used in SA.4.  It has the
orthogonal Boolean `+4` eigenvectors `1` and

```text
v_0=(-,-,-,+; -,-,+,-; +,-,+,+; -,+,+,+).
```

For every `j>=1`, set

```math
H_j=H_{16}^{\otimes j},
\qquad n_j=16^j,
\qquad r_j=4^j,
\qquad
a_j=\mathbf1_{n_j},
\qquad
b_j=v_0\otimes\mathbf1_{16^{j-1}}.                  \tag{BG.18}
```

Then (BG.1)--(BG.2) hold and `m_j=r_j/2` is an integer.  Therefore

```math
{\mathcal S_{m_j}-\mathcal B_{m_j}\over n_j^{3/2}}
={\sqrt2-1\over2}                                   \tag{BG.19}
```

for every `j`.

The full parent has `N_j=n_j+2m_j=n_j+r_j` vertices.  Hence the gap in
total-order units tends to the same constant:

```math
{\mathcal S_{m_j}-\mathcal B_{m_j}\over N_j^{3/2}}
\longrightarrow {\sqrt2-1\over2}.                  \tag{BG.20}
```

## 5. Exact-sign completion does not remove the gap

To obtain a complete signing, fill all pairs among the `2m` auxiliary
vertices by any public hollow sign matrix `C`.  Let `B_C,S_C` denote the
exact Boolean response and the old-spin spherical relaxation with this
same auxiliary term retained.  Uniform cap Lipschitzness gives

```math
|B_C-\mathcal B_m|\le Q(C),
\qquad
|S_C-\mathcal S_m|\le Q(C).                          \tag{BG.21}
```

Since

```math
Q(C)\le {2m\choose2}=O(r^2)=O(n),                   \tag{BG.22}
```

the completed gap obeys

```math
S_C-B_C
\ge(\sqrt2-1)mn-O(n).                               \tag{BG.23}
```

At `m=r/2`, this is still
`((sqrt(2)-1)/2+o(1))n^(3/2)`.  No special auxiliary signing or optimizer
claim is required.

## 6. What this falsifies

1. There is no uniform `o(1)` normalized theorem rounding the
   regular-Hadamard **spherical trust value** to the exact Boolean response
   under the sole condition `mp/r=O(1)`.
2. Low rank does not help by itself: the exposed field lies in a
   two-dimensional top eigenspace.
3. A trust margin does not address this obstruction.  The gap is not caused
   by instability at `alpha=1/2`; it is the fixed `l_1/l_2` geometry of a
   continuous mixture of two Boolean poles.
4. Any positive recovery theorem must add a checkable Boolean-net property
   strong enough to approximate normalized sums of relevant eigenspace
   poles by actual Boolean vectors while preserving the quadratic channel.
   Such a property is strictly extra information beyond `(G,R)`.

This does not show that every Gram state has a gap, that `(G,R)` is
information-insufficient for every possible Boolean-response functional, or
that no smaller non-spherical Boolean quotient exists.  Such a claim would
require two equal-Gram states with separated Boolean responses.  What is
ruled out is the current spherical trust value as a uniform Boolean
surrogate, even in the bounded-total-mass regime where its own metric entropy
and gluing laws are favorable.

## 7. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_boolean_spherical_gap.py
```

The verifier exhausts all old Boolean spins at order 16 for `m=r/2` and
`m=r`, checks the exact formulas, and verifies the tensor eigenvector,
orthogonality, `l_1/l_2`, and normalized-gap identities through order 4096.
