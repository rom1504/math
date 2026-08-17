# Relative involution algebra controls graph-composition energy

Status: **rigorous theorem, independently audited after the quadratic-loss
repair**.  This abstracts
the Walsh commutation example into a quantitative theorem for any symmetric
orthogonal child/bridge pair.  It is a spectral/Boolean composition bound,
not a complete response quotient.

## 1. Model

Let `C,F` be real symmetric `N by N` involutions:

```math
C^2=F^2=I_N.
```

For a graph `G` on `k` vertices, with adjacency matrix `A_G`, put

```math
mathcal M_G(C,F)=I_k tensor C+A_G tensor F.               \tag{RIC.1}
```

For Boolean block spins `X=(x_1,...,x_k)`, define, at scale `lambda>0`,

```math
E_G(X)={lambda\over2}X^Tmathcal M_G(C,F)X.               \tag{RIC.2}
```

The child contribution is `lambda x_i^TCx_i/2`, and every undirected edge
contributes `lambda x_i^TFx_j` once.

## 2. A robust two-sided law

### Theorem RIC.1 (commutator section versus anticommutator ceiling)

Write

```math
rho_G=||A_G||_(2->2),
\qquad
eta_+=||CF+FC||_(2->2).                                  \tag{RIC.3}
```

Then every graph obeys the Boolean spectral ceiling

```math
\max_(X in {+-1}^{kN})E_G(X)
\le {lambda kN\over2}
\sqrt{1+rho_G^2+rho_G eta_+}.                            \tag{RIC.4}
```

Suppose in addition that `G` is bipartite with classes `L,R`, and there is a
Boolean vector `s` such that

```math
Cs=s,
\qquad Fs in {+-1}^N.                                    \tag{RIC.5}
```

Put `eta_-=||CF-FC||_(2->2)`.  Then

```math
\max_XE_G(X)
\ge lambda N\left{
 {k\over2}+|E(G)|-{eta_-^2\over4}\sum_(H in cc(G))\min(|L\cap H|,|R\cap H|)
 \right}.                                               \tag{RIC.6}
```

In particular, exact commutation (`eta_-=0`) gives the exact termwise value

```math
\max_XE_G(X)=lambda N\left({k\over2}+|E(G)|\right).       \tag{RIC.7}
```

If `G` is `r`-regular and bipartite with `r>0`, comparison with a second pair
having anticommutator norm at most `eta` yields an extensive energy gap of at
least

```math
{lambda kN\over2}
\left[1+r-\sqrt{1+r^2+r eta}\right].                    \tag{RIC.8}
```

whenever the first pair has an exact commuting Boolean section.  The bracket
is positive for `eta<2`; thus the extensive separation is stable under a
fixed anticommutation error.

### Proof

Squaring (RIC.1) gives the exact identity

```math
mathcal M_G(C,F)^2
=(I_k+A_G^2) tensor I_N
 +A_G tensor(CF+FC).                                     \tag{RIC.9}
```

Therefore

```math
||mathcal M_G(C,F)||^2
\le1+rho_G^2+rho_G eta_+.                                \tag{RIC.10}
```

In fact the two terms in (RIC.9) commute, so the exact spectral formula is

```math
||mathcal M_G(C,F)||^2
=\max_(alpha in spec(A_G),sigma in spec(CF+FC))
 (1+alpha^2+alpha sigma).                                \tag{RIC.10a}
```

For a regular bipartite graph the adjacency spectrum contains both
`+-rho_G`; hence (RIC.10) is exact, not a triangle-loss estimate, for the
operator ceiling used in (RIC.8).

Every Boolean `X` has squared Euclidean norm `kN`, so the Rayleigh bound in
(RIC.2) proves (RIC.4).

For the lower bound, assign `s` to one color class and `Fs` to the other.
Every bridge is saturated because

```math
s^TF(Fs)=N.                                              \tag{RIC.11}
```

The `s` child has value `s^TCs=N`. Put `t=Fs`. Moreover

```math
Ct-t=(CF-FC)s=[C,F]s.
```

Symmetry and `C^2=I` give the exact sectionwise identity

```math
N-t^TCt={1\over2}||Ct-t||_2^2
={1\over2}||[C,F]s||_2^2.                               \tag{RIC.12}
```

Thus one imperfect child loses at most `eta_-^2N/2` in its quadratic form,
or `lambda eta_-^2N/4` in the energy. Orient every connected component so
that `Fs` lies on its smaller color class and sum the terms, proving
(RIC.6). When `eta_-=0`, the separate bounds
`x^TCx<=N` and `x^TFy<=N` show that the constructed value is also the global
upper bound, proving (RIC.7).  Finally, an `r`-regular graph has
`rho_G=r` and `|E|=kr/2`; subtract (RIC.4), with `eta_+<=eta`, from (RIC.7)
to obtain (RIC.8). `square`

## 3. Interpretation and limits

The pair of numbers

```math
(||[C,F]||,||\{C,F\}||)
```

is not a complete state.  It supplies two robust certificates:

- a near-common optimizer section when the commutator is small and a Boolean
  pole survives transport;
- a global composition ceiling when the anticommutator is small.

This cleanly generalizes the exact Walsh parity bit. There `eta_-=0` in the
even class and `eta_+=0` in the odd class. In that Walsh specialization the
isolated spectra are identical while the relative algebra changes the
extensive optimum; matching isolated spectra are not asserted for arbitrary
pairs in Theorem RIC.1.

The quantitative lower bound is conditional on the exact Booleanity of
`Fs`. That property is brittle under generic matrix perturbation. Thus the
theorem is stable in the relative algebra *given a transported Boolean
pole*, not a general perturbative realization theorem.

The result does not compute the Boolean maximum in the intermediate-angle
regime, provide a finite response update, or prove that commutator data alone
is sufficient under varying words `C_1,...,C_k`. A pair of families agreeing
in a proposed relative-algebra state but having separated Boolean graph
responses is the relevant test of that state's incompleteness; it would not
contradict RIC.1.
