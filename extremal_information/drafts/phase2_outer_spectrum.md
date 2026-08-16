# Phase 2 director derivation: outer spectra as a root-averaged quotient

Status: theorem draft for independent coding-theory audit.  This concerns the
covering-radius model, not the original signing convergence problem.

## 1. The object forced by the counterexample

Let `(X,d)` be a finite metric space and let `A` be a nonempty subset.  Write

```math
\delta_A(x)=d(x,A)=\min_{a\in A}d(x,a).
```

The full labeled map `x -> delta_A(x)` reconstructs `A` from its zero set.
For the root-averaged pressure experiment below, retain only the outer
distance polynomial

```math
O_A(z)=\sum_{x\in X}z^{\delta_A(x)}.                 \tag{OS.1}
```

In the binary Hamming cube `Q_n`, it has only `n+1` integer coefficients.
This is the distance-layer distribution of the code in the ambient space,
not its inner codeword distance enumerator.

For `beta in R`, define the outer root pressure

```math
P_A(\beta)=\log\sum_{x\in X}e^{\beta\delta_A(x)}
=\log O_A(e^\beta).                                  \tag{OS.2}
```

## 2. Exact quotient and product algebra

### Theorem OS.1 (outer-spectrum response theorem)

Fix a finite integer-valued metric space of diameter at most `n`.

1. The polynomial `O_A` is the coarsest exact deterministic quotient for the
   complete pressure experiment `(P_A(beta))_(beta in R)`, up to a one-to-one
   recoding.
2. The covering radius is its degree and is the zero-temperature limit

   ```math
   \rho(A)=\deg O_A
   =\lim_{\beta\to\infty}{P_A(\beta)\over\beta}.      \tag{OS.3}
   ```

   Quantitatively,

   ```math
   \rho(A)
   \le {P_A(\beta)\over\beta}
   \le \rho(A)+{\log|X|\over\beta}
   \qquad(\beta>0).                                  \tag{OS.4}
   ```
3. For `l_1` metric products,

   ```math
   O_{A\times B}(z)=O_A(z)O_B(z),
   \qquad
   P_{A\times B}(\beta)=P_A(\beta)+P_B(\beta),       \tag{OS.5}
   ```

   and `rho(A times B)=rho(A)+rho(B)`.

#### Proof

The pressure is obtained from the polynomial by (OS.2).  Conversely, equality
of pressure on any interval implies equality of `O_A(e^beta)` there; two
real polynomials equal on an interval have the same coefficients.  Thus the
pressure answer function and `O_A` determine one another, proving exact
minimality in the quotient sense.

At least one root has distance `rho(A)` and every root has distance at most
`rho(A)`.  Bounding the sum in (OS.2) between `e^(beta rho(A))` and
`|X|e^(beta rho(A))` proves (OS.3)--(OS.4).

In the product metric,

```math
\delta_{A\times B}(x,y)=\delta_A(x)+\delta_B(y).
```

Expanding (OS.1) gives polynomial multiplication, and the other identities
follow. `square`

For `Q_n`, exact storage of `O_A` uses at most `O(n^2)` bits: there are
`n+1` coefficients, each between zero and `2^n`.  This is exponentially
smaller than the complete labeled distance map.  It is not claimed minimal in
bit length for a restricted code family.

### Linear-code interpretation

If `C <= Q_n` is linear, `delta_C` is constant on each coset and equals the
minimum weight of that coset.  Therefore

```math
O_C(z)=|C|\sum_{U\in Q_n/C}z^{\operatorname{wt}_{\min}(U)}. \tag{OS.6}
```

Thus `O_C/|C|` is the classical coset-leader weight distribution.  Its
largest occupied degree is the classical maximum coset-leader weight, namely
the covering radius.  The theorem adds an operational query-minimality and
composition interpretation; it does not rename this coding invariant as a
new polynomial.

## 3. It genuinely forgets the code

The outer spectrum is not a disguised labeled code table.  In `Q_4`, set

```math
A=\{0000,0001,0010,0011\},
\qquad
B=\{0000,0001,0010,0101\}.
```

Both ambient distance-layer counts are

```math
O_A(z)=O_B(z)=4+8z+4z^2.                              \tag{OS.7}
```

For `A`, this follows because it is a two-dimensional coordinate face and
distance is the weight of the other two coordinates.  Directly partitioning
the sixteen roots for `B` gives the same counts.

The ordered inner distance distributions are different:

```math
I_A(z)=4+8z+4z^2,
\qquad
I_B(z)=4+6z+4z^2+2z^3.                               \tag{OS.8}
```

Hamming isometries preserve inner distances, so the codes are not isometric.
Hence the state in (OS.1) forgets genuine code geometry while retaining every
query in (OS.2).

Cartesian powers preserve the collision and keep the codes nonisometric:
their common outer spectrum is `(4+8z+4z^2)^m`, while their inner enumerators
remain different powers of the polynomials in (OS.8).

There is a complementary separation from the previous phase.  Theorem 3.3
in `theorems.md` gives, for every fixed replica order `k`, code pairs whose
complete ambient unrooted membership/distance censuses agree through `k`
points while their covering radii differ.  Their outer spectra must therefore
differ in degree.  Thus an `O(n^2)`-bit outer state can retain information
missed by every fixed unrooted replica hierarchy.  The two states are
incomparable: the exact `Q_4` collision above shows that the outer spectrum
can forget inner data, while the parity hierarchy shows that fixed
inner-replica data can forget the outermost distance layer.  The obstruction
is not simply that every successful state must be enormous; different
experiments expose different coordinates of code geometry.

## 4. Exact scope boundary

The outer spectrum is sufficient for symmetric root pressure, covering
radius, and Cartesian-product composition.  It is not sufficient for:

- a labeled nearest-code query, which recovers `delta_A(x)` and hence `A`;
- puncturing or shortening at a named coordinate;
- unions, because `delta_(A union B)=min(delta_A,delta_B)` needs the joint
  distribution of the two distance maps; or
- an arbitrary bridge coupling to codewords rather than roots.

Thus the correct conclusion is not “rooted information compresses.”  It is:

> quotient the root labels only when the allowed environment is invariant
> under that quotient, and retain the outer spectrum when Cartesian products
> are the declared composition.

## 5. Literature coordinates

For linear codes, coset leaders and their maximum weight are classical.  A
primary modern algorithmic reference explicitly treating the coset-leader
weight distribution is Borges-Quintana, Borges-Trenard, Marquez-Corbella, and
Martinez-Moro,
[*Computing coset leaders and leader codewords of binary codes*](https://arxiv.org/abs/1211.5568).
The earlier external-distance framework is Delsarte,
[*Four fundamental parameters of a code and their combinatorial significance*](https://doi.org/10.1016/S0019-9958(73)80007-5).

The algebra (OS.5) is elementary and the zero-temperature estimate (OS.4) is
the usual log-sum-exp sandwich.  For this program the result supplies a
precise strict quotient for one declared pressure experiment between the
failed unrooted replica hierarchy and the lossless labeled distance map.  It
is classical generating-function structure, not by itself a new coding
mechanism.

## 6. Director checkpoint

This is not a new coding invariant, and it is not the minimal state for
covering radius under Cartesian products: the scalar `rho(A)` already
composes additively and answers that zero-temperature query.  It is a useful
Level-2 taxonomy result for three narrower reasons:

1. the previous phase had only the two extremes—failed inner/replica data and
   a labeled rooted map that reconstructs the code;
2. `O_A` is a strict intermediate state with a proved exact
   *complete-pressure* experiment;
3. it has an exact product algebra and a zero-temperature continuity bound.

The strongest next coding theorem is to identify a comparably strict quotient
for a richer composition than Cartesian product, such as bounded-state
trellis gluing.  The boundary response kernel from the independent feature-
growth report is the natural candidate; for linear trellises its realizable
kernel semigroup may be much smaller than the universal max-plus matrix cube.
