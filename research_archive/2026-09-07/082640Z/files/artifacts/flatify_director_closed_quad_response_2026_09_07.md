# Exact Boolean response of a constant-diagonal closed-quad trade

Date: 2026-09-07. Status: elementary theorem, root reconstruction;
finite replay accompanies the proof. This is a local identity, not a descent
or convergence theorem.

Let S be symmetric with every entry a sign, including its diagonal, and let
I be four indices. Suppose the product of the four entries in every column
of S[I,:] is +1, and S[I,I] has constant diagonal d. Put R=J_4/2-I_4,
extended by the identity outside I. Suppose R S R is again a sign matrix.
Then S[I,I]=d vv^T for an even-parity sign vector v.

Proof: after multiplying the principal block by d, its diagonal is +1 and
the graph of negative off-diagonal entries has even degrees. On four
vertices it is empty, a triangle, or a four-cycle. The triangle case has
total matrix sum 4 (up to sign), and gives entries in {0,+2,-2} in R S R,
so is excluded. The empty and four-cycle cases are exactly vv^T with
even-parity v. Equivalently this can be checked on the 64 possible
off-diagonal principal signings. Notice that this classification does NOT
apply to the mixed-diagonal paired trades in H_2 tensor S.

Fix outside spins z. Write h=S[I,I^c]z and let e be their internal hollow
energy. For x in the four-spin cube the full hollow energy is

```math
e+h^T x+(d/2)[(v^T x)^2-4].
```

The even-parity four-spin vectors are permuted by R (the all-equal pair is
fixed and the other six vectors are negated); their internal energy is
also preserved because Rv=+v or -v. The odd-parity vectors are mapped
bijectively onto the eight vectors +/-2e_i. The internal energy of either
an odd-parity vector or +/-2e_i is zero. Thus, defining

```math
C(e,h)=max_{x: product_i x_i=+1}
 |e+h^T x+(d/2)[(v^T x)^2-4]|,
p_odd(h)=max_{x: product_i x_i=-1} h^T x,
```

the conditional caps before and after the trade are exactly

```math
max{C(e,h), |e|+p_odd(h)},
max{C(e,h), |e|+2 ||h||_infinity},                 (1)
```

respectively. If h has no zero coordinate, p_odd(h) is ||h||_1 when the
coordinate-sign product is -1, and ||h||_1-2 min_i |h_i| otherwise; the same
formula extends by setting the correction to zero when a coordinate is zero.
Maximizing (1) over outside z gives the exact global caps.

The closed-quad hypothesis constrains exterior column patterns but does not
force the favorable inequality 2||h||_infinity <= p_odd(h). Both directions
occur, for example h=(4,0,0,0) and h=(2,2,2,-2). Both are integer sums of
even-parity sign columns. Also R is an involution. Therefore no universal
cap descent follows simply from being an allowed trade. The unresolved
constructive question is simultaneous control of these competing branches
over outside spins for a sequence of trades.

The exact replay also supplies a 12-vertex full-sign example with cap 24
before the trade and 26 afterwards. Both full matrices are stored; the
reverse trade gives strict descent. This verifies actual global cap change,
not merely a conditional-field possibility. It is finite evidence about
trade direction, not a scalable obstruction for optimal children.

## Literature scope

[Crnkovic--Egan--Svob, arXiv:2511.07020](https://arxiv.org/html/2511.07020v1),
Theorems 6.1 and 6.3, give sufficient constant-column/zero-sum or orthogonal
restricted-column conditions for sign-preserving Hadamard switches. Their
Theorem 7.3 proves that a nontrivial complex Hadamard trade changes at least
the order-many entries. These statements concern Hadamard orthogonality,
not Boolean-cap descent or optimal signing transfer. None supplies the
missing inequality in (1). The local response calculation here also does
not assume that S is Hadamard.

Replay: `computations/flatify_director_closed_quad_response_2026_09_07.py`.
