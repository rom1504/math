# Same-spin flat phases from projective nondisplaceability

Date: 2026-09-07. Status: complete consequence of the primary geometric
theorem, with the discrete losses stated separately. This is not an
original-problem convergence theorem or a Boolean half-floor.

## 1. The imported theorem is stronger than unitary Sinkhorn scaling

For every projective dimension d, the Clifford torus
`T={ [w]: |w_0|=...=|w_d| }` intersects every Hamiltonian image of
`RP^d` in `CP^d`. The following primary sources were checked directly.

- Entov--Polterovich, *Rigid subsets of symplectic manifolds*,
  [arXiv:0704.0105](https://arxiv.org/pdf/0704.0105): Theorem 1.2(iii),
  printed p.11, says a superheavy set meets a heavy set for the same
  idempotent. Example 1.22, pp.23--24, makes `RP^d` superheavy over `Z2`
  for the fundamental idempotent. Theorem 1.9, p.15, makes the special
  torus fibre superheavy for every nonzero idempotent. In projective
  space this fibre is the Clifford torus. Thus the coefficient field and
  idempotent agree; no cross-coefficient inference is used.
- Tamarkin, *Microlocal condition for non-displaceability*,
  [arXiv:0809.1584](https://arxiv.org/pdf/0809.1584), abstract and
  introduction p.1, states this mutual nondisplaceability in all
  dimensions and gives an independent microlocal proof.
- Alston, *Lagrangian Floer homology of the Clifford torus and real
  projective space in odd dimensions*,
  [arXiv:0902.0197](https://arxiv.org/pdf/0902.0197), Theorem 1.2 p.1,
  independently supplies nonempty intersection in every odd projective
  dimension. This alone covers every even matrix order. Its Floer
  computation is explicitly only in odd projective dimension; the
  all-dimensional claim is not inferred from that computation.

A projective unitary transformation is Hamiltonian: write its unitary
representative as `exp(iK)` with K Hermitian and use the projective
quadratic Hamiltonian associated to K. The arbitrary normalization of
the Fubini--Study form only rescales that Hamiltonian.

By contrast, ordinary unitary Sinkhorn/biunimodular results only supply
flat z and Uz. They do not by themselves impose `Uz=bar z`.

## 2. Exact same-spin theorem

For every complex symmetric unitary U of order n, there is a vector
`z in T^n` such that

```math
Uz=\bar z,\qquad z^TUz=n.
```

Proof. Write `U=VV^T` with V unitary. This special Takagi factorization
can be seen without singular-value theory: the antiunitary map
`w -> U bar w` is an involution; an orthonormal real basis of its fixed
space forms V. Nondisplaceability gives a projective point of
`V RP^(n-1)` on the Clifford torus. Choose its representative
`w=Vr`, with r real, rescaled so every coordinate of w has modulus one.
Then `U bar w=VV^T bar V r=Vr=w`. Set `z=bar w`.
The spectral upper bound is `|z^TUz|<=||z||^2=n`, so the result is exact.

In particular, if H is any real symmetric Hadamard of order n,

```math
\max_{|z_i|=1}|z^THz|=n^{3/2}.
```

For its hollowing `A=H-diag(H)`, the complex full quadratic is at least
`n^(3/2)-n`. More generally, for real hollow symmetric A and any complex
symmetric unitary U,

```math
\max_{|z_i|=1}|z^TAz|\ge nr-n\|A-rU\|_{op}\quad(r\ge0).
```

This approximate-unitary statement has no cap-only hypothesis: small
Boolean cap is not asserted to control its operator error.

## 3. Four-phase and Boolean losses

Put `C(A)=max_(|z_i|=1)|z^TAz|` and retain the original convention
`Q(A)=max_(x_i=+-1)|x^TAx|/2`. For real hollow symmetric A,

```math
Q(A)\ge\tfrac14 C(A).                              (1)
```

Indeed rotate a maximizing z by a common phase so `z^TAz=C(A)` is
real. Write `z=u+iv`. Both u and v lie in the real cube, and
`C(A)=u^TAu-v^TAv`. A hollow quadratic is affine in each real coordinate,
so its absolute value on the cube is at most `2Q(A)`. This proves (1).
Thus the geometric theorem gives only
`Q(H-diag H)>=n^(3/2)/4-n/4`, weaker than the already banked universal
lower bound. No claim that the factor in (1) is sharp is made.

For four phases `D4={1,i,-1,-i}`, its convex hull contains the disc of
radius `1/sqrt(2)`. Round coordinates independently with
`E q_i=z_i/sqrt(2)`. Hollowness gives the exact identity

```math
\mathbb E(q^TAq)=\tfrac12z^TAz,
\qquad \max_{q\in D4^n}|q^TAq|\ge\tfrac12 C(A).      (2)
```

This is a quadratic half-loss, not a lossless conversion to four phases.
The diagonal must first be removed, or its expectation treated separately.

## 4. Exact relation to the canonical H2 lift

Let `L(A)=[[A,A+I],[A+I,-A]]`. For q=a+ib in `D4^n`, set
`x=a+b`, `y=a-b`; x and y are Boolean and the supports of a,b partition
the coordinates. Direct expansion gives

```math
q_{L(A)}(x,y)
=\operatorname{Re}(q^TAq)+\operatorname{Im}(q^TAq)
  +\sum_i(a_i^2-b_i^2),                            (3)
```

where `q_L(w)=w^TLw/2`. Conjugating q reverses the imaginary part, so

```math
Q(L(A))\ge\max_{q\in D4^n}|q^TAq|-n
          \ge\tfrac12 C(A)-n.                      (4)
```

For a hollowed Hadamard this is only
`Q(L(A))>=(1/2)n^(3/2)-(3/2)n`, which after normalization by
`(2n)^(3/2)` is `1/(4sqrt(2))-o(1)`. Applying the geometric theorem
directly to an H2 lift that is itself within O(1) operator distance of a
scaled symmetric unitary gives the stronger but still insufficient
`1/4-o(1)` normalized bound. Neither proves a half-floor for selected H2
lifts. An arbitrary bounded-op signing has no asserted unitary completion
with vanishing relative error.

The real restriction is substantive even at order two:
`H2=[[1,1],[1,-1]]` has complex maximum `2sqrt(2)` but real Boolean full
quadratic maximum 2. Its flat conjugate witness is
`(exp(i pi/8),exp(-3i pi/8))`.

## 5. Products do not close the original order problem

Flat conjugate witnesses tensorize: if `Uz=bar z` and `Vw=bar w`, then
`(U tensor V)(z tensor w)=overline(z tensor w)`. The product theorem for
heavy sets (Entov--Polterovich Theorem 1.5) is consistent with this elementary
fact. It preserves continuous phases, not the discrete cube. It neither
removes the fixed rounding losses above nor supplies flat sign fillings of
the zero cross blocks in direct sums. Therefore no all-order Boolean
subadditivity or cap-preserving construction is obtained here.
