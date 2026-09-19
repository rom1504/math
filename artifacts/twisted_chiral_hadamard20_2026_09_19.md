# A short arithmetic certificate for the order-20 witness

Status: proved; independently reconstructed by the adversarial agent.
This is a finite witness certificate, not a proof of the global value M20.
The underlying excess bound is classical: M. R. Best, *The excess of a
Hadamard matrix*, CWI ZW87/76 (1976), Section6, Theorem5, proves that the
maximum excess of an order-20 Hadamard matrix is80. The argument below
uses symmetry to give a self-contained certificate for our construction.
Primary source: https://ir.cwi.nl/pub/6858/6858D.pdf . No external novelty
is asserted.

## 1. Symmetric Hadamard quadratic excess

If H is a symmetric order-20 sign matrix with H²=20I, then

```math
 \max_{z\in\{\pm1\}^{20}}|z^THz|\le80.             \tag{1}
```

Fix z and set K=diag(z)Hdiag(z), r=K1. Then K²=20I and
sum r_i²=400. Every r_i is even, and all r_i have the SAME residue
modulo4. Indeed, any two Hadamard rows differ in exactly10 positions;
their sums differ by twice a sum of10 signs, hence by0 modulo4.

If r_i=2 modulo4, each (r_i-2)(r_i-6)>=0. Summing gives
400-8 sum r_i+240>=0, hence sum r_i<=80.

If r_i=0 modulo4, let q=r/4. Then q is integral, sum q_i²=25,
and Cauchy gives sum q_i<=sqrt(500)<23. Its parity is odd because
sum q_i=sum q_i² modulo2. Thus a sum greater than20 would have
to be21. In that case

```math
 \sum_i(q_i-1)^2=3,\qquad \sum_i(q_i-1)=1,
```

forcing q=1+e_a+e_b-e_c for three distinct coordinates. But Kr=20·1
then gives

```math
 K(e_a+e_b-e_c)=5\,\mathbf1-r.
```

At coordinate c the right side is5, while the left side is a signed
sum of three unit entries and has absolute value at most3. Contradiction.
So sum q_i<=19 in this residue class, in particular sum r_i<=80.
Apply the same reasoning to -H for the negative side. This proves(1).

## 2. The conference double

Let C be any symmetric hollow conference matrix of order10:
C has off-diagonal signs and C²=9I. Set

```math
 H=\begin{pmatrix}C-I&C+I\\ C+I&-C+I\end{pmatrix},
 \qquad
 D=\begin{pmatrix}C&C+I\\ C+I&-C\end{pmatrix}.
```

The commuting blocks give H²=20I. Its diagonal is(-1_10,+1_10),
whose sum is zero, so z^TDz=z^THz for every Boolean z.
Consequently Q(D)<=40. A stored spin attaining40 certifies equality
for the explicit conference witness. This is exactly the requested
twisted family with B=C and d=1, without any relaxation of the signs.

Thus M20<=40. The separate exhaustive fixed-child computation gives
F(A)=44 for BOTH globally minimizing order-10 child classes, each of
cap13. Our conference child has cap15. This establishes a finite failure
of exact optimal-child inheritance, not a global lower bound at order20.

## 3. A reusable row-feedback restriction

The same elementary mechanism works at any symmetric Hadamard order N.
For K=diag(z)Hdiag(z), r=K1 and any real k, put e=r-k1. Then

```math
 K e=(N-k^2)\mathbf1-k e,
 \qquad
 \max_i|N-k^2-k e_i|\le\|e\|_1.                   \tag{2}
```

Together with sum r_i²=N² and the row-sum lattice, (2) can exclude
formally feasible nearly regular row profiles. It is only a necessary
constraint. No fixed asymptotic discount from N^(3/2) follows here.

## 4. Why this finite discount does not tensor automatically

If H is ANY symmetric Hadamard matrix of order N, the sign vector
v=vec(H) satisfies

```math
 (H\otimes H)v=Nv.
```

Therefore v^T(H⊗H)v=N³, attaining the operator bound at order N².
If tr(H)=0, hollowing H⊗H changes every Boolean energy by
tr(H⊗H)=tr(H)²=0. Its hollow signing thus has EXACT cap N³/2,
or normalized cap1/2. In particular the cap40 order-20 example does
not retain its strict normalized discount under this tensor square.
This is a limitation of that operation, not of every recursive chiral lift.

## Dependencies and status

The proof of the upper certificate requires only conference orthogonality,
integer row sums, symmetry and the displayed algebra. Equality uses the
explicit spin certificate in the accompanying computation. The full-family
44 statements have separate exhaustive programs and census certificates.
Neither global M20=40 nor an asymptotic convergence claim is made.
