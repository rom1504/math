# Equatorial switching and balanced weight-two codewords

2026-09-06. Independently derived finite theorem; director audit PASS in
`decisive_director_equatorial_audit_2026_09_06.md`.
This is an intrinsic all-order regularity statement, not convergence.

## 1. A near-zero Boolean value for every hollow quadratic form

Let A be any real symmetric zero-diagonal n-by-n matrix, n>=2, and set

```
H_A(x)=sum_(i<j) A_ij x_i x_j,
v=max_i sum_(j!=i) A_ij^2,
L=sqrt(2 v log(7200 n)).
```

There exists a Boolean x with |H_A(x)|<=L. If v=0 this is immediate.

First let sigma^2=E H_A(X)^2=sum_(i<j) A_ij^2 for a uniform Boolean X.
The fourth moment satisfies

```
E H_A(X)^4 <= 15 sigma^4.                              (1)
```

One way to verify (1) is to compare with the centered Gaussian quadratic
form of the same coefficients. All fourth-moment terms whose Gaussian
and Boolean moments differ have positive products (one edge four times,
or two doubled edges). Four distinct edges contribute only when they
form a four-cycle, with the same weight in both input laws. Thus the
Boolean fourth moment is at most the Gaussian fourth moment. Diagonalizing
the latter gives 3 sigma^4+3 tr(A^4)<=15 sigma^4, since
tr(A^2)=2 sigma^2. This argument permits signed real coefficients.

Interpolation of L1,L2,L4 gives E|H_A|>=sigma/sqrt(15). Since E H_A=0,
E(H_A)_+=E(H_A)_- >=sigma/(2 sqrt(15)). Cauchy--Schwarz therefore shows

```
P(H_A>0)>=1/60,  P(H_A<0)>=1/60.                       (2)
```

Take independent uniform Boolean vectors X,Y and the coordinate path
Z^(t)=(Y_1,...,Y_t,X_(t+1),...,X_n). Its i-th energy increment is

```
H_A(Z^(i))-H_A(Z^(i-1))=(Y_i-X_i) F_i,
F_i=sum_(j<i) A_ij Y_j + sum_(j>i) A_ij X_j.
```

Each F_i is a sum of independent signs of squared coefficient sum at
most v. No independence BETWEEN the F_i is claimed. Hoeffding and a
union bound give

```
P(max_i |F_i|>L)<=2n exp(-L^2/(2v))=1/3600.
```

By (2), the independent endpoint energies have opposite signs with
probability at least 2/60^2=1/1800. Therefore with probability at least
1/3600 the endpoints have opposite signs and every jump has magnitude
at most 2L. At a zero-crossing jump at least one endpoint has absolute
energy at most L, proving the theorem. The construction is a finite
randomized search with success probability bounded below independently
of the matrix and n; it makes no hypothesis about near-minimality.

## 2. Balanced signings have the same asymptotic minimum at every order

Now A is a hollow signing and N=binom(n,2). Let M_n^bal be the minimum
of Q(A)=max_x |H_A(x)| over signings with

```
|sum_(i<j) A_ij|<=1.
```

Every n admits such a signing. Then, with
L_n=sqrt(2(n-1) log(7200n)),

```
M_n <= M_n^bal <= M_n+L_n.                             (3)
```

Indeed start with an exact minimizer. The switching gauge
A'_ij=A_ij x_i x_j preserves Q and, by Section 1, has total edge sum
S with |S|<=L_n. If N is even, flip |S|/2 edges of its majority sign;
if N is odd, flip (|S|-1)/2. The resulting total is zero or +/-1.
An edge flip changes the cap by at most two, so the total cap increase
is at most |S|<=L_n. This proves (3), including the parity normalization.

In particular the normalized balanced and unrestricted sequences differ
by O(sqrt(log n)/n). They have identical liminf, limsup and convergence
status. The conclusion is all-order, not confined to powers of two.

The square-root scale cannot be removed for arbitrary switching classes.
For the all-positive signing and n=4r(r+1),
H_A(x)=((sum_i x_i)^2-n)/2. The sum is even, and the nearest permitted
squares to n are (2r)^2 and (2r+2)^2. Thus min_x|H_A(x)|=2r.
Here N is even, so every exactly balanced representative requires at
least r edge edits after switching. This does not concern near-minimizers;
it only shows that a universal exact-equator gauge is false and that
the proved universal estimate is sharp up to a square-root logarithm.

## 3. Exact code interpretation and its limits

The augmented cut code is

```
C_n^+={(c+b_i+b_j)_(i<j):c,b_i in F_2}
      =RM(1,n) restricted to {e_i+e_j:i<j}.
```

For the sign word a_ij=(-1)^f_ij,
d(f,C_n^+)=(N-Q(A))/2. Let rho_n be its unrestricted covering radius
and rho_n^bal the maximum distance to the code among words of weight
floor(N/2) or ceil(N/2). Equation (3) is exactly

```
rho_n-L_n/2 <= rho_n^bal <= rho_n.                     (4)
```

More generally EVERY coset has a switching representative within L_n/2
of the equatorial weight band, and hence within L_n/2 edge edits of a
balanced word. This is not row-sum regularity, spectral regularity,
ground-state isotropy, or an order-transfer construction. It does not
assert that exact minimizers themselves admit an exactly balanced gauge
without the additional edge edits.

The code mapping was already present in
`critical_scale_code_audit_2026_09_05.md`; the new part here is the
quantitative near-equator switching/edit estimate. Thus balanced
weight-two truth-table constructions may be compared to the original
optimum without a leading asymptotic loss, but their recursions must
still control actual cap to yield convergence.

## 4. Initial architecture screen

The initial independent ideas considered were: cut-code shortening;
dual weight-enumerator bounds; sparse spherical Reed--Muller recursions;
equatorial representative regularity; quadratic Boolean Sidon constants;
quadratic unimodular polynomial constructions; graph-design composition;
rare-selector large deviations; exact deep-hole exchange; and stabilized
bilinear norms. The three prioritized tests were spherical coding,
equatorial regularity and exact-constant Banach-space literature.

Shortening retains a boundary completion problem rather than automatically
giving a powered recurrence. Generic enumerator sufficiency is already
falsified in the older code audit. Graph-design composition with fixed
seed blocks enters a Gaussian first-moment regime and has not retained
the actual seed cap. Exact deep-hole exchange overlaps the director's
current convex-geometric work and is not duplicated here.

The present theorem is the first concrete positive output of this screen.
It is not, by itself, a primary convergence mechanism.

## 5. Finite verification

`computations/decisive_audit_equatorial_switching_checks_2026_09_06.py`
passes 210 weighted/sign fourth-moment cases, 105 exact switching/edit
cap checks and 36 certified coordinate paths at orders 24,48,96.
It is self-contained and writes no generated artifacts. These finite
checks accompany, and do not replace, the finite theorem's proof.
