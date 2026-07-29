# Analytic audit of the dependent compressed four-lifts

## Status

Two order-\(12\) compressed lifts, one for each switching class of
triangle signings, have been independently verified:

\[
Q(B_-)=Q(B_\triangle)=40<48=4^{3/2}Q(C_3).
\]

The finite contraction is real, but it does not furnish an iterative
amplification mechanism. The decisive facts are:

1. both witnesses are chiral/self-complementary;
2. any exact compressed lift of a chiral seed has normalized norm at
   least that of its seed;
3. more generally, a fixed-fibre compressed lift can only gain, at
   leading conditional-mean level, by recentering the midpoint of the
   seed's energy interval;
4. the available diagonal-fibre shift is \(O(ns^2)\), which is
   \(o((ns)^{3/2})\) for fixed \(s\) and \(n\to\infty\);
5. exhaustive search over every seed-independent uniform
   \(4\)-fibre gadget finds no common contraction of the two triangle
   profiles.

Thus the \(Q=40\) examples demonstrate genuinely dependent finite
absorption, but their improvement is a small-order centering effect.
The exact-compression composition law remains useful only if one can
prove normalized **nonincrease** on already centered, chiral profiles.

## 1. Independent exact verification

The first witness is the all-negative triangle matrix in
`dependent_profile_recovery_witness.json`. The second witness has seed
edge signs \((-1,-1,+1)\); both matrices are included in
`audit_dependent_4lift.py`.

Enumerating \(2^{11}\) spins after fixing the global sign gives:

| seed class | cross-block sums | \(Q\) | antipodal maximizers |
|---|---:|---:|---:|
| \((-1,-1,-1)\) | \((-8,-8,-8)\) | \(40\) | \(18\) |
| \((-1,-1,+1)\) | \((-8,-8,+8)\) | \(40\) | \(12\) |

For both witnesses the energy support is exactly

\[
\{-40,-36,-32,\ldots,32,36,40\},
\]

and the complete histogram is symmetric under \(E\mapsto-E\).

The cross blocks are not common templates and are not Hadamard blocks.
For the first witness their ranks/determinants are

\[
(3,0),\quad(4,8),\quad(4,8),
\]

while all three blocks of the second witness have rank four and
determinant of absolute value eight.

## 2. Exact chiral form

The first witness has the signed antiautomorphism already recorded in
the source JSON. For the second witness an independent switching
isomorphism search gives

\[
\begin{aligned}
p={}&(3,9,8,0,11,6,5,10,2,1,7,4),\\
\varepsilon={}&(1,1,1,-1,-1,-1,1,-1,-1,-1,1,1).
\end{aligned}
\]

If \(S e_i=\varepsilon_i e_{p(i)}\), exact multiplication gives, for
both matrices,

\[
\boxed{
S^2=-I,\qquad S^\top BS=-B,\qquad SB=-BS.
}
\tag{2.1}
\]

After a signed permutation of coordinates,

\[
S=
\begin{pmatrix}
0&-I_6\\
I_6&0
\end{pmatrix}.
\]

Every real symmetric matrix anticommuting with this \(S\) has the
form

\[
\boxed{
B=
\begin{pmatrix}
A&C\\
C&-A
\end{pmatrix},
}
\tag{2.2}
\]

where \(A,C\) are symmetric. For a zero-diagonal signing \(B\),
\(A\) is a zero-diagonal signing and \(C\) is a full symmetric sign
matrix.

Writing \(z=x+iy\), with \(x,y\in\{\pm1\}^6\), and

\[
M=A-iC,
\]

one has the exact complex representation

\[
\boxed{
\begin{pmatrix}x\\y\end{pmatrix}^{\!\top}
B
\begin{pmatrix}x\\y\end{pmatrix}
=\operatorname{Re}(z^\top Mz).
}
\tag{2.3}
\]

Multiplication \(z\mapsto iz\) negates (2.3), explaining the exact
energy symmetry.

This is a genuine analytic pattern, but not a recursively contracting
one: chiral symmetry centers the energy range, eliminating the only
conditional-mean source of strict compressed-lift improvement.

## 3. The exact energy-range obstruction

Let \(A\) be an order-\(n\) signing and let \(B\) be an exact
compressed \(s\)-lift, meaning that its \(r=s\) cross blocks satisfy

\[
\mathbf 1^\top B_{ij}\mathbf 1=s^{3/2}a_{ij}.
\tag{3.1}
\]

Write

\[
U(A)=\max_xx^\top Ax,\qquad
L(A)=\min_xx^\top Ax,
\]

and define the midpoint and half-width

\[
\mu(A)=\frac{U(A)+L(A)}2,\qquad
w(A)=\frac{U(A)-L(A)}2.
\]

Since the uniform average of \(x^\top Ax\) is zero,

\[
L(A)\le0\le U(A)
\]

and therefore

\[
Q(A)=w(A)+|\mu(A)|.
\tag{3.2}
\]

For a macro spin \(x\), repeat \(x_i\) on all \(s\) vertices of the
\(i\)-th fiber. Its lifted energy is

\[
\boxed{
E_B(x\otimes\mathbf1_s)
=s^{3/2}x^\top Ax+d_B,
}
\tag{3.3}
\]

where

\[
d_B=\sum_{i=1}^n\mathbf1^\top D_i\mathbf1
\]

is the total diagonal-fibre energy. Because each \(D_i\) is a
zero-diagonal signing,

\[
|d_B|\le ns(s-1).
\tag{3.4}
\]

The maximum over the two endpoints of the seed energy interval gives

\[
\begin{aligned}
Q(B)
&\ge
\max\left\{
\left|s^{3/2}U(A)+d_B\right|,
\left|s^{3/2}L(A)+d_B\right|
\right\}\\
&=
s^{3/2}w(A)
+\left|d_B+s^{3/2}\mu(A)\right|.
\end{aligned}
\tag{3.5}
\]

Optimizing the last term subject to (3.4) proves:

\[
\boxed{
Q(B)
\ge
s^{3/2}w(A)
+\left(s^{3/2}|\mu(A)|-ns(s-1)\right)_+.
}
\tag{3.6}
\]

Equivalently,

\[
\boxed{
Q(B)
\ge
s^{3/2}Q(A)
-\min\left\{
s^{3/2}|\mu(A)|,\,
ns(s-1)
\right\}.
}
\tag{3.7}
\]

This improves the cruder conditional-mean estimate by identifying the
only possible gain: shifting the midpoint of an asymmetric energy
interval.

### Fixed-fibre consequence

For fixed \(s\),

\[
\frac{ns(s-1)}{(ns)^{3/2}}
=O\!\left(\sqrt{\frac{s}{n}}\right)\to0.
\]

Thus no fixed-\(s\) exact compressed lift can improve the normalized
seed constant by a fixed amount as \(n\to\infty\).

### Chiral consequence

If \(A\) is switching-isomorphic to \(-A\), then

\[
U(A)=-L(A),\qquad \mu(A)=0.
\]

Equation (3.6) becomes

\[
\boxed{
Q(B)\ge s^{3/2}Q(A)
}
\tag{3.8}
\]

for every exact compressed \(s\)-lift \(B\).

Both order-\(12\) witnesses satisfy this hypothesis. Consequently,
every exact compressed four-lift \(C\) of either witness obeys

\[
\boxed{Q(C)\ge8\cdot40=320.}
\tag{3.9}
\]

The strict \(48\to40\) contraction therefore cannot repeat.

## 4. Why the triangle contraction is finite-size recentering

For either triangle switching class the energy interval has width
eight:

\[
U-L=8,\qquad Q=6,\qquad |\mu|=2.
\]

At \(s=4\), the scaled half-width is

\[
s^{3/2}w=8\cdot4=32.
\]

The diagonal blocks can shift the midpoint by up to

\[
ns(s-1)=3\cdot4\cdot3=36,
\]

so complete recentering is feasible at this tiny order. The universal
constant-fibre lower bound is therefore only \(32\). The observed
value \(40\) consists of this centered baseline plus an eight-unit
microscopic penalty.

At large \(n\) with \(s=4\), the available centering shift is only
\(12n\), whereas a leading energy midpoint may be of order
\(n^{3/2}\). The mechanism responsible for \(48\to40\) is therefore
not scale-preserving.

## 5. Audit of seed-independent uniform gadgets

To test whether the two witnesses came from a hidden common local
operator, every uniform four-fibre lift

\[
\mathcal T_{R,D}(A)=A\otimes R+I_n\otimes D
\tag{5.1}
\]

was enumerated, where

- \(R\) runs over all \(2^{10}\) symmetric full sign matrices of
  order four;
- \(D\) runs over all \(2^6\) zero-diagonal symmetric sign matrices
  of order four;
- the exact-compression condition is
  \(\mathbf1^\top R\mathbf1=8\).

For each pair \((R,D)\), all \(2^{12}\) Boolean states were evaluated
for both triangle classes.

The exact results are:

| optimization target | minimum lifted \(Q\) |
|---|---:|
| first class alone | \(36\) |
| second class alone | \(36\) |
| one common \((R,D)\) for both classes | \(\boxed{48}\) |

Thus each triangle profile can be recentered by a class-dependent
uniform gadget, but no seed-independent common tensor gadget contracts
both. The edge-dependent \(Q=40\) witnesses are not disguised common
block lifts.

## 6. Composition verdict

The exact compressed classes do compose:

\[
B\in\mathcal L_s(A),\quad C\in\mathcal L_t(B)
\quad\Longrightarrow\quad
C\in\mathcal L_{st}(A).
\]

However, the range theorem shows what an iterative proof would now
have to establish:

\[
G_t(B)\le t^{3/2}Q(B)+o((Nt)^{3/2})
\]

for already centered, typically chiral, intermediate profiles.
There is no further midpoint gain available. All nonconstant fiber
states must be controlled at essentially exact equality.

The standard common gadget tested on the first order-\(12\) witness
already gives a heuristic lower value above the required \(320\), so
it does not supply this equality step. The chiral complex form (2.2)
suggests a possible recursively constrained search, but at present
there is no composable analytic rule.

## 7. Conclusion

The compressed-lift discovery survives verification but changes
interpretation:

\[
\boxed{
\text{dependent finite absorption exists, but strict contraction is
entirely a recentering resource and stops on the first chiral output.}
}
\]

It therefore does not prove the desired amplification inequality.
Its remaining value is narrower and concrete: search for a chiral
lift operator attaining normalized **equality** on centered profiles.
Anything weaker cannot compose.

