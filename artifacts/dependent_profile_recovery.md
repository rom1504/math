# Dependent profile recovery: exact compressed-lift checkpoint

## Status

For a seed signing \(A\) of order \(n\) and a square fibre size \(s\),
let \(\mathcal L_s(A)\) be the set of order-\(ns\) signings \(B\),
partitioned into \(n\) fibres of size \(s\), whose cross-block sums obey
\[
\mathbf1^\top B_{ij}\mathbf1
=a_{ij}s^{3/2}\qquad(i\ne j). \tag{1}
\]
The diagonal fibre blocks are arbitrary zero-diagonal signings.  Define
\[
G_s(A)=\min_{B\in\mathcal L_s(A)}Q(B). \tag{2}
\]
This is the most permissive exact-compression version of dependent
microblock recovery: every microblock may depend jointly on the entire
seed.

The first exact finite computation gives a positive result.  For the
all-negative triangle \(C^-\),
\[
Q(C^-)=6,\qquad
\boxed{G_4(C^-)\le40<48=4^{3/2}Q(C^-).} \tag{3}
\]
Thus the independent-Hadamard absorption wall does not extend to
arbitrary jointly dependent microblocks.  The value \(40\) is a
certified upper witness, not a certified optimum; the current mixed
integer lower bound is \(32\).

## 1. Exact witness

The matrix is stored without abbreviation in
`dependent_profile_recovery_witness.json`.  Exhaustive enumeration of
the \(2^{11}\) configurations modulo \(x\sim-x\) proves
\[
Q(B)=40. \tag{4}
\]
Its three \(4\times4\) cross blocks all have sum \(-8\), exactly
\[
-4^{3/2},
\]
so it satisfies (1).  The normalized values are
\[
\frac{Q(B)}{12^{3/2}}
=0.962250448649\ldots,\qquad
\frac{Q(B)}{2\cdot12^{3/2}}
=0.481125224324\ldots . \tag{5}
\]
The latter is the original half-energy normalization and is below
\(1/2\).

This is a genuinely global construction.  The three cross blocks are
different and non-Hadamard: one has rank three, while the other two
have determinant \(8\).  Hence it is not a common tensor gadget in
disguise.

## 2. Exact chiral symmetry

Although the witness has no nontrivial signed automorphism, it has one
signed antiautomorphism.  Let \(S\) be the signed permutation specified
in the JSON file.  Exact multiplication gives
\[
S^2=-I,\qquad S^\top BS=-B,\qquad SB=-BS. \tag{6}
\]
Consequently the spectrum is symmetric about zero.  Numerically it is
\[
\pm\{4.83759104,\ 4.29540995,\ 3.38243991,\ 2.75976426,\
\sqrt5,\ 0.29994629\}.
\]
The exact characteristic polynomial is
\[
\lambda^{12}-66\lambda^{10}+1627\lambda^8-18604\lambda^6
+98663\lambda^4-196850\lambda^2+16925. \tag{7}
\]
This chiral structure was not imposed in the optimization.  It is the
only clear algebraic pattern currently visible, and may be useful for
a recursively constrained search.

## 3. Universal conditional-mean bound

There is a simple theorem valid for every \(B\in\mathcal L_s(A)\).
Fix fibre magnetizations
\[
m_i\in\{-1,-1+2/s,\ldots,1\}
\]
and sample the spins in fibre \(i\) independently and uniformly subject
to that magnetization.  For \(i\ne j\),
\[
\mathbb E[x_i^\top B_{ij}x_j]
=m_im_j\,\mathbf1^\top B_{ij}\mathbf1
=s^{3/2}a_{ij}m_im_j. \tag{8}
\]
The total expected diagonal-fibre contribution has absolute value at
most \(ns(s-1)\).  Since a maximum absolute value dominates the
absolute mean,
\[
\boxed{
Q(B)\ge
s^{3/2}\max_{m\ \mathrm{on\ the\ grid}}|m^\top Am|
-ns(s-1).
}
\tag{9}
\]
In particular, taking a Boolean maximizer of \(A\),
\[
G_s(A)\ge s^{3/2}Q(A)-ns(s-1). \tag{10}
\]
For fixed \(s\) and \(n\to\infty\), (10) says that exact compression
cannot improve the normalized seed value by a fixed amount.

On the other hand,
\[
\max_{m\in[-1,1]^n}|m^\top Am|=Q(A), \tag{11}
\]
because the zero diagonal makes the quadratic form affine in each
coordinate separately, so its maximum and minimum over the cube occur
at vertices.  Therefore conditional averaging reaches precisely the
desired baseline and supplies no positive gap above it.  A universal
impossibility theorem must use fluctuations or joint entropy, not just
conditional means.

The order-12 witness shows why the \(O(ns^2)\) term in (10) cannot be
ignored at small \(n\): dependent diagonal and cross blocks can center
the finite system enough to reduce \(48\) to \(40\).

## 4. Exact composition law

The constrained lift classes compose.  If
\[
B\in\mathcal L_s(A),\qquad C\in\mathcal L_t(B),
\]
then regrouping the \(nst\) vertices into the original \(n\) fibres of
size \(st\) gives
\[
C\in\mathcal L_{st}(A). \tag{12}
\]
Indeed, for an original macro edge \(i,j\), its aggregate block sum is
\[
\sum_{u\in i,\ v\in j}
\mathbf1^\top C_{uv}\mathbf1
=t^{3/2}\sum_{u\in i,\ v\in j}b_{uv}
=(st)^{3/2}a_{ij}. \tag{13}
\]
Consequently,
\[
\boxed{
G_{st}(A)
\le
\min_{B\in\mathcal L_s(A)}G_t(B).
}
\tag{14}
\]
This is a genuine profile-valued composition inequality.  It is not a
closed scalar inequality, because the right side depends on the full
intermediate signing \(B\), not only on \(Q(B)\).

If one could prove the uniform absorption estimate
\[
G_t(B)\le t^{3/2}Q(B)+o((Nt)^{3/2}) \tag{15}
\]
for every near-minimizing order-\(N\) signing \(B\), then (14) would
iterate and provide the amplification theorem needed for convergence.
Equation (15) is exactly the unresolved step.

## 5. What the finite witness does and does not prove

The witness proves:

- joint dependence can beat every fixed common-block tensor rule;
- the independent residual-covariance wall is not universal;
- exact compression admits nontrivial finite absorption;
- constrained lift classes possess an exact composition law.

It does not yet prove:

- that \(40\) is the optimum constrained value;
- a construction for arbitrary seeds;
- that the chiral pattern iterates;
- or the asymptotic estimate (15).

The next concrete experiment is to impose the chiral relation (6) in a
cutting-plane search for a compressed lift of the order-12 witness.
The theoretical target is a conditional-discrepancy theorem controlling
all Boolean fibre configurations at once while retaining the exact
compression constraints.

## Reproducibility

Run
```text
python3 verify_dependent_profile_recovery.py
```
to verify the signing, all block sums, all \(2^{11}\) quadratic
energies, both normalizations, and the signed antiautomorphism.
