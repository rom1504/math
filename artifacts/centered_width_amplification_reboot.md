# Centered-width amplification reboot

## Status

For a symmetric zero-diagonal signing \(A\), use the one-copy
Hamiltonian
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j
\]
and define
\[
P(A)=\max_xH_A(x),\qquad
m(A)=\min_xH_A(x),\qquad
W(A)=\frac{P(A)-m(A)}2.
\tag{1}
\]
Thus \(W(A)\) is the centered half-range.  If
\[
G_n=\min_AW(A),\qquad
F_n=\min_A\max_x|H_A(x)|,
\]
then \(G_n\le F_n\).

This note audits two possible steps toward convergence:

1. scale-preserving amplification of \(G_n\);
2. asymptotic midpoint balancing \(F_n-G_n=o(n^{3/2})\).

The first new result is negative but exact: compressed lifting is
monotone in normalized centered width, and equality already fails for
the four-lift of the all-negative triangle.

## 1. Exact cut identity

If \(x,y\) are Boolean vectors and
\[
S=\{i:x_i\ne y_i\},
\]
then
\[
H_A(x)-H_A(y)
=2\sum_{\substack{i\in S\\j\notin S}}a_{ij}x_ix_j.
\]
Consequently
\[
\boxed{
W(A)=
\max_{S\subset[n]}\|A_{S,S^c}\|_{\infty\to1}.
}
\tag{2}
\]
This is the minimum signed hereditary rectangular cut norm, rather
than the ordinary same-sign quadratic norm.

## 2. Exact compressed lifts can never contract centered width

For a square fibre size \(s\), let \(\mathcal L_s(A)\) be the class of
order-\(ns\) signings \(B\), partitioned into \(n\) fibres of size
\(s\), satisfying
\[
\mathbf1^\top B_{ij}\mathbf1
=a_{ij}s^{3/2}\qquad(i\ne j).
\tag{3}
\]
The diagonal fibre blocks are arbitrary signings.  Put
\[
\Gamma_s(A)=\min_{B\in\mathcal L_s(A)}W(B).
\tag{4}
\]

For a constant-fibre spin \(X=x\otimes\mathbf1_s\),
\[
H_B(X)=s^{3/2}H_A(x)+d_B,
\tag{5}
\]
where
\[
d_B=\sum_i\sum_{u<v\text{ in fibre }i}b_{uv}
\]
does not depend on \(x\).  Taking the range over \(x\) proves
\[
\boxed{
\Gamma_s(A)\ge s^{3/2}W(A).
}
\tag{6}
\]
Equivalently,
\[
\frac{W(B)}{(ns)^{3/2}}
\ge
\frac{W(A)}{n^{3/2}}
\qquad(B\in\mathcal L_s(A)).
\tag{7}
\]

This is stronger than the corresponding absolute-norm range theorem:
finite recentering disappears completely because translations cancel
from a range.  Exact compressed lifting can only preserve or increase
the normalized centered width.

The lift classes compose:
\[
B\in\mathcal L_s(A),\quad
C\in\mathcal L_t(B)
\Longrightarrow
C\in\mathcal L_{st}(A).
\tag{8}
\]
Thus an amplification proof in this framework requires asymptotic
equality in (6); no strict contraction can be accumulated at
intermediate scales.

## 3. General equality constraints

Suppose \(B\in\mathcal L_s(A)\) attains equality in (6).  Then the
constant-fibre states in (5) already span the full energy interval of
\(B\).  Hence
\[
\boxed{
\begin{aligned}
P(B)&=s^{3/2}P(A)+d_B,\\
m(B)&=s^{3/2}m(A)+d_B,
\end{aligned}}
\tag{9}
\]
and every microstate \(Y\in\{\pm1\}^{ns}\) satisfies
\[
s^{3/2}m(A)+d_B
\le H_B(Y)\le
s^{3/2}P(A)+d_B.
\tag{10}
\]

Every repeated seed maximizer is therefore a global maximizer of
\(B\), and every repeated seed minimizer is a global minimizer.
After switching by any such repeated state \(X\), all cut sums obey
\[
\sum_{\substack{u\in T\\v\notin T}}b_{uv}X_uX_v
\begin{cases}
\ge0,&X\text{ is a maximum},\\
\le0,&X\text{ is a minimum},
\end{cases}
\qquad T\subset[ns].
\tag{11}
\]
In particular, because every microvertex has odd degree when \(ns\)
is even,
\[
X_u(BX)_u
\begin{cases}
\ge1,&X\text{ is a maximum},\\
\le-1,&X\text{ is a minimum}.
\end{cases}
\tag{12}
\]
Equations (10)--(12) are exact necessary conditions for equality,
not relaxations.

## 4. Equality four-lift of the negative triangle is impossible

Let \(C^-\) be the all-negative triangle.  Its one-copy energy
interval is
\[
[m(C^-),P(C^-)]=[-3,1],\qquad W(C^-)=2.
\tag{13}
\]
For \(s=4\), (6) gives the equality target
\[
\Gamma_4(C^-)\ge4^{3/2}W(C^-)=16.
\tag{14}
\]

Assume equality.  Write, for a microvertex \(u\) in fibre \(i\),
\[
d_u=\text{its internal-fibre row sum},\qquad
r_{uj}=\text{its row sum into fibre }j.
\]
The repeated all-equal macro spin is a global minimum.  The three
macro spins with one exceptional sign are global maxima.  Applying
(12) to a vertex in fibre \(1\) gives
\[
\begin{aligned}
d_u+r_{u2}+r_{u3}&\le-1,\\
d_u-r_{u2}-r_{u3}&\ge1,\\
d_u-r_{u2}+r_{u3}&\ge1,\\
d_u+r_{u2}-r_{u3}&\ge1.
\end{aligned}
\tag{15}
\]
Here
\[
d_u\in\{-3,-1,1,3\},\qquad
r_{uj}\in\{-4,-2,0,2,4\}.
\]
The finite integer system (15) forces
\[
r_{u2},r_{u3}\in\{-4,-2\}.
\tag{16}
\]
But every cross block has total sum \(-8\).  Summing its four row
sums in (16) therefore forces every row sum to equal \(-2\).
Applying the same argument from the other fibre forces every column
sum to equal \(-2\).  Consequently every cross block is
\[
\boxed{-J_4+2P}
\tag{17}
\]
for a permutation matrix \(P\).

Substitution in (15) also gives \(d_u\in\{1,3\}\).  Hence the negative
edges in every diagonal fibre form a matching.  Independent
permutations of the three fibres reduce the three cross-block
permutations to
\[
I,\quad I,\quad P.
\tag{18}
\]
There are \(24\) choices for \(P\) and \(10\) matchings, including the
empty matching, in each diagonal fibre: exactly
\[
24\cdot10^3=24{,}000
\tag{19}
\]
reduced equality candidates.

Exact enumeration of all \(2^{11}\) antipodal spin states for every
candidate gives
\[
\boxed{
\min_{\text{24,000 equality candidates}}W(B)=20>16.
}
\tag{20}
\]
The minimum is attained, for example, by \(P=I\) and three empty
matchings; its energy interval is
\[
[-14,26].
\tag{21}
\]
Thus:
\[
\boxed{
\text{No exact compressed four-lift of }C^-
\text{ attains the centered-width range bound.}
}
\tag{22}
\]

Run
```text
python3 verify_centered_width_equality_no_go.py
```
for the pure-integer enumeration certificate.

## 5. Consequence for an amplification proof

The order-\(12\), \(Q=40\) dependent lift found previously has
one-copy interval \([-20,20]\), so its centered width is \(20\).
It contracts the *absolute* normalized seed value by finite
recentering, but (6) shows why it cannot contract centered width.

The no-go (22) does not exclude
\[
\Gamma_s(A)
\le s^{3/2}W(A)+o((ns)^{3/2})
\tag{23}
\]
for large near-minimizing seeds \(A\).  It does show that (23) cannot
come from a universal exact-equality gadget, even on the smallest
nontrivial seed.  A surviving centered-width proof must do at least
one of the following:

1. prove near-equality only for an asymptotically extremal class whose
   complete endpoint cut profiles satisfy additional structure;
2. relax exact block compression and control the quotient error
   together with the compulsory microscopic action;
3. use a nonlocal multiscale recovery in which intermediate objects
   are profiles rather than signings and only the final scale is
   Boolean-realized.

The full all-level action profile is closed under composition, but its
\(d\)-th iterate reads tuple level \(s^d\); compactness alone does not
give the uniform absorption required in (23).

## 6. Midpoint-balancing target

Encode the negative edges of \(A\) by a word \(a\), let \(C_n\) be the
cut code, and put
\[
f(a)=d(a,C_n),\qquad
g(a)=d(a,\mathbf1+C_n),\qquad
\rho=\rho(C_n\cup(\mathbf1+C_n)).
\]
Then
\[
F_n=\binom n2-2\rho,\qquad
G_n=\binom n2-\max_a(f(a)+g(a)),
\tag{24}
\]
so
\[
\boxed{
F_n-G_n
=\max_a(f(a)+g(a))-2\rho.
}
\tag{25}
\]
Therefore midpoint balancing is exactly the asymptotic normality
inequality
\[
\max_a(f(a)+g(a))
\le2\rho+o(n^{3/2}).
\tag{26}
\]
The stronger conjectured finite bound has \(+1\) in place of the
\(o(n^{3/2})\) term and is verified through order ten.  A
deepest-hole plateau-connectivity proof is impossible: at several
small orders every deepest hole is isolated even when perfectly
balanced.  The remaining midpoint work must therefore be a metric or
counting inequality in the edge-plus-star Cayley graph, not a path
argument.

## 7. Current verdict

The two desired convergence lemmas remain logically sufficient:

\[
\begin{array}{ll}
\text{amplification:}&
G_{sn}\le s^{3/2}G_n+o((sn)^{3/2}),\\[2mm]
\text{midpoint:}&
F_n-G_n=o(n^{3/2}).
\end{array}
\]

This audit proves that exact compressed lifting supplies the opposite
inequality automatically and that equality has already failed at the
first four-lift test.  The viable amplification target is therefore
asymptotic profile recovery for large near-minimizers, not an exact
finite block gadget.  The midpoint target remains (26).

Related sidecars:

- `centered_width_rectangular_system.md`
- `midpoint_centering.md`
- `augmented_cut_code_midpoint_audit.md`
- `dependent_profile_recovery.md`
- `profile_renormalization_semigroup.md`
- `pythagorean_centered_width_block.md`
