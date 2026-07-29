# Midpoint balancing for the augmented cut code of \(K_n\)

## 1. Exact rooted graph model

Let \(N=\binom n2\), let \(C_n\) be the cut code of \(K_n\), and put
\[
D_n=C_n\cup({\bf1}+C_n).
\]
For an edge word \(a\), write
\[
f(a)=d(a,C_n),\qquad g(a)=d(a,{\bf1}+C_n),\qquad
\rho=\rho(D_n).
\]

Gauge-fix vertex \(n\): switch so that all \(n-1\) edges incident with it
are positive.  A switching class is then represented uniquely by a graph
\(G\) on \(m=n-1\) vertices.  If \(\delta S\) denotes the cut of \(S\) in
\(K_m\), then
\[
\boxed{
 f(G)=\min_{S\subseteq[m]}
 \left(|E(G)\mathbin\triangle\delta S|+|S|\right),\qquad
 g(G)=f(\overline G).
}
\tag{1}
\]
The term \(|S|\) counts the root edges made negative by the switch.

Equivalently, in the quotient
\[
\mathbb F_2^{E(K_m)}
\]
the Cayley generators are all single edges and all vertex stars.  The
antipodal center is the all-one edge word.  This is the extra structure
which an arbitrary binary code does not have.

If
\[
P(A)=\max_x H_A(x),\qquad Q(A)=-\min_xH_A(x),
\]
then
\[
\boxed{P=N-2f,\qquad Q=N-2g.}
\tag{2}
\]
Thus
\[
|f-g|\le1\quad\Longleftrightarrow\quad |P-Q|\le2.
\tag{3}
\]

## 2. A sufficient normal-coordinate inequality

The desired theorem would follow at once from
\[
\boxed{
\max_a\{f(a)+g(a)\}\le2\rho(D_n)+1.
}
\tag{4}
\]
Indeed, if \(a\) is any deepest hole, then
\(\min(f(a),g(a))=\rho\).  Inequality (4) forces the other distance to
be at most \(\rho+1\).

There is an exact connection with Graham--Sloane normality.  Form the
length-\(N+1\) graph code
\[
\widetilde D_n=
\{(c,0):c\in C_n\}\cup
\{({\bf1}+c,1):c\in C_n\}.
\]
For its added coordinate, the coordinate norm is
\[
\boxed{
N^{(\mathrm{new})}(\widetilde D_n)
=1+\max_a(f(a)+g(a)).
}
\tag{5}
\]
So (4) is precisely the sharp acceptability-type estimate needed at this
distinguished coordinate.  General normal-code theory does not prove it:
the covering radius of \(\widetilde D_n\) can be \(\rho(D_n)+1\), losing
exactly the two units that matter here.

In energy language, (4) is equivalent to
\[
\frac{P(A)+Q(A)}2\ge F(n)-1
\qquad\hbox{for every signing }A.
\tag{6}
\]
It therefore also says that the minimum centered width differs from
\(F(n)\) by at most one.

## 3. New exact cut-code computation

The program `_cut_midpoint_audit.cpp` enumerates one gauge-fixed
representative of every switching class, evaluates all spin states, and
computes both endpoints.  Its \(n=8\) calculation covers all
\[
2^{\binom72}=2,097,152
\]
switching classes.

| \(n\) | \(F(n)\) | number of optimal classes | optimal endpoint profiles |
|---:|---:|---:|---|
| 3 | 3 | 2 | \((3,1),(1,3)\) |
| 4 | 4 | 6 | \((4,4)\) |
| 5 | 4 | 12 | \((4,4)\) |
| 6 | 5 | 12 | \((5,5)\) |
| 7 | 9 | 3,240 | \((9,7),(7,9)\) |
| 8 | 10 | 4,200 | **all \((10,10)\)** |

Thus there is no cut-code counterexample at \(n=8\); in fact every
optimal switching class is exactly centered.  Independently found
known-optimum witnesses at \(n=9\) and \(n=10\) have exact profiles
\[
(P,Q)=(12,12),\qquad (13,13),
\]
respectively.  These last two statements verify witnesses, while their
global optimality uses the campaign's independent certificates for
\(F(9)\) and \(F(10)\).

There is also an exact, independently formulated centered-width
certificate at \(n=9\).  The script `_cut_width_milp.py` gauge-fixes the
root, uses one binary variable for each of the remaining edge signs, and
imposes the upper and lower endpoint inequalities for all \(2^{n-1}\)
spin states.  HiGHS closed the branch-and-bound tree with zero gap and
proved
\[
\boxed{\min_A(P(A)+Q(A))=24,\qquad G_9=12=F(9).}
\tag{7}
\]
Equivalently,
\[
\max_a(f(a)+g(a))=36-12=24=2\rho(D_9),
\]
so (4) holds with one unit to spare at \(n=9\).  The same MILP
formulation reproduces the exhaustive \(G_8=10\) result.

The same model also closed at \(n=10\) (51,759 branch-and-bound nodes):
\[
\boxed{\min_A(P(A)+Q(A))=26,\qquad G_{10}=13=F(10).}
\tag{8}
\]
Thus the centered-width and absolute objectives agree exactly through
order ten except for the unavoidable one-unit lattice loss at orders
three and seven.

The same exhaustive computation gives the minimum half-range:

| \(n\) | \(\min_A(P+Q)/2\) | \(F(n)\) |
|---:|---:|---:|
| 3 | 2 | 3 |
| 4 | 4 | 4 |
| 5 | 4 | 4 |
| 6 | 5 | 5 |
| 7 | 8 | 9 |
| 8 | 10 | 10 |

Hence the sharp possible one-unit loss in (6) already occurs at orders
3 and 7.

## 4. A plateau-connectivity proof is impossible

A tempting proof is to connect a deepest hole to its antipode while
remaining in the deepest-hole layer, and use the discrete intermediate
value theorem on \(f-g\).  The exact quotient computation rules this out.

Using the full quotient Cayley generating set (all single edges and all
rooted stars), the deepest-hole plateau has the following connectivity:

| \(n\) | components | largest component |
|---:|---:|---:|
| 5 | 12 | 1 |
| 6 | 12 | 1 |
| 8 | 4,200 | 1 |

Thus every deepest hole is isolated at these orders, despite being
perfectly centered.  Balance must come from a metric/counting inequality,
not from a path inside the optimal plateau.

## 5. Star-flip local structure

If a counterexample has
\[
f(a)=\rho,\qquad g(a)\ge\rho+2,
\]
then every coordinate neighbor has distance at most \(\rho\) from
\(C_n\), because its distance to the complement half remains at least
\(\rho+1\).  In the rooted model this means:

* toggling any internal edge cannot increase \(f(G)\);
* toggling any full vertex star cannot increase \(f(G)\).

Using (1), a coset leader \(E(G)\triangle\delta S\) of size \(\rho-|S|\)
must consequently possess near-leaders covering every edge and every
rooted star direction.  This is substantially stronger than local
maximality for an arbitrary code, but no counting contradiction with
\(g\ge\rho+2\) has yet been completed.

## 6. Exact vertex-deletion recurrence and its loss

The rooted formula has a useful exact recursion, but it also exposes why
naive induction misses the sharp result.  Delete a reduced vertex \(v\),
write \(H=G-v\), and let \(r\in\mathbb F_2^{m-1}\) be its incidence row.
For \(T\subseteq V(H)\), put
\[
k_T=|r\mathbin\triangle{\bf1}_T|.
\]
Comparing switches \(T\) and \(T\cup\{v\}\) gives
\[
\boxed{
f_m(G)=
\min_T\left\{
|E(H)\mathbin\triangle\delta T|+|T|
+\min(k_T,m-k_T)
\right\}.
}
\tag{9}
\]
For the complement, the last penalty is
\[
\min(m-1-k_T,k_T+1).
\tag{10}
\]
In particular,
\[
\max_G(f_m(G)+f_m(\overline G))
\le
\max_H(f_{m-1}(H)+f_{m-1}(\overline H))
+2\lfloor m/2\rfloor.
\tag{11}
\]
The added term is sharp for the pointwise extension penalty.  It is too
large for a proof of (4): the corresponding increment of \(2\rho(D_n)\)
is smaller by order \(\sqrt n\).  Averaging over deleted vertices does
not remove this loss without joint information about the two optimizing
switch sets.  Thus a successful induction must couple those switch sets;
separate extension bounds cannot work.

## 7. Why the finite MILP did not yield a symbolic dual

The continuous root relaxation is far from the integer answer: at
\(n=9\) its width bound is \(16\), versus the exact integer value \(24\);
at \(n=10\) it starts at \(18\), versus \(26\).  The final certificates
use extensive branching, conflict cuts, and restarts.  There is no small
LP dual supported on orbit-averaged spin constraints to read off.

This is structurally unsurprising: after root gauge fixing, the internal
edge signs are completely free Boolean variables.  Triangle relations
describe the quotient presentation, but impose no linear constraints on
these gauge variables.  Any symbolic certificate must exploit Boolean
integrality (or a genuinely nonlinear triangle/star inequality), not
just an averaged linear combination of the energy slabs.

## 8. Current verdict

The graph-specific target has been reduced to the concrete inequality
(4), and all cut-code data through \(n=10\) support it.  The exhaustive
\(n=8\) result is new evidence strong enough to kill both a small
counterexample and the plateau-connectivity proof strategy.

The remaining proof problem is:

> Show, using the edge-plus-star metric (1), that
> \[
> f(G)+f(\overline G)\le2\rho(D_n)+1
> \]
> for every graph \(G\) on \(n-1\) vertices.

This is a distinguished-coordinate normality statement for a very
structured evaluation code; it is not presently implied by the general
theory of normal binary codes.
