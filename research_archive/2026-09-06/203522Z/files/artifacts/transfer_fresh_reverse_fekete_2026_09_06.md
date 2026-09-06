# Lower powered aggregation: a frozen target and exact convergence implication

Date: 2026-09-06. This is a proof target, not a proved original-parameter
recurrence. Write `Q(A)=max_x |x^T A x|/2`, `M_n=min_A Q(A)`,
`b_n=M_n^(2/3)`.

## 1. Target frozen before inspecting this archive direction

The genuinely lower-direction target is: there are absolute C and n_0
such that, whenever m,n>=n_0 and `1/2<=m/n<=2`,

`b_(m+n) >= b_m+b_n-C sqrt(m+n)`.                        (RF)

This is deliberately quantitative. A merely qualitative o(m+n) defect
does not suffice: slow oscillations can obey such a defect. The target
is also for the ACTUAL minima M_m,M_n, not caps of arbitrary restrictions.

The naked cut-code formulation gives the same obligation. If D_n is
`binom(n,2)` and rho_n is the covering radius of the cut code augmented
by its all-one coset, then `M_n=D_n-2rho_n`. Thus (RF) seeks a nonlinear
upper bound on the parent covering radius, with the full bridge and
the common-spin decoding constraint retained. Ordinary direct-sum
covering-radius addition controls a different code; independent child
decodings do not by themselves establish (RF).

## 2. Exact implication of (RF)

Principal restriction gives M_n, and hence b_n, monotone. The known
original scale bound gives `b_n=O(n)`. Under these two facts, (RF)
would prove convergence of `b_n/n` and hence of `M_n/n^(3/2)`.

Fix k>=n_0. Combine L identical blocks of order k using a binary tree
which at every node splits the number s of leaves into
`floor(s/2),ceil(s/2)`. Thus every use of (RF) has child-order ratio
between 1/2 and 2, including the three-leaf split 1+2.

Let T(L) be the sum, over internal nodes, of the square roots of their
numbers of leaves. Then

`T(1)=0`, `T(L)=sqrt(L)+T(floor(L/2))+T(ceil(L/2))`.

The elementary induction bound

`T(L)<=3(L-sqrt(L))`                                    (2)

holds. Indeed, if a+b=L and 1/2<=a/b<=2, then
`sqrt(a)+sqrt(b)>=(4/3)sqrt(L)`, so the induction step reads

`sqrt(L)+3L-3(sqrt(a)+sqrt(b))<=3(L-sqrt(L))`.

Applying (RF) at every node gives

`b_(Lk)>=L b_k-C sqrt(k) T(L)>=L b_k-3C L sqrt(k)`.       (3)

For arbitrary N use L=floor(N/k) and monotonicity to obtain
`b_N>=b_(Lk)`. At fixed k, letting N tend to infinity in (3) yields

`liminf_N b_N/N >= b_k/k-3C/sqrt(k)`.

Finally choose k tending to infinity through a limsup subsequence.
This proves convergence, with no unfilled dyadic phases or unstated
integer allocation requirement. Equal-halves-only transfer would not
justify this argument: the three-leaf split is used essentially.

## 3. First aggregation obstacle

A sufficient pointwise strengthening would have to control a genuine
hollow parent together with its two principal children, but one cannot
simply replace the ACTUAL minima in (RF) by arbitrary child caps.
The known small-retention theorem makes this especially dangerous:
random small principal children of bounded-cap parents can have a
larger normalized cap than the parent. Any proof must exploit the
minimum-order child benchmark, rather than silently treating those
large actual child caps as values of M_m.

No proof of (RF) is claimed. The next task is to test exact lower
aggregation mechanisms or a code-radius inequality with the correct
common-spin constraints, not to rename (RF) as an established lemma.

## 4. A concrete pointwise powered-aggregation falsifier

The unrestricted pointwise version fails at a leading linear b-scale,
even with both child normalized caps arbitrarily close to 1/2 from
ABOVE. This does not falsify (RF), whose child benchmarks are M_n.

Fix delta>0. Take `m=4^a`,
`H_m=(J_4-2I_4)^(tensor a)`, and
`k=floor(sqrt(2 delta) m^(3/4))`. Thus `H_m^2=mI` and
`H_m 1=sqrt(m)1`. Use the hollowed core matrix

`L_m = hollow([[H_m,H_m],[H_m,-H_m]])`.

Its cap is at most `sqrt(2) m^(3/2)+O(m)` by the operator bound.
Its two principal m-vertex blocks have caps
`m^(3/2)/2+O(m)`; the all-one vectors already give the lower estimates.

Adjoin k new vertices to each child. Make the first new k-set a
positive clique and the second a negative clique. Fill every remaining
missing edge with signs so that their combined quadratic cap is
`o(m^(3/2))`. Such a filling exists by independent signs and a union
bound: there are `O(mk+k^2)` missing edges and `O(m)` vertices, so the
cap of the random-edge remainder is `O(m sqrt(k)+k sqrt(m))`, which is
`o(m^(3/2))` for this k. Principal remainders obey the same bound.

Let A_m,D_m be the two children, of equal order `n=m+k`, and P_m the
completed parent. For the first child, the all-one core and clique
states align at the positive orientation. They align at the negative
orientation for the second child. Consequently

`Q(A_m)=Q(D_m)=(1/2+delta+o(1)) m^(3/2)`.

More precisely both equalities here are asymptotic equalities, not an
assertion of equal finite caps. For the parent, the *difference* of
the positive- and negative-clique energies has absolute maximum
`k^2/2+O(k)=delta m^(3/2)+o(m^(3/2))`, not twice that quantity.
Therefore

`Q(P_m)<=(sqrt(2)+delta+o(1)) m^(3/2)`.

It follows that

`Q(A_m)^(2/3)+Q(D_m)^(2/3)-Q(P_m)^(2/3)`
` >= [2(1/2+delta)^(2/3)-(sqrt(2)+delta)^(2/3)-o(1)] m`.  (4)

The displayed coefficient is strictly positive for delta>0: after
raising positive quantities to the power 3/2, this is precisely
`2 sqrt(2)(1/2+delta)>sqrt(2)+delta`.

Thus no O(sqrt(n)) error repairs universal pointwise powered
superadditivity in the bounded-cap class. The child constants converge
to `1/2+delta`; the parent constant is at most
`1/2+delta/(2 sqrt(2))`. Taking delta arbitrarily small does NOT move
the example into the original sub-1/2 minimizing class.

This construction exploits the difference between two-sided cap and
centered range: the one-sided coherent additions align inside each
child but pay only one copy at the parent. It is compatible with the
exact global-reversal block identity, and assumes no cancellation of
the random bridge with its children.
