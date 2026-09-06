# Exact gap-two geometry of minimizing signings

Date: 2026-09-06. Status: director's proposed theorem independently
reconstructed; all steps pass. This gives a spatial separation inside
the smallest nontrivial energy window. It does not prove a Gaussian
width bound, an insertion theorem, or original convergence.

Let A be a hollow symmetric signing of order n. Define

    H_A(x)=sum_(i<j) a_ij x_i x_j,
    Q=max_x |H_A(x)|.

Assume A is an exact minimizer, or more generally that no single
coefficient flip can decrease Q. An oriented state is (sigma,x), with
energy sigma H_A(x); its gap is Q-sigma H_A(x). Spins are identified
under global negation. The projective Hamming distance between spins
is min(d_H(x,y),d_H(x,-y)); the orientation sigma is not counted in it.

## 1. Exact finite theorem

Fix ANY oriented ground state. For every integer 1<=r<=n/2 satisfying

    n(2r-1)<Q,
    Q>r(n-r)+1,                                      (1)

there exists an oriented state with gap at most 2 whose projective
Hamming distance from that ground state is greater than r.

The quantifier is from EVERY ground state, not only from one selected
ground. No normalized operator bound or random-matrix assumption is used.

## 2. Switch to a positive all-one ground

Switch vertices by the chosen ground spin and multiply the whole matrix
by its orientation. These operations preserve the cap and single-edge
local minimality. Call the resulting signing B. Then

    H_B(1)=Q, ell_i=sum_(j!=i) B_ij>=0,
    sum_i ell_i=2Q.

Nonnegativity follows because flipping one spin changes its positive
ground energy from Q to Q-2ell_i, which cannot exceed Q.

Suppose for contradiction that every gap-at-most-two oriented state is
within projective radius r. Choose its representative by flipping a set
S from 1, with |S|<=r. Write cut_B(S)=sum_(i in S,j notin S)B_ij. Then

    H_B(1^S)=Q-2 cut_B(S),
    |cut_B(S)|<=|S|(n-|S|)<=r(n-r).

If this state has NEGATIVE orientation and gap at most 2, then
H_B(1^S)<=-Q+2. The two displayed bounds would imply
Q<=r(n-r)+1, contradicting (1). Thus every such near-active state has
positive orientation.

## 3. Large-field vertices are fixed throughout the near-active window

Write the gap of a positive-oriented near-active state as d, with
0<=d<=2. At that state every signed local field ell_i^S satisfies

    Q-d-2ell_i^S<=Q, hence ell_i^S>=-d/2>=-1.

This uses the cap at the spin obtained by flipping i from the near-state,
not an assumption that the near-state is itself a local maximum. For
each flipped vertex i in S, the exact field transformation is

    ell_i^S=-ell_i+2 sum_(j in S excluding i)B_ij.

Consequently

    ell_i<=2(|S|-1)+1<=2r-1=:h.                       (2)

For comparison, summing cut energies alone only gives the weaker bound

    gap=2 cut_B(S)
       =2 sum_(i in S)ell_i-4 sum_(i<j in S)B_ij<=2.

The pointwise estimate (2) shows directly that no vertex in

    T={i:ell_i>h}

can belong to S. Every oriented near-active state therefore has both
orientation + and spin +1 on all of T.

There must be a POSITIVE coefficient on some edge inside T. Indeed,

    sum_(i in T)ell_i>=2Q-nh>Q.

If every edge inside T were negative, then

    sum_(i in T)ell_i
      =2 sum_(i<j in T)B_ij+cut_B(T)<=cut_B(T)<=Q.

The last inequality is exact: flipping T gives energy
Q-2cut_B(T), and the absolute cap forces it to be at least -Q.
This contradicts the preceding strict inequality.

## 4. A single coefficient flip contradicts minimality

Choose a positive edge inside T and flip its coefficient to negative.
Every gap-at-most-two oriented state's energy decreases by exactly 2,
because its orientation and both endpoint spins equal +1. Every other
oriented state's old gap was greater than 2, and changing one coefficient
can increase an oriented energy by at most 2. It too remains below Q.

The new signing therefore has cap strictly less than Q, contradicting
single-edge local minimality. This proves the theorem.

All oriented energies have the same integer parity; consequently the
new cap is actually at most Q-2. The strict-gap argument already suffices
and does not need this final integrality observation.

## 5. Asymptotic consequence and limits

Suppose exact minimizers satisfy Q>=c n^(3/2) for a fixed c>0. For every
epsilon>0 choose

    r=floor((c/2-epsilon)sqrt(n)),

where 0<epsilon<c/2.

Both conditions (1) hold for all sufficiently large n: the first has
a strict leading-coefficient margin, and the leading coefficient of
r(n-r) is c/2-epsilon<c. Thus, from every ground state, another oriented state in the
GAP-TWO window has projective distance at least

    (c/2-o(1))sqrt(n).                                 (3)

The currently audited bound permits c=.4333221116640807, giving the
separation coefficient .21666105583204035. Statement (3)
uses only the positive universal lower constant; it does not require
the upper coefficient 1/2.

The conclusion is stronger spatial information than merely counting
near-active states. It is still sublinear separation and only asserts
one distant state from each ground. It does not give an extensive
packing, an overlap distribution, or the discrepancy control required
for a fixed-block insertion estimate. It also does not automatically
extend to signings with an o(n^(3/2)) optimality error: their cap could
tolerate the two-unit improvement used in the contradiction.
