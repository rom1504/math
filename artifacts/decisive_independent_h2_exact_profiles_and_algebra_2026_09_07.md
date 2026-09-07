# H2 lift: exact finite profiles and the clique-flip identity

Status: exact finite computations plus elementary identities. No asymptotic
selected-minimizer estimate is established.

Define `L(A)=[[A,A+I],[A+I,-A]]`. This is a hollow signing of order 2n.
The following table gives exact exhaustive caps; the parent rows are exact
minimizers at their stated orders, with the stored proof/certificate scopes
recorded in the earlier nesting replay.

| Parent order | Parent cap | Selected parent | Lift cap |
|---:|---:|---|---:|
| 5 | 4 | gauged code 13 | 13 |
| 6 | 5 | gauged code 220 | 18 |
| 7 | 9 | gauged code 828 | 31 |
| 7 | 9 | gauged code 826 | 25 |
| 8 | 10 | gauged code 53014 | 40 |
| 9 | 12 | gauged codes 898008, 6737136 | 41 |
| 10 | 13 | exact_m10 | 50 |
| 10 | 13 | first10 of nested_10_in_11 | 48 |
| 11 | 17 | nested_10_in_11 | 49 |
| 11 | 17 | heuristic_m11, subsequently certified optimal | 53 |
| 12 | 18 | extension_nested_m11_to_12 | 56 |

Repeated lifts of the n5 and n6 matrices give respectively `4->13->48`
and `5->18->56`. The first lifts happen to be exact minimizers at their
new orders. This finite observation is not a hereditary asymptotic theorem.

Exact source:
`computations/decisive_independent_h2_lift_caps_2026_09_07.cpp` and `.py`;
stored matrices, maximizing spin codes, support energies, and results in
the matching `.json`. The C++ scan fixes one spin and visits every remaining
spin by Gray code with exact integer field/energy updates. The driver checks
the maximizing witness by an independent quadratic identity.

## 1. Cross-support reduction

Given lift spins x,y, put `p=(x+y)/2`, `r=(x-y)/2`. Their supports partition
the old vertices. Write `u=q_A(p)`, `v=q_A(r)`, `w=p^TAr`. Then

```math
q_{L(A)}(x,y)=2(u-v+w)+x\cdot y.                       (1)
```

In particular, the n8 witness has support sizes 4,4 and
`(u,v,w)=(6,-6,8)`, yielding cap40 with zero diagonal correction.
The bad stored n10 witness has sizes5,5 and `(10,-10,5)`, yielding50.
These witnesses exploit oppositely coherent principal blocks. Such blocks
at a fixed positive fraction of arbitrarily large orders would have cap
of order n² and hence are not by themselves asymptotic-minimizer examples.

Product rounding on the original cube gives the additional exact constraint

```math
|u a^2+v b^2+wab|\le Q(A)\qquad(a,b\in[-1,1]).          (2)
```

These constraints alone have only the sharp bound

```math
|u-v+w|\le(1+\sqrt2)Q(A).                              (3)
```

Indeed let `t=sqrt2-1`. Subtract the two values at `(a,b)=(1,t)` and
`(t,-1)`; the difference is `2t(u-v+w)`, so its absolute value is at most
2Q. The coherent two-block gadget in the accompanying planting artifact
asymptotically attains this ratio. Thus (2) alone cannot yield the desired
`sqrt2 Q+O(n)` control of (1).

## 2. Exact clique-flip orbit interpretation

Let T be the support of r and z=x. Let A^T denote the signing obtained by
reversing exactly the edges whose BOTH endpoints lie in T. Then

```math
q_{L(A)}(x,y)=2q_{A^T}(z)+n-2|T|.                      (4)
```

Consequently

```math
\left|Q(L(A))-2\max_{T\subseteq[n]}Q(A^T)\right|\le n.  (5)
```

This relates the lift to the worst member of a clique-flip orbit. If A is
an exact minimizer, each orbit member has cap at least Q(A); that lower
statement does not upper-bound the worst orbit member. A proposed sharp
H2 theorem would therefore need genuine uniform control of this orbit for
appropriately SELECTED exact minimizers.

## 3. A mixed row identity (not a deterministic row selection)

For every original spin x let `d=Q(A)-|q_A(x)|`, and choose s so
`s q_A(x)=|q_A(x)|`. Comparing x to each one-coordinate flip gives
`s x_i(Ax)_i>=-d/2`. Since their sum is `2|q_A(x)|`,

```math
\|Ax\|_1\le2Q(A)+(n-2)d.                              (6)
```

If the diagonal hole of each row is filled by a sign, the uniform
distribution over these n rows consequently satisfies, simultaneously
for all x,

```math
\mathbb E_i[|a_i\cdot x|-d]\le2Q(A)/n+1-2d/n.          (7)
```

Equation (7) is a MIXED-row certificate. Interchanging its expectation
and maximum, or purifying it to one signing row at comparable cost, is
not justified. The coefficient 2Q/n also differs from the sharp drift
coefficient 3Q/(2n) required by the one-row convergence architecture.

## 4. Convergence scope

An H2 recurrence for selected exact minimizers with summable normalized
errors would control dyadic rays, but does not alone prove full convergence.
A corresponding recurrence at a multiplicatively independent scale such
as12 would supply a complete route: reoptimize at each stage, sum the
geometrically decaying errors, and use density of the semigroup generated
by2 and12 plus principal restriction. None of these sharp recurrences is
proved by the identities or finite scans here.
