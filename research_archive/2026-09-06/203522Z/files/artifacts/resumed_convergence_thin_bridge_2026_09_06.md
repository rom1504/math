# Sharp thin bridges and the actual-minimizer response obstruction

Date: 2026-09-06. This is a proved original-problem cross-order estimate,
not a convergence proof. Its lower obstruction applies to iid bridges,
not to arbitrary designed bridges or to the original minimizing sequence.

**Archive reconciliation:** Sections 1--2 independently rediscover the
already proved `mesoscopic_completion_2026_09_05.md`, and are not new
campaign progress. The response form and iid comparison are retained for
their exact scope, without claiming a new research route. Work on this
line was stopped after the root identified the duplication.

Write

`H_A(x)=x^T A x/2`, `Q(A)=max_x |H_A(x)|`, `M_n=min_A Q(A)`.

All child matrices below are hollow symmetric signings. A rectangular
bridge has arbitrary signs and no symmetry condition. Let

`mu_r=E |epsilon_1+...+epsilon_r|`,

where the epsilons are independent uniform signs. Thus

`mu_(2k)=2k binom(2k,k)/2^(2k)` and
`mu_(2k+1)=(2k+1) binom(2k,k)/2^(2k)`.

In particular `mu_r/sqrt(r) -> sqrt(2/pi)`.

## 1. An exact rectangular estimate, sharp in the thin regime

Let `b(n,h)` be the minimum of
`beta(B)=max_(x,y) |x^T B y|` over `n by h` sign matrices. Then

```
n mu_h <= b(n,h)
        <= n mu_h + sqrt(2 n h (h-1) log 2).                 (1)
```

For the lower bound, average
`sum_i |sum_j b_ij y_j|` over a uniform sign vector y. Every row has
expected absolute sum `mu_h`, independently of its pattern. The maximum
over y is at least this average.

For the upper bound take all entries of B independently uniform. For a
fixed y, put `Z_y=sum_i |sum_j b_ij y_j|`. Its mean is `n mu_h`.
Changing one of the `nh` independent signs changes `Z_y` by at most 2.
The conditional-expectation martingale, followed by Hoeffding's elementary
bounded-interval exponential estimate at each step, therefore gives

`E exp(t(Z_y-E Z_y)) <= exp(nh t^2/2)`.

There are `2^(h-1)` distinct queries modulo `y -> -y`. Consequently

`E max_y Z_y <= n mu_h + log(2^(h-1))/t + nh t/2`.

Optimize over positive t. For h=1 there is only one query and
`beta(B)=n`, so (1) also holds at that endpoint. An outcome no larger
than the expectation proves existence.

Thus, for `h -> infinity` and `h/n -> 0`,

```
b(n,h)/(n sqrt(h)) -> sqrt(2/pi).                           (2)
```

This is an exact asymptotic for the minimum rectangular norm in this
regime, not merely a random-matrix estimate.

## 2. Original-problem cross-order consequence

For every n,h at least one,

```
M_(n+h) <= M_n + M_h + n mu_h
                    + sqrt(2 n h (h-1) log 2).             (3)
```

Indeed choose exact child minimizers A,D and a bridge satisfying (1).
The exact bridge reversal identity is

`Q([[A,B],[B^T,D]]) = max_(x,y)(|H_A(x)+H_D(y)|+|x^T B y|)`.

It is bounded above by `M_n+M_h+beta(B)`. No cancellation between
optimizing children and bridge is assumed.

Let `c_n=M_n/n^(3/2)`. For `h -> infinity`, `delta=h/n -> 0`, (3)
and the known all-order boundedness of c_h give

```
c_(n+h) <= c_n/(1+delta)^(3/2)
  + [sqrt(2/pi)+o(1)] sqrt(delta)/(1+delta)^(3/2)
  + sqrt(2 log 2) delta/(1+delta)^(3/2)
  + O(delta^(3/2)).                                       (4)
```

The leading bridge coefficient `sqrt(2/pi)` is optimal for a proof
that separately pays the entire rectangular norm, by (2). Equation (4)
improves the elementary union-over-both-spin-families coefficient, but
retains a square-root modulus and does not imply convergence.

## 3. An energy-aware, exact response refinement

For a fixed child A define its absolute external-field response

`F_A(z)=max_x (|H_A(x)|+z dot x)`.

Global spin reversal gives `F_A(-z)=F_A(z)`. If `S_h^(n)` is a vector
of n independent copies of a sum of h independent signs, then a random
bridge as above gives

```
M_(n+h) <= M_h + E F_A(S_h^(n))
                    + sqrt(2 n h (h-1) log 2)              (5)
```

for every exact order-n minimizing child A. More generally replace M_n
implicitly present in `F_A` by the cap of any selected child.

For each fixed y, `By` has the law `S_h^(n)`. Flipping one entry of B
changes `F_A(By)` by at most 2, since it changes one field coordinate
by 2 and all spins have magnitude 1. The same martingale argument and
the same `2^(h-1)` query quotient bound `E max_y F_A(By)` by the
right side of (5), without its M_h term. Finally, the exact block cap
is at most `M_h+max_y F_A(By)`.

The crude pointwise inequality
`F_A(z)<=Q(A)+||z||_1` recovers (3). Thus (5) retains actual labelled
energy geometry that (3) discards. It does not replace this geometry
by a shell-cardinality or scalar entropy hypothesis.

## 4. Iid bridges miss the required derivative coefficient

There is a complementary lower statement on the same actual children.
Fix arbitrary A,D of orders n,h and let B have iid sign entries.
Choose an absolute ground state x of A and its orientation sigma,
so `sigma H_A(x)=Q(A)`. Maximizing the resulting linear field in y
and paying the most adverse possible internal D energy proves

```
Q([[A,B],[B^T,D]]) >= Q(A)+||B^T x||_1-Q(D).               (6)
```

Here `||B^T x||_1` is a sum of h independent copies of `|S_n|`.
Its mean is `h mu_n` and its variance is at most hn. Hence, whenever
`n -> infinity` and `h -> infinity`, Chebyshev's inequality gives

`||B^T x||_1 = [sqrt(2/pi)+o_P(1)] h sqrt(n)`.

Suppose now that A,D are actual minimizers, `h=o(n)`, and
`c_n<=1/2+o(1)`. The known upper bound makes
`M_h=O(h^(3/2))=o(h sqrt(n))`. Equation (6) becomes

```
Q(parent) >= M_n+[sqrt(2/pi)-o_P(1)]h sqrt(n).              (7)
```

In contrast, preserving the order-n normalized cap at order n+h requires

`(1+h/n)^(3/2) M_n
 = M_n+[3 c_n/2+o(1)]h sqrt(n)`.

The coefficient in this target is at most 3/4+o(1), whereas

`sqrt(2/pi)-3/4 = 0.047884560802865... > 0`.

Therefore iid bridges exceed this local scaling target with probability
tending to one, by a positive multiple of `h sqrt(n)`. This conclusion
already uses a single ground state of the actual child. It neither
requires a large ground-state layer nor appeals to a surrogate cap.

The implication is limited but concrete: an upper-preserving local
construction cannot draw its bridge independently of all child ground
states. Correlating or balancing the bridge against those states is
not covered by (7), and is not ruled out here.

## 5. Scope

Equations (1), (3), and (5) are positive finite-order results. Equation
(2) identifies the best thin-block bridge coefficient under separate
norm payment. Equation (7) excludes iid bridges from the needed local
scaling construction, not arbitrary bridges or original convergence.
No new value of liminf or limsup is claimed.

Only elementary martingale exponential bounds, Chebyshev, and the
central-binomial asymptotic are used in the displayed derivations.
