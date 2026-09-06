# Wave 44: hard-column switching Euler and the incidence-mass cliff

## Status and verdict

**Verified.**  Specializing the old cube Euler identities
(10.991)--(10.992) to the *hard complement-incidence column* gives an exact
coupling between row descent and changes of incidence mass.  At a minimizer of
the exact scalar column objective, it yields a clean conditional theorem: a
geometric-mean switching Harnack bound at cost `G` forces the parent row to be
within `G/lambda` of the uniform baseline.

**Falsified as a generic lemma.**  Exact enumeration supplies an actual
order-nine minimizer column above the baseline row for which every
row-descending one-vertex switch has zero incidence mass.  Hence exact
minimality does not imply even one-step incidence retention, and the natural
Euler summation cannot close without a genuinely new orbit-expansion/Harnack
property.

**Open.**  No project-scale column and no scalable actual-minimizer
obstruction is proved.  The finite cliff isolates the missing property; it is
not an asymptotic falsifier.  Conflict, agreement, decoding, and shared
priorities play no role here.

## 1. Hard complement columns

Fix an exact minimizer `A`, put `Q(A)=q`, and let `d=(sigma,x)` be a full
oriented parent cut.  Gauge its unordered edges and degrees by

```math
s_{ij}=\sigma a_{ij}x_ix_j,qquad
r_i=\sum_{j\ne i}s_{ij},qquad
E=\sum_i r_i=\langle A,d\rangle,qquad
R=\sum_i r_i^2=R_2(d).
\tag{R44.1}
```

For an `m`-selector `S`, write

```math
c_S(d)=2\sum_{\{i,j\}\subset S}s_{ij},qquad
L_S(d)=2c_S(d)-E.
```

The hard complement column and its surprise are

```math
\mathcal I_d=\{S:L_S(d)\ge q\},qquad
\alpha_d=U_m(\mathcal I_d),qquad
h_d=-\log\alpha_d,
\tag{R44.2}
```

with `h_d=+infinity` for an empty column.  By (10.1045), every member of
this column is already favorable.

## 2. Exact block switching specialized to incidence

For `U subseteq [n]`, let `d^U` be obtained by flipping the physical spin on
`U`.  Define the parent and internal crossing weights

```math
C_d(U)=\sum_{i\in U,j\notin U}s_{ij},qquad
C_d^S(U)=\sum_{i\in U\cap S,j\in S\setminus U}s_{ij}.
```

Direct expansion gives

```math
\boxed{
L_S(d^U)-L_S(d)=4C_d(U)-8C_d^S(U).}
\tag{R44.3}
```

This is the hard-threshold complement-incidence specialization of the payoff
Euler equation (10.991).  It gives the exact column-mass transport identity

```math
\boxed{
\alpha_{d^U}
=\mathbb E_{S\sim U_m}
\mathbf 1\{L_S(d)+4C_d(U)-8C_d^S(U)\ge q\}.}
\tag{R44.4}
```

For one vertex, put

```math
a_i^S=\mathbf1_{\{i\in S\}}\sum_{j\in S\setminus\{i\}}s_{ij},qquad
t_i=\sum_{j\ne i}s_{ij}r_j=x_i(A^2x)_i.
```

Then (R44.3) and the old row Euler identity become

```math
\boxed{
L_S(d^i)-L_S(d)=4r_i-8a_i^S,qquad
R_2(d^i)-R_2(d)=4[(n-1)-t_i].}
\tag{R44.5}
```

Since `sum_i t_i=R`, summing the row increments recovers (10.992):

```math
\boxed{
\sum_{i=1}^n [R_2(d^i)-R_2(d)]
=4[n(n-1)-R].}
\tag{R44.6}
```

There is also a useful fixed-block average.  If `U` is uniform among the
`k`-sets and

```math
\theta_k=1-\frac{4k(n-k)}{n(n-1)},
```

then pair counting in the row quadratic form gives

```math
\boxed{
\mathbb E_{|U|=k}R_2(d^U)
=n(n-1)+\theta_k[R-n(n-1)].}
\tag{R44.7}
```

Equations (R44.3)--(R44.7) were checked coefficient by coefficient on every
oriented cut of `A_5,A_6,A_8,A_9`.

## 3. Exact scalar Euler consequence and the missing Harnack theorem

Let `d_lambda` minimize the exact Wave-43 scalar objective

```math
F_\lambda(d)=h_d+\lambda R_2(d),qquad \lambda>0,
\tag{R44.8}
```

over nonempty columns.  Scalar optimality says for every block `U`

```math
h_{d_\lambda^U}-h_{d_\lambda}
\ge-\lambda[R_2(d_\lambda^U)-R_2(d_\lambda)].
\tag{R44.9}
```

Average (R44.9) over `|U|=k` and use (R44.7).  If `R=R_2(d_lambda)` is
above the baseline, then

```math
\boxed{
\mathbb E_{|U|=k}[h_{d_\lambda^U}-h_{d_\lambda}]
\ge
\lambda(1-\theta_k)[R-n(n-1)].}
\tag{R44.10}
```

Consequently the following conditional statement is exact:

```math
\boxed{
\mathbb E_{|U|=k}[h_{d_\lambda^U}-h_{d_\lambda}]\le G_k
\quad\Longrightarrow\quad
R_2(d_\lambda)
\le n(n-1)+\frac{G_k}{\lambda(1-\theta_k)}.}
\tag{R44.11}
```

For `k=Theta(n)`, `1-theta_k=Theta(1)`.  At the project price

```math
\lambda\asymp
\frac{n^{3/4-c}}{n^{9/4-c}}=n^{-3/2},
```

a geometric-mean orbit bound `G_k=O(n^{3/4-c})` would force
`R_2=O(n^{9/4-c})`.  To finish the leading route one would still have to
show `h_{d_lambda}=O(n^{3/4-c})`, or equivalently bound the scalar value at
the project scale.  Thus (R44.11) regularizes row but does not manufacture
the missing incidence mass.

The exact missing global-minimality property is now explicit: one needs
control of the **geometric mean** of the neighboring hard-column masses,
not merely pointwise nonemptiness of each selector fiber.  Equation
(10.1044) says only `max_d L_S(d)>=q` separately for each `S`.  It gives no
upper bound on the left side of (R44.10), which may be infinite when one
neighboring column is empty.

## 4. Exact finite cliff

The failure is visible on actual minimizers, not only on an abstract set
system.

- On `A_6,m=5`, the unique Pareto type has `(R,|I_d|,E)=(30,3,6)`.
  All three incidences have zero margin `L_S-q=0`, and none of the six
  one-vertex switches remains incident to the same selector.  Each complete
  threshold fiber nevertheless has 12 oriented projective states: fiber
  nonemptiness/size does not imply local retention around the chosen column.
- More sharply, on `A_9,m=7` there is an active column with `R=80>72=n(n-1)`,
  `|I_d|=2`, and `E=8`, for which **every** row-descending one-vertex switch
  has `|I_{d^i}|=0`.  Thus no universal assertion of the form “a high-row
  active column has an active row-descending neighbor” follows from exact
  minimality.  One certificate is
  `sigma=-1`, `x=(1,-1,-1,1,-1,-1,-1,1,1)`; its descending vertices
  `4,6,8` change the row by `-40,-32,-24`, respectively, and all three
  resulting columns are empty.
- Pareto columns also show mass loss under descent: on `A_8,m=6`, the two
  exact Pareto types are `(40,2,0)` and `(72,6,0)`; four of the eight
  row-descending neighbors of a representative of the latter are empty.

The second item is a finite, scoped falsifier of neighbor retention.  It is
not a scalable obstruction to the target column because its rows are only at
the natural `Theta(n^2)` scale.

## 5. Independent audit of the selector second moment

For completeness, put `p_j=(m)_j/(n)_j`.  Edge-pair counting gives

```math
\mathbb E c_S=p_2E,
```

and

```math
\boxed{
\operatorname{Var}(c_S)
=(p_4-p_2^2)E^2
+4(p_3-p_4)R
+2n(n-1)(p_2-2p_3+p_4).}
\tag{R44.12}
```

This is an independent rederivation of the existing identity (10.1123), not
a new theorem.  Incidence requires

```math
c_S\ge\frac{q+E}{2},
```

whose displacement above the mean is

```math
\frac q2+(1/2-p_2)E\ge(1-p_2)q
\tag{R44.13}
```

when `p_2>=1/2`.  The formula correctly exposes row as the variance-bearing
term, but a second moment supplies an upper tail bound, not the reverse tail
or orbit Harnack needed in (R44.11).  All coefficients in (R44.12) were
independently checked by exhaustive selector averaging for every oriented cut
of the four finite minimizers.

## 6. Exact finite scalar frontiers

Writing a Pareto point as `(R,|I_d|,E)`, exhaustive enumeration gives

```text
A5, m=4: (16,2,0)
A6, m=5: (30,3,6)
A8, m=6: (40,2,0), (72,6,0)
A9, m=7: (56,1,0), (72,2,0), (88,3,0)
```

The transitions are exact: for `A_8` they occur at
`lambda=log(3)/32`; for `A_9` at `log(3/2)/16` and `log(2)/16`.
These data confirm that the scalar objective genuinely trades incidence mass
against row; neither coordinate monotonically controls the other.

## 7. Surviving target

A positive switching-Euler continuation must prove a minimizer-specific
hard-column expansion estimate strong enough to make (R44.11) finite at
project cost, together with the project scalar-value bound.  A theorem about
arithmetic mean fiber size, first moments, or the existence of one witness per
selector is insufficient.  Alternatively one must attack (10.1149) directly:

```math
\min_d\{-\log U_m(\mathcal I_d)+\lambda R_2(d)\}
=O(n^{3/4-c})
\quad\text{for }\lambda\asymp n^{-3/2}.
```

The present wave neither proves nor falsifies this estimate.

## Verification

`tmp/r44_direct_incidence_check.py` verifies all switching, row, incidence,
and variance identities; exhausts the scalar Pareto frontiers; and checks the
finite neighbor cliffs.  It ends with

```text
PASS r44_direct_incidence_check
```
