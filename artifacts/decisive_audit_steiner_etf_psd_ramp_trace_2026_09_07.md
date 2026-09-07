# Independent audit: Steiner full signs force linear PSD ramp trace

Date: 2026-09-07. Status: independent proof reconstruction PASS.

Audited canonical theorem:
`artifacts/decisive_independent_steiner_etf_ramp_trace_obstruction_2026_09_07.md`.
This is a counterexample to universal ramp availability for bounded-cap
parents. Neither width-minimality nor a failure of the conditional sign
surgery is asserted.

## 1. The exact imported existence statement

I read the primary [Fickus--Mixon--Tremain paper](https://arxiv.org/pdf/1009.5730)
directly. Theorem 1, printed page 4, constructs a real ETF from a Steiner
system S(2,k,v) and a real Hadamard matrix of order
`1+(v-1)/(k-1)`. Section 3.1.4, printed page 9, states the existence of
S(2,5,v) for v congruent to 1 or 5 modulo 20. Alternatively the fixed-k,
sufficiently-large admissible existence statement preceding and used in
Theorem 2, printed pages 11--13, suffices here.

For r=16^j-1, v=4r+1, both r and vr/5 are integers, v is 1 modulo 20,
and r+1 is an unconditional Sylvester Hadamard order. No unproved
Hadamard existence, growing-k existence theorem, or all-order density
conclusion enters this audit.

## 2. Independent frame and sign normalization

The design has d=vr/5 blocks, replication r, and n=v(r+1) frame columns.
At each point, put the r nonconstant rows of the normalized-sign
Hadamard matrix into its incident block positions, then scale by
1/sqrt(r). This produces F of size d by n.

Distinct design blocks share at most one point; at that point they use
distinct Hadamard rows. Thus frame rows are orthogonal, each with
squared norm 5(r+1)/r. A column has r entries of squared magnitude 1/r.
For columns at the same point, deleting the constant Hadamard row gives
inner product -1/r. Columns at distinct points overlap at exactly one
design block, giving inner product +/-1/r. Consequently

```math
G=F^T F,qquad G_{ii}=1,qquad
A=r(G-I_n)in\{0,\pm1\}^{n\times n}
```

is hollow and has no zero off-diagonal entry. Its eigenvalues are exactly

```math
4r+5\quad(d\text{ times}),\qquad
-r\quad(n-d\text{ times}),\qquad
n=4r^2+5r+1.
```

For P=max H_A, R=max(-H_A), w=(P+R)/2, I=(P-R)/2, the spectral bounds give

```math
P\le n(4r+5)/2,\quad R\le nr/2,\quad
w\le n(5r+5)/4,\quad I\le n(4r+5)/4.                 \tag{1}
```

The last bound uses R>=0, which follows from the uniform-spin mean of
the hollow quadratic being zero. In particular Q(A)/n^(3/2)<=1+o(1)
and ||A||op/sqrt(n)->2. Nothing here says that these signs minimize Q
or width.

## 3. Exact Boolean covariance, including the singular case

Let Z be centered Gaussian with covariance G and X=sign(Z). Although G
is singular, every coordinate variance is one, so the signs have no
tie ambiguity and are individually centered. The bivariate Gaussian
angle formula follows by representing the two unit Gaussian directions
in their two-dimensional span. It applies to the present correlations
without any nondegeneracy assumption on the full Gaussian vector.

Set a_r=(2/pi)arcsin(1/r). Oddness of arcsin and exact flatness of the
off-diagonal Gram magnitudes imply the matrix identity

```math
K=\mathbb E XX^T=I_n+a_r A,qquad
\|K\|_{\rm op}=1+a_r(4r+5).                           \tag{2}
```

This is an exact covariance, not an entrywise approximation promoted to
an operator statement. Since Tr(A^2)=n(n-1),

```math
\mathbb E H_A(X)=\tfrac12\operatorname{Tr}(AK)
 ={a_r n(n-1)\over2}
 =\left({2\over\pi}+o(1)\right)n^{3/2}.               \tag{3}
```

Combining (3) with R<=nr/2 shows I>0 eventually. Combining (3) with
the bound on w gives a positive gap because 2/pi>5/8.

## 4. Finite trace inequality and asymptotic constant

Suppose T>=0 is arbitrary and satisfies the literal all-spin ramp

```math
H_A(x)\le w+{I\over n}x^T T x+e_n
\quad\text{for every }x\in\{\pm1\}^n.                 \tag{4}
```

No upper bound T<=I, rank condition, or sign compatibility is needed.
From (2), PSD ordering gives Tr(TK)<=||K||op Tr(T), regardless of
whether T and K commute. Average (4), then use the upper bounds in (1).
Whenever the numerator is positive, the resulting exact finite bound is

```math
{\operatorname{Tr}T\over n}
\ge
{2a_r(n-1)-5(r+1)-4e_n/n
 \over (4r+5)(1+a_r(4r+5))}.                           \tag{5}
```

For e_n=o(n^(3/2)), divide numerator and denominator by r. Since
n/r^2->4 and r a_r->2/pi, (5) yields

```math
\liminf_j {\operatorname{Tr}T\over n}
\ge {16-5\pi\over4(\pi+8)}>{1\over156}.               \tag{6}
```

For the strict rational bound, f(z)=(16-5z)/(4(z+8)) has derivative
-14/(z+8)^2<0 and f(22/7)=1/156. Hence pi<22/7 proves (6), with no
floating-point positivity test.

In particular all o(n)-trace PSD contractions are excluded, including
every sign-compatible one. This does not exclude a different sign
repair, nor prove anything about the existence of a useful ramp for a
selected sequence of exact or asymptotic width minimizers.

## 5. Reproducible finite check

I read and reran
`computations/decisive_independent_steiner_etf_ramp_check_2026_09_07.py`
with `.venv/bin/python`. It passed all assertions. The three cyclic
base blocks have each nonzero ordered difference exactly once; their
183 translates cover all 1830 point pairs once. The integer frame has
shape 183 by 976 and exact row Gram 80 I, and its column Gram minus
15 I is a hollow full-sign matrix. Rank, trace, and squared-trace
identities give eigenvalues 65 and -15 without numerical eigensolvers.

The exact rational replay returned

```text
E H_A(X) >= 222040/11
w <= 19520
||K||op <= 86/21
zero-error ramp Tr(T) >= 61488/6149.
```

The last line is a supplementary finite instance, not the proof of the
linear asymptotic lower bound. The only correction requested of the
canonical note was bibliographic: Theorem 1 is on printed page 4, not
printed page 3. No mathematical correction was needed.
