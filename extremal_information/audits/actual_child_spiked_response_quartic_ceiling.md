# The quartic spiked response is submacroscopic for actual thermal children

Status: **rigorous actual-child perturbative ceiling, not a solution of the
physical-scale response problem**.  The degree-two spiked carrier exposes an
exponential switching-orbit response table.  This note computes its first
nonconstant temperature coefficient exactly.  For contracted-temperature
minimizing children in a balanced split, that coefficient contributes only
`O(sqrt(N))` at the physical scale.  Therefore a linear excursion of the
spiked response, if it exists, is necessarily a genuinely collective
higher-order phenomenon; it cannot be certified by the quartic/sector--Gram
tangent.

The result uses the actual minimizing children.  It neither bounds the full
physical response nor constructs a scalable collision, so bounded-row-degree
cross-row synchronization remains open.

## 1. Setup

Let `A,D` be hollow signings of orders `m,n`, let `N=m+n`, and put

```math
 H_{A,D,B}^{\epsilon}(x,z)
 =H_A(x)+\epsilon H_D(z)+x^{\mathsf T}Bz.           \tag{SQ.1}
```

For an auxiliary real inverse temperature `s`, define

```math
 L_s(B)=\log E_{x,z}\cosh\bigl(sH_{A,D,B}^{\epsilon}(x,z)\bigr).
                                                               \tag{SQ.2}
```

Fix `y in {+-1}^n`.  With the exact degree-two row density from Theorem
37.40,

```math
 q_{v_i,y}(b)={1\over2}
 \left(1+{v_i\langle y,b\rangle\over\sqrt n}\right)^2,
 \qquad
 P_{v,y}=\bigotimes_{i=1}^m q_{v_i,y}U_n,           \tag{SQ.3}
```

write

```math
 R_s(v)=E_{P_{v,y}}L_s(B).                          \tag{SQ.4}
```

The children below may be selected at the physical temperature
`t=beta/sqrt(N)`; the auxiliary path `s -> R_s` keeps those children fixed.

## 2. Exact first response coefficient

**Theorem SQ.1 (quartic switching-orbit formula).**  For any two query words
`v,w in {+-1}^m`,

```math
 \boxed{
 [s^4]\{R_s(v)-R_s(w)\}
 ={1\over2}\bigl(\|Av\|_2^2-\|Aw\|_2^2\bigr)
 +{2\epsilon H_D(y)\over n}
       \bigl(H_A(v)-H_A(w)\bigr).}                 \tag{SQ.5}
```

Equivalently, the fourth derivative at zero is `4!` times the right-hand
side.  The quadratic coefficient is independent of `v`, and all odd
coefficients vanish.

*Proof.*  If `H=H_(A,D,B)^epsilon`, then

```math
 \log E\cosh(sH)
 ={s^2\over2}EH^2+{s^4\over24}
     \{EH^4-3(EH^2)^2\}+O(s^6).                   \tag{SQ.6}
```

The second moment is the number of edges and is independent of every sign.
For a Rademacher quadratic form with edge coefficients `c_e`, its fourth
cumulant is

```math
 -2\sum_ec_e^4
 +24\sum_{C_4}\prod_{e\in C_4}c_e,                \tag{SQ.7}
```

where the last sum is over unoriented simple four-cycles.  This follows by
classifying the even four-edge multigraphs: two repeated pairs cancel the
Gaussian subtraction, one edge used four times contributes `-2`, and four
distinct edges contribute precisely when they form a four-cycle.

Under (SQ.3), for distinct coordinates `j,k` and distinct rows `a,c`,

```math
 E B_{aj}={v_ay_j\over\sqrt n},\qquad
 E(B_{aj}B_{ak})={y_jy_k\over n},\qquad
 E(B_{aj}B_{ck})={v_av_cy_jy_k\over n}.            \tag{SQ.8}
```

Four-cycles on four left vertices or four right vertices are independent
of `v`.  So are alternating bridge four-cycles.  A cycle on three left
vertices and one right vertex contributes, after summing the right vertex,

```math
 {1\over2}\sum_b
 \left\{\left(\sum_{a\ne b}A_{ab}v_a\right)^2-(m-1)\right\}
 ={1\over2}\{\|Av\|_2^2-m(m-1)\}.                 \tag{SQ.9}
```

The analogous one-left/three-right term is independent of `v`.  Finally,
on two left and two right vertices, the two cycles containing one `A` edge,
one `epsilon D` edge, and two bridge edges sum to

```math
 {2\epsilon\over n}H_A(v)H_D(y).                  \tag{SQ.10}
```

Subtracting the formulas for `v` and `w`, and using the factor `24` in
(SQ.7) against `1/24` in (SQ.6), proves (SQ.5). `square`

## 3. Actual-minimizer scale

Let

```math
 p_A(t)=\log E_x\cosh(tH_A(x)),\qquad
 Q(A)=\max_x|H_A(x)|.                              \tag{SQ.11}
```

If `A` minimizes `p_A(t)` among order-`m` signings, randomizing all edge
signs and averaging the unlogged partition function gives

```math
 p_A(t)\le {m\choose2}\log\cosh t.                \tag{SQ.12}
```

On the other hand, retaining one extremizing projective spin gives

```math
 p_A(t)\ge tQ(A)-m\log2.                           \tag{SQ.13}
```

Consequently

```math
 Q(A)\le
 {{m\choose2}\log\cosh t+m\log2\over t}.         \tag{SQ.14}
```

For a symmetric hollow matrix, elementary polarization gives
`max_(r,v)|r^TAv|<=4Q(A)`.  Therefore, for every Boolean `v`,

```math
 \|Av\|_2^2
 \le\|Av\|_\infty\|Av\|_1
 \le4(m-1)Q(A).                                   \tag{SQ.15}
```

**Corollary SQ.2 (submacroscopic quartic range).**  Fix `beta>0` and a
balanced window `theta N<=m,n<=(1-theta)N`.  If `A,D` are the actual
contracted-temperature minimizers at `t=beta/sqrt(N)`, then

```math
 \boxed{
 t^4\max_{v,w}|C_4(v,w)|
 =O_{\beta,\theta}(\sqrt N)=o(N).}                 \tag{SQ.16}
```

Here `C_4(v,w)` denotes the right-hand side of (SQ.5).

More explicitly, the absolute value of the bracket in (SQ.5) is at most

```math
 4(m-1)Q(A)+{4Q(A)Q(D)\over n}.                    \tag{SQ.17}
```

Equations (SQ.14), `log cosh t<=t^2/2`, and `t=beta/sqrt(N)` turn the first
term after multiplication by `t^4` into `O(sqrt(N))` and the second into
`O(1)`.

Equation (SQ.16) evaluates only the quartic coefficient at the physical
scale.  It does not evaluate the complete analytic series at `s=t`.

## 4. Consequence and limitation

The exact small actual-child computations show a response dominated by row
order two, but they cannot establish its asymptotic range.  Corollary SQ.2
explains why the first such mode cannot itself be extensive: the child
minimization inequality already reduces it to `O(sqrt(N))`.

Thus a scalable fixed-density excursion of `R_t` cannot arise from the
quartic sector--Gram tangent alone.  It must first enter at sixth or higher
order; SQ.2 does not distinguish a fixed higher diagram from a genuinely
collective resummation.  In particular, SQ.2 gives no remainder estimate at
the physical point.  It does not prove `max R_t-min R_t=o(N)`, decide the
product phase, or change the bounded-row-degree cross-row SML.
