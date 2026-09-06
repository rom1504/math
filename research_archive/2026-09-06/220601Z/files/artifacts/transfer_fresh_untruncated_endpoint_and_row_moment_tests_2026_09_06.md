# Exact untruncated endpoint and first-marked quantitative tests

2026-09-06. Follow-up to the certified principal-deletion loss barrier.
All results in Sections 1--4 are finite statements for original hollow
signings. No bounded operator norm is assumed. They do not establish an
unrestricted feedback gain or replace the remaining Stein brackets.

Normalize `B=A/sqrt(n-1)`, and put

```
Lambda(B)=max_(x Boolean) |x^T B x|/(2n) <= C,
beta(B)=max_(x,y Boolean) |x^T B y| <= 4 C n.
```

An original bound `Q(A)<=C0 n^(3/2)` gives
`Lambda(B)<=C0 sqrt(n/(n-1))`; any strict asymptotic upper slack may be
absorbed into a fixed C. Write KG for any valid real Grothendieck constant.

## 1. A cap-only randomized spin-ascent cost

Let u be any possibly random Boolean vector, let sigma be either energy
orientation, and set

```
e_B(u)=u^T B u/(2n),
d=(1/n) E sum_i [-sigma u_i (Bu)_i]_+.
```

For each sample u, put `v_i=u_i 1_{sigma u_i(Bu)_i<0}`. Independently flip
each selected spin with probability p. Since B is hollow, the exact
conditional expected change in oriented energy is

```
2p (1/n) sum_i [-sigma u_i(Bu)_i]_+
       +(2 sigma p^2/n) v^T B v.
```

Every `v in [-1,1]^n` satisfies `|v^T B v|<=2Cn`: the hollow quadratic
is affine in each coordinate, or equivalently independent Boolean rounding
of v preserves its expected quadratic energy. Consequently

```
E[sigma e_B(u_new)-sigma e_B(u)] >= 2 p d - 4 C p^2.
```

Also `d<=E||Bu||1/n<=beta(B)/n<=4C`. Thus the admissible choice
`p=d/(4C)<=1` gives the exact improvement

```
E sigma e_B(u_new) >= E sigma e_B(u) + d^2/(4C).       (1)
```

The displacement is `u_new=u-2v` on flipped coordinates; it is important
that the coefficient 2 is allowed. Artificially restricting displacement
to coefficient at most 1 would unnecessarily clip this optimization.

For first-marked endpoints with means `E e_B(u_+)=b+j` and
`E e_B(u_-)=b-j`, if the positive instability at u+ and the negative
instability at u- are both at least d0, then (1) gives

```
Lambda(B) >= j + |b| + d0^2/(4C).                    (2)
```

Thus the deterministic endpoint cost itself does not require an operator
cap. Obtaining the positive instability is a separate obligation.

## 2. Cap-only averaged root standard deviations

For ANY random vector X on ANY probability space, suppose
`sup_j E X_j^2<=v^2`. Then

```
(1/n) sum_i sqrt(E(BX)_i^2) <= KG beta(B) v/n
                            <=4 KG C v.             (3)
```

Proof: regard `h_j=X_j/v` as unit-ball vectors of the real Hilbert space
L2. For `r_i=sum_j B_ij h_j`, choose the Hilbert unit vector
`u_i=r_i/||r_i||` when nonzero, and zero otherwise. The exact bilinear
vector objective is

```
sum_(i,j) B_ij <u_i,h_j> = sum_i ||r_i||.
```

Grothendieck's inequality proves (3). No independence or centering of X is
needed. For v=0 the assertion is immediate.

In particular, the rows of B have Euclidean norm exactly one, so
`E G_i^2=1` for `G=BS`. Applying (3) to X=G gives

```
(1/n) sum_i sqrt((B^4)_ii) <=4 KG C,                 (4)
```

because `BG=QS`, `Q=B^2`, and `E(QS)_i^2=(B^4)_ii`.
The equivalent deterministic vector construction uses the rows of B as
the h_j and normalized rows of B squared as the u_i.

This is a bound on the average STANDARD DEVIATION, not on the average
variance or the operator norm.

## 3. The full old literal list has cap-only root moment caps

Keep exactly

```
D_i=S_i h2(G_i),   Y=BD,
W=(S,G,Y,QS,QD),   F=f(G,Y),
H=1-F^2,   C_field=H sign(BF),
T=BF,   K=BC_field.
```

Here F is ternary, so both F and C_field lie pointwise in [-1,1]. The
letter C in (1)--(4) is the numerical cap bound, not C_field.

The exact covariance, already banked in the marked-return archive, is

```
Cov(D)=[1-3/(n-1)] I + [2/(n-1)] Q.                 (5)
```

Indeed `D_i=sqrt(2) sum_(j<k, j,k!=i) B_ij B_ik S_i S_j S_k`.
Its variance is `1-1/(n-1)`, while for i!=j the only matching triples
give `E D_i D_j=2 Q_ij/(n-1)`. Since `||Q||op<=n-1` for every signing,
equation (5) gives `Cov(D)<=3I`, hence `E Y_i^2<=3` at EVERY root.

Apply (3) to X=G,Y,F,C_field. Together with (4), this gives

```
avg_i sqrt(E(QS)_i^2) <=4 KG C,
avg_i sqrt(E(QD)_i^2) <=4 sqrt(3) KG C,
avg_i sqrt(E T_i^2)   <=4 KG C,
avg_i sqrt(E K_i^2)   <=4 KG C.                     (6)
```

Define the deterministic root statistic

```
m_i=sqrt(E[(QS)_i^2+(QD)_i^2+T_i^2+K_i^2]).
```

Then `avg_i m_i<=4KG C(3+sqrt(3))`. For every fixed R>0, therefore,

```
#{i:m_i>R}/n <=4KG C(3+sqrt(3))/R.                  (7)
```

On all remaining roots the four displayed fields have second moment at
most R squared; S,G,Y already have uniform second moments 1,1,3. This
simultaneous root cutoff does not change B or delete any vertices from
the original optimization problem.

## 4. An untruncated application of the cubic chain theorem

The coefficient of `S_a S_b S_c`, for distinct a,b,c, in Y_i is exactly

```
sqrt(2) [B_ia B_ab B_ac + B_ib B_ab B_bc
                         + B_ic B_ac B_bc].
```

It has absolute value at most `3sqrt(2)/(n-1)^(3/2)`. Consequently

```
max_a E(Delta_a Y_i)^2 <=9(n-2)/(n-1)^2,
sum_a E(Delta_a Y_i)^2 =3 E Y_i^2 <=9.              (8)
```

Y is a PURE cubic Fourier polynomial; there is no hidden linear part in
these formulas. For the two-dimensional query `V_i=(G_i,Y_i)`, equations
(8) and the exact coefficients of G give

```
sup_i sum_a I_a(V_i)<=10,
max_(i,a) I_a(V_i)<=10/(n-1).                       (9)
```

Let Z_j be pure degree-K probe polynomials with
`max_(j,a) E(Delta_a Z_j)^2<=c/n`, and let the smooth tests f_i satisfy
`sup_i ||D3 f_i||op,infty<=M`. The director's independently audited cubic
chain theorem now applies with d=3, I0=10, and (9), with NO operator cap:

```
R_ij=E[Z_j f_i(G_i,Y_i)]
 -(1/K) sum_l E[partial_l f_i(G_i,Y_i)
                     sum_a Delta_a Z_j Delta_a V_il],

||R||_*/n <= [10 M sqrt(c)/(3K)] 3^((K+5)/2)
                                      sqrt(10/(n-1)).           (10)
```

This is a quantitative original-signing application, not a proposed
conjecture. It controls only the CHAIN ERROR for fixed smooth old tests.
The brackets remain actual and must still be analyzed. Adding the literal
own sign S_i or the coherent QS,QD queries does not preserve (9): those
can have order-one coordinate influences. No hard-threshold passage or
Gaussian law for Y without an operator cap is claimed here.

## 5. The exact remaining quantitative checkpoint

The current feedback proof extracts a useful positive return covariance
on roots satisfying simultaneous query/probe caps. Section 3 now supplies
the query caps from the original Boolean cap alone. It does NOT show that
the return covariance survives that cutoff. For a proposed probe E_i, a
concrete sufficient checkpoint is to find fixed R,M,delta,a>0 such that

```
liminf_n (1/n) sum_(i: m_i<=R, delta<=E E_i^2<=M)
                         E[(BC_field)_i E_i] >= a.              (11)
```

One must separately supply the joint literal regressed-query comparison
on these roots, including the actual endpoints, T, and
`U_i=K_i-[E K_i E_i/E E_i^2]E_i`. Positive covariance or independent
one-dimensional probe marginals alone do not supply that premise.

For clarity there is an explicit conditional constant. Set

```
rho=a/(4R sqrt(M)),   s=a/(4 sqrt(M)),
h=E(sN-2R)_+,
d0=min(4C,rho h/4),   gain=d0^2/(4C)>0.              (12)
```

For all sufficiently large n the covariance in (11) is at least an/2.
Discarding its roots with covariance less than a/4 leaves at least rho n
roots, since each covariance is at most R sqrt(M). The Gaussian probe
component of K there has standard deviation at least s; its regression
coefficient is at most R/sqrt(delta). The raw second moment identity
`E U_i^2=E K_i^2-(E K_i E_i)^2/E E_i^2` needs only that E_i is centered.
Thus `E(plus/minus T_i+U_i)^2<=2R^2`, so the absolute coherent shift is
at most 2R with probability at least one half. Independent probe noise
then gives each endpoint instability density at least rho h/2 in the
limit. The harmless extra factor two in d0 absorbs approximation slack.
All relevant row second moments are fixed on these roots, licensing the
same linear-growth truncation passage as in the endpoint audit.

Consequently (2) gives gain (12), with no principal retention factor.
Equation (11) is a sufficient target for this argument,
not asserted necessary for every possible proof.

Bounded average root standard deviation alone does not imply (11).
For example, in abstract scalar arrays take a vanishing fraction theta
of roots with `K_i=theta^-1 N_i`, `E_i=N_i`, and take both fields zero
elsewhere. The average standard deviation of K and the average covariance
are both one; probe variances are at most one. Yet every fixed K-moment
cutoff removes the entire covariance. This example is not presented as
an original signing or as a counterexample to its additional structure.
It shows exactly why (3)--(7) cannot by themselves complete the argument.

No terminal-response optimization is performed, and no new unrestricted
decimal lower bound or convergence statement follows from this note.

## 6. Provenance

The exact D covariance and its operator bound were independently rederived
here but are NOT new: see `continued_director_marked_return_regression_2026_09_06.md`,
Section 2, and `continued_feedback_boundary_graph_and_cubic_return_audit_2026_09_06.md`,
Section 4. The chain theorem is
`transfer_director_boolean_cubic_chain_nuclear_remainder_2026_09_06.md`,
independently checked in `transfer_fresh_boolean_cubic_chain_remainder_audit_2026_09_06.md`.
The cap-only Hilbert-space consequences, exact cap-only endpoint constants,
and the explicit all-signing first-marked chain application are the scoped
conclusions recorded here.

Finite integer replay: `python3 computations/transfer_fresh_untruncated_first_marked_verify.py`
passed 120 signings at orders 2 through 16, including positive cliques.
The small cases enumerate 4,064 cube states and 8,128 oriented ascent
checks. The checker verifies the exact cubic coefficients, D covariance,
variance/influence bounds, and rational cap-only ascent constant. This
supports the algebra only; it is not an asymptotic feedback test.
