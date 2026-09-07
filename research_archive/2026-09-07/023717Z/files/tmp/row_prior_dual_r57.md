# Wave 57: the exact row--incidence dual and why high row alone does not migrate into favorability

Status: the finite-law identities, LP duals, and the sector-positive migration
consequence below are verified.  The stored-minimizer statements were checked
by exact enumeration of all selectors and oriented projective cuts (apart from
the harmless floating evaluation of `(m/n)^(3/2)`, with a positive boundary
margin).  The hub-avoidance construction is a sharp scalable **abstract
incidence obstruction**, not a claimed sequence of exact-minimal signings.  No
convergence theorem or asymptotic signing counterexample is obtained.

## 1. Setup and the KL projection that exposes the real row certificate

Fix the bare favorable incidence

```math
F=F_t\subseteq\mathcal X_m\times\mathcal D,
\qquad
u_d=U_m(F^d),
\qquad
h_d=\log\frac1{u_d},
```

and discard the empty columns (`u_d=0`) when taking minima.  Write
`R_d=x^{\mathsf T}A^2x`.  If `P` is any joint selector--cut law supported on
`F`, let `pi=P_D` and let `U_d=U_m(.|F^d)`.  The chain rule and KL projection
onto a set give the exact Pythagorean identity

```math
\boxed{
K(P):=I_P(S;D)+D(P_S\Vert U_m)
=\sum_d\pi_d\left\{h_d+
D(P_{S|d}\Vert U_d)\right\}.}
\tag{R57.1}
```

Thus replacing every posterior by `U_d` preserves the output law and its row
cost and can only decrease `K`.  The lower attainable information--row region
is exactly the convex hull of the column points

```math
(h_d,R_d).
\tag{R57.2}
```

This is sharper than a statement about each selector having a witness and
automatically respects the matching obstruction (10.1357): in a perfect
matching every `h_d=log |X_m|`, so the entire convex hull has that information
coordinate.

For `H>=min_d h_d`, define the least mean row at information budget `H` by

```math
\mathcal R(H)
=\inf\{\mathbb E_PR_D:P(F)=1,\ K(P)\le H\}.
```

Projection in (R57.1) makes this an ordinary finite LP, and LP duality gives

```math
\boxed{
\begin{aligned}
\mathcal R(H)
&=\min_{\pi\in\Delta(\mathcal D)}
\left\{\sum_d\pi_dR_d:\sum_d\pi_dh_d\le H\right\}\\
&=\sup_{\lambda\ge0}
\left[\min_d\{R_d+\lambda h_d\}-\lambda H\right].
\end{aligned}}
\tag{R57.3}
```

An optimizer uses at most two columns.  More generally, a joint law with
`K(P)<=H` and `E R_D<=R_0` exists if and only if the convex hull in (R57.2)
meets the southwest rectangle at `(H,R_0)`.  Its exact obstruction certificate
is

```math
\boxed{
\exists\alpha,\beta\ge0,\ (\alpha,\beta)\ne(0,0):
\quad
\min_d\{\alpha(h_d-H)+\beta(R_d-R_0)\}>0.}
\tag{R57.4}
```

When `alpha>0`, this is the scalar gap

```math
\min_d\{h_d+\lambda R_d\}>H+\lambda R_0,
\qquad \lambda=\beta/\alpha;
\tag{R57.5}
```

the case `alpha=0` is simply `min_d R_d>R_0`.  Consequently, “every column
with `h_d=O(H)` has large row” is only a coarse shadow of the true
certificate.  It can miss a two-column tradeoff and it does not specify the
separating row price.

## 2. The exact selected-prior incidence LP

There is a second, slightly different, exact LP for priors.  Let `nu` be an
arbitrary selector-independent cut prior and put

```math
q_d=\nu_du_d,
\qquad Z_\nu=\sum_dq_d.
```

Then `sum_d q_d/u_d=1`, and the captured output law is `pi_d=q_d/Z_nu`.
The largest selected-prior incidence whose captured mean row is at most `R_0`
is therefore

```math
\boxed{
\begin{aligned}
Z_{\rm sel}(R_0)=\max_{q_d\ge0}\quad&\sum_dq_d\\
\text{subject to}\quad
&\sum_d\frac{q_d}{u_d}=1,\\
&\sum_dq_d(R_d-R_0)\le0.
\end{aligned}}
\tag{R57.6}
```

Its direct LP dual is

```math
\boxed{
Z_{\rm sel}(R_0)=
\min_{y\in\mathbb R,\ \gamma\ge0}
\left\{y:
\frac{y}{u_d}+\gamma(R_d-R_0)\ge1\ \text{for every }d
\right\}.}
\tag{R57.7}
```

Equivalently, parametrizing by the captured law gives

```math
\boxed{
Z_{\rm sel}(R_0)^{-1}
=\min_{\pi:\,\mathbb E_\pi R\le R_0}
\mathbb E_\pi e^{h_D},
\qquad
-\log Z_\nu=\log\mathbb E_\pi e^{h_D}.}
\tag{R57.8}
```

Again an optimum uses at most two columns.  Notice the distinction between
the two exact costs:

```math
\min K(P)=\mathbb E_\pi h_D,
\qquad
-\log Z_\nu=\log\mathbb E_\pi e^{h_D}.
\tag{R57.9}
```

The latter includes the output-prior KL term in (10.1354).  Jensen orders the
two but does not identify them.

For a **fixed** reference prior `nu`, row penalization has the exact KL dual

```math
\boxed{
\begin{aligned}
Z_{\nu,\gamma}
&=\sum_d\nu_du_de^{-\gamma R_d},\\
-\log Z_{\nu,\gamma}
&=\inf_{P:\,P(F)=1}
\{K(P)+D(P_D\Vert\nu)+\gamma\mathbb E_PR_D\}.
\end{aligned}}
\tag{R57.10}
```

The optimizing law is proportional to
`U_m(S)nu(d)1_F(S,d)e^(-gamma R_d)`.  Optimizing the reference prior itself
collapses back to one column:

```math
\boxed{
-\log\sup_\nu Z_{\nu,\gamma}
=\min_d\{h_d+\gamma R_d\}.}
\tag{R57.11}
```

Equations (R57.3)--(R57.11) are the requested exact row-penalized
LP/KL dual.  They do **not** say that KL to an arbitrary `nu` controls row.
Such a statement is false: `nu` may be singular and supported entirely on
high-row cuts.  Row control must be imposed in the LP, obtained from a
reference measure with a proved row moment bound, or proved by signing
geometry.

## 3. What high row does give in the correct energy sector

There is a clean but insufficient migration consequence.  For an oriented
parent state `d=(tau,x)`, gauge the edge signs and fields by

```math
s_{ij}(d)=\tau a_{ij}x_ix_j,
\qquad
r_i(d)=\sum_{j\ne i}s_{ij}(d).
```

Then

```math
\sum_i r_i=q-\Delta(d),
\qquad
R_d=\sum_i r_i^2,
\qquad
|r_i|\le n-1.
\tag{R57.12}
```

Suppose `Delta(d)<=q`, so the oriented energy is nonnegative.  With
`P_+=sum_i(r_i)_+` and `P_-=sum_i(-r_i)_+`, one has `P_-<=P_+`; hence

```math
\boxed{R_d\le2(n-1)P_+.}
\tag{R57.13}
```

For every integer `L<=P_+`, choose integers `0<=k_i<=r_i^+` summing to `L`.
At a positive-field vertex there are `(n-1+r_i)/2>=r_i` incident positive
signed edges.  Select `k_i` of them and take the union `E`.  Since an edge is
counted at most twice,

```math
L/2\le |E|\le L,
\qquad s_e(d)=+1\quad(e\in E).
\tag{R57.14}
```

Take `0<c<1/4`, `L asymp n^(5/4-c)`, and assume
`R_d>=2(n-1)L`.  Apply the verified localized migration theorem (10.1323) to
`E` with its parameter `lambda_0=1/2`.  Here

```math
\eta_E=\sqrt{8|E|n\log2}
=O(n^{9/8-c/2})=o(T_n),
\qquad \eta_E/|E|=o(1).
```

Taking response tolerance `5 eta_E` gives a parent state `omega` with

```math
\boxed{
\Delta(\omega)=O(n^{9/8-c/2})=o(T_n),
\qquad
|\{e\in E:s_e(\omega)=-1\}|>|E|/4.}
\tag{R57.15}
```

Thus a near-side cut with row at or above the project threshold really does
force a quantitatively separated near-ground response.  If
`R_d/n^(9/4-c)->infinity`, the block size can also be multiplied by a slowly
diverging factor.  This is a one-response theorem, not an abundance theorem:
it controls neither `R(omega)` nor `u_t(omega)`, and many input columns may
share the same response.

## 4. Two sharp walls to turning (R57.15) into favorability

### 4.1 Exact orientation wall on a stored minimizer

For the stored exact minimizer `A_8`, at `m=5` and `t=0`, exhaustive
enumeration gives `q=20` and 56 selectors.  All eight maximum bare-incidence
columns have

```math
u_t(d)=40/56,
\qquad h_d=0.3364722366\ldots,
\qquad R_d=64,
\qquad \Delta(d)=40=2q.
\tag{R57.16}
```

Every signed field is strictly negative; its positive-sector mass `P_+` is
zero.  Reversing the orientation makes the same physical word an exact parent
ground, keeps row 64, and lowers its coverage to `30/56`.  Therefore even on
an exact minimizer, maximal bare favorability plus high row need not expose
any positive field block: the parent-deficit subsidy in (10.1352) can select
the opposite energy sector.  This is an exact finite, scale-free wall.  It is
not an asymptotic counterexample to a theorem that adds `Delta<=q` or an
equivalent sector hypothesis.

The finite Pareto curve is genuinely two-dimensional.  At captured row
budget 40, (R57.6) on `A_8` uses captured weights `1/4,3/4` on

```math
(u,R)=(40/56,64),\qquad (39/56,32),
```

and has `Z_sel=260/371=0.7008086253...`.  This is larger than the best
single-column coverage under the strict row cap 40, namely `39/56`.
Randomization only convexifies the exact tradeoff; it does not make row a
consequence of arbitrary-prior KL.

### 4.2 Scalable hub-avoidance wall

There is a sharper scale obstruction to any attempt to infer fibre overlap
from the entropy budget alone.  Let

```math
k=\omega_n n^{1/4-c},
\qquad \omega_n\longrightarrow\infty,
\qquad \omega_n=o(n^{1/4+c}),
```

and let `K` be a fixed set of `k` vertices.  The selector family avoiding it
has exact mass

```math
u_K=\frac{\binom{n-k}{m}}{\binom nm}
=\frac{(n-m)_k}{(n)_k},
```

so at fixed density

```math
-\log u_K=k\log\frac1{1-p}+O(k^2/n)=o(H).
\tag{R57.17}
```

Concentrating fields of order `n` on these `k` vertices contributes

```math
R\asymp kn^2
=\omega_n n^{9/4-c}\gg n^{9/4-c}.
\tag{R57.18}
```

The positive-edge certificate extracted from those fields can be supported
entirely on edges incident to `K`.  Every selector in the low-information
family above sees none of those edges internally.  Exact-minimizer migration
may return a state that reverses a constant fraction of the certificate, but
it gives no control on that state's signs away from the certificate.  Hence
it supplies no implication for `F_t` on the avoiding family.

This can be made into an exact abstract incidence/response model: declare
`F^d={S:S cap K=emptyset}`, give `d` row (R57.18), define its local deficits
so that (10.1352) cuts out exactly this fibre, and attach the near-ground
response states required by (R57.15) with their favorable fibres empty.  The
KL projection identity, the matching-safe collision bound, the bare-deficit
identity, and the stated migration conclusion all remain true.  What is
absent is precisely a theorem coupling the response signs outside the chosen
edge block to the original favorable fibre.  This is an incidence-level
logical obstruction, not a realization by an exact-minimal signing.

The exponents are sharp for the proposed conversion: the row threshold can
live on only `R/n^2=n^(1/4-c)` hub vertices, whereas the allowed information
budget is `H=n^(3/4-c)`, larger by `n^(1/2)`.  Entropy therefore cannot force
the favourable selectors to visit the high-field support.  The upper
condition on `omega_n` also keeps the hubs' total positive field
`k n=o(n^(3/2))`, so this profile is compatible with the parent energy scale.

## 5. What additional hypothesis would actually suffice

The high-row certificate alone cannot be converted into a contradiction or
useful abundance theorem using the currently proved perturbation statements.
Two independent additions are needed:

1. **Correct-sector control.**  Require `Delta(d)<=q` (preferably
   `Delta=o(T_n)`) for the relevant low-information columns, or directly
   require the positive-sector field mass in (R57.13).  This excludes the
   exact `A_8` opposite-orientation wall and activates (R57.15).

2. **Fibre-stable response or row descent.**  A response theorem must control
   more than one selected edge block.  One exact sufficient version is: for
   each separating price `lambda` in (R57.5), every row-high favorable column
   `d` in the correct sector has an exact-minimality response `omega` such
   that

   ```math
   h_\omega+\lambda R_\omega
   <h_d+\lambda R_d.
   \tag{R57.19}
   ```

   Applied to a minimizer of the scalar column objective, this contradicts
   the obstruction certificate.  A more structural sufficient substitute is
   a saved overlap

   ```math
   U_m(F^d\cap F^\omega)\ge e^{-O(H)}
   \tag{R57.20}
   ```

   together with a proved row descent to `R_omega=O(n^(9/4-c))`.  Merely
   showing that favorable selectors contain a positive fraction of the
   migrated block is not enough unless the uncontrolled outside-edge change
   is also bounded.

A delocalization/visibility condition ruling out the hub-avoidance family is
necessary for any entropy-only implementation of (R57.20), but by itself it
does not give row descent.  Similarly, iterating (R57.15) gives no packing:
without a bounded-preimage or separation hypothesis, one response can serve
many high-row columns.

## 6. Resulting judgment

The exact obstruction is the scalar convex-hull gap (R57.4)--(R57.5), not the
verbal statement that all low-information fibres have high row.  Exact
minimality does turn a **correct-sector** high-row column into a separated
near-ground state at deficit `o(T_n)`, but current migration has no mechanism
that preserves favorable incidence or lowers project row.  Opposite
orientation and hub avoidance independently show why.  The next viable row
attack must prove a fibre-stable, row-descending response such as
(R57.19), or construct the selected prior directly through (R57.6).  No claim
that arbitrary-prior KL controls row survives this audit.

Reproducible audit: `tmp/row_prior_dual_r57_check.py`.
