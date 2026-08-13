# Asymmetric Krivine outputs: an exact same-spin recoupling theorem

## Status

**Proved structural bridge; no universal bound improvement is claimed.**
For a pair of genuinely asymmetric Boolean outputs \(x,y\), the ordinary
polarization loss is not intrinsic unless the quadratic energies on their
agreement and disagreement shores have opposite signs. In the opposite-sign
case, the entire loss is exactly localized to one *aligned one-sided cap* of
one principal block.

That exact cap is still a Boolean optimization. The new usable result is a
polynomial-time replacement: each one-sided cap has a rigorous lower
certificate from the nuclear norm, or a stronger sign-specific spectral
projector certificate. Combining this certificate with one fixed shore gives
a same-spin bound at the \(n^{3/2}\) scale on flat/conference matrices.

This route strictly escapes both prior algebraic no-gos:

- it does not use the same Gaussian response at the two original endpoints;
- it does not separately round or bound the two full vectors \(x\) and \(y\);
- the cross block is retained and shown never to hurt, rather than discarded;
  and
- the statistic used by the polynomial theorem is an eigendecomposition of a
  specified principal block, not the unknown value \(Q(A)\).

The result does not yet improve the project lower bound. To do so, an
asymmetric rounding must also be shown to have a subleading recoupling defect.

## 1. Normalization

Throughout the theorem and the checker, use doubled energy

```math
Q(A)=\max_{z\in\{\pm1\}^n}|z^{\mathsf T}Az|
```

and bilinear cap

```math
B(A)=\max_{x,y\in\{\pm1\}^n}|x^{\mathsf T}Ay|.
```

The project's undoubled cap is

```math
H(A)=\max_z\left|\sum_{i<j}a_{ij}z_iz_j\right|
=\frac12Q(A).
```

Thus ordinary polarization is \(Q(A)\ge B(A)/2\), or
\(H(A)\ge B(A)/4\). Every displayed recoupling certificate below is in
doubled normalization; divide it by two for the project normalization.

## 2. Exact agreement/disagreement algebra

Fix a real symmetric zero-diagonal matrix \(A\) and Boolean vectors \(x,y\).
Let

```math
I=\{i:x_i=y_i\},\qquad J=\{i:x_i=-y_i\},
```

and let

```math
p=x_I=y_I,\qquad q=x_J=-y_J.
```

Write

```math
P=p^{\mathsf T}A[I]p,\qquad
R=q^{\mathsf T}A[J]q.                            \tag{2.1}
```

Equivalently, after extending \(p,q\) by zero to the full vertex set,
\(x=p+q\) and \(y=p-q\). Symmetry cancels the cross terms:

```math
\boxed{x^{\mathsf T}Ay=P-R.}                    \tag{2.2}
```

For a symmetric zero-diagonal block \(D\), define its two doubled one-sided
caps

```math
C_+(D)=\max_{r\in\{\pm1\}^{|D|}}r^{\mathsf T}Dr,
\qquad
C_-(D)=-\min_{r\in\{\pm1\}^{|D|}}r^{\mathsf T}Dr. \tag{2.3}
```

Both are nonnegative because the average Boolean energy is zero. Write
\(C_{\operatorname{sgn}t}(D)=C_+(D)\) for \(t>0\) and \(C_-(D)\)
for \(t<0\).

### Theorem 2.1 (exact two-shore recoupling)

With the notation above:

1. If \(PR\ge0\), including \(P=0\) or \(R=0\), then

   ```math
   \boxed{Q(A)\ge|x^{\mathsf T}Ay|.}             \tag{2.4}
   ```

2. If \(PR<0\), then

   ```math
   \boxed{
   Q(A)\ge
   \max\left\{
   |P|+C_{\operatorname{sgn}P}(A[J]),
   |R|+C_{\operatorname{sgn}R}(A[I])
   \right\}.}                                    \tag{2.5}
   ```

Consequently, recoupling is completely lossless whenever either

```math
C_{\operatorname{sgn}P}(A[J])\ge|R|
\quad\text{or}\quad
C_{\operatorname{sgn}R}(A[I])\ge|P|.            \tag{2.6}
```

The second condition in (2.6), for example, asks the agreement block for an
energy with the sign of the disagreement-block witness. It does not ask for
another full \(n\)-vertex optimizer.

#### Proof

Fix any \(r\in\{\pm1\}^J\). The two full common-spin vectors
\((p,r)\) and \((p,-r)\) have energies

```math
P+r^{\mathsf T}A[J]r
\ \pm\ 2p^{\mathsf T}A[I,J]r.
```

The elementary identity

```math
\max\{|a+b|,|a-b|\}=|a|+|b|
```

therefore gives

```math
Q(A)\ge
|P+r^{\mathsf T}A[J]r|
+2|p^{\mathsf T}A[I,J]r|.                       \tag{2.7}
```

Choose \(r\) whose principal-block energy has the sign of \(P\). This proves

```math
Q(A)\ge|P|+C_{\operatorname{sgn}P}(A[J]).        \tag{2.8}
```

Interchanging \(I\) and \(J\) proves the other term in (2.5). In particular,
the cross block in (2.7) is never a defect: a global sign choice on the
recoupled shore makes it helpful.

The same argument, or random completion of a partial spin, gives
\(Q(A)\ge\max\{|P|,|R|\}\). If \(PR\ge0\), then
\(|P-R|\le\max\{|P|,|R|\}\), proving (2.4). If \(PR<0\), (2.2) gives
\(|x^{\mathsf T}Ay|=|P|+|R|\), so (2.6) follows from (2.5).

### Exact defect form

When \(PR<0\), define

```math
\begin{aligned}
\delta_J(x,y)
&=\left[|R|-C_{\operatorname{sgn}P}(A[J])\right]_+,\\
\delta_I(x,y)
&=\left[|P|-C_{\operatorname{sgn}R}(A[I])\right]_+.
\end{aligned}                                    \tag{2.9}
```

Then

```math
\boxed{
Q(A)\ge |x^{\mathsf T}Ay|
-\min\{\delta_I(x,y),\delta_J(x,y)\}.}           \tag{2.10}
```

The theorem also has the clean averaged form requested in the audit:

```math
\boxed{
Q(A)\ge
\frac{|x^{\mathsf T}Ay|
+C_{\operatorname{sgn}R}(A[I])
+C_{\operatorname{sgn}P}(A[J])}{2}}
\qquad(PR<0).                                    \tag{2.11}
```

Indeed, the maximum of the two terms in (2.5) is at least their average.

## 3. Polynomial one-sided spectral certificate

The exact one-sided caps in (2.5) isolate the right obstruction, but computing
them is still a Boolean optimization. This section replaces them by a
polynomial-time statistic.

The assumptions are explicit. Let \(D\) be an \(m\times m\) real symmetric
matrix such that

```math
d_{ii}=0,\qquad |d_{ij}|\le1.                    \tag{3.1}
```

All principal blocks of a signing satisfy (3.1). For
\(\sigma\in\{+1,-1\}\), let

```math
s_\sigma(D)=\operatorname{tr}((\sigma D)_+),
\qquad
r_\sigma(D)=\operatorname{rank}((\sigma D)_+),  \tag{3.2}
```

where \((\sigma D)_+\) is the positive spectral part. Set

```math
\kappa=\frac\pi2-1
```

and define the sign-specific projector certificate

```math
\Gamma_\sigma(D)=
\max_{0\le\theta\le1}
\frac2\pi\left(
\theta s_\sigma(D)-\kappa\theta^2r_\sigma(D)
\right).                                        \tag{3.3}
```

The scalar maximizer is explicit:

```math
\theta_*=
\min\left\{1,\frac{s_\sigma(D)}
{2\kappa r_\sigma(D)}\right\}                   \tag{3.4}
```

when \(r_\sigma(D)>0\).

### Theorem 3.1 (one-sided projector rounding)

For every \(D\) satisfying (3.1),

```math
\boxed{C_\sigma(D)\ge\Gamma_\sigma(D).}          \tag{3.5}
```

Because \(\operatorname{tr}D=0\),

```math
s_+(D)=s_-(D)=\frac12\|D\|_*.
```

Taking \(\theta=1\) and \(r_\sigma(D)\le m\) gives the simpler common
nuclear certificate

```math
\boxed{
C_+(D),C_-(D)\ge
S(D):=
\left[\frac{\|D\|_*}{\pi}
-\left(1-\frac2\pi\right)m\right]_+.}           \tag{3.6}
```

Here \(\|D\|_*\) is the nuclear norm of the principal block \(D\), not the
unknown Boolean cap of the full parent.

#### Proof

It is enough to prove the positive case for \(\sigma D\), so write
\(E=\sigma D\). Let \(\Pi\) be the orthogonal projector onto the positive
eigenspace of \(E\). For \(0\le\theta\le1\), define

```math
K_\theta=
\theta\Pi+\operatorname{diag}
\left(1-\theta\Pi_{11},\ldots,1-\theta\Pi_{mm}\right). \tag{3.7}
```

Both summands are positive semidefinite, and \(K_\theta\) has diagonal one.
It is therefore a Gaussian correlation matrix. If
\(G\sim N(0,K_\theta)\) and \(z_i=\operatorname{sgn}G_i\), the Gaussian
sign identity gives, for \(i\ne j\),

```math
\mathbb E[z_iz_j]
=\frac2\pi\arcsin(\theta\Pi_{ij}).               \tag{3.8}
```

The elementary sharp remainder bound

```math
|\arcsin u-u|\le
\left(\frac\pi2-1\right)u^2
=\kappa u^2,\qquad |u|\le1,                     \tag{3.9}
```

has the following direct derivative proof. For \(0<u<1\), put

```math
r(u)=\frac{\arcsin u-u}{u^2}.
```

Then

```math
u^3r'(u)
=u\left(\frac1{\sqrt{1-u^2}}-1\right)
-2(\arcsin u-u).                                \tag{3.9a}
```

Writing \(u=\sin\vartheta\), the right side becomes

```math
F(\vartheta)=\tan\vartheta+\sin\vartheta-2\vartheta.
```

Now \(F(0)=0\), and, with \(c=\cos\vartheta\in(0,1]\),

```math
F'(\vartheta)=\sec^2\vartheta+\cos\vartheta-2
=\frac1{c^2}+c-2\ge0.                           \tag{3.9b}
```

Indeed, \(h(c)=c^{-2}+c-2\) satisfies
\(h'(c)=1-2c^{-3}<0\) on \((0,1]\) and \(h(1)=0\), so
\(h(c)\ge0\) there. Thus \(r'(u)\ge0\). Since
\(r(0+)=0\) and \(r(1)=\pi/2-1\), (3.9) follows for \(u\ge0\);
oddness of \(\arcsin u-u\) gives the absolute-value statement for
\(u<0\).

Using \(|e_{ij}|\le1\), (3.8), and (3.9),

```math
\begin{aligned}
\mathbb E[z^{\mathsf T}Ez]
&=\frac2\pi\sum_{i\ne j}e_{ij}
  \arcsin(\theta\Pi_{ij})\\
&\ge\frac2\pi\left(
\theta\operatorname{tr}(E\Pi)
-\kappa\theta^2\sum_{i\ne j}\Pi_{ij}^2
\right)\\
&\ge\frac2\pi\left(
\theta s_\sigma(D)-\kappa\theta^2r_\sigma(D)
\right),
\end{aligned}                                    \tag{3.10}
```

because

```math
\sum_{i\ne j}\Pi_{ij}^2
\le\|\Pi\|_F^2=\operatorname{rank}\Pi.
```

Some Boolean realization is at least its expectation. Optimizing \(\theta\)
proves (3.5). Equation (3.6) follows as stated.

### Corollary 3.2 (polynomial recoupling)

If \(PR<0\), Theorems 2.1 and 3.1 give

```math
\boxed{
Q(A)\ge
\max\left\{
|P|+\Gamma_{\operatorname{sgn}P}(A[J]),
|R|+\Gamma_{\operatorname{sgn}R}(A[I])
\right\}.}                                      \tag{3.11}
```

In particular,

```math
\boxed{
Q(A)\ge
\max\{|P|+S(A[J]),\,|R|+S(A[I])\}.}             \tag{3.12}
```

The averaged nuclear form is

```math
\boxed{
Q(A)\ge
\frac{|x^{\mathsf T}Ay|+S(A[I])+S(A[J])}{2}}
\qquad(PR<0).                                    \tag{3.13}
```

The polynomial sufficient condition for no-loss recoupling is

```math
S(A[J])\ge|R|
\quad\text{or}\quad
S(A[I])\ge|P|,                                  \tag{3.14}
```

with the stronger sign-specific version obtained by replacing \(S\) with
the appropriate \(\Gamma\).

All data in (3.11)--(3.14) are computable from \(A,x,y\) by two principal
eigendecompositions. This is strictly weaker than computing \(Q(A)\).

### Anchored-shore strengthening: retain the cross field

The preceding certificate throws away the nonnegative cross term in (2.7).
It can instead be retained while collapsing one whole shore to a single
weighted vertex. Fix the anchor \(p\) on \(I\), and put

```math
h_J=A[J,I]p,
\qquad
E_J(p)=
\begin{pmatrix}
A[J]&h_J\\
h_J^{\mathsf T}&0
\end{pmatrix}.                                  \tag{3.15}
```

For \(r\in\{\pm1\}^J\) and \(t\in\{\pm1\}\), the full common spin
\((tp,r)\) has energy

```math
(tp,r)^{\mathsf T}A(tp,r)
=P+(r,t)^{\mathsf T}E_J(p)(r,t).                \tag{3.16}
```

Consequently the exact recoupling theorem strengthens to

```math
\boxed{
Q(A)\ge |P|+C_{\operatorname{sgn}P}(E_J(p)).}  \tag{3.17}
```

Equation (3.17) is asserted for `P != 0`.  If `P=0`, there is no aligned
sign to choose; the valid replacement is simply

```math
Q(A)\ge \max\{C_+(E_J(p)),C_-(E_J(p))\}.
```

There is a symmetric bound anchored on \(q\):

```math
E_I(q)=
\begin{pmatrix}
A[I]&A[I,J]q\\
q^{\mathsf T}A[J,I]&0
\end{pmatrix},
\qquad
Q(A)\ge |R|+C_{\operatorname{sgn}R}(E_I(q)).   \tag{3.18}
```

Likewise, (3.18) is for `R != 0`, with the maximum of the two one-sided
caps used when `R=0`.

Here \(C_\sigma(E)\) is the one-sided Boolean cap of the weighted matrix
\(E\). Equations (3.17)--(3.18) remain strictly weaker than maximizing over
the full parent: they optimize only \(|J|+1\) or \(|I|+1\) spins, with
the anchor shore compressed to one coordinate. They dominate (2.5), since
fixing the new coordinate and averaging over its two signs recovers the
shore-only cap.

The augmented entries \(h_J\) need not have magnitude at most one. The
projector theorem therefore needs a weighted remainder. For an arbitrary
real symmetric zero-diagonal matrix \(E\), let \(\Pi_\sigma\) project
onto the positive eigenspace of \(\sigma E\), and define

```math
\mathcal W_E(\Pi_\sigma)
=\sum_{a\ne b}|e_{ab}|(\Pi_\sigma)_{ab}^2.    \tag{3.19}
```

Exactly the proof of Theorem 3.1, without replacing \(|e_{ab}|\) by one,
gives

```math
\boxed{
C_\sigma(E)\ge
\widetilde\Gamma_\sigma(E):=
\max_{0\le\theta\le1}\frac2\pi
\left[
\theta\operatorname{tr}((\sigma E)_+)
-\kappa\theta^2\mathcal W_E(\Pi_\sigma)
\right].}                                      \tag{3.20}
```

If \(\mathcal W_E(\Pi_\sigma)>0\), its maximizer is

```math
\theta_*=
\min\left\{
1,\frac{\operatorname{tr}((\sigma E)_+)}
{2\kappa\mathcal W_E(\Pi_\sigma)}
\right\}.                                     \tag{3.21}
```

No entry bound is assumed in (3.19)--(3.21). Thus

```math
\boxed{
Q(A)\ge\max\left\{
|P|+\widetilde\Gamma_{\operatorname{sgn}P}(E_J(p)),
|R|+\widetilde\Gamma_{\operatorname{sgn}R}(E_I(q))
\right\}.}                                    \tag{3.22}
```

When \(PR<0\), define the anchored defect

```math
\widetilde\Delta_\Gamma(x,y)=
\min\left\{
\left[|R|-\widetilde\Gamma_{\operatorname{sgn}P}(E_J(p))\right]_+,
\left[|P|-\widetilde\Gamma_{\operatorname{sgn}R}(E_I(q))\right]_+
\right\}.                                     \tag{3.23}
```

Then

```math
Q(A)\ge |x^{\mathsf T}Ay|-\widetilde\Delta_\Gamma(x,y). \tag{3.24}
```

This is the strongest polynomial recoupling interface proved here. It uses
the principal-shore spectrum and the complete anchored cross field, but
never optimizes a full parent spin.

### Interface with a random asymmetric scheme

The inequalities are samplewise. For any jointly distributed asymmetric
outputs \((X,Y)\), define the polynomial defect to be zero when \(PR\ge0\),
and when \(PR<0\) define

```math
\Delta_\Gamma(X,Y)=
\min\left\{
\left[|R|-\Gamma_{\operatorname{sgn}P}(A[J])\right]_+,
\left[|P|-\Gamma_{\operatorname{sgn}R}(A[I])\right]_+
\right\}.                                       \tag{3.25}
```

Then

```math
\boxed{
Q(A)\ge
\left|\mathbb E[X^{\mathsf T}AY]\right|
-\mathbb E\Delta_\Gamma(X,Y).}                  \tag{3.26}
```

This is the minimal rigorous interface for importing a genuinely asymmetric
Krivine output. A large bilinear expectation is useful without polarization
exactly when the *joint* agreement/disagreement law makes the defect in
(3.25), or the stronger anchored defect (3.23), subleading. No separate
estimate on \(X^{\mathsf T}AX\) or
\(Y^{\mathsf T}AY\) appears.

## 4. Conference/flat specialization

If \(D\) is a full sign block of order \(m\), then

```math
\|D\|_*
\ge\frac{\|D\|_F^2}{\|D\|_{\mathrm{op}}}
=\frac{m(m-1)}{\|D\|_{\mathrm{op}}}.             \tag{4.1}
```

Thus

```math
C_+(D),C_-(D)\ge
\left[
\frac{m(m-1)}{\pi\|D\|_{\mathrm{op}}}
-\left(1-\frac2\pi\right)m
\right]_+.                                      \tag{4.2}
```

Let \(C\) be a symmetric conference signing of order \(n\), so

```math
C^2=(n-1)I.
```

Every principal compression satisfies
\(\|C[S]\|_{\mathrm{op}}\le\sqrt{n-1}\). Hence, for \(s=|S|\),

```math
\boxed{
C_+(C[S]),C_-(C[S])\ge
F_n(s):=
\left[
\frac{s(s-1)}{\pi\sqrt{n-1}}
-\left(1-\frac2\pi\right)s
\right]_+.}                                     \tag{4.3}
```

For an asymmetric pair on \(C\) with \(PR<0\), \(i=|I|\), \(j=|J|\),

```math
\boxed{
Q(C)\ge\max\{|P|+F_n(j),\,|R|+F_n(i)\}.}         \tag{4.4}
```

Equivalently,

```math
Q(C)\ge
\frac{|x^{\mathsf T}Cy|+F_n(i)+F_n(j)}2.         \tag{4.5}
```

If \(i=\alpha n+o(n)\) and \(j=(1-\alpha)n+o(n)\), the extra term in
(4.5) is

```math
\frac{\alpha^2+(1-\alpha)^2}{2\pi}n^{3/2}
+o(n^{3/2}).                                    \tag{4.6}
```

This is a genuine project-scale contribution. It does not prove a better
constant by itself, because the asymmetric scheme must simultaneously
control \(P,R\) and its bilinear response. It does show that the
same-spin conversion need not pay the universal factor two on flat inputs.

The cross-block row fields need no centering hypothesis in this theorem:
equation (2.7) makes them nonnegative after choosing the global sign of the
recoupled shore. The only remaining hypothesis is a one-sided spectral
condition on a diagonal shore.

## 5. Why a structure-free phase recoupling is impossible

### 5.1 Positive-semidefinite obstruction to sign reversal

After switching by \(q\), reversing the sign of every pair interaction on a
disagreement shore of size \(m\) would require a random sign vector \(r\)
with

```math
\mathbb E[r_ir_j]\le-\gamma\qquad(i\ne j)
```

for a constant \(\gamma>0\). But

```math
0\le\mathbb E\left(\sum_i r_i\right)^2
\le m-\gamma m(m-1),
```

so

```math
\boxed{\gamma\le\frac1{m-1}.}                   \tag{5.1}
```

For even \(m\), equality is attained by a uniform balanced sign vector.
Thus a matrix-independent correlated phase rule can reverse the disagreement
block only at a vanishing \(1/m\) coefficient. A constant-scale theorem must
use a statistic of the actual block \(A[J]\); (3.2) supplies one.

### 5.2 Four-phase encoding does not evade the obstruction

The complex vector

```math
w=p+iq
```

satisfies

```math
w^{\mathsf T}Aw
=P-R+2i\,p^{\mathsf T}A[I,J]q.                  \tag{5.2}
```

So the desired bilinear value is the real part of a four-phase quadratic
form. However, any deterministic antipodal coordinate map

```math
f:\{\pm1\}^2\to\{\pm1\},
\qquad f(-a,-b)=-f(a,b),
```

is one of \(\pm x\) or \(\pm y\): its two values at \((1,1)\) and
\((1,-1)\) determine it, and they either agree or disagree. Random global
phase halfspace maps are mixtures of exactly these four maps. They therefore
see only the original full energies, not \(P-R\).

Non-antipodal four-state maps require additional block statistics and cannot
be guaranteed from \(x^{\mathsf T}Ay\) alone. The one-sided shore statistic
in Theorem 2.1 is the exact missing information.

### 5.3 Ordinary polarization is asymptotically sharp without flatness

There is also a full-sign-matrix obstruction, not merely a correlation-matrix
argument. Let \(H_d=J_d-I_d\), choose a \(d\times d\) sign matrix \(R_d\)
with

```math
\|R_d\|_{\infty\to1}=O(d^{3/2}),
```

and form the symmetric zero-diagonal signing

```math
A_d=
\begin{pmatrix}
H_d&R_d\\
R_d^{\mathsf T}&-H_d
\end{pmatrix}.                                  \tag{5.3}
```

Such \(R_d\) exist by a random-sign union bound. With

```math
x=(\mathbf1,\mathbf1),\qquad
y=(\mathbf1,-\mathbf1),
```

the two cross terms cancel and

```math
x^{\mathsf T}A_dy=2d(d-1).                      \tag{5.4}
```

For any common spin \((u,v)\),

```math
(u,v)^{\mathsf T}A_d(u,v)
=\left(\sum_i u_i\right)^2
-\left(\sum_i v_i\right)^2
+2u^{\mathsf T}R_dv,
```

and therefore

```math
Q(A_d)\le d^2+O(d^{3/2}).                        \tag{5.5}
```

Together with the universal \(B(A)\le2Q(A)\), (5.4)--(5.5) imply

```math
\frac{B(A_d)}{Q(A_d)}\longrightarrow2.           \tag{5.6}
```

These matrices have \(Q(A_d)=\Theta(n^2)\), not the near-minimal
\(n^{3/2}\) scale. They prove that universal lossless recoupling is false
and that an extra flatness/statistical hypothesis is logically necessary.
The conference inequality (4.4) is the corresponding correct-scale theorem.

## 6. Reproducible finite audit

The checker is
[check_same_spin_recoupling.py](../computations/check_same_spin_recoupling.py).
Run it in the repository environment:

```text
.venv/bin/python computations/check_same_spin_recoupling.py
```

It exhausts every bilinear maximizing pair, including all choices at zero
row fields, for:

- every saved exact minimizer orbit representative of orders \(3\) through
  \(8\);
- the saved exact order-\(10\) minimizer; and
- saved symmetric conference matrices of orders \(6,10,14,18\).

The columns below all use doubled normalization. “Exact face” uses the actual
one-sided caps in (2.5). “Projector” uses (3.3). “Nuclear” uses the coarser
(3.6). The last column is the agreement/disagreement decomposition for a
bilinear maximizing pair producing the exact-face entry.

| case | \(n\) | \(Q\) | \(B\) | \(B/2\) | exact face | projector | nuclear | \((|I|,|J|,P,R)\) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| min-n3-orbit0 | 3 | 6 | 6 | 3.000 | 6.000 | 6.000 | 6.000 | (3, 0, 6, 0) |
| min-n4-orbit0 | 4 | 8 | 8 | 4.000 | 8.000 | 8.000 | 8.000 | (4, 0, 8, 0) |
| min-n5-orbit0 | 5 | 8 | 8 | 4.000 | 8.000 | 8.000 | 8.000 | (5, 0, 8, 0) |
| min-n6-orbit0 | 6 | 10 | 12 | 6.000 | 8.000 | 6.558 | 6.183 | (3, 3, 6, -6) |
| min-n7-orbit0 | 7 | 18 | 18 | 9.000 | 18.000 | 18.000 | 18.000 | (2, 5, 2, -16) |
| min-n7-orbit1 | 7 | 18 | 18 | 9.000 | 18.000 | 18.000 | 18.000 | (7, 0, 18, 0) |
| min-n7-orbit2 | 7 | 18 | 18 | 9.000 | 18.000 | 18.000 | 18.000 | (1, 6, 0, -18) |
| min-n8-orbit0 | 8 | 20 | 24 | 12.000 | 16.000 | 12.836 | 12.456 | (4, 4, 12, -12) |
| min-n8-orbit1 | 8 | 20 | 24 | 12.000 | 16.000 | 12.836 | 12.456 | (4, 4, 12, -12) |
| exact-n10 | 10 | 26 | 40 | 20.000 | 24.000 | 21.115 | 20.730 | (5, 5, 20, -20) |
| conference-n6 | 6 | 10 | 12 | 6.000 | 8.000 | 6.558 | 6.183 | (3, 3, 6, -6) |
| conference-n10 | 10 | 30 | 30 | 15.000 | 30.000 | 30.000 | 30.000 | (10, 0, 30, 0) |
| conference-n14 | 14 | 42 | 46 | 23.000 | 42.000 | 31.093 | 29.639 | (8, 6, 28, -18) |
| conference-n18 | 18 | 66 | 72 | 36.000 | 64.000 | 41.879 | 40.062 | (9, 9, 36, -36) |

The checker asserts for every tested pair

```math
\text{nuclear}\le\text{projector}\le
\text{exact face}\le Q(A).
```

The polynomial certificates strictly exceed ordinary \(B/2\) on each listed
nontrivial \(B>Q\) representative. This is finite verification of a real
algebraic escape, not evidence for a uniform asymptotic improvement.

## 7. Precise next theorem

Given a genuinely asymmetric Krivine construction \((X,Y)\), the remaining
target is no longer “convert bilinear to quadratic.” It is the following
joint-law statement:

> Prove uniformly for every near-minimizing signing that the defect
> \(\mathbb E\Delta_\Gamma(X,Y)\) in (3.25) is \(o(n^{3/2})\), or is smaller
> than the gain of the asymmetric bilinear response over the current doubled
> constant \(c_*\).

The polynomial version is strictly weaker than computing \(Q(A)\): it uses
only the two random shores, their witnessed energies, and two polynomial
spectral profiles.  The exact augmented version retains, for each branch,
the free principal block `D` as well as `(P,h)`; it is a restricted Boolean
maximum rather than a bounded-complexity statistic.  Computing both branches
costs `O(2^max(|I|,|J|))`, while deliberately using only the branch with the
smaller free shore costs `O(2^min(|I|,|J|))` and yields a weaker valid bound.
Either formulation is falsifiable on finite and conference data. If it
fails, the failure identifies which induced shore has a spectrally invisible
one-sided energy, rather than returning the undifferentiated polarization
factor two.
