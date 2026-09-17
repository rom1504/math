# Wave 29 Route 2: joint aligned payoff and low row square

## Result

There is no generic scalar-to-joint transfer, even at the target exponents
and even if the missing scalar lower tail is granted. An abstract quadratic
pair satisfying all coarse constraints currently known for \(A^2\) has
aligned mass \(\exp\{-O(n^{3/4-c})\}\) but zero aligned mass below a
\(\Theta(n^{9/4-c})\) row cap. An exact finite \(A_8\) selector exhibits
the same positive-correlation mechanism in the actual
\(H_S=A\circ(\xi_S\xi_S^{\mathsf T}-p_2\mathbf1\mathbf1^{\mathsf T})\)
geometry.

The useful extra hypothesis is a **conditional row-Laplace lower bound**
under the aligned event. Equivalently, it is an integrated bound on the
positive aligned-event/row-square covariance along a row tilt. This gives
an exact transfer without subtracting unconditional row-bad mass.

## 1. Exact row-tilt transfer theorem

Let \(R\ge0\) on a finite probability space and let \(E\) have probability
\(p>0\). In the application,

```math
E=E_T=\{|x^{\mathsf T}H_Sx|\ge T\},
\qquad R=x^{\mathsf T}A^2x.
```

For \(\theta>0\), define

```math
Z(\theta)=\mathbb E e^{-\theta R},\qquad
d\mathbb P_\theta=Z(\theta)^{-1}e^{-\theta R}\,d\mathbb P,
\qquad p_\theta=\mathbb P_\theta(E).
```

The conditional Laplace transform factors exactly as

```math
L_E(\theta):=\mathbb E[e^{-\theta R}\mid E]
=Z(\theta)\frac{p_\theta}{p}.
```

If \(r_C=\mathbb P(R\le C\mid E)\), splitting at the cap gives

```math
L_E(\theta)\le r_C+(1-r_C)e^{-\theta C}.
```

Therefore

```math
\boxed{
\mathbb P(E,\ R\le C)
\ge
p\,\frac{L_E(\theta)-e^{-\theta C}}
{1-e^{-\theta C}}.
}
\tag{JT.1}
```

This uses conditional row mass and is not the failed unconditional Markov
subtraction.

There is an exact derivative identity. Put

```math
\kappa(\theta)
=\int_0^\theta
\{\mathbb E_s[R\mid E]-\mathbb E_sR\}\,ds.
```

Then

```math
\boxed{
\log\frac{p_\theta}{p}=-\kappa(\theta),\qquad
\frac{d}{d\theta}\log p_\theta
=\mathbb E_\theta R-\mathbb E_\theta[R\mid E]
=-\frac{\operatorname{Cov}_\theta(\mathbf1_E,R)}{p_\theta}.
}
\tag{JT.2}
```

Writing \(\mu=\mathbb ER\), Jensen gives
\(Z(\theta)\ge e^{-\theta\mu}\). Consequently, if

```math
\theta(C-\mu)-\kappa(\theta)\ge\log2,
\tag{JT.3}
```

then (JT.1) implies

```math
\boxed{
\mathbb P(E,\ R\le C)
\ge \frac p2\exp\{-\theta\mu-\kappa(\theta)\}.
}
\tag{JT.4}
```

Thus, with \(L=n^{3/4-c}\), a **still unproved** scalar lower bound
\(p\ge e^{-K_HL}\) transfers to the same exponent if some \(\theta\)
satisfies (JT.3) and

```math
\theta\mu+\kappa(\theta)=O(L).
\tag{JT.5}
```

Equivalently, the weakest direct sufficient bound on the tilted hit ratio is

```math
\boxed{
\log\frac{p_\theta}{p}
\ge-\theta(C-\mu)+\log2,
}
\tag{JT.6}
```

together with (JT.5). For a signing \(\mu=n(n-1)\). At the target scales,
take

```math
\theta=\frac{L}{C}=\Theta(n^{-3/2}),\qquad
C=\Theta(n^{9/4-c}).
```

Then \(\theta C=L\) and
\(\theta\mu=\Theta(n^{1/2})=o(L)\). For example, any fixed-margin bound

```math
\kappa(\theta)\le(1-\varepsilon)L
```

implies joint mass \(\exp\{-O(L)\}\), conditional on the unproved scalar
mass \(p\ge e^{-O(L)}\). More generally one can multiply \(\theta\) by a
fixed constant to absorb the corresponding uniform constant in
\(\kappa/L\).

The transparent stronger condition

```math
\operatorname{Cov}_s(\mathbf1_E,R)\le0
\quad(0\le s\le1/\mu)
```

has \(\kappa(1/\mu)\le0\). Since
\(C/\mu=\Theta(n^{1/4-c})\), it retains a constant fraction of the scalar
tail. Pointwise negative covariance is not needed; (JT.3)--(JT.5) are the
actual integrated conditions.

## 2. Abstract target-scale counterexample

Fix \(0<c<1/4\), let \(L=n^{3/4-c}\), and take an unbounded sequence of
fourth-power orders. Let \(u=\mathbf1\), \(a=n^{-1/4}\), and
\(b=n^{1/2}\). Define

```math
H=a(uu^{\mathsf T}-I),\qquad
B=(n-1-b)I+buu^{\mathsf T}
=(n-1)I+\frac baH.
```

The matrix \(H\) is symmetric and zero diagonal. The row kernel \(B\)
obeys exactly

```math
B\succeq0,\qquad B_{ii}=n-1,\qquad
\operatorname{tr}B=n(n-1),\qquad
\lVert B\rVert_{\rm op}=\Theta(n^{3/2}).
```

Also \(\lVert H\rVert_{\rm op}=\Theta(n^{3/4})\) and
\(\lVert H\rVert_F^2=\Theta(n^{3/2})\). These match the coarse scales
available for an exact minimizer, although this abstract \(B\) is not
asserted to equal the square of a complete signing.

For uniform \(x\in\{\pm1\}^n\), put \(Z=u^{\mathsf T}x\). Then

```math
x^{\mathsf T}Hx=a(Z^2-n),\qquad
R=x^{\mathsf T}Bx
=n(n-1)+b(Z^2-n)
=n(n-1)+n^{3/4}x^{\mathsf T}Hx.
\tag{JT.7}
```

Choose an integer \(r\), of the parity of \(n\), with
\(r=2\sqrt{nL}+O(1)\), and set

```math
T=a(r^2-n)=\Theta(n^{3/2-c}),\qquad
C=n(n-1)+b(r^2-n)-1=\Theta(n^{9/4-c}).
```

For large \(n\), \(T>an\), so the negative side of \(a(Z^2-n)\) cannot
reach absolute threshold \(T\). Hence

```math
\{|x^{\mathsf T}Hx|\ge T\}=\{|Z|\ge r\}.
```

An elementary Stirling bound for one binomial atom gives

```math
\Pr\{|Z|\ge r\}
\ge2^{-n}{n\choose(n+r)/2}
\ge\exp\{-O(r^2/n)-O(\log n)\}
=\exp\{-O(L)\}.
```

But (JT.7) gives \(R>C\) everywhere on this event. Thus

```math
\boxed{\Pr\{|x^{\mathsf T}Hx|\ge T,\ R\le C\}=0.}
\tag{JT.8}
```

No theorem using only positivity, diagonal, trace, operator norm, and a
separate scalar lower tail can prove joint retention. A deterministic repair
preserving the same threshold is impossible in this example.

## 3. Exact finite signing counterexample

Use the audited exact order-eight minimizer \(A_8\), take \(m=4\), and use
the zero-based selector

```text
S = {2,3,5,7}.
```

Here

```math
q_8=20,\qquad Q(A_8[S])=12,\qquad p_2=\frac3{14},
\qquad T=Y_A(S)=12-5\sqrt2.
```

For

```math
H_S=A_8\circ
\left(\xi_S\xi_S^{\mathsf T}-\frac3{14}\mathbf1\mathbf1^{\mathsf T}\right),
\qquad R=x^{\mathsf T}A_8^2x,
```

exact enumeration of the 128 projective spins gives

```math
\begin{array}{c|ccc}
R\text{ on the aligned event}&56&64&72\\ \hline
\#\text{ spins}&8&12&6.
\end{array}
```

The largest non-hit payoff is \(34/7\), the smallest hit payoff is \(40/7\),
and

```math
\frac{34}{7}<12-5\sqrt2<\frac{40}{7}.
```

Therefore

```math
\Pr(E_T)=\frac{26}{128}=\frac{13}{64},\qquad
\Pr(R\le40)=\frac{30}{128}=\frac{15}{64},\qquad
\boxed{\Pr(E_T,\ R\le40)=0.}
\tag{JT.9}
```

The obstruction is visible in the row-tilt derivative:

```math
\mathbb ER=56,\qquad
\mathbb E[R\mid E_T]=\frac{824}{13},\qquad
\left.\frac d{d\theta}\log p_\theta\right|_{\theta=0}
=-\frac{96}{13}<0.
\tag{JT.10}
```

Thus exact signing minimality does not force even local negative
row/alignment correlation. This finite cap lies below the asymptotic target
ratio \(C/\mathbb ER\to\infty\), so (JT.9) is a mechanism wall, not an
asymptotic falsifier of (10.795). The abstract construction is the
target-scale falsifier to arguments using only coarse matrix constraints.

The census is checked by:

```bash
.venv/bin/python tmp/joint_tail_r29.py
```

## 4. Open minimizer-specific input

This route does **not** prove the scalar lower tail. Even if that separate
bound were proved, (JT.8) shows that joint retention needs extra structure.
The most economical exact extra input is, for some target-scale \(\theta\),

```math
\boxed{
\mathbb E[e^{-\theta R_2(x)}\mid
|x^{\mathsf T}H_Sx|\ge T]
\ge e^{-K_Rn^{3/4-c}},
\qquad
\theta C\ge K_Rn^{3/4-c}+\log2.
}
\tag{JT.11}
```

Together with scalar mass \(e^{-K_Hn^{3/4-c}}\), (JT.1) gives joint mass
at least \(\frac12e^{-(K_H+K_R)n^{3/4-c}}\).

This formulation is also necessary up to constants at the exponent scale.
Indeed, if
\(\mathbb P(R\le C\mid E)\ge e^{-K n^{3/4-c}}\), then

```math
L_E(\theta)\ge
e^{-\theta C}\mathbb P(R\le C\mid E).
```

For \(\theta C=O(n^{3/4-c})\), this is (JT.11) with a changed constant.
Thus the conditional row-Laplace criterion neither hides an unconditional
subtraction nor imposes a stronger exponential-scale target.

At the proof-mechanism level, (JT.11) should be attacked through the
integrated covariance in (JT.2), not unconditional row pruning. The finite
\(A_8\) sign reversal means a viable theorem must allow some positive
covariance and control only its integral/scale. No such minimizer-specific
bound follows from the current ledger, and no asymptotic complete-signing
counterexample to it is known.
