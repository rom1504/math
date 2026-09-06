# Ordinary good-signing entropy: two exact structures and the order gap

Date: 2026-09-05. Status: new same-order entropy inequalities and an exact
monotone-family reduction are proved below. No all-order entropy limit or
original optimum limit is proved. The failed imports are scoped precisely;
they do not rule out correlated-disorder entropy methods.

Write N=binom(n,2), H_A(x)=sum_(i<j) a_ij x_i x_j,
q(A)=max_x |H_A(x)|, G_n(c)={A:q(A)<=c n^(3/2)}, and
Z_n(c)=|G_n(c)|. All logarithms in the formulas below are natural.

## 1. An entropy inequality for the entire good-signing family

For a nonempty G_n(c), set

```math
s_n(c)=\frac{\log Z_n(c)}{N}\in[0,\log2].
```

Let h(u)=-u log u-(1-u)log(1-u), with h^(-1) taking values in [0,1/2].
Fix p in (0,1/2), and choose
t>2 sqrt(p(1-p) log 2). Define

```math
c'=(1-2p)c+t,
\qquad
\Psi_p(s)=h\bigl(p+(1-2p)h^{-1}(s)\bigr).
```

Then uniformly over every nonempty cap G_n(c),

```math
\boxed{\quad s_n(c')\ge \Psi_p(s_n(c))-o(1).\quad}
```

This strengthens the one-seed estimate s_n(c')>=h(p)-o(1): it retains
the entropy of **all** good signings at the input.

### Proof, with the uniform finite-n error

Take A uniformly from G_n(c). Independently reverse each edge with
probability p to obtain B. Conditional on each A, Bernstein gives

```math
\Pr\{q(B)>c'n^{3/2}\mid A\}
\le r_n:=2^n\exp\left\{-\frac{t^2n^3}
 {8p(1-p)N+(4/3)t n^{3/2}}\right\}.
```

For our fixed p,t, r_n decays exponentially in n. Thus the same bound
holds after averaging over the arbitrary correlated law of A.
The binary-noise entropy inequality gives

```math
H(B)\ge N\Psi_p\bigl(H(A)/N\bigr)
=N\Psi_p(s_n(c)).
```

If epsilon_n=Pr(B not in G_n(c')), splitting its entropy according to
that event yields

```math
H(B)\le h(\epsilon_n)+(1-\epsilon_n)\log Z_n(c')
+\epsilon_n N\log2.
```

For sufficiently large n, r_n<=1/2, and consequently

```math
s_n(c')\ge\Psi_p(s_n(c))
-r_n\log2-\frac{h(r_n)}N.
```

In particular the error is independent of the seed family and of c.
Nonemptiness is essential: substituting s=0 for an **empty** input would
incorrectly manufacture signings.

### Primary theorem actually imported, including its stronger form

Alex Samorodnitsky, *On the entropy of a noisy function*,
[arXiv:1508.01464v4](https://arxiv.org/pdf/1508.01464), Theorem 1.12,
equation (6), proves the following stronger inequality. Put
lambda=(1-2p)^2 and let T include each input coordinate independently
with probability lambda. For an arbitrary binary random vector A and
independent p-noise B,

```math
H(B)\ge N h\left(p+(1-2p)h^{-1}
 \left(\frac{\mathbb E_T H(A_T)}{\lambda N}\right)\right).
```

The paper uses bits; the displayed natural-log version is obtained by
multiplication by log 2. There is no independence assumption on A.
Theorem 1.12 is the author's own strengthening; ordinary Mrs. Gerber's
lemma is also stated in its Theorem 1.6. The latter was originally proved
by Wyner and Ziv, *A Theorem on the Entropy of Certain Binary Sequences and
Applications: Part I*, IEEE Transactions on Information Theory 19 (1973),
769--772, [DOI](https://doi.org/10.1109/TIT.1973.1055107).

For completeness, the needed projection inequality
E_T H(A_T)>=lambda H(A) follows directly from the entropy chain rule:
condition each retained coordinate on only the earlier retained coordinates,
which gives at least the entropy conditional on all earlier coordinates,
then average its retention indicator. Monotonicity of Psi_p gives the
simpler boxed inequality above. Thus the application also yields the
stronger projected-entropy version by exactly the same conditioning proof.

### Consequence and exact missing state

For p>0, Psi_p(s)>s whenever 0<=s<log 2. Given c<d, p can be chosen
small enough that (1-2p)c+2 sqrt(p(1-p) log 2)<d. Therefore any
pointwise subsequential entropy profile has strictly increasing values
between c and d whenever its value at c is positive and below log 2.
The same statement applies at a zero-entropy cap if its events are known
to be eventually nonempty. There can be no positive, non-full plateau
in such a profile.

This constrains genuine entropy profiles but does **not** identify their
values across subsequences. The strengthened theorem's missing state is
the entropy of random **edge projections** A_T. These are not uniformly
distributed smaller-order good signings, and their tests remain tied to
n vertices. Replacing them by Z_m at some effective m is not justified
by Shearer, exchangeability, or the noise theorem.

## 2. Exact reduction to one monotone family of signings

Define a coordinatewise increasing family F_n by

```math
F_n=\left\{B:\
\sum_{e\in\delta(S)}b_e\ge0\ \hbox{and}\
\sum_{e\notin\delta(S)}b_e\ge0
\quad\hbox{for every }S\subseteq[n]\right\}.
```

Both are majority constraints on a fixed subset of edge coordinates;
thus changing any -1 edge to +1 preserves membership. Let
T(B)=sum_e b_e. Since

```math
T(B)-H_B(x)=2\sum_{e\in\delta(S_x)}b_e,
\qquad
T(B)+H_B(x)=2\sum_{e\notin\delta(S_x)}b_e,
```

membership in F_n is equivalent to q(B)=T(B): the all-ones spin
attains the positive absolute maximum.

The augmented switching group consists of vertex switches and global
coefficient complementation. For n>=3 its action is free with orbit size
2^n. Every orbit has at least one representative in F_n: switch an
absolute maximizing spin to all ones and choose the coefficient sign
so its value is positive. Define

```math
L_n(c)=\#\{B\in F_n:T(B)\le cn^{3/2}\}.
```

The orbit count gives

```math
\boxed{\quad L_n(c)\le Z_n(c)\le 2^n L_n(c).\quad}
```

The inequalities also hold for empty events. Hence their nonnegative
ordinary entropies differ by at most n log 2/n^2. Our entropy question
is exactly the leading entropy of this monotone cut/co-cut-majority family
at a positive-edge excess of order n^(3/2).

There is an analogous pressure identity. For beta>0, put

```math
\mathcal P_n(\beta)=\sum_A e^{-\beta\sqrt n\,q(A)},
\qquad
\mathcal L_n(\beta)=\sum_{B\in F_n}e^{-\beta\sqrt n\,T(B)}.
```

Then Lcal_n<=Pcal_n<=2^n Lcal_n. Equivalently,

```math
\mathcal L_n(\beta)
=(2\cosh(\beta\sqrt n))^N\Pr_{p_n}(F_n),
\qquad p_n=(1+e^{2\beta\sqrt n})^{-1}.
```

This is a scalar external-field representation retaining all original
constraints. It is not an ordinary fixed-parameter dense random-graph
partition function: the displayed Bernoulli bias is exponentially small
in sqrt(n), and the leading order n^(5/2) terms must cancel before the
desired order n^2 pressure is visible. No dimensional superadditivity
for this particular F_n is established here.

### A small exact check on potential combinatorial imports

Writing negative edges as a down-set, F_5 allows every set of at most
two negative edges, while its three-negative-edge members are exactly
the 60 labeled copies of P4. It is not a matroid: the bases
{12,23,34} and {12,24,45} fail basis exchange upon removing 34 from
the first. Adding 24 produces a star, and adding 45 produces P3 plus K2;
both are forbidden. Thus a matroid-base/Lorentzian support theorem
cannot be imported directly.

Independent exhaustive enumeration gave the following positive-edge
layer counts (only nonzero layers shown):

```text
n=3: k=3:1
n=4: k=5:6, 6:1
n=5: k=7:60, 8:45, 9:10, 10:1
n=6: k=10:72, 11:360, 12:395, 13:105, 14:15, 15:1
n=7: k=15:7140, 16:14742, 17:5880, 18:1330, 19:210, 20:21, 21:1
```

They are logconcave and ultra-logconcave in these small cases. No general
logconcavity theorem, and no cross-order consequence, is asserted.

## 3. Exact audit of the standard nonlinear-LDP bound

Let f_n(A)=-beta sqrt(n) q(A), with the N edge signs as coordinates.
Ronen Eldan, *Gaussian-width gradient complexity, reverse log-Sobolev
inequalities and nonlinear large deviations*,
[arXiv:1612.04346v6](https://arxiv.org/pdf/1612.04346), Corollary 2,
bounds the product-measure variational error by

```math
64\operatorname{Lip}(f_n)^{2/3}D(f_n)^{1/3}N^{2/3},
```

where the gradient is the discrete half-difference, and D is its Gaussian
width (Definitions (2)--(3)). Here Lip(f_n)=beta sqrt(n). At a rank-one
signing a_ij=x_i x_j, the unique projective absolute maximizing spin
remains maximizing after any single-edge flip for all sufficiently large
n. Therefore the discrete-gradient set contains
{-beta sqrt(n) v_x:x in {+/-1}^n}, and coefficient complementation gives
the opposite set. It follows that

```math
D(f_n)\ge\beta\sqrt n\,\mathbb E q(G)
\ge\left(\frac{\beta}{2\sqrt\pi}+o(1)\right)n^2,
```

where G has independent standard Gaussian edge entries. The last inequality
needs no spin-glass theorem: fix spins on half the vertices and optimize
the other half against their independent cross-edge fields; the internal
edge contribution has mean zero. The upper Gaussian union bound gives
E q(G)=O(n^(3/2)). Thus even the explicit lower bound on the theorem's
error expression is of order beta n^(7/3), not o(n^2). The small number
exp(O(n)) of cut normals alone is not a sufficient complexity estimate.

The product variational expression in Eldan's theorem uses
E_product[f_n(A)], **not** f_n(EA). These must not be confused.
The latter, naive mean-matrix surrogate is demonstrably wrong here:

```math
\sup_{m\in[-1,1]^N}
\left[-\beta\sqrt n\,q(m)+\sum_e h((1+m_e)/2)\right]
=N\log2,
```

whereas log Pcal_n<=N log2-beta sqrt(n) M_n, a fixed order-n^2 gap
using any positive lower bound on M_n/n^(3/2). This rules out that
particular surrogate, not all correlated or product-expectation
variational principles. Clipping or localization could change the global
gradient complexity and would require a fresh proof; no no-go claim for
those modifications is made.

## 4. What remains after the independent inequality attack

The true row identity is

```math
q\begin{pmatrix}A&y\\y^T&0\end{pmatrix}
=\max_x\bigl(|H_A(x)|+|y\cdot x|\bigr).
```

It requires the locations and joint geometry of the seed's energy levels,
not just its scalar norm or entropy. The archived exact cavity audit
already records this obstruction. Block gluing similarly leaves a leading
cross-edge contribution; neither monotonicity of F_n nor the new entropy
noise inequality removes it.

A precise sufficient next theorem would be entropy-preserving size
transfer: for every c and epsilon>0, show

```math
\liminf_{n\to\infty}S_n(c+\varepsilon)
\ge\limsup_{m\to\infty}S_m(c),
\qquad S_n(c)=n^{-2}\log(1+Z_n(c)).
```

This statement compares entropy, not merely a minimum or one fixed tensor
seed. Together with uniform positive-entropy thickening it forces the
original minimum to converge: a threshold c_0 between a putative liminf
and limsup has good seeds infinitely often; thicken below another threshold
c_1, giving positive limsup entropy there; transfer to a slightly larger
threshold still below the limsup, contradicting infinitely many empty
events. The inequality has not been proved. A direct all-order monotone-
family count theorem for F_n, or a correlated-disorder pressure theorem
at the displayed n-dependent field, would be alternative landings.

For a short scope check only: scalar deletion, symmetry, edge Lipschitzness,
and even the exact biased-noise cloud do not by themselves force an
all-order limit. On the same sign cubes, max{q(A),b_n n^(3/2)} with
b_n=.60+.02 sin(log n) obeys those scalar properties while its minimum
oscillates (the original all-order upper bound is <=.5+o(1)). The floor
b_n n^(3/2) increases with n, and the noise cap stays above it because
b_n<sqrt(log2). This artificial example is not the original quadratic
norm and is not a barrier to using its additional cut-code structure.
