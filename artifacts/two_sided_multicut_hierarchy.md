# Two-sided multicut hierarchies: exact identities and entropy barriers

## 1. Setup

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad a_{ij}\in\{\pm1\},
\]

and switch and globally negate \(A\), if necessary, so that

\[
H_A(\mathbf 1)=M=\max_x|H_A(x)|.
\]

For the cut represented by \(x\), put

\[
C_A(x)=\frac{M-H_A(x)}2.
\]

The full ground-state condition is exactly

\[
\boxed{0\le C_A(x)\le M\quad\text{for every }x\in\{\pm1\}^n.}
\tag{1.1}
\]

The point of this note is to retain both inequalities.  One-sided
switching-minimality is not used.

## 2. A \(k\)-bit partition is an exact Fourier quotient

Let \(\Gamma\le\{\pm1\}^n\) be a subgroup under coordinatewise
multiplication, of rank \(k\).  Equivalently, assign a label
\(\lambda_i\in\mathbb F_2^k\) to every vertex and write

\[
g_t(i)=(-1)^{t\cdot\lambda_i},
\qquad t\in\mathbb F_2^k.
\]

Give edge \(ij\) the difference label

\[
d_{ij}=\lambda_i+\lambda_j
\]

and define the signed fiber totals

\[
W_d=\sum_{\substack{i<j\\d_{ij}=d}}a_{ij}.
\]

Then the energies on the entire symmetric-difference closure of the
\(k\) generating cuts are the Fourier transform of \(W\):

\[
\boxed{
H_A(g_t)=\sum_{d\in\mathbb F_2^k}W_d(-1)^{t\cdot d}.
}
\tag{2.1}
\]

Consequently Parseval gives

\[
\boxed{
2^{-k}\sum_tH_A(g_t)^2=\sum_dW_d^2.
}
\tag{2.2}
\]

Thus (1.1), restricted to \(\Gamma\), implies

\[
\sum_dW_d^2\le M^2.
\tag{2.3}
\]

This formulation keeps every two-cut identity.  For example, if
\(g_s,g_t\in\Gamma\), then \(g_sg_t\in\Gamma\), and

\[
D(s,t)=\frac{C(g_s)+C(g_t)-C(g_sg_t)}2
\]

is automatically the signed sum of edges crossing both cuts, with

\[
-\frac M2\le D(s,t)\le M.
\]

All triangle identities among the selected cuts are therefore already
encoded in the Fourier quotient.

## 3. Random labels recover flat second moments, but only the \(O(n)\) floor

Let \(N=\binom n2\), and now choose the labels
\(\lambda_1,\ldots,\lambda_n\) independently and uniformly in
\(\mathbb F_2^k\).  For two distinct edges \(e,f\),

\[
\Pr(d_e=d_f)=2^{-k},
\]

whether or not the edges share a vertex.  Expanding the right side of
(2.2) therefore gives the exact identity

\[
\boxed{
\mathbb E_\lambda\sum_dW_d^2
=(1-2^{-k})N+2^{-k}H_A(\mathbf1)^2.
}
\tag{3.1}
\]

If the two-sided cap holds on every character cut for every label
realization, (2.3) and (3.1) yield only

\[
M^2\ge
(1-2^{-k})N+2^{-k}M^2,
\qquad\text{hence}\qquad
\boxed{M^2\ge N.}
\tag{3.2}
\]

This remains only linear in \(n\), independently of \(k\).  Random
partitioning does recover the flat edge variance before contraction;
the loss is not a missing diagonal term.  The loss occurs because
second moments do not see the exponentially rare spin configurations
responsible for the \(n^{3/2}\) maximum.

There is an exact all-moments version.  For every test function
\(\Phi\),

\[
\boxed{
\mathbb E_\lambda\,2^{-k}\sum_t\Phi(H_A(g_t))
=2^{-k}\Phi(H_A(\mathbf1))
 (1-2^{-k})\mathbb E_x\Phi(H_A(x)),
}
\tag{3.3}
\]

where \(x\) is a uniform independent sign vector.  Indeed, \(t=0\)
gives \(\mathbf1\), while for every fixed \(t\ne0\) the variables
\((-1)^{t\cdot\lambda_i}\) are independent signs.  Thus random
\(k\)-bit quotients do not create a new moment inequality: after
averaging, they reproduce the ordinary Rademacher-chaos moments.

## 4. A flat-signing no-go theorem for every sublinear generated multicut hierarchy

The preceding failure is not peculiar to averaging.  There are actual
flat signings which satisfy the full two-sided cap on any prescribed
sublinear-rank cut group at \(o(n^{3/2})\) scale.

### Theorem 4.1

Let \(\Gamma\le\{\pm1\}^n\) have size \(L=2^k\).  There is a signing
\(B\) and a number \(M_\Gamma\) such that

\[
H_B(\mathbf1)=M_\Gamma,\qquad
|H_B(g)|\le M_\Gamma\quad(g\in\Gamma),
\tag{4.1}
\]

and

\[
\boxed{
M_\Gamma\le \sqrt{2N\log(48L)}.
}
\tag{4.2}
\]

Moreover, in the Fourier quotient (2.1) the same signing can be chosen
so that

\[
\boxed{
\sum_dW_d^2\ge\frac N2.
}
\tag{4.3}
\]

In particular, every tested cut obeys

\[
0\le C_B(g)\le M_\Gamma,
\]

all symmetric-difference and triangle identities among the cuts hold
exactly, every original edge still has magnitude one, and a constant
fraction of the flat second moment survives contraction.

#### Proof

Choose \(A\) with independent Rademacher edge signs.  For each fixed
\(g\), \(H_A(g)\) is a sum of \(N\) independent signs, so

\[
\Pr\{|H_A(g)|>t\}\le2e^{-t^2/(2N)}.
\]

With

\[
t=\sqrt{2N\log(48L)},
\]

a union bound shows

\[
\Pr\left\{\max_{g\in\Gamma}|H_A(g)|>t\right\}\le\frac1{24}.
\tag{4.4}
\]

For fixed edge-label fibers, let

\[
Z=\sum_dW_d^2.
\]

The fiber sums are independent, and if a fiber contains \(m\) edges,
its Rademacher sum \(S_m\) satisfies

\[
\mathbb ES_m^2=m,\qquad
\mathbb ES_m^4=3m^2-2m.
\]

It follows that

\[
\mathbb EZ=N,\qquad
\mathbb EZ^2
=N^2+2\sum_dm_d^2-2N
\le3N^2.
\]

Paley--Zygmund therefore gives

\[
\Pr\{Z\ge N/2\}\ge\frac1{12}.
\tag{4.5}
\]

The events (4.4)--(4.5) have positive-probability intersection.  Fix
such an \(A\), choose \(g_0\in\Gamma\) attaining the maximum absolute
energy, and put

\[
B_{ij}
=\operatorname{sgn}(H_A(g_0))\,a_{ij}g_0(i)g_0(j).
\]

Because \(\Gamma\) is a group,

\[
H_B(g)
=\operatorname{sgn}(H_A(g_0))H_A(g_0g),
\]

which proves (4.1)--(4.2).  Switching by \(g_0\) multiplies every
fiber total \(W_d\) by a character sign, so it leaves \(Z\) unchanged.
This proves (4.3).  \(\square\)

For \(k=o(n)\), (4.2) becomes

\[
\boxed{
\frac{M_\Gamma}{n^{3/2}}
=O\!\left(\sqrt{\frac{k+1}{n}}\right)=o(1).
}
\tag{4.6}
\]

Hence no argument which uses only \(k=o(n)\) generating cuts, their
entire symmetric-difference closure, the exact two-cut/triangle
identities inside that closure, and flat second moments can force even
the correct order \(n^{3/2}\).  Linear cut entropy is necessary.

The same statement covers a partition into \(r\) cells and all
cell-union cuts: those cuts form a group of rank at most \(r\), so
\(r=o(n)\) is again insufficient.  If “\(k\)-bit partition” instead
means \(2^k\) cells and all unions of cells, the obstruction applies
whenever \(2^k=o(n)\).

## 5. Even all local cuts through radius \(o(\sqrt n)\) are insufficient

There is a complementary deterministic obstruction which includes
all singleton cuts and many higher local cuts.

Suppose first that \(n\) is even.  Choose a balanced
\(u\in\{\pm1\}^n\) and set

\[
a_{ij}=-u_iu_j.
\]

Then

\[
H_A(x)=\frac{n-(u\cdot x)^2}{2},
\qquad
H_A(\mathbf1)=\frac n2.
\tag{5.1}
\]

For the spin vector obtained by flipping \(S\) from \(\mathbf1\),

\[
u\cdot x^S=-2\sum_{i\in S}u_i
\]

and hence

\[
\boxed{
C_A(S)=\left(\sum_{i\in S}u_i\right)^2.
}
\tag{5.2}
\]

Therefore

\[
0\le C_A(S)\le\frac n2
\qquad\text{whenever}\qquad
|S|\le\sqrt{\frac n2}.
\tag{5.3}
\]

By cut complementation the same is true for
\(|S|\ge n-\sqrt{n/2}\).  Thus a completely flat signing with only
\(M=n/2\) satisfies the exact two-sided inequalities for every cut in
a Hamming ball of radius \(\lfloor\sqrt{n/2}\rfloor\), including all
singletons and pairs.  It also satisfies every multi-cut identity
whose participating cuts and symmetric differences remain in that
ball.

For odd \(n\), take \(|\sum_i u_i|=1\).  The same calculation gives
\(H_A(\mathbf1)=(n-1)/2\), and the two-sided condition holds for every
\(S\) satisfying

\[
\left|\sum_i u_i-2\sum_{i\in S}u_i\right|^2\le2n-1.
\]

Again every Hamming ball of radius \(o(\sqrt n)\) is eventually
admissible.

This proves that bounded-size, fixed-degree, and even
\(o(\sqrt n)\)-radius local cut tests cannot detect the
\(n^{3/2}\) rigidity.  The obstruction is precisely the remote cut
\(x=u\), where \(H_A(u)=-\binom n2\).

## 6. Exact deletion and restriction second moments

The deletion/repair lemma starts from the switched row sums

\[
r_i=\sum_{j\ne i}a_{ij}=C_A(\{i\}),
\qquad
0\le r_i\le M,
\qquad
\sum_i r_i=2M.
\tag{6.1}
\]

Deleting a uniformly random vertex \(V\) leaves the bare signed total

\[
Y_{n-1}=M-r_V,
\]

so

\[
\boxed{
\mathbb EY_{n-1}=M\left(1-\frac2n\right),
\qquad
\operatorname{Var}(Y_{n-1})
=\frac1n\sum_i r_i^2-\frac{4M^2}{n^2}.
}
\tag{6.2}
\]

Re-minimizing the induced switching class replenishes an amount
\(2\Delta_V\) with

\[
0\le2\Delta_V\le r_V.
\tag{6.3}
\]

Thus the scalar deletion martingale has no forced loss beyond (6.2);
the allowed replenishment can cancel the entire one-step decrement.

There is also an exact formula at every restriction size.  Let \(T\)
be a uniformly random \(t\)-subset, put

\[
Y_T=\sum_{\{i,j\}\subseteq T}a_{ij},
\qquad
p_s=\frac{(t)_s}{(n)_s},
\]

and write \(N=\binom n2\).  Splitting two-edge products according to
whether the two edges coincide, meet in one vertex, or are disjoint
gives

\[
\boxed{
\mathbb EY_T^2
=p_4M^2+(p_2-p_4)N
 (p_3-p_4)\left(\sum_i r_i^2-n(n-1)\right).
}
\tag{6.4}
\]

The negative term in the last parenthesis is real: flat edge variance
does not force large row-sum variance because off-diagonal
correlations can cancel the diagonal contribution.  In the balanced
rank-one construction of Section 5, \(r_i=1\) for every \(i\), so
\(\sum_i r_i^2=n\) despite \(\sum_{i<j}a_{ij}^2=N\).

Equations (6.2)--(6.4) show exactly where a deletion argument must add
new information.  Tracking only the current total, row-sum second
moment, and flat Frobenius norm remains compatible with the rank-one
remote-cut obstruction.  A successful deletion martingale would have
to control the joint distribution of repair sets across many deletion
steps, or an equivalent linear-entropy family of remote cuts.

## 7. Consequence for the proposed route

The random-\(k\)-bit/multicut route has a sharp entropy threshold:

* second-moment averaging, even with exact flatness, gives only
  \(M\ge\sqrt{\binom n2}\);
* \(k=o(n)\) correlated cut generators and their full triangle closure
  admit actual flat signings with cap \(o(n^{3/2})\);
* all cuts within radius \(o(\sqrt n)\), including singleton and pair
  constraints, admit a flat rank-one signing with cap \(O(n)\);
* scalar deletion-repair data permit complete replenishment and do not
  prevent the same remote-cut obstruction.

Thus a two-sided multicut proof of the correct scale must retain a
family of \(\exp(\Omega(n))\) genuinely remote cuts, or use additional
structure equivalent in strength to such a family.  The natural next
positive target is not another finite quotient identity; it is a
linear-entropy theorem about the joint energy profile of remote cuts.
