# Dependent two-step AMP / Onsager rounding

## Status

This note separates three claims:

1. the scalar two-step variational problem and its numerical value are
   exact consequences of the stated population law;
2. the fixed-dither Hermite tail reduction and finite-degree diagram
   lemma prove the resulting lower bound for every sequence of exact
   symmetric conference matrices;
3. extension from conference matrices to every competing signing
   remains the proof obligation identified below.

No universal \(0.78338\ldots\) theorem is claimed yet.

## 1. The paired population law

Let \(U=U^\top\), \(U^2=I\), \(U_{ii}=0\), and
\(|U_{ij}|=(n-1)^{-1/2}\).  The intended application is
\(U=A/\sqrt{n-1}\) for a conference signing.

At one coordinate let
\[
S\in\{\pm1\},\qquad G,W\sim N(0,1)
\]
be independent.  For a threshold \(t\), put
\[
F=\operatorname{sign}(G+tS),
\]
\[
a=\mathbb E(GF)=2\phi(t),\qquad
b=\mathbb E(SF)=2\Phi(t)-1,
\]
\[
s^2=1-a^2-b^2,\qquad
R=\frac{F-aG-bS}{s}.
\tag{1}
\]
Thus \(R\) has mean zero and variance one and is orthogonal to both
\(S\) and \(G\).

The matrix action pairs the scalar directions
\[
US=G,\quad UG=S,\qquad UR=W,\quad UW=R.
\tag{2}
\]
For a second Boolean response
\[
Y=g(S,G,R,W),
\]
write
\[
\alpha=\mathbb E(SY),\quad
\beta=\mathbb E(GY),\quad
\gamma=\mathbb E(RY),\quad
\delta=\mathbb E(WY).
\]
The paired state-evolution energy is
\[
\boxed{
\mathcal E(g)=2(\alpha\beta+\gamma\delta).
}
\tag{3}
\]
This is genuinely dependent: \(R\) is the nonlinear residual of the
first response, and \(W=UR\) is its matrix backtracking partner.  It is
not an independent extra probe.

The first variation of (3) shows that every nondegenerate stationary
response obeys
\[
\boxed{
Y=\operatorname{sign}
\left(\beta S+\alpha G+\delta R+\gamma W\right).
}
\tag{4}
\]
This gives a four-dimensional self-consistency iteration rather than an
unstructured search over Boolean boundaries.

## 2. Optimized two-step constant

Splitting the one-dimensional Gaussian integral at \(G=-tS\), solving
(4), and then optimizing \(t\) gives
\[
\boxed{
t_2=0.8414699114\ldots,\qquad
c_2=0.783387533648\ldots .
}
\tag{5}
\]
At the optimum, one choice of score coefficients in the order
\((S,G,R,W)\) is
\[
(p,q,r,d)=
(0.5859761744,\ 0.6179560304,\
 0.2396817825,\ 0.4661704739),
\tag{6}
\]
and the corresponding moments are
\[
(\alpha,\beta,\gamma,\delta)=
(0.5618430437,\ 0.5327670921,\
 0.4238402492,\ 0.2179176762).
\tag{7}
\]
Substitution in (3) gives (5).

For an exact conference matrix, the same rule can be written without
explicit residual variables:
\[
\boxed{
Y=\operatorname{sign}\left(
-0.1225227631\,S
-0.1063660230\,US
+0.4194582166\,F
+0.8158276928\,UF
\right).
}
\tag{8}
\]
The negative \(S,US\) terms are the Onsager/backtracking correction.
Removing them is precisely the invalid “fresh residual” approximation.

The reproducible population calculation is
`_two_step_amp_variation.py`.  Whole-line Gauss--Hermite quadrature must
not be used without splitting at \(G=-tS\): doing so biased an earlier
value in the fourth decimal place.

Direct simulations of (8) on symmetric Paley conference matrices give:

| order | mean \(Y^\top UY/n\) (100 starts) |
|---:|---:|
| 102 | \(0.7776\) |
| 194 | \(0.7854\) |
| 402 | \(0.7861\) |
| 602 | \(0.7834\) |

This is strong numerical confirmation of (5), not a substitute for the
uniform conference proof.

## 3. Fixed-dither formulation

For a proof, keep fixed \(\tau_0,\tau_1>0\) and use the soft first
message
\[
F_{\tau_0}
=\psi_{\tau_0}(G+tS),\qquad
\psi_\tau(z)=2\Phi(z/\tau)-1.
\tag{9}
\]
Let
\[
a_\tau=\mathbb E(GF_\tau),\quad
b_\tau=\mathbb E(SF_\tau),\quad
s_\tau^2=\mathbb EF_\tau^2-a_\tau^2-b_\tau^2,
\]
\[
R_\tau=(F_\tau-a_\tau G-b_\tau S)/s_\tau.
\tag{10}
\]
The final Boolean output is randomized with conditional mean
\[
\Psi=
\psi_{\tau_1}(pS+qG+rR_\tau+dW).
\tag{11}
\]
For distinct coordinates, conditional independence of the final
dithers gives
\[
\mathbb E(Y_iY_j)=\mathbb E(\Psi_i\Psi_j).
\]
One first takes \(n\to\infty\), then sends
\(\tau_0,\tau_1\downarrow0\).

## 4. Rank-three Onsager cancellation

Write
\[
F_\tau(g,s)=u_\tau(g)+s\,v_\tau(g).
\]
Oddness under \((g,s)\mapsto(-g,-s)\) makes \(u_\tau\) odd and
\(v_\tau\) even.  In the Gaussian Hermite expansion, subtracting
\(a_\tau g+b_\tau s\) removes exactly
\[
H_1(g),\qquad sH_0(g).
\]
Consequently
\[
\boxed{
R_\tau
=\sum_{\substack{\ell\ge3\\\ell\ {\rm odd}}}
u_\ell H_\ell(G)
+S\sum_{\substack{\ell\ge2\\\ell\ {\rm even}}}
v_\ell H_\ell(G).
}
\tag{12}
\]
Every surviving term has total hybrid rank at least three.  This is the
Onsager cancellation responsible for a fresh second field.  Without the
two subtracted linear components, rank-one diagrams survive and the
claimed state evolution is false.

## 5. Rigorous Hermite-tail reduction

Let \(R_{\tau,L}\) be (12) truncated at degree \(L\), and set
\[
W_\tau=UR_\tau,\qquad W_{\tau,L}=UR_{\tau,L}.
\]
For each conference row, \(G_i\) is a normalized Rademacher sum
independent of \(S_i\).  Scalar CLT and polynomial approximation give
\[
\lim_{n\to\infty}
\frac1n\mathbb E\|R_\tau-R_{\tau,L}\|_2^2
=\operatorname{Tail}_\tau(L),
\tag{13}
\]
where \(\operatorname{Tail}_\tau(L)\downarrow0\).

Conference orthogonality gives the exact isometry
\[
\boxed{
\|W_\tau-W_{\tau,L}\|_2
=\|R_\tau-R_{\tau,L}\|_2.
}
\tag{14}
\]
Since
\[
\operatorname{Lip}(\psi_{\tau_1})
=\sqrt{2/\pi}\,\tau_1^{-1},
\]
the squared empirical difference between the final soft messages is at
most
\[
C_{\tau_1,r,d}\,
\frac{\|R_\tau-R_{\tau,L}\|_2^2+
\|W_\tau-W_{\tau,L}\|_2^2}{n}.
\tag{15}
\]
Finally, \(\|U\|_{\rm op}=1\) and \(\|\Psi\|_2\le\sqrt n\) imply
\[
\left|
\frac1n\mathbb E\Psi^\top U\Psi
-\frac1n\mathbb E\Psi_L^\top U\Psi_L
\right|
\le
C_{\tau_1,r,d}\sqrt{\operatorname{Tail}_\tau(L)}+o_n(1).
\tag{16}
\]
Thus the order
\[
n\to\infty,\qquad L\to\infty,\qquad
\tau_0,\tau_1\downarrow0
\]
is legitimate.  The infinite Hermite tail is not a remaining gap.

## 6. Finite-degree conference lemma

The needed finite polynomial statement follows from homogeneous-sum
universality and an explicit contraction calculation.

> **Flat-involution diagram lemma.**  Fix \(L\).  Substitute
> \(G=US\) into a hybrid polynomial of the form (12), let
> \(W=UR_L\), and test it with a fixed smooth function of
> \((S,G,R_L,W)\).  Then the empirical one-site law converges to
> independent \((S,G,R,W)\), and the oriented two-site response has
> only the two first-order channels
> \[
> S_i\leftrightarrow G_j,\qquad
> R_i\leftrightarrow W_j.
> \]
> Consequently
> \[
> \frac1n\mathbb E Y^\top UY
> =2(\alpha\beta+\gamma\delta)
> +O_{\tau_0,\tau_1,L}(n^{-c_L}).
> \tag{17}
> \]
> Here \(c_L>0\); no uniformity in \(L\) is asserted or needed.

Here are the contraction details.  Write \(u_k\) for row \(k\) of
\(U\).  After multilinearizing the Rademacher inputs, the leading
order-\(\ell\) coefficient tensor of a pure \(H_\ell(G_k)\) term in
\(W_i\) is
\[
T_i(j_1,\ldots,j_\ell)
=\sum_k U_{ik}\prod_{r=1}^\ell U_{kj_r}.
\tag{18}
\]
Orthogonality gives, for every input coordinate \(j\),
\[
\sum_{j_2,\ldots,j_\ell}
T_i(j,j_2,\ldots,j_\ell)^2
=\sum_kU_{ik}^2U_{kj}^2
=\frac1{n-1}.
\tag{19}
\]
For a hybrid \(S_kH_{\ell-1}(G_k)\) term the local tensor is the
symmetrization of
\[
e_k\otimes u_k^{\otimes(\ell-1)}.
\tag{20}
\]
The assumptions \(U_{kk}=0\), \(u_k\perp u_l\) for \(k\ne l\), and
\(|\langle e_k,u_l\rangle|=(n-1)^{-1/2}\) imply the same
\(O_L(1/n)\) influence bound.

More generally, in each odd total chaos degree \(\ell\ge3\), the local
tensor is a linear combination of
\[
u_k^{\otimes\ell},
\qquad
\operatorname{Sym}
\left(e_k\otimes u_k^{\otimes(\ell-1)}\right).
\tag{21}
\]
These local tensors are mutually orthogonal as \(k\) varies.  A direct
expansion of every nontrivial \(r\)-contraction, \(1\le r<\ell\), gives
\[
\boxed{
\|T_i\otimes_rT_i\|_{\rm HS}^2
\le \frac{C_L}{n-1}.
}
\tag{22}
\]
For the pure tensors this is exactly
\(\sum_kU_{ik}^4=1/(n-1)\).  For the hybrid tensors, every cross term
either contains \(\langle u_k,u_l\rangle=0\), or contains at least two
cross factors
\(\langle e_k,u_l\rangle\langle u_k,e_l\rangle=U_{kl}^2\);
summation and orthogonality give the stated \(C_L/(n-1)\) bound.

Repeated Rademacher indices do not create a new leading diagram.
Multilinearizing \(H_\ell(u_k\cdot S)\) leaves the Wick-ordered
distinct-index chaos plus lower-order leakage whose total \(L^2\)
mass is \(O_L((n-1)^{-1})\).  In (20), terms in which the distinguished
index \(k\) reappears among the \(G_k\) indices have the same bound,
using \(U_{kk}=0\).

The homogeneous-sum invariance theorem now replaces the Rademachers by
Gaussians with error \(O_L(n^{-c_L})\).  Equations (19) and (22), via
the multivariate fourth-moment theorem, give a standard Gaussian
second field independent of the local variables.  The same contraction
calculation after deleting two coordinates gives the two-site cavity
law.

Finally set \(\varepsilon=U_{ij}\).  With the direct edge removed,
\[
G_i=G_i^{(ij)}+\varepsilon S_j,
\qquad
W_i=W_i^{(ij)}+\varepsilon R_j+O(\varepsilon^2),
\tag{23}
\]
and symmetrically for \(j\).  Taylor expansion of the fixed-dither
response and Gaussian integration by parts give
\[
\mathbb E(\Psi_i\Psi_j)
=\text{orientation-even baseline}
+2\varepsilon(\alpha\beta+\gamma\delta)
+O_{\tau_0,\tau_1,L}(\varepsilon^2)
+O_{\tau_0,\tau_1,L}(n^{-c_L}|\varepsilon|).
\tag{24}
\]
Average the two matrix orientations to cancel the baseline, multiply
by \(U_{ij}\), and sum.  Since
\[
\frac1n\sum_{i\ne j}U_{ij}^2=1,
\qquad
\frac1n\sum_{i\ne j}|U_{ij}|^3=O(n^{-1/2}),
\]
equation (17) follows.

Combining the finite-degree lemma with the tail passage in Section 5
and then sending both dithers to zero proves the conference theorem:
\[
\boxed{
\liminf_{n\to\infty}
\frac{Q(C_n)}{n\sqrt{n-1}}
\ge c_2=0.783387533648\ldots
}
\tag{17a}
\]
for every sequence of exact symmetric conference matrices \(C_n\).
The decimal is the numerical evaluation of the explicitly stated
one-dimensional population variational problem; the theorem itself
uses that variational value, not a simulation fit.

## 7. Extension beyond conference matrices

For a general normalized signing
\[
B=A/\sqrt{n-1},\qquad
q_{ij}=(B^2)_{ij}\quad(i\ne j),
\]
the spectral estimate
\[
\|A\|_{\rm op}^2\le2Q(A)
\]
implies, on every \(Q(A)=O(n^{3/2})\) competing sequence,
\[
\frac1{n^2}\sum_{i\ne j}q_{ij}^2=o(1).
\tag{20}
\]
Rank three improves the relevant obstruction: the first residual
correlation is cubic in \(q_{ij}\), so the variance correction to
\(BR\) begins at
\[
\boxed{
\Lambda_4(B)
=\frac1n\sum_{i\ne j}q_{ij}^4.
}
\tag{21}
\]
If \(\Lambda_4(B)\to0\), the same diagram bounds are expected to give
the conference law (17), hence the \(c_2\) witness.

There is now a rigorous spectral purification before this dichotomy.
The weighted same-sign Grothendieck bound gives
\[
Q(A)\ge\frac{\operatorname{tr}|A|^3}{2K_G(n-1)}.
\tag{21a}
\]
Here the factor \(2K_G\) is obtained by a mixed-sign bilinear
assignment: the left vectors have spectral coordinates
\(|\lambda_r|u_r(i)\), and the right vectors have coordinates
\(\operatorname{sgn}(\lambda_r)|\lambda_r|u_r(i)\).
Combining this with
\[
\operatorname{tr}|B|^3
\ge\frac{\operatorname{tr}B^4}{\|B\|_{\rm op}},
\qquad
\|A\|_{\rm op}^2\le2Q(A),
\]
shows that every sequence capable of violating the target \(c_2\)
must obey
\[
\frac1n\operatorname{tr}B^4
\le(3.49519+o(1))n^{1/4}.
\tag{21b}
\]
Equivalently,
\[
\sum_{i\ne j}q_{ij}^2=O(n^{5/4}).
\tag{21c}
\]
This rules out the formerly problematic dense profile
\(q_{ij}\asymp n^{-1/4}\) on \(\Theta(n^2)\) pairs.  A critical
profile not yet excluded is
\[
|q_{ij}|\asymp n^{-1/8}
\quad\text{on }\Theta(n^{3/2})\text{ pairs},
\tag{21d}
\]
which simultaneously has (21c) and \(\Lambda_4(B)=\Theta(1)\).

In fact, \(\Lambda_4=o(1)\) cannot be forced even after deleting
\(o(n)\) vertices.  Here is an explicit flat-sign obstruction.  Let
\(C_k\) be a symmetric conference matrix and put
\[
R=J_3-2I_3,\qquad D=J_3-I_3,
\]
\[
\boxed{
A_k=C_k\otimes R+I_k\otimes D.
}
\tag{21e}
\]
This is a symmetric zero-diagonal signing of order \(n=3k\).  Since
\[
R^2=4I_3-J_3,\qquad D^2=I_3+J_3,\qquad RD+DR=4I_3,
\]
and \(C_k^2=(k-1)I_k\),
\[
A_k^2
=I_k\otimes\big((k-1)R^2+D^2\big)+4C_k\otimes I_3.
\tag{21f}
\]
Consequently, within each three-vertex fibre,
\[
q_{(u,a),(u,b)}
=\frac{-k+2}{3k-1}\longrightarrow-\frac13
\qquad(a\ne b),
\]
whereas between distinct fibres the only nonzero correlations are
\[
q_{(u,a),(v,a)}
=\frac{4(C_k)_{uv}}{3k-1}=O(k^{-1}).
\]
It follows exactly that
\[
\boxed{
\Lambda_4(A_k/\sqrt{3k-1})\longrightarrow\frac2{81}>0.
}
\tag{21g}
\]
Deleting \(o(n)\) vertices damages only \(o(k)\) of the disjoint
three-vertex fibres, so the same positive limit remains on every
\((1-o(1))n\)-vertex principal submatrix.

This family still has \(Q(A_k)=O(n^{3/2})\): indeed
\[
\|A_k\|_{\rm op}
\le \sqrt{k-1}\,\|R\|_{\rm op}+\|D\|_{\rm op}
=2\sqrt{k-1}+2,
\]
and \(Q(A_k)\le n\|A_k\|_{\rm op}\).
Thus neither the \(O(n^{3/2})\) hypothesis nor \(o(n)\)-vertex
purification can make the conference diagram law universal.  The
missing theorem must handle a finite-type correlated cavity law,
not merely eliminate it.

If \(\Lambda_4\not\to0\), mean-square pseudorthogonality is
insufficient.  Thresholding the graph of large \(|q_{ij}|\) yields a
precise combinatorial dichotomy:

* if its large-correlation edges have an \(o(n)\)-vertex cover, delete
  those vertices and apply monotonicity on the remaining principal
  submatrix;
* otherwise a linear-sized matching (at a suitable threshold scale)
  of correlated row pairs survives.

The missing universal lemma is to convert the second branch at the
refined scale (21d) into a Boolean quadratic witness of size at least
\[
(c_2-o(1))n^{3/2}.
\tag{22}
\]
Fixed-degree moment/SOS bounds do not do this: a constant amount of
\(\Lambda_4\) changes fixed moments only below the \(n^{3/2}\) scale.
An orientation-even nonlinear rounding or a principal-subset
purification theorem is required.

## 8. Current conclusion

The dependent second step is quantitatively real and beats every
independent-probe rule:
\[
0.6729867289\longrightarrow0.7833875336
\]
in doubled normalization on the paired flat-involution state evolution.
The exact implementable rule is (8), and the infinite Hermite tail is
controlled by (13)--(16).

One gap remains before this becomes a universal lower bound:

1. prove that every non-conference finite-type cavity law, including
   the explicit three-fibre obstruction (21e), supplies at least the
   same \(c_2\) lower bound.
