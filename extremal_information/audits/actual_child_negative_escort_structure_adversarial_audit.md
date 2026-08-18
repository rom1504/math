# Adversarial audit: actual-child negative escort structure

**Object audited:**
[`drafts/actual_child_negative_escort_structure.md`](../drafts/actual_child_negative_escort_structure.md)

**Verdict:** **PASS, with an essential scope warning already stated in the
draft.**  I found no sign, factor, marginalization, or variational error in
AC.1--AC.33.  In particular, the note correctly avoids the tempting but false
conclusion that a linear escort gain by itself forces linear row total
correlation.

## 1. Escort and orientation identities

For fixed orientation, the forward output formula gives

```math
p_\epsilon(B)=\frac{d\Pi_\epsilon}{dU_B}(B)
=c_\epsilon\exp L_\epsilon(B).
```

The factor `c_epsilon` disappears after normalizing a negative power, so
AC.6 is exact.  Conditioning the joint escort on `epsilon` also gives AC.3;
there is no hidden sector-weight factor.  The KL peel AC.4 is the ordinary
chain rule, and its orientation contribution is at most `log 2`.

## 2. Marginalization in Lemma AC.1

This is the main technical point where a false proof would be easy.  Let `C`
be fixed coordinates, `S` retained coordinates, and `T` marginalized
coordinates.  Up to normalization, the retained density is

```math
h(s)=\mathbb E_{U_T}\exp[-\lambda F(s,C,T)].
```

If `s` and `s'` differ only in coordinate `j`, pairing equal `T` values gives

```math
e^{-\lambda c_j}\le h(s)/h(s')\le e^{\lambda c_j}.
```

Thus `g=log h` really does have coordinate oscillation at most
`lambda c_j` after both conditioning and marginalization.  Reveal the
retained bits sequentially and let `a_j` be their predictable half-log-odds.
Then `|a_j|<=lambda c_j/2`, and the likelihood ratio is

```math
\ell=\prod_j e^{a_jb_j}/\cosh a_j.
```

The product `prod_j e^(2a_jb_j)/cosh(2a_j)` is a mean-one fair-bit
martingale.  Hence

```math
\mathbb E_U\ell^2
\le\prod_j\{1+\tanh^2(\lambda c_j/2)\},
```

which verifies the sharp first inequality in AC.7; the quadratic bound
follows from `log(1+tanh^2u)<=u^2`.  The same half-log-odds bound gives AC.8
and all factors in AC.10--AC.12.

## 3. Sign check for AC.11

Writing

```math
\phi(\lambda)=-\log\mathbb E_Ue^{-\lambda L}
```

gives `phi'(lambda)=E_q L` and

```math
D(q\Vert U)=\phi(\lambda)-\lambda\phi'(\lambda).
```

Thus the sign in AC.11a is correct.  Summing the bitwise KL bound with
`c_j=2t` gives

```math
D(q\Vert U_B)\le mn\,\kappa(\lambda t)
\le \lambda^2t^2mn/2.
```

The total-correlation identity AC.13b follows by subtracting the marginal
row divergences from the divergence to uniform product measure.  It is an
upper bound only; nothing here supplies an extensive lower bound.

## 4. Directed product projection

For every bridge law `p`, direct substitution of
`dq/dU=e^{-lambda L}/E_Ue^{-lambda L}` gives

```math
D(p\Vert q)=\lambda\left[
\mathbb E_pL+\lambda^{-1}D(p\Vert U_B)-V_\lambda
\right].
```

Therefore minimizing over independent-row laws proves AC.16 exactly.  At a
global product minimizer, variation of one row while fixing the other rows
has a unique entropic minimizer, yielding AC.17.  Averaging over the other
rows preserves the `2t` one-bit oscillation, so AC.18 follows from AC.1.
The decomposition AC.24 is consequently an identity, not an inequality.

The two directed projections must not be interchanged:

```math
\inf_{p\ \text{ product}}D(p\Vert q)
\quad\hbox{and}\quad
\inf_{p\ \text{ product}}D(q\Vert p)=\operatorname{TC}(q)
```

need not control one another in dimension-free fashion.

## 5. Effective support

If `d_j<=kappa`, then

```math
(\sum_j\sqrt{d_j})^2/\sum_jd_j
\ge (\sum_jd_j)/\kappa,
```

which verifies AC.26.  The uniform cube law has subgaussian proxy
`t^2mn`, so entropy duality gives AC.27.  Concavity of `phi` gives
`E_qL<=V_lambda`, hence a gain `G_lambda>=eta N` lower-bounds the actual
uniform-to-escort energy displacement by the same amount.  The stated
`Omega(N)` KL and `Omega(N^2)` effective-support conclusions are therefore
correct for comparable splits.

## 6. Exact non-product falsifier

Condition on `epsilon`.  The latent child law retains full support at every
finite `t`.  For distinct rows `i,i'` and columns `j,j'`, its rank-one bridge
word satisfies

```math
Q_{ij}Q_{ij'}=y_jy_{j'},
\qquad
Q_{ij}Q_{ij'}Q_{i'j}Q_{i'j'}=1.
```

Passing through the independent channel of bias `rho=tanh t` gives exactly
AC.29.  Full support implies
`|E[y_jy_{j'}\mid epsilon]|<1`, so row independence would equate `rho^4`
with a strictly smaller `rho^4c_{jj'}^2`.  This proves AC.30.

Under the canonical child-spin posterior, multiplication of the escort
density by the posterior density gives

```math
q(B\mid Z,\epsilon)
\propto \Pi(B\mid Z,\epsilon)p_\epsilon(B)^{-(\lambda+1)}.
```

The first factor is a strictly positive row product.  If the result were a
row product, division would make `p_epsilon` a row product, contradicting
the rectangle calculation.  AC.31 is therefore a valid exact falsifier of
the obvious Gibbs-latent decomposition.

## 7. Optimizer inequalities

Flipping a child edge set `S` changes its Hamiltonian by
`-2 sum_(e in S) a_e x_e`.  The ratio between the flipped and original
augmented partition functions is precisely the expectation in AC.32.
Minimality makes it at least one.  For one edge,

```math
\cosh(2t)-a_e\mathbb E[\tau x_ux_v]\sinh(2t)\ge1,
```

which is equivalent to AC.33.  This is the only part of the theorem note
that exploits optimizer status beyond selecting the children, and the draft
correctly does not claim that it controls AC.24.

## 8. Sharp scope warning and falsifier

Uniform row-filtration `D_2` is **not** a common-latent row-product
decomposition.  A filtration kernel may depend on the complete preceding
row history, so treating that history as a reusable latent state can encode
the full bridge law.

This is not merely a missing proof.  Tight autoregressive components do not
imply the hypothesis of the archived latent-iid no-gain theorem.  Let the row
alphabet be a finite cube `X`, let `U_X` be uniform, and choose a balanced
character `h:X->\{+-1\}`.  For `m>=3` and `0<delta<1`, put

```math
 {dq_m\over dU_X^{\otimes m}}(r_1,\ldots,r_m)
 =1+\delta\prod_{i=1}^m h(r_i).                       \tag{A.1}
```

In every row reveal order, the first `m-1` conditional rows are exactly
uniform, while the last has density

```math
 1+\delta\left(\prod_{i<m}h(r_i)\right)h(r_m).
```

Consequently every conditional component satisfies

```math
 D_2(q_m(R_i\mid R_{<i})\Vert U_X)
 \le \log(1+\delta^2),                                \tag{A.2}
```

uniformly in `m`.  Nevertheless `q_m` has no representation
`int nu_z^{\otimes m} pi(dz)`.  Indeed, such a representation would give

```math
 E h(R_1)h(R_2)=\int(E_{\nu_z}h)^2\,\pi(dz).
```

The left side is zero under (A.1), so `E_{nu_z}h=0` almost surely.  It would
then force `E prod_i h(R_i)=0`, whereas (A.1) gives exactly `delta`.  Taking
`h` to be an even two-coordinate character also makes the law invariant
under simultaneous global sign reversal of all bridge bits.  Thus tight
row-filtration `D_2`, row exchangeability, and the relevant sign symmetry
still do not produce a latent-iid decomposition.  The no-gain theorem
cannot be invoked without an additional de Finetti/extendibility or
approximate-product transfer theorem.  The actual child law is generally
not even row-exchangeable, since the child signing labels its rows.
This example is a falsifier of the claimed hypothesis transfer, not of the
no-gain conclusion: its total correlation is only `O_delta(1)`, so a
different archived theorem may still control its particular pressure.

The barrier persists at extensive total correlation and even at the actual
escort's `O(n^(-1/2))` one-bit oscillation scale.  Take rows of length `3r`.
On a first block of `r` bits, choose a common fair latent sign `Z` and,
conditional on `Z`, draw the rows iid from bit-products of bias
`Z tanh(c/sqrt(r))`.  Group the remaining `2r` bits into disjoint pairs and
let `h_j` be the product in pair `j`.  Independently of the first block,
tilt the fair parity coordinates by

```math
 \prod_{j=1}^r\left(1+{d\over\sqrt r}
                  \prod_{i=1}^m h_j(R_i)\right).
```

Posterior mixtures of the two weakly biased first-block laws have row
Renyi-two constant at most `e^(c^2)`.  For every row order the first `m-1`
parity conditionals are uniform, and the last has Renyi-two constant
`(1+d^2/r)^r<=e^(d^2)`.  Every one-bit log-density oscillation is
`O((c+d)/sqrt(r))`.  Thus every autoregressive row component is uniformly
tight at the same local scale as AC.10.  The common-latent sign symmetry and
the even characters also make the joint law invariant under global bridge
reversal.

The common weak latent contributes `Omega(m)` total correlation by the
Gaussian row-sum limit.  For every `j`, however, the pair moment
`E h_j(R_1)h_j(R_2)` is zero while the `m`-row moment is `d/sqrt(r)>0`; the
same square-moment argument still forbids *every* latent-iid representation.
The parity factor is not even `o(1)`-close to its fair counterpart: its `r`
independent parity products have bias `d/sqrt(r)`, and the central limit
theorem gives a nonzero limiting total-variation separation.  Hence tight
filtration complexity can coexist simultaneously with the two exact
loopholes left by Theorems 37.14 and 37.17: linear row dependence and failure
of common-latent iid disintegration.  This remains a structural scope
falsifier, not a favorable-pressure construction.

There is a second independent scope barrier: Theorems 37.13--37.17 prove
no-gain for the **conference-child pressure** `f_r`, using its projected
coupling structure.  They are not law-only inequalities for an arbitrary
parent landscape.  Even an exact tight-component latent-iid representation
of the actual escort would not by itself transfer those theorems to
`L_epsilon` for arbitrary optimized children.  One also needs a uniform
actual-child analogue of the component no-gain estimate (RT.2), or a proved
comparison of the actual-child shortfall to the conference observable.

The exact transfer that *would* be legitimate is correspondingly strong.
If there are laws

```math
 \widetilde q_m=\int\nu_{z,m}^{\otimes m}\,\pi_m(dz)
```

whose component `D_2(nu_(z,m)||U_X)` is tight in the sense of RT.5 and
`TV(q_m,\widetilde q_m)=o(1)`, then for the bounded shortfall observable
`0<=S_m<=h_beta`,

```math
 E_{q_m}S_m
 \le E_{\widetilde q_m}S_m+h_\beta TV(q_m,\widetilde q_m)=o(1). \tag{A.3}
```

This follows directly from RT.1 **when `S_m` is the conference shortfall**
and is the minimal generic route from an approximate latent-iid
representation in that benchmark.  AC.10 supplies neither the mixture nor
the required `o(1)` law-level error, and it does not replace the observable.
In particular an `o(m)` KL or total-correlation estimate is far too weak:
Pinsker only helps when the relevant KL itself is `o(1)`.

Nor can bounded differences alone turn a linear gain into extensive row
dependence.  On an `m by n` Boolean bridge, consider the shifted linear toy
potential

```math
L_0(B)=C-t\sum_{ij}B_{ij}.
```

It has the same `2t` one-bit oscillation as the parent log partition.  Its
negative escort is the independent-bit law of bias `tanh(lambda t)`, so both
row total correlation and the reverse product-projection cost vanish
exactly.  Nevertheless

```math
\mathbb E_UL_0-V_\lambda
=\frac{mn}{\lambda}\log\cosh(\lambda t)
=\Theta(N)
```

for comparable splits and `t=beta/sqrt(N)`.  The example is not a signing
partition function; it proves that an optimizer-specific theorem about the
row-product term in AC.24 is indispensable.

There is a second, stronger warning: even **linear** total correlation does
not certify irreducible dependence under the regularity proved in AC.2.
Let the bridge have `m` rows of length `n`, put `a=c/sqrt(n)`, and let
`P_+` and `P_-` be the bit-product laws with respective biases
`+tanh(a)` and `-tanh(a)`.  For fixed `0<delta<1`, set

```math
Q=(1-\delta)U+{\delta\over2}(P_++P_-).                 \tag{A.4}
```

This law is invariant under global sign reversal.  Its log density has
one-bit oscillation at most `2a`: every component likelihood changes by a
factor in `[e^(-2a),e^(2a)]`, and the same interval contains their mixture
ratio.  After every row prefix, the next-row kernel is a posterior mixture
of the three product rows, so

```math
D_2(Q(R_i\mid R_{<i})\Vert U_n)
\le n\log(1+\tanh^2a)\le c^2.                         \tag{A.5}
```

With the three-valued mixture label `Z`, rows are iid conditionally on `Z`,
and

```math
\operatorname{TC}(Q)
=mI(Z;R_1)-I(Z;R_1,\ldots,R_m).                        \tag{A.6}
```

The central-limit theorem for the normalized row sum shows
`I(Z;R_1)->i(c,delta)>0`: the limiting row-sum laws are a standard Gaussian
and the symmetric mixture of translates by `+-c`.  Since the second term in
(A.6) is at most `H(Z)`, `TC(Q)=Omega(m)`.  Nevertheless `U` is row product
and `dQ/dU>=1-delta`, so

```math
\inf_{p\ \text{ row product}}D(p\Vert Q)
\le D(U\Vert Q)\le-\log(1-\delta)=O(1).               \tag{A.7}
```

Thus global-sign symmetry, weak single-bit oscillation, uniformly tight
filtration `D_2`, and extensive total correlation can all coexist with a
constant reverse product-projection cost.  The actual-child alternative must
be phrased using the directed variational charge in AC.16 (or an equally
strong common-latent complexity), not total correlation alone.

Accordingly the new SML is genuinely narrower than “understand the actual
child law,” but it is not yet a Level-5-to-6 bridge: determine, using AC.32
or another minimizing-child identity, whether the linear part of AC.24 lies
in the bounded-`D_2` product shadow or in the directed dependence charge.
Any claim replacing that charge by ordinary total correlation requires an
additional theorem.

There is also an information-complexity caveat.  Each factor in the product
shadow is still an arbitrary law on `2^n` rows, and its fixed-point equation
AC.17 queries the expectation of the full parent landscape.  More sharply,
the complete best-response oracle

```math
(r,p_{-i})\longmapsto E_{p_{-i}}L(r,R_{-i})
```

reconstructs the whole table of `L`: choose every factor of `p_-i` to be a
point mass.  AC.17 only evaluates this oracle at a minimizing fixed point,
so the fixed-point data do **not** logically reconstruct `L`; they contain
only `m(2^n-1)` row-law parameters rather than `2^(mn)` landscape values.
But no theorem finds or controls that fixed point from a smaller child
statistic, and a generic implementation of the response map has access to
the full landscape.  Moreover, as `lambda` tends to infinity, row-product
laws contain every point mass `delta_B`, so the row-product variational
minimum tends to the full deterministic bridge minimum.

Thus AC.24 is a strict fixed-`lambda` structural decomposition, but it is
not by itself an information-theoretically strict reduction.  Even the SML
“decide which term is linear” is presently an exhaustive classification,
not a demonstrably easier theorem: it could require evaluating the same
full response landscape.  A successful next lemma must control one term
uniformly from a stated low-information statistic of the minimizing
children, without solving AC.17 by enumerating the parent bridge landscape.
