# Power-saving near-order regularization by exact sign gadgets

2026-09-17. **Verified derivation:** root construction, independently
reconstructed by the Bernoulli and localization researchers. External
novelty is unestablished. This is an unconditional transformation of
EVERY actual full signing, not a favorable-profile assumption.

## 1. The theorem

For every sufficiently large n and every hollow symmetric full signing
A of order n, there is a hollow symmetric full signing W containing A
as an unchanged principal block, of order

```math
N=n+q,\qquad q=\lfloor n^{2/3}\rfloor\lfloor n^{1/24}\rfloor,
```

with the following simultaneous properties, for universal constants:

```math
\begin{split}
Q(W)&\le Q(A)+O(n^{11/8}),\\
w\bigl(\{z\in\{\pm1\}^N:Q(W)-|H_W(z)|\le n^{5/4}\}\bigr)
 &\le O(n^{15/16}),\\
\log\bigl|\{z:Q(W)-|H_W(z)|\le n^{5/4}\}\bigr|
 &\le O(n^{23/24}\log n).
\end{split}                                                   \tag{1}
```

Here w(C)=E max_(z in C) g^Tz is ordinary Gaussian width, not a
factorization norm. All constants are uniform in A; no spectral
condition or prior near-ground geometry is used. Both polarities and
EVERY new spin are included. This is an existence theorem, not a
polynomial-time optimizer.

Consequently any liminf-realizing sequence may be replaced by actual
near-order full signings with the same limiting normalized cap and
subexponentially many states even in a polynomially large near-ground
window. This is NOT all-order recovery and does not prove convergence.
The controlled normalized window n^(-1/4) shrinks with n.

## 2. The exact physical gadget

Put K=floor(n^(2/3)), r=floor(n^(1/24)), q=Kr. Choose an n by K
array h of independent fair signs. The bridge has r IDENTICAL copies
of its a-th column h^a for every a<=K. Let D be any full signing of
order q with Q(D)<=q^(3/2), available for all sufficiently large q
by the elementary independent-sign union bound. Set

```math
W=\begin{pmatrix}A&C\\C^T&D\end{pmatrix},\qquad
F_h(x)=\frac{|H_A(x)|}{\sqrt n}
       +\frac r{\sqrt n}\sum_{a=1}^K|h^a\cdot x|,
\qquad M(h)=\max_x F_h(x).
```

Ignoring D, maximizing the new spins gives exactly sqrt(n) M(h).
Indeed the sum of the r spins in group a ranges in[-r,r], and the
maximum of an affine function is attained at its endpoints. The two
polarities can be absorbed into the endpoint pattern sigma in{+-1}^K.
Equivalently

```math
M(h)=\max_{x,\sigma}
 \left[|H_A(x)|/\sqrt n+
       (r/\sqrt n)\sum_{i,a}h_{ia}x_i\sigma_a\right].       \tag{2}
```

Restoring D changes the cap by at most Q(D). For ANY near-extreme
parent word (x,y) of tolerance T, the projection x satisfies

```math
F_h(x)\ge M(h)-(T+2Q(D))/\sqrt n.                         \tag{3}
```

This follows by comparing the actual energy to its endpoint maximum,
not by assuming that every optimizing group is constant in the presence
of D. Thus no hidden restriction on the new spins is made.

## 3. Gaussian smoothing, while keeping the hidden maximum

Write H0(x)=|H_A(x)|/sqrt(n), delta=r sqrt(K/n), and

```math
G(t)=\mathbb E_g\max_x[H0(x)+\sqrt t\,g\cdot x].
```

For all t,s>0, coupling the SAME Gaussian vector at the two scales gives

```math
G(t+s)-G(t)\le\sqrt{2/\pi}\,n(\sqrt{t+s}-\sqrt t)
 \le\sqrt{2/\pi}\,ns/(2\sqrt t).                         \tag{4}
```

No differentiability or random-child hypothesis is required.

Let M(h;u) denote(2) with an added field u dot x, and define

```math
\Delta_s(h)=\mathbb E_{g'}M(h;\sqrt s\,g')-M(h).
```

For the joint old-spin nearcode C_b(h)={x:F_h(x)>=M(h)-b},
choosing its best x under the secondary field proves pointwise

```math
\sqrt s\,w(C_b(h))\le b+\Delta_s(h).                      \tag{5}
```

This avoids the invalid use of driver-dependent optimum offsets in a
replacement theorem. We compare the expectations of TWO separate
maxima, each of which has deterministic offsets (conditioning on g').

The coefficients of every driver h_ia in(2) have magnitude r/sqrt(n).
There are nK drivers, at most2^(n+K) queries, and fourth budget

```math
W_4=Kr^4/n,\qquad L=(n+K+1)\log2,\qquad
U=8W_4^{1/4}L^{3/4}.
```

The independently reconstructed fourth-order sign/Gaussian maximum
comparison bounds either endpoint replacement error by U, uniformly
in ALL offsets and in the secondary field. In this gadget the driver
magnitudes are query independent; even the separate-coordinate symmetric
fourth-order proof suffices. We do not credit a nonexistent extra gain
to the common-Gibbs refinement here.

At the Gaussian endpoint each fixed sigma has an iid Gaussian field
of variance delta^2, despite dependence between different sigma fields.
After the secondary field its variance is delta^2+s. Gaussian
concentration plus log-sum-exp over the 2^K patterns yields

```math
\begin{split}
\mathbb E M(h)&\le Q(A)/\sqrt n+sqrt{2/\pi}\,\delta n
 +\sqrt{2\delta^2nK\log2}+U,\\
\mathbb E\Delta_s(h)&\le
 \sqrt{2/\pi}\,ns/(2\delta)
 +\sqrt{2(\delta^2+s)nK\log2}+2U.                         \tag{6}
\end{split}
```

For the second line, the unperturbed Gaussian JOINT maximum is bounded
below by one fixed-pattern mean G(delta^2); the perturbed joint maximum
is bounded above by G(delta^2+s) plus its displayed pattern overhead.
No independence of pattern maxima is used.

## 4. Concentration and all exponents

Changing one h_ia changes M by at most2r/sqrt(n), and Delta_s by at
most4r/sqrt(n). Bounded differences therefore gives centered MGF
proxies delta^2 n and4delta^2 n, respectively. With probability tending
to one, BOTH(6) hold with an additional n^(3/4) on their right sides;
failure probability is at most2 exp(-c n^(3/4)).

Take s=n^(-1/4) delta. Integer floors cause no difficulty, because
r,K tend to infinity. The exact scales are

```math
\begin{array}{c|c}
q&\Theta(n^{17/24})\\
\delta&\Theta(n^{-1/8})\\
s&\Theta(n^{-3/8})\\
U,\ \sqrt{(\delta^2+s)nK}&O(n^{17/24})\\
ns/\delta&n^{3/4}\\
Q(D)&O(n^{17/16}).
\end{array}
```

Thus sqrt(n) times the first line of(6), plus Q(D), proves the cap
bound in(1). For T=n^(5/4), the b in(3) is n^(3/4)+O(n^(9/16)).
Equations(5)--(6) therefore give w(C_b)<=O(n^(15/16)). The width of
the full parent code is at most the width of its old-coordinate
projection plus E||g_new||_1, the latter O(q). This proves the width
bound in(1), simultaneously for one realization of h.

## 5. Elementary conversion from width to entropy

For any nonempty C subset{+-1}^d, let X be uniform on C, independent
of Gaussian g, and observe Y=X+a g. A nearest-codeword decoder Xhat
satisfies

```math
\|\widehat X-X\|^2\le2a\,g\cdot(\widehat X-X),\qquad
\mathbb E d_H(\widehat X,X)\le a\,w(C)/2.
```

Gaussian channel capacity and coordinatewise binary entropy give,
whenever a w(C)/(2d)<=1/2,

```math
\frac{\log|C|}{d}\le\frac12\log(1+a^{-2})
             +h\left(\frac{a w(C)}{2d}\right).            \tag{7}
```

The capacity bound follows from trace Cov(X)<=d and Gaussian maximum
entropy. The decoding bound remains valid with correlated coordinates
of X and with a nonzero mean. Choosing a=(w(C)/d)^(-1/3) yields
O((w(C)/d)^(2/3)log(d/w(C))) for small relative width. With d=N
and(1)'s width estimate this proves its entropy estimate. No approximate
enumeration of the actual nearcode is assumed.

## 6. Exact implication and remaining gap

Restriction to an unchanged principal block implies Q(A)<=Q(W): average
the new spins, retaining each old word, then use convexity of absolute
value. Hence, when Q(A_n)=O(n^(3/2)), N/n->1 and(1) preserve every
subsequential normalized cap limit. The transformation is valid in
particular on exact minimizers, though its outputs need not be exact
minimizers at their new orders.

This supplies an actual-sign, power-saving geometry regularization
theorem, rather than assuming a small nearcode. It does NOT supply
cheap isotropic response centers, a favorable fixed-eta hierarchy, or
all-order recovery. Sparse random codes can have small Gaussian width
and still obstruct small isotropic absolute response. Also even a pure
Gaussian-field landscape has small-field coordinates whose nearlevel
Hamming radius is of order sqrt(eta/delta), so the center-cover slope
used elsewhere need not be finite. These are genuine remaining gaps,
not reasons to omit the unconditional transformation.

## 7. Dependencies and verification

- Fourth-order all-offset maximum comparison:
  paper_symmetric_frame_universality_2026_09_17.md, Sections1--3.
- Gaussian-field and information mechanisms independently reconstructed
  in paper_localization_external_field_regularization_2026_09_17.md.
- Exact sign gadget, all-query projection, limit orders and exponents:
  reconstructed separately by root, Bernoulli and localization tracks.

Gaussian smoothing and perturb-and-MAP are established mechanisms; the
specific complete-sign compilation is presented as a proved result in
this repository, NOT an established external-priority claim. No change
to the original reported asymptotic interval follows from(1).

## 8. Same-order regularization: no order-selection qualification

The construction also regularizes EVERY signing at its ORIGINAL order.
For a given A_N choose q=ceil(N^(17/24)), m=N-q, and take any
principal m by m submatrix B. Averaging the omitted spins proves
Q(B)<=Q(A_N). Put K=floor(m^(2/3)) and divide q into K positive
integer group sizes r_a differing by at most1. They all satisfy
r_a~m^(1/24). Repeat the independent column h^a exactly r_a times.

Every proof above goes through verbatim with

```math
\delta^2=\frac1m\sum_{a=1}^K r_a^2,\qquad
W_4=\frac1m\sum_{a=1}^K r_a^4.
```

These formulas preserve all exponents. In particular every Gaussian
pattern still has the SAME variance delta^2, even though the integer
group lengths differ. Bounded differences sums their squared lengths;
there is no rounding error paid separately per spin query. A tolerance
N^(5/4) is at most2m^(5/4), harmless in the finite bound(5).

We therefore obtain a full signing Atilde_N of EXACTLY order N with

```math
\begin{split}
Q(\widetilde A_N)&\le Q(A_N)+C N^{11/8},\\
w(E_{\widetilde A_N}(N^{5/4}))&\le C N^{15/16},\\
\log|E_{\widetilde A_N}(N^{5/4})|&\le C N^{23/24}\log N.
\end{split}                                                   \tag{8}
```

Here E_A(T) means the full ABSOLUTE nearlevel set at unnormalized
tolerance T. Only edges incident to q=O(N^(17/24)) designated vertices
are changed, so at most O(N^(41/24)) edges are rewritten. The large
principal block is left exactly unchanged.

Define R_N by the two geometry bounds in(8), with the universal constant
supplied by this proof, independently of the unknown M_N. Then

```math
M_N\le\min_{A\in R_N}Q(A)\le M_N+C N^{11/8}.             \tag{9}
```

Thus optimization may be restricted to this mesoscopically regularized
class at power-saving cost, for EVERY order. This is a constructive
restriction, not a definition referring to the target optimum. We do
not claim it is sufficient for convergence: fixed normalized windows,
all future response laws, and cross-order value transfer remain separate.
An explicit small-width Hadamard subcode still has a universal isotropic
response obstruction, as recorded in the companion audit artifacts.
