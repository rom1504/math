# A uniform correlated improvement over product-cluster extraction

Status: complete finite argument, pending independent audit. This is an
actual-signing pressure/cap inequality, not an original convergence theorem.
It strengthens the product-flip lemma in
`transfer_reconstruction_gibbs_cluster_extraction_2026_09_06.md`; that file's
trial law has independent coordinates, whereas the present law has small
signed correlations on a principal spectral core. No tuning is attempted.

## Statement

Let A be a hollow symmetric sign matrix of order n with
`Q(A)<=C n^(3/2)`. Fix beta>0 and `0<delta<1/2`, and put

```math
r=1-2\delta,\quad z=\Phi^{-1}(\delta),\quad a=\phi(z)^2,
\quad K=8K_G C,\quad
k=\min\{\beta a,1/(2K),1\},
\qquad K_G=\pi/(2\operatorname{arsinh}1).
```

Then, writing `Z_++Z_-=sum_(s,x) exp(s beta H_A(x)/sqrt(n))`,

```math
\log(Z_++Z_-)
\ge n h(\delta)+{\beta r^2Q(A)\over\sqrt n}
 +\left({n\over4}-{1\over2}\right){3\beta ak\over2}
 -\beta C(1-r^2)k^2.                                      (1)
```

In particular the product-cluster lower bound has a uniform strictly
positive pressure correction at every fixed beta and interior delta:

```math
\liminf_n\left\{{1\over n}\log(Z_++Z_-)
 -h(\delta)-\beta r^2{Q(A)\over n^{3/2}}\right\}
\ge {3\over8}\beta a k>0.                                (2)
```

The statement concerns any sequence satisfying the indicated cap bound.
It does not require minimizing A or any spectral bound on the whole matrix.

## Spectral core, with no vertices removed from the trial law

Gauge and reverse the matrix so that `H_B(1)=Q(A)` and `Q(B)=Q(A)`.
The simultaneous diagonal-majorant lemma gives a set I of size
`m>=n/2` with

```math
\|B_I\|_{op}\le8K_G Q(A)/n\le K\sqrt n.                  (3)
```

For its self-contained Grothendieck-rounding and SDP proof see Section 1
of `resumed_bound_audit_minimal_proof_2026_09_06.md`: first
`max_(x,y)|x^TBy|<=4Q(B)`, then a diagonal D majorizes both B and -B
with `Tr D<=4K_G Q(B)`. Retain coordinates
`D_ii<=8K_G Q(B)/n`. This yields (3).

Choose a Gaussian vector G with covariance equal to
`Sigma_I=I_m+k B_I/sqrt(n)` on I, and identity covariance on its complement,
independent between the two parts. Equation (3) and `kK<=1/2` ensure
positive definiteness. Set `Y_i=1` if `G_i>z`, and `Y_i=-1` otherwise.
Every marginal has mean r; all n coordinates remain present.

## Entropy bound

Relative entropy decreases under coordinatewise thresholding. Since the
reference independent Gaussian thresholds to the independent bias-r law,

```math
nh(\delta)-H(Y)
=D(\mathcal L(Y)\Vert\mathcal L(Y_1)^{\otimes n})
\le D(N(0,\Sigma_I)\Vert N(0,I_m))
=-\tfrac12\log\det\Sigma_I.
```

The trace of `B_I` is zero. Applying
`-log(1+u)+u<=u^2/[2(1-kK)]` to its normalized eigenvalues gives

```math
H(Y)\ge nh(\delta)
-{k^2m(m-1)\over4n(1-kK)}.                              (4)
```

## Energy bound and the crucial flat-sign cancellation

Let f=r plus its orthonormal Hermite expansion be the biased threshold
function. Its first coefficient is `2 phi(z)`, and its total nonconstant
squared coefficient mass is `1-r^2`. Consequently for a Gaussian pair
of correlation u,

```math
F(u):=\mathbb E f(G_1)f(G_2)
=r^2+\sum_{j\ge1}c_j^2u^j,\qquad c_1^2=4a.
```

For `v=k/sqrt(n)`, split the series into odd and even parts O(v), E(v).
For `0<=v<=1`,

```math
O(v)\ge4av,\qquad 0\le E(v)\le(1-r^2)v^2.
```

Every core edge has correlation `v B_ij`. Flatness `B_ij=+-1` therefore
gives the EXACT identity

```math
\mathbb E H_B(Y)=r^2Q(A)
 +{m(m-1)\over2}O(v)+H_{B_I}(1)E(v).
```

Outside-core and crossing edges have independent marginals and are already
included in the first term. By principal monotonicity,
`|H_(B_I)(1)|<=Q(A)<=C n^(3/2)`. Thus

```math
{\beta\over\sqrt n}\mathbb E H_B(Y)
\ge {\beta r^2Q(A)\over\sqrt n}
 +{m(m-1)\over n}2\beta a k
 -\beta C(1-r^2)k^2.                                    (5)
```

In particular the even Hermite terms cost only O(1) in log partition,
not O(n): their signed edge sum is cap-controlled. This is where an
absolute edgewise error estimate would incorrectly lose the improvement.

The finite Gibbs variational inequality applied to this trial law, together
with (4)--(5), yields the stronger finite bound with correction

```math
{m(m-1)\over n}
\left[2\beta a k-{k^2\over4(1-kK)}\right]
-\beta C(1-r^2)k^2.
```

Since `k<=beta a` and `kK<=1/2`, the bracket is at least
`(3/2) beta a k`. Also `m(m-1)/n>=n/4-1/2`, proving (1).
The additional condition k<=1 ensures v<=1 at every order.

## Exact cap consequence and scope

Any all-order selected family with normalized log-partition upper bound P
and any fixed cap bound C satisfies

```math
\limsup_n Q(A_n)/n^{3/2}
\le {P-h(\delta)-(3/8)\beta a k\over\beta(1-2\delta)^2}.
```

Thus every fixed interior product-cluster extraction can be strictly
improved for cap-bounded actual signings. A cap bound is available from the
bare extraction before invoking this theorem. This does not by itself show
strict improvement below the infimum of ALL product-corrected parameter
choices: that would additionally require control of their joint boundary
regimes. No such compactness or universality assertion is made here.

Optimizing arbitrary correlated trial laws closes at the exact Gibbs
variational principle; this proof does not supply a finite-dimensional
closure or stationarity theorem that identifies the original minimum.
