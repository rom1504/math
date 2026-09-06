# Exact Gibbs-cluster extraction for both signs

Date: 2026-09-06. This note changes the extraction of a quadratic cap from
an annealed partition upper bound. It does not change the row potential,
assert a second-moment theorem, or invalidate the uncorrected-partition
floors proved elsewhere. All logarithms are natural.

## 1. Finite lemma, without an imported variational theorem

Let `A` be any real symmetric hollow matrix of order `N`, and put

```math
H_A(x)=\frac12x^TAx,\qquad Q(A)=\max_x|H_A(x)|,
\qquad Z_\sigma(b)=\sum_{x\in\{-1,1\}^N}e^{\sigma bH_A(x)}.
```

For `0<=delta<1/2`, set `q=1-2delta` and
`h(delta)=-delta log delta-(1-delta)log(1-delta)`. Then, for every `b>0`,

```math
Z_+(b)+Z_-(b)
\ge 2e^{Nh(\delta)}\cosh(bq^2Q(A))
\ge e^{Nh(\delta)+bq^2Q(A)}.                 \tag{1}
```

Proof. Choose `x_*` with `|H_A(x_*)|=Q(A)`. Let independent signs
`xi_i` have probability `delta` of being negative, and put
`Y_i=(x_*)_i xi_i`. Its product law `mu` has entropy `Nh(delta)`.
Hollowness gives the exact identity

```math
\mathbb E H_A(Y)=q^2 H_A(x_*).               \tag{2}
```

For `0<delta<1/2`, Jensen's inequality applied under `mu` gives

```math
\log Z_\sigma(b)
=\log\mathbb E_\mu\frac{e^{\sigma bH_A(Y)}}{\mu(Y)}
\ge \sigma bq^2H_A(x_*)+Nh(\delta).
```

Sum the two resulting inequalities. The endpoint `delta=0` follows by
continuity or by taking the two terms at `x_*`. No local-maximality,
positivity of the matrix, or distributional assumption is used. In
particular, the same product law may be used for the two objective signs.

Consequently,

```math
Q(A)\le
\frac{\log(Z_+(b)+Z_-(b))-Nh(\delta)}{b(1-2\delta)^2}. \tag{3}
```

If a probability distribution on hollow sign matrices satisfies
`E[Z_+(b)+Z_-(b)]<=B`, there exists an actual matrix in its support with
`Z_+(b)+Z_-(b)<=B`. Equation (3) therefore gives an existence bound with
`log B`. Finite ensembles suffice below, so no support or essential-infimum
qualification is needed. A union bound over extremizers is unnecessary.

## 2. Exact weave normalization and diagonal cost

Use the notation of Sections 5--7 of
`transfer_reconstruction_standalone_2026_09_06.md`. The full symmetric
weave `W=W_T` has order `N=mk`, and its hollow signing is
`A=W-diag(W)`. Set `D_sigma(x)=2(m^2 k-sigma x^TWx)`. For every realization,

```math
\sum_xe^{\sigma(t/k)x^TAx}
=e^{tm^2-\sigma(t/k)\operatorname{tr}W}
  \sum_x e^{-tD_\sigma(x)/(2k)}.             \tag{4}
```

Here `|tr W|<=mk`, so the diagonal factor costs at most `e^{tm}`,
not `e^{tm^2}`. Notice that the inverse temperature in (1) is
`b=2t/k`, because `H_A=x^TAx/2`.

The all-spin Finner calculation before Markov's inequality gives, for
each objective sign,

```math
\mathbb E\sum_xe^{-tD_\sigma(x)/(2k)}
\le(\mathbb EZ_{\mathrm{row}})^m.
```

Fibre matrices are independent; this is a power of an expectation, not
the expectation of a power. Therefore

```math
\mathbb E[Z_+(2t/k)+Z_-(2t/k)]
\le 2e^{tm^2+tm}(\mathbb EZ_{\mathrm{row}})^m. \tag{5}
```

At any fixed recursion depth `r`, the established finite-type estimate is

```math
\limsup_{m\to\infty}\frac1m\log\mathbb EZ_{\mathrm{row}}
\le p\log2+(\mathcal B^r\Phi_t)(\nu_p),
\qquad k/m\longrightarrow p.
```

Define `A_r=p log2+t+(B^r Phi_t)(nu_p)`. Select a realization by (5)
and apply (3). Since

```math
\frac{(2t/k)(mk)^{3/2}}{m^2}=2t\sqrt{k/m},
```

the resulting bound along the admissible orders is

```math
\limsup\frac{Q(A)}{N^{3/2}}
\le \frac{A_r-p h(\delta)}{2t\sqrt p(1-2\delta)^2}. \tag{6}
```

This checks the two potentially hidden factors of two, and incorporates
the physical diagonal before invoking the hollow quadratic lemma.
The factors `2` and `e^{tm}` in (5) contribute `o(m^2)` only.

## 3. Fixed-depth quantifiers and the existing exact certificate

The already reconstructed conditional-variance certificate gives

```math
p=\frac{31}{32},\quad t=4,\quad
a=\frac{91470529542342299}{20460000000000000000}>0,
```

and, for every fixed `a'<a`, a finite depth `r` such that
`A_r<=t sqrt p-a'`. Fix `delta` first, then `a'`, then such a depth.
Take the order limit at this fixed depth. The admissible terminal orders
`2^u12^v` and principal restriction, exactly as in Section 7 of the
standalone reconstruction, extend (6) to every sufficiently large order.
Only then let `a'` increase to `a`. Thus

```math
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le C_\delta:=(1-2\delta)^{-2}
 \left[\frac12-\frac{a+p h(\delta)}{8\sqrt p}\right]. \tag{7}
```

There is no estimate uniform in a growing recursion depth or alphabet,
and no exchange of the finite-depth and order limits.

Let `C_0=1/2-a/(8 sqrt p)`. The improvement is equivalent to

```math
p h(\delta)>4\delta(1-\delta)(4\sqrt p-a).   \tag{8}
```

Since `h(delta)/delta` diverges as `delta` decreases to zero, every fixed
finite pressure upper bound has a strictly improved extraction at some
sufficiently small fixed positive `delta`. For the explicit rational
choice `delta=2^{-24}`, the accompanying outward-rational checker proves

```math
C_\delta<0.499432211<C_0.
```

This modest numerical gain comes from a different counting implication:
an extremizer forces an entire positive-entropy Gibbs-weighted cluster.
The permanent and typical-profile floors still apply to the uncorrected
partition criterion; they do not apply to (3) with its entropy credit.
The independently falsified near-identical Gaussian replica inequality
is consistent with, but is not needed by, this finite proof.

## 4. Reproduction and precise scope

Run

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/transfer_reconstruction_gibbs_cluster_exact_checks_2026_09_06.py
```

The checker uses outward rational logarithm and square-root intervals,
checks (8) and (7) at the fixed certificate parameters, and enumerates
every hollow sign matrix of orders two through four. At three rational
flip probabilities it verifies the exact entropy/energy identity and
the both-sign lower bound at inverse temperature `log 2`; all partition
sums in this finite test are rational. These tests complement, not
replace, the Jensen proof. The analytic closure and original root
certificate are dependencies of (7), not of the finite lemma (1).
