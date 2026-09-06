# Independent reconstruction of the improved all-order signing upper bound

Date: 2026-09-06. Status: verified consequence of the new director precision
supersolution, the independently replayed exact E certificate, and the
existing reconstructed orbital/type/weave realization. Entropy extraction
is omitted: it is unnecessary for this minimal proof.

Let `M_n=min_A max_x |sum_(i<j) a_ij x_i x_j|` over hollow signings.
Put

```math
p={31\over32},\quad t=4,\quad
a_E={19678127864847\over800000000000000}.
```

Then

```math
\limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}-{a_E\over8\sqrt{31/32}}
<0.496876095.
```

This is an all-order upper bound on the ORIGINAL signing problem. It does
not prove convergence or identify actual minimizing signings with the
recursive construction.

## 1. Analytic certificate and limit order

The new precision-side Schur proof gives `B E_t<=E_t`; conditional copies
give the opposite inequality. The independently audited Gaussian-boundary
replacement therefore identifies

```math
\lim_r B^r\Phi_t(\nu)=E_t(\nu).
```

The old exact finite-channel certificate, replayed in this campaign,
gives

```math
p\log2+E_4(\nu_p)+4(1-\sqrt p)\le-a_E.
```

Fix any a' with 0<a'<a_E. By the decreasing Bellman limit there is a
FINITE depth r such that the same expression with `B^r Phi_4` is at most
-a'. Next fix `0<eta<a'/4`. These parameters remain fixed throughout the
matrix-order limit. No effective depth bound or growing-depth type theorem
is needed.

## 2. One-row orbital estimate and independent fibre assembly

For m divisible by 2^r, build the exact depth-r randomized Hadamard bases
from binary Hadamard gates, independent child bases, and fresh uniform
signed input permutations. Every resulting U_m is orthogonal with entries
of absolute value 1/sqrt(m). Terminal Hadamards may be arbitrary of their
specified order.

Set `k=floor(pm)` and use the ternary input vector with k entries
`+-sqrt(m/k)` and the remaining entries zero. Its symmetrized empirical
law tends to nu_p. The finite-depth type bound and uniform orbital bound
give

```math
\limsup_m{1\over m}\log\mathbb E_H Z_T(4)
\le p\log2+(B^r\Phi_4)(\nu_p).
```

Here the physical fibre Hadamard is `H=sqrt(m) U_m^T`, so its normalized
row spectrum is EXACTLY `U_m(1_T x)/sqrt(k/m)`. The extra diagonal deletion
in the orbital permanent norm costs at most a factor sqrt(2m), negligible
on the m scale. Fixed-depth factorial and type errors are o(m), uniformly
over the terminal Hadamards; root and reachable alphabets remain finite.

Choose a fresh independent recursive base in every fibre. Independently
sample the outer symmetric edge signs and output-column permutations.
The PSD-kernel graph-contraction bound then factors as `(E Z_T)^m`, not
`E[Z_T^m]`. This independence is essential.

## 3. Exact defect, event bound, and cap normalization

The restricted weave is the symmetric FULL sign matrix K of order N=mk
with entries

```math
K_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i),\quad a\in T_i,b\in T_j.
```

For a spin vector x, let `h_i=H_i[T_i,:]^T x_i`. Orthogonality gives
`sum_i ||h_i||²=m²k`, while the quadratic form is
`x^T Kx=sum_(i,j) S_ij h_i(j)h_j(i)`. For each sigma=+-1,

```math
D_\sigma(x)=\sum_{i,j}(h_i(j)-\sigma S_{ij}h_j(i))^2
=2(m^2k-\sigma x^TKx).
```

With `gamma_m=1-sqrt(k/m)+eta`, the exponential Markov bound and complete
spin/energy-sign union bound give

```math
\Pr\{\exists x,\sigma:D_\sigma(x)\le2\gamma_m m^2k\}
\le2\exp(4\gamma_m m^2)(\mathbb E Z_T(4))^m.
```

Its logarithm divided by m² has limsup at most `-a'+4eta<0`. Hence an
actual K exists outside this event for every sufficiently large allowed m.
It satisfies

```math
\max_x|x^TKx|\le m^2k(\sqrt{k/m}-\eta).
```

Delete only the diagonal to obtain a hollow signing A. Its Hamiltonian is
one half of its matrix quadratic form, and diagonal deletion costs at most
N/2. Thus

```math
{Q(A)\over N^{3/2}}
\le{1\over2}-{\eta\over2\sqrt{k/m}}+{1\over2\sqrt N}.
```

The factor 1/2 and the square root of the retention fraction have both
been checked directly from N=mk.

## 4. All orders, not just powers of two

Use the explicit Hadamards H_2 and H_12. The latter follows, for example,
from the elementary Paley skew matrix over the eleven-element field:
`R^T=-R`, `RR^T=11I`, so H_12=I+R. Kronecker products supply every order
`s=2^a 12^b`, a,b nonnegative integers. The terminal matrix need not be
symmetric.

These orders are multiplicatively relatively dense. The ratio
log(12)/log(2) is irrational, since otherwise some positive powers of
12 and 2 would agree, contradicting unique factorization. For any
epsilon>0, a FINITE set of nonnegative multiples of log(12) is an
epsilon-net modulo log(2). This is the elementary irrational-rotation
density statement, provable by pigeonholing to obtain a sufficiently
small nonzero rotation step and taking a finite number of its iterates.
For any sufficiently large target log X, choose one of those finitely
many b values and then a>=0 so that
`log X<=a log2+b log12<=log X+epsilon`. The lower bound on X ensures
that the integer a is nonnegative. Thus the ordered available s_j satisfy
`s_(j+1)/s_j->1`.

No prime number theorem is needed. I separately checked the alternative
fixed-progression PNT argument in the earlier artifact, but it is not a
dependency of this minimal reconstruction.

Take `m_j=2^r s_j`, `N_j=m_j floor(p m_j)`. Then
`N_(j+1)/N_j->1`. The uniform terminal theorem makes the construction valid
at every sufficiently large N_j.

For arbitrary n, choose the least N_j>=n and take a principal n-by-n
restriction. Its cap does not increase: average all removed spins
independently to recover any prescribed restricted Hamiltonian. Since
N_j/n->1, this proves the all-order limsup bound
`1/2-eta/(2sqrt p)`. Finally send eta up to a'/4 and a' up to a_E.
These are limits AFTER the fixed-depth matrix-order limit.

## 5. Exact constant check and scope warning

The exact replay outputs and constant conversion are

- `computations/decisive_bridge_supersolution_ternary_E_replay_2026_09_06.json`;
- `computations/decisive_bridge_improved_cap_constant_2026_09_06.py` and JSON.

The outward upper square-root endpoint is 246062746063/250000000000.
It yields the rational upper bound

```math
{1\over2}-{a_E\over8\sqrt p}
\le {3129925021741553\over6299206299212800}<0.496876095.
```

The lower/upper interval displays for the constant are approximately
0.4968760940775465 and 0.49687609407754973. No floating decimal enters
the certificate.

**Important distinction:** H=E is an exact identification of the
GAUSSIAN-BOUNDARY BELLMAN CERTIFICATE. The actual one-row annealed pressure
is only bounded ABOVE by the Bellman iterates; concatenation/projection
and terminal orbital inequalities may have slack. The outer Finner and
union bounds may also have slack. No equality with actual ensemble
pressure, best ensemble cap, or original minimizing cap has been proved.
