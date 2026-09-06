# Arbitrary seed mixers: an exact tilted weave and a Gram-cost obstruction

Date: 2026-09-06. Seed-transfer track. This gives an exact sign realization
for every full sign seed, identifies its new pressure obligation, and
proves that the natural diagonal-temperature Fock bound cannot benefit
from a seed's smaller Boolean cap. The obstruction is to that certificate,
not to the actual ensemble or to convergence.

## 1. Nonorthogonal sign bases are allowed if the row energy is retained

Let `H_i`, `1<=i<=m`, be arbitrary full sign matrices of order `m`, not
necessarily Hadamard. Independently permute their columns. Let `S` be a
symmetric outer full sign matrix with independent fair off-diagonal signs.
Select `k` rows in each fibre, and form the symmetric full sign matrix

```math
K_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i),\qquad N=mk.
```

For a fibre spin `x_i`, put `v_i=H_i[T_i,:]^T x_i/sqrt(k)`. Define

```math
L_t^+(v)=\left(\mathbb E_{g\in G_m}
                     e^{2t\langle v,gv\rangle}\right)^{1/2}
=e^{t\|v\|^2}L_t(v),
\qquad
Z_i^+(t)=\sum_{x_i\in\{\pm1\}^k}L_t^+(v_i).
```

The corresponding positive-semidefinite kernel is
`C_t(a,b)=cosh(2tab)`. The exact extension of the soft weave estimate is

```math
\Pr\{\exists x:\ |x^T Kx|\ge 2cN^{3/2}\}
\le 2(2m)^{m/2}
       e^{-2ctm^2\sqrt{k/m}}\prod_i Z_i^+(t).                (1)
```

To prove this, average the off-diagonal signs in
`exp(t sigma x^T Kx/k)` and condition on the diagonal coordinate chosen
by each column permutation. Each edge gives `cosh(2t v_i(j)v_j(i))`.
Apply the same PSD-kernel graph contraction as in the existing weave
proof. The diagonal term is at most `exp(t v_i(i)^2)`. In the permanent,
the terms fixing that coordinate give

```math
(L_t^+(v))^2\ge\frac{\cosh(2tv_i^2)}m
                      (L_t^+(v\setminus v_i))^2
\ge\frac{e^{2tv_i^2}}{2m}
                      (L_t^+(v\setminus v_i))^2.
```

This pays the polynomial factor in (1). Markov, the spin union bound,
and the two objective signs complete the proof. No constant-norm identity
or spectral upper bound was used.

For independent identically distributed fibre bases, the product in (1)
averages to `(E Z^+)^m`. Thus a sufficient normalized hollow-cap bound is

```math
\limsup_m\frac1m\log\mathbb E Z^+(t)<2ct\sqrt p,
\qquad k/m\longrightarrow p.                              (2)
```

Hollowing pays only `N/2`. When the bases are Hadamard, `||v||^2=m` and
this is exactly the old criterion with the extra constant `t`, not a
different normalization.

## 2. Every actual full sign seed has an exact flat realization

Fix a full sign seed `B` of order `d`, and set `M=B/sqrt(d)`. For a
terminal Hadamard order `q`, choose independent normalized Hadamard bases
`U_j`, output signed permutations `g_j`, and an input signed permutation
`g_0`. Set

```math
V=(M\otimes I_q)\operatorname{diag}(g_1U_1,\ldots,g_dU_d)g_0,
\qquad H=\sqrt{dq}\,V^T.                                  (3)
```

Every entry of `H` is a sign: there is only one nonzero summand for a
given physical input block and output block. The base need not be
symmetric because the outer weave is symmetric by its own formula.
The transpose convention makes the weave row spectrum `H^T x/sqrt(k)`
equal to `V x/sqrt(k/m)`. Its exact input-energy Gram matrix is

```math
V^TV=g_0^T\operatorname{diag}(U_j^Tg_j^T)
        (B^TB/d\otimes I_q)
        \operatorname{diag}(g_jU_j)g_0.                     (4)
```

Thus (1)--(2) apply to an exact ensemble which accepts an ARBITRARY seed,
not only a conference/Hadamard seed. Relatively dense terminal Hadamard
orders give relatively dense parent orders `m floor(pm)`, where `m=dq`
(or `d 2^r q` at fixed child depth). Principal restriction fills all
orders, as in the existing all-order proof.

What remains is the seed-dependent pressure comparison in (2). Its
target for a hollow seed `A`, completed by a sign diagonal to `B`, is
`c <= Q(A)/d^(3/2)+epsilon_d`, with `epsilon_d -> 0` along a near-optimal
seed sequence. The exact realization itself supplies no such comparison.

## 3. A multitemperature nonorthogonal Fock bound

Let `v=(v_1,...,v_d)` with block length `q`, `m=dq`, and
`||v||^2<=Cm`. Let the independent `g_j` be uniform signed permutations.
For fixed `d,t,C,M` choose positive scalar temperatures `tau_j` such that

```math
\operatorname{diag}(\tau_1,\ldots,\tau_d)\succeq tM^TM.       (5)
```

Then, uniformly in `v`,

```math
\mathbb E_g L_t^+((M\otimes I_q)gv)
\le\exp(O_{d,t,C,M,\tau}(\sqrt m))
                           \prod_{j=1}^d L_{\tau_j}^+(v_j). (6)
```

Here is a complete operator proof. Use the UNNORMALIZED symmetric-Fock
feature

```math
\phi^+(z)=\bigoplus_{r\ge0}\frac{2^{r/2}}{\sqrt{r!}}z^{\otimes r},
\qquad\langle\phi^+(z),\phi^+(w)\rangle=e^{2\langle z,w\rangle}.
```

Write `D=diag(sqrt(tau_j)) otimes I_q` and
`C_0=sqrt(t)(M otimes I_q)D^(-1)`. By (5), `||C_0||op<=1`; its
second quantization is contractive on every degree. Since `D` commutes
with the product signed-permutation group,
`phi^+(sqrt(t)Mgv)=Gamma(C_0)phi^+(gDv)`.

The orbit covariance of `phi^+(gDv)` is a tensor product over blocks.
Each factor has operator norm exactly `(L_(tau_j)^+(v_j))^2`: its finite
orbit Gram matrix has positive entries and constant row sum. Hence the
product covariance has norm `prod_j(L_(tau_j)^+(v_j))^2`.

Project to the global `G_m` invariant subspace through degree `R=O(m)`.
Its rank is `exp(O(sqrt(m)))` by the partition bound. The trace estimate
for the projected contractive image of the covariance is at most this
rank times its operator norm. Choose the degree cutoff sufficiently
large that the unnormalized exponential-series tail, with parameter
`2||Dv||^2=O(m)`, is at most one. This is possible uniformly with a
linear cutoff by the same Poisson Chernoff estimate as before. The
discarded norm remains at most one after the contraction. Finally
`L_tau^+(w)>=1`, from the constant kernel term. Taking square roots and
Jensen proves (6).

The zero-temperature or zero-diagonal endpoint, if desired, is obtained
by a decreasing positive regularization of the `tau_j`; it is not
needed for fixed positive `t` and full sign `B`.

## 4. Why this certificate cannot exploit a better Boolean seed

Define the tilted finite-type boundary

```math
\Psi_t(\nu)=t\,m_2(\nu)+\Phi_t(\nu).
```

It is increasing in `t`. Indeed `F_t` is an infimum of affine functions
of `t`, hence concave, with right derivative at zero `2 Var(nu)`.
For symmetric `nu`, `Var(nu)=m_2(nu)`, so `Psi_t` is increasing.
Equivalently this follows directly from the positive, increasing
cosh-kernel permanent. The same monotonicity holds for every fixed-depth
tilted orthogonal recursion
`t m_2(nu)+(B^r Phi_t)(nu)`: the orthogonal recursion preserves the
average second moment, and each terminal tilted reward is increasing.

For a full sign seed, every diagonal entry of `M^TM` is exactly one.
Consequently every feasible temperature in (5) obeys `tau_j>=t`.
Thus (6), followed by any fixed-depth child orbital/type bounds, is
NEVER sharper than the corresponding bound with all child temperatures
equal to `t` and an orthogonal flat mixer. This holds termwise for every
input split/type, and therefore also after taking their supremum and
paying their exact type entropies.

The seed's only role in this certificate is an extra positive Gram
payment. Its smaller `Q(B)` cannot lower the bound. Optimizing the
diagonal temperatures does not fix this: all of them have the same
lower bound `t`, and the normalized target denominator remains the
ORIGINAL `t` in (2).

This is a scalable obstruction to the diagonal-temperature Fock
contraction certificate. It does not establish the pressure of the
actual ensemble, and it does not prohibit a covariance-sensitive bound
which retains more than the orbit covariance's largest eigenvalue.
In particular, replacing the variable row energy by its spectral
maximum is a strictly coarser version of the same obstruction.

## 5. Boolean near-optimality does not give uniform Gram control

There is a simple additional warning about universal seed claims.
Starting from ANY near-optimal hollow signing `A_n`, replace the edges
inside a set of `s_n=floor(n^(3/4)/log n)` vertices by all plus signs,
and call the result `B_n`. Then

```math
|Q(B_n)-Q(A_n)|\le s_n(s_n-1)=o(n^{3/2}).                  (7)
```

So the modified sequence remains asymptotically near-optimal. Its
planted principal clique gives `||B_n||op>=s_n-1`, hence
`||B_n||op/sqrt(n) -> infinity`. More specifically, let `z` equal one
on the clique and zero elsewhere. Averaging independent fair spins
outside the clique gives

```math
\max_{x\in\{\pm1\}^n}\|B_nx\|^2
\ge\|B_nz\|^2\ge s_n(s_n-1)^2.
```

Therefore `max_x||B_nx||^2/n^2 -> infinity` as well. This proves that
neither a bounded normalized operator norm nor a bounded uniform output
second moment follows from the Boolean cap of an arbitrary near-optimal
seed. It does not exclude selecting a different well-conditioned
near-optimal sequence. Nor does it obstruct (1), which deliberately
keeps the fluctuating row energy inside its tilted partition sum.
