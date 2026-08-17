# Fractional balance forced by near-minimality

**Status:** proof draft for independent audit.  This concerns genuine
near-minimizers at arbitrary orders.  It is not yet a composable response
state and does not imply convergence.

Let `E=binom(n,2)`, identify a hollow signing with
`a in {+-1}^E`, and let

```math
Z_n^+=\{(sigma x_ix_j)_(i<j):sigma,x_i in\{+-1\}\}.
```

Then `Q(a)=max_(z in Z_n^+)<a,z>`.  For `u>0`, define the positive near-top
shell

```math
S_u(a)=\{z in Z_n^+:<a,z>\ge Q(a)-u n^{3/2}\}.              \tag{FB.1}
```

For a probability measure `mu` on a shell put

```math
m_e(mu)=E_mu[a_ez_e],
\qquad
V(mu)={1\over E}\sum_e(m_e(mu))_+.                         \tag{FB.2}
```

Thus `V` measures how much common sign-aligned edge mass remains after
mixing near-top signed cuts.

## Theorem FB.1 (near-minimizers have a fractionally balanced active shell)

There is an absolute constant `C` with the following property.  Let

```math
Q(a)\le M_n+epsilon n^{3/2},
\qquad 0<kappa<1/4,
\qquad q={kappa\over\sqrt n},                                \tag{FB.3}
```

and put

```math
eta_n(kappa)=C\left(\sqrt{kappa}\,n^{-1/4}+n^{-1/2}\right).
                                                                    \tag{FB.4}
```

If `epsilon+eta_n(kappa)<kappa/2`, then there is a probability measure `mu`
supported on `S_(2kappa)(a)` such that

```math
\boxed{
V(mu)\le
 {epsilon+eta_n(kappa)\over kappa(1-1/n)}.}                  \tag{FB.5}
```

If also `Q(a)>2kappa n^(3/2)`, then the same measure obeys

```math
{1\over E}\sum_e|m_e(mu)|
\le {2(epsilon+eta_n(kappa))\over kappa(1-1/n)}.             \tag{FB.6}
```

Consequently, for fixed sufficiently small `epsilon` and then `n ->
infinity`, choosing `kappa=sqrt(epsilon)` gives a shell of normalized width
`2sqrt(epsilon)` and total signed-edge bias `O(sqrt(epsilon))`; the known
positive lower bound on `M_n/n^(3/2)` supplies the additional hypothesis in
(FB.6).  For exact
minimizers, the explicit choice `kappa=n^(-1/6)` gives both shell width and
normalized total bias `O(n^(-1/6))`.

### Proof

For a nonempty finite shell `S`, finite minimax gives

```math
\begin{split}
v(S)
&=\max_{w in[0,1]^E}\min_{z in S}{1\over E}
       \sum_e w_ea_ez_e\\
&=\min_{mu in Delta(S)}{1\over E}
       \sum_e\big(E_mu[a_ez_e]\big)_+ .                    \tag{FB.7}
\end{split}
```

Indeed, after exchanging max and min, the maximizing coordinate choice is
`w_e=1` on positive mean coordinates and zero on negative ones.  Apply this
to `S=S_(2kappa)(a)` and choose a maximizing `w`.

Independently flip signing coordinate `e` with probability `qw_e`, obtaining
an exact signing `a^F`.  For every `z in Z_n^+`,

```math
<a^F,z>=<a,z>-2\sum_e xi_ea_ez_e,
\qquad
E<a^F,z>=<a,z>-2q\sum_ew_ea_ez_e.             \tag{FB.8}
```

There are at most `2^(n+1)` signed cuts.  Bernstein's inequality and a union
bound therefore give, with positive probability,

```math
\max_{z in Z_n^+}|<a^F,z>-E<a^F,z>|
\le C\left(\sqrt{kappa}\,n^{5/4}+n\right)
=eta_n(kappa)n^{3/2}.                            \tag{FB.9}
```

For `z in S`, (FB.7)--(FB.8) give

```math
E<a^F,z>
\le Q(a)-2qEv(S)
=Q(a)-kappa(1-1/n)v(S)n^{3/2}.                  \tag{FB.10}
```

For `z notin S`, use `sum_e w_ea_ez_e>=-E` to obtain

```math
E<a^F,z>
<Q(a)-2kappa n^{3/2}+2qE
\le Q(a)-kappa n^{3/2}.                         \tag{FB.11}
```

Equations (FB.9)--(FB.11) show that if the right side of (FB.5) were
strictly smaller than `v(S)`, then

```math
Q(a^F)<Q(a)-epsilon n^{3/2}\le M_n,
```

contradicting the definition of `M_n`.  This proves (FB.5).

For a minimizing measure in (FB.7), write
`P=sum_e(m_e)_+` and `N=sum_e(-m_e)_+`.  Since the measure is supported on
the shell,

```math
P-N=\sum_em_e=E_mu<a,z>\ge Q(a)-2kappa n^{3/2}>0.             \tag{FB.12}
```

Hence `N<=P`, so `sum_e|m_e|=P+N<=2P`; (FB.6) follows from
(FB.5). `square`

### Corollary FB.2 (a short fractional witness list)

If a measure in Theorem FB.1 satisfies

```math
{1\over E}\sum_e|m_e(mu)|\le delta,
```

then for every integer `K>=1` there are `K` (not necessarily distinct)
signed cuts `z_1,...,z_K` in the same shell such that

```math
{1\over E}\sum_e\left|{1\over K}\sum_(r=1)^K a_e(z_r)_e\right|
\le delta+K^{-1/2}.                              \tag{FB.13}
```

Indeed sample the cuts independently from `mu`.  For each edge, Jensen and
the variance bound give expected absolute empirical error at most
`K^(-1/2)`; summing over edges proves existence.  Thus fixed normalized
near-minimizer accuracy has a list of `O(delta^(-2))` near-top witnesses
whose average cancels almost every edge.  The list still stores
`O(Kn)` Boolean bits and is not by itself a subextensive response carrier.

## Theorem FB.3 (exact minimizers have a stretched-exponential near-top shell)

Assume `M_n>=c_0n^(3/2)` eventually for some fixed `c_0>0`.  There is
`c_1=c_1(c_0)>0` such that the following holds for all sufficiently large
`n`.  If `a` is an exact minimizer and `r=r_n` is an integer sequence with

```math
n\le r,
\qquad {r\over n^{3/2}}\longrightarrow0,
```

then

```math
\boxed{
 |S_(2r/n^{3/2})(a)|\ge\exp(c_1r/n).}           \tag{FB.14}
```

In particular, for every `L_n -> infinity` with `L_n<=sqrt(n)`, take
`r=floor(n^(3/2)/L_n)` and decrease `c_1` by an absolute factor to absorb
rounding.  Every exact minimizer has at least

```math
\exp\left(c_1{\sqrt n\over L_n}\right)
```

signed cuts of energy at least
`M_n-2n^(3/2)/L_n`.  Taking `L_n=log n` gives a vanishing-width shell with
logarithmic cardinality `Omega(sqrt(n)/log n)`.

### Proof

For `z in S_(2r/n^(3/2))(a)`, put `g_e(z)=a_ez_e`.  Since
`M_n>=c_0n^(3/2)` and `r=o(n^(3/2))`, after decreasing `c_0` by an absolute
factor if necessary,

```math
{1\over E}\sum_eg_e(z)
={<a,z>\over E}\ge {c_0\over\sqrt n}             \tag{FB.15}
```

uniformly on the shell.  Choose an `r`-element edge set `F` uniformly.
Hoeffding--Serfling for sampling without replacement from the sign population
`(g_e(z))_e` gives

```math
\Pr\left\{\sum_(e in F)g_e(z)\le0\right\}
\le\exp(-c_1r/n).                               \tag{FB.16}
```

If the shell had fewer than `exp(c_1r/n)` members (with `c_1` decreased once
more), a union bound would produce one `F` for which every shell member has
positive sum.

Flip exactly the edges in this `F`.  Since `a` minimizes over all signings,
there is a signed cut `z_*` with

```math
M_n=Q(a)\le<a^F,z_*>
=Q(a)-d_a(z_*)-2\sum_(e in F)a_e(z_*)_e,        \tag{FB.17}
```

where `d_a(z)=Q(a)-<a,z>`.  Hence

```math
d_a(z_*)+2\sum_(e in F)a_e(z_*)_e\le0.          \tag{FB.18}
```

The sum is at least `-r`, so `d_a(z_*)<=2r`; thus `z_*` belongs to the shell.
But (FB.18) also says its sum on `F` is nonpositive, contradicting the choice
of `F`. `square`

The proof is a genuine exact-minimizer consequence rather than finite
evidence.  Cardinality alone is not contextual packing: the near-top cuts
could still be generated by a small algebra or symmetry orbit.  The theorem
rules out a bounded list that explicitly enumerates the whole shell, but not
a bounded algebraic, symmetric, or generative description.

## Interpretation and limitations

1. This is a **PROVES AN ARROW** result at benchmark Level 5: genuine
   near-minimality forces a collective active-shell equilibrium at arbitrary
   orders.
2. The conclusion is weaker than the original optimization.  It asserts the
   existence of one mixture with `E` first marginals; it neither gives the
   maximum, the full energy histogram, nor a maximizing-spin table.
3. It is not yet the desired small reusable state.  The measure may require
   many atoms, and first edge marginals need not determine response under a
   dense continuation.
4. The theorem explains a necessary repair mechanism.  A sparse direction
   cannot lower every near-top signed cut by a leading amount; otherwise its
   exact random realization would beat `M_n`.
5. The same proof applies to any `exp(O(n))` family of sign queries, so it is
   not tied to Walsh or conference structure.
6. Exact minimizers necessarily have a vanishing-width shell containing
   `exp(n^(1/2-o(1)))` signed cuts, but this is a count, not a response-metric
   packing or an entropy-rate lower bound under composition.

The next discriminating question is whether augmented-cut geometry turns
the fractional balance in (FB.6) into a low-information collective response
state, or whether exact minimizers can realize (FB.6) with an extensive
contextual packing.
