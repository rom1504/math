# Wave 23: a regular-parent-ground codebook improves the rate exponent

## Statement

Use the ordered pairing of the ledger, let `A` be an exact order-`n`
minimizer with `Q(A)=q_n`, and keep a fixed selector density bounded away
from zero and one.  For a full oriented cut `d=(sigma,x)`, write

```math
w_{ij}=\sigma a_{ij}x_ix_j,
\qquad
r_i(d)=\sum_{j\ne i}w_{ij}.
```

Suppose that, for some `0<c<1/4`, a finite codebook
`\mathcal C\subseteq\mathcal G(A)` of exact **parent** grounds satisfies

```math
\begin{aligned}
\mathbb E_{S\sim U_m}
\left[Q(A[S])-\max_{d\in\mathcal C}c_A(S,d)\right]
&=O(n^{3/2-c}),\\
\log|\mathcal C|&=O(n^{3/4-c}),\\
\max_{d\in\mathcal C}\max_i r_i(d)&=O(n^{3/4-c}).
\end{aligned}
\tag{R23.1}
```

Then

```math
\boxed{
V_{\rm ad}(A,m)-p_2q_n=O(n^{3/2-c}).
}
\tag{R23.2}
```

Thus (R23.1) is another sufficient form of the power-saving restriction
lemma.  Compared with (10.714), it permits codebook log-size
`O(n^{3/4-c})` rather than `O(n^{1/2-2c})`, at the price of supporting the
codebook on row-regular parent grounds.

## Proof

Every parent ground has nonnegative row fields.  Indeed, flipping spin `i`
changes its ordered payoff from `q_n` to `q_n-4r_i`; ground optimality gives
`r_i\ge0`.  Moreover,

```math
\sum_i r_i=q_n,
\qquad
R_2(d):=\sum_i r_i^2
\le R_\infty(d)\sum_i r_i
=R_\infty(d)q_n.
\tag{R23.3}
```

The known universal upper bound `q_n=O(n^{3/2})` and (R23.1) therefore give

```math
R_\infty=O(n^{3/4-c}),
\qquad
R_2=O(n^{9/4-c})
\tag{R23.4}
```

uniformly on the codebook.

Choose a maximizing codeword `D=D(S)` with fixed tie breaking, and let `nu`
be its output law.  Before taking the supremum over all cuts in (10.671),
the entropy inequality gives the sharper support-dependent form

```math
\mathbb E[c_A(S,D)-p_2\langle A,D\rangle]
\le
\frac{I(S;D)+
\log\mathbb E_{d\sim\nu}e^{\Lambda_d^U(\lambda)}}{\lambda}.
\tag{R23.5}
```

Here every output is a parent ground, so
`\langle A,D\rangle=q_n`, and `I(S;D)\le\log|\mathcal C|`.
Equation (10.691), restricted to the support and using (R23.4), gives

```math
\log\mathbb E_\nu e^{\Lambda_d^U(\lambda)}
\le
\chi_{n,m}+\lambda\epsilon_{n,m}q_n
+O\!\left(\lambda^2
[n^{9/4-c}+n^2]\right)
\tag{R23.6}
```

provided the two denominators stay bounded away from zero.

The exact-minimizer spectral bound is
`\lVert A\rVert_{\rm op}^2\le2q_n=O(n^{3/2})`.  Choose

```math
\lambda=\eta n^{-3/4}
\tag{R23.7}
```

with a sufficiently small fixed `eta`.  The operator denominator in
(10.691) is then uniformly positive, while
`\lambda R_\infty=O(n^{-c})`.  Substitution of (R23.6)--(R23.7) into
(R23.5), followed by the distortion bound, yields

```math
\begin{aligned}
V_{\rm ad}(A,m)-p_2q_n
\le{}&
O(n^{3/2-c})+\epsilon_{n,m}q_n\\
&+O\!\left(
n^{3/4}\log|\mathcal C|
+n^{3/4}\chi_{n,m}
+n^{-3/4}[n^{9/4-c}+n^2]
\right).
\end{aligned}
\tag{R23.8}
```

At fixed density,
`\epsilon_{n,m}q_n=O(n^{1/2})`,
`\chi_{n,m}=O(\log n)`, and the last row is

```math
O(n^{3/2-c})+O(n^{3/4}\log n)+O(n^{5/4}).
```

For `c<1/4` this proves (R23.2).

## Interpretation and scope

The stronger entropy allowance comes from stopping at the
operator-controlled boundary `lambda=Theta(n^{-3/4})` rather than using the
uniform worst-row optimizer `lambda=Theta(n^{-1-c})`.  Parent-ground local
optimality converts one row-field cap into the required square-field bound.

This does not prove that such a codebook exists.  Its exact remaining
combinatorial content is:

1. cover most optimized principal sections to deficit `O(n^{3/2-c})` using
   at most `exp(O(n^{3/4-c}))` parent grounds; and
2. choose those grounds with maximum row field
   `O(n^{3/4-c})`.

The second condition is much weaker than square-root row regularity but is
not automatic: `sum_i r_i=O(n^{3/2})` permits a migrating set of heavier
coordinates.  The edge-cube results in Waves 22--23 explain why unweighted
flip witnesses alone do not rule that out.

There is also a critical exponent barrier to obtaining the cap by elementary
trimming.  For a parent ground and threshold `T`, the heavy set

```math
H_T=\{i:r_i>T\}
```

has only `|H_T|\le q_n/T`.  Conditioning a fixed-density selector to include
or exclude this known set costs `O(|H_T|)` relative entropy, and leaves
light-coordinate coefficients bounded only at the scale
`T+|H_T|` without further cancellation.  The entropy term at the
corresponding Bernstein boundary is therefore controlled only by

```math
|H_T|\,[T+|H_T|]
\le \frac{q_n}{T}\left(T+\frac{q_n}{T}\right).
\tag{R23.9}
```

Balancing at `T\asymp\sqrt{q_n}\asymp n^{3/4}` reproduces the full
`O(n^{3/2})` scale, with no power saving.  This is a limitation of the
Markov/conditioning certificate, not a lower bound on the actual best
ground.  A proof of (R23.1) must obtain genuine row regularity or exploit
cancellation beyond the total-field identity.
