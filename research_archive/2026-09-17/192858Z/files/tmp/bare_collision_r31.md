# Wave 31 Route 2: bare mesoscopic collision moments

Status: exact abstract incidence results, with a scoped extremal wall.  None
of the statements below supplies the still-missing exact-minimizer input.

## 1. Setup and the diagonal moment

Let `X` be the selector set, let `A` be a family of `M` candidate hit sets
`H_a subset X`, and let `w` be an arbitrary selector law.  Write

```math
p_a=w(H_a),\qquad h_w=\max_a p_a,\qquad
S_s(w)=\sum_a p_a^s.
```

For an iid batch `X_1,...,X_r`, put

```math
N_a=\sum_{i=1}^r\mathbf1_{\{X_i\in H_a\}},\qquad
J_{r,s}=\{\max_aN_a\ge s\}.
```

Thus `J_(r,r)` is the common-coset event of (10.889).  Define the
nonnegative integer

```math
Z_{r,s}=\sum_a\binom{N_a}{s}.
```

Expanding over the `s`-subsets of batch coordinates gives the exact first
moment

```math
\mathbb E_w Z_{r,s}=\binom rs S_s(w).
```

Moreover `J_(r,s)={Z_(r,s)>0}` and, on this event,
`1<=Z_(r,s)<=M binom(r,s)`.  Consequently

```math
\boxed{
\frac{S_s(w)}M\le \Pr_w(J_{r,s})
\le \binom rs S_s(w).
}
\tag{R31.1}
```

This audits and slightly refines the direct maximum-degree sandwich

```math
h_w^s\le\Pr_w(J_{r,s})\le M\binom rs h_w^s,
```

because `S_s/M<=h_w^s<=S_s`.

At the project scale `log M<=rL`, take `s=ceil(alpha r)` for fixed
`alpha>0`.  Then both `log M` and `log binom(r,s)` are `O(rL)=O(sL)`.
It follows that, uniformly in `w`,

```math
\Pr_w(J_{r,s})\ge e^{-O(rL)}
\quad\Longleftrightarrow_{\rm exponent}\quad
S_s(w)\ge e^{-O(rL)}.
\tag{R31.2}
```

Here “equivalent at exponent” means that the constants in the `O(rL)`
may change by the candidate-count and binomial terms.  It is only a
diagnostic reduction, not a proof of the needed lower bound.

There is a useful power-saving tradeoff beyond fixed `alpha`.  Suppose the
construction starts with `L_0=n^(3/4-c_0)` and its corresponding batch size
`r`, and put `s=ceil(r/T)` for `1<=T<=r`.  The maximum-degree sandwich,
uniformly in `w`, gives

```math
\delta_*
\ge [M\binom rs]^{-1/s}
\inf_w\Pr_w(J_{r,s})^{1/s}.
\tag{R31.2a}
```

If

```math
\inf_w\Pr_w(J_{r,s})\ge e^{-K rL_0},
```

then `r/s<=T`, `log M<=rL_0`, and
`log binom(r,s)<=s log(eT)` imply the exact estimate

```math
\boxed{
\log\tau_*
\le (K+1)T L_0+\log(eT).
}
\tag{R31.2b}
```

Consequently `T<=n^eta` with `eta<c_0` still gives a power-saving cover:
with `c'=c_0-eta>0`, the right side is
`O(n^(3/4-c'))`.  The original row and tolerance bounds at exponent `c_0`
are stronger than those required at `c'`, so the usual fixed-cut argument
still proves convergence.  Thus it is enough to compress a likely batch
into `n^eta`, `eta<c_0`, hit cosets: pigeonhole gives `J_(r,ceil(r/T))`.
Pointwise planting only gives the deterministic choice `T=r`, for which
`T L_0` is of order `n log n`, and remains too expensive.

The partial-event power sum remains the exact scalar diagnostic throughout
this tradeoff.  Since `log binom(r,s)=O(r)=o(rL_0)`, (R31.1) says

```math
\Pr(J_{r,s})\ge e^{-O(rL_0)}
\quad\Longleftrightarrow_{\rm exponent}\quad
S_s(w)\ge e^{-O(rL_0)}.
\tag{R31.2c}
```

The loss `T` appears only when the `s`th root is taken to recover
`delta_*`.

More generally, for any real `1<=t<=r`, power means and the fixed first
`t` batch coordinates give

```math
\boxed{
\left(\frac{S_t(w)}M\right)^{r/t}
\le h_w^r\le\Pr_w(J_{r,r})\le S_r(w).
}
\tag{R31.3}
```

Thus a sufficient `t`-moment hypothesis is

```math
\log\frac{M}{S_t(w)}=O(tL).
\tag{R31.4}
```

Among this hierarchy, the `t=r` diagonal moment is the weakest scalar
incidence datum at the desired exponent: `Pr(J_(r,r))<=S_r` makes its
target-scale lower bound necessary, while (R31.1) makes it sufficient up
to the already affordable factor `M`.  Cross-codegrees are not the first
missing datum.

There is an exact entropy reformulation.  If `S_1>0` and
`eta_a=p_a/S_1`, then, with Renyi entropy `H_s(eta)`,

```math
\log S_s=s\log S_1-(s-1)H_s(\eta).
\tag{R31.5}
```

The missing statistic is therefore concentration of incidence edge mass
on candidate endpoints, after crediting the total incidence mass `S_1`.
Low Renyi entropy by itself is sufficient but not necessary when `S_1` is
large; the exact quantity is the left side of (R31.5).

## 2. What higher codegrees and Holder can do

For `J_(r,r)`, put

```math
Z=\sum_a\prod_{i=1}^r\mathbf1_{\{X_i\in H_a\}},
\qquad
p_{a_1\ldots a_q}=w(H_{a_1}\cap\cdots\cap H_{a_q}).
```

Then

```math
\mathbb EZ=S_r,
\qquad
T_{q,r}:=\mathbb EZ^q
=\sum_{a_1,\ldots,a_q}p_{a_1\ldots a_q}^{\,r}.
```

Holder applied to `Z 1_(Z>0)` gives, for every real `q>1`,

```math
\boxed{
\Pr(J_{r,r})\ge
\frac{S_r^{q/(q-1)}}{(\mathbb EZ^q)^{1/(q-1)}}.
}
\tag{R31.6}
```

For integer `q`, the denominator is the displayed `q`-fold codegree sum.
At `q=2` this is Paley--Zygmund,

```math
\Pr(J_{r,r})\ge
\frac{S_r^2}{\sum_{a,b}w(H_a\cap H_b)^r}.
\tag{R31.7}
```

The local de Caen form is also exact and sometimes sharper:

```math
\Pr(J_{r,r})\ge
\sum_a\frac{p_a^{2r}}
{\sum_b w(H_a\cap H_b)^r}.
\tag{R31.8}
```

Terms with `p_a=0` are omitted (equivalently assigned value zero) in
(R31.8).

These formulas identify valid possible inputs, but at the present
candidate entropy they are not needed after the diagonal moment is known.
Indeed `Z<=M 1_(Z>0)` implies

```math
\mathbb EZ^q\le M^{q-1}S_r,
```

and (R31.6) already recovers `Pr(J_(r,r))>=S_r/M`.  Conversely
`Pr(J_(r,r))<=S_r`, so no favorable cross-codegree estimate can compensate
for `S_r=e^{-omega(rL)}`.

Ordinary Finner inequalities do not introduce a hidden gain here: all
candidate events `H_a^r` depend on the same full set of `r` independent
coordinates, so their dependency hyperedges coincide and a fractional
matching has total weight at most one.  Any improvement needs a new
incidence property forced by exact minimization, not merely the product
cylinder form.

## 3. Pointwise degree and KKT audit

Let

```math
D(x)=|\{a:x\in H_a\}|,\qquad D_{\min}=\min_xD(x).
```

Since `S_1=E_w D(X)>=D_min`, power means give

```math
S_s\ge M^{1-s}D_{\min}^s,
\qquad
\Pr(J_{r,s})\ge(D_{\min}/M)^s.
\tag{R31.9}
```

Hence the useful normalized-degree input is

```math
D_{\min}/M\ge e^{-O(L)}.
\tag{R31.10}
```

Pointwise planting gives only `D_min>=1`, hence in the worst case only
`M^{-s}=e^{-O(srL)}` when `s=Theta(r)`, a factor `r` too costly in the
exponent.  A huge raw degree is likewise irrelevant unless it is compared
with `M`.

Now let

```math
\delta_*=\min_w\max_ap_a
=\max_{\mu\in\Delta(A)}\min_x
g_\mu(x),\qquad
g_\mu(x)=\sum_a\mu_a\mathbf1_{H_a}(x).
```

At a saddle point `(w,mu)`, complementary slackness says

```math
p_a=\delta_*\quad(\mu_a>0),\qquad
g_\mu(x)=\delta_*\quad(w_x>0).
```

Thus for every `s>=1`, exactly

```math
\sum_a\mu_ap_a^s=\delta_*^s.
\tag{R31.11}
```

The weighted common-candidate indicator for the first `s` samples is at
most `1_(J_(r,s))`, so (R31.11) recovers
`Pr(J_(r,s))>=delta_*^s`.  This is the lower half of the collision
sandwich, not a numerical estimate on `delta_*`.  Finite-game optimality
alone therefore adds no missing concentration.

If a batch is coverable by at most `T` hit candidates, one of them hits at
least `ceil(r/T)` batch members.  Therefore a probability `e^{-O(rL_0)}`
of such a batch cover is a sufficient (stronger) statistic for the event in
(R31.2a).  The growing-group tradeoff shows that even
`T<=n^eta`, `eta<c_0`, suffices.  Pointwise planting gives only an
`r`-candidate cover and does not imply this compression.

## 4. A symmetric extremal wall, including positive-fraction collisions

This construction is abstract; its candidates need not be row-good cosets
of an exact signing minimizer.  It cleanly separates incidence/KKT facts
from the missing minimizer correlation.

Let `|X|=N` and take as candidates all `t`-element subsets `B subset X`,
with `H_B=B`, where `t` is fixed (for example `t=3`).  Then

```math
M=\binom Nt,
\qquad
D(x)=\binom{N-1}{t-1},
\qquad
D(x)/M=t/N.
```

Every tuple of at most `t` selectors shares a candidate, including tuples
with repetitions.  The system is transitive and incidence-regular.

For uniform `w`, every candidate has `p_B=t/N`.  For arbitrary `w`, the
sum of the `t` largest atom masses is at least `t/N`, and that sum is the
mass of one candidate.  Hence

```math
\delta_*=t/N.
\tag{R31.12}
```

The uniform candidate law includes every selector with probability `t/N`,
so it is a primal optimum; all dual and primal constraints are active and
the complete KKT conditions hold.

Nevertheless, `J_(r,r)` occurs exactly when the sample has at most `t`
distinct values, and

```math
\Pr(J_{r,r})
=N^{-r}\sum_{u=1}^t(N)_u\,{r\brace u}
\le\binom Nt(t/N)^r.
\tag{R31.13}
```

More strongly, for every `s<=r`, a union bound over the candidate and the
`s` captured batch coordinates gives

```math
\boxed{
\Pr(J_{r,s})
\le\binom Nt\binom rs(t/N)^s.
}
\tag{R31.14}
```

Choose `log N=Theta(n)` (as for a fixed-density selector slice), fixed
`t`, and first take `s=ceil(alpha r)`.  Since `r` grows and `L=o(n)`,

```math
\log\Pr(J_{r,s})
\le-s\log N+O(n+r)
=-\Theta(rn)=-\omega(rL).
```

Equivalently, writing it without the deliberately expanded middle term,
the leading contribution is `-s log N=-Theta(rn)`, while
`log M=O(n)` and `log binom(r,s)=O(r)` are lower order.  Thus pointwise
planting, enormous regular raw degree, perfect common cover through every
fixed order `t`, and exact minimax/KKT optimality all coexist with failure
of every fixed-positive-fraction collision bound.

The same wall audits the growing-group tradeoff.  With
`s=ceil(r/T)`, (R31.14) has leading logarithm `-Theta(rn/T)`.  It is
`-omega(rL_0)` whenever

```math
T=o(n/L_0)=o(n^{1/4+c_0}).
```

In particular it defeats the entire sufficient range
`T<=n^eta`, `eta<c_0`.  Again this is an abstract warning that even a
highly symmetric, fully optimal incidence game need not have the required
batch compression; it is not a signing counterexample.

This wall is sharp for the diagonal diagnostic: here

```math
S_s=\binom Nt(t/N)^s=e^{-\Theta(sn)}
```

at mesoscopic `s`, which is already below the needed scale.  It does not
falsify the signing project because it imports no exact-minimizer
structure.

The exact `A_9` data fit this diagnosis but do not resolve it.  Their
`J_(2,2)=J_(3,3)=1` statements force `S_2,S_3>=1` by (R31.1), for every
selector law.  They give no growing-order moment, and the scoped
`t`-subset wall proves that no abstract interpolation from finitely many
fixed moments to `s=Theta(r)` is possible.

## 5. Resulting open target

The bare route can now be stated without higher-overlap ambiguity.  At the
original exponent `c_0`, it is enough to prove for every adversarial
selector law `w`, with either fixed `alpha>0` or more generally
`s=ceil(r/T)` and `T<=n^eta`, `eta<c_0`, that

```math
\boxed{
\sum_{a\in\mathscr A_C}w(H_a)^s
\ge e^{-O(rL)}.
}
\tag{R31.15}
```

For fixed `alpha` this is necessary and sufficient at the same cover
exponent; for growing `T` it yields the degraded but still positive power
`c'=c_0-eta`.  This is weaker than prescribing coherent completions and
stronger than pointwise planting.  The next proof input must show that
exact-minimizer hit incidence cannot spread its edge mass across candidates
as in the `t`-subset wall.  Neither candidate counting, fixed-order
codegrees, ordinary Holder/Finner, nor saddle-point optimality establishes
(R31.15).
