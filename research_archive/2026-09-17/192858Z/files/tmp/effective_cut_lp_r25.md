# Wave 25 root audit: all-cut effective-loss reduction

Wave 24 specialized the global-overlap construction to exact parent grounds and discarded the favorable deficit between an output cut and the parent optimum. Neither restriction is required by the information/mgf argument. Retaining the deficit gives a strictly weaker all-cut criterion, and the row-square integrality extraction found independently in Wave 25 reduces its mixture LP to one low-cost cut up to a polynomial factor.

The analytic identities below are **Verified** from (10.759)--(10.760). The finite `A_6,A_8,A_9` tables are exact exhaustive calculations checked by `effective_cut_lp_r25.py`. The resulting all-cut coverage lemma is an **Open target**; no convergence proof is claimed.

## 1. Retain the parent deficit

Let `A` be an exact order-`n` minimizer, `q_n=Q(A)`, `S` an `m`-selector, and `d` any oriented full cut. Put

```math
\Delta_A(d)=q_n-\langle A,d\rangle\ge0,
\qquad
\ell(S,d)=Q(A[S])-c_A(S,d),
```

Let

```math
B_{n,m}=\left[\left(\frac mn\right)^{3/2}-p_2\right]q_n>0
```

be the deterministic coefficient slack from (10.646), and define the centered effective loss

```math
\boxed{
\widehat\ell(S,d)
=\ell(S,d)-p_2\Delta_A(d)-B_{n,m}
=Q(A[S])-p_2q_n-
\bigl[c_A(S,d)-p_2\langle A,d\rangle\bigr]-B_{n,m}.
}
\tag{E.1}
```

For an arbitrary selector--cut channel `P`, the proof of (10.761) before dropping `-p_2 E Delta_A(D)` gives

```math
\boxed{
V_{\rm ad}(A,m)-\left(\frac mn\right)^{3/2}q_n
\le
\mathbb E_P\widehat\ell(S,D)
+\epsilon_{n,m}q_n
+\frac{J}{\lambda}
+\lambda\left[p^2\mathbb E_PR_2(D)+\frac{n^2}{2}\right],
}
\tag{E.2}
```

on the same halved spectral domain as (10.761), where

```math
J=I_{\rm Sh}(S;D)+D_{\rm KL}(P_S\Vert U_m)+\chi_{n,m}.
```

Indeed, write `Q(A[S])=ell(S,D)+c_A(S,D)`, retain
`p_2 E<A,D>-p_2q_n=-p_2 E Delta_A(D)`, subtract `B_{n,m}`, and then apply the exact output-reference optimization (10.760) and the arbitrary-cut slice mgf (10.759). The latter depends on signed row fields only through `R_2`; nonnegativity of rows and the parent-ground property are never used. The bound `|<A,d>|<=q_n` holds for every oriented cut by definition of `Q(A)`.

For every fixed cut, the uniform-selector mean is cut-independent:

```math
\boxed{
\mathbb E_{S\sim U_m}\widehat\ell(S,d)
=\mathbb E_{S\sim U_m}Q(A[S])
-\left(\frac mn\right)^{3/2}q_n.
}
\tag{E.3}
```

Thus this is a tail/compression reformulation, not a new scalar averaging inequality.

## 2. Global conditioning over all cuts

Let `nu` be any law on the entire cut set `D_n`, and put

```math
G=\{(S,d):\widehat\ell(S,d)\le t\},
\qquad Z=(U_m\otimes\nu)(G).
```

Conditioning the product law on `G` gives `E hat ell<=t` and

```math
I(S;D)+D_{\rm KL}(P_S\Vert U_m)\le-\log Z,
\qquad
\mathbb E_PR_2(D)=\mathbb E[R_2(D)\mid G].
```

Consequently, for fixed `0<c<1/4`, it suffices to prove

```math
\boxed{
t=O(n^{3/2-c}),\qquad
-\log Z=O(n^{3/4-c}),\qquad
\mathbb E[R_2(D)\mid G]=O(n^{9/4-c}).
}
\tag{E.4}
```

Taking `lambda` proportional to `n^{-3/4}` in (E.2) then yields the required power-saving principal-restriction edge. Exact parent grounds are only the special case `Delta_A(d)=0`. Even there the event is `ell(S,d)<=B_{n,m}+t`; (10.764) unnecessarily discarded the available leading slack by requiring `ell<=t`.

For

```math
u_d=U_m\{S:\widehat\ell(S,d)\le t\},
\qquad c_d=R_2(d),
```

the same one-constraint LP and dual as (10.765)--(10.766) hold, now over all cuts:

```math
Z_*^{\rm all}(C)
=\max_{\nu\in\Delta(\mathcal D_n\cup\{0\})}
\left\{\sum_d\nu_du_d:
\sum_d\nu_du_d(c_d-C)\le0\right\}
=\min_{\theta\ge0}\max_d
\bigl[u_d(1+\theta(C-c_d))\bigr]_+.
\tag{E.5}
```

## 3. Polynomial extraction makes one cut enough

Take the budget `C` to be a nonnegative integer; replacing a real asymptotic budget by its ceiling loses nothing because every `c_d` is an integer. Let `mu_d=nu_d u_d/Z` be the captured output law of any feasible positive-overlap solution. Since `E_mu c_d<=C` and every cost above `C` is at least `C+1`,

```math
\mu\{c_d\le C\}\ge\frac1{C+1}.
```

In general

```math
\mathbb E_\mu\frac1{u_D}
=\frac{\nu\{d:u_d>0\}}Z\le\frac1Z;
```

equality holds after deleting dummy and zero-coverage mass, as one may for an optimizer. This inequality already forces some captured cut with

```math
\boxed{
c_d\le C,
\qquad
u_d\ge\frac{Z}{C+1}.
}
\tag{E.6}
```

Otherwise the contribution of the low-cost set alone would exceed `1/Z`, contradicting the displayed upper bound. Conversely one such cut is a feasible singleton prior. Since `log(C+1)=O(log n)`, (E.4) is therefore equivalent at its stretched-exponential scale to the following simpler **Open target**:

> For one active fixed-ratio window and every requested target order, choose a target-specific exact minimizer `A` and find one full cut `d` with `R_2(d)=O(n^{9/4-c})` and
> `U_m{hat ell(S,d)<=O(n^{3/2-c})}>=exp{-O(n^{3/4-c})}`.

Mixtures of two cuts can improve exact finite LP values, but cannot repair a stretched-exponential failure of every under-budget singleton.

## 4. Exact one-deletion audit

At `m=n-1`, exhaustive enumeration at the stronger threshold `ell(S,d)-p_2 Delta_A(d)<=0` gives these nondominated all-cut `(u_d,R_2(d))` pairs, and hence lower bounds the centered-effective-loss coverage at `t=0`:

| minimizer | nondominated all-cut pairs |
|---|---|
| `A_6` | `(5/6,30)` |
| `A_8` | `(3/8,32)`, `(1/2,64)` |
| `A_9` | `(1/9,24)`, `(2/9,32)`, `(1/3,48)`, `(4/9,112)` |

The exact-ground front begins at row-square `30,64,80`, respectively. Thus the all-cut formulation agrees for `A_6` but is strictly better for `A_8` and `A_9`; for example, `A_9` already has positive zero-threshold coverage at row-square `24`, whereas every positive-coverage parent ground has row-square at least `80`. The uncentered means of `ell-p_2 Delta` are `4/3`, `3`, and `16/3`. The centered means in (E.3) are, respectively,

```math
8-10(5/6)^{3/2},
\qquad 18-20(7/8)^{3/2},
\qquad 24-24(8/9)^{3/2}.
```

These examples validate that the deficit term is operational rather than cosmetic. They do not supply an asymptotic coverage theorem.

## 5. Consequence for route selection

The strongest current implementation of optimized restriction should use the singleton all-cut effective-loss target (E.6), not the parent-ground two-column target (10.764). Exact-ground exchange geometry may still be useful for constructing a cut, but it is no longer a necessary part of the sufficient lemma. Any proof must create a fixed-cut upper tail for the centered selector payoff while controlling that cut's row-square; (E.3) shows why scalar means alone remain circular.
