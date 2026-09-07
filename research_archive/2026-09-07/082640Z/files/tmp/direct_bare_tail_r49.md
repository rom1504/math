# Wave 49B: an annealed bare-incidence reduction

## Status

- **Verified:** a lower bound on the *total* favorable incidence of arbitrary
  cuts automatically supplies both the project row bound and one fixed column
  satisfying (10.795).  The row need not be imposed in the incidence event.
- **Verified:** the total incidence has exact principal-norm, quadratic-chaos,
  and completion-deficit profile formulas.  These retain the two allowances in
  `widehat ell` and do not impose `L_S(d)>=q_n` or complement incidence.
- **Open:** the required annealed lower bound.  No generic concentration
  argument proves its saved exponent across a genuine normalized restriction
  gap.
- **Verified scoped wall:** robust coordinate-cylinder certificates at the
  project entropy already force the desired recurrence by conditional
  averaging.  They are therefore not an independent mechanism.
- **Exact finite enumeration:** at threshold zero on `A6,A8,A9`, the annealed
  favorable mass is between `0.2367` and `0.5`.  In the high-ratio cases,
  `86%--95%` of favorable incidences use the retained full deficit and lie
  outside the older local-deficit incidence.  This is mechanism evidence only.

## 1. Annealed favorable incidence makes the row automatic

Let `Pi_n` be uniform on the `2^n` oriented projective cuts
`d=(sigma,x)`.  For a fixed exact minimizer `A` with `Q(A)=q_n`, define

```math
F_t=\{(S,d):\widehat\ell(S,d)\le t\},\qquad
Z_t=(U_m\otimes\Pi_n)(F_t),\qquad L_t=-\log Z_t.
```

Assume `Z_t>0`, and condition the product law on `F_t`; call the resulting
law `P`.  Since its density is `1/Z_t` on `F_t`, the KL chain rule gives the
exact identity

```math
L_t=D(P\Vert U_m\otimes\Pi_n)
=D(P_D\Vert\Pi_n)+D(P_S\Vert U_m)+I_P(S;D).
\tag{R49.1}
```

In particular, `D(P_D||Pi_n)<=L_t`.  The row cost is independent of the
orientation and is

```math
R_2(d)=x^{\mathsf T}A^2x
=n(n-1)+x^{\mathsf T}Cx,\qquad
C=A^2-(n-1)I.
```

The standard Rademacher Hanson--Wright mgf and entropy duality imply, for
every cut law `mu` with `D(mu||Pi_n)<=L`,

```math
\mathbb E_\mu R_2
\le n(n-1)+C_0\bigl(\lVert C\rVert_F\sqrt L
+\lVert C\rVert_{\rm op}L\bigr).
\tag{R49.2}
```

For an exact minimizer, the established `||A||_op^2<=2q_n` gives

```math
\lVert C\rVert_F^2
=\operatorname{tr}A^4-n(n-1)^2
\le2q_n n(n-1)=O(n^{7/2}),
\qquad
\lVert C\rVert_{\rm op}\le2q_n+n=O(n^{3/2}).
\tag{R49.3}
```

Consequently

```math
\overline R:=\mathbb E_P R_2(D)
\le n(n-1)+O\bigl(n^{7/4}\sqrt{L_t}+n^{3/2}L_t\bigr).
\tag{R49.4}
```

For a fixed output put

```math
u_t(d)=U_m\{S:\widehat\ell(S,d)\le t\}.
```

The posterior under the conditioned product law is exactly `U_m` conditioned
on this set.  Hence, as useful bookkeeping,

```math
\mathbb E_{P_D}\log\frac1{u_t(D)}
=D(P_S\Vert U_m)+I_P(S;D)
=L_t-D(P_D\Vert\Pi_n)\le L_t.
\tag{R49.5}
```

There is a sharper extraction than averaging the two costs.  The conditioned
output law is exactly

```math
P_D(d)=\frac{\Pi_n(d)u_t(d)}{Z_t}.
```

Let `G={d:R_2(d)<=2\overline R}`.  Markov gives `P_D(G)>=1/2`, and therefore

```math
\sum_{d\in G}\Pi_n(d)u_t(d)=Z_tP_D(G)\ge Z_t/2.
```

Since `Pi_n(G)<=1`, one cut `d in G` satisfies

```math
\boxed{u_t(d)\ge Z_t/2=\tfrac12e^{-L_t},\qquad
R_2(d)\le2\overline R.}
\tag{R49.6}
```

This proves the following exact sufficient interface:

```math
\boxed{
t=O(n^{3/2-c}),\qquad
-\log Z_t=O(n^{3/4-c})
\quad\Longrightarrow\quad (10.795),
}
\tag{R49.7}
```

uniformly on a compact ratio window, for every fixed `0<c<1/4`.  Indeed,
at `L_t=O(n^(3/4-c))`, the two errors in (R49.4) are respectively
`O(n^(17/8-c/2))` and `O(n^(9/4-c))`; the first is no larger than the second
precisely for `c<=1/4`.  Thus (R49.6) has the project row and mass exponents.

This is not yet a proof of convergence: (R49.7)'s annealed mass inequality is
open.  Its concrete gain is that arbitrary-cut construction no longer has to
carry a row constraint.  A scalar two-variable counting estimate suffices,
and entropy of the conditioned output enforces the row afterward.

## 2. Exact norm-profile and completion-profile targets

For a physical projective spin `x`, put

```math
W_S=P_SAP_S-p_2A,\qquad
X_S(x)=x^{\mathsf T}W_Sx,\qquad
Y_S=Q(A[S])-p^{3/2}q_n.
```

For the two orientations of the same physical cut,

```math
\boxed{\widehat\ell(S,(\sigma,x))=Y_S-\sigma X_S(x).}
\tag{R49.8}
```

Writing `a_S=Y_S-t`, the conditional favorable density is therefore exactly

```math
u_t(S):=\Pi_n\{d:\widehat\ell(S,d)\le t\}
=\begin{cases}
\frac12\Pr_x\{|X_S(x)|\ge a_S\},&a_S>0,\\[2mm]
\frac12\bigl(1+\Pr_x\{|X_S(x)|\le-a_S\}\bigr),&a_S\le0.
\end{cases}
\tag{R49.9}
```

Thus `Z_t=E_S u_t(S)`.  Since `Q(A[S])` takes only `O(n^2)` integer values,
a concrete profile-conditioned route is to find one value `r` for which,
when `r-p^(3/2)q_n-t>0`,

```math
\boxed{
U_m\{Q(A[S])=r\}\,
\Pr_{\substack{S\mid Q(A[S])=r\\x\sim\mathrm{Unif}}}
\left\{|x^{\mathsf T}(P_SAP_S-p_2A)x|
\ge r-p^{3/2}q_n-t\right\}
\ge e^{-O(n^{3/4-c})}.}
\tag{R49.10}
```

The factor `1/2` and the polynomial number of norm profiles cost only
`O(log n)`.  Equation (R49.10) is a genuine quadratic-chaos/counting target,
not a fixed-cut slice tail and not complement incidence.

There is a second exact form which displays the retained deficit more
directly.  Let `a=(sigma,y)` be a uniform oriented projective local state on
`S`, let

```math
\delta_S(a)=Q(A[S])-\sigma y^{\mathsf T}A[S]y,\qquad
B=(p^{3/2}-p_2)q_n,
```

and let `E_{S,a}(w)` be the full oriented energy after a uniform outside
completion `w`.  Define its conditional CDF
`F_{S,a}(r)=Pr_w{E_{S,a}(w)<=r}`.  Direct rearrangement of `widehat ell<=t`
gives

```math
\boxed{
u_t(S)=\mathbb E_a
F_{S,a}\!\left(q_n+\frac{B+t-\delta_S(a)}{p_2}\right).}
\tag{R49.11}
```

If `delta_S(a)<=B+t`, the argument of the CDF is at least `q_n`, so this
term is one.  This recovers the older approximate-child incidence.  But
(R49.11) also retains every state with larger local deficit whose outside
completion has enough full-energy deficit.  This second population is the
new part that the complement and local-ground formulations discard.

The exact open lemma exposed by this wave is (R49.7), equivalently an
`e^{-O(n^(3/4-c))}` lower bound on the average of (R49.9) or (R49.11), for a
target-specific exact minimizer at every active pair.  A fixed proposed set
of constants is falsified by an unbounded active family on which
`-log Z_t` exceeds that budget.  Such a failure only falsifies the annealed
implementation, not the bare tail itself.

## 3. Why robust coordinate profiles are not the mechanism

One tempting way to prove (R49.10) is to find, for each useful selector, a
coordinate subcube of favorable cuts.  It has an exact wall.

Suppose a coordinate cylinder is entirely contained in `F_t(S)` and has
total chart codimension `K`: it fixes the orientation bit and `K-1`
independent tail-spin bits.  Equivalently, after including the always-fixed
projective gauge anchor, its physical spin values are fixed on a set `J` of
`K` vertices.  (If `Y_S-t>0`, the orientation cannot remain free.)  Its
product mass is exactly `2^{-K}`.  Averaging (R49.8) over the free physical
bits kills every coefficient except those with both endpoints in `J`.
Since every off-diagonal coefficient of `W_S` is either `(1-p_2)a_ij` or
`-p_2a_ij`, it has absolute value at most one.  Therefore

```math
\boxed{Y_S\le t+K(K-1).}
\tag{R49.12}
```

At `K=O(n^(3/4-c))` and `t=O(n^(3/2-c))`, (R49.12), together with
`Q(A[S])>=q_m`, already gives the desired restriction recurrence.  Thus a
project-entropy robust cylinder cannot explain a rare upper tail across a
genuine normalized gap.  The live version must estimate the curved
quadratic-tail mass in (R49.9)/(R49.11), not replace it by an all-completions
cell.

Equivalently, if `r` denotes only the additional fixed tail spins, then
`K=r+1`, the mass is `2^{-(r+1)}`, and the right side of (R49.12) is
`t+r(r+1)`.  Conditional uniform completion also has

```math
\mathbb E R_2
=\lVert A[:,J]u_J\rVert_2^2+(n-|J|)(n-1)
\le2q_n|J|+n(n-1),
```

so such a cylinder does give (10.795) when `K` has the project exponent.
Equation (R49.12) shows why this correct theorem is circular rather than a
new route.

## 4. Finite audit

`tmp/direct_bare_tail_r49_check.py` enumerates all oriented projective cuts,
all selectors, and all coordinate cylinders for `A6,A8,A9`, at `t=0`.  The
comparisons with `p^(3/2)q_n` are made exactly by squaring nonnegative
rationals; the printed floating decision margin is at least `0.0717`.
It also verifies (R49.11), both KL data-processing inequalities, and the
cylinder assertions.  Selected outputs are:

| case | `Z_0` | strong local mass | retained-only share of `Z_0` | `-log Z_0` | best fixed-column coverage |
|:--|--:|--:|--:|--:|--:|
| `A6,m=3` | `0.265625` | `0.125000` | `52.94%` | `1.32567` | `0.35000` |
| `A8,m=4` | `0.367857` | `0.116071` | `68.45%` | `1.00006` | `0.57143` |
| `A8,m=6` | `0.400670` | `0.055804` | `86.07%` | `0.91462` | `0.64286` |
| `A8,m=7` | `0.390625` | `0.023438` | `94.00%` | `0.94001` | `0.62500` |
| `A9,m=5` | `0.348230` | `0.147569` | `57.62%` | `1.05489` | `0.67460` |
| `A9,m=7` | `0.236708` | `0.018229` | `92.30%` | `1.44093` | `0.44444` |
| `A9,m=8` | `0.252604` | `0.012153` | `95.19%` | `1.37593` | `0.44444` |

Across all eleven audited `(A,m)` pairs, `Z_0` lies in
`[0.2367079,0.5]`.  The annealed output mean row lies between `30` and
`74.29`; it is at the ordinary `Theta(n^2)` scale.  Minimum robust-cylinder
codimensions range from four to six, comparable to the unsaved small-order
scale and carrying no asymptotic conclusion.

The sharp finite message is the retained-only column: at high ratios, almost
all favorable pairs would be erased by the stronger local-deficit incidence.
Thus (R49.7)/(R49.11) is not a renamed version of that retired route.  The
asymptotic saved-exponent lower bound remains completely open and should be
tested next either by a minimizer-specific lower-tail theorem for the
completion CDFs in (R49.11), or by searching for an exact-minimizer family
whose `-log Z_t` has the unsaved `Theta(n^(3/4))` scale.
