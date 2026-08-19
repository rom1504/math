# Labelled Fourier quotients need at least square-root dimension

Status: **proved exact labelled-Fourier reduction and scalable no-go for
small Fourier quotients of the actual switch-convolution certificate**.  The
Fourier coefficients of both actual kernels are derived explicitly.  An
`ell`-dimensional character quotient gives a rigorous joint
child--bridge cancellation certificate, but its negative log gain is at
most

```math
2t\max\{m,n\}\,\ell.                               \tag{0.1}
```

On comparable splits, this is `O_beta(ell sqrt(N))`.  Hence a positive
linear pre-cancellation defect cannot be removed by any
`ell=o(sqrt(N))` labelled Fourier quotient, regardless of the children and
regardless of the chosen sign bridge.  Because the actual child kernel has
no odd--odd Fourier modes, the `(N-2)`-dimensional even--even quotient is
already exactly the original fixed-child switch/bridge optimization.  The
theorem leaves a real intermediate window, but rules out a bounded or
slowly growing list of negative Fourier modes.

## 1. Characters of the rank-one switching group

Let

```math
\mathcal G_{m,n}=\{xy^{\mathsf T}:x\in\{\pm1\}^m,
 y\in\{\pm1\}^n\}.
```

It is an elementary abelian group of order `2^(N-1)`, `N=m+n`.  Its
characters are indexed by pairs `(S,T)` with

```math
S\subseteq[m],\quad T\subseteq[n],\quad
|S|+|T|\equiv0\pmod2,
```

through

```math
\chi_{S,T}(xy^{\mathsf T})=x_Sy_T,
\qquad x_S=\prod_{i\in S}x_i.                      \tag{1.1}
```

The parity condition makes (1.1) invariant under the simultaneous gauge
`(x,y)->(-x,-y)`, and the resulting `2^(N-1)` characters are distinct.
Use normalized Fourier coefficients

```math
\widehat f(S,T)=\mathbb E_{Q\in\mathcal G}f(Q)
 \chi_{S,T}(Q).                                    \tag{1.2}
```

For child signings `A,D`, orientation `epsilon`, sign bridge `B`, and
`t=beta/sqrt(N)`, the exact switch-convolution kernels are

```math
w_\epsilon(xy^{\mathsf T})
=\cosh t\{H_A(x)+\epsilon H_D(y)\},
\qquad
k_B(xy^{\mathsf T})=e^{t x^{\mathsf T}By}.         \tag{1.3}
```

## 2. Exact Fourier coefficients of the actual kernels

For an even subset `S`, define

```math
C_A(S)=\mathbb E_x x_S\cosh(tH_A(x)),
\qquad
S_A(S)=\mathbb E_x x_S\sinh(tH_A(x)),              \tag{2.1}
```

and analogously for `D`.  Since `H_A(-x)=H_A(x)`, both quantities vanish
when `|S|` is odd.  Expanding the addition formula for `cosh` proves

```math
\boxed{
\widehat w_\epsilon(S,T)=
\begin{cases}
C_A(S)C_D(T)+\epsilon S_A(S)S_D(T),
   & |S|,|T|\text{ even},\\
0, & |S|,|T|\text{ odd}.
\end{cases}}                                       \tag{2.2}
```

Thus the internal kernel deletes the entire odd--odd character sector; it
does not make the surviving even--even coefficients scalar or unsigned.

Put `theta=tanh t`.  The edgewise high-temperature expansion gives, for
every allowed `(S,T)`,

```math
\boxed{
\widehat k_B(S,T)
=(\cosh t)^{mn}
\sum_{\substack{F\subseteq K_{m,n}:\\
 \partial_LF=S,\ \partial_RF=T}}
 \theta^{|F|}\prod_{ia\in F}B_{ia}.}              \tag{2.3}
```

Here the two boundaries are the odd-degree vertex sets of `F`.  Formula
(2.3) is a signed boundary enumerator, not a function merely of `|S|,|T|`
or of the singular values of `B`.

Normalize

```math
a={w_\epsilon\over\widehat w_\epsilon(\varnothing,\varnothing)},
\qquad
b={k_B\over\widehat k_B(\varnothing,\varnothing)}. \tag{2.4}
```

Both have Haar mean one.  Fourier inversion of the exact convolution is

```math
\boxed{
(a*b)(g)=\sum_{S,T}
 \widehat a(S,T)\widehat b(S,T)\chi_{S,T}(g).}      \tag{2.5}
```

Only the even--even terms from (2.2) survive.  A bridge switch multiplies
each coefficient in (2.3) by the corresponding character sign, so (2.5)
is precisely the labelled joint geometry available to switching.

## 3. Every Fourier subspace is an exact coarse certificate

Let `V` be an `ell`-dimensional subspace of the character group.  Choose a
basis and let

```math
\pi_V:\mathcal G\longrightarrow\mathbb F_2^\ell
```

record its `ell` character values in additive `0/1` notation.  For a
density `f` on `mathcal G`, define its quotient density

```math
f_V(z)=\mathbb E[f(Q)\mid\pi_V(Q)=z].               \tag{3.1}
```

It has mean one on the uniform quotient.  Conditional expectation commutes
with group convolution, so

```math
\boxed{(a*b)_V=a_V*b_V.}                            \tag{3.2}
```

Equivalently, (3.2) is the Fourier polynomial (2.5) truncated to the
subspace `V`; unlike an arbitrary truncation, it remains nonnegative.
Since every quotient value is an average over one fibre,

```math
\min_{g\in\mathcal G}(a*b)(g)
\le\min_{z\in\mathbb F_2^\ell}(a_V*b_V)(z).        \tag{3.3}
```

Let

```math
\mathcal C_{\epsilon,B}
=\log\mathbb Ew_\epsilon+\log\mathbb Ek_B
-P_m(\beta)-P_n(\beta).                            \tag{3.4}
```

For exact own-scale minimizing children, (3.3) gives the required direct
defect arrow

```math
\boxed{
E_{m,n}(\beta)
\le \mathcal C_{\epsilon,B}
 +\log\min_z(a_V*b_V)(z).}                         \tag{3.5}
```

In particular, with no omitted normalization or comparison step,

```math
\boxed{
\mathcal C_{\epsilon,B}+\log\min_z(a_V*b_V)(z)
\le C_\beta N^{1-\delta}
\quad\Longrightarrow\quad
E_{m,n}(\beta)\le C_\beta N^{1-\delta}.}           \tag{3.5a}
```

For `ell=1`, if `V` is generated by a nontrivial character `chi`, (3.5)
reduces to the explicit negative-mode certificate

```math
\boxed{
E_{m,n}(\beta)
\le\mathcal C_{\epsilon,B}
 +\log\{1-|\widehat a(\chi)\widehat b(\chi)|\}.} \tag{3.6}
```

Indeed the two quotient values have mean one and difference twice the
displayed Fourier product.

At the other endpoint, let `W` be the even--even character subspace.  It has
dimension `N-2`, and (2.2) says that every Fourier coefficient of `a*b`
outside `W` vanishes.  Hence `a*b` is constant on cosets of the annihilator
of `W`; conditional expectation onto `W` loses no information, (3.3) is
equality, and (3.5) is exactly

```math
\min_gL_{\epsilon,g}(B)-P_m-P_n.                   \tag{3.7}
```

After minimizing `B`, (3.7) is the complete fixed-child bridge
optimization.  Taking all `N-1` characters is also exact but is redundant
by one dimension.  Thus the quotient dimension measures, without a change
of problem, how much of the exact bridge search has been reconstructed.

## 4. Gibbs smoothness bounds every small quotient

Parameterize `mathcal G` by the `N-1` independent vertex switches

```math
x_1,\ldots,x_m,y_2,\ldots,y_n,
```

fixing `y_1=1`.  Flipping one `x_i` changes `x^TBy` by at most `2n`, and
flipping one free `y_a` changes it by at most `2m`.  Consequently

```math
e^{-2tn}\le{k_B(Qe_i)\over k_B(Q)}\le e^{2tn},
\qquad
e^{-2tm}\le{k_B(Qe_a)\over k_B(Q)}\le e^{2tm}.    \tag{4.1}
```

Represent the quotient map `pi_V` by an `ell` by `(N-1)` binary matrix.
Choose `ell` pivot columns.  Their quotient increments form a basis, so any
two quotient fibres can be connected by at most `ell` corresponding vertex
flips.  Summing (4.1) along that path shows that their unnormalized bridge
masses have ratio at most

```math
\exp\{2t\max(m,n)\ell\}.                            \tag{4.2}
```

The same holds for their normalized quotient densities.  Since `b_V` has
uniform mean one, (4.2) implies

```math
\boxed{
\min_z b_V(z)\ge e^{-2t\max(m,n)\ell}.}            \tag{4.3}
```

Convolution with the nonnegative mean-one density `a_V` cannot go below
the minimum of `b_V`.  Therefore

```math
\boxed{
\log\min_z(a_V*b_V)(z)
\ge-2t\max(m,n)\ell.}                              \tag{4.4}
```

This is uniform over every child pair, every orientation, every sign bridge,
and every labelled `ell`-dimensional Fourier subspace.

A slightly sharper intrinsic version is available.  If the chosen pivot
columns contain `p` left and `q` right vertex switches (`p+q=ell`), then
the right side of (4.4) can be replaced by

```math
-2t(np+mq).                                         \tag{4.5}
```

Minimizing this cost over all pivot bases of the quotient gives its exact
elementary Lipschitz bound.

## 5. Scalable no-go and the surviving window

Suppose that for some bridge/orientation family the pre-cancellation term
has a positive linear floor

```math
\mathcal C_{\epsilon,B}\ge cN-o(N),\qquad c>0.     \tag{5.1}
```

This is exactly the regime in which a joint Fourier gain is needed.  The
numerical value of the quotient certificate on the right side of (3.5),
not merely our estimate of it, satisfies by (4.4)

```math
\boxed{
\mathcal C_{\epsilon,B}+\log\min(a_V*b_V)
\ge cN-o(N)-2\beta{\max(m,n)\over\sqrt N}\ell.}   \tag{5.2}
```

For comparable splits, every `ell=o(sqrt(N))` quotient therefore retains a
positive linear floor.  More quantitatively, if
`max(m,n)<=rho N`, cancellation of `cN+o(N)` requires

```math
\boxed{
\ell\ge {c\over2\beta\rho}\sqrt N-o(\sqrt N).}    \tag{5.3}
```

Thus one negative mode, any bounded collection of modes, and every quotient
with `2^o(sqrt(N))` response sectors are rigorously too small.  This is not
the old scalar-channel loss: (3.2)--(3.6) preserve child--bridge
cancellation exactly inside the selected labelled quotient.

The quantifier in this conclusion is conditional on (5.1).  In particular,
if the scalar orientation factor inside `mathcal C_(epsilon,B)` already
supplies a negative linear correction, (5.1) need not hold and this theorem
does not obstruct that separate escape.  What (5.2) says uniformly is that,
whenever a linear amount remains to be cancelled by labelled switch modes,
an `o(sqrt(N))` quotient cannot do it.

The theorem does **not** say that full landscape information is necessary,
and it does not rule out a useful growing quotient in the intermediate
window.  It leaves the interval

```math
\sqrt N\lesssim\ell<N-2                             \tag{5.4}
```

open.  At `ell=N-2` the even--even state is exactly the full switch
optimization by (3.7).  A successful labelled-Fourier architecture must
therefore either

1. find a structured quotient of at least square-root dimension whose
   `2^ell` sector table closes algebraically, or
2. use an inequality which accesses many labelled modes implicitly without
   materializing their quotient table.

## 6. Exact finite falsification of the one-mode route

[`audit_switch_fourier_actual_children.py`](../computations/audit_switch_fourier_actual_children.py)
enumerates exact child pressure minimizers and all bridge switching classes
through equal child order five.  It verifies (2.5) against direct parent
enumeration and compares (3.6) with the full quotient (3.7).

At `beta=4`, the best certificates over both orientations and every bridge
class are:

| child order | best one-character | best full switch convolution |
|---:|---:|---:|
| 3 | `0.3420077422` | `-3.2037699582` |
| 4 | `1.9498729156` | `-1.0903249513` |
| 5 | `5.6921143951` | `1.3628812578` |

At `beta=8` they are:

| child order | best one-character | best full switch convolution |
|---:|---:|---:|
| 3 | `0.9975037338` | `-8.8961345428` |
| 4 | `5.5346941698` | `-3.0151911658` |
| 5 | `15.1460190415` | `3.3623410788` |

The finite data are a falsifier, not an asymptotic theorem: even where the
full labelled convolution finds a negative certificate, its best single
mode remains positive.  Frozen JSON outputs accompany the script.  The
all-order content is Theorem (4.4), not these small values.

### Exhaustive subspace profile at child order four

[`audit_switch_coordinate_quotient_profile.py`](../computations/audit_switch_coordinate_quotient_profile.py)
goes beyond the one-mode audit.  For a fixed exact rare bridge it enumerates
**every** binary linear character subspace when the switch dimension is at
most seven.  Thus its dimension profiles optimize over all quotient
orientations, not merely over coordinate subsets.

For the best order-four bridge at `beta=4`, the best certificate at quotient
dimensions `0,...,7` is

```text
4.83849, 3.03183, 1.49940, 0.63887,
0.25106, -0.05651, -1.09032, -1.09032.
```

At `beta=8` it is

```text
14.07255, 9.13620, 5.53467, 1.47959,
0.97526, 0.57787, -3.01519, -3.01519.
```

There are `127`, `2667`, `11811`, `11811`, `2667`, and `127` nontrivial
subspaces in dimensions one through six; all were checked.  Thus at
`beta=4` a negative certificate first appears at dimension five, and at
`beta=8` it first appears at dimension six.  Dimension six is `N-2` and
therefore exactly recovers the even--even convolution, explaining the
equality of the last two entries.  At child order three, a dimension-two
quotient already gives a negative certificate at `beta=4,8`; the required
fraction is not monotone at these tiny orders.

This is a finite falsifier, not evidence for a linear asymptotic dimension
lower bound.  Its useful conclusion is narrower: even after optimizing over
all labelled modes, the first actual example beyond the one-mode audit does
not reveal a small fixed-dimensional closure, and the low-temperature
order-four case needs the entire `N-2` effective character space.

## 7. Verdict

The labelled Fourier route does not collapse to a scalar entropy or a
separately paid channel.  It gives an exact hierarchy of genuinely joint
certificates.  But bridge Gibbs smoothness proves that a sub-square-root
quotient cannot cancel a linear defect, and the endpoint of the hierarchy
is the original bridge optimization itself.  No theorem currently closes
the intermediate window (5.4) for actual optimizing children.  Therefore
the negative-mode idea does not produce `E=o(N)` in this campaign; it gives
a scalable quantitative obstruction and identifies the precise amount of
labelled state any continuation of this route must carry.
