# Fair-base replicated smoothing for the full carrier certificate

Status: **rigorous carrier-independent interpolation identity**.  The
fair-base capped overlap `rho_N(S)` in CT.7 is, up to the explicit cap tail
and the `O(N^(-1/2))` cavity correction, exactly one normalized secant of a
scalar Gaussian-smoothing derivative.  The continuous value at replica
number one is a mixed derivative.  This extends CI.1 from the spiked
marginal to the complete bounded-`L^2` square carrier without introducing a
carrier prior, net, or response table.

Existing contracted-temperature minimality proves the tail estimate needed
to make the identity uniform on every declared compact tilt window.  It does
not prove the derivative rigidity, because the smoothing changes parent
cross edges and the replica-one derivative cancels the overlap identically.

## 1. Fair-base cap tail on a compact tilt window

Retain the setup of CT.1.  Thus `A,D` are actual
contracted-temperature minimizing children of orders `m,n`, `N=m+n`,
`t=beta/sqrt(N)`, and `L` is their exact bridge pressure.  Under the fair
bridge law `U=U_(mn)`, the exact fixed-replica estimate RT.3 says, for every
positive integer `k`,

```math
 \log E_Ue^{kL}
 \le k\{\log2+p_A(t)+p_D(t)\}
      +{k^2t^2mn\over2}
 \le A_{\beta,k}N,                                \tag{FB.1}
```

where, for example, one may take

```math
 A_{\beta,k}=k\log2+{\beta^2k\over4}
                         +{\beta^2k^2\over8}.     \tag{FB.2}
```

Fix a compact tilt window `S<infinity`, choose an integer `k>S`, and choose

```math
 C>\max\left\{{3\beta^2\over8}+2\log K,
               {A_{\beta,k}\over k-S}\right\}.  \tag{FB.3}
```

The first condition is the full-carrier response-truncation threshold from
CT.1; the second controls the retained mass uniformly on the tilt window.
Increasing the cap to satisfy both only decreases the truncation error.  Put

```math
 F=L\wedge CN,
 \qquad
 {d\Pi_{s,F}\over dU}={e^{sF}\over E_Ue^{sF}}.    \tag{FB.4}
```

**Lemma FB.1 (uniform retained-mass estimate).**  There is
`c=c(beta,k,S,C)>0` such that

```math
 \boxed{
 \sup_{|s|\le S}\Pi_{s,F}\{L\ge CN\}\le e^{-cN}.}              \tag{FB.5}
```

*Proof.*  Markov's inequality and (FB.1) give

```math
 U\{L\ge CN\}\le e^{(A_{\beta,k}-kC)N}.          \tag{FB.6}
```

For `0<=s<=S`, the numerator of the capped tilted mass is
`e^(sCN)U{L>=CN}`, while the denominator is at least one because `F>=0`.
For `-S<=s<0`, the denominator is at least `e^(sCN)`, so the tilted mass is
at most the fair mass.  Condition (FB.3) proves (FB.5). `square`

The result is uniform over the actual children and orientation.  Its only
optimizer input is the annealed child pressure bound used in RT.3.

## 2. One scalar smoothing curve

Extend the bridge coefficients to real matrices `h` and write

```math
 L(h)=\log E_{x,z}\cosh\left(t\{H_A(x)+\epsilon H_D(z)
                                  +x^{\mathsf T}hz\}\right).    \tag{FB.7}
```

Let `G` have iid standard Gaussian entries and be independent of the fair
sign bridge `B`.  Choose `CN` outside the finite set of sign-bridge pressure
values; an arbitrarily small increase of `C` ensures this.  For `s ne0`,
define

```math
 \mathscr H_{s,F}(u)
 ={1\over s}\log E_{B,G}
 \exp\left\{s\min(L(B+\sqrt uG),CN)\right\},
 \qquad u\ge0,                                    \tag{FB.8}
```

and use the continuous expectation-valued extension at `s=0`.

Let

```math
 m_a(B)=E_{\nu_B}[\tau X_iZ_j],
 \qquad
 \mathcal O(B)={1\over mn}\sum_a m_a(B)^2.       \tag{FB.9}
```

Thus

```math
 \mathcal O(B)
 =E_{\nu_B^{\otimes2}}[\tau^1\tau^2R_XR_Z].      \tag{FB.10}
```

Define two scalar functions

```math
 A_N(s)={2\mathscr H_{s,F}'(0+)\over t^2mn},
 \qquad
 p_N(s)=\Pi_{s,F}\{L<CN\}.                       \tag{FB.11}
```

**Theorem FB.2 (exact fair-base smoothing identity).**  For every real
`s`, with the formula interpreted continuously at `s=0`,

```math
 \boxed{
 A_N(s)
 =E_{\Pi_{s,F}}
 \left[1_{\{L<CN\}}\{1+(s-1)\mathcal O(B)\}\right].}          \tag{FB.12}
```

Consequently the continuous secant

```math
 \mathcal Q_N(s)=
 \begin{cases}
 \displaystyle {A_N(s)-p_N(s)\over s-1},&s\ne1,\\[6pt]
 \displaystyle \partial_s\{A_N(s)-p_N(s)\}\big|_{s=1},&s=1
 \end{cases}                                      \tag{FB.13}
```

is exactly the retained overlap:

```math
 \boxed{
 \mathcal Q_N(s)
 =E_{\Pi_{s,F}}[1_{\{L<CN\}}\mathcal O(B)].}      \tag{FB.14}
```

*Proof.*  The heat-semigroup derivative at `u=0+` is one half of the
Euclidean bridge Laplacian.  At `L(B)<CN`, Gibbs differentiation gives

```math
 \partial_{h_a}L=t m_a,
 \qquad
 \partial_{h_a}^2L=t^2(1-m_a^2).                 \tag{FB.15}
```

Therefore

```math
 \Delta_he^{sL}
 =st^2e^{sL}\left\{mn+(s-1)\sum_am_a^2\right\}. \tag{FB.16}
```

At `L(B)>CN`, the capped function is locally constant, so its Laplacian
vanishes.  Averaging (FB.16), dividing by `s`, the partition function, and
`t^2mn/2` proves (FB.12).  Subtracting `p_N(s)` and dividing by `s-1`
proves (FB.14) away from one.  Differentiation at one gives its continuous
value and proves the second line of (FB.13). `square`

The formula at replica number one deserves emphasis:

```math
 A_N(1)=p_N(1),
 \qquad
 \mathcal Q_N(1)=
 \partial_s\{A_N(s)-p_N(s)\}|_{s=1}.             \tag{FB.17}
```

The overlap cancels completely from the one-replica smoothing derivative
and reappears only in the mixed replica/smoothing derivative.

## 3. Exact equivalence with the CT overlap target

Let `r_a(B_-a)` be the edge-cavity response and let

```math
 \rho_N^{\rm CT}(S)
 ={1\over mn}\sup_{|s|\le S}
 E_{\Pi_{s,F}}\sum_ar_a^2                        \tag{FB.18}
```

be exactly CT.7.  Edge insertion FI.17 gives, pointwise in `B`,

```math
 \left|{1\over mn}\sum_ar_a^2-\mathcal O(B)\right|
 \le {2\tanh t\over1-\tanh t}=O_\beta(N^{-1/2}).              \tag{FB.19}
```

Equations (FB.5), (FB.14), and (FB.19) yield

```math
 \boxed{
 \left|\rho_N^{\rm CT}(S)
       -\sup_{|s|\le S}\mathcal Q_N(s)\right|
 \le O_\beta(N^{-1/2})+e^{-cN}.}                 \tag{FB.20}
```

Thus, on every declared compact window,

```math
 \boxed{
 \rho_N^{\rm CT}(S)=o(1)
 \quad\Longleftrightarrow\quad
 \sup_{|s|\le S}\mathcal Q_N(s)=o(1).}          \tag{FB.21}
```

Equivalently, the carrier-independent derivative-rigidity target is

```math
 \boxed{
 \begin{aligned}
 &\sup_{\substack{|s|\le S\\s\ne1}}
 {\left|A_N(s)-p_N(s)\right|\over|s-1|}=o(1),\\
 &\left|\partial_s\{A_N(s)-p_N(s)\}|_{s=1}\right|=o(1).
 \end{aligned}}                                   \tag{FB.22}
```

This includes `s=1` without dividing by a vanishing secant.  It consists of
one scalar replicated smoothing curve; `p_N(s)=1+O(e^(-cN))` is an
independently certified tail correction.  No carrier labels or coefficients
occur.

For a growing window `S_N`, the same theorem remains valid if one chooses
`k_N>S_N` and a cap coefficient satisfying the carrier threshold and

```math
 \delta_N:=(k_N-S_N)C_N-A_{\beta,k_N},
 \qquad \delta_NN\longrightarrow\infty.          \tag{FB.23}
```

For an exponential tail one may require `delta_N>=c_0>0`; for example take
`k_N=2S_N+O(1)` and
`C_N>=2A_(beta,k_N)/(k_N-S_N)`, together with
`C_N>3beta^2/8+2log K`.  This permits
`C_N=O_beta(1+S_N)`.  Enlarging the cap does not harm full-carrier response
truncation.  With a fixed cap and an unbounded positive tilt window,
however, FB.5 need not hold: the capped tilt can concentrate on the plateau
`{L>=CN}`.  Therefore the fixed-cap identity is an exact compact-window
equivalence, while recurrence-scale growing windows require the declared
adaptive cap with quantitative slack (or separate plateau control).  This
qualification is essential.

## 4. What minimality does and does not prove

Contracted-temperature minimality proves RT.2 and hence the fixed-replica
moment/tail input FB.1--FB.5.  Exact Gibbs edge insertion, independently of
minimality, gives the cavity/full-Gibbs comparison FB.19.  Minimality does
**not** prove FB.22.

The reason is visible in the exact identity.  At `s=1`, the overlap cancels
from `A_N(1)` for every pair of children, so the one-replica annealed
minimality bound is algebraically blind to it.  For `s ne1`, FB.22 asks for
`o(N)` precision in the unnormalized derivative
`mathscr H_(s,F)'(0+)` around its `Theta(N)` baseline.  Present child
minimality supplies only `O(N)` pressure and moment bounds.  Finally, the
Gaussian smoothing in FB.8 perturbs parent cross edges, whereas child
minimality compares internal signings at the two child orders.  There is no
proved variational inequality between those operations.

Accordingly FB.22 is a genuine carrier-independent scalar SML, but it is
not a consequence of current minimality.  Proving it requires a replicated
external-field comparison or a synchronization principle transferring
internal child rigidity to cross-edge overlap.  Repeating one-replica
annealing, internal sign flips, or fixed-cap all-tilt concentration cannot
supply the missing derivative precision.
