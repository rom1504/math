# Wave 35 main attack: direct arbitrary-cut hinge free energy

Status: the variational and cube identities below are **Verified**, including
independent exhaustive checks on `A_6,A_8,A_9`.  The resulting entropic-hinge
condition is an **Open sufficient lemma**.  The linear row-penalty conclusion
is a scoped mechanism wall, not a falsifier of (10.795).

## 1. Exact existential hinge dual

Fix `A,m`, the uniform selector law `U=U_m`, a row cap `C`, and

```math
h_d(S)=\widehat\ell(S,d),
\qquad \phi_t(u)=[u-t]_+,
\qquad
K_{\lambda,t}(d)=\mathbb E_U e^{-\lambda\phi_t(h_d(S))}.
```

Let `D_C={d:R_2(d)<=C}`.  The finite Gibbs variational principle, followed
only by commuting two maxima, gives

```math
\boxed{
\log\max_{d\in D_C}K_{\lambda,t}(d)
=\sup_{w\in\Delta(\binom{[n]}m)}
\left\{-D(w\Vert U)
-\lambda\min_{d\in D_C}\mathbb E_w\phi_t(h_d(S))\right\}.
}
\tag{M35.1}
```

For fixed `d`, the maximizing selector law is exactly

```math
\pi_d(S)=\frac{U(S)e^{-\lambda\phi_t(h_d(S))}}
                    {K_{\lambda,t}(d)},
\qquad
-\log K_{\lambda,t}(d)
=D(\pi_d\Vert U)+\lambda\mathbb E_{\pi_d}\phi_t(h_d).
\tag{M35.2}
```

This is existential, not adversarial: only one low-entropy selector tilt and
one responding low-row cut are needed.  For every buffer `u>0`,

```math
\boxed{
K_{\lambda,t}(d)
\le U\{S:h_d(S)\le t+u\}+e^{-\lambda u}.
}
\tag{M35.3}
```

Consequently, if for every relevant `(A,n,m)` there exist `w,d,t,u`, with
uniform constants `0<A<Lambda`, satisfying the package

```math
D(w\Vert U)+\lambda\mathbb E_w\phi_t(h_d)\le A T L_0,
\qquad
\lambda u\ge\Lambda T L_0,
\qquad d\in D_C
\tag{M35.4}
```

give

```math
U\{h_d\le t+u\}
\ge e^{-A T L_0}-e^{-\Lambda T L_0}
=e^{-(A+o(1))T L_0}.
\tag{M35.5}
```

At the Wave 34 scales one may take `t+u=O(n^(3/2-c'))`,
`c'=c_0-eta`, and `lambda=Theta(n^(-3/4))`.  With
`C=2n(n-1)`, (M35.5) is (10.795) and proves convergence.

The condition is exponent-equivalent to the direct tail, not a proof of it:
if a cut covers a set `G`, conditioning `U` on `G` gives
`D(U(.|G)||U)=log(1/U(G))` and zero hinge on `G`.  Its value is that it
isolates a possible softer interface for a conditional replacement theorem:
an **existential** entropy budget plus expected positive-part loss.

## 2. Why the unclipped exponential is too lossy

Since `h_d=Y_A-X_d`, the same calculation without `phi_t` is completely
linear:

```math
\boxed{
\log\max_{d\in D_C}\mathbb E_U e^{-\lambda h_d(S)}
=\sup_w\left\{
\lambda\left[\max_{d\in D_C}\langle H_w,d\rangle-Y_w\right]
-D(w\Vert U)
\right\}.
}
\tag{M35.6}
```

At a joint optimizer, `d` is a constrained maximizing cut of `H_w`, with
`w(S) proportional to exp(-lambda h_d(S))`.  Thus the direct soft optimizer
has an exact self-consistent selector-response interpretation through
(10.821)--(10.822).

However,

```math
h_d(S)=\delta_S(d)+p_2\langle A,d\rangle-p^{3/2}q_n
\ge-(p_2+p^{3/2})q_n.
\tag{M35.7}
```

Converting the raw exponential back to a hard upper tail can therefore pay
`Theta(lambda n^(3/2))=Theta(n^(3/4))`, larger than
`T L_0=o(n^(3/4))`.  The positive-part clipping in (M35.1) is essential.

## 3. Exact block increments and Euler system

Write `d=sigma xx^T`.  For `H subset [n]`, flip the physical spins in `H`
and put

```math
W_H=\sum_{i\in H,j\notin H}A_{ij}x_ix_j,
\qquad
W_{H,S}=\sum_{i\in H\cap S,j\in S\setminus H}A_{ij}x_ix_j.
```

Direct expansion gives

```math
\boxed{
h_{d^H}(S)-h_d(S)=4\sigma[W_{H,S}-p_2W_H],
}
\tag{M35.8}
```

and

```math
\boxed{
R_2(d^H)-R_2(d)
=-4\sum_{i\in H}x_i(A^2x)_i
+4\lVert A[:,H]x_H\rVert_2^2.
}
\tag{M35.9}
```

For the `n` auxiliary one-vertex involutions (not the coordinate basis),

```math
\sum_i[h_{d^i}(S)-h_d(S)]=4X_d(S),
\qquad
\sum_i[R_2(d^i)-R_2(d)]=4[n(n-1)-R_2(d)].
\tag{M35.10}
```

Orientation reversal leaves `R_2` unchanged and changes the loss by
`h_{-d}-h_d=2X_d`.

For the linearly row-penalized objective

```math
\Psi_{\lambda,\gamma}(d)
=-\log K_{\lambda,t}(d)+\gamma R_2(d),
```

every global minimizer, with posterior `pi_d` from (M35.2), satisfies for
every block `H`

```math
\boxed{
\log\mathbb E_{\pi_d}
\exp\{-\lambda[\phi_t(h_{d^H})-\phi_t(h_d)]\}
\le\gamma[R_2(d^H)-R_2(d)].
}
\tag{M35.11}
```

The analogous orientation inequality has zero right side.  This is exact,
but summing it does not bound the selector KL in (M35.2); the hinge also
destroys the raw linear closure in (M35.10).

## 4. Smooth row enforcement is off scale

For a one-vertex flip,

```math
|h_{d^i}(S)-h_d(S)|
\le L_i:=4[(m-1)+p_2(n-1)].
\tag{M35.12}
```

Thus `|-log K(d^i)+log K(d)|<=lambda L_i`.  If
`R_2(d)>2n(n-1)`, (M35.10) supplies an `i` with
`Delta_i R_2<-4(n-1)`.  Therefore the black-box bit-descent choice

```math
\gamma\ge\frac{\lambda L_i}{4(n-1)}=Theta(\lambda)
\tag{M35.13}
```

forces a linear-penalty optimizer into `D_C`.  At
`lambda=Theta(n^(-3/4))`, however, it charges any ordinary
`R_2=Theta(n^2)` center by `Theta(n^(5/4))`, whereas the target logarithmic
budget is `T L_0=o(n^(3/4))`.  A target-scale bound on the whole penalized
objective would instead have to produce the much stronger
`R_2=O(TL_0/lambda)=O(Tn^(3/2-c_0))`; generic bit descent does not supply
such a center.  This is a scale mismatch for the black-box linear-penalty
proof, not an impossibility theorem.  A hard or shifted barrier avoids paying
inside `D_C` but then reduces to constrained optimization and loses
informative outgoing Euler directions.  Extra signing-minimizer structure is
required.

## 5. Finite audit

`tmp/direct_cut_free_energy_r35.py` exhaustively constructs the cut/selector
data and verifies (M35.2)--(M35.3) and all one-bit/orientation instances of
(M35.11) on `A_6,A_8,A_9`; it also checks (M35.8)--(M35.9) on 300
reproducibly sampled block flips.  The algebraic derivation covers every
block, and the checker passes.

The data also delimit surrogate shortcuts.  At `A_8,m=4,t=0,lambda=.05`,
the unpenalized soft optimizer covers `0.257143` of the slice while the best
hard cut covers `0.571429`.  At `A_9,m=5,lambda=.2`, adding normalized linear
row penalty `gamma=.1` moves from `(R_2,coverage)=(112,.674603)` to
`(16,.126984)`.  All `A_9` cuts in this audit already lie under the Wave 34
cap `2n(n-1)=144`.  These are finite numerical mechanism walls only: small
soft temperature and gratuitous smooth row pressure need not preserve the
hard-tail optimizer.

## 6. Frontier contribution

The direct-cut route now has a clean nonlinear response target: prove (M35.4)
for one low-entropy selector tilt and one genuinely low-row cut, with a strict
soft-to-hard margin.  Uniform first moments remain circular, and the bare
cube Euler system plus a generic linear row penalty cannot prove it.  The
most plausible sources of the missing selector entropy compression are the
projective agreement mechanism and the conditional multi-selector
replacement game, not another average-loss calculation.
