# Adaptive physical response from a quadratic high-energy constraint

2026-09-17. Director extension of the Bernoulli track's signed-ground
covariance argument. **Verified by director derivation and independent
reconstruction by the Bernoulli and localization researchers.** The
third independent adversarial review is in progress.
This concerns genuine hollow full signings. It assumes a bounded spectral
norm of the ACTUAL signing; it does not claim that exact minimizers have
that property. No external novelty claim is made.

## 1. Statement

Write kappa=sqrt(2/pi). Fix L>0 and 0<c<=L/2. Let A be a symmetric
hollow full signing of order n with ||A||op<=L sqrt(n), and let

```math
 C=\{x\in\{\pm1\}^n:|H_A(x)|\ge c n^{3/2}\}.
```

If C is nonempty, there is a centered EXACTLY isotropic physical sign
law nu such that uniformly for x in C,

```math
 E_\nu |h\cdot x|/\sqrt n
 \le \kappa-\delta(c,L)+O(n^{-1/6}\sqrt{\log n}).       \tag{1}
```

The error constant is absolute. An explicit positive choice is

```math
 a=(c/(4L))^2,\quad b=1/2-\pi/12,\quad
 t=\min\{1/12,c/\sqrt{432}\},
 \delta(c,L)=\min\{\kappa^5 b^2a^2/8,
                         \kappa^5t^2c^4/32\}.          \tag{2}
```

Moreover nu is a mixture of centered Gaussian-sign laws whose latent
correlation matrices have spectrum in[1/2,3/2]. Hence every linear
response is subGaussian with variance proxy(3/2)||w||_2^2. No support-size
bound is claimed. The law need not be radial in A or a scalar energy tilt.

## 2. Reusable scalar ingredient and exact mixture identity

For a correlation R with (1/2)I<=R<=(3/2)I, h=sign(G), G~N(0,R),
the archived Gaussian-boundary replacement theorem gives UNIFORMLY for
Boolean x

```math
 E|h\cdot x|/\sqrt n
 =\kappa\sqrt{x^T S_R x/n}+O(n^{-1/6}\sqrt{\log n}),
 \qquad S_R=(2/\pi)\arcsin[R].                       \tag{3}
```

Arcsine is entrywise, including diagonal pi/2. This theorem does not
require small off-diagonal correlations. The precise reconstruction and
parameter selection are in
[the Gaussian response audit](paper_localization_gaussian_sign_response_2026_09_17.md).
Gaussian Holder gives the stated linear subGaussian proxy. Independent
audit must check these imported hypotheses, not just a previous verdict.

For a symmetric hollow H with ||tH||op<=1/2, take the equal mixture of
the two Gaussian-sign laws R_plus/minus=I plus/minus tH. This law is
centered and EXACTLY isotropic, because the off-diagonal arcsines cancel.
Put V=(2/pi)arcsin[tH], whose diagonal is zero. Its two sign covariance
matrices are I plus/minus V, so |x^T Vx/n|<=1. Formula(3) yields

```math
 E_\nu |h\cdot x|/\sqrt n
 =\kappa f(v_x)+O(n^{-1/6}\sqrt{\log n}),
 f(v)=(\sqrt{1+v}+\sqrt{1-v})/2\le1-v^2/8.          \tag{4}
```

The elementary inequality follows by squaring, or from the even Taylor
series and continuity at the endpoints. Errors are uniform in H, x and
all probability measures subsequently used.

## 3. Dual distribution and the high-covariance branch

Fix ANY probability mu on C. Let

```math
 \Sigma=E_\mu xx^T,\quad \sigma(x)=\operatorname{sign}H_A(x),
 K=E_\mu\sigma(x)xx^T,\quad s=E_\mu\sigma,\quad K_0=K-sI.
```

Then diag(Sigma)=1, diag(K_0)=0, -Sigma<=K<=Sigma, and
<A,K_0>=2E_mu|H_A|>=2c n^(3/2).
Let Q be the spectral projector of Sigma on eigenvalues STRICTLY above2,
P=I-Q, r=rank Q, and m=tr(Sigma Q). Thus r<=m/2.

If m>=an, use H=Q-diag(Q) and coefficient1/2 in(4).
Its norm is at most1: both 0<=Q<=I and 0<=diag(Q)<=I, so
-I<=Q-diag(Q)<=I. For V=kappa^2 offdiag(arcsin[Q/2]),

```math
 E_\mu v_x
 =\frac{\kappa^2}{n}
    \left(\operatorname{tr}(\Sigma\arcsin[Q/2])
                   -\sum_i\arcsin(Q_{ii}/2)\right)
 \ge\frac{\kappa^2}{n}(m/2-(\pi/6)r)
 \ge\kappa^2 b a.                                   \tag{5}
```

Every odd entrywise power of the PSD matrix Q is PSD by the Schur
product theorem, and all arcsine coefficients are nonnegative. Hence
the first trace is at least m/2. Convexity on[0,1] gives
arcsin(u/2)<=u arcsin(1/2), bounding the diagonal sum by(pi/6)r.
Using E v^2>=(E v)^2 in(4) gives the first discount in(2).

## 4. Low-covariance-mass branch: signed energy survives truncation

Suppose m<an. Put B=PK_0P. Because -P Sigma P<=PKP<=P Sigma P<=2P,
we have ||B||op<=3. Furthermore

```math
 |\langle A,K_0-B\rangle|
 \le \|A\|_{op}(2\sqrt{nm}+r)
 \le L n^{3/2}(2\sqrt a+a/2)<c n^{3/2}.             \tag{6}
```

For detail, K-PKP is the average of sigma times
(Qx)x^T+(Px)(Qx)^T. Each term's pairing with A is bounded by
||A||op ||Qx|| ||x||, and Cauchy--Schwarz gives2||A||op sqrt(nm).
The remaining signed diagonal term is -sQ and costs at most||A||op r.
For a=(c/(4L))^2 and c/L<=1/2, the factor in(6) is strictly below c.
Consequently <A,B>>=c n^(3/2), and since ||A||F^2=n(n-1),

```math
 \|B\|_F^2\ge c^2 n.                               \tag{7}
```

Set H=B-diag(B). Then ||H||op<=6 and

```math
 \langle K_0,H\rangle=\operatorname{tr}(K_0PK_0P)
                     =\|B\|_F^2.                  \tag{8}
```

The removed diagonal pairs to zero with K_0. For t in(2), ||tH||<=1/2.
The compression identity X circle Y=V^*(X tensor Y)V gives
||X circle Y||op<=||X||op||Y||op for any real symmetric matrices.
The absolutely convergent arcsine series therefore gives

```math
 R=\arcsin[tH]-tH,\qquad
 \|R\|_{op}\le\arcsin(6t)-6t\le216t^3.             \tag{9}
```

R is hollow. Since K is a signed mixture of xx^T,
|<K_0,R>|=|<K,R>|<=n||R||op. Combining(7)--(9),

```math
 E_\mu\sigma(x)v_x
 =\frac{\kappa^2}{n}\langle K_0,\arcsin[tH]\rangle
 \ge\kappa^2(tc^2-216t^3)\ge\kappa^2tc^2/2.        \tag{10}
```

Cauchy--Schwarz, sigma^2=1 and(4) now give the second discount in(2).
In particular the arcsine remainder is paid in OPERATOR norm; an entrywise
sign agreement after projection is neither true nor assumed.

## 5. Minimax produces ONE physical law for the entire high-energy code

Let V_n be the convex hull of the paired Gaussian-sign laws in(4) with
hollow tH and ||tH||op<=1/2. This is a compact convex subset of the
finite-dimensional probability simplex on the physical sign cube.
All its elements are centered, exactly isotropic, and have the common
linear subGaussian proxy3/2. Sections3--4 prove that for EVERY mu on C,
there exists nu in V_n with average response at most the right side of(1).
Finite-dimensional minimax therefore gives

```math
 \min_{\nu\in V_n}\max_{x\in C}E_\nu|h\cdot x|/\sqrt n
 =\max_{\mu\in\mathcal P(C)}\min_{\nu\in V_n}
          E_{\mu,\nu}|h\cdot x|/\sqrt n
 \le\kappa-\delta(c,L)+O(n^{-1/6}\sqrt{\log n}).
```

This is an existence theorem, not an efficient recovery algorithm. It
does not encode or presume the unknown optimum M_n. Its bound is uniform
on an entire macroscopic high-energy set, but the positive discount in(2)
is small and need not meet the slope required for original convergence.

## 6. Exact limits of the conclusion

- The signing norm assumption is substantive and is NOT proved for exact
  minimizers. The bounded-covariance dual lemma alone did not require it.
- Isotropy is exact; physical signs and the scalar absolute response are
  retained throughout, not replaced by vector values or separate channels.
- Minimax produces a single law but no polynomial support-size bound.
- Sampling many bridge columns still incurs old-code selection costs.
  Microscopic preparation can control those costs, but not automatically
  exclude lower-energy old words.
- The original interval and convergence question are unchanged.
