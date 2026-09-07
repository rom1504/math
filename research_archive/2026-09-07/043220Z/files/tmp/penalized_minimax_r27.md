# Wave 27 root memo: every linear row-penalized learner is circular

## Status

The constrained and priced identities below are **proved** and are checked
on `A_6,A_8,A_9`, at one and two deletions, by
`tmp/penalized_minimax_r27.py`.  They falsify the proposed *linear-mean*
penalized all-cut route as an independent way to prove (10.795).  They do
not rule out a nonlinear coverage, positive-part, or Laplace theorem.

## Exact centered response matrix

Let `w` be any law on the `m`-slice and put

```math
r_{ij}(w)=\Pr_{S\sim w}\{i,j\in S\},
\qquad
H_w=A\circ(R_w-p_2\mathbf 1_{\ne}),
```

where the diagonal is zero.  Also put

```math
Y_w=\mathbb E_{S\sim w}Q(A[S])-p^{3/2}q_n.
```

For every oriented full cut `d`, (10.806) gives the exact identity

```math
\boxed{
\mathbb E_w\widehat\ell(S,d)
=Y_w-\langle H_w,d\rangle.
}
\tag{R27.1}
```

Represent `d=\sigma xx^T` off the diagonal.  Since orientation does not
affect row square,

```math
\langle H_w,d\rangle=\sigma x^TH_wx,
\qquad
R_2(d)=\lVert Ax\rVert_2^2=x^TA^2x.
\tag{R27.2}
```

Thus the exact pure constrained value is

```math
\boxed{
\min_{d:R_2(d)\le C}\mathbb E_w\widehat\ell(S,d)
=Y_w-max_{x:x^TA^2x\le C}|x^TH_wx|.
}
\tag{R27.3}
```

This is the complete all-cut response; no parent-ground or sign condition is
hidden in it.

## Mixed constraint and price dual

Allow a law `nu` on cuts and assume the row budget is feasible.  Finite LP
duality gives

```math
\boxed{
\begin{aligned}
&\min_{\nu:\mathbb E_\nu R_2\le C}
\mathbb E_{w\otimes\nu}\widehat\ell(S,D)\\
&\quad=Y_w-inf_{\theta\ge0}
\left\{\theta C+max_x
\left[|x^TH_wx|-\theta x^TA^2x\right]\right\}.
\end{aligned}
}
\tag{R27.4}
```

Equivalently, the priced pure problem is

```math
\boxed{
\min_d\{\mathbb E_w\widehat\ell(S,d)+\theta(R_2(d)-C)\}
=Y_w-\theta C-max_x
\{|x^TH_wx|-\theta x^TA^2x\}.
}
\tag{R27.5}
```

These formulas expose exactly what a spectral or SDP relaxation of the
penalized learner would have to estimate.

## Uniform-selector collapse

For `w=U_m`, every off-diagonal inclusion probability is `p_2`, so

```math
H_{U_m}=0.
```

Therefore every fixed cut, every low-row cut, and every mixture have the
same mean:

```math
\boxed{
\mathbb E_{U_m}\widehat\ell(S,d)
=\mathbb E_{U_m}Q(A[S])-p^{3/2}q_n.
}
\tag{R27.6}
```

If `C` is at least the minimum cut row square, the mixed value (R27.4) is
exactly the right side of (R27.6); its dual chooses `theta=0`.  Consequently
any statement of the form

> for every selector law `w`, some cut or row-feasible mixture has expected
> effective loss `O(n^(3/2-c))`

already assumes the stronger uniform mean restriction lemma at `w=U_m`.
Adding `R_2`, a Lagrange price, an SDP relaxation, or arbitrary cuts cannot
change this conclusion.  Apparent gain from fixing a positive price and
using unused row budget disappears when the exact constraint dual optimizes
the price.

## Route disposition

The linear penalized minimax route should be stopped.  A viable all-cut
functional must be nonlinear in selector loss--for example a hard coverage
event, positive-part Laplace transform, or another quantity that detects a
stretched-exponential favorable tail despite the cut-independent mean.
Equations (R27.3)--(R27.5) can still be used as diagnostics or relaxations of
such a nonlinear theorem, but cannot by themselves prove it.
