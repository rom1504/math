# Wave 42: mixed principal-norm pressure and replica collapse

## 1. Outcome

The replica clause in (10.1125) is unnecessary.  Put

```math
D_S=Q(A[S])-\sigma C_S\ge0,
\qquad
M_S=H-D_S,
\qquad H=B_{n,m}+t.
\tag{R42.1}
```

Since `M_S<=H`, the first mixed pressure controls both its replica gap and
the hard favorable mass pointwise.  If `H>0`, `theta H<=bL_0`, and
`K_z(theta)>=aL_0` for `0<a<b`, then

```math
U_m\{\delta_S^{abs}(z)\le H\}
\ge\exp\{-(b-a+o(1))L_0\}.
\tag{R42.2}
```

No Paley--Zygmund or separate second-replica hypothesis is needed.  At fixed
density `H=Theta(n^(3/2))`, the natural pressure temperature is therefore

```math
\theta=\frac{bL_0}{H}=\Theta(n^{-3/4-c}),
\tag{R42.3}
```

not the larger `Theta(n^(-3/4))` initially proposed in (10.1125).

The remaining one-pressure input is still open and is essentially the strict
soft exceptional-center target.  Exact selector exchange makes
`S -> Q(A[S])` a one-sided self-bounding function, but its audited fixed-slice
variance proxy is `O(nq_n)=O(n^(5/2))`.  Even an optimistic quadratic entropy
conversion at (R42.3) costs `n^(1-2c)`, larger than
`L_0=n^(3/4-c)` by `n^(1/4-c)`.  More importantly, such an upper fluctuation
bound cannot manufacture positive pressure without the forbidden first-mean
input.  The full oriented deficit is not self-bounding: exact finite
minimizers already violate the tempting constant-two deletion bound.

## 2. Exact two-replica identities

Let

```math
Z(\theta)=\mathbb E_{U_m}e^{-\theta D_S},
\qquad
L(\theta)=\log Z(\theta),
\qquad
K(\theta)=\theta H+L(\theta),
\tag{R42.4}
```

and let `pi_theta(S)=U_m(S)e^(-theta D_S)/Z(theta)`.  If
`N=binom(n,m)`, direct calculation gives

```math
\boxed{
\begin{aligned}
K(2\theta)-2K(\theta)
&=L(2\theta)-2L(\theta)\\
&=D_2(\pi_\theta\Vert U_m)\\
&=\log\!\left(N\sum_S\pi_\theta(S)^2\right).
\end{aligned}
}
\tag{R42.5}
```

Thus the proposed replica gap is exactly the Renyi collision localization of
the `Q_S-sigma C_S`-weighted selector law.  Since
`L''(u)=Var_(pi_u)(D_S)`, it also has the exact curvature representation

```math
\boxed{
\begin{aligned}
K(2\theta)-2K(\theta)
={}&\int_0^\theta u\,\operatorname{Var}_{\pi_u}(D_S)\,du\\
&+\int_\theta^{2\theta}(2\theta-u)
\operatorname{Var}_{\pi_u}(D_S)\,du.
\end{aligned}
}
\tag{R42.6}
```

These identities are exact but do not themselves control localization.

## 3. Pointwise cap removes the replica condition

Let `W=e^(theta M_S)`.  From `D_S>=0`,

```math
0<W\le e^{\theta H},
\qquad
W^2\le e^{\theta H}W.
\tag{R42.7}
```

The second inequality and (R42.5) immediately give

```math
\boxed{
K(2\theta)-2K(\theta)
\le\theta H-K(\theta).
}
\tag{R42.8}
```

Hence `K(theta)>=aL_0` and `theta H<=bL_0` automatically give the old replica
bound with constant `b-a`.

There is a sharper direct conversion.  On `{M_S<0}`, `W<1`, while everywhere
`W<=e^(theta H)`.  If `p=P{M_S>=0}`, then

```math
e^{K(\theta)}=\mathbb EW
\le1+p(e^{\theta H}-1).
```

Therefore

```math
\boxed{
P\{M_S\ge0\}
\ge\frac{e^{K(\theta)}-1}{e^{\theta H}-1}.
}
\tag{R42.9}
```

Because `|C_S|>=sigma C_S`, `M_S>=0` implies
`delta_S^abs=Q_S-|C_S|<=H`.  Equations (R42.2)--(R42.3) follow.  This is the
exact strict soft-to-hard argument in its simplest form; it shows that the
only unresolved part of (10.1125) was always the first pressure.

If `H<0`, the target event `delta_S^abs<=H` is empty.  If `H=0`, positive
pressure is impossible because `M_S=-D_S<=0`; zero-deficit incidence must be
handled by an atom or an infinite-temperature limit, and the fixed orientation
can miss grounds in the opposite sector.  If `H` is positive but much smaller
than its fixed-density `Theta(n^(3/2))` scale, (R42.9) remains exact but
`theta=L_0/H` may be too large for perturbative selector estimates.

## 4. Exact selector-exchange formula

Let adjacent selectors be `S=U union {a}` and `T=U union {b}`.  Define

```math
\delta_a(S)=Q(A[S])-Q(A[U]),
\qquad
c_a(U)=2\sum_{i\in U}A_{ai}z_az_i,
\tag{R42.10}
```

and similarly for `b,T`.  Principal cap is monotone under adding vertices,
so both deletion losses are nonnegative.  Since
`C_S=C_U+c_a(U)`, one has the exact port decomposition

```math
\boxed{
\begin{aligned}
D_S&=D_U+\delta_a(S)-\sigma c_a(U),\\
D_T-D_S
&=[\delta_b(T)-\sigma c_b(U)]
-[\delta_a(S)-\sigma c_a(U)],\\
\frac{\pi_\theta(T)}{\pi_\theta(S)}
&=\exp\{-\theta(D_T-D_S)\}.
\end{aligned}
}
\tag{R42.11}
```

Thus every tilted Johnson edge is a difference of two core ports.  The
deletion term is one-sided and stable; the center-field term has no matching
sign.  Existing global minimality does not compare those two terms.

## 5. Principal norm is deletion-self-bounding

Let a sector `epsilon` and spin `y` attain `Q(A[S])`, and put

```math
\ell_i=\epsilon y_i(A[S]y)_i.
```

Single-spin stability of an absolute ground and summation of its local fields
give

```math
\ell_i\ge0,
\qquad
\sum_{i\in S}\ell_i=Q(A[S]).
\tag{R42.12}
```

Restricting the same ground after deleting `i` proves

```math
\boxed{
0\le\delta_i(S)\le2\ell_i,
\qquad
\sum_i\delta_i(S)\le2Q_S,
\qquad
\sum_i\delta_i(S)^2\le4(m-1)Q_S.
}
\tag{R42.13}
```

For an adjacent swap through `U`,

```math
(Q_S-Q_T)_+\le Q_S-Q_U=\delta_a(S).
\tag{R42.14}
```

Average directed Johnson edges uniformly.  Reversal symmetry and (R42.13)
then give

```math
\boxed{
\mathbb E_{S\to T}(Q_S-Q_T)^2
\le8\,\mathbb E_SQ_S.
}
\tag{R42.15}
```

The exact Johnson spectral gap is `n/[m(n-m)]`; hence

```math
\boxed{
\operatorname{Var}_{U_m}(Q_S)
\le\frac{4m(n-m)}n\,\mathbb E Q_S
=O(nq_n)=O(n^{5/2}).
}
\tag{R42.16}
```

This is a genuine minimizer-compatible structural estimate, but it is too
large for the project pressure.  A fixed-slice entropy/Herbst argument using
only the same one-sided carré du champ has quadratic proxy `nq_n`.  At
`theta=Theta(n^(-3/4-c))`, its scale is

```math
\theta^2nq_n=O(n^{1-2c})
=n^{1/4-c}L_0.
\tag{R42.17}
```

At the older `theta=Theta(n^(-3/4))`, the miss is `n^(1/4+c)L_0`.
Moreover, (R42.16)--(R42.17) are upper fluctuation estimates.  Integrating
curvature cannot lower-bound `K(theta)` above zero unless one supplies a
boundary value; linearizing that boundary is precisely the circular
`K'(0)=t-G-p_2Delta` mean.

An improvement of the deletion-square statistic to the rough tilted scale

```math
\mathbb E_{\pi_u}\sum_{i\in S}\delta_i(S)^2
=O(H^2/L_0)=O(n^{9/4+c}),
\qquad0\le u\le\theta,
\tag{R42.18}
```

would repair the quadratic entropy exponent at the new temperature.  It
still would not create positive pressure by itself; a noncircular lower
boundary or mixed alignment theorem would remain necessary.

## 6. The full deficit is not self-bounding

It is tempting to replace `Q_S` in (R42.13) by
`D_S=Q_S-sigma C_S`.  This is false even inside the audited exact minimizers.
For `A_8`, take

```text
z=(1,-1,-1,-1,-1,-1,-1,-1),
sigma=-1,
S=(0,1,3,4,5,6).
```

Then `Q_S=18`, `C_S=-14`, and `D_S=4`, while the six deletion increments
`D_S-D_(S\setminus i)` are

```text
[-4, 4, -4, 4, 4, 4].
```

Thus

```math
\sum_i(D_S-D_{S\setminus i})_+=16=4D_S.
\tag{R42.19}
```

`A_9` contains the same ratio.  The negative increments show exactly why
the signed center ports in (R42.11) defeat the monotone deletion proof.
This is a finite mechanism wall, not an asymptotic counterexample; it rules
out the natural constant-two self-bound already on exact minimizers.

## 7. Conclusion and verification

`tmp/mixed_pressure_r42_check.py` verifies:

- (R42.12)--(R42.16) for every selector size in `A_6,A_8,A_9`;
- the Renyi, cap, and direct soft-to-hard identities numerically;
- every adjacent-selector port identity (R42.11) on `A_8,m=4`; and
- the exact deficit-deletion wall (R42.19).

The useful Wave 42 correction is that mixed pressure is a **one-input**
criterion: delete the replica clause from (10.1125) and use (R42.9).  This
does not yet advance the lower bound on that pressure.  Principal-cap
self-bounding misses the exponent and controls the wrong direction, while
the mixed deficit loses monotonicity.  A successor must directly prove a
strict lower bound on

```math
\mathbb E_{U_m}\exp\{-\theta[Q(A[S])-\sigma C_S]\}
```

for one low-row center, using a nonlinear actual-minimizer alignment not
encoded by uniform means or generic selector smoothness.
