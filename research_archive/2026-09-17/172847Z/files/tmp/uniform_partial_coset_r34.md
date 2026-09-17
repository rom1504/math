# Wave 34 root route: uniform partial-coset reduction

## Status

The reductions below are proved from (10.795), (10.831)--(10.832),
(10.853)--(10.854),
(10.873), (10.904), (10.925)--(10.926), and (10.947).  They do not prove
the missing signing estimate.  Their point is that the adversarial-selector
power sum and its common-weight capacity certificate are stronger than the
restriction theorem actually needs.

## 1. One coset may cover only a rare uniform-selector set

For an eligible row-good block coset `a`, let `C_a` be its oriented cut
support and retain the hit set `H_a` from (10.873).  Define

```math
E_d(t)=\{S:\widehat\ell(S,d)\le t\}.
```

By the definition of the compressed norm,

```math
H_a=\bigcup_{d\in C_a}E_d(t).
```

Consequently, for the uniform selector law `U_m`,

```math
\boxed{
\max_{d\in C_a}U_m(E_d(t))
\ge \frac{U_m(H_a)}{|C_a|}.
}
\tag{R34.1}
```

If `a` has `k` blocks, `|C_a|<=2^k`.  At the project scale
`k=O(L_0/log n)`, this costs `o(L_0)` in the logarithm.  Thus, if for some
`T<=n^eta`, `eta<c_0`,

```math
\boxed{
\max_{a\in\mathscr A_C}U_m(H_a)\ge e^{-O(TL_0)},
}
\tag{R34.2}
```

then (R34.1) gives a cut satisfying (10.795) with saving
`c'=c_0-eta>0`; the tolerance and row cap constructed at `c_0` are stronger
than those needed at `c'`.  Hence (R34.2), uniformly in the target window,
proves convergence.

Conversely, a cut satisfying (10.795) can be planted as the all-ones word of
the balanced row-good coset supplied by (10.853)--(10.854).  Its coset hit
set contains `E_d(t)`.  Therefore (R34.2) with `T=1` is exponent-equivalent
to the arbitrary-cut lemma (10.795), up to the already harmless codebook
factor.

The fractional value

```math
\delta_*=\inf_w\max_a w(H_a)
```

is at most `max_a U_m(H_a)`.  (The exact `A_9,m=5` table already separates
them: `delta_*=839/995`, while the best uniform one-coset coverage is
`124/126`.)  A lower bound on `delta_*`, and hence the
all-law power sum (10.907), is sufficient for (R34.2) but is not necessary.
This does not solve (R34.2); it removes an adversarial quantifier from the
sharp restriction target.

## 2. Only the uniform center moment is needed

Let `r,s,T,L_0,k_0` be as in (10.906), and put

```math
h_U=\max_a U_m(H_a),\qquad M=|\mathscr A_C|.
```

If `J_(r,s)` is the partial-collision event for `r` iid **uniform**
selectors, then (10.904) gives

```math
\Pr(J_{r,s})
\le \binom rs\sum_aU_m(H_a)^s
\le \binom rs M h_U^s.
\tag{R34.3}
```

Since `log M=O(n log k)=O(rL_0)` and `s=ceil(r/T)`, any uniform-law bound

```math
\Pr(J_{r,s})\ge e^{-O(rL_0)}
```

implies `h_U>=e^{-O(TL_0)}`, hence (R34.2) and convergence.  The infimum over
all selector laws in (10.907), (10.926), and (10.927) is unnecessary for
this implication.

In particular, it is enough to establish the center-star event in (10.926)
only for `w=U_m`:

```math
\boxed{
\mathbb E_{z\sim\nu_2}
\Pr_{S_1,\ldots,S_s\sim U_m}
\left\{\sum_{j=1}^sa_z(S_j)\le D\right\}
\ge e^{-O(rL_0)}.
}
\tag{R34.4}
```

## 3. The uniform high moment is exponent-equivalent to one center degree

For the hard equal-radius event, set

```math
p_d(z)=U_m\{S:a_z(S)\le d\},
\qquad p_*=\max_{z\in\mathcal C_2}p_d(z).
```

Because `nu_2` is uniform and `|C_2|<=2^n`, exactly

```math
\boxed{
2^{-n}p_*^s
\le \mathbb E_{\nu_2}p_d(z)^s
\le p_*^s.
}
\tag{R34.5}
```

Moreover `n=O(rL_0)` and `s=ceil(r/T)`.  Therefore the uniform version of
(10.927) is exponent-equivalent to the single-center condition

```math
\boxed{
\max_{z\in\mathcal C_2}
U_m\{S:a_z(S)\le d\}\ge e^{-O(TL_0)}.
}
\tag{R34.6}
```

One center may cost its full atom `2^{-n}`: this is still within the allowed
`e^{-O(rL_0)}` batch probability.  Thus no abundance of good centers and no
minimum over selectors is required.

The soft statement has the same form.  Put

```math
q_U(z)=\mathbb E_{S\sim U_m}e^{-\lambda a_z(S)},
\qquad q_*=\max_{z\in\mathcal C_2}q_U(z),
\qquad L_\lambda(U_m)=\mathbb E_{\nu_2}q_U(z)^s.
```

Then `2^(-n)q_*^s<=L_lambda(U_m)<=q_*^s`.  With
`lambda=Lambda rL_0/(D+1)`, it is enough to prove, for uniform
`0<A<Lambda`,

```math
\boxed{q_*\ge e^{-ATL_0}.}
\tag{R34.7}
```

Indeed the main Laplace term is at least
`e^(-n log 2-A sT L_0)`, while the subtraction in (10.947) is
`e^(-Lambda rL_0)`.  Since `n/(rL_0)=O(1/log k)` and `sT/r->1`, the fixed
gap `Lambda>A` gives the strict tail margin eventually.  Equations
(10.947), (R34.3), and (R34.1) then prove convergence.  Conversely, a soft
moment of size `e^{-O(rL_0)}` forces `q_*>=e^{-O(TL_0)}`.  Thus the soft
high moment and the one-center soft degree are exponent-equivalent.

## 4. A scalar first moment is a stronger uniform-law certificate

Let

```math
q_U(z)=\mathbb E_{S\sim U_m}e^{-\lambda a_z(S)},
\qquad L_\lambda(U_m)=\mathbb E_{\nu_2}q_U(z)^s.
```

The exact uniform soft target is the single high moment `L_lambda(U_m)`,
with no infimum over `w`.  Convexity also gives the stronger scalar
certificate

```math
\boxed{
L_\lambda(U_m)
\ge
\left(\mathbb E_{z\sim\nu_2,S\sim U_m}
e^{-\lambda a_z(S)}\right)^s.
}
\tag{R34.8}
```

Take `lambda=Lambda rL_0/(D+1)`.  For uniform constants `0<A<Lambda`, the
following first-moment statement therefore suffices:

```math
\boxed{
\mathbb E_{z\sim\nu_2,S\sim U_m}e^{-\lambda a_z(S)}
\ge e^{-ATL_0},
\qquad \Lambda r>A sT.
}
\tag{R34.9}
```

Indeed, (R34.8) gives `L_lambda(U_m)>=e^{-AsTL_0}`.  The sharp subtraction
in (10.947) is `e^{-Lambda rL_0}`, so (R34.9) proves (R34.4), then (R34.2),
then convergence.  The common dual weight with a minimum over every selector
in (10.950) is a sufficient certificate for all adversarial laws, but it is
not required for this uniform reduction.  Conversely, (R34.9) is not
necessary for large `L_lambda(U_m)`: rare centers with large `q_U(z)` can
make the high moment large while the first moment remains small, as quantified
by the Renyi decomposition (10.949).

More explicitly, put `Z_U=E_(nu_2)q_U` and let `P_U` have density `q_U/Z_U`
relative to `nu_2`.  Then the exact fixed-law identity

```math
\frac1s\log L_\lambda(U_m)
=\log Z_U+\frac{s-1}{s}D_s(P_U\Vert\nu_2)
\tag{R34.9a}
```

shows the two possible mechanisms.  The exact uniform soft criterion is

```math
-\log Z_U-\frac{s-1}{s}D_s(P_U\Vert\nu_2)\le O(TL_0)
\tag{R34.9b}
```

with a strict constant advantage over the tail subtraction.  Equation
(R34.9) uses only the first term; a rare-center proof may instead use the
positive Renyi correction.

The same statement has a hard-radius version.  With `d=floor(D/s)`, it is
enough that

```math
\mathbb E_{z\sim\nu_2,S\sim U_m}
\mathbf1_{\{a_z(S)\le d\}}\ge e^{-O(TL_0)},
\tag{R34.10}
```

because Jensen applied to the center degrees makes `s` iid selectors share
one center with probability `e^{-O(rL_0)}`.

## 5. Exact remaining wall

This is a logical weakening, not a proof of the stronger first-moment
certificate.  If `N_S=|F_S|`, the same union count as (10.929), averaged only
over uniform `S`, gives

```math
\mathbb E_{S\sim U_m}g_S(d)
\le \frac{4V(m,d)}{2^m}\mathbb E_{S\sim U_m}N_S.
\tag{R34.11}
```

Since `log V(m,d)=o(L_0)` and `TL_0=o(n)`, (R34.10) still requires

```math
\mathbb E_{S\sim U_m}N_S\ge2^{m-o(n)}.
\tag{R34.12}
```

Unlike (10.929), this is an **average** favorable-label multiplicity, not a
uniform lower bound for every selector.  Proving (R34.12), or bypassing the
center star while proving (R34.2), remains minimizer-specific.  The abstract
independent-label wall (10.930) also defeats the uniform version and is still
relevant; it is not realized by an exact signing minimizer.
