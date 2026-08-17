# A fixed-rate near-top tail for every bounded-cap complete signing

**Status.** Task-local theorem draft, frozen for independent audit.  The
proof combines three previously separate ingredients: one-sided block
superadditivity, the archived Bollobas--Scott one-sided energy-product
bound, and the factorized conditional concentration theorem PC.1.  The
conclusion proves the proposed exact-minimizer lemma `L_tail`, but proves no
cross-order recurrence and no convergence theorem.

Throughout,

```math
H_A(x)=\frac12x^TAx=\sum_{i<j}a_{ij}x_ix_j,
\qquad
P(A)=\max_xH_A(x),
\qquad
N(A)=-\min_xH_A(x),
\qquad
Q(A)=\max\{P(A),N(A)\}.
```

All matrices in the main theorem are hollow complete signings.  The
normalization is the half-quadratic, unordered-edge normalization used in
the current project.

## 1. Two deterministic inputs

### Lemma BT.1 (oriented block superadditivity)

For every principal partition `[n]=T\sqcup R`,

```math
\boxed{P(A)\ge P(A[T])+P(A[R]),
\qquad N(A)\ge N(A[T])+N(A[R]).}                 \tag{BT.1}
```

Indeed, choose positive ground states `u,v` of the two diagonal blocks.
The two full spins `(u,v)` and `(-u,v)` have the same internal energy and
opposite cross energies, so one has energy at least the sum of the two
positive block caps.  Apply this to `-A` for the second inequality.

This is the half-energy form of ledger (10.13); it is not new.

### Lemma BT.2 (one-sided cap retained by a linear principal block)

Fix `C>0` and `0<epsilon<1`.  There is `gamma=gamma(C,epsilon)>0`
such that, for all sufficiently large `n`, the following holds.  If
`Q(A)<=Cn^(3/2)` and `R` is a principal set of size
`m>=(1-epsilon)n`, then

```math
\boxed{P(A[R])\ge\gamma n^{3/2}.}                \tag{BT.2}
```

For the constants used below, put

```math
C_0=\max\{C,1\},\qquad \epsilon=\frac14,
\qquad \gamma_C={1\over100000C_0}.              \tag{BT.3}
```

Then (BT.2) holds with `gamma=gamma_C`.

#### Proof and normalization check

Put `B=A[R]`, `p=P(B)`, `q=N(B)`, and `M=Q(A)`.  Principal
monotonicity gives `Q(B)<=M`.  The archived Bollobas--Scott translation is
stated for the doubled energy

```math
P_2(B)=2p,\qquad N_2(B)=2q,\qquad R_2(B)=2(p+q).
```

It gives

```math
P_2(B)R_2(B)
\ge {\{1-r_B^2\}m^3\over1600},
\qquad
r_B={P_2(B)\over m(m-1)},                       \tag{BT.4}
```

provided the associated graph density satisfies the theorem's harmless
condition `p_G(1-p_G)>=1/m`.  Here

```math
0\le r_B\le {2M\over m(m-1)}=O_C(n^{-1/2}),
```

so that condition holds and `1-r_B^2>=1/2` for all sufficiently large
`n`.  Since `p+q<=2Q(B)<=2M`, (BT.4) yields

```math
p\ge{(1-r_B^2)m^3\over12800M}.                  \tag{BT.5}
```

For `m>=3n/4` and `M<=Cn^(3/2)<=C_0n^(3/2)`,

```math
p\ge {27\over1638400C_0}n^{3/2}
   > {1\over100000C_0}n^{3/2}.                  \tag{BT.6}
```

This proves the explicit assertion.  The factor `12800`, rather than
`3200`, is the required conversion from the archived doubled-energy
product theorem to the present half-energy one-sided lower bound.

## 2. Uniform near-top entropy deficit

### Theorem BT.3 (bounded cap forces a fixed-rate two-sided thin tail)

Fix `C>0`.  There are constants `d_C,kappa_C>0` such that every
sufficiently large hollow complete signing satisfying

```math
Q(A)\le Cn^{3/2}                                  \tag{BT.7}
```

obeys, for each `sigma in {+-1}`, on writing

```math
P_\sigma(A)=\max_x\sigma H_A(x),
```

```math
\boxed{
\#\{x:P_\sigma(A)-\sigma H_A(x)<d_Cn^{3/2}\}
\le\exp\{(\log2-\kappa_C)n\}.}                 \tag{BT.8}
```

Consequently, after decreasing `kappa_C` if necessary to absorb a union of
two events,

```math
\#\{x:Q(A)-|H_A(x)|<d_Cn^{3/2}\}
\le\exp\{(\log2-\kappa_C)n\}.                  \tag{BT.8a}
```

One may take

```math
\gamma_C={1\over100000\max\{C,1\}},
\qquad d_C={\gamma_C\over2}.                    \tag{BT.9}
```

The value of `kappa_C` is ineffective only to the extent that the absolute
constants in Grothendieck--Pietsch and Hanson--Wright are left symbolic.

#### Proof

It is enough to prove the positive statement `sigma=+1`; applying it to
`-A` proves the negative statement.  Apply the Grothendieck--Pietsch
construction in PC.1 with fixed
`epsilon=1/4`.  It gives a set `T`, its complement `R`, and, for each
fixed `x_T`, a cross field `h=A_(R,T)x_T` satisfying

```math
|T|<\frac n4,
\qquad
\|A[R]\|_(2\to2)\le32K_GC_0\sqrt n,
\qquad
\|h\|_2\le8\sqrt2K_GC_0n.                     \tag{BT.10}
```

Lemma BT.2 and then Lemma BT.1 give

```math
P(A[R])\ge\gamma_Cn^{3/2},
\qquad
P(A[T])\le P(A)-\gamma_Cn^{3/2}.                \tag{BT.11}
```

Conditionally on `x_T`, with `X_R` uniform,

```math
H_A(x_T,X_R)
=H_(A[T])(x_T)+h^TX_R+H_(A[R])(X_R).             \tag{BT.12}
```

Therefore the event in (BT.8) is contained in

```math
h^TX_R+H_(A[R])(X_R)>{\gamma_C\over2}n^{3/2}.   \tag{BT.13}
```

This is the important uniformity point.  In the notation of PC.1 the
threshold is allowed to vary with the signing,

```math
t_A={P(A)\over n^{3/2}}-d_C,
```

but its proof uses only the fixed gap
`eta=gamma_C-d_C=gamma_C/2`, not the value of `t_A`.  Equation (BT.13)
is the direct event inclusion, so no uniformity in an unrecorded threshold
parameter is being assumed.

Split (BT.13) between its two summands.  The Rademacher linear tail and
(BT.10) give

```math
\Pr\left\{h^TX_R>{\gamma_C\over4}n^{3/2}\right\}
\le\exp\left\{-{\gamma_C^2\over
4096K_G^2C_0^2}n\right\}.                       \tag{BT.14}
```

Since `||A[R]||_F^2<=n^2`, Hanson--Wright gives

```math
\Pr\left\{H_(A[R])(X_R)>{\gamma_C\over4}n^{3/2}\right\}
\le2\exp\left[-c_(HW)n\min\left\{
{\gamma_C^2\over4},{\gamma_C\over64K_GC_0}
\right\}\right].                               \tag{BT.15}
```

Both estimates are uniform in the frozen value of `x_T`.  Average over
`x_T` and put

```math
\kappa_0=\min\left\{
{\gamma_C^2\over4096K_G^2C_0^2},
c_(HW)\min\left({\gamma_C^2\over4},
                 {\gamma_C\over64K_GC_0}\right)
\right\}.                                       \tag{BT.16}
```

The uniform probability is at most `3exp(-kappa_0n)`.  For all sufficiently
large `n` this is at most `exp(-kappa_0n/2)`.  Multiplication by `2^n`
proves the positive half of (BT.8) with `kappa_C=kappa_0/2`.  Repeating for
`-A` proves its negative half.  The absolute near-ground event in (BT.8a)
is contained in the union of those two one-sided events; one more harmless
halving of `kappa_C` absorbs its factor two.  `square`

### Corollary BT.4 (`L_tail` is proved)

Every exact minimizer has `Q(A)=M_n`.  The archived random-sign estimate

```math
M_n\le\sqrt{(\log2)(n^3-n)}<n^{3/2}
```

allows Theorem BT.3 with `C=1`.  Hence fixed constants

```math
d_0={1\over200000},\qquad \kappa>0              \tag{BT.17}
```

satisfy the exact-minimizer upper-tail statement `L_tail` for every
sufficiently large order.  For the application of Theorem 21.8, globally
orient the larger one-sided cap so that `P(A)=Q(A)`.

The conclusion is actually uniform over the full bounded-cap class
`Q(A)<=Cn^(3/2)`; exact coefficient minimality is unnecessary.

## 3. Anatomy of any zero-density positive core

The following consequence addresses an arbitrary core, rather than the
particular fixed-density Pietsch set used above.

### Theorem BT.5 (opposite-polarity core/complement separation)

Let `A_n` be exact minimizers, globally oriented so that
`P(A_n)=M_n`, and let `[n]=T_n\sqcup R_n`, where `k_n=|T_n|=o(n)` and
`m_n=n-k_n`.  Put

```math
\Delta_n=M_n-M_(m_n).
```

Then

```math
0\le\Delta_n
\le M_(k_n)+\sqrt{2(\log2)m_nk_n n}=o(n^{3/2}).  \tag{BT.18}
```

If, for fixed `t>0`,

```math
P(A_n[T_n])\ge(t-o(1))n^{3/2},                  \tag{BT.19}
```

then

```math
\begin{aligned}
N(A_n[T_n])&\le\Delta_n=o(n^{3/2}),\\
N(A_n[R_n])&=Q(A_n[R_n])=M_n-o(n^{3/2}),\\
P(A_n[R_n])&\le M_n-P(A_n[T_n]),\\
N(A_n)&=M_n-o(n^{3/2}).                          \tag{BT.20}
\end{aligned}
```

In particular, the two blocks have opposite macroscopic polarities:

```math
P(A_n[T_n])-N(A_n[T_n])\ge(t-o(1))n^{3/2},
```

```math
N(A_n[R_n])-P(A_n[R_n])\ge(t-o(1))n^{3/2}.      \tag{BT.21}
```

Moreover, the one-sided product theorem forces

```math
P(A_n[R_n])
\ge\left({1\over12800c_n}-o(1)\right)n^{3/2},
\qquad c_n={M_n\over n^{3/2}},                  \tag{BT.22}
```

and hence

```math
P(A_n[T_n])
\le\left(c_n-{1\over12800c_n}+o(1)\right)n^{3/2}.
                                                               \tag{BT.23}
```

Using the rigorous `limsup c_n<=1/2`, the retained complement cap in
(BT.22) is at least `(1/6400-o(1))n^(3/2)`.

#### Proof

Principal monotonicity and the definition of `M_(m_n)` give

```math
M_(m_n)\le Q(A_n[R_n])\le M_n.                  \tag{BT.24}
```

The random-bridge near-order inequality gives (BT.18).  Lemma BT.1 gives

```math
P(A_n[R_n])\le M_n-P(A_n[T_n]).                 \tag{BT.25}
```

For all large `n`, (BT.19) makes the right side strictly below
`M_n-\Delta_n<=Q(A_n[R_n])`.  Thus the absolute cap of the complement is
its negative cap.  Applying the negative half of BT.1 now yields

```math
N(A_n[T_n])
\le M_n-N(A_n[R_n])
\le M_n-M_(m_n)=\Delta_n.                       \tag{BT.26}
```

Equations (BT.20)--(BT.21) follow.  Finally apply the one-sided estimate
in equation (BT.5), with
`m_n/n->1`, `Q(A_n[R_n])<=M_n`, and `r_B=o(1)`, to obtain
(BT.22); combine it with (BT.25) for (BT.23).  `square`

The trivial internal edge count also gives

```math
|T_n|\ge(\sqrt{2t}-o(1))n^{3/4}.                \tag{BT.27}
```

Thus exact minimality does not presently exclude every fixed positive
core.  It forces any surviving core into a sharply oriented two-block
configuration, and Theorem BT.3 excludes precisely the near-ground regime
where the core leaves less than a fixed positive one-sided cap for its
linear-sized complement.  Excluding smaller fixed-level cores would require
a new lemma about this opposite-polarity interface; it is no longer needed
for `L_tail`.

## 4. Consequences and archive comparison

1. **Physical contextual consequence.**  Corollary BT.4 supplies the sole
   tail hypothesis in Theorem 21.8, so exact-minimizer switchings have an
   exponential matched-roof boundary-response code with operator-`O(sqrt n)`
   bridge and `Theta(n^(3/2))` directed gaps.  Combined with the independently
   audited BR.2, this becomes an `Omega(n)`-bit
   all-spins-free scalar physical contextual packing of complete order-`2n`
   signings, all with cap `O(n^(3/2))`.  This is an information-heaviness
   conclusion, not convergence.
2. **Archive collisions.**  BT.1 is ledger (10.13).  BT.2 is the
   half-energy specialization of ledger (10.148)--(10.151), and the
   near-order estimate in BT.5 is Theorem 36.15.  PC.1--PC.3 already proved
   that a zero-rate fixed upper tail creates a zero-density principal core.
   The qualitative opposite-polarity conclusion in BT.5 was independently
   derived as OC.2 in `exact_minimizer_oriented_core_separation.md`; only
   the explicit retained-complement estimates (BT.22)--(BT.23) sharpen that
   parallel synthesis.
   The new theorem-level increment is the synthesis: apply the one-sided
   product theorem to the *linear complement of the Pietsch heavy set* and
   feed its fixed positive cap back through block superadditivity before
   conditional concentration.  No searched archive statement made this
   feedback step or concluded the uniform bounded-cap thin-tail theorem.
3. **Why the old core formulation looked open.**  PC.3 replaced the moving
   near-ground level by a fixed lower threshold `t_0`.  At that lower level,
   an opposite-polarity core is not excluded by current identities.  The
   proof above retains the actual moving level `P_sigma(A)-d_Cn^(3/2)`.  The
   complement's universal positive cap then supplies exactly the missing
   fixed margin.
4. **No convergence claim.**  BT.3 is a same-order entropy theorem.  It
   proves neither a cross-order composition inequality nor a recovery
   theorem and does not alter the rigorous interval for `M_n/n^(3/2)`.
