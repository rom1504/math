# Exact deletion--extension response and the canonical erased-row mismatch

Status: **rigorous actual-minimizer theorem and exact all-actual-child
falsifier**.  Every row of an exact pressure-minimizing signing minimizes a
one-vertex external-field response of its own deletion.  This gives an exact
dynamic recurrence for the minimized pressure and identifies a narrow
inverse-escort mode statement.  The canonical erased-row law for a general
two-child bridge is a different, sector-biased extension response.  An
order-three/order-two example, in which every signing used is an exact
thermal minimizer, shows that the actual deleted row need not even be a mode
of that canonical law.

No conference or Paley object is used.

## 1. Sector partitions and one-vertex response

For a signing `C` of order `d`, define

```math
Z_C^\sigma(t)=2^{-d}\sum_{x\in\{\pm1\}^d}
 e^{\sigma tH_C(x)},
\qquad \sigma\in\{\pm1\},                           \tag{DER.1}
```

and

```math
\overline Z_C(t)={Z_C^+(t)+Z_C^-(t)\over2}
=2^{-d}\sum_x\cosh(tH_C(x)).                        \tag{DER.2}
```

Let `mu_(C,sigma,t)` be the sector Gibbs law proportional to
`exp(sigma t H_C(x))`.  For a proposed new row
`b in {+-1}^d`, internal temperature `t`, and star amplitude `u`, put

```math
R_C^\sigma(b;t,u)
=E_{X\sim\mu_{C,\sigma,t}}
  \cosh(u\langle b,X\rangle).                       \tag{DER.3}
```

The augmented sector weights and response are

```math
\omega_C^\sigma(t)
={Z_C^\sigma(t)\over Z_C^+(t)+Z_C^-(t)},
\qquad
R_C^{\rm aug}(b;t,u)
=\sum_{\sigma=\pm1}\omega_C^\sigma(t)
 R_C^\sigma(b;t,u).                                 \tag{DER.4}
```

Write `C oplus b` for the order-`d+1` signing obtained by adjoining one
vertex with incident row `b`.

**Lemma DER.1 (exact extension identity).**  At the same star and internal
temperature,

```math
\boxed{
\overline Z_{C\oplus b}(t)
=\overline Z_C(t)R_C^{\rm aug}(b;t,t).}              \tag{DER.5}
```

More generally, if the old edges have amplitude `t` and the new star has
amplitude `u`, the right side is
`bar Z_C(t) R_C^aug(b;t,u)`.

*Proof.*  If `z` is the new spin, then sectorwise

```math
\begin{aligned}
Z_{C\oplus b}^\sigma(t)
 &=2^{-(d+1)}\sum_{x,z}
 e^{\sigma t\{H_C(x)+z\langle b,x\rangle\}}\\
 &=2^{-d}\sum_xe^{\sigma tH_C(x)}
       \cosh(t\langle b,x\rangle)\\
 &=Z_C^\sigma(t)R_C^\sigma(b;t,t).
                                                               \tag{DER.6}
\end{aligned}
```

Average the two sectors.  Replacing the star occurrence of `t` by `u`
proves the anisotropic statement. `square`

## 2. Every row of an exact minimizer is an optimal reinsertion

Let

```math
F_n(t)=\min_{A\in\mathcal S_n}\log\overline Z_A(t), \tag{DER.7}
```

where `S_n` is the set of hollow signings of order `n`.  For an order-`n`
signing `A`, delete vertex `i`, write `C_i=A[V\setminus\{i\}]`, and let
`b_i=(a_(ij))_(j ne i)` be its deleted row.

**Theorem DER.2 (actual-row cavity optimality).**  If `A` attains `F_n(t)`,
then for every vertex `i`,

```math
\boxed{
b_i\in\mathop{\rm argmin}_{b\in\{\pm1\}^{n-1}}
 R_{C_i}^{\rm aug}(b;t,t).}                         \tag{DER.8}
```

Moreover `C_i` attains the combined dynamic objective:

```math
F_n(t)
=\log\overline Z_{C_i}(t)
 +\min_b\log R_{C_i}^{\rm aug}(b;t,t).              \tag{DER.9}
```

*Proof.*  Replacing only row `i` by any `b` produces another order-`n`
signing.  Minimality of `A`, followed by DER.1, proves (DER.8).  The left
side of (DER.9) is attained by `(C_i,b_i)`.  No other pair `(C,b)` can have
a smaller value, since its extension would contradict the definition of
`F_n`. `square`

The star-flip form is

```math
E_{\nu_A}\exp\left{
-2t\tau x_i\sum_{j\in S}a_{ij}x_j\right}\ge1
\quad(S\subseteq V\setminus\{i\}).                 \tag{DER.10}
```

Thus DER.2 is a useful cavity interpretation of actual optimality, but
(DER.10) is precisely the star-subset specialization of AC.32; it is not an
additional independent constraint.

## 3. Exact dynamic recurrence

Define the deletion-state value

```math
V_C(t)=\min_b\log R_C^{\rm aug}(b;t,t).             \tag{DER.11}
```

Every order-`n` signing has a unique representation after declaring its
last vertex as `C oplus b`.  Therefore

```math
\boxed{
F_n(t)=\min_{C\in\mathcal S_{n-1}}
 \{\log\overline Z_C(t)+V_C(t)\}.}                 \tag{DER.12}
```

This is the exact one-step Bellman recurrence.  Importantly, its optimizing
`C` need not attain `F_(n-1)(t)`: a higher-pressure deletion can compensate
with a better extension response.  Replacing (DER.12) by a recurrence only
over order-`n-1` pressure minimizers is unjustified.

Two sharp elementary bounds follow.  Since `cosh>=1`, `V_C>=0`.  Also

```math
E_{b\sim U_d}R_C^{\rm aug}(b;t,t)
=(\cosh t)^d,                                       \tag{DER.13}
```

because the uniform row average of
`cosh(t sum_j b_jx_j)` is `(cosh t)^d` for every `x`.  Hence

```math
\boxed{
F_{n-1}(t)\le F_n(t)
\le F_{n-1}(t)+(n-1)\log\cosh t.}                  \tag{DER.14}
```

For every deletion of an exact order-`n` minimizer,

```math
0\le V_{C_i}(t)
=F_n(t)-\log\overline Z_{C_i}(t)
\le F_n(t)-F_{n-1}(t).                              \tag{DER.15}
```

At `t=beta/sqrt(N)`, DER.14 is an `O_beta(1)` one-vertex increment.  It is
exact but does not by itself control a macroscopic deletion chain, because
the Bellman state includes the complete function `b -> R_C^aug(b;t,t)`.

## 4. The narrow inverse-escort consequence

Normalize the same-temperature augmented response by its uniform average:

```math
z_C^{\rm aug}(b;t)
={R_C^{\rm aug}(b;t,t)\over(\cosh t)^d},
\qquad E_{U_d}z_C^{\rm aug}=1.                      \tag{DER.16}
```

For `lambda>0`, define its inverse escort

```math
{dr_{C,\lambda}^{\rm aug}\over dU_d}(b)
={z_C^{\rm aug}(b;t)^{-\lambda}
  \over E_{U_d}z_C^{\rm aug}(b;t)^{-\lambda}}.     \tag{DER.17}
```

DER.2 proves the exact support statement

```math
\boxed{
\text{every deleted row }b_i\text{ is a mode of }
r_{C_i,\lambda}^{\rm aug}.}                         \tag{DER.18}
```

It also gives `z_(C_i)^aug(b_i;t)<=1`.  This is the strongest universal
mass statement available from minimization alone: the escort has full
support, and a mode can still have only its trivial uniform mass
`2^{-(n-1)}`.  No basin-size or effective-support lower bound follows from
DER.18.

## 5. What the canonical erased-row law actually escorts

Now take the actual two-child setup of CR.0.  Let the left child be `A`, the
right child be `D` of order `d`, fix relative orientation `epsilon`, and keep
the internal child amplitude `t`.  Reindex the right sector by
`a=epsilon s`.  The exact weights in the erased-row channel are

```math
\widehat\omega_a^{A,D,\epsilon}(t)
={Z_A^{\epsilon a}(t)Z_D^a(t)
  \over\sum_{c=\pm1}Z_A^{\epsilon c}(t)Z_D^c(t)}.   \tag{DER.19}
```

Consequently Proposition CR.0 is exactly

```math
\boxed{
z_{A\to D}^{\epsilon}(b;t,u)
={1\over(\cosh u)^d}
 \sum_{a=\pm1}\widehat\omega_a^{A,D,\epsilon}(t)
 R_D^a(b;t,u).}                                    \tag{DER.20}
```

Indeed each sector law is invariant under `Y -> -Y`, so
`E exp(u<b,Y>)=E cosh(u<b,Y>)`; the binary-channel identity
`prod_j(1+tanh(u)b_jY_j)=exp(u<b,Y>)/(cosh u)^d` supplies exactly the
displayed normalization.

The canonical erased-row law is

```math
{dr_{\rm row,u}\over dU_d}(b)
\propto z_{A\to D}^{\epsilon}(b;t,u)^{-\lambda}.   \tag{DER.21}
```

Thus it is indeed an exact inverse escort of a one-vertex **anisotropic,
sector-biased** extension response.  It equals the normalized augmented
extension response in DER.16 when

1. `u=t` (the new star has the same amplitude as the old edges);
2. the left child is sector-neutral,
   `Z_A^+(t)=Z_A^-(t)`, so that
   `widehat omega_a=omega_D^a`.

To transfer the actual-row mode conclusion DER.18, one additionally needs
`D` to be the particular deletion of a minimizing parent whose deleted row
is being tested.

In the macroscopic two-child bridge, none of these alignments is automatic.
Even if the left weights are neutral, minimality of a right child `D` says
that its existing rows optimally reinsert into `D\setminus i`; it says
nothing about a new row extending the full `D` to order `d+1`.

## 6. Exact all-actual-child falsifier to sector-weight transfer

The sector-bias mismatch occurs at the first nontrivial orders.  Let `T` be
the all-positive triangle and `E` the positive single-edge signing.

Every order-three signing is an exact augmented-pressure minimizer for every
`t>0`: according to its triangle product, its projective energy multiset is
either

```math
\{3,-1,-1,-1\}\quad\text{or}\quad\{-3,1,1,1\},
```

and the absolute-energy multisets coincide.  The order-two signing `E` is
also an exact minimizer.  Deleting one vertex of `T` leaves `E` and the
actual deleted row

```math
b_+=(1,1).                                           \tag{DER.22}
```

For `E`, one has `Z_E^+=Z_E^-=cosh t`.  Put
`b_times=(1,-1)`.  Direct sector summation gives

```math
\begin{aligned}
R_E^+(b_+;t,t)
 &={e^t\cosh(2t)+e^{-t}\over2\cosh t},\\
R_E^-(b_+;t,t)
 &={e^{-t}\cosh(2t)+e^t\over2\cosh t},              \tag{DER.23}\\
R_E^+(b_\times;t,t)&=R_E^-(b_+;t,t),\\
R_E^-(b_\times;t,t)&=R_E^+(b_+;t,t).
\end{aligned}
```

Their augmented weights are equal, so `b_+` and `b_times` tie, as DER.2
requires.  But use a second copy of the exact minimizing triangle `T` as the
left child in the canonical erased-row law, take `epsilon=+1`, and retain
`E` as the right child.  Its sector partitions are

```math
Z_T^+={e^{3t}+3e^{-t}\over4},
\qquad
Z_T^-={e^{-3t}+3e^t\over4},
\qquad Z_T^+-Z_T^-=2\sinh^3t>0.                    \tag{DER.24}
```

Hence the canonical right-sector weight

```math
w={Z_T^+\over Z_T^++Z_T^-}
```

is strictly larger than `1/2`.  Moreover

```math
R_E^+(b_+)-R_E^-(b_+)
={(e^t-e^{-t})(\cosh2t-1)\over2\cosh t}>0.          \tag{DER.25}
```

Combining DER.20 and DER.23--DER.25 at `u=t` gives

```math
\boxed{
z_{T\to E}^{+}(b_+;t,t)
>z_{T\to E}^{+}(b_\times;t,t).}                    \tag{DER.26}
```

Therefore the actual deleted row `b_+`, although a minimizer of the exact
augmented reinsertion response, has **strictly smaller** inverse-escort
density than `b_times` in the canonical erased-row law.  All three roles
here—the minimizing parent `T`, the left child `T`, and the right child
`E`—are occupied by exact pressure minimizers at the common raw temperature.

Changing orientation to `epsilon=-1` reverses the left-sector bias and
reverses which row is preferred; it does not restore a sectorwise theorem.

## 7. Consequence for the actual-child campaign

Deletion--reinsertion gives a real optimizer-specific statement, but its
state does not coincide with the canonical erased-row state:

```text
actual row optimality
  = neutral-sector, same-temperature response of its own deletion;

canonical erased row
  = other-child-biased, possibly anisotropic response extending the full
    right child.
```

DER.26 is a precise no-go to transferring the mode property across that
gap.  Closing it would require a new theorem of one of the following forms:

- **sector synchronization:** the two functions `R_D^+` and `R_D^-` have
  asymptotically the same low-response set at physical scale;
- **orientation compensation:** the two biased escorts can be combined
  before paying an absolute value or KL cost;
- **adjacent-order stability:** a pressure minimizer `D` is also a near-best
  Bellman deletion state for order `d+1`, uniformly over the relevant sector
  weights.

None follows from exact pressure minimization, and the finite falsifier
already kills an exact version of sector synchronization.  The Bellman
recurrence DER.12 is exact but retains the complete extension-response
table; without one of these new inputs it is a reformulation, not a strict
compression of the actual child law.
