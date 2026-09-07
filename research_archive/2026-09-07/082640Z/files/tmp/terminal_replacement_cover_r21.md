# Wave 21 main audit: terminal replacement is a max-plus translate cover

## Status

The identities below are exact in the ordered `Q` normalization.  They turn
global minimality under every switched optimal block replacement into a
translate-cover condition.  The order-nine hybrid wall shows that the cover
can hold with room to spare while the terminal restriction has positive
excess.  Thus the cover is a useful interface, but it does not by itself
control terminal refresh or prove convergence.

## Exact response profile

Let `A` be an exact order-`n` minimizer, let `T` have order `m`, and put
`B=A[T]`.  Let `D_T` be the group of oriented projective cuts on `T` under
entrywise multiplication.  It has order `2^m`.  For `d in D_T`, define the
best external response

```math
Z(d)=\max\{\langle A,d'\rangle-\langle B,d\rangle:
             d'\in\mathcal D_n,\ d'[T]=d\}.
```

Then

```math
Q(A)=\max_d\{Z(d)+\langle B,d\rangle\}.
```

Fix any exact order-`m` minimizer `C`.  For every `h in D_T`, replacing `B`
by the switched/complemented signing `C odot h` gives another order-`n`
signing.  Exact global minimality therefore implies

```math
\max_d\{Z(d)+\langle C,hd\rangle\}\ge Q(A)
\qquad(h\in\mathcal D_T).
```

Put

```math
q=Q(A),\quad q_B=Q(B),\quad q_m=Q(C),\quad
\varepsilon=q_B-q_m,
```

```math
s(d)=Z(d)-(q-q_B),\qquad
b_B(d)=q_B-\langle B,d\rangle,\qquad
b_C(d)=q_m-\langle C,d\rangle.
```

The parent cap and its attainment say

```math
s(d)\le b_B(d),\qquad \max_d[s(d)-b_B(d)]=0.
```

All switched optimal replacements give the exact covering condition

```math
\boxed{
\forall h\in\mathcal D_T\quad
\max_d\{s(d)-b_C(hd)\}\ge\varepsilon.
}
```

Thus every translate of the optimal terminal deficit profile must meet enough
external surplus to pay the full terminal excess.

## Layer and pressure consequences

For every realized deficit level `t` of `C`, put

```math
N_t=\{u:b_C(u)\le t\},\qquad
R_t=\{d:s(d)\ge\varepsilon+t\}.
```

Choosing a maximizing witness in the boxed condition gives

```math
\boxed{
\mathcal D_T=\bigcup_t N_tR_t^{-1},
\qquad
2^m\le\sum_t|N_t||R_t|.
}
```

There is also an exact soft consequence.  Write

```math
D_\lambda(C)=\sum_u e^{-\lambda b_C(u)}.
```

Exponentiate the boxed maximum, sum over `d`, and average over the translate
`h`.  Translation invariance gives

```math
\boxed{
e^{\lambda\varepsilon}
\le
\left(2^{-m}\sum_de^{\lambda s(d)}\right)D_\lambda(C)
\qquad(\lambda>0).
}
```

This is a genuine global-minimality constraint, but a useful excess bound
still requires an upper tail estimate for the external surplus.  The parent
cap only gives `s<=b_B`, whose exponential moment can have leading scale.

## Exact order-nine wall

Take the displayed `A_9`, let `T={3,4,5,6,7,8}`, and use the exact order-six
minimizer `C=A_6` from the verifier.  Then

```math
q=24,\qquad q_B=14,\qquad q_6=10,\qquad\varepsilon=4.
```

Exhaustive enumeration gives

```math
\boxed{s(d)\le b_B(d),\qquad
       \max_d[s(d)-b_B(d)]=0.}
```

Equality holds at 17 of the 64 terminal cuts.  The parent deficit
`b_B(d)-s(d)` has exact histogram
`0^{17},8^{20},16^{18},24^8,32^1`; the external surplus itself has range
`[-4,20]`.  Thus the response does not refund every terminal deficit, but it
has a substantial exact-ground contact layer.  Moreover

```math
\min_h\max_d[s(d)-b_C(hd)]=8,
```

so every switched/complemented optimal replacement has full norm at least
`24-4+8=28`; this recovers the exact hybrid wall (10.551).  The layer sizes
`(|N_t|,|R_t|)` at `t=0,4,16,20` are respectively

```text
(12,35), (32,17), (52,1), (64,0).
```

Hence access to every optimal switched replacement does not force descent;
the cover can be strictly stronger than its required excess at every
translate.  A positive theorem must bound the external surplus profile using
additional minimizer structure.  The translate cover alone is an exact
restatement of the compatibility gap already visible in the hybrid wall.
