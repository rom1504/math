# Wave 15 route 3: the exact bicriteria problem and its scalar wall

Status: all algebraic claims and the finite order-nine census below are
verified.  This memo makes no tracked-file edit.

## 1. The nondegenerate Pareto problem

Let `A` be an order-`n` global minimizer, so `Q(A)=q_n`.  For a partition
`P=(V_i)` with `m_i=|V_i|`, write

```math
A=C+\bigoplus_iD_i,
\qquad D_i=A[V_i],
```

and define

```math
S(P)=\sum_iQ(D_i),\qquad
B(P)=\sum_iq_{m_i},\qquad
X(P)=S(P)-B(P),
```

```math
d(P)=Q(C)-q_n,
\qquad O(P)=[d(P)]_+.
```

Here `X` is the captured internal optimality excess and `O` is the cross-only
overshoot.  For a fixed, genuinely nondegenerate size profile `m=(m_i)`, the
pure Pareto envelope is

```math
\boxed{
X^*_{A,m}(\eta)
=\max\{X(P): P\text{ has profile }m,\ O(P)\le\eta\}.
}
\tag{R15.1}
```

Equivalently one can minimize `O-lambda X` for `lambda>=0`.  Since the pure
frontier need not be convex, the exact randomized relaxation is worth stating.
For distributions `mu` on partitions of profile `m`,

```math
\boxed{
X^{\rm mix}_{A,m}(\eta)
=\max_{\mu:\ \mathbb E_\mu O\le\eta}\mathbb E_\mu X
=\inf_{\lambda\ge0}
\left\{\lambda\eta+
\max_{P}\bigl(X(P)-\lambda O(P)\bigr)\right\}.
}
\tag{R15.2}
```

This is finite-dimensional LP duality.  Fixing a nondegenerate profile is
essential.  If arbitrary profiles are allowed, both the all-singleton
partition and the one-block partition have `(O,X)=(0,0)`; cross control alone
then says nothing.

## 2. Exact insertion into the common-mosaic identity

Choose arbitrary local minimizers `G_i`, put

```math
K_G=\bigoplus_iG_i,
\qquad H_G=C+K_G,
\qquad E_G=Q(H_G)-q_n\ge0.
```

Because `Q(K_G)<=B(P)`, triangle and reverse triangle give

```math
|Q(H_G)-Q(C)|\le B(P).
```

Consequently every tuple, not merely the best tuple, satisfies

```math
\boxed{
\max\{0,d(P)-B(P)\}
\le E_G\le d(P)+B(P).
}
\tag{R15.3}
```

The right endpoint is nonnegative because global minimality gives
`q_n<=Q(H_G)<=Q(C)+B(P)`.  In terms of positive overshoot only, the slightly
weaker but symmetric statement is

```math
\boxed{
[O(P)-B(P)]_+\le E_G\le O(P)+B(P).
}
\tag{R15.4}
```

The responses in (10.524) have the elementary closed forms

```math
\Phi_D=q_n-S(P),
\qquad
\Phi_G=Q(H_G)-B(P).
```

Thus (10.524) becomes

```math
\boxed{
\Phi_G-\Phi_D=X(P)+E_G.
}
\tag{R15.5}
```

Combining (R15.4) and (R15.5),

```math
\boxed{
X(P)+[O(P)-B(P)]_+
\le\Phi_G-\Phi_D
\le X(P)+O(P)+B(P).
}
\tag{R15.6}
```

For the selector-adjusted Pareto value

```math
R^*_{A,m}(\eta)
=\max_{P:\ O(P)\le\eta}\min_G(\Phi_G-\Phi_D),
```

the fixed-profile constant `B_m=sum_i q_{m_i}` therefore gives

```math
\boxed{
X^*_{A,m}(\eta)
\le R^*_{A,m}(\eta)
\le X^*_{A,m}(\eta)+\eta+B_m.
}
\tag{R15.7}
```

This is the exact useful meaning of the bicriterion.  When `B_m` and `eta`
are subleading, maximizing captured excess is asymptotically the same as
maximizing the forced increase of the common-mosaic response.  The excess is
not a credit that lowers the hybrid norm: it enters (R15.5) with a plus sign.

## 3. What scalar comparison actually follows

The only direct scalar statement furnished by global minimality and the
triangle bound is

```math
\boxed{
[q_n-Q(C)]_+\le B(P)=\sum_iq_{m_i}.
}
\tag{R15.8}
```

It is a lower bound on the sum of child minima, not the upper comparison
needed for the directed convergence argument.  Notice that `X` cancels
completely.  Knowing only `O=[Q(C)-q_n]_+` is weaker still: when `O=0`, it
does not retain the signed cross deficit `q_n-Q(C)`.

To display the desired recurrence exactly, set

```math
a_r=\frac{q_r}{r^{3/2}},
\qquad
\beta_m=\sum_i\left(\frac{m_i}{n}\right)^{3/2}.
```

Then

```math
\boxed{
B(P)-\beta_m q_n
=\sum_i m_i^{3/2}(a_{m_i}-a_n).
}
\tag{R15.9}
```

For `k=n/m` equal blocks this is

```math
\boxed{
B(P)-\sqrt{\frac mn}\,q_n
=k\left[q_m-\left(\frac mn\right)^{3/2}q_n\right].
}
\tag{R15.10}
```

Since `B=S-X`, a directed comparison

```math
q_m\le\left(\frac mn\right)^{3/2}q_n+E_{n,m}
```

is equivalent to

```math
\boxed{
X(P)\ge
S(P)-\sqrt{\frac mn}\,q_n-kE_{n,m}.
}
\tag{R15.11}
```

But for a fixed equal-block profile, `S(P)-X(P)=kq_m` for every partition.
Therefore (R15.11) is independent of the partition and is exactly the desired
recurrence rewritten.  An existence theorem phrased only as (R15.11) is
circular.  Cross overshoot `O` does not occur in (R15.9)--(R15.11), and
(R15.8) cannot supply the missing upper bound on `B`.

Accordingly, a noncircular positive result needs an additional input beyond
the `(O,X)` Pareto pair.  Examples would be an independently proved upper
bound `S(P)<=T(P)` together with a structural lower bound
`X(P)>=T(P)-beta_m q_n-error`, or a later transport theorem that converts the
large response increase in (R15.7) into a separately signed scalar credit.
The common-mosaic identity alone provides neither.

## 4. Exact order-nine wall

Use `A_9` from (10.550) and the zero-based partition

```math
V_1=\{0,1,2,4,5,6\},
\qquad V_2=\{3\},\quad V_3=\{7\},\quad V_4=\{8\}.
```

Exact Boolean enumeration gives

```math
\boxed{
q_9=24,quad Q(C)=22,quad O=0,quad
Q(D_1)=22,quad B=q_6=10,quad X=12.
}
\tag{R15.12}
```

Enumeration of all 384 labelled order-six minimizers gives

```math
\boxed{
\min_{G_1\in\mathcal M_6}Q(C\oplus G_1)=24,
}
\tag{R15.13}
```

with 40 minimizers attaining 24.  Hence the best hybrid has `E_G=0`, while
the common-mosaic response rises by exactly `X=12`.  This is a sharp finite
demonstration of the sign in (R15.5).  In particular, the tempting
"captured-excess credit" inequality is false:

```math
q_9=24
>Q(C)+B-X=22+10-12=20.
\tag{R15.14}
```

The checker enumerates all 21,147 set partitions of `A_9`.  The global Pareto
frontier has `O=0` and maximal `X=12`, realized by several profiles including
`6+1+1+1`, `3+6`, and `4+5`.  In contrast, the all-singleton and one-block
partitions both give `(O,X)=(0,0)`, and every `3+3+3` partition with minimal
overshoot also has `X=0`.  Thus perfect cross control can coexist with zero,
moderate, or large captured excess even inside one fixed global minimizer.

Run the audit with

```bash
.venv/bin/python tmp/check_bicriteria_recurrence_r15.py
```

## 5. Conclusion

The Pareto/minimax problem is meaningful for finding common-mosaic response,
and (R15.7) is its exact reduction.  It does not itself yield a recurrence for
`q_n`.  At a fixed size profile, the desired scalar upper comparison is
already encoded in the partition-invariant number `B=sum_i q_{m_i}`; captured
excess enters the response with the wrong sign for direct cancellation.  A
future partition theorem must therefore add an independent estimate tying
`S-X=B` to the normalized parent scale, or connect the response increase to a
new temporal/scalar mechanism.  Merely requiring small cross overshoot and
large internal excess, even with an exactly minimizing hybrid, is insufficient.
