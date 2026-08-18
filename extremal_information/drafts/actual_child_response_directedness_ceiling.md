# Why extensive leave-one-out response is not yet directed row dependence

Status: **rigorous scalable method ceiling**.  Corollary RR.3 gives an exact
within-row/cross-row split for the actual inverse escort.  The construction
below shows that even extensive cross-row erasure information, together with
the exact weak-coordinate and tight-row-Renyi scales of the actual law, does
not force an extensive reverse row-product projection.  Any Level-6 use of
the response identity must retain its KL direction or add a functional
inequality specific to the optimized-child channel.

This is a falsifier of an inference, not a child-induced counterexample.

## The weak common-latent construction

Fix `epsilon in (0,1)` and a sufficiently small constant `c>0`.  Let there
be `r` rows of `r` fair-sign coordinates, put `d=r^2` and
`a_r=c/sqrt(r)`, and introduce a latent variable

```math
\Pr(W=0)=1-\epsilon,
\qquad \Pr(W=1)=\Pr(W=-1)=\epsilon/2.                 \tag{DC.1}
```

Conditionally on `W`, all `d` bits are independent with mean `a_rW`.  Denote
their joint law by

```math
q_r=(1-\epsilon)U_d
 +{\epsilon\over2}\nu_{+,r}^{\otimes d}
 +{\epsilon\over2}\nu_{-,r}^{\otimes d}.             \tag{DC.2}
```

It is centrally symmetric and has full support.

**Proposition DC.1 (directionality ceiling).**  The laws (DC.2) have all of
the following properties.

1. Every full one-bit conditional mean has magnitude at most
   `a_r=O(r^(-1/2))`.
2. After arbitrary conditioning outside a row and arbitrary marginalization
   inside it, every retained row component obeys

   ```math
   D_2(q_{r,S}\Vert U_S)\le |S|\log(1+a_r^2)\le c^2. \tag{DC.3}
   ```

3. The reverse row-product projection is tight:

   ```math
   \inf_{p=\otimes_i p_i}D(p\Vert q_r)
   \le D(U_d\Vert q_r)\le-\log(1-\epsilon)=O(1).      \tag{DC.4}
   ```

4. Nevertheless the irreducible cross-row erasure information from
   (RR.23a) is extensive:

   ```math
   \boxed{
   \sum_{i,j}I(B_{ij};R_{-i}\mid R_{i,-j})
   \ge\kappa_{c,\epsilon}r}                           \tag{DC.5}
   ```

   for all sufficiently large `r`, with
   `kappa_(c,epsilon)>0`.

*Proof.*  Given all other bits, the posterior is some law on `W`; hence the
conditional mean of the omitted bit is `a_r E[W|rest]` and has magnitude at
most `a_r`.

Given anything outside a retained set `S` in one row, its conditional law is
a mixture of the fair product law and the two product laws with means
`+-a_r`.  If `f_0,f_+,f_-` are their densities relative to `U_S`, convexity
of the squared `L^2` norm gives

```math
\left\|w_0f_0+w_+f_++w_-f_-\right\|_2^2
\le w_0+w_+(1+a_r^2)^{|S|}+w_-(1+a_r^2)^{|S|}
\le e^{c^2}.                                          \tag{DC.6}
```

This proves (DC.3).  The pointwise domination
`q_r>=(1-epsilon)U_d` proves (DC.4).

It remains to prove (DC.5).  Fix a bit `e=(i,j)`, let `A=R_(i,-j)` be the
other bits in its row, and let `C=R_(-i)` be all other rows.  Put

```math
M_A=\mathbb E[W\mid A],
\qquad M_{AC}=\mathbb E[W\mid A,C].                   \tag{DC.7}
```

The two conditional laws of `B_e` have means `a_rM_A` and `a_rM_(AC)`.
Pinsker's inequality for a fair-sign bit and the tower property therefore
give

```math
I(B_e;C\mid A)
\ge {a_r^2\over2}\mathbb E(M_{AC}-M_A)^2.             \tag{DC.8}
```

The `r-1` observations in `A` carry only constant signal-to-noise.  More
quantitatively,

```math
I(W;A)
\le\epsilon(r-1)D(\nu_{+,r}\Vert U_1)
\le C\epsilon c^2.                                    \tag{DC.9}
```

Choosing `c` as a sufficiently small fixed constant makes the joint law of
`(W,A)` a fixed positive distance from perfect recovery.  Quantitatively,
write `M_A=E[W|A]`.  Pinsker applied to the joint law and the product of its
marginals gives

```math
\mathbb E M_A^2
=\mathbb E[WM_A]
\le 2\|P_{W,A}-P_WP_A\|_{\rm TV}
\le\sqrt{2I(W;A)}.                                   \tag{DC.9a}
```

Since `E W^2=epsilon`, choose the fixed `c` so that the last term is at most
`epsilon/2`.  Hence, uniformly in `r`,

```math
\mathbb E\operatorname{Var}(W\mid A)
=\epsilon-\mathbb E M_A^2
\ge\delta_{c,\epsilon}:={\epsilon\over2}>0.           \tag{DC.10}
```

In contrast, `C` contains `r(r-1)` observations.  Their empirical sum has
means `0,+-a_rr(r-1)` in the three latent sectors and standard deviation at
most `r`.  The estimator which thresholds this sum at half the nonzero mean
has error at most

```math
2\exp\{-a_r^2r(r-1)/8\}=e^{-\Omega_{c}(r)}.           \tag{DC.10a}
```

Conditional expectation minimizes squared error, so
`E Var(W|A,C)=o(1)`.  Martingale orthogonality now gives the exact identity

```math
\mathbb E(M_{AC}-M_A)^2
=\mathbb E\operatorname{Var}(W\mid A)
 -\mathbb E\operatorname{Var}(W\mid A,C)
\ge\delta_{c,\epsilon}-o(1).                         \tag{DC.11}
```

All `r^2` bits are symmetric.  Summing (DC.8) and using
`a_r^2=c^2/r` proves (DC.5). `square`

## Consequence for the actual-child campaign

The counterexample simultaneously has:

- the actual law's `O(r^(-1/2))` one-bit conditional scale;
- uniformly tight conditional row `D_2`;
- `Omega(r)` cross-row leave-one-out dependence;
- only `O(1)` reverse distance to a row product.

Thus neither the response norm in SH.0 nor the cross-row term in RR.23a is,
by itself, the directed resource in (AC.16).  The useful next theorem must be
one of:

1. a child-specific exclusion of a constant-weight product background such
   as (DC.2);
2. a row functional inequality comparing cross erasure information to the
   reverse projection under an additional actual-channel hypothesis; or
3. a coarse row feature whose **reverse** product projection is certified
   directly by data processing.

The third alternative is implemented finitely by GC.2; its scalable form is
now the narrowest evidence-backed target.
