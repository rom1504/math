# The sharp interface scale for quadratic orientation visibility

**Status.** Rigorous task-local theorem.  This sharpens the positive boundary
in BH.3.  The sign of an even old landscape is hidden by the outer absolute
value until the continuation has enough **internal cap**.  For unrestricted
exact-sign quadratic continuations the first possible width is
`n^(3/4)`, and a biased flat bridge attains that exponent while the complete
parent remains at the `n^(3/2)` scale.

## 1. The continuation-cap law

Let `H` be any even landscape on `{+-1}^n`, let
`B in {+-1}^{n times m}`, and let `K` be any even landscape on
`{+-1}^m`.  Define

```math
R_\sigma(B,K)
=\max_{x,y}|\sigma H(x)+x^TBy+K(y)|,
\qquad \sigma\in\{+-1\}.                         \tag{OV.1}
```

### Theorem OV.1 (orientation visibility is paid by internal future cap)

For every `H,B,K`,

```math
|R_+(B,K)-R_-(B,K)|\le2||K||_\infty.             \tag{OV.2}
```

In particular, if `K=H_C` is a quadratic signing on `m` spins, then

```math
|R_+-R_-|\le2Q(C)\le m(m-1).                     \tag{OV.3}
```

More generally, on a continuation class with `Q(C)<=K_0m^alpha`, a response
gap at least `epsilon n^(3/2)` requires

```math
m\ge\left({\epsilon\over2K_0}\right)^{1/\alpha}
        n^{3/(2\alpha)}.                          \tag{OV.4}
```

Thus arbitrary complete sign continuations have the width barrier
`Omega(n^(3/4))`.  If the continuation block is itself bounded-cap,
`Q(C)=O(m^(3/2))`, a target-scale orientation gap requires `m=Omega(n)`.

#### Proof

First omit `K`.  Inverting all new spins and then negating the expression
gives

```math
\max_{x,y}|-H(x)+x^TBy|
=\max_{x,y}|H(x)+x^TBy|.                          \tag{OV.5}
```

Adding `K` changes either cap by at most `||K||_infinity`; the triangle
inequality proves (OV.2).  Equations (OV.3)--(OV.4) follow immediately.
`square`

This result applies to every even old landscape; no spectral or Walsh
hypothesis enters the lower bound.  It identifies the actual budget as
internal future cap, not the number of old--new edges.

## 2. Exact signs attain the `n^(3/4)` threshold

Let `n=4^j`, write `q=sqrt(n)`, and take a regular symmetric Hadamard signing
`mathcal H` satisfying

```math
\mathcal H^2=nI,
\qquad \mathcal H\mathbf1=q\mathbf1,
\qquad \operatorname{tr}\mathcal H=0.            \tag{OV.6}
```

Put `A=mathcal H-diag(mathcal H)` and

```math
S=qn=n^(3/2),
\qquad m=\lfloor n^(3/4)\rfloor,
\qquad a={m\over n}.                              \tag{OV.7}
```

### Lemma OV.2 (a biased flat exact-sign bridge)

There is an exact sign matrix `B in {+-1}^{n times m}` such that

```math
B=aJ+E,
\qquad ||E||_(2->2)\le C\sqrt n                 \tag{OV.8}
```

for an absolute constant `C`.

#### Proof

Choose the entries of `B` independently with mean `a`; this is possible for
`0<=a<=1`.  The centred entries of `E=B-aJ` are independent, mean zero, and
uniformly bounded.  The standard rectangular subgaussian norm bound gives

```math
||E||\le C(\sqrt n+\sqrt m)\le2C\sqrt n
```

with positive probability.  Equivalently, this follows directly from
Hoeffding on fixed unit vectors and `1/4`-nets of the two Euclidean spheres.
Absorb the factor two into `C`. `square`

Let `C_m=J_m-I_m` be the positive clique and define the two complete
exact-sign parents of order `N=n+m` by

```math
P_\sigma(x,y)
=\sigma {1\over2}x^T\mathcal Hx+x^TBy
 +{(\mathbf1^Ty)^2-m\over2}.                      \tag{OV.9}
```

The diagonal deletion in the old block is harmless because of (OV.6).

### Theorem OV.3 (sharp orientation exposure at sublinear width)

The parents (OV.9) satisfy

```math
Q(P_+)\ge(2-o(1))S,
\qquad
Q(P_-)\le\left({5\over4}+o(1)\right)S.           \tag{OV.10}
```

Consequently

```math
Q(P_+)-Q(P_-)
\ge\left({3\over4}-o(1)\right)n^(3/2).           \tag{OV.11}
```

Both parents have cap `O(N^(3/2))`.  Together with Theorem OV.1 this proves
that `n^(3/4)` is the sharp exponent for exposing orientation by an
arbitrary exact-sign quadratic continuation.

#### Proof

The error bridge obeys uniformly

```math
|x^TEy|
\le||E||\sqrt{nm}=O(n\sqrt m)=O(n^(11/8))=o(S).  \tag{OV.12}
```

At `x=1,y=1`, the positive orientation has energy

```math
{S\over2}+m^2+{m^2-m\over2}+O(n^(11/8))
=(2-o(1))S,                                       \tag{OV.13}
```

because `m^2/S` tends to one.  This proves the first half of (OV.10).

For the negative orientation put

```math
p={\mathbf1^Tx\over n},
\qquad s={\mathbf1^Ty\over m},
\qquad \lambda={m^2\over S}\le1.                \tag{OV.14}
```

Decompose `x=p1+x_perp`.  Regularity and the spectral bound give

```math
-{1\over2}x^T\mathcal Hx
\le S\left({1\over2}-p^2\right).                 \tag{OV.15}
```

For the positive outer channel, after dividing by `S` and absorbing
(OV.12), the remaining scalar envelope is

```math
{1\over2}-p^2+\lambda ps+{\lambda\over2}s^2
\le {1\over2}+left({\lambda^2\over4}
                    +{\lambda\over2}\right)s^2
\le {5\over4}.                                   \tag{OV.16}
```

For the negative outer channel, the old term is at most `S/2`, and

```math
\lambda(-ps-s^2/2)
\le\lambda(|s|-s^2/2)\le {1\over2}.              \tag{OV.17}
```

This channel is at most `(1+o(1))S`.  Equations (OV.16)--(OV.17) prove the
second half of (OV.10).  Finally, the child, biased bridge, error bridge, and
clique have caps at most `S/2`, `m^2`, `o(S)`, and `m^2/2`, respectively,
so both complete parents have `O(S)=O(N^(3/2))` cap. `square`

## 3. Interpretation

There are now two distinct continuation thresholds.

1. If every component must itself have natural `m^(3/2)` cap, a sublinear
   component cannot expose orientation at the old target scale.
2. If only the final parent must have `O(n^(3/2))` cap, a smaller component
   may spend quadratic internal cap.  Width `n^(3/4)` is then necessary and
   sufficient.

Thus continuation width alone is not the dynamic resource.  The correct
budget is the pair `(width, internal cap)`.  The construction is a strict
near-original benchmark: it uses complete exact signs and a bounded-cap
parent, but does not concern minimizers or prove a recurrence for `M_n`.
