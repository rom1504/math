# Independent audit: PC.3 labelled sparse-flip response gap

This audits Theorem 21.64 / DS.6 independently of its derivation.

## Verdict

**PASS.**  The operator normalization, two trust channels, hollow transfer,
gap constant, and labelled-versus-unlabelled scope all check.

## Normalized operator and field

Write `r=sqrt(N)`, `T=H/r`, `e=x/sqrt(N)`, and `f=h/||h||_2`.  The exact
sparse-flip mean and matrix Bernstein give, on the same positive-probability
event used for the pole estimates,

```math
{H'\over r}=T-kappa ee^T+o_{op}(1).
```

PC.3 product closure gives `Te=e,Tf=f`.  The exact row law gives

```math
{||h||_2^2\over N}={7j\over4}+O(1),
\qquad
e^Tf={||h||_1\over\sqrt N||h||_2}longrightarrow
rho=\sqrt{2/\mathop{\rm pi}}.
```

For `m~lambda sqrt(N/j)`, the normalized field strength is
`b=m||h||_2/N -> lambda sqrt(7)/2`.

## Positive channel

For `u=y/sqrt N`, decompose under the `+-1` eigenspaces of `T` and put
`a=e^Tu`.  Writing `s=sqrt(1-rho^2)`, every positive-channel spherical
competitor satisfies

```math
J_+(u)le {1\over2}+bs+brho|a|-{kappa\over2}a^2+o(1)
\le {1\over2}+bs+{b^2rho^2\over2kappa}+o(1).
```

This permits every direction in the positive eigenspace; no Boolean or
spectral competitor is omitted.

## Negative channel

If `c=||P_+u||_2`, then

```math
J_-(u)le {1\over2}-(1-kappa/2)c^2+bc+o(1)
\le {1\over2}+{b^2\over2(2-kappa)}+o(1).
```

The unflipped Boolean vector `x=sgn(h)` attains at least
`1/2+brho+o(1)`.  Thus the minimum of the two displayed channel gaps is
positive under (21.378).  At `kappa=1/2,lambda=1/10`, the smaller gap is
approximately `0.014665>0.0146`.

The diagonal of `H'` is unchanged and has trace zero.  Hollowing therefore
preserves every Boolean quadratic energy exactly and changes the spherical
operator by only `o(1)` after normalization.

## Scope

The maximum includes every old-child Boolean spin and both quadratic signs,
so this is a genuine full response separation for the fixed labelled field.
It is not yet a separation between unconstrained appended parents: free
shore spins also optimize the endpoint label and may select another field.

## Free-shore corollary

Corollary 21.65 / DS.7 also passes.  Here `E X=1/2,E Y=0`, so the
all-positive endpoint has `||W_j1||_1/N_j>=1+j/2`.  There is no missing
factor two in the cross energy.  At multiplicity
`m_j~lambda sqrt(N_j/j)` this is
`(lambda/2+o(1))N_j^(3/2)sqrt(j)`, whereas the child contributes only
`O(N_j^(3/2))` and every internal exact-sign shore contributes `O(N_jj)`.
The no-go applies to every internal completion of this direct unbalanced
lift, but not to a redesigned balanced or cancelled cross representation.

## Microcanonical compiler

Theorem 21.66 / DS.8 passes.  For a row conditioned to have sum `g_i` and
an endpoint of coordinate sum `c`, the exact variance is

```math
\left(1-(g_i/s)^2\right){s^2-c^2\over s-1}\le2s.
```

The stated Jensen bound follows.  Hoeffding at deviation
`C s^(3/2)sqrt(N)` has tail `exp(-2C^2s)`, so a union bound over `2^s`
endpoints succeeds for an absolute `C`.  For `g=m_jh_j`, magnitude and
parity are automatic; the two excess terms are respectively
`N_j^(5/4)j^(1/4)` and `N_j^(5/4)j^(3/4)`, both `o(N_j^(3/2))`.

The scalar cap conclusion controls the whole cross-response roof and realizes
the target field exactly, but by itself does not prove endpoint stability.
The strengthened affine endpoint conclusion is audited next.

## Unconstrained-parent lift

Theorem 21.67 passes.  For a fixed shore endpoint and global absolute sign
`tau`, replacing the old spin by `tau y` turns the parent objective into one
of the two trust channels plus a shore term.  Hence the flipped parent is at
most the maximum labelled trust response plus `Q(C)`, while the target
endpoint gives the unflipped trust response minus `Q(C)`.

The uniform affine compiler estimate, field-`l_1` Lipschitzness, and evenness
reduce every endpoint to strength `|a|b`.  Both spherical upper bounds are
increasing on `0<=|a|<=1`, so the target-strength gap is uniform.  Finally
`2Q(C)=O(s_j^2)=o(N_j^(3/2))` and `N_j+s_j=(1+o(1))N_j`.  The direction and
normalization in (21.391) are therefore correct.  This compares two
constructed exact-sign parents, not their minima over signings.
