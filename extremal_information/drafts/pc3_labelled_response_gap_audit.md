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
