# Independent audit of `favorable_gram_r38.md`

## Verdict

The main algebraic bridge is **PASS after scope/wording corrections**.  The
Parseval factors, Bhatia--Davis inequality, complement-deficit bound and its
`0<=Delta<=2q_n` range, external-Gram exponent, and independent-noise
formulas are correct.  The checker reruns successfully.

Four statements should be narrowed before ledger integration:

1. The endpoint bound with `q_n` applies to an exact minimizer `A`.  For an
   arbitrary signing, replace `q_n` by `Q(A)`; only the Parseval and
   Bhatia--Davis parts are signing-universal as written.
2. From an **average** deficit hypothesis one obtains average cap
   persistence, not pointwise persistence until Markov thinning is applied.
3. Completion data do not *determine* the internal Gram term and Parseval
   does not bound it at project scale.  The broader phrase that no estimate
   based on `(w,a,b,e)` can control it is too strong: exact-ground stability
   itself gives `||A[S]y||^2<=(m-1)e`.  The finite examples establish
   non-determination, not impossibility of every additional inequality.
4. Independent noise certainly loses near-ground/near-parent energy at the
   displayed scaling.  It does not, without specifying the lifted deficit,
   falsify every notion of a favorable label, because the favorable loss also
   contains `-p_2 Delta-B_(n,m)`.  State the obstruction for preservation of
   `Theta(n^(3/2))` child energy or of a project-deficit complement
   certificate.

Also define `beta=U_v(G)` in the sufficient package and retain the ambient
route assumption `0<c'<1/4`, under which the noise baseline `Theta(n^2)`
fits the target.

## Derivation audit

Let `T=S^c`, `|T|=t`, and

```math
f(x)=x^TA[T]x+2y^TA[S,T]x.
```

The quadratic Walsh coefficients are `2A_ij`, one for every unordered edge
of `T`, and the linear coefficients are the entries of `2A[T,S]y`.
Orthogonality therefore gives exactly

```math
E f^2=4 binom(t,2)+4||A[T,S]y||^2
       =2t(t-1)+4||A[T,S]y||^2.
```

Writing the range as `[-Z_-,Z_+]`, mean zero and
`(Z_+-f)(f+Z_-)>=0` give

```math
E f^2<=Z_+Z_-=w^2-a^2.
```

For an exact minimizer, every completion has energy `e+f(x)` in
`[-q_n,q_n]`.  Thus `Z_+<=q_n-e`, `Z_-<=q_n+e`, and

```math
w^2-a^2<=q_n^2-e^2.
```

The same endpoints give `|a+e|<=b`, so

```math
w^2-a^2<=(q_n-b)^2-(|e|-b)_+^2.
```

All signs and factors agree with (R38.1), (R38.4)--(R38.6).

For an oriented complement certificate, `E_d=q_n-Delta` and

```math
e>=q_n-Delta/2.
```

Since `|E_d|<=q_n`, exactly `0<=Delta<=2q_n`; hence the lower bound on `e`
is nonnegative and may be squared.  This gives

```math
4||A[S^c,S]y||^2+2t(t-1)
<=q_n^2-e^2
<=q_n Delta-Delta^2/4
<=q_n Delta.
```

Thus the factor `1/4` in

```math
E ||A[S^c,S]y^S||^2 <= (q_n/4) E Delta_S
```

is correct.  Also `e<=Q(A[S])` gives the pointwise cap-persistence inequality

```math
Delta>=2(q_n-Q(A[S])).
```

Under `E Delta=O(n^(3/4-c'))`, its immediate consequence is

```math
E[q_n-Q(A[S])]=O(n^(3/4-c')).
```

A constant-loss Markov thinning then makes the same scale pointwise on the
retained subfamily.

## Favorable-support and internal-Gram scope

For `p_2>=1/2`, every complement certificate is favorable by (10.1045), so
the labels used in the positive bridge are on the required support.  The
finite exact child grounds are also favorable: after positive orientation,
`Q(A[S])-e=0`, so their favorable loss is
`-p_2 Delta-B_(n,m)<=0` for any full extension.

The internal/external decomposition is exact, and the completion polynomial
contains the internal block only through the scalar `e`.  The finite data
correctly show non-determination.  In particular, on `A_8,m=5` there are
records with the identical tuple

```text
(w,a,b,e,external)=(8,0,12,12,3)
```

but internal Gram `32` or `40`.  The sharper `A_9` group with
`(w,a,b,external)=(4,-2,20,0)` has internal Gram `46` or `86`, but its two
internal values correspond to different energies (`e=14` and `e=22`).  It
therefore supports the claim for `(w,a,b)`, not by itself for the full tuple
including `e`.

## Noise audit

For `Y_i=y_i z_i`, independent with `E z_i=theta`, the Gram matrix of the
selected full columns has diagonal `n-1`.  Therefore exactly

```math
E||A[:,S]Y||^2
=theta^2||A[:,S]y||^2+(1-theta^2)m(n-1),
E[Y^TA[S]Y]=theta^2e.
```

Starting from the worst allowed exact-ground scale `Theta(n^(5/2))`, a
target `O(n^(9/4-c'))` forces
`theta^2=O(n^(-1/4-c'))`; for `0<c'<1/4` the diagonal baseline is affordable.
A `Theta(n^(3/2))` energy then falls to `O(n^(5/4-c'))`.  This rigorously
blocks generic independent noise from preserving a near-ground or
near-parent complement certificate while guaranteeing the Gram target.  It
does not alone rule out a different favorable lift with a large deficit.

## Checker

```text
.venv/bin/python tmp/favorable_gram_r38_check.py
PASS favorable_gram_r38_check
```

All reported finite values use integer enumeration; no floating-point LP is
involved.
