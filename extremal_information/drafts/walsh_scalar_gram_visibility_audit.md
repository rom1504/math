# Adversarial audit: scalar visibility of a Walsh Gram flux

Audit target:
[`walsh_scalar_gram_visibility.md`](walsh_scalar_gram_visibility.md).

Verdict: **PASS**.  The proof gives a genuine unrooted scalar separation.
The exact query is one unweighted triangle, the compared tuples agree in
self-pairings, relation kernel, and characteristic-root fibre, and the
normalizations yield the stated leading-scale gap.

## 1. State collision

For

```math
(a,b,c)=(e_1+e_2,e_3+e_4,a+b)
```

and

```math
(a',b',c')=(e_1+e_2,e_1+e_3,a'+b'),
```

both first pairs are independent and the third label is their sum.  Hence
both kernels are exactly `{000,111}`.  Every label has even weight.  The
first Gram matrix is zero and the second has zero diagonal and unit
off-diagonal.  At `m>=5`, neither rank-two span contains `omega`; both root
fibres are therefore empty.  Thus no omitted rooted distinction can explain
the scalar gap.

## 2. Weyl signs

With `J_a=M_aT_aF`, direct multiplication gives

```math
J_aJ_b=(-1)^{b\cdot b+a\cdot b}M_{a+b}T_{a+b}.
```

For even `a,b,c=a+b`, multiplying by `J_c` gives

```math
J_aJ_bJ_c=(-1)^{a\cdot b}F.
```

All three children commute with one another and with `F`, because their
self-parities vanish.  Thus the simultaneous-eigenspace reduction is valid;
it is not merely a trace or averaged statement.

## 3. Boolean saturation is proved, not inferred from spectrum

The good-state lower bound requires one Boolean vector common to four
`+1` eigenspaces.  The draft constructs it as follows.

1. The standard chirp `x_0(u,v)=(-1)^(u dot v)` is self-dual.
2. For `f_p=(p,p)` with even `p`, `D_(f_p)x_0` is also self-dual, so
   `C_(f_p)x_0=x_0`.
3. The two source labels and two target labels are independent totally
   isotropic pairs, and both spans avoid the characteristic vector.
   Characteristic-rooted Witt extension therefore gives an ambient
   orthogonal coordinate permutation between them.
4. Transporting `x_0` back gives the required Boolean section.

Witt extension is used to construct a coordinate permutation, not as a
query.  The query itself has no root, field, or pinned spin.  The verifier
also constructs the transport explicitly at `m=5` as the product of the
three isotropic transvections with axes `0x44,0x47,0x125`, then checks the
Walsh eigen-identities by exact integer FWHT.  Higher `m` follow directly
from the abstract extension argument (or by tensoring with standard
two-coordinate chirps).

## 4. Spectral ceiling and scaling

On a joint internal eigenspace, write `f` for the `F_E` sign and `lambda_i`
for the child signs.  Unit flux gives

```math
(f\lambda_1)(f\lambda_2)(f\lambda_3)=-1.
```

After multiplying the sector matrix by `f`, the only cases are:

- three negative diagonal signs, with norm `2`;
- one negative diagonal sign, with eigenvalues
  `0,(1+sqrt(17))/2,(1-sqrt(17))/2`.

Therefore the normalized block operator norm is at most
`(1+sqrt(17))/2`.  A Boolean three-block vector has norm squared `3n`, and
the original energy is `q/2` times its normalized Rayleigh quadratic.  The
bad upper coefficient is consequently

```math
{3\over2}{1+\sqrt{17}\over2}
={3(1+\sqrt{17})\over4}.
```

The good common section saturates three child terms of coefficient `1/2`
and three bridge terms of coefficient `1`, giving `9/2`.  Their difference
is `3(5-sqrt(17))/4>0`.  All coefficients multiply
`qn=n^(3/2)`.

## 5. Scope

The result proves that **some** off-diagonal Gram information survives in
the minimal scalar semantic quotient.  It does not prove entrywise recovery
of a general Gram matrix or full scalar minimality of `(G,R)`.  The exposed
quantity is naturally the bilinear flux on a rank-two span presented by the
triangle relation.  The draft states this limitation explicitly.
