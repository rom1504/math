# Audit of the all-parity rooted Walsh orbit extension

**Verdict: PASS.**  I found no mathematical counterexample.  The affine-
symplectic parameterization, the extension argument in both parities, and the
rooted-response normalization are correct.  The result should still be scoped
as an orbit classification and a rooted-response quotient; it does not prove
that the rooted fibre is necessary for every unrooted graph maximum.

## 1. Even-dimensional coordinates

For even `m`, the proposed splitting

```math
V=\langle e,\omega\rangle\perp W,qquad
B(e,e)=B(e,\omega)=1
```

exists.  The restriction to `W` is nondegenerate alternating.  For

```math
T\omega=\omega,qquad
Tw=Sw+B(t,Sw)\omega,qquad
Te=e+t+c\omega,
```

the potentially nontrivial pairing is

```math
B(Te,Tw)=B(t,Sw)+B(t,Sw)=0.
```

The self-pairings and the pairings with `omega` are also preserved.  In the
converse direction, fixing `omega` forces the coefficient of `e` in `Te` to
be one; the induced map on `omega^perp/<omega>` is symplectic, and
orthogonality to `Te` uniquely forces the displayed `omega` coefficient of
`Tw`.  Hence `(S,t,c)` is both sufficient and unique.  The formula remains
valid at `m=2`, where `W=0` and the two choices of `c` are exactly the two
coordinate permutations.

## 2. Extension theorem and edge cases

Equality of relation kernels makes the label map a well-defined linear
isomorphism, Gram equality makes it an isometry, and equality of rooted
fibres is exactly the condition that it respect the characteristic vector
when that vector lies in the presented span.

The odd-dimensional projection has kernel `<omega>`; the draft explicitly
checks the only possible ambiguity.  This also covers `m=1`, where the
alternating complement is zero-dimensional.

For even `m`, adjoining `omega` when it is absent is legitimate because

```math
B(omega,u)=B(u,u),
```

which the partial isometry preserves.  Once `omega` is present, write the
even part as `<omega> direct-sum P`.  The induced map on `P` extends by
symplectic Witt.  If the source span contains an odd vector
`z=e+w_0+b omega`, the identity

```math
ell(w)=B(w_0,w)+B(w_0',S_0w)
```

is precisely the compatibility needed for the shear
`t=w_0'+Sw_0`; the stated value of `c` then sends `z` to `phi(z)` exactly.
No nondegeneracy of the restricted form on `P` is being assumed.  This
settles the potentially dangerous degenerate-subspace case.

Exact enumeration supplied with the draft passes.  I additionally checked:

- all tuples of lengths at most five for `m=1` and `m=2`;
- all tuples of lengths at most four for `m=3`;
- **all 65,536 ordered four-tuples for `m=4`**, where the invariant gives
  exactly 2,076 orbit classes under the 48-element orthogonal group.

Thus the small exceptional dimensions and every possible subspace type at
`m=4` were exercised; no collision between the rooted state and the true
orbit partition occurred.

## 3. Rooted-fibre necessity

For even `m>=4`, `a=omega` and `b=e_1+e_2` are distinct nonzero isotropic
vectors.  Therefore their singleton Gram matrices and relation kernels
agree, while their rooted fibres differ.  Since every orthogonal map fixes
the characteristic vector, they cannot be orbit-equivalent.  This is a
scalable orbit obstruction.  The draft also proves semantic necessity for
the declared **rooted Walsh response** family; it correctly does not claim
necessity for every unrooted composition query.

The restriction `m>=4` in that proposition is essential: at `m=2`, the
suggested `b` equals `omega`, and there is no second nonzero isotropic
singleton.  The all-parity classification theorem itself remains valid at
`m=2`.

## 4. Independent response-normalization check

Write `q=2^m`, `n=q^2`, `F=W/q`, and

```math
H_c(x)={q\over2}x^T\widehat C_cx,qquad
(P_WH_c)(y)=\max_x\{H_c(x)+q x^TFy\}.
```

For `s_c(u,v)=(-1)^(u dot v+c dot v)`, both `s_c` and `y_c=Fs_c` are
Boolean.  At the matched query, `x=s_c` gives

```math
{q\over2}n+qn={3\over2}n^{3/2}.
```

At the crossed query, modulation reduces the objective divided by `q` to

```math
{1\over2}u^TFu+w^Tu,qquad \|u\|_2^2=n,qquad w^TFw=0.
```

Since `(2I-F)^(-1)=(2I+F)/3`, completing the square gives the upper bound
`4n/3`, hence crossed response at most `4n^(3/2)/3`.  The two reciprocal
queries make the response difference at least `+n^(3/2)/6` and at most
`-n^(3/2)/6`; half its oscillation is therefore at least
`n^(3/2)/6`, exactly as stated.  There is no missing factor of `q`, `1/2`,
or two in the projective normalization.

## 5. Scope and minor presentation notes

- The `O(k^2)` state is globally coordinated.  It is not an independently
  composable summary, because cross-Gram entries and new cross-relations are
  absent from separately stored component states.
- “Minimal extra datum” is valid in the stated orbit-theoretic sense.  It
  should not be upgraded to minimality among all possible semantic response
  summaries without a separate lower-bound theorem.
- In (AP.7), `perp` denotes that the two-dimensional plane is orthogonal to
  `W`; the plane itself is not alternating because `B(e,e)=1`.  This is
  mathematically consistent, though `direct-sum` with an explanatory phrase
  would be slightly harder to misread.

Subject to those already-respected scope boundaries, the draft is ready for
promotion.
