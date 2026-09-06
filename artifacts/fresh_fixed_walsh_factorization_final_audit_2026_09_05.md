# Fixed Walsh amplification: exact factorization boundary

Date: 2026-09-05. This targeted audit neither proves `R=T` nor exhibits
a strict `R<T` seed. It records the exact norm convention, the primary
theorem boundary, and the finite realization step still required.

## 1. The two norms

For a symmetric real seed `B`, use

```
R(B) = (1/2) sup_s beta(H_s tensor B)/s^(3/2),
T(B) = (1/2) min_(P>=B,P>=-B) max_(x signs) x^T P x.
```

The generators in `H_s` are the fixed regular symmetric order-four and
order-144 Hadamards. The Walsh-only subfamily uses order-four powers.
The already proved inequalities are `R<=T`; the exact Hilbert
factorization identification is

```
2T(B) = Gamma_2(B : ell_infinity^n -> ell_1^n).
```

This is not the entrywise factorization norm for `ell_1->ell_infinity`,
and not its SDP dual (the ordinary quantum XOR bias). The five-cycle
Seidel example in the earlier audit already separates `R` from half the
ordinary quantum bias: `R<=12sqrt(5)/5<5sqrt(5)/2`.
Consequently a complex-phase/quantum rounding import at constant one
cannot be silently substituted for the required real operator norm.

## 2. Exact primary Fourier-type check

Aicke Hinrichs, *Hilbert space factorization and Fourier type of
operators*, Studia Mathematica 145 (2001), 199--212,
[publisher PDF](https://www.impan.pl/shop/en/publication/transaction/download/product/89984).
Theorems 1.1--1.3 and 6.1 were inspected directly.

Theorems 1.1--1.2 compare Fourier and all-orthogonal amplification on
vector-valued **L2-to-L2** spaces; the latter characterizes Hilbert
factorization. Theorem 1.3 gives nonuniformity for a prescribed matrix
at the same dimension. Theorem 6.1 includes Walsh gradations. These are
not the endpoint **L-infinity-to-L1** norms defining `R`. The finite-
dimension nonuniformity also does not refute a supremum over arbitrarily
large catalyst dimensions. Thus no equality or strict fixed-seed
separation for the present `R` is imported from this paper.

The source's counterexample operators act from `ell_1^n` to
`ell_infinity^n`, not in our reversed convention, and its standing
Banach-space convention is complex. Its historical open-
problem statement is not asserted here to describe current literature.

## 3. A finite exact dual description of the missing coupling

Let `Cut_n=conv{xx^T : x in {+/-1}^n}`. Direct finite-dimensional SDP
duality gives

```
2T(B) = max_(C in Cut_n) ||C^(1/2) B C^(1/2)||_* .       (1)
```

Indeed the primal for `2T` minimizes `t` with `P>=+/-B` and
`x^T P x<=t`. Let the dual PSD matrices be `X,Y` and the nonnegative
cube multipliers be `p_x`. Stationarity is
`sum p_x=1`, `X+Y=sum p_x xx^T=C`. The objective is `tr B(X-Y)`.
For fixed `C`, write `X-Y=C^(1/2) K C^(1/2)` with `-I<=K<=I` on
the support of `C`; maximizing is the nuclear norm in (1).
Strict primal feasibility is obtained by a sufficiently large scalar
`P` and then a strictly larger `t`, so no duality gap is used implicitly.

For singular `C`, positivity of `X,Y` and `X+Y=C` forces both to vanish
on `ker(C)`. The contraction `K` is defined using the inverse square
root only on `range(C)`; extend it arbitrarily by a contraction on the
kernel. Conversely every such `K` gives
`X=C^(1/2)(I+K)C^(1/2)/2`, `Y=C^(1/2)(I-K)C^(1/2)/2`.
The maximizing `K` is the sign of the symmetric compressed matrix,
with any orthogonal completion on its kernel. All feasible dual
variables are bounded (`tr C=n`, `0<=X,Y<=C`), so the dual maximum
is attained. On a bounded primal sublevel, `P>=0` and `tr P<=t`
(average the cube inequalities), giving primal attainment as well.

The cut matrix `C` needs at most `n(n-1)/2+1` Boolean atoms, by
Caratheodory in the fixed-diagonal affine space. Approximate their
weights by rational numbers and form a large uniform-row Boolean
matrix `F` with `F^T F/N` close to `C`. For this fixed `F`, allowing
**every** orthogonal outer matrix `U` gives

```
sup_(U orthogonal) |tr[U F B F^T]|/N
 = ||F B F^T||_*/N
 = ||(F^T F/N)^(1/2) B (F^T F/N)^(1/2)||_*.
```

The optimum can even use a symmetric orthogonal polar completion.
This explains exactly why unrestricted orthogonal amplification gives
`T`. One sufficient operational route to `R=T` is to realize, for each
seed `B`, an optimizing polar coupling using the prescribed family of
flat Hadamard outers and Boolean rows, to arbitrarily small objective
error. Merely representing `C` by finite Boolean atoms does not give
that step. Equality of the norms would not require every individual
cut-covariance/contraction pair to be realized: a different optimizer
or convex approximation of the relevant support functional could
suffice. No such realization or approximation is proved here.

The approximation is also valid when `C` is singular: square root is
continuous on the positive semidefinite cone in finite dimension and
the nuclear norm is continuous. The atom weights may be rounded with
denominator any sufficiently large allowed outer order, including
`4^a`; their total rounding error tends to zero. None of this reduces
the seed dimension `n` or supplies the missing prescribed-outer
realization. It only makes that obligation mathematically explicit.

As a sanity check, if `B>=0`, taking `P=B` gives `T(B)=q(B)`;
regularity of the actual outers gives `R(B)>=q(B)`, while `R<=T`
gives equality throughout. The same is true for `B<=0`, by changing
the global sign.

## 3a. The full-sign floor and its equality case

Formula (1) gives a short exact interpretation of the certificate
floor. For any symmetric full sign `B`, its `n` columns are Boolean,
so `C=B^2/n` belongs to `Cut_n`. Hence

```
2T(B) >= tr |B|^3/n >= n^(3/2),
```

where the last step is the power-mean inequality and
`tr B^2=n^2`. Equality in that step requires every singular value
to equal `sqrt(n)`. Conversely a symmetric Hadamard `B` has
`P=sqrt(n)I>=+/-B` and `2T(B)<=n^(3/2)`. Thus

```
T(B)>=n^(3/2)/2,
equality if and only if B^2=nI.
```

An absolute-PSD-majorant upper certificate for a full sign seed can
therefore never certify a normalized value below `1/2`. A proof of
`R=T` would make this the exact regularized seed floor; it would not
provide the separate upper-recovery relation needed for convergence
of the original minima.

The hollow analogue is exact as well. For a symmetric hollow signing
`A` of order `n>=2`, complete the zero in each of its columns by `+1`
and by `-1`, and average the `2n` resulting Boolean outer products.
This gives the feasible cut covariance `C=(A^2+I)/n`. Formula (1) yields

```
2T(A) >= (1/n) sum_j (s_j^3+s_j),
```

where `s_j` are the singular values of `A`. Set `a=sqrt(n-1)>=1`.
The following exact scalar identity supplies the required supporting
line, including near zero:

```
s^3+s - (a^3+a) - [(3a^2+1)/(2a)](s^2-a^2)
 = (s-a)^2 [s+(a^2-1)/(2a)] >=0.
```

Since `sum s_j^2=n(n-1)=n a^2`, summing gives

```
T(A)>=n sqrt(n-1)/2.
```

Equality forces all singular values to equal `sqrt(n-1)`; at `n=2`
the same conclusion follows using their fixed sum of squares. Conversely
a symmetric conference matrix attains the bound with majorant `aI`.
Thus equality holds precisely when `A^2=(n-1)I`. Again this concerns
`T`, an upper relaxation of `R` and `q`, so it must not be read as the
same inequality for those smaller norms.

## 4. What Schmidt supplies

The separately audited primary theorem of Schmidt gives asymptotically
flat scalar Walsh spectra, hence normalized maximal odd-Walsh excess
tending to one and `R(H2)=sqrt(2)`. See
`fresh_schmidt_odd_walsh_regularization_2026_09_05.md` for the exact
same-spin lift and source normalization.

In its actual coset construction, Fourier transport on the selected
balanced coset functions approaches inversion of coset labels. Applying
this coordinatewise to Boolean vector labels therefore transports the
same label distribution; it does not automatically realize the polar
coupling in Section 3. A vector-valued extension with precisely that
additional conclusion has not been obtained in this audit.

Even the weaker full-sign landing statement
`R(B)>=n^(3/2)/2` for every symmetric full sign seed would be useful:
it would identify the minimum regularized limit as `1/2`, without
requiring `R=T` for arbitrary real seeds. No theorem establishing this
weaker landing was found either. Neither result alone is an original
convergence proof without the already identified near-minimizer
recovery/tensor-restriction step.

## 5. Finite nonregularizability is not tensor separation

The new primary order-36 Hadamard examples with bilinear excess 204 or
208, below spectral 216, refute exact regularizability at that order.
They do not supply a positive defect uniform under all fixed Walsh
amplifications. The current real-symmetric equal-square examples do
prove a different statement: their regularized caps differ, but no
strict `R<T` value for either seed is established. These distinctions
must be preserved when using the finite certificates.

There is also an elementary reason not to replace fixed catalysts by
self-tensor asymptotics. For every symmetric full sign seed `B`, the
Boolean vector `vec(B)` satisfies

```
q_full(B tensor B) >= tr(B^4)/2 >= n^3/2.
```

The last inequality uses `tr B^2=n^2`. Thus self-squaring always
reaches normalized at least `1/2`, regardless of a possible smaller
normalized cap of the original seed. It cannot serve as an
upper-preserving recovery operation for such an advantage.
