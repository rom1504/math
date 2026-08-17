# Proof audit: the fixed-seed parity-pole no-go

**Audit target.** `fixed_seed_amplification_no_go.md`.

**Verdict.** The theorem is correct for its declared common-factor,
odd-monomial amplification class. The exposed full-product pole supplies a
genuine no-go, and recursive selector products make that pole worse. The
statement must remain scoped exactly as written: it is not a theorem about
arbitrary exact-sign lifts or full trust optima.

## 1. Seed normalization

For the order-16 regular Hadamard seed, `H_0^2=16I`; hence
`T_0=H_0/4` is a self-adjoint contraction. Direct exact evaluation gives

```math
{z_S^TH_0z_S\over4\cdot16}={52\over64}={13\over16}
```

for all sixteen odd seed products and value `32/64=1/2` for the prescribed
majority witness. Thus every odd pole, including whatever odd pole parity
selects after padding, has the same ratio. The proof does not silently
assume that the five generators alone exhaust the odd shell.

## 2. Parity and Fourier support

If every lifted port is an odd seed monomial, exponent vectors live in
`F_2^5`. The symmetric difference of an odd number of odd vectors has odd
parity. The common Boolean multiplier occurs `p` times and therefore remains
`u`. This proves PA.4 even when ports are duplicated, negated, or linearly
dependent.

For odd `p`, the top-degree majority coefficient is

```math
(-1)^((p-1)/2)2^{-(p-1)}\binom{p-1}{(p-1)/2},
```

so it is nonzero. Linear dependence can make the resulting product column
duplicate another column, but cannot remove the query from the Fourier
identity. Consequently the maximum over all declared individual deficits
does inspect this pole.

The selector-preservation premise PA.3 is not used to prove parity. It is
used only to identify the joint witness's base Rayleigh ratio as `1/2`.
The class is nonempty: the original five ports belong to it, as do odd
replication and arbitrary cancelling pairs of odd seed monomials.

There is also a stronger one-block argument on the concrete seed support.
The eight top-energy exponent masks form the subgroup
`K={0,5,9,12,17,20,24,29}`. Seed rows `3` and `10` agree on every
`K`-monomial but the prescribed selector has opposite values there. Hence
arbitrary monomial ports `z_(A_j) tensor u_j` which reproduce that selector
must use at least one `A_j notin K`. A singleton majority Fourier
coefficient is nonzero, so that port is active; tensor factorization and
the enumerated next shell bound its Rayleigh ratio by `13/16`. This
validates PA.1b without a common auxiliary factor or a formal full-cube
oddness assumption.

Averaging removes the possibility of rotating this non-top port label
between mixture blocks. Quotient the exponent map by `K`. It is nonzero on
every block. Some nonzero scalar functional of that quotient is either
constant one or balanced on the affine hyperplane of odd Fourier subsets
`B`. Thus at least half the active product channels leave `K` on each
block. Their deficits are at least `3/16`; all remaining deficits are
nonnegative. Fubini's interchange of the finite channel average and the
block mixture proves a common global channel with deficit at least `3/32`.
This validates PA.1c even when monomial classes, auxiliary factors, and
port labels vary by block.

## 3. Tensor arithmetic and positivity

Let `s=<u,Su>`. Since `S` is a self-adjoint contraction and `u` is a unit
Boolean function, `s in[-1,1]`. Tensor factorization gives exactly

```math
r_X={s\over2},\qquad r_Z={13s\over16}.
```

Therefore

```math
(1-r_Z)-{3\over8}(1-r_X)={5\over8}(1-s)\ge0.
```

No absolute value is being inserted: these are positive-roof Rayleigh
deficits, exactly the quantities in the robust synchronization premise.
The conclusion would need reformulation for an absolute-Rayleigh query.

The tensor child remains a contraction even if `S` is indefinite, so
`I-T_0 tensor S` is PSD. Equality occurs at `s=1`; hence `3/8` cannot be
improved within this class.

## 4. Recursive selectors

In a composition `F(g_1,...,g_L)` on disjoint leaf blocks, the coefficient
of the union of all leaf variables can only arise by choosing the full
outer Fourier monomial and the full monomial of every inner function.
It is therefore

```math
\widehat F([L])\prod_t\widehat g_t([5]).
```

For `g_t=Maj_5`, this is nonzero whenever the stipulated outer coefficient
is nonzero. Tensor Rayleigh values multiply, yielding `(13/16)^L`.
This proves the recursive claim without assuming anything about the
quadratic value of the composed selector itself.

## 5. Mixtures and perturbations

On a weighted disjoint union, normalized inner products are convex
combinations. Every seed-generated block satisfies the same linear
inequality, so their mixture does. A synchronized background block with
both deficits zero also satisfies it. This is a contraction-level direct
sum; it is not silently identified with a dense sign matrix.

If the two Rayleigh values have errors `e_X,e_Z`, each of magnitude at most
`eta`, then

```math
d_Z-{3\over8}d_X
={5\over8}(1-s)-e_Z+{3\over8}e_X
\ge-{11\over8}\eta.
```

The constant `11/8` is the correct worst-case triangle constant. An
operator-norm perturbation gives the required two quadratic estimates;
an edge-count or Frobenius estimate alone is not claimed to do so.

## 6. Exact-sign and hollowing audit

If `K` is entrywise signed then `H_0 tensor K` is entrywise signed. Because
`tr H_0=0`, its trace is zero regardless of `K`. For every Boolean `y`,

```math
y^T\operatorname {diag}(H_0\otimes K)y
=\operatorname {tr}(H_0\otimes K)=0,
```

so hollowing preserves all Boolean energies exactly. Removing the diagonal
changes operator norm by at most one. If the auxiliary Boolean `u` is a
top eigenvector, the tensor spectral roof is attained by a Boolean vector;
otherwise roof tightness is not asserted.

## 7. Claims deliberately not made

The proof does not cover:

1. a lift whose ports are not monomials of one retained seed factor;
2. a dense exact-sign cross-block completion with leading operator effect;
3. a growing exact-sign seed whose odd shell itself becomes near-top; or
4. the full maximum over Boolean spins in a completed trust parent.

These exclusions are substantive escape routes, not technical omissions.
Accordingly the correct research conclusion is: **stop trying to amplify
the fixed five-port seed by natural tensor/replica/dilution operations; a
vanishing-marginal exact-sign example requires a nonlocal lift or a new
growing seed.**
