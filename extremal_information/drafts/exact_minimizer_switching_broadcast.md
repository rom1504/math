# Exact minimizers have a target-scale switching broadcast before physical exposure

**Status:** rigorous task-local theorem; independently audited after the
boundary-profile/scalar-cap scope repair.

This note separates the two clauses in the broadcast target left by the
mesoscopic affine-child no-go.  The coefficient-metric clause is automatic
and much stronger than previously required: every exact minimizer has an
exponential switching orbit at mutual Boolean distance `Omega(n^(3/2))`.
What remains is precisely bounded-cap physical exposure of that orbit.

## 1. Notation

For a hollow symmetric matrix `A`, put

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j,
\qquad
\|A\|_{\rm B}=\max_{x\in\{\pm1\}^n}|H_A(x)|.
```

For `u in {+-1}^n`, let

```math
A^u=D_uAD_u.
```

The labels `u` and `-u` give the same switched matrix.  Switching preserves
the complete Boolean energy multiset and in particular
`Q(A^u)=Q(A)`.

## 2. A bipartite Boolean-norm lower bound

### Lemma SB.1

If `R in {+-1}^{k times ell}`, then

```math
\max_{p\in\{\pm1\}^k,q\in\{\pm1\}^\ell}|p^TRq|
\ge {k\sqrt\ell\over\sqrt2}.                    \tag{SB.1}
```

The same argument after transposition gives the maximum of this bound and
`ell sqrt(k/2)`.

#### Proof

Choose `q` uniformly.  The sharp `p=1` Khintchine inequality gives, for
every row `i`,

```math
\mathbb E_q\left|\sum_jR_{ij}q_j\right|
\ge\sqrt{\ell/2}.
```

Hence some `q` has `sum_i |(Rq)_i|>=k sqrt(ell/2)`.  Choose
`p_i=sign((Rq)_i)`. `square`

Only the displayed elementary Rademacher inequality is used; no
Grothendieck or polarization loss enters.

## 3. The switching orbit is exponentially separated

### Theorem SB.2 (universal switching broadcast in old-block norm)

There is an absolute `c>0` such that, for every sufficiently large `n`,
every hollow order-`n` signing `A`
has a projective label set

```math
U\subseteq\{\pm1\}^n/\{\pm1\},
\qquad |U|\ge\exp(cn),                            \tag{SB.2}
```

for which

```math
\boxed{\|A^u-A^v\|_{\rm B}\ge {1\over4}n^{3/2}}
\qquad(u\ne v).                                  \tag{SB.3}
```

All children have exactly the same cap and the same unlabelled energy
landscape.  In particular, when `A` is an exact minimizer, every member of
the packing is an exact minimizer at the same order.

#### Proof

A standard random-code argument gives `U` of size `exp(cn)` with

```math
n/4\le d_{\rm P}(u,v)\le n/2                    \tag{SB.4}
```

for distinct labels.  For example, Hoeffding gives
`Pr{|u dot v|>n/2}<=2exp(-n/8)`, so any fixed `c<1/32` works by a union
bound.

Fix `u,v` and replace `v` by `-v` if needed.  Put `s=u odot v` and
`S={i:s_i=-1}`.  Then `k=|S|` lies in `[n/4,n/2]` and
`ell=n-k>=n/2`.  Switching by `v` preserves Boolean norm, while

```math
D_v(A^u-A^v)D_v=A^s-A.
```

This difference vanishes inside `S` and its complement and equals `-2A`
on the complete bipartite block between them.  Lemma SB.1 therefore yields

```math
\|A^u-A^v\|_{\rm B}
\ge2k\sqrt{\ell/2}
\ge {1\over4}n^{3/2}.                            \tag{SB.5}
```

The switching action on complete signings is free modulo the global label:
a nonconstant `s` changes every cross-edge coefficient.  Thus all members
are distinct. `square`

The constant is deliberately coarse.  The important point is the fixed
target scale, uniformly over the old signing.

## 4. Exact pinned exposure and its physical ceiling

For a child `C` and `r in {+-1}^n`, define

```math
R_C(nr)=\max_x\{H_C(x)+n r^Tx\}.
```

The elementary pin identity is

```math
R_C(nr)=n^2+H_C(r),                              \tag{SB.6}
```

because changing `d` coordinates can improve the quadratic part by at most
`2d(n-d)<=2nd`, exactly the field loss.

For every pair in Theorem SB.2, choose a displayed Boolean witness `r_(uv)`
and an ordering of `(u,v)` so that

```math
H_(A^u-A^v)(r_(uv))\ge n^{3/2}/4.
```

Then (SB.6) gives

```math
R_(A^u)(nr_(uv))-R_(A^v)(nr_(uv))
\ge n^{3/2}/4.                                  \tag{SB.7}
```

Thus every exact minimizer class has `Omega(n)` bits of exact pinned
response entropy without the perturbations used in Theorem 36.3.  This is
not yet a bounded-cap all-spins-free contextual packing: the literal
amplitude-`n` pin has a quadratic calibration, and the universal-pin barrier
prevents compiling it uniformly over arbitrary children.

Theorem SA.1 of the mesoscopic no-go makes the separation (SB.3) necessary
for any public physical exposure.  Theorem SB.2 proves that necessary
old-block clause in the strongest possible qualitative form.  The remaining
clause is a restricted anti-pin theorem for the switching orbit.

## 5. A precise near-minimizer lemma for boundary-response exposure

Theorem 21.8 supplies the boundary-response compiler under one tail
hypothesis.  After multiplying `A` globally by `-1` if necessary, put
`P=max_x H_A(x)=Q(A)`.  The following statement is therefore sufficient:

> **Uniform exact-minimizer tail deficit (`L_tail`).**  There exist fixed
> `d_0,kappa>0` such that every sufficiently large exact minimizer obeys
> ```math
> #\{x:P-H_A(x)<d_0n^{3/2}\}
> \le\exp\{(\log2-\kappa)n\}.                    \tag{SB.8}
> ```

Combining (SB.8) with Theorem 21.8 produces one exact-sign bridge of
operator norm `O(sqrt n)` and an `exp(gamma n)` subfamily of exact-minimizer
switchings whose **conditional boundary-response profiles** are pairwise
separated by `Omega(n^(3/2))`.  Filling the new--new block by any one common
bounded-cap signing gives complete exact-sign parents of order `2n` and cap
`O(n^(3/2))`; the same conditional profiles remain separated because the
fill adds one common function of the boundary spin.

This is not yet a scalar all-spins-free parent-cap packing.  Optimizing the
boundary spin can erase a sup-norm difference between two profiles.  A
further low-cap boundary selector/restricted anti-pin is required to expose
one of the separated coordinates without the quadratic calibration of the
amplitude-`n` pin.

`L_tail` is strictly less information than the original optimization: it
asks only for an exponential upper bound on one fixed-width upper level set
of each already-minimizing landscape.  It neither identifies an optimizer,
computes `M_n`, compares orders, nor determines the response outside that
level set.  Its truth would be a **negative** conclusion for contextual
compression: exact minimizers themselves would retain a linear bounded-cap
response rate.

A convenient stronger sufficient condition for (SB.8) is

```math
\|A\|_{2\to2}\le C\sqrt n                       \tag{SB.9}
```

uniformly over exact minimizers.  Invoke the known uniform lower bound
`M_n>=c_*n^(3/2)` and choose `d_0<c_*`.  Since
`||A||_F^2=n(n-1)`, the Hanson--Wright inequality at the resulting fixed
positive fraction of `n^(3/2)` then gives (SB.8).  Condition (SB.9) is not
asserted here; sparse-edit examples show that an edit-robust analogue over
arbitrary vanishing near-minimizers would need qualification.

## 6. Frontier classification

- **PROVES / weakens an arrow:** the old-block norm-packing clause of
  `L_broadcast` is universal and no longer missing.
- **Remaining SML:** for boundary-response incompressibility, prove
  `L_tail`.  For all-spins-free physical incompressibility, one must
  additionally construct a low-cap selector/restricted anti-pin for the
  resulting switching profiles.
- **Benchmark:** Level 5 for exact-minimizer coefficient and pinned response;
  physical bounded-cap response remains conditional.
- **No convergence implication:** the theorem establishes information
  heaviness, not cross-order transfer.

This is compatible with Theorem 36.3, which concerns an exponential family
inside every positive halo and does not use switching copies, and with BCX,
which proves bounded-cap anti-pin exposure for one special structured
switching family.  It does not assume the mesoscopic shell packing and does
not repair the AO affine children ruled out by SA.3.
