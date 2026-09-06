# Actual-seed losses in the rank-one Hadamard weave

Date: 2026-09-06. Independent reconstruction for the renewed convergence
campaign. These are scoped transfer obstructions, not nonconvergence of the
original minima. The strict upper theorem is used only in Section 4; Sections
1--3 are unconditional finite identities.

Write `Q(A)=max_x |x^T A x|/2` for hollow symmetric signings. For a full
symmetric sign matrix, distinguish its uncentered quadratic cap from the cap
after deleting its diagonal.

## 1. The full weave has an exact hollow cap, for every seed

Let `m` be even, let every `H_i` be a real Hadamard matrix of order `m`, and
let `S` be ANY symmetric full sign matrix. Define

```math
K_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i).
```

Then `K^2=m^2 I` and `tr K=m tr S`. A fixed-point-free involution `pi`
gives Boolean eigenvectors of BOTH signs: choose `u_i` on one vertex of
each matched pair, set `u_{pi(i)}=sigma S_{i,pi(i)}u_i`, and put
`x_{i,a}=u_i H_i(a,pi(i))`. Direct column orthogonality gives `Kx=sigma m x`.
Consequently, for `A_K=K-diag K`,

```math
Q(A_K)=\frac{m^3+m|\operatorname{tr}S|}{2}.
```

The spectral upper bound and the two Boolean witnesses agree exactly.
In particular, choosing a balanced seed diagonal gives normalized hollow
cap exactly `1/2`, independently of the seed's original Boolean cap,
spectral norm, or near-optimality. This is stronger than merely observing
that a spectral certificate fails to improve `1/2`: these are actual Boolean
witnesses in the resulting actual signing.

## 2. Arbitrary balanced fibre restriction retains a Boolean floor

Retain arbitrary sets `T_i` of the same size `k`. The resulting full matrix
has order `N=mk`. For the same matching and the restricted column vectors,
choose one independent uniform sign `u` per matched pair. In the quadratic
Hamiltonian, edges joining different matched pairs have mean zero. Each
matched pair contributes exactly `sigma k^2`. The contribution inside
individual fibres is independent of `sigma`, because changing the relative
pair sign only changes an entire fibre's spin vector by global reversal.
After hollowing, denote this deterministic contribution by `c_pi`. Thus

```math
\mathbb E_u H_{A_{K,T}}(x^{(\sigma,u)})
 =c_\pi+\sigma\frac{mk^2}{2}.
```

Taking the larger absolute value over the two signs and using
`max |H| >= |E H|` proves the exact bound

```math
Q(A_{K,T})\ge \frac{mk^2}{2}+|c_\pi|
\ge\frac{mk^2}{2},\qquad
\frac{Q(A_{K,T})}{(mk)^{3/2}}\ge\frac12\sqrt{\frac{k}{m}}.
```

There is no diagonal error in this lower bound. It uses global reversal
correctly: the internal contribution cannot cancel both energy signs.
It holds for arbitrary Hadamard choices, arbitrary seed signs and arbitrary
selectors. It is only a floor for this rank-one weave, not for all signings.

An asymptotically cap-preserving transfer of seeds of normalized cap `c`
through balanced restrictions therefore requires `limsup k/m <=4c^2`.
At `p=31/32` the floor is `sqrt(31/32)/2`, approximately `0.4921254921257382`.
This does not conflict with the proved upper construction.

## 3. Where seed information is erased

Fix two full symmetric seeds `S,S'` with the same diagonal. Choose column
sign matrices `D_i` so that

```math
D_i(j,j)D_j(i,i)=S_{ij}S'_{ij}\quad (i\ne j).
```

This is always possible independently on each unordered pair. Replacing
`H_i` by `H_i D_i` makes the weave with seed `S'` IDENTICAL entrywise to
the weave with seed `S`. Selectors are unchanged. Therefore the class of
all outputs, when arbitrary Hadamard bases are allowed, does not depend on
the off-diagonal seed at all.

Likewise, if independent uniform signs are supplied to every Hadamard column,
the output distribution is identical for every deterministic seed with the
same diagonal. One may add these signs without changing any absolute row
spectrum used in the permanent certificate. Thus that seed-blind ensemble
can prove a uniform construction but cannot explain an improvement through
the small value of a particular input seed's `Q`.

Scope: a restricted recursive distribution need not itself be invariant under
arbitrary output signs until those signs are added. This statement does not
rule out a new, correlated seed-adapted ensemble, nor does it assert that
every norm or every construction called weaving is seed-blind.

## 4. The strict upper theorem makes the loss apply to actual minimizers

Assume the all-order bound in
`continued_director_strict_all_order_upper_2026_09_06.md`, with constant
`c_*<0.499432220485404`. Along even Hadamard orders, take ACTUAL minimizing
hollow seeds, and complete their diagonals with balanced signs. Section 1
then shows a normalized actual-output loss at least

```math
\frac12-c_*-o(1)>0.000567779514596-o(1).
```

More generally any attempted balanced restricted weave transfer preserving
all these near-minimizers must delete a fraction at least
`1-4c_*^2-o(1)`, approximately `0.0022698285640787`.

The same new upper theorem upgrades two previously non-near-minimizer-only
warnings to actual minimizing seeds:

1. The absolute-PSD-majorant norm satisfies
   `T(A)>=n sqrt(n-1)/2`, whereas actual minimizers satisfy
   `Q(A)<= (c_*+o(1))n^(3/2)`. Hence its relaxation gap on actual minimizers
   is at least `(1/2-c_*-o(1))n^(3/2)`. This is a gap for `T`, NOT a claimed
   gap for the smaller stabilized Hadamard norm `R`.
2. Complete an actual minimizing seed to any symmetric full sign matrix `B`
   of order `n`, and hollow `B tensor B`. The Boolean vector `vec(B)` has
   full energy `tr B^4`, while hollowing subtracts `(tr B)^2`. Therefore
   `Q(hollow(B tensor B)) >= (n^3-n^2)/2`, because
   `tr B^4 >= (tr B^2)^2/n=n^3`. This is again an actual Boolean cap loss
   at least `1/2-c_*-o(1)`, not merely an upper-certificate failure.

Neither observation establishes the value of `R`, forbids tensoring with a
different carefully chosen catalyst, or rules out a nonlocal sign-preserving
recovery operation.

### Reconstruction of the precise majorant inequality

The convention here is

```math
2T(A)=\min_{P\succeq A,\ P\succeq-A}\ \max_{x\in\{-1,1\}^n}x^TPx.
```

Complete the zero entry in each of the `n` columns of hollow `A` by both
possible signs, and average their `2n` Boolean outer products. The resulting
cut covariance is `C=(A^2+I)/n`. For any feasible `P`, its maximum Boolean
quadratic value is at least `tr PC`. Let `J=sign(A)` (assign any value in
`[-1,1]` on the zero eigenspace), and set
`X=C^(1/2)(I+J)C^(1/2)/2`, `Y=C^(1/2)(I-J)C^(1/2)/2`.
Both are PSD and `X+Y=C`. The two majorant inequalities give

```math
\operatorname{tr}PC\ge\operatorname{tr}A(X-Y)
=\frac1n\sum_j(s_j^3+s_j),
```

where `s_j` are the singular values of `A`. Thus no appeal to unproved
`R=T`, nor to `tr|A|` alone, is involved. Put `a=sqrt(n-1)>=1`. The identity

```math
s^3+s-(a^3+a)-\frac{3a^2+1}{2a}(s^2-a^2)
=(s-a)^2\left(s+\frac{a^2-1}{2a}\right)\ge0
```

and `sum_j s_j^2=n(n-1)=na^2` imply
`sum_j(s_j^3+s_j)>=n(a^3+a)=n^2 sqrt(n-1)`.
Taking the minimum over `P` proves the stated lower bound on `T`.
This independently reconstructs the argument in
`fresh_fixed_walsh_factorization_final_audit_2026_09_05.md`, Section 3a.

## Verification

`computations/transfer_adversary_seed_loss_2026_09_06.py` independently checks
the finite identities in integer arithmetic. It uses an archived actual
order-eight minimizing seed only as a finite illustration, and reports its
enumerated Boolean cap directly rather than trusting the stored value.
The asymptotic near-minimizer conclusion is the analytic argument above.
