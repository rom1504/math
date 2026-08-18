# Independent audit of macroscopic child-marginal spread

Status: **passed**.  This note independently verifies
[`actual_child_macroscopic_marginal_spread.md`](actual_child_macroscopic_marginal_spread.md).
The theorem correctly distinguishes sector conditioning on exterior child
spins from an unconditional vertex marginal.  Its conclusions do not extend
to a negative-bridge posterior without an additional argument.

## 1. Subset-flip identities

For a uniform `r`-subset of a `k`-vertex block, an internal edge is cut with
probability `2r(k-r)/(k(k-1))`, so its character has mean multiplier
`1-4r(k-r)/(k(k-1))`.  A cross edge is cut with probability `r/k`, giving
multiplier `1-2r/k`.  This proves MS.3 with the constants in MS.2.

Uniform averaging over the exterior annihilates both its quadratic
Hamiltonian and every cross monomial, hence

```math
H_U(u)=E_{v'}H_A(u,v').                              \tag{VMS.1}
```

Also `H_V(-v)=H_V(v)` while the cross term reverses sign, proving the second
identity in MS.4.  Both affected pieces consequently have absolute value at
most the global cap `K_A`.  Multiplication by either sector sign gives MS.6;
no positivity assumption on either piece is being used.

## 2. Conditional versus marginal law

For fixed exterior configuration `v`, all points `u^R` on the radius-`r`
Hamming sphere are distinct.  Restricting the conditional partition sum to
that sphere and applying Jensen to MS.6 gives

```math
\max_u\mu_{A,s}(X_U=u\mid X_V=v)
\le {k\choose r}^{-1}e^{t(a_{k,r}+b_{k,r})K_A}.      \tag{VMS.2}
```

Thus MS.8 is genuinely uniform in `v`.  The marginal is a convex mixture of
these conditional laws, so each marginal atom obeys the same ceiling.  This
proves MS.9 directly.  The alternative log-sum-exp proof in MS.10 has the
correct Jensen direction: log-sum-exp is convex and coordinatewise
increasing.

## 3. Constants and uniform exponential rate

From `tK_A/m<=C_beta` and `k>=theta m`,

```math
{tK_A\over k}(a_{k,r}+b_{k,r})
\le {C_\beta\over\theta}
 \left(4q_k(1-q_k){k\over k-1}+2q_k\right),          \tag{VMS.3}
```

which, combined with the standard type-class lower bound, is exactly
MS.12.  For any fixed `q<1/2`, `r=floor(qk)` converges uniformly along
`k>=theta m`; taking the supremum after this uniform lower bound justifies
MS.13 and MS.15.

Put `D=C_beta/theta` and `q_*=e^(-6D)`.  Since

```math
h(q)\ge q\log(1/q)+q(1-q),
```

one obtains

```math
h(q_*)-D(6q_*-4q_*^2)
\ge q_*(1-q_*)+4Dq_*^2
\ge q_*(1-q_*).                                     \tag{VMS.4}
```

Thus the explicit constant in MS.14 is valid.  Taking half this positive
margin absorbs the floor, type-class polynomial factor, and finite
`k/(k-1)` correction uniformly, proving MS.16--MS.17.

## 4. Common-block consequences

If two independent block laws each have maximum atom at most `e^(-eta k)`,
then their coincidence probability is at most that maximum atom.  A union
over equality and antipodal equality gives the factor two in MS.18.
Conditioning on the second block word and covering the two radius-`delta k`
balls gives at most

```math
2(k+1)e^{k h(\delta)}e^{-\underline\eta k},          \tag{VMS.5}
```

which is MS.19.  The stated catalogue consequence is the same atom bound
followed by a union bound.

## 5. Scope and source correction

The proof uses actual-child optimality only through the scalar cap
contraction.  It establishes spread for every positive-density coordinate
marginal and even after exact conditioning on all complementary child
spins.  It does not imply independence, correlation decay, stability under
the negative bridge tilt, or a bound on row-product lifetime.  The source's
scope statement is therefore accurate.

The draft initially contained one rendering/control-character defect in the
subscript of MS.15.  It was replaced during integration by the ordinary
condition `A\ {\rm actual}`.  This was purely typographical and did not alter
the quantified statement.

## 6. Verdict

All mathematical identities, conditional-to-marginal implications,
constants, rates, and scope claims pass independent verification.  No
normalization or logical correction is required.
