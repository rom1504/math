# Independent audit: regular-Hadamard amplification limit

**Verdict: PASS WITH MINOR SCOPE AND VERIFIER REPAIRS.**  The monotone
embedding, the factor `1/2`, the operator-norm ceiling, and the compact
response-carrier argument are all correct.  The Walsh and hollow-sign
corollaries are also correct for the stated zero-extended, scalar-onsite
Walsh language.  Before canonical promotion, the draft should (i) make that
continuation boundary explicit, (ii) distinguish the positive-eigenvalue
assumption needed by the upper functional from the weaker assumption needed
by the absolute functional, and (iii) either verify or separately cite the
stated `104/64` wind-tunnel value.  The result is a valid near-original
restricted-model theorem, but it is an elementary exact tensor-direct-limit
mechanism rather than an all-order or optimizer-transfer theorem for `M_n`.

## 1. Normalization and monotonicity

Let `s_r=h^r`, `N_r=ds_r`, and `B_r=B tensor H_r`.  If

```math
Hu=(\sqrt h)u,\qquad u\in\{\pm1\}^h,
```

then

```math
u^THu=(\sqrt h)||u||_2^2=h^{3/2}.
```

Consequently, for every Boolean `x` at depth `r`,

```math
(x\otimes u)^T(B_r\otimes H)(x\otimes u)
=h^{3/2}x^TB_rx.                                      \tag{A.1}
```

The denominator obeys

```math
2N_{r+1}^{3/2}=2(hN_r)^{3/2}=h^{3/2}(2N_r^{3/2}),      \tag{A.2}
```

so every *signed* normalized Rayleigh value at level `r` occurs exactly at
level `r+1`.  Taking either the maximum or the maximum absolute value proves
both monotonicities with precisely the draft's `1/2` convention.  In
particular, when the amplified matrix is hollow,

```math
{1\over2}\max_x|x^TAx|
=\max_x\left|\sum_{i<j}A_{ij}x_ix_j\right|=Q(A),       \tag{A.3}
```

so the absolute normalization agrees with the main project.

The common ceiling is also exact as an inequality:

```math
||B_r||_{2\to2}=||B||_{2\to2}h^{r/2},
\qquad ||x||_2^2=dh^r,
```

and hence

```math
{ |x^TB_rx|\over2(dh^r)^{3/2}}
\le {||B||_{2\to2}\over2\sqrt d}.                    \tag{A.4}
```

This proves the upper and absolute bounds, including when the un-absolute
maximum happens to be negative at the initial level.

There is one important scope distinction.  A Boolean eigenvector with
eigenvalue `+sqrt(h)` is an explicit hypothesis, not a consequence of
symmetry alone.  It also forces `sqrt(h)` to be an integer, since `Hu` is
integral.  Thus no hidden nonsquare order is being used.  If only

```math
Hu=-(\sqrt h)u
```

is available, (A.1) reverses the signed value: the absolute sequence still
embeds monotonically, while the stated one-step monotonicity for `q_r^+`
does not follow.  Replacing `H` by `-H` restores the positive case, but that
is a different generator.  This should be mentioned because the word
"regular" is sometimes used without fixing the sign of the row sum.

## 2. Compact response-carrier argument

For every pair of Boolean blocks,

```math
|x_i^TH_rx_j|
\le ||H_r||_{2\to2}||x_i||_2||x_j||_2
=s_r^{3/2},                                            \tag{A.5}
```

so `K_r^(d)` really lies in the stated symmetric cube.  Amplifying *all*
blocks by the same `u` gives

```math
{(x_i\otimes u)^T(H_r\otimes H)(x_j\otimes u)
 \over (s_rh)^{3/2}}
= {x_i^TH_rx_j\over s_r^{3/2}},                        \tag{A.6}
```

including for `i != j`.  This proves nesting before convexification and
therefore after convexification.

The remaining topological claims need no extra realizability assumption.
An increasing family of nonempty compact subsets of the fixed compact cube
converges in Hausdorff distance to the closure of its union: a finite cover
of that closure by balls centered in the union is eventually contained in
one common level.  The nested union is convex, so its closure is compact and
convex.  Strict convex-geometric terminology would call it a compact convex
**set** rather than necessarily a "body," because nonempty interior in all
`d(d+1)/2` coordinates is not proved (and fails for degenerate generators
such as the order-one example).

For a block vector `x=(x_1,...,x_d)`, the identity

```math
{x^T(B\otimes H_r)x\over2(ds_r)^{3/2}}
={1\over2d^{3/2}}\sum_{i,j}B_{ij}
  {x_i^TH_rx_j\over s_r^{3/2}}                         \tag{A.7}
```

has the correct off-diagonal multiplicity: both sides sum over ordered
pairs.  Linear maximization is unchanged by convexification.  For the
absolute functional, the same fact follows from

```math
\left|L\left(\sum_t\lambda_tK_t\right)\right|
\le\sum_t\lambda_t|L(K_t)|\le\max_t|L(K_t)|.          \tag{A.8}
```

Hausdorff convergence therefore gives both support-function limits.  In
finite dimension it is uniform on every coefficient-bounded family because
support functions are Lipschitz with constant given by the corresponding
dual norm of `B`.

Finally, a cube grid gives an external `l_infinity` epsilon-net with at most

```math
(1+2/\epsilon)^{d(d+1)/2}                              \tag{A.9}
```

points (the bound is deliberately loose).  This is an external net; an
internal-net statement with exactly the same constant has not been proved.
To turn (A.9) into a numerical response-error bound one must also multiply
by the `l_1` size of the queried coefficient matrix, but the draft only
claims metric entropy of the carrier itself.

## 3. Exact boundary of the Walsh tensor identity

The asserted identity is valid for the standard linear-label Walsh graph
language, but the proof should be displayed rather than left implicit.  Let
`E_m=F_2^m direct-sum F_2^m`.  After grouping old and new coordinates,

```math
W_{E_{m_0+r}}=W_{E_{m_0}}\otimes W_{E_r}
=W_{E_{m_0}}\otimes W_4^{\otimes r}.                  \tag{A.10}
```

For a label extended as `(a,0^r)`, its modulation is

```math
D_{(a,0^r)}=D_a\otimes I_{4^r}.                        \tag{A.11}
```

Thus every child block satisfies

```math
C_{(a,0^r)}^{(m_0+r)}
=C_a^{(m_0)}\otimes W_4^{\otimes r},                  \tag{A.12}
```

and the common Walsh bridge satisfies the same identity.  If the graph,
scalar onsite coefficients, and scalar bridge coefficients are frozen,
reordering block and Walsh coordinates gives

```math
M_{m_0+r}=M_{m_0}\otimes W_4^{\otimes r}.             \tag{A.13}
```

The corollary therefore follows from HA.1--HA.2 with
`B=M_(m_0)` and `d=dim M_(m_0)`.

The declaration matters.  Equation (A.13) need not hold for a newly
appended ambient-coordinate label, a coordinate-dependent onsite field, a
bridge whose microscopic coefficients change with `r`, or an external pole
that singles out the new coordinates.  "Scalar onsite weights" and
zero-extension of every exposed label are sufficient and should remain in
the formal statement.

The order-four generator in the draft is correct:

```math
W_4(1,1,1,-1)^T=2(1,1,1,-1)^T.                       \tag{A.14}
```

Hence it supplies the required *positive* regular Boolean eigenvector.

## 4. Hollow sign corollary

Suppose `B` is symmetric with every entry in `{+-1}` and `tr(B)=0`.  Then
`B tensor W_4^(tensor r)` has sign entries, and hollowing changes a Boolean
quadratic by the constant

```math
x^T(B\otimes W_4^{\otimes r})x-x^TA_rx
=tr(B\otimes W_4^{\otimes r})
=tr(B)tr(W_4)^r=0.                                    \tag{A.15}
```

Thus `A_r` is a valid symmetric hollow signing and its full and hollow
quadratic energies agree pointwise, not merely at the optimum.  HA.1 then
proves convergence of

```math
Q(A_r)/(d4^r)^{3/2}.                                  \tag{A.16}
```

For this particular generator `tr(W_4)=0`, so the condition `tr(B)=0` is
actually stronger than necessary once `r>=1`; it is a clean sufficient
condition that also covers the base level `r=0`.  Complete signed Walsh
block graphs with signed scalar onsite and bridge coefficients meet the
entrywise-sign condition.  Programs with missing graph edges have zero
off-diagonal blocks and are not complete signings, although their response
limits remain covered by HA.1.

## 5. Verifier audit

I reran

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_amplification_limit.py
```

and obtained

```text
regular-Hadamard amplification checks passed: 127
```

The script correctly checks:

- the positive regular eigenvector of `W_4`;
- the exact factor `8=4^(3/2)` under one amplification;
- upper and absolute monotonicity on three small outer matrices;
- the operator-norm ceiling with the factor `1/2`;
- and pointwise energy invariance after hollowing one zero-trace sign
  template.

It is a useful finite diagnostic, but its comments slightly exceed its
coverage.  It checks whole quadratic values for random witnesses, not every
pair entry of `K_r`; it does not check (A.13), the Hausdorff claim, or the
stated strict wind-tunnel transition.  I independently enumerated the
`m=1` two-block program and confirmed upper energy `12`.  The `m=2` value
`104` is consistent with the repository's existing Walsh response data, but
the accompanying amplification verifier never constructs that program.
Before promotion, either add an exact specialized check returning
`12 -> 104` or cite the existing exact computation that certifies `104`.
The limit proof itself does not depend on strictness.

The random tests of (A.15) are harmless but unnecessary: once the template
is checked to be sign-valued and trace zero, (A.15) proves the identity for
all Boolean vectors.

## 6. Novelty and near-original significance

The theorem is rigorous and generative in the limited sense requested by
the theory program:

- it gives an exact scale-preserving recovery map rather than a leading-loss
  bridge estimate;
- it yields a fixed-dimensional carrier complete for all fixed outer
  quadratic queries;
- and it proves a whole-sequence limit on a nontrivial dense hollow-sign
  hierarchy, not just convergence of a selected numerical subsequence.

Its mathematical engine is nevertheless elementary tensor
supermultiplicativity plus finite-dimensional compactness.  It neither
computes `K_infinity^(d)` nor gives an effective finite presentation of it.
It covers the geometric orders `d4^r` of one fixed tensor-generated family;
it does not optimize over that family, transfer near-minimizers between
orders, tolerate nonsummable perturbations, realize arbitrary dense bridges,
or compare its limit with `M_n/n^(3/2)`.

Accordingly, **"near-original structured benchmark" is accurate**, while
"progress on the original convergence theorem" would be too strong.  The
first genuinely new extension would be the one already identified in the
draft: prove that an approximately regular, nonexact amplification gives
Hausdorff inclusions with summable error, or identify a larger family closed
under an equally lossless recovery map.
