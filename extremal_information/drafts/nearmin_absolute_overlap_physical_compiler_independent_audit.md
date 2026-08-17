# Independent audit: absolute-overlap physical compiler

Status: **PASS after repair.**

Recheck note: all four repairs requested in the first audit have been applied.
The repaired verifier passes `408996` checks.  The added two-cap dichotomy
(AO.21a)--(AO.21c) was also checked: its sign alignment follows from Hamming
triangle inequalities and its final strict inequality follows correctly from
(AO.20).

This audit was performed independently against
`nearmin_absolute_overlap_physical_compiler.md` and
`verify_nearmin_absolute_overlap_compiler.py`.  It checks the exact
normalizations, concentration scale, free-shore realization, and the
unconditional corollary.  It does not assert that exact-minimizer shells have
a growing projective packing.

## Audited derivation

Let `E=n(n-1)/2`.  For augmented cuts,

```math
\langle \sigma c(u),\tau c(v)\rangle
=\sigma\tau\frac{\langle u,v\rangle^2-n}{2}.
```

Consequently

```math
\frac{|\langle z^u,z^v\rangle|}{E}\le 1-\gamma
\quad\Longrightarrow\quad
1-\frac{\langle u,v\rangle^2}{n^2}
\ge \gamma(1-1/n).
```

The constant `1/4` in the spherical lemma is valid.  Indeed the loss from
the two maxima is at least

```math
\frac12\min\{\alpha,\lambda\}
\bigl[(1-|V\cdot Y|^2)+(1-|U\cdot Y|^2)\bigr],
```

and the bracket is at least `1-|U dot V|`, hence at least `theta/2`.

For a sparse-flip child the coordinate identity

```math
\mathbb E b^u=(1-p)a+pz^u
```

is exact.  The variance proxy is `O(pE)=O(alpha n^(3/2))` and the summands
are bounded.  Bernstein at deviation

```math
C(\sqrt\alpha n^{5/4}+n)
```

has exponent `Omega(n)`, uniformly in the augmented cut, so a sufficiently
large absolute `C` pays for `exp(O(n))` children times `exp(O(n))` cuts.
For the simultaneous edit-count event one additionally uses
`pE=Theta(alpha n^(3/2))`.  It pays for exponentially many children when
`alpha` is fixed, and more generally when `alpha sqrt(n) -> infinity`.

The free-shore identity is also correct despite the absolute cap.  With the
shore Hamiltonian omitted, choosing the all-equal shore endpoint with its
sign matched to the outer absolute value gives exactly

```math
\max_{y,\tau}\{\tau H_b(y)+h|u\cdot y|\}.
```

Adding a common shore signing perturbs this maximum by at most `Q(C_h)` in
each direction.  Thus comparisons cost `2Q(C_h)`, not a hidden endpoint or
pinning term.

For the cross child, the valid one-sided estimate is

```math
\langle z^v,\tau c(y)\rangle
\le |\langle z^v,\tau c(y)\rangle|
\le \frac{\langle v,y\rangle^2+n}{2}.
```

After spherical relaxation, subtraction from the target lower bound leaves

```math
pE-\frac{pn}{2}-\frac\alpha2n^{3/2}=-\alpha n^{1/2}.
```

This verifies (AO.13), including its `alpha/n` coefficient.  The parent cap
is at most

```math
Q(a)+2\alpha n^{3/2}+nh+Q(C_h)=O(n^{3/2})
```

under the stated bounded-cap hypothesis.  No separately paid scalar
channels, full pinning, or unbalanced shore endpoint are hidden here: the
single rank-one shore is optimized jointly with the child before the
absolute value is taken.

## Unconditional corollary

The mathematical conclusion of Corollary AO.3 is valid.  Choose `v` with
`|u dot v| <= 1`.  Then

```math
\frac{|\langle z^u,z^v\rangle|}{E}
=\frac{|\langle u,v\rangle^2-n|}{n(n-1)}
\le\frac1{n-1}.
```

Only the target `z^u` is used in the lower bound; the cross-child upper bound
uses no energy condition on `z^v`.  Repeating the ordered-pair proof with
`d=0` therefore gives the advertised gap.  This is a common all-spins-free
context: the two children are fixed before evaluating the same rank-one
query.  It proves one physical response bit, not a reusable growing family.

The same conclusion has a useful vanishing-strength form.  If

```math
\alpha=\alpha_n\to0,
\qquad \alpha_n\sqrt n\to\infty,
\qquad h=\lfloor\alpha_n\sqrt n\rfloor,
```

then the simultaneous concentration and edit-count union bounds remain
valid, and

```math
Q(P^{u|u})-Q(P^{v|u})
\ge \bigl(\alpha_n/4-o(\alpha_n)\bigr)n^{3/2}.
```

Indeed

```math
\frac{\rho_n}{\alpha_n n^{3/2}}
=O\left((\alpha_n\sqrt n)^{-1/2}
       +(\alpha_n\sqrt n)^{-1}\right)=o(1),
```

while `Q(C_h)/(alpha_n n^(3/2))=O(alpha_n/sqrt(n))=o(1)` and the floor in
`h` is relatively negligible.  The threshold `alpha_n >> n^(-1/2)` is
therefore the correct one for this proof.

## Repairs applied and rechecked

1. **Parameter regime for the exponential simultaneous
   edit-count event.**  The text's fixed-`alpha` theorem is fine.  If `alpha`
   is allowed to depend on `n`, the assertion that (AO.7) can be intersected
   for all exponentially many children requires at least
   `alpha sqrt(n) -> infinity` (or a correspondingly smaller child family).
   This qualification and the vanishing-strength corollary are now present.

2. **Absolute value in (AO.19).**  At finite orders the cross-pair
   overlap is

   ```math
   \left|1-\frac{2k(n-k)}E\right|,
   ```

   not always `1-2k(n-k)/E`; the latter is eventually positive and has the
   stated asymptotic, but is negative for some small `n` (for example
   `n=10,k=4`).  The formula is now corrected.

3. **One-target extension.**  Formally Theorem AO.2 asks
   every member of `U` to be near-top, whereas Corollary AO.3 uses an
   arbitrary decoy `v`.  The proof already establishes the stronger
   one-target/decoy statement.  State it as a sentence in the theorem or say
   in the corollary that its proof repeats (AO.15)--(AO.16), rather than
   literally invoking all hypotheses of AO.2.  The theorem now states this
   extension explicitly.

4. **Hardened overlap verifier.**  In
   `check_edge_overlap_identity`, the edge vectors are `int8`, so NumPy's
   matrix product can overflow once `E>127`.  Cast both operands to `int64`
   before `@`, and add deterministic aligned, orthogonal, and one-flip pairs;
   random pairs almost never exercise the overflow.  The current successful
   run (`408930` checks) supports the other identities but does not cure this
   test weakness.  The operands are now cast to `int64`, and aligned,
   one-flip, and balanced deterministic cases have been added.

With these applied repairs, the verdict is **PASS** for AO.1--AO.16,
AO.20--AO.21, the physical/contextual interpretation, and the vanishing
`alpha_n >> n^(-1/2)` extension.  The projective-packing premise remains a
conditional structural obligation, exactly as the draft says.
