# Adversarial self-audit of `gauge_reset_necessity_report.md`

**Verdict.** The algebraic statements and constants in the report survive
rechecking, subject to two explicit qualifications. Theorem 3.1 is correct
for `r>=2`, a factorial declared language, and genuinely arbitrary
per-step projective disturbances. The `2 by 2` clamp is a decisive
counterexample to the *narrow* completeness claim “entrywise kernel gauge or
small full-image reset,” but it is **not** a decisive counterexample to a
literal, maximally broad phrase “zero-holonomy recurrent dynamics plus
resets”: after one transient step its paired dynamics are recurrent with zero
increment. Finite-semigroup absorption should therefore be presented as a
finite algebraic realization/refinement of zero-holonomy recurrence, or as a
third mechanism relative to the repository's currently proved
entrywise/additive criteria, not as a refutation of every possible broad
gauge--reset decomposition.

## 1. Audit of Theorem 3.1

### 1.1 Product orientation and the suffix argument

The recursion is

```math
e_t=P_t e_{t-1}+\eta_t.
```

Thus

```math
e_T=P_T\cdots P_1e_0+
\sum_{s=1}^T P_TP_{T-1}\cdots P_{s+1}\eta_s.       \tag{A.1}
```

This matches the report's orientation. Suppose the last length-`L` window
contains a reset factor occupying times `a,...,b`, so

```math
P_bP_{b-1}\cdots P_a=0
```

on projective space. Every disturbance with `s<a` traverses the whole factor
`P_b...P_a` in (A.1), and is therefore killed. A disturbance inserted at
time `a` occurs *after* `P_a`, so it need not be killed; this is why the
correct remaining count is `T-a+1`, not merely the tail after `b`. Since
`a>=T-L+1`, at most `L` disturbances remain. The upper bound `L epsilon` is
therefore correct.

The factorial-language hypothesis is doing real work: it guarantees that
the final window and every suffix used below are allowed words. For a
non-factorial language, the condition should instead be imposed directly on
the prefix closure / actual allowed trajectories.

### 1.2 Adversarial lower-bound construction

Let a length-`T` word contain no zero-product factor and put

```math
A_s=P_T\cdots P_{s+1}.
```

For `s<T`, `A_s` is the selector product of a nonempty contiguous factor,
so it is nonzero. For `s=T`, it is the identity, also nonzero provided
`r>=2`. Each nonzero selector has two output coordinates `(j_s,k_s)` mapped
to distinct input coordinates. There are at most `r(r-1)` ordered output
pairs, hence a common pair `(j,k)` occurs for at least
`ceil(T/[r(r-1)])` indices (the report's weaker floor is valid).

For each chosen `s`, if `A_s` sends `j,k` to `alpha_s,beta_s`, choose

```math
\eta_s(\alpha_s)=+\epsilon,
\qquad
\eta_s(\beta_s)=-\epsilon,
```

and zero elsewhere. Then

```math
{1\over2}\operatorname {osc}(\eta_s)=\epsilon,
\qquad
(A_s\eta_s)_j-(A_s\eta_s)_k=2\epsilon.
```

All chosen summands have the same sign on the *same final coordinate pair*,
so no cancellation is possible there. The final Hilbert norm is at least
the number of chosen indices times `epsilon`. This verifies (3.4).

### 1.3 Exact scope of the equivalence

The theorem should explicitly assume `r>=2`. Its equivalence is quantitative
but not constant-for-constant:

- reset gap `L` implies robust constant at most `L`;
- robust constant `C` implies no reset-free word of length exceeding roughly
  `C r(r-1)`.

This is exactly what “equivalent up to the displayed constants” says. It is
not an equivalence between one fixed `L` and the same robustness constant.

### 1.4 Gauge corollary

For

```math
e_t=P_te_{t-1}+h_t-P_th_{t-1}+\eta_t,
```

the substitution `\bar e_t=e_t-h_t` indeed gives

```math
\bar e_t=P_t\bar e_{t-1}+\eta_t.                 \tag{A.2}
```

One small qualification is needed: `\bar e_0=e_0-h_0` is generally not
zero. A completed tangent reset kills this initial term, so after the first
reset the report's bound is valid. Before the first reset one must also pay
`\|P_t\cdots P_1(e_0-h_0)\|_H`. At the endpoint, a bound on `e_t` rather
than `\bar e_t` also pays `\|h_t\|_H`. This is consistent with the report's
phrase “bounded endpoint gauges,” but an integrated theorem should display
the endpoint term explicitly.

**Conclusion on Theorem 3.1:** verified, with the harmless missing
qualification `r>=2` and the just-noted initial/terminal gauge terms.

## 2. Audit of the max-plus clamp calculations

The convention is

```math
(F_Su)_b=\max_a(u_a+S_{ab}),
\qquad z=u_2-u_1,
\qquad z'=(F_Su)_2-(F_Su)_1.                     \tag{A.3}
```

### 2.1 Exact clamp

For

```math
S_0=\begin{pmatrix}0&0\\-1&0\end{pmatrix},
```

one gets

```math
z'=\max(0,z)-\max(0,z-1)=\operatorname {clip}(z,0,1).
```

For

```math
S_\delta=\begin{pmatrix}0&\delta\\-1&0\end{pmatrix},
```

one gets

```math
z'=\max(\delta,z)-\max(0,z-1)
  =\operatorname {clip}(z,\delta,1).
```

Both formulas, including all breakpoint cases, are correct. Both maps are
idempotent.

For representatives `(0,z)` and `(0,z')`,

```math
d_H([u],[v])={1\over2}\operatorname {osc}(u-v)
            ={1\over2}|z-z'|.                   \tag{A.4}
```

Hence

```math
\sup_z d_H(P_0z,P_\delta z)=\delta/2,
```

and the image diameters are exactly `1/2` and `(1-delta)/2`. The report's
half-Hilbert constants are correct.

The difference matrix is

```math
E=S_\delta-S_0=
\begin{pmatrix}0&\delta\\0&0\end{pmatrix},
```

whose alternating rectangle is `-delta`; it is not row-plus-column
separable. This check is also correct.

### 2.2 Drifting comparator

For

```math
\widehat S_\delta=
\begin{pmatrix}0&0\\-1+\delta&\delta\end{pmatrix},
```

equation (A.3) yields

```math
z'=\max(0,z+\delta)-\max(0,z-1+\delta)
  =\operatorname {clip}(z+\delta,0,1).
```

Thus `\widehat P_delta^t(0)=min(t delta,1)` in the `z` coordinate, and its
Hilbert distance from `P_0^t(0)` is `min(t delta,1)/2`. The report says that
the coordinate error reaches order one, not that the Hilbert distance equals
one, so there is no constant error.

Moreover

```math
\widehat S_\delta-S_0=
\begin{pmatrix}0&0\\\delta&\delta\end{pmatrix}
```

is a row potential. Repeating it fails adjacent-interface compatibility,
exactly as claimed.

**Conclusion on the clamp computations:** fully verified. The formulas were
also evaluated numerically at all three affine regions for several values of
`delta`.

## 3. Audit of finite-semigroup absorption

### 3.1 Proof

If both actions factor exactly through the same finite semigroup `S`, a word
`w` and a shortest representative `v` of its element give exactly the same
map under each action:

```math
F_w=F_v,
\qquad G_w=G_v.
```

Hybridizing the at most `L` factors of `v` costs at most `epsilon` per
factor because all prefixes supply the same input to the factor being
replaced and all suffixes are nonexpansive. The `L epsilon` bound is correct.
No hidden commutativity or invertibility assumption is used.

The clamp pair supplies two actions of the monoid `{1,p}` with `p^2=p`, so
the theorem applies with `L=1`.

### 3.2 Is this genuinely distinct from the repository's mechanisms?

The strongest objection is valid in part.

- It is **genuinely distinct from the proved entrywise gauge criterion**:
  the nonzero rectangle in (2.6) is decisive.
- It is **genuinely distinct from a small full-image reset** at the target
  `O(delta)` scale: every positive power has image diameter asymptotic to
  `1/2`.
- It is **not necessarily distinct from the broad English phrase
  “zero-holonomy recurrent dynamics.”** In the clamp example, after one
  transition the paired point `(P_0z,P_delta z)` is fixed. The incremental
  error on that recurrent paired orbit is zero. Thus a decomposition which
  is allowed an arbitrary finite transient and a state-dependent paired
  recurrence classifies the example as zero-holonomy recurrence, without a
  reset.

The current repository theorem about repeatable holonomy is narrower: it
uses additive labels on a finite context graph, while the clamp's relevant
state is the exact semigroup action / paired orbit and is not recovered from
the kernel rectangle labels. Consequently the honest claim is:

> Finite-semigroup absorption is a new **finite certificate and algebraic
> realization** of stable recurrence relative to the current proved
> criteria. It refutes completeness of “kernel gauge + small full-image
> reset,” but by itself does not refute a decomposition whose
> “zero-holonomy recurrent” clause is defined broadly enough to include all
> exact semigroup relations.

This downgrade should be made explicit in any repository integration. Calling
it an unconditional “third mechanism beyond zero holonomy” would overstate
the result.

### 3.3 Strength and novelty of the hypothesis

Requiring both families to factor through the same finite semigroup is
strong. It does not follow from one-step closeness. It is nevertheless a
checkable structural hypothesis: an idempotent, band, or finite automaton
transition algebra supplies it. In the current theory's language it is close
to **quotient compatibility / bounded query-feature algebra**. The theorem
is a quantitative consequence (`L epsilon` uniformly over all words), not a
new universal abstraction.

Thus it is useful benchmark progress, but should not be advertised as a
classification of all stable coherent perturbations. Near conjugacies,
invariant graphs, and infinite compact transition semigroups remain outside
it.

## 4. Audit of the affine-selector cycle theorem

For `A(e)=P_sigma e+b`,

```math
A^k(e)=P_\sigma^ke+\sum_{t=0}^{k-1}P_\sigma^tb.  \tag{A.5}
```

Coordinate `j` samples `b` along the forward orbit
`j,sigma(j),sigma^2(j),...`. Two basins have linearly separating sums exactly
when their recurrent cycles have different means. If all cycle means equal
`beta`, the equations

```math
p_j-p_{\sigma(j)}=b_j-\beta
```

are consistent on every cycle and then recursively on every in-tree. Hence

```math
b=p-P_\sigma p+\beta1,
\qquad
A^k(e)=p+P_\sigma^k(e-p)+k\beta1.
```

All signs and orientations are correct.

One wording should be sharpened. If cycles `C,C'` have means
`beta_C,beta_C'`, the coordinate difference grows as
`k(beta_C-beta_C')+O(1)`, while the Hilbert norm is at least
`k|beta_C-beta_C'|/2-O(1)`. Saying it grows “linearly at the difference” is
qualitatively correct, but an exact projective constant includes the factor
`1/2` and may involve the range over all cycle means.

This theorem is stronger and cleaner than the clamp as evidence for the
needed refinement: selector transport converts ordinary holonomy into a
twisted cocycle, and transient selector trees provide **relative/tangent**
resets without collapsing the entire state image.

## 5. Strongest defensible synthesis after audit

The report's conclusions should be ranked as follows.

1. **Strong theorem:** Theorem 3.1 is a sharp gauge--tangent-reset converse
   for adversarial residuals on a fixed finite selector language.
2. **Strong theorem:** Theorem 4.1 is the exact one-cycle twisted-holonomy
   criterion.
3. **Valid new sufficient mechanism:** common finite-semigroup factorization
   gives a depth-uniform bounded-normal-form estimate.
4. **Sharp falsifier, with limited target:** the clamp disproves completeness
   of entrywise gauge plus small full-image reset for coherent exact maps.
5. **Not established:** a counterexample to every broadly formulated
   zero-holonomy-recurrence plus reset decomposition for nonlinear max-plus
   systems.

The most reliable revision to the current theory is therefore not simply to
add an unrelated third box. It is to refine the two existing boxes:

```math
\begin{array}{c}
\text{ordinary endpoint gauge}
\longrightarrow
\text{selector-transported (twisted) cocycle},\\[2mm]
\text{small full image}
\longrightarrow
\text{small/zero image on the residual-error bundle},
\end{array}                                       \tag{A.6}
```

and then record finite-semigroup relations as a checkable way for a coherent
pair to have zero residual on its recurrent action algebra.

## 6. Recommended integration language

A theorem-level integration can safely say:

> On a fixed coordinate-selector language, after subtracting a transported
> endpoint gauge, uniform robustness against arbitrary local residuals is
> equivalent (up to an `r(r-1)` constant) to syndetic rank-one reset of the
> residual directions. For a recurrent affine selector, the exact obstruction
> is unequal cycle mean; equality is equivalent to a twisted coboundary.

A counterexample-level integration can safely say:

> For coherent fixed continuations, arbitrary residual robustness is too
> strong a model. Nearby idempotent max-plus clamps stay `O(delta)`-close at
> all depths despite nonzero kernel rectangle defect and no `O(delta)`
> full-image reset. More generally, exact actions of a common finite
> semigroup are stable with constant equal to its bounded normal-form
> diameter.

It should **not** say without qualification:

> Gauge and reset are false as a complete decomposition even when
> zero-holonomy recurrence is allowed in its broadest possible sense.

That stronger statement is not proved by the clamp.
