# Independent audit: relative-involution composition stability

Audited draft:
[`relative_involution_composition_stability.md`](relative_involution_composition_stability.md).

## Verdict

**PROMOTE WITH REPAIR.**  The displayed theorem is correct as written except
for one edge-case overstatement: the positivity claim following (RIC.8) needs
`r>0`.  All scale factors in the block energy, square identity, Rayleigh
ceiling, bipartite construction, exact commuting equality, and regular-graph
comparison check out.

There is, however, an exact identity that makes the approximate-commuting
half of the theorem substantially stronger than the draft states.  If
`t=Fs`, then

```math
N-t^TCt={1\over2}\|[C,F]s\|_2^2.                         \tag{A.1}
```

Consequently the loss is quadratic in the commutator, not linear.  This
should be incorporated before promotion.  The result is a genuine
conditional robustification of the Walsh commute/anticommute argument, but
not a robust realization theorem, a response quotient, or a proof that the
two operator norms form a sufficient state.

## 1. Block-energy normalization

Put

```math
M=I_k\otimes C+A_G\otimes F.
```

Because `A_G` is the adjacency matrix of a simple undirected graph,

```math
X^TMX
=\sum_i x_i^TCx_i
2\sum_{\{i,j\}\in E(G)}x_i^TFx_j.
```

Thus multiplication by `lambda/2` gives exactly

```math
{\lambda\over2}\sum_i x_i^TCx_i
+\lambda\sum_{\{i,j\}\in E(G)}x_i^TFx_j.
```

There is no missing factor of two in (RIC.2), (RIC.6), or (RIC.7).
In the Walsh specialization, `lambda=q=sqrt(N)`, so the saturated child and
edge values are respectively `N^(3/2)/2` and `N^(3/2)`, as in Theorem 21.10.

## 2. Squaring and the operator ceiling

Using `C^2=F^2=I` gives

```math
M^2=(I_k+A_G^2)\otimes I_N
    +A_G\otimes(CF+FC),
```

so (RIC.9) is exact.  Since `M` is symmetric,
`||M||^2=||M^2||`; moreover

```math
\|I_k+A_G^2\|=1+\rho_G^2,
\qquad
\|A_G\otimes(CF+FC)\|=\rho_G\eta_+.
```

The triangle inequality therefore proves (RIC.10), and
`||X||_2^2=kN` then gives exactly the coefficient in (RIC.4).

There is a useful stronger spectral identity.  The two summands of `M^2`
commute, so if `S=CF+FC`, then

```math
\|M\|^2
=\max_{a\in\operatorname{spec}(A_G),\,\sigma\in\operatorname{spec}(S)}
  (1+a^2+a\sigma).                                      \tag{A.2}
```

In particular, when `G` is `r`-regular and bipartite, both `r` and `-r`
belong to the adjacency spectrum.  Pairing the appropriate one with an
eigenvalue of `S` of absolute value `eta_+` shows

```math
\|M\|^2=1+r^2+r\eta_+.                                  \tag{A.3}
```

Thus the anticommutator ceiling used in (RIC.8) has no hidden spectral
triangle loss in the regular bipartite case.  It may still have the intended
Boolean Rayleigh loss.

## 3. Bipartite Boolean section and the missing quadratic identity

Let `t=Fs`.  The hypotheses imply `s,t` are Boolean, `Cs=s`, and

```math
s^TFt=s^TF^2s=N.
```

Hence assigning `s` and `t` to opposite color classes saturates every edge,
independently of orientation.  Also

```math
Ct-t=(CF-FC)s=[C,F]s.                                   \tag{A.4}
```

The draft bounds `t^TCt` by applying Cauchy--Schwarz to this difference.
That proves its stated, valid estimate

```math
t^TCt\ge N(1-\eta_-).
```

But orthogonality of the involution `C` gives the exact identity

```math
\|Ct-t\|_2^2
=2N-2t^TCt,
```

which combined with (A.4) proves (A.1).  Define the sectionwise defect

```math
\delta_s={\|[C,F]s\|_2\over\sqrt N}\le\eta_-.
```

If `b` blocks receive `t`, the precise constructed energy is

```math
\lambda N\left({k\over2}+|E(G)|-{\delta_s^2\over4}b\right).
                                                                    \tag{A.5}
```

For a fixed global bipartition one may take
`b=min(|L|,|R|)`.  More sharply, connected components may be oriented
independently, giving

```math
b=\sum_{H\in\operatorname{cc}(G)}
  \min(|L\cap H|,|R\cap H|).                            \tag{A.6}
```

Thus the recommended replacement for (RIC.6) is

```math
\max_XE_G(X)
\ge\lambda N\left\{
 {k\over2}+|E(G)|
 -{\eta_-^2\over4}
  \sum_H\min(|L\cap H|,|R\cap H|)
\right\}.                                               \tag{A.7}
```

The original linear-loss formula follows because `eta_-<=2` for two
orthogonal involutions and hence `eta_-^2/4<=eta_-/2`.  Therefore this is a
strengthening, not a repair of a false inequality.

The assumption `Fs in {+-1}^N` is doing indispensable work.  A small
operator perturbation of `F` need not preserve it.  Accordingly (A.7) is
robust in the **relative algebra conditional on an exact Boolean transported
pole**; it is not yet robustness of the Walsh construction under arbitrary
matrix perturbations.

## 4. Exact commuting equality

When `[C,F]=0`, (A.4) says `Ct=t`.  The bipartite assignment therefore
saturates all terms.  Conversely, for Boolean vectors,

```math
x^TCx\le\|x\|_2^2=N,
\qquad
x^TFy\le\|x\|_2\|y\|_2=N.
```

Summing these termwise upper bounds proves (RIC.7), with precisely the stated
factor.  In fact the same triangle estimates show that the equality also
holds if the objective is `max_X |E_G(X)|`, a potentially useful connection
to the motivating absolute quadratic cap.

No claim of equality is justified for small nonzero commutator: (A.7) is
only a constructed lower section, as the draft correctly states.

## 5. The regular comparison and positivity threshold

For an `r`-regular graph,

```math
\rho_G=r,
\qquad
|E(G)|={kr\over2}.
```

The commuting pair therefore has value

```math
{\lambda kN\over2}(1+r),
```

while (RIC.4) bounds the second pair by

```math
{\lambda kN\over2}\sqrt{1+r^2+r\eta}.
```

Their difference is exactly (RIC.8).  For `r>0`, its bracket has the useful
rationalized form

```math
1+r-\sqrt{1+r^2+r\eta}
={r(2-\eta)\over
  1+r+\sqrt{1+r^2+r\eta}}.                              \tag{A.8}
```

It is strictly positive exactly when `eta<2`.  For `r=0` (the edgeless
graph), it is zero even when `eta<2`; the draft must add `r>0` to its
positivity sentence.  At `eta=2` it is also zero, and `eta>2` is impossible
because `||CF+FC||<=2`.

The quantity displayed in (RIC.8) is a **raw extensive energy gap**, not a
normalized gap.  Dividing by `kN` gives

```math
{\lambda\over2}
\left[1+r-\sqrt{1+r^2+r\eta}\right].
```

The terminology should be changed or the normalization explicitly stated.

## 6. Novelty and exact scope

Relative to the committed Walsh theorem, this draft does more than merely
rename the parity bit:

1. it isolates the algebraic proof for arbitrary real symmetric involutions;
2. it allows every bipartite composition graph;
3. it gives a continuous anticommutator ceiling; and
4. after (A.1), it gives a quadratic stability law for a transported Boolean
   maximizing section.

Those are genuine reusable statements.  Substituting `eta_-=0` for the even
Walsh class and `eta_+=0` for the odd class recovers the mechanism of Theorem
21.10 exactly.  The new content is the quantitative neighborhood of those
two algebraic poles.

The result nevertheless has deliberately limited theoretical reach:

- it constructs no non-Walsh family with an exact Boolean transported pole;
- it does not make that pole stable under arbitrary perturbations;
- it gives certificates, not a finite response update or compositional
  quotient;
- the two norms do not determine the intermediate Boolean optimum;
- it proves neither necessity nor sufficiency of commutator data under
  varying operator words; and
- the general comparison does not require the two pairs to have matching
  isolated spectra or matching static response summaries.  Only the Walsh
  specialization currently supplies that controlled contextual comparison.

Therefore the interpretation sentence that “isolated spectra ... are
identical in the two cases” must explicitly be restricted to the Walsh
example (or matching-spectrum hypotheses must be added).  It is not a
consequence of RIC.1 for arbitrary pairs.  Likewise, a pair of families with
the same two norm values but separated responses would demonstrate that the
proposed data are incomplete; it would not falsify any theorem presently
claimed.  The final “falsifier” sentence should be relabeled as an
incompleteness test.

## 7. Required repairs before canonical promotion

1. Replace the linear commutator loss by the exact sectionwise identity
   (A.1) and the quadratic bound (A.7).
2. Add `r>0` to the strict-positivity assertion after (RIC.8).
3. Rename (RIC.8) an extensive energy gap, or display its normalization.
4. Restrict the identical-isolated-spectrum interpretation to the Walsh
   specialization.
5. Describe the Boolean-section assumption as brittle under generic matrix
   perturbation; do not market the result as full perturbative robustness.
6. Relabel the proposed same-state separated pair as a test of state
   incompleteness, not a falsifier of RIC.1.

With those changes, this is theorem-level progress: a concise robust spectral
certificate that cleanly abstracts and quantitatively strengthens the Walsh
holonomy mechanism, while remaining honest that it is not yet a dynamic
compression theorem.
