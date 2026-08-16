# Independent audit: response geometry and robust tropical rank

**Status.** Independent verifier report.  This audits Sections 2--6 of
`phase2_response_geometry_composition.md` and all of
`phase2_robust_tropical_rank.md`.  It does not edit a surface theory file.

## 1. Verdict

All six mathematical cores requested for audit are correct:

1. the response-separation polytope has the stated exact orthogonal-product
   and scaling law;
2. identical child polytopes can give either total same-space cancellation or
   an orthogonal box at the parent;
3. every latent coordinate of a universal max-plus product is killed by an
   exact one-bit collision once there are at least two factors, including a
   coordinate in the first or last factor;
4. a four-cell gap `G` forces min-plus factor rank at least the fooling-set
   size under strict uniform error below `G/4`, and `G/4` is sharp in general;
5. the zero-diagonal/one-off-diagonal family has exact rank `r` but rank-one
   normalized mean-square error `1/r`; and
6. Sheshadri's transversal block therefore makes the optimal approximate
   min-plus rank of a linear-code conditional distance table exactly `2^s`
   for raw uniform additive error strictly below `1/2`.

There are no constant or sign errors in the proofs.  Three qualifications
should accompany promotion:

- “Every approximation retains the exact `2^s` count” should be replaced by
  “no approximation can use fewer than `2^s` terms; the optimal approximate
  rank remains `2^s`.”  An individual approximating matrix can have rank
  greater than `2^s`.
- The code result is robust at the **unscaled Hamming lattice scale**.  If the
  landscape is normalized by `m`, its hypothesis is error `<1/(2m)`, not a
  fixed normalized error.  It gives no control at error `epsilon m`, at
  relative error, or under average error.
- The exact source is a July 2026 arXiv preprint, not a peer-reviewed
  publication.  Its Remark 7(ii) publicly asks the approximate-factorization
  question and does not state this corollary.  A targeted search found no
  prior statement of the precise robust code corollary, but the stability
  argument is elementary enough that external novelty should remain
  unclaimed absent a fuller literature review.

Subject to these qualifications, the drafts are safe to use as rigorous
internal results.

## 2. Orthogonal product and scaling

Let

```math
T_{a,b}=\alpha R_a\oplus\beta S_b,
```

with the standard direct-sum Hilbert norm and `alpha,beta>0`.  For every pair
of parent inputs,

```math
\|T_{a,b}-T_{a',b'}\|^2
=\alpha^2\|R_a-R_{a'}\|^2
 +\beta^2\|S_b-S_{b'}\|^2.                         \tag{V.1}
```

If `gamma in Gamma(R)` and `eta in Gamma(S)`, multiplying their inequalities
by `alpha^2` and `beta^2` and adding proves

```math
(\alpha^2\gamma,\beta^2\eta)\in\Gamma(T).
```

Conversely, in a parent certificate set `b=b'`; this gives
`gamma/alpha^2 in Gamma(R)`.  Setting `a=a'` gives the other factor.  Hence

```math
\Gamma(T)=\alpha^2\Gamma(R)\times\beta^2\Gamma(S). \tag{V.2}
```

The inverse-Hamming quotient for a parent pair is a weighted average of the
two scaled child quotients, so it is bounded below by their minimum.  Varying
only a minimizing child attains that bound:

```math
\kappa(T)=\min\{\alpha^2\kappa(R),
                 \beta^2\kappa(S)\}.              \tag{V.3}
```

Thus the `t`-fold product with scale `1/sqrt(t)` has `kappa(R)/t`, as stated.
For complete formal scope either assume both latent cubes have positive
dimension or state a convention for the minimum over an empty set.  This is
only a degenerate-domain issue.  The zero-scale convention in the draft is
also correct: single-bit comparisons force every weight belonging to the
zero response to vanish.

## 3. Same-space cancellation and the joint repair

For `R_a=a e_1` and `S_b=b e_1`, each child one-bit separation polytope is
`[0,4]`.  But

```math
T_{1,-1}=T_{-1,1}=0.
```

The two inputs differ in both latent bits, so the defining inequality gives
`gamma_1+gamma_2<=0`.  Nonnegativity gives

```math
\Gamma(T^\parallel)=\{(0,0)\}.
```

For `S_b=b e_2`, squared response distance is exactly four times parent
Hamming distance, proving

```math
\Gamma(T^\perp)=[0,4]^2.
```

This is a conclusive counterexample to composition from the two child
polytopes alone.  The proposed repair is also correct.  Expanding

```math
\|\Delta R+\Delta S\|^2
=\|\Delta R\|^2+\|\Delta S\|^2
 +2\langle\Delta R,\Delta S\rangle
```

and subtracting the coordinate budgets `c,d` proves that
`(gamma-c,eta-d)` is a parent certificate whenever these differences are
coordinatewise nonnegative and the displayed cross-Gram lower bound holds.
The draft correctly does not infer this Hilbert-space condition from its
separate profile-synchronization result.

## 4. Max-plus collapse at every layer

For one binary `Q by Q` kernel, normalized Frobenius distance squared is
Hamming distance divided by `Q^2`.  Single-coordinate comparisons and their
sum give exactly

```math
\Gamma(F_1)=[0,Q^{-2}]^{Q^2},
\qquad \kappa(F_1)=Q^{-2}.                          \tag{V.4}
```

Fix a coordinate `(s,u,v)` of a `t`-factor product, with `t>=2`.  Put every
other factor equal to zero and compare `J` with `J-E_(u,v)` in factor `s`.
Both parent kernels are the all-one kernel:

- if `s=1`, then for fixed left endpoint choose the next state different
  from `v` whenever the left endpoint is `u`;
- if `s=t`, then for fixed right endpoint choose the preceding state
  different from `u` whenever the right endpoint is `v`;
- if `1<s<t`, both states around the exceptional transition are free, so
  choose any transition other than `(u,v)`.

The choice exists for every endpoint pair because `Q>=2`.  The two inputs
differ in exactly the selected bit but have identical output.  Its certificate
weight must therefore be zero.  Since this works for every factor and every
entry,

```math
\Gamma(F_t)=\{0\},\qquad \kappa(F_t)=0
\quad(t\ge2).                                      \tag{V.5}
```

As a finite check, the construction was exhaustively evaluated for all
layers and endpoints at `Q=2,3` and `t=2,3,4`; every pair of products was the
all-one kernel.  The analytic argument covers all `Q>=2,t>=2`.

This collision concerns factor-to-product identifiability, not the
sufficiency of the resulting boundary kernel for future endpoint contexts.
The draft keeps those two statements separate.

## 5. Robust four-cell crossing

Let `epsilon=||M-Mtilde||_infty`.  In a min-plus factorization of `Mtilde`,
every separable term majorizes `Mtilde`, and some term is tight at each finite
cell.  If one term were tight at two distinguished cells `i,j`, separability
would imply

```math
\widetilde M(x_i,y_i)+\widetilde M(x_j,y_j)
\ge
\widetilde M(x_i,y_j)+\widetilde M(x_j,y_i).        \tag{V.6}
```

The hypothesis and four perturbation errors instead give

```math
\widetilde M(x_i,y_j)+\widetilde M(x_j,y_i)
-\widetilde M(x_i,y_i)-\widetilde M(x_j,y_j)
\ge G-4\epsilon.                                   \tag{V.7}
```

This is positive when `epsilon<G/4`, a contradiction.  Distinct
distinguished cells therefore require distinct tight terms.  The constant
and strict inequality are both correct.

They are sharp in the class stated.  For

```math
M=\begin{pmatrix}0&1\\1&0\end{pmatrix}
```

the gap is `G=2`; the constant matrix with entries `1/2` is at uniform
distance `G/4=1/2` and has min-plus factor rank one.  Consequently no theorem
with the same hypotheses can include equality or replace `G/4` by a larger
universal threshold.

The source paper permits `+infinity` in factor functions, whereas the drafts
write finite real factors.  This creates no mismatch for finite matrices:
the lower-bound argument already permits infinities, and unused infinite
values can be replaced by sufficiently large finite values without changing
the represented finite table.

## 6. The average-MSE obstruction

For `D_r(i,i)=0` and `D_r(i,j)=1` when `i!=j`, its `r` diagonal cells satisfy
the crossing theorem with `G=2`, so its min-plus factor rank is at least `r`.
The explicit `r`-term representation

```math
u_k(i)=\mathbf 1_{i\ne k},
\qquad
v_k(j)=\mathbf 1_{j\ne k}
```

has minimum zero on the diagonal and one off it, proving equality.  The
rank-one all-one matrix differs on exactly the `r` diagonal entries, hence

```math
{1\over r^2}\|D_r-J_r\|_F^2={1\over r}.            \tag{V.8}
```

The conclusion is exactly the one claimed: an unweighted fooling set can
have vanishing mass under the uniform entry distribution, so fooling-set
cardinality and gap alone cannot imply an average-MSE rank lower bound.
Calling (V.8) “normalized mean-square error” is less ambiguous than
“normalized Frobenius distortion,” since the latter can also mean
`||D_r-J_r||_F/r=1/sqrt(r)`.

The repository's reproducibility script independently returns MSE `1/r` for
`r=2,4,8,16,32`.

## 7. Linear-code corollary and source audit

The primary source is Karthik Sheshadri,
[*Trellis State Complexity as an Exact Tropical Factorization Rank*](https://arxiv.org/abs/2607.23471),
arXiv:2607.23471v1 (26 July 2026).  Its definitions and normalization match
the draft:

```math
W(x_L,x_R)=d((x_L,x_R),C),
\qquad
s=\dim C-\dim C_L-\dim C_R.
```

Theorem 1 proves min-plus factor rank, tropical rank, and Kapranov rank all
equal `2^s`.  Lemma 2 selects one lifted codeword from each class of
`P_R(C)/C_R`; the corresponding block has diagonal zero and off-diagonal
entries at least one.  For two distinct selected diagonal cells, the two
cross entries are each at least one.  The robust four-cell hypothesis
therefore holds with `G=2`.  It follows that every `Mtilde` satisfying

```math
\|\widetilde W-W\|_\infty<1/2
```

has factor rank at least `2^s`.  Choosing `Mtilde=W` and invoking the exact
factor-rank upper bound proves

```math
\min_{\|\widetilde W-W\|_\infty\le\epsilon}
\operatorname{rank}_{\min,+}(\widetilde W)=2^s
\quad(0\le\epsilon<1/2).                           \tag{V.9}
```

Thus the corollary and its strict threshold are fully justified.  Notice
that (V.9) states the **minimum** approximate rank.  It does not assert that
every admissible approximant itself has rank exactly `2^s`; the proved
pointwise assertion is rank at least `2^s`.

Remark 7(ii) of the source asks how the rank degrades under additive
approximation and whether there is a robust state-complexity analogue.  It
does not contain (V.9).  Targeted searches of the foundational tropical-rank
literature and work on approximate tropical/subtropical factorization found
factor-rank definitions, exact fooling-set/crossing arguments, and numerical
Frobenius-loss methods, but not this exact uniform-gap corollary.  This is not
an exhaustive novelty search.  The safest claim is therefore:

> (V.9) is a correct elementary robust deduction from Sheshadri's explicit
> transversal, absent from that preprint, and a valid partial answer to its
> approximate-factorization question in the raw sub-half-unit regime.

It should not be advertised as a macroscopic robustness theorem.  For the
normalized distance landscape `W/m`, (V.9) permits error only

```math
\epsilon_{\rm normalized}<1/(2m),                 \tag{V.10}
```

which tends to zero at the lattice spacing.  Errors `epsilon m` in raw
distance, constant relative error, and average error remain completely open;
the `D_r` example explains why the last of these requires a mass/exposure
hypothesis.

## 8. Director-facing assessment

The response-separation results are valuable negative structure: they prove
that `Gamma` is an exact posterior-width certificate for a fixed response
embedding but not a closed max-plus feature algebra.  The robust tropical
theorem is a correct generative deduction and supplies a sharp growth law at
the integer query scale.  Its strongest defensible interpretation is
**lattice-robust min-plus channel incompressibility**, not response
compression at fixed normalized distortion.

The proposed next target—a positive-measure or query-exposed tropical
fooling-set theorem under average loss—is logically forced rather than a
cosmetic extension.  Any such theorem must add a quantitative mass condition;
Proposition TR.3 rules out a distribution-free version.
