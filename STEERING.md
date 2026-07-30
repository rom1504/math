# Strategic steering

Evidence cutoff: Wave 34, §10.87 (2026-07-30). Regenerate at the Wave 39
boundary, or earlier after a decisive change.

## User-stated research objective

Determine whether $\lim_{n\to\infty}M_n/n^{3/2}$ exists. The conjectural value $1/2$ is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains $0.336493364431\ldots\le\liminf\le\limsup\le1/2$.

The leading route is **adaptive optimized principal restriction**. Its sharp bare target is the arbitrary-cut tail (10.795), not an adversarial-selector fractional cover. Wave 34 shows that one eligible row-good coset need hit only an $\exp\{-O(TL_0)\}$ fraction of the **uniform** selector slice. The all-selector power sum (10.907) and common-weight capacity remain sufficient but stronger certificates.

The best structured implementation is now the **exceptional-center degree**
(10.967). It reduces the batch collision to one low-row center that is close
to favorable child fibers for a moderately rare set of uniform selectors.
This is a real quantifier improvement, but it does not itself prove a new
signing estimate.

The constant-shortfall pressure criterion (10.617) is **not** currently the
strongest route. It has no mechanism producing its required square-root
deletion reward, whereas the restriction and harmonic routes have exact,
weaker, and more locally testable remaining lemmas.

## Exact sufficient lemma and convergence chain

Fix $c_0\in(0,1/4)$ and a ratio window $m/n\in[\rho,1)$. Set

$$
L_0=n^{3/4-c_0},\quad k_0=\Theta(L_0/\log n),\quad r=\left\lceil n\log(2k_0)/L_0\right\rceil .
$$

Choose an integer $T=T(n,m)$ with $1\le T\le n^\eta$, fixed $\eta<c_0$, and
put $s=\lceil r/T\rceil$, $D=\Theta(k_0)$. Uniformly for every relevant
exact-minimizer target pair, it would suffice to prove

$$
\max_{z\in\mathcal C_2}U_m\{S:a_z(S)\le\lfloor D/s\rfloor\}\ge \exp\{-O(TL_0)\}. \tag{10.967}
$$

Here $\mathcal C_2$ is the admissible low-row center set and $a_z(S)$ is the
shortest projective distance from the restricted center to a favorable child
label. This is the concrete equal-radius sufficient lemma now sought; the
total-cost star event (10.926) is logically weaker.

Why it proves convergence:

1. One such center makes $s$ iid uniform selectors share a favorable star of total projective length at most $D$.
2. The star produces an eligible row-good block coset hitting those selectors.
3. Paying the candidate-coset count and at most $2^k$ support cuts yields one fixed cut with hit mass $\exp\{-O(TL_0)\}$.
4. With $c'=c_0-\eta>0$, that cut satisfies the arbitrary-cut lemma (10.795):
   $R_2(d)=O(n^{9/4-c'})$, favorable mass
   $\exp\{-O(n^{3/4-c'})\}$, and tolerance $O(n^{3/2-c'})$.
5. The restriction inequality, geometric-window summability, and exact landing prove convergence.

The logically bare target remains (10.795). Equivalently, it is enough to
prove directly that some eligible row-good coset has
$U_m(H_a)\ge\exp\{-O(TL_0)\}$. A proof need not pass through centers.

## Known obstructions and falsification criteria

- A center may have ambient mass as small as $2^{-n}$ without harming the
  project exponent; average center degree and Jensen first moments are
  therefore unnecessarily strong.
- The one-sided distance-deficit bound (10.969) supplies a necessary
  positive-part-Laplace obstruction. The exact parent mean (10.970) is the
  unknown endpoint excess and survives every selector-independent
  slack/replacement mixture, so first-moment reweighting is circular.
- A strict soft bound $q_*\ge e^{-ATL_0}$ with
  $\lambda=\Lambda rL_0/(D+1)$ and $A<\Lambda$ implies (10.967). The hard
  degree directly proves the star event but need not recover the same strict
  soft margin; hard and soft are not interchangeable.
- Global block replacement does not condition its witness into a favorable
  child fiber. The exact $A_9$ game has zero favorable mass in every optimal
  law, and every selector of size three, four, or five on $A_6$ has $V_B=V_0$
  at zero target slack. A viable project-scale route must add structure beyond
  unconditional minimax and, downstream, establish cross-selector overlap.
- Uniform degrees, generic incidence/KKT arguments, Johnson-entry trees,
  small-radius noise smoothing, and selector-independent mean tilts have
  already met recorded exponent or mechanism walls.

For one fixed sufficient parameter choice, its (10.967) certificate is
falsified by an unbounded actual-minimizer family satisfying

$$
\max_{z\in\mathcal C_2}U_m\{S:a_z(S)\le\lfloor D/s\rfloor\}=\exp\{-\omega(TL_0)\}.
$$

Falsifying the implementation as a whole requires this failure for every
admissible fixed parameter choice; even that would not falsify (10.795). The bare route is falsified only by ruling out every cut meeting the row,
tail, and tolerance bounds of (10.795), for every fixed saving choice, along
an unbounded actual-minimizer family. Abstract label systems and finite
zero-slack examples do not meet either falsification standard.

## Strongest independent alternative

The harmonic parent-Gibbs route has the exact identity (10.980): weighted
Dirichlet energy equals endpoint conditional Bernoulli KL minus an
edge-context migration covariance. Its present sufficient package is
(10.983): uniformly positive restoring scale and correlation,

$$
\widehat\kappa_t\ge\kappa_0,\qquad \rho_t\ge\rho_0,
$$

together with endpoint KL and adverse migration totaling
$O(n^{1/2-2c})$. This gives the exact parent entropy target through (10.900).
Adjacent-selector Hellinger control is still separate. Vertex crossings force
posterior selector inclusion, but the orientation edge has no such anchor.
Finite data show that resonance, covariance sign, and restoring constants
cannot simply be assumed benign.

## Constant-shortfall audit

Criterion (10.617) would suffice if every target pair admitted an exact root
and deletion path with

$$
\kappa_{\beta,i}(B_r)\ge
\frac{q_n}{n^{3/2}}\bigl[r^{3/2}-(r-1)^{3/2}\bigr]-K
$$

at every step, for fixed parameters. No such path theorem is known. A local
falsifier would have $\max_i\kappa_{\beta,i}=o(\sqrt r)$ on an unbounded
exact-minimizer family; the path version needs an unavoidable deletion
cutset. Neither has been established.

## Ranked routes

1. **Exceptional-center nonlinear tail.** Prove or falsify (10.967)/(10.969)
   directly for exact minimizers under uniform selectors.
2. **Direct arbitrary-cut tail.** Prove (10.795), or the one-coset uniform
   target (10.964), without imposing a center representation.
3. **Harmonic endpoint/migration control.** Prove (10.983) and the separate
   selector-Hellinger estimate.
4. **Soft conditional replacement.** Seek a project-scale penalized response
   gap plus cross-selector overlap that survives the $A_6/A_9$ walls.
5. **Exact terminal min-cut / quadratic susceptibility.** Retain as geometric
   and discrete-minimality fallbacks without Johnson-entry costs.
6. **Constant shortfall.** Keep as an exact fallback, but below the three
   sharper routes until a square-root deletion mechanism appears.

## Next-wave decision

Wave 35 should attack, independently: the nonlinear exceptional-center tail;
a direct uniform one-coset/arbitrary-cut argument; and harmonic endpoint KL plus migration. Do not restore an adversarial all-selector quantifier that the restriction theorem no longer needs.
