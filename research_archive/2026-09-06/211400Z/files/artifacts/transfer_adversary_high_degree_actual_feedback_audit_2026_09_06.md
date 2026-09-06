# Independent audit: the high-degree actual feedback probe

2026-09-06. Status: **PASS at the fixed first-marked-frame, fixed-operator-
norm scope stated below.** The finite Gaussian rigidity theorem was audited
separately in `transfer_adversary_gaussian_return_rigidity_audit_2026_09_06.md`.
This audit reconstructs the additional Boolean mixed-covariance and local
query argument of
`transfer_seed_high_degree_actual_feedback_probe_2026_09_06.md` rather than
inferring it from Gaussian marginal laws or from positive covariance.

The proof uses the already banked exact source-cut, local surgery,
coherent boundary-graph and finite-cube Stein modules. The relevant
sections were read directly in the marked local-noise separation, first
marked energy projection, open-mark Hall, boundary-graph/next-return,
two-factor flat transport, and Boolean stable-noise audit artifacts.
The conclusion is NOT an arbitrary-history theorem, a comparison across
matrix orders, or a proof of convergence of the original minimum.

## 1. Normalization, coefficients, and harmless approximation errors

Here `B=A/sqrt(n-1)`, `Q=B^2`, and every row of B has norm one.
The actual old fields are exactly

```
G=BS, D_i=S_i h_2(G_i), Y=BD, W=(S,G,Y,QS,QD).
```

Neither QS nor QD is assigned a Gaussian law. For the fixed residual
original-degree block P>3, the source and transported Schur mains are

```
R_P=sum_(a+3b=P) r_ab^2 Q^{circ(a+b)},
T_P=B R_P B, q_i=(T_P)_ii.
```

Same-original-degree local monomials are handled by the banked covariance
comparison, not by an unsupported orthogonality assertion. The odd local
degrees `a+b` give the exact subset lower bound for the q_i. In particular
the deterministic cutoff removes at most epsilon*n rows. The operator
and row-variance estimates in (1) follow from the unit diagonal of Q and
the Schur multiplier bound.

For ORTHONORMAL h_k, the correct probe source is

```
U_ki=q_i^(k/2) h_k(Z_Pi/sqrt(q_i))
    =He_k(Z_Pi;q_i)/sqrt(k!).
```

This is the normalization in the current source. Thus the Gaussian
covariance is `T_P^{circ k}`, not `k! T_P^{circ k}`. Conditional Gaussian
Hermite projection gives the mixed coefficient
`d_i=1_J alpha_ik/q_i^(k/2)` on the LEFT root. Consequently the proposed
mixed matrix is `D_k T_P^{circ k}`; the ensuing trace with Q is correctly
`tr(Q D_k T_P^{circ k} D_k)`.

At fixed P,k,epsilon, all d_i are bounded. Normal-ordering the k noise
copies and projecting their distinct original seed labels gives their
exact degree-kP main. The local surgery module controls the difference in
averaged L2. These errors are sufficient here: for arbitrary random
vectors X,Y,

```
||E X Y^T||_* <= sqrt(E||X||^2 E||Y||^2).
```

Thus an averaged L2 error o(1), against an averaged bounded second
moment, is o(n) in the required nuclear norm. This statement must be
used with fixed bounded diagonal coefficients; no unbounded cutoff
removal is made on a raw covariance error. Exact channel covariance
powers can be replaced by the Schur mains using bounded Gram Schur
multipliers and the averaged diagonal variance error.

## 2. Why the mixed covariance has the claimed stronger norm

After the local surgery, the right side contains k separate original-P
noise vertices, all of whose original seed slots are distinct. The left
response is expanded in noise Hermites with coefficients in the literal
W. Internal coherent collisions remain exact. Every coherent primitive
has at most three original slots, whereas every right noise has P>3.

For a fixed equality diagram, first remove inter-block restrictions by
finite inclusion-exclusion. In the resulting unrestricted diagram a
connected isolated whole-noise pair contributes its exact cross-root
covariance. Every non-isolated merge has a small proper-cut or retained-
label influence bound on a noise end. The only way that a noise could
be wholly consumed without such a gain is an isolated equal-degree
noise/noise pair: P>3 excludes the coherent primitive exception.

There are precisely the following sufficient dispositions.

1. Two vertex-disjoint nontrivial merges give two factors
   `n^(-1/2) polylog(n)`. Subsequent merges are Hilbert-contractive,
   including shared retained labels. Their matrix is entrywise
   `O(n^-1 polylog(n))`, hence has Frobenius norm polylogarithmic and
   nuclear norm o(n).
2. A single nontrivial star, with an unaffected isolated whole pair,
   gives an entrywise-small star factor times a covariance matrix whose
   Frobenius norm is `O(sqrt(n) polylog(n))`. The product again has
   polylogarithmic Frobenius norm. Its coherent coefficients are literal
   moments; no Gaussianization of them is involved.
3. If no unaffected whole pair remains, a lone nontrivial star must
   consume all k right vertices and have a left noise as its center.
   A right-centered star cannot contain all k>=3 right vertices. A
   coherent left primitive cannot be the center because it has at most
   three slots. Pull the left noise through its one B transport to its
   old G/Y source. A set of k' right noises has P*k' distinct labels;
   source primitives have at most three slots each. Hall matching
   therefore supplies k distinct primitive neighbors. Each matched
   pair gains `n^-1/2 polylog(n)`, since P>3. The source scalar is
   `O(n^(-k/2) polylog(n))`; summing the pulled-back B row costs at most
   sqrt(n). For k>=3 the final entry is O(n^-1 polylog(n)).

The graph alternative is elementary: a connected bipartite graph which
is neither an edge nor a star contains a four-vertex path, and its outer
edges are vertex-disjoint. Several nontrivial components also give two
disjoint merges. The slot-capacity condition, not graph terminology
alone, excludes the coherent-centered full star.

The exact-distinctness qualification is important. A whole covariance
factor is extracted ONLY in an unrestricted component diagram. If an
inclusion-exclusion equality joins that component to another component,
its whole-pair merge retains the shared label and gains a fixed-slot
influence factor. It can be paired with the previous star merge, or
with the retained merge from a second joined whole pair. If another
whole component remains unaffected, that component instead supplies the
Frobenius sqrt(n) bound. A repeated slot inside one squarefree primitive
is zero. Hence no raw distinct-label factorization is assumed. This
qualification was queried independently and has been made explicit in
the source's Section 2.1.

The sole surviving terms pair all k right noises wholly with left
original-P noises. What remains of the coherent coefficient is its
ACTUAL mean. Their sum is exactly `D_k T_P^{circ k}`. For `Cov(U_k)`,
both sides have k equal-degree branches, and a lone all-consuming star
is impossible; the same argument leaves `T_P^{circ k}`. This proves
both normalized nuclear comparisons in source (2), using the preceding
ordered raw replacements.

## 3. Full-return positivity and all constants in that step

The scalar g_i is odd and increasing: global reversal preserves H and
negates the literal V; auxiliary Gaussian symmetry handles the rest.
Its derivative projection gives

```
<gbar,h_1> = average_i 2 q_i E[H_i phi(V_i/sigma_i)/sigma_i].
```

At least half the rows have `q_i>tau_P^2/2`. The old marginal yields
`E H_i=mu+o(1)` uniformly, and the actual averaged V second moment is
at most L^2+o(1). The proposed cutoff `R0^2=8 L^2/mu` therefore leaves
enough H mass for the conservative positive constant c0 in (4).
The finite Gaussian rigidity audit gives one odd k in the fixed finite
set up to K0, with the claimed coefficient lower bound.

Removing q_i<=epsilon*tau_P^2 costs at most
`tau_P epsilon^(3/2)` in any averaged Hermite coefficient. Applying the
signed rank identity to T_P and using
`Q >= T_P/(L^2 tau_P^2)` proves the retained correlation lower bound.
The return is then restored to the FULL C by

```
average |E[(B(C-C^J))_i E_ki]|
 <= L sqrt(average E|C-C^J|^2) sqrt(average E E_ki^2)
 <= L^2 sqrt(epsilon)+o(1).
```

This does not use a coefficient supremum growing with the variance
cutoff. The source ideal variances of `D_k U_k` are at most one. Both
inequalities in (5) therefore give the stated `a_*>0` for the actual,
unmodified BC. The finite choice of k may depend on n without changing
the order of limits: all comparison errors are uniform over that fixed
finite set.

## 4. The literal regressed query comparison passes

This is the step which cannot be inferred from Section 3. The probe's
exact main has degree p=kP and proper-cut norm
`O(n^-1/2 polylog(n))`, by two-factor flat transport. At a fixed
polynomial stage, inspect every homogeneous degree q in each literal
query. For q<p, all derivative contractions are proper for the probe.
For q=p, retain the exact covariance (and remove it by the stated
regression for BC). For q>p, the full-p contraction needs a separate
bound:

* An old Z channel has small proper cuts on its right tensor; this
  controls its full-p contraction. At q=p its covariance with the probe
  is small by the source Hall argument below, so no old-channel
  regression is needed.
* The q>p part of `Bc0(W)` has small proper cuts by the exact coherent
  boundary-graph theorem, since q>p>3. Its internal Boolean equality
  diagrams are retained by that theorem.
* Each main term of `B[(A-a)Z+R_{>=2}]` contains at least two positive
  rooted factors, and two-factor flat transport gives the required
  small right cut.
* For `BD_aB r(G,Y)`, pull BOTH root transports outside the contraction.
  The left source is exactly k original-P noises and has NO additional
  coherent factor. With all p left labels contracted and the q-p right
  labels left OPEN, Hall matching again gives k distinct old primitive
  neighbors. Every selected merge is proper or retains an influence
  label. The source tensor satisfies

  ```
  sup_(a,b) ||H_(a,b,free)||_F = O(n^(-k/2) polylog(n)),
  ||H||_F = O(n^(1-k/2) polylog(n)).
  ```

  Applying B and `BD_aB` on the separate root axes is bounded; selecting
  the diagonal output roots is a coordinate projection, of norm one.
  The resulting squared norm divided by n tends to zero. There is no
  root-l1 loss here. Exact exclusions touching a free mark retain that
  mark and use the influence bound; it is not summed as a scalar.

This is a direct k-branch strengthening of the old two-branch open-mark
Hall proof. Its former centered-one-noise obstruction is absent because
the present probe source has k>=3 noise factors and no coherent factor.
It also proves the small equal-degree old-channel covariance. Its
coefficient/root transports remain bounded at each fixed approximation
and variance-cutoff stage.

These are precisely the derivative-contraction conditions in the
finite-cube characteristic-function proof. Apply that proof to an
ARBITRARY JOINT characteristic function of the finite query list, not
separately to each coordinate. The resulting comparison retains the
actual joint law of

```
(W_i, the retained old Z channels, BC_i-beta_i E_ki)
```

and makes the probe an independent Gaussian with its actual variance.
No Gaussian law is assigned to the regressed query. Positive covariance
or separate independence tests would not have established this claim.

For bounded/hard responses first keep P,k,epsilon, the response
approximation, an upper root second-moment cap, and a positive probe
variance cutoff fixed. These make the regression coefficients bounded.
The old finite-catalog approximation applies because this SPECIFIC W
has the audited exponential moments. Averaged L2 errors and the root
caps transfer the fixed-degree contraction conditions. Then take n to
infinity and remove approximation errors in the ordered way. Old-channel
degree tails restore full Z, hence BF=V+Z, jointly with the regressed
return. The old small-ball estimate handles sign(BF) and the endpoint
spins. Root cutoffs can either be removed in the existing ordered-gain
fashion or simply kept on a fixed positive fraction of useful roots.

## 5. Cap consequence and its limits

The independent endpoint proof
`transfer_fresh_conditional_actual_feedback_endpoint_gain_2026_09_06.md`
uses exactly the joint comparison established in Section 4 above. Its
conditional premise is therefore available at this scope. Its fixed
cutoffs give a positive fraction of roots with nonzero probe regression,
and the independent Gaussian makes both required endpoint orientations
have a positive density of negative oriented local fields. Exact
independently thinned spin flips, or the equivalent multilinear cube
move, yield a strict expected oriented energy gain. In particular the
source's qualitative conclusion

```
liminf_n [Lambda(B_n)-j_n] > 0
```

passes for each FIXED admissible f and FIXED L. The endpoint normalizing
factor is `Lambda=Q_abs(A)/(n sqrt(n-1))`. Both actual endpoint fields
are retained: `Bu_+=BF+BC`, `Bu_-=-BF+BC`. Negative-energy ascent uses
the reversed orientation at u_-; there is no artificial bridge
cancellation or global-reversal shortcut.

This does NOT show that the gain dominates a spectral-core loss of
order 1/L, give a new unrestricted decimal lower bound, control an
arbitrary longer history, transfer actual minimizers across orders, or
decide the original convergence question.

## 6. Replayed finite checks and their exact evidentiary role

Read and reran
`computations/transfer_seed_high_degree_feedback_checks_2026_09_06.py`
using `.venv/bin/python`. Output:

```
signed_rank_cases 196
gaussian_mixed_normalization_cases 3240
bipartite_graph_cases 2452
actual_spin_ascent_cases 72 strict_finite_cases 70
PASS: finite algebra only; no asymptotic mixed closure certified by computation
```

The mixed Gaussian checks use an explicit non-Gaussian coherent shift;
the spin checks enumerate actual Boolean cubes. Neither proves the
asymptotic covariance or local comparison. Those conclusions rest on
the contraction arguments above, with the exact scope and ordered
limits retained.
