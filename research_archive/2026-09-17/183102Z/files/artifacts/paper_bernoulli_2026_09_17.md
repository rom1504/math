# Bernoulli decomposition and margin-preserving sign completion

2026-09-17. Initial primary mapping frozen before consulting historical
gamma2/augmentation artifacts. Status: derived comparison theorems,
counterexamples, and stronger landscape/extension criteria. Sections
11–12 received independent mathematical PASS; individual audits are
recorded below. No improvement to the original cap interval or proof
of convergence is claimed.

**Applicability correction (18:30 UTC).** A new universal Hamming-sphere
argument proves that every bounded-cap sequence has full near-level
entropy at least order `eta log(1/eta)`. Therefore the finite full-code
entropy-slope criteria audited below are VACUOUS in the original bounded
cap class. The positive-e0, finite-K_es regime in Section 12 is likewise
impossible there. A further universal response bound rules out whole-code
`mu(eta)^2=o(eta)` for every isotropic sign law. These corrections do not
invalidate the finite comparison/realization lemmas, nor the ANCHORED
Hamming-cover criteria, whose residual variance pays the forced local
entropy. See Section 24 for the exact quantifiers and proof.

## 1. Primary mechanism and exact scope

Primary sources: Bednorz--Latala,
[On the boundedness of Bernoulli processes](https://annals.math.princeton.edu/2014/180-3/p08),
Annals of Mathematics 180 (2014), 1167--1203;
[primary preprint](https://arxiv.org/html/1305.4292).
The journal PDF is preserved under
`tmp/paper_portfolio_2026_09_17/bernoulli/bednorz_latala_2014.pdf`.
The decisive proof components being reconstructed are Proposition 2.10,
Theorem 3.1, the chopped-process Corollary 5.3, Proposition 6.2, and
the final assembly in Section 7.

For deterministic finite T in R^d, define
`b(T)=E_epsilon sup_(t in T)<epsilon,t>` and
`g(T)=E_g sup_(t in T)<g,t>`, with independent standard signs and
standard Gaussian coordinates, respectively. Their theorem supplies a
pointwise decomposition `t=a(t)+z(t)` with

```
sup_t ||a(t)||1 <= C_BL b(T),
g({z(t):t in T}) <= C_BL b(T).                       (1)
```

The constructed z(t) freezes each coordinate at the first large jump
of a suitable admissible chain. It need not belong to span(T), preserve
rank-one features, or preserve coefficient signs. The construction is
not a direct decomposition of a quadratic chaos. Its use below is only
inside an upper bound for a LINEAR random process; it does not replace
the actual sign matrix or its Boolean maximum.

## 2. Derived comparison for bounded subGaussian rounding

Let X in R^d be centered, with `||X||infinity<=B` almost surely, and
assume every linear form obeys

```
E exp(lambda <X,u>) <= exp(lambda^2 L^2 ||u||2^2/2)
                         for all u and lambda in R. (2)
```

Then for every deterministic finite T,

```
E sup_(t in T)<X,t> <= C (B+L) b(T).                 (3)
```

Proof: use (1). The a-part is pointwise at most
`B sup||a||1`. The z-part is a centered subGaussian process with its
increments bounded by the L-scaled Euclidean metric. The generic
chaining upper bound and Gaussian majorizing-measure equivalence give
`E sup_t<X,z(t)> <= C L g(z(T))`. Combine the bounds. Translation of
T is harmless because X and the reference signs are centered.

Every assumption may instead hold conditionally on a sigma-field G,
provided the finite index set T and its margins are G-measurable. Then
(3) holds conditionally; no independence of the coordinates of X is
claimed or needed. This is the precise point where subGaussian sign
rounding can replace a Gaussian-width certificate by a Bernoulli-width
certificate. Historical collision and actual-value benefit remain to
be tested.

## 3. Exact scalar-margin completion implication

Let E be a specified collection of undecided edges of an actual sign
quadratic, with all other edges fixed. Let r0 in [-1,1]^E be a fractional
center and let

```
H_r(x)=h0(x)+sum_(e={i,j} in E) r_e x_i x_j.
```

Fix a proposed full scalar cap U and suppose the signed margins

```
m_(x,sigma)=U-sigma H_(r0)(x)>0
                    for every Boolean x and sigma in {+1,-1}.
```

Define the FINITE deterministic margin-normalized index set

```
T_U={0} union
 { (sigma x_i x_j/m_(x,sigma))_(e in E) : x,sigma }.
```

If a rounding law R in {+1,-1}^E has centered displacement X=R-r0
satisfying (2) and the coordinate bound B, and

```
C(B+L) b(T_U)<1,                                    (4)
```

then some outcome is a full sign completion with `||H_R||infinity<U`.
Indeed Z=sup_(t in T_U)<X,t> is nonnegative; (3)--(4) give an outcome
with Z<1, and its individual inequalities are exactly
`sigma H_R(x)<U`. Thus no supremum is interchanged with an infimum,
and no averaging over physical Boolean inputs replaces their maximum.
Zero margins require exact protection or separate treatment; they cannot
be divided by zero or silently discarded.

For one-vertex augmentation this becomes
`H_new(x,s)=H_A(x)+s r dot x` and
`||H_new||infinity=max_x(|H_A(x)|+|r dot x|)`.
At fractional center zero, write `D_x=Q(A)-|H_A(x)|`.
The exact target `Q(A)+delta` is governed by
`T_delta={0} union {+/-x/(D_x+delta)}`.
Obtaining a useful original recurrence still requires a bound on this
actual near-level Bernoulli geometry for selectable exact minimizers.

Important baseline: with no extra constraints, independent sign rounding
already satisfies `E sup_T<epsilon,t>=b(T)` with constant one. Therefore
(3) by itself cannot improve that independent-rounding argument. Its
potential benefit is retaining additional constraints or anisotropic
control furnished by a dependent rounding law, while avoiding an
unnecessary Gaussian-width overpayment.

## 4. Chain-rule source under examination

Chu--Raginsky,
[A Chain Rule for the Expected Suprema of Bernoulli Processes](https://arxiv.org/html/2304.14474),
2023, applies the BL decomposition to a Lipschitz composite index class.
Its Gaussian part is a different set S, not necessarily the original T;
the cost of the uniformly Lipschitz function class must also be paid.
For a general vector-to-vector composition both l1 and l2 Lipschitz
control are relevant. This is not yet an application to the discontinuous
sign map or an automatic license to flatten quadratic-response sets.

## 5. Primary proof reconstruction: what the decomposition actually does

The decisive BL mechanism is not merely a Gaussian comparison. Proposition
2.10 separates coordinates J from their complement. Independent copies on
J let the argument compare a local Bernoulli supremum with the full one;
Sudakov separation on the complementary coordinates either supplies a
small covering or pays a definite decrease in the local functional. The
independence is that of the REFERENCE Bernoulli coordinates. One cannot
repeat the coordinate-removal step on an arbitrary dependent rounding law.

Theorem 3.1 converts an admissible multiscale chain into the decomposition.
For each coordinate, freeze the chain representative at its first jump
larger than the current scale. The sum of these first-large-jump errors is
paid in l1; the frozen process is paid by a Gaussian chaining functional.
The relevant budget is `sup_t sum_n 2^n r^(−j_n(t))`. Under the theorem's
capped-square increment hypothesis, the two costs are respectively at
most `37 M` and `C sqrt(M)` times this budget. Chopping maps isolate the
large-coordinate increments without destroying the contraction property
of the independent Bernoulli process. Sections 5–7 construct partitions
whose functional drops telescope and bound the budget by `C b(T)`.

This explains both the power and the scope limitation of (1): coordinate
freezing is nonlinear in t. The Gaussian piece is not required to remain
inside any prescribed linear image of the original feature set.

## 6. Fixed-law decomposition: the 2026 source and its decisive proof

Liu--Zadik, [A Bayesian Proof of the Bernoulli Theorem,
arXiv:2608.11031v2](https://arxiv.org/html/2608.11031v2), supplies a stronger
distributional formulation. For a finitely supported law mu on R^d put

```
B(mu) = sup E<epsilon,U>,  U~mu, epsilon~uniform independent-coordinate signs;
G(mu) = sup E<G0,U>,       U~mu, G0~N(0,I).
```

Each supremum is over ALL couplings of the two displayed marginals.
The coordinates of epsilon are independent of each other, but epsilon
and U are generally dependent. Define

```
delta(mu) = inf_(a:T->R^d) [ E_mu||a(U)||1 + G((id−a)#mu) ].
```

Their Theorem 2.1 gives `delta(mu) <= C_LZ B(mu)` for every prescribed mu.
The explicit constants in Lemmas 5.1 and 6.1 give the coarse admissible
value `C_LZ=2400 pi`. No claim of a useful sharp constant is made.

Here is the proof mechanism of the imported implication, independently
reconstructed from Sections 5.3 and 6. For scale t let

```
phi_t(x,u) = sum_i min(1, |x_i−u_i|^2/t^2),
RD_mu(t) = inf_(X,U both~mu) [ I(X;U) + E phi_t(X,U) ].
```

Take dyadic scales r_k=2^(−k) and optimal self-couplings (X,U_k), made
conditionally independent given X. Start with U_k0 independent of X and
end with U_k1=X. Set

```
Delta_k = clip_(3r_k/2)(U_(k+1)−U_k),
Vtilde = U_k0 + sum_k Delta_k,
V(x) = E[Vtilde | X=x],  A(x)=x−V(x).
```

For each coordinate, telescope across maximal runs for which
`|X_i−U_(k,i)|>r_k`. Endpoint terms plus clipped increments cost at most
`(17/2)r_p` for a run starting at p. Consequently

```
E||A(X)||1 <= (17/2) sum_(interior k) r_k E phi_(r_k)(X,U_k).
```

For any coupling of X with a standard Gaussian G0, extend the U_k kernels
conditionally independently of G0. The initial U_k0 term has zero mean.
Writing d_k=E phi_(r_k)(X,U_k), q_k=I(X;U_k), one has

```
E||Delta_k||2^2 <= (9/4)r_k^2(d_k+d_(k+1)),
I(G0; U_k,U_(k+1)) <= q_k+q_(k+1).
```

The Gaussian mean-information inequality
`E||E[G0|W]||2^2 <= 2 I(G0;W)`, followed by Cauchy--Schwarz, pays the
Gaussian coupling objective by
`(9/(2sqrt(2))) sum_k r_k(d_k+q_k)`.
The endpoints vanish; monotonicity of RD gives
`sum_k r_k RD_mu(r_k) <= 2 integral_0^infinity RD_mu(t)dt`.
Thus `delta(mu)<=30 integral RD_mu`.

For the reverse comparison, take the exact type class of N-tuples with
empirical law mu. Its entropy deficit is at most
`(|supp(mu)|−1)log(N+1)`. Conditional entropy subadditivity gives
`RD_(type law)(t)>=N RD_mu(t)−deficit`. Integrate only to
`sqrt(N/(deficit+1))`; the independent self-coupling bounds the remaining
tail by `diam2(T)^2/t^2`. Apply the paper's set-level Cauchy-channel bound
`integral RD <=80pi b(T)`, divide by N, and use the exact-type identity
`N^−1 b(type class) -> B(mu)`. Rational laws are dense, and maximal-coupling
bounds make both B and the RD integral continuous in mu. This proves
`integral RD_mu<=80pi B(mu)` for every finite law.

The fixed-law result is genuinely useful when a maximizing index has a
prescribed marginal, since one must keep the offsets associated with
that same index distribution. It does not identify that distribution
for the original extremal quadratic problem.

## 7. All affine margins: bounded subGaussian convex comparison

**Theorem.** There is an absolute C such that every centered X satisfying
(2) and `||X||infinity<=B` obeys

```
X <=_(convex order) C(B+L) epsilon.                  (5)
```

Equivalently, for every finite family of slopes t and arbitrary real
offsets m_t,

```
E max_t [m_t+<X,t>]
    <= E max_t [m_t+C(B+L)<epsilon,t>].             (6)
```

In particular, the offsets are NOT multiplied by C. The theorem does
not follow simply by applying a homogeneous comparison to each level
set and summing: that would change the margin accounting.

There are two proofs, which clarify attribution.

**Proof A (BL plus tensorization).** Independent concatenations of X
retain the same B and L. Apply (3), with a dimension-independent constant,
to the exact-type slope set in R^(dN), divide by N, and let N tend to
infinity. Liu's tensorization principle, as presented in Proposition 3.1
of van Handel, [On the subgaussian comparison theorem,
arXiv:2512.18588v2](https://arxiv.org/html/2512.18588v2), gives

```
F_X(mu):=sup_(U~mu coupled to X) E<X,U>
                 <= C(B+L) B(mu).                  (7)
```

First do this for rational mu and then use continuity in the finite-law
probabilities. Finally use the exact identity

```
E max_t[m_t+<X,t>] = sup_mu [ F_X(mu)+E_mu m_U ].
```

This identity follows by choosing a maximizing index, not by exchanging
a minimization and a maximization. Equation (6) follows, and monotone
approximation of convex functions by finite maxima of affine functions
gives (5). Subtract an affine lower bound before using monotone convergence.

**Proof B (fixed-law decomposition plus Gaussian convex comparison).**
Theorem 1.1 and Corollary 1.2 of van Handel give a coupling with
`X=c_sg L E[G0|X]`. For any coupling (X,U), extend the Gaussian coupling
conditionally on X, independently of U. For every a:T->R^d,

```
E<X,U> <= B E||a(U)||1
             + c_sg L G((id−a)#mu).
```

Optimize a and use Section 6 to obtain
`F_X(mu)<=2400pi max(B,c_sg L) B(mu)`. Keep the same mu in the offset
identity to conclude (6). This proof provides an explicit coarse factor
apart from the Gaussian comparison constant, but the new Cauchy proof
is NOT essential: Proof A already combines BL with tensorization.

Conditional versions hold by applying the theorem to each conditional
law, with the slopes and offsets measurable before the current rounding.
Conditioning later on the selected rounding outcome does not preserve
these hypotheses automatically.

For the actual sign-completion setting of Section 3, take
`m_(x,sigma)=sigma H_(r0)(x)` and slopes
`t_(x,sigma)=(sigma x_i x_j)_(e in E)`. Thus (6) bounds the expected FULL
scalar cap of the actual rounded matrix, retaining every existing energy
offset. However, independent rounding already corresponds to factor one.
Without extra constraints or a useful anisotropic refinement, the
comparison cannot improve the independent-rounding construction.

## 8. A valid anisotropic certificate, and a false shortcut

Suppose instead that
`E exp(<u,X>)<=exp(u*Q u/2)` for a positive semidefinite Q, while
`||X||infinity<=B`. Gaussian convex comparison, after restricting to the
range of Q, gives the rigorous all-margin upper bound

```
E max_t[m_t+<X,t>]
 <= inf_(a_t) E_G max_t[
       m_t+B||a_t||1+c_sg <G,Q^(1/2)(t−a_t)> ].     (8)
```

Indeed the a-part is bounded pointwise, and van Handel's affine theorem
handles the remaining slopes. This is an upper certificate, not a proved
sharp characterization. It preserves the actual Q and every m_t.

A tempting replacement is FALSE: for arbitrary A, the assumptions

```
X centered, X in A[-1,1]^k a.s.,
E exp(<u,X>) <= exp(||A*u||2^2/2)
```

do not imply `X <=_(convex order) C A epsilon` for any universal C.
This is already false in physical dimension one.

For an integer m>=2 let A_m be the row with one entry 1 and m^3 entries
equal to m^(−2). Its zonotope is `[−m−1,m+1]` and
`A_m A_m*=1+1/m`. Let `X_m=clip_m(G)`, with G standard Gaussian.
Symmetry and `|X_m|<=|G|` give
`E exp(lambda X_m)=E cosh(lambda X_m)<=exp(lambda^2/2)`.
All the proposed hypotheses hold. But

```
A_m epsilon = epsilon_0+Z_m,
E Z_m=0, E Z_m^2=1/m.
```

For a fixed C>0 test the convex function `f_C(x)=(|x|−2C)_+`.
The pointwise inequality `(s−1)_+<=s^2/4`, s>=0, gives

```
E f_C(C A_m epsilon) <= C E(|Z_m|−1)_+ <= C/(4m).
```

On the other hand, as m tends to infinity,
`E f_C(X_m)->E(|G|−2C)_+>0`. Hence no fixed C works.
Writing `h(a)=2(phi(a)−a Phi(−a))`, the exact left side for m>2C is
`h(2C)−h(m)`, giving a finite, non-asymptotic witness whenever this
exceeds `C/(4m)`.

Even requiring X_m to be the image of an ACTUAL sign vector does not fix
this counterexample. For even m, the image `A_m{−1,+1}^(1+m^3)` contains
every point in the lattice `(2/m^2)Z` between −m−1 and m+1. Replace X_m
by its symmetric rounding toward zero to this lattice. The result is
still 1-subGaussian and centered and has an actual (dependent) sign lift.
The left-side hinge expectation decreases by at most `2/m^2`, so the
same contradiction follows. This lift is NOT claimed to obey an ambient
dimension-free subGaussian bound on the full sign vector; that stronger
condition is relevant to identity-augmented GS constructions.

The discarded information is the intermediate-scale tail profile of the
weighted independent sign sum. Its support is of order m and its variance
is near one, but most of that variance lies in a single bounded coordinate;
the remaining large support has vanishing variance. The two crude data
do not constrain X to have that same intermediate-scale behavior.

BL Problem 1.3 separately asks for a Banach-space-valued decomposition.
Its stated span-preservation obstruction explains why one cannot pull
an arbitrary BL decomposition of A*T back through A*. It is NOT needed
to prove the elementary counterexample above, and no claim about the
current open/closed status of that broader problem is made here.

## 9. Historical collision and current value status

After freezing Sections 1–3, the following archive was inspected:

- `artifacts/correlated_internal_fill_rounding.md` already contains the
  actual affine cap margins, exact-protection kernel, a hypothesized
  constrained subGaussian sign law, and a generic-Gaussian-chaining
  sufficient criterion. Sections 2 and 7 replace its Gaussian-width
  payment by a Bernoulli comparison when the relevant law exists.
  They do not supply that law or show its cap certificate is small.
- `artifacts/principle_invent_2026_09_07_hierarchy_complexity_dilution.md`
  already uses GS plus a multilevel near-extreme response hierarchy.
  Its gamma2 is a MATRIX FACTORIZATION norm, not the generic-chaining
  gamma2 in BL. The two notions must not be identified. The all-affine
  theorem avoids peeling as a formal comparison, but does not estimate
  the resulting penalized Bernoulli process for actual exact minimizers.
- `artifacts/retrieval_panel_2026_08/discrepancy_toolkit.md` already lists
  BL and its random-linear-process versus deterministic-minimax boundary.
  This campaign's additional outputs are the fixed-offset comparison,
  its exact source dependence, and the explicit zonotope falsifier—not
  a claim that BL itself is a new discovery.

No improvement to the original lower constant, upper constant, or
convergence problem follows yet. To affect those values one still needs
a quantitatively useful sign law with retained anisotropic constraints,
and an estimate of its actual penalized response process.

## 10. Independent audit of the cumulative-precision GS theorem

Independently read
`artifacts/paper_discrepancy_2026_09_17.md`, its complete replay script,
the relevant archived hierarchy Section 7, and the primary JASA
supplement PDF pages 46–51 and 106. The primary proof establishes the
variance parameter one for arbitrary initial bias; the random phase
compensator must not be replaced by its mean in an exponential bound.

The cumulative-precision augmentation theorem in that artifact passes
this audit, including the coefficient `9/8`. Details checked:

1. For a query in nested level j, its ancestral factorization vectors
   represent the same original sign word. Minimizing the squared norm
   of their weighted direct-sum representation gives
   `n/(sum_(k<=j) a_k^2/g_k)`. No independence of the feature blocks is
   assumed; GS supplies one joint Euclidean subGaussian bound.
2. A sign code of factorization square at most ng has VC dimension at
   most ng. On a shattered d-coordinate set, use the factorization row
   for each labeling epsilon and average
   `d^2 <= gamma_fac^2 ||sum_i epsilon_i v_i||^2`.
   This gives d<=gamma_fac^2 and the stated Sauer entropy bound.
3. With required precisions b_j, the running maximum B_j is enough.
   Since g decreases down the hierarchy,
   `sum g_j(B_j−B_(j−1))` is bounded by the weighted positive variation
   integral, including its initial boundary term. The possible terminal
   jump costs at most `g_J b_J=o(1)`.
4. If `g h(g)=kappa eta`, differentiating gives
   `g (h(g))'=kappa rho(g)`, where rho tends to 1/2. Scaling
   `eta=R tau epsilon u` leaves the integral
   `integral_0^infinity [2u/(1+u)^3−1/(2(1+u)^2)]_+ du=9/16`.
   The positivity interval is u>1/3. The subGaussian tail factor two
   therefore gives `9(1+zeta)kappa R/(8tau)`.
5. The new-word entropy contribution is O(sqrt(epsilon)); the initial
   boundary term is O(epsilon). A fixed identity-feature weight pays
   responses outside the first near-level code. Every shell has enough
   old-energy deficit to pay its corresponding bridge allowance.
6. A full-sign new principal block with cap at most q^(3/2) exists for
   all large q by the elementary independent-edge union bound. The old
   block is preserved exactly. Both polarities and every new word are
   included, so the scalar cap accounting is valid.
7. The hierarchy is finite for each fixed epsilon; n tends to infinity
   only after all hierarchy parameters are fixed. Hence the limsup
   factorization majorants apply simultaneously. Letting epsilon tend
   to zero last is legitimate.

Consequently the sufficient coefficient `tau>9K/8` and necessary
liminf-landscape inequality `K>=4c_*/3` are verified. This is a genuine
quantitative strengthening of the archived landscape obstruction, not
a proof of a new cap constant or convergence. The reproduction script
passed all its deterministic checks when independently rerun.

The later actual-entropy Section 9 was also independently checked. Its
criterion `2K_es/tau+2(log2)e0/tau^2<1` and mixed necessary condition
pass. In particular, its new-word precision increments correctly charge
the terminal plateau only once. The extension below improves that
part of the argument by treating the new-spin maximization exactly.

## 11. Exact new-spin composition: an improved mixed augmentation criterion

This section's finite lemmas and formal extension received an independent
mathematical PASS. Its positive-e0/finite-K regime is now ruled out for
bounded-cap families by Section 24; only the zero-e0 weighted-entropy
regime remains an unexcluded conditional route. Its
mechanism is GS plus the exact absolute-value composition and a
Gaussian-integration bound for squares. BL and the 2023 chain-rule
paper motivated examining the composition, but neither is needed as
an imported premise of the following elementary improvement.

### 11.1 Finite independent-column lemma

Let Z_1,...,Z_q be independent centered real variables with
`E exp(lambda Z_j)<=exp(lambda^2 sigma^2/2)`. Then

```
max_(y in {−1,1}^q) |sum_j y_j Z_j| = sum_j |Z_j|.
```

A first bound is

```
P(sum_j |Z_j| > q sigma+2 sigma sqrt(q u)) <= exp(−u).
```

Indeed `E|Z_j|<=sigma`; Jensen and an independent copy Z_j' give

```
E exp(lambda(|Z_j|−E|Z_j|))
 <= E exp(lambda(|Z_j|−|Z_j'|))
  = E cosh(lambda(|Z_j|−|Z_j'|))
 <= E cosh(lambda(Z_j−Z_j')) <= exp(lambda^2 sigma^2).
```

No symmetry of Z_j is required: exchangeability of the pair supplies
the cosh identity, and the reverse triangle inequality supplies the
pointwise contraction. Independence across j then gives the displayed
tail bound.

There is a stronger bound that keeps the leading old-code entropy
coefficient unchanged. For `0<=t<1/2`, Gaussian integration gives

```
E exp(t sum_j Z_j^2/sigma^2)
 = E_G E_Z exp(sqrt(2t) sum_j G_j Z_j/sigma)
 <= E_G exp(t ||G||2^2) = (1−2t)^(−q/2).
```

Consequently, for w>=1,

```
P(sum_j |Z_j| >= q sigma w)
 <= exp(−q J(w)),
J(w)=(w^2−1−2log w)/2.                             (9)
```

Here Cauchy--Schwarz reduces the event to
`sum Z_j^2>=q sigma^2 w^2`, and the optimizing Chernoff parameter is
`t=(w^2−1)/(2w^2)`. A convenient explicit consequence is

```
P(sum_j |Z_j| > sigma sqrt(q[q+2sqrt(q u)+2u]))
                                                  <= exp(−u). (10)
```

For example, (10) follows by minimizing the square-MGF bound, or from
`−log(1−2t)−2t <= 2t^2/(1−2t)` and optimizing t.
There is no separate `2^q` union cost in (9) or (10).

### 11.2 Cumulative-precision theorem

Use the exact definitions e(eta), s(eta), e0, and K_es from the audited
GS artifact. Write K=K_es here. Suppose K<infinity. Then the actual
full-sign extension conclusion

```
limsup_n [Q(B_(n+floor(epsilon n)))−Q(A_n)]/n^(3/2)
                  <= tau epsilon+epsilon^(3/2)      (11)
```

holds for every positive tau satisfying

```
2K/tau + e0/tau^2
       + (pi/2) sqrt(K e0)/tau^(3/2) < 1.           (12)
```

As usual, epsilon is sufficiently small and FIXED before n tends to
infinity; the old block is preserved and every new edge is a sign.
One may take the better of (12), the audited all-new-word union
criterion, and the factorization-only `tau>9K_gamma/8` criterion.

Proof: keep the same nested finite hierarchy and GS feature construction.
For x in level j, the q independent columns have scalar responses with
proxy `sigma^2=n/B_j`. Applying (10) and unioning only over x in E_j
requires the precision

```
b_j=(1+zeta)
 [epsilon^2+2epsilon^(3/2)sqrt(s_j)+2epsilon s_j]
 /(eta_j/R+tau epsilon)^2,                          (13)
```

with terminal denominator `(tau epsilon)^2`. Fixed finite union costs
and all `o(n)` entropy errors are absorbed by the strict factor.
Allocate features by the same running maximum B_j=max_(k<=j)b_k.

Split (13) into old-entropy, mean, and cross terms. With
`g_j s_j<=kappa eta_j`, the old term has resource at most
`2(1+zeta)kappa R/tau+o(1)`, exactly as in the audited proof. The mean
term increases down the hierarchy. Its increments, including the final
one, telescope after splitting at a fixed delta>0, giving resource at
most `(1+zeta)e0/tau^2+o(1)`.

For the cross term, monotonicity of s gives, at nonterminal levels,

```
(sqrt(s_j)f_j−sqrt(s_(j−1))f_(j−1))_+
       <= sqrt(s_j)(f_j−f_(j−1)),
f_j=(eta_j/R+tau epsilon)^(−2).
```

Use `g_j sqrt(s_j)<=sqrt(g_j kappa eta_j)`. Above a fixed delta, the
cross resource tends to zero. Below delta, `g_j<=e(delta)+o(1)`.
The remaining integral is bounded by

```
(4(1+zeta)epsilon^(3/2)/R) sqrt(kappa[e(delta)+o(1)])
       integral_0^infinity sqrt(eta)/(eta/R+tau epsilon)^3 d eta
 = (pi/2)(1+zeta) sqrt(kappa R[e(delta)+o(1)])/tau^(3/2),
```

using `integral_0^infinity sqrt(u)/(1+u)^3 du=pi/8`.
The initial cross boundary term vanishes. At `eta_J` of order
epsilon^2, its terminal contribution is O(sqrt(epsilon)), since

```
g_J [2epsilon^(3/2)sqrt(s_J)/(tau epsilon)^2]
 <= (2/tau^2) sqrt(g_J kappa eta_J/epsilon).
```

Let epsilon tend to zero, then delta to zero. Choosing kappa>K,
R>1, and zeta>0 sufficiently close to their limiting values leaves
a positive identity-feature budget under (12). The same outside-code
bound, shell-deficit payment, and new principal-block filling prove
(11). No stochastic independence between different old queries is used.

For a liminf-realizing family, put T=3c_*/2. Necessarily

```
2K/T + e0/T^2 + (pi/2)sqrt(K e0)/T^(3/2) >= 1.     (14)
```

In particular, if K=0 then

```
e0 >= 9c_*^2/4,                                    (15)
```

improving the all-new-word union consequence
`e0>=9c_*^2/(8log2)`. At e0=0 the coefficient remains `tau>2K`:
this composition argument does not claim an improvement in that limit.

## 12. Positive inner complexity: a sharper one-factorization bound

**Applicability correction:** the hypotheses e0>0 and finite K_es cannot
hold for full near-level sets of any bounded-cap sequence, by Section 24.
The following calculation is mathematically valid but is NOT a viable
original-problem criterion. It is preserved as an audited formal estimate.

When e0>0, one can avoid paying cumulative precision altogether in the
innermost asymptotic regime. Pick a single small fixed outer code
E(eta0), whose factorization density is arbitrarily close to e0.
Its one GS feature block represents EVERY inner query. Retain a small
fixed identity block for the queries outside E(eta0). Because
`e(eta)>=e0` and `e(eta)s(eta)<=kappa eta`,

```
s(eta)<=kappa eta/e0   (0<eta<=eta0).               (16)
```

Thus one common variance proxy suffices, and the only remaining
optimization is scalar. Put r=sqrt(e0), k=K/r. Define

```
w(k)=(k+sqrt(k^2+4))/2,
F_chi(K,e0)=r [ w(k)/2 + log(w(k))/k ]  (K>0),
F_chi(0,e0)=r.
```

Then (11) holds for every `tau>F_chi(K,e0)`.

To verify the constant, in the ideal limiting parameters, a code at
window `eta=epsilon r z` has entropy rate at most `epsilon k z`.
The proposed bridge allowance is `epsilon r (t+z)`, where t=tau/r.
Equation (9) and the code union therefore require

```
J(t+z)>kz  for every z>=0.
```

This is equivalent to

```
t > sup_(w>=1) [w−J(w)/k].                         (17)
```

The derivative vanishes uniquely at `w−1/w=k`; the maximum is
`w/2+log(w)/k`, since `w^2−1=kw`. At K=0 the threshold is simply t>1.
This proves the displayed formula.

For complete finite-parameter accounting, use an outer density g0>e(eta0)
and feature weight 1−theta, giving `r0=sqrt(g0/(1−theta))`. A geometric
shell ratio R replaces k in (17) by `kappa R r0/e0`, and the threshold
by `r0 sup[w−J(w)/(kappa R r0/e0)]`. Strict tau above the limiting
formula permits choosing kappa, eta0, g0, theta, and R in that order
so this remains below tau. The terminal window of order epsilon^2
has entropy O(epsilon^2), while t>1 supplies a strictly positive
tail exponent of order epsilon n. The finite hierarchy and all outer
queries are then paid exactly as before. Fix these parameters and
epsilon before taking n to infinity.

For comparison, the old all-new-word union bound on this SAME single
factorization gives another improved sufficient threshold. Put
`L=sqrt(2(log2)e0)` and define

```
F_union(K,e0)=L                         if 0<=K<=L,
F_union(K,e0)=(K^2+L^2)/(2K)            if K>L.
```

Indeed it is precisely
`sup_(u>=0)[sqrt(2Ku+L^2)−u]`. Thus one may take

```
tau > min(F_chi(K,e0), F_union(K,e0))     (e0>0).   (18)
```

The common-factorization observation was independently suggested by the
discrepancy researcher while this proof was being derived. It is distinct
from the absolute-value/square-MGF improvement. For K/r near zero,
`F_chi/r=1+(K/r)/4+O((K/r)^2)`; for large K/r its leading term is K/2,
but F_union can be smaller. These are extension/landscape criteria, not
new global cap constants. For liminf-realizing signings with e0>0,
both displayed thresholds must be at least `3c_*/2`.

Sections 11–12 received an independent PASS from the discrepancy
researcher, including the square-MGF calculation, pi/2 resource,
terminal accounting, finite-parameter proxy, and both scalar thresholds.
The replay `computations/paper_bernoulli_2026_09_17_abs_columns.py` checks
600 tail inequalities, the pi/8 integral, 12 scalar maximizations, and
finite hierarchy resources. This supplements, rather than replaces,
the proofs above.

## 13. A sharp finite loss example for bounded-subGaussian convex comparison

The general convex comparison in Section 7 cannot have factor one even
when `B=L=1`. Let X be uniform on the even-parity code in `{−1,1}^5`:
`product_i X_i=1`. Its coordinates are centered and bounded by one, and
X is exactly 1-subGaussian:

```
E exp(<theta,X>)
 = product_i cosh(theta_i)+product_i sinh(theta_i)
 <= exp(||theta||2^2/2).                            (19)
```

Here is an elementary all-direction proof. Put `u_i=tanh|theta_i|`
and

```
g(u)=atanh(u)^2/2+log(1−u^2)/2.
```

All coefficients in the power series of g from degree four onward are
positive. The first three give

```
g(u)>=u^4/12+4u^6/45+71u^8/840
     >=(u^4+u^6+u^8)/12 >=u^5/5,    0<=u<1.
```

For the last inequality, its numerator is
`p(u)=5u^4+5u^2−12u+5`. Since p''>=10, Taylor's lower bound at 2/3 gives

```
p(u)>=17/81+(16/27)(u−2/3)+5(u−2/3)^2
     >=701/3645>0.
```

AM–GM now implies `sum_i g(u_i)>=product_i u_i`, which is at least
`log(1+product_i u_i)`. This proves (19), including the worst positive
sign of the product term; a negative product only reduces the MGF.

Nevertheless the SHARP dilation for convex domination by five
independent signs epsilon is

```
X <=_(convex order) (5/4)epsilon,
and no smaller dilation works.                     (20)
```

For the upper bound, given X, set epsilon=X with probability one half;
otherwise flip one uniformly chosen coordinate. The resulting epsilon
is uniform on the entire cube and `E[epsilon|X]=(4/5)X`. Jensen proves
(20). For sharpness use the convex homogeneous support function
`f(z)=max_(x in even parity code)<x,z>`.
One has `E f(X)=5`, whereas `E f(epsilon)=4`: even words give five,
odd words give three. Thus any dilation must be at least 5/4.

The parent independently supplied the simpler dimension-six proof and
the coupling; the dimension-five extension above improves that exact
benchmark. It shows a real loss in ambient convex comparison, not a
counterexample to the theorem or a new extremal quadratic construction.

The exact coupling and sharp support-function ratio, 1001 rational
polynomial checks, and 10000 numerical MGF directions are replayed by
`computations/paper_bernoulli_2026_09_17_parity.py`. The all-direction
proof is the power-series argument, not the numerical check.

## 14. Dense independent-reference universality: an archived scope boundary

An initially fresh calculation here collided with a stronger existing
archive theorem, so it is recorded only as an application. For any
finite T in `[-1,1]^m` with |T|=N and arbitrary deterministic offsets a_t,
ordinary one-coordinate Lindeberg replacement gives

```
| E max_t[a_t+<epsilon,t>] − E max_t[a_t+<G,t>] |
                    <= C m^(1/3)(log N)^(2/3).      (21)
```

Indeed the log-sum-exp approximation at inverse temperature beta has
error at most `log N/beta`. Every third coordinate derivative is bounded
by `8 beta^2`, independently of the offsets. Matching the first two
moments and summing the third-order Taylor errors costs `C m beta^2`;
optimize beta. The case N=1 is an exact equality by centering.

For the actual dense bridge, m=nq with q=floor(epsilon n), and
`T={(sigma x_i y_j):x in a specified code,y in {−1,1}^q,sigma=±1}`.
There are at most `2^(n+q+1)` queries. Hence (21) is `O_epsilon(n^(4/3))`,
which is `o(n^(3/2))`, uniformly over every fixed old-energy offset.

Thus an isotropic replacement of independent physical-edge Gaussian
noise by independent Bernoulli noise has no leading-order cap benefit
for a macroscopic bridge. This does NOT Gaussianize the dependent GS
law, its anisotropic feature process, or a one-vertex extension.

Archive collision:
`artifacts/flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`
already proves a more general all-offset bounded-feature comparison;
`artifacts/flatify_independent_2026_09_07_iid_bridge_control_audit.md`
contains the direct third-derivative replacement. No novelty is claimed
for (21). The new-spin composition gains in Sections 11–12 avoid this
boundary because they use the actual GS column laws and their targeted
variances, not merely this independent reference swap.

## 15. Sparse-in-basis isotropic full-sign columns: a stronger cover theorem

New application proposed by the parent and reconstructed here; independent
audit pending. This uses an explicit full-sign law, not an abstract
subGaussian-existence assumption. Its advantage is in mean absolute
responses, not in covariance contraction.

Fix k admitting a real sign Hadamard matrix H_k. Along orders n=kp,
suppose H_p is also a sign Hadamard matrix. Let h_a and v_b denote their
rows, and use the ANTIPODAL center family

```
F_n={ ±h_a tensor v_b : 1<=a<=k, 1<=b<=p }.
```

The antipodal inclusion is essential because the absolute-energy codes
E_n(eta) are antipodal. Generate one random full-sign column by

```
Z=h_L tensor g,
L uniform in {1,...,k},  g uniform in {−1,1}^p,
L and g independent.
```

Then EZ=0 and `Cov(Z)=I_n`. Indeed, for coordinates (a,j),(b,l),
independence of g gives zero when j differs from l, and Hadamard row
orthogonality gives `E_L h_(L,a)h_(L,b)=1_(a=b)`.
For any deterministic u,

```
E[exp(<u,Z>)|L]
 <=exp((1/2)sum_j(sum_a u_(a,j)h_(L,a))^2)
 <=exp(k||u||2^2/2).
```

Thus Z has subGaussian variance proxy k in ambient Euclidean coordinates.
For every center f, however,

```
<f,Z> = ±k 1_(L=a) sum_j (v_b)_j g_j,
E|<f,Z>|=a_p:=E|sum_(j=1)^p g_j|
          =(sqrt(2/pi)+o(1))sqrt(p).                (22)
```

The large response is activated with probability 1/k. Its variance is
still n; covariance alone does not detect the gain in (22).

### 15.1 Cover hypothesis and exact extension statement

For a sequence of actual full-sign matrices A_n with bounded normalized
caps, define the limiting cover radius

```
r(eta)=limsup_n n^(-1)
       max_(x in E_n(eta)) min_(f in F_n) d_H(x,f).
```

Assume r(eta)<1/2 for sufficiently small positive eta and

```
K_r=limsup_(eta down to0) r(eta) h(r(eta))/eta < infinity,
```

where h is binary entropy with natural logarithms. Then, for every

```
tau > sqrt(2/pi)/sqrt(k) + 2k K_r,                  (23)
```

and every sufficiently small fixed epsilon>0, there exist actual
full-sign extensions preserving A_n as an exact principal block and
satisfying the extension conclusion (11).

In particular `r(eta)h(r(eta))=o(eta)` permits every
`tau>sqrt(2/pi)/sqrt(k)`. The condition holds, for example, when
`r(eta)=O(eta)`. This is a sufficient cover hypothesis; it is not claimed
that a current minimizer family satisfies it.

### 15.2 Proof, paying the centers only once

Let q=floor(epsilon n) and use q independent columns with the law (22)
to form C. For every center f,

```
max_y |f^T C y| = sum_(j=1)^q |<f,C_j>|.
```

The centered absolute-value argument of Section 11.1 gives variance
proxy `2kn` per summand. A union over only 2n centers therefore shows,
with probability tending to one,

```
max_(f in F_n,y) |f^TCy|
    <=q a_p+O(sqrt(qkn log n)),
limsup_n n^(−3/2) max_(f,y)|f^TCy|
    <=sqrt(2/pi) epsilon/sqrt(k).                  (24)
```

This is one simultaneous center event, used in every energy shell.
It is never multiplied by the number of levels.

For a word x within Hamming radius rho n of f, write w=x−f.
Then `||w||2^2<=4rho n`. For every fixed new word y,

```
E exp(lambda w^TCy)
                  <=exp(lambda^2 4kq rho n/2).
```

The number of pairs (f,x) at that radius is at most
`2n sum_(l<=rho n) binom(n,l)`, and all 2^q new words are included.
A two-sided union bound, for each fixed rho<1/2, therefore gives

```
n^(−3/2) max_(f,x,y) |(x−f)^TCy|
 <=sqrt(8k epsilon rho[h(rho)+epsilon log2])+o_n(1). (25)
```

The uniform entropy bound is valid for the finite hierarchy below;
if rho=0 the residual is identically zero.

Choose K_r<kappa and a small fixed eta0. Use geometric levels
`eta_j=eta0 R^(−j)`, R>1, down to `eta_J` of order epsilon^2.
For these finitely many levels choose strict monotone radius majorants
rho_j>r(eta_j) with `rho_j h(rho_j)<=kappa eta_j`.
They can be chosen with arbitrarily small extra gaps, preserving all
fixed-level limsup cover statements. In particular rho0 can be made
arbitrarily small by shrinking eta0.

Using `sqrt(a+b)<=sqrt(a)+sqrt(b)`, (25) is at most

```
sqrt(8k epsilon kappa eta_j)
       +sqrt(8k log2) epsilon sqrt(rho0)+o_n(1).    (26)
```

For x in `E_n(eta_j)\E_n(eta_(j+1))`, its old scalar deficit exceeds
`eta_j/R`. The elementary optimization

```
sup_(eta>=0)[sqrt(8k epsilon kappa eta)−eta/R]
                                      =2k kappa R epsilon
```

pays every nonterminal shell. At the terminal level the first term
in (26) is O(epsilon^(3/2)), so no old deficit is needed. Every shell
is combined with the SAME center event (24).

Outside E_n(eta0), the ambient k-subGaussian bound and a union over all
old and new words give bridge cap
`sqrt(2k epsilon(1+epsilon)log2)n^(3/2)+o_n(n^(3/2))`.
For small enough fixed epsilon this is below eta0 n^(3/2), so the
old deficit pays it without invoking a center.

Finally fill the new principal block with a full-sign block of cap at
most q^(3/2). Every bridge entry is already a sign. The simultaneous
events in (24)–(26) and the outside-code event have positive probability
by choosing their finite union failure bounds to sum below one.
Choose kappa near K_r, R near one, and eta0 small so that the coefficient
in (24) plus `2k kappa R+sqrt(8k log2)sqrt(rho0)` is below tau.
Then choose small fixed epsilon and its finite hierarchy, and only then
let n tend to infinity. The terminal O(epsilon^(3/2)) residual can be
absorbed into the strict tau margin for small epsilon; it does not
alter the stipulated coefficient one on the new-block term in (11).
This proves (23).

### 15.3 Why the certificate is different from low-factorization centers

Uniform measure on the orthogonal basis rows has covariance I_n.
Thus `gamma_fac(F_n)^2=n` exactly, by the covariance lower bound and
the identity-factorization upper bound. A low-factorization-center
certificate cannot certify this full center family as cheap. Nevertheless
the explicit law (22) makes its mean absolute response smaller by
1/sqrt(k), while retaining exact covariance I_n and usable tails.

This does not assert that the actual near-level code contains every
center, nor that old methods fail for every subset of these centers.
It establishes a new sufficient cover criterion using this given full
center family. If a liminf-realizing sequence satisfies the cover
hypothesis, necessarily

```
sqrt(2/pi)/sqrt(k)+2k K_r >= 3c_*/2.               (27)
```

In the zero-slope case, any admissible fixed k with
`sqrt(2/pi)/sqrt(k)<3c_*/2` rules out that cover geometry.
For the reported lower bound on c_*, already k=2 meets this numerical
inequality. This is a landscape obstruction and an actual extension
construction under the cover hypothesis, not a new cap constant.

### 15.4 All-order, subexponential-dictionary generalization

The parent observed that H_p is unnecessary: for EVERY Boolean word v
of length p, `E|sum_j v_j g_j|=a_p`. Consequently the same theorem holds
for any nonempty center code G_p in `{−1,1}^p` with
`log|G_p|=o(n)`, using

```
F_n={ ±h_a tensor v : 1<=a<=k, v in G_p }.
```

The center union has subexponential, rather than necessarily polynomial,
size. Its concentration error is
`O(sqrt(qkn[log|G_p|+log n]))=o(n^(3/2))` for fixed k and epsilon.
The residual count becomes `exp(n h(rho)+o(n))`, so every constant in
(23) is unchanged. G_p may be selected from the actual old matrix before
the bridge columns are sampled; the law Z itself does not depend on v.

There is also no divisibility restriction. For arbitrary n put
`p=floor(n/k)` and `ell=n−kp<k`. Append ell independent fair coordinates
xi to Z, and include arbitrary signs on these leftover coordinates in
the center family:

```
Z=(h_L tensor g, xi),
F_n={ (±h_a tensor v,u): a<=k,v in G_p,u in {−1,1}^ell }.
```

The appended vector is still centered, has covariance I_n, and has
ambient subGaussian proxy k. Every center's mean absolute response is
at most `a_p+ell=a_p+O(k)` by the triangle inequality. The extra bridge
center cost `q O(k)=o(n^(3/2))`; the center count gains only a fixed
factor 2^ell. Since p/n tends to 1/k, the leading coefficient remains
`sqrt(2/pi)/sqrt(k)`. Thus (23) is an ALL-ORDER theorem for such
selectable structured center codes. Only the fixed small H_k is needed.
A common coordinate permutation and common coordinate sign-switching,
selected before sampling, may also be applied to both dictionary and
column law without changing any assertion.

The structured-cover premise remains substantive: arbitrary
subexponential Boolean codes are not asserted to admit this common
repeated-block/mode representation. Nor is a near-level cover inferred
from covariance, cardinality, or a few prescribed ground states.

### 15.5 Audit, archive distinction, and a quantitative Walsh corollary

Both other researchers independently audited Sections 15.1–15.3 and
returned PASS; the localization researcher also reran the finite-law
script successfully. The all-order generalization above is an explicit
extension of the same proof. The replay is
`computations/paper_bernoulli_2026_09_17_basis_columns.py`.

Archive checks found adjacent but different mechanisms:

- `principle_invent_2026_09_07_hierarchy_complexity_dilution.md` pays
  low-factorization center covers. It does not make the full orthogonal
  dictionary cheap in mean absolute response.
- `flatify_independent_2026_09_07_mixed_paired_bridge_sector.md` pins a
  fixed physical center and pays specified paired-noise sectors, not
  every center in this orthogonal dictionary and its near-level cover.
- `decisive_bridge_walsh_latent_spike_routing_2026_09_07.md` routes latent
  Walsh spikes in a different ensemble to prove a pressure LOWER
  obstruction. It does not give the current one-column sampling upper
  construction or its extension slope.
- `flatify_construct_2026_09_07_high_rank_template_embedding.md` uses
  exact isotropic sign orbits for adversarial lower witnesses. It does
  not assert that isotropy forces large mean absolute responses.

The block dependence in Z is essential. Its single shared mode L affects
kp physical coordinates, so this is not an independent physical-edge
reference swap covered by Section 14's universality boundary.

For a fixed Sylvester-Walsh dictionary at dyadic orders, the SAME cover
radius r works for each fixed dyadic k dividing the orders eventually.
Applying (27) for each such k gives

```
K_r >= sup_(k dyadic) [(3c_*/2−sqrt(2/pi)/sqrt(k))_+/(2k)].
```

In particular k=4 gives
`K_r>=(3c_*/2−sqrt(2/pi)/2)_+/8`, approximately 0.03138 when the
reported lower bound for c_* is inserted. This is a statement along
liminf-realizing dyadic orders if such a subsequence is under discussion;
it does not assume the global liminf is realized at those orders.
The continuous-k maximization has value `pi c_*^3/8`, but that is only
a relaxation and is NOT asserted to be attained by available Hadamard
orders. The corollary was independently derived by the localization
researcher as well.

### 15.6 A self-contained exclusion independent of the reported lower bound

The zero-slope exclusion already follows at k=8 from an elementary
universal lower bound; it does not need the reported numerical c_*.
For a symmetric hollow full-sign A, let `M(A)=max_x |H_A(x)|` and
`B(A)=max_(x,y) |x^T A y|`. For independent fair y, optimizing x gives

```
B(A)>=E_y sum_i |sum_(j!=i) a_ij y_j|
     =n E|S_(n−1)|=(sqrt(2/pi)+o(1)) n^(3/2).
```

For Boolean x,y put z=(x+y)/2 and w=(x−y)/2. Both lie in the real
cube, and polarization gives
`x^T A y=2[H_A(z)−H_A(w)]`. A multilinear polynomial's absolute
maximum on that cube is attained on its Boolean vertices, so
`B(A)<=4M(A)`. Consequently

```
c_*:=liminf_n M_n/n^(3/2) >=sqrt(2/pi)/4.
```

The fixed Sylvester matrix H_8 exists, and
`sqrt(2/pi)/sqrt(8)<(3/2)sqrt(2/pi)/4`. Therefore a liminf-realizing
sequence cannot have K_r=0 for any of the all-order structured center
families of Section 15.4 with k=8. More quantitatively it must satisfy

```
K_r >=sqrt(2/pi) [3/8−1/sqrt(8)]/16 >0.
```

This remains a conditional restriction on the geometry of an actual
near-extreme set. It does not establish that any minimizer has the
structured cover, nor improve the universal lower bound used in this
paragraph. The all-order finite replay additionally checks arbitrary
short codes and leftover coordinates at (n,k)=(7,2),(11,4),(13,4).

## 16. One mixed column law protects several incompatible frame dictionaries

This extends Section 15 to a finite union of different structured
dictionaries, without separately paying their scalar response bounds.
The dictionaries may use unrelated coordinate permutations and switches.
Their number and mode orders are fixed before taking n to infinity.
Section 17 below improves the residual coefficient by balancing the
sampled modes; the independent-mixture theorem here remains valid but
is not the strongest final certificate.

### 16.1 Uniform mean response away from the protected dictionary

Put `kappa_0=sqrt(2/pi)`. The elementary weighted-Rademacher estimate

```
E|sum_j a_j epsilon_j|
 <=kappa_0 (sum_j a_j^2)^(1/2)
      +C (sum_j |a_j|^3)^(1/3)                  (28)
```

holds for a universal C. Here is a self-contained proof of the required
uniform error. Smooth the absolute value to
`F_delta(t)=sqrt(t^2+delta^2)`. Its third derivative has magnitude at
most `3/delta^2`. Replace each independent sign in turn by an independent
standard Gaussian. Taylor expansion to order two has cancelling first
and second moments, and the total replacement error is at most
`(1+2kappa_0) sum_j|a_j|^3/(2delta^2)`. Since
`|t|<=F_delta(t)<=|t|+delta`, the desired one-sided estimate follows by
choosing `delta=(sum_j|a_j|^3)^(1/3)`; the zero-sum case is immediate.

For ANY Boolean word f, not only a protected center, condition the
Section 15.4 tensor law on its mode L. The coefficients of the remaining
independent signs have magnitude at most k. Their conditional squared
sum is `sigma_L^2<=kn`, and exact isotropy gives
`E_L sigma_L^2=n`. Consequently (28) and Jensen give the UNIFORM bound

```
E|f^T Z|<=kappa_0 sqrt(n)+O_k(n^(1/3)).          (29)
```

The error is uniform in f, the chosen coordinate switching/permutation,
and the short code. Independent leftover coordinates have coefficient
magnitude one and are already included in these estimates. For a
protected center f in its own dictionary F, the sharper earlier bound
is instead

```
E|f^T Z|<=[kappa_0/sqrt(k)+o_n(1)]sqrt(n).       (30)
```

### 16.2 Explicit finite-mixture extension theorem

Fix R dictionaries `F_(r,n)` satisfying Section 15.4, possibly with
different fixed mode orders k_r. Let Z_r be their explicit full-sign
column laws. Draw one label r with probabilities pi_r and then sample
Z=Z_r; use independent samples of this mixed law for the q bridge
columns. This gives an ACTUAL sign-valued column, with exact covariance
I and subGaussian proxy `K=max_r k_r`.

For the union `F_n=union_r F_(r,n)`, define the guaranteed membership
discount

```
d_n(f)=sum_(r:f in F_(r,n)) pi_r (1−1/sqrt(k_r)),
delta=liminf_n min_(f in F_n) d_n(f).
```

Equations (29)–(30) give, uniformly on this entire union,

```
E|f^T Z|<=[kappa_0(1−delta)+o_n(1)]sqrt(n).     (31)
```

There are only subexponentially many centers. The centered-absolute
concentration argument of Section 11, with scalar proxy Kn, therefore
pays the center cost in (31) in ONE event. In particular this is not a
sum of R separately maximized scalar responses.

Define r(eta) using the distance to the UNION F_n and let
`K_r=limsup_(eta down to0) r(eta)h(r(eta))/eta<infinity` as before.
The residual proof in (25)–(26) uses proxy K and the union's still
subexponential size, and is otherwise unchanged. Thus every

```
tau>kappa_0(1−delta)+2K K_r                    (32)
```

gives the actual full-sign extension conclusion (11), with the old
block unchanged and the same unit coefficient on epsilon^(3/2).
All dictionaries and mixing probabilities are fixed from the old matrix
before sampling. The finite hierarchy, outside-code deficit, terminal
absorption, and all-order limits are exactly those of Section 15.2.
For fixed n, choosing pi to maximize the least membership discount is
a finite fractional-cover linear program; this describes an optional
optimization of an already explicit column law, not an assumed rounding
oracle.

In particular R arbitrary dictionaries of the same mode order k, mixed
equally, always have `delta>=(1−1/sqrt(k))/R`. Hence

```
tau>kappa_0[1−(1−1/sqrt(k))/R]+2k K_r.         (33)
```

For TWO possibly incompatible k=4 dictionaries, the zero-cover-slope
coefficient is `3kappa_0/4=0.598413...`. This is below `3c_*/2` at the
reported lower bound and therefore rules out such a zero-slope union
cover on a liminf-realizing sequence. More centers are genuinely allowed
than in a single shared grouping; no common representation of the union
is required. Three k=8 dictionaries have coefficient about 0.62595,
while four k=16 or five k=256 dictionaries have coefficient
`13kappa_0/16=0.648281...`. These numerical consequences use the reported
lower bound, unlike the single-frame k=8 conclusion of Section 15.6.

### 16.3 Scope and a ground-landscape obstruction

The new mixture theorem still needs a cover of the FULL near-extreme
code. Its premise is not inferred from low cardinality, covariance, or
a selected ground-state family. The independent discrepancy-track
artifact, Section 15, constructs a recursive Boolean eigenvector law
for `(J_4−2I_4)^(tensor a)` and proves that every fixed repeated-block
dictionary misses a near-ground word by a positive fraction of n.
We independently audited its complete moment proof and replayed its
exact computations. The key bound for ANY signed matching M of p edges
is

```
E S_M(X)^2 <=(11p^2+32p)/27,
S_M(X)=sum_({i,j} in M) sigma_ij X_i X_j.
```

A partial matching between any two mode blocks therefore gives distance
at least `[1−sqrt(11/27)]n/(2k)−o(n)` from the entire fixed-k family,
even when the inner code is the whole cube and the coordinates are
arbitrarily switched/permuted. This is an actual application obstruction
for the near-half tensor family, not a defect in (32).

For TWO fixed mode geometries this same witness law still gives a
simultaneous obstruction: choose any t with `sqrt(22/27)<t<1`.
Markov's inequality and a union bound show positive probability that
both matching statistics obey `|S_(M_r)|<=t p_r`, since the total bad
probability is at most `22/(27t^2)+o(1)<1`. The selected near-ground word
is then at distance at least `(1−t)n/(2 max(k_1,k_2))−o(n)` from their
UNION. Thus even the explicit two-frame gain does not supply an
unconditional improvement of that familiar near-half construction.
For three or more frames this second-moment argument alone is
inconclusive; no stronger exclusion is asserted here.

## 17. Balanced mode composition removes every mode-order residual loss

The localization-track researcher supplied the following strengthening,
which we independently checked and combine here with the finite-mixture
mean-response estimate. It changes the column sampling rule, not the
geometric cover assumption or scalar center payment.

Allocate deterministic column counts q_r to the R frames, with
`sum_r q_r=q` and `|q_r−pi_r q|<=1`. Within frame r assign deterministic
mode counts q_(r,a), each floor or ceiling of q_r/k_r. The remaining
short-word signs, and all leftover signs, are sampled independently
across columns. Thus the columns are independent but not identically
distributed. Every physical bridge entry is still a sign.

For any fixed old vector w and new Boolean word y, the subGaussian
variance proxy of its full bridge response is the SUM of the actual
mode proxies. Orthogonality of H_(k_r) gives, in frame r's own grouped
coordinates,

```
sum_a q_(r,a) ||sum_b h_(a,b) w_b||2^2
 <=(q_r+k_r) ||w_core||2^2.
```

The independent leftover coordinates contribute exactly
`q_r ||w_leftover||2^2`. A signed coordinate permutation preserves both
norms. Consequently, writing `D=sum_r k_r`,

```
E exp(lambda w^T C y)
 <=exp(lambda^2 (q+D)||w||2^2/2).              (34)
```

There is no requirement that the different frames' groupings agree.
The signs y_j disappear upon squaring their coefficients. No averaging
of a random subGaussian compensator has occurred: (34) is a direct
deterministic sum of the independent short-sign MGF bounds.

For a fixed center, the sum of the scalar column proxies is similarly
at most `(q+D)n`. Centered absolute-value concentration from Section 11
therefore controls the sum of its absolute column responses, and a
single union over the subexponential center family costs o(n^(3/2)).
Its mean still satisfies (31): replacing equal mode proportions by
integer-balanced counts changes each frame's mean sum by only
`O_(k_r)(sqrt(n))`, by conditional Cauchy--Schwarz. The O(n^(1/3))
uniform error in (29), multiplied by q, is also o(n^(3/2)).

Equation (34) replaces the factor k, or K, in the residual variance
calculation by `1+D/q=1+o_n(1)`. Thus the center cost is unchanged,
while the shell cost becomes `2 K_r`. The final combined certificate is

```
tau>sqrt(2/pi)(1−delta)+2K_r,                 (35)
```

with exactly the same actual full-sign extension conclusion and order
of limits as before. For R arbitrary dictionaries of a common mode
order k, the explicit equal allocation gives

```
tau>sqrt(2/pi)[1−(1−1/sqrt(k))/R]+2K_r.        (36)
```

For one dictionary this is
`tau>sqrt(2/pi)/sqrt(k)+2K_r`, improving Section 15's `2k K_r` term.
For two unrelated k=4 dictionaries it is `tau>0.598413...+2K_r`.
This is a joint mode-balancing construction with one scalar center
event, not separately paid channels. Its lack of a currently established
cover hypothesis for an actual minimizing sequence is unchanged.

If the center dictionary is the fixed Sylvester-Walsh basis at dyadic
orders, it has the required representation for EVERY fixed dyadic k
eventually. Therefore (35), first with k fixed and only then k tending
to infinity, strengthens Section 15.5's conditional landscape bound to

```
K_r>=3c_*/4.                                 (37)
```

As before this concerns a liminf-realizing dyadic subsequence if one
is available; it does not assume such a subsequence exists. No growing-k
uniformity, arithmetic interpolation, or unconditional cap improvement
is inferred. The mixture-law and aggregate-Gram reproductions are in
`computations/paper_bernoulli_2026_09_17_mixture_columns.py`.

## 18. Independent audits of the combined positive theorem and its scope

The finite lemmas in this section remain valid. The full near-level
finite-K_s application is now known to be universally impossible for
bounded-cap sequences; Section 24 supersedes any earlier suggestion
that this was merely an unestablished optimizer hypothesis. The separate
anchored Hamming-cover theorem is not ruled out by that obstruction.

The director's `paper_director_spectral_entropy_response_2026_09_17.md`
has been independently reconstructed in full here and passes. Its
Gaussian-observation entropy bound retains the actual Boolean set;
balanced scalar-sign bounded differences give variance proxy `(q+k)n`
for the EXACT absolute new-spin response. Consequently the direct
near-level hypotheses `mu(eta)->0` and
`K_s=limsup s(eta)/eta<infinity` yield extension slope `tau>K_s/2`.
The terminal contribution is absorbed by strict tau, and the outer
deficit and all-order leftovers are correctly accounted for. On a
liminf-realizing family the necessary implication is `K_s>=3c_*`.
No optimizer hypothesis is established by this audit.

The localization artifact `paper_localization_balanced_modes_2026_09_17.md`
also passes independent reconstruction: its entire low-effective-mode
code `R_eff<=u=o(k)` has subexponential size for `k->infinity, k=o(n)`,
maximal squared factorization complexity n, and vanishing scalar mean
response under balanced modes. The Hamming-cover extension slope is
`tau>2K_H`. The finite Gaussian entropy proof and its separate net
cross-check both have the stated constants and limit order.

The actual near-half Hadamard family is nevertheless an obstruction,
not a demonstrated application. Sections 17–18 of the discrepancy
artifact construct a random-involution Boolean eigenvector law and
prove the following, also independently audited here:

- For `n=p^2`, every signed matching of m pairs satisfies
  `E S_M^2<=c_p m^2+(1−c_p)mp/2`, where
  `c_p=(2p−1)/[(p−1)(p−3)]=O(1/p)`.
  Hence every fixed finite union of fixed-order repeated-block
  dictionaries misses a near-ground word by a positive linear distance.
- For any switched/permuted Sylvester frame, the normalized mode-energy
  defect `D_k=k sum_a (||P_aY||^2/n)^2−1` satisfies
  `E D_k<=(k−1)[c_p+(1−c_p)/p]`.
  Thus for `k=o(sqrt(n))`, `u=o(k)`, the full low-mode code is at
  asymptotic Hamming distance n/2 from some near-ground word. Finite
  unions, and suitable summed defect budgets, have the same obstruction.

The Fourier kernel in that eigenvector construction is explicitly
alternating; no ordinary-dot-product self-duality classification is
silently imported. We replayed its exact 1680-atom law, 10560 tested
fourth moments, 100 signed matchings, and exact mode Parseval/defect
checks at k=2,4,8,16. These pass. The remaining mode range from order
sqrt(n) up to o(n) is not excluded by that proof.

## 19. Independent product-mixture proof of the fourth-moment transfer

The canonical Gaussian-transfer theorem is maintained by the director
and localization track in `paper_symmetric_frame_universality_2026_09_17.md`;
its full proof, including the finite 19/12 constant, independently passes
our audit. The following independent proof confirms the
same sharper rate and records exactly why a balanced FOURTH-moment
budget can be summed before taking a query maximum. It is not a
separate claim of novelty for Lindeberg replacement.

Consider finitely many affine queries
`L_t(z)=b_t+sum_i a_(t,i) z_i`, with arbitrary deterministic offsets.
Let N be the number of queries, `B=max_(t,i)|a_(t,i)|`, and
`S_4=max_t sum_i |a_(t,i)|^4`. For inverse temperature lambda set

```
F(z)=lambda^(-1) log sum_t exp(lambda L_t(z)).
```

The fourth coordinate derivative is lambda^3 times the fourth Gibbs
cumulant of a_(t,i). Since its centered fourth moment is at most
16 times its uncentered fourth moment, and its squared variance is at
most its fourth moment,

```
|partial_i^4 F(z)|
 <=19 lambda^3 <|a_(t,i)|^4>_z.                (38)
```

Interpolate the INDEPENDENT driver law by the product of
`nu_s=(1−s)Rad+s Gaussian`, 0<=s<=1. The derivative of E F is the
sum of single-coordinate Rad-to-Gaussian replacements with every
other coordinate distributed according to the SAME product nu_s law.
The first three moments agree. Taylor expansion at coordinate zero
and (38) bound the i-th replacement by

```
(19/24) lambda^3 m_i
 [E|Rad|^4 exp(2lambda B|Rad|)
    +E|G|^4 exp(2lambda B|G|)],
```

where m_i is the Gibbs expectation of `|a_(t,i)|^4` with that coordinate
deleted. This follows because changing a coordinate by z changes any
Gibbs expectation of a nonnegative observable by factors between
`exp(−2lambda B|z|)` and `exp(2lambda B|z|)`.

Assume `lambda B<=1`. The bracket is bounded by a universal constant.
Crucially, reinserting a fresh `Z~nu_s` gives

```
E_Z <|a_(t,i)|^4>_(coordinate Z)
 >=m_i E exp(−2lambda B|Z|)
 >=e^(−2) m_i,
```

by Jensen and `E|Z|<=1`. Thus each deleted-coordinate expectation is
controlled by the expectation under one COMMON complete product law.
The sum over coordinates is now at most `S_4`, inside that common
Gibbs expectation. Integrating the product-mixture derivative proves

```
|E F(Rad)−E F(Gaussian)|<=C lambda^3 S_4.       (39)
```

Softmax differs from the actual maximum by between zero and
`log N/lambda`. Therefore

```
|E max_t L_t(Rad)−E max_t L_t(Gaussian)|
 <=log N/lambda+C lambda^3 S_4,
                    provided lambda B<=1.    (40)
```

For the balanced tensor bridge, each scalar driver has coefficient
at most k, while EVERY Boolean query satisfies
`sum_i a_(t,i)^2<=(q+k)n`. Hence the JOINT bound is
`S_4<=k^2(q+k)n`, rather than the sum of the separate coordinate maxima.
With q=floor(epsilon n), take `lambda=n^(−1/4)k^(−1/2)`; for
`k=o(sqrt(n))`, lambda k tends to zero and (40) gives

```
|E Q(sign bridge)−E Q(mode-conditioned Gaussian bridge)|
       <=C_epsilon n^(5/4) sqrt(k)=o(n^(3/2)). (41)
```

All child offsets are retained exactly. The Gaussian reference keeps
the actual singular bundle orientations and balanced aggregate Gram;
it is NOT the independent-edge Gaussian model. The canonical proof
uses a smart Gaussian/Rademacher path and an exact trapezoid remainder,
which gives the same joint budget by a different argument. Both proofs
explain why sequential separate-coordinate maxima lose an unnecessary
factor of k. The discrepancy-track researcher independently confirmed
the product-mixture comparison and audited the canonical smart path.

The director subsequently removed the coefficient-cap restriction by
an independently audited truncation step. With `L=log(2N)` and
`B_0=(S_4/L)^(1/4)`, delete, separately in each query, the coefficients
larger than B_0. The Rademacher deleted response is deterministically
at most `S_4/B_0^3`; the Gaussian expected absolute deleted maximum is
at most `sqrt(2L S_4/B_0^2)`. Apply (40) to the truncated queries at
`lambda=1/B_0`. All three terms are bounded by
`C S_4^(1/4) L^(3/4)`. If S_4=0 the two maxima have identical
expectations. Thus the all-offset comparison has this universal rate
WITHOUT a coefficient-cap assumption. This is a proof decomposition,
not a modification of physical signs in the final construction. In
particular a critical-scale query class with `S_4<=theta n^3` and
`log N=O(n)` has error `O(theta^(1/4)n^(3/2))`. The director and
localization track maintain its canonical statement.

## 20. Audit of the bounded-coefficient fourth-mass Laplace envelope

The director's focused fourth-paper branch uses Gao--Qian,
[*Fourth-Moment Geometry of Rademacher Sums*, Section 4.4](https://arxiv.org/html/2608.17802#S4.SS4), as the motivation for retaining
fourth-power mass in a scalar Laplace bound. We read and reconstructed
Sections 3.1 and 4.1--4.4. The following coefficient-cap refinement has
an independent elementary proof, so it does not rely on the broader
moment-extremizer theorem imported by that paper. No external novelty
claim is made.

For a finite sign sum with variance v, fourth mass w>0, and coefficient
bound B>0, set `m=floor(w/B^4)` and `r=w−mB^4`. The proposed bound is

```
log E exp(t sum_i a_i epsilon_i)
 <=(t^2/2)(v−mB^2−sqrt(r))
       +m log cosh(tB)+log cosh(t r^(1/4)).     (42)
```

It passes audit. Indeed put
`g_t(z)=t^2 sqrt(z)/2−log cosh(t z^(1/4))`. This is increasing and
concave on z>=0, with g_t(0)=0. For x=|t|z^(1/4), the sign of its
second derivative is the sign of
`3tanh(x)/x−sech(x)^2−2`. To prove this is nonpositive, define

```
N(x)=2x+x sech(x)^2−3tanh(x).
N'(x)=2tanh(x)[tanh(x)−x sech(x)^2].
```

The bracket starts at zero and has derivative
`2x sech(x)^2 tanh(x)>=0`. Therefore N>=0. Concavity now minimizes
`sum_i g_t(a_i^4)`, at fixed fourth mass and cap B^4, by filling m
coordinates to the cap and one to r. Subtract that deficit from
`t^2v/2` to obtain (42). The same packing argument applied to sqrt(z)
shows `mB^2+sqrt(r)<=v`, so the displayed residual variance is valid.

The bound is sharp in the finite-sum closure at feasible parameters:
retain the capped spikes and the residual spike, and realize remaining
variance by increasingly many small signs. To preserve fourth mass
exactly, decrease the residual spike slightly; if r=0, decrease one
capped spike instead. The two scalar second/fourth moment equations
have the same one-spike-plus-flat-tail solution used in the paper,
and the adjusted spike converges from below. At zero remaining variance
the packed vector itself attains equality. Actual-max versus upper-cap
conventions cause no exception for feasible nonzero coefficient data.

Concavity also gives the simpler chord consequence

```
log MGF <=t^2v/2
 −(w/B^4)[t^2B^2/2−log cosh(tB)].              (43)
```

At the critical balanced scale `B~sqrt(n)`, `w~n^3`,
`t=beta/sqrt(n)`, this deficit is of order n for fixed beta>0.
The cap-free one-spike envelope has deficit only of order sqrt(n)
there. This is a scalar coefficient-sensitive improvement; combining
it with a supremum still requires all query entropies and offsets to be
paid. A pressure correction alone is not an upper bound for the
original minimum over full signings.

## 21. Operator-norm obstruction to every Boolean rank-one tensor dictionary

This strengthens the mode-specific stress test, using the independently
audited random-involution law from discrepancy Sections 17--18. The
discrepancy researcher independently checked the argument below.
Let n=p^2=4^a and let Y be that law after its fixed switch into a
Boolean eigenvector of `(J_4−2I_4)^(tensor a)`. Every output is in the
absolute near-level code at every fixed positive window eventually.
For any signed matching of m coordinates it satisfies

```
E S_M^2<=c_p m^2+(1−c_p)mp/2,
c_p=(2p−1)/[(p−1)(p−3)].
```

Choose ANY deterministic coordinate permutation and switching, and
reshape n0=k m<=n selected coordinates into a k by m matrix X, leaving
ell=n−n0 coordinates aside. The full Boolean rank-one dictionary is

```
F={ (c tensor v,z): c in {+-1}^k,
                    v in {+-1}^m,z in {+-1}^ell }.
```

It contains every pattern dictionary C tensor G, regardless of the
cardinalities of C and G or whether C is a Hadamard system. Let
`A=X X^T−mI_k`. Its diagonal is zero. Each off-diagonal entry is a
signed matching statistic on m disjoint pairs of physical coordinates.
The moment theorem therefore yields the exact normalized bound

```
E ||A||_F^2/n0^2
 <=(1−1/k)c_p+(k−1)(1−c_p)p/(2n0).           (44)
```

Also, for EVERY c,v in the dictionary,

```
|c^T X v|/n0 <=||X||op/sqrt(n0),
||X||op^2/n0 <=1/k+||A||_F/n0.                (45)
```

Suppose `k->infinity`, `k=o(sqrt(n))`, and `ell=o(n)`, with
n0/n->1. The right side of (44) tends to zero. Markov's inequality
selects one near-ground output with `||A||_F/n0=o(1)`. Equations
(45) give correlation o(n) with EVERY center in F; arbitrary leftover
signs add at most ell=o(n). Consequently

```
liminf_n n^(-1) max_(x in E_n(eta))
                       min_(f in F) d_H(x,f)>=1/2     (46)
```

for every fixed eta>0. This holds uniformly over the chosen grouping,
permutation, and switches. A fixed finite union is covered by the same
argument and a union bound; more generally it suffices that the sum of
the bounds in (44) tends to zero. There is no Sylvester-mode hypothesis.
Every more finely separable Boolean tensor dictionary is a subset of
such a rank-one dictionary after flattening, so it is excluded too in
this dimensional range.

This does not exclude all possible isotropic sign-column laws. It
specifically obstructs using a near-level cover by these very large
separable tensor families on the explicit near-half signing. It says
nothing unconditional about unknown minimizing sequences. The exact
finite Gram identity and moment bound are replayed in
`computations/paper_bernoulli_2026_09_17_rankone_ground.py`.

A further support-independent consequence applies in the smaller range
`k=o(n^(1/4))`. Put `S=XX^T/n0`, so tr S=1. Eigenvalue-wise,

```
1−tr(sqrt(S))/sqrt(k)
 =0.5 sum_i (sqrt(lambda_i(S))−1/sqrt(k))^2
 <=(k/2)||S−I/k||F^2.                               (46a)
```

By (44), its expectation is `O(k^2/sqrt(n))=o(1)`.
Consequently one near-ground word has
`tr(sqrt(S))/sqrt(k)=1−o(1)`. For EVERY isotropic sign law nu,
the nuclear inequality proved in discrepancy Section 19 gives

```
Phi_nu(X)>=tr(sqrt(S))/sqrt(k)=1−o(1).                (46b)
```

The selected word works simultaneously for all such laws on the chosen
grouping; no support-size or mass condition appears. Thus arbitrary-law
small-response hypotheses also fail on this explicit near-half family
when k=o(n^(1/4)). This stronger law-uniform statement is not asserted
for the remaining n^(1/4)-to-sqrt(n) range.

## 22. Support-free deterministic realization of any isotropic sign-mode law

This closes the atom-count qualification in discrepancy Section 19.
The director independently reconstructed the argument, and the full
localization-track audit returned PASS, including the finite constants,
uniform quantifier, leftovers, and Gaussian-transfer corollary. The result is existential
and quantitative; no efficient method for testing its uniform response
condition is claimed.

### 22.1 Finite simultaneous covariance and response theorem

Let nu be ANY probability law on `{+-1}^k` with `E_nu hh^T=I_k`.
For every integer q>=1 there exist deterministic sign vectors
`h_1,...,h_q` such that

```
||q^(-1) sum_j h_j h_j^T−I_k||op <=4k/sqrt(q),          (47)
```

and SIMULTANEOUSLY for every positive semidefinite k by k matrix S
with trace at most one,

```
q^(-1) sum_j sqrt(h_j^T S h_j)
 <=E_nu sqrt(h^T S h)+8sqrt(pi/2) k/sqrt(q).           (48)
```

There is no restriction on nu's support, its smallest atom, or its
mean vector. The same deterministic label list handles every S.

Proof: initially sample H_1,...,H_q independently from nu. Since
`||H_j||2^2=k`, the exact Frobenius second moment is

```
E ||q^(-1) sum_j H_j H_j^T−I_k||F^2=(k^2−k)/q.         (49)
```

Thus the expected operator error is at most k/sqrt(q). For the
response class write `f_T(h)=||Th||2`, allowing ALL real k by k
matrices T with Frobenius norm at most one. This class is precisely
the class in (48), since `S=T^T T`. Symmetrization gives

```
E sup_T [q^(-1) sum_j f_T(H_j)−E_nu f_T(H)]
 <=(2/q) E_(H,epsilon) sup_T sum_j epsilon_j ||T H_j||.
```

Condition on the sampled H_j. Writing a standard Gaussian as its
independent sign times magnitude, convexity of the supremum and
`E|G|=sqrt(2/pi)` compare the last Rademacher supremum to
`sqrt(pi/2)` times its Gaussian counterpart. The Gaussian processes

```
X_T=sum_j g_j ||T H_j||,
Y_T=sum_(j,l) g_(j,l) (T H_j)_l
```

satisfy
`E(X_T−X_U)^2<=E(Y_T−Y_U)^2`, by the reverse triangle inequality.
Gaussian comparison therefore gives `E sup X_T<=E sup Y_T`.
For completeness, on a finite index set this follows by Gaussian
interpolation of softmax: its derivative is a nonnegative Gibbs
average of the differences of squared increment distances. Approximate
the compact Frobenius unit ball by finite nets and pass to the limit.

The latter supremum is the Frobenius norm of the Gaussian matrix
`sum_j g_j H_j^T`, where now each g_j is a standard k-vector. Its
expected squared Frobenius norm is q k^2. Consequently

```
E sup_(S>=0,trS<=1)
 [q^(-1) sum_j sqrt(H_j^T S H_j)−E_nu sqrt(H^T S H)]
 <=2sqrt(pi/2) k/sqrt(q).                              (50)
```

The supremum is nonnegative because S=0 is included. Markov's
inequality with four times each expectation bound shows that (47)
and (48) hold together with probability at least one half. Fix one
such realization. This proves the deterministic existence claim.

### 22.2 Actual full-sign bridge, including all-order quantitative errors

For n0=kp use those fixed labels and independent fair p-vectors g_j
to form the actual sign columns `C_j=h_j tensor g_j`. For an old
Boolean k by p array X define

```
Phi_nu(X)=E_nu ||h^T X||2/sqrt(n0).
```

The mean of the EXACT new-spin response obeys, uniformly for every X,

```
E max_y |x^T C y|=E sum_j |h_j^T X g_j|
 <=q sqrt(n0) Phi_nu(X)
        +8sqrt(pi/2) k sqrt(q n0).                    (51)
```

This uses conditional Cauchy--Schwarz and (48) with
`S=XX^T/n0`; no Gaussian approximation is needed. For every real old
vector w and every Boolean new word y, the independent scalar signs
give the deterministic MGF certificate

```
E exp(t w^T C y)
 <=exp[t^2(q+4k sqrt(q))||w||2^2/2].                  (52)
```

The same proxy applies to the CENTERED exact response
`sum_j|w^T C_j|` by bounded differences in the independent scalar
signs. In particular no random covariance has been replaced by its
expectation: the labels satisfying (47) are fixed before any g_j is
sampled, and their actual aggregate covariance is bounded in (52).

For arbitrary n take p=floor(n/k), n0=kp, ell=n−n0<k and append
independent fair leftover signs to each column. The mean bound gains
only q sqrt(ell), while (52) holds with the full norm of w unchanged.
For fixed epsilon>0, `q=floor(epsilon n)` and `k=o(sqrt(n))`,

```
q+4k sqrt(q)=q+o(n),
8sqrt(pi/2) k sqrt(q n0)+q sqrt(ell)=o(n^(3/2)).        (53)
```

Therefore EVERY sequence of isotropic sign-mode laws with
`k_n^2=o(n)` has the required actual-column realization, regardless
of support size. The laws may depend on the old signing; the selected
labels must be fixed before inner-sign sampling. Combining (51)--(53)
with the support-independent entropy theorem removes its extra
finite-atom realization premise in precisely the same dimension range.
The sharper anisotropic-noise entropy refinement is maintained by the
director; the realization theorem applies to it unchanged.

Section 24 rules out the previously proposed FULL-near-level finite
entropy-slope application. The realization theorem itself is unaffected
and may instead be used on a low-response CENTER code, with the forced
near-level local variation paid by anchored residual increments.

The same construction has driver coefficient bound k and common
quadratic budget `(q+4k sqrt(q))n`. Hence the canonical symmetric-frame
universality theorem also transfers it, with all child offsets, to
its ACTUAL mode-conditioned Gaussian law throughout `k=o(sqrt(n))`.
It does not transfer to independent physical-edge Gaussian disorder.

## 23. Exact critical-scale response gap on an actual near-ground word

The coherent sector in the director's capped-fourth-mass split cannot
be given a vanishing sign/Gaussian comparison error in general. There
is an exact witness within the actual tensor-column construction.

Let n=k^2, let H_k be a sign Hadamard matrix, and form a k by k old
array whose columns are a signed permutation of its rows:
`X_(b,t)=eta_t h_(pi(t),b)`. For a column with ANY fixed mode a,

```
<X,h_a tensor g_j>
   =k eta_(pi^(-1)(a)) g_(j,pi^(-1)(a)).              (54)
```

Thus for q independent columns, regardless of balanced mode counts,

```
max_y |x^T C y|=qk                 for the sign bridge,
E max_y |x^T C^G y|=sqrt(2/pi) qk   for its Gaussian reference.
```

The exact gap is `(1−sqrt(2/pi))qk`. Every fixed-new-word query has
variance `qk^2=qn` and fourth coefficient mass `qk^4`. At
`q~epsilon n`, `k=sqrt(n)`, the gap is of order n^(3/2), and the
fourth mass is of order n^3. Moreover
`W_4^(1/4)[log(2^{q+1})]^(3/4)` has exactly the same order qk.
This checks the necessity of a coherent-sector loss in the universal
comparison; the coefficient bound is not a technical nuisance.

For Sylvester H_k, choose the array
`X_(u,v)=(-1)^(u dot v)` on F_2^a times itself. It is a +k eigenvector
of the alternating Fourier matrix. Apply the fixed coordinate switch D
from discrepancy Section 17 to both the array and the column law.
The resulting old word is a Boolean eigenword of `(J_4−2I_4)^(tensor a)`
with absolute deficit exactly n after hollowing. Thus this response
witness can lie in an ACTUAL full-sign child's near-extreme code.

This is deliberately a FIXED-old-word bridge response statement.
It does not assert the same gap for the whole parent cap: other old
words and both child offsets can change that maximum. It also does not
claim every selectable frame has this obstruction. Pointwise improved
Laplace bounds remain valid, but they do not by themselves compare
maxima of correlated processes with Gaussian leading coefficient one.

## 24. Universal near-level entropy and response rigidity: applicability correction

The discrepancy researcher found the entropy obstruction below; the
director then derived the response strengthening. Both have been
independently reconstructed here. They change the applicability verdict,
not the correctness of the finite information or rounding lemmas.

### 24.1 Forced local entropy makes finite full-code slopes impossible

Suppose actual sign matrices satisfy `Q(A_n)/n^(3/2)->c`, with
`0<c<infinity`. Choose a Boolean maximizer x and a polarity sigma with
`sigma H_A(x)=Q(A)`. Flip a uniformly random set of exactly r coordinates,
producing Y. Every pair product has the same multiplier, so EXACTLY

```
E[Q(A)−sigma H_A(Y)]
 =4r(n−r)Q(A)/[n(n−1)].                              (55)
```

The random deficit on the left is nonnegative. Fix a small positive
eta and `rho<eta/(4c)`, and take r=floor(rho n). Markov's inequality
leaves a positive fraction, independent of n, of that Hamming sphere
inside the absolute near-level set E_n(eta). Therefore

```
liminf_n n^(-1) log|E_n(eta)|>=h(rho).
liminf_(eta down to0)
 [liminf_n n^(-1)log|E_n(eta)|]/[eta log(1/eta)]
                                                >=1/(4c).  (56)
```

To obtain the second line, use
`rho=(1−delta)eta/(4c)` and send delta to zero last.
For bounded but nonconvergent normalized caps, replace c by their
finite limsup. Consequently `K_s=limsup s(eta)/eta` is ALWAYS infinite
for the full near-level code in this class. The direct finite-entropy-
slope augmentation criterion, although formally proved, has NO original-
class application. It is not merely awaiting a favorable optimizer.

Likewise, if a factorization profile has e0>0, then
`limsup e(eta)s(eta)/eta=infinity`. This rules out the positive-inner-
complexity regimes of Sections 11--12, including their proposed
zero-K_es boundary, for bounded-cap families.

### 24.2 Every isotropic law sees square-root response growth

Let nu_n be ANY isotropic sign law in dimension n, and define
`f_n(z)=E_(h~nu_n)|h dot z|`. This function is convex and 1-Lipschitz
in Euclidean norm, by Cauchy--Schwarz and isotropy. Independently flip
each coordinate of the same ground word with probability rho, where
0<rho<1/2. The exact expected energy deficit is now
`4rho(1−rho)Q(A_n)`.

For every product distribution on the Boolean cube,

```
Var f_n(Y)<=4.                                        (57)
```

Indeed the one-sided Efron--Stein inequality bounds variance by the
sum of squared positive coordinate-resampling differences. Choose one
subgradient g at Y; ||g||<=1. Convexity gives
`(f(Y)−f(Y with bit i flipped))_+<=2|g_i|`.
Their squared sum is at most four, even before averaging. Resampling
can also leave the bit unchanged and does not enlarge this bound.

Uniformly for EVERY fixed h, the centered random sum h dot Y has
variance `4rho(1−rho)n`. Its mean may depend on h. Smooth the absolute
value by `sqrt((a+z)^2+delta^2)` and replace the centered biased signs
by matching Gaussians. The third-derivative bound is uniform in the
shift a, and the total third absolute moments are O(n). Taking
delta=n^(1/3) yields a uniform O(n^(1/3)) comparison error. A shifted
Gaussian has smallest mean absolute value at zero shift. Consequently

```
E_Y f_n(Y)/sqrt(n)
 >=2sqrt(2/pi)sqrt(rho(1−rho))−o(1).                   (58)
```

By (57), the response is within o(sqrt(n)) of its mean with probability
tending to one. If `4rho(1−rho)c<eta`, the near-energy event has a
fixed positive probability by Markov. The events therefore intersect.
Optimizing rho, and then taking arbitrarily small strict slacks, proves

```
liminf_n max_(y in E_n(eta)) f_n(y)/sqrt(n)
 >=sqrt(2/pi) sqrt(eta/c),          0<eta<c.            (59)
```

The laws may be chosen from the old signing and may vary arbitrarily
with n. The estimate is uniform in them. It follows that whole-code
`mu(eta)^2=o(eta)` is also universally impossible. For a tensor-mode
law, its scalar absolute response is at most the array response Phi;
thus the corresponding whole-code small-Phi profiles are obstructed
as well. In particular the proposed
`mu(eta)^2 log(1/mu(eta)^2)=O(eta)` application cannot occur.

### 24.3 The viable remaining scope is anchored increments

The forced entropy comes from a local cloud around a ground word.
Its radius is order eta, and a bridge INCREMENT over that cloud has
variance proportional to that radius, rather than the whole-word
variance qn. Therefore `r h(r)/eta` is of order
`eta log(1/eta)` when r=O(eta), and tends to zero. The anchored
center-cover theorems of Sections 15--17 are NOT made vacuous by
(55)--(59). They still require an actual cover of the complete near-
level set; the explicit Hadamard stress tests show that several natural
candidate centers fail. The new low-response realization theorem can
be applied to CENTER codes, with all residual local variation paid
separately by its genuine small increment variance. It must not be
advertised as a rescued full-code entropy-slope theorem.
