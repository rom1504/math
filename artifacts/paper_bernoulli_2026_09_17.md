# Bernoulli decomposition and margin-preserving sign completion

2026-09-17. Initial primary mapping frozen before consulting historical
gamma2/augmentation artifacts. Status: derived comparison theorems,
counterexamples, and stronger landscape/extension criteria. Sections
11–12 received independent mathematical PASS; individual audits are
recorded below. No improvement to the original cap interval or proof
of convergence is claimed.

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

This section is a derived extension, pending independent audit. Its
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
