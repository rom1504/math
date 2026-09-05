# Ordinary entropy support criterion: independent audit and primary LDP scope

Date: 2026-09-05. Status: the support criteria and uniform noise-cloud estimate
are verified. No lower-tail LDP or entropy limit for the original signing
problem is imported or proved.

This report supplements the director's
[proof and archive correction](fresh_entropy_support_audit_2026_09_05.md).
The present researcher read that proof in full and independently reconstructed
its atom-floor and Bernstein/shell arguments. The old Section 5 claim in
`good_signing_entropy_threshold.md` was too broad and has now been corrected
by the director. No edit to that old file was made by this researcher.

## 1. Exact criteria and constants checked

Let N=binom(n,2), let A have independent uniform signs above its zero diagonal,
and write

```math
q(A)=\max_x\left|\sum_{i<j}a_{ij}x_ix_j\right|,
\quad c_n=\min_A q(A)/n^{3/2},
\quad Z_n(c)=\#\{A:q(A)\le cn^{3/2}\}.
```

Every nonempty event has probability at least 2^(-N). Therefore

```math
I_n(c)=-n^{-2}\log\Pr(q(A)\le cn^{3/2})
\in[0,(N/n^2)\log2]\cup\{+\infty\}.
```

An empty event has infinite rate, not the maximal finite rate (log 2)/2.
If liminf c_n < limsup c_n, every threshold strictly between them alternates
infinitely often between empty and nonempty events. Consequently an exact
extended limit of I_n(c) on a dense set of thresholds forces c_n to converge.
This does not require rate differentiability, uniform convergence, a good
rate function, or exclusion of a phase containing only one switching orbit.

For the stronger ordinary nonnegative-entropy criterion, set
S_n(c)=n^(-2) log(1+Z_n(c)). Given any seed A_0 satisfying
q(A_0)<=c_0 n^(3/2), reverse each edge independently with fixed probability
p in (0,1/2), obtaining B. Then

```math
H_B(x)=(1-2p)H_{A_0}(x)+W_x,
\qquad \operatorname{Var}(W_x)=4p(1-p)N.
```

The centered summands of W_x are independent and bounded in absolute value
by 2. There are 2^(n-1) projective spin vectors and two objective signs.
Bernstein therefore gives exactly

```math
\Pr\left(\max_x|W_x|\ge t n^{3/2}\right)
\le 2^n\exp\left\{-\frac{t^2 n^3}
 {8p(1-p)N+(4/3)t n^{3/2}}\right\}.
```

For fixed t>2 sqrt(p(1-p) log 2), this is exponentially small in n,
uniformly over all seeds and all c_0. The shell of flip counts satisfying
|K-pN|<=N^(3/4) has probability tending to one. A particular outcome in
this shell has probability at most

```math
\exp\{-Nh(p)+N^{3/4}|\log(p/(1-p))|\}.
```

The good-event/shell intersection has probability at least 1/2 for all
sufficiently large n, uniformly over seeds. Thus

```math
\log Z_n((1-2p)c_0+t)
\ge Nh(p)-N^{3/4}|\log(p/(1-p))|-\log2.
```

This verifies the constant in Theorem 2.1 of the archived
`entropic_franz_parisi_bernoulli.md`; its uniformity is adequate for
arbitrary seed sequences. For a completely explicit fixed relaxation
epsilon>0, take

```math
p=\min\{1/4,\varepsilon^2/(64\log2)\},\qquad t=\varepsilon/2.
```

Then t>2 sqrt(p(1-p) log 2), and for c_0>=0 the new cap is at most
c_0+epsilon. In particular, eventually and uniformly in c_0>=0,

```math
Z_n(c_0)>0\quad\Longrightarrow\quad
S_n(c_0+\varepsilon)\ge h(p)/4>0.
```

If S_n has a limit on a dense set and c_n oscillates, choose
liminf c_n<c_0<c_1<limsup c_n with c_1 in that dense set. The preceding
implication gives positive limsup S_n(c_1), while infinitely many empty
events give liminf S_n(c_1)=0. This is the contradiction. The ordinary
n^2 scale is sufficient; the previously sufficient n^(3/2) log n scale
is not required by the support argument.

The noise parameter is fixed before n tends to infinity. No claim here
is uniform as p tends to zero with n. The argument creates entropy above
the seed cap; it does not improve the seed's value.

## 2. Even a lower-interval LDP is enough

A weaker imported theorem than a full compact-space LDP would suffice.
It is enough to have one extended function J such that the laws of
q(A)/n^(3/2) satisfy the usual LDP upper bounds on closed lower intervals
and lower bounds on open lower intervals, at speed n^2.

Indeed, if a=liminf c_n<b=limsup c_n, choose a<u<v<b. The event
(-infinity,u] is nonempty infinitely often, so the atom floor and the
LDP upper bound give inf_{z<=u} J(z)<infinity. The same infimum over
(-infinity,v) is finite. The LDP lower bound then forces that larger
event to be nonempty eventually, contradicting limsup c_n>v.

This version needs neither compactness nor attainment of an infimum.
Alternatively clip the observable above any fixed K>1/2 and ask for a
standard LDP on [0,K]. The clipping retains the eventual support minimum.
Do not confuse a pointwise CDF-rate limit with an arbitrary LDP at a
discontinuity of its rate; the open/closed-interval proof avoids that issue.

## 3. Primary results actually checked

### 3.1 Gaussian Ising ground states: a very recent upper-tail theorem

Hong-Bin Chen, Alice Guionnet, Justin Ko, Bertrand Lacroix-a-Chez-Toine,
and Jean-Christophe Mourrat, *One-sided large deviations for the ground-state
energy of spin glasses*, March 6, 2026,
[arXiv:2603.06368](https://arxiv.org/abs/2603.06368),
[author manuscript](https://perso.ens-lyon.fr/jean-christophe.mourrat/gs_ldp.pdf).

The model is a **Gaussian** mixed p-spin field with covariance
n xi(sigma dot tau/n), plus deterministic external field. Its ground-state
variable L_n is the maximum over Ising spins divided by n, not an absolute
maximum. Theorem 1.2 proves, for r>=the typical ground-state value,

```math
\lim_n -n^{-1}\log\Pr(L_n\ge r)=\Lambda^*(r).
```

This is the upper tail at speed n. The introduction explicitly separates
the lower tail and says its zero-field speed is expected to be n^2,
citing the spherical result below. This is not our required lower-tail
Rademacher LDP. In particular, convergence of the usual spin free energy
or its fixed positive fractional moments does not supply the signing-space
lower-tail rate in Section 1.

### 3.2 Spherical models: an n^2 lower-tail bound, not the desired exact rate

Brice Huang and Mark Sellke, *A Constructive Proof of the Spherical Parisi
Formula*, [arXiv:2311.15495v3](https://arxiv.org/pdf/2311.15495), May 1, 2024.

The spins range over the radius-sqrt(n) sphere, and the disorder is Gaussian.
Theorem 5.8, equation (5.10), proves that a fixed lower deviation below
the zero-external-field ground-state value has probability at most
exp(-C(xi,epsilon)n^2), asymptotically. It does **not** assert existence
of an exact n^2 lower-tail rate there. Its equation (5.9) concerns finite
speed-n cost above that threshold. The exact rate function elsewhere in
the paper is for the upper tail of 1RSB models.

Thus neither the spin domain, disorder law, nor exact lower-tail-limit
hypothesis matches our problem. An exponential upper bound alone does
not trigger the finite-atom convergence criterion.

### 3.3 Rademacher largest eigenvalue: exact LDP at the other speed

Alice Guionnet and Jonathan Husson, *Large deviations for the largest
eigenvalue of Rademacher matrices*,
[final author manuscript](https://perso.ens-lyon.fr/aguionne/LDPRad.pdf),
[arXiv:1810.01188](https://arxiv.org/abs/1810.01188).

The final author's Theorem 1.4 assumes sharp sub-Gaussian moment-generating
functions and bounded support or a uniform log-Sobolev inequality. It
gives a speed-n LDP for lambda_max(A/sqrt(n)), with rate infinity below 2
and rate

```math
I(r)=\frac12\int_2^r\sqrt{x^2-4}\,dx\quad(r\ge2).
```

It also treats -lambda_min. Rademacher entries satisfy the assumptions;
the stated diagonal variance can be arranged with bounded diagonal entries
and then deleted at a deterministically vanishing operator-norm cost.

The rate being infinite below 2 means superexponential decay at speed n,
not an empty event. A signing atom costs exp(-Theta(n^2)), so the finite
atom floor is not bounded at this theorem's speed. No n^2 lower-tail
formula for the Boolean norm follows.

### 3.4 Exact small operator norm: GOE, with a demonstrable discrete mismatch

Antoine Maillard, *Average-case matrix discrepancy: satisfiability bounds*,
[arXiv:2410.17887](https://arxiv.org/pdf/2410.17887),
[published article](https://doi.org/10.1002/rsa.70033).

Proposition 2.1 gives the exact GOE lower-tail rate

```math
\lim_{d\to\infty}-d^{-2}\log\Pr(\|W\|_{op}\le\kappa)
=-\frac{\kappa^4}{128}+\frac{\kappa^2}{8}
-\frac12\log(\kappa/2)-\frac38
\quad(0<\kappa\le2),
```

and zero for kappa>2. Here W is normalized GOE with limiting spectrum
[-2,2]. The proof uses the constrained Gaussian beta-ensemble and its
spectral-measure variational problem. The paper's discrepancy model uses
independent GOE matrices W_i and chooses signs in their sum; this is not
the deterministic cut-incidence matrix of our feasibility system.

There is an elementary obstruction to importing this Gaussian rate for
sign matrices: for every hollow sign A,

```math
\|A/\sqrt n\|_{op}^2\ge n^{-1}\operatorname{Tr}((A/\sqrt n)^2)
=1-1/n.
```

Every fixed kappa<1 is therefore eventually impossible for Rademacher
sign matrices, while the displayed Gaussian rate is finite. Universality
of the upper tail at speed n cannot repair this n^2-scale support mismatch.
Also q(A)<=n||A||op/2 is only one inequality: even a sign-matrix spectral
LDP would need further Boolean information to imply a q-LDP.

### 3.5 Recent spectral-measure LDPs remain in a sparse regime

Fanny Augeri, *Large deviations of the empirical spectral measure of
supercritical sparse Wigner matrices*, Advances in Mathematics 466 (2025),
110156, [arXiv:2401.11925](https://arxiv.org/pdf/2401.11925),
[publication](https://doi.org/10.1016/j.aim.2025.110156).

Theorem 1.4 studies bounded centered entries multiplied by independent
Bernoulli(p) sparsifiers, normalized by sqrt(np), under the hypotheses
p->0 and np/log n->infinity. It proves a weak spectral-measure LDP at
speed n^2 p, with a rate obtained from quadratic-vector-equation measures
of nonnegative kernels. Theorem 1.5 gives the analogous sparse
Erdos-Renyi result. The assumption p->0 excludes dense Rademacher Wigner.
Moreover weak spectral-measure information is not automatically a small
operator-norm event or a Boolean ground-state observable.

For historical calibration only, Guionnet's
[Bernoulli random matrices survey](https://perso.ens-lyon.fr/aguionne/ECM_proc3.pdf),
Section 4.2, explicitly listed dense Bernoulli spectral-measure LDPs as
open at the time. This old survey is not a certificate of current status;
the directly checked 2025 theorem above supplies the precise newer scope.

### 3.6 Discrete feasibility entropy: an exact theorem for independent rows

Emmanuel Abbe, Shuangping Li, and Allan Sly, *Proof of the Contiguity
Conjecture and Lognormal Limit for the Symmetric Perceptron*,
[arXiv:2102.13069](https://arxiv.org/pdf/2102.13069).

The model has an m by d matrix G with i.i.d. uniform sign entries, m=floor(alpha d),
and constraints |(Gx)_j|<=kappa sqrt(d). Theorem 2.1, for
0<alpha<alpha_c(kappa), proves that its number of feasible x divided by its
expectation converges to a positive lognormal random variable. In particular
the normalized logarithmic count has a limit in probability in that regime.

The relevant missing mapping is substantial: our d=binom(n,2) sign
coordinates are tested by **all** 2^(n-1) deterministic rows
(x_i x_j)_(i<j), with their exact overlap geometry. They are not alpha d
independent random rows, and our cap divided by sqrt(d) grows as sqrt(n).
The theorem is useful evidence that discrete feasibility entropy can control
support, but it cannot be applied by treating these cut rows as independent.

## 4. What the search does and does not establish

No checked primary theorem supplies a genuine n^2 lower-tail LDP or dense
threshold entropy limit for the Rademacher Boolean absolute ground state.
This is a scoped negative search result, not a theorem that such an LDP
cannot exist, nor a claim that ordinary entropy is inherently inadequate.

The clean noncircular landing target is now the labeled entropy limit
S(c)=lim_n n^(-2) log(1+Z_n(c)) on a dense set. It is unnecessary to
compute S, determine its critical value in advance, show strict positivity
at the critical point, or classify exceptional switching orbits. If the
limits exist, the noise-cloud estimate identifies the limiting optimum as
the boundary between zero and positive entropy, using thresholds on either
side of that boundary. Alternatively, the lower-interval LDP in Section 2
is already sufficient without the noise cloud.

Potential methods still have to retain the discrete support and the full
correlation geometry of the cut tests. The demonstrated failures above are
failures of particular Gaussian, spectral, or independent-row mappings;
they are not general no-go theorems for an entropy/LDP approach.
