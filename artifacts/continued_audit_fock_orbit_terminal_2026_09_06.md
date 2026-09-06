# Independent audit: Gaussian-Fock orbit terminal estimate

Date: 2026-09-06. Independent reconstruction by the audit agent.
Status: the finite theorem below passes. Its application to a particular
recursive construction still requires the conditional-randomness and
normalization checks listed in Section 5. No original cap improvement is
asserted here.

## 1. Exact normalization

Let G be the full group of signed permutation matrices on R^m. For t>=0
use the Gaussian feature vector in the symmetric tensor Fock space

```math
\phi_t(v)=e^{-t\|v\|^2}\bigoplus_{d\ge0}
                  \frac{(2t)^{d/2}}{\sqrt{d!}}v^{\otimes d}.
```

It has norm one and inner product
`<phi_t(v),phi_t(w)>=exp(-t||v-w||^2)`. Orthogonal matrices act
unitarily on every degree and carry `phi_t(v)` to `phi_t(Uv)`.
Write P_G for the orthogonal projection onto the G-invariant subspace and

```math
P_t(v)=\|P_G\phi_t(v)\|^2
       =\mathbb E_{g\in G}e^{-t\|v-gv\|^2},
\qquad L_t(v)=\sqrt{P_t(v)}.
```

If a_i=|v_i| and
`K_t(a,b)=exp[-t(a^2+b^2)] cosh(2tab)`, averaging the independent signs
and then the permutation gives exactly

```math
P_t(v)=\frac{\operatorname{per}(K_t(a_i,a_j))_{i,j=1}^m}{m!}.
```

In particular there is no missing `2^m`, and the factorial is `m!`.
Also `P_t(v)>=exp(-4t||v||^2)` simply by the diameter of the orbit.
The stronger lower bound `exp(-2t||v||^2)` follows from `cosh>=1`,
but is not needed below.

## 2. Exact orbit-covariance operator norm

Set

```math
C_v=\mathbb E_{g\in G}
          |\phi_t(gv)\rangle\langle\phi_t(gv)|.
```

Its nonzero eigenvalues are the eigenvalues of the finite orbit Gram
matrix `exp[-t||gv-hv||^2]/|G|`. This remains true when orbit points repeat.
The Gram matrix is entrywise nonnegative, symmetric, and every row sum
is P_t(v). Its all-ones eigenvector and its row-sum norm therefore give

```math
\|C_v\|_{\mathrm{op}}=P_t(v).
```

This equality, rather than a trace-one bound, is the key estimate.

## 3. A uniform finite-dimensional truncation

Let Pi_D project onto Fock degrees at most D. Both rotations and P_G
commute with Pi_D. Thus, for every orthogonal U,

```math
\mathbb E_g L_t(Ugv)^2
 =\operatorname{Tr}(P_G\rho(U)C_v\rho(U)^*)
 \le r_{m,D}P_t(v)
       +\Pr\{\operatorname{Pois}(2t\|v\|^2)>D\},
```

where `r_{m,D}=rank(P_G Pi_D)`. To justify the tail, split the positive
trace into its orthogonal low- and high-degree blocks. The low block is
bounded by its rank times `||C_v||`; the high block by the full
high-degree mass of C_v. Cross-degree blocks contribute zero to this
trace. The squared degree masses of phi_t(v) have the displayed Poisson
law.

A degree-d signed-permutation invariant symmetric tensor corresponds to
a homogeneous invariant polynomial. Odd degrees have dimension zero;
degree 2j consists of symmetric polynomials in `v_1^2,...,v_m^2` and has
dimension p_m(j), the number of partitions of j with at most m parts.
Consequently

```math
r_{m,D}=\sum_{j\le D/2}p_m(j)
 \le\sum_{j\le D/2}p(j)
 \le\exp\{\pi\sqrt{D/3}\}.
```

The last bound follows, without an asymptotic partition formula, from
the partition generating function and
`log prod_{k>=1}(1-e^{-sk})^{-1}<=pi^2/(6s)`, optimizing s. The D=0
case is interpreted as r=1.

Fix t,C>0 and assume `||v||^2<=Cm`. Put `D=ceil(8tCm)` and
`lambda_max=2tCm`. The Poisson Chernoff bound at exponent log 4 gives

```math
\Pr\{\operatorname{Pois}(2t\|v\|^2)>D\}
 \le e^{-(4\log4-3)\lambda_{\max}}
 \le e^{-2\lambda_{\max}}
 \le P_t(v).
```

Therefore the explicit uniform result is

```math
\boxed{\quad
\mathbb E_g L_t(Ugv)
 \le \sqrt{1+e^{\pi\sqrt{\lceil8tCm\rceil/3}}}\ L_t(v).
\quad}                                                     (1)
```

The first moment follows by Cauchy--Schwarz only after the relative
second-moment estimate has been established. For fixed t,C its factor
is `exp(O_{t,C}(sqrt(m)))`. The t=0 or C=0 cases are immediate.
This is not the invalid exponential-scale replacement of
`E sqrt(P)` by `sqrt(E P)` used in a global annealed first moment.

## 4. Removing a diagonal coordinate

The diagonal of the even Gaussian kernel is
`K_t(a,a)=(1+exp(-4ta^2))/2>=1/2`. Restrict the full permanent sum to
permutations fixing coordinate i. If P_del,i denotes the normalized
permanent after deleting that row and column, then exactly

```math
P_t(v)\ge \frac{K_t(a_i,a_i)}m P_{\mathrm{del},i}
          \ge\frac1{2m}P_{\mathrm{del},i}.
```

Hence `max_i sqrt(P_del,i)<=sqrt(2m)L_t(v)`. This applies to the
one-coordinate deletion in the weave estimate; no sum over the m
deletions is needed.

## 5. Exact scope of recursive use

For a sequence `v_{j+1}=U_j g_j v_j`, (1) iterates conditionally if,
given the past, g_j is fresh uniform on the FULL signed-permutation
group, and U_j is fixed or is sampled independently of g_j. Orthogonal
steps preserve the norm bound. A fixed number of full-size steps costs
only `exp(O(sqrt(m)))`; more generally the accumulated costs must be
summed explicitly at the dimensions of the steps.

The theorem does not by itself cover a subgroup randomization, a U_j
chosen after observing g_j, a transform with changed amplitude
normalization, or exponentially many separately charged blocks. A
construction using paired blocks must specify how the full orbit law
appears conditionally and how blockwise partition sums are assembled.
It also does not count input selectors or spin rows. Those entropy
factors remain part of the recursive variational calculation.

The finite PSD permanent product inequality and its finite-type limit
are compatible with this theorem, but are separate statements. In
particular, a product of coordinatewise even Gaussian kernels is not
the same kernel as averaging a single shared global sign of a vector.
An application must declare which sign group is being averaged.

## 6. Independently verified permanent tensorization and type limit

For any two m-by-m positive-semidefinite Gram matrices A and B,

```math
\frac{\operatorname{per}(A\circ B)}{m!}
 \ge\frac{\operatorname{per}A}{m!}
       \frac{\operatorname{per}B}{m!}.                  (2)
```

Indeed write `A_ij=<v_i,v_j>`, `B_ij=<w_i,w_j>`, and let u and w
be their m-fold tensor products. On the regrouped tensor space the
simultaneous permutation average P_diag is an orthogonal projection.
Its range contains the range of the product of the two separate
symmetrizing projections. Evaluate this projection inequality on
`u tensor w`. The three quadratic forms are precisely the three
normalized permanents in (2). This proof does not assert the stronger,
generally different unnormalized inequality without the factorial.

For a probability law mu on a finite subset of R^s define

```math
F_t(\mu)=\inf_{\pi\in\Pi(\mu,\mu)}
 \left\{D(\pi\Vert\mu\otimes\mu)
                    +t\int\|x-y\|^2\,d\pi(x,y)\right\}.
```

Apply (2) repeatedly to the coordinate Gaussian Gram matrices of a
word of joint types with empirical distribution approaching mu. Their
Hadamard product is the Gaussian Gram matrix of the joint vector.
The exact finite-table permanent formula and fixed-type Stirling limit
give normalized log permanents `-F_t(mu)` and `-F_t(mu_r)` respectively.
Therefore, with no missing entropy or factor of two,

```math
F_t(\mu)\le\sum_{r=1}^sF_t(\mu_r).                  (3)
```

This extends to every finite-second-moment probability law by finite
W2 quantization. To check the passage directly, couple mu to nu by a
transport kernel and apply that kernel independently to the endpoints
of an almost-optimal self-coupling of mu. Relative entropy contracts.
If the endpoint displacement has L2 norm delta=W2(mu,nu), the L2 norm
of the difference of the endpoint differences is at most 2delta, so
the quadratic-cost change is bounded by
`4delta sqrt(cost_mu)+4delta^2`. Since
`sqrt(cost_mu)<=2sqrt(integral ||x||^2 dmu)`, reversing mu and nu proves
continuity of F_t in W2. Coordinate marginals converge in W2 as well.

The same finite-type proof applies to a product of coordinatewise even
Gaussian kernels, using their positive semidefiniteness; its functional
uses the sum of the negative logarithms of those kernels. It must not
be silently relabeled as the shared-global-sign-folded vector kernel.

Finally, (3) is an upper bound on F, hence a lower bound on the associated
normalized permanent. Using it as an upper bound on a partition sum
requires a separate algebraic rearrangement in the specified recursion.
The inequality alone cannot supply that rearrangement.

## 7. Full-read audit of the specified finite-depth recursion

The subsequently written
`continued_convergence_recursive_orbit_bound_2026_09_06.md` supplies that
separate rearrangement. Its finite-depth upper recursion and soft
supersolution reconstruct, with the weave interface clarifications below.

At a node of size s=2q the exact matrix is

```math
U_s=\operatorname{diag}(U_q^{(1)},U_q^{(2)})
       2^{-1/2}\begin{pmatrix}I&I\\I&-I\end{pmatrix}g.
```

Both child matrices have entries of magnitude `1/sqrt(q)`, so the node
has entries of magnitude `1/sqrt(s)` and is orthogonal. Fresh uniform
signed permutation g precedes the pairing and is independent of both
children, precisely as required by Section 5.

For a symmetric empirical input law nu, an admissible signed pair law pi
prescribes the average ABSOLUTE marginal, not each signed marginal.
If N_a are the absolute multiplicities and n_ab the ordered pair counts,
the orbit has `s! 2^(s-N_0)/prod_a N_a!` distinct signed words, while the
pair table has `q!/prod_ab n_ab!` realizations. The logarithmic probability
divided by s is

```math
\tfrac12H(\pi)-H(\nu)
   =-\tfrac12D(\pi\Vert\nu\otimes\nu)+o(1).
```

The equality holds for unequal signed marginals because `log nu(a)`
depends only on |a|. This checks the entropy normalization in the Bellman
operator. Concatenation of output lists gives
`L(v_+,v_-)<=L(v_+)L(v_-)` by inclusion of the full-group invariant
space into the block-group invariant space. Thus the two child
contributions have weights 1/2, while the pair entropy cost also has
weight 1/2. Independence of the child bases gives the required product
of conditional expectations.

At any fixed depth r the alphabet and number of type coordinates are
finite; all type-count errors are polynomial in s. There are only `2^r`
terminal nodes, each of dimension `s/2^r` and with norm squared at most
the original root norm squared. Their total orbit-theorem loss is
`O_{r,t,C}(sqrt(s))`. Consequently the limit is indeed taken as
`s -> infinity` with r,t fixed. The claimed bound does not license a
growing-depth substitution without additional estimates.

For the supersolution, let `C_t(mu)=2H(mu)-F_t(mu)`. Mixing feasible
self-couplings proves concavity of C_t, while coupling entropy proves
its coordinate subadditivity. Concatenation proves concavity of raw
Gaussian F_t as well: here use the UNSIGNED permutation-group analogue
of the projection argument, so it applies also to nonsymmetric laws.
Both functionals are reflection invariant; symmetrizing a law can only
increase either one. Combined with (3), an admissible pair law pi and
its orthogonal output coordinates C,D therefore satisfy

```math
F_t(\nu_+)+F_t(\nu_-)\ge F_t(\pi),
\qquad C_t(\pi)\le2C_t(\nu).
```

The second inequality uses subadditivity at the two possibly unequal
signed marginals, reflection symmetrization, and then concavity at their
average, which is exactly nu. Since
`D(pi||nu tensor nu)=2H(nu)-H(pi)`, the difference between Phi_t(nu)
and the proposed Bellman summand is exactly

```math
\frac14\{F_t(\nu_+)+F_t(\nu_-)-F_t(\pi)
                  +2C_t(\nu)-C_t(\pi)\}\ge0.
```

This checks both the direction and every factor of two. The operator is
monotone, hence the successive Bellman upper bounds decrease. Selected
test couplings give LOWER bounds on these suprema, not certified upper
evaluations of them.

At t=infinity and nu uniform on {-1,+1}, put e=Pr(A=B). The output
entropy sum is `2h(e)+log 2`, while the input KL cost is at least
`log 2-h(e)`, with equality after balancing the two signs within each
category. The first Bellman gain is therefore exactly `(log 2)/4`.
This is an infinite-tilt calculation, not a finite-t cap certificate.

## 8. Necessary weave-interface choices

The existing weave definition has row spectrum `H[T,:]^T x`. To apply
the recursion to a supported input spin, choose the base matrix
`H=sqrt(m) U_m^T`. With input values `0,+/-1/sqrt(p_m)`, where
`p_m=k/m`, its normalized spectrum is exactly U_m applied to that input.
Using `H=sqrt(m) U_m` without transposing would reverse this interface.

Furthermore, choose an independent recursive base in EACH fibre, then
apply fresh independent output column permutations for the Finner proof.
That proof works for distinct deterministic fibre bases, so averaging
over these independent bases yields the product of their expected
one-row sums. Sharing a single random base across all fibres instead
would yield a moment of the one-row sum, not the m-th power of its mean.

With these choices the stated sufficient criterion
`p log 2 + B^r Phi_t(nu_p) + t(1-sqrt(p)) < 0` is correctly normalized.
No strict finite variational certificate has been supplied in the audited
file. The matching-row lower obstruction must still be respected, and
the square root remains inside the base expectation throughout.

## 9. A valid first-iterate symmetry reduction

There is an additional concavity useful for finite evaluation:

```math
R_t(\mu):=H(\mu)-F_t(\mu)
 =\sup_{\eta\in\Pi(\mu,\mu)}
       \{H_\eta(Y\mid X)-t\mathbb E_\eta(X-Y)^2\}
```

is concave in mu. Conditional entropy is jointly concave in the joint
law, and feasible self-couplings mix with their marginals. Therefore
for every mixture the Jensen gaps satisfy `Delta F_t<=Delta H`.

Fix the input nu and mix feasible pair laws pi. The symmetrized output
laws are fixed stochastic pushforwards of pi. The entropy mixture gap
of either output is at most that of pi, by data processing for the
mixture label. Consequently the first Bellman summand has Jensen gap

```math
\tfrac12\Delta H(\pi)
 -\tfrac14[\Delta F_t(\nu_+)+\Delta F_t(\nu_-)]\ge0.
```

It is thus a concave function of the pair law. Averaging over independent
input sign flips and input-coordinate swap preserves its value on each
group image and cannot reduce its value after averaging. An optimizer
may therefore be taken to have a symmetric coupling theta of the input
magnitudes, with conditionally independent uniform signs. The magnitude
marginals of theta are both the prescribed magnitude law.

For ternary input `(1-p)delta_0+(p/2)(delta_{-a}+delta_a)`, write
`theta=[[1-p-u,u],[u,p-u]]`, `0<=u<=min(p,1-p)`. Both symmetrized child
laws have magnitudes `0,a/sqrt(2),sqrt(2)a` with probabilities

```math
(1-p/2-3u/2, 2u, (p-u)/2).
```

The input pair KL is exactly `D(theta||alpha tensor alpha)`, where
`alpha=(1-p,p)`. This justifies the one-dimensional FIRST-iterate
reduction in the convergence agent's diagnostic script. A rigorous
upper evaluation still requires controlled functional evaluations or
an analytic bound; a floating optimizer alone is not that certificate.

No preservation of the concavity of `H+2f` under the Bellman operator
has been proved here. Thus this symmetry reduction must not be imported
automatically into deeper iterates with f different from Phi_t.
