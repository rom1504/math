# A renormalized spectral-spike state for dense quadratic landscapes

**Status:** scoped Level-3 theorem benchmark for the rare-event portfolio.
The landscape is deterministic after one predeclared generic sequence is
frozen, the interaction may be made dense, and a zero-density extremal state
has an asymptotically exact generic composition law.  It is **not** a uniform
theorem for adversarial spike directions, and it has no present implication
for the Boolean signing problem.

## 1. Model and the information that weak bulk limits lose

Let

```math
H_C(x)=x^{\mathsf T}Cx,\qquad x\in S^{N-1},
```

for a real symmetric `N x N` matrix `C`.  Its extremal response is

```math
\operatorname{OPT}(C)=\max_{x\in S^{N-1}}H_C(x)=\lambda_{\max}(C).
                                                                    \tag{1}
```

Suppose that `B_N` is deterministic, `\sup_N\|B_N\|<\infty`, its empirical
spectral measure

```math
\mu_{B_N}=N^{-1}\sum_{i=1}^N\delta_{\lambda_i(B_N)}
```

converges weakly to a compactly supported probability measure `\mu`, and

```math
\lambda_{\max}(B_N)\longrightarrow b:=\max\operatorname{supp}\mu. \tag{2}
```

The usual bulk state `\mu` assigns mass zero to every fixed-rank component.
To retain such rare components in a declared bulk-plus-perturbation
presentation, augment it by an **un-normalized positive spike multiset**
`\Theta=(\theta_1,\ldots,\theta_r) subset (0,infinity)^r`.  The bulk measure
uses mass `1/N` per eigenvalue while each presented finite-rank perturbation
is retained as a unit mark.  This is not, in general, an intrinsic
renormalization of the total matrix's spectrum.  The candidate presented
rare-event state is therefore

```math
\mathfrak S(B_N,P_N)=(\mu,\Theta),\qquad
P_N=U_N\operatorname{diag}(\Theta)U_N^{\mathsf T}.          \tag{3}
```

Here `r` is fixed and the columns of `U_N` are Haar orthonormal directions,
independent of `B_N`.  Once an almost-surely good sequence of directions is
chosen, all matrices in the theorem below are deterministic.  A generic
common orthogonal conjugation makes both the bulk and perturbation dense
almost surely without changing the landscape optimum.

Write, for `z>b`,

```math
G_\mu(z)=\int {1\over z-t}\,d\mu(t),\qquad
G_\mu(b+)=\lim_{z\downarrow b}G_\mu(z).                    \tag{4}
```

## 2. Exact finite-`N` carrier and asymptotic collapse

For one spike `P_N=\theta u_Nu_N^{\mathsf T}`, the matrix determinant lemma
gives the exact equation

```math
\det(zI-B_N-P_N)
=\det(zI-B_N)\left(1-\theta G_{N,u_N}(z)\right),            \tag{5}
```

where

```math
G_{N,u}(z)
=u^{\mathsf T}(zI-B_N)^{-1}u
=\int {1\over z-t}\,d\nu_{N,u}(t),
\quad
\nu_{N,u}=\sum_i|\langle u,v_i\rangle|^2\delta_{\lambda_i(B_N)}. \tag{6}
```

Thus the exact rooted carrier is the weighted spectral measure `\nu_{N,u}`,
not merely the empirical spectrum.  For a Haar direction, concentration of
quadratic forms gives, almost surely,

```math
\nu_{N,u_N}\Rightarrow\mu,
\qquad G_{N,u_N}\longrightarrow G_\mu
```

locally uniformly on `(b,\infty)`.  Generic orientation therefore
**synchronizes the rooted carrier to the unrooted bulk state**.  This is the
mechanism that makes (3) smaller than the full spherical energy landscape.

The following is the rank-one specialization of the finite-rank deformation
theorem of
[Benaych-Georges and Nadakuditi](https://arxiv.org/abs/0910.2120); in the
present form it also follows directly from (5)--(6).

### Theorem 1 (bulk-plus-spike response composition)

Under (2), fix `\theta>0` and take a Haar direction `u_N` independent of
`B_N`.  Almost surely,

```math
\lambda_{\max}(B_N+\theta u_Nu_N^{\mathsf T})
\longrightarrow \mathcal R(\mu,\theta),                    \tag{7}
```

where

```math
\mathcal R(\mu,\theta)=
\begin{cases}
\rho, & \theta G_\mu(b+)>1,
          \quad G_\mu(\rho)=\theta^{-1},\ \rho>b,\\
b, & \theta G_\mu(b+)\le1.
\end{cases}                                                \tag{8}
```

The root in the first line is unique.  For every fixed positive spike
multiset carried by one jointly Haar-generic orthonormal frame, the
corresponding finite-rank theorem gives a finite list of possible upper
outliers determined by the same Cauchy-transform rule (with multiplicities).
Adding further independent Haar-generic finite-rank spikes therefore updates
the presented state asymptotically by multiset union,

```math
(\mu,\Theta)\star\Theta'=(\mu,\Theta\uplus\Theta'),         \tag{9}
```

and all upper outliers from these fixed positive spikes are recovered
asymptotically from this size-independent state.  The bulk component remains
`\mu` because finite-rank edits do not change its weak limit.  The measure
`\mu` need not itself be finite-dimensional unless the bulk class is fixed
or finitely parametrized.

#### Proof in the rank-one case

For `z>\lambda_{\max}(B_N)`, `G_{N,u_N}(z)` is positive and strictly
decreasing.  The only perturbed eigenvalue above the bulk edge is the unique
solution of `\theta G_{N,u_N}(z)=1`, when such a solution remains separated
from the edge.  Local uniform convergence in (6) sends that solution to the
unique `\rho` in (8).  If `\theta G_\mu(b+)\le1`, then for every
`\varepsilon>0`,

```math
\theta G_\mu(b+\varepsilon)<1.
```

Uniform convergence excludes a perturbed eigenvalue above
`b+\varepsilon` for all large `N`; interlacing and (2) give the matching
lower bound.  This proves (7).  The rank-`r` statement is precisely the
finite-rank determinant/secular-matrix extension proved in the cited paper.

The theorem is probabilistic only as an existence device.  For a countable,
predeclared family of deterministic bulk sequences and generic direction
frames, intersecting the probability-one resolvent events and freezing one
realization yields deterministic dense matrix sequences obeying (7)--(9)
simultaneously for all positive strengths.  This does not make (9) uniform
over directions selected adversarially after the realization is frozen.

## 3. A sharp fixed-depth insufficiency theorem

The rare component is not merely inconvenient for normalized moments: it is
topologically absent from every statistic continuous in the weak empirical
spectral law.

### Theorem 2 (fixed-depth bulk-indistinguishability, macroscopic response separation)

Assume `\theta G_\mu(b+)>1`, and let `\rho>b` solve
`G_\mu(\rho)=1/\theta`.  Put

```math
C_N^{(0)}=B_N,
\qquad C_N^{(1)}=B_N+\theta u_Nu_N^{\mathsf T}.             \tag{10}
```

Then:

1. the empirical spectral laws of `C_N^{(0)}` and `C_N^{(1)}` have the same
   weak limit `\mu`;
2. for every fixed `L`,

   ```math
   {1\over N}\operatorname{tr}(C_N^{(1)})^k
   -{1\over N}\operatorname{tr}(C_N^{(0)})^k\longrightarrow0,
   \qquad 1\le k\le L;                                    \tag{11}
   ```

3. nevertheless their extremal responses remain separated:

   ```math
   \operatorname{OPT}(C_N^{(1)})-\operatorname{OPT}(C_N^{(0)})
   \longrightarrow\rho-b>0.                               \tag{12}
   ```

Consequently no fixed-depth normalized trace/closed-walk summary, and more
generally no summary whose limiting value is continuous only in the weak
empirical spectral law of the already composed matrix, is uniformly
sufficient for the declared dense-spike continuation queries.

#### Proof

A rank-one perturbation changes the empirical distribution function by at
most `1/N` (rank inequality).  Uniform operator-norm bounds then imply (11)
for every fixed polynomial.  Equation (12) is Theorem 1 together with (2).

This conclusion has deliberately narrow quantifiers.  It does **not** say
that every bounded-dimensional statistic fails: the spike strength itself is
a bounded-dimensional statistic and is exactly the needed renormalized mark.
It says that the entire natural hierarchy of fixed-depth normalized bulk
observables loses a component whose mass is `1/N` but whose optimized response
is order one.

## 4. Why this is an actual rare-event composition theorem

The two scales are explicit:

| component | ordinary normalization | retained normalization | role |
|---|---:|---:|---|
| bulk eigenvalue cloud | mass `1/N` | probability measure `\mu` | supplies the response medium through `G_\mu` |
| each finite-rank spike | mass `1/N\to0` | unit mark `\theta_j` | creates a possible order-one outlier |

The response is neither a function of the bulk nor of the spike alone.  It is
the nonlinear composition `G_\mu(\rho)=1/\theta`.  At finite `N`, the rooted
weighted measure (6) is needed.  Generic orientation supplies a deterministic
synchronization theorem in the limit, collapsing that rooted information to
`\mu`.  Thus this is more than the tautology “record the maximizer”: it gives
a strict presented-state reduction and an exact checkable asymptotic response
law under generic contexts.  It is not asserted minimal: for a top-eigenvalue
only query, `(mu,current top)` is a smaller state than the full outlier
multiset.

If one also composes independent Haar-conjugated bulk matrices and the
deterministic bulk tuple has a strong limiting distribution (including the
required no-outlier/edge control), strong asymptotic freeness identifies the
limiting bulk by free additive convolution and controls the spectral edge;
see
[Collins and Male](https://arxiv.org/abs/1105.4345).  That extension is not
needed for Theorems 1--2 and is not claimed here for arbitrary deterministic
relative orientations.

## 5. Novelty audit against the Gaussian tangent state

This benchmark does not repackage the repository's Gaussian
Morse/tangent-mass composition.

1. **Different rare object.**  The Gaussian theorem retains saddle amplitude,
   Hessian and covariance data of a positive array.  Here the retained object
   is a zero-density spectral atom whose ordinary empirical mass vanishes.
2. **Different algebra.**  Gaussian tangent data compose by local convolution
   near a saddle.  Spectral spikes compose through a resolvent/secular
   equation and multiset union.
3. **Different synchronization.**  The Gaussian collapse comes from a local
   limit theorem.  Here Haar isotropy makes a rooted spectral measure converge
   to the unrooted empirical law.
4. **Different falsifier.**  The obstruction is finite-rank instability of
   weak spectral topology: all fixed normalized moments agree while an
   outlying maximum differs by a constant.

The imported random-matrix outlier theorem is classical.  The contribution
of this note is not a new random-matrix theorem; it is the exact placement of
that theorem in the extremal-information hierarchy and the rigorous
incompressibility comparison (Theorem 2).

## 6. Boundary of the result

This model is a Level-3 positive benchmark in the following scoped sense:

- matrices and landscapes are dense and can be frozen deterministically;
- for each fixed regular bulk sequence, a countable family of generic
  direction frames is declared before freezing, while spike strengths may
  range over positive values;
- the relative spike directions are generic, not adversarial;
- the number of spikes is fixed as `N\to\infty`;
- the configuration space is the sphere, not the Boolean cube.

If the direction is chosen adversarially or correlated with a prior spike,
the rooted measure `\nu_{N,u}` need not approach `\mu`; it may concentrate on
a single eigenspace, and the compact state (3) fails.  For example, two unit
spikes in the same direction have top eigenvalue two over zero bulk, while
two orthogonal generic unit spikes have top eigenvalue one; the formal
multiset `{1,1}` alone cannot distinguish them.  This is a useful exact
boundary: **presented rare marks compose only after a synchronization
mechanism identifies how they couple to the bulk.**  The analogous missing
theorem in harder deterministic landscapes would have to replace Haar
isotropy by a structural synchronization law.

No canonical files should be changed on the strength of this draft alone.
