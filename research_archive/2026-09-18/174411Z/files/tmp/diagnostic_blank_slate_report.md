# Diagnostic 2/3: blank-slate convergence attack

Finalized: 2026-07-30 (UTC)

Pre-ledger candidates: `tmp/diagnostic_blank_slate_scratch.md`

Finite checker: `tmp/diagnostic_blank_slate_checker.py`

## Verdict

None of the three candidates supplies a new convergence theorem.  Candidate 2
does yield a decisive falsifier of the proposed *mechanism*: in the quantifier
order needed for its convergence proof, far-down principal restriction is
locally universal and collapses exactly to the desired scalar inequality.
Candidates 1 and 3 require a summable power saving, not an unnamed little-oh;
the available block and interpolation identities have an extensive defect.

The useful new item in this diagnostic is the cut-quasirandomness/local-
universality proof below.  It explains why taking the large parent arbitrarily
far above the target order cannot make optimized restriction easier.

Throughout,

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
D(A)=\max_x|H_A(x)|,\qquad c_n=M_n/n^{3/2}.
\]

## 1. Candidate 2: far-down restriction is the scalar conclusion in disguise

### 1.1 Exact proposed quantifiers and convergence implication

The pre-ledger proposal, with a general positive exponent, is:

> There are constants `C,gamma>0` such that, for all large `n`, all `N>=n`,
> and every competitive order-`N` signing `A` (say
> `D(A)<=N^(3/2)`), some `n`-set `S` satisfies
> \[
> \frac{D(A[S])}{n^{3/2}}
> \le \frac{D(A)}{N^{3/2}}
> +C\left(n^{-\gamma}+\sqrt{n/N}\right).                 \tag{R_far}
> \]

Only exact minimizers are needed, but the stronger competitive-signing
quantifier was the candidate as stated.  Its convergence implication is
immediate and should be kept explicit.  Let

\[
\ell=\liminf_N c_N.
\]

Fix `n`, take `N_k/n -> infinity` with `c_(N_k) -> ell`, and apply
`(R_far)` to an optimal `A_(N_k)`.  Since every induced order-`n` signing has
norm at least `M_n`,

\[
c_n\le c_{N_k}+C\left(n^{-\gamma}+\sqrt{n/N_k}\right).
\]

Sending `k` to infinity gives

\[
c_n\le\ell+Cn^{-\gamma}.
\]

Then `limsup c_n <= ell = liminf c_n`, so the normalized minima converge.

### 1.2 Cut-quasirandomness lemma

The apparent structural content disappears in exactly that order of limits.

**Lemma.**  Let `A_N` be any sequence of signings with
`D(A_N)=O(N^(3/2))`.  For every fixed finite signing `P`, every sufficiently
large `A_N` contains a principal induced copy of `P`.

**Proof.**  For disjoint vertex sets `U,V`, write

\[
w_A(U,V)=\sum_{i\in U,j\in V}a_{ij},
\qquad
C_A(S)=w_A(S,S^c).
\]

If `x^S` is `-1` on `S` and `+1` off `S`, then

\[
H_A(\mathbf1)-H_A(x^S)=2C_A(S),
\]

so

\[
|C_A(S)|\le D(A).                                          \tag{1.1}
\]

For disjoint `U,V`, with `R` the remaining vertices,

\[
C_A(U)+C_A(V)-C_A(U\cup V)=2w_A(U,V),
\]

and hence

\[
|w_A(U,V)|\le\frac32D(A).                                 \tag{1.2}
\]

Also `|sum_(i<j in U) a_ij|<=D(A)`: restrict to `U`, use the
conditional-extension inequality `D(A[U])<=D(A)`, and evaluate its all-one
spin.  Decomposing arbitrary (possibly overlapping) rectangles into three
disjoint rectangles and one internal square now gives, for the ordered
zero-diagonal matrix,

\[
\max_{X,Y}\left|\sum_{i\in X,j\in Y}a_{ij}\right|
\le\frac{13}{2}D(A).                                       \tag{1.3}
\]

The constant is irrelevant.  Let `G_A` be the graph whose edge `ij` is
present when `a_ij=+1`.  Its step graphon therefore obeys

\[
\left\|W_{G_A}-\frac12\right\|_\square
\le \frac{13D(A)}{4N^2}+O(N^{-1})=O(N^{-1/2}).              \tag{1.4}
\]

(The `O(N^-1)` term only handles diagonal cells.)  Thus `G_A` converges in
cut norm to the constant-one-half graphon.

Let `P=(p_ij)` have fixed order `k`.  The indicator that an injective tuple
`phi:[k]->[N]` realizes `P` is

\[
\prod_{i<j}\frac{1+p_{ij}a_{\phi(i)\phi(j)}}2.
\]

The ordinary graphon counting lemma, applied after expanding the absent-edge
factors, and (1.4) show that the proportion of such injective tuples is

\[
2^{-\binom k2}+o(1).                                       \tag{1.5}
\]

Collisions cost only `O_k(1/N)`.  The limit in (1.5) is positive, proving
the claim.  QED.

### 1.3 Exact collapse of the candidate

Fix `n` and choose one optimal order-`n` signing `P_n`.  The lemma gives, for
every competitive sequence `A_N`,

\[
\boxed{
\min_{|S|=n}D(A_N[S])=M_n
}
\qquad\text{for all sufficiently large `N`.}               \tag{1.6}
\]

The lower bound is the definition of `M_n`; the induced copy of `P_n` gives
the upper bound.  Therefore, in the only quantifier order used in Section
1.1 (fix `n`, then send the liminf parent order to infinity), `(R_far)` says
precisely

\[
c_n\le\ell+Cn^{-\gamma}.
\]

It has not converted parent structure into scale transfer.  It has asked for
the scalar conclusion after quasirandomness has made the parent contain every
fixed local pattern.

This also explains why fixed graphons and fixed subgraph statistics are blind
to the problem: all competitive signings have the same local limit, even if
their `n^(3/2)` fluctuation norms differ.

### 1.4 Coupled scales return to the known adaptive-selector obligation

To avoid (1.6), the child order must grow with the parent.  At fixed density
`m/N -> alpha`, define in one-copy normalization

\[
V_{\rm ad}(A,m)
=\min_{|S|=m}\max_{\sigma=\pm1,x}
 \sigma H_{A[S]}(x_S)
=\min_{|S|=m}D(A[S]).                                      \tag{1.7}
\]

If the selector is hidden from the maximizing spin, uniform `S` retains each
edge with probability

\[
p_2=\frac{(m)_2}{(N)_2},
\]

so the hidden minimax value satisfies

\[
V_{\rm hid}(A,m)\le p_2D(A).                               \tag{1.8}
\]

Writing `I(A,m)=V_ad-V_hid`, a sufficient power-saving rounding statement is

\[
I(A,m)
\le\left[(m/N)^{3/2}-p_2\right]D(A)+O(N^{3/2-\delta}).      \tag{1.9}
\]

Equations (1.7)--(1.9) are exactly the hidden-versus-revealed selector game
recorded in ledger (10.643)--(10.646).  A fixed-density recurrence, uniform
over a compact ratio window,

\[
M_m\le(m/N)^{3/2}M_N+O(N^{3/2-\delta})                     \tag{1.10}
\]

can be iterated along geometrically decreasing orders; the normalized errors
sum to `O(n^-delta)` at terminal order `n`, and convergence follows.  Thus the
non-tautological form of Candidate 2 is not weaker than the existing
power-saving restriction target.  Its exact missing content is the adaptive
maximum in (1.9).

For `m/N -> 0`, the hidden term is only

\[
p_2D(A)=O(m^2/\sqrt N)=o(m^{3/2}),
\]

so essentially the whole desired child bound sits in the adaptivity gap.  A
hidden-selector average supplies no far-down shortcut.

## 2. Candidate 1: block composition needs a summable defect

### 2.1 What is unconditional

For optimal children `A,B` of orders `m,r`, form

\[
G_C=\begin{pmatrix}A&C\\C^T&B\end{pmatrix}.
\]

An independent sign bridge and a union bound give a deterministic `C` with

\[
M_{m+r}\le M_m+M_r+
\sqrt{2mr(m+r+2)\log2}.                                    \tag{2.1}
\]

This proves useful local continuity when `r=o(m)`, but the bridge in (2.1) is
of leading order for a balanced split.

### 2.2 Correct approximate-Fekete threshold

Put `b_n=M_n^(2/3)`.  A genuine sufficient theorem would be

\[
b_{m+r}\le b_m+b_r+e(m+r),                                 \tag{2.2}
\]

with a Hammersley-summable modulus, for example

\[
e(N)=O(N^{1-\delta})
\quad\text{or more generally}\quad
\sum_{N\ge1}\frac{e(N)}{N^2}<\infty.                       \tag{2.3}
\]

Then the almost-subadditive lemma forces `b_n/n`, hence `c_n`, to converge.
Uniform `e(N)=o(N)` is not enough.  The explicit countermodel

\[
b_n=n\left[1+\varepsilon\sin(\log\log(n+n_0))\right]       \tag{2.4}
\]

has no normalized limit but satisfies

\[
\sup_{m+r=N}\frac{|b_N-b_m-b_r|}{N}\longrightarrow0.       \tag{2.5}
\]

For (2.5), split at `min(m,r)<=eta N`.  That part contributes
`O(eta N)+o(N)`; when both parts exceed `eta N`, slow variation is uniform on
`[eta N,N]`.  Send `N` to infinity and then `eta` to zero.

### 2.3 Why the block map does not approach (2.3)

The exact relative-sign identity is

\[
D(G_C)=\max_{x,y}
\left(|H_A(x)+H_B(y)|+|x^TCy|\right).                       \tag{2.6}
\]

Thus cross edges cannot cancel internal energy.  With exact energy layers
`X_A(p),X_B(q)`, (2.6) is

\[
D(G_C)=\max_{p,q}
\left\{|p+q|+
\max_{x\in X_A(p),y\in X_B(q)}|x^TCy|\right\}.             \tag{2.7}
\]

If the target energy is

\[
T=(M_m^{2/3}+M_r^{2/3})^{3/2},
\]

then an `e(N)=O(N^(1-delta))` defect in (2.2) requires the right side of
(2.7) to be at most

\[
T+O(N^{3/2-\delta}).                                       \tag{2.8}
\]

The scalar child values do not control the restricted bilinear norms in
(2.7).  An iid bridge has a `Theta(N^(3/2))` bulk-layer floor, and (2.1) pays
the same leading scale.  The ledger leaves only a structured flat bridge with
uniform anti-alignment on every high child-energy layer; no candidate map or
estimate here supplies the power saving in (2.8).

The finite exact value `M_10=13` also makes the zero-defect `5+5` inequality
fail:

\[
13^{2/3}-2\,4^{2/3}=0.4890906\ldots>0.
\]

This is only finite evidence, not an asymptotic falsifier.

## 3. Candidate 3: finite-temperature interpolation has extensive curvature

### 3.1 Exact zero-temperature map

Define

\[
Z_A(\beta)=\mathbb E_x
 \exp\left(\frac{\beta|H_A(x)|}{\sqrt n}\right),
\qquad
P_n(\beta)=\frac1n\min_A\log Z_A(\beta).
\]

For every `n,beta`,

\[
c_n-\frac{\log2}{\beta}
\le\frac{P_n(\beta)}\beta\le c_n.                          \tag{3.1}
\]

Hence convergence of `P_n(beta)` at every fixed `beta`, followed by
`beta -> infinity`, would prove convergence of `c_n`.

### 3.2 Correct additivity threshold

Writing `a_n(beta)=nP_n(beta)`, a sufficient interpolation theorem is

\[
|a_{m+r}(\beta)-a_m(\beta)-a_r(\beta)|\le e_\beta(m+r),     \tag{3.2}
\]

where

\[
\sum_N e_\beta(N)/N^2<\infty;
\]

`e_beta(N)=O_beta(N^(1-delta))` is enough.  Again, a bare uniform `o(N)`
does not suffice, by the slowly oscillating model (2.4).

### 3.3 Exact annealed bridge bound and the scale wall

Let `N=m+r`, `t=beta/sqrt(N)`, and use a random sign bridge `C`.  The bounds

\[
e^{|z|}\le2\cosh z,
\qquad
\cosh(u+v)\le2\cosh u\cosh v
\]

give, after averaging `C` and then choosing one realization,

\[
\begin{aligned}
NP_N(\beta)\le{}&
mP_m\!\left(\beta\sqrt{m/N}\right)
+rP_r\!\left(\beta\sqrt{r/N}\right)\\
&+mr\log\cosh(\beta/\sqrt N)+\log4.                       \tag{3.3}
\end{aligned}
\]

For a balanced split, the explicit bridge term is

\[
mr\log\cosh(\beta/\sqrt N)
=\frac{\beta^2mr}{2N}+O_\beta(1)=\Theta_\beta(N).           \tag{3.4}
\]

Moreover, (3.3) evaluates each child at the smaller normalized temperature
`beta sqrt(m/N)`, not at `beta`.  Reaching (3.2) requires a structural
cancellation between two extensive terms.  Convexity or scalar pressure
bounds do not provide it.

Ordinary spin-glass covariance interpolation is also blind here.  If a fixed
signing is randomly switched by independent vertex signs `g_i`, then

\[
\mathbb E_g H_{A^g}(x)H_{A^g}(y)
=\sum_{i<j}x_ix_jy_iy_j
=\frac{(x\cdot y)^2-n}{2},                                 \tag{3.5}
\]

independently of `A`.  Yet Boolean ground norms range from order `n^(3/2)` to
order `n^2`.  An interpolation using only (3.5) cannot distinguish the
structures being minimized.  Higher cumulants/growing replicas or an
equivalent full energy-layer state would be needed to reduce the extensive
curvature in (3.4); this diagnostic supplies no such estimate.

## 4. Finite checker and disposition

The checker exhausts all switching classes through order seven and records:

- the exact hard minima and minimizing finite-temperature pressures;
- finite `2/3`-power defects;
- best principal-child normalized gaps among hard optimizers; and
- same-temperature pressure additivity defects.

These data confirm finite obstructions (including optimizer principal gaps and
sign-changing pressure defects), but they are deliberately treated only as
evidence.  They do not amplify to separated asymptotic subsequences.

Final disposition:

1. Candidate 1 is a scalar/block restatement until a structured bridge proves
   the power-saving layer estimate (2.8).
2. Candidate 2 is tautological at far-down scales by (1.6); at coupled scales
   it is exactly the existing adaptive-selector/power-saving recurrence.
3. Candidate 3 is a valid sufficient language, but the available interpolation
   has an extensive, non-summable defect and covariance universality erases the
   signing.

Accordingly this blank-slate attack finds no new route with leverage.  It
recommends against treating far-down sampling, bare `o(N)` composition, or
ordinary covariance interpolation as independent convergence mechanisms.
