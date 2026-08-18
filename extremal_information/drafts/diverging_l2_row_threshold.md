# Diverging row-L2 complexity: the bulk regime, an STP obstruction, and a sharp-edge counterexample

**Status.** Task-local threshold report.  This does not claim a pressure
theorem for diverging row-density complexity.  It proves that every
deterministic ingredient except the sharp upper edge remains subcritical
when `K_r=o(r)`, proves Marchenko--Pastur bulk convergence in that regime,
shows that the existing Strong Tail Projection argument can fail for
arbitrarily slowly diverging `K_r`, and gives a separate Hadamard-cluster
construction where the sharp edge itself fails at a larger subexponential
`K_r` scale.

## 1. Uniform deterministic reduction through `K_r=o(r)`

Let `mu_r` be an exact-sign row law with density `g_r=dmu_r/dU_r`, and put

```math
K_r=\mathbb E_{U_r}g_r^2.
\tag{DT.1}
```

No central symmetry is assumed.  Let

```math
m_r=\mathbb E R_r,
\qquad
\Sigma_r=\mathbb E R_rR_r^T.
\tag{DT.2}
```

Parseval gives

```math
\|m_r\|_2^2+{1\over2}\|\Sigma_r-I\|_F^2\le K_r-1.
\tag{DT.3}
```

Assume for this section that `K_r=o(r)` and choose

```math
\delta_r=(K_r/r)^{1/4}.
\tag{DT.4}
```

Since `K_r>=1`, this tends to zero.  Let `P_r^0` be the spectral projection
of `Sigma_r` outside `[1-delta_r,1+delta_r]`, join its range with
`span(m_r)`, call the joined projection `P_r`, and set `V_r=I-P_r`.  Then

```math
k_r:=\operatorname{rank}P_r
\le {2(K_r-1)\over\delta_r^2}+1
=O(\sqrt{rK_r})=o(r).
\tag{DT.5}
```

As in the bounded-`L2` theorem, noncommutation of the joined projection is
harmless: `ran(V_r)` is a subspace of the regular spectral space, so the
compression

```math
S_r=V_r\Sigma_rV_r|_{\operatorname{ran}V_r}
\tag{DT.6}
```

has spectrum in `[1-delta_r,1+delta_r]`, and `V_rm_r=0`.

### Proposition DT.1 (all non-edge costs remain subcritical)

Let `B_r` have independent rows with law `mu_r`, and let `W_r` be an iid
Rademacher bridge.  There is a rowwise coupling such that

```math
\begin{aligned}
\mathbb E\|(B_r-W_r)V_r\|_F
&=O\big(r^{3/4}(\log K_r)^{1/4}\big)=o(r),\\
\mathbb E\|B_rP_r\|_*&=o(r^{3/2}),\\
\mathbb E\|W_rP_r\|_*&=o(r^{3/2}).
\end{aligned}
\tag{DT.7}
```

**Proof.**  Jensen under `mu_r` gives

```math
D(\mu_r\|U_r)\le\log K_r.
\tag{DT.8}
```

The cube `T1` inequality therefore couples one row to a uniform row with
expected Hamming cost at most `sqrt(r log(K_r)/2)`.  Independent row copies
and the exact-sign mismatch factor four give the first line of (DT.7).

For the conditioned bridge, use only the trace identity
`tr(Sigma_r)=r`, not the growing operator bound:

```math
\begin{aligned}
\mathbb E\|B_rP_r\|_*
&\le\sqrt{k_r}
  \big(r\operatorname{tr}(P_r\Sigma_r)\big)^{1/2}\\
&\le r\sqrt{k_r}=o(r^{3/2}).
\end{aligned}
\tag{DT.9}
```

For the iid bridge,

```math
\mathbb E\|W_rP_r\|_*\le k_r\sqrt r=o(r^{3/2}).
\tag{DT.10}
```

This proves the proposition. `square`

The existence of a vanishing spectral window with uniformly `o(r)` peel
rank follows from (DT.3) exactly when `K_r=o(r)`, at the level of this
worst-case bound.  Transport alone would allow the larger condition
`log K_r=o(r)`; the Fourier rank is the first deterministic bottleneck.

## 2. The Marchenko--Pastur bulk still survives

Define the centered isotropic projected row

```math
X_r=S_r^{-1/2}V_rR_r\in\mathbb R^{d_r},
\qquad d_r=r-k_r\sim r.
\tag{DT.11}
```

For every cube event `A`, Cauchy--Schwarz gives

```math
\mu_r(A)\le\sqrt{K_r}\,U_r(A)^{1/2}.
\tag{DT.12}
```

If `A_r` is a positive semidefinite `d_r by d_r` matrix with uniformly
bounded operator norm, the usual Hanson--Wright calculation for the lifted
quadratic form gives, for every fixed `epsilon>0`,

```math
\Pr\{|X_r^TA_rX_r-\operatorname{tr}A_r|>\epsilon d_r\}
\le C\sqrt{K_r}\,e^{-c_\epsilon d_r}+o(1).
\tag{DT.13}
```

The `o(1)` records the deterministic centering displacement
`O(delta_r d_r)=o(d_r)`.  Since `K_r=o(r)` in particular has
`log K_r=o(r)`, the right side tends to zero.  Yaskov's weak quadratic-form
criterion applies.

### Proposition DT.2 (bulk and lower edge)

If `K_r=o(r)`, the empirical spectral distribution of

```math
r^{-1}\mathbb X_r^T\mathbb X_r
\tag{DT.14}
```

for `r` independent copies of `X_r` converges to the parameter-one
Marchenko--Pastur law.  Consequently

```math
\liminf_{r\to\infty}{\|B_rV_r\|_{op}\over\sqrt r}\ge2
\quad\hbox{in probability}.
\tag{DT.15}
```

**Proof.**  Equation (DT.13) is the required quadratic-form condition, and
`d_r/r->1`.  Marchenko--Pastur convergence puts positive empirical mass in
every fixed interval immediately below four.  Hence the largest sample
eigenvalue is at least `4-o_Pr(1)`.  Multiplication by `S_r^(1/2)`, whose
smallest singular value tends to one, gives (DT.15). `square`

This is a lower edge only.  It supplies no upper bound and therefore cannot
be inserted into the projected pressure theorem.

## 3. STP can fail for arbitrarily slow divergence

The square-root likelihood loss in (DT.12) is harmless for macroscopic
quadratic deviations, but not uniformly for low-rank projection tails.  The
following construction makes the obstruction explicit.

Let

```math
T_r=r^{-1/2}\sum_{i=1}^rW_i,
\qquad
a_r=\sqrt{\log\log r}.
\tag{DT.16}
```

Along even orders, let

```math
A_r=\{a_r\le|T_r|\le a_r+1\},
\qquad
L_r=\{|T_r|\le1/2\}.
\tag{DT.17}
```

Write `nu_A`, `nu_L` for the corresponding uniform conditional laws and
let `v_A=E_(nu_A)T_r^2`, `v_L=E_(nu_L)T_r^2`.  Put

```math
p_r={1-v_L\over v_A-v_L},
\qquad
\mu_r=p_r\nu_A+(1-p_r)\nu_L.
\tag{DT.18}
```

For large `r`, `0<p_r<1`, `p_r=Theta(a_r^(-2))`, and by construction
`E_mu T_r^2=1`.  The law is centrally symmetric and permutation invariant.
Therefore its mean is zero and

```math
\Sigma_r=I.
\tag{DT.19}
```

Indeed, permutation invariance makes every off-diagonal second moment equal,
and `E(sum R_i)^2=r` forces that common value to vanish.

Stirling's formula in the moderate-deviation range
`a_r=o(r^(1/6))` gives

```math
U_r(A_r)=\exp\{-a_r^2/2+O(a_r+\log a_r)\},
\qquad U_r(L_r)=\Theta(1).
\tag{DT.20}
```

The two mixture supports are disjoint, so

```math
K_r={p_r^2\over U_r(A_r)}
    +{(1-p_r)^2\over U_r(L_r)}
=(\log r)^{1/2+o(1)}(\log\log r)^{-2}+O(1).
\tag{DT.21}
```

Thus `K_r->infinity` but `K_r=r^{o(1)}=o(r)`.

Let `Q_r` be projection onto the all-ones unit vector.  On `A_r`,
`||Q_rR||_2^2=T_r^2>=a_r^2`, with probability `p_r`.  Were the array to
satisfy the Chafaï--Tikhomirov STP inequality, its rank-one instance at
`t_r=a_r^2-1` would require

```math
p_r\le {g(1)\over t_r^2}.
\tag{DT.22}
```

But `p_rt_r^2=Theta(a_r^2)->infinity`, so no finite dimension-independent
value `g(1)` exists.  Hence:

### Proposition DT.3 (the current upper-edge proof has no diverging-`K` regime)

There are centered, isotropic exact-sign row laws with `K_r=r^(o(1))` for
which Strong Tail Projection fails.  Therefore the bounded-`K` STP proof
does not extend uniformly to **any** condition merely asserting
`K_r->infinity` slowly.

This proposition does not show that the sample edge for (DT.18) exceeds two.
It is a falsifier of the imported sufficient hypothesis, not of the desired
sharp edge.

## 4. A genuine sharp-edge counterexample from Hadamard latent clusters

The sharp edge itself eventually becomes false.  The construction is
robust against every deterministic sublinear-rank column peel.

Take Hadamard orders `r=2^m`, let `h_1,...,h_r in {+-1}^r` be the rows of a
Walsh--Hadamard matrix, and put `u_j=h_j/sqrt(r)`.  Thus `(u_j)` is an
orthonormal basis.  Fix a constant `L>2` and set

```math
c_r=L\sqrt{{\log\log r\over\log r}}.
\tag{DT.23}
```

To sample a row `R`, choose `(J,S)` uniformly from
`[r] times {+-1}`, and conditional on `(J,S)`, choose coordinates
independently with

```math
\Pr\{R_\ell=x\mid J=j,S=s\}
={1+s c_r h_{j\ell}x\over2}.
\tag{DT.24}
```

Call the resulting row law `mu_r`.

### Lemma DT.4 (exact covariance and density complexity)

The law `mu_r` is centrally symmetric, has mean zero, and satisfies

```math
\Sigma_r=I.
\tag{DT.25}
```

Its squared `L2(U_r)` density is exactly

```math
\boxed{
K_r=
{(1+c_r^2)^r+(1-c_r^2)^r\over2r}
+{r-1\over r}(1-c_r^4)^{r/2}.}
\tag{DT.26}
```

In particular,

```math
\log K_r
=(L^2+o(1)){r\log\log r\over\log r}.
\tag{DT.27}
```

**Proof.**  Averaging over `S` gives zero mean and central symmetry.
Conditional off-diagonal second moments are
`c_r^2 h_(j,a)h_(j,b)`.  Averaging over the Hadamard basis makes them zero;
the diagonal is one.

The component density relative to `U_r` is

```math
g_{j,s}(x)=\prod_{\ell=1}^r(1+s c_rh_{j\ell}x_\ell).
\tag{DT.28}
```

For two components,

```math
\langle g_{j,s},g_{k,t}\rangle_U
=\prod_{\ell=1}^r
 (1+st c_r^2h_{j\ell}h_{k\ell}).
\tag{DT.29}
```

If `j=k`, this is `(1+st c_r^2)^r`.  If `j!=k`, orthogonality gives equally
many agreeing and disagreeing coordinates, so it is
`(1-c_r^4)^(r/2)`.  Sum all ordered pairs in
`g=(2r)^(-1)sum_(j,s)g_(j,s)` to obtain (DT.26).  Expanding logarithms in
(DT.26) proves (DT.27). `square`

### Lemma DT.5 (restricted occupancy)

Throw `r` independent balls uniformly into the `2r` signed bins
`(j,s)`.  For every deterministic set of bins `G_r` with
`|G_r|=2r-o(r)`,

```math
\max_{(j,s)\in G_r}N_{j,s}
=(1+o_\Pr(1)){\log r\over\log\log r}.
\tag{DT.30}
```

The statement is uniform over the choice of the deterministic set.

**Proof sketch.**  The upper bound is the standard union bound using
`Pr(Bin(r,1/(2r))>=k)`.  For the lower bound, Poissonize the number of balls.
The counts in the retained bins become independent Poisson variables with
mean `1/2`.  At

```math
k=(1-\epsilon){\log r\over\log\log r},
```

Stirling gives `Pr(Pois(1/2)>=k)=r^{-(1-epsilon)+o(1)}`.  Since there are
`2r-o(r)` retained bins, the expected number reaching `k` is
`r^(epsilon+o(1))`, and independence makes the probability of none tend to
zero.  Standard monotone de-Poissonization returns to exactly `r` balls.
All estimates depend on the retained set only through its cardinality.
`square`

### Theorem DT.6 (every deterministic `o(r)` peel leaves an excess edge)

Let `B_r` have `r` independent rows with law (DT.24).  For every
deterministic orthogonal projection `P_r` with rank `d_r=o(r)`, put
`Q_r=I-P_r`.  Then

```math
\boxed{
\liminf_{r\to\infty}
{\|B_rQ_r\|_{op}\over\sqrt r}\ge L
\quad\hbox{in probability}.}
\tag{DT.31}
```

Consequently, for every fixed `c_0<L-2`,

```math
\Pr\{\|B_rQ_r\|_{op}>(2+c_0)\sqrt r\}\longrightarrow1.
\tag{DT.32}
```

**Proof.**  Write

```math
a_j^2=\|Q_ru_j\|_2^2=1-\|P_ru_j\|_2^2.
\tag{DT.33}
```

Since `(u_j)` is an orthonormal basis,

```math
\sum_{j=1}^r\|P_ru_j\|_2^2=\operatorname{tr}P_r=d_r.
\tag{DT.34}
```

Choose any deterministic `eta_r downarrow0` with `d_r/eta_r=o(r)`, for
example `eta_r=sqrt(d_r/r)` when `d_r>0`.  All but `o(r)` indices satisfy
`a_j^2>=1-eta_r`.  Retain the two signed bins over these good indices.
Lemma DT.5 supplies a bin `(j,s)` with

```math
N_{j,s}=(1-o_\Pr(1)){\log r\over\log\log r}.
\tag{DT.35}
```

For a good `j`, define the unit vector

```math
w_j={Q_ru_j\over a_j}.
\tag{DT.36}
```

Conditional on a row belonging to latent bin `(j,s)`, independence in
(DT.24) gives

```math
\mathbb E\langle R,w_j\rangle
=s c_r\sqrt r\,a_j,
\qquad
\operatorname{Var}(\langle R,w_j\rangle)=1-c_r^2.
\tag{DT.37}
```

Moreover the centered variable is uniformly subgaussian by Hoeffding,
because `||w_j||_2=1`.  Since

```math
c_r^2r\asymp {r\log\log r\over\log r}\gg\log r,
\tag{DT.38}
```

a union bound over all `r` sampled rows shows that, simultaneously for every
row in its own good latent bin,

```math
s\langle R,w_j\rangle
=(1-o(1))c_r\sqrt r\,a_j.
\tag{DT.39}
```

Because `Q_rw_j=w_j`, restrict `B_rQ_rw_j=B_rw_j` to the rows in the
maximally occupied bin.  Equations (DT.35) and (DT.39) give

```math
{\|B_rQ_r\|_{op}\over\sqrt r}
\ge(1-o_\Pr(1))c_ra_j\sqrt{N_{j,s}}
=(1-o_\Pr(1))L.
\tag{DT.40}
```

This proves the theorem. `square`

The projection is allowed to know the row law and the Hadamard basis, but it
must be deterministic rather than chosen after seeing the sampled latent
occupancies.  This is exactly the setting of a covariance/mean peel.  Since
`Sigma_r=I` and `m_r=0`, the canonical peel is in fact zero.

## 5. Exact threshold picture currently proved

The present evidence separates three statements.

1. If `K_r=o(r)`, one can always choose a vanishing spectral window with
   `o(r)` exceptional rank; Hamming transport and both nuclear costs are
   subcritical, and the projected empirical spectral **bulk** is
   Marchenko--Pastur.

2. Even for `K_r=r^(o(1))`, Strong Tail Projection can fail.  Thus the
   bounded-`K` sharp-edge proof supplies no theorem in a general diverging
   regime.  A different upper-edge argument would be required.

3. At

   ```math
   \log K_r=\Theta(r\log\log r/\log r),
   ```

   the sharp edge itself can fail by a fixed amount, after every
   deterministic `o(r)`-rank peel.

There remains a genuine gap between the slowly diverging STP obstruction
and the Hadamard excess-edge construction.  In particular this report does
**not** decide whether the sharp edge remains two uniformly under the sole
condition `K_r=o(r)`.  No pressure no-gain theorem is claimed in that open
regime.
