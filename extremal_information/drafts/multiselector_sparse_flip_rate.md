# Packing diffuse PC.3 sparse-flip channels

**Status.**  Phase-I theorem and scope audit.  There is a polynomial-in-order
packing of sparse-flipped PC.3 children which is pairwise separated by one
common bank of selector-indexed, balanced-compiler queries.  The exponent
proved here is deliberately crude.  The result is a one-hot message packing,
not a coordinate-decodable Boolean edit cube and not one physical parent with
all query shores simultaneously free.  Two elementary ceilings are recorded
for those stronger models.

Only Theorems 21.62--21.67 and their cited PC.3 drafts/verifier are used.

## 1. The metric comes before the construction

Let `A(u)` be hollow exact-sign children of the common PC.3 order
`N_j=16^j`.  A labelled query `ell` is an integer field `g(ell)`, and its
child response is

```math
 R_j(u\mid\ell)=\mathcal B_{A(u)}(g(\ell)).       \tag{MS.1}
```

After balanced compilation, let `B(ell)` be an `N_j times s_j` exact-sign
cross block, let `C_j` be one common hollow exact-sign shore, and put

```math
 P_j(u\mid\ell)=
 \begin{pmatrix}A(u)&B(\ell)\\B(\ell)^T&C_j\end{pmatrix},
 \qquad
 \widetilde R_j(u\mid\ell)=Q(P_j(u\mid\ell)).    \tag{MS.2}
```

For a fixed query bank `U_j`, the pairwise contextual metric is

```math
 d_j(u,v)={1\over (N_j+s_j)^{3/2}}
 \max_{\ell\in U_j}
 |\widetilde R_j(u\mid\ell)-\widetilde R_j(v\mid\ell)|.
                                                               \tag{MS.3}
```

A `delta`-packing means `d_j(u,v)>=delta+o(1)` for every distinct pair,
with the `o(1)` uniform over the family.  A stronger oriented one-hot code
has

```math
 \widetilde R_j(v\mid u)-\widetilde R_j(u\mid u)
 \ge(\delta+o(1))N_j^{3/2}
 \qquad(u\ne v).                                 \tag{MS.4}
```

This definition separates three statements which otherwise look similar.

1. A list of one-bit constructions only asserts, separately for every `u`,
   that `(A_0,A(u))` is separated by query `u`.  It says nothing about
   `(A(u),A(v))`.
2. Equations (MS.2)--(MS.4) use one common family of child states and one
   common query bank.  The context `ell` is externally selected; all child
   states are tested against the same compiled block `B(ell)`.
3. A Boolean edit cube would instead have states `A_b`,
   `b in {0,1}^{k_j}`, and require coordinate-wise decoding whenever
   `b_ell ne b'_ell`.  The theorem below does not construct that cube.

Likewise, (MS.2) is not a claim that the blocks `B(ell)` can all be placed
as free shores in one matrix.  Free shores add their fields rather than
externally selecting one of them; Corollary 21.65 is already a warning that
this distinction changes the leading scale.

## 2. Counterexample first: diffuse one-bit channels need not pack

Write the projectivized PC.3 endpoint as

```math
 \epsilon_0=1,
 \qquad (\alpha_t,\beta_t)
   =(\epsilon_{2t+1},\epsilon_{2t+2})\in\{\mathord\pm1\}^2.
                                                               \tag{MS.5}
```

After cancelling the common base pole, its row field and selector are

```math
 L_{\alpha,\beta}=1+\sum_{t<j}(\alpha_tX_t+\beta_tY_t),
 \qquad f_{\alpha,\beta}=\operatorname {sgn}L_{\alpha,\beta}, \tag{MS.6}
```

under the independent three-atom PC.3 law

```math
 (X,Y)=(1,1),(1,-1),(-1,1)
 \quad\hbox{with probabilities}\quad {1\over4},{1\over2},{1\over4}.
                                                               \tag{MS.7}
```

Take two endpoints with the same `beta`, and obtain `alpha'` from a balanced
`alpha` by swapping one `+1` and one `-1`.  Both have
`|sum_t alpha_t|=0`, but their fields differ pointwise by at most four.
Uniform anti-concentration for the remaining `Theta(j)` independent
increments gives

```math
 \Pr\{f_{\alpha,\beta}\ne f_{\alpha',\beta}\}=O(j^{-1/2}),
 \qquad
 {x_{\alpha,\beta}^Tx_{\alpha',\beta}\over N_j}
       =1-O(j^{-1/2}).                           \tag{MS.8}
```

Both selectors are diffuse in the sense proved in Lemma MS.1 below.  Yet
their two sparse-flip means satisfy

```math
 {H^{(\alpha,\beta)}-H^{(\alpha',\beta)}\over\sqrt {N_j}}
 =-\kappa(e e^T-e'e'^T)+o_{op}(1)=o_{op}(1),    \tag{MS.9}
```

because the norm of the difference of the two rank-one projections is
`sqrt(1-(e^Te')^2)=O(j^{-1/4})`.  Hence every labelled response with a
common field differs by `o(N_j^{3/2})`: for Boolean `y`, the quadratic
difference is at most `N_j\|H-H'\|_{op}/2`.  The same argument applies after
any one common balanced compiler, because the two parent blocks then differ
only in this child block.  Thus even a large collection of individually valid,
diffuse base-versus-flip channels can collapse to one point in (MS.3).
Pairwise selector geometry is indispensable.

## 3. A uniformly diffuse endpoint class

### Lemma MS.1 (balanced first coordinates are uniformly diffuse)

There is an absolute `C` such that every endpoint in (MS.5) satisfying

```math
 \left|\sum_{t<j}\alpha_t\right|\le1            \tag{MS.10}
```

obeys

```math
 \max_{z\in Z_j}{|x_{\alpha,\beta}^Tz|\over N_j}
 \le {C\over\sqrt j},                           \tag{MS.11}
```

uniformly over the arbitrary signs `beta`.

#### Proof

For a local score `A_t=alpha_tX_t+beta_tY_t`,

```math
 \mathbb EA_t={\alpha_t\over2},
 \qquad
 \operatorname {Var}(A_t)={7\over4}-\alpha_t\beta_t
       \in\left\{{3\over4},{11\over4}\right\}. \tag{MS.12}
```

Every active product, after cancelling the base pole, is a tensor product
of local characters in `{1,X,Y,XY}`, whose means are respectively
`1,1/2,0,-1/2`.  The proof of Theorem 21.62/DS.4 now applies verbatim and
uniformly: after deleting any `k` local factors, the remaining score has
variance comparable to `j-k`, while its mean has modulus at most
`(k+1)/2`.  Berry--Esseen supplies uniform `O(m^{-1/2})`
anti-concentration for `m` remaining factors.  Eliminating nonconstant
characters produces coefficients at most `2^{-k}`.  More precisely, after
deleting `k` local factors the remaining score mean has modulus at most
`(k+3)/2`; this harmless additive constant replaces the `(k+1)/2` bound in
the periodic-endpoint proof.  Splitting the resulting geometric error
series at `k=j/2` gives (MS.11).  The signs `beta_t` affect
neither the mean bound nor the uniform variance bounds. `square`

The same moment calculation gives, uniformly over this class,

```math
 V_{\alpha,\beta}:=\sum_{t<j}
       \left({7\over4}-\alpha_t\beta_t\right),
 \qquad {3j\over4}\le V_{\alpha,\beta}\le {11j\over4},       \tag{MS.13}
```

and

```math
 {\|h_{\alpha,\beta}\|_1\over N_j}
   =\sqrt{2\over\pi}\sqrt{V_{\alpha,\beta}}+O(1),
 \qquad
 {\|h_{\alpha,\beta}\|_2^2\over N_j}
   =V_{\alpha,\beta}+O(1).                     \tag{MS.14}
```

Indeed Wasserstein Berry--Esseen applies uniformly to the absolute-value
function, and the field mean is bounded under (MS.10).  In particular,

```math
 {x_{\alpha,\beta}^Th_{\alpha,\beta}
  \over\sqrt{N_j}\|h_{\alpha,\beta}\|_2}
 =\sqrt{2\over\pi}+O(j^{-1/2}).                 \tag{MS.15}
```

Thus all these endpoints have the same limiting response geometry as the
periodic endpoint in Theorem 21.64, even though their field variances need
not be equal.

## 4. Exponentially many endpoints with controlled mutual overlap

For two balanced endpoints write `C((alpha,beta),(alpha',beta'))` for the
covariance of their centred fields.  From (MS.7),

```math
 C=\sum_{t<j}\left{
 {3\over4}\alpha_t\alpha'_t+\beta_t\beta'_t
 -{1\over2}(\alpha_t\beta'_t+\beta_t\alpha'_t)
 \right\}.                                      \tag{MS.16}
```

### Lemma MS.2 (positive-rate endpoint code)

Fix `0<t<1/2`.  For all sufficiently large `j`, there is a set `U_j` of
endpoints satisfying (MS.10), of size

```math
 |U_j|\ge \exp(t^2j/100),                       \tag{MS.17}
```

such that every distinct pair has `|C|<=tj`.  Consequently,

```math
 \max_{u\ne v\in U_j}{|x_u^Tx_v|\over N_j}
 \le \mu_t+O(j^{-1/2}),
 \qquad
 \mu_t={2\over\pi}\arcsin {4t\over3}.         \tag{MS.18}
```

#### Proof

For even `j`, draw `alpha'` uniformly from balanced sign words and draw
`beta'` uniformly and independently.  Conditional on a fixed `(alpha,beta)`,
write (MS.16) as

```math
 C=\sum_t a_t\alpha'_t+\sum_t b_t\beta'_t,
 \quad
 a_t={3\over4}\alpha_t-{1\over2}\beta_t,
 \quad b_t=\beta_t-{1\over2}\alpha_t.          \tag{MS.19}
```

Here `a_t in {+-1/4,+-5/4}` and
`b_t in {+-1/2,+-3/2}`.  Hoeffding for the independent `beta'` signs and
Hoeffding for sampling without replacement for the balanced `alpha'` signs
give

```math
 \Pr\{|C|>tj\}\le4\exp(-t^2j/25).              \tag{MS.20}
```

For odd `j`, balance the first `j-1` coordinates and fix the last one; its
bounded contribution changes `25` to, for example, `30` for all sufficiently
large `j`.  Sampling `exp(t^2j/100)` endpoints and taking a
union bound over their pairs proves the existence assertion.  Repeated
endpoints may be excluded; their probability is negligible compared with
the `Theta(4^j/sqrt j)` available balanced words.

The two field variances are at least `3j/4`, so (MS.16) bounds the absolute
correlation of the centred score pair by `4t/3`.  Its two-dimensional
covariance has least eigenvalue at least `(3/4-t)j`.  The fixed-dimensional
multivariate Berry--Esseen theorem is therefore uniform.  Applying it to
the four quadrants, allowing for the bounded nonzero field means, and using
the Gaussian sign identity

```math
 \mathbb E[\operatorname {sgn}G\operatorname {sgn}G']
 ={2\over\pi}\arcsin(\operatorname {Corr}(G,G'))
```

proves (MS.18). `square`

The exponent in (MS.17) is not optimized.  Its purpose is to distinguish a
genuine positive packing rate in tensor depth from a merely fixed or
logarithmic number of examples.

## 5. Common-family sparse-flip packing

Fix the parameters `kappa=1/2` and `lambda=1/10` from Theorem 21.64, and
write

```math
 \rho=\sqrt{2/\pi},\qquad s=\sqrt{1-\rho^2},
 \qquad b={\lambda\sqrt7\over2},
```

```math
 \delta=
 \min\left\{b(\rho-s)-{b^2\rho^2\over2\kappa},
 b\rho-{b^2\over2(2-\kappa)}\right\}
 =b(\rho-s)-{b^2\rho^2\over2\kappa}>0.0146.    \tag{MS.21}
```

For each endpoint `u`, choose an even integer `m_u` such
that

```math
 {m_u\|h_u\|_2\over N_j}=b+o(1),
 \qquad m_u=\Theta(\sqrt{N_j/j}),               \tag{MS.22}
```

and set `g(u)=m_uh_u`.  Rounding to the nearest even integer costs only `o(1)`
in (MS.22).

For each `u`, independently as an existence argument, apply the sparse-flip
construction to the selector `x_u`.  Choose a good realization `H(u)` and
hollow it to `A(u)`.  Lemma MS.1 and Theorem 21.62 give, uniformly over all
active products for each member,

```math
 \max_{z\in Z_j}d_z(A(u))=O(j^{-1})+o(1),       \tag{MS.23}
```

while matrix concentration gives

```math
 {H(u)\over\sqrt {N_j}}
 =T_j-\kappa e_ue_u^T+o_{op}(1),
 \qquad e_u={x_u\over\sqrt {N_j}}.             \tag{MS.24}
```

There is no aggregate edge budget across (MS.24): these are different
members of one family, and the same matrix coordinate may be reused in
different members.

### Theorem MS.3 (polynomial-order one-hot packing)

Take `t=1/10` in Lemma MS.2 and construct `A(u)` as above.  Then the common
labelled query bank satisfies, uniformly for `u ne v`,

```math
 \mathcal B_{A(v)}(g(u))-\mathcal B_{A(u)}(g(u))
 \ge(\delta_*+o(1))N_j^{3/2},                  \tag{MS.25}
```

where

```math
 \delta_*=\delta-{\kappa\over2}\mu_{1/10}^2
 >0.012,
 \qquad
 \mu_{1/10}={2\over\pi}\arcsin {2\over15}.    \tag{MS.26}
```

In particular the number of pairwise separated states obeys

```math
 k_j:=|U_j|\ge e^{j/10000}
 =N_j^{1/(10000\log16)}\longrightarrow\infty,  \tag{MS.27}
```

where `log` is the natural logarithm.

#### Proof

Equations (MS.14)--(MS.15) make the proof of Theorem 21.64 uniform over the
balanced endpoint class.  The target-flipped member therefore has

```math
 {\mathcal B_{A(u)}(g(u))\over N_j^{3/2}}
 \le {1\over2}+b\rho-\delta+o(1).              \tag{MS.28}
```

For the cross member `A(v)`, use the explicit Boolean competitor `x_u` in
the positive trust channel.  Equations (MS.24) and (MS.15) give

```math
 {\mathcal B_{A(v)}(g(u))\over N_j^{3/2}}
 \ge {1\over2}-{\kappa\over2}(e_u^Te_v)^2
       +b\rho+o(1).                             \tag{MS.29}
```

Subtract (MS.28) from (MS.29), then use (MS.18).  At `t=1/10`,
`mu_t<0.086`, so (MS.26) follows. `square`

This is stronger than separate base-versus-flip assertions: query `u`
separates state `u` from every other member of the same family, with a
uniform orientation.

## 6. Balanced compilation preserves the whole packing

Choose one common even integer `s_j=Theta(sqrt(N_jj))` with

```math
 s_j\ge\max_{u\in U_j}m_u(2j+1).                \tag{MS.30}
```

For each query `u`, Theorem 21.66 supplies an `N_j times s_j` block `B(u)`
and endpoint `eta_*(u)` satisfying `B(u)eta_*(u)=g(u)` and the uniform
affine error (21.386a).  The error is still

```math
 N_j\sqrt{s_j}+s_j^{3/2}\sqrt {N_j}=o(N_j^{3/2}).             \tag{MS.31}
```

Use any one common hollow exact-sign `C_j` of order `s_j` in (MS.2).

### Corollary MS.4 (unconstrained parents in every fixed context)

For every distinct `u,v in U_j`,

```math
 Q(P_j(v\mid u))-Q(P_j(u\mid u))
 \ge(\delta_*+o(1))N_j^{3/2}.                  \tag{MS.32}
```

Consequently (MS.3) makes `{u:u in U_j}` a
`delta_*`-packing, up to an arbitrarily small fixed reduction of
`delta_*` for all sufficiently large `j`.

#### Proof

For the target-flipped child `A(u)`, the two spherical response bounds from
Theorem 21.64 are increasing with field strength `|a|b` on `0<=|a|<=1`.
The affine compiler estimate and field-`l_1` Lipschitzness therefore give

```math
 Q(P_j(u\mid u))
 \le \left({1\over2}+b\rho-\delta+o(1)\right)N_j^{3/2}
       +Q(C_j).                                  \tag{MS.33}
```

For `P_j(v|u)`, select `eta_*(u)` and use the cross-child witness in
(MS.29).  The shore term is bounded below by `-Q(C_j)`, so

```math
 Q(P_j(v\mid u))
 \ge \left({1\over2}-{\kappa\over2}(e_u^Te_v)^2
             +b\rho+o(1)\right)N_j^{3/2}-Q(C_j).              \tag{MS.34}
```

Now `Q(C_j)=O(s_j^2)=o(N_j^{3/2})`; subtraction proves (MS.32). `square`

## 7. Rates and ceilings for stronger meanings of "simultaneous"

### 7.1 Query-count ceiling for the proved model

There are at most

```math
 2^{(2j+1)-1}=4^j=\sqrt {N_j}                  \tag{MS.35}
```

projective PC.3 endpoints.  Hence any packing which assigns a distinct
native selector query to every one-hot state has `k_j<=sqrt(N_j)`.  Together
with (MS.27), the rigorously established, non-sharp rate window is

```math
 N_j^{1/(10000\log16)}\le k_j\le N_j^{1/2}.     \tag{MS.36}
```

Determining the sharp exponent in this window is a finite-alphabet coding
problem for the covariance kernel (MS.16); it is not decided here.

### 7.2 Disjoint edit-layer edge budget

Suppose a Boolean cube is implemented by fixed, pairwise disjoint unordered
edge layers `F_1,...,F_k`, and toggling coordinate `ell` must change its
selector quadratic energy by at least `cN_j^{3/2}`.  One flipped unordered
edge changes the normalized energy `sum_(i<j)a_ijx_ix_j` by at most two.
Therefore

```math
 |F_\ell|\ge {c\over2}N_j^{3/2},
 \qquad
 k\le {1+o(1)\over c}\sqrt {N_j}.              \tag{MS.37}
```

This is a genuine edge-budget ceiling for disjoint layers.  It is not a
ceiling for the one-hot family above, whose different matrices may reuse
edges, nor for a general cube with overlapping parity layers.

### 7.3 Rank ceiling for canonical simultaneous rank-one superposition

All PC.3 selectors lie in the span of the active products.  The local raw
character Gram matrix on the three-atom support has rank three, so this
span has dimension exactly `3^j`.  Suppose an all-on cumulative sparse-edit
state has the canonical first-order form

```math
 {H_{\rm all}\over\sqrt {N_j}}
 =T_j-\kappa\sum_{\ell=1}^k e_\ell e_\ell^T+o_{op}(1),
 \qquad {\|H_{\rm all}\|_{op}\over\sqrt {N_j}}\le1+o(1).     \tag{MS.38}
```

On the selector span, `T_j=I`.  If
`S=sum_ell e_ell e_ell^T`, the roof in (MS.38) forces

```math
 \lambda_{\max}(S)\le {2\over\kappa}+o(1).     \tag{MS.39}
```

Since `tr(S)=k` and `rank(S)<=3^j`,

```math
 \boxed{k\le\left({2\over\kappa}+o(1)\right)3^j
 =O\left(N_j^{\log_{16}3}\right).}             \tag{MS.40}
```

This is stronger than the disjoint-edge ceiling in its stated canonical
superposition model.  It does not rule out a different nonlinear edit cube
whose all-on mean is not the sum in (MS.38).

## 8. Conclusion and unresolved boundary

The phase-I question has a positive answer in the exact pairwise metric
(MS.3): diffuse sparse-flip channels admit `k_j>=exp(cj)` pairwise states,
and rowwise microcanonical compilation preserves their parent-scale gaps in
every externally selected context.  The construction is counterexample-led:
balanced diffuseness alone is insufficient, while a covariance code makes
the sparse-flip directions genuinely distinguishable.

What remains open is the stronger physical multiplexing problem: construct
one `2^{k_j}` exact-sign edit cube with coordinate-wise decoding, or place
many query compilers into one all-spins-free parent with a selector gadget
whose own energy remains subleading.  Equations (MS.37) and (MS.40) are
ceilings only for the two explicitly stated composition rules and should
not be quoted as general impossibility theorems.
