# Hard feedback thresholds without a pointwise variance floor

Date: 2026-09-06. Status: proof extension, submitted for independent audit.
The subset-variance inequality in Section 1 was proposed by the director;
the audit agent sharpened it to the inverse-moment bound below. This note
checks both and records the precise threshold consequences.

## 1. An exact subset lower bound

Let `B` be real symmetric with every row of Euclidean norm one, put
`Q=B²`, and let `p=2r−1` be a positive odd integer. For every subset
`J` of the row indices,

```math
\sum_{i\in J}(B Q^{\circ p}B)_{ii}\ \ge\ \frac{|J|^2}{n}.       (1)
```

This needs neither flat entries nor a bounded operator norm. Write `b_j`
for row `j`, set `w_j=b_j^{\otimes r}`, and define the positive operator
`L=\sum_j w_jw_j^{\mathsf T}` on the `r`-fold tensor space. Let
`P=P_J\otimes I^{\otimes(r−1)}` be the orthogonal projection restricting
the first tensor coordinate to `J`. Symmetry and unit row norms give

```math
\operatorname{Tr}(PLP)
 =\sum_j\sum_{i\in J}B_{ji}^2\|b_j\|^{2r-2}=|J|.
```

Moreover, direct expansion yields

```math
\|PL\|_F^2
 =\sum_{j,k}\langle Pw_j,Pw_k\rangle\langle w_j,w_k\rangle
 =\sum_{i\in J}\sum_{j,k}B_{ij}B_{ik}Q_{jk}^{2r-1}.
```

Since `PLP` has rank at most `n`, is positive semidefinite, and
`\|PLP\|_F\le\|PL\|_F`, Cauchy--Schwarz for its nonzero eigenvalues
proves (1).

Consequently, for nonnegative weights `w_p` of total mass `tau²>0`,

```math
T=B\left[\sum_{p\ {m odd}}w_pQ^{\circ p}\right]B,
\qquad
\sum_{i\in J}T_{ii}\ge\tau^2|J|^2/n.                         (2)
```

Finite mixtures follow by addition and infinite mixtures by monotone
convergence of the nonnegative quadratic forms. In particular,

```math
T_{ii}\ge\tau^2/n>0,\qquad
\#\{i:T_{ii}\le\varepsilon\tau^2\}\le\varepsilon n.          (3)
```

The bound `tau²/n` is not a fixed positive floor. The apex signing in
`continued_feedback_apex_vanishing_variance_2026_09_06.md` is compatible
with it. The important new fact is that small variances occupy a small
fraction of the roots.

There is a stronger exact conclusion. Let `Pi` be the orthogonal
projection onto the range of `L`, and `P_i` the singleton version of
`P`. Then `L Pi=L`, and

```math
1=\langle P_i\Pi,P_iL\rangle_F
 \le\|P_i\Pi\|_F\|P_iL\|_F
 =\sqrt{d_i(BQ^{\circ p}B)_{ii}},\qquad
d_i=\operatorname{Tr}(P_i\Pi).
```

Since the `P_i` partition the identity, `sum_i d_i=rank(L)<=n`.
Therefore `sum_i 1/(B Q^{circ p}B)_ii<=n`. Convexity of inverse,
applied to the normalized mixture weights, gives

```math
\frac1n\sum_i\frac{\tau^2}{T_{ii}}\le1,\qquad
\frac1n\sum_i\sigma_i^{-1}\le\tau^{-1},\qquad
\frac1n\sum_i\sigma_i\ge\tau.                              (3a)
```

For an infinite mixture, inverse convexity is applied directly to its
probability measure. All terms are positive by (1). In particular, for
any shifts `v_i`, deterministic or random and independent of the scalar
Gaussian,

```math
\frac1n\sum_i\mathbb P(|v_i+\sigma_iN|\le\eta)
 \le\sqrt{2/\pi}\,\eta/\tau.                               (3b)
```

This improves the cutoff-based small-ball estimate used below; the
cutoff is still useful when weighting the RAW actual channel.

## 2. Zero-first masked feedback: full threshold covariance

Use the hypotheses and notation of
`continued_feedback_zero_first_masked_covariance_2026_09_06.md`, with
`tau²=E f(N)²>0`, but now set `psi=sign`, with `sign(0)=0`. Thus

```math
Y=Bf(BS),\qquad C_i=H((BS)_i)\operatorname{sign}(Y_i),
\quad\sigma_i^2=T_{ii}.
```

Let independent Gaussians `U~N(0,Q), Z~N(0,T)` and put

```math
K_{ij}=\mathbb E[H(U_i)H(U_j)]\,
        \mathbb E[\operatorname{sign}(Z_i)\operatorname{sign}(Z_j)].
```

Then the previously proved smooth closure extends to

```math
\|\mathbb E CC^{\mathsf T}-K\|_*/n\longrightarrow0.             (4)
```

Here is the ordered approximation, including the small-variance rows.
For fixed `epsilon>0`, discard `J_epsilon={i:sigma_i²<=epsilon tau²}`.
Their number is at most `epsilon n`. Replacing a uniformly bounded
random vector on those coordinates changes its covariance by
`O(n sqrt(epsilon))` in nuclear norm. On the remaining roots, use a
fixed smooth odd function `psi_eta` bounded by one and equal to `sign`
outside `[-eta,eta]`. The one-root limit from the smooth proof gives
the Gaussian small-ball bound

```math
\limsup_{n\to\infty}\frac1n\sum_{i\notin J_\varepsilon}
 \mathbb P(|Y_i|\le\eta)
 \le C\eta/(\tau\sqrt\varepsilon),
```

where a slightly larger smooth interval majorant can be used at the
endpoints. Thus the averaged squared response error is at most
`C(epsilon+eta/(tau sqrt(epsilon)))` after `n→infinity`. The Gaussian
comparison has the same bound. First send `n→infinity`, next
`eta→0` at fixed `epsilon`, and finally `epsilon→0`. This proves (4).
Alternatively (3b) gives a direct `O(eta/tau)` averaged Gaussian
small-ball bound, and the smooth local comparison transfers its smooth
majorant without first discarding rows.

The resulting cutoff-free energy formula is

```math
\frac{\mathbb E C^{\mathsf T}BC}{2n}
 =\frac{\mu^2}{\pi n}
   \operatorname{Tr}(B D_{1/\sigma}T D_{1/\sigma})+o(1),
\qquad \mu=\mathbb E H(N).                                   (5)
```

To justify the unbounded diagonal notation in (5), first prove it with
`D_{1/sigma}` set to zero on `J_epsilon`. The full normalized Gaussian
vector has covariance `D_{1/sigma}T D_{1/sigma}` with diagonal exactly
one. Removing `epsilon n` coordinates therefore changes this matrix
by `O(n sqrt(epsilon))` in nuclear norm. This does NOT require a bound
on `||D_{1/sigma}||op`.

For feasible endpoints `H>=0, |f|+H<=1`, (5) and the local gain give

```math
\Lambda(B)\ge\sqrt{2/\pi}\,\mu\,\frac1n\sum_i\sigma_i
 +\frac{\mu^2}{\pi n}
  \left|\operatorname{Tr}(B D_{1/\sigma}T D_{1/\sigma})\right|-o(1).
                                                                    (6)
```

Thus the earlier threshold formula is valid without a pointwise
variance-floor assumption. It still need not improve the universal
lower constant: the common energy may vanish.

## 3. Nonzero-first response: a cutoff-safe energy identity

Use the hypotheses and notation of
`continued_feedback_coherent_return_energy_projection_2026_09_06.md`:
`G=BS`, `V=QS`, `f=b h_1+r`, `Z=Br(G)`, and
`T=B[sum_{p>=3 odd} f_p² Q^{circ p}]B`. Suppose `tau²=E r(N)²>0`.
Define the actual hard feedback and its literal coherent smoothing by

```math
C_i=H(G_i)\operatorname{sign}(bV_i+Z_i),
\qquad
c_i^0=H(G_i)\left[2\Phi(bV_i/\sigma_i)-1\right].
```

Let

```math
a_i=\sqrt{2/\pi}\,\frac1{\sigma_i}
       \mathbb E_S\!\left[H(G_i)
       \exp\{-b^2V_i^2/(2\sigma_i^2)\}\right],
\qquad
a_i^\varepsilon=a_i\mathbf1_{\{\sigma_i^2>\varepsilon\tau^2\}}.
```

Write `e_n(W)=E W^T B W/(2n)`. For every fixed `epsilon>0`,

```math
\limsup_{n\to\infty}\left|
 e_n(C)-e_n(c^0)
 -\frac{\mathbb E(c^0)^{\mathsf T}B D_{a^\varepsilon}Z}{n}
 -\frac{\operatorname{Tr}(B D_{a^\varepsilon}T
                                    D_{a^\varepsilon})}{2n}
 \right|\le C_{L,H}\sqrt\varepsilon.                         (7)
```

The trace in (7) can equivalently use untruncated `a`: indeed
`|a_i|sigma_i<=sqrt(2/pi)||H||infinity`, so its Gaussian covariance
changes by `O(n sqrt(epsilon))`. We deliberately do NOT remove the
cutoff from the cross with the RAW actual vector `Z`. A normalized
nuclear covariance error for `Z` cannot silently be multiplied by an
unbounded `D_a`. Equation (7), with `n→infinity` before
`epsilon→0`, is the precise general energy conclusion.

For completeness, apply the smooth theorem to `psi_eta` as in Section
2. Conditional on the actual coherent pair `(G_i,V_i)`, the comparison
Gaussian has density bounded by `1/(sqrt(2 pi) sigma_i)`, so the same
small-ball estimate works on the retained roots regardless of the
possibly non-Gaussian law of `V_i`. The corresponding coherent
response and derivative converge uniformly there as `eta→0`.
On discarded roots, the smooth derivative coefficient `a_{i,eta}`
satisfies

```math
|a_{i,\eta}|\sigma_i
 =|\mathbb E_{S,N}H(G_i)N\psi_\eta(bV_i+\sigma_iN)|
 \le\sqrt{2/\pi}\|H\|_\infty.                               (8)
```

For each fixed `eta`, its magnitude is bounded independently of `n`.
The already proved normalized nuclear covariance comparison for raw
`Z` therefore transfers the weighted squared-norm bound to actual
`D_{a_eta}Z`. After `n→infinity`, the discarded coordinates contribute
at most `C epsilon` to that normalized squared norm. Cauchy--Schwarz
and `||B||op<=L` give their `O(sqrt(epsilon))` energy and cross costs.
This proves (7) by the indicated order of limits.

The local gain has the simpler cutoff-free form

```math
\frac1n\mathbb E\sum_i H(G_i)|bV_i+Z_i|
 =\frac1n\sum_i\mathbb E_{S,N}H(G_i)|bV_i+\sigma_iN|+o(1).      (9)
```

It follows directly by uniform integrability and the local comparison;
absolute value is continuous and Lipschitz. The exact feasible-endpoint
identity combines (7), (9), and the old self-energy `b² m_3/2`.
The same proof works for the fixed-dimensional independent-colored-seed
extension, replacing `bV` by its literal summed coherent return. It
does not yet prove a theorem for genuine dependent old-tree fields.

## 4. Scope

No pointwise fixed positive variance floor has been asserted. The exact
subset inequality replaces that unnecessary hypothesis. The conclusions
are feedback covariance/energy identities in the already audited scalar
and independent-colored-seed classes. They establish neither spectral
flatness of minimizing signings nor a new universal bound nor convergence
of the original extremal sequence.
