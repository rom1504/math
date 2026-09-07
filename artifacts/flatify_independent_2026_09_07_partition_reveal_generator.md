# Hidden balanced-partition reveal: exact profiles, generator, and finite walls

Date: 2026-09-07. Status: analytic derivation plus exhaustive-signing numerical
diagnostics. No original convergence or nonconvergence theorem follows yet.

The root proposed revealing one initially hidden balanced partition, avoiding
the invalid step of freely reselecting a deterministic interpolation's target
partition midway. This note derives that specific operation before imposing
any martingale sign assumption.

## 1. A row-regular variance martingale with exact child endpoints

Let N=m+n, with m,n>=2, and select an m-set S uniformly. The terminal
variance profile is

```math
V_{ij}(S)=\begin{cases}
(N-1)/(m-1),&i,j\in S,\\
(N-1)/(n-1),&i,j\notin S,\\
0,&\text{otherwise}.
\end{cases}                                               (1)
```

Every profile has off-diagonal row sum N-1, and E V_ij=1. Reveal membership
labels in a fixed vertex order. The posterior mean variance v_ij is a
martingale, remains exactly row regular, and is bounded by a fixed constant
for comparable child sizes. It starts at the full uniform profile and ends
at one genuine two-child profile, without changing the selected partition.

For the explicit posterior formula, let k labels have been revealed, of
which r are positive; write u=N-k, R=m-r, a=(N-1)/(m-1), b=(N-1)/(n-1).
Known equal-sign pairs have value a or b; known opposite pairs have value
zero. Known-positive/unknown and known-negative/unknown pairs have values
aR/u and b(u-R)/u. Two distinct unknown vertices have value

```math
{aR(R-1)+b(u-R)(u-R-1)\over u(u-1)}.                      (2)
```

Terms involving fewer than two unknown vertices are simply absent. These
profiles and their row sums are checked with exact Fractions in the replay.

Use the natural normalization

```math
\Psi_N(\beta)=\min_A {1\over2}\log(Z_A^+Z_A^-),\qquad
Z_A^\sigma=\sum_x\exp\left({\sigma\beta H_A(x)\over\sqrt{N-1}}\right).
```

At a posterior profile v the same definition uses coefficients
A_ij sqrt(v_ij). The initial optimized pressure is Psi_N(beta); the final
one is EXACTLY Psi_m(beta)+Psi_n(beta). No endpoint temperature correction
is needed with sqrt(N-1) normalization. Consequently a one-sided o(N)
integrated reveal comparison at every fixed beta would give the corresponding
almost-additive pressure route. It is not established here.

## 2. The exact one-step innovation is low-rank modulo its diagonal

Choose one unknown vertex w to reveal. Let U be the other u-1 unknown
vertices, and define q_i=a on known positive vertices, q_i=-b on known
negative vertices, and q_i=0 elsewhere. Put

```math
h={a(R-1)-b(u-R-1)\over u-1}.
```

Let delta=v^(positive)-v^(negative). Its nonzero off-diagonal entries are:

- delta_wi=q_i for known i;
- delta_wj=h for j in U;
- delta_ij=-q_i/(u-1) for known i and j in U;
- delta_jl=-2h/(u-2) for distinct j,l in U, when u>2.

Thus delta is a sum of matrices supported on the span of e_w,q,1_U,
apart from a diagonal correction irrelevant to the quadratic Hamiltonian.
It has rank at most three modulo that diagonal. Each row sum is exactly
zero; in particular sum_i q_i=-(u-1)h. The reveal probabilities are R/u
and 1-R/u, so its conditional covariance is
(R/u)(1-R/u) delta tensor delta.

This structure is stronger than an arbitrary row-regular perturbation.
However multiplying it entrywise by the current optimizing signing does
not preserve matrix rank. No operator-rank bound for the actual physical
quadratic observable is inferred.

## 3. Continuous posterior reveal and the exact branch generator

For any continuous posterior realization, write dv_e=sum_k Gamma_ek dW_k.
Localize away from zero variance, and let A be an active optimizing signing.
Put lambda_e=beta sqrt(v_e)/sqrt(N-1), B_e(x)=A_e x_i x_j, and

```math
d_e={E_+ B_e-E_- B_e\over2},\qquad
K_{ef}={\operatorname{Cov}_+(B_e,B_f)+
                 \operatorname{Cov}_-(B_e,B_f)\over2}.
```

K is PSD. For the smooth branch F_A(v), differentiation gives

```math
\partial_e F_A={\beta d_e\over2\sqrt{N-1}\sqrt{v_e}},
\partial_{ef}F_A=
 {\beta^2 K_{ef}\over4(N-1)\sqrt{v_ev_f}}
 -\mathbf1_{e=f}{\beta d_e\over4\sqrt{N-1}v_e^{3/2}}.     (3)
```

Define the physical quadratic observables

```math
L_k(x)=\sum_e {\Gamma_{ek}\over\sqrt{v_e}}B_e(x).
```

The smooth branch's Ito drift is exactly

```math
\mathcal D_A=
 {\beta^2\over16(N-1)}\sum_{\sigma=\pm,k}
                    \operatorname{Var}_\sigma(L_k)
 -{\beta\over8\sqrt{N-1}}
       \sum_e{d_e\|\Gamma_{e,:}\|^2\over v_e^{3/2}}.       (4)
```

Taking the minimum of finitely many smooth branches additionally introduces
a NONPOSITIVE switching local-time term. This follows by repeated use of
min(X,Y)=(X+Y-|X-Y|)/2 and Tanaka's formula. Therefore a putative lower
Jensen/submartingale bound must control this switching loss as well as (4).
For an upper/supermartingale comparison the loss helps, but (4) itself has
no established sign. At an actual edgewise optimum d_e<=tanh(lambda_e),
yet this bound does not compare the physical variances in (4) with the
weighted radial term.

The elementary identity

```math
E\int {\|\Gamma_{e,:}\|^2\over v_e}\,dt
 =2[E(v_e\log v_e)_{\rm final}-(v_e\log v_e)_{\rm initial}]
```

is an entropy-type budget for the posterior variance process. Summed over
edges it has order N^2 for a nontrivial balanced terminal split. On its own,
the prefactor beta^2/(N-1) in (4) leaves an order-beta^2 N budget, not the
required subleading one. It does not supply the missing comparison.

## 4. Exact N=4: local drifts have both signs, and their first terms cancel

Take m=n=2, put s=beta/sqrt(3), and subtract the common 4log2. Initially,
exhausting the eight switching classes gives

```math
F_0=2\log\cosh s+\log\cosh(2s).
```

The first revealed label carries no partition information. At the next
reveal, equal labels have probability 1/3 and determine the final matching;
opposite labels have probability 2/3 and leave a K_(2,2) profile of variance
3/2 on its four edges. The exact optimized pressures are

```math
F_{\rm end}=2\log\cosh(\sqrt3 s),\qquad
F_{\rm mid}=\log\cosh(\sqrt6 s).                         (5)
```

For the cycle, its negative edge-product minimizes the paired partition
function, and 4logcosh(y)+log(1-tanh(y)^4)=logcosh(2y) proves (5).

The first informative mean drift is
(F_end+2F_mid)/3-F_0<0, whereas the later conditional drift
F_end-F_mid>0, for every beta>0. To verify the first sign, apply concavity
of z -> logcosh(sqrt(z)) to the majorization

```math
(3/2,3/2,3/4,3/4,0,0,0,0,0)
 \succ (1,1,1,1/4,1/4,1/4,1/4,1/4,1/4),
```

after rescaling s. The second sign is the same concavity applied to two
equal summands. Strictness holds for beta>0.

Their leading expansions in beta are respectively -beta^4/9 and
+beta^4/6; the latter branch has probability 2/3. Hence the fourth-order
terms cancel in the full endpoint comparison. The remaining total is

```math
F_{\rm end}-F_0=-{4\over405}\beta^6+O(\beta^8)<0
```

near zero. In fact it is negative for every positive beta: use the product
cosh(z)=product_(j>=0)(1+4z^2/[pi^2(2j+1)^2]) and, for each positive u,

```math
(1+3u)^2 < (1+u)^2(1+4u),
```

whose right side minus left side equals 4u^3. Thus literal local
submartingale and supermartingale hypotheses are both false already on
actual finite global minimizers. This is not a scalable obstruction to a
sublinear total error: the cancellation itself is a reason not to bound
every adverse step separately.

## 5. Numerical exhaustive diagnostic and remaining obligation

`computations/flatify_independent_2026_09_07_partition_reveal.py` exhausts all
switching-gauged signings at N=4,6,7 and five temperatures, and checks every
posterior row sum and state probability exactly. Its optimized pressures
are floating point, not transcendental certificates. Results are preserved
in `computations/results/flatify_independent_2026_09_07_partition_reveal.json`.

At N=6 and N=7 the total endpoint difference is positive at all five tested
temperatures, unlike N=4. Several individual expected reveal drifts change
sign along the same path. For example at N=6,beta=3 the nonzero mean drifts
are approximately .65307,.84513,-.19290,.36457. Nothing here supports a
universal finite endpoint sign.

The genuine surviving target is a cumulative, optimizer-specific estimate
for (4) together with switching local time, or a discrete equivalent that
retains the cancellations in Section 4. Merely revealing a fixed hidden
partition removes the target-reselection gap but does not by itself prove
that estimate. Width-to-absolute comparison must also be kept explicit in
any eventual application to the original sequence.
