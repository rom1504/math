# Spectral cores with row-variance completion

2026-09-07. **Proved lemma; independent audit requested.** This is a
weighted-matrix regularization operation, not sign flatification. The
consequence for the marked lower bound depends on the separately audited
bounded-operator weighted extension.

## Statement

Let W_N be symmetric and hollow, with max entry at most a fixed K, and
row squared norms N-1+o(N), uniformly. Suppose its quadratic cap is at
most C N^(3/2). Fix epsilon>0 sufficiently small in terms of K.

There is a principal set I of size k>=(1-epsilon)N and a symmetric hollow
weighted matrix V_k such that, as N tends to infinity at fixed epsilon,

```math
\begin{aligned}
\max_{ij}|V_{ij}|&\le K+O_K(\sqrt\epsilon),\\
\|V\|_{\rm op}&\le O_{C,K,\epsilon}(\sqrt k),\\
\sum_j V_{ij}^2&=k-1+O_K(\sqrt{\epsilon N\log N})\quad\text{uniformly},\\
\mathcal W(V)&\le(1+O_K(\epsilon))\mathcal W(W_N)
                  +O_K(\sqrt\epsilon\,N^{3/2}),
\end{aligned}
```

where the displayed row-error constant is for sufficiently large N at
fixed epsilon, and W denotes energy half-width. The same upper estimate
holds with cap Q in place of half-width. No constraint on the shape of
the entrywise variance profile is required.

The construction preserves the desired lower-bound limit order:
first N tends to infinity at fixed epsilon and fixed operator bound,
then epsilon tends to zero. In particular a universal fixed-operator
lower bound for this weighted row-regular class extends to its entire
bounded-cap class, with no uniformity in the operator constant required.

## 1. Spectral deletion

For every hollow symmetric W, beta(W)<=4Q(W) by cube polarization.
The finite SDP/Grothendieck diagonal-majorant proof in
`resumed_bound_audit_minimal_proof_2026_09_06.md`, Section 1, applies to
arbitrary real coefficients: D>=W,-W and Tr D<=K_G beta(W).
Delete coordinates above K_G beta(W)/(epsilon N). The retained W_0 has

```math
\|W_0\|_{\rm op}\le4K_G C\epsilon^{-1}\sqrt N.
```

Write q_i=sum_j (W_0)_{ij}^2. Bounded entries and deletion give
q_i=N-1+O_K(epsilon N)+o(N), uniformly. Set

```math
q_+=\max_iq_i,\quad q_-=\min_iq_i,\quad
D_0=2(q_+-q_-)+\epsilon N,\quad
\gamma^2=\frac{k-1}{q_++D_0}.
```

Thus gamma=1+O_K(epsilon)+o(1). Define d_i=k-1-gamma^2 q_i and
d_0=(k-1)D_0/(q_++D_0). Exactly,

```math
d_0\le d_i\le\tfrac32d_0,\qquad
d_0=O_K(\epsilon N)
```

for all sufficiently large N at fixed epsilon. The uniform original
o(N) row error is absorbed at this fixed epsilon; no rate is assumed.

## 2. An exact nonnegative symmetric deficit realization

For k>=4 and S=sum_i d_i put, for i!=j,

```math
u_{ij}=\frac{d_i+d_j}{k-2}-\frac{S}{(k-1)(k-2)},\qquad u_{ii}=0.
```

Then sum_j u_ij=d_i exactly. Moreover
2d_0>=S/(k-1), since S<=3kd_0/2 and k>=4. Hence u_ij>=0 and
max u_ij=O_K(epsilon). This is the reason for adding a common positive
variance deficit before attempting completion. Arbitrary nonnegative
deficits alone need not have a symmetric realization.

Choose independent signs on unordered pairs and set R_ij=sign_ij sqrt(u_ij).
Let V=gamma W_0+R. The diagonal is zero exactly. Its squared row norm is

```math
\sum_jV_{ij}^2=k-1+2\gamma\sum_j(W_0)_{ij}R_{ij}.
```

The cross sum has mean zero and subgaussian variance proxy at most
(max u_ij)q_i=O_K(epsilon N). A union bound over rows proves the asserted
uniform row error with probability tending to one.

## 3. Simultaneous operator and objective control

Let U=max u_ij. For every real unit x,
x^T R x=2 sum_{i<j}sqrt(u_ij)sign_ij x_i x_j is subgaussian with
variance proxy at most 2U. A 1/4-net of the unit sphere of size at most
9^k and the symmetric net inequality therefore give

```math
\|R\|_{\rm op}=O(\sqrt{Uk})=O_K(\sqrt{\epsilon N})
```

with probability at least 1-exp(-c k), by choosing a sufficiently large
absolute constant. This event and the row event have nonempty intersection.
Consequently Q(R)<=k||R||op/2=O_K(sqrt(epsilon)N^(3/2)).

Half-width is a seminorm on hollow quadratic forms. Thus
W(V)<=gamma W(W_0)+W(R), and principal averaging gives W(W_0)<=W(W_N).
Also W(R)<=Q(R). This proves the theorem. The cap version is identical.

## 4. Lower-bound consequence and nonclaims

Suppose a constant c is a proved lower bound on W(V_k)/k^(3/2) for every
fixed-entry-bound, fixed-operator-bound sequence with uniform row variance
k-1+o(k). Applying that theorem at fixed epsilon and using Section 3 gives

```math
\liminf_N\frac{\mathcal W(W_N)}{N^{3/2}}
\ge c(1-\epsilon)^{3/2}-O_K(\sqrt\epsilon).
```

Sending epsilon to zero gives c for the original bounded-amplitude class.
If the half-width is potentially below c, Q<=2W supplies the bounded-cap
hypothesis. Thus unbounded Q is not an omitted case.

The output is weighted and has only asymptotically equal row variance.
It is not a full sign matrix, does not flatten amplitudes, and does not
transfer the unknown optimized cap between orders. The error is small
only after choosing epsilon small, not a power-saving sign realization.
