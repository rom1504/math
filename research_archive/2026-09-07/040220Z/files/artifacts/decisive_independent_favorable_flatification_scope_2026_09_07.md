# Favorable flatification: exact scope and a sufficient original-value theorem

Status: the elementary statements and conditional implication below are proved.
The favorable comparison at fixed oversaturation is **not proved**. This note
does not claim original convergence or a new lower endpoint.

For a hollow real symmetric matrix B write

```math
Q(B)=\max_{x\in\{-1,1\}^n}\left|\sum_{i<j}B_{ij}x_ix_j\right|,
\qquad q_n(B)=Q(B)/n^{3/2}.
```

Let M_n be the minimum over full off-diagonal signs and m_n=M_n/n^{3/2}.
A row-regular weighted matrix here means sum_{j!=i} B_ij^2=n-1 for every i.

## 1. Exact favorable flatification already fails at order four

Let H_2=[[1,1],[1,-1]] and set

```math
B_4=\sqrt{3/2}\begin{pmatrix}0&H_2\\H_2^T&0\end{pmatrix}.
```

Every row has squared norm 3, all entries have magnitude at most sqrt(3/2),
and Q(B_4)=sqrt(6), since the infinity-to-one norm of H_2 is 2.
But M_4=4: every flat signing has E_x q_A(x)^2=6, while its cap is an even
integer, so M_4>=4. A signing of K_4 with exactly one negative edge attains 4.
Thus an exact finite statement M_n<=Q(B) is false even at bounded amplitude
and exact row regularity. This does not refute an o(n^{3/2}) error theorem.

There is also a finite example built from two actual optimal children.
The order-five pentagon Seidel signing has cap M_5=4. Consequently
B_10=(3/2)diag(A_5,A_5) has exact row variance 9 and Q(B_10)=12, whereas the
independently stored exact value is M_10=13. This second example uses the
archived exact M_10 certificate; the order-four example is self-contained.

For completeness, order four also admits an exact arbitrary-amplitude
identity. Group its six magnitudes into the three opposite-edge pairs,
and write h_j and l_j for each pair's maximum and minimum. Then

```math
\min_{A_e=\pm1}Q(A\circ w)
=\sum_{j=1}^3 h_j+
 \left|2\max_j l_j-\sum_{j=1}^3l_j\right|.
```

To check it, fix the product p of the four spins. The maximum absolute
energy with this product is the sum over opposite pairs of
|a_e+p a_opposite|. Each pair contributes h_j+p r_j l_j, where its relative
edge sign r_j can be chosen independently. Maximizing over p and minimizing
the three r_j gives the formula. This is a finite diagnostic, not an
asymptotic Gaussian-amplitude comparison.

## 2. Vanishing oversaturation is already covered

Suppose B is row regular and max|B_ij|<=L, where L>=1. Apply independent
unbiased sign rounding to B/L. The proved strong-variance smooth maximum
bound gives, with a_n=(n-1)log(2)/n,

```math
m_n\le {q_n(B)\over L}
       +\sqrt{a_n(1-L^{-2})}.
```

Indeed the total variance of the rounding noise is
n(n-1)(1-L^{-2})/2, and there are 2^n signed spin tests. In particular
L_n=1+o(1) gives m_n<=q_n(B)+o(1). The outstanding issue is **fixed**
oversaturation, already L approximately sqrt(2) for two equal children.

This proof is a direct corollary of
`decisive_independent_soft_flatness_variational_bridge_2026_09_07.md`;
it is not a new universality assertion.

## 3. A strong rectangular consequence of unrestricted flatification

For any full sign m by m matrix C put N=2m and

```math
B=\sqrt{(2m-1)/m}\begin{pmatrix}0&C\\C^T&0\end{pmatrix}.
```

This has exact row variance N-1, maximum amplitude below sqrt(2), and

```math
q_{2m}(B)={\sqrt{1-1/(2m)}\over2}
              {\|C\|_{\infty\to1}\over m^{3/2}}.
```

Thus a uniform favorable theorem m_N<=q_N(B)+o(1) for this class alone
would imply

```math
\liminf_m {\min_C\|C\|_{\infty\to1}\over m^{3/2}}
\ge 2\liminf_m m_{2m}\ge .8666442233281614.
```

The last inequality uses the currently certified original lower endpoint.
Primary references checked directly state the classical asymptotic lower
sqrt(2/pi), and an all-large-order upper 1+epsilon. The exact finite value
G_12=36 gives sqrt(3)/2 at order 12, not an infinite family.

- Dumitrescu, arXiv:2412.16994v1, abstract and original-game discussion:
  <https://arxiv.org/html/2412.16994v1>.
- Pellegrino--Raposo, arXiv:2111.00445v3, introduction and Table 3:
  <https://arxiv.org/html/2111.00445v3>.

This is an implication, not a new lower-bound proof from flatification.
The separate campaign proof
`decisive_audit_bipartite_marked_response_lower_2026_09_07.md`, independently
checked in `decisive_bridge_bipartite_marked_lower_audit_2026_09_07.md`, now
establishes this rectangular lower bound directly from the marked mechanism.
That proof does not need favorable flatification. No claim is made here
about arbitrary weighted-array extension of that argument.

## 4. Only a two-child theorem with summable error is sufficient

The following is a precise sufficient theorem for the original limit, not
an established hypothesis. It needs neither arbitrary weighted parents nor
two-sided amplitude universality.

**Hypothesis.** There are C<infinity and delta in (0,1) such that, whenever
N=m+l is sufficiently large and m,l lie in [N/4,3N/4], one can choose actual
optimal flat signings A_m,A_l for which the row-regular matrix

```math
B=\operatorname{diag}\left(
 \sqrt{(N-1)/(m-1)}A_m,
 \sqrt{(N-1)/(l-1)}A_l\right)
```

satisfies M_N<=Q(B)+C N^{3/2-delta}.

**Conclusion.** m_n converges.

**Proof.** Put u_n=M_n/sqrt(n-1), n>=2. The triangle inequality for the two
block energies, which requires no assertion about their polarity balance,
gives

```math
u_{m+l}\le u_m+u_l+C' (m+l)^{1-\delta}.
```

Fix a sufficiently large integer k. For N>=k^2 write N=qk+r, where
q=floor(N/k), 0<=r<k<=q. Partition N into q leaf sizes k or k+1, exactly r
of the latter. Merge the leaves using a balanced binary tree on their count.
Every pair of child vertex counts has ratio between 1/4 and 3/4: their leaf
counts differ by at most one, and each leaf size is k or k+1. This includes
the worst count ratio 1:2; k>=2 is enough. Thus the recurrence applies to
all internal nodes once k is above its size threshold.

At tree depth j there are at most 2^j nodes, with total size at most N.
Concavity of t^(1-delta) bounds the total defect on that level by
C' 2^(j delta) N^(1-delta). The height is at most ceil(log_2 q), so summing
the geometric series bounds the total defect by C_delta N k^(-delta).
Therefore

```math
\limsup_{N\to\infty}{u_N\over N}
\le \max\left\{{u_k\over k},{u_{k+1}\over k+1}\right\}
       +C_\delta k^{-\delta}.
```

Appending an arbitrary row to an optimal order-k signing gives
M_{k+1}<=M_k+k. Since (k+1)sqrt(k)>k sqrt(k-1), this implies
u_{k+1}/(k+1)<=u_k/k+k/((k+1)sqrt(k)). Choose k along a liminf subsequence
and send k to infinity. Finally u_n/n=m_n sqrt(n/(n-1)), so the limits are
the same. QED.

A nonquantitative o(N^{3/2}) local defect is not automatically summable
under this proof. A power saving, or an appropriate dyadic Dini bound, is
a substantive part of the hypothesis.

## 5. Exact remaining gap

No theorem here produces the favorable comparison at fixed L>1, including
the narrower two-optimal-child class in section 4. The sparse-amplitude
lower >=2/pi from the companion universality obstruction is compatible
with this favorable direction and does not refute it. The finite examples
in section 1 refute only zero-error versions. A norm bound on A-B cannot
settle the two-block claim: every flat filling has a bridge discrepancy of
order N^{3/2}, so genuine cancellation with the child landscapes is needed.

The original convergence problem is still unresolved by this route.
