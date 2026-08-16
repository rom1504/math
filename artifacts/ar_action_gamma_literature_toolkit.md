# Action convergence, rescaled graph limits, and Gamma recovery toolkit

Date: 2026-08-16.

Status: primary-source literature packet.  It records exact quantifier and
normalization boundaries; it does not propose a recovery theorem.

## Verdict

No located theorem gives directed one-profile recovery of a selected signed
action limit by exact symmetric hollow signings on an upward ratio-dense set
of orders.

- action convergence supplies subsequential compactness, not prescribed-order
  inverse realization;
- general finite graphops are not dense even subsequentially in the full
  topology;
- bounded and \(L^p\) graphon sampling work at every order only for a fixed
  norm-controlled kernel, which fails for the \(A/\sqrt n\) scaling; and
- graphon Gamma-limsup results recover vertex states for a kernel sequence
  already assumed to converge, not the kernels themselves.

The closest exact-sign statement is Backhausz--Szegedy's deterministic
representative of the iid full sign model at every order.  Even there the
paper extracts convergence only on subsequences and explicitly leaves the
all-natural-orders question open.  The matrices are nonsymmetric and include
the diagonal.

## Primary theorem cards

1. **Action compactness.** Backhausz--Szegedy, *Action convergence of
   operators and graphs*, Lemma 2.6 gives action-convergent subsequences under
   a uniform \(\infty\to1\) bound.  Theorems 2.9--2.10 give a \(P\)-operator
   limit and compact weak-equivalence-class space under uniform
   \(p\to q\) bounds with \(p<\infty\), \(q>1\).
   [Primary paper](https://doi.org/10.4153/S0008414X2000070X),
   [arXiv](https://arxiv.org/abs/1811.00626).

2. **Graphops are a positive subclass.** Theorem 6.3 of the same paper
   represents a graphop by a symmetric finite positive measure with absolutely
   continuous marginals.  A signed self-adjoint action limit need not be a
   graphop, so graphop density results do not automatically apply.

3. **Closest random-sign result.** In Section 11, \(H_n\) has independent
   \(\pm n^{-1/2}\) entries.  Lemma 11.1 proves action-distance concentration.
   Lemma 11.2 constructs deterministic full sign matrices \(M_n\) at every
   order with \(\|M_n\|_{2\to2}\le3\) and
   \(d_M(M_n,H_n)\to0\) in probability.  Proposition 11.1 still extracts a
   convergent subsubsequence from every infinite order set.  The paper leaves
   open whether all natural numbers are a good sequence.

4. **Finite graphops are not universally dense.** Backhausz--Szegedy Remark
   3.4 cites finite-approximability failures.  Kun--Thom Theorem 1.3 constructs
   almost-free probability-measure-preserving actions of certain
   \(\Gamma\times\Delta\) that weakly contain no finite labelled-graph
   sequence.  This obstructs full-profile density, not necessarily the
   selected one-profile condition here.
   [Kun--Thom](https://arxiv.org/abs/1901.03963).

5. **Endpoint pathology.** Hrusková Theorem 1.4 gives graph sequences bounded
   only in \((\infty,1)\) with noncanonical weakly equivalent action-limit
   representatives, including representatives that need not be self-adjoint
   or positivity preserving.  Stronger \((p,q)\) control is essential.
   [arXiv:2210.10720](https://arxiv.org/abs/2210.10720).

6. **Bounded graphon sampling.** Borgs--Chayes--Lovasz--Sos--Vesztergombi
   Theorem 4.7 gives, for every sample size \(k\),

   ```math
   \delta_\square(U,\mathbf H(k,U))
   \le \frac{10\|U\|_\infty}{\sqrt{\log_2 k}}
   ```

   with probability at least
   \(1-\exp[-k^2/(2\log_2k)]\), together with the simple-graph form for
   \([0,1]\)-valued graphons.
   [Primary paper](https://doi.org/10.1016/j.aim.2008.07.008),
   [arXiv](https://arxiv.org/abs/math/0702004).

7. **Unbounded sampling.** Fekete--Kunszenti-Kovacs Theorems 2.1 and 2.3
   give sampling errors for a fixed \(L^p\) graphon with constants depending
   on \(\|U\|_p\); Corollary 2.5 gives almost-sure convergence for \(p>4\).
   [arXiv:2203.07581](https://arxiv.org/abs/2203.07581).

8. **Sparse \(L^p\) graphons.** Borgs--Chayes--Cohn--Zhao Theorem 2.8 is
   subsequential for upper-regular graph sequences, while Theorems 2.13--2.14
   give compactness and an all-order random recovery for a fixed \(L^1\)
   graphon after sparse \(\rho_n^{-1}\) scaling.  Those realizers contain
   zeros and are not fully supported signs.
   [Primary paper](https://doi.org/10.1090/tran/7543),
   [arXiv](https://arxiv.org/abs/1401.2906).

9. **Dense Ising limits.** Borgs et al. prove graphon limits of dense
   ground-state and free energies for uniformly bounded weights and no
   dominant vertices.  Their Remark 2.16 notes that iid signed couplings are
   trivial under \(n^{-2}\) normalization; the nontrivial spin-glass scale is
   \(n^{-3/2}\).
   [Annals paper](https://doi.org/10.4007/annals.2012.176.1.2).

10. **Gamma recovery of spins, not kernels.** Braides--Cermelli--Dovetta
    Theorems 4.3--4.4 assume supplied dense graph kernels \(W_n\) converge in
    cut norm and then construct finite-label recovery fields.  Proposition
    4.5 transfers constrained minima; Remark 4.2 shows weak \(L^1\) kernel
    convergence is insufficient.
    [Primary PDF](https://www.numdam.org/item/10.1051/cocv/2019029.pdf).

11. **2025 Gamma extension.** Zhang--Scott--Du--Porter Theorem 6.1 and its
    corollaries handle bounded real states and cut-convergent \(L^p\) kernels.
    The kernel sequence remains an input, and the sharp-interface scaling for
    singular graphons is left for future work.
    [arXiv:2408.00422](https://arxiv.org/abs/2408.00422).

12. **Regular weighted discretization.** Le--Jegelka Theorem 2 discretizes a
    regular \(L^2\)-operator at every order with restricted Lipschitz-profile
    error of order \(n^{-1/2}\); Theorem 4 recovers action limits when the
    allowed profile Lipschitz constant grows slowly.  The hypotheses are
    strong regularity assumptions and the finite operators are arbitrary
    weighted matrices, not hollow signs.
    [arXiv:2306.04495](https://arxiv.org/abs/2306.04495).

## Scale calculation

For \(T_A=A/\sqrt n\) acting on the uniform \(n\)-point probability space,
the equivalent step kernel is

```math
K_A=nT_A=\sqrt n A.
```

Thus \(\|K_A\|_p\asymp\sqrt n\) for every fixed \(p\).  Applying bounded
sampling at sample size \(n\) gives a formal error of order
\(\sqrt{n/\log n}\), not \(o(1)\).  More conceptually, raw iid sign noise has
ordinary dense cut norm of order \(n^{-1/2}\); multiplying by \(\sqrt n\)
leaves order-one action noise.  The fluctuation is the target object, not a
rounding error.

The Gamma quantifiers are also reversed.  Existing theorems say

```math
\forall\text{ target vertex state }u\ \exists u_n
\text{ recovering its energy for an already convergent }W_n,
```

whereas directed recovery requires

```math
\forall\text{ finite-model state }f_n\ \exists f\text{ in the limit}
\text{ matching the law of }(f_n,T_nf_n).
```

No theorem in this packet bridges that reversal.

## Exact-sign one-profile boundary

For every hollow signing \(A\), choose a random Rademacher vector \(x\).  Each
coordinate of \(T_Ax\) has second moment \((n-1)/n\) and fourth moment at
most three times the squared second moment.  Paley--Zygmund implies that some
deterministic \(x\) satisfies

```math
\frac1n\#\{i:|(T_Ax)_i|\ge1/2\}\ge\frac1{12}.
```

Every one-profile of the zero operator has output coordinate zero, so the
directed one-profile distance from \(T_A\) to zero is at least \(1/12\).
Exact sign models are therefore not universally dense even in directed
one-profile topology.  The zero operator is not asserted to be a signed
action cluster; the point is that any positive recovery theorem must use the
special structure of the selected cluster, not only abstract operator bounds.
