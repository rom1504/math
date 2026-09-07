# Exact Gaussian phase: all-order consequence and director audit

Status: proved, using the dependencies listed below. This is an upper bound,
not a proof of convergence. Date: 2026-09-07.

Put p=31/32 and rho=(sqrt(257)-1)/16. Then

```math
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le\frac{4+p\log 2-4(1-\rho)+\tfrac14\log(1-\rho^2)}{8\sqrt p}
\le\frac{7787631971809}{15748015748016}<0.494515125.
```

## Proof dependencies, including the new phase theorem

1. The precision/Schur supersolution and Gaussian-boundary replacement
   identify the deep binary Bellman certificate H with E. See
   `decisive_director_precision_schur_supersolution_2026_09_06.md` and its
   independent scalar/matrix audits. This is not an identification with
   actual ensemble pressure.
2. `decisive_bridge_exact_ternary_gaussian_phase_2026_09_07.md` proves
   E_4(nu_p)=g_4(1). Its analytic reduction does not assume a three-atom
   reproduction alphabet: the Gaussian-mixture KKT derivative has at most
   one positive local maximum. Endpoints cause no exception: at the
   positive source endpoint the zero and negative source terms make the
   derivative strictly negative, and outside the source hull every term
   decreases. The quotient-of-positive-series monotonicity proves the
   claimed derivative shape. Optimizing the three-atom weight is then exact.
3. For lambda<=1/2, unit subgaussianity gives J_lambda=lambda; its conjugate
   optimum has lambda_star<1/2. For lambda in [1/2,4], the directed-interval
   certificate covers every reproduction amplitude in [0,1/sqrt(p)] and
   excludes this entire region by a strict gap. The director independently
   inspected the support proof, mixture-weight elimination, branch tests,
   monotone rectangle bound, and exact comparison constants, then replayed
   the full certificate. Floating arithmetic chooses the subdivision only;
   every accepted rectangle is interval-verified.
4. The all-order weave reconstruction in
   `decisive_bridge_improved_all_order_cap_audit_2026_09_06.md` applies with
   a= -[p log2+g_4(1)+4(1-sqrt(p))]>0. Its cap is 1/2-a/(8sqrt(p)), which
   simplifies to the expression displayed above. Independence gives
   (E Z)^m, not E Z^m. Full-spin counting includes both energy signs;
   deleting the diagonal costs N/2. H2/H12 terminal orders have consecutive
   ratios tending to one, and principal restriction fills all orders.

The limit order is fixed strict margin, then a finite recursive depth,
then all sufficiently large matrix orders, and only finally margin to zero.
There is no interchange requiring a uniform growing-depth theorem.

## Reproducibility and limitations

Replay `computations/decisive_bridge_gaussian_phase_rectangles_2026_09_07.py
--verify` (212505 accepted rectangles, 425009 tree nodes). Replay the exact
rational constant conversion with
`computations/decisive_director_gaussian_phase_cap_2026_09_07.py`.

The new result determines a nontrivial uninformative-channel phase of the
certificate, not merely its numerical optimum over a guessed channel set.
It improves the original universal upper bound. It does not identify the
original minimum, supply a vanishing-loss seed transfer, or settle whether
liminf equals limsup. The lower endpoint 0.4333221116640807 is retained from
its separate proof and certificate, not newly derived here.
