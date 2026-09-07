# Independent audit of sign-compatible, trace-budget midpoint surgery

Date: 2026-09-07. Status: PASS. I read and reconstructed Sections 2.1--2.2
of `decisive_audit_low_rank_midpoint_sign_surgery_2026_09_07.md`, including
the extension from projectors to positive contractions. No correction found.

For 0<=T<=Id satisfying the stated ramp, write r=Tr T, t=2I/n and
kappa=max off-diagonal |T_ij|. Sign compatibility and t*kappa<=2 make
p_ij=t A_ij T_ij/2 valid probabilities. The exact mean is
A-t offdiag(T), not a surrogate with an unpaid Schur term. Since A is
hollow,

```math
S=\sum_{i<j}p_{ij}={t\over4}\operatorname{Tr}(AT).
```

Compatibility turns this trace into the sum of the off-diagonal absolute
entries of T. Cauchy--Schwarz and Tr(T^2)<=Tr T give
Tr(AT)<=n sqrt(r), hence S<=I sqrt(r)/2. This does not need an operator
bound on A. The sharper bound Tr(AT)<=||A||op Tr T is also valid.

The centered edge noise has total variance at most 4S and individual
absolute bound 2. Bernstein and the full-cube union bound therefore give
the quoted sqrt(8Sa)+(4/3)a error, a=(n+2)log2. Markov controls the edit
count by 4S; if S=0, no flips occur and the separate zero case is immediate.
The forbidden diagonal costs Ir/n. Substitution yields exactly the
normalized term
2sqrt((I/n^(3/2))(1+2/n)log2)*(r/n)^(1/4).

For the more general masked construction, replacing projection norms by
the norms of T^(1/2)x preserves the masking estimate, because ||T^(1/2)||<=1.
Spectral-square-root columns have total squared norm Tr T, giving the same
Schur and flip-count bounds. Thus the trace-budget extension is genuine.
Entrywise probability feasibility remains an independent assumption.

The new Steiner-family obstruction in
`decisive_independent_steiner_etf_ramp_trace_obstruction_2026_09_07.md`
does not contradict this conditional theorem. It shows that arbitrary
bounded-cap/full-sign parents need not supply ANY small-trace PSD ramp,
before compatibility or entrywise feasibility is imposed. Selection for
actual width minimizers remains open.
