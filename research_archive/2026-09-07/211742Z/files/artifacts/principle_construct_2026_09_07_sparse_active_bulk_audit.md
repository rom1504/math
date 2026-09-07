# Independent audit: sparse-active balanced bulk

2026-09-07. **PASS**, complete theorem in
`principle_synthesis_2026_09_07_balanced_bulk_sparse_active_bound.md`.

Balanced columns annihilate every constant fibre. For ell nonconstant
fibres only their directed mutual ports survive. The signed reciprocal
swap is a contraction, so the quadratic energy is at most half the
squared port norm. Each physical row map uses at most ell-1 columns;
its Gram matrix has diagonal at most q and at most ell-2 off-diagonal
entries in each row, each bounded by mu. Gershgorin therefore gives
q+(ell-2)mu exactly. This proves the deterministic bound for EVERY seed
and mask, including zero energy when ell<=1.

Conditioned on the full transform and prescribed partition, the exactly-k
choices in distinct groups remain independent. Selected column products
have expectation p times their full inner product, hence zero for distinct
Hadamard columns. Group contribution ranges have length at most2L,
giving2 exp[-u^2/(2mL)]. The union over at most m^3 pairs and columns
supplies the stated O_L(sqrt(m log m)) coherence/excess simultaneously.
Any exact balancing changes a pairwise inner product by at most twice
the sum of its two flip counts; no independence of repair locations is
needed. Thus the coherence premise remains valid after actual sign repair.

Finally q/sqrt(N)=sqrt(p), so ell=o(sqrt(m/log m)) gives a RELATIVE
variance coefficient sqrt(p)/2+o(1). The explicit larger finite range
(ell-2)mu <=(2b-sqrt(p))sqrt(N) is also normalized correctly. The result
is simultaneous over outer seeds and all relevant physical words,
including microscopic variance and arbitrary minority-density scales.

A possible later extension would replace pairwise coherence by a
restricted-isometry upper bound for the actual stratified selector.
Uniform-row Hadamard RIP cannot simply be conditioned on stratification,
whose probability has a leading exponential cost. No such stronger range
is being asserted here.
