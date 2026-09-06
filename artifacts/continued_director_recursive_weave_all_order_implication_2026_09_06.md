# A strict recursive-weave certificate would improve the all-order upper bound

**Final campaign update:** the finite variational hypothesis is now proved
for some fixed finite depth. Its unconditional numerical consequence is in
`continued_director_strict_all_order_upper_2026_09_06.md`.

Date: 2026-09-06. Director reconstruction. The finite variational hypothesis
below remains OPEN. This file removes only a potential order-density gap
from its consequence; it is not a convergence theorem.

Use Phi_t and the finite-depth operator B from
`continued_convergence_recursive_orbit_bound_2026_09_06.md`. Suppose, for
some fixed p in (0,1), t>0, integer r>=0, and a>0, one proves

```math
 p\log2+(\mathcal B^r\Phi_t)(\nu_p)+t(1-\sqrt p)\le-a,
 \quad \nu_p=(1-p)\delta_0+\frac p2\delta_{-1/\sqrt p}
                               +\frac p2\delta_{1/\sqrt p}.    (1)
```

Then the original minimax problem satisfies an ALL-ORDER conclusion

```math
 \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
       \le\frac12-\frac{\eta}{2\sqrt p},                (2)
```

for any fixed 0<eta<min(a/t,sqrt(p)). The restriction on eta only avoids
an impossible negative cap bound. In particular, a single strict finite
certificate would lower the current all-order upper bound, not merely
produce an isolated subsequence improvement.

## 1. Uniformity of the finite-depth construction

The terminal orbital estimate is uniform over every real Hadamard matrix
of its terminal order. For fixed depth and initial alphabet, all type and
factorial remainders are o(m), uniformly in those terminal choices. Thus
the proved upper recursion works along ANY unbounded Hadamard-order
sequence s, with m=2^r s and k=floor(pm). Rounding the initial counts is
harmless in the fixed finite-type limit: at finite t the kernel entries
are positive and continuous on the bounded finite alphabet, entropy is
continuous at zero masses, and the relevant compact marginal polytopes
can be rounded with O(1) count changes. No claim uniform in growing r
or t is being made.

With gamma_m=1-sqrt(k/m)+eta, (1) and the exact soft weave theorem give
an exponentially vanishing failure probability, because its exponent
divided by m² is at most -a+t eta+o(1)<0. Hence an actual full sign
matrix K of order N=mk exists with

```math
 \max_x |x^T Kx|\le m^2 k(\sqrt{k/m}-\eta).
```

After deleting its diagonal, A=K-diag(K) is hollow with exact sign
off-diagonal entries, and

```math
 \frac{Q(A)}{N^{3/2}}
 \le\frac12-\frac{\eta}{2\sqrt{k/m}}+\frac1{2\sqrt N}.
                                                               (3)
```

This checks the factor 1/2 between the matrix quadratic form and Q.
The terminal matrices need not be symmetric: the outer weave is symmetric
by its own exact entry formula.

## 2. A relatively dense supply of terminal orders

For each prime q congruent to3 modulo4, let chi be the quadratic character
of F_q, extended by chi(0)=0, and C_xy=chi(x-y). Then
C^T=-C, C1=0, and CC^T=qI-J. The off-diagonal character sum is -1:
for a nonzero shift, count solutions of u²-v²=constant using the bijection
(u-v,u+v); there are q-1 such pairs, giving sum_z chi(z(z+c))=-1.

Consequently

```math
 R=\begin{pmatrix}0&1^T\\-1&C\end{pmatrix},\qquad H=I+R
```

satisfies R^T=-R, RR^T=qI, and HH^T=(q+1)I. Its entries are signs.
This directly verifies the needed Paley terminal construction; no Boolean
extremal property of Paley matrices is used.

The prime number theorem in the FIXED progression3 modulo4 gives
pi(x;4,3)~x/(2 log x). Its hypotheses gcd(3,4)=1 and fixed modulus are
exactly satisfied. In particular consecutive such primes q_j obey
q_(j+1)/q_j->1: otherwise some fixed relative interval would eventually
contradict the asymptotic difference of the two prime counts.

Primary source, Theorem2.1:
[Soprounov, A short proof of the prime number theorem for arithmetic progressions](https://academic.csuohio.edu/soprunov-ivan/wp-content/uploads/sites/93/2023/02/primes.pdf).
Only this qualitative fixed-modulus theorem is imported; no unproved
short-interval estimate or variable-modulus uniformity is required.

Choose s_j=q_j+1, m_j=2^r s_j and N_j=m_j floor(p m_j). Then
N_(j+1)/N_j->1. Equation(3) holds at every sufficiently large such order.
This use of Paley is solely an order-realization device behind a uniform
terminal theorem, not an assumption about optimizing children or their
spectra.

## 3. Principal restriction fills every remaining order

For ANY hollow signing A on N vertices, and any principal n-vertex
restriction A_I, Q(A_I)<=Q(A). Indeed, for a fixed spin assignment on I,
average the outside spins independently with mean zero. The expected
full Hamiltonian is exactly the restricted Hamiltonian; its absolute
value is at most the maximum absolute full Hamiltonian.

For arbitrary n, take the least N_j>=n. Relative density gives N_j/n->1.
Restrict the signing supplied by (3). It follows that

```math
 M_n/n^{3/2}\le [Q(A_{N_j})/N_j^{3/2}](N_j/n)^{3/2},
```

which proves (2). There is no accumulated insertion/bridge defect and no
spectral bound on the restriction is needed.

## 4. What remains

Hypothesis(1) is not proved, and shallow or selected-policy calculations
do not prove it: B is a supremum over all admissible pair laws. Moreover,
even an improved limsup would not show that liminf equals limsup. The
original convergence/nonconvergence question remains open. The significance
here is limited but precise: all-order realization is NOT an additional
open hypothesis for this particular finite-depth recursive certificate.
