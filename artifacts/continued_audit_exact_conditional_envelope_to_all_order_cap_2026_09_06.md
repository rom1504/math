# Fresh independent audit: exact conditional envelope to an all-order signing cap

Date: 2026-09-06. The final verification continued for a few minutes after
the campaign checkpoint at the director's explicit request.

**Verdict: PASS.** Independently reconstructing the numerical certificate
and the actual signing implication gives

`limsup_n M_n/n^(3/2) <= 1/2-a/(8 sqrt(31/32))`,

where

`a=91470529542342299/20460000000000000000>0`.

The right side is approximately `0.49943222048540312`. This is a strict
all-order upper-bound improvement. It does not establish convergence of
`M_n/n^(3/2)` and does not change the preserved universal lower bound.

## 1. Files read and exact independent replay

The audit freshly read the relevant proofs, rather than relying only on
earlier PASS messages:

- `continued_convergence_conditional_variance_supersolution_2026_09_06.md`;
- `continued_feedback_conditional_variance_exact_certificate_2026_09_06.md`;
- `continued_convergence_restricted_weave_2026_09_06.md`, Sections 1--4;
- `continued_convergence_recursive_orbit_bound_2026_09_06.md`, Sections 1--5
  and its safe symmetry restriction;
- `continued_director_recursive_weave_all_order_implication_2026_09_06.md`;
- the previously independently reconstructed terminal-gap and unbounded
  Sinkhorn-rigidity arguments, with their finite-depth Gaussian-boundary
  consequence.

Both the complete exact certificate script and its imported elementary
logarithm, entropy, square-root and hull routines were read. An independent
run was performed with

```sh
.venv/bin/python computations/continued_feedback_conditional_variance_exact_certificate_2026_09_06.py --output tmp/continued_audit_conditional_variance_exact_certificate_2026_09_06.json
```

The independently obtained exact outputs were

`grid upper = -13307450618189/1056000000000000`,

`continuity payment = 630156538579809/77500000000000000`,

`total upper = -91470529542342299/20460000000000000000`.

No optimizer or floating transcendental value entered this verification.

## 2. Exact numerical certificate reconstruction

For `p=31/32`, `t=4`, the full conditional-variance posterior envelope
has offset

`-h(p)+t(1-sqrt(p))+sup_(E z=p) E q(z,s)`,

where

`q(z,s)=h(z)+z h((1+s)/2)+g_t((z-z^2s^2)/p)`.

Reflection symmetrization leaves precisely the one barycentre constraint
`E z=p`; it does not impose an unwarranted fixed posterior bias. The
least concave majorant over z therefore evaluates the full posterior
optimization after maximizing over s.

The checker uses all 2501-by-2501 mesh pairs. For variance
`v=(i/N-i^2j^2/N^4)/p`, its lookup bin is exactly

`floor[32000(i N^3-i^2 j^2)/(31 N^4)]`.

The Gaussian potential decreases with variance, so using the lower
variance bin gives an upper value. Its potential is enclosed from the
exact formula

`r=(sqrt(1+16t^2v^2)-1)/(4tv)`,
`g_t(v)=-tv(1-r)+(1/4)log(1-r^2)`.

The upper r endpoint is used in the linear term, and the lower endpoint
in the logarithmic term. Both signs are correct. Zero variance is
handled separately and exactly. The atanh-series logarithm remainder,
integer-square-root bounds, and upward integer rounding all reconstruct.

The largest possible array numerator is bounded by
`32*1000*N^4<2^63`; the entropy multiplication is bounded by
`N*10^12<2^63`. All intermediate signs, indices and divisions were
checked. The decreasing-slope hull stack is the upper, not lower, hull.

To pass to every continuous posterior, round z randomly to adjacent
mesh points, preserving its mean, and s to its nearest point. The entropy
payment is
`h(1/N)+(log 2)/N+h(1/(4N))`. The variance derivatives obey
`|partial_z v|<=1/p`, `|partial_s v|<=2/p`, and `|g_t'|<=t`; thus their
combined payment is `2t/(pN)`. These are exactly the four terms paid in
the output above. No extra variance-bin payment is required, since each
finite grid value already uses the optimistic bin potential. Conversely,
that optimistic lookup was not mistaken for an exact point evaluation.

## 3. From the envelope to a strict finite-depth Bellman certificate

Write `T_t(nu)=sup_L[E g_t(Var(X|L))-I(X;L)]`. The independent finite-source
proof of `B T_t<=T_t` is in
`continued_audit_temperature_alignment_boundary_2026_09_06.md`, Section 5.
It retains the correct sequential labels, Schur-complement comparison,
log-variance concavity, eigenvalue majorization, and all information costs.
Finite labels suffice on every state reachable from the ternary source.

The trivial label gives `G<=T_t`, with `G(nu)=g_t(m_2(nu))`. Monotonicity
and the supersolution inequality give `B^r G<=T_t` at every finite depth.
The independently proved Gaussian-boundary replacement says

`lim_r B^r Phi_t = lim_r B^r G`.

Hence the limiting Phi-boundary Bellman offset is at most `-a`. For each
strict `a'<a`, some **fixed finite depth** r therefore satisfies

`p log 2+(B^r Phi_t)(nu_p)+t(1-sqrt(p)) <= -a'`.

This deduction does not require identifying the Bellman limit with the
old latent envelope, making the old envelope a supersolution, or bounding
the required depth effectively. It is essential not to replace `a'` by a
claimed attained finite-depth margin a. The endpoint a enters only after
the final limiting argument.

## 4. Fresh audit of the actual weave and its soft probability bound

The full matrix has exact sign entries

`K_(i,a),(j,b)=S_ij H_i(a,j)H_j(b,i)`.

It is symmetric even when the individual terminal Hadamards are not.
After keeping k coordinates in each of m fibres, its order is `N=mk`.
For every Boolean block assignment, define
`h_i=H_i[T_i,:]^T x_i`. Orthogonality gives

`sum_i ||h_i||^2=m^2 k`,
`E=x^T K x=sum_ij S_ij h_i(j)h_j(i)`,
`D_sigma=2(m^2k-sigma E)`.

The exponential Markov variable is exactly
`exp[-t D_sigma/(2k)]`. Each undirected off-diagonal edge occurs twice
in D, so averaging its independent sign gives precisely the kernel

`[exp(-t(a-b)^2/k)+exp(-t(a+b)^2/k)]/2`.

There is no missing factor two. Dropping the diagonal terms increases
this exponential. After conditioning on the diagonal coordinate, the
remaining row spectra are independent uniform multiset permutations.
Factoring the PSD kernel and applying graph Cauchy--Schwarz gives the
product of square roots of normalized deleted permanents. The squared
vertex norm is exactly the permutation-average permanent, with divisor
`(m-1)!`, including repeated-coordinate multiplicities.

Summing over **all** spins and the two energy signs yields

`P(failure)<=2 exp(t gamma m^2) product_i Z_i(t)`.

The randomized recursive basis is sampled independently in each fibre,
and its fresh output column permutation and edge signs are independent
of that basis. Conditioning first on all bases proves the same inequality
for arbitrary different bases. Averaging then gives `(E Z)^m`, not
`E Z^m`. The square root stays inside the one-row expectation throughout;
the invalid annealed-square-root shortcut is not used.

## 5. Terminal orbital bound, pair types and order of limits

The Gaussian-Fock orbit covariance has operator norm exactly the orbit
kernel average, because its finite orbit Gram matrix is nonnegative with
constant row sum. Degree truncation preserves both orthogonal and signed-
permutation actions. The invariant rank is bounded by the partition
generating function and is `exp(O(sqrt(m)))`; the omitted Poisson tail
is dominated by the full orbital norm's elementary lower bound. Thus the
terminal orbital loss is subexponential, uniformly in the terminal
Hadamard. The independent fresh signed permutation is on its input side.

For concatenated child vectors, the full invariant subspace is contained
in the product-block invariant subspace. This gives
`L(v_+,v_-)<=L(v_+)L(v_-)` in the needed direction. Fixing a diagonal
permutation coordinate gives the deletion payment `sqrt(2m)`.

The signed input-word count is
`m! 2^(#nonzero)/product_abs count!`; the ordered-pair-table count is
`(m/2)!/product_pair count!`. Their logarithmic ratio is exactly
`-m D(pi||nu tensor nu)/2+o(m)`, even before safe symmetry averaging.
The absolute-marginal constraint makes the reference-log term equal to
`2H(nu)`. The permanent square root has exponent `Phi=-F/2`.

At fixed r, every descendant alphabet is finite. The child bases and
fresh node permutations are independent conditional on the input word.
The total terminal loss is `O_r(sqrt(m))`, and the type errors are
`o_r(m)`. Initial rounding `k=floor(pm)` only changes a convergent finite
alphabet and compact type constraints. Therefore the one-row exponent is
bounded above by

`p log 2+(B^r Phi_t)(nu_p)`.

This is an upper bound on the supremum over every pair policy, not a
selected-policy lower witness. The basis orientation is also correct:
choose `H_i=sqrt(m) U_i^T`, so the restricted row spectrum is
`U_i(1_T x)/sqrt(k/m)` after normalization.

Choose `0<eta<a'/t` and `gamma=1-sqrt(k/m)+eta`. The failure exponent
divided by `m^2` is at most `-a'+t eta+o(1)<0`. Thus an actual full
sign matrix exists for all sufficiently large allowed orders with

`max_x |x^T K x| <= m^2 k (sqrt(k/m)-eta)`.

Deleting its diagonal pays at most N in the quadratic form, or N/2 in
the problem's Hamiltonian Q. Since `N=mk`, its normalized cap is

`Q(A)/N^(3/2) <= 1/2-eta/(2 sqrt(k/m))+1/(2 sqrt(N))`.

## 6. Dense terminal orders and the all-order endpoint

For a prime q congruent to 3 modulo 4, the quadratic-character matrix C
has `C^T=-C`, `C1=0`, and `CC^T=qI-J`. The bordered skew matrix R satisfies
`RR^T=qI`; hence `H=I+R` is a sign Hadamard of order `q+1`. This uses no
symmetry or Boolean extremal property of H.

The fixed-modulus prime number theorem applies with modulus 4 and residue
3, so consecutive such primes have ratio tending to one. I freshly read
Theorem 2.1 and its hypotheses in
[Soprounov's primary exposition](https://academic.csuohio.edu/soprunov-ivan/wp-content/uploads/sites/93/2023/02/primes.pdf).
Only qualitative fixed-progression relative density is needed; no short-
interval rate or uniformity in a growing modulus is imported.

For fixed depth r, set `m_j=2^r(q_j+1)` and
`N_j=m_j floor(p m_j)`. Then `N_(j+1)/N_j->1`. Principal restriction of
a hollow signing cannot increase Q: average the outside spins with mean
zero to recover each fixed restricted Hamiltonian. Restricting from the
least `N_j>=n` therefore extends the bound to every sufficiently large n.

Finally let eta approach `a'/t`, then let `a'` approach a. Depth is fixed
before each matrix-order limit, and may vary only between these outer
margin choices. This proves

`limsup M_n/n^(3/2) <= 1/2-a/[2t sqrt(p)]`
`                       =1/2-a/[8 sqrt(31/32)]`.

There is no growing-depth terminal estimate, minimizer transfer, or
unproved dense-order assumption hidden in this conclusion. The matching
lower-limit equality needed for the original convergence question remains
unproved.

## 7. Director synthesis and outward decimal

The complete final director synthesis
`continued_director_strict_all_order_upper_2026_09_06.md` was then read.
Its statement, strict finite-depth margin, order of limits, and explicit
convergence disclaimer agree with this reconstruction. The displayed
decimal was independently checked using the exact square-root upper
endpoint: the resulting rational upper endpoint is

`80459630021641337701/161102201102367360000`

and is strictly less than `499432220485404/10^15`. Thus the convenient
outward decimal `0.499432220485404` is safe; the longer approximate
decimal should not be substituted as an uncertified outward rounding.
