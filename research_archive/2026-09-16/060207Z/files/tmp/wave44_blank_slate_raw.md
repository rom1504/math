# Wave 44 blank-slate abstraction audit — generation phase

Generation constraint: use only

~~~math
M_n=\min_{a_{ij}\in\{\pm1\}}\max_{x_i\in\{\pm1\}}
\left|\sum_{i<j}a_{ij}x_ix_j\right|
~~~

and the question whether M_n/n^(3/2) converges. The comparisons with the
research ledger are intentionally deferred to a separate phase below.

These are main-agent research hypotheses, not user directives.

## Candidate A: discrepancy-preserving decorated lifts

Map an order-n signing to an order-N signing by replacing each old vertex
with an almost equal fiber and decorating every complete bipartite fiber pair
by a balanced sign matrix. A Boolean state on the lift is then compressed
not just to one fiber magnetization, but to its correlations with the chosen
decorations. The desired decorations would make all nonconstant fiber modes
nearly orthogonal while the constant modes reproduce the old quadratic form.

The concrete theorem is the following uniform lift inequality. For every
order n, every minimizing signing A_n, and every N>=n, construct a signing
B_N such that

~~~math
\max_{y\in\{\pm1\}^N}|H_{B_N}(y)|
\le \left(\frac Nn\right)^{3/2}M_n+\varepsilon_nN^{3/2},
\qquad \varepsilon_n\longrightarrow0.
\tag{A}
~~~

Taking n along a subsequence realizing the liminf and then allowing every
larger N in (A) gives limsup<=liminf, hence convergence. The mathematical
mechanism sought is an orthogonal-array or discrepancy-corrected composition
which controls every Boolean choice inside every fiber, not only constant
fiber states.

## Candidate B: a convergent sum-of-squares surrogate hierarchy

For a signing A, let S_(n,r)(A) be the degree-2r sum-of-squares upper
bound for

~~~math
\max_{x_i^2=1}|H_A(x)|.
~~~

Put V_(n,r)=min_A S_(n,r)(A). This maps the original problem to a nested
sequence of finite semidefinite moment problems, with V_(n,r)>=M_n and
equality when the hierarchy reaches full degree.

A concrete two-part theorem would prove convergence:

~~~math
\text{for every fixed }r,\quad
\lim_{n\to\infty}\frac{V_{n,r}}{n^{3/2}}=v_r,
\tag{B1}
~~~

and

~~~math
\lim_{r\to\infty}\ \limsup_{n\to\infty}
\frac{V_{n,r}-M_n}{n^{3/2}}=0.
\tag{B2}
~~~

Indeed, for large fixed r, (B2) places the original normalized sequence
uniformly within o_r(1) of the convergent sequence in (B1), making it
Cauchy. Candidate machinery is exchangeability/compactness of fixed-degree
pseudomoment arrays for (B1), plus a signing-uniform integrality-gap theorem
for (B2).

## Candidate C: asymptotically exact coset-weight linear programming

Encode an edge signing as a binary word. Adding a cut or its complement
forms the augmented cut code

~~~math
\mathcal C_n=\{(t+z_i+z_j)_{i<j}:t,z_i\in\mathbb F_2\}.
~~~

For an edge word u, the closest-codeword distance is determined by its
coset weight enumerator. Thus the original minimum-maximum problem is the
deficit of the largest coset minimum weight from half the block length.

Let LP_(n,r) be the order-r truncated MacWilliams/Krawtchouk moment
program over feasible coset weight distributions of C_n, optimized for
the largest possible coset minimum weight, and let

~~~math
D_{n,r}=\frac12\binom n2-LP_{n,r}.
~~~

The concrete theorem is

~~~math
\lim_{n\to\infty}\frac{D_{n,r}}{n^{3/2}}=d_r
\quad\text{for every fixed }r,
\tag{C1}
~~~

together with asymptotic tightness

~~~math
\lim_{r\to\infty}\ \limsup_{n\to\infty}
\frac{|D_{n,r}-M_n/2|}{n^{3/2}}=0.
\tag{C2}
~~~

Krawtchouk-to-Hermite asymptotics are a candidate compact limiting theory for
(C1); (C2) is the exact nontrivial requirement that bounded-order coset
moments determine the covering-radius deficit at the n^(3/2) scale.
Together they force convergence of M_n/n^(3/2).

## Ledger-comparison phase

Performed only after Candidates A--C above were fixed.

### Candidate A comparison

This is exactly the uniform-amplification objective in ledger §6.1, with the
same nonconstant fiber-state issue already isolated in §§3.4, 3.8 and
10.31--10.32. The ledger goes substantially beyond the raw proposal:
multi-channel Fourier lifts relax to the spectral ceiling, exact dependent
lifts require all tuple levels, compulsory orthogonal-fiber variance survives
static filters, and centered width cannot contract. Candidate A supplies no
new decoration theorem that evades those walls. Reject it as a renamed
amplification route, while retaining inequality (A) as an exact sufficient
benchmark.

### Candidate B comparison

Ledger §§3.14--3.16 and 3.18 already show that fixed-degree cut-cone/SOS,
fixed moments, and fixed-replica profiles cannot see a zero-entropy planted
extremizer; a rank-one obstruction survives every local test of radius
o(sqrt(n)). Sections 10.31 and 10.54 further show that compact fixed profiles
do not give all-order realization and that SDP dual slack can push near-ground
mass into thin curvature kernels. Thus (B2), not (B1), contains essentially
all the missing mathematics. It would require a genuinely growing-order,
signing-uniform integrality-gap theorem, at which point the hierarchy again
encodes the remote Boolean extrema. Candidate B is testable but not a new
mechanism at present; keep it dormant rather than displace sharper routes.

### Candidate C comparison

Ledger §3.15 gives the exact signed Eulerian/Krawtchouk enumerator and proves
that the full signed coset weight enumerator is invertibly the original
energy histogram. Fixed dual degree sees only fixed moments and misses sparse
resonant maxima; growing degree returns the original problem. Section 10.50
adds that even the complete one-point Johnson/Delsarte algebra is blind to
the distinction between quadratic- and 3/2-scale regular signings. Hence
(C2) is precisely the unavailable resonance-sensitive high-degree theorem,
not a consequence of MacWilliams compactness. Reject Candidate C as a coding
reformulation without a stronger separating statistic.

### Audit judgment

All three candidates give exact convergence implications, but the ledger
comparison finds no genuinely stronger mechanism at the required scale.
Candidate A is the known amplification problem, Candidate C is the known
signed-histogram dual, and Candidate B hides the original remote-extremum
problem in a uniform hierarchy gap. None changes the leading strategy.
