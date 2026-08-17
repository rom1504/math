# Independent audit: PC.3 sparse-flip coherence obstruction

**Verdict: PASS, with exact scope.**  Proposition DS.3 and Theorem DS.4 in
`pc3_diagonal_switching_coherence.md` rigorously produce symmetric exact-sign
children whose active-product Rayleigh deficits vanish while one prescribed
same-field selector defect stays constant.  This is not a separation of the
full Boolean trust optima.

## 1. Sparse-flip constants

Let `H^2=NI`, `r=sqrt(N)`, and let both `x` and every `z in Z` be Boolean
positive `r`-eigenvectors.  Eligible unordered edges obey
`H_ij x_i x_j=1`, and are flipped independently with probability
`q=kappa/r`.

For a Boolean test `y`, flipping one eligible unordered edge changes
`y^THy` by `-4H_ijy_iy_j`.  For `z`,

```math
\begin{aligned}
4\sum_{i<j}1_{H_{ij}x_ix_j=1}H_{ij}z_iz_j
 &=2\sum_{i<j}H_{ij}z_iz_j
   +2\sum_{i<j}x_ix_jz_iz_j\\
 &=rN+(x^Tz)^2-N.
\end{aligned}
```

This confirms (DS.30), including its factors of two.  Division by
`rN=N^(3/2)` and substitution `q=kappa/r` leave
`kappa(x^Tz/N)^2+o(1)`.

For `y=x`, the eligible-edge count is

```math
{1\over2}\left({N(N-1)\over2}+{rN\over2}\right)
={N(N-1)+rN\over4},
```

so the normalized expected loss is `kappa+o(1)`.  Scalar Bernstein has
variance `O(qN^2)=O(N^(3/2))`; the assumption
`log|Z|=o(N^(3/2))` indeed permits a uniform `o(N^(3/2))` fluctuation.

For the operator norm, off the diagonal

```math
\mathbb EH'=(1-q)H-qxx^T.
```

The diagonal correction has norm at most `2q`.  Since `Hx=rx`, the rank-one
term commutes with `H`; for fixed `0<kappa<1`, all eigenvalues of the mean
have modulus at most `r+o(r)`.  Matrix Bernstein contributes only
`O(N^(1/4)sqrt(log N)+log N)=o(r)`.  Finally the full sign matrix has
Frobenius norm `N`, forcing `||H'||_op>=sqrt(N)`.  Thus the use of the actual
operator norm in (DS.24)--(DS.26) is valid rather than an unproved choice of
normalization.

Deleting the diagonal preserves every Boolean energy because the trace is
unchanged and zero; it changes operator norm by at most one.  Hence the
hollow matrix `A'=H'-diag(H')`, with explicit roof `||A'||_op<=r'+1`, has
the same asymptotic conclusions.

## 2. Diffuse-selector recurrence

The PC.3 row law for the two relative generators is

```math
(X,Y)=(1,1),(1,-1),(-1,1)
\quad\hbox{with probabilities}\quad1/4,1/2,1/4.
```

For endpoint (DS.35), the score increment is
`A_t=(-1)^tX_t+Y_t`.  It is even and lies in `[-2,2]`, so
`1+sum_t A_t` is always odd: there is no hidden tie convention.  Its mean
is `(-1)^t/2`; the full mean is bounded.  Both parity laws have variance
bounded below (`3/4` and `11/4`, respectively), and all third moments are
uniformly bounded.  Berry--Esseen therefore applies uniformly after
deleting any subset of factors.

An active product, after cancelling the common base pole, is exactly a raw
product of local characters in `{1,X,Y,XY}`.  Their means are
`1,1/2,0,-1/2`; hence every nonconstant factor has mean modulus at most
one half.  Conditioning on such a factor differs from replacing it by its
mean only when the remaining score is within two of the threshold.  The
localized error is controlled by the uniform `O(m^(-1/2))`
anti-concentration estimate.

After `ell` eliminations, the accumulated coefficient is at most `2^-ell`.
Splitting the resulting geometric series at `ell=j/2` gives
`O(j^(-1/2))`.  The terminal term is also valid:

* if at least `j/2` factors were removed, `2^-k` is negligible;
* otherwise the remaining variance is `Theta(j)`, its mean is at most
  `(k+1)/2`, and Berry--Esseen gives
  `O((k+1)/sqrt(j))`, which remains `O(j^(-1/2))` after multiplication by
  `2^-k`.

This confirms uniformity over all `4^j` products.  There is no missing union
bound and no assumption that the raw characters are orthogonal under the
three-atom law.

## 3. Scope boundary

The result proves a coherent **certificate** obstruction:

```math
\max_{z\in Z_j}d_z=o(1)
\quad\hbox{but}\quad d_{x_\epsilon}=\kappa+o(1).
```

It does not say the selector is the optimal Boolean witness, nor that the
complete trust response loses `kappa`.  Another Boolean spin may repair the
full optimum.  It also does not contradict DS.2: targeted edge flips are
not diagonal conjugations of the PC.3 child.
