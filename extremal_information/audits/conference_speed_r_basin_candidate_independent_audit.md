# Independent audit: speed-`r` conference basin candidates

**Frozen source:**
`extremal_information/drafts/conference_speed_r_basin_candidate_audit.md`

**SHA-256:**
`fba0a166c67337776aa7c2ecac42a00a3acc3830000992957ce26ab6fb8a7f47`

**Verdict:** **PASS WITH ONE NOTATIONAL REPAIR.**  CB.1--CB.4 and their
stated scope are mathematically valid.  In (CB.24)--(CB.25), the variance
profile must be frozen before the convolution is regarded as a function of
its translation argument.  This is the evident intended proof and requires
no change to the theorem statement or bound.

## 1. Edit modulus and Hamming retraction

For a real bridge increment `Delta`, every Boltzmann exponent changes by at
most `t||Delta||_1`.  The inequality
`cosh(u+v)<=exp(|v|)cosh(u)` therefore gives (CB.7), and a sign flip costs
`2t`, giving (CB.8).

If `pi_r(U_r)` is uniform on `F_r`, coupling it to `U_r` and using
`s_r=o(r^(3/2))` gives

```math
{2ts_r\over r}
={\sqrt2\beta s_r\over r^{3/2}}=o(1).
```

Thus CB.1 follows directly from the archived convergence in probability.
No cardinality or injectivity hypothesis on the retraction is missing.

For a rank-`q_r` binary parity system, choosing pivot columns of a full-rank
parity-check matrix gives one pivot assignment for every free assignment.
The overwrite map changes at most `q_r` coordinates and every fibre output
has exactly `2^(q_r)` preimages.  Hence its pushforward is uniform.  The
row/column-product system has the single familiar dependency and rank
`2r-1`, so its dimension is `(r-1)^2`.  CB.2 is correct.

## 2. Lindeberg and convex-even comparison

Representing `cosh` by an auxiliary sign `sigma` makes differentiation in
one real bridge coordinate exact:

```math
\partial_e^3f=t^3\kappa_3(\sigma x_i y_j).
```

The observable is sign-valued, so its third cumulant has absolute value at
most `2` (indeed a smaller universal bound is possible).  The biased sign
`B_e` and `Y_e=m_e+sqrt(1-m_e^2)V_e` match their first two raw moments, and
both replacement variables are uniformly bounded.  Taylor replacement
therefore costs `O(t^3)` per coordinate and

```math
r^2t^3=O_\beta(\sqrt r).
```

The comparison direction is correct, subject to the following notational
repair.  For the already fixed mean vector `M`, set

```math
a_e=\sqrt{1-m_e^2},\qquad W_M=(a_eV_e),
\qquad G_M(Z)=\mathbb E_V f(Z+W_M).
```

Now `W_M` is frozen as `Z` varies.  Convexity of `f` makes `G_M` convex;
symmetry of `W_M` and evenness `f(-B)=f(B)` make `G_M` even.  Therefore

```math
\mathbb E f(M+W_M)=G_M(M)\ge G_M(0)=\mathbb E f(W_M).
```

Writing simply `G(M)` while allowing `a_e` to appear to vary with its
argument would not justify convexity.  Freezing the profile, as above,
repairs (CB.24)--(CB.25) without altering any conclusion.

Finally,
`1-sqrt(1-u^2)<=u^2` and (CB.7) give the asserted `O(sqrt r)` comparison
between `W_M` and a uniform sign bridge.  The argument proves the lower
bound direction in CB.3; it does not accidentally reverse the desired
inequality.

## 3. Product concentration

Under the biased product law, one sign-coordinate oscillation is at most
`2t`.  Hence

```math
\sum_e c_e^2=r^2(2t)^2=2\beta^2r.
```

McDiarmid gives
`P{f<=Ef-u}<=exp(-u^2/(beta^2r))`.  Since CB.19 implies
`Ef>=(h_beta-eta/2)r` for all large `r`, taking `u=eta r/2` yields exactly
the exponent displayed in CB.20.

## 4. Exact type shell

For `N_r=theta r^2+O(r)`, `m_r=c/sqrt(r)`, and
`k_r=floor((1+m_r)N_r/2)`, Stirling expansion gives

```math
\log {N_r\choose k_r}
=N_r\log2-{N_rm_r^2\over2}+O(\log r)
=N_r\log2-{\theta c^2\over2}r+O(\log r).
```

This proves CB.32.  Under the associated product bias, the agreement count
is binomial with its mean at `k_r+O(1)`, so the conditioning point has
probability `Theta(N_r^(-1/2))=Theta(1/r)`.  Conditional on it, the law is
uniform on the shell.  Dividing CB.20 by this polynomial probability
preserves an `exp(-c r)` bound, proving CB.4.  The finitely-many-block
extension similarly incurs only a fixed polynomial conditioning cost.

## 5. Scope checks

Pairing `y` and `-y` proves the exact product identity (CB.36).  A single
rank-one witness then gives (CB.37), including its `beta r^(3/2)/sqrt(2)`
leading term and the `2r log 2` normalization charge.

The report correctly does **not** infer a superexponential lower-pressure
tail from CB.3 or CB.4.  An `exp(-Theta(r))` exceptional subset can remain
inside a speed-`r` shell, and product bounded differences operates only at
speed `r`.  The conclusions exclude the proposed affine and independent
weak-bias mechanisms as whole basins; they do not solve the original
pressure-selected lower-deviation problem.

## Required source repair

Replace the two sentences around (CB.24)--(CB.25) by the frozen-profile
definition `G_M(Z)=E_V f(Z+W_M)` above.  No other correction is required.
