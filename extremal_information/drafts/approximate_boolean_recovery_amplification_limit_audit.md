# Independent audit: approximate Boolean recovery amplification limit

**Verdict: PASS.**  The directed-excess lemma has the correct orientation,
the matrix and tensor certificates have the correct normalizations, the
slow logarithmic drift is a valid square-summable-step/nonconvergent
falsifier, and the perfect-matching construction gives a genuine dense sign
hierarchy with perturbation norm exactly `2`.  No theorem-level repair is
required.  Two optional precision/coverage improvements identified during
the audit were subsequently applied: the reverse-rate wording was sharpened,
and the verifier now samples both drift subsequences.

## 1. Lemma AR.1: orientation, compactness, and the absent reverse rate

The direction in (AR.4) is correct.  Given `x_r in K_r`, compactness of
`K_(r+1)` makes the distance from `x_r` to `K_(r+1)` attainable, so one may
choose a forward chain with

```math
d(x_s,x_{s+1})\leq\epsilon_s.
```

It is Cauchy, and if its limit is `z`, then

```math
d(x_r,z)\leq\sum_{s\geq r}\epsilon_s.
```

Since this `z` belongs to the proposed limit set, taking the supremum over
`x_r` proves exactly

```math
e(K_r,K_\infty)\leq\sum_{s\geq r}\epsilon_s.
```

This is the old-to-limit direction, not the reverse direction.

The proof of reverse convergence is also sound.  The set `E` of forward-chain
limits is nonempty and its closure is compact because `X` is compact.  A
finite net for `closure(E)` may be chosen with centers in `E`.  Every such
center has one chain point in every sufficiently late `K_r`, at distance at
most the same tail.  Taking the largest of the finitely many starting levels
therefore proves `e(K_infinity,K_r)->0`.  No convexity is needed.

The claimed lack of a reverse defect-only rate is correctly witnessed by
(AR.34)--(AR.35).  The appended Boolean pair `(1,1)` annihilates every block

```math
b_j\begin{pmatrix}-1&1\\1&-1\end{pmatrix},
```

so every forward defect is exactly zero.  On self-responses, each block
contributes either `0` or `-4b_j`, hence convexification gives

```math
K_r^{(1)}=\left[-4\sum_{j\leq r}b_j,0\right],
\qquad
K_\infty^{(1)}=[-4B,0].
```

Thus the reverse excess and the Hausdorff distance are
`4 sum_(j>r)b_j`, even though the displayed forward-defect tail is identically
zero.  The sequence `(b_j)` can make this convergence arbitrarily slow.

The only suggested wording repair is semantic, not mathematical.  A compact
ambient space always supplies the trivial bound `diam(X)`, so the sentence
"There is no bound" is most literally stated as:

> There is no universal bound on `e(K_infinity,K_r)` that tends to zero as a
> function only of `sum_(s>=r) epsilon_s`; the tail may even vanish
> identically while the reverse excess is positive.

## 2. Theorem AR.3: compressed kernels and operator factors

For `C_r=A_r/a_r` and a replication matrix
`T_r : R^(n_r) -> R^(n_(r+1))`, the lifted pair error is exactly

```math
\kappa_{r+1}(T_rx,T_ry)-\kappa_r(x,y)
=x^T(T_r^TC_{r+1}T_r-C_r)y=x^TD_r^\uparrow y.
```

Maximizing over the two independent Boolean vectors proves (AR.17), with no
factor of two.  The notation in (AR.16) agrees with the induced norm because

```math
\max_{\|y\|_\infty\leq1}\|Dy\|_1
=\max_{x,y\in\{\pm1\}^n}|x^TDy|.
```

The two advertised certificates are correctly normalized:

```math
|x^TDy|\leq\sum_{ij}|D_{ij}|,
\qquad
|x^TDy|\leq(\sqrt n)\|D\|_{2\to2}(\sqrt n)
=n\|D\|_{2\to2}.
```

For the backward map, `D_r^down` acts on `R^(n_(r+1))`, so the analogous
operator factor is `n_(r+1)`, exactly as the draft states.  Requiring one
nonzero sign in each row of `T_r` (respectively `S_r`) is sufficient to send
every Boolean input to a Boolean output; neither injectivity nor balanced
column multiplicities are silently used.

## 3. Corollary AR.4: tensor normalization and outer support

Put `n=n_r`, `h=h_r`, `n_(r+1)=hn`, and
`T=I_n tensor u_r`.  Directly,

```math
T^T(A_r\otimes H_r)T
=A_r(u_r^TH_ru_r)=\rho_r h^{3/2}A_r.
```

Dividing by `(hn)^(3/2)` gives `rho_r A_r/n^(3/2)`, so subtraction of
`C_r` yields precisely (AR.27).  For the first summand,

```math
\left\|{A_r\over n^{3/2}}\right\|_{\infty\to1}
\leq n{\|A_r\|_{2\to2}\over n^{3/2}}=M_r.
```

For the perturbation, `Tx` and `Ty` are Boolean vectors of length `hn`, so

```math
{ |(Tx)^TE_r(Ty)|\over(hn)^{3/2}}
\leq{(\sqrt{hn})\|E_r\|_{2\to2}(\sqrt{hn})
       \over(hn)^{3/2}}
={\|E_r\|_{2\to2}\over\sqrt{hn}}=e_r.
```

This proves (AR.28) with no lost factor of `h`, `n`, or `2`.  The intrinsic
kernel bound follows from the same calculation:

```math
{|x^TA_ry|\over n_r^{3/2}}\leq M_r.
```

The vector-residual estimate is also exact as an inequality:

```math
|\rho_r-1|
={|u_r^T(H_ru_r-\sqrt{h_r}u_r)|\over h_r^{3/2}}
\leq{\|H_ru_r-\sqrt{h_r}u_r\|_2\over h_r}.
```

Finally,

```math
M_{r+1}
\leq M_r{\|H_r\|_{2\to2}\over\sqrt{h_r}}+e_r
\leq(1+\sigma_r)M_r+e_r,
```

so (AR.29) and the stated finite-product bound are valid.

For the outer query, writing a Boolean vector as blocks
`x=(x_1,...,x_d)` gives

```math
x^T(B\otimes A_r)x
=\sum_iB_{ii}x_i^TA_rx_i
+2\sum_{i<j}B_{ij}x_i^TA_rx_j.
```

Consequently `theta_B` must have diagonal entries `B_ii` and off-diagonal
entries `2B_ij`, as in the draft, and

```math
{1\over2(dn_r)^{3/2}}\max_x x^T(B\otimes A_r)x
={1\over2d^{3/2}}h_{K_r^{(d)}}(\theta_B).
```

The outer factor `1/2` and the factor `d^(-3/2)` are therefore correct.
Convexification does not change either the signed maximum or the maximum
absolute linear response.

## 4. The logarithmic-drift falsifier

For (AR.30), the normalized kernel matrix is

```math
C_r={A_r\over n_r^{3/2}}={c_r\over n_r}I_{n_r}.
```

With the duplication matrix `T=I_(n_r) tensor (1,1)^T`,

```math
T^TC_{r+1}T-C_r
={c_{r+1}-c_r\over n_r}I_{n_r}.
```

Its Boolean bilinear norm is exactly `|c_(r+1)-c_r|`, proving (AR.31) for
the all-pairs distortion, not merely for self-responses.  The step errors
satisfy

```math
\alpha_r
\leq\log{r+3\over r+2}
=\log\left(1+{1\over r+2}\right)
\leq{1\over r+2},
```

and hence `sum_r alpha_r^2<infinity`.

Every Boolean `x` has

```math
{x^TA_rx\over n_r^{3/2}}=c_r,
```

so `K_r^(1)={c_r}`.  To verify genuine nonconvergence, let `m_k` be the
nearest integer to `exp(pi/2+2pi k)` and put `r_k=m_k-2`.  Relative rounding
error tends to zero, hence `log(r_k+2)-(pi/2+2pi k)->0` and
`c_(r_k)->1`.  Repeating this with `exp(3pi/2+2pi k)` gives a subsequence
tending to `-1`.  Thus the doubling lift has vanishing and square-summable
per-step distortion while its scalar carrier does not converge.  Since
finite total variation would force `(c_r)` to converge, the step errors are
necessarily nonsummable, as claimed by the subsection title.

## 5. Corollary AR.6: perfect-matching flips and hollowing

A Hadamard matrix of order `h>1` has even order: two distinct orthogonal
sign rows agree and disagree equally often.  Thus `N_(r+1)=hN_r` is even
and admits a fixed-point-free perfect matching at every step, regardless of
the parity of `N_0`.

Reversing an off-diagonal sign changes it by `+2` or `-2`.  Because the
changed edges form a perfect matching, a simultaneous row/column
permutation puts the symmetric perturbation into blocks

```math
\begin{pmatrix}0&\pm2\\\pm2&0\end{pmatrix}.
```

Each block has eigenvalues `+2` and `-2`; their direct sum therefore has

```math
\|E_r\|_{2\to2}=2,
\qquad e_r={2\over\sqrt{N_{r+1}}}.
```

Since `N_r=N_0h^r` and `h>1`, these errors are summable.  The Boolean
eigenvector assumption gives
`u^THu=h^(3/2)`, hence `rho_r=1`; Hadamard orthogonality gives
`||H||_(2->2)=sqrt(h)`, hence `sigma_r=0`.  Recurrence (AR.29) supplies the
uniform `M_r` bound.  All hypotheses of AR.4 are therefore met.  The
construction changes exactly `N_(r+1)/2`, or `Theta(N_r)`, undirected edges
at each level and retains symmetric full sign entries.

For `A_r=C_r^circ`, Booleanity makes the deleted diagonal a constant:

```math
x^TC_rx-x^TA_rx=\sum_i(C_r)_{ii}x_i^2=\operatorname{tr}(C_r).
```

Using the project convention
`Q(A)=(1/2)max_x|x^TAx|` for a hollow symmetric signing,

```math
\left|{Q(A_r)\over N_r^{3/2}}-p_r^{abs}(C_r)\right|
\leq{|\operatorname{tr}(C_r)|\over2N_r^{3/2}}
\leq{1\over2\sqrt{N_r}}.
```

The correction vanishes because `h>1`, so the hollow normalized energies
converge.  The factors `2` and `1/2` in (AR.29e)--(AR.29h) are correct.

## 6. Verifier audit and rerun

I reran

```bash
python3 \
  extremal_information/experiments/verify_approximate_boolean_recovery_amplification_limit.py
```

and obtained

```text
approximate Boolean recovery checks passed: 66693
```

The verifier correctly checks the exact regular-Hadamard compression, the
compressed perturbation norm against exhaustive all-pairs distortion, the
`n||D||_(2->2)` and entrywise bounds, an exact signed left inverse and its
reverse distortion, both finite falsifiers, matching-flip symmetry/signs
and operator norm `2`, the quasi-monotone scalar inequalities, and the
hollowing correction.

Its scope is appropriately finite.  The two optional audit refinements were
applied after the first run:

- The draft now says precisely that no reverse-excess rate tending to zero
  can depend only on the forward-defect tail.
- The drift block now samples rounded positive and negative logarithmic-phase
  subsequences.  Its self-response comment now accurately identifies the
  tested all-ones witness; the all-Boolean statement follows algebraically
  from `x^Tx=n`.

These refinements do not change the PASS verdict.
