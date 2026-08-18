# Independent audit: constant-mass row-magnitude fibres

**Frozen source:**
`extremal_information/drafts/conference_row_magnitude_fibres.md`

**SHA-256:**
`6e896153ebe6eaec1ecaea67d21731aa42d74e0efbf90079887c09b18ec60150`

**Verdict:** **PASS.**  The family mass, conditional covariance, layer
coupling, dependent-row operator tail, high-temperature normalization,
nuclear/Frobenius conversion, probability and `L^1` transfer, and stated
scope are correct.  No mathematical repair is required.

## Exact mass and covariance

After gauging `v_r` to the all-ones vector, membership in `E_r` depends only
on the Hamming layer and is invariant under both coordinate permutations and
global row negation.  Consequently a conditioned row is exchangeable and
centered.  The product family has exactly

```math
|E_r|^r=(2^rp_r)^r=2^{r^2}p_r^r,
```

so its log-density in the full bridge cube is exactly `r log p_r`.  Under
the standing lower bound `p_r >= p_0`, this is an `O(r)` entropy cost.

Writing `S=sum_j R_j`, exchangeability gives a common off-diagonal
covariance `rho_r`.  Since

```math
E[S^2\mid E_r]=r+r(r-1)\rho_r=r\alpha_r,
```

one has `rho_r=(alpha_r-1)/(r-1)`.  The covariance eigenvalues are therefore
exactly `alpha_r` along the all-ones direction and
`(r-alpha_r)/(r-1)` on its orthogonal complement.  Moreover

```math
\alpha_r={1\over r}E[S^2\mid E_r]
\le {1\over rp_r}E S^2={1\over p_r}\le {1\over p_0}.
```

Thus the finite-rank covariance conclusion is correct, as is the warning
that this conclusion alone says nothing sufficient about Boolean pressure.

## Layer coupling and Frobenius cost

Let `K` have the unconditional binomial layer law and let `K'` have that law
conditioned on the allowed absolute layers.  Even the independent coupling
obeys

```math
E|K-K'|
\le {1\over2}E|2K-r|+{1\over2}E|2K'-r|
\le {1+p_0^{-1}\over2}\sqrt r.
```

For fixed `k,k'`, choosing the smaller plus set uniformly and then extending
it uniformly to the larger layer gives uniform marginals on both Hamming
layers and distance exactly `|k-k'|`.  This verifies RM.1 without an implicit
loss of uniformity.

With independent copies across rows, if `D=R_r-W_r`, then every changed sign
contributes `4` to `||D||_F^2`, so

```math
E\|D\|_F^2=4E d_H(R_r,W_r)=O_{p_0}(r^{3/2}),
\qquad
E\|D\|_F=O_{p_0}(r^{3/4}).
```

The second estimate follows by Cauchy--Schwarz/Jensen and is indeed `o(r)`,
which is the scale needed after the nuclear/Frobenius conversion.

## Uniform dependent-row operator tail

The absolute-magnitude event is centrally symmetric, hence the conditioned
row has mean zero.  For every deterministic unit vector `z`, conditioning
the ordinary Rademacher Hoeffding bound costs at most `p_0^{-1}`:

```math
P(|\langle R,z\rangle|>u\mid E_r)
\le {2\over p_0}e^{-u^2/2}.
```

This is a uniform `psi_2` bound depending only on `p_0`.  Therefore
`<R_i,z>^2` has uniformly bounded `psi_1` norm.  Bernstein's inequality for
the independent rows, at a sufficiently large constant multiple of `r`,
has exponent whose constant can be made larger than `log 9`; a `1/4`-net
of size at most `9^r` then gives

```math
P(\|R_r\|_{op}>L_{p_0}\sqrt r)\le 2e^{-c_{p_0}r}.
```

No independence among coordinates inside a row is used.  The constants are
uniform over the magnitude set and the gauge vector.

As a separate robustness check, one can also obtain an adequate conditioned
operator tail directly from the full bridge measure: choose a sufficiently
large `L=L(p_0)` so that the rectangular Rademacher norm tail is at most
`exp[-(log(1/p_0)+c)r]`, and divide by the fibre probability
`p_r^r >= p_0^r`.  This gives an exponential in-fibre tail independently of
RM.2 and confirms that no coupling-dependent norm claim is hidden in RM.3.

## High-temperature threshold and stability constants

For a conference signing, `||A_r||_op=sqrt(r-1)`.  Block triangle inequality
therefore gives, for either orientation,

```math
\left\|{\beta\over\sqrt{2r}}
\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix}\right\|_{op}
\le {\beta\over\sqrt2}
\left(\sqrt{1-1/r}+{\|B\|_{op}\over\sqrt r}\right).
```

Thus RM.16 places both coupled parents in the same closed operator ball of
radius `kappa<1/2`.  The operator ball is convex, so the whole interpolation
segment required by the pressure-stability proof remains in that ball.

For `D=R_r-W_r`, the off-diagonal symmetric dilation has singular values of
`D` twice and hence

```math
\left\|\begin{pmatrix}0&D\\D^T&0\end{pmatrix}\right\|_*
=2\|D\|_*.
```

With `t=beta/sqrt(2r)` and `||D||_*<=sqrt(r)||D||_F`, the archived stability
bound becomes exactly

```math
{K_\kappa\over2}(2t\|D\|_*)
\le {K_\kappa\beta\over\sqrt2}\|D\|_F.
```

There is no missing factor of two or power of `r`.  Taking expectations on
the joint regular event produces `O(r^(3/4))=o(r)` pressure error.

## Exceptional event, probability, and `L^1`

The two marginal operator failures have probabilities `e^{-Omega(r)}`;
their dependence under the coupling is irrelevant by a union bound.  Every
sign parent has

```math
0\le f_{\epsilon,r}(B)\le C_\beta r^{3/2},
```

so the exceptional contribution is `O(r^(3/2)e^{-cr})=o(r)`.  Hence

```math
E|f(R_r)-f(W_r)|=o(r).
```

The archived uniform-conference pressure theorem then transfers convergence
in probability.  It also transfers `L^1`: the uniform bridge input has the
needed mean/probability conclusion (and its normalized pressure is
nonnegative), while the displayed expected coupling error is `o(r)`.

There is an independent conditional proof of the probability statement:
on a fixed operator ball the archived conference lower/upper tail is
`e^{-Omega(r^2)}`; dividing by the row-product fibre mass, at worst
`p_0^r`, leaves `e^{-Omega(r^2)}`.  Together with the conditioned operator
tail, this yields the same conclusion and shows that the result is genuinely
a quenched pressure theorem rather than only a covariance or norm statement.

## Scope and one non-substantive wording point

The proved interval is nonempty but may be much smaller than the full
campaign interval because `L_{p_0}` is not a sharp spectral-edge constant.
The theorem makes no claim when `p_r` tends to zero.  Both limitations are
stated accurately.  The result applies to central bands, fixed normalized
tails, finite unions of magnitude layers, and parity-compatible layer sets
whenever their total row mass stays bounded below.

The sentence saying that every “fixed nontrivial” constraint has *precisely*
speed-`r` entropy loss should be read with the usual additional meaning that
`p_r` stays bounded away from both zero and one.  Under only RM.2, the exact
statement is the stronger and fully correct formula RM.4, with entropy loss
`r log(1/p_r)`, which can be `o(r)` when `p_r -> 1`.  This prose point is not
used by RM.1--RM.3 and does not require a mathematical repair.

## Corrections

None required.
