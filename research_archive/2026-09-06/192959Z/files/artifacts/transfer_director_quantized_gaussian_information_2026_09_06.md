# Quantized Gaussian information through random principal restrictions

2026-09-06. **Proved; independently audited.** This develops the
OU/resampling mechanism found by the seed researcher into a weighted,
non-flat information inequality. It is not a convergence theorem. All
norms below are Euclidean operator norms and all logarithms are natural.

## 1. A finite-dimensional, channel-sensitive inequality

Let `B` be any real symmetric hollow `N` by `N` matrix, `||B||<=C`.
Choose real `s` and `r` with `|s|C<r<1`, and put `t=s/r`. Under
`G~N(0,I+sB)`, observe `Y_i=h_i(G_i)`, where each `h_i` is a measurable
finite-valued quantizer with nonzero cell probabilities under the
standard normal law `gamma`. Let `mu` be the product of these marginals.

On `L2(gamma)`, let `P_i` be conditional expectation onto the quantizer,
`Pi` projection onto constants, and `T_a` the Gaussian OU operator
`T_a f(x)=E f(ax+sqrt(1-a^2)Z)`. Define

```math
\kappa_i=\|T_{\sqrt r}(P_i-\Pi)T_{\sqrt r}\|_{2\to2}\in[0,1].
```

For independent selectors `Z_i~Ber(kappa_i)` and their support `S`,

```math
 D(\mathcal L(Y)\Vert\mu)
 \le \log(1+\chi^2(\mathcal L(Y),\mu))
 \le \log\mathbb E_S\det(I-t^2 B_S^2)^{-1/2}.       \tag{1}
```

The empty determinant is one. The finite matrix determining `kappa_i`
has dimension at most the number of nonempty quantizer cells minus one;
it does not encode an `N`-spin response table.

### Proof

Let `L_u` be the likelihood of `N(0,I+uB)` relative to `gamma^N`.
Gaussian convolution, or characteristic functions, gives the exact
identity `L_s=T_sqrt(r)^{tensor N} L_t`. The assumption `|t|C<1`
implies both positive covariance and `L_t in L2(gamma^N)`.

The output likelihood is conditional expectation of `L_s`; hence its
squared `L2` norm is the quadratic form of
`tensor_i [T_sqrt(r) P_i T_sqrt(r)]` at `L_t`. Each factor is positive
and is bounded in positive-semidefinite operator order by
`Pi+kappa_i(I-Pi)`. Tensoring preserves this order: telescope the
difference of the two tensor products into positive tensor products.
Expand each bounding factor as `(1-kappa_i)Pi+kappa_i I`.
The resulting quadratic form is the average squared norm of the
Gaussian marginal likelihood on `S`. Direct Gaussian integration gives

```math
 \|\mathbb E[L_t\mid G_S]\|_2^2
 =\det(I-t^2B_S^2)^{-1/2}.
```

Finally, `D(p||mu)<=log E_p[p/mu]` is Jensen's inequality. This proves
(1), without a Taylor remainder, Gaussian universality, or assumptions
on eigenvectors.

## 2. Weighted bound requiring no entrywise flatness

Write `w_ij=B_ij^2`,
`D_*=max_i sum_(j!=i) w_ij`, and

```math
 c=\frac{t^2}{2(1-t^2 C^2)}.
```

If `D_*=0`, the outputs are independent and their divergence is zero.
Otherwise (1) implies the entirely explicit estimate

```math
 D(\mathcal L(Y)\Vert\mu)
 \le \sum_{i<j}\frac{w_{ij}}{D_*}
       \log\{1+\kappa_i\kappa_j(e^{2cD_*}-1)\}.       \tag{2}
```

In particular, if all `kappa_i<=kappa`,

```math
 D(\mathcal L(Y)\Vert\mu)
 \le \frac{\operatorname{tr}B^2}{2D_*}
       \kappa^2(e^{2cD_*}-1).                        \tag{3}
```

### Proof of the weighted step

The scalar inequality `-log(1-u)<=u/(1-t^2 C^2)` for
`0<=u<=t^2 C^2` bounds the logarithm of the determinant in (1) by
`c tr B_S^2=2c sum_(i<j) w_ij Z_i Z_j`.
Set `theta_ij=w_ij/D_*` on positive-weight edges. The sum of these
weights at each vertex is at most one. For independent coordinates,
the fractional Holder inequality is

```math
 \mathbb E\prod_e f_e(Z_e)^{\theta_e}
 \le\prod_e(\mathbb E f_e(Z_e))^{\theta_e}.          \tag{4}
```

Here is a self-contained proof for the finite setting needed here:
integrate one coordinate at a time, apply Holder with exponents
`1/theta_e` to the incident functions, and insert the constant function
one for any unused exponent mass. Each incident function is replaced
by its integral over that coordinate, retaining its exponent. Repeat
until every coordinate has been integrated. All functions in the
present application are bounded and positive, so no integrability
qualification is hidden. This is the classical fractional Holder
mechanism, not a newly claimed inequality.

Apply (4) with `f_ij=exp(2cD_* Z_i Z_j)`. Its expectation is
`1+kappa_i kappa_j(e^{2cD_*}-1)`. Taking logarithms proves (2), and
`log(1+u)<=u` proves (3).

## 3. Binary thresholds and the rare-event scale

For identical thresholds with rare probability `delta`, let
`z=Phi^{-1}(1-delta)`, `m=1-2delta`, `v=1-m^2`,
`a=2 phi(z)`, and `eta=a^2/v`. For example use
`h(g)=1-2*1{g>=z}`; reversing all threshold labels does not change
`kappa`. The centered quantizer subspace has rank one, so

```math
 \kappa=\eta\int_0^r
       \frac{\exp\{z^2u/(1+u)\}}{\sqrt{1-u^2}}\,du
 \le \eta r\,\frac{e^{z^2r}-1}{z^2r\sqrt{1-r^2}}.  \tag{5}
```

The quotient at `z=0` is interpreted continuously. One may derive the
identity by differentiating the bivariate threshold probability with
respect to its Gaussian correlation, or by integrating its elementary
two-variable density derivative. Thus (3) has a genuine `eta^2`
factor, with a fully specified threshold-dependent constant. This
does not contradict the independent-pair counterexample as
`delta->0` at fixed nonzero `s`: its factor in (5) diverges.

For flat `B=A/sqrt N`, a sharper estimate keeps the cardinality of
`S` and uses a binomial large-deviation bound. That is the seed
researcher's separate theorem. Inequality (2) instead allows arbitrary
sparse or uneven variance geometry and nonidentical finite quantizers.
Its constant is generally worse; it is not claimed to improve the
flat theorem's numerical pressure credit.

## 4. Immediate quantitative use and scope

For an actual hollow signing `A`, take `B=A/sqrt N` after switching
and global sign reversal so that the all-plus configuration has
energy `Q(A)`. A threshold-Gaussian trial law with mean `m` has
first-Hermite covariance `a^2 s B_ij`. The remaining Hermite terms
have absolute value at most `(v-a^2)s^2/N` on each edge. The Gibbs
variational principle therefore gives

```math
 \log\sum_x2\cosh\!\left(\frac{\beta H_A(x)}{\sqrt N}\right)
 \ge Nh(\delta)+\frac{\beta m^2 Q(A)}{\sqrt N}
  +\frac{\beta a^2s}{2}(N-1)
  -\frac{\beta(v-a^2)s^2}{2}\frac{N-1}{\sqrt N}
  -\frac{N}{2}\kappa^2(e^{2c}-1).                  \tag{6}
```

In (6) use `beta>=0`, `s>=0`, the same operator bound `C`, and
`D_*=(N-1)/N<=1`; the function `(e^{2cD}-1)/D` is increasing,
which permits the displayed slightly weaker final term. All limits
fix `delta,s,r,C,beta` first. The remainder is `O(sqrt N)`.
An available partition-function upper certificate can therefore use
the explicit credit `beta a^2s/2-kappa^2(e^{2c}-1)/2` directly.
Without such a certificate, (6) alone supplies no cross-order theorem
or cap upper bound. No assumption that actual minimizers have bounded
operator norm is made; the construction to which (6) is applied does.

## 5. Relation to earlier failed ideas

Raw Gaussian KL paid for radial information lost by quantization.
Here smoothing first allows the quantizer projection to be dominated
by coordinate resampling, so two surviving endpoints pay two factors
of `kappa`. No Gaussian equality for the discrete entropy is asserted.
The proof retains the entire likelihood only in its analytic derivation;
the resulting bound needs the covariance norm, squared-entry geometry,
and one finite channel operator per coordinate. It does not reconstruct
Boolean maximization. Its immediate application is a pressure bound,
not resolution of the signing convergence problem.
