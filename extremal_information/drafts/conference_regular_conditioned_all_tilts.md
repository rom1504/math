# All fixed bridge tilts are annealed on the operator-regular conference class

**Status.** Task-local theorem draft, separate from the two frozen conference
sources.  The theorem proves a speed-`r^2` lower-pressure tail after
conditioning on the standard high-temperature operator event.  Consequently
every fixed entropy tilt, and in fact every tilt `lambda_r=o(sqrt(r))`, has
the typical conference pressure rate on that conditioned class.

This does not settle the broader cut-cap-conditioned or unconditioned laws.
It localizes any finite-tilt phase there: it must be carried by the
exponentially rare operator-irregular bridges.

## 1. Setup

Use the conference parent and pressure

```math
S_{\epsilon,B}
=\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix},
\qquad
t={\beta\over\sqrt{2r}},
\qquad
f_{\epsilon,r}(B)
=\log\overline Z_{2r}(S_{\epsilon,B},t),
\tag{RC.1}
```

with

```math
0<\beta<{\sqrt2\over6}.
\tag{RC.2}
```

Fix constants `delta>0` and `kappa<1/2` satisfying

```math
{\beta(3+\delta)\over\sqrt2}<\kappa<\frac12.
\tag{RC.3}
```

This is possible in (RC.2), after taking `delta` small.  For each
orientation define

```math
\mathcal K_{\epsilon,r}
=\left\{B:\|tS_{\epsilon,B}\|_{2\to2}\le\kappa\right\},
\tag{RC.4}
```

and let `K_r` be the corresponding event in the joint `(epsilon,B)` output
space.

The rectangular Rademacher norm tail and

```math
\|tS_{\epsilon,B}\|_{2\to2}
\le{\beta\over\sqrt2}\left(
 \sqrt{1-1/r}+{\|B\|_{2\to2}\over\sqrt r}
 \right)
\tag{RC.5}
```

give constants `c_(beta,kappa)>0` and `r_0` such that

```math
\Pr(\mathcal K_r^c)\le2e^{-c_{\beta,\kappa}r}
\qquad(r\ge r_0).
\tag{RC.6}
```

## 2. Dimension-free convex concentration on the regular event

The high-temperature covariance theorem used in Section 6 of
`artifacts/two_temperature_bridge_audit.md` gives, on the convex real set
`K_(epsilon,r)`,

```math
\|\nabla_B f_{\epsilon,r}(B)\|_F
\le L_{\beta,\kappa},
\tag{RC.7}
```

where `L_(beta,kappa)` is independent of `r`.  The function `f` is convex
in the real bridge entries.  Taking the supremum of its supporting affine
functions over `K_(epsilon,r)` gives a convex
`L_(beta,kappa)`-Lipschitz extension `g_(epsilon,r)` to the full cube that
agrees with `f` on `K_(epsilon,r)`.

Talagrand convex-Lipschitz concentration supplies constants
`C_(beta,kappa),c_(beta,kappa)>0`, independent of `r`, for which

```math
\log\mathbb E_B
 e^{-s(g_{\epsilon,r}-\mathbb Eg_{\epsilon,r})}
\le C_{\beta,\kappa}s^2
\qquad(s\ge0),
\tag{RC.8}
```

and consequently

```math
\Pr_B\{g_{\epsilon,r}-\mathbb Eg_{\epsilon,r}\le-u\}
\le e^{-c_{\beta,\kappa}u^2}
\qquad(u\ge0).
\tag{RC.9}
```

The archived proof also gives `|g|=O_(beta,kappa)(r)` on the cube, while
the elementary cap bound gives `0<=f=O_beta(r^(3/2))`.  Equation (RC.6)
therefore implies

```math
|\mathbb Eg_{\epsilon,r}-\mathbb Ef_{\epsilon,r}|=o(1).
\tag{RC.10}
```

The audited conference pressure theorem, separately for both orientations,
then gives

```math
\mathbb Eg_{\epsilon,r}=h_\beta r+o(r),
\qquad
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}.
\tag{RC.11}
```

## 3. Main theorem

Let `U_r^K` be the uniform joint output law conditioned on `K_r`, and put

```math
\mathcal R^K_{\lambda,r}
=-{1\over\lambda}\log
 \mathbb E_{U_r^K}e^{-\lambda f_{\epsilon,r}(B)}.
\tag{RC.12}
```

### Theorem RC.1 (speed-`r^2` lower tail and all-fixed-tilt annealing)

For every fixed `eta>0`, there is `c_(beta,kappa,eta)>0` such that, for all
large `r`,

```math
\boxed{
U_r^K\{f_{\epsilon,r}(B)\le(h_\beta-\eta)r\}
\le\exp\{-c_{\beta,\kappa,\eta}r^2\}.}
\tag{RC.13}
```

Consequently, for every fixed `lambda>0`,

```math
\boxed{
{\mathcal R^K_{\lambda,r}\over r}\longrightarrow h_\beta.}
\tag{RC.14}
```

More generally, (RC.14) holds for every positive sequence
`lambda_r=o(sqrt(r))`.

**Proof.**  For either orientation, (RC.11) implies, for all large `r`,

```math
(h_\beta-\eta)r
\le\mathbb Eg_{\epsilon,r}-{\eta r\over2}.
\tag{RC.15}
```

On `K_(epsilon,r)`, `f=g`.  Hence (RC.9), (RC.6), and conditioning give

```math
\begin{aligned}
U_r^K\{f\le(h_\beta-\eta)r\}
&\le {\sum_{\epsilon=+-1}
 \Pr_B\{g_{\epsilon,r}-\mathbb Eg_{\epsilon,r}
               \le-\eta r/2\}
 \over2\Pr(\mathcal K_r)}\\
&\le \exp\{-c_{\beta,\kappa,\eta}r^2\},
\end{aligned}
\tag{RC.16}
```

after reducing the constant.  This proves (RC.13).

For the moment statement, center at `Eg_(epsilon,r)`.  Equation (RC.8)
gives

```math
1\le\mathbb E_Be^{-\lambda(g-\mathbb Eg)}
\le e^{C_{\beta,\kappa}\lambda^2}.
\tag{RC.17}
```

The omitted part of this centered moment is, by Cauchy--Schwarz, (RC.6),
and (RC.8) at `2lambda`, at most

```math
\begin{aligned}
\mathbb E_B[1_{\mathcal K_{\epsilon,r}^c}
 e^{-\lambda(g-\mathbb Eg)}]
&\le\Pr(\mathcal K_{\epsilon,r}^c)^{1/2}
 \left(\mathbb Ee^{-2\lambda(g-\mathbb Eg)}\right)^{1/2}\\
&\le\exp\{-c_{\beta,\kappa}r/2
          +2C_{\beta,\kappa}\lambda^2+O(1)\}.
\end{aligned}
\tag{RC.18}
```

For fixed `lambda`, this is `o(1)`.  Thus restriction to `K`, normalization
by `Pr(K)=1-o(1)`, and averaging the two orientations change the logarithm
of the centered negative moment by only `O_(beta,kappa,lambda)(1)`.
Together with (RC.11),

```math
\mathcal R^K_{\lambda,r}
=h_\beta r+o(r)+O_{\beta,\kappa}(\lambda).
\tag{RC.19}
```

This proves (RC.14).  For an arbitrary positive sequence
`lambda_r=o(sqrt(r))`, one must also handle exponentially small tilts rather
than divide an unspecified normalization error by `lambda_r`.  If
`lambda_r>=e^(-c_1r)`, with `c_1` smaller than the exponent in (RC.18), the
same calculation gives

```math
\mathcal R^K_{\lambda_r,r}
=h_\beta r+o(r)+O_{\beta,\kappa}(\lambda_r).
\tag{RC.20}
```

If `0<lambda_r<e^(-c_1r)`, then `0<=f<=kappa r` on `K`, and the
bounded-range exponential lemma under the conditioned law gives

```math
0\le
 \mathbb E[f\mid\mathcal K]
 +{1\over\lambda_r}\log
   \mathbb E[e^{-\lambda_r f}\mid\mathcal K]
\le {\lambda_r\kappa^2r^2\over8}=o(r).
\tag{RC.21}
```

The conditioned mean is `h_beta r+o(r)`, proving the claimed full range.
`square`

## 4. Every regular output also has the bounded-cap thin tail

The conditioning is not a purely analytic relaxation.  Every retained
matrix is still a complete hollow sign signing.  From (RC.4), with order
`N=2r`,

```math
\begin{aligned}
Q(S_{\epsilon,B})
&\le {N\over2}\|S_{\epsilon,B}\|_{2\to2}\\
&\le r\,{\kappa\sqrt{2r}\over\beta}
={\kappa\over2\beta}N^{3/2}.
\end{aligned}
\tag{RC.22}
```

Thus Theorem 36.26 applies uniformly with the fixed constant

```math
C_{\beta,\kappa}={\kappa\over2\beta}.
\tag{RC.23}
```

Every output in `K_r` therefore has a common fixed-rate two-sided spin tail.
The theorem is stronger than merely saying that thin-tail outputs dominate:
it identifies an operator-regular thin-tail subclass on which bridge-pressure
lower deviations already have the conjectural speed `r^2`.

## 5. Exact localization of a possible finite-tilt phase

Let `U_r` denote the full bridge law.  Decompose its negative moment as

```math
\mathbb E_{U_r}e^{-\lambda f}
=U_r(\mathcal K_r)\,
  \mathbb E_{U_r^K}e^{-\lambda f}
 +\mathbb E_{U_r}[1_{\mathcal K_r^c}e^{-\lambda f}].
\tag{RC.24}
```

The first term has exponential rate `-lambda h_beta r` for every fixed
`lambda`, by RC.1.  Therefore any departure of the unconditioned (or the
broader cut-cap-conditioned) soft minimum from the typical branch must be
caused entirely by `K_r^c`:

```math
\boxed{
\text{finite-tilt phase}
\quad\Longrightarrow\quad
\mathbb E[1_{\mathcal K_r^c}e^{-\lambda f}]
\text{ has the competing exponential rate}.}
\tag{RC.25}
```

Since `U_r(K_r^c)=e^{-Omega(r)}`, this requires a joint speed-`r` tradeoff
between operator irregularity and pressure reduction.  Merely changing the
empirical spectral law on the regular class cannot create the phase.

This is a strict sharpening of the lower-LDP target in
`conference_negative_moment_phase_boundary.md`: the only still-uncontrolled
part of that LDP is the operator-irregular sector.

## 6. Archive comparison and scope

1. Section 6 of `artifacts/two_temperature_bridge_audit.md` constructs the
   same extension `g` and uses it to prove an unconditioned exact-rate
   theorem only on a small fixed-`lambda` interval; the exponentially rare
   complement could dominate larger negative moments.  RC.1 asks a
   different, newly legitimate question--condition on the regular event--
   and obtains the speed-`r^2` tail and all-fixed-tilt result there.
2. The frozen thin-tail-conditioned theorem uses a broader cut-cap event.
   RC.1 does not silently replace that theorem: it identifies a stricter
   overwhelming subclass and proves that any unresolved phase is supported
   on the difference.
3. The result does not prove the full bridge lower-LDP, a same-temperature
   recurrence, or pressure-minimizer compensation.  Conference children are
   still not known contracted-temperature pressure minimizers.
4. The March 2026 Gaussian spin-glass paper suggests speed-`N^2` lower
   deviations at zero field.  RC.13 proves precisely that speed on the
   operator-regular conference sector by a different elementary mechanism;
   it does not transfer the Gaussian Parisi theorem.

**Research judgment.**  Do not spend further effort on regular bridge
outputs: their entire finite-disorder-temperature hierarchy is now
controlled.  The exact next falsifier/theorem is whether
`K_r^c` contains a speed-`r` family whose pressure gain overcomes its entropy
cost.  That is a joint operator-irregularity/pressure large-deviation
problem, not a generic bridge or scalar thin-tail problem.
