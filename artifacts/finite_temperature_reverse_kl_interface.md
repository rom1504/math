# Joint reverse-KL compensation for signing free energy

Status: **exact theorem and finite exhaustive diagnostic**.  The theorem
identifies a signing-specific finite-temperature condition that would turn
the contracted annealed recurrence into same-temperature almost
subadditivity.  The condition is not proved uniformly.  It is strictly more
information than the scalar pressures, but it does not optimize a parent
ground state or select a spin configuration.

## 1. One normalization throughout

For a signing `A` of order `n`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad E_n={n\choose2},
```

and use the normalized partition function

```math
\overline Z_n(A,t)
=2^{-n}\sum_{x\in\{-1,1\}^n}\cosh(tH_A(x)),
\qquad
F_n(t)=\min_A\log\overline Z_n(A,t).               \tag{1.1}
```

The pressure in the proposed finite-temperature formulation is

```math
\Phi_n(\beta)
=\left(1+{1\over n}\right)\log2+{1\over n}P_n(\beta),
\qquad
P_n(\beta)=F_n\!\left({\beta\over\sqrt n}\right). \tag{1.2}
```

The ground-state squeeze is

```math
\beta {M_n\over n^{3/2}}-\log2
\le {P_n(\beta)\over n}
\le \beta {M_n\over n^{3/2}}.                     \tag{1.3}
```

Indeed, for a fixed `A` with cap `K`,
`tK-n log 2 <= log Zbar_n(A,t) <= tK`; the lower bound uses the two
maximizers `x,-x`.  Equation (1.3) gives the stated uniform squeeze for
`Phi_n`, with the slightly stronger lower entropy correction
`(log 2)/(n beta)`.

Write

```math
\ell(t)=\log\cosh t.                               \tag{1.4}
```

Uniformly random edge signs give

```math
0\le F_n(t)\le E_n\ell(t).                         \tag{1.5}
```

Thus `0 <= P_n(beta) <= beta^2(n-1)/4`, which supplies the linear
boundedness needed below.

## 2. The joint orientation-and-bridge output law

Let `N=m+n`, `L=mn`, and `t=beta/sqrt(N)`.  Choose child signings `A,D`
and give their augmented Gibbs variables the independent laws

```math
\nu_A(x,\tau_1)
=\frac{\exp(t\tau_1H_A(x))}
       {2^{m+1}\overline Z_m(A,t)},
\qquad
\nu_D(y,\tau_2)
=\frac{\exp(t\tau_2H_D(y))}
       {2^{n+1}\overline Z_n(D,t)}.                \tag{2.1}
```

From one sample define the relative child orientation and latent bridge word

```math
\epsilon=\tau_1\tau_2,
\qquad q_{ij}=\tau_1x_iy_j.                        \tag{2.2}
```

Pass every coordinate of `q` independently through the binary channel

```math
\Pr\{B_{ij}=b\mid q_{ij}\}
=\frac{e^{tbq_{ij}}}{2\cosh t},
\qquad b\in\{-1,1\}.                              \tag{2.3}
```

Let `Pi_(A,D,t)` be the resulting joint law of `(epsilon,B)`.  For the
block signing with internal blocks `A,epsilon D` and bridge `B`, direct
substitution of `tau_2=tau_1 epsilon` proves the exact formula

```math
\boxed{
\Pi_{A,D,t}(\epsilon,B)
=\frac{\overline Z_N(A,\epsilon D,B;t)}
 {2^{L+1}(\cosh t)^L
  \overline Z_m(A,t)\overline Z_n(D,t)}.}          \tag{2.4}
```

Summing (2.4) over all outputs proves normalization and recovers the
annealed factorization

```math
\mathbb E_{\epsilon,B}\overline Z_N(A,\epsilon D,B;t)
=(\cosh t)^L\overline Z_m(A,t)\overline Z_n(D,t), \tag{2.5}
```

where the expectation in (2.5) is uniform.  More importantly, if `U` is
uniform on the same `2^(L+1)` outputs, taking logarithms in (2.4) gives the
**quenched identity**

```math
\boxed{
\mathbb E_U\log\overline Z_N(A,\epsilon D,B;t)
=\log\overline Z_m(A,t)+\log\overline Z_n(D,t)
 +L\ell(t)-D(U\|\Pi_{A,D,t}).}                    \tag{2.6}
```

In particular some deterministic joint choice `(epsilon,B)` is no larger
than the right side.  The reverse relative entropy in (2.6) is exactly the
joint cancellation discarded when Jensen is applied to (2.5); it is not a
sum of separately paid left and right channels.

## 3. Exact same-temperature compensation theorem

Set

```math
\theta={m\over N},
\qquad \beta_m=\beta\sqrt\theta,
\qquad \beta_n=\beta\sqrt{1-\theta}.              \tag{3.1}
```

Let `A,D` range over minimizers defining `P_m(beta_m),P_n(beta_n)` and put

```math
\begin{aligned}
\mathcal D_{m,n}(\beta)
&=\max_{A,D}D(U\|\Pi_{A,D,\beta/\sqrt N}),\\
\mathcal T_{m,n}(\beta)
&=P_m(\beta)-P_m(\beta_m)
 +P_n(\beta)-P_n(\beta_n).                        \tag{3.2}
\end{aligned}
```

Both terms are nonnegative.  Applying (2.6) to a maximizing child pair and
then using the minimum over all parent signings proves

```math
\boxed{
P_N(\beta)\le P_m(\beta)+P_n(\beta)
 +L\ell\!\left({\beta\over\sqrt N}\right)
 -\mathcal T_{m,n}(\beta)-\mathcal D_{m,n}(\beta).}
                                                               \tag{3.3}
```

This separates three effects without changing temperature on the final
right side:

1. `L ell(beta/sqrt(N))` is the annealed bridge cost;
2. `T` is the cost already paid while heating both children from their
   contracted parameters back to `beta`;
3. `D` is joint bridge/orientation cancellation before any absolute value.

The following is consequently an exact positive interface.

**Theorem 3.1 (all-splits compensation criterion).**  Fix `beta>0`.
Suppose there are `delta_beta>0`, `C_beta<infinity`, and `N_beta` such that,
for every `m+n=N>=N_beta`,

```math
\boxed{
\mathcal T_{m,n}(\beta)+\mathcal D_{m,n}(\beta)
\ge mn\ell\!\left({\beta\over\sqrt N}\right)
 -C_\beta N^{1-\delta_\beta}.}                    \tag{3.4}
```

Then `P_n(beta)/n`, and hence `Phi_n(beta)`, converges.

**Proof.**  Equations (3.3)--(3.4) give

```math
P_{m+n}(\beta)
\le P_m(\beta)+P_n(\beta)+C_\beta(m+n)^{1-\delta_\beta}. \tag{3.5}
```

The sequence is nonnegative and `O_beta(n)` by (1.5).  The standard
balanced-tree proof of the almost-subadditive lemma applies: starting with
blocks of size `k`, the total error per vertex at merging level `j` is
`O_beta((2^j k)^(-delta_beta))`.  Its sum is `O_beta(k^(-delta_beta))`.
Choosing `k` along a subsequence attaining `liminf P_k(beta)/k`, and then
letting `k` tend to infinity, gives `limsup <= liminf`.  This proves the
claim.  Equivalently, the defect in (3.5) has Hammersley-summable density.

An unspecified `o(N)` in (3.4) is not enough: slowly oscillating linear
sequences can have such a two-block defect.  A power saving is convenient;
the proof only needs the corresponding dyadic error densities to be
summable.

## 4. Comparable splits suffice

The all-splits quantifier can be reduced to a fixed balanced window.  First
note the exact raw-temperature cavity bounds

```math
F_n(t)\le F_{n+1}(t)\le F_n(t)+n\ell(t).           \tag{4.1}
```

For the lower bound, average the new spin and use
`cosh(tH) cosh(tL) >= cosh(tH)`.  For the upper bound, start with an
`F_n(t)` minimizer and average the new signed row; the average partition is
`Zbar_n(cosh t)^n`.

At fixed scaled `beta`, (4.1), the entropy squeeze, and
`sqrt(1+1/n)-1 <= 1/(2n)` give the uniform adjacent estimate

```math
-{1\over2}\left(\log2+{\beta^2\over4}\right)
\le P_{n+1}(\beta)-P_n(\beta)
\le {\beta^2\over2}.                              \tag{4.2}
```

For the nontrivial direction, let `A` minimize `F_n(beta/sqrt(n+1))` and
write `K=max_x|H_A(x)|`.  Then

```math
K\le {n(\log2+\beta^2/4)\over \beta/\sqrt{n+1}},
```

and changing the raw parameter from `beta/sqrt(n+1)` to `beta/sqrt n`
costs at most `K` times their difference, which is at most the constant in
(4.2).

It follows that Theorem 3.1 remains true if (3.4) is assumed only for

```math
{N\over4}\le m,n\le{3N\over4}.                    \tag{4.3}
```

Indeed, partition a large order into leaves of sizes `k` and `k+1` and
merge them in a balanced binary tree.  Every internal split lies in (4.3),
the accumulated error per vertex is `O_beta(k^(-delta_beta))`, and (4.2)
makes the two leaf ratios asymptotically equal.  Choosing `k` along a
`liminf` subsequence completes the same proof.

## 5. Replica form and comparison with Gaussian interpolation

Order the output bits as the orientation followed by the `L` bridge bits.
The chain rule for reverse relative entropy rewrites (2.6) as

```math
\mathcal D_{m,n}(\beta)
=\frac12\sum_{j=0}^{L}
 \mathbb E_{U_{<j}}[-\log(1-s_j^2)],               \tag{5.1}
```

where `s_j` is the conditional bias of output bit `j` under `Pi`.  For a
bridge bit,

```math
s_j=\tanh(t)r_j,
\qquad
r_j=\mathbb E[q_j\mid (\epsilon,B)_{<j}].          \tag{5.2}
```

Consequently the concrete squared-response condition

```math
\mathcal T_{m,n}(\beta)
 +{\tanh^2(\beta/\sqrt N)\over2}
   \sum_{j=1}^{mn}\mathbb E_U r_j^2
\ge mn\ell(\beta/\sqrt N)-C_\beta N^{1-\delta_\beta} \tag{5.3}
```

is sufficient for (3.4).  Unlike independently bounded response channels,
all `r_j` are posterior means in one evolving joint Gibbs law and their
cancellation is taken before the logarithm.

This is the exact point at which the Guerra--Toninelli mechanism does not
transfer automatically.  Gaussian stability under interpolation and
integration by parts turn the corresponding derivative into an overlap
square with a forced sign.  Fixed-modulus Bernoulli signings are not closed
under that rescaling.  Random vertex switching of any fixed signing has
mean-zero, identity edge covariance while leaving its partition function
unchanged, so covariance-level Gaussian interpolation cannot distinguish
different signing pressures.  Here the missing information survives as the
higher-order output divergence in (5.1), but no universal lower bound for it
is known.  The comparison is with Guerra and Toninelli's original
[Gaussian interpolation theorem](https://arxiv.org/abs/cond-mat/0204280),
whose thermodynamic-limit conclusion uses same-temperature subadditivity of
the quenched free energy.

## 6. Exhaustive small-order diagnostic

The script
`computations/check_finite_temperature_reverse_kl.py` exhausts one
representative of every switching orbit of every child and parent signing,
all low-parameter minimizing child pairs, both relative orientations, and
every bridge through total order seven.  The stored result is
`computations/results/finite_temperature_reverse_kl_small.json`.  It is
reproduced by

```text
.venv/bin/python computations/check_finite_temperature_reverse_kl.py \
  --max-n 7 \
  --output computations/results/finite_temperature_reverse_kl_small.json
```

Its SHA-256 digest is
`951b14b968636c0394fc5d0abf0f40095ab7b568919b0f3474c1195b118fc6e3`.
The enumeration is exhaustive; the transcendental evaluations use ordinary
double precision and are therefore classified as **exhaustive numerical**.

For the most balanced split at each order, define the uncancelled interface
defect

```math
G_{m,n}(\beta)
=mn\ell(\beta/\sqrt N)
 -\mathcal T_{m,n}(\beta)-\mathcal D_{m,n}(\beta). \tag{6.1}
```

The computed values are:

| `beta` | `N=4` | `N=5` | `N=6` | `N=7` |
|---:|---:|---:|---:|---:|
| 0.25 | 0.015665 | 0.015684 | 0.015705 | 0.015875 |
| 0.5 | 0.063080 | 0.063362 | 0.063690 | 0.066232 |
| 1 | 0.254574 | 0.255973 | 0.258160 | 0.289807 |
| 2 | 0.912557 | 0.849614 | 0.797231 | 1.037759 |
| 4 | 2.497477 | 1.871166 | 1.492197 | 2.452983 |
| 8 | 5.820336 | 3.547893 | 2.550662 | 5.105247 |

Thus exact zero-defect compensation is not automatic.  It already fails
exactly at `1+1`: both `T` and `D` vanish while
`G_(1,1)(beta)=ell(beta/sqrt 2)>0`.  Every tested nontrivial balanced case
also has positive `G`.  On the other hand, these orders do not distinguish a
bounded or other power-saving defect from a linear one.  For `beta<=1`, the
nearly order-independent values are compatible with `G=O_beta(1)`; the
larger-`beta` data have substantial parity/order variation and support no
scaling claim.

The test also computes the actual same-temperature scalar defect

```math
P_N(\beta)-P_m(\beta)-P_n(\beta).                  \tag{6.2}
```

At the balanced order-six split it is respectively
`0.014143, 0.041401, 0.015572, -0.645080, -3.203770, -8.896135`
for the six displayed temperatures.  In particular, at `beta>=2` an
exceptional parent already gives same-temperature subadditivity while the
uniform-output interface (6.1) still has a positive defect.  This proves at
finite order that the reverse-KL criterion is sufficient but strictly
stronger than merely finding the best parent bridge.

## 7. Research judgment and falsification criterion

The weakest scalar addition needed by this architecture is the
Hammersley-summable same-temperature defect (3.5).  Stating that bound alone
does not explain why complete quadratic signings should satisfy it.  The
joint law (2.4) supplies a noncircular signing-specific mechanism: `D` is
computed entirely from contracted child minimizers and a fixed noisy
rank-one channel, before the parent minimum is taken.

It is therefore genuinely weaker than full parent ground-state
maximization at each fixed finite `beta`: it controls an averaged logarithm,
uses no extremal spin, and its constants may deteriorate arbitrarily as
`beta` tends to infinity.  It is not equivalent to parent minimization.
However, it is stronger than the desired scalar same-temperature inequality,
because a rare good bridge can beat the uniform logarithmic average.  The
finite audit exhibits precisely that separation.

The sharp next theorem, if this route is retained, is (3.4) on the balanced
window (4.3), or the squared-response version (5.3).  A sequence of balanced
orders for which

```math
G_{m,n}(\beta)\ge c_\beta N
```

for some fixed `beta,c_beta>0` would falsify this uniform-output mechanism,
though not finite-temperature convergence itself.  Conversely any uniform
`O_beta(N^(1-delta))` bound would be primary progress: it would prove the
thermodynamic limit at that `beta`, and doing so for every fixed `beta`
would prove convergence of `M_n/n^(3/2)` through (1.3).
