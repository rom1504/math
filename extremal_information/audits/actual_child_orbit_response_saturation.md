# The finite actual-child posterior saturates its symmetry quotient

Status: **rigorous exact finite actual-minimizer theorem and reproducible
modular certificate**.  This note tests the orbit quotient of OQ.1 on the
two certified order-eight pressure-minimizer classes.  In this instance the
quotient is not merely sufficient: the scalar averaged-posterior response
distinguishes every cell of the largest simultaneous-symmetry quotient.

This is a finite theorem at raw temperature `t=3`.  It proves neither an
all-order orbit bound nor a bound on canonical row-product regret.

## 1. Actual latent experiment

Let `A` be either of the two order-eight matrices `A_0,A_1` in the certified
minimizer classification, and let `D` be the unique order-two child.  Both
`A_0,A_1` are exact pressure minimizers at every raw `t>=3`; `D` is the exact
order-two minimizer.  Use projective spins, represented by first coordinate
`+1`, and put

```math
Q=xy^{\mathsf T},
\qquad E(Q)=H_A(x)+H_D(y).                           \tag{OR.1}
```

For an indeterminate `z`, write

```math
C_k(z)={z^k+z^{-k}\over2},
\qquad
w_Q(z)=C_{E(Q)}(z),
\qquad
K_Q(B;z)=C_{\langle B,Q\rangle}(z).                \tag{OR.2}
```

The unnormalized output likelihood is

```math
O_B(z)=\sum_Qw_Q(z)K_Q(B;z).                        \tag{OR.3}
```

Channel constants common to all `Q` and `B` have been omitted; they cancel
from every posterior and disorder escort below.  At negative-disorder
exponent `lambda=1`, direct Bayes substitution gives

```math
\boxed{
{\bar\mu_z(Q)\over\mu_z(Q)}
=\left(\sum_Rw_R(z)\right)
 {\displaystyle\sum_B K_Q(B;z)O_B(z)^{-2}
  \over\displaystyle\sum_BO_B(z)^{-1}}.}          \tag{OR.4}
```

Thus the exact minimal scalar response partition is the level-set partition
of

```math
N_Q(z)=\sum_B{K_Q(B;z)\over O_B(z)^2}.              \tag{OR.5}
```

Every term is a rational function in `z` over `Q`.

## 2. The correct actual-child symmetry

For signed permutation matrices `U,V`, define the simultaneous similarity
group

```math
\mathcal G^{\pm}(A,D)
=\{(U,V):\text{ for some }\delta\in\{\pm1\},
 U^{\mathsf T}AU=\delta A,
 V^{\mathsf T}DV=\delta D\}.                       \tag{OR.6}
```

It acts projectively on `(x,y)`, hence on `Q`.  If `(U,V)` has sign `delta`,
then

```math
E(Ux,Vy)=\delta E(x,y),
\qquad
\langle UBV^{\mathsf T},UQV^{\mathsf T}\rangle
=\langle B,Q\rangle.                               \tag{OR.7}
```

Since `C_k=C_(-k)` and the fair bridge law is invariant, (OR.2)--(OR.5)
show that `N_Q` is constant on every `\mathcal G^{\pm}(A,D)` orbit.  This is the
actual signed/anti-signed version of the abstract quotient OQ.1.  The
anti-isomorphisms matter: ordinary factor automorphisms give `38` and `44`
joint cells below, while simultaneous anti-isomorphisms pair them into the
true `19` and `22` response cells.

For comparison, define the rooted child profile

```math
\mathcal P_A(x)
=\operatorname {hist}_{z\in\{\pm1\}^8/\{\pm1\}}
       \bigl(H_A(z),|\langle x,z\rangle|\bigr).     \tag{OR.8}
```

It is a polynomial-bin histogram for one declared root, not a table of
bridge pressures.  In both actual order-eight classes, exhaustive exact
enumeration shows that its level sets are exactly the signed-automorphism
orbits of projective child spins.  Therefore
`(\mathcal P_A(x),H_D(y))`, followed by the simultaneous-complement
identification in (OR.6), is an exact finite posterior quotient here.

## 3. Exact saturation theorem

**Theorem OR.1 (actual order-eight orbit-response saturation).**  Pair each
certified order-eight pressure minimizer with the unique order-two child,
take orientation `+`, raw temperature `t=3`, and negative-disorder exponent
`lambda=1`.  Then:

| left class | projective signed automorphisms | left spin/profile cells | ordinary joint cells | simultaneous `\mathcal G^{\pm}` cells | distinct values of `\bar\mu/\mu` |
|---:|---:|---:|---:|---:|---:|
| `A_0` | `24` | `19` | `38` | `19` | `19` |
| `A_1` | `16` | `22` | `44` | `22` | `22` |

In particular, the scalar posterior response separates **every**
simultaneous-similarity orbit.  The combined-energy label has only `12`
values.  It is strictly too coarse: individual energy shells contain as
many as four response values for `A_0` and six for `A_1`.

*Exact certificate.*  Evaluate the rational functions at `z=2` in the
prime field

```math
\mathbb F_p,\qquad p=1,000,003.                     \tag{OR.9}
```

Here

```math
C_k(2)=\frac{2^k+2^{-k}}2
```

is well-defined modulo `p`.  Exhaustive enumeration of all `2^16` bridges
checks that every `O_B(2)` is nonzero modulo `p`; the inverse-escort
normalizer residues are respectively

```text
431332, 833738.
```

It then evaluates (OR.5) exactly in `F_p`.  The residues are constant on
the simultaneous-similarity orbits and are pairwise distinct across all
`19`, respectively `22`, orbits.  Distinct residues imply distinct rational
numbers at `z=2`, hence distinct rational functions in `Q(z)`.  If two such
functions agreed at the actual `z=e^3`, their nonzero difference numerator,
a Laurent polynomial over `Q`, would have the transcendental number `e^3`
as a root.  This is impossible.  Therefore the response values are also
pairwise distinct at the actual raw temperature `t=3`. `square`

The stored SHA-256 digests of the sorted orbit-response residues are

```text
A_0: 2865ec1d90548756777b363a16cf926423e59dfb410c48b62ca7cf1ea0700b3c
A_1: bc9f892f51f3bd4738c86e83dfe3a87acbfc57a716cdd89fb776cf62615f43f5
```

## 4. What this decides

The theorem gives a positive and a negative answer at the finite actual
frontier.

1. **Positive:** a rooted energy--overlap feature plus the exact simultaneous
   signed-similarity congruence captures all averaged posterior retuning for
   these actual children.  This is the optimizer-specific analogue of the
   magnetization orbit count in the diffuse BSC phase.
2. **Negative:** the response uses every cell of that quotient.  Neither the
   combined energy, the unrooted child spectra, nor symmetry alone supplies
   a further collapse in this experiment.

The result does not show that rooted profiles are orbit-complete at other
orders.  Even when they are, their number can be exponential; storing one
profile label per latent word would then reconstruct the full landscape.
The exact next all-order obligation is therefore not to add another scalar
shell statistic, but to prove one of:

```math
\log |\{\mathcal G^\pm\text{-orbits}\}|=o(N)
```

for actual optimizing children, an approximate quotient with `o(N)`
symmetry-breaking KL as in (OQ.9), or a theorem showing that the actual
response coalesces across exponentially many orbit cells.  Nothing here
controls `J-I^\leftarrow`, target reach, or a Level-6 recurrence.

## 5. Reproduction

From the repository root:

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_orbit_posterior_quotient.py
```

The machine-readable certificate is
[`../../computations/results/actual_child_orbit_posterior_quotient.json`](../../computations/results/actual_child_orbit_posterior_quotient.json).
