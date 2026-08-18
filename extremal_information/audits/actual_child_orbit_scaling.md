# Exact scaling of the actual-child orbit quotient through order eight

Status: **rigorous exact finite actual-minimizer theorem with a reproducible
modular certificate**.  This audit extends the order-eight experiment OR.1 to
every order for which the repository contains an exhaustive cap-minimizer
classification.  It uses no conference or Paley surrogate.

The conclusion is deliberately finite.  It gives no all-order orbit bound,
row-product-regret theorem, or recurrence.

## 1. Which stored children are actual pressure minimizers?

Write

```math
Z_A(t)=2^{-(n-1)}
       \sum_{x\in\{\pm1\}^n/\{\pm1\}}\cosh(tH_A(x)).       \tag{OS.1}
```

The stored files `m3_minimizer_orbits.json` through
`m8_minimizer_orbits.json` are exhaustive classifications of all cap
minimizers modulo signed permutation and global complement.  A signing which
is not a cap minimizer has cap at least `M_n+2`.  Hence, for `3<=n<=8` and
`t>=3`,

```math
{\cosh((M_n+2)t)\over 2^{n-1}\cosh(M_nt)}
\ge {e^{2t}\over2^n}>1.                               \tag{OS.2}
```

Thus no non-cap-minimizer can minimize (OS.1).  There is one classified cap
class at orders three through six.  The two order-eight classes have the same
absolute-energy histogram.  At order seven the three projective
absolute-energy histograms are

| class | `|H|=1` | `3` | `5` | `7` | `9` |
|---:|---:|---:|---:|---:|---:|
| `0` | 15 | 21 | 15 | 9 | 4 |
| `1` | 21 | 21 | 7 | 8 | 7 |
| `2` | 21 | 13 | 15 | 12 | 3 |

The exact differences between the unnormalized partition sums of classes
zero/one and class two are

```math
\begin{aligned}
S_0(t)-S_2(t)
 &=\cosh(9t)-3\cosh(7t)+8\cosh(3t)-6\cosh(t),\\
S_1(t)-S_2(t)
 &=4\cosh(9t)-4\cosh(7t)-8\cosh(5t)+8\cosh(3t).
                                                               \tag{OS.3}
\end{aligned}
```

Both are positive for `t>=3`: indeed
`cosh(9t)/cosh(7t)>=e^(2t)/2>9` and
`cosh(7t)>=cosh(5t)>=cosh(t)`.  Therefore the actual pressure-minimizer
classes at `t=3` are

```text
n=2: unique;  n=3,4,5,6: class 0;
n=7: class 2; n=8: classes 0 and 1.                 (OS.4)
```

## 2. Exact quotient and response calculation

Pair each eligible left child with the unique order-two child, use projective
spins, raw temperature `t=3`, and negative-disorder exponent `lambda=1`.
The definitions of the signed-similarity quotient and the inverse-escort
posterior response are exactly (OR.1)--(OR.6) in
`actual_child_orbit_response_saturation.md`.

For each left child, enumerate:

1. projective signed automorphisms;
2. their orbits on left spins;
3. rooted energy--absolute-overlap profiles;
4. simultaneous signed/anti-signed orbits of the left/right latent pair;
5. the inverse-escort posterior-response rational functions.

The last functions are evaluated in `F_p` at `z=2`, where
`p=1,000,003`.  Distinct residues imply distinct rational functions over
`Q(z)`.  A nonzero rational function over `Q` cannot vanish at the
transcendental actual point `z=e^3`; therefore this proves distinction at the
actual temperature.

## 3. Finite scaling theorem

**Theorem OS.1 (actual-child orbit scaling through the classified
frontier).**  For the pressure minimizers (OS.4), the exact counts are:

| left order/class | projective spins | signed-aut group | spin orbits | rooted-profile cells | simultaneous cells | posterior-response cells |
|---:|---:|---:|---:|---:|---:|---:|
| `2/0` | 2 | 2 | 2 | 2 | 2 | 2 |
| `3/0` | 4 | 6 | 2 | 2 | 4 | 3 |
| `4/0` | 8 | 4 | 5 | 5 | 5 | 5 |
| `5/0` | 16 | 10 | 4 | 4 | 4 | 4 |
| `6/0` | 32 | 60 | 4 | 4 | 4 | 4 |
| `7/2` | 64 | 12 | 13 | 13 | 26 | 26 |
| `8/0` | 128 | 24 | 19 | 19 | 19 | 19 |
| `8/1` | 128 | 16 | 22 | 22 | 22 | 22 |

At every eligible class, the rooted energy--absolute-overlap partition is
exactly the signed-automorphism orbit partition of the left child.  Except at
order three, the negative posterior response separates every cell of the
largest simultaneous signed-similarity quotient.

At order three there are exactly three response values at `t=3`.  More
strongly, two of the four simultaneous-cell rational functions agree
identically, not just modulo `p`.  (The remaining functions need not be
distinct at every temperature; all responses coincide at `t=0`.)  To
certify the identity, group the 64 bridges by the exact histogram

```math
\operatorname{hist}_Q\bigl(|E(Q)|,|\langle B,Q\rangle|\bigr),              \tag{OS.5}
```

which determines the output denominator as a rational function of `z`.
There are six denominator types.  Within each type, the coefficient of every
kernel `C_k(z)` in the difference of the two colliding response numerators is
exactly zero.  Thus their rational functions agree identically. `square`

## 4. What the data decide

The small symmetric orders do have a strict quotient, but there is no stable
finite-state trend.  The number of rooted cells rises from four at order six
to 13 at order seven and 19 or 22 at order eight.  Moreover, aside from the
single order-three identity, the scalar negative posterior uses every
available simultaneous cell.  Thus the finite actual data provide no
response coalescence beyond exact symmetry.

The repository contains a certified order-nine **cap** minimizer but no
exhaustive classification of thermal-pressure minimizers at that order.  It
cannot be admitted as an actual child without a new broad classification.
As an explicitly ineligible diagnostic, that representative has 256
projective spins, signed-automorphism group order two, 144 spin orbits, and
116 rooted-profile cells; the rooted-profile and automorphism partitions are
already unequal.

Therefore an all-order actual-child quotient still requires a genuinely new
theorem: subexponential orbit/profile growth, approximate response
coalescence across exponentially many cells, or a different coordinate mark
whose symmetry-breaking information is `o(n)`.  Finite rooted profiles do not
supply such a theorem by themselves.

## 5. Reproduction

From the repository root:

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_orbit_scaling.py
```

The exact certificate is written to
[`../../computations/results/actual_child_orbit_scaling.json`](../../computations/results/actual_child_orbit_scaling.json).
