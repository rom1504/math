# The composition law for exposed Boolean flatness

**Status.** Rigorous exact chain rule, hierarchical recovery theorem, and a
pumpable tensor benchmark.  The carrier applies to a **specified exposed
spherical optimizer**.  It does not select that optimizer and does not claim
to compress a switching family of all spherical response channels.

The composable state is elementary but nontrivial:

```math
(E,L)=(\|u\|_2^2,\|u\|_1).
```

Both coordinates add under disjoint block composition.  Their normalized
ratio separates two distinct recovery resources: within-block Boolean
flatness and allocation of Euclidean amplitude among blocks.

## 1. Exact one-level chain rule

Partition a vector `u in R^N` into blocks `u_i in R^(n_i)`, with

```math
N=\sum_i n_i,
\qquad \|u\|_2^2=N,
\qquad \lambda_i={n_i\over N}.                       \tag{FC.1}
```

Define the block RMS amplitude

```math
\rho_i={\|u_i\|_2\over\sqrt{n_i}},
\qquad \sum_i\lambda_i\rho_i^2=1.                   \tag{FC.2}
```

If `rho_i>0`, normalize `v_i=u_i/rho_i`, so `||v_i||_2^2=n_i`, and put

```math
\phi_i=1-{\|v_i\|_1\over n_i}.                       \tag{FC.3}
```

For a zero block take `phi_i=0`; it is multiplied by `rho_i`.  The global
flatness is

```math
\phi(u)=1-{\|u\|_1\over N}.                          \tag{FC.4}
```

### Theorem FC.1 (flatness chain rule)

```math
\boxed{
\phi(u)
={1\over2}\sum_i\lambda_i(\rho_i-1)^2
 +\sum_i\lambda_i\rho_i\phi_i.}                    \tag{FC.5}
```

The first term is the **amplitude-allocation defect**

```math
\mathcal A(\lambda,\rho)
:={1\over2}\sum_i\lambda_i(\rho_i-1)^2.             \tag{FC.6}
```

It vanishes exactly when every nonempty block receives its dimension-
proportional RMS amplitude `rho_i=1`.  The second term is the transported
within-block flatness.

#### Proof

Additivity of `l_1` gives

```math
\phi(u)=1-\sum_i\lambda_i\rho_i(1-\phi_i)
=1-\sum_i\lambda_i\rho_i+
  \sum_i\lambda_i\rho_i\phi_i.                      \tag{FC.7}
```

Using `sum lambda_i=sum lambda_i rho_i^2=1`,

```math
1-\sum_i\lambda_i\rho_i
={1\over2}\sum_i\lambda_i(\rho_i-1)^2.              \tag{FC.8}
```

Substitution proves (FC.5). `square`

Thus local Booleanity alone is insufficient.  Even if every `phi_i=0`, an
unequal Euclidean amplitude allocation creates global rounding loss.

## 2. Exact hierarchical expansion

Let a rooted block tree partition the coordinates recursively.  For a node
`v`, let `N_v` be its block size and write

```math
R_v={\|u_v\|_2/\sqrt{N_v}\over
          \|u_{root}\|_2/\sqrt N}.                   \tag{FC.9}
```

The root has `R_root=1`.  At each internal node apply FC.1 to its normalized
restriction and call the local allocation defect `A_v`.  Define the
transport weight

```math
\omega_v={N_v\over N}R_v.                            \tag{FC.10}
```

If `u_v=0`, set its local defect and every descendant contribution to zero;
then `omega_v=0`, so the convention does not change any identity.

### Theorem FC.2 (tree chain rule)

For any finite partition tree,

```math
\boxed{
\phi(u)=
\sum_{v\text{ internal}}\omega_v\mathcal A_v
+\sum_{\ell\text{ leaf}}\omega_\ell\phi_\ell.}     \tag{FC.11}
```

At every fixed depth,

```math
\sum_{v\text{ at that depth}}\omega_v\le1.       \tag{FC.12}
```

#### Proof

Iterate (FC.5).  Along a root-to-`v` path, the products of relative size
fractions and relative RMS amplitudes telescope to `(N_v/N)R_v`, proving
(FC.11).  At a fixed depth the blocks partition the root and

```math
\sum_v{N_v\over N}R_v^2=1.
```

Cauchy--Schwarz with `sum_v N_v/N=1` proves (FC.12). `square`

In particular, if a depth-`D` hierarchy obeys

```math
\mathcal A_v\le a_d\quad(v\text{ at internal depth }d),
\qquad \phi_\ell\le b,                               \tag{FC.13}
```

then

```math
\phi(u)\le b+\sum_{d=0}^{D-1}a_d.                   \tag{FC.14}
```

This is the explicit accumulation law.  Allocation errors add by level
inside flatness before the final square-root rounding loss; they are not
separately paid in the energy.

## 3. Boolean recovery under repeated composition

The following formulation is independent of the Hadamard notation.  Let

```math
F(u)={1\over2}u^TMu+h^Tu,
\qquad \|M\|_{op}\le\Lambda,                         \tag{FC.15}
```

and compare its maximum on the sphere `||u||^2=N` with its maximum on the
Boolean cube.  Put

```math
\kappa={\|h\|_2\over\Lambda\sqrt N}                 \tag{FC.16}
```

when `Lambda>0`.

### Theorem FC.3 (hierarchical exposed recovery)

Assume `Lambda>0`.  Suppose an exposed spherical vector `u` has value at least
`S-xi Lambda N` and has a partition tree satisfying (FC.13).  Then

```math
\boxed{
0\le S-B
\le\Lambda N\left[
\xi+(1+\kappa)
\sqrt{2\left(b+\sum_{d=0}^{D-1}a_d\right)}
\right].}                                             \tag{FC.17}
```

#### Proof

The tree chain rule and (FC.14) bound the global flatness.  Coordinatewise
sign rounding has Euclidean distortion

```math
\|u-\operatorname{sgn}u\|_2
=\sqrt{2N\phi(u)}.                                   \tag{FC.18}
```

The quadratic loss is at most `Lambda N sqrt(2phi)`, and the linear loss is
at most `||h|| sqrt(N) sqrt(2phi)`.  Add the exposure error `xi Lambda N`.
This proves (FC.17). `square`

For the Hadamard trust response, `Lambda=r` and
`kappa<=mp/r=c`, recovering BS.1.  A growing-depth composition has
`o(Lambda N)` Boolean loss whenever

```math
\xi+b+\sum_d a_d=o(1)                                \tag{FC.19}
```

and `kappa=O(1)`.

Equation (FC.19) is severe but exact.  Merely requiring every level defect
to tend to zero is not enough when depth diverges; the defects must be
summable on the realized exposed path distribution.

## 4. Pumpable amplitude-localization obstruction

Take an equal binary split at every node and fix `delta in (0,1)`.  Give the
two children relative RMS amplitudes

```math
\rho_+=\sqrt{1+\delta},
\qquad \rho_-=\sqrt{1-\delta}.                       \tag{FC.20}
```

They satisfy `(rho_+^2+rho_-^2)/2=1`.  Put

```math
s={\rho_++\rho_-\over2}<1,
\qquad \mathcal A=1-s>0.                             \tag{FC.21}
```

At depth `D`, take every scalar leaf to be positive Boolean and let the
amplitude of word `omega in {+,-}^D` be

```math
u_D(\omega)=\prod_{t=1}^D\rho_{\omega_t}.             \tag{FC.22}
```

Then `N=2^D`, `||u_D||_2^2=N`, all leaf flatnesses vanish, but

```math
{\|u_D\|_1\over N}=s^D,
\qquad
\boxed{\phi(u_D)=1-s^D.}                             \tag{FC.23}
```

The chain expansion is exact level by level:

```math
\phi(u_D)=\sum_{j=0}^{D-1}s^j(1-s)=1-s^D.            \tag{FC.24}
```

Thus any fixed nonzero amplitude imbalance is pumpable: repeated
composition drives the exposed optimizer maximally far from Boolean
flatness despite perfect leaves.

For small level-dependent imbalances `delta_j`,

```math
s_j={\sqrt{1+\delta_j}+\sqrt{1-\delta_j}\over2}
=1-{\delta_j^2\over8}+O(\delta_j^4).                 \tag{FC.25}
```

Hence the product remains flat only in the regime where the accumulated
squared imbalance is small; a divergent sum drives the `l_1/l_2` ratio to
zero.

## 5. Benchmark: a rank-one tensor linear landscape

The obstruction is not specific to Hadamard matrices.  Consider the pure
linear landscape on `N=2^D` coordinates

```math
F_D(x)=u_D^Tx.                                        \tag{FC.26}
```

Its spherical and Boolean maxima are exactly

```math
S_D=\sqrt N\|u_D\|_2=N,
\qquad
B_D=\|u_D\|_1=Ns^D.                                 \tag{FC.27}
```

Therefore

```math
{S_D-B_D\over N}=1-s^D\longrightarrow1.             \tag{FC.28}
```

This is a one-state multiplicative weighted-language/tensor-network
benchmark.  Every local factor and every leaf is exactly known, yet the
composition-created amplitude allocation produces an extensive extremal
gap.  The two-scalar exposed carrier `(E,L)` predicts it exactly.

## 6. What the state does and does not accomplish

1. **Exact compositional state for one witness.**  `(E,L)` adds under block
   union, and FC.5 is its normalized gluing law.
2. **Two independent failure modes.**  Local coordinate nonflatness and
   block amplitude imbalance are distinct and both nonnegative.
3. **Dynamic stopping criterion.**  Repeated Boolean recovery requires a
   summable allocation defect, not merely good leaf rounding.
4. **No optimizer-selection claim.**  If composition switches to a
   different spherical optimizer, its `(E,L)` state must be supplied or
   certified separately.  The theorem does not turn one exposed witness
   into a carrier for the full response roof.
5. **Relation to synchronization.**  The close-pole Walsh family in BS.4
   has one exposed direction with vanishing flatness.  FC.1 identifies the
   additional condition needed to preserve that property when such systems
   are assembled: their RMS amplitudes must also synchronize.

## 7. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_exposed_flatness_composition.py
```

The verifier checks the exact one-level and tree identities on random block
vectors, the hierarchical recovery inequality, and the pumpable tensor
formulas through depth 12.
