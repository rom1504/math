# Common-pole synchronization and its tensor algebra

**Status.** Rigorous task-local Boolean-recovery theorem and exact
composition law.  This is a sufficient condition for the regular-Hadamard
trust interface.  It does not claim that every recoverable response has a
common pole.

## 1. A one-witness recovery certificate

Let `H` be symmetric with

```math
H^2=r^2I,
\qquad \operatorname{tr}H=0,                       \tag{CS.1}
```

and let `w_1,...,w_p in {+-1}^n` be Boolean ports, each repeated `m`
times.  Define the Boolean and spherical responses

```math
\mathcal B=\max_{x,\sigma,\epsilon}
 \left\{{\sigma\over2}x^THx
       +m\sum_{i=1}^p\epsilon_iw_i^Tx\right\},      \tag{CS.2}
```

where `x in {+-1}^n`, and define `mathcal S` by replacing the cube with
`||u||_2^2=n`.  Suppose there is a Boolean top pole

```math
x_0\in\{+-1\}^n,
\qquad Hx_0=rx_0.                                   \tag{CS.3}
```

Its average port synchronization deficit is

```math
\delta(x_0;W)
=1-{1\over pn}\sum_{i=1}^p|w_i^Tx_0|.              \tag{CS.4}
```

### Theorem CS.1 (common-pole Boolean recovery)

With `c=mp/r`,

```math
\boxed{
0\le\mathcal S-\mathcal B
\le c\,\delta(x_0;W)\,rn.}                         \tag{CS.5}
```

Consequently, bounded total port mass and one top pole with
`delta=o(1)` give `mathcal S-mathcal B=o(rn)`.  If all auxiliary pairs are
then completed by an arbitrary exact signing `C` on `pm=O(r)` vertices,
the completed gap is at most

```math
c\delta rn+2Q(C)=c\delta rn+O(r^2),                 \tag{CS.6}
```

so the same conclusion holds when `r^2=O(n)` and `rn` is the leading scale.

#### Proof

The spectral and triangle inequalities give

```math
\mathcal S\le {rn\over2}+mpn.                      \tag{CS.7}
```

In the Boolean problem use `x=x_0`, the positive quadratic channel, and
choose every endpoint sign to match `w_i^Tx_0`.  Then

```math
\mathcal B\ge {rn\over2}
 +m\sum_i|w_i^Tx_0|
={rn\over2}+mpn(1-\delta).                         \tag{CS.8}
```

The cube is contained in the sphere, so subtraction proves (CS.5).
Adding the same auxiliary Hamiltonian changes either optimum by at most its
cap, which proves (CS.6). `square`

The certificate is strictly smaller than a Boolean response table: it is
one Boolean top pole and `p` correlations.  It is not necessary.  A response
may round through a channel-dependent witness even when no common pole sees
all ports.

## 2. Exact tensor composition

For `j=1,2`, let `(H_j,x_j,W_j)` satisfy (CS.1)--(CS.4), with orders `n_j`,
port counts `p_j`, and deficits `delta_j`.  Form

```math
H=H_1\otimes H_2,
\qquad x_0=x_1\otimes x_2,
\qquad
W=\{w\otimes v:w\in W_1,\ v\in W_2\},             \tag{CS.9}
```

where the Cartesian tensor family is counted with multiplicity if two
presentations yield the same vector.

### Theorem CS.2 (multiplicative synchronization quality)

The tensor family has

```math
\boxed{1-\delta=(1-\delta_1)(1-\delta_2),}
\qquad
\delta=\delta_1+\delta_2-\delta_1\delta_2.          \tag{CS.10}
```

For `L` factors,

```math
1-\delta_{[L]}=\prod_{j=1}^L(1-\delta_j).           \tag{CS.11}
```

Thus common-pole recovery is reusable under tensor composition exactly when
the response-weighted loss `c_[L] delta_[L]` vanishes at the required scale
(and the completion term is lower order), where `c_[L]` is the total port
mass of the composed response.  Deficit alone is sufficient only when
`c_[L]=O(1)`.  In particular,

```math
\delta_{[L]}\le\sum_{j=1}^L\delta_j,                \tag{CS.12}
```

whereas repeating any fixed positive deficit drives `delta_[L]` to one.

#### Proof

The pole in (CS.9) remains Boolean and has top eigenvalue `r_1r_2`.
Moreover

```math
{1\over p_1p_2n_1n_2}
\sum_{w,v}|(w\otimes v)^T(x_1\otimes x_2)|
=\left({1\over p_1n_1}\sum_w|w^Tx_1|\right)
 \left({1\over p_2n_2}\sum_v|v^Tx_2|\right).       \tag{CS.13}
```

This proves (CS.10), and iteration gives (CS.11)--(CS.12). `square`

The scalar `1-delta` is therefore a genuine, exact composition coordinate
for this restricted recovery certificate, but not a full contextual state:
equal deficits can have different Boolean responses.  It explains two opposite
examples.  Identical or `o(n)`-perturbed poles retain quality near one;
orthogonal poles have zero correlation with any chosen member and incur a
fixed loss.  It does not encode coefficient-side cycle holonomy or certify
families with optimizer switching.

## 3. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_common_pole_synchronization.py
```

The verifier exhausts finite Boolean responses for regular Walsh orders,
checks (CS.5) on random port families, verifies completion Lipschitzness,
and checks the tensor deficit identity exactly.
