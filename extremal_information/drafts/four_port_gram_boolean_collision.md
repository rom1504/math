# A four-port equal-Gram Boolean support collision

**Status.** Rigorous exact benchmark counterexample with finite verifier.
It proves that a port Gram matrix does not determine the exact Boolean joint
linear response.  The construction has no quadratic child, so it is a
falsifier for Gram-only Boolean decoding in general, not yet an equal
Gram--Rayleigh collision inside the regular-Hadamard family.

## 1. Two orthogonal four-port systems

Fix `L>=1` and put `n=16L`.  Let `U` have one row equal to every sign pattern
in `{+-1}^4`, each repeated `L` times.  Let `V` have one row equal to every
even-parity sign pattern, each repeated `2L` times.  Regard their four
columns as labelled Boolean ports.

### Theorem GB.1 (equal Gram, separated Boolean support)

Both port systems have the same normalized Gram matrix,

```math
{1\over n}U^TU={1\over n}V^TV=I_4.                 \tag{GB.1}
```

Nevertheless their joint Boolean field responses

```math
\mathcal R(W)
=\max_{x\in\{+-1\}^n}\sum_{i=1}^4|w_i^Tx|          \tag{GB.2}
```

are

```math
\boxed{\mathcal R(U)={3n\over2},\qquad
       \mathcal R(V)=2n.}                           \tag{GB.3}
```

Thus equal exact Gram data can hide a response gap `n/2` already at four
ports.

#### Proof

Both row distributions have unbiased coordinates and zero pair
correlations, proving (GB.1).  Boolean duality gives

```math
\mathcal R(W)
=\max_{\epsilon\in\{+-1\}^4}
 ||W\epsilon||_1.                                  \tag{GB.4}
```

For the uniform row distribution, multiplication by `epsilon` preserves the
distribution and

```math
\mathbb E|S_1+S_2+S_3+S_4|={3\over2}.              \tag{GB.5}
```

For the even-parity distribution this expectation is one when `epsilon`
has even parity and two when it has odd parity.  Maximization therefore
gives (GB.3). `square`

Three ports cannot realize this moment-counting mechanism.  Modulo a global
row sign there are only four projective row patterns, and their total mass
plus three pair correlations determine all four frequencies.  At four
ports there are eight projective patterns but only seven such constraints;
the unresolved fourth-order parity mass is exactly what (GB.3) exposes.

## 2. Theory consequence and scope

The spherical response of both systems is identical because it depends only
on (GB.1):

```math
\max_{||x||_2^2=n}\sum_i|w_i^Tx|
=\sqrt n\max_\epsilon||W\epsilon||_2=2n.           \tag{GB.6}
```

For `V` this relaxation is exact; for `U` it has gap `n/2`.  Therefore there
is no universal decoder from a Gram matrix to exact Boolean joint support.
The missing coordinate is a fourth-order row-pattern statistic.

This does not yet prove information insufficiency of the full `(G,R)` state
for regular-Hadamard trust responses.  A common sign-matrix quadratic term
can distinguish the two row-pattern systems, and realizing both inside one
Hadamard spectral geometry is an additional algebraic obligation.  The
example is the mandatory finite linear benchmark for any proposed
Gram-only Boolean recovery theorem.

## 3. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_four_port_gram_boolean_collision.py
```

The verifier constructs both port systems for several multiplicities,
checks their exact Gram matrices, exhausts all endpoint words and old Boolean
spins at the base order, and verifies the response formulas by duality at
larger orders.
