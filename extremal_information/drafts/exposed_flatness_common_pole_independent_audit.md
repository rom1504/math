# Independent audit: exposed flatness and common-pole synchronization

**Verdict.**  The two main algebraic results are correct.  FC.1--FC.3 and
the tensor benchmark pass; CS.1 and the exact tensor deficit law pass.  The
drafts provide useful restricted certificates, but neither scalar summary is
a complete reusable contextual state.  Three wording repairs are needed to
keep the scope exact.

The independent verifier is
`experiments/verify_flatness_common_pole_independent_audit.py`.  It does not
import either primary verifier.  It checks an irregular partition tree with
a zero block, exposed rank-one landscapes with exact spherical and Boolean
optima, all one- and two-port multisets at regular Hadamard order four, exact
rational tensor identities, and explicit scope counterexamples.

## 1. Exposed-flatness composition law

### Exact chain identity: pass

Write `r_i=||u_i||_2/sqrt(n_i)` and `lambda_i=n_i/N`.  Under
`sum_i lambda_i r_i^2=1`,

```math
1-\sum_i\lambda_i r_i
=\frac12\sum_i\lambda_i(r_i-1)^2.
```

Since `||u_i||_1=n_i r_i(1-phi_i)`, this gives FC.1 exactly.  There
is no separately paid scalar channel: the allocation and transported local
terms are nonnegative pieces of one exact `l_1/l_2` identity.

At a nonzero internal node `v`, the child coefficient produced by one more
application is

```math
\frac{N_v}{N}\frac{r_v}{r_{root}}
\frac{N_w}{N_v}\frac{r_w}{r_v}
=\frac{N_w}{N}\frac{r_w}{r_{root}}=\omega_w.
```

Thus the transport weights telescope exactly, including on irregular trees.
At any depth,

```math
\sum_v\frac{N_v}{N}R_v^2\le 1,
\qquad
\sum_v\omega_v\le1,
```

where equality holds in the first display when the nodes at that depth still
partition the root.  Cauchy--Schwarz proves the second display.  The leaf
blocks also form a disjoint partition, so their total transport mass is at
most one.  Consequently the depth estimate FC.14 is valid.

### Repair FC-R1: state the zero-node convention

For an internal block with zero Euclidean mass, relative child RMS values are
formally `0/0`, so “apply FC.1 to its normalized restriction” is undefined.
This is harmless but should be explicit:

> If `||u_v||_2=0`, set `A_v=0` and all descendant normalized flatnesses
> arbitrarily.  Their transport weights are zero, so they make no contribution.

The independent verifier exercises this case.

### Recovery constants: pass

For `s=sgn(u)`, with arbitrary signs at zero coordinates,

```math
\|u-s\|_2^2=2N\phi(u),
\qquad
\|u+s\|_2\le2\sqrt N.
```

Therefore

```math
\left|\frac12u^TMu-\frac12s^TMs\right|
\le \Lambda N\sqrt{2\phi(u)},
```

and

```math
|h^T(u-s)|
\le\kappa\Lambda N\sqrt{2\phi(u)}.
```

Adding the exposure loss gives FC.17 with exactly the claimed constants.
The theorem should say explicitly that `Lambda>0` (the draft defines `kappa`
only in that case); the purely linear `Lambda=0` case needs a separate
normalization.

### Pumpable benchmark: pass

For `rho_+-` as in FC.20 and `s=(rho_++rho_-)/2`, tensorization gives

```math
\frac{\|u_D\|_2^2}{2^D}
=\left(\frac{\rho_+^2+\rho_-^2}{2}\right)^D=1,
\qquad
\frac{\|u_D\|_1}{2^D}=s^D.
```

The total transport mass at level `j` is `s^j`; hence the exact contribution
is `s^j(1-s)` and the sum is `1-s^D`.  The pure linear benchmark then has
the exact extensive gap claimed in FC.28.

### Research value and limitation

The state `(E,L)` is a genuine exact carrier for the **flatness of one supplied
exposed witness** and yields a pumpable stopping criterion.  It is not a
contextual response carrier: two vectors can have the same `(E,L)` and respond
differently to a fixed linear query.  The verifier records a four-coordinate
example.  Thus this result is more than a tautology, but its reusable content
is conditional on preserving or separately selecting the witness through
composition.  The draft already acknowledges this limitation accurately.

## 2. Common-pole synchronization algebra

### One-witness response bound: pass

From `H^2=r^2I`, symmetry gives `||H||_op=r`.  Therefore every spherical
candidate is at most

```math
\frac12rn+mpn.
```

The Boolean pole `x_0` gives

```math
\frac12rn+m\sum_i|w_i^Tx_0|
=\frac12rn+mpn(1-\delta).
```

The difference is exactly `mpn delta = c delta rn`.  This proves CS.1.
Adding the same bounded auxiliary Hamiltonian to two optima changes their
gap by at most twice its cap, proving CS.6.  When `pm=O(r)` and `r^2=O(n)`,
that completion term is `o(rn)`.

### Tensor deficit law: pass

Every tensor-port correlation factors before taking the absolute value:

```math
|(w\otimes v)^T(x_1\otimes x_2)|
=|w^Tx_1|\,|v^Tx_2|.
```

Averaging proves

```math
1-\delta=(1-\delta_1)(1-\delta_2)
```

exactly.  The pole and involution hypotheses also tensor exactly.

### Repair CS-R1: carry port mass as a separate recovery coordinate

The sentence “common-pole recovery is reusable ... exactly when the
accumulated synchronization loss vanishes” is too strong without a uniform
port-mass hypothesis.  CS.1 controls the normalized gap by

```math
c_{[L]}\delta_{[L]},
\qquad c_{[L]}=\frac{m_{[L]}p_{[L]}}{r_{[L]}}.
```

Thus the exact criterion supplied by the theorem is
`c_[L] delta_[L]=o(1)`, together with a negligible normalized completion
term.  If repetitions tensorize, then `c_[L]=prod_j c_j`; if instead a new
common multiplicity is imposed after tensorization, `c_[L]` must be recomputed.
The draft should specify which operation is intended.  Deficit quality alone
does not control recovery when port mass grows.

### Repair CS-R2: distinguish a compositional observable from a full state

The scalar `1-delta` is an exact multiplicative observable of a **presented
common-pole certificate**, not by itself a complete response state.  At order
four there are two two-port families with the same pole deficit `delta=1/2`
but Boolean responses 8 and 10.  In lexicographic cube order they are

```math
W_A=\{(-1,-1,-1,-1),(-1,-1,1,1)\},
```

and

```math
W_B=\{(-1,-1,-1,1),(-1,-1,-1,1)\}.
```

Accordingly, a reusable tensor certificate consists at least of the common
pole presentation, its tensor-compatible port family, synchronization quality,
and port mass.  The scalar quality is enough for the stated upper bound only
after those structural data are supplied.

## 3. Overall classification

| Claim | Audit | Mathematical status |
|---|---:|---|
| FC.1 one-level chain rule | pass | exact identity |
| FC.2 tree weights and depth bound | pass with zero-node convention | exact hierarchy |
| FC.3 recovery constants | pass with `Lambda>0` stated | rigorous sufficient theorem |
| pumpable rank-one benchmark | pass | scalable obstruction |
| CS.1 common-pole bound | pass | rigorous sufficient theorem |
| CS.2 tensor deficit law | pass | exact multiplicative observable |
| `(E,L)` as full response state | not claimed; false | explicit fixed-query counterexample |
| `delta` alone as full response state | false | same-deficit responses 8 and 10 |
| `delta_[L]->0` alone as recovery criterion | repair | require `c_[L]delta_[L]->0` |

The flatness law supplies a real but witness-conditional reusable state.  The
common-pole law supplies a real restricted synchronization algebra, provided
the certificate presentation and port-mass coordinate are retained.  Neither
result solves optimizer switching or general contextual compression, and the
drafts should not be promoted beyond that scope.

## 4. Reproduction

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_flatness_common_pole_independent_audit.py
```

