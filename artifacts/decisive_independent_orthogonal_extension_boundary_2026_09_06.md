# Orthogonal signing insertion: a calibrated boundary theorem

Status: elementary proof, pending independent audit. This is not a selectable
near-minimizer theorem, and its leading constant does not prove convergence.

Let C be a symmetric hollow conference signing, `C^2=lambda^2 I`,
`lambda=sqrt(n-1)`, and put `S=n lambda/2`, `M=M(C)`, `D=S-M>=0`.
Completing any row of C at its diagonal with either sign gives a valid new
sign row a. Then

```math
E(C)\le\lambda+1+
\begin{cases}
\lambda+D,&D\le\lambda,\\
2\sqrt{\lambda D},&D\ge\lambda.
\end{cases}                                                \tag{1}
```

Proof. For any cube x choose s so that `s q_C(x)=|q_C(x)|`. Orthogonality gives

```math
\|Cx-s\lambda x\|_2^2=4\lambda(S-|q_C(x)|).
```

For any completed row a, therefore,

```math
|a\cdot x|\le\lambda+1+2\sqrt\lambda\sqrt{S-|q_C(x)|}.
```

Writing `v=M-|q_C(x)|>=0`, maximize
`-v+2sqrt(lambda)*sqrt(D+v)` over v>=0. Its maximum is lambda+D when
D<=lambda, and `2sqrt(lambda D)` when D>=lambda. The permitted interval of v
can only reduce the maximum. This proves (1).

For a symmetric full Hadamard H with `H^2=nI` and `tr H=0`, the same proof
applies to the hollow signing `A=H-diag(H)`, since q_A=q_H. A full row of H
is already a valid new sign row, so the additive 1 disappears and
`lambda=sqrt(n)`. In particular exact spectral saturation gives
`E(A)<=2sqrt(n)`.

Consequences and limitations:

- Spectrally saturated symmetric Hadamard matrices are NOT examples with
  `E/sqrt(n)->infinity`, even if their Boolean ground layers are very large.
- Spectral deficit D=O(sqrt(n)) still gives E=O(sqrt(n)).
- A fixed fractional spectral deficit D=Theta(n^(3/2)) only yields O(n), not
  the desired small-row theorem. Thus the sub-half regime does not become
  easier under this argument.
- Isotropic absolute-ground laws give the opposite inequality E>=sqrt(n).
  This already exceeds `(3/2)c sqrt(n)` when c<1/2. A sharp selecting-row
  theorem must exploit quantitative non-isotropy of selected exact minimizers.

## Attempted scalable obstructions

Constant block replication of a signing multiplies its normalized cap by the
square root of the replication factor, leaving the sub-half range. Hadamard
tensor multiplication produces a vector-spin relaxation, not preservation of
the original scalar cap. Sparse zero-coupling models can have extensive row
discrepancy, but filling their missing edges at magnitude one introduces a
leading n^(3/2) objective; no sign-preserving reduction was found.

An intuitive encoding obstruction was identified, not proved as a theorem:
forcing r independent logical spins by ferromagnetic clusters against dense
cross fields appears to require cluster sizes at least r. Balanced extension
rows can then cancel within clusters except for O(1) per logical spin, giving
only r=O(sqrt(total dimension)) discrepancy. This does not exclude a more
efficient encoding and is not used in any proof.
