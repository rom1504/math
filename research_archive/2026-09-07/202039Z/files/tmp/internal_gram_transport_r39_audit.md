# Independent audit of `internal_gram_transport_r39.md`

## Verdict

The concentration, conditioning, entropy transport, joint-KL averaging, and
derandomization arguments pass.  The checker passes unchanged.  There is one
substantive scope correction: the new full-column estimate makes parent
deficit/cap persistence optional, so it must not appear in the claimed exact
remaining lemma or its falsification criterion.

## Concentration and scales

For Bernoulli selector indicators `xi_i` of mean `p`, with
`G=D_z A^2 D_z`,

```math
K_\xi=(p\mathbf1+\eta)^TG(p\mathbf1+\eta)
```

has the baseline and centered expansion stated in (R39.G5)--(R39.G6).
Bounded-variable Hanson--Wright applies to the centered quadratic part; the
linear subgaussian mgf applies to `2p<G1,eta>`.  Cauchy--Schwarz combines the
two without an independence assumption.  The spectral inputs are exactly

```math
\lVert G\rVert_F^2=\operatorname{tr}A^4\le2q_n n(n-1),
\qquad
\lVert G\mathbf1\rVert^2=z^TA^4z\le2q_nR_2(z).
```

Conditioning on `sum xi_i=m` costs at most `log(n+1)`, because `m` is a
binomial mode.  Entropy duality and the constrained `lambda` optimization
produce the square-root term plus `q_n H`; no pair-marginal transfer is used.

Conditioning (R39.G9) on `D` is valid.  The chain rule gives

```math
\mathbb E_DD(P_{S|D}\Vert U_m)
=D(P_S\Vert U_m)+I(S;D),
```

and Cauchy--Schwarz gives (R39.G1).  At the project scales, the square-root
term has exponent

```text
(3/2 + (9/4-c) + (3/4-c))/2 = 9/4-c,
```

while `q_n H` has the same exponent.  The `n^2` and `log n` contributions
are smaller for fixed `c<1/4`.

## Support and derandomization

In the high-ratio window `p_2>=1/2`, complement-incidence support makes each
restricted projective label favorable by (10.1045).  Sampling one kernel
certificate independently for each anchored selector preserves the expected
row, full Gram, and label conflict.  On the diagonal `S=T`, a deterministic
assignment has zero self-conflict, so its expected corrected conflict is no
larger than the kernel conflict.  Applying the probabilistic method to the
normalized nonnegative costs gives one deterministic assignment up to an
absolute factor.  This does not assume the desired collision; low conflict
remains an explicit hypothesis.

The fractional-cover identity also passes:

```math
\mathbb E_SD(\pi_S\Vert r)=\mathbb E_S\log(W/Z_S)\le\log W,
```

and its decomposition implies `I(S;D)<=log W`.

## Required correction: deficit is optional

Once (R39.G1) supplies

```math
\mathbb E\lVert A[:,S]z_S^D\rVert^2=O(n^{9/4-c}),
```

the external Parseval estimate and cap persistence are not needed for the
agreement route.  Complement-incidence support already supplies favorable
labels for every value of the parent deficit.  The actual sufficient kernel
package is therefore:

```math
D(P_S\Vert U_m)+I(S;D)=O(n^{3/4-c}),
\qquad
\mathbb E R_2(D)=O(n^{9/4-c}),
```

together with anchored density/affordability, complement-incidence support,
and project-scale anchored conflict.  A small `E Delta_D` may be recorded as
extra cap-persistence information, but it is not a required clause.

Accordingly:

- remove the deficit line from (R39.G3), (R39.G13), (R39.G16), and the
  surviving fractional-cover target, or mark it explicitly optional;
- do not say the third line is needed to prove the favorable Gram package;
- replace the falsifier “row high or deficit high.”  Large deficit alone does
  not obstruct this implementation.  A correct falsifier must rule out every
  affordable complement kernel having both project row/information and
  project conflict (or force the row above scale).

The statement should also explicitly attach favorability to the high-ratio
condition `p_2>=1/2`.  There is a harmless missing backslash in `lVert` in
(R39.G7).
