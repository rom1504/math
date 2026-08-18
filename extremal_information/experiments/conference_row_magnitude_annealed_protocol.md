# Exact annealed audit for row-magnitude conference fibres

**Classification.** Exact finite row/spin enumeration and exact annealed
bridge pressure.  This is a falsification experiment, not a quenched theorem
or asymptotic proof.

## Question

Fix a conference child `A_r`, a sign direction `v`, and a row event

```math
E=\{R:|\langle R,v\rangle|\in I\}.
```

Making all bridge rows independent and uniform on `E` produces a fibre of
relative mass `p_E^r=exp(-Theta(r))` whenever `p_E` is fixed.  Does even its
annealed pressure fall below the uniform-bridge value?

## Exact reduction

The conditioned row law is centrally symmetric.  With
`t=beta/sqrt(2r)` and

```math
M_v(y,t)=\mathbb E_{R\mid E}e^{t\langle R,y\rangle},
```

independence of the rows gives

```math
\mathbb E_B\cosh\{t(U+x^TBy)\}
=\cosh(tU)M_v(y,t)^r.
```

Thus `E_B Zbar` is computed by exhaustive sums over only one row and the
two child spin sets.  The program evaluates all `M_v(y,t)` by exact XOR
convolution using the Walsh--Hadamard transform; it never samples a bridge.

## Preregistered finite grid and outcome

The run exhausts every lower and upper absolute-row-sum threshold for:

- conference orders `r=6,10,14`;
- `beta=0.1,0.2,0.5`;
- both child orientations;
- the all-one direction, all distinct universal-double row directions,
  eight seeded random directions, and every projective direction at `r=6`.

There are respectively `1152`, `1080`, and `1848` comparisons.  No
conditioned-minus-uniform annealed difference is below `-10^(-8)` per row.
Restricting to row-event probabilities in `[0.1,0.9]`, the smallest
differences per row are

```text
r=6:  2.1509790698e-7
r=10: 6.6548221656e-8
r=14: 3.1332534844e-8
```

At `beta=0.5`, the all-one central event `|<R,1>|<=2` has total annealed
increments `0.00139061`, `0.00211452`, and `0.00248289` at orders
`6,10,14`; these are `O(1)` on the observed scale, not a negative linear
gain.  A raw minimum `-2.36e-9` per row at order fourteen occurs only in a
nearly full event and is below the stated floating tolerance.

The experiment therefore falsifies no part of the row-magnitude no-gain
hypothesis and supplied the right target for the subsequent rigorous
high-temperature coupling theorem.  It does **not** prove equality of
quenched pressures: `log E Z` and `E log Z` are different quantities.

## Reproduction

From the repository root:

```bash
.venv/bin/python \
  extremal_information/experiments/conference_row_magnitude_annealed_audit.py \
  --output \
  extremal_information/experiments/results/conference_row_magnitude_annealed_audit.json
```

Source SHA-256:
`246243fc3518ecbab83e6a17706e5970c4fca599f159c6030f61fa5127832d44`.
Result SHA-256:
`e93cbb9bed65c15d7b2f1b882db6bc98bb6cfafbd965f17cbb774b94f6bef537`.
No temporary directory is used.
