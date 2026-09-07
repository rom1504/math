# Finite checker record

2026-09-07. Command:

```text
.venv/bin/python computations/principle_synthesis_2026_09_07_logarithmic_local_check.py
```

PASS. Every nonempty one-edge and two-distinct-edge Fourier coefficient was
exactly zero. All row chi-square identities and their rational upper bounds
passed. All full-law numerical KL values were below the proved rational
upper bound.

| Full sign input | n | k | Support size | Full KL (numerical) | KL upper (exact) | TV (numerical) |
|---|---:|---:|---:|---:|---:|---:|
| Symmetric H4 | 4 | 4 | 64 | 0.08301074146762076 | 1 | 0.1875 |
| Pentagon + I | 5 | 4 | 56 | 0.5270500293604091 | 36/25 | 0.395 |
| Arbitrary symmetric order 3 | 3 | 4 | 56 | 0.5961424232033783 | 44/27 | 0.430555555555555 |

Mean row chi-square at t=1,2,3, respectively:

- H4: 0, 1/4, 3/4.
- Pentagon + I: 0, 9/25, 27/25.
- Arbitrary order 3: 0, 11/27, 11/9.

At these row lengths every nonempty even subset has size two, so the row upper
inequality is equality. The final order-3 case deliberately has k>n: it tests
only the with-replacement statement (2), which has no k<=n restriction.

The initial run failed before executing mathematics because repository Python
3.9 has no `int.bit_count()`. It was replaced by `bin(value).count("1")` and
the complete checker rerun. This is an environment compatibility correction,
not a changed mathematical test or discarded failure.
