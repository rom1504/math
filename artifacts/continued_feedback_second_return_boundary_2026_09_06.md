# Second actual feedback return: exact boundary and focused falsifier

Date: 2026-09-06. This note records a **conditional next-query criterion**
and a finite actual-signing test, not a proved second-iteration theorem.
The first marked-history theorem remains proved at its stated scope in
`continued_feedback_first_marked_history_energy_projection_2026_09_06.md`.

## 1. The actual next field, without deleting its coherent return

Use that theorem's notation, including literal
`V=b0 QS+b1 QD`, raw residual `Z=Br(G,Y)`,
`C=H(G,Y) psi(V+Z)`, conditional mean `c0`, and deterministic coefficient
`a`. Write the exact identity

```math
BC=L+\eta,\qquad
L=B(c_0+D_aZ),\qquad
\eta=B(C-c_0-D_aZ).                                      (1)
```

At a fixed polynomial approximation stage, the proof already separates
the centered-coefficient and higher-noise pieces of eta. Their surviving
original Boolean degrees are at least five, and every proper local tensor
cut is `n^(-1/2) polylog(n)`. The same proof gives local separation from
the OLD coherent primitives `(S,G,Y,QS,QD)` and the old noise `Z`.
It does **not** give separation from the NEW `L`: its polynomial degree
can exceed the degree of a noise component of eta.

Here is the exact remaining condition at each fixed polynomial stage.
Let `K_{i,p}` be eta's exact degree-p Walsh tensor, and `J_{i,q}` that of
L. For every `q>p`, ask for

```math
\frac1n\sum_i\left\|
 K_{i,p}\mathbin{\star_p}J_{i,q}\right\|_F^2\longrightarrow0.       (2)
```

The star includes the exact distinct-label restrictions, not an
unrestricted pre-Walsh source tensor. Proper eta contractions are already
small; (2) is the only new full-noise contraction in the exact Boolean
Stein derivative product. At equal degree the full contraction is a
scalar covariance and can be retained by Gaussian regression. This
criterion is falsifiable and is strictly more than another root-map bound.

For clarity, small proper cuts of K alone do NOT imply (2). A Gaussian
variable N and the coherent response `L=N^3` give the elementary abstract
obstruction: after linear regression, N and `N^3-3N` still have nonzero
mixed fourth moments. This is not asserted to arise from (1).

The exact cube calculation proving sufficiency of the mixed derivative
condition is in `continued_audit_boolean_stable_noise_stein_2026_09_06.md`,
Sections 2--5. With (2), regress L against the finite vector of eta's
homogeneous components, retaining their deterministic covariance and the
literal regression residual R. At fixed variance cutoffs, the same Stein
calculation makes the eta vector independent Gaussian noise relative to
`(R,S,G,Y,QS,QD,Z)`. Low eta-variance components may be left unregressed;
their contribution to the next field itself is small in averaged L2.
Passing from polynomial stages to bounded responses would still require
the ordered approximation and variance-truncation bookkeeping; this note
does not silently claim that passage.

## 2. Operational consequence if this criterion is proved

Assume feasible `|f|+H<=1` and `|psi|<=1`, and set

```math
u=f+C,\qquad J=1-|u|,\qquad
j_2=\frac1n\mathbb E\sum_iJ_i|(Bu)_i|.
```

These are actual fields and `u` is in the cube. Exactly,

```math
Bu=V+Z+L+\eta.                                             (3)
```

If the above comparison is established, write at the retained finite
stage `L=R+sum_p beta_p eta_p`, and let `v_next,i` be the Gaussian
variance of `sum_p(1+beta_{i,p}) eta_{i,p}`. The local comparison would give

```math
j_2=\frac1n\sum_i\mathbb E\left[
 J_i\,\mathbb E_N|V_i+Z_i+R_i+\sqrt{v_{\rm next,i}}N|\right]+o(1)
\ \ge\sqrt{2/\pi}\,\frac1n\sum_i\mathbb E[J_i]
                    \sqrt{v_{\rm next,i}}+o(1).              (4)
```

The first line retains the actual coherent bias; the inequality uses
that a centered Gaussian has minimum expected absolute shift at zero.
There is no proved positive lower bound on `v_next`: regression can
cancel the innovation. The two feasible endpoints
`+u+J sign(Bu)` and `-u+J sign(Bu)` give the usual exact certificate
`j2+|e(u)+e(J sign(Bu))|`. Thus (2) has a concrete next-query consequence,
but neither (2) nor a useful uniform improvement has been established.

## 3. Focused actual-signing falsification attempt

The independent three-stage pilot/test script
`computations/continued_feedback_second_return_probe_2026_09_06.py`
uses

```math
f=\sin((G+Y)/\sqrt2),\quad H=1-|f|,\quad\psi(t)=\tanh(2t).
```

The first pilot estimates a, the second estimates the per-root
regression `R=L-beta eta`, and the third tests mixed moments and the
actual gain in (3). Conditional standard errors exclude BOTH pilots and
48-node Gaussian quadrature error. These are falsifiers, not certificates.

| Actual family | n | normalized Cov(eta²,R²) | normalized eta fourth cumulant | actual j2 |
| --- | ---: | ---: | ---: | ---: |
| Steiner | 120 | .06541 (.03526) | 15.8017 (1.1538) | .103104 (.000311) |
| Steiner | 496 | .06228 (.00569) | 2.0104 (.1410) | .102362 (.000103) |
| Steiner | 2016 | .003838 (.000647) | .14812 (.00630) | .103675 (.000063) |
| hidden B³ return | 64 | -.06915 (.01824) | 6.3783 (.9547) | .152581 (.000501) |
| hidden B³ return | 256 | -.03658 (.00217) | .78570 (.04721) | .154093 (.000272) |
| hidden B³ return | 1024 | -.009418 (.000807) | .13311 (.00486) | .153493 (.000135) |

The square covariance is normalized by the pilot average of
`Var(eta_i) Var(R_i)`; the fourth cumulant by the average of
`Var(eta_i)^2`. Steiner uses 10000 samples per stage at n120,496 and 6000
at n2016; hidden-B³ uses 6000 per stage. The hidden family is the genuine
hollow signing construction in `continued_feedback_hidden_third_return_2026_09_06.md`.
No persistent normalized obstruction is visible in these sizes, but the
finite cumulants are nonzero and no asymptotic inference is justified.
The measured second-query gains are well below the preserved universal
bound and do not supply a new constant.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_feedback_second_return_probe_2026_09_06.py --kind steiner --sizes 4 5 --samples 10000
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_feedback_second_return_probe_2026_09_06.py --kind steiner --sizes 6 --samples 6000
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_feedback_second_return_probe_2026_09_06.py --kind hidden --sizes 16 64 256 --samples 6000
```

The focused test is now banked. Further feedback iteration is paused at
(2), rather than advertised as an unlimited old-history closure.
