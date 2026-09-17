# Wave 43: integrated context-mass control of harmonic migration

## Status

The statements through (R43.16) below are exact.  They give a global,
cross-coordinate sufficient norm for signed vertex migration and two
baseline-free context-mass majorants.  They do **not** prove the required
project estimate for exact minimizing signings.

The main conclusion is that the three-factor cost from (10.1144) is useful
locally but its raw magnitude cannot be the project target: it contains a
large constant omission baseline even when the likelihood is identically
one.  The intrinsic object is instead the variance of the **summed actual
binary KL load**.  The exact `A9` edge from (10.1146) is harmless for this
integrated formulation because its context mass is exponentially small,
whereas the abstract Wave-36 transient-mode example is still detected.

The numerical audits are in
`/home/math/quadra/tmp/harmonic_integrated_context_r43_check.py`.

## 1. Context normalization and the global covariance

Let `nu` be the positive base law, let `f>0` satisfy
`E_nu f=1`, put `g=log f`, and interpolate by

```math
\mu_s(d)=\frac{\nu(d)e^{s g(d)}}{\mathbb E_\nu e^{s g}},
\qquad 0\le s\le1.
```

Use the oriented-cut chart and let `V_*` be its `n-1` nonreference vertex
coordinates.  For `i in V_*`, a context `e=D_{-i}` is an unordered edge
`{x,y}`.  Define

```math
M_{s,i}(e)=\mu_s(x)+\mu_s(y),
\qquad
A_{i,e}(s)=\log\mathbb E_{\nu(D_i\mid e)}e^{s g},
```

```math
X_{i,e}(s)=A'_{i,e}(s)
=\mathbb E_{\mu_s}[g\mid D_{-i}=e],
```

and the conditional binary KL

```math
k_{i,e}(s)=sA'_{i,e}(s)-A_{i,e}(s)
=D(\mu_s(D_i\mid e)\Vert\nu(D_i\mid e)).
```

For every fixed coordinate, `M_(s,i)` is a probability law on contexts:
`sum_e M_(s,i)(e)=1`.  If `bar g_s=E_(mu_s)g`, direct differentiation gives

```math
\boxed{
M'_{s,i}(e)=M_{s,i}(e)\{X_{i,e}(s)-\bar g_s\}.
}
\tag{R43.1}
```

There is no missing context normalization.  In particular,

```math
\boxed{
\operatorname{Cov}_{M_{s,i}}(X_{i,e},k_{i,e})
=\sum_e M'_{s,i}(e)k_{i,e}(s).
}
\tag{R43.2}
```

Put the total vertex load on a state

```math
L_s(d)=\sum_{i\in V_*}k_{i,d_{-i}}(s).
```

Conditioning `g` on each context and summing **before** taking any absolute
value gives the already identified collapse, now with its normalization
explicit:

```math
\boxed{
\sum_{i\in V_*}
\operatorname{Cov}_{M_{s,i}}(X_{i,e},k_{i,e})
=\operatorname{Cov}_{\mu_s}(g,L_s).
}
\tag{R43.3}
```

Thus cross-coordinate cancellation is retained exactly.

## 2. The three-factor cost is finite and controls score variation

Assume finite temperature and `0<m<n`.  Every completion weight and every
selector posterior is then positive, and every coordinate has both
containing and omitting selectors.  Hence, on every vertex context edge,

```math
0<r_i(x),r_i(y)<1.
```

The matched exclusion identity and the cost in (10.1144) are

```math
\chi_{i,e}=g(y)-g(x)=\log\frac{r_i(x)}{r_i(y)},
\qquad
C_{i,e}=-\frac12\log\{r_i(x)r_i(y)\}.
```

Writing `u_x=-log r_i(x)>=0` and `u_y=-log r_i(y)>=0` proves, on **every**
vertex edge rather than only on a crossing edge,

```math
\boxed{
C_{i,e}=\frac{u_x+u_y}{2}
\ge\frac{|u_x-u_y|}{2}
=\frac{|\chi_{i,e}|}{2}.
}
\tag{R43.4}
```

This is consistent with the factorization

```math
C_{i,e}=C_{\rm lev}+C_{\rm aff}+C_{\rm rev}.
```

The total is nonnegative, as (R43.4) shows, while `C_aff` and `C_rev` are
nonnegative and `C_lev` can be negative.  No argument may discard that
negative cross-level term before reconstructing the total.

There are two elementary binary-KL bounds.  The tilted log-density ratio has
range `s|chi|`, while `A''(s)<=chi^2/4`; therefore

```math
\boxed{
0\le k_{i,e}(s)
\le\min\left\{s|\chi_{i,e}|,
\frac{s^2\chi_{i,e}^2}{8}\right\}
\le2sC_{i,e}.
}
\tag{R43.5}
```

The middle expression is baseline-free.  It is convenient to set

```math
w_{i,e}(s)
=\min\left\{|\chi_{i,e}|,\frac{s\chi_{i,e}^2}{8}\right\},
\qquad k_{i,e}(s)\le s w_{i,e}(s).
\tag{R43.6}
```

## 3. Intrinsic global load-variance bound

Let

```math
V_s=\operatorname{Var}_{\mu_s}(g),
\qquad
\mathscr H=\operatorname{Ent}_\nu(f)
=\int_0^1sV_s\,ds,
```

and define the intrinsic total-load norm

```math
\boxed{
\mathcal J_L^2
=\int_0^1\frac{\operatorname{Var}_{\mu_s}(L_s)}{s}\,ds.
}
\tag{R43.7}
```

At `s=0`, use the continuous value zero.  Indeed `k_(i,e)(s)=O(s^2)`, so
the integrand is `O(s^3)`.

For the globally signed adverse vertex migration

```math
\mathcal A_V
=\left[-\int_0^1\operatorname{Cov}_{\mu_s}(g,L_s)\,ds\right]_+,
```

Cauchy--Schwarz in state and then in time gives

```math
\boxed{
\mathcal A_V
\le\int_0^1\sqrt{V_s\operatorname{Var}_{\mu_s}(L_s)}\,ds
\le\sqrt{\mathscr H}\,\mathcal J_L.
}
\tag{R43.8}
```

This is one Cauchy step applied **after** summing the coordinate loads.  It
has no accidental `sqrt(n)` or `n` from coordinatewise Cauchy--Schwarz.

A concrete baseline-free score majorant follows from (R43.5).  Put

```math
W_s(d)=\sum_{i\in V_*}w_{i,d_{-i}}(s).
```

Since `0<=L_s(d)<=sW_s(d)`, one has

```math
\boxed{
\mathcal J_L^2
\le\int_0^1s\,\mathbb E_{\mu_s}W_s^2\,ds
\le\int_0^1s\,\mathbb E_{\mu_s}
\left(\sum_{i\in V_*}|\chi_{i,D_{-i}}|\right)^2ds.
}
\tag{R43.9}
```

The first quantity in (R43.9) vanishes for a flat likelihood and weights a
large score edge by its actual state mass.

There is also a genuinely signed context-mass version.  Define

```math
\Phi_s
=\sum_{i\in V_*}\sum_e
w_{i,e}(s)[-M'_{s,i}(e)]_+.
\tag{R43.10}
```

Equations (R43.2) and (R43.5), with `k>=0`, give

```math
\boxed{
\mathcal A_V
\le\int_0^1\sum_{i,e}
k_{i,e}(s)[-M'_{s,i}(e)]_+\,ds
\le\int_0^1s\Phi_s\,ds.
}
\tag{R43.11}
```

If `V_s=0`, then every `M'_(s,i)=0`; set the following quotient to zero.
Another time Cauchy inequality yields

```math
\boxed{
\mathcal A_V\le\sqrt{\mathscr H}\,\mathcal J_\Phi,
\qquad
\mathcal J_\Phi^2
=\int_0^1s\frac{\Phi_s^2}{V_s}\,ds.
}
\tag{R43.12}
```

This is the cost-weighted, direction-sensitive form: only contexts losing
mass are charged.  The looser substitution `w_(i,e)<=2C_(i,e)` is always
valid, but the next section shows why raw `C` should not define the main
project target.

## 4. The raw-cost baseline wall

Take the completely flat matched endpoint `f=1` at fixed selector density
`rho=m/n<1`.  The posterior selector law is uniform and

```math
r_i(d)=1-\rho,
\quad \chi_{i,e}=0,
\quad k_{i,e}(s)=0,
\quad L_s=0.
```

Thus there is no entropy, energy, or migration.  Nevertheless

```math
C_{i,e}=-\log(1-\rho)>0.
```

Consequently, with `n-1` chart vertices,

```math
\int_0^1s\,\mathbb E_{\mu_s}
\left(\sum_iC_{i,D_{-i}}\right)^2ds
=\frac12(n-1)^2\log^2\frac1{1-\rho}
=\Theta(n^2).
\tag{R43.13}
```

Demanding an `O(n^(1/2-2c))` bound on this raw square would fail on the
trivial endpoint.  This is not a counterexample to (R43.8), (R43.9), or
(R43.12), all of which are zero.  It shows that reveal, affinity, and
cross-level magnitudes must be baseline-corrected or used through actual
score variation and signed context movement.

## 5. Closed bootstrap, with orientation separate

Let

```math
I_V=\int_0^1s\mathcal E_{s,V}(g)\,ds,
\qquad
I_0=\int_0^1s\mathcal E_{s,0}(g)\,ds
```

be the vertex and global-orientation energies.  The vertex endpoint identity
is

```math
I_V=C_V(1)-\int_0^1\operatorname{Cov}_{\mu_s}(g,L_s)\,ds
\le C_V(1)+\mathcal A_V.
\tag{R43.14}
```

The orientation coordinate has its own exact identity

```math
I_0=C_0(1)-\int_0^1
\operatorname{Cov}_{M_{s,0}}(A'_{0,e},k_{0,e})\,ds,
\tag{R43.15}
```

but no selector-exclusion probability and no version of (R43.4).  It must be
bounded independently, either directly as `I_0=O(a_n)` or by separately
controlling its endpoint cost and adverse migration.

Set `a_n=n^(1/2-2c)`.  Suppose uniformly on the project family that

```math
C_V(1)\le K_Ea_n,
\qquad I_0\le K_Oa_n,
\qquad \mathcal J_L^2\le K_Ja_n,
```

and the restoring comparison is

```math
V_s\le R\{\mathcal E_{s,V}(g)+\mathcal E_{s,0}(g)\}.
```

Writing `h=sqrt(H)` and using (R43.8) gives the noncircular quadratic

```math
h^2\le R(K_E+K_O)a_n+R\sqrt{K_Ja_n}\,h.
```

Hence

```math
\boxed{
h\le
\frac{R\sqrt{K_J}+\sqrt{R^2K_J+4R(K_E+K_O)}}{2}\sqrt{a_n},
\qquad
\mathscr H=O(a_n).
}
\tag{R43.16}
```

The alternative hypothesis `J_Phi^2=O(a_n)` closes the same bootstrap by
(R43.12).  Endpoint vertex cost, restoring, orientation, and the separate
adjacent-selector Hellinger estimate remain independent inputs.

## 6. Exact and numerical stress tests

### The `A9,m=4` pointwise wall passes the integrated rarity test

On the edge from (10.1146), at `beta=2`, the exact factor costs are

```text
C_rev   =  6.6033335388
C_aff   =  6.7679672345
C_lev   = -7.0790096615
C_total =  6.2922911118
chi     = 11.3924395965
```

Thus `C_lev` is genuinely negative, but
`C_total>=|chi|/2` and the three terms reconstruct the total.  Direct
quadrature gives the following values for the same context:

| `beta` | maximum `M_s(e)` | `integral s M_s(e) C_e^2 ds` |
|---:|---:|---:|
| 1 | `3.768e-15` | `6.35e-15` |
| 2 | `2.490e-27` | `1.38e-26` |
| 4 | `1.190e-51` | `1.62e-50` |

The mass exponent between the displayed temperatures is numerically `28`,
matching (10.1147).  Large pointwise reveal/affinity cost on this edge is
therefore negligible after actual context-mass weighting.

For the full `A9,m=4,beta=2` endpoint, quadrature gives

```text
H                    = 1.6675784e-2
signed adverse A_V   = 8.2014154e-5
J_L^2                = 1.7199750e-5
sqrt(H J_L^2)        = 5.3555515e-4
```

The exact endpoint identity error is below `2.3e-15`.

### The abstract transient-mode wall remains visible

The strictly positive canonical Wave-36 example has fixed-size selector
omission and a separate orientation bit, but is not a quadratic-signing
endpoint.  At `M=10000`, the vertex quantities are

```text
endpoint vertex cost = 8.49450e-6
signed adverse A_V   = 1.24004e-3
J_L^2                = 3.03300e-6
sqrt(H J_L^2)        = 1.44841e-3
exact decreasing-context flux majorant = 1.24004e-3
```

Thus endpoint rarity alone misses migration by a factor about `146`, while
the intrinsic global norm and exact context flux see it.  Broad pointwise
majorants are deliberately much looser on this example.  This confirms that
the surviving theorem must be signing/minimizer-specific and integrated in
time; the `A9` rare-edge mechanism by itself is not enough.

## 7. Ledger-facing conclusion

The meaningful new sufficient harmonic target is

```math
\boxed{
\mathcal J_L^2
=\int_0^1\frac{\operatorname{Var}_{\mu_s}
\left(\sum_{i\in V_*}k_{i,D_{-i}}(s)\right)}{s}\,ds
=O(n^{1/2-2c}).
}
\tag{R43.17}
```

It preserves cancellation inside the total local load and across vertex
coordinates, is zero for the flat endpoint, passes the exact `A9` rare-edge
test, and rejects the known abstract transient-mode wall.  The
direction-sensitive alternative is `J_Phi^2=O(n^(1/2-2c))` from (R43.12).

What remains open is an exact-minimizer estimate for one of these two norms.
The reveal/affinity factorization proves the local domination (R43.4), but it
does not supply that estimate: raw factor costs have the baseline wall
(R43.13), and endpoint rarity alone has the Wave-36 transient wall.
