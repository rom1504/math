# Audit of the Wave 40 abstract harmonic wall

## Verdict

The construction and all four quantitative claims are correct:

1. `V_s=(2v_s,v_s,v_s,0)`;
2. the weighted two-level isotonic residual is
   `w_s(1-w_s)(4v_s/3)^2`;
3. `K_L/sqrt(L)` converges to the displayed positive finite constant; and
4. endpoint conditional KL and normalized parent entropy converge to
   `2 log(2)/3` and `log(4/3)`.

The supplied checker passes.  Its curvature and residual formulas are coded
in already-reduced form rather than independently reconstructed from the four
states, but the direct audit below verifies those reductions.

Two wording corrections are advisable.  First, identify the normalized
parent likelihood explicitly as `f_L=e^(g_L)/Z_L` (equivalently subtract
`log Z_L` from `g_L`).  Second, “nor one using only endpoint cost and parent
entropy” should be scoped to a right side which remains bounded when those
two quantities remain bounded, such as `O(C_V(1)+H)`, or to a locally bounded
function of the endpoint pair.  An arbitrary pathological function of the
exact endpoint values could recover `L` from their exponentially small
deviations.  The example also shows only that **additional target structure**
is necessary; it does not logically establish that exact minimality is the
unique possible additional input.

## Direct checks

For either bit, conditioning on the other bit equal to zero leaves scores
`0,L`, and conditioning on it equal to one leaves scores `L,L`.  Hence the
two conditional curvatures are respectively

```math
v_s=L^2\frac{e^{sL}}{(1+e^{sL})^2},\qquad0.
```

Adding the two bit curvatures in state order `00,01,10,11` gives exactly
`(2v_s,v_s,v_s,0)`.  Conditional on score zero there is only state `00`;
conditional on score `L`, the three high states are equiprobable at every
time.  Therefore

```math
m_s(0)=2v_s,\qquad m_s(L)=\frac{v_s+v_s+0}{3}=\frac23v_s.
```

The nondecreasing two-point regression must pool this decreasing pair.  For
low-level mixture mass `w_s`, the weighted variance around the pooled mean is

```math
w_s(1-w_s)\left(2v_s-\frac23v_s\right)^2,
```

which is (W1).

With `u=sL`, write `h(u)=e^u/(1+e^u)^2`.  Then

```math
\frac{K_L}{\sqrt L}
=\frac43\int_0^L h(u)
\sqrt{u(1-u/L)w_{u/L}(1-w_{u/L})}\,du.
```

For fixed `u`, `w_(u/L)` tends to
`1/[2(1+3e^u)]`, yielding (W2).  Since `a_1<=a_(u/L)`, the integrand is
dominated by a constant times `sqrt(u)e^(-3u/2)` after a harmless adjustment
on a bounded interval.  This proves dominated convergence and positivity.
Independent quadrature gives

```math
c_*=0.12341114780933377\ldots,
\qquad c_*^2=0.01523031140361723\ldots,
```

consistent with the checker's increasing values of `K_L^2/L`.

For one endpoint bit, the active-context probability is

```math
\frac{1+e^L}{1+3e^L}\longrightarrow\frac13,
```

and its conditional KL tends to `log 2`.  Summing two bits proves (W3).
Finally, with `Z_L=(1+3e^L)/4` and normalized likelihood
`f_L=e^(g_L)/Z_L`,

```math
\operatorname{Ent}_\nu(f_L)
=\frac{3e^L}{1+3e^L}L-\log Z_L
\longrightarrow\log\frac43.
```

Thus `K_L^2` diverges linearly while both endpoint quantities stay bounded.
This rules out `K^2=O(C_V(1))`, `K^2=O(C_V(1)+H)`, and any endpoint-only
bound with a locally bounded right side near the limiting pair.

## Scope

The wall is an arbitrary two-bit exponential family with an unbounded score
gap.  It is not shown to arise from the matched selector/restriction
likelihood, an exact signing minimizer, restoring, or the orientation and
adjacent-selector package.  It therefore does not falsify (10.1096) for the
research problem.  It correctly falsifies only a structure-free attempt to
deduce the isotonic defect from endpoint conditional KL and parent entropy.
