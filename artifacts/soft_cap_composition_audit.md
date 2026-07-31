# Scaled soft-cap composition: exact theorem and obstruction

Status: blank-slate derivation, subsequently compared with the fixed-pressure
candidate in ledger Section 10.113.2. The exact cosh formulation sharpens
that diagnostic but does not prove convergence.

## 1. Soft absolute cap and ground-state squeeze

For a symmetric signing `A` of order `n`, write

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,qquad
K(A)=\max_x|H_A(x)|,
```

and define the normalized cosh partition function

```math
Z_n(A;\gamma)=2^{-n}\sum_{x\in\{\pm1\}^n}
                  \cosh(\gamma H_A(x)).             \tag{1}
```

Let

```math
F_n(\gamma)=\min_A\log Z_n(A;\gamma),qquad
p_n(t)={1\over n}F_n(t/\sqrt n),qquad
g_n={M_n\over n^{3/2}}.                              \tag{2}
```

The two states `x,-x` at an absolute maximizer and the elementary bounds
`e^u/2<=cosh(u)<=e^{|u|}` give, uniformly in `A`,

```math
K(A)-{n\log2\over\gamma}
\le {1\over\gamma}\log Z_n(A;\gamma)
\le K(A).                                            \tag{3}
```

Minimizing preserves both inequalities. At the project scaling
`gamma=t/sqrt(n)`,

```math
g_n-\frac{\log2}{t}
\le \frac{p_n(t)}t
\le g_n.                                             \tag{4}
```

Consequently, if `p_n(t)` converges as `n` tends to infinity for every fixed
`t>0`, then `g_n` converges. Indeed, for every `t`,

```math
\liminf g_n\ge {p(t)\over t},\qquad
\limsup g_n\le {p(t)\over t}+{\log2\over t},         \tag{5}
```

and then `t` can tend to infinity. This needs no interchange of the limits
in `n` and `t`; (4) is uniform in `n`.

At a fixed unscaled inverse temperature `gamma>0`, the conclusion is even
more circular:

```math
{F_n(\gamma)\over\gamma n^{3/2}}-g_n\longrightarrow0. \tag{6}
```

Thus convergence of the fixed-`gamma` free-energy minima at the
`n^(3/2)` scale is equivalent to ground-state convergence, not a weaker
surrogate theorem.

## 2. Exact bridge composition and parameter rescaling

Take signings `A,B` of orders `m,n`, put `N=m+n`, and form a block signing
with bridge `C`:

```math
H_S(x,y)=H_A(x)+\epsilon H_B(y)+x^{\mathsf T}Cy,
\qquad \epsilon\in\{\pm1\}.                         \tag{7}
```

Choose the `mn` entries of `C` independently and uniformly, and also average
over `epsilon`. For each fixed `(x,y)`,

```math
\mathbb E_C\cosh\!\left(\gamma(H_A+\epsilon H_B+x^TCy)\right)
=(\cosh\gamma)^{mn}\cosh\!\left(\gamma(H_A+\epsilon H_B)\right).
```

The identity
`[cosh(a+b)+cosh(a-b)]/2=cosh(a)cosh(b)` then gives the exact factorization

```math
\mathbb E_{\epsilon,C}Z_N(S;\gamma)
=(\cosh\gamma)^{mn}Z_m(A;\gamma)Z_n(B;\gamma).       \tag{8}
```

Since a finite collection contains an element no larger than its average,
(8) proves, without exchanging a minimum and an expectation,

```math
F_N(\gamma)\le F_m(\gamma)+F_n(\gamma)
                    +mn\log\cosh\gamma.             \tag{9}
```

For a fixed total scaled parameter `t`, the same physical inverse
temperature `gamma=t/sqrt(N)` is the child scaled parameter

```math
t_m=t\sqrt{m/N},\qquad t_n=t\sqrt{n/N}.              \tag{10}
```

Writing `theta=m/N`, (9) is exactly

```math
p_N(t)\le
\theta p_m(t\sqrt\theta)
+(1-\theta)p_n(t\sqrt{1-\theta})
+{mn\over N}\log\cosh{t\over\sqrt N}.              \tag{11}
```

Equivalently, with `phi_n(t)=p_n(t)/t`,

```math
\phi_N(t)\le
\theta^{3/2}\phi_m(t\sqrt\theta)
+(1-\theta)^{3/2}\phi_n(t\sqrt{1-\theta})
+{mn\over tN}\log\cosh{t\over\sqrt N}.             \tag{12}
```

The bridge term tends to
`[t/2]theta(1-theta)` in (12), so it is leading order, not a summable
composition defect.

There is a useful exact renormalization. Put

```math
q_n(t)=p_n(t)-{t^2\over4}.                           \tag{13}
```

Using `log cosh u<=u^2/2`, the quadratic terms in (11) cancel exactly and
give

```math
q_N(t)\le
\theta q_m(t\sqrt\theta)
+(1-\theta)q_n(t\sqrt{1-\theta}).                   \tag{14}
```

This is defect-free, but it is not ordinary subadditivity: every balanced
split sends `t` to `t/sqrt(2)`. Iterating (14) drives the leaves to
infinite-temperature `t=0`, where `q_n(0)=0`. It recovers the annealed bound
`p_N(t)<=t^2/4` along recursively composable orders, rather than convergence
of `p_N(t)` at one fixed nonzero parameter.

## 3. Ordinary subadditivity subtracts the wrong scale

At fixed `gamma`, define

```math
G_n(\gamma)=F_n(\gamma)-{n\choose2}\log\cosh\gamma.  \tag{15}
```

Equation (9) says exactly

```math
G_{m+n}(\gamma)\le G_m(\gamma)+G_n(\gamma).          \tag{16}
```

Fekete's lemma therefore applies, but only to a quantity with asymptotics

```math
G_n(\gamma)
=-{n\choose2}\log\cosh\gamma+O_\gamma(n^{3/2}).     \tag{17}
```

In particular `G_n(gamma)/n` tends to minus infinity. The subadditivity is
dominated by its artificial order-`n^2` term and gives no control of the
order-`n^(3/2)` remainder. Restoring that remainder is exactly the
fixed-temperature problem (6), hence exactly the original convergence
problem at the required scale.

## 4. Exact quenched bridge interpolation

Jensen's annealed step can be refined without any illegal min/expectation
exchange. Reveal the bridge signs one at a time. Represent the cosh by an
auxiliary sign:

```math
\cosh(\gamma H)=\mathbb E_{\tau\in\{\pm1\}}e^{\gamma\tau H}.
```

Before revealing a bridge edge `e=(i,j)`, let `mu` be the current Gibbs law
on `(x,y,tau)` and set

```math
r_e=\mathbb E_\mu[\tau x_i y_j].                    \tag{18}
```

For its fresh sign `epsilon_e`, the exact partition-function ratio is

```math
{Z_{\epsilon_e}\over Z_0}
=\cosh\gamma+\epsilon_e r_e\sinh\gamma.
```

Therefore

```math
\mathbb E_{\epsilon_e}\log Z_{\epsilon_e}-\log Z_0
=\log\cosh\gamma
+{1\over2}\log(1-r_e^2\tanh^2\gamma).              \tag{19}
```

Summing (19), conditional on previously revealed signs, gives a genuine
quenched interpolation. Since `min_C log Z_C<=E_C log Z_C`, it can be used
constructively. It improves the annealed bridge only through the negative
correlation corrections in (19).

At `gamma=t/sqrt(N)`, the net leading bridge cost is heuristically and, after
Taylor control, quantitatively governed by

```math
{t^2\over2N}\sum_{e\text{ bridge}}(1-r_e^2).         \tag{20}
```

Thus making the balanced bridge contribution subextensive by this mechanism
requires `r_e^2=1-o(1)` for almost every edge at its reveal time. That is a
near-polarization theorem for the joint Gibbs state, not a consequence of
the parent free energies. No uniform such statement is available; at an
unpolarized reveal step `r_e=0`, (19) is exactly the full annealed cosh cost.

## 5. Resulting judgment after ledger comparison

Ledger Section 10.113.2 already identified the extensive bridge term and
temperature mismatch for the exponential absolute pressure. The cosh audit
adds three exact facts:

1. randomizing the orientation of one child gives the factorization (8)
   with no extra factor of two;
2. subtracting `t^2/4` gives the defect-free but parameter-contracting law
   (14);
3. the quenched reveal formula (19) identifies the precise missing datum as
   accumulated squared Gibbs correlations.

None closes convergence. Fixed unscaled temperature is asymptotically
equivalent to the ground problem; scaled temperature has the useful squeeze
(4), but its composition contracts the parameter. A successful continuation
would need either a uniqueness theorem for the entire parameter-rescaling
system (14), with matching lower control, or a uniform polarization/overlap
theorem making (20) subextensive. Proving either from soft free-energy values
alone would be circular, because those scalar values do not determine the
edge correlations in (18).
