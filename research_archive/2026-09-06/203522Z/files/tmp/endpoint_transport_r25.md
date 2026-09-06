### Wave 25, Route 3: a Hellinger endpoint reduction and the missing parent-cut theorem

**Status.** Sections 1--5 and the finite audits in Sections 7--8 are
**Verified** by `tmp/endpoint_transport_r25.py`.  The two estimates isolated
in Section 6 are **Proposed / open**.  In particular, this route does not yet
prove the restriction theorem.

#### 1. Lifted endpoint measures and exact local formulas

Retain the notation of (10.772)--(10.775).  Thus

```math
h_S(d)=e^{-\ell_S(d[S])},\qquad
a(d)=\mathbb E_{S\sim U_m}h_S(d),\qquad
C={Z_A(\beta)\over Z_0}.
\tag{R25.1}
```

Let

```math
q_S(d)=\mu_{\gamma,S}(d[S])
       \nu_\beta(d\mid d[S]),
\qquad
w_S={Z_S(\gamma)\over Z_0}.
\tag{R25.2}
```

Here every `q_S` is a normalized full-cut channel and
`\mathbb E_Uw_S=1`.  Direct cancellation of the child and conditional
partition functions gives the useful finite-measure identity

```math
\boxed{
w_Sq_S(d)=C\nu_\beta(d)h_S(d).
}
\tag{R25.3}
```

Indeed,

```math
q_S(d)={e^{\gamma c_S(d[S])}\over Z_S(\gamma)}
{e^{\beta(\langle A,d\rangle-c_S(d[S]))}\over K_{\beta,S}(d[S])},
```

and `\ell_S=\log K_{\beta,S}+(\beta-\gamma)c_S`.

There is also a completely local expression for an adjacent edge.  If
`S=C_0\cup\{u\}` and `T=C_0\cup\{v\}`, then (10.745)'s notation gives

```math
\boxed{
\ell_S-\ell_T
=\log{\nu_\beta(D_u\mid D_{C_0})
          \over\nu_\beta(D_v\mid D_{C_0})}
-\gamma(c_S-c_T).
}
\tag{R25.4}
```

This follows from
`G_S=\beta c_S+F_S=\log Z_A+\log\nu_{\beta,S}` and
`\ell_S=G_S-\gamma c_S`.  Equivalently, if
`r=h_T/h_S=e^{\ell_S-\ell_T}`, then the two endpoint edge integrands are

```math
h_S(r-1)\log r
\quad\hbox{and}\quad
h_S(\sqrt r-1)^2.
\tag{R25.5}
```

The first is the modified-LSI/Jeffreys integrand in (10.775); the second is
the ordinary-LSI/Hellinger integrand below.  Formula (R25.4) is exact, but no
available exact-minimizer inequality controls its exponentially tilted
average at the required scale.

#### 2. Finite-measure ordinary log-Sobolev inequality

Use throughout the convention

```math
H^2(\mu,\eta)
={1\over2}\sum_x(\sqrt{\mu(x)}-\sqrt{\eta(x)})^2
={\mu(1)+\eta(1)\over2}-\sum_x\sqrt{\mu(x)\eta(x)}.
\tag{R25.6}
```

This definition applies to finite measures of unequal mass.  For a positive
function `h` on the Johnson slice, put `a=\mathbb E_Uh` and
`\rho(S)=U(S)h(S)/a`.  Then

```math
aD_{\rm KL}(\rho\Vert U)=\operatorname{Ent}_U(h).
\tag{R25.7}
```

The ordinary Johnson log-Sobolev inequality, with precisely the
normalization of (10.728)--(10.729), says

```math
\operatorname{Ent}_U(h)
\le \tau_{\rm ls}{m(n-m)\over2n}
\mathbb E_{S\sim S'}(\sqrt{h_{S'}}-\sqrt{h_S})^2,
\tag{R25.8}
```

where

```math
\tau_{\rm ls}
\le {2\over\log2}\log{n^2\over m(n-m)}.
\tag{R25.9}
```

Apply (R25.8) with `h=h_\bullet(d)` for every fixed `d`, multiply by
`C\nu_\beta(d)`, and sum over `d`.  Since

```math
\mathbb E_{d\sim Ca\nu_\beta}D_{\rm KL}(\rho_d\Vert U)
=C\mathbb E_{d\sim\nu_\beta}\operatorname{Ent}_U(h_\bullet(d)),
```

one obtains the exact endpoint bound

```math
\boxed{
\mathbb E_{Ca\nu_\beta}D_{\rm KL}(\rho_d\Vert U)
\le \mathcal H_{\rm sel},
}
\tag{R25.10}
```

where

```math
\boxed{
\mathcal H_{\rm sel}
=\tau_{\rm ls}{m(n-m)\over2n}C
\mathbb E_{\nu_\beta,\,S\sim S'}
(\sqrt{h_{S'}}-\sqrt{h_S})^2.
}
\tag{R25.11}
```

There is no missing factor of two in (R25.11): (R25.6) contains the
`1/2`.  In terms of the finite measures in (R25.3),

```math
\boxed{
\mathcal H_{\rm sel}
=\tau_{\rm ls}{m(n-m)\over n}
\mathbb E_{S\sim S'}H^2(w_Sq_S,w_{S'}q_{S'}).
}
\tag{R25.12}
```

For clarity, its two sources of variation split exactly as

```math
\boxed{
H^2(w_Sq_S,w_Tq_T)
={1\over2}(\sqrt{w_S}-\sqrt{w_T})^2
+\sqrt{w_Sw_T}\,H^2(q_S,q_T).
}
\tag{R25.13}
```

Thus the local target simultaneously asks for smooth child partition
weights and Hellinger overlap of the adjacent lifted channels.  It is not
merely a conditional-channel estimate.

#### 3. Exact comparison with the exponential Jeffreys form

For every `x,y>0`,

```math
\boxed{
(x-y)(\log x-\log y)
\ge4(\sqrt x-\sqrt y)^2.
}
\tag{R25.14}
```

For example, set `t=\sqrt{x/y}\ge1`.  After cancelling the positive
factor `2(t-1)`, the claim is exactly
`(t+1)\log t\ge2(t-1)`, which follows from
`\log t\ge2(t-1)/(t+1)`.  Symmetry handles `x<y`.

Consequently, if the constants `\tau` are omitted but the common Johnson
coefficient is retained, the raw functionals obey

```math
\boxed{
\mathcal D_{\rm sel}^{\rm raw}
\ge4\mathcal H_{\rm sel}^{\rm raw}.
}
\tag{R25.15}
```

This is a pointwise statement, not just an averaged one.  It explains why
the square-root functional is much less sensitive to exponentially tiny
endpoint weights.  The certified entropy bounds are respectively
`\tau_{\rm mls}\mathcal D_{\rm sel}^{\rm raw}` and
`\tau_{\rm ls}\mathcal H_{\rm sel}^{\rm raw}`; their numerical ordering is
not implied by (R25.15), because their log-Sobolev constants differ.

#### 4. Exact endpoint sufficient lemma

Combining (10.774) with (R25.10) proves the following without any
approximation:

```math
\boxed{
D_{\rm KL}(M_0\Vert M_1)
\le \operatorname{Ent}_{\nu_\beta}(Ca)
+\tau_{\rm ls}{m(n-m)\over n}
\mathbb E_{S\sim S'}H^2(w_Sq_S,w_{S'}q_{S'}).
}
\tag{R25.16}
```

At a fixed selector density `m/n\to\rho\in(0,1)`, (R25.9) is `O_\rho(1)`.
It is therefore enough to prove the two separate estimates

```math
\boxed{
\begin{aligned}
\operatorname{Ent}_{\nu_\beta}(Ca)
&=O(n^{1/2-2c}),\\
\mathbb E_{S\sim S'}H^2(w_Sq_S,w_{S'}q_{S'})
&=O(n^{-1/2-2c}).
\end{aligned}}
\tag{R25.17}
```

Together they give `D_{\rm KL}(M_0\Vert M_1)=O(n^{1/2-2c})`, hence close
the soft full-interpolation alternative (10.773) at the same information
exponent required in (10.771).  They do not by themselves prove the newer
row-square-conditioned target (R24.8).

#### 5. What the local theorem does and does not control

The second line of (R25.17) is the sharp local theorem suggested by this
route.  It is strictly more tail-robust at the integrand level than (10.775),
and (R25.13) exposes its two components.  A proof would still need an
exact-minimizer competitor inequality that survives the endpoint weights;
(10.746)--(10.750) provide neither component at this scale.

More importantly, no selector-edge inequality alone controls the first line
of (R25.17).  The **minimal separate parent-cut theorem** is exactly

```math
\boxed{
\operatorname{Ent}_{\nu_\beta}
\left(
{Z_A(\beta)\over Z_0}
\mathbb E_{S\sim U_m}e^{-\ell_S(D[S])}
\right)
=O(n^{1/2-2c})
}
\tag{R25.18}
```

for exact minimizing signings, the target temperature regime, and fixed
retention density.  Calling (R25.18) a theorem is a specification, not a
proof.  A useful stronger formulation could bound this entropy by a parent
Gibbs flip/transport form, but the currently available conditional KL budget
(10.746) is under the wrong marginal and does not do so.

#### 6. Open theorem package and falsification criteria

**Proposed.** The endpoint route closes if both (R25.18) and the second line
of (R25.17) hold uniformly for exact minimizers at the parameters used by the
restriction argument.  It is falsified in this form by either of the
following asymptotic families:

1. exact minimizers for which
   `\operatorname{Ent}_{\nu_\beta}(Ca)/n^{1/2-2c}\to\infty`;
2. exact minimizers for which the average finite Hellinger edge in (R25.17)
   is `\omega(n^{-1/2-2c})`.

Finite examples can falsify proposed *generic comparison inequalities* but
cannot falsify these asymptotic scaling statements.  In particular, the
calculation below rules out absorbing the parent-cut term into the selector
Dirichlet form; it does not rule out a distinct exact-minimizer theorem such
as (R25.18).

#### 7. Exact-minimizer obstruction to absorption

Consider

```math
A_4=
\begin{pmatrix}
0&-1&1&1\\
-1&0&1&-1\\
1&1&0&-1\\
1&-1&-1&0
\end{pmatrix}.
```

Exhaustive enumeration gives `Q(A_4)=8=q_4`, so this is an exact minimizer.
For one deletion and `(\beta,\gamma)=(0.5,0.1)`, the exact-enumeration
values are

| quantity | value |
|---|---:|
| `\operatorname{Ent}_{\nu_\beta}(Ca)` | `0.627697022731` |
| actual conditional selector entropy | `0.147161088198` |
| `\mathcal D_{\rm sel}^{\rm raw}` | `0.363652934734` |
| `\mathcal H_{\rm sel}^{\rm raw}` | `0.080437866343` |
| Lee--Yau `\tau_{\rm ls}\mathcal H_{\rm sel}^{\rm raw}` | `0.388520927162` |
| total endpoint KL | `0.774858110929` |

Here `\tau_{\rm ls}=4.83007499856` is the explicit value from (R25.9).
Thus the cut entropy exceeds both the raw modified-LSI functional and the
full explicit ordinary-LSI certificate.  Even exact minimizers do not obey
either tempting absorption inequality.  A separate theorem of the form
(R25.18) is genuinely necessary for this decomposition.

#### 8. Low-temperature audit on `A_6,A_8,A_9`

For matched temperatures `\beta=\gamma`, one deletion, and the same
normalization as above, selected rows of the exhaustive audit are:

| matrix | `\beta` | cut entropy | selector entropy | Jeffreys raw | Hellinger raw | `\tau_{\rm ls}` times Hellinger raw |
|---|---:|---:|---:|---:|---:|---:|
| `A_6` | `0.5` | `0.008845` | `0.180890` | `0.637437` | `0.126814` | `0.722331` |
| `A_6` | `4` | `0.000000` | `0.182322` | `5.333334` | `0.166667` | `0.949332` |
| `A_6` | `8` | `0.000000` | `0.182322` | `10.666667` | `0.166667` | `0.949332` |
| `A_8` | `0.5` | `0.011454` | `0.412514` | `1.078315` | `0.227771` | `1.454383` |
| `A_8` | `4` | `0.000000` | `0.980827` | `11.999997` | `0.624664` | `3.988664` |
| `A_8` | `8` | `0.000000` | `0.980829` | `24.000000` | `0.625000` | `3.990806` |
| `A_9` | `0.5` | `0.012776` | `0.588663` | `1.689913` | `0.322284` | `2.152763` |
| `A_9` | `4` | `0.035687` | `1.355057` | `20.824471` | `0.733818` | `4.901685` |
| `A_9` | `8` | `0.035687` | `1.355062` | `42.157809` | `0.734127` | `4.903747` |

The verified numerical phenomenon is clear: from `\beta=4` to `8`, the raw
Jeffreys functional approximately doubles (`2.000000`, `2.000001`, and
`2.024436` respectively), while the raw Hellinger functional and the actual
selector entropy change by at most `0.054%` and `0.0004%`.  Thus the
Jeffreys form has a linear low-temperature blow-up on these examples, while
the Hellinger form and the entropy saturate.  This supports retaining the
ordinary-LSI reduction, but it supplies no asymptotic exact-minimizer bound.

#### 9. Disposition

This route yields a sharper exact endpoint formulation but not a closure.
The actionable theorem package is (R25.17): a local finite-measure
Hellinger estimate plus the independent parent-cut entropy theorem
(R25.18).  The `A_4` audit rigorously forbids replacing the latter by the
former.  Future work on this route should therefore either prove a parent
Gibbs transport inequality for the specific density `Ca`, or abandon the
endpoint chain rule in favor of a formulation in which the parent-cut
marginal is fixed by construction.
