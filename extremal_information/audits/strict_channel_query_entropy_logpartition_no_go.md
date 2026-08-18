# Query entropy defeats strict row-channel contraction

Status: **rigorous scalable no-go theorem with precise scope**.  A strict
contraction in every row channel, global bridge-sign symmetry, physical
`O(N^{-1/2})` bit Lipschitzness, and even exact zero-field Gibbs
log-partition structure do not force sublinear range over the `2^m` query
words.  There is an explicit family with all of these properties and
`Theta(N)` response range.  Under the natural channel marginal its pressure
variance is also `Theta(N)`, so the sublinear-MGF hypothesis in SC.15 fails
at the sharp scale.

The construction is a frozen rank-one two-block Gibbs system, not an actual
contracted-temperature minimizing child.  Its role is solely adversarial:
it proves that channel contraction plus generic physical structure cannot
close the bounded-row-degree cross-row problem.  The analogous rank-one
ferromagnetic sign children are quantitatively excluded by actual
minimality.  No claim is made that the construction is an actual-child
counterexample.

## 1. The exact channel and pressure

Let `m+n=N`, fix words

```math
 u\in\{\pm1\}^m,\qquad y\in\{\pm1\}^n,
```

and, for a bridge row `b in {+-1}^n`, put

```math
 z(b)={\langle y,b\rangle\over\sqrt n},\qquad
 q_v(b)={1\over2}\{1+v z(b)\}^2\quad(v\in\{\pm1\}).          \tag{QE.1}
```

Thus `q_v=e+vz`, where `e=(1+z^2)/2`.  As in SC.1--SC.7, if

```math
 d\mu=e\,dU_n,\qquad a={z\over e},
```

then `q_vdU_n=(1+va)dmu`, and

```math
 \kappa_n^2=E_\mu a^2
 =E_{U_n}{2z^2\over1+z^2}<1.                                  \tag{QE.2}
```

Moreover `kappa_n^2` tends to

```math
 \kappa_\infty^2
 =E_{G\sim N(0,1)}{2G^2\over1+G^2}
 =0.6886409151\ldots<1.                                       \tag{QE.3}
```

For a query word `v in {+-1}^m`, let `P_v` be the product of the row
laws `q_(v_i)dU_n`.  Define

```math
 S(B)=u^{\mathsf T}By,
 \qquad
 L_t(B)=\log\cosh\{tS(B)\},
 \qquad
 R_t(v)=E_{P_v}L_t(B).                                         \tag{QE.4}
```

This is a literal zero-field Gibbs log partition.  Namely, constrain the
left and right spins to the two-point sets `{+-u}` and `{+-y}` and give the
four pairs uniform mass.  Then

```math
 \log E_{x\in\{\pm u\},\,r\in\{\pm y\}}
       \exp\{t x^{\mathsf T}Br\}
 =\log\cosh\{t u^{\mathsf T}By\}=L_t(B).                       \tag{QE.5}
```

It has the two generic physical properties used by the spiked-response
audit:

```math
 L_t(-B)=L_t(B),                                                \tag{QE.6}
```

and flipping one bridge bit changes `L_t` by at most `2t`, because
`log cosh` is one-Lipschitz.  Consequently `R_t(-v)=R_t(v)` and all odd
query-Walsh coefficients vanish.

## 2. A pumpable linear response range

**Theorem QE.1 (strict contraction with extensive rare-query response).**
Suppose `m/N -> alpha in (0,1)`, `n/N -> 1-alpha`, and set

```math
 t={\beta\over\sqrt N},\qquad \beta>0.                         \tag{QE.7}
```

Then the response (QE.4) satisfies

```math
 \boxed{
 \operatorname {range}_{v\in\{\pm1\}^m}R_t(v)
 \ge {\beta m\sqrt n\over\sqrt N}
      -{\beta\sqrt{m(n-1)+n}\over\sqrt N}
      -\log2
 =\beta\alpha\sqrt{1-\alpha}\,N-O_\beta(\sqrt N).}           \tag{QE.8}
```

In fact the two displayed witness words have the sharp limits

```math
 {R_t(u)\over N}\longrightarrow
 \beta\alpha\sqrt{1-\alpha},
 \qquad
 {R_t(w)\over N}\longrightarrow0.                            \tag{QE.8a}
```

At the same time the exact weighted-Parseval inequality holds:

```math
 \sum_{S\ne\varnothing}\kappa_n^{-2|S|}\widehat R_t(S)^2
 \le \operatorname {Var}_{\mu^{\otimes m}}L_t.                \tag{QE.9}
```

Thus every individual row-order-`k` channel is attenuated by
`kappa_n^k`, uniformly away from one, while the exponentially large query
family still contains an extensive excursion.

*Proof.*  Under one row law in (QE.1), symmetry and the elementary moments
of a normalized Rademacher sum give

```math
 E_{q_vU_n}z=v,
 \qquad
 E_{q_vU_n}z^2={E z^2+E z^4\over2}=2-{1\over n},
 \qquad
 \operatorname {Var}_{q_vU_n}z=1-{1\over n}.                  \tag{QE.10}
```

Since `S=sqrt(n) sum_i u_i z(B_i)`, row independence yields

```math
 E_{P_v}S=\sqrt n\sum_i u_iv_i,
 \qquad
 \operatorname {Var}_{P_v}S=m(n-1).                            \tag{QE.11}
```

For the aligned query `v=u`, the inequalities
`log cosh r >= |r|-log 2` and `E|S|>=|ES|` give

```math
 R_t(u)\ge tm\sqrt n-\log2.                                   \tag{QE.12}
```

Choose `w` with `|sum_i u_iw_i|<=1`.  Since
`log cosh r<=|r|`, Cauchy--Schwarz and (QE.11) give

```math
 R_t(w)
 \le t E_{P_w}|S|
 \le t\sqrt{m(n-1)+n}.                                        \tag{QE.13}
```

Subtracting proves (QE.8).  Equation (QE.9) is the exact conditional-
expectation/Bessel identity SC.11--SC.12 applied to this `L_t`; it does not
use any special feature of `L_t`.  Finally,

```math
 E_{P_u}|S|\le m\sqrt n+\sqrt{m(n-1)},
```

so `R_t(u)<=tm sqrt(n)+O(sqrt(N))`; (QE.12) gives the matching lower
bound.  Equation (QE.13) gives `R_t(w)=O(sqrt(N))`, while `R_t(w)>=0`.
This proves (QE.8a). `square`

The construction is pumpable in the literal sense: taking any fixed
`alpha` and increasing `N` repeats the same one-row binary channel, while
the number of independent query inputs grows linearly and the favorable
word aligns all of them.  Pairwise contraction suppresses a prescribed
high-order Fourier channel; it does not pay for the entropy of choosing the
rare aligned word.

## 3. The natural pressure variance is sharply linear

The example also shows that the missing superconcentration condition cannot
be obtained from the generic log-partition geometry.

**Theorem QE.2 (no sublinear natural-MGF proxy).**  Under the hypotheses of
Theorem QE.1,

```math
 {1\over N}\operatorname {Var}_{\mu^{\otimes m}}L_t
 \longrightarrow
 2\beta^2\alpha(1-\alpha)\left(1-{2\over\pi}\right)>0.        \tag{QE.14}
```

In particular, any number `sigma_N^2` satisfying

```math
 \log E_{\mu^{\otimes m}}
 \exp\{s(L_t-EL_t)\}\le {s^2\sigma_N^2\over2}
 \quad(s\in\mathbb R)                                         \tag{QE.15}
```

must obey `sigma_N^2=Omega(N)`.

*Proof.*  Under `mu`, one row variable `z` is symmetric and

```math
 E_\mu z^2=2-{1\over n}.                                      \tag{QE.16}
```

Its fourth moment is uniformly bounded, because it is one half the sum of
the fourth and sixth moments of a normalized Rademacher sum.  The triangular-
array Lindeberg theorem therefore gives

```math
 {tS\over\sqrt N}
 ={\beta S\over N}
 \Longrightarrow
 G\sim N\bigl(0,2\beta^2\alpha(1-\alpha)\bigr),               \tag{QE.17}
```

with convergence of second moments.  Since

```math
 0\le |r|-\log\cosh r\le\log2,
```

`L_t/sqrt(N)` converges in distribution and second moment to `|G|`.
The variance of the absolute value of a centered Gaussian with variance
`s_G^2` is `s_G^2(1-2/pi)`, proving (QE.14).  Differentiating (QE.15) twice
at zero gives `Var(L_t)<=sigma_N^2`. `square`

Equivalently, the entropy union bound in SC.16 is scale-sharp on this
family: `2^m` strict-noise queries convert `Theta(N)` natural variance into
a `Theta(N)` rare minimum/range.

## 4. Exactly what actual child optimality excludes

The constrained spins in (QE.5) are frozen rank-one children.  The closest
ordinary hollow sign-matrix realization is the switched ferromagnet

```math
 A_{ij}=u_iu_j\quad(i\ne j).                                  \tag{QE.18}
```

It is not a contracted-temperature minimizing child.  If

```math
 p_A(t)=\log E_x\cosh\{tH_A(x)\},
```

then retaining `x=u` gives

```math
 p_A(t)\ge t{m\choose2}-m\log2.                               \tag{QE.19}
```

On the other hand, averaging the unlogged partition function over all edge
signs shows that an actual minimizing signing `A_*` satisfies

```math
 p_{A_*}(t)\le {m\choose2}\log\cosh t.                         \tag{QE.20}
```

At `t=beta/sqrt(N)` and `m asymp N`, (QE.19) is `Theta(N^(3/2))`, whereas
(QE.20) is `O_beta(N)`.  Thus the scalar optimizer pressure already excludes
the frozen ferromagnetic mechanism by a macroscopic margin.

This exclusion is narrower than the desired response theorem.  It does not
show that an actual minimizing child has sublinear external-disorder
variance, nor that exponentially many delocalized switching representatives
cannot collectively imitate (QE.8).  The exact optimizer identities in the
repository are zero-field or gauge-invariant identities; the spiked query
is an externally gauge-fixed observation.  SQ.2 additionally removes the
quartic/sector--Gram tangent as the source of an extensive response, but it
does not control the resummed sixth-and-higher response.

Therefore Theorems QE.1--QE.2 isolate the genuinely optimizer-specific
missing input:

```math
 \boxed{
 \text{actual-child external-disorder superconcentration,
 or an equivalent rare-query synchronization theorem}.}      \tag{QE.21}
```

Strict row-channel contraction, global evenness, physical bit Lipschitzness,
and generic Gibbs pressure structure cannot substitute for (QE.21).

## 5. Scope decision

This is a scalable general no-go, not a collision for actual children.  It
rejects any proposed proof that combines only the exact channel coefficient
`kappa<1`, a generic bounded-difference or log-partition estimate, and a
union bound over query words.  It leaves open precisely whether contracted-
temperature minimality supplies the extra superconcentration or
synchronization needed to beat the constructed entropy mechanism.
