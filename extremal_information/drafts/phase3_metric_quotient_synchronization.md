# Metric-quotient synchronization of presented carriers

**Status.** The statements below are proved.  They complement the carrier-
capacity lower theorem with a sufficient condition for strict compression.
The condition is finite, deterministic, and stated only in the metric seen by
the endpoint queries.

Carrier Hausdorff entropy can force full response information, but raw carrier
size need not.  The two-scale counterexample suggests the correct upper
mechanism: a coarse metric quotient may collapse a large carrier while every
discarded fibre has only subscale diameter.  The theorem below makes that
mechanism exact and shows that it is stable under min-plus continuation.

## 1. Synchronizing metric quotients

Let `(X,d_X)` and `(Y,d_Y)` be finite metric spaces and let
`varpi:X->Y` be onto, with `a,b>=0`.  Call `varpi` an **`(a,b)` metric synchronization** if

1. `varpi` is one-Lipschitz;
2. every fibre has diameter at most `a`; and
3. for every `x in X` and `y in Y`, some `z in varpi^(-1)(y)` obeys

   ```math
   d_X(x,z)\le d_Y(\varpi x,y)+b.              \tag{MQ.1}
   ```

Condition (MQ.1) says that quotient distances can be lifted with additive
defect `b`.  An exact metric submetry has `b=0`; unlike an abstract linkage
assumption, all three clauses can be checked directly on finite objects.

Let `C subseteq X` be nonempty, let `p>=0`, let `alpha:C->[0,p]`, and define its presented
carrier response

```math
F_{C,\alpha}(x)=\min_{c\in C}
 \{d_X(x,c)+\alpha(c)\}.                       \tag{MQ.2}
```

The proposed quotient decoder is

```math
G_{\varpi C}(x)=d_Y(\varpi x,\varpi C).        \tag{MQ.3}
```

### Theorem MQ.1 (deterministic metric synchronization)

For every presented carrier,

```math
\boxed{
0\le F_{C,\alpha}(x)-G_{\varpi C}(x)
\le a+b+p\qquad(x\in X).}                     \tag{MQ.4}
```

Consequently the quotient carrier `varpi C`, rather than the full carrier or
its presentation, answers every endpoint query to uniform error `a+b+p`.

#### Proof

For every `c in C`, one-Lipschitzness gives

```math
d_Y(\varpi x,\varpi c)\le d_X(x,c)
\le d_X(x,c)+\alpha(c).
```

Taking minima proves the lower inequality in (MQ.4).  Conversely choose
`y in varpi C` nearest to `varpi x`, choose `c in C` with `varpi c=y`, and
use (MQ.1) to find `z` in the same fibre as `c` such that

```math
d_X(x,z)\le d_Y(\varpi x,y)+b.
```

The fibre bound gives `d_X(z,c)<=a`, so

```math
d_X(x,C)\le d_Y(\varpi x,\varpi C)+a+b.
```

Finally `F_(C,alpha)<=d_X(x,C)+p`, proving the upper bound. `square`

The fibre term and presentation term are both necessary.  For
`X=Y times Z` with the `ell_1` product metric, projection to `Y` is an
`(diam Z,0)` synchronization.  A singleton carrier at one end of a diameter
of `Z`, with presentation cost `p`, attains error `diam Z+p` at the opposite
end.

### Corollary MQ.2 (response complexity upper bound)

Let `mathcal C` be any family of presented carriers with common bounds
`a,b,p`.  If their quotient carriers have a Hausdorff `eta`-net of size `N`
in `Y`, then all responses have a uniform `(a+b+p+eta)`-net of size at most
`N`.  Equivalently,

```math
R_{a+b+p+\eta}^{\rm det}
\le\log_2\operatorname{Cov}
 \bigl(\{\varpi C:C\in\mathcal C\},d_H^Y,\eta\bigr).       \tag{MQ.5}
```

#### Proof

Distance-to-set functions on `Y` are isometric to the Hausdorff metric.
Choose a nearest quotient carrier from the net and combine that `eta` error
with (MQ.4). `square`

This is a strict quotient whenever the projected-carrier family has much
smaller metric entropy than the original carrier family.  Merely exhibiting
a projection is not enough: `a+b+p` must be below the declared response
scale.

## 2. Stability under future min-plus contexts

A min-plus context is any operator on real-valued profiles

```math
(T_Kf)(x)=\min_{z\in X}\{f(z)+K(x,z)\},         \tag{MQ.6}
```

where `K:X times X->R union {+infinity}` is fixed and every displayed output
is finite.  Using order inequalities and infima gives the corresponding
extended-real formulation without any indeterminate infinity differences.

### Proposition MQ.3 (no continuation amplification)

If `||f-g||_infinity<=epsilon`, then

```math
\|T_Kf-T_Kg\|_\infty\le\epsilon.              \tag{MQ.7}
```

The same is true after taking a maximum or minimum over the output.  Hence
the error `a+b+p` in MQ.1 survives any sequence of future min-plus contexts
without amplification, provided the quotient state is decoded once and is
not repeatedly rounded.

#### Proof

The inequalities `g-epsilon<=f<=g+epsilon` remain true after adding
`K(x,z)` and minimizing over `z`.  Iterate, and use the same monotonicity for
the terminal maximum or minimum. `square`

This is stronger than a one-query approximation, but it is not an automatic
closed algebra: the family of quotient carriers must itself update under the
declared composition.  The next proposition supplies one important case.

### Proposition MQ.4 (closed additive carrier algebra)

Assume `X,Y` are abelian groups and `varpi` is a homomorphism.  Compose two
presented carriers by

```math
C\star D=C+D,
```

with presentation cost

```math
(\alpha\square\beta)(w)
=\min_{c+d=w}\{\alpha(c)+\beta(d)\}.           \tag{MQ.8}
```

Then

```math
\varpi(C\star D)=\varpi C+\varpi D.           \tag{MQ.9}
```

Thus projected carriers form an exact compositional feature algebra.  If the
children have presentation radii `p_C,p_D`, the parent radius is at most
`p_C+p_D`.  The maintained approximate state is therefore the pair consisting
of the projected carrier and one scalar radius certificate.  After repeated
composition, the theorem guarantees submacroscopic error whenever the actual
or tracked presentation radius, together with `a+b`, is submacroscopic.  This
is a sufficient certificate, not a converse: infimal convolution can reduce
the true radius, and a large bound need not cause a large error.

#### Proof

Equation (MQ.9) is homomorphism of `varpi`.  Formula (MQ.8) is associative
infimal convolution, and every summand is at most `p_C+p_D`. `square`

The set algebra above needs no metric invariance.  If `d_X` is
translation-invariant, it also gives the exact response update

```math
F_{C+D,\alpha\square\beta}(x)
=\min_{d\in D}\{F_{C,\alpha}(x-d)+\beta(d)\}.   \tag{MQ.9a}
```

Translation invariance of `d_Y` gives the analogous interpretation for the
projected decoder.  Without it, (MQ.9) remains an exact carrier-set update but
need not be a min-plus identity for distance responses.

For linear-subspace carriers, the update is simply quotient-space join.
This is a genuine bounded or subexponential feature algebra whenever the
projected dimension stays small; it does not reconstruct the points in the
discarded metric fibres.

## 3. Three model tests

### Example MQ.5 (two-scale finite-field carrier)

Let `X=F_q^D`, let `varpi:X->F_q^r` be linear and onto, and put

```math
d_X(x,x')=s\mathbf1_{\varpi x\ne\varpi x'}
          +\mathbf1_{x\ne x'},
\qquad
d_Y(y,y')=(s+1)\mathbf1_{y\ne y'}.             \tag{MQ.10}
```

This is a `(1,0)` metric synchronization.  A target quotient fibre can be
reached from `x` at exactly the quotient distance by retaining its hidden
coordinate.  Therefore a mixed-holonomy carrier with circuit toll at most
`2k` is decoded from the subspace `varpi(im V)` with error at most `2k+1`.
For fixed `q,r` there are only finitely many such subspaces, although the exact
gluing map has `Dk` field coordinates and `diam(X)=s+1` can be linear in
`D`.

The specialized fragment argument is stronger for complete weighted kernel
alphabets: retaining the labeled map `varpi V` gives `q^(rk)` states and
all-future error one, because replacing every lifted letter changes the
aggregate kernel endpoint only inside one diameter-one fibre.  This sharper
statement is proved in
[`phase3_carrier_relation_response_law.md`](phase3_carrier_relation_response_law.md).

### Example MQ.6 (rank-metric row projection)

Let `X=M_D(F_q)` with rank distance, let `Y=M_{r,D}(F_q)` with rank
distance, and let `varpi` retain the first `r` rows.  The map is one-Lipschitz,
every fibre has diameter at most `D-r`, and a prescribed top block can be
reached while leaving the bottom rows unchanged.  Thus it is a
`(D-r,0)` metric synchronization.  Every rank-metric presented carrier is
decoded from its projected carrier with error

```math
D-r+p.                                         \tag{MQ.11}
```

When `D-r=o(D)` and `p=o(D)`, this is a non-Hamming submacroscopic factor.
It becomes a strict information quotient only for a carrier class whose
projected Hausdorff covering rate is demonstrably smaller.  It complements
the rank-metric capacity example: large separated projected carriers force
information, while geometry confined to the discarded rows does not.

### Example MQ.7 (punctured Hamming carriers)

Let `X=F_q^{r+h}` with Hamming distance and let `varpi` puncture the last
`h` coordinates.  It is an `(h,0)` metric synchronization.  Consequently
every code/coset-distance carrier with presentation radius `p` is recovered
from its punctured carrier to error `h+p`, stably under subsequent fixed
min-plus queries.  This is a small-error factor; strict information
compression additionally requires a smaller projected-carrier entropy, and
closed updating requires a composition, such as Minkowski addition, that
descends through puncturing.  Puncturing a fixed positive fraction does not
compress at vanishing normalized distortion.

These examples show that synchronization is neither field-specific nor a
claim about averages.  It is a uniform statement that hidden fibres have
small diameter in the exact query metric.

## 4. Boundary of the theorem

Small fibres alone do not suffice if quotient distances cannot be lifted;
small quotient state alone does not suffice if `a+b` is macroscopic; and a
small endpoint decoder does not imply closure under a different composition
that fails (MQ.9) or (MQ.7).  Conversely, MQ.1 does not require probabilistic
overlap identities, ultrametricity, or an assumed relation between hidden
labels.  It gives a finite, falsifiable route to deterministic
synchronization:

> find a response-metric submetry whose fibres have subscale diameter, then
> prove that the projected carriers close under composition.

Together with the carrier-capacity lower theorem, this gives a two-sided law
inside the presented-carrier class.  Hausdorff-rich quotient carriers force
information; low-entropy quotient carriers with small fibres compress it.
The word “compress” here includes the four strictness checks: subscale error,
smaller projected covering entropy, a descended composition law, and a
controlled presentation-radius certificate.
