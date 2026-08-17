# A biased boundary fill converts matched response roofs into scalar caps

**Status.** Task-local theorem draft.  The argument below is self-contained
apart from the standard net proof for the operator norm of a centred random
sign matrix.  It closes the *physical selector* step for the particular
matched-roof family constructed in Theorem 21.8.  It does not claim that an
arbitrary projectively separated collection of response functions can be
scalarized.

The main point is that the boundary does not need to be pinned exactly.
An exact-sign fill with rank-one bias `lambda/sqrt n` confines every relevant
optimizer to a small projective ball.  The bridge operator bound then makes
the response deficit stable inside that ball.  The same biased fill has a
positive `Theta(n^(3/2))` value at its pole, which prevents the outer
negative channel in the absolute cap from masking the desired positive
response.

The construction is query-owned in the usual contextual sense: for each
declared pole `q` one freezes one fill `D_q`, and that same fill is applied
to every child.  A finite (even exponentially large) query bank causes no
simultaneous-probability issue, since the existence statement can be applied
separately to each pole and the resulting fills then frozen.

## 1. Profiles and the matched roof

Let `C` be a hollow sign matrix on `n` old spins and let
`B in {+-1}^{n times n}`.  Put

```math
 H_C(x)=\sum_{i<j}C_{ij}x_ix_j,
 \qquad
 f_C^\sigma(y)=\max_x \sigma\{H_C(x)+x^TBy\},
 \quad \sigma\in\{+-1\}.                         \tag{BR.1}
```

Suppose `Q(C)<=P`.  Both signed profiles are even and have the common upper
roof

```math
 f_C^\sigma(y)\le U(y):=P+\|By\|_1.              \tag{BR.2}
```

If `\|B\|_(2 to2)<=L sqrt n` and `d=d_P(y,q)` is projective Hamming
distance, then

```math
 |f_C^\sigma(y)-f_C^\sigma(q)|
 \le2Ln\sqrt d,
 \qquad
 |U(y)-U(q)|\le2Ln\sqrt d.                       \tag{BR.3}
```

Indeed, replace `y` by `-y` if necessary and compare the two maxima.  The
difference of the linear fields is at most

```math
 \max_x|x^TB(y-q)|
 \le\sqrt n\,\|B\|\,2\sqrt d.
```

The special input needed below is one target child `C_*`, one pole `q`, and
one decoy `C`, satisfying

```math
 f_(C_*)^+(q)=U(q),
 \qquad
 f_C^+(q)\le U(q)-\delta n^{3/2}.                 \tag{BR.4}
```

This is a *directed matched-roof deficit*.  It is stronger than merely
knowing that two profiles have large projective distance, but it is exactly
the certificate produced in the proof of Theorem 21.8.

## 2. An exact-sign projective roof selector

For `q in {+-1}^n`, let `R_q` be the hollow rank-one signing

```math
 (R_q)_{ij}=q_iq_j\quad(i\ne j),
 \qquad
 H_(R_q)(y)={(q^Ty)^2-n\over2}.                   \tag{BR.5}
```

### Lemma BR.1 (biased exact-sign fill with spectral error)

There is an absolute constant `K_0` with the following property.  For every
fixed `lambda>0`, all sufficiently large `n`, and every `q`, there is a
hollow exact signing `D_q` such that, with `p=lambda/sqrt n`,

```math
 D_q=pR_q+E_q,
 \qquad
 \|E_q\|_(2 to2)\le K_0\sqrt n.                  \tag{BR.6}
```

Consequently

```math
 Q(D_q)\le{\lambda+K_0+o(1)\over2}n^{3/2}.       \tag{BR.7}
```

#### Proof

Choose the upper-triangular entries independently with

```math
 \mathbb P\{(D_q)_{ij}=q_iq_j\}={1+p\over2}.
```

Then `E_q=D_q-pR_q` is symmetric, centred, has independent bounded
upper-triangular entries, and a fixed-vector Hoeffding estimate followed by
a constant Euclidean net gives
`\|E_q\|<=K_0 sqrt n` with positive probability.  (The quadratic form on
one unit vector has a subgaussian tail because
`sum_(i<j)v_i^2v_j^2<=1/2`; a `1/4`-net has `9^n` points.)  This proves
(BR.6).  Since `\|R_q\|=n-1`,

```math
 Q(D_q)\le{n\over2}\|D_q\|
 \le{n\over2}\{p(n-1)+K_0\sqrt n\},
```

which is (BR.7). `square`

For later use, if `d=d_P(y,q)<=n/2` and `h=d(n-d)`, (BR.5)--(BR.6) give

```math
 H_(D_q)(q)-H_(D_q)(y)
 \ge2ph-2K_0\sqrt{nh},                            \tag{BR.8}
```

and

```math
 H_(D_q)(q)+H_(D_q)(y)
 \ge p\{n(n-1)-2h\}-K_0n^{3/2}
 \ge\left({\lambda\over2}-K_0-o(1)\right)n^{3/2}.
                                                               \tag{BR.9}
```

The first error estimate is sharper than a generic quadratic-form
Lipschitz bound.  After gauging `q=1`, the difference contains only the
`d by (n-d)` crossing block, so its absolute value is at most
`2\|E_q\|sqrt{d(n-d)}`.  For (BR.9), use
`h<=n^2/4` and `|H_(E_q)(z)|<=K_0n^(3/2)/2`.
Notice the asymmetry which makes the absolute channel manageable: the
rank-one mean is `+(lambda/2+o(1))n^(3/2)` at `q`, whereas its most negative
value over Boolean `y` is only `-O(lambda sqrt n)`.  Equivalently, the sum
of its values at `q` and any `y` has the leading lower bound in (BR.9).

## 3. Scalarization theorem

### Theorem BR.2 (bounded-cap matched-roof anti-pin)

Fix `L,delta,C_P>0`.  There is a constant
`C=C(L,delta,C_P)` such that the following holds for all sufficiently large
`n`.

Let `C_*` and `C` be hollow order-`n` signings with

```math
 Q(C_*),Q(C)\le P\le C_Pn^{3/2},
 \qquad \|B\|_(2 to2)\le L\sqrt n,                \tag{BR.10}
```

and suppose (BR.4) holds at `q`.  There is a hollow exact signing `D_q`,
depending on `q,L,delta` but **not** on which child is attached, for which
the complete exact-sign parents

```math
 \mathcal P_J^(q)=
 \begin{pmatrix}J&B\\B^T&D_q\end{pmatrix},
 \qquad J\in\{C_*,C\},                            \tag{BR.11}
```

have order `2n`, cap at most `Cn^(3/2)`, and

```math
 \boxed{
 Q(\mathcal P_(C_*)^(q))-Q(\mathcal P_C^(q))
 \ge {\delta\over2}n^{3/2}.}                     \tag{BR.12}
```

#### Proof

Take `D_q` from Lemma BR.1.  First choose
`theta in (0,1/4)` so small that

```math
 2(L+K_0)\sqrt\theta\le\delta/4,                 \tag{BR.13}
```

and then choose the fixed `lambda` so large that

```math
 2\lambda\theta(1-\theta)
 \ge2L+2K_0+\delta,
 \qquad
 {\lambda\over2}\ge2L+K_0+\delta+1.             \tag{BR.14}
```

Put

```math
 T=U(q)+H_(D_q)(q).
```

The matched equality in (BR.4) gives an actual positive-channel parent
configuration of value `T`, so

```math
 Q(\mathcal P_(C_*)^(q))\ge T.                   \tag{BR.15}
```

We bound both absolute-value channels of the decoy.  In the positive
channel, if `d=d_P(y,q)<=theta n`, (BR.3), the fact that the rank-one mean
in (BR.6) decreases away from `q`, and the spectral error estimate give

```math
 f_C^+(y)+H_(D_q)(y)
 \le T-\delta n^{3/2}
       +2Ln\sqrt d+2K_0\sqrt{nd(n-d)}
 \le T-{3\delta\over4}n^{3/2}.                   \tag{BR.16}
```

If instead `d>=theta n`, then (BR.2), (BR.3), and (BR.8) give

```math
 f_C^+(y)+H_(D_q)(y)
 \le T+{2L+2K_0-2\lambda\theta(1-\theta)}n^{3/2}
 \le T-\delta n^{3/2}.                           \tag{BR.17}
```

Here and below the harmless `o(1)` terms are absorbed by the spare `1` in
(BR.14).

For the negative channel, use `f_C^-(y)<=U(y)`, the coarse consequence
`U(y)-U(q)<=2Ln^(3/2)` of (BR.3), and (BR.9):

```math
 f_C^-(y)-H_(D_q)(y)
 \le T-{\lambda/2-K_0-2L-o(1)}n^{3/2}
 \le T-\delta n^{3/2}.                           \tag{BR.18}
```

Equations (BR.16)--(BR.18) show

```math
 Q(\mathcal P_C^(q))
 =\max_{\sigma,y}\{f_C^\sigma(y)+\sigma H_(D_q)(y)\}
 \le T-{3\delta\over4}n^{3/2}.                  \tag{BR.19}
```

Together with (BR.15) this proves the weaker displayed constant (BR.12).
Finally,

```math
 Q(\mathcal P_J^(q))
 \le Q(J)+\|B\|_(infinity to1)+Q(D_q)
 \le P+Ln^{3/2}+Q(D_q)=O(n^{3/2}),               \tag{BR.20}
```

which proves the cap claim. `square`

### Corollary BR.3 (a matched-roof code becomes a scalar contextual code)

Let `{C_i:i in I}` be a family satisfying (BR.10), and suppose there are
boundary words `{q_i:i in I}` such that

```math
 f_(C_i)^+(q_i)=U(q_i),
 \qquad
 f_(C_j)^+(q_i)\le U(q_i)-\delta n^{3/2}
 \quad(i\ne j).                                  \tag{BR.21}
```

For every `i`, choose the public query fill `D_(q_i)` supplied by BR.2.
Then the scalar complete-parent response table obeys

```math
 Q(\mathcal P_(C_i)^(q_i))-Q(\mathcal P_(C_j)^(q_i))
 \ge{\delta\over2}n^{3/2}
 \quad(i\ne j).                                  \tag{BR.22}
```

Thus the family is an `ell_infinity` scalar contextual packing of size
`|I|` under all-spins-free exact-sign contexts, with every parent of order
`2n` and cap `O(n^(3/2))`.

## 4. Application to Theorem 21.8

Orient a hollow sign quadratic globally so that

```math
 P=\max_xH_A(x)=Q(A).
```

(This is exactly the orientation used in the exact-minimizer `L_tail`
application.)  When Theorem 21.8 is applied to this `H=H_A`, its proof
supplies more than the displayed projective profile separation.  Every
switch `A^(s_y)` still has absolute cap `P`, as BR.2 requires.  With

```math
 C_y=A^(s_y),
 \qquad s_y=u_*\mathbin\odot\operatorname{sign}(By),
```

Lemma 21.7 gives the matched equality

```math
 f_(C_y)^+(y)=P+\|By\|_1,                        \tag{BR.23}
```

and the probabilistic event in that proof gives, for every distinct code
pair `y,z`, the two directed deficits

```math
 f_(C_z)^+(y)\le P+\|By\|_1-dn^{3/2},
 \qquad
 f_(C_y)^+(z)\le P+\|Bz\|_1-dn^{3/2}.           \tag{BR.24}
```

Its bridge already satisfies `\|B\|=O(sqrt n)`.  Corollary BR.3 therefore
turns all `exp(gamma n)` switched children into a scalar complete-parent
packing with gap `Omega(n^(3/2))`, exact signs, order `2n`, and cap
`O(n^(3/2))`.

Consequently, for a globally oriented bounded-cap sign quadratic with
`P=Q(A)`, the tail hypothesis (21.38) is sufficient for an `Omega(n)`-bit
**all-spins-free scalar physical contextual packing**.  In particular this
applies to the proposed exact-minimizer hypothesis `L_tail`.  The
boundary-selector clause previously left in that implication chain is no
longer a separate missing lemma.

## 5. Scope, information content, and archive comparison

1. **Not an arbitrary profile scalarizer.**  BR.2 uses a common upper roof,
   exact equality for the target at its named pole, a fixed directed roof
   deficit for the decoy, and a spectral Lipschitz bound on the bridge.  A
   bare projective response distance need not supply these data.
2. **Not a universal coordinate pin.**  The biased fill localizes only the
   profiles generated through one `O(sqrt n)`-operator bridge.  It does not
   make `q` optimal against every complete sign child.  Hence it does not
   satisfy the hypothesis of the universal-pin barrier UP.1.  Exact pinning
   is replaced by localization plus stability of the directed deficit.
3. **No full optimization is reconstructed.**  The query stores `q`, one
   scalar bias `lambda/sqrt n`, and a spectrally controlled exact-sign
   error.  It neither stores a child energy table nor computes its optimizer.
   The target pole and directed deficit are already outputs of Theorem
   21.8's code construction.
4. **Relation to AO.2 and BCX.**  Those results use a repeated rank-one
   *cross shore* to expose fixed-scale child geometry.  BR.2 instead uses a
   weak rank-one bias in the compulsory boundary--boundary fill to scalarize
   an already separated conditional profile.  The proof mechanism is
   related, but the matched-roof scalarization and the suppression of both
   absolute channels do not appear in those theorems.  BCX already gives
   scalar physical separation for one special regular-Hadamard switching
   family; BR.2 is the general compiler for any matched-roof family meeting
   BR.10--BR.21.
5. **Relation to the exact coordinate compiler and orientation ceiling.**
   Theorem 21.29 locks an arbitrary coordinate exactly but pays a quadratic
   fill.  BR.2 keeps the fill at `O(n^(3/2))` by asking only for projective
   localization of a spectrally Lipschitz matched-roof family.  The
   orientation law OV.1 is respected: the query-owned fill itself has
   `Theta(n^(3/2))` internal cap, exactly enough to distinguish the two
   absolute channels at the target scale.  The biased-rounding device in
   BR.1 is the square symmetric analogue of the biased rectangular bridge
   in OV.2; the new content is what that device proves for matched roofs.
6. **The query language is not compressed.**  One freezes a query fill for
   each pole.  A generic exact error matrix `E_q` may require
   `Theta(n^2)` public description bits.  The conclusion is a lower bound on
   child-state information with query complexity uncharged, as in the
   contextual packing definition; it is not a small-description query-bank
   theorem.
7. **What remains.**  For exact minimizers, the single remaining structural
   input on this negative route is the uniform upper-tail deficit
   `L_tail`.  A uniform `O(sqrt n)` operator bound would imply it, while the
   separate tail-or-spectral-spike dichotomy describes the alternative.

The theorem proves information heaviness, not convergence or a cross-order
recurrence.  Its role is to decide whether the contextual-compression route
can possibly discard a linear rate of near-ground information.
