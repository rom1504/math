# Equal Gram--Rayleigh states with separated Boolean trust responses

**Status.**  Rigorous scalable counterexample.  Two four-port systems in the
same dense regular-Hadamard family have exactly the same full Gram--Rayleigh
state `(G,R)` but Boolean trust responses separated by at least `rn/8` at
bounded total port mass.  The separation survives tensor powers and arbitrary
public exact-sign completion at leading order.

This closes the information-sufficiency question left open by the spherical
Gram carrier: `(G,R)` determines the spherical trust response, but it does not
determine the exact Boolean response, even up to `o(rn)`.

## 1. Response and state

Let `H` be a symmetric trace-zero sign matrix with

```math
H^2=r^2I.
```

For an ordered Boolean port tuple `W=(w_1,...,w_p)` and a common integer
multiplicity `m`, define the old-spin Boolean trust response

```math
\mathcal B_m(H;W)
=\max_{x\in\{\pm1\}^n}
 \left\{\left|\frac12x^THx\right|
 +m\sum_{i=1}^p|w_i^Tx|\right\}.                 \tag{EG.1}
```

Its normalized Gram--Rayleigh state is

```math
G_{ij}=\frac{w_i^Tw_j}{n},
\qquad
R_{ij}=\frac{w_i^THw_j}{rn}.                     \tag{EG.2}
```

The collective pure-field support is

```math
L(W):=\max_x\sum_i|w_i^Tx|
=\max_{\epsilon\in\{\pm1\}^p}
 \left\|\sum_i\epsilon_iw_i\right\|_1.          \tag{EG.3}
```

Unlike the spherical support, `L(W)` is not determined by `G`.

## 2. The order-16 seed

Use the regular symmetric Walsh matrix `H_0` returned by
`regular_hadamard(2)` in the project verifier.  Thus

```math
n_0=16,\qquad r_0=4,\qquad H_0^2=16I,
\qquad H_0\mathbf1=4\mathbf1,
\qquad \operatorname{tr}H_0=0.                   \tag{EG.4}
```

Write `+` and `-` for `+1` and `-1`.  Define five Boolean top eigenvectors:

```text
w0 = ----------------
w1 = -----++--++-----
w2 = ---+--+-+-++-+++
w3 = ---+++-+-+---+++
w4 = --+----+-++++-++
```

Every displayed word obeys `H_0 w_i=4w_i`.  Consider

```math
W^A=(w_0,w_1,w_2,w_3),
\qquad
W^B=(w_0,w_1,w_2,w_4).                            \tag{EG.5}
```

### Lemma EG.1 (exact equal state, unequal support)

The two tuples have the same exact normalized state

```math
G^A=G^B=R^A=R^B=
\begin{pmatrix}
1&1/2&0&0\\
1/2&1&0&0\\
0&0&1&0\\
0&0&0&1
\end{pmatrix},                                   \tag{EG.6}
```

but

```math
\boxed{L(W^A)=32,\qquad L(W^B)=28.}               \tag{EG.7}
```

Moreover the Boolean word

```text
x_A = ------+---+-----
```

satisfies

```math
\frac12x_A^TH_0x_A=24,
\qquad
\sum_{i=1}^4|w_i^{A\,T}x_A|=32.                  \tag{EG.8}
```

#### Proof

Direct inner products give the unnormalized Gram matrix

```math
\begin{pmatrix}
16&8&0&0\\
8&16&0&0\\
0&0&16&0\\
0&0&0&16
\end{pmatrix}
```

for both tuples.  Since every port is a top eigenvector, its Rayleigh matrix
is four times this matrix, proving (EG.6).

There are only sixteen endpoint sign words in (EG.3).  Exact evaluation gives
the possible `l_1` supports `{16,24,32}` for `W^A` and `{20,28}` for `W^B`,
which proves (EG.7).  Substitution of the displayed `x_A` gives (EG.8).
All integer checks are included in the verifier. `square`

For orientation, exact enumeration of all `2^16` old spins gives

```math
\mathcal B_2(H_0;W^A)=88,
\qquad
\mathcal B_2(H_0;W^B)=82,                         \tag{EG.9}
```

and

```math
\mathcal B_4(H_0;W^A)=152,
\qquad
\mathcal B_4(H_0;W^B)=138.                       \tag{EG.10}
```

These finite values motivated the scalable argument, but the theorem below
does not rely on extrapolating them.

## 3. Dense tensor amplification

For `j>=1`, put

```math
H_j=H_0^{\otimes j},
\quad n_j=16^j,
\quad r_j=4^j,
\quad N_j=16^{j-1}.                               \tag{EG.11}
```

For every seed port set

```math
w_{i,j}=w_i\otimes\mathbf1_{N_j},
\qquad
W_j^A=(w_{0,j},w_{1,j},w_{2,j},w_{3,j}),
\qquad
W_j^B=(w_{0,j},w_{1,j},w_{2,j},w_{4,j}).          \tag{EG.12}
```

These are Boolean ports in one dense entrywise-sign common factor.  Choose

```math
m_j=r_j.                                           \tag{EG.13}
```

The collective port-mass parameter remains bounded:

```math
c={m_jp\over r_j}=4.                               \tag{EG.14}
```

### Theorem EG.2 (scalable equal-state Boolean separation)

For every `j>=1`,

```math
(G_j^A,R_j^A)=(G_j^B,R_j^B),                       \tag{EG.15}
```

exactly, but

```math
\boxed{
\mathcal B_{r_j}(H_j;W_j^A)
-\mathcal B_{r_j}(H_j;W_j^B)
\ge {r_jn_j\over8}.}                              \tag{EG.16}
```

In this family `r_j=sqrt(n_j)`, so the gap is `n_j^(3/2)/8`.

#### Proof

The all-one vector is a top pole of the last `j-1` tensor factors.  Hence
all inner products in (EG.2) multiply by `N_j`, while the Rayleigh numerators
multiply by `r_j/4` in addition.  Normalization cancels these factors, so
the common state remains (EG.6).  This proves (EG.15).

The pure-field support tensorizes exactly:

```math
L(W_j^A)=32N_j,
\qquad
L(W_j^B)=28N_j.                                   \tag{EG.17}
```

Indeed every signed port sum is the corresponding seed sum tensored with
`mathbf1`, so its `l_1` norm is multiplied by `N_j`.

For the `A` response, evaluate at

```math
x_{A,j}=x_A\otimes\mathbf1_{N_j}.
```

Its quadratic term is

```math
\left|\frac12x_{A,j}^TH_jx_{A,j}\right|
=24\,{r_j\over4}N_j=6r_jN_j,                     \tag{EG.18}
```

and its field term at `m_j=r_j` is `32r_jN_j`.  Thus

```math
\mathcal B_{r_j}(H_j;W_j^A)\ge38r_jN_j.           \tag{EG.19}
```

For every Boolean `x`, the spectral bound and (EG.17) give

```math
\mathcal B_{r_j}(H_j;W_j^B)
\le {r_jn_j\over2}+28r_jN_j
=36r_jN_j,                                        \tag{EG.20}
```

because `n_j=16N_j`.  Subtraction gives

```math
2r_jN_j={r_jn_j\over8},
```

proving (EG.16). `square`

The proof uses neither a separately paid scalar decomposition nor the
spherical relaxation.  It exhibits one Boolean witness on the high side and
uses a state-invisible fourth-order support statistic on the low side.

## 4. Exact-sign completion

The response (EG.1) is **exactly** the cap of the incomplete old--new parent,
not a further relaxation.  Delete the diagonal of `H_j` and call the resulting
hollow signing `A_j`.  Since `tr(H_j)=0`, for every Boolean old spin

```math
\sum_{a<b}(A_j)_{ab}x_ax_b
=\frac12\left(x^TH_jx-\operatorname{tr}H_j\right)
=\frac12x^TH_jx.                                  \tag{EG.21}
```

Append four shores, each of width `m_j=r_j`, and give every vertex in shore
`i` the old--new column `w_{i,j}`.  For fixed `x`, independently choosing the
four shore signs gives

```math
\max_y\left|
 \frac12x^TH_jx+
 \sum_{i=1}^4\sum_{a=1}^{m_j}y_{ia}w_{i,j}^Tx
 \right|
=\left|\frac12x^TH_jx\right|
 +m_j\sum_{i=1}^4|w_{i,j}^Tx|.                   \tag{EG.22}
```

Indeed choose every endpoint sign so that its bridge term has the sign of
the old quadratic term (either sign may be used when that term vanishes).
Maximizing over `x` proves the exact identification with (EG.1).

Now fill all pairs among the `4r_j` new vertices by the same arbitrary hollow
sign matrix `C_j`.  Pointwise,

```math
\bigl||E(x,y)+H_{C_j}(y)|-|E(x,y)|\bigr|
\le |H_{C_j}(y)|\le Q(C_j).                       \tag{EG.23}
```

Taking maxima proves that either completed cap differs from its incomplete
cap by at most `Q(C_j)`.  Thus, if the two completed responses are
`B^A_C,B^B_C`, then

```math
B^A_C-B^B_C
\ge {r_jn_j\over8}-2Q(C_j).                       \tag{EG.24}
```

Since

```math
Q(C_j)\le {4r_j\choose2}=O(r_j^2)=O(n_j),         \tag{EG.25}
```

while `r_jn_j=n_j^(3/2)`, every public completion retains the leading gap:

```math
B^A_C-B^B_C
\ge\left({1\over8}-o(1)\right)n_j^{3/2}.          \tag{EG.26}
```

Every off-diagonal entry of the completed parent is an exact sign: the old
block uses the off-diagonal entries of `H_j`, the bridge uses Boolean ports,
and the auxiliary block uses `C_j`.  Its total order is
`n_j+4r_j=n_j+4sqrt(n_j)`, so (EG.26) has the same leading normalization in
total-parent units.

## 5. Consequence for the Gram carrier

The spherical trust formula is a function of `(G,R)`, so the two systems have
identical spherical responses in every endpoint channel.  More strongly,
Theorem EG.2 is an information collision for **any** Boolean decoder from this
state.  Let a proposed decoder receive `(G,R,n,r,m,p)` and even the common
public completion `C_j`.  It must return the same number on the two inputs.
Since their completed caps differ by `(1/8-o(1))rn`, at least one decoding
error is `(1/16-o(1))rn`.  Therefore no `(G,R)`-only decoder can have uniform
`o(rn)` Boolean error on this class.  This conclusion is not limited to
rounding the particular spherical value.

The missing datum is already visible as the signed-sum support body

```math
\epsilon\longmapsto
\left\|\sum_i\epsilon_iw_i\right\|_1,             \tag{EG.27}
```

or, equivalently here, a fourth-order row-pattern statistic.  This statistic
is strictly smaller than the full Boolean landscape, but it is not encoded by
second-order Gram--Rayleigh geometry.

The conclusion is information-theoretic, not merely an integrality gap:
equal carrier states have different exact outputs.  It does not say that the
full support body is always necessary, nor does it rule out a different
compressed Boolean carrier carrying suitable higher-order information.

## 6. Secondary one-port collision and flatness

There is also a minimal finite collision at one port.  At order 16 let

```text
u- = ----+--++--+----
u+ = +--+---++-------
```

(the strings have length 16 as spelled explicitly in the verifier).  They
satisfy

```math
G^- =G^+ =[1],\qquad R^-=R^+=[0],                 \tag{EG.28}
```

but at `m=r=4`,

```math
\mathcal B_4(H_0;u^-)=64,
\qquad
\mathcal B_4(H_0;u^+)=78.                         \tag{EG.29}
```

The projective Hamming-shell maxima of the quadratic child energy, for
distances `0,...,8`, are respectively

```text
u-: 0, 6, 16, 22, 32, 26, 24, 26, 32
u+: 0,22, 24, 26, 32, 26, 32, 26, 32.
```

Since

```math
\mathcal B_m(H;u)=mn+max_{0\le d\le n/2}
 \{E_u(d)-2md\},                                  \tag{EG.30}
```

these tables are a compact exhaustive certificate for (EG.29).

This one-port collision also has a dense scalable exposed-flatness version.
For any zero-Rayleigh port put `J=H/r` and `z=Jw`.  At `m=r`, the unique
optimizer in quadratic channel `sigma` is

```math
v_\sigma={\sqrt3w+\sigma z\over2},                \tag{EG.31}
```

with common spherical value `3sqrt(3)rn/4`.  For `u^-`, both channels have
`l_1` norm `8sqrt(3)`, while for `u^+` the positive channel has `l_1` norm
`3+7sqrt(3)`.  Thus their minimum exposed flatnesses differ by

```math
{3-\sqrt3\over16}.                                 \tag{EG.32}
```

Tensoring each port with a common all-one top pole repeats (EG.31)
coordinatewise, so this flatness separation persists through the dense
regular-Hadamard tensor family.  This secondary fact is a recovery-state
collision; Theorem EG.2 is the stronger exact-response collision.

## 7. Audit of the pure-linear four-port seed

The separate draft `four_port_gram_boolean_collision.md` is correct.  Its
uniform four-bit row distribution and duplicated even-parity distribution
both have Gram matrix `I_4`, while their signed-sum `l_1` supports are
`3n/2` and `2n`.  The claimed minimality of four ports for that pure-linear
moment mechanism is also correct: after quotienting row patterns by global
sign, three ports leave four frequencies, and total mass plus the three pair
moments form an invertible four-character Walsh system.  At four ports there
are eight projective frequencies but only seven mass/pair constraints; the
fourth-order parity character spans the missing direction.

The Walsh tuples in EG.5 are not literally the uniform/even-parity seed, but
they implement the same lesson inside equal full `(G,R)` states: pair moments
agree while a higher row-pattern statistic changes the Boolean support.

## 8. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_equal_gram_rayleigh_collision.py
```

The verifier checks every displayed inner product, exhausts all endpoint and
old-spin words at order 16, verifies the tensor formulas, performs the full
one-port Rayleigh-class search, and checks the exposed-flatness identities.
