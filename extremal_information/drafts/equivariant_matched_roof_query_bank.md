# One switched boundary fill gives the complete matched-roof query bank

**Status.** Task-local theorem draft awaiting independent audit.  This note
removes the quadratic *per-query* description caveat from BR.2 and records
the exact response-information rate of every bounded-cap switching orbit.
It does not charge the one-time public description of the common bridge,
base fill, or codebook.

## 1. Equivariant exact-sign fills

For `q in {+-1}^n`, put `S_q=diag(q)` and

```math
R_q=qq^T-I.
```

### Lemma EQ.1 (one base fill serves every pole)

Fix `lambda>0` and put `p=lambda/sqrt n`.  For all sufficiently large `n`
there is one hollow exact signing `D_1` with

```math
D_1=pR_1+E_1,
\qquad \|E_1\|_{2\to2}\le K_0\sqrt n.           \tag{EQ.1}
```

For every pole define

```math
\boxed{D_q=S_qD_1S_q.}                           \tag{EQ.2}
```

Then, simultaneously for every `q`,

```math
D_q=pR_q+E_q,
\qquad E_q=S_qE_1S_q,
\qquad \|E_q\|_{2\to2}\le K_0\sqrt n,           \tag{EQ.3}
```

and `D_q` is a hollow exact signing with

```math
Q(D_q)\le{\lambda+K_0+o(1)\over2}n^{3/2}.       \tag{EQ.4}
```

#### Proof

Lemma BR.1, applied only at the all-positive pole, supplies `D_1`.  Since

```math
S_qR_1S_q=S_q(J-I)S_q=qq^T-I=R_q,
```

equations (EQ.2)--(EQ.3) follow.  Diagonal switching preserves hollowness,
exact signs, operator norm, and Boolean cap, proving (EQ.4). `square`

Thus the random-existence argument is used once, not once per pole.  The
entire fill bank is an exact switching orbit of one public base signing.

### Corollary EQ.2 (simultaneous equivariant BR selector)

Fix the constants `L,delta,C_P` in BR.2.  Choose `theta,lambda` there once,
and choose `D_1` from EQ.1.  Then every directed matched-roof certificate at
every pole `q` is scalarized by the fill (EQ.2), with the same response gap
and parent-cap constants as BR.2.

#### Proof

The proof of BR.2 uses no entry of `D_q` beyond the decomposition
`D_q=pR_q+E_q` and the uniform operator bound on `E_q`.  Equations
(BR.8)--(BR.9), the projective near/far split, and suppression of the
negative absolute channel therefore hold verbatim for every switched fill
in (EQ.2).  No simultaneous union bound is required. `square`

Once `D_1` is public, a query fill is specified by its projective pole using
at most `n-1` varying bits (or by an index of length `ceil(log_2|I|)` in a
public selected codebook).  The full fill orbit has exactly `2^(n-1)`
members: if `S_rD_1S_r=D_1`, completeness and nonzero entries force
`r_ir_j=1` for every `i ne j`, hence `r=+-1`.  Storing an arbitrary public
base uses at most `binom(n,2)` bits; no subquadratic explicit construction is
proved here.  The response model does not charge the one-time public base or
codebook.

## 2. Exact switching-orbit response rate

Fix a hollow complete signing `A` and its diagonal switching orbit

```math
\mathcal O(A)=\{S_sAS_s:s\in\{+-1\}^n/\{+-1\}\}.
```

For a fixed public cross block `B` and a public query fill `D`, let the
scalar response of a child `C` be the cap of the complete parent

```math
\operatorname{Resp}_{B,D}(C)
=Q\begin{pmatrix}C&B\\B^T&D\end{pmatrix}.        \tag{EQ.5}
```

### Theorem EQ.3 (every bounded-cap switching orbit has `Theta(n)` response bits)

Fix `C_0>0`.  There are constants `gamma,c,C_1>0` such that every
sufficiently large order-`n` complete signing satisfying

```math
Q(A)\le C_0n^{3/2}                               \tag{EQ.6}
```

has a suborbit `\{A_i:i\in I\}` of size

```math
|I|\ge\exp(\gamma n),                            \tag{EQ.7}
```

one common exact-sign bridge `B` with `\|B\|_{2\to2}=O(\sqrt n)`, and an
equivariant query bank `D_i=S_{q_i}D_1S_{q_i}` such that

```math
\operatorname{Resp}_{B,D_i}(A_i)
-\operatorname{Resp}_{B,D_i}(A_j)
\ge c n^{3/2}
\quad(i\ne j),                                  \tag{EQ.8}
```

while every displayed parent has order `2n`, exact signs, and cap at most
`C_1n^(3/2)`.

Consequently any encoder/decoder which approximates all responses in this
public bank with uniform error less than `cn^(3/2)/3` needs at least
`exp(gamma n)` states, hence `Omega(n)` bits.  Conversely `n-1` bits always
suffice to identify the switch and reconstruct the child exactly.  Thus the
worst-case fixed-accuracy response-description complexity of this switching
orbit is

```math
\boxed{\Theta(n)\text{ bits}.}                   \tag{EQ.9}
```

#### Proof

Choose `sigma in {+-1}` so that `P(sigma A)=Q(A)`.  Theorem BT.3 gives a
uniform positive near-top entropy deficit depending only on `C_0`.
Theorem 21.8 applied to `sigma A` produces one bounded-operator exact-sign
bridge `B`, an exponential switching subcode, and directed matched-roof
deficits at the named code poles.  Corollary EQ.2 scalarizes them using the
bank generated from `D_1`.  If `sigma=-1`, globally negate every resulting
parent.  Its child block becomes

```math
sigma(\sigma A)^{s_i}=A^{s_i},
```

while the common bridge and fill bank become `sigma B` and
`sigma D_i=S_(q_i)(sigma D_1)S_(q_i)`.  Absolute caps and all response gaps
are unchanged.  Thus in either orientation the children in
(EQ.7)--(EQ.8) belong to the switching orbit of the originally declared
`A`, and the parent-cap assertion follows.

If two children shared one approximate state, their decoded response to
query `i` would be the same, while the two true responses differ by at least
`cn^(3/2)`; two errors smaller than `cn^(3/2)/3` cannot cover the gap.
Hence the encoder is injective on the subcode.  For the upper bound, encode
the projective switch label `s`, using at most `n-1` bits.  Given the public
base child, that label reconstructs `S_sAS_s` and therefore every declared
future response exactly, with no computational-efficiency claim. `square`

This theorem applies to every bounded-cap signing, not only exact or near
minimizers.  In particular it proves the campaign's physical multi-selector
lower bound at a linear information rate.  It rules out bounded,
poly(`n`)-state, and more generally `exp(o(n))`-state response summaries for
this future language.  It does **not** rule out an extensive `O(n)`-bit
state, since the switch label itself is one; it does not assert that the
order-`2n` parents are near-minimal; and it supplies no cross-order update or
recurrence.

## 3. Archive comparison

- BR.1--BR.3 prove the biased matched-roof selector, but choose each fill
  separately and explicitly leave query description uncharged.  EQ.1 is
  the equivariant simultaneous refinement and identifies the varying query
  label as linear rather than quadratic.
- BCL.1 proves `Theta(n)` response bits for a regular-Walsh switching orbit.
  EQ.3 extends the rate statement to **every** bounded-cap complete signing,
  using BT.3 instead of Walsh spectral flatness.
- Switching itself and the `n-1`-bit exact upper bound are elementary, and
  the linear lower-rate implication is already latent in the audited
  BT.3--Theorem 21.8--BR chain.  The genuinely new increment here is that
  one common switched physical query architecture realizes the whole bank;
  EQ.3 records the resulting exact `Theta(n)` classification cleanly.

Classification:

```text
PROVES / CLOSES:
  the physical multi-selector information rate is Theta(n) bits;
  query fills can be generated equivariantly from one common base.

DOES NOT PROVE:
  near-minimality of query parents, subextensive compositional state,
  cross-order transfer, or convergence.
```
