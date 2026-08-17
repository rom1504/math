# Rank-one common-factor amplification of the one-port Gram collision

**Status.** Rigorous exact lifting theorem, scalable one-port response
collision, exact-sign completion, and a sharp scope boundary for the obvious
regular-Hadamard tensor quotient.  The amplification uses a rank-one
ferromagnetic common factor.  It deliberately leaves the two-eigenvalue
Hadamard class; the separate four-port construction in
`regular_hadamard_equal_gram_rayleigh_collision.md` retains that stronger
class.

## 1. A general exact blow-up law

Let `H` be a symmetric entrywise sign matrix of order `n` with

```math
\operatorname{tr}H=0,                              \tag{RA.1}
```

and let `w in {+-1}^n`.  For an integer `m>=0`, define the one-port Boolean
trust response

```math
\mathcal B_m(H;w)
=\max_{x\in\{+-1\}^n}
 \left\{\left|\frac12x^THx\right|+m|w^Tx|\right\}.
                                                               \tag{RA.2}
```

Its spherical relaxation is

```math
\mathcal S_m(H;w)
=\max_{\|u\|_2^2=n}
 \left\{\left|\frac12u^THu\right|+m|w^Tu|\right\}.
                                                              \tag{RA.2a}
```

For `k>=1`, let

```math
H^{[k]}=H\otimes J_k,
\qquad w^{[k]}=w\otimes\mathbf1_k,
\qquad m^{[k]}=km,                                  \tag{RA.3}
```

and put

```math
\Gamma_k=\{-1,-1+2/k,\ldots,1-2/k,1\}.             \tag{RA.4}
```

The matrix `H^[k]` is still symmetric and entrywise sign, and its trace is
zero.  Hence deleting its diagonal gives a hollow complete signing with the
same Boolean quadratic energy.

### Theorem RA.1 (rank-one common-factor lifting law)

For every `k>=1`,

```math
\boxed{
\mathcal B_{km}(H\otimes J_k;w\otimes\mathbf1_k)
=k^2\max_{y\in\Gamma_k^n}
 \left\{\left|\frac12y^THy\right|+m|w^Ty|\right\}.}
                                                               \tag{RA.5}
```

Moreover, with

```math
D(H):=\frac12\sum_{i=1}^n|H_{ii}|,                 \tag{RA.6}
```

one has the uniform comparison

```math
\boxed{
k^2\mathcal B_m(H;w)
\le \mathcal B_{km}(H\otimes J_k;w\otimes\mathbf1_k)
\le k^2\bigl(\mathcal B_m(H;w)+D(H)\bigr).}       \tag{RA.7}
```

For an entrywise sign matrix, `D(H)=n/2`.

#### Proof

Write a Boolean vector on the product as fibers
`X=(X_1,...,X_n)`, with `X_i in {+-1}^k`, and set

```math
s_i=\mathbf1_k^TX_i,
\qquad y_i=s_i/k\in\Gamma_k.                       \tag{RA.8}
```

Every point of `Gamma_k^n` is realized by some collection of fibers.  Since
`J_k=1_k1_k^T`, exact contraction gives

```math
X^T(H\otimes J_k)X=s^THs=k^2y^THy,                 \tag{RA.9}
```

and

```math
(w\otimes\mathbf1_k)^TX=w^Ts=kw^Ty.                \tag{RA.10}
```

Substitution proves (RA.5).  Since the grid contains every Boolean vertex,
the lower bound in (RA.7) is immediate.

For the upper bound, fix `y in [-1,1]^n` and signs
`sigma,epsilon in {+-1}`.  Independently round `y_i` to a sign `Z_i` with
`E Z_i=y_i`.  Put

```math
f_{\sigma,\epsilon}(z)
=\frac\sigma2z^THz+m\epsilon w^Tz.                 \tag{RA.11}
```

Independence and `tr H=0` give

```math
f_{\sigma,\epsilon}(y)
=\mathbb E f_{\sigma,\epsilon}(Z)
 +\frac\sigma2\sum_iH_{ii}y_i^2.                  \tag{RA.12}
```

The expectation is at most the maximum of the same channel on the cube,
which is at most `B_m(H;w)`, while the last term is at most `D(H)`.
Maximizing over the four channels and over the box proves the upper bound in
(RA.7), hence also the grid bound. `square`

The diagonal allowance in (RA.7) is not a numerical artifact.  Blow-up
magnetizations are fractional, so the trace-zero cancellation that is exact
at Boolean vertices need not hold coordinate by coordinate.  The theorem
isolates its complete worst-case cost.

### Corollary RA.2 (seed-gap amplification criterion)

If two ports satisfy

```math
\mathcal B_m(H;w^+)-\mathcal B_m(H;w^-)>D(H),       \tag{RA.13}
```

then their rank-one blow-ups obey

```math
\mathcal B_{km}(H^{[k]};w^{+,[k]})
-\mathcal B_{km}(H^{[k]};w^{-,[k]})
\ge k^2\left[
\mathcal B_m(H;w^+)-\mathcal B_m(H;w^-)-D(H)
\right].                                           \tag{RA.14}
```

This is an exact scalable separation, not an extrapolation from finite
values.

## 2. Applying the order-16 zero-Rayleigh collision

Let `H_0` be the order-16 regular Walsh signing from
`regular_hadamard_equal_gram_rayleigh_collision.md`.  Thus

```math
H_0^2=16I,
\quad H_0\mathbf1=4\mathbf1,
\quad \operatorname{tr}H_0=0.                     \tag{RA.15}
```

Use the frozen one-port words

```text
w- = ----+--++--+----
w+ = +--+---++-------
```

Both have length 16 and satisfy

```math
(w^-)^TH_0w^-=(w^+)^TH_0w^+=0,                    \tag{RA.16}
```

but exact hypercube enumeration gives

```math
\mathcal B_4(H_0;w^-)=64,
\qquad
\mathcal B_4(H_0;w^+)=78.                         \tag{RA.17}
```

For every `k>=1`, put

```math
H_k=H_0\otimes J_k,
\quad n_k=16k,
\quad r_k=\|H_k\|_{op}=4k,
\quad m_k=4k,
\quad w_k^\pm=w^\pm\otimes\mathbf1_k.             \tag{RA.18}
```

Then `H_k` is regular with row sum `r_k`.  The total one-port mass remains

```math
\frac{m_k}{r_k}=1.                                 \tag{RA.19}
```

### Theorem RA.3 (scalable one-port equal-state separation)

The two lifted ports have exactly the same singleton Gram--Rayleigh state,

```math
G_k^-=G_k^+=[1],
\qquad R_k^-=R_k^+=[0],                            \tag{RA.20}
```

but for every `k>=1`,

```math
\boxed{
\mathcal B_{m_k}(H_k;w_k^+)
-\mathcal B_{m_k}(H_k;w_k^-)
\ge6k^2={3\over32}r_kn_k.}                        \tag{RA.21}
```

#### Proof

The Gram identity is immediate.  The Rayleigh numerator factors as

```math
(w^\pm)^TH_0w^\pm
\;\mathbf1_k^TJ_k\mathbf1_k=0,                    \tag{RA.22}
```

which proves (RA.20).  Here `D(H_0)=8`.  Apply the lower side of (RA.7) to
`w+` and the upper side to `w-`:

```math
\mathcal B_{4k}(H_k;w_k^+)\ge78k^2,
\qquad
\mathcal B_{4k}(H_k;w_k^-)\le(64+8)k^2.           \tag{RA.23}
```

Their difference is at least `6k^2`.  Since `r_kn_k=64k^2`, (RA.21)
follows. `square`

The spherical trust responses are nevertheless equal.  The range of `J_k`
reduces the nonzero part of `H_k` to `kH_0`, both ports have normalized
Rayleigh coordinate zero, and the kernel of `J_k` affects neither the field
nor the quadratic term.  Scaling shows that an optimizer uses the full norm
inside the range.  The two common spherical values are therefore

```math
\mathcal S_{m_k}(H_k;w_k^-)
=\mathcal S_{m_k}(H_k;w_k^+)
={3\sqrt3\over4}r_kn_k.                            \tag{RA.24}
```

Thus RA.3 is an equal-state, equal-spherical-value, unequal-Boolean-response
collision.

## 3. Completing the amplified response to exact signings

Delete the diagonal of `H_k` to obtain the old hollow signing.  Append
`m_k=4k` spins connected to the old block by the repeated port `w_k^\pm`.
Before filling the auxiliary shore, endpoint optimization is exactly
(RA.2).

There exists a hollow signing `C_M` on every sufficiently large order `M`
with

```math
Q(C_M)\le2M^{3/2}.                                 \tag{RA.25}
```

Indeed, choose its edge signs independently.  For a fixed spin the energy
is a sum of `L=M(M-1)/2` independent signs, so Hoeffding and a union bound
over `2^M` spins give

```math
\Pr\{Q(C_M)>2M^{3/2}\}
\le2^{M+1}e^{-4M}<1.                               \tag{RA.26}
```

Use the same `C_{4k}` for both parents.  Uniform cap Lipschitzness and
(RA.21) give

```math
B^+_{C}-B^-_{C}
\ge6k^2-2Q(C_{4k})
\ge6k^2-32k^{3/2}=\Theta(k^2).                    \tag{RA.27}
```

Hence the collision survives an exact complete-signing realization for all
large `k`.

This completed family has `r_k=Theta(n_k)`, not `sqrt(n_k)`.  It is a
theorem about scalable response information, not a near-original
`n^(3/2)` construction.  The four-port regular-Hadamard collision retains
`r=sqrt(n)` and supplies that stronger scale.

## 4. Why the obvious regular-Hadamard magnetization lift does not close

The proof of RA.1 uses the rank-one identity `J_k=1_k1_k^T`; it annihilates
every fiber mode invisible to the magnetization `y`.  This closure fails at
full leading order for an involutive regular-Hadamard factor.

Let `K` be a regular Hadamard signing of order `N`, with

```math
K^2=s^2I,
\qquad Ka=sa,
\qquad Kb=sb,\qquad a^Tb=0,                       \tag{RA.28}
```

for orthogonal Boolean top poles `a,b`.  The regular Walsh tensor family has
such a pair.  Form `H_0 tensor K` and use the port `w tensor a`.  The natural
top-pole magnetization of a fiber `X_i` is `a^TX_i/N`.

### Proposition RA.4 (leading hidden-mode obstruction)

The two Boolean product spins

```math
X_1=\mathbf1_{16}\otimes b,
\qquad X_2=w^-\otimes b                           \tag{RA.29}
```

have the same zero top-pole magnetization in every fiber and zero port field,
but their quadratic energies are respectively

```math
\left|\frac12X_1^T(H_0\otimes K)X_1\right|
={r'n'\over2},
\qquad
\frac12X_2^T(H_0\otimes K)X_2=0,                 \tag{RA.30}
```

where `n'=16N` and `r'=4s`.

#### Proof

Both magnetizations and fields vanish because `a^Tb=0`.  Product quadratic
forms factor.  Since `H_0 1=4 1` and `Kb=sb`,

```math
\frac12X_1^T(H_0\otimes K)X_1
=\frac12(64)(sN)=32sN={r'n'\over2}.               \tag{RA.31}
```

Equation (RA.16) makes the second product zero. `square`

Thus no quotient retaining only the top-pole fiber magnetizations can extend
the exact rank-one lifting identity to the regular-Hadamard tensor.  It loses
a full leading `rn/2` hidden-mode variable.  Proposition RA.4 is a no-go for
the **obvious quotient proof**, not a proof that the actual one-port tensor
responses become equal or remain separated.  The p=1 strict involutive
response amplification remains unresolved by this argument; the p=4 support
construction avoids the issue and already gives the required strict-class
collision.

## 5. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_rank_one_common_factor_amplification.py
```

The verifier exhausts the order-16 seed, checks the exact blow-up contraction
on random product spins, verifies the equal normalized states and product
witnesses for several `k`, and checks the leading hidden-mode obstruction in
the first regular-Hadamard tensor extension.
