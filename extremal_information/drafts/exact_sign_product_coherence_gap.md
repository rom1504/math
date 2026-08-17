# Exact-sign product coherence can fail at five ports

**Status.**  Rigorous finite seed, scale-preserving exact-sign tensor
amplification, and finite verifier.  The result separates *joint selector
defect* from every individual product Rayleigh deficit.  It does **not**
separate the full Boolean trust optima, and its marginal deficit is the fixed
constant `3/16`, not a quantity tending to zero.

## 1. The defect tested here

Let `H` be symmetric with `H^2=r^2I`, let
`w_1,...,w_p in {+-1}^n`, and put

```math
z_S=\bigodot_{i\in S}w_i
```

for every odd set in the Fourier support of odd majority.  The individual
positive-Rayleigh deficit and the majority-selector deficit at endpoint
`epsilon` are

```math
\begin{aligned}
d_S&=1-{z_S^THz_S\over rn},\\
d_{\rm maj}(\epsilon)
 &=1-{x_\epsilon^THx_\epsilon\over rn},\\
x_\epsilon(j)&=\operatorname {Maj}
 (\epsilon_1w_1(j),\ldots,\epsilon_pw_p(j)).       \tag{ES.1}
\end{aligned}
```

Writing `a^epsilon` for the majority Fourier coefficient vector and
`D=Z^T(I-H/r)Z/n`, the second quantity is exactly

```math
d_{\rm maj}(\epsilon)=(a^\epsilon)^TD a^\epsilon. \tag{ES.2}
```

Thus (ES.1) asks whether separately good product poles force the *one same
Boolean spin* which pays all port fields jointly to have a good child
energy.  The weighted-contraction example CB.1 says no with vanishing
marginal deficits.  The theorem below is the first exact-sign version, at a
fixed nonzero marginal deficit.

## 2. Why the existing three-port seed is rigid

Let `H_0` be the order-16 regular Hadamard matrix from TC.2 and let its
product-closed triple be

```text
a = +-----+--+-----+
b = +--+--------+--+
c = +++-+-++++-+-+++
t = +++++--++--+++++
```

where `t=a odot b odot c`.  Put

```math
x_\epsilon={\epsilon_1a+\epsilon_2b+\epsilon_3c
 -\epsilon_1\epsilon_2\epsilon_3t\over2}.          \tag{ES.3}
```

### Lemma ES.1 (projective selector--pole rigidity at three ports)

For this particular seed,

```math
\boxed{
\{[x_\epsilon]:\epsilon\in\{+-1\}^3\}
=\{[a],[b],[c],[t]\}.}                             \tag{ES.4}
```

More precisely, in lexicographic endpoint order the eight witnesses are

```text
-a, -b, -c, +t, -t, +c, +b, +a.
```

Consequently, for **every** real symmetric matrix `K`, not only for a
Hadamard child or a diagonal switching of one,

```math
\max_\epsilon x_\epsilon^TKx_\epsilon
=\max_{z\in\{a,b,c,t\}}z^TKz,                    \tag{ES.5}
```

and the same holds with minima, absolute values, or any scalar function of
the quadratic value.  This explains the exhaustive diagonal-switching
rigidity seen at `p=3`: the selector and product test sets are literally the
same projective set.

#### Proof

Substitution in (ES.3) gives the displayed eight-word list.  Quadratic
forms are invariant under a global sign, proving (ES.4)--(ES.5). `square`

This is seed-specific rather than a general three-port theorem.  It also
does not persist in the growing PC.3 presentation.  At tensor depth two,
the five PC.3 ports have 16 projective odd products and 16 projective
majority witnesses, but the two sets have only one common element.  Exact
product closure makes both sets top for the unperturbed child; it does not
make them setwise identical.

## 3. An exact-sign five-port gap

Use the same regular Hadamard matrix `H_0`, so

```math
H_0^2=16I,
\qquad H_0\mathbf1=4\mathbf1,
\qquad \operatorname {tr}H_0=0.                  \tag{ES.6}
```

Take the following five Boolean ports:

```text
w1 = -++-+--++-++-++-
w2 = +++-+-+++--+-+++
w3 = ----------+-----
w4 = +-----+--++----+
w5 = +--+-++--+--+--+
```

For five-input majority every odd subset is active.  Its normalized Fourier
coefficients are `3/8` in degrees one and five and `-1/8` in degree three.

### Theorem ES.2 (exact-sign marginal/coherent separation)

For all 16 odd subsets `S subseteq [5]`,

```math
\boxed{z_S^TH_0z_S=52,
\qquad d_S={3\over16}.}                            \tag{ES.7}
```

At the all-positive endpoint, however,

```text
x = +-----+---+----+
```

and

```math
\boxed{x^TH_0x=32,
\qquad d_{\rm maj}(\mathbf1)={1\over2}.}          \tag{ES.8}
```

Hence

```math
d_{\rm maj}(\mathbf1)-\max_Sd_S={5\over16}.       \tag{ES.9}
```

In the robust-selector normalization, every separately inspected product
witness loses `3rn/32` from the positive child bound, whereas the one
majority witness loses `rn/4`.  The additional coherent loss is `5rn/32`.

The value `3/16` is also the smallest nonzero positive-Rayleigh deficit of
a Boolean vector for this `H_0`: exhaustive exact evaluation gives top
quadratic value `64` and next value `52`.  Thus the construction uses the
closest non-top Boolean shell available in this seed.

#### Proof

Multiply each of the 16 displayed odd coordinate products by `H_0`; exact
integer inner products give (ES.7).  The pointwise majority of the five
displayed rows is the word in (ES.8), and another integer multiplication
gives its quadratic value `32`.  Dividing by `rn=64` proves the deficits.

Equivalently, the exact majority identity

```math
x={1\over8}\left(
3\sum_{|S|=1}z_S-\sum_{|S|=3}z_S+3z_{[5]}
\right)                                             \tag{ES.10}
```

shows directly that the joint value is a coherent quadratic combination of
the product poles, rather than an independently paid sum of their diagonal
values.  The verifier exhausts all `2^16` Boolean vectors for the final
shell claim. `square`

## 4. Scale-preserving tensor amplification

For `j>=1`, put

```math
H_j=H_0^{\otimes j},\qquad
n_j=16^j,\qquad r_j=4^j,\qquad
w_{i,j}=w_i\otimes\mathbf1_{16^{j-1}}.             \tag{ES.11}
```

### Corollary ES.3 (scalable exact-sign coherence gap)

Every `H_j` is symmetric and entrywise signed, satisfies

```math
H_j^2=n_jI,
\qquad \operatorname {tr}H_j=0,
\qquad r_j=\sqrt {n_j},                            \tag{ES.12}
```

and all ratios in (ES.7)--(ES.9) remain exact at order `n_j`.

Indeed `z_{S,j}=z_S tensor 1` and `x_j=x tensor 1`, while the tail all-one
vector is a positive top pole.  Therefore both Rayleigh numerators acquire
the same factor `r_{j-1}n_{j-1}` as their normalizing denominator.  The
separation is consequently a fixed `5/16` at the natural
`r_jn_j=n_j^(3/2)` scale.

Deleting the diagonal of `H_j` gives a hollow exact signing.  Since its
trace is zero, every Boolean quadratic energy is unchanged:

```math
y^T(H_j-\operatorname {diag}H_j)y=y^TH_jy.         \tag{ES.13}
```

Thus this is not merely a weighted contraction construction.  The exact
deficit ratios above use the diagonal-completed Hadamard roof `r_j=sqrt n_j`.
For the hollow child `A_j=H_j-diag(H_j)`, the safe contraction scale
`r'_j=sqrt(n_j)+1` gives deficits

```math
1-{13\over16}{\sqrt{n_j}\over\sqrt{n_j}+1},
\qquad
1-{1\over2}{\sqrt{n_j}\over\sqrt{n_j}+1},        \tag{ES.14}
```

whose excess tends to `5/16`.  Hollow Boolean energies are exact, but the
hollow and completed operator roofs are not identical.

## 5. Exact scope and next boundary

The theorem proves that exact signs do **not** impose the three-port
setwise rigidity at five ports.  It rules out any inference of joint
selector quality from the statement `d_S<=3/16` for every active product,
even in a regular-Hadamard family and at every tensor scale.

It does not yet match the asymptotic strength of CB.1:

1. `p=5` is fixed;
2. the marginal deficit remains `3/16` rather than tending to zero;
3. (ES.8) is the loss of the prescribed same-field selector witness, not a
   proof that the full Boolean trust optimum has the same loss.  A different
   Boolean spin may repair some of it.

The sharp remaining exact-sign question is therefore whether there are
growing generated port algebras with

```math
\max_S d_S=o(1)
\quad\hbox{but}\quad
\max_\epsilon d_{\rm maj}(\epsilon)\ge c>0,       \tag{ES.15}
```

or whether regular-Hadamard/sign integrality forces a genuine coherence
theorem once all marginal deficits vanish.

## 6. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_exact_sign_product_coherence_gap.py
```

The verifier reconstructs `H_0`, checks every Fourier coefficient and odd
product, exhausts the order-16 Rayleigh shell, verifies the tensor identities
physically at depth two, and checks both the three-port rigidity and its
failure for the five-port PC.3 tensor presentation.
