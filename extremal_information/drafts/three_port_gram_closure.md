# Gram closure through three Boolean top ports, and its four-port boundary

**Status.**  Rigorous algebraic closure theorem, scalable exact-sign
regular-Hadamard counterexample at four ports, and finite verifier.  The
positive theorem identifies a nontrivial class on which the Gram--Rayleigh
state determines the exact Boolean trust response.  The counterexample
shows that the conclusion already fails at four ports, even when all ports
are Boolean top eigenvectors of one common regular Hadamard matrix, so that
`R=G` exactly.

## 1. Trust responses

Let `H` be symmetric with

```math
H^2=r^2I,
\qquad \operatorname {tr}H=0,                         \tag{TC.1}
```

and let `w_1,...,w_p in {+-1}^n` be Boolean `+r`
eigenvectors.  For an integer port width `m>=0`, put

```math
\begin{aligned}
\mathcal B_m(H;w_1,\ldots,w_p)
 &=\max_{x\in\{+-1\}^n}
   \left\{ {1\over2}|x^THx|+m\sum_{i=1}^p|w_i^Tx|\right\},\\
\mathcal S_m(H;w_1,\ldots,w_p)
 &=\max_{\|u\|_2^2=n}
   \left\{ {1\over2}|u^THu|+m\sum_{i=1}^p|w_i^Tu|\right\}.
                                                               \tag{TC.2}
\end{aligned}
```

These are respectively the Boolean and spherical versions of the complete
channel maximum over the outer quadratic sign and all endpoint signs.  Write

```math
G_{ij}={w_i^Tw_j\over n},
\qquad
R_{ij}={w_i^THw_j\over rn}.                          \tag{TC.3}
```

Here `R=G`, since all ports lie in the positive top eigenspace.

## 2. Two ports close automatically

### Proposition TC.1 (two-port exact Gram formula)

For two Boolean top ports with `rho=w_1^Tw_2/n`,

```math
\boxed{
\mathcal B_m={rn\over2}+mn(1+|\rho|),
\qquad
\mathcal S_m={rn\over2}+mn\sqrt{2(1+|\rho|)}.}       \tag{TC.4}
```

#### Proof

For `epsilon_i in {+-1}`, Boolean duality gives

```math
\max_x\sum_i|w_i^Tx|
=\max_\epsilon\|\epsilon_1w_1+\epsilon_2w_2\|_1
=n(1+|\rho|).                                      \tag{TC.5}
```

For the maximizing endpoint sign, one of the signed ports itself attains
the field support and the positive child bound `rn/2`.  The separate child
and field bounds prove Boolean equality.  The analogous Euclidean field
support is `n sqrt(2(1+|rho|))`; its normalized port sum remains in the
positive top eigenspace, so it simultaneously attains the spherical child
bound. `square`

## 3. A three-port algebraic closure theorem

For three ports define

```math
T(G)=\max_{\epsilon\in\{+-1\}^3}
       \sum_{1\le i<j\le3}G_{ij}\epsilon_i\epsilon_j. \tag{TC.6}
```

### Theorem TC.2 (triple-product closure)

Suppose, in addition to (TC.1), that

```math
t=w_1\odot w_2\odot w_3
\quad\hbox{satisfies}\quad Ht=rt.                  \tag{TC.7}
```

Then the exact Boolean and spherical responses are determined by the
three-port Gram matrix:

```math
\boxed{
\begin{aligned}
\mathcal B_m
 &= {rn\over2}+{mn\over2}\bigl(3+T(G)\bigr),\\
\mathcal S_m
 &= {rn\over2}+mn\sqrt{3+2T(G)}.
\end{aligned}}                                      \tag{TC.8}
```

#### Proof

For three scalar signs,

```math
\operatorname {sgn}(a+b+c)={a+b+c-abc\over2},
\qquad
|a+b+c|={3+ab+ac+bc\over2}.                        \tag{TC.9}
```

There are no ties.  Apply the first identity coordinatewise to
`a=epsilon_1w_1`, `b=epsilon_2w_2`, `c=epsilon_3w_3`.  The resulting
Boolean sign vector is

```math
x_\epsilon={\epsilon_1w_1+\epsilon_2w_2+
 \epsilon_3w_3-\epsilon_1\epsilon_2\epsilon_3t\over2}.
                                                               \tag{TC.10}
```

Hypothesis (TC.7) puts every term, hence `x_epsilon`, in the positive top
eigenspace.  Thus `x_epsilon` simultaneously attains the Boolean field
support and the child value `rn/2`.  The second identity in (TC.9) gives

```math
\left\|\sum_i\epsilon_iw_i\right\|_1
={n\over2}\left(3+\sum_{i<j}G_{ij}\epsilon_i\epsilon_j\right). \tag{TC.11}
```

The separate child and field bounds prove the first formula in (TC.8).
For the sphere,

```math
\left\|\sum_i\epsilon_iw_i\right\|_2^2
=n\left(3+2\sum_{i<j}G_{ij}\epsilon_i\epsilon_j\right),       \tag{TC.12}
```

and the normalized sum is again a positive top eigenvector.  It therefore
attains the Euclidean field and child bounds simultaneously, proving the
second formula. `square`

The new hypothesis is genuinely algebraic: it asks that the Boolean poles
be closed under the one higher-order character needed by three-input
majority.  It does not store a Boolean response table or optimize over old
spins.  It is not a theorem from `(G,R)` alone; rather, it identifies a
strict structured class on which `(G,R)` is complete.

### A scalable nontrivial closed family

Let `H_16` be the regularized order-16 Walsh matrix defined in Section 4.
It has positive top eigenvectors

```text
w1 = +-----+--+-----+
w2 = +--+--------+--+
w3 = +++-+-++++-+-+++
t  = +++++--++--+++++
```

with `t=w1 odot w2 odot w3`, and normalized Gram matrix

```math
G=R=
\begin{pmatrix}
1&1/2&0\\
1/2&1&-1/2\\
0&-1/2&1
\end{pmatrix}.                                      \tag{TC.13}
```

Here `T(G)=1`, so

```math
\mathcal B_m={rn\over2}+2mn,
\qquad
\mathcal S_m={rn\over2}+\sqrt5,mn.                \tag{TC.14}
```

Tensoring `H_16`, all four displayed vectors in the first factor, and the
all-one positive top pole in every later factor preserves (TC.7) and
(TC.13) at orders `16^j`.  Thus the closure is not a one-off finite
coincidence.  For `m=r/4` it operates at bounded total port mass
`3m/r=3/4`, although its spherical relaxation still has the fixed normalized
gap `(sqrt(5)-2)/4`.

## 4. Four ports: a common-Hadamard equal-`(G,R)` collision

Index the order-16 Walsh matrix by `F_2^4`, put

```math
b(u,v)=(-1)^{u\cdot v},\qquad u,v\in\mathbb F_2^2,
\qquad H_{16}=D_bW_{16}D_b.                         \tag{TC.15}
```

Then `H_16` is a symmetric regular Hadamard matrix,

```math
H_{16}^2=16I,\qquad H_{16}\mathbf1=4\mathbf1,
\qquad \operatorname {tr}H_{16}=0.                 \tag{TC.16}
```

Consider the five Boolean `+4` eigenvectors

```text
w0 = ----------------
w1 = -----++--++-----
w2 = ---+--+-+-++-+++
w4 = ---+++-+-+---+++
w5 = --+----+-++++-++
```

and the two four-port tuples

```math
P^+=(w_0,w_1,w_2,w_4),
\qquad P^-=(w_0,w_1,w_2,w_5).                       \tag{TC.17}
```

### Theorem TC.3 (scalable regular-Hadamard Gram collision)

The two tuples have exactly the same normalized Gram--Rayleigh data,

```math
G^+=G^-=R^+=R^-=
\begin{pmatrix}
1&1/2&0&0\\
1/2&1&0&0\\
0&0&1&0\\
0&0&0&1
\end{pmatrix},                                      \tag{TC.18}
```

but their joint Boolean field supports are

```math
J(P):=\max_{x\in\{+-1\}^{16}}\sum_{i=1}^4|w_i^Tx|,
\qquad J(P^+)=32,quad J(P^-)=28.                  \tag{TC.19}
```

Moreover this collision has a leading-scale exact-sign lift.  For `j>=1`,
put

```math
H_j=H_{16}^{\otimes j},\quad n_j=16^j,\quad r_j=4^j,
\quad
w_{i,j}=w_i\otimes\mathbf1_{16^{j-1}}.             \tag{TC.20}
```

Let `P_j^+` and `P_j^-` be the corresponding tuples and take `m_j=r_j`.
Then they still have identical `(G,R)`, while

```math
\boxed{
\mathcal B_{m_j}(H_j;P_j^+)
-\mathcal B_{m_j}(H_j;P_j^-)
\ge {1\over8}r_jn_j.}                              \tag{TC.21}
```

Thus Gram--Rayleigh data does not determine the exact Boolean trust
response, even for a common regular Hadamard child, Boolean top-eigenvector
ports, fixed arity four, and bounded total port mass `4m_j/r_j=4`.

#### Proof

Direct multiplication proves the eigenvector claims and (TC.18).  By
Boolean duality,

```math
J(P)=\max_{\epsilon\in\{+-1\}^4}
       \left\|\sum_i\epsilon_iw_i\right\|_1.        \tag{TC.22}
```

The 16 endpoint values for `P^+` have multiset

```math
\{16^{(4)},24^{(8)},32^{(4)}\},                    \tag{TC.23}
```

whereas those for `P^-` have multiset

```math
\{20^{(8)},28^{(8)}\}.                             \tag{TC.24}
```

This proves (TC.19) by a finite exact certificate, with every operation an
integer sum.  If `N=16^{j-1}`, grouping an arbitrary Boolean vector into
16 fibres gives fibre sums in `[-N,N]`.  The joint support is a convex
function of these sums, hence its maximum over that box occurs at a vertex.
The vertices are realized by fibre-constant Boolean vectors.  Consequently

```math
J(P_j^+)=2n_j,
\qquad J(P_j^-)={7\over4}n_j.                       \tag{TC.25}
```

For `P^+`, the base Boolean word

```text
x+ = ------+---+-----
```

has port correlations `(12,12,4,-4)`, hence joint field support 32, and
`|x_+^TH_16x_+|/2=24`.  Tensoring it with the all-one pole gives

```math
\mathcal B_m(H_j;P_j^+)
\ge {3\over8}r_jn_j+2mn_j.                          \tag{TC.26}
```

For the other tuple, the Rayleigh and field-support bounds give

```math
\mathcal B_m(H_j;P_j^-)
\le {1\over2}r_jn_j+{7\over4}mn_j.                 \tag{TC.27}
```

Subtract (TC.27) from (TC.26) and set `m=r_j` to obtain
(TC.21). `square`

At the finite seed, exhaustive Boolean evaluation additionally gives

```math
\begin{array}{c|ccc}
m&1&2&4\\ \hline
\mathcal B_m(H_{16};P^+)&56&88&152\\
\mathcal B_m(H_{16};P^-)&56&82&138.
\end{array}                                        \tag{TC.28}
```

These exact seed values are not used in the tensor lower bound.

The collision also locates the extra coordinate in the exposed-flatness
carrier.  The spherical-maximizing endpoint words are exactly those with
the first two endpoint signs equal, because their squared field norm is
`5n`.  Among these channels the best normalized `l_1` ratios are

```math
{1\over n}\left\|{\sum_i\epsilon_iw_i^+\over\sqrt5}\right\|_1
={2\over\sqrt5},
\qquad
{1\over n}\left\|{\sum_i\epsilon_iw_i^-\over\sqrt5}\right\|_1
={7\over4\sqrt5}.                                  \tag{TC.29}
```

Thus the best exposed flatness deficits are respectively
`1-2/sqrt(5)` and `1-7/(4sqrt(5))`.  Equal `(G,R)` determines the common
spherical value but not the `L=||u||_1` coordinate of the exposed carrier.
That carrier correctly distinguishes the systems; its extra coordinate is
precisely higher-order row-pattern information absent from the Gram state.

## 5. Exact-sign completion and the structural boundary

Each tuple represents four repeated shores of width `m_j`.  Fill all pairs
among the `4m_j` auxiliary vertices by the same arbitrary hollow signing in
both systems.  Each completed cap differs from its incomplete trust response
by at most

```math
{4m_j\choose2}=O(r_j^2)=O(n_j).                    \tag{TC.30}
```

Therefore the completed exact-sign caps remain separated by
`r_jn_j/8-O(n_j)=Theta(n_j^(3/2))`.  The collision is not an artifact of a
weighted or incomplete parent.

The three-port theorem and the four-port collision isolate the first hidden
coordinate.  Three-input majority has a Walsh expansion using the ports and
their triple product, while its absolute support has degree at most two and
is therefore Gram-visible.  At four ports the unresolved fourth-order row
pattern changes `l_1` support by a fixed fraction.  A sufficient carrier
beyond three ports must retain an appropriate higher-order product/row-type
algebra, or prove a synchronization law making those higher features
functions of the pairwise state.  Pairwise `(G,R)` alone cannot do so.

This is consistent with exposed-flatness recovery: that theorem controls a
specific spherical optimizer after one knows its coordinate geometry.  The
two tuples above have identical spherical Gram response, but different
Boolean joint support because `(G,R)` omits the fourth-order coordinate
pattern.

## 6. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_three_port_gram_closure.py
```

The verifier checks the two- and three-port formulas by exhaustive cube
evaluation, checks the product-closed triple, certifies every endpoint sum
and trust value in the four-port collision, and verifies the tensor support,
eigenvector, Gram, and leading-separation identities at the next order.
