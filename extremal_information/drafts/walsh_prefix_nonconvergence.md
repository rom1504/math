# A single infinite Walsh signing has two separated prefix phases

Status: rigorous near-original benchmark theorem. This is a nonconvergence
result for one explicit all-order sequence of dense hollow signings, **not**
for the minimizing values `M_n`.

## 1. The coherent infinite signing

Let

```math
H=
\begin{pmatrix}
1&1&1&1\\
1&-1&1&-1\\
1&1&-1&-1\\
1&-1&-1&1
\end{pmatrix},
\qquad
u=(1,1,1,-1)^T.                                      \tag{WP.1}
```

Then `H^2=4I` and `Hu=2u`. Put `H_r=H^(tensor r)`. Since the top-left
entry of `H` is one, `H_r` is the top-left principal block of `H_(r+1)`.
There is therefore one infinite symmetric sign matrix `S` whose leading
`4^r by 4^r` block is `H_r` for every `r`.

Let `C_n` be the leading `n by n` principal block of `S`, and let

```math
A_n=C_n^circ
```

be obtained by zeroing its diagonal. Each `A_n` is a valid dense hollow
signing. Define its normalized Boolean energy

```math
q_n={1\over2n^(3/2)}\max_(x in {+-1}^n)|x^TA_nx|
={Q(A_n)\over n^(3/2)}.                               \tag{WP.2}
```

## 2. Two persistent phases

Let `B` be the leading `3 by 3` block of `H`:

```math
B=\begin{pmatrix}1&1&1\\1&-1&1\\1&1&-1\end{pmatrix}. \tag{WP.3}
```

### Theorem WP.1 (Walsh-prefix nonconvergence)

The sequence `(q_n)` does not converge. More precisely,

```math
q_(4^r)={1\over2}\quad(r>=1),                         \tag{WP.4}
```

whereas

```math
q_(3*4^r)>= {89\over96sqrt3}
=0.535251812061\ldots\quad(r>=2).                    \tag{WP.5}
```

Consequently

```math
liminf_n q_n<=1/2
< {89\over96sqrt3}
<=limsup_n q_n.                                       \tag{WP.6}
```

#### Proof

At order `n=4^r`, the full matrix is `C_n=H_r`. Its operator norm is
`sqrt(n)`, while `u^(tensor r)` is a Boolean eigenvector of eigenvalue
`sqrt(n)`. Hence

```math
\max_x|x^TH_rx|=n^(3/2).                              \tag{WP.7}
```

Also `tr(H_r)=tr(H)^r=0` for `r>=1`, so hollowing does not change any
Boolean quadratic value. This proves (WP.4).

At order `3*4^r`, the coherent principal block is

```math
C_(3*4^r)=B tensor H_r.                               \tag{WP.8}
```

At `r=2`, the following three length-16 blocks form a Boolean vector `z`:

```text
(+ + + -  + + - +  + - + +  - + + -)
(+ - - -  + + + -  + - + +  + + - -)
(+ + + +  - + - +  - + + -  - - + -)
```

Direct integer multiplication gives

```math
z^T(B tensor H_2)z=356.                               \tag{WP.9}
```

For `r=2+s`, tensor `z` with `u^(tensor s)`. This remains Boolean and the
regular-Hadamard identity gives

```math
(z tensor u^(tensor s))^T(B tensor H_(2+s))
 (z tensor u^(tensor s))
=356*8^s.                                             \tag{WP.10}
```

Since `(3*4^(2+s))^(3/2)=48^(3/2)8^s`, (WP.9)--(WP.10) give

```math
q_(3*4^r)>= {356\over2*48^(3/2)}
={89\over96sqrt3}.                                   \tag{WP.11}
```

Finally `tr(B tensor H_r)=tr(B)tr(H_r)=0` for `r>=1`, so hollowing again
changes nothing. The strict inequality in (WP.6) is equivalent to
`89>48sqrt(3)`. `square`

## 3. The complete continuous phase profile

The two subsequences are part of a deterministic log-periodic limit law.
For `R_r=4^r` and `1<=t<=4`, define

```math
F_r(t)={Q(A_(floor(tR_r)))\over R_r^(3/2)}.            \tag{WP.12}
```

### Theorem WP.2 (continuous Walsh prefix-phase law)

The functions `F_r` converge uniformly on `[1,4]` to a continuous
nondecreasing function `F`. Consequently, for every `t in [1,4]`,

```math
{Q(A_(floor(t4^r)))\over floor(t4^r)^(3/2)}
\longrightarrow L(t):={F(t)\over t^(3/2)}.            \tag{WP.13}
```

The normalized phase profile `L` is continuous, satisfies

```math
L(1)=L(4)=1/2,
\qquad L(3)>=89/(96sqrt3)>1/2,                        \tag{WP.14}
```

and is therefore nonconstant.

#### Proof

First consider the dense set of base-four rationals

```math
t={p\over4^k}\in[1,4].                                \tag{WP.15}
```

For `r>=k`, put `n=t4^r=p4^(r-k)`. Associativity of the Kronecker powers
shows that the leading `n by n` block is

```math
C_n=B_(p,k) tensor H_(r-k),                            \tag{WP.16}
```

where `B_(p,k)` is the leading `p by p` block of `H_(k+1)`. The outer
template is fixed as `r` grows. Regular-Hadamard amplification, with the
`O(n^(-1/2))` diagonal-removal error, proves convergence of the normalized
Boolean maximum and hence of `F_r(t)` at every `t` in (WP.15).

It remains to make this dense-set convergence uniform. Let

```math
1<=t<=s<=4,
\qquad n=floor(tR_r),\quad m=floor(sR_r),\quad h=m-n.
```

Write the hollow prefix at order `m` in blocks as

```math
A_m=\begin{pmatrix}A_n&R\\R^T&D\end{pmatrix}.         \tag{WP.17}
```

Both `R` and the unhollowed matrix underlying `D` are submatrices of
`H_(r+1)`, whose operator norm is `2sqrt(R_r)`. Therefore

```math
||R||<=2sqrt(R_r),
\qquad ||D||<=2sqrt(R_r)+1.                           \tag{WP.18}
```

For Boolean `x,y`, the newly exposed energy obeys

```math
|x^TRy+{1\over2}y^TDy|
<=2sqrt(R_rnh)+{h\over2}(2sqrt(R_r)+1).               \tag{WP.19}
```

Principal deletion is lossless, `Q(A_n)<=Q(A_m)`: average the extended
quadratic over independent random missing spins and choose the appropriate
sign. Combining this with (WP.19) gives

```math
0<=F_r(s)-F_r(t)
<=2sqrt{(n/R_r)(h/R_r)}+h/R_r+{h\over2R_r^(3/2)}.     \tag{WP.20}
```

The right side is bounded uniformly by a modulus tending to zero with
`s-t` (plus `o_r(1)`). Thus `(F_r)` is asymptotically equicontinuous and
uniformly bounded. Convergence on the dense set (WP.15), followed by a
finite-net argument, makes it uniformly Cauchy and yields a continuous
limit `F`. Monotonicity follows from principal deletion. Dividing by
`t^(3/2)` proves (WP.13). Equations (WP.4)--(WP.5) give (WP.14). `square`

Equivalently, if

```math
t_n={n\over4^(floor(log_4 n))}\in[1,4),
```

then `q_n-L(t_n)->0`. The obstruction is a continuous scale phase, not an
isolated arithmetic anomaly.

## 4. What the example teaches

Every fixed prefix phase `d*4^r`, with `d` fixed, has a limiting response by
regular-Hadamard amplification. The failure is **between** phases: the
all-order direct limit remembers the leading base-four prefix, and different
outer response templates have separated limits. Thus subsequential carrier
compactness plus exact realization at every order does not imply that the
realizers approach one common extremal state.

The state needed by this automatic hierarchy is small: a finite outer
cross-correlation carrier together with the scale phase. Yet the phase does
not synchronize. This is a concrete all-order realization obstruction that
does not come from parity, unavailable Hadamard orders, or sampling noise.

It does not imply anything about nonconvergence of `M_n`. Minimization may
choose unrelated signings at the two order classes, and proving
nonconvergence of `M_n` would require universal lower and constructive upper
separation. The result instead validates the theory's distinction between
compact response images and a reusable all-order congruence.

## 5. Reproducibility

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_walsh_prefix_nonconvergence.py
```

The verifier checks nesting, Hadamard regularity, the exact 48-coordinate
certificate (WP.9), the tensor amplification identity through two further
levels, trace removal, and the strict rational-radical gap.
