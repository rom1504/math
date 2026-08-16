# Scalar-closed multichannel holonomy over finite fields

**Status.** The statements below are proved.  The normal form, local shear,
distance-transform bounds, and response geometry are exhaustively checked over
`F_3` in
[`verify_phase3_qary_multichannel_holonomy.py`](../experiments/verify_phase3_qary_multichannel_holonomy.py).

This note extends the binary multichannel packing to every finite field.  The
only new issue is that a quotient channel has several nonzero scalar letters.
Scalar closure does **not** create a cheaper hidden relation: after the letters
in each fibre are added, their two total coefficients must be opposite, so a
nonzero channel still costs exactly two letters.

## 1. Exact scalar-closed normal form

Let `q` be a prime power,

```math
W=\mathbb F_q^D,\qquad Q=\mathbb F_q^k,
```

and let `b_1,...,b_D` and `q_1,...,q_k` be their coordinate bases.  Give the
kernel the scalar-closed Hamming alphabet

```math
K^\times=\{(\alpha b_i,0):\alpha\in\mathbb F_q^\times,
                         1\le i\le D\}.
```

For an ordered independent tuple `V=(v_1,...,v_k)` in `W`, define

```math
P=K^\times\cup
  \{(0,\alpha q_j):\alpha\ne0,\ 1\le j\le k\},
```

```math
R_V=K^\times\cup
  \{(\alpha v_j,\alpha q_j):\alpha\ne0,\ 1\le j\le k\}.
                                                        \tag{QMP.1}
```

Word length is in the additive group `W\oplus Q`.  All alphabets are
symmetric because they are scalar closed.  Write `wt` for Hamming weight and
let `Vz=\sum_j z_jv_j`.

### Proposition QMP.1 (exact kernel profile and local gauge triviality)

For every `u\in W`,

```math
F_V(u):=\ell_{P\cup R_V}(u,0)
 =\min_{z\in\mathbb F_q^k}
   \bigl(2\operatorname{wt}(z)+\operatorname{wt}(u+Vz)\bigr). \tag{QMP.2}
```

Moreover each fragment separately is gauge-equivalent to `P` by an
automorphism fixing `W` pointwise.  Thus all dependence on `V` is created by
composition, not visible in either fragment alone.

#### Proof

Consider any word for `(u,0)`.  In quotient channel `j`, let `z_j` be the sum
of the coefficients of its `R_V` letters.  Quotient cancellation forces the
sum of the coefficients of its `P` letters to be `-z_j`.  If `z_j\ne0`, at
least one letter of each kind was used, at cost at least two.  If `z_j=0`,
all letters in that channel contribute zero in both `W` and `Q` after
consolidation and may be deleted.  The total `W` contribution of the retained
`R_V` letters is `Vz`; the kernel letters therefore cost at least
`wt(u-Vz)`.  Replacing `z` by `-z` gives the lower bound in (QMP.2).

Conversely, for every nonzero coordinate `z_j`, use the two letters
`(z_jv_j,z_jq_j)` and `(0,-z_jq_j)`, then represent `u-Vz` with one scalar
kernel letter per nonzero coordinate.  This attains the lower bound.

Finally let `L_V:Q\to W` be the linear map `L_V(q_j)=v_j`.  The shear

```math
(w,x)\longmapsto(w-L_Vx,x)
```

fixes `W` and sends `R_V` exactly to `P`.  This proves local triviality.
`\square`

The proof also pinpoints why using all nonzero scalar generators is harmless:
several letters in one fibre can only change their summed coefficient.  A
zero sum has zero holonomy; a nonzero sum can already be represented by one
letter.

## 2. Response geometry

Put `C_V=\operatorname{im}V` and use ambient Hamming distance.  Equation
(QMP.2) gives

```math
d(u,C_V)\le F_V(u)\le d(u,C_V)+2k.             \tag{QMP.3}
```

Indeed a nearest codeword has a coordinate vector with at most `k` nonzero
coordinates.  Hence for two `k`-subspaces `C,C'\le W`, equipped with arbitrary
ordered bases `V,V'`,

```math
\left|\,\|F_V-F_{V'}\|_\infty-d_H(C,C')\,\right|\le2k,       \tag{QMP.4}
```

where `d_H` is Hausdorff distance in the Hamming cube.  This uses the general
identity
`\|d(\cdot,C)-d(\cdot,C')\|_\infty=d_H(C,C')` and (QMP.3).

## 3. Macroscopic `q`-ary packing

### Lemma QMP.2 (uniform good linear host)

For every prime power `q` and all sufficiently large `D`, there is a
`q`-ary linear code `C_0\le\mathbb F_q^D` with

```math
\dim C_0=\lfloor D/4\rfloor,
\qquad d(C_0)>D/8.                              \tag{QMP.5}
```

#### Proof

Choose an `r=\lfloor D/4\rfloor` dimensional subspace uniformly.  A fixed
nonzero vector belongs to it with probability at most `q^{r-D+1}`.  Also

```math
\sum_{i\le D/8}{D\choose i}(q-1)^i
 \le 2^{H_2(1/8)D}q^{D/8}
 \le q^{(H_2(1/8)+1/8)D}.
```

Here the last inequality uses `\log_2q\ge1`, and
`H_2(1/8)+1/8=0.66856\ldots<3/4`.  The expected number of nonzero words of
weight at most `D/8` is therefore `o(1)`. `\square`

### Theorem QMP.3 (full macroscopic response rate over `F_q`)

For every prime power `q`, all sufficiently large `D`, and

```math
1\le k\le\lfloor D/32\rfloor,
```

there is a family `\mathcal V_{q,D,k}` of independent ordered `k`-tuples in
`\mathbb F_q^D` such that

```math
|\mathcal V_{q,D,k}|\ge q^{3Dk/16}             \tag{QMP.6}
```

and distinct members obey

```math
\boxed{\ \|F_V-F_{V'}\|_\infty>D/16.\ }        \tag{QMP.7}
```

Consequently any deterministic summary answering all kernel-endpoint queries
to error `\varepsilon D`, for fixed `\varepsilon<1/32`, has at least

```math
\frac{3}{16}Dk\log_2q                           \tag{QMP.8}
```

bits on this family.  If `V` is uniform on the family and a randomized
message reconstructs the whole profile to this accuracy with probability at
least `1-\eta`, then

```math
I(V;S)\ge(1-\eta)\frac{3}{16}Dk\log_2q-H_2(\eta). \tag{QMP.9}
```

#### Proof

Take `C_0` from Lemma QMP.2 and, for every `k`-subspace `C\le C_0`, choose
one ordered basis `V_C`.  The Gaussian binomial coefficient satisfies

```math
{r\brack k}_q
=\prod_{i=0}^{k-1}\frac{q^r-q^i}{q^k-q^i}
\ge q^{k(r-k)}\ge q^{3Dk/16}.                  \tag{QMP.10}
```

For distinct `C,C'`, choose `c\in C\setminus C'`.  Every `c-c'`, with
`c'\in C'`, is a nonzero word of `C_0`, so `d(c,C')>D/8`.  At query `c`,
(QMP.2)--(QMP.3) give `F_{V_C}(c)\le2k` and
`F_{V_{C'}}(c)>D/8`.  Since `2k\le D/16`, (QMP.7) follows.  The deterministic
claim is the packing argument, and the randomized claim is nearest-profile
decoding followed by Fano's inequality. `\square`

More generally, a `q`-ary `[D,r,d]` host yields

```math
\log_2|\mathcal V|\ge k(r-k)\log_2q,
\qquad
\|F_V-F_{V'}\|_\infty\ge d-2k.                \tag{QMP.11}
```

## 4. What this adds

The lower bound is not a binary parity artifact.  Composition exposes a
linear map in `\operatorname{Hom}(\mathbb F_q^k,\mathbb F_q^D)`, and a
constant fraction of its `Dk\log_2q` raw information remains visible at
macroscopic response accuracy.  The proof applies simultaneously to
scalar-closed Cayley landscapes and to `q`-ary systematic-code coset-leader
profiles.  It therefore supplies the same information-growth mechanism in
two model languages.

This does not claim that every multichannel family is incompressible.
Synchronization or special decoding structure may still give a strict
quotient.  It proves that field enlargement and unlabelled queries do not, by
themselves, remove the macroscopic mixed-holonomy charge.
