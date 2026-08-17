# An off-diagonal Walsh Gram bit is visible to an unrooted scalar maximum

Status: rigorous task-local theorem with an exact finite verifier.  The
result answers a minimality question left open by the unrooted ambient-orbit
collapse: although the characteristic-root fibre disappears, the binary
Gram form does not disappear altogether.  No rooted field, external pole,
or arbitrary lookup query is used.

## 1. The two algebraic states

Put

```math
q=2^m,\qquad n=q^2,
```

and let `F_E` be the normalized Walsh involution on
`E=F_2^m\oplus F_2^m`.  For `a in F_2^m`, write

```math
\ell_a=(0,a),\qquad
D_a(z)=(-1)^{\ell_a\cdot z},\qquad
C_a=D_aF_ED_a.                                           \tag{SG.1}
```

Thus `F_E` and every `C_a` are symmetric orthogonal involutions.  The
corresponding unrooted scalar query on the ordinary triangle is

```math
\mathcal E_{\mathbf a}(x_1,x_2,x_3)
 ={q\over2}\sum_{i=1}^3x_i^TC_{a_i}x_i
 +q\sum_{1\le i<j\le3}x_i^TF_Ex_j.                     \tag{SG.2}
```

Every `x_i` is Boolean of length `n`.  This is the unweighted `K_3` member
of the declared unrooted Walsh-graph language.

For every `m>=5`, use

```math
\begin{aligned}
\mathbf a^0&=(e_1+e_2,\ e_3+e_4,\ e_1+e_2+e_3+e_4),\\
\mathbf a^1&=(e_1+e_2,\ e_1+e_3,\ e_2+e_3).
\end{aligned}                                           \tag{SG.3}
```

In each tuple all three labels are nonzero and even, the first two are
independent, and the third is their sum.  Hence both have the same relation
kernel

```math
\mathcal R=\{000,111\}                                  \tag{SG.4}
```

and the same three self-pairings, all zero.  Their off-diagonal Gram values
are different:

```math
a_i^0\cdot a_j^0=0,
\qquad
a_i^1\cdot a_j^1=1\quad(i\ne j).                       \tag{SG.5}
```

Because `m>=5`, neither two-dimensional label span contains
`omega=(1,...,1)`.  Thus both characteristic-root fibres are empty as well.
The only orbit datum being changed is the off-diagonal bilinear form.

## 2. The flux identity

Let `F` be the normalized order-`q` Walsh matrix, and on functions on
`F_2^m` let `M_a` and `T_a` denote modulation and translation.  The reduced
child involution is

```math
J_a=M_aFM_a=M_aT_aF.                                    \tag{SG.6}
```

If `a,b` are even and `c=a+b`, Walsh commutation gives

```math
J_aJ_bJ_c=(-1)^{a\cdot b}F.                             \tag{SG.7}
```

Indeed,

```math
J_aJ_b=(-1)^{a\cdot b}M_cT_c,
\qquad (M_cT_c)^2=I,
```

which proves (SG.7).  The identical identity holds for the full involutions
`C_a,C_b,C_c,F_E`.  Moreover all four involutions commute, because every
label is even.  Consequently every common eigenspace has signs

```math
F_E\mapsto f,\qquad C_{a_i}\mapsto\lambda_i,
\qquad
\lambda_1\lambda_2\lambda_3=(-1)^{a\cdot b}f.          \tag{SG.8}
```

The bit in (SG.7) is a scalar, gauge-invariant triangle flux.  It is not a
rooted coordinate.

## 3. A common Boolean section in the zero-flux state

We give a dimension-uniform construction rather than infer Boolean
saturation from a spectral eigenvector.

On `E` put

```math
x_0(u,v)=(-1)^{u\cdot v}.                               \tag{SG.9}
```

Then `F_Ex_0=x_0`.  If `p` is even and `f_p=(p,p)`, direct Walsh summation
also gives

```math
F_E(D_{f_p}x_0)=D_{f_p}x_0,
\qquad C_{f_p}x_0=x_0.                                  \tag{SG.10}
```

Choose independent even vectors `p_1=e_1+e_2` and
`p_2=e_1+e_3`.  The source pair

```math
(\ell_{e_1+e_2},\ell_{e_3+e_4})
```

and target pair `(f_{p_1},f_{p_2})` are two independent totally isotropic
pairs in `E`.  Both spans avoid the characteristic vector
`Omega_E=(omega,omega)`.  The characteristic-rooted Witt extension lemma
therefore supplies `O in O(E)` carrying the source pair to the target pair.
The induced coordinate permutation commutes with `F_E` and conjugates the
corresponding children.  Transporting `x_0` back through `O` gives one
Boolean vector `x` satisfying

```math
F_Ex=x,
\qquad C_{a_i^0}x=x\quad(i=1,2,3).                      \tag{SG.11}
```

The third identity follows either by transport or because the third label
is the sum of the first two.  This use of Witt extension is only an
existence proof for a Boolean section; it does not invoke a rooted query.

## 4. Exact separation

### Theorem SG.1 (off-diagonal Gram visibility)

For every `m>=5`, the tuples in (SG.3) have identical self-pairings,
relation kernel, and empty characteristic-root fibre, but the unrooted
unweighted triangle query satisfies

```math
\max_X\mathcal E_{\mathbf a^0}(X)
 ={9\over2}n^{3/2},                                    \tag{SG.12}
```

whereas

```math
\max_X\mathcal E_{\mathbf a^1}(X)
 \le {3(1+\sqrt{17})\over4}n^{3/2}.                   \tag{SG.13}
```

Thus the two scalar maxima differ by at least

```math
{3(5-\sqrt{17})\over4}n^{3/2}
=0.6576707807\ldots\,n^{3/2}.                          \tag{SG.14}
```

#### Proof

For the zero-flux tuple, put the vector from (SG.11) in all three blocks.
Every child contributes `qn/2` and every edge contributes `qn`.  The three
child and three edge operator-norm bounds give the reverse inequality, so
(SG.12) is exact.

For the unit-flux tuple, decompose the common internal space into joint
eigenspaces of `F_E,C_{a_1},C_{a_2},C_{a_3}`.  The normalized three-block
coefficient operator restricts on such an eigenspace to

```math
B=\operatorname{diag}(\lambda_1,\lambda_2,\lambda_3)
  +fA(K_3).                                             \tag{SG.15}
```

Multiplying by `f` does not change its norm.  If
`s_i=f\lambda_i`, then (SG.8) says `s_1s_2s_3=-1`.
When all three signs are negative, `||diag(s)+A(K_3)||=2`.  When exactly
one is negative, permutation symmetry reduces the matrix to

```math
\begin{pmatrix}-1&1&1\\1&1&1\\1&1&1\end{pmatrix},
```

whose eigenvalues are `0` and `(1+-sqrt(17))/2`.  Therefore the full
normalized block operator has norm at most `(1+sqrt(17))/2`.

Every Boolean three-block vector has squared Euclidean norm `3n`.
Multiplying the Rayleigh bound by the prefactor `q/2` in (SG.2) proves
(SG.13). `square`

## 5. Semantic consequence

The ambient collapse proves that `(Gram,relation kernel)` is sufficient for
the whole unrooted weighted graph landscape.  Theorem SG.1 proves a genuine
partial converse: after fixing self-pairings and the complete relation
kernel, an off-diagonal Gram/flux bit can still change one ordinary scalar
Boolean maximum at the full `n^(3/2)` scale.  Hence the minimal scalar
semantic quotient cannot in general collapse all of the Gram form.

The result does **not** prove that every individual Gram entry is separately
recoverable, nor that `(Gram,R)` is scalar-minimal for arbitrary tuples.
What the triangle exposes is the bilinear form on a rank-two presented span,
through the relation `a_1+a_2+a_3=0`.  This suggests that the semantic
resource is a relation-cycle flux rather than a list of raw pair entries.

## 6. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_walsh_scalar_gram_visibility.py
```

The verifier checks the two algebraic states, the Weyl product identity,
the exact common Boolean section at `m=5` using an explicit product of three
orthogonal transvections, the triangle spectra, and the stated constants.
