# Independent audit: robust product synchronization

**Status.**  Independent proof and adversarial audit.  The proposed robust
selector theorem and its tensor estimate are correct under the stated
operator-contraction hypothesis.  The potentially dangerous Schur multiplier
is indefinite, but it is nevertheless contractive on Hermitian matrices for
a structural reason recorded below.  The constants in both estimates are
first-order sharp.  No exact-sign or hollowness assumption is used by this
abstract theorem.

## 1. Setup and the exact defect identity

Let `H` be real symmetric, let `r>0`, and assume

```math
\|H\|_{op}\le r.                                      \tag{RPS.1}
```

Let `w_1,...,w_p in {+-1}^n`.  Fix an antipodally odd Boolean selector
`tau:{+-1}^p->{+-1}` which agrees with majority off its zero layer.  Write

```math
\tau(a)=\sum_{S\in\mathcal A}\widehat\tau(S)
                 \prod_{i\in S}a_i,                  \tag{RPS.2}
```

where `A` is its nonzero Walsh support; every `S in A` has odd cardinality.
Put `v_S=odot_(i in S)w_i`, and let `V` be the matrix with these product
ports as columns.  Define the normalized product-port Gram and Rayleigh
matrices

```math
G={1\over n}V^TV,
\qquad
R={1\over rn}V^THV,
\qquad
D=G-R.                                               \tag{RPS.3}
```

For an endpoint word `epsilon in {+-1}^p`, define

```math
c_\epsilon(S)=\widehat\tau(S)\prod_{i\in S}\epsilon_i,
\qquad
x_\epsilon=Vc_\epsilon.                             \tag{RPS.4}
```

Fourier inversion says that `x_epsilon` is the coordinatewise selector
`tau(epsilon_1w_1,...,epsilon_pw_p)`, and hence is Boolean.  Parseval gives

```math
\|c_\epsilon\|_2^2=1.                               \tag{RPS.5}
```

Since `||x_epsilon||_2^2=n`, one also has the exact, sometimes useful,
identity

```math
c_\epsilon^TGc_\epsilon=1.                          \tag{RPS.6}
```

Consequently

```math
rn-x_\epsilon^THx_\epsilon
=rn\,c_\epsilon^TDc_\epsilon.                     \tag{RPS.7}
```

This is the whole source of the robust response estimate; no product
channel is paid separately.

## 2. Positivity and robust Boolean recovery

### Proposition RPS.1 (positive selector defect)

Under (RPS.1),

```math
\boxed{D={1\over n}V^T(I-H/r)V\succeq0.}           \tag{RPS.8}
```

For `m>=0`, put `z_epsilon=sum_i epsilon_iw_i` and

```math
B_\epsilon=
\max_{x\in\{+-1\}^n,\ \sigma\in\{+-1\}}
\left\{{\sigma\over2}x^THx+mz_\epsilon^Tx\right\}.
                                                               \tag{RPS.9}
```

Then every endpoint obeys the sharper querywise estimate

```math
0\le {rn\over2}+m\|z_\epsilon\|_1-B_\epsilon
\le {rn\over2}c_\epsilon^TDc_\epsilon
\le {rn\over2}\|D\|_{op}.                         \tag{RPS.10}
```

#### Proof

The upper spectral bound `H<=rI` proves (RPS.8).  The norm assumption in
(RPS.1) bounds either signed quadratic channel by `rn/2`, while Holder
bounds the field by `m||z_epsilon||_1`; this proves the left inequality in
(RPS.10).  The majority property gives
`z_epsilon^Tx_epsilon=||z_epsilon||_1`, including ties because their field
coordinate is zero.  Evaluate the positive quadratic channel at this one
Boolean selector, and use (RPS.7), (RPS.5), and (RPS.8).  `square`

Only the upper bound `H<=rI` is needed for positivity and for the positive
witness.  The two-sided norm bound is needed for the displayed roof because
the response also maximizes the negative quadratic channel.

## 3. The indefinite contraction-kernel Schur lemma

The tempting tensor decomposition

```math
D_{12}=D_1\circ G_2+R_1\circ D_2                \tag{RPS.11}
```

contains an indefinite Schur symbol `R_1`.  Entrywise boundedness of `R_1`
would **not** by itself imply operator-norm contraction.  Here contraction
does hold, because `R_1` is a Rayleigh kernel dominated by its Gram kernel.

### Lemma RPS.2 (Hermitian Schur contraction)

Let `G` be a real correlation matrix and let `R` be symmetric with

```math
G+R\succeq0,
\qquad G-R\succeq0.                                \tag{RPS.12}
```

Then for every real symmetric `X`,

```math
\boxed{\|R\circ X\|_{op}\le\|X\|_{op}.}          \tag{RPS.13}
```

#### Proof

Put `P=(G+R)/2` and `N=(G-R)/2`.  The Schur maps `Phi_P` and `Phi_N` are
positive, and their sum `Phi_G` is unital because `diag(G)=1`.  If
`-I<=X<=I`, positivity gives

```math
-\Phi_P(I)\preceq\Phi_P(X)\preceq\Phi_P(I),
\qquad
-\Phi_N(I)\preceq\Phi_N(X)\preceq\Phi_N(I).
```

Subtracting the two maps and using
`Phi_P(I)+Phi_N(I)=Phi_G(I)=I` yields

```math
-I\preceq (P-N)\circ X=R\circ X\preceq I.
```

Scale by `||X||_op`.  `square`

For (RPS.3), condition (RPS.12) follows directly from

```math
G\mathbin\pm R={1\over n}V^T(I\mathbin\pm H/r)V\succeq0.   \tag{RPS.14}
```

Thus the indefinite step in (RPS.11) is valid.  This lemma is restricted to
Hermitian inputs; that is all the tensor proof requires.  One should not
replace (RPS.12) merely by `|R_ST|<=1`.

## 4. Tensor subadditivity

Consider two systems with the same selector and active-set labels.  Let
their data be `(H_i,r_i,n_i,V_i,G_i,R_i,D_i)`.  Give the tensor child
`H_1 tensor H_2` the corresponding ports `w_j^(1) tensor w_j^(2)`.  Its
active product-port columns are `v_S^(1) tensor v_S^(2)`.  Therefore

```math
G_{12}=G_1\circ G_2,
\qquad R_{12}=R_1\circ R_2,                        \tag{RPS.15}
```

and

```math
\begin{aligned}
D_{12}
 &=D_1\circ G_2+R_1\circ D_2                       \tag{RPS.16a}\\
 &=D_1\circ G_2+G_1\circ D_2-D_1\circ D_2.         \tag{RPS.16b}
\end{aligned}
```

### Theorem RPS.3 (robust product synchronization)

Under (RPS.1) in both factors,

```math
D_{12}\succeq0,
\qquad
\boxed{\|D_{12}\|_{op}
       \le\|D_1\|_{op}+\|D_2\|_{op}.}             \tag{RPS.17}
```

Consequently an `L`-factor tensor has selector defect at most the sum of
the `L` factor defects, and Proposition RPS.1 gives the corresponding
joint Boolean-response guarantee at every endpoint.

#### First proof

Positivity is (RPS.8) applied to the tensor contraction.  Schur
multiplication by a correlation matrix is a unital positive map and hence
contracts Hermitian operator norm.  Apply this to the first term of
(RPS.16a), and apply Lemma RPS.2 to the second:

```math
\|D_{12}\|_{op}
\le\|D_1\circ G_2\|_{op}+\|R_1\circ D_2\|_{op}
\le\|D_1\|_{op}+\|D_2\|_{op}.                    \tag{RPS.18}
```

#### Second proof, avoiding the indefinite multiplier

The Schur product theorem makes all three matrices
`D_1 circ G_2`, `G_1 circ D_2`, and `D_1 circ D_2` positive semidefinite.
Equation (RPS.16b) therefore gives the Loewner bound

```math
0\preceq D_{12}
\preceq D_1\circ G_2+G_1\circ D_2.                \tag{RPS.19}
```

Correlation Schur contraction on the two terms proves (RPS.17).  This
second proof is the safer route if one does not wish to invoke RPS.2.

Iteration of (RPS.17) proves the `L`-factor statement.  `square`

## 5. Sharpness and scope

### 5.1 The growing-arity boundary: raw operator defect is not stable

The fixed-label tensor theorem uses corresponding ports and hence Hadamard
products.  The PC.3 affine-coset construction is different: when two
product algebras are combined, **all pairs** of active products occur.  If
`V_i` lists those products, then

```math
V_{12}=V_1\otimes V_2,
\qquad
G_{12}=G_1\otimes G_2,
\qquad
R_{12}=R_1\otimes R_2.                              \tag{RPS.20}
```

Thus (RPS.17) does not extend with the raw norm `||G-R||_op`.

This failure is exact in the PC.3 regular-Hadamard seed.  In the order
`(a,b,c,abc)`, the four active poles have Gram matrix

```math
G_0=
\begin{pmatrix}
1&1/2&0&-1/2\\
1/2&1&-1/2&0\\
0&-1/2&1&1/2\\
-1/2&0&1/2&1
\end{pmatrix},
\qquad \|G_0\|_{op}=2.                              \tag{RPS.21}
```

Let `T_0=H_16/4`, let `U` be the span of these poles, and let `P_U` be its
orthogonal projection.  For `0<eta<=2`,

```math
T_\eta=T_0-\eta P_U                                \tag{RPS.22}
```

is still a symmetric contraction.  Its active state has
`R_0=(1-eta)G_0` and `D_0=eta G_0`.  Tensor it with `j-1` exact seed
factors.  The level-`j` PC.3 active algebra has

```math
G_j=G_0^{\otimes j},
\qquad D_j=\eta G_j,
\qquad
\boxed{\|D_j\|_{op}=\eta 2^j}.                    \tag{RPS.23}
```

Here the port arity is `p_j=2j+1`.  Exact factors of zero defect therefore
multiply the raw defect by two at each step.  This is a genuine
arity/dimension loss, not a loose estimate.

### 5.2 A dimension-free relative defect repairs the direct product

The raw norm is also stronger than the selector response needs.  Define the
intrinsic Gram-relative defect

```math
\delta(V,T)=
\sup_{c:\ c^TGc>0}{c^T(G-R)c\over c^TGc}
=1-\lambda_{min}(P_UT|_U),                         \tag{RPS.24}
```

where `T=H/r` and `U=range(V)`.  The equality follows by writing `u=Vc`;
it remains valid for redundant columns by using the generalized Rayleigh
quotient.  Since every selector coefficient obeys `c_epsilon^TG
c_epsilon=1`, Proposition RPS.1 sharpens to

```math
{rn\over2}+m\|z_\epsilon\|_1-B_\epsilon
\le {rn\over2}\delta(V,T).                         \tag{RPS.25}
```

### Proposition RPS.4 (Cartesian-product relative synchronization)

For the full Cartesian product algebra in (RPS.20),

```math
\boxed{\delta(V_1\otimes V_2,T_1\otimes T_2)
       \le\delta(V_1,T_1)+\delta(V_2,T_2).}         \tag{RPS.26}
```

If both compressed contractions `P_(U_i)T_i|_(U_i)` are positive
semidefinite, then the exact formula is

```math
\delta_{12}=\delta_1+\delta_2-\delta_1\delta_2.     \tag{RPS.27}
```

#### Proof

Let `A_i=P_(U_i)T_i|_(U_i)`.  The compression on `U_1 tensor U_2` is
exactly `A_1 tensor A_2`.  Every eigenvalue `a` of `A_1` lies in
`[1-delta_1,1]`, and similarly `b in [1-delta_2,1]`.  Hence

```math
1-ab=(1-a)+(1-b)-(1-a)(1-b)
\le\delta_1+\delta_2.                              \tag{RPS.28}
```

Taking the largest defect proves (RPS.26).  When both spectra are
nonnegative, the smallest product eigenvalue is
`(1-delta_1)(1-delta_2)`, proving (RPS.27).  `square`

The same upper bound holds if the next active span is only a subspace of
`U_1 tensor U_2`, as for a corresponding-port tensor: restriction cannot
increase the largest defect.  Equality of the active span with the full
tensor product is what gives the spectral calculation and exact formula
(RPS.27).  In the PC.3 generator enlargement, the affine coset in the new
factor is multiplied independently by the old coset, so its active columns
are exactly all `v_s tensor v_t`; hence this equality hypothesis holds.

For (RPS.22)--(RPS.23), the relative defect is exactly `eta` at every
level, and every selector coefficient has `c^TD_jc=eta`.  Thus PC.3 has no
actual response loss from growing arity; only the **wrong raw norm** has the
factor `2^j`.  The correct boundary is:

* corresponding-port tensors admit the ordinary defect norm RPS.3;
* expanding affine-coset product algebras require the generalized
  `D preceq delta G` defect RPS.4.

This relative certificate still stores the declared active product span;
it does not reduce an exponentially supported selector to pairwise port
data.  It removes an artificial conditioning loss once that product algebra
has already been selected.

### 5.3 The endpoint defect table is not a tensor congruence

One might try to compress `D` further to the finitely queried numbers

```math
q_\epsilon(D)=c_\epsilon^TDc_\epsilon.             \tag{RPS.29}
```

That table gives every one-step loss in RPS.1, but it is not reusable under
tensoring, already for three-port majority.

Take the four coordinates indexed by the active monomials
`(w_1,w_2,w_3,w_1w_2w_3)`.  Uniform projective three-bit rows make their
Gram matrix `G=I_4`.  Up to antipodes, the endpoint Fourier vectors are the
columns of the orthogonal matrix

```math
C={1\over2}
\begin{pmatrix}
1&1&1&1\\
1&1&-1&-1\\
1&-1&1&-1\\
-1&1&1&-1
\end{pmatrix}.                                     \tag{RPS.30}
```

Set

```math
D_0={1\over2}I_4,
\qquad
D_1=C\left[{1\over2}I_4+{1\over4}(E_{12}+E_{21})\right]C^T.
                                                               \tag{RPS.31}
```

Both are admissible Gram--Rayleigh defects: `D_i>=0` and
`2I-D_i>=0`, so `R_i=I-D_i` is the compression of a symmetric contraction.
Indeed, if `V/2` is the orthogonal normalized active-port matrix, take
`T_i=(V/2)R_i(V/2)^T`.  Every initial endpoint sees exactly the same defect,

```math
c_\epsilon^TD_0c_\epsilon
=c_\epsilon^TD_1c_\epsilon={1\over2}.              \tag{RPS.32}
```

For the corresponding-port self tensor, however,

```math
D_i^{(2)}=I-(I-D_i)\circ(I-D_i).                   \tag{RPS.33}
```

Direct calculation gives

```math
c_\epsilon^TD_0^{(2)}c_\epsilon={48\over64},
\qquad
c_\epsilon^TD_1^{(2)}c_\epsilon={47\over64}       \tag{RPS.34}
```

for every endpoint.  Thus identical one-step defect tables can evolve to
different tables.  The hidden datum is off-query coherence in `D`; the full
matrix law (RPS.16), or a genuine congruent quotient of it, cannot in
general be replaced by the exposed diagonal values (RPS.29).

More generally, with `D_0=aI` and
`C^TD_1C=aI+t(E_12+E_21)`, the initial tables agree, while the **average**
second-step table differs by `t^2/4`.  This survives throughout the open
admissible range `0<t<a` and `a+t<2`, so the example is not a boundary
degeneracy.

### 5.4 Sharp constants for fixed-label tensors

Both constants are first-order sharp already for one active port.  Let
`w in {+-1}^n`, `0<=d<=1`, and

```math
H=I-d{ww^T\over n},
\qquad r=1.                                         \tag{RPS.35}
```

Then `||H||_op=1`, `G=1`, `R=1-d`, and `D=d`.  If `m>=d`, direct
one-variable maximization in `t=w^Tx/n` shows that the positive selector
`x=w` is globally optimal in (RPS.9), and

```math
{n\over2}+mn-B={dn\over2}.                         \tag{RPS.36}
```

Thus the factor `1/2` in (RPS.10) cannot be reduced.  Tensoring two such
systems gives

```math
D_{12}=d_1+d_2-d_1d_2.                             \tag{RPS.37}
```

As `d_1,d_2` tend to zero, the ratio of (RPS.37) to `d_1+d_2` tends to
one, so the additive coefficient in (RPS.17) is also sharp.

The theorem controls a **declared selector** through the product ports in
its Fourier support.  It does not say that a small defect follows from the
original port Gram--Rayleigh table, and it does not optimize over selectors.
If the active support has exponential size, `D` can itself be an
exponentially large state.  Its content is instead a robust, jointly paid
version of exact product-algebra synchronization and a clean tensor law for
that certificate.

## 6. Adversarial checks

The accompanying verifier:

1. enumerates exact rational symmetric contractions obtained as averages of
   signed involutions and Boolean port tuples;
2. checks `G+-R>=0`, `D>=0`, the selector Fourier identities, and (RPS.10)
   by exhaustive Boolean maximization;
3. checks (RPS.15)--(RPS.17) on tensor pairs;
4. samples Hermitian inputs for the indefinite Schur contraction RPS.2; and
5. checks the sharp one-port family.

See
[`../experiments/verify_robust_product_synchronization.py`](../experiments/verify_robust_product_synchronization.py).
