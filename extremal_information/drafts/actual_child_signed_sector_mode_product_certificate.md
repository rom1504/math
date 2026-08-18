# Signed sector modes as product-phase certificates

Status: **rigorous actual-child theorem and scope audit**.  The theorem below
does not replace the balanced product-phase lemma.  It identifies the extra
signed information which the scalar sector--Gram mass discards, gives an exact
reverse-product certificate in the spectrally stable regime, and gives a
one-parameter coherent product retuning certificate in the unstable regime.
The latter can be tested from a polynomial-size sector mode and scalar
expectations along one declared product path; it does not solve the full
row-product variational problem.

## 1. Exact product-gain identity

Put two contracted-temperature optimizing children in the balanced row
direction and orientation of Theorem 37.32.  Let

```math
r=r_{\rm row}^{\otimes m},
\qquad {dq\over dr}={e^{-\lambda h}\over E_re^{-\lambda h}},
```

where `h` is the exact full-versus-row-erased bridge interaction.  Write

```math
\mathcal J=D(r\Vert q),
\qquad
\mathcal I^{\leftarrow}
=\inf_{P=\otimes_iP_i}D(P\Vert q).
```

For a row product `P`, define its gain over the canonical product by

```math
\mathscr G_h(P)
=\lambda\{E_rh-E_Ph\}-D(P\Vert r).                 \tag{SM.1}
```

Then the Gibbs density gives the exact identity

```math
\boxed{
\mathcal J-\mathcal I^{\leftarrow}
=\sup_{P\ {\rm row\ product}}\mathscr G_h(P).}     \tag{SM.2}
```

Indeed

```math
D(P\Vert q)=D(P\Vert r)+\lambda E_Ph+\log E_re^{-\lambda h},
```

and subtracting this expression from its value at `P=r` proves (SM.2).
Unlike a coordinate best-response equation, (SM.1) gives a rigorous lower
certificate from any single declared product law.

## 2. The signed quadratic carrier

Let `t=beta/sqrt(N)`.  With the notation of Theorems 37.30--37.33, put

```math
H_2(B)=\sum_{i<k}B_i^{\mathsf T}\Gamma_{ik}B_k
      ={1\over2}B^{\mathsf T}MB,                    \tag{SM.3}
```

where `M` is the symmetric `mn by mn` block matrix with diagonal blocks zero
and off-diagonal blocks `M_(ik)=Gamma_(ik)`.  Thus

```math
\|M\|_F^2=2K_\epsilon.                              \tag{SM.4}
```

The matrix does not require an `O(N^4)` table.  If `C_A^a,C_D^a` are the
two sector correlation matrices and
`\widehat C_A^a=C_A^a-I_m`, then the exact sector disintegration gives

```math
\boxed{
M=\sum_{a=\pm1}\pi_a^\epsilon
  \widehat C_A^a\otimes C_D^{\epsilon a}.}          \tag{SM.5}
```

Consequently a matrix-vector product with `M` uses four child covariance
matrices, containing `O(m^2+n^2)` real entries.  A negative edge eigenvalue
can therefore be computed without storing the bridge pressure, an
external-field table, or even the `O(N^4)` entries of `M`.

Define the exact physical remainder, with no series assumption, by

```math
h=t^2H_2+R+c,
\qquad \Omega=\operatorname {osc}R.                 \tag{SM.6}
```

Under the convergence premise of Theorem 37.33,

```math
\Omega\le2\mathfrak C_{\ge4}(t).                   \tag{SM.7}
```

The point of retaining `M`, rather than only `K_epsilon`, is that product
retuning sees the **negative spectral edge**

```math
\rho_-(M)=\max\{0,-\lambda_{\min}(M)\},             \tag{SM.8}
```

whereas `K_epsilon` is only half the squared Frobenius norm.

## 3. A global stable-mode certificate for reverse-product information

The balanced optimizer theorem gives

```math
{dr_{\rm row}\over dU_n}\le e^{C_0},
\qquad C_0=\lambda(\beta^2/2+\log2),                \tag{SM.9}
```

and `r_row` is centrally symmetric.

**Theorem SM.1 (signed spectral stability).**  If

```math
\boxed{\lambda t^2\rho_-(M)\le e^{-C_0},}           \tag{SM.10}
```

then

```math
\boxed{
0\le\mathcal J-\mathcal I^{\leftarrow}
\le\lambda\Omega.}                                \tag{SM.11}
```

In particular, along any actual-minimizer subsequence on which
`Omega=o(N)` and `J>=eta N`, one has

```math
\mathcal I^{\leftarrow}\ge\eta N-o(N).             \tag{SM.12}
```

*Proof.*  Central symmetry and (SM.9) imply, exactly as in Lemma SP.2,

```math
\log E_{r_{\rm row}}e^{\langle v,B\rangle}
\le {e^{C_0}\over2}\|v\|_2^2.                     \tag{SM.13}
```

If `P_i` is any row law and `m_i=E_(P_i)B_i`, entropy duality applied to
(SM.13) gives

```math
D(P_i\Vert r_{\rm row})
\ge {e^{-C_0}\over2}\|m_i\|_2^2.                  \tag{SM.14}
```

For `P=\otimes_iP_i`, independence between distinct rows and the zero
diagonal blocks of `M` give

```math
E_PH_2={1\over2}m^{\mathsf T}Mm,
\qquad E_rH_2=0.                                   \tag{SM.15}
```

Substitute (SM.6), (SM.14), and (SM.15) into (SM.1):

```math
\mathscr G_h(P)
\le {1\over2}\{\lambda t^2\rho_-(M)-e^{-C_0}\}
       \|m\|_2^2+\lambda\Omega.                   \tag{SM.16}
```

Under (SM.10) this is at most `lambda Omega`.  Take the supremum in
(SM.2).  Nonnegativity follows by taking `P=r`. `square`

This is a global product certificate.  It does not assume that `r` is a
coordinate fixed point, and it covers arbitrary within-row changes, not
only linear-field tilts.

## 4. One coherent binary retuning direction

There is a complementary certificate which never optimizes over all row
factors.  For a set of rows `S`, choose odd binary features

```math
\phi_i:\{-1,1\}^n\longrightarrow\{-1,1\},
\qquad \phi_i(-b)=-\phi_i(b),                       \tag{SM.17}
```

and put

```math
w_i=E_{r_{\rm row}}[B\phi_i(B)].                    \tag{SM.18}
```

Set `w_i=0` outside `S`, concatenate the `w_i` into `w`, and, for one common
amplitude `a`, define the product path

```math
{dP_{i,a}\over dr_{\rm row}}(b)
=\begin{cases}
 e^{a\phi_i(b)}/\cosh a,&i\in S,\\
 1,&i\notin S.
\end{cases}                                         \tag{SM.19}
```

This is one scalar path, not `m` coupled best responses.  Since an odd bit
is fair under a centrally symmetric law,

```math
E_{P_{i,a}}B=\tanh(a)w_i,
\qquad
D(P_{i,a}\Vert r_{\rm row})
=d(a):=a\tanh a-\log\cosh a.                       \tag{SM.20}
```

**Theorem SM.2 (binary sector-mode retuning certificate).**  If `k=|S|`,
then the exact quadratic score of (SM.19) is

```math
\mathscr G_{t^2H_2}(P_a)
=-{\lambda t^2\over2}\tanh^2(a)w^{\mathsf T}Mw
 -k d(a).                                           \tag{SM.21}
```

For the full actual interaction,

```math
\boxed{
\mathscr G_h(P_a)
\ge-{\lambda t^2\over2}\tanh^2(a)w^{\mathsf T}Mw
    -k d(a)-\lambda\Omega.}                        \tag{SM.22}
```

Suppose `0<delta<=1` and

```math
-\lambda t^2w^{\mathsf T}Mw\ge(1+\delta)k.         \tag{SM.23}
```

With

```math
a^2={3\delta\over4(1+\delta)},                     \tag{SM.24}
```

one obtains

```math
\boxed{
\mathcal J-\mathcal I^{\leftarrow}
\ge {3\delta^2\over16(1+\delta)}k-\lambda\Omega.} \tag{SM.25}
```

Thus, if `k>=alpha N` and `Omega=o(N)`, a rounded negative sector mode
satisfying (SM.23) proves the coherent-retuning branch
`J-I^leftarrow=Omega(N)`.  By Theorem 37.27, the resulting linear optimal
gain also forces order-one regular retuning on a positive density of rows.

*Proof.*  Equations (SM.20), row independence, and (SM.15) prove (SM.21).
The remainder expectation changes by at most `Omega`, proving (SM.22).
For `0<=a<=1`,

```math
d(a)\le a^2/2,
\qquad \tanh^2a\ge a^2(1-2a^2/3).                  \tag{SM.26}
```

The first inequality follows by integrating
`d'(a)=a sech^2(a)<=a`; the second follows from
`tanh a>=a-a^3/3`.  Use (SM.23) in (SM.21), substitute (SM.24), and obtain

```math
\mathscr G_{t^2H_2}(P_a)
\ge {a^2k\over2}
 \{(1+\delta)(1-2a^2/3)-1\}
={3\delta^2\over16(1+\delta)}k.
```

Now apply (SM.2) and (SM.22). `square`

A canonical way to declare the features is to take a negative eigenvector
`v` of the matrix in (SM.5), split it into row blocks, and use

```math
\phi_i(b)=\operatorname {sgn}_*\langle v_i,b\rangle, \tag{SM.27}
```

where ties are resolved by any fixed odd rule.  Condition (SM.23) is a
single signed Rayleigh check after this hyperplane rounding.  It is not
implied by the Frobenius mass (SM.4); that is precisely the information
lost by the sector--Gram scalar.  Such a hyperplane bit can be genuinely
aggregate and high-transport, but that geometric property is not automatic
and is not needed for the product-trial certificate.

## 5. Exact directional observable without a cluster premise

The remainder is unnecessary if one evaluates the actual interaction along
the declared path.  Define

```math
\boxed{
\mathscr D_{\phi}(a)
=\lambda\{E_rh-E_{P_a}h\}-k d(a).}                 \tag{SM.28}
```

Then (SM.1)--(SM.2) give, with no cumulant convergence or approximation,

```math
\boxed{
\mathscr D_{\phi}(a)\ge cN
\quad\Longrightarrow\quad
\mathcal J-\mathcal I^{\leftarrow}\ge cN.}        \tag{SM.29}
```

This observable is operationally smaller than the forbidden product
best-response oracle:

1. its direction is specified by `O(mn)` mode coefficients and one odd bit
   per row;
2. it has one scalar parameter and queries only two scalar expectations;
3. it never asks for an effective row potential as a function of an
   arbitrary point mass, so its queries cannot reconstruct `B -> h(B)`;
4. exact one-row preprocessing uses at most `2^n` row likelihood values,
   rather than `2^(mn)` bridge values.

There is also a sampling formulation.  Both `r` and `P_a` are row products.
Changing a complete bridge row changes `h` by at most `4tn`, because each
of its `n` bit changes moves `log p` and its erased-row term by at most
`2t` each.  Hence bounded differences under either product law gives a
variance proxy at most

```math
{1\over8}m(4tn)^2=2mt^2n^2=O_{\beta}(N^2)           \tag{SM.30}
```

at comparable splits.  Under sampling access to the balanced one-row law,
`O_(beta)(epsilon^(-2)log(1/zeta))` independent bridge point evaluations
estimate each expectation in (SM.28) to error `epsilon N` with failure
probability at most `zeta`.  Thus a fixed-density margin in (SM.29) has a
finite-query statistical certificate; no Gibbs table is hidden in the
definition.

## 6. What this does and does not reset

The signed negative edge and the rounded-mode profile are genuinely new
information beyond the four-coordinate sector--Gram state:

```text
sector--Gram mass K = half the Frobenius mass of M;
reverse-product stability = the negative spectral edge of M;
coherent retuning = a signed rounded Rayleigh profile of M and the actual h.
```

The child covariance matrices in (SM.5) are polynomial-size and are
strictly smaller than the bridge/Gibbs landscape.  Theorem SM.1 therefore
decides the reverse-dependence branch whenever the physical remainder is
sublinear and the negative edge is below the explicit threshold.  Theorem
SM.2 decides the retuning branch when one extensive rounded mode crosses the
explicit threshold, and (SM.28) is an exact low-query version which does not
need the remainder bound.

What is not proved is equally important.  Neither linear Frobenius mass nor
linear **absolute** connected-cluster mass forces (SM.10), (SM.23), or a
linear value of (SM.28).  Absolute cluster mass has no sign, and a Frobenius
norm has no spectral-edge or eigenvector information.  No theorem here says
that every positive-density optimal retuning is caught by the single
sector-mode rounding.  Therefore the result supplies concrete directional
certificates but does not verify one of them asymptotically for all actual
minimizers and does not remove `L_balanced-product-phase`.
