# Three independent Gauss phases give a genuine new Walsh catalyst family

Date: 2026-09-06. Exact arithmetic and asymptotic realization proved below.
A bounded label search does NOT improve the weighted seed's known value.
Neither `R=T` nor a strict `R<T` seed is proved.

This differs from the archived `7^e` construction: its available phases
are tied to one angle. Here three independent angles are realized by actual
even-dimensional Walsh transforms, not inserted as free numerical choices.

## 1. Base field and Gauss sums

Let `F32=F2[g]/(g^5+g^2+1)`, let `omega=exp(2pi i/31)`, and define

```math
G_j=\sum_{t=0}^{30}(-1)^{\operatorname{Tr}(g^t)}\omega^{jt},
\qquad j\ne0\pmod {31}.
```

The integer checker
`computations/transfer_adversary_gauss31_phase_2026_09_06.py`
verifies that the powers of `g` are all 31 nonzero elements and return to
one. Thus the quotient is a field and `g` is primitive; this is not an
assumed primitive-polynomial table. The additive character is nontrivial,
since its value at one is `(-1)^5=-1`.

The usual two-sum orthogonality proof gives `|G_j|^2=32`.
Substitution by squaring gives `G_(2j)=G_j`, and the real additive
character gives `G_(-j)=conjugate(G_j)`. The six Frobenius orbits have
representatives `1,3,5,-1,-3,-5`:

```text
  1: {1,2,4,8,16}       -1: {15,23,27,29,30}
  3: {3,6,12,17,24}     -3: {7,14,19,25,28}
  5: {5,9,10,18,20}     -5: {11,13,21,22,26}.
```

## 2. Exact independence, using only a small local-ring calculation

Put `L=Q2[alpha]/(alpha^5+alpha^2+1)`. The reduction modulo two is
the field just verified. On its integral ring the minimum 2-adic
valuation of the five coefficients is multiplicative: after dividing out
the minimum powers of two, the two nonzero residue classes have nonzero
product in `F32`. This gives a discrete valuation `v` with `v(2)=1`.

In the corresponding ring modulo 32, the checker computes

```math
z=\alpha^{32}=(4,9,26,20,2),\qquad z^{31}=1,\qquad z\bmod2=g.
```

The derivative of `T^31-1` at `z` is a unit. Successive Newton/Hensel
lifting therefore gives an actual root `zeta` congruent to `z` modulo32;
its reduction has order31, so `zeta` has order31. For each `k=1,3,5`
there is an embedding `omega -> zeta^k` of the cyclotomic field into `L`.
More explicitly, start with the displayed polynomial representative of
`z` in the complete finite free ring `Z2[alpha]` and iterate
`z <- z-(z^31-1)/(31 z^30)`. The initial residual is divisible by `2^5`,
the denominator is a unit, and the residual valuation at least doubles
at each step. This preserves the displayed residue modulo32. Since
`zeta != 1` and `zeta^31=1`, it satisfies the 31st cyclotomic polynomial;
substitution therefore defines the asserted injective field map from
`Q(omega)`, even though that field has degree30 over `Q` and `L` has
degree5 over `Q2`.

Evaluating `G_j` in these three embeddings gives the following residues
modulo32. Each displayed integer means that scalar coefficient, with
the four other coefficients zero modulo32:

| embedding `k` / character `j` | 1 | 3 | 5 |
| --- | ---: | ---: | ---: |
| 1 | 16 | 8 | 24 |
| 3 | 8 | 24 | 26 |
| 5 | 24 | 26 | 4 |

Every entry is nonzero modulo32, so these determine the exact valuations.
Consequently the valuation matrix of the three cyclotomic units-in-modulus
`z_j=G_j^2/32`, `j=1,3,5`, is

```math
V=\begin{pmatrix}3&1&1\\1&1&-3\\1&-3&-1\end{pmatrix},
\qquad \det V=-36.
```

Here 'unit-in-modulus' means complex modulus one, NOT a 2-adic unit.
If `z_1^a z_3^b z_5^c=1` for integers `a,b,c`, applying these three
valuations yields `V(a,b,c)^T=0`; therefore all three integers vanish.
The same argument excludes a nontrivial root of unity on the right.
This proves multiplicative independence without invoking Stickelberger's
formula or making an unverified orientation choice.

It follows that `(z_1^r,z_3^r,z_5^r)`, `r>=1`, is dense in the entire
three-torus. One direct proof averages each nonconstant torus character:
the average is a geometric sum with ratio different from one, hence
tends to zero. Approximation of continuous functions by trigonometric
polynomials gives equidistribution, in particular density and infinitely
many visits to every open neighborhood.

## 3. Why the phases are supplied by actual even Walsh orders

Lift each multiplicative character by the norm from `F_(32^s)` to `F32`,
and lift the additive character by the trace. Their Gauss sums satisfy

```math
G_j^{(s)}=(-1)^{s-1}G_j^s.                              (1)
```

For completeness, (1) has a short elementary generating-function proof.
For a monic polynomial `P(T)=T^d+a_1T^(d-1)+...+a_d` over the base field,
put `w(P)=chi((-1)^d a_d) psi(-a_1)`, with the weight zero when `a_d=0`.
Weights multiply under polynomial multiplication. Their sum over monic
polynomials of degree one is `G(chi,psi)`, and for every degree at least
two the independent sum over `a_1` vanishes. Thus

```math
\sum_{P\ \mathrm{monic}}w(P)T^{\deg P}=1+G(\chi,\psi)T.
```

The Euler product over monic irreducibles and its formal logarithm
identify the coefficient of `T^s` as `G^(s)/s`: every nonzero element of
the extension is a root of an irreducible polynomial whose degree divides
`s`, and its lifted norm and trace give the corresponding powered weight.
Comparing with the logarithm of `1+GT` proves (1).
Indeed, an irreducible of degree `d | s` contributes its `d` roots to
the extension sum, each with weight `w(P)^(s/d)`. Thus that sum is
`sum_(d|s) d sum_(P irreducible, deg P=d) w(P)^(s/d)`, exactly `s`
times the coefficient of the formal Euler-product logarithm. The
constant monic polynomial has weight one. These are identities of
formal power series; no convergence or logarithm branch is assumed.

Take `s=2r`. Then `q=32^(2r)=2^(10r)=4^(5r)` is an allowed even Walsh
order, and

```math
G_j^{(2r)}/\sqrt q=-z_j^r.
```

The harmless common minus sign and the density from Section2 therefore
realize independently prescribed limits `exp(i theta_1)`,
`exp(i theta_3)`, `exp(i theta_5)` on the three positive orbits, with
conjugate limits on the opposite orbits.

## 4. The resulting lower tests for the prescribed norm R

The convention here is
`R(B)=(1/2) sup_h beta(H4^(tensor h) tensor B)/(4^h)^(3/2)`,
where `beta(M)=max_(x,y Boolean)|x^T M y|`. There is no additional
normalization by the fixed seed order in this definition.

Given any three real angles, set `lambda_0=0`, set
`lambda_j=exp(i theta_l)` on the orbit of `l in {1,3,5}`, and use its
complex conjugate on the opposite orbit. Define the real symmetric matrix

```math
K_\theta(a,b)=\frac1{31}\sum_{j=1}^{30}
                   \lambda_j e^{2\pi i j(a+b)/31}.
```

Character orthogonality gives exactly
`K_theta 1=0` and `K_theta^2=I-J/31`. No constant component is restored.
The opposite Fourier-sign convention only replaces all three angles by
their negatives and gives the same available family.

Let `B` be any fixed symmetric real seed and `F` any Boolean `31-by-d`
label array. On each nonzero element of `F_q`, use the corresponding
multiplicative-coset row of `F`, and put zero at additive zero. Expanding
these functions in multiplicative characters and applying (1) gives,
along the chosen subsequence, the limit

```math
\frac1q\|U_q f_q B\|_1\longrightarrow
\frac1{31}\|K_\theta F B\|_1.
```

Here is the exact sign convention. If `chi_q(y)=omega^b`, use row `F_b`
and put `Fhat_j=(1/31) sum_b F_b omega^(-jb)`. For nonzero `x` with
`chi_q(x)=omega^a`, additive Fourier transformation gives

```math
(U_q f_q)(x)=-q^{-1/2}\widehat F_0+
 \sum_{j=1}^{30}\widehat F_j\,omega^{-ja}\,G_j^{(2r)}/\sqrt q.
```

Choose the dense phase subsequence so that `G_(-j)^(2r)/sqrt(q)` tends
to `lambda_j`. Replacing `j` by `-j` in this formula gives precisely
the plus-sign kernel `K_theta(a,b)` above. Each of the 31 nonzero cosets
has `(q-1)/31` elements, establishing the normalized `L1` limit.

The nonzero-point error from the trivial character is `O(q^(-1/2))`.
The additive-zero Fourier coefficient is `O(sqrt q)`, so that single
coordinate also contributes only `O(q^(-1/2))` to normalized `L1`.
Replacing the zero input by any Boolean row changes every transformed
coordinate by `O(q^(-1/2))`, so the same limit holds for genuine Boolean
input arrays. The constants in these statements may depend on fixed
`B` and its order, but not on `q`. Standard even Walsh matrices
are independently sign/permutation equivalent to the prescribed regular
order-four powers, as in the previously audited bilinear stabilization.
Thus the actual asymptotic lower family is

```math
R(B)\ge\frac1{62}\max_{F\in\{\pm1\}^{31\times d}}
                      \|K_\theta F B\|_1.                (2)
```

For fixed Boolean `F,G`, its bilinear objective is a sum of three terms
`alpha_l cos(theta_l)+beta_l sin(theta_l)`. Hence those three angles can
be eliminated exactly as `sum_l sqrt(alpha_l^2+beta_l^2)`. This is not
valid for the archived one-angle family, where the angles are tied.

## 5. Bounded weighted-seed test and scope

For `B=[[-1,2],[2,4]]`, the search with 256 restarts, 300 iterations per
restart, seed260906, found only the zero-phase baseline

```math
\frac72\left(1-\frac1{31^2}\right)=3.496357960457856\ldots,
```

not a value above the already known `R(B)>=7/2`. The floating search is
only a lower-witness search; it proves no upper bound or nonexistence.
The full exact certificate, search parameters, and final Boolean witness
are preserved in
`computations/results/transfer_adversary_gauss31_phase_2026_09_06.json`.
The new result is the exact three-dimensional asymptotic phase family (2)
and its finite arithmetic independence certificate. Larger primes,
tensor products, and a global phase/label optimum have not been explored
or asserted here. No further constant improvement is banked from this test.
