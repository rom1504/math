# Stable four-point recovery for nearly flat Walsh spectra

Date: 2026-09-06. Director derivation; full independent audit and rational
certificate replay passed.
This is an actual supported-spin count, not a repair-to-bent assumption.

Use m=2^d, k/m=p in (0,1), q=1-p, and the epsilon-flatness definition
from `continued_feedback_stable_flat_profile_count_2026_09_06.md`:

```math
  \sum_y (|(H\xi)(y)|-\sqrt{k})^2\le\epsilon^2km,
  \qquad \xi\in\{0,\pm1\}^m,\quad |\operatorname{supp}\xi|=k.
                                                               (1)
```

Let G_m be the selector event that every nonconstant normalized Walsh
coefficient of 1_T is at most m^(-1/4) in absolute value. Then
Pr(G_m)->1 for a uniform k-selector. In fact sampling-without-replacement
bounded-difference concentration, followed by a union bound over m-1
characters, gives a bound of the form 2m exp(-c sqrt(m)).

Let X=(X_1,...,X_4) have independent coordinates with probabilities
P(0)=q, P(+1)=P(-1)=p/2. Define

```math
 Y=(X_1+X_2+X_3+X_4,
    X_1+X_2-X_3-X_4,
    X_1-X_2+X_3-X_4).
```

The conclusion is

```math
 \lim_{\epsilon\downarrow0}\limsup_m\frac1m
 \log E_T[1_{G_m}N_T(\epsilon)]
 \le c_4(p):=(p+3/4)\log2-\tfrac14 H(Y).                (2)
```

Here N_T counts the actual epsilon-flat rows on T. More quantitatively,
for small epsilon the right side can be increased by

```math
 C_p\epsilon+h(e)+e\log13,
 \qquad e=64p\epsilon^2<13/14,                          (3)
```

where a fixed finite C_p suffices on each compact subinterval of (0,1).
For example C_p=10000 log(1/min(q,p/2)) is more than sufficient in the
proof below. These loose error constants are not optimized.

## 1. Flatness supplies the four-point information actually needed

For normalized Fourier coefficients \hat f=m^(-1)Hf, write
U_2(f)^4=sum_y |\hat f(y)|^4. This is a norm and satisfies the exact
four-point Fourier identity on F_2^d. If h=H xi and
h=sqrt(k) z+v with z Boolean and ||v||_2²<=epsilon²km, then

```math
 U_2(\xi)^4\le8p^2(m^{-1}+\epsilon^4).                 (4)
```

Indeed (a+b)^4<=8(a^4+b^4), and sum v_y^4<=||v||_2^4.
For T in G_m, Parseval and the maximum nonconstant coefficient give
U_2(1_T-p)<=m^(-1/8). The centered indicator of each of the three
symbols is a linear combination of 1_T-p and xi:

```math
 1_{\xi=0}-q=-(1_T-p),\qquad
 1_{\xi=\pm1}-p/2=((1_T-p)\pm\xi)/2.
```

Consequently all these centered indicators have U_2 norm
O(epsilon+m^(-1/8)). Expanding a four-symbol pattern at
(x,x+a,x+b,x+a+b), its average over x,a,b differs from the product-P
probability by O(epsilon+m^(-1/8)). To verify the bound without an
imported mixing theorem, expand the product of the four Fourier series:
the surviving terms have a common frequency. Holder's inequality with
four exponents 4 bounds every term containing a centered indicator by
its U_2 norm, since the other bounded indicators have U_2 norm at most 1.
There are at most 15 nonconstant expansion terms. Restricting to linearly
independent a,b changes the estimate by O(1/m).

Thus every one of the 81 ordered ternary patterns on a random affine
two-plane has its iid-P frequency up to O(epsilon+m^(-1/8)). This uses
only four-point information. It makes no claim about higher cube laws.

## 2. A short conditional recovery description on some plane direction

Fix independent a,b. On each four-point coset of span(a,b), omit the
parity character chi=(1,-1,-1,1). The orthogonal projection Q onto the
other three characters has kernel span(chi), so QX and Y determine each
other. Under the product-P reference law let

```math
 c(x)=-\log P^{\otimes4}(X=x\mid Y=Y(x)).
```

This lies in [0,4 log(1/min(q,p/2))]. Its value is unchanged when the
four positions are translated within their affine two-plane: the product
prior is permutation-invariant and chi changes only by a sign. Hence
the coset cost is independent of the choice of its representative.

By Section 1, averaging the sum of these costs over independent a,b
gives at most

```math
 (m/4)H(X\mid Y)+C_p\epsilon m+o(m).
```

Some direction pair therefore has cost no larger. Encode that pair
using O(log m) bits. Given all the exact coset projections, the product
conditional probabilities on their fibres sum to one. The number of
preimages with cost at most L is at most exp(L). This is a counting
argument, not an assumption that the actual rows follow the reference law.

## 3. Three quarters of the spectrum approximately give the projections

The union of the frequency hyperplanes a-perp and b-perp has size 3m/4.
The inverse Walsh transform after retaining exactly those coefficients
is the block projection Q just described. Encode their signs (3m/4 bits),
replace their magnitudes by sqrt(k), and invert. Orthogonality and (1)
bound the squared error of the resulting approximate projection by
epsilon² k.

Every true projected coordinate is a multiple of 1/4 in [-3/2,3/2].
Round to the nearest of these 13 values. A wrong rounded coordinate
costs at least 1/64 in squared error; at most 64epsilon² k coordinates
are wrong. Encode their locations and exact values. This uses at most

```math
 \sum_{j\le em}\binom mj13^j
   =\exp\{m[h(e)+e\log13]+o(m)\}
```

choices. The true projection is now known, and Section 2 bounds its
possible low-cost preimages. The number of eligible ternary xi is at most

```math
 \exp\{m[\tfrac34\log2+\tfrac14H(X\mid Y)
             +C_p\epsilon+h(e)+e\log13]+o(m)\}.
```

Divide by binomial(m,k). Since H(X)=4[h(p)+p log2], this proves (2)-(3).
Conditioning the selector law on G_m, or intersecting with the previously
proved high-probability selector filters, leaves the rate unchanged.

## 4. Exact numerical comparison and the limitation

At p=15/16 the finite 81-pattern distribution gives

```text
 H(Y) = 3.4556375720844534964125897411693377...
 c_4  = 0.3057764741737943355384317696683385...
```

The pair-sum constant is 0.3168851694355724..., so this is a strict
improvement. Apply the already proved uniform permanent perturbation
bound to obtain

```math
 \limsup_m\frac1m\log E_T\left[1_{G_m}
       \sum_{x\text{ epsilon-flat}}L_t(|H[T,:]^Tx|)\right]
 \le c_4(p)+\tfrac12\log[(1+e^{-4t})/2]
          +C_p\epsilon+h(e)+e\log13+4t\epsilon.          (5)
```

At p15/16 and t near 1.0317185, the limiting expression (5) plus
t(1-sqrt(p)) is approximately -0.00003352074. A reproducible rational
interval calculation at t=33/32 accompanies this file. Therefore for
sufficiently small fixed epsilon this declared dense profile class has
a negative tilted exponent, even though the pair-sum estimate did not.

This does NOT give an upper bound for the full signing. At this tilt
the typical-Gaussian profile already contributes a positive exponent.
At the tilts roughly 4--7 not excluded by the full-sum lower obstructions,
the flat-class bound (5) remains positive. One cannot choose separate
tilts for separately paid row classes inside the existing Finner product
and thereby infer the full first moment. A new combined count or a valid
multi-class graph argument would be required.

The useful statement is the stable recovery/counting theorem itself:
near-flat information, plus a cheaply imposed selector property, supplies
precisely enough four-point data for a better description without any
exact Walsh divisibility or an assumed repair to a bent function.
