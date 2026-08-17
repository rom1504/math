# Independent audit: Boolean product-algebra trust closure

**Verdict: PASS.**  PC.1--PC.3 are rigorous with the stated scope.  The
endpoint selector proof pays the quadratic and linear channels jointly at
one Boolean top eigenvector; the affine-coset equivalence and both
composition laws are exact.  The dense tensor family really has
`p=(log_2 n)/2+1`, `sqrt(n)` histogram bins, and bounded total port mass.

## 1. Selector and witness check

For an even number of inputs, the zero-sum layer has no fixed point under
the antipode.  Hence an odd tie selector exists.  The concrete rule
`tau(a)=a_1` on ties is antipodally odd and agrees with the sign of the sum
off ties.  Its Fourier support contains only odd sets.

For every endpoint word `epsilon`, PC.10 fixes one coordinatewise selector
witness `x_epsilon`; it does not choose a different old spin after looking
only at the aggregate support.  If every active Fourier product is a
positive top eigenvector, Fourier synthesis gives

```math
Hx_\epsilon=rx_\epsilon,
\qquad
z_\epsilon^Tx_\epsilon=\|z_\epsilon\|_1.           \tag{PCA.1}
```

The latter identity remains true on the even-arity tie layer because both
sides contribute zero there.  Thus the same switch attains the child and
field upper bounds.  The separate assumption that all singleton ports are
top eigenvectors is needed and correctly retained for the spherical
formula, even when a chosen selector has a vanishing singleton coefficient.

## 2. Coset and composition check

With `c_i=w_1 odot w_i`, multiplication by
`w_1 product_(i in T)c_i` gives the odd port product on `T` or on
`{1} union T` according to the parity of `|T|`.  This is a bijection of
sets, modulo any relations among the generators, so

```math
\{w_S:|S|\text{ odd}\}=w_1\langle c_2,\ldots,c_p\rangle. \tag{PCA.2}
```

Block concatenation preserves every product identity.  For tensor ports,
the `S`-product is `w_S tensor v_S`, proving spectral closure.  A row pair
has projective type equal to the group product of its two marginal types;
therefore the exact state law is the unnormalized convolution

```math
\mu_{W\boxtimes V}=\mu_W*_G\mu_V.                  \tag{PCA.3}
```

This makes PC.2 a state algebra, not merely preservation of a certificate.
The verifier checks (PCA.3) exactly and checks common-pole amplification as
convolution with a point mass.

## 3. Independent TC/EG seed cross-check

Using the independent lexicographic construction in TC.3, the order-16
regular Hadamard has the product-closed triple `(a,b,c)` and fourth pole
`a odot b odot c`.  The relative words `a odot b` and `a odot c` generate
four distinct elements, and their affine coset is exactly those four top
poles.  This independently corroborates the seed used by PC.3.

For comparison, both four-port equal-`(G,R)` tuples of EG.2/TC.3 fail the
active `p=4` selector condition: at least one required triple product is not
a positive top eigenvector.  Hence PC.1 does not assume away the collision
silently; it identifies precisely the missing higher-product
synchronization.  The independent TC tensor calculation also reproduces
EG.2's support values `2n` and `7n/4` and leading trust separation `rn/8`.

## 4. Growing-family arithmetic

At tensor depth `j`, the two base relative generators in each factor are
independent, giving a group of order `4^j`.  The base pole plus `2j`
generator-translates has

```math
p_j=2j+1,
\qquad
2^{p_j-1}=4^j=\sqrt{16^j}.                         \tag{PCA.4}
```

All odd products form the affine tensor coset and are Boolean positive top
poles.  Storing `sqrt(n)` public histogram counts needs at most
`sqrt(n) ceil(log_2(n+1))` bits.  With `m=floor(r/p)`, one has
`pm<=r=sqrt(n)`, so arbitrary auxiliary exact-sign completion costs
`O(n)=o(n^(3/2))`.  These are strict sub-landscape and bounded-mass
statements, but only for the structured affine-pole interface; the draft
does not claim compression of arbitrary dense bridges.

## 5. Reproduction

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_port_product_algebra_closure.py
```

The verifier exhausts the abstract antipodal benchmark, checks the selector
Fourier expansion through six ports, and checks the dense `j=1,2` tensor
identities without numerical tolerances.
