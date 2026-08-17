# Independent audit: triple response linkage

**Verdict:** PASS, with one editorial scope caution.

**Audited source:**
`extremal_information/drafts/exact_ground_endpoint_linkage.md`, SHA-256
`6be4d518be53c071db5b101127267d6cd1fede8736a0ca5dabbd035351bc3ea4`.

This audit checks TL.1--TL.4 independently, including the conversion from
the GL.7 code shell to positive energy responses, both uses of MI.1, and the
order-five exact computation.  No canonical file is edited here.

## 1. TL.1: PASS

Put

```math
b_e=a_e(z_0)_e,\qquad w_i=z_i z_0\quad(i=1,2).
```

On the cell `w_1=w_2=+1`, direct Fourier expansion gives

```math
\sum_e b_e{1+w_{1,e}\over2}{1+w_{2,e}\over2}
={H_0+H_1+H_2+\langle a,z_0z_1z_2\rangle\over4}.
```

The product of three augmented cuts is an augmented cut, so the last
response is at least `-M`.  The left side is the number of all-correct
edges minus the number of all-wrong edges in that cell.  The all-correct
count therefore dominates it.  Hence

```math
|R(z_0,z_1,z_2)|
\ge {H_0+H_1+H_2-M\over4}
={M\over2}-{d_0+d_1+d_2\over4}.
```

If `d_i<=2s_i`, this is exactly

```math
|R(z_0,z_1,z_2)|\ge{M-s_0-s_1-s_2\over2}.
```

The pair identity TL.6 follows from the same expansion with one relative
word:

```math
|R(z_0,z_1)|\ge{H_0+H_1\over2}
=M-{d_0+d_1\over2}.
```

No opposite-lift, exact-ground, or pairwise-distance hypothesis is hidden
in either identity.

As a separate finite check, exhaustive enumeration of all signings and all
augmented-cut triples at order four tested 37,440 applicable oriented
triples and found no violation of TL.1.

## 2. GL.7-to-energy conversion: PASS

For the augmented cut code,

```math
N={n\choose2},\qquad k_n=\Theta(n),qquad
\rho={N-M_n\over2},\qquad t_n={M_n\over2}.
```

With

```math
u_n=\lceil\sqrt{k_nt_n}\rceil,
```

the known `M_n=Theta(n^(3/2))` bounds give

```math
k_n=o(u_n),\qquad u_n=o(t_n),qquad
u_n=\Theta(n^{5/4}).
```

Thus GL.7 applies uniformly to every deepest coset and supplies two members
of `mathcal L_(u_n)` at projective distance

```math
(2-o(1))t_n=(1-o(1))M_n.
```

A codeword at distance at most `rho+u_n` has response at least

```math
N-2(\rho+u_n)=M_n-2u_n>0
```

eventually.  Hence these are indeed positive words of energy deficit at
most `2u_n`; no orientation or factor-two loss is missing.

Projective distance is the quotient metric modulo complementation.  For any
positive exact ground `g`, its triangle inequality implies that at least one
endpoint `u` of the GL.7 pair satisfies

```math
d_P(g,u)\ge(1/2-o(1))M_n.
```

This step does not require the GL.7 pair itself to contain an exact ground.

## 3. First MI.1 application: PASS

The deficits of `g,u` are at most `0,2u_n`.  TL.6 supplies

```math
|R(g,u)|\ge M_n-u_n=(1-o(1))M_n.
```

Take

```math
\theta_n={1\over2}-{1\over\log n},
\qquad
r_n=\left\lceil{(n\log2+1)\log^2n\over2}\right\rceil.
```

Then

```math
2(1/2-\theta_n)^2r_n\ge n\log2+1,
```

so the MI.1 sampling expression is at most `e^(-1)<1`.  Moreover,

```math
r_n\le |R(g,u)|,qquad 2r_n<M_n,qquad
u_n+r_n<M_n
```

for all sufficiently large `n`.  Applying MI.1 with common shell parameter
`s=u_n` produces a positive word `v` of deficit at most `2r_n` and

```math
d_P(v,g),d_P(v,u)
\ge\min\{\theta_n(M_n-u_n),M_n-u_n-r_n\}
=(1/2-o(1))M_n.
```

The complementary-distance term is asymptotically `M_n`, so it does not
set the displayed constant.

## 4. Second MI.1 application: PASS

TL.1 applied to `g,u,v` gives

```math
|R(g,u,v)|
\ge {M_n-u_n-r_n\over2}
=(1/2-o(1))M_n.
```

For the second MI.1 use, take

```math
s_*=\max\{u_n,r_n\}.
```

Then, eventually,

```math
r_n\le|R(g,u,v)|,qquad 2r_n<M_n,qquad
s_*+r_n<M_n,
```

and the same sampling inequality applies.  The resulting positive word
`w` has deficit at most `2r_n` and satisfies, for each `z` in `{g,u,v}`,

```math
d_P(w,z)
\ge\min\{\theta_n|R(g,u,v)|,M_n-s_*-r_n\}
\ge(1/4-o(1))M_n.
```

Together with the earlier pairwise bounds, this proves TL.8.  The common
energy-deficit width

```math
D_n=2\max\{u_n,r_n\}=o(M_n)
```

has the correct factor two.  Each MI witness comes from a fresh perturbation
of the same exact minimizer, so deficits do not accumulate across the two
applications.

## 5. TL.3 four-anchor ceiling: PASS

After gauging by `z_0`, the indicator that all four words agree is

```math
2^{-3}(1+w_1)(1+w_2)(1+w_3).
```

Multiplication by `az_0` produces exactly four singleton/anchor responses
and four triple-product responses.  Thus TL.18 has all eight terms with
coefficient `1/8` and no missing orientation sign.  Every triple product is
again an augmented cut.  If the anchor responses are `M-d_i`, bounding the
four remaining responses below by `-M` gives precisely

```math
\sum_{e\in\mathcal C_4}a_e(z_0)_e
\ge-{d_0+d_1+d_2+d_3\over8}.
```

In particular the cap-only lower bound is exactly zero, not positive, when
all four anchors are grounds.  This proves the stated ceiling for this
Fourier/cell argument; it does not assert that every four-ground common
reservoir is empty.

The holonomy rewrite is also exact.  With

```math
h=z_0z_1z_2z_3,
```

multiplication by `h` pairs each anchor with the complementary triple
product, so TL.18 becomes

```math
{1\over8}\sum_{i=0}^3
\big(\langle a,z_i\rangle+\langle a,hz_i\rangle\big).
```

If `h=\mathbf1`, then `hz_i=z_i`; the signed common-cell mass is therefore

```math
{1\over4}\sum_{i=0}^3(M-d_i)
=M-{d_0+d_1+d_2+d_3\over4}.
```

The common-correct reservoir is the positive part of that common cell and
hence has at least this size.  The positive-trivial-holonomy corollary is
therefore valid with the displayed constant.

The final jointly paid corollary also has the stated sharp constant.  If
`d_i=o(M)` and

```math
\sum_{i=0}^3\langle a,hz_i\rangle\ge-(4-\eta)M,
```

then TL.18 gives

```math
\sum_{e\in\mathcal C_4}a_e(z_0)_e
\ge {4M-\sum_i d_i-(4-\eta)M\over8}
={\eta M\over8}-o(M).
```

The reservoir, being the positive part of this cell, satisfies the same
lower bound.  Conversely, if this certificate fails to give an
`Omega(M)` reservoir for every fixed `eta>0`, the translated-response sum
must be `-4M+o(M)`.  Since each of its four summands is at least `-M`, all
four are then individually `-M+o(M)`.  Thus “simultaneous near-saturation”
is exactly the remaining failure mode **of this certificate**; it is not
asserted to classify all possible fourwise-collapse mechanisms.

## 6. TL.4 and the order-five computation: PASS

If `z_0z_1z_2z_3=-\mathbf 1`, an edge correct for all four would obey both

```math
\prod_{i=0}^3 a_e(z_i)_e=+1
```

and

```math
\prod_{i=0}^3 a_e(z_i)_e
=a_e^4\prod_{i=0}^3(z_i)_e=-1,
```

a contradiction.  Thus the fourwise reservoir is empty.

Independent exact enumeration of the displayed order-five matrix gives the
sixteen energies

```text
4, 4, 4, 0, 0, 4, -4, -4, 0, -4, 4, -4, 0, 0, 0, -4.
```

Exhausting all `2^10` order-five signings gives the cap histogram

```text
Q=4: 192,   Q=6: 480,   Q=8: 320,   Q=10: 32,
```

so `M_5=4`.  This also follows from the draft's shorter argument:
orthogonality of the ten degree-two characters gives
`E_x H_A(x)^2=10`, and all energies are even.

For masks

```text
0000, 0001, 1010, 1011
```

(decimal `0,1,10,11`), the raw energies are `4,4,4,-4`.  Positive
orientation gives product `-\mathbf 1`; all six projective distances equal
four.  The four triple-reservoir sizes are exactly `2,3,2,2`, while the
fourwise reservoir is empty.

## 7. Scope and one editorial caution

| Claim | Audit |
|---|---|
| TL.1 universal triple reservoir | PASS |
| TL.6 pair reservoir | PASS |
| GL.7 shell/energy normalization | PASS |
| ground-triangle selection | PASS |
| first MI.1 constants and hypotheses | PASS |
| second MI.1 constants and hypotheses | PASS |
| TL.2 common shell width | PASS |
| TL.3 exact four-anchor formula/ceiling | PASS |
| TL.3 jointly paid `eta M/8-o(M)` corollary | PASS |
| TL.4 pointwise holonomy obstruction | PASS |
| order-five exact minimizer | PASS |
| Example 159 reclassification | PASS |
| no implication to fixed-scale `L_projective` | PASS |

TL.2 guarantees only four words and only separation `Omega(M_n)`.  Since

```math
{M_n\over\binom n2}=\Theta(n^{-1/2}),
```

it supplies neither a fixed-edge-scale pair nor the growing fixed-scale
packing required by `L_projective`.  The order-five example is finite and
already has a fixed-normalized four-point packing, so it is not a scalable
counterexample to that lemma.

Example 159 is correctly removed as the depth-three obstruction: it did not
certify its displayed responses to be within `o(M)` of the actual cap, while
TL.1 now handles every genuine positive thin-shell triple.

The only scope caution concerns the proposed next-step language.  TL.4
exhibits **one exact mechanism** for fourwise collapse, but it does not prove
that every small fourwise reservoir contains a negative affine
parallelogram.  Thus “orientation holonomy is the known next obstruction”
is justified; a claim that it classifies all depth-four collapse would not
be.  The selection-or-holonomy statement at the end of the source should be
read as a candidate next lemma, not as an established exhaustive dichotomy.
