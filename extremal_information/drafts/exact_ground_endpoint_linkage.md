# Triple response linkage: four universal energy-scale witnesses and the holonomy ceiling

**Status:** proved draft; no canonical edits.

This note closes one exact gap left by GL.7, MB.2, and MI.2.  Three positive
near-ground augmented cuts always retain a large common-correct reservoir.
The proof is a one-line odd-product/cell identity, but it was absent from the
previous iteration: a fourth Fourier coefficient is itself a legal signed
cut response and is bounded below by `-Q(A)`.

Combined with low-dimensional deep-hole diffusion and exact edge-flip
minimality, the identity proves that **every exact minimizer** has four
mutually `Omega(M_n)`-separated words in an `o(M_n)` positive shell.  This
strictly improves the previous conditional three-word opposite-lift result
and removes Example 159 as a depth-three obstruction.

The gain remains at the energy scale, not at a fixed fraction of
`binom(n,2)`.  An explicit order-five exact minimizer shows a sharp
obstruction at the next step: four exact positive grounds can have negative orientation
holonomy, hence empty fourwise common-correct reservoir.  The next theorem
must select against or exploit that holonomy; an unqualified four-anchor
Helly statement is false.

## 1. Correlation cells and the exact triple identity

Let `A` be an order-`n` signing, let

```math
E={n\choose2},\qquad M=Q(A),
```

and let `z_0,z_1,z_2` be positive augmented-cut words.  Write

```math
H_i=\langle a,z_i\rangle=M-d_i,
\qquad d_i\ge0,                                   \tag{TL.1}
```

For any finite list of augmented-cut words, define its common-correct
reservoir by

```math
R(z_0,\ldots,z_\ell)
=\{e:a_e(z_i)_e=+1\text{ for }0\le i\le\ell\}.    \tag{TL.2}
```

### Theorem TL.1 (universal triple common-correct mass)

For every signing and every three positive augmented-cut words,

```math
\boxed{
|R(z_0,z_1,z_2)|
\ge {H_0+H_1+H_2-M\over4}
= {M\over2}-{d_0+d_1+d_2\over4}.}                 \tag{TL.3}
```

In particular, if `d_i<=2s_i`, then

```math
|R(z_0,z_1,z_2)|
\ge {M-s_0-s_1-s_2\over2}.                        \tag{TL.4}
```

No pairwise-distance or opposite-lift hypothesis is required.

#### Proof

Gauge by `z_0`: put

```math
b_e=a_e(z_0)_e,
\qquad w_1=z_1z_0,
\qquad w_2=z_2z_0.
```

The signed `b`-mass of the cell `w_1=w_2=+1` is

```math
\begin{aligned}
\sum_e b_e{1+w_{1,e}\over2}{1+w_{2,e}\over2}
={1\over4}\big(&\langle a,z_0\rangle
                 +\langle a,z_1\rangle
                 +\langle a,z_2\rangle\\
                &+\langle a,z_0z_1z_2\rangle\big).
                                                               \tag{TL.5}
\end{aligned}
```

The odd product `z_0z_1z_2` is again an augmented-cut word, so its response
is at least `-M`.  Thus (TL.5) is at least the right side of (TL.3).

On this cell, a positive `b_e` is exactly an edge correct for all three
original words; a negative `b_e` is an edge wrong for all three.  The number
of positive edges is at least positive count minus negative count, i.e. at
least the signed mass in (TL.5).  This proves (TL.3), and (TL.4) follows
from (TL.1). `square`

### Pair version

The same calculation with one relative word gives

```math
|R(z_0,z_1)|
\ge {H_0+H_1\over2}
=M-{d_0+d_1\over2}.                               \tag{TL.6}
```

This recovers the `Omega(M)` pair reservoir without invoking the PP.4
bipartite shore.  PP.4 remains useful for the anatomical description of an
opposite lift, but not for the reservoir size used below.

### Cut-endpoint form

For comparison with the previous mesoscopic-shore language, if `z_0` is an
exact ground, switch to `b=a z_0` and put

```math
C(S)=\sum_{e\in\delta(S)}b_e.
```

Then `0<=C(S)<=M`.  A same-lift positive word `c(T)` has `C(T)` near zero;
an opposite-lift positive word `-c(T)` has `C(T)` near `M`.  The two cases
of TL.1 reduce respectively to

```math
\sum_{\delta(S)\setminus\delta(T)}b
={C(S)-C(T)+C(S\triangle T)\over2},
```

and

```math
\sum_{\delta(S)\cap\delta(T)}b
={C(S)+C(T)-C(S\triangle T)\over2}.                \tag{TL.7}
```

The gauge-free identity (TL.5) contains both cases and also allows the base
word to have `o(M)` deficit.

## 2. Every exact minimizer has four thin-shell energy-scale directions

The next theorem combines three already audited inputs with TL.1:

1. the augmented cut code has dimension `Theta(n)`, length `E=Theta(n^2)`,
   and radius deficit `t=M_n/2=Theta(n^(3/2))`;
2. GL.7 forces diameter `(1-o(1))M_n` in a code shell of width
   `u_n=o(M_n)`;
3. MI.1 turns an `Omega(M_n)` common-correct reservoir into one new jointly
   separated word using `O(n log^2 n)=o(M_n)` exact edge flips.

### Theorem TL.2 (unconditional four-witness theorem)

There is a deterministic sequence `D_n=o(M_n)` such that the positive
`D_n`-deficit shell

```math
\{z:\langle a,z\rangle\ge M_n-D_n\}
```

of **every** exact order-`n` minimizer contains four
words `g,u,v,w` satisfying

```math
\boxed{
d_{\rm P}(z,z')\ge(1/4-o(1))M_n
\quad\text{for every distinct }z,z'\in\{g,u,v,w\}.}             \tag{TL.8}
```

One may take the code-shell excess

```math
u_n=\left\lceil\sqrt{k_nt_n}\right\rceil
=Theta(n^{5/4}),                                  \tag{TL.9}
```

where `k_n=dim C_n=Theta(n)` and `t_n=M_n/2`, together with localized-flip
sample size

```math
r_n=\left\lceil{(n\log2+1)\log^2n\over2}\right\rceil
=O(n\log^2n),                                     \tag{TL.10}
```

and `D_n=2max{u_n,r_n}`.

#### Proof

Fix an exact minimizer `A`, and let `g` be any positive ground word.  The
code shell `mathcal L_(u_n)` is the positive response shell of deficit
`2u_n`.  By GL.7 it contains two words whose projective distance is at least

```math
(2-o(1))t_n=(1-o(1))M_n.                           \tag{TL.11}
```

Projective distance is a metric, so one member, call it `u`, satisfies

```math
d_{\rm P}(g,u)\ge(1/2-o(1))M_n.                   \tag{TL.12}
```

Both `g` and `u` lie in the deficit-`2u_n` shell.  By (TL.6), their
common-correct reservoir has size at least

```math
M_n-u_n=(1-o(1))M_n.                               \tag{TL.13}
```

Take

```math
\theta_n={1\over2}-{1\over\log n}.
```

Equation (TL.10) gives

```math
2^n\exp\{-2(1/2-\theta_n)^2r_n\}<1.               \tag{TL.14}
```

Apply the two-anchor case of MI.1 in the reservoir (TL.13).  Exact
minimality of `A` supplies a positive word `v` of deficit at most `2r_n`
such that

```math
d_{\rm P}(v,g),\ d_{\rm P}(v,u)
\ge(1/2-o(1))M_n.                                 \tag{TL.15}
```

The complementary-distance term in MI.1 is
`M_n-o(M_n)` and therefore is not active here.

More explicitly, this first application uses shell parameter `s=u_n`.
For all sufficiently large `n`, (TL.13) and `u_n,r_n=o(M_n)` give

```math
r_n\le M_n-u_n,\qquad 2r_n<M_n,\qquad u_n+r_n<M_n,
```

while (TL.14) is MI.1's sampling condition.

Now apply TL.1 to `g,u,v`.  Their deficits are at most
`0,2u_n,2r_n`, respectively, so

```math
|R(g,u,v)|
\ge {M_n-u_n-r_n\over2}
=(1/2-o(1))M_n.                                   \tag{TL.16}
```

Run MI.1 once more, with a fresh `r_n`-set inside this triple reservoir.
This time use the common shell parameter
`s_*=max{u_n,r_n}`.  Equation (TL.16) gives `r_n<=|R(g,u,v)|` eventually,
and `2r_n<M_n`, `s_*+r_n<M_n`, and (TL.14) verify the remaining hypotheses.
The theorem returns a positive word `w` of deficit at most `2r_n` and,
simultaneously for `z=g,u,v`,

```math
d_{\rm P}(w,z)
\ge\min\{\theta_n|R(g,u,v)|, M_n-o(M_n)\}
\ge(1/4-o(1))M_n.                                 \tag{TL.17}
```

Equations (TL.12), (TL.15), and (TL.17) prove (TL.8).  Since
`u_n,r_n=o(M_n)`, the declared common shell width `D_n` is `o(M_n)`.
`square`

### Quantifier and normalization check

The theorem is uniform over every exact minimizer.  GL.7 is applied to the
augmented cut code of the chosen signing's deepest coset; it supplies the
pair before `g` is used.  The triangle inequality is the only step selecting
which endpoint of that pair becomes `u`.  Code distance excess `u_n`
corresponds to energy deficit `2u_n`, which accounts for the factor two in
`D_n`.  The lower bound `M_n>=c_0n^(3/2)` makes both (TL.9) and (TL.10)
`o(M_n)` and verifies all positivity conditions in MI.1.

## 3. Why the same argument stops at four

For four anchors, the analogue of (TL.5) has eight Fourier terms.  The four
terms not fixed by the individual near-ground energies may all have the
unfavorable sign, so cap boundedness gives no positive lower bound.  This is
not merely a weakness of the estimate.

### Proposition TL.3 (exact four-anchor cell formula and cap ceiling)

For four positive augmented-cut words `z_0,z_1,z_2,z_3`, write
`\langle a,z_i\rangle=M-d_i` and let

```math
\mathcal C_4=\{e:(z_0)_e=(z_1)_e=(z_2)_e=(z_3)_e\}.
```

Then their signed common-cell mass is exactly

```math
\begin{aligned}
\sum_{e\in\mathcal C_4}a_e(z_0)_e={1\over8}\big(&
 \langle a,z_0\rangle+\langle a,z_1\rangle
 +\langle a,z_2\rangle+\langle a,z_3\rangle\\
 &+\langle a,z_0z_1z_2\rangle
  +\langle a,z_0z_1z_3\rangle\\
 &+\langle a,z_0z_2z_3\rangle
  +\langle a,z_1z_2z_3\rangle\big).              \tag{TL.18}
\end{aligned}
```

If `h=z_0z_1z_2z_3` is the four-anchor holonomy, the same identity can be
written more symmetrically as

```math
\sum_{e\in\mathcal C_4}a_e(z_0)_e
={1\over8}\sum_{i=0}^3
 \big(\langle a,z_i\rangle+\langle a,hz_i\rangle\big).
```

The reservoir `R(z_0,z_1,z_2,z_3)` is precisely the positive part of
`\mathcal C_4`, so its size is at least this signed mass whenever the latter
is positive.  In particular, positive trivial holonomy `h=\mathbf1` gives

```math
|R(z_0,z_1,z_2,z_3)|
\ge M-{d_0+d_1+d_2+d_3\over4}.
```

Consequently, if the four deficits are `d_0,d_1,d_2,d_3`, cap boundedness
alone gives only

```math
\sum_{e\in\mathcal C_4}a_e(z_0)_e
\ge-{d_0+d_1+d_2+d_3\over8}.                     \tag{TL.19}
```

#### Proof

Gauge by `z_0`, put `w_i=z_iz_0`, and expand
`2^{-3}\prod_{i=1}^3(1+w_i)` exactly as in (TL.5).  The four terms indexed
by the empty set and the three singletons are the four displayed anchor
responses.  The other four are the displayed triple-product responses.
Every triple product is an augmented-cut word and hence has response at
least `-M`.  Substituting `\langle a,z_i\rangle=M-d_i` proves (TL.19).
`square`

At zero deficit the right side of (TL.19) is exactly zero.  Thus even the
complete four-anchor Fourier expansion cannot certify one common-correct
edge from cap boundedness.  More quantitatively, for `d_i=o(M)` the identity
does certify an `Omega(M)` reservoir whenever, for some fixed `eta>0`,

```math
\sum_{i=0}^3\langle a,hz_i\rangle\ge-(4-\eta)M.
```

This is a genuinely joint four-channel condition: the four translated
responses are paid only through their sum.  Therefore the precise remaining
failure mode for this certificate is simultaneous near-saturation of that
sum at `-4M`, not an individual scalar-channel estimate.  The following
obstruction realizes the zero ceiling pointwise.

### Proposition TL.4 (negative holonomy kills the common reservoir)

For any four augmented-cut words satisfying

```math
z_0z_1z_2z_3=-\mathbf1,                            \tag{TL.20}
```

and any signing `A`,

```math
R(z_0,z_1,z_2,z_3)=\varnothing.                   \tag{TL.21}
```

#### Proof

If an edge were correct for all four words, then all four numbers
`a_e(z_i)_e` would equal `+1`, so their product would be `+1`.  But that
product is `a_e^4 prod_i(z_i)_e=-1` by (TL.20), a contradiction. `square`

This obstruction occurs inside a genuine exact minimizer.  Consider

```math
A=\begin{pmatrix}
0& 1& 1& 1& 1\\
1& 0&-1& 1&-1\\
1&-1& 0&-1& 1\\
1& 1&-1& 0& 1\\
1&-1& 1& 1& 0
\end{pmatrix}.                                    \tag{TL.22}
```

Over the sixteen projective spins (`x_1=1`), its energies in binary-mask
order are

```text
 4, 4, 4, 0, 0, 4,-4,-4, 0,-4, 4,-4, 0, 0, 0,-4.
```

Thus `Q(A)=4`.  Conversely every order-five signing has

```math
\mathbb E_x H_A(x)^2={5\choose2}=10.
```

Its energies are even integers, so their maximum absolute value is at least
`4`.  Hence `M_5=4` and (TL.22) is an exact minimizer.

Take projective spin masks

```math
(0000)_2,\quad(0001)_2,\quad(1010)_2,\quad(1011)_2
\qquad\text{(decimal }0,1,10,11\text{)}.
```

Their energies are `4,4,4,-4`.  Orient each signed-cut word positively, so
the fourth receives a minus sign.  The masks XOR to zero but the four
orientation signs multiply to `-1`; therefore (TL.20) holds.  All four words
are exact positive grounds, all six projective edge distances equal
`4=M_5`, every triple common-correct reservoir has size at least two (the
bound in TL.1), yet their fourwise reservoir is empty.

The matrix and enumeration are also contained in the exact-minimizer
record
`extremal_information/experiments/nearmin_shell_parallelogram_results.json`;
the displayed table makes the verification independent of that file.

This is **not** a scalable counterexample to `L_projective`: the four words
already form the desired fixed-normalized finite-order packing, and the
example says nothing about large `n`.  Its rigorous role is to falsify
unqualified four-anchor common-correct persistence, even under exact
minimality and zero shell deficit.

## 4. Revision of the archived obstruction and the SML

1. **Example 159 is no longer the depth-three obstruction.**  Its displayed
   words were not certified to lie within `o(M)` of their signing's actual
   cap.  TL.1 proves that a genuine positive thin-shell triple has
   `(1/2-o(1))M` common-correct mass, regardless of pairwise shore geometry.
2. **The MB/MI iteration advances unconditionally once more.**  GL.7, TL.6,
   MI.1, TL.1, and MI.1 yield the four-point theorem TL.2 with no cumulative
   shell loss.
3. **The revised obstruction is orientation holonomy at depth four.**  A
   negative affine parallelogram can recycle the fourth response while
   annihilating the common reservoir identically.  The order-five example
   shows that exact minimization alone does not prohibit this finite pattern.
4. **`L_projective` does not change.**  Since
   `M_n/binom(n,2)=Theta(n^(-1/2))`, the separation in (TL.8) remains
   mesoscopic.  It neither gives a fixed-edge-scale pair nor a growing
   fixed-scale packing.

The smallest honest next lemma is a **selection-or-holonomy dichotomy**:
for the four words produced in TL.2, either select the fourth response so an
`Omega(M_n)` fourwise reservoir survives, or prove that negative orientation
holonomy forces a genuinely new projective direction at scale
`Theta(n^2)`.  Proposition TL.4 shows that the first alternative cannot hold
for every possible active quadruple; witness selection or holonomy
amplification is essential.
