# Walsh/bent stability and counting audit for the regular-Hadamard basin

## Scope and conclusion

This note audits primary literature for a theorem that could turn a
linear-sized diffuse unmatched core in the regular-Hadamard example into
quantitative bent-like or plateaued structure.  The normalization and the
logical interface with the greedy problem are made explicit below.

The useful result is a **counting theorem**, not a stability theorem.  The
explicit bad endpoint is an affine-support ("trivial") 2-plateaued Boolean
function.  Potapov's 2024 bent-function bound and the exact affine-support
classification imply that this entire class, and every sufficiently small
fixed relative Hamming neighborhood of it, is exponentially sparse.  I did
not find a theorem proving that a large diffuse core, a Walsh fourth-moment
excess, or a Walsh `L1` deficit places the *input* in such a neighborhood.
Nor does sparsity of terminal states bound the basins of a deterministic
greedy map.  Thus the literature does not currently prove decay of the bad
basin probability.

## 1. Exact Walsh coordinates

Let

```math
V=\mathbb F_2^k\times\mathbb F_2^k,\qquad
d=2k,\qquad m=|V|=2^d=s^2,
```

and put `q(a,b)=a\mathbin\cdot b` and
`f(u)=(-1)^{q(u)}`.  For the polar form

```math
\beta((a,b),(c,e))=a\mathbin\cdot e+c\mathbin\cdot b,
```

write

```math
H_{u,v}=(-1)^{\beta(u,v)},\qquad D=\operatorname{diag}(f),
\qquad K=DHD,\qquad C=K-I.
```

For a spin vector `X\in\{\pm1\}^V`, set `g=DX=fX` and use the
unnormalized Walsh transform

```math
\widehat g=Hg.
```

Then, with no loss or inequality,

```math
KX=D\widehat g,
\qquad X_u(CX)_u=g_u\widehat g(u)-1,
\qquad X^{\mathsf T}CX=g^{\mathsf T}\widehat g-m.       \tag{1}
```

The regularizing input `X=\mathbf1` becomes `g=f`.  The quadratic Gauss-sum
identity is precisely

```math
\widehat f=sf,                                           \tag{2}
```

so this input is a self-dual bent vector in the standard normalized
Hadamard convention `m^{-1/2}Hf=f`.

For the codimension-two subspace `S` and the bad terminal
`r=\mathbf1-2\mathbf1_S` from
[`regular_hadamard_local_stability_obstruction.md`](regular_hadamard_local_stability_obstruction.md),
put `h=Dr=fr`.  The identity `Kr=2s\mathbf1_S` is equivalent to

```math
\widehat h=2s\,f\mathbf1_S.                              \tag{3}
```

Thus `h` is exactly 2-plateaued: its Walsh values are `0` or `\pm2s`, and
its Walsh support is an affine (here linear) subspace of codimension two.
Moreover

```math
h_u\widehat h(u)=
\begin{cases}
-2s,&u\in S,\\
0,&u\notin S,
\end{cases}                                             \tag{4}
```

with the displayed signs after the identification in the obstruction.
The terminal inequality for the objective with signing `-C` is therefore
the pointwise Walsh condition `h_u\widehat h(u)\le0`; the diagonal `-I`
makes every one-flip inequality strict even where (4) vanishes.

The energy bookkeeping also has a direct spectral form.  If
`d_\xi=|\widehat g(\xi)|` for the initial regularizing vector and
`c_\xi=|\widehat z(\xi)|` for an anti-regularizing terminal, Parseval gives
`\sum d_\xi^2=\sum c_\xi^2=m^2`, and

```math
B(g):=ms-\sum_\xi d_\xi
 ={1\over2s}\sum_\xi(d_\xi-s)^2,                        \tag{5}
```

```math
A(z):=ms-\sum_\xi c_\xi
 ={1\over2s}\sum_\xi(c_\xi-s)^2.                        \tag{6}
```

In the full-shore branch the clipped defect is

```math
\Delta=[A(z)-B(g)-2m]_+.                                \tag{7}
```

Consequently the dangerous event is not merely Walsh non-flatness.  It is a
large *increase* in non-flatness along the particular greedy dynamics.

For autocorrelation
`R_z(t)=\sum_u z(u)z(u+t)`, Wiener--Khinchin gives

```math
\sum_\xi|\widehat z(\xi)|^4=m\sum_tR_z(t)^2.             \tag{8}
```

Since

```math
(c_\xi^2-m)^2=(c_\xi-s)^2(c_\xi+s)^2
\ge m(c_\xi-s)^2,
```

(6)--(8) imply the rigorous one-sided estimate

```math
\sum_{t\ne0}R_z(t)^2\ge2sA(z).                          \tag{9}
```

This is not a stability theorem: it supplies no inverse conclusion or
closeness to a plateaued class.

## 2. Exact self-dual and plateaued results

### Self-dual bent vectors

Carlet, Danielsen, Parker and Solé define the normalized Sylvester transform
`U=2^{-d/2}H`.  A Boolean sign vector `F` is bent exactly when `UF` is
Boolean, self-dual when `UF=F`, and anti-self-dual when `UF=-F`.  Their
Rayleigh-quotient theorem says

```math
|F^{\mathsf T}HF|\le2^{3d/2},                            \tag{10}
```

with equality precisely for a self-dual or anti-self-dual bent vector
(with the sign of the quotient selecting which).  This gives the exact
standard vocabulary for (2), but no neighborhood stability or asymptotic
count.  They enumerate, for example, `2`, `20`, and `42,896` self-dual bent
functions in dimensions `2`, `4`, and `6`.

Source: C. Carlet, L. E. Danielsen, M. G. Parker and P. Solé,
["Self-dual bent functions"](https://doi.org/10.1504/IJICOT.2010.032864),
*International Journal of Information and Coding Theory* 1 (2010),
384--399; [author manuscript](https://ii.uib.no/~larsed/sdbent.pdf),
especially Section 3.

Shi, Li, Cheng, Crnković, Krotov and Solé extend this language to an
arbitrary Hadamard matrix `H`: a self-dual Hadamard bent sequence satisfies
`HX=\sqrt v X`.  Their strong automorphism group is the signed-monomial
centralizer

```math
\operatorname{SAut}(H)=\{P:PH=HP\}.                     \tag{11}
```

For the Sylvester matrix they prove that the extended affine transformation
`T_{L,b,e,c}` is strong exactly under the displayed orthogonality and
translation restrictions in their theorem (`L^{\mathsf T}=L^{-1}`, equal
translations, and even translation weight).  This is useful for quotienting
the *Hadamard landscape*.  It does not make a lexicographically tie-broken
greedy map equivariant.

Source: M. Shi et al.,
["Self-dual Hadamard bent sequences"](https://doi.org/10.1007/s11424-023-2276-8),
*Journal of Systems Science and Complexity* 36 (2023), 894--908;
[arXiv:2203.16439](https://arxiv.org/abs/2203.16439).

Kutsenko determines the full group among Hamming isometries preserving
self-duality (for even `d\ge4`).  In his notation it consists of transforms

```math
F(x)\longmapsto F(L(x+c))+\langle c,x\rangle+e,
\qquad LL^{\mathsf T}=I,
```

with the stated even-weight restriction on `c`.  This prevents silently
using a larger Sylvester symmetry group.  It is not, by itself, a theorem
about the full automorphism group of the quadratically switched matrix `K`.

Source: A. Kutsenko,
["The group of automorphisms of the set of self-dual bent functions"](https://doi.org/10.1007/s12095-020-00438-y),
*Cryptography and Communications* 12 (2020), 881--898;
[IACR ePrint 2019/1408](https://eprint.iacr.org/2019/1408), Theorems 6 and 8.

### Affine-support 2-plateaued vectors

Hodžić, Pasalic, Wei and Zhang prove the precise spectral construction
needed here.  An `r`-plateaued `d`-variable Boolean function has Walsh values
in

```math
\{0,\ \pm2^{(d+r)/2}\}.
```

When its support is an affine subspace of dimension `d-r`, the signs on that
support, after any fixed affine coordinatization, must be the sign sequence
of a bent function on `d-r` variables, and this condition is sufficient.
They call these affine-support examples trivial plateaued functions and
identify them with partially bent functions.  Their linear-structure result
also gives an `r`-dimensional space of linear structures in this case.
Equation (3) is exactly their `r=2` case, not merely an analogy.

Source: S. Hodžić, E. Pasalic, Y. Wei and F. Zhang,
["Designing plateaued Boolean functions in spectral domain and their classification"](https://doi.org/10.1109/TIT.2019.2909910),
*IEEE Transactions on Information Theory* 65 (2019), 5865--5879;
[arXiv:1811.04171](https://arxiv.org/abs/1811.04171), especially the affine
support theorem and the following linear-structure characterization.

## 3. Counting consequences, with constants

Let `B(d)` be the number of bent Boolean functions in even dimension `d`.
Potapov proves

```math
\log_2 B(d)\le {11\over32}2^d(1+o(1)).                  \tag{12}
```

He also gives an upper bound for all fixed-order plateaued functions.  For
`r=2`, substitution in his Theorem 1 yields

```math
\log_2 \#\{\text{2-plateaued functions on }V\}
\le(0.5114896483\ldots+o(1))m.                          \tag{13}
```

Indeed the coefficient is

```math
{1+{3\over8}\log_2 6\over8}
+{H_2(1/4)+1/4\over4}=0.5114896483\ldots .              \tag{14}
```

Source: V. N. Potapov,
["Upper bounds on the numbers of binary plateaued and bent functions"](https://doi.org/10.1007/s12095-024-00766-3),
*Cryptography and Communications* 16 (2024), 1347--1356;
[arXiv:2303.16547](https://arxiv.org/abs/2303.16547), Theorems 1 and 2.

For the smaller class that contains the explicit obstruction, the affine
support classification gives a sharper count directly.  There are

```math
4{d\brack2}_2={2(m-1)(m-2)\over3}                       \tag{15}
```

affine codimension-two subspaces of `V`.  For each fixed support, the
spectral sign patterns are in bijection with bent functions on `d-2`
variables.  Therefore

```math
N_{\rm triv,2}(d)
={2(m-1)(m-2)\over3}\,B(d-2),                           \tag{16}
```

and (12) gives

```math
\log_2N_{\rm triv,2}(d)
\le\left({11\over128}+o(1)\right)m.                     \tag{17}
```

Thus exact affine-support 2-plateaued vectors occupy at most
`2^{-(117/128-o(1))m}` of the Boolean cube.

The same argument gives a genuine neighborhood statement.  The union of
relative Hamming balls of radius `\delta m` around this class has size at
most

```math
2^{[11/128+H_2(\delta)+o(1)]m}.                          \tag{18}
```

It is exponentially sparse whenever

```math
H_2(\delta)<117/128,
\quad\text{equivalently}\quad
\delta<0.3291567\ldots .                                \tag{19}
```

Using (13) instead, the neighborhood of *all* exact 2-plateaued functions is
exponentially sparse at least when

```math
H_2(\delta)<1-0.5114896483\ldots,
\quad\text{i.e.}\quad \delta<0.1062533\ldots .          \tag{20}
```

Equations (17)--(20) are useful only after an inverse theorem or basin
preimage estimate has been proved; they do not supply either one.

Agievich's two-row bent-rectangle method gives exact finite upper bounds on
`B(d)` that are useful at small dimensions but does not improve the
asymptotic exponent (12), and it contains no near-flat stability result.

Source: S. Agievich,
["Upper bounding the number of bent functions using 2-row bent rectangles"](https://eprint.iacr.org/2023/497),
IACR ePrint 2023/497.

## 4. Why `L1` or fourth moment alone cannot identify the bad class

For a uniformly random Boolean sign vector `g`, every Walsh coefficient has
the law of a sum `S_m` of `m` independent signs.  Hence

```math
\mathbb E\sum_\xi|\widehat g(\xi)|
=m\,\mathbb E|S_m|
\sim\sqrt{2/\pi}\,m^{3/2},                              \tag{21}
```

so

```math
{\mathbb EB(g)\over ms}\longrightarrow
1-\sqrt{2/\pi}=0.202115\ldots .                         \tag{22}
```

Also, by the fourth moment of a Rademacher sum,

```math
\mathbb E\sum_\xi|\widehat g(\xi)|^4
=m(3m^2-2m)=3m^3-2m^2.                                 \tag{23}
```

By comparison, a bent vector has Walsh `L1=ms` and fourth moment `m^3`,
whereas the explicit 2-plateaued endpoint has Walsh `L1=ms/2`,
`A=ms/2`, and fourth moment `4m^3`.  Thus an order-`m^{3/2}` `L1` deficit
and an order-`m^3` fourth-moment excess are generic-scale events.  The bad
endpoint's fourth moment is only asymptotically `4/3` times the random mean.
No inference of plateaued structure follows from their magnitude alone.

Schmidt proves the sharper generic peak statement

```math
{\max_\xi|\widehat g(\xi)|\over\sqrt{2m\log m}}
\longrightarrow1
```

almost surely.  Generic spectra therefore have a logarithmically larger
peak than bent spectra, but this still does not control the greedy endpoint
or its basin.

Source: K.-U. Schmidt,
["Nonlinearity measures of random Boolean functions"](https://doi.org/10.1007/s12095-015-0164-3),
*Cryptography and Communications* 8 (2016), 637--645;
[arXiv:1308.3112](https://arxiv.org/abs/1308.3112), Theorem 1 specialized to
first-order Reed--Muller nonlinearity.

## 5. Reed--Muller covering radius: exact map, no stability payoff

For the Boolean function whose sign vector is `g`, the distance to the
first-order Reed--Muller code is exactly

```math
\operatorname{dist}(g,\mathrm{RM}(1,d))
={m-\max_\xi|\widehat g(\xi)|\over2}.                   \tag{24}
```

Bent functions are exact deep holes when `d` is even.  Schmidt proves the
covering-radius asymptotic (in his notation)

```math
2^{d/2}-{\rho_d\over2^{d/2-1}}\longrightarrow1.         \tag{25}
```

He also constructs asymptotically optimal functions in odd dimensions,
where bent functions do not exist.  Thus near-optimal covering radius does
not by itself force exact bent algebra, and (24) concerns the Walsh maximum,
not the `L1` increase in (7).

Source: K.-U. Schmidt,
["Asymptotically optimal Boolean functions"](https://doi.org/10.1016/j.jcta.2018.12.005),
*Journal of Combinatorial Theory, Series A* 164 (2019), 50--59;
[arXiv:1711.08215](https://arxiv.org/abs/1711.08215).

Metric-complement results for `RM(1,d)` likewise identify exact bent deep
holes but provide neither a quantitative neighborhood inverse theorem nor a
greedy-basin estimate.

Source: V. Oblaukhov,
["On metric regularity of Reed--Muller codes"](https://arxiv.org/abs/1912.10811),
arXiv:1912.10811.

## 6. Automorphisms and the deterministic-basin gap

For the quadratically switched Walsh matrix `K`, affine isometries of `q`
give an explicit objective-preserving subgroup (translations combined with
the appropriate quadratic correction, and the linear group
`O^+(2k,2)`).  The exact full group can instead be obtained from the
strong-automorphism/digraph machinery of Shi et al.  Either can quotient
terminal states and the unoriented local-search landscape.

There is a decisive algorithmic caveat.  A rule such as "flip a largest
positive gain and break ties by the smallest coordinate" is not invariant
under those permutations.  Orbit quotienting is exact for basin
probabilities only if

- the chosen gain is unique at every step;
- the tie-breaking is replaced by an equivariant randomized rule; or
- the quotient retains the stabilizer of the imposed coordinate order and
  all tie information.

Consequently the automorphism results above do not by themselves determine
the basin probability of the explicit stable core.

There is one general fact about terminal *states*.  Spink proves that a real
quadratic polynomial on the `m`-cube has at most
`\binom m{\lfloor m/2\rfloor}` strict local maxima, and this is sharp.  It
only bounds the fraction of states that are strict terminals by
`O(m^{-1/2})`; one strict terminal may still have an exponentially large
basin.

Source: H. Spink,
["Local maxima of quadratic Boolean functions"](https://doi.org/10.1017/S0963548315000322),
*Combinatorics, Probability and Computing* 25 (2016), 633--640;
[arXiv:1310.1570](https://arxiv.org/abs/1310.1570).

## 7. Precise missing theorem

The counting results would become relevant if one could prove any one of
the following genuinely new statements.

1. **Initial inverse theorem.**  If the prescribed dynamics creates a
   linear diffuse core, or satisfies
   `A(z)-B(g)\ge\eta ms`, then its initial vector `g` lies within
   `\delta(\eta)m` of the affine-support 2-plateaued class, with
   `\delta(\eta)<0.3291567`; alternatively, closeness to the full
   2-plateaued class with `\delta(\eta)<0.1062533` would suffice.

2. **Basin-entropy theorem.**  The union of preimages of bad affine-support
   2-plateaued terminals has size `2^{(1-\varepsilon)m}`.  From (17), it
   would suffice to bound each such basin by
   `2^{(117/128-\varepsilon)m}` (or prove the union bound directly).

3. **Direct dynamic tail theorem.**  Bound the probability of
   `A(z)-B(g)\ge\eta ms` for the actual deterministic tie rule, without
   classifying the endpoint.

No checked source proves one of these statements.  In particular, there is
no located theorem turning the Walsh fourth moment or `L1` deficit alone
into Hamming closeness to a bent or plateaued function.  The current exact
identities therefore support the following research judgment:

> Large diffuse cores should not be called bent-like merely because their
> Walsh fourth moment is large.  Continue the bent/plateaued route only if
> computation reveals a dynamics-specific correlation with affine-support
> plateaued structure or supplies evidence for one of the three inverse or
> basin statements above.  If large diffuse cores occur for generic Walsh
> spectra, the available literature gives no leverage and this route should
> be stopped.

## 8. Recent papers checked but not promoted to a bridge

- Kharaghani et al.,
  ["On Regular Quaternary Hadamard Matrices"](https://doi.org/10.1002/jcd.70006)
  (2026), study vectors whose transformed coordinate magnitudes are
  constant.  In the real Boolean Walsh setting Parseval makes this exactly
  the bent condition; their finite quaternary constructions do not give
  approximate stability or basin counts.
- Kölsch and Polujan,
  ["The Combinatorial Structure and Value Distributions of Plateaued Functions"](https://doi.org/10.1007/s00145-026-09578-5)
  (2026), give exact plateaued-function structure and distributions, not an
  approximate-plateaued inverse theorem for (7).
- Work titled "stability of Walsh spectrum with biased inputs" concerns
  robustness under input bias, not Hamming or moment stability toward the
  bent/plateaued classes needed here.
