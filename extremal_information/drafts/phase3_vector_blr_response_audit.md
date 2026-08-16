# Independent audit: vector-Hamming synchronization and future responses

**Status.**  Independently reconstructed and proved.  The finite searches in
[`verify_phase3_vector_blr_response_audit.py`](../experiments/verify_phase3_vector_blr_response_audit.py)
check the scalar constants, the simultaneous vector conclusion, all
cycle-contractive maps on `F_2^3` up to irrelevant linear coordinates, and
the all-context word-metric corollary for every spanning binary support
through ambient dimension three.

This note deliberately does not use the proof in
`phase3_geodesic_synchronization.md`.  Its purpose is to decide whether the
proposed synchronization and the constant `11` may safely be used.

## 1. The vector-Hamming statement

Let `Q=F_2^k`, and give `F_2^D` Hamming weight `|.|`.  For

```math
f:Q\longrightarrow F_2^D,qquad f(0)=0,
```

put

```math
\partial f(a,b)=f(a)+f(b)+f(a+b).                \tag{VA.1}
```

### Theorem VA.1 (joint BLR synchronization)

If

```math
|\partial f(a,b)|\le\delta\qquad(a,b\in Q),      \tag{VA.2}
```

then there is a linear map `L:Q to F_2^D` for which

```math
\mathbb E_x|f(x)+L(x)|\le\delta                  \tag{VA.3}
```

and, for the **same** `L`,

```math
\max_x|f(x)+L(x)|\le3\delta.                    \tag{VA.4}
```

The average constant one is asymptotically sharp.  No claim of sharpness is
made for the uniform constant three.

### Scalar BLR constant, with normalization checked

First let `h:Q to F_2`, `h(0)=0`, and write

```math
\delta_h=\Pr_{a,b}[h(a)+h(b)+h(a+b)=1].          \tag{VA.5}
```

Set `g=(-1)^h` and normalize the Fourier transform by

```math
\widehat g(\alpha)=\mathbb E_x g(x)(-1)^{\alpha\cdot x}.
```

Orthogonality and Parseval give

```math
1-2\delta_h
=\mathbb E_{a,b}g(a)g(b)g(a+b)
=\sum_\alpha\widehat g(\alpha)^3,
\qquad
\sum_\alpha\widehat g(\alpha)^2=1.             \tag{VA.6}
```

Let `M=max_alpha widehat g(alpha)`.  Since
`sum_alpha widehat g(alpha)=g(0)=1`, one has `M>0`.  For every real
coefficient `c=widehat g(alpha)`, `c^3<=M c^2`; hence (VA.6) implies

```math
1-2\delta_h\le M.                               \tag{VA.7}
```

Choose a maximizing character `chi_alpha`.  The relative Hamming distance
between `h` and the corresponding linear Boolean function is

```math
\Pr_x[h(x)\ne\alpha\cdot x]={1-M\over2}\le\delta_h. \tag{VA.8}
```

This also covers `delta_h>=1/2`; no sign assumption on the left side of
(VA.7) was used.  The condition `h(0)=0` is what guarantees a positive
Fourier coefficient and hence a *linear*, rather than complemented affine,
approximant.

The constant in (VA.8) cannot be reduced uniformly.  If `k` is even and
`h` is a nondegenerate quadratic bent function with `h(0)=0`, then every
Fourier coefficient has modulus `2^{-k/2}`.  Consequently

```math
\operatorname{dist}(h,\mathrm{Lin})
 ={1-2^{-k/2}\over2},
\qquad
\delta_h={1-2^{-k}\over2},                      \tag{VA.9}
```

and their ratio tends to one.

The joint vector average constant in (VA.3), measured against the
**pointwise** defect budget, is also sharp.  Let `phi` be the same bent
quadratic and put `c(q)=1_(q ne 0)`.  For

```math
f(q)=(\phi(q),c(q)+\phi(q))\in F_2^2,
```

every nondegenerate additive triangle has defect weight exactly one and
every degenerate triangle has defect zero, so (VA.2) holds with `delta=1`.
The two coordinate minimizations separate.  Their nearest linear distances
are respectively

```math
2^{k-1}-2^{k/2-1}
\quad\hbox{and}\quad
2^{k-1}-2^{k/2-1}-1.
```

Consequently every linear vector map has average error at least

```math
1-2^{-k/2}-2^{-k},
```

which tends to the right side `delta=1` of (VA.3).

### Proof of Theorem VA.1

Apply (VA.8) separately to the coordinate functions `f_j`.  If

```math
\delta_j=\Pr_{a,b}[\partial f_j(a,b)=1],
```

the selected linear coordinate `L_j` has error probability at most
`delta_j`.  Keeping the coordinates together when summing gives

```math
\mathbb E_x|f(x)+L(x)|
 \le\sum_j\delta_j
 =\mathbb E_{a,b}|\partial f(a,b)|
 \le\delta,                                     \tag{VA.10}
```

which proves (VA.3).

Put `e=f+L`.  Since `L` is linear, `partial e=partial f`.  For fixed `x`
and every `y`, the Hamming triangle inequality gives

```math
|e(x)|
\le |\partial e(x,y)|+|e(y)|+|e(x+y)|.           \tag{VA.11}
```

Average over `y`, use (VA.2), (VA.3), and translation invariance of the
uniform measure.  This yields `|e(x)|<=delta+2delta=3delta` for every `x`.
The linear map in (VA.4) is therefore exactly the one already selected for
(VA.3).  `square`

The proof is best viewed as an `l_infinity(l_1)` stability statement for a
binary one-cochain.  It is stronger than applying a scalar worst-case bound
to each output channel, because (VA.10) uses the joint pointwise budget in
(VA.2).

## 2. From a diametral chart to a linear graph support

Let `G=F_2^w`, let `S subseteq G\{0}` span `G`, and suppose

```math
B=\{b_1,\ldots,b_D\}\subseteq S,
\qquad W=\operatorname{span}B,
\qquad Q=G/W                                      \tag{VA.12}
```

is a geodesic chart in the following exact sense: whenever
`R subseteq S\B` has projected sum zero,

```math
\left|\sum_{s\in R}s\right|_B\le |R|.           \tag{VA.13}
```

Here `|.|_B` is Hamming weight in the basis `B`.  Assume also that every
nonzero element of `Q` occurs in `pi(S\B)`.

Choose one `s_q in S cap pi^{-1}(q)` for every nonzero `q`, put `s_0=0`,
and choose an arbitrary linear section `L_0:Q to G`.  In `B` coordinates,

```math
f(q)=s_q+L_0(q)\in W.                            \tag{VA.14}
```

For nonzero distinct `a,b`, the three selected generators over
`a,b,a+b` form a projected-zero set, so (VA.13) gives
`|partial f(a,b)|<=3`.  Degenerate pairs have zero defect.  Theorem VA.1
therefore supplies a linear section `L=L_0+A` satisfying

```math
|s_q+L(q)|_B\le9\qquad(q\in Q).                 \tag{VA.15}
```

Define the synchronized graph support

```math
\Gamma_L=B\cup\{L(q):q\in Q\setminus\{0\}\}.   \tag{VA.16}
```

For a spanning support `C`, write `lambda_C(x)` for its Cayley word length.

### Theorem VA.2 (uniform future-response comparison)

For every arbitrary appended support `U subseteq G\{0}`,

```math
\lambda_{S\cup U}(x)-9
\le\lambda_{\Gamma_L\cup U}(x)
\le\lambda_{S\cup U}(x)+11                     \tag{VA.17}
```

for every `x in G`.  In particular,

```math
\sup_U\|\lambda_{S\cup U}-\lambda_{\Gamma_L\cup U}\|_\infty
\le11,                                          \tag{VA.18}
```

and the covering-radius responses differ by at most `11` under every future
append operation.

#### Proof

Fix a shortest representation of `x` from `S union U`, leave all generators
assigned to `U` unchanged, and let `R subseteq S` be the remaining part.
Write `R=I disjoint-union R_0`, where `I subseteq B` and
`R_0 subseteq S\B`, and put `q=sum_{s in R_0}pi(s)`.

If `q=0`, (VA.13) represents `sum R` using at most
`|I|+|R_0|=|R|` members of `B subseteq Gamma_L`.  If `q ne0`, apply
(VA.13) to the symmetric difference `R_0 triangle {s_q}` and then use
(VA.15):

```math
\left|\sum_{s\in R}s+L(q)\right|_B
\le |I|+(|R_0|+1)+9
=|R|+10.                                         \tag{VA.19}
```

One graph generator `L(q)` plus those basis generators replaces `R` at a
cost increase at most `11`.  This proves the right inequality in (VA.17).

Conversely, in any representation from `Gamma_L union U`, all graph
generators not assigned to `U` can be consolidated: their sum is either
zero or one element `L(q)`, because `L` is linear and (VA.16) contains every
nonzero quotient label.  In the latter case replace `L(q)` by `s_q` and at
most nine members of `B`, using (VA.15).  The cost increases by at most
nine.  This proves the left inequality.  The argument never assumes that
`U` spans or is disjoint from either support. `square`

The two additive units in `11=3*3+2` are real steps in this proof: one
closes an arbitrary quotient sum by `s_q`, and one pays for the consolidated
graph generator.  The audit found no justification for deleting either
unit from the theorem.  It also found no small example attaining `11`, so
this is a certified bound, not a sharpness assertion.

## 3. Adversarial finite search

The verifier performs five independent checks.

1. It enumerates every scalar truth table with `f(0)=0` through `k=4` and
   checks the integer form of (VA.8).
2. Modulo addition of linear coordinates, it enumerates all vector maps for
   `(k,D)=(2,D)`, `D<=4`, and `(k,D)=(3,D)`, `D<=3`.  It constructs the
   Fourier-maximizing `L` and checks (VA.3)--(VA.4) simultaneously.
3. For `Q=F_2^3`, it enumerates every multiset of nonlinear scalar classes
   allowed by all additive-cycle inequalities.  Linear coordinates are
   irrelevant and may be added freely.  There are `2,636` admissible
   multisets; their largest optimal uniform distance from a linear map is
   two.  This does **not** prove a better general constant, but it found no
   hidden small counterexample.
4. It enumerates every spanning support, every geodesic chart with complete
   quotient projection, every transversal, every appended support, and
   every target through `w=3`.  It constructs the section by the audited BLR
   rule and verifies both directional constants in (VA.17).
5. It repeats the all-future enumeration for incomplete quotient
   projections, verifying `2h`, `10h+1`, and `10h-1`, and separately checks
   the sharp selector family through `h=5`.

The finite data suggest that `3delta` and `11` are conservative.  Improving
them is a separate `l_infinity(l_1)` cochain-stability problem; finite data
alone are not grounds for changing the rigorous constants.

## 4. Judgment

The proposed theorem and response corollary are valid with the stated
constants.  The substantive new fact is not scalar BLR itself.  It is the
chain

```text
one joint Hamming budget on every additive triangle
    -> one linear section with bounded pointwise vector error
    -> a constant-error all-future word-profile quotient.
```

The resulting state is strictly smaller than the original transversal: it
retains a linear section rather than one unrestricted point in every affine
fibre.  Theorem VA.2 **does** erase arbitrary membership choices inside the
complete fibres: the whole raw support, not merely the chosen transversal,
is within all-context distance `11` of the linear graph.  What is not proved
is an exact summary-only composition law, or a constant-error result when
the quotient projection is sparse and has growing diameter.

## 5. Audit of fibre stripping and incomplete projection

The originating note also proposes two extensions.  They are valid after
making their endpoint hypotheses explicit.

### Theorem VA.3 (fibre stripping in an abelian extension)

Let

```math
0\longrightarrow K\longrightarrow G\overset\pi\longrightarrow Q
\longrightarrow0
```

be an extension of finite abelian groups.  Let `B=-B` generate `K`, let
`A=-A subseteq G\K`, and assume that

```math
\ell_B(a_1+\cdots+a_r)\le r                    \tag{VA.20}
```

for every word in `A` whose projected sum is zero.  Assume
`P=pi(A)` generates `Q`, let its Cayley diameter be `h`, and select lifts
`s_p in A` with `s_{-p}=-s_p`.  Put `A^tr={s_p:p in P}`.  Then, for every
symmetric future generator set `T`,

```math
0\le\ell_{A^{tr}\cup B\cup T}(x)
       -\ell_{A\cup B\cup T}(x)\le2h           \tag{VA.21}
```

for every `x in G`.  The factor two is attained for every `h>=1` in a
binary extension.

#### Proof audit

Take the `r` source letters in an `A union B union T` word, with quotient
sum `q`.  Write `q` using `m<=h` labels of `P`.  The original letters
together with the negatives of the selected lifts form a projected-zero
`A`-word of length `r+m`; sign compatibility is exactly what ensures that
the negative lifts are still in `A`.  By (VA.20), its kernel sum costs at
most `r+m` letters of `B`.  Hence the original source block is replaced by
`m+(r+m)<=r+2h` letters.  The lower inequality is inclusion.

For sharpness, take `Q=F_2^h`, `K=F_2^{2h}`, kernel basis
`e_1,...,e_(2h)`, and two lifts over each quotient basis vector,

```math
c_i=(0,q_i),\qquad c_i'=c_i+e_{2i-1}+e_{2i}.
```

Select `c_i`.  The point `sum_i c_i'` costs `h` with all lifts and exactly
`3h` after stripping.  Quotient independence forces the `h` selected lifts,
and then all `2h` kernel coordinates.  Thus the gap is `2h`.

More generally, include the perturbed lift `c_i'` only for `i in J`.  At the
root `x_P=sum_(i in P)c_i'`, quotient independence and disjoint kernel
coordinate pairs give

```math
\ell_{S_J}(x_P)=|P|+2|P\setminus J|.
```

It follows by choosing `P` in either directional set difference that the
all-root response distance between `S_J` and `S_K` is exactly

```math
2\max\{|J\setminus K|,|K\setminus J|\}.
```

Thus the sharp pair belongs to an `h`-bit selector family; it is not an
isolated endpoint construction.

No commutativity or sign issue is hidden in this proof: commutativity is used
to group source letters, symmetry permits inverse words, and the compatible
choice of lifts is a genuine hypothesis.  In exponent two it is automatic.
For the binary specialization `A=S\setminus B` is indeed disjoint from the
kernel: applying (VA.13) to a singleton in `S cap W` forces its `B`-weight
to be one, hence that element was already in `B`.

### Theorem VA.4 (partial-projection linear summary)

In the binary setup (VA.12)--(VA.13), let

```math
P=\pi(S\setminus B),\qquad h=\max_{q\in Q}\ell_P(q).
```

Assume `Q ne 0`, so `h>=1`.  Select one source lift for every `p in P`.
For each `q`, take a shortest selected word `C_q` of length at most `h`,
take `C_0` empty, and let `g(q)=sum C_q`.  There is a linear section `L`
such that

```math
\max_q|g(q)+L(q)|_B\le9h.                       \tag{VA.22}
```

For `Gamma_L` as in (VA.16), every future `U` satisfies

```math
\lambda_{\Gamma_L\cup U}(x)
 \le\lambda_{S\cup U}(x)+10h+1,
\qquad
\lambda_{S\cup U}(x)
 \le\lambda_{\Gamma_L\cup U}(x)+10h-1.         \tag{VA.23}
```

If `Q=0`, cycle contraction forces `S=B`, and the exact statement has zero
error instead.

#### Proof audit

Choose an arbitrary linear section `L_0` and put
`f(q)=g(q)+L_0(q) in W`.  The symmetric difference of
`C_x,C_y,C_(x+y)` is a projected-zero source set of size at most `3h`, so
`|partial f(x,y)|<=3h`.  Theorem VA.1 gives a linear `A:Q to W` with uniform
error at most `9h`; `L=L_0+A` proves (VA.22).  This explicit `L_0` step is
needed to type the BLR application and to guarantee that `L` is a section.

To convert a source block `R` of quotient sum `q ne0`, compare it with
`C_q`.  Cycle contraction and (VA.22) leave at most

```math
|R|+h+9h
```

kernel letters, and the one graph generator costs the final `+1`.  This is
the first inequality in (VA.23).  Conversely, consolidate any nonempty
graph block to one `L(q)` and replace it by the at most `h` selected source
letters in `C_q` plus at most `9h` kernel letters.  Relative to the at least
one graph letter removed, the overhead is at most `10h-1`.  A zero quotient
block is deleted.  These replacements leave the arbitrary future word
unchanged.

The conditions `P` generates `Q`, `C_0` is empty, and either `Q ne0` or a
separate zero-error endpoint are essential for a literal statement.  In
particular, applying `10h-1` at `h=0` would give the false bound `-1`; the
originating draft was corrected during this audit.
