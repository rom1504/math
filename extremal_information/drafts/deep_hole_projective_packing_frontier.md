# Deep-hole projective packing: exact floor, generic countermodel, and the remaining cut-code lemma

**Status:** theorem draft for independent audit.  This note makes no
canonical edit.  Statements PP.1--PP.4 are proved below.  They do **not**
prove a fixed-scale projective packing for exact minimizers.

## Executive verdict

For the augmented cut code, deepest-coset stability does force more than one
or two isolated nearest representatives:

* the radius-`rho+1` shell of every deep hole contains at least three
  distinct projective codewords;
* these three are automatically separated by the projective minimum distance
  `n-1`, i.e. only by normalized scale `4/n`;
* if three of their error supports already cover all coordinates, those three
  are separated by at least `M_n-2`, i.e. normalized scale
  `Theta(n^(-1/2))`.

Neither scale survives the error term in the AO.2 physical compiler.  The
desired fixed-`gamma` conclusion remains open.

There is also a sharp warning against importing a theorem for generic linear
codes.  An explicit family of antipodal linear codes has length `N`,
covering-radius deficit `Theta(N^(3/4))` (the analogue of `Theta(n^(3/2))`
when `N=Theta(n^2)`), and deep-hole shells of every
`o(N^(3/4))` excess contained in two projective caps of total diameter
`o(N)`.  That construction has dimension `Theta(N^(3/4))`, not the
augmented cut code's `Theta(sqrt N)` dimension.  A sphere-covering argument
proves that the same separable two-cap construction cannot have both
dimension `O(sqrt N)` and deficit `Theta(N^(3/4))` unless the cap support is
already `Theta(N)`.  Thus a successful argument must use at least the low
dimension or the weight-two cut geometry; generic deep-hole theory alone
cannot settle the question.

## 1. Code notation

Let `C <= F_2^N` be a binary linear code containing the all-one word
`mathbf 1`.  For a received word `y`, put

```math
rho=d(y,C),
\qquad
\mathcal L_s(y)=\{c\in C:d(y,c)\le rho+s\}.
```

The projective distance is

```math
d_{\rm pr}(c,c')=\min\{d(c,c'),N-d(c,c')\}.
```

For the augmented cut code `C_n`, `N=binom(n,2)`,

```math
rho={N-M_n\over2},
```

and `mathcal L_s(y)` is the positive energy shell of deficit `2s`.  Its
nonzero projective codeword weights are

```math
\min\{k(n-k),N-k(n-k)\},\qquad 1\le k\le n-1,
```

so, for `n>=5`, the projective minimum distance is exactly `n-1`.

## 2. One-coordinate stability forces three projective shell points

### Theorem PP.1 (the exact radius-one projective floor)

Let `y` be a deep hole of a binary linear code `C <= F_2^N`, so `rho` is the
covering radius of `C`.  Then

```math
|\mathcal L_1(y)|
\ge
\left\lceil{N\over rho+1}\right\rceil.             \tag{PP.1}
```

If `mathbf 1 in C` and `rho+1<N/2`, no two members counted in (PP.1) are
projectively identical.  In particular, for the augmented cut code, whenever
`n>=5` and `M_n>2`, the positive deficit-two shell of every exact minimizer contains at
least three projective words, pairwise separated by at least `n-1` edges.
Equivalently it is a three-point AO.1 packing with

```math
gamma_n={2(n-1)\over\binom n2}={4\over n}.          \tag{PP.2}
```

There is a useful strengthened alternative.  If three members
`c_1,c_2,c_3 in mathcal L_1(y)` have error supports

```math
N_i=\operatorname{supp}(y+c_i)
```

which cover all `N` coordinates, then

```math
d_{\rm pr}(c_i,c_j)\ge N-2(rho+1)
\quad(i\ne j).                                    \tag{PP.3}
```

For `C_n`, the right side is `M_n-2`.  Hence a three-member subcover would
give a three-point packing at normalized AO scale

```math
gamma_n={2(M_n-2)\over N}=Theta(n^{-1/2}).         \tag{PP.4}
```

#### Proof

Fix a coordinate `e`.  Since `y` is deepest,

```math
d(y+e,C)\le rho.
```

Choose `c_e` attaining this distance and write
`N_e=supp(y+c_e)`.  The coordinate `e` must belong to `N_e`: otherwise
`d(y,c_e)=d(y+e,c_e)-1<rho`, contradicting the definition of `rho`.
Moreover

```math
|N_e|-1=d(y+e,c_e)\le rho,
```

so `c_e in mathcal L_1(y)`.  Thus the error supports of
`mathcal L_1(y)` cover all coordinates.  Every such support has size at most
`rho+1`, proving (PP.1).

If `c' = c+mathbf 1`, then

```math
d(y,c')=N-d(y,c)\ge N-(rho+1)>rho+1.
```

Thus at most one member of each projective pair lies in the shell.  For
`C_n`, `rho+1<N/2` is exactly `M_n>2`; then
`N/(rho+1)>2`, so (PP.1) gives at least three distinct projective words.
The projective minimum-distance calculation gives (PP.2).

For (PP.3), fix a pair `{i,j}` and let `k` be the remaining index.  Since
the three supports cover,

```math
N
\le |N_i\cup N_j|+|N_k|
\le |N_i|+|N_i\mathbin\triangle N_j|+|N_k|
\le 2(rho+1)+d(c_i,c_j).
```

Therefore `d(c_i,c_j)>=N-2(rho+1)`.  In the other direction,

```math
d(c_i,c_j)=|N_i\mathbin\triangle N_j|
\le |N_i|+|N_j|\le2(rho+1),
```

so also `N-d(c_i,c_j)>=N-2(rho+1)`.  Taking the smaller of the two proves
(PP.3).  For the cut code,
`N-2(rho+1)=M_n-2`. `square`

### What PP.1 does and does not say

The first conclusion is a generic deepest-coset fact; only its projective
minimum distance uses the augmented cut code.  It is an exact improvement
over mere nonemptiness, but it is far below the frozen target:

```math
gamma_n=4/n
```

in (PP.2), and even the three-subcover alternative gives only
`gamma_n=Theta(n^(-1/2))`.  The sparse-flip/free-shore theorem AO.2 has an
`O(n^(5/4))` simultaneous cut-control remainder.  Its signal at (PP.4) is
only `Theta(n)`, so neither result yields a physical contextual packing.

The proof also explains why a growing shell count alone is insufficient.
Single-coordinate stability says that the near-nearest error supports cover
the coordinate set, but many codewords lying in two narrow projective caps
can perform that cover.

## 3. A scalable two-cap family for generic antipodal linear codes

The next construction is deliberately **not** the augmented cut code.  It
shows exactly which generic coding hypotheses fail to force projective
packing.

### Theorem PP.2 (generic deep holes may have collapsed thin shells)

There is an infinite family of binary linear codes `C_m <= F_2^{N_m}`
containing `mathbf 1` with

```math
N_m\longrightarrow\infty,
\qquad
rho(C_m)={N_m\over2}-Theta(N_m^{3/4}),              \tag{PP.5}
```

and a deep hole `y_m` such that, for every `s_m=o(N_m^{3/4})`, the whole
shell `mathcal L_(s_m)(y_m)` has projective diameter `o(N_m)`.  In
particular it has no two-point, hence no three-point, fixed-`gamma`
projective packing.

The dimensions are

```math
\dim C_m=Theta(N_m^{3/4}),                           \tag{PP.6}
```

so this does not falsify the augmented cut-code statement.

#### Construction and proof

Let `m` range through the positive even integers, and put

```math
N_m=m^4,\qquad D_m=m^3,\qquad L_m=N_m-D_m.
```

Both `D_m` and `L_m` are even.  Put

```math
C_m=\operatorname{Rep}_{L_m}\oplus\mathbb F_2^{D_m}.
                                                               \tag{PP.7}
```

The second factor is the full code.  Covering radii add under direct sums,
hence

```math
rho(C_m)=L_m/2=N_m/2-D_m/2.
```

Since `D_m=N_m^(3/4)`, this proves (PP.5), and
`dim C_m=D_m+1` proves (PP.6).  The code contains the global all-one word.

Take a balanced word `y_A in F_2^(L_m)` and any
`y_B in F_2^(D_m)`; put `y=(y_A,y_B)`.  This is a deep hole, and its two
nearest `C_m` words are

```math
(0_{L_m},y_B),\qquad (\mathbf 1_{L_m},y_B).          \tag{PP.8}
```

They are distinct projective lines at projective distance
`D_m=o(N_m)`.

More generally, the shell is exactly

```math
\mathcal L_s(y)=
\{(r,b):r\in\{0_{L_m},\mathbf1_{L_m}\},
                 d(b,y_B)\le s\}.                   \tag{PP.9}
```

Words with the same repetition coordinate are therefore at mutual distance
at most `2s`.  Words with opposite repetition coordinates have projective
distance at most `D_m`.  Consequently

```math
\operatorname{diam}_{\rm pr}\mathcal L_s(y)
\le\max\{D_m,2s\}.                                  \tag{PP.10}
```

For `s=o(N_m^(3/4))`, the right side is `o(N_m)`, proving the claim.
`square`

The error supports have a particularly transparent form.  Partition the
coordinates into `U,V,H`, where `U` and `V` are the two halves of the
repetition block and `H` is the `D_m`-coordinate full-code block.  The
radius-`rho+r` error supports are exactly

```math
U\cup T\quad\hbox{and}\quad V\cup T,
\qquad T\subseteq H,\quad |T|\le r.                \tag{PP.10a}
```

For every `r`-set `F`, choose `T=F\cap H`.  If `h=|T|`, one of the two
supports in (PP.10a) meets `F` in at least

```math
h+{r-h\over2}={r+h\over2}.
```

Its excess radius is exactly `h`, so this is precisely the complete
deep-hole certificate
`|F cap N_c| >= (r+d(y,c)-rho)/2`.  The two-cap family therefore satisfies
all exact simultaneous `r`-flip majority-cover constraints, at every radius,
not merely the one-coordinate cover used in PP.1.

This construction has all the superficial features that a generic coding
argument might use: linearity, an antipodal code, a deepest coset, and the
correct `N^(3/4)` covering-radius deficit.  Yet its near-nearest landscape is
projectively collapsed.  Therefore any successful theorem must use at least
the augmented cut code's low dimension or its special cut-weight algebra.

## 4. Why the direct-sum countermodel cannot simply be made low-dimensional

### Theorem PP.3 (sphere-covering barrier for separable two-cap models)

Let

```math
C=\operatorname{Rep}_L\oplus B,
```

where `B` is a binary `[D,k]` code with covering radius

```math
rho(B)=D/2-t.
```

Assume `k>=1` and `t>=0`.  Then necessarily

```math
\boxed{D\ge {2t^2\over k\log2}.}                   \tag{PP.11}
```

Consequently, in the augmented-cut scaling

```math
N=Theta(n^2),\qquad t=Theta(n^{3/2}),\qquad k=O(n),
```

every such separable two-cap construction has `D=Omega(n^2)=Omega(N)`.
Thus it cannot inherit PP.2's particular `D=o(N)` projective-diameter
certificate.  This does not by itself exclude a different source of
projective shell collapse inside the `B` factor.

#### Proof

The radius-`rho(B)` balls about the `2^k` codewords cover the `D`-cube, so

```math
2^k\sum_{j\le D/2-t}\binom Dj\ge2^D.               \tag{PP.12}
```

For `X~Bin(D,1/2)`, Hoeffding's inequality gives

```math
2^{-D}\sum_{j\le D/2-t}\binom Dj
=\Pr\{X\le D/2-t\}
\le\exp(-2t^2/D).                                   \tag{PP.13}
```

Combining (PP.12)--(PP.13) yields
`k log2>=2t^2/D`, which is (PP.11). `square`

The theorem rules out only a **separable** localization mechanism.  The
augmented cut code has all `n-1` independent cut directions supported across
the same complete edge set; a hierarchical two-cap shell could exploit that
overlap and is not reduced to a code on a small coordinate block.  Thus
PP.3 is evidence for the plausibility of the cut-specific lemma, not a proof
of it.

## 5. Exact anatomy of the surviving two-cap obstruction

The augmented cut code does add one rigid fact absent from PP.2: in the
asymptotic `D=o(N)` regime, a pair of projectively close but oppositely
oriented positive shell words can agree only on a complete bipartite cut,
and that cut carries the whole energy bias.  The quantitative hypothesis
which excludes the finite balanced-cut alternative is stated below.

### Proposition PP.4 (bipartite agreement core)

Let `a` be an exact minimizer, `Q(a)=M_n`, let `y` be its binary received
word under the code dictionary, and let `z_0,z_1` be two positive oriented
shell words of deficits at most `2s`.  Let their corresponding codewords be
`c_0,c_1`, and suppose their **actual** code distance is the large lift

```math
d(c_0,c_1)=N-D,
\qquad D=d_{\rm pr}(c_0,c_1)<N/2.                    \tag{PP.14}
```

Then

```math
D\ge M_n-2s.                                         \tag{PP.15}
```

More precisely, write `N_i=supp(y+c_i)` and split their `D` agreement
coordinates into

```math
B=N_0\cap N_1,
\qquad
Z=[N]\setminus(N_0\cup N_1).
```

If `w_i=|N_i|`, then exactly

```math
|B|={w_0+w_1-N+D\over2},
\qquad
|Z|={N-w_0-w_1+D\over2},
\qquad
|Z|-|B|=N-w_0-w_1\ge M_n-2s.                       \tag{PP.16}
```

For the augmented cut code, if additionally

```math
D<N-\lfloor n^2/4\rfloor,                            \tag{PP.16a}
```

then the agreement set `B union Z` is exactly an edge cut `delta(S)` of
`K_n`, with

```math
D=|S|(n-|S|).                                        \tag{PP.17}
```

After switching so that the first oriented word is `z_0=mathbf1`, the
second is `z_1=-c(v_S)`, and the signing bias on this bipartite core is

```math
\sum_{e\in\delta(S)}a_e= {\langle a,z_0\rangle+
                              \langle a,z_1\rangle\over2}
\ge M_n-2s.                                          \tag{PP.18}
```

In particular, an asymptotic two-cap obstruction with `D=o(N)` must have

```math
|S|=o(n)
\quad\hbox{but}\quad
|S|\ge {M_n-2s\over n}=Omega(\sqrt n)               \tag{PP.19}
```

when `s=o(n^(3/2))`.

#### Proof

The two error supports have symmetric difference `N-D`.  Solving

```math
w_0+w_1=2|B|+(N-D),
\qquad
D=|B|+|Z|
```

gives (PP.16).  Since `w_i<=rho+s` and
`2rho=N-M_n`, its last inequality and hence (PP.15) follow.

The difference of two augmented-cut codewords is either a cut or its
complement.  A cut of `K_n` has weight at most `floor(n^2/4)`.  Under
(PP.14) and (PP.16a), the actual difference weight `N-D` is therefore too
large to be a cut and must be the complement of one.  The coordinates on
which the corresponding oriented sign words agree form `delta(S)`, proving
(PP.17).  In the displayed gauge both words
equal `+1` on `delta(S)` and are opposite off it.  Adding their two energy
inequalities therefore gives (PP.18).  Finally
`|S|(n-|S|)=o(n^2)` implies, after replacing `S` by its complement,
`|S|=o(n)`, while `|S|n>=D>=M_n-2s` gives (PP.19).
`square`

Condition (PP.16a) is automatic in the intended asymptotic two-cap regime
`D=o(N)`; it is not automatic from `D<N/2` at small or balanced-cut scales.

PP.4 is the cleanest exact description of what a failure of projective
diffusion would have to look like: a positive and a negative near-ground
configuration differ on only `o(n)` vertices, and an
`Omega(n^(3/2))` signed bias is concentrated on the resulting
`K_(|S|,n-|S|)` interface.  It is not itself a contradiction.  Replacing or
optimizing that interface while preserving cancellation against all other
spins is again a joint bridge problem; paying its channels separately would
repeat an archived obstruction.

## 6. Exact remaining lemma and implication to AO.2

The weakest fixed-scale statement currently visible is:

> **Cut-DH(3).** There are a constant `gamma>0` and a sequence
> `s_n=o(n^(3/2))` such that, for every deep hole `y` of the augmented cut
> code `C_n`, the shell `mathcal L_(s_n)(y)` contains three words
> `c_1,c_2,c_3` with
>
> ```math
> d_{\rm pr}(c_i,c_j)\ge {gamma\over2}\binom n2
> \qquad(i\ne j).
> ```

The factor `gamma/2` is chosen so that the signed edge overlap obeys the AO.1
condition `|R|<=1-gamma`.  A growing version replaces three by
`K_n -> infinity`; an information-rate version asks for
`K_n>=exp(cn)`.

The complete proved implication is

```text
Cut-DH(K_n) with shell excess s_n=o(n^(3/2)) and fixed gamma
  --[exact code/energy dictionary]-->
(d_n,gamma) positive-shell packing, d_n=2s_n=o(n^(3/2))
  --[audited Theorem AO.2, fixed alpha=lambda]-->
K_n all-spins-free exact-sign children at order n+O(sqrt n),
common low-cap query bank, pairwise response gap Theta(n^(3/2)),
total parent caps O(n^(3/2))
  --[response packing]-->
at least log_2 K_n reusable bits for that declared physical language.
```

Thus `K_n=3` supplies a genuine three-state low-cap obstruction,
`K_n -> infinity` supplies unbounded information, and
`K_n=exp(cn)` supplies `Omega(n)` bits.  This is an incompressibility result,
not a convergence theorem.

**Assumption-distance assessment.**  Cut-DH(3) is syntactically and
informationally much weaker than the original minimization: once three
witnesses are exhibited, their shell energies and three pairwise overlaps
are checked in polynomial time, and they retain none of the remaining
`2^(n-1)` energy landscape.  There is no reverse implication from the
packing to the value of `M_n`.  However, present methods prove only the
vanishing-scale PP.1/DH.1 conclusions, and PP.2 shows that deepest-coset
minimality alone cannot prove it.  Its proof must exploit a cut-specific
low-dimensional or weight-two identity not yet identified.  It should
therefore be called a **strictly stated reduction target**, not a proved
strict reduction or a high-confidence lemma.

## 7. Director recommendation

1. Retain Cut-DH(3) as the frozen structural target only if the campaign is
   willing to seek a genuinely cut-specific theorem.  Do not cite generic
   deep-hole, UPWS, or nearest-leader multiplicity results as evidence for it.
2. Record PP.1 as the exact unconditional floor.  Its normalized separation
   is too small for AO.2, so it does not reset the physical-transfer strike.
3. Use PP.2 as the scalable two-cap falsifier for generic coding arguments,
   and PP.3 as the precise reason it does not yet falsify the augmented code.
4. The next discriminating attack should either prove that low-dimensional
   cut directions cannot realize the PP.2 localization nonseparably, or
   construct such a hierarchical augmented-cut countermodel.  More shell
   cardinality estimates at the FB.3 scale will not decide fixed `gamma`.
