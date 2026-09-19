# Chiral follow-up: universal padding, proof-class audit, and remaining gap

This continues `twisted_chiral_adversary_2026_09_18.md` after the usage
interruption. Throughout, `Q(A)=max_x |sum_{i<j} A_ij x_i x_j|`, and
`beta(A)=max_{x,y in {+1,-1}^N} x^T A y`. Matrices are symmetric and
hollow unless explicitly stated. No numerical failure is a lower bound.

## 1. A universal low-cap padding theorem

The earlier opposite-clique/Hadamard construction does not need a
Hadamard baseline. The following applies to **every** bounded-cap sequence.

**Theorem.** Fix `C<infinity` and `epsilon>0`. For every sufficiently
large order-N full signing A with `Q(A)<=C N^(3/2)`, there is a full
signing A' of the same order, obtained by replacing two disjoint principal
blocks of size `m=floor(sqrt(2 epsilon) N^(3/4))` by opposite cliques, with

```math
Q(A')\le Q(A)+\epsilon N^{3/2}+O_{C,\epsilon}(N^{11/8}),
\qquad
\beta(A')\ge\beta(A)+4\epsilon N^{3/2}
 -O_{C,\epsilon}(N^{5/4}).                       \tag{1}
```

The constants are uniform in A. This is an existence theorem, not an
efficient search claim and not a statement that A' is a minimizer.

### Fourth moment of a random principal block

Every column of A belongs to the real cube. Consequently

```math
|(A^3)_{ij}|=|(Ae_i)^T A(Ae_j)|\le\beta(A),
\quad
\operatorname{tr}(A^4)
=\sum_{ij}A_{ij}(A^3)_{ji}\le N(N-1)\beta(A).      \tag{2}
```

Let S be a uniformly random m-subset. The closed length-four walks with
repeated vertices have total weight
`r_m=m(m-1)(2m-3)`: their weights are all +1. The other walks use four
distinct vertices. Therefore, exactly,

```math
\mathbb E\operatorname{tr}(A[S]^4)
=\frac{(m)_4}{(N)_4}\bigl(\operatorname{tr}(A^4)-r_N\bigr)+r_m.
```

In particular, since `beta(A)<=4Q(A)`, this expectation is at most

```math
T:=\frac{m^4}{N^2}\beta(A)+2m^3
\le4m^4Q(A)/N^2+2m^3.                           \tag{3}
```

The trace is nonnegative. On the event `tr(A[S]^4)<=4T`,

```math
Q(A[S])\le\frac m2\|A[S]\|_{op}
\le\frac m2(4T)^{1/4}.                          \tag{4}
```

### Simultaneous retention of a bilinear witness

Fix a maximizing pair x,y with `x^T A y=beta(A)`. Choose ordered
disjoint uniform m-subsets S1,S2, and write C0 for their complement.
Put `F=x[C0]^T A[C0] y[C0]`. Extending the restricted signs by zero
places both arguments in the cube, so `F<=beta(A)` pointwise. Hence
the deficit is nonnegative and

```math
\mathbb E[\beta(A)-F]
=\left(1-\frac{(N-2m)_2}{(N)_2}\right)\beta(A)=:L.
```

Each S_i is a uniform m-subset. Markov's inequality with factor four,
applied to both traces and to this nonnegative deficit, shows that with
probability at least 1/4 all three inequalities hold:

```math
\operatorname{tr}(A[S_i]^4)\le4T\ (i=1,2),
\qquad F\ge\beta(A)-4L.                         \tag{5}
```

Zero expectation cases hold identically. Fix such a choice. This argument
does not presume concentration or independent restrictions.

### Effect of overwriting the blocks

Let A' replace the two old principal blocks by `J_m-I_m` and
`-(J_m-I_m)`. The new opposite-clique energy has absolute maximum
`floor(m^2/2)`: one block is constant and the other as balanced as
parity permits. Therefore

```math
Q(A')\le Q(A)+\lfloor m^2/2\rfloor+m(4T)^{1/4}.  \tag{6}
```

The bilinear cube norm is superadditive over disjoint principal blocks.
To see this, choose a maximizing x/y pair on each block, then multiply
both arguments on each block by an independent common sign. The internal
terms stay fixed and the cross terms average to zero. Some choice of
block signs therefore makes the cross contribution nonnegative. Thus

```math
\beta(A')\ge\beta(A[C0])+2m(m-1)
\ge\beta(A)+2m(m-1)-4L.                         \tag{7}
```

For `m=Theta(N^(3/4))`, (3) gives `T=O(N^(5/2))`, so the error in
(6) is `O(N^(11/8))`. Also `L=O(m beta(A)/N)=O(N^(5/4))`.
Rounding m only costs `O(N^(3/4))`. This proves (1).

The exact trace identity, witness-retention expectation, simultaneous
Markov selection, and finite overwrite inequalities are exhaustively
checked on seven small matrices (including both order-ten minimizing
classes) by `computations/twisted_chiral_padding_audit_2026_09_19.py`.
Its JSON retains the exact rational expectations. These finite checks
support the identities, not the asymptotic quantifiers.

### Optional stronger cap error by decoupling row restrictions

The cap error in (1) can in fact be improved to `O(N^(9/8))`. The
fourth-moment proof above is retained because it is particularly short
and independent of this refinement.

For independent Bernoulli row and column subsets U,V of rates a,b,
respectively, standard ghost-sample symmetrization followed by the
no-outer-absolute Rademacher contraction inequality gives

```math
\mathbb E\beta(P_U A P_V)
\le ab\beta(A)+bN\sqrt{2aN}+aN\sqrt{2bN}.        \tag{10}
```

Here is the reduction, specifying the constants. Conditional on V,
center the selected row indicators and introduce an independent ghost
sample. Their differences have the law `epsilon_i zeta_i`, where
`zeta_i` has mean `2a(1-a)`. Contracting the absolute row-field functions
to their underlying linear fields leaves
`sum_{j in V}|sum_i epsilon_i zeta_i A_ij|`, with expectation at most
`|V|sqrt(2aN)`. The uncentered term is `a beta(A P_V)`. Apply the same
one-sided calculation to V and average to obtain (10). Contraction
has constant one here because there is no absolute value outside the
supremum of the centered sum.

For a principal Bernoulli-p set S, randomly color each of its vertices
into U or V with equal probabilities. For a fixed maximizing x,y on S,
the expected cross evaluation is `beta(A[S])/4`, since A is hollow.
Taking the cross supremum therefore gives

```math
\beta(A[S])\le4\mathbb E_{\rm color}\beta(P_U A P_V).
```

Unconditionally U has rate p/2, and conditional on U the vertices outside
U enter V independently with rate `p/(2-p)`. Extend V by sampling the
vertices inside U at the same rate. The extended V' is independent of
U, and enlarging either rectangular index set cannot decrease its
bilinear cube norm. Both rates are at most p. Formula (10) now implies

```math
\mathbb E\beta(A[S])
\le4p^2\beta(A)+8\sqrt2\,(pN)^{3/2}.            \tag{11}
```

To pass to a uniform m-subset, use a Bernoulli reservoir of rate
`p=2m/N`. For sufficiently large m, the probability that the reservoir
has at least m vertices is at least 1/2. Conditional on its size,
selecting a uniform m-subset preserves the uniform distribution, and
principal restriction does not increase beta. Hence

```math
\mathbb E_{|S|=m}Q(A[S])
\le16m^2\beta(A)/N^2+32m^{3/2}.                 \tag{12}
```

One uses `Q<=beta/2` after comparing the two expectations. Apply the
same factor-four Markov selection to the two principal caps and the
complement witness deficit. For `m=Theta(N^(3/4))`, their summed cap
payment is `O(N^(9/8))`, while the complement deficit remains
`O(N^(5/4))`. Thus the sharper version of (1) is

```math
Q(A')\le Q(A)+\epsilon N^{3/2}+O_{C,\epsilon}(N^{9/8}),
\qquad
\beta(A')\ge\beta(A)+4\epsilon N^{3/2}
 -O_{C,\epsilon}(N^{5/4}).                      \tag{13}
```

In particular the resulting all-twist gap has an explicit
`O(N^(5/4))` error. The argument uses ordinary row/column norm
monotonicity, not a false entrywise monotonicity of signed quadratic
forms.

## 2. Consequences for every signed-permutation twist

For every twisted double built from A', every matching d satisfies

```math
Q(D)\ge\beta(A')-N.
```

This follows from exchanging the two parent spin blocks and then from
`beta(B+diag(d))>=beta(B)-N=beta(A')-N`. Combining with (1), if
`Q(A_N)/N^(3/2)->c` and `beta(A_N)/N^(3/2)->b`, gives the quantitative
curve

```math
\liminf\frac{Q(D_N)-2\sqrt2 Q(A'_N)}{N^{3/2}}
\ge b-2\sqrt2 c+(4-2\sqrt2)\epsilon.             \tag{8}
```

In particular `b>=2c`. More uniformly, suppose a baseline sequence has
`Q(A_N)<= (U+o(1))N^(3/2)`. For any fixed `epsilon>U/sqrt(2)`, (1)
and `beta(A_N)>=2Q(A_N)` force a positive leading gap in (8), without
requiring the baseline normalized caps to converge. The modified sequence
has cap at most `(U+epsilon+o(1))N^(3/2)`.

Consequently the uniform twisted-doubling statement with any
`o(N^(3/2))` defect fails on the class `Q(A)<=C N^(3/2)` whenever

```math
C>(1+1/\sqrt2)U.                                \tag{9}
```

Using the archived all-order upper bound `U=0.493608094` gives a threshold
below 0.84265. The imported upper theorem is recorded in
`flatify_adversary_2026_09_07_ternary_upper_reconstruction.md`; no new
reverification of its analytic chain is claimed here. One may use the
explicit Hadamard upper U=1/2 instead for a completely elementary
baseline, recovering the earlier threshold 0.853553....

This still does **not** refute selection among exact minimizers, nor a
uniform theorem restricted to caps near the present 0.43--0.50 range.
Even starting from exact minimizers, overwriting blocks does not preserve
minimality. Formula (1), rather than the small improvement of the numerical
threshold, is the main new structural statement.

## 3. Independent audit of the factor-four certificate obstruction

The section “Linear cube-pullback certificates have an asymptotic
factor-four floor” in `twisted_chiral_uniform_2026_09_18.md` passes.
For an order-2n twisted core D0, its off-diagonal entrywise mass is
`4n(n-1)`. An order-n full signing pulled back by a real cube contraction
R has total entrywise mass at most `n(n-1)`, including any new diagonal.
Thus an off-diagonal identity

```math
D_0=\sum_r\lambda_r R_r^T A_r R_r+E
```

forces `4n(n-1)<=n(n-1) sum|lambda_r|+|E|_1`. For hollow symmetric E,
fourth-moment interpolation gives
`E|row(E) dot signs|>=||row(E)||_2/sqrt(3)`, and polarization gives
`beta(E)<=4Q(E)`. Hence `|E|_1<=4sqrt(6n)Q(E)` and a remainder
`Q(E)=o(n^(3/2))` requires `sum|lambda_r|>=4-o(1)`.

The hollow condition on E matters: diagonal entries are not controlled
by the off-diagonal Q. The stated theorem includes this condition.
Arbitrary diagonal cancellation in the pullbacks is harmless after
projecting onto off-diagonal entries. This is a limitation of absolute
triangle/coefficient certificates, not of nonlinear maps, joint
cancellation estimates, or the actual twisted matrices.

## 4. Order-ten completion audit

Both complete order-ten censuses cover exactly all `2^36` root-gauged
signings in two disjoint shards. The direct popcount formula and the
independent signed half-edge-table recurrence are correct. Every surviving
matrix is checked on all 512 projective spins. Both runs find minimum 13
and exactly two switching/permutation/global-negation classes; each has
root orbit size and enumerated count 362880.

Both fixed-child family outputs are `complete_family`, with cap 44,
185794560 distinct B matrices, and 9292808 canonical B representatives.
The second class has separate complete normal- and extreme-order runs.
The quotient, stabilizer weights, and matching-bitset exclusions were
audited in the earlier adversarial artifact; the extreme-order version
only reorders the full x domain and translates the full y Gray traversal.

Fresh direct reconstruction and all 524288 projective parent-spin checks
for both final witnesses reproduce cap 44; the output is
`computations/results/twisted_chiral_adversary_n10_recheck_2026_09_19.json`.
Thus **every exact-minimizing order-ten child has twisted-family optimum
44**. This is not a global lower bound for order-twenty signings.

## 5. Global chirality does not yet yield a half-constant theorem

For `D=[[A,C],[C,-A]]`, put `M=A-iC` and
`u=(x+iy)/(1+i)` in `{1,i,-1,-i}^n`. Then exactly

```math
H_D(x,y)=-\operatorname{Im}(u^T M u).
```

This is a single real projection on a four-phase alphabet, not the
modulus supremum of a polynomial on the continuous complex torus.
Bounds for the latter cannot be imported without paying the actual
rounding/projection loss. Spectral symmetry follows from anticommutation
with the canonical quarter-turn, but does not identify Boolean maximizers.

The known strict-upper restricted-weave proof also does not immediately
transfer. Enforcing a macro-pair chiral involution ties paired fiber
bases and row selectors. Its scalar proof needs independent fresh
fiber data, specifically `(E Z)^m`, and does not justify replacing this
by paired dependent averages. A symmetric complex-Hadamard weave naturally
realifies to a chiral sign matrix, but has a different four-phase source
and a real-projection counting problem; no required entropy theorem has
been proved here.

Finite chiral caps below one half, such as Q16=30, do not disprove an
asymptotic half lower bound. Conversely, their tensor amplification is
not automatic: tensoring flat orthogonal factors can create new Boolean
spectral saturators. No asymptotic half lower bound or chiral strict-upper
counterexample is asserted by this follow-up.

## 6. Independent bipolar-model audit

The arbitrary-proportion bipolar theorem in
`twisted_chiral_uniform_2026_09_18.md` passes the independent algebra
audit. For masses p>=q and `T=[[1,t],[t,-1]]`, balancing the switching
inside both blocks gives the exact continuous core cap
`max{p^2+tpq+q|q-tp|, p^2+q^2/2+tpq}`. The four diamond faces have
the claimed maxima, including unequal p,q: expanding the bridge terms
cancels the mixed linear seed terms. On the zero/zero face the branch
choice follows from

```math
p\min(p-c,tq)\ge(q-a)\min(t(p-c),q).
```

The seed formula, factor-four comparison, and `4N||E||+O(N)` finite
rounding error are also correct. This is a positive selection theorem
for a specified quadratic-scale model, not a transfer theorem for
near-minimizing `N^(3/2)`-scale children.

## 7. Balanced bilinear maximizers do not improve the diagonal loss

The universal matching-invariant obstruction is `gamma(A)>=beta(A)-N`.
One might hope that balanced bilinear extremizers improve the loss to
`O(sqrt(N))`. Balance alone does not suffice: the full loss N is attained
on an explicit infinite family whose maximizing products are all balanced.

Let `m=4^k>=64` and

```math
R=(J_4-2I_4)^{\otimes k},\qquad
R^T=R,\quad R^2=mI,\quad R1=\sqrt m\,1.
```

For `a in {1,2}`, define the symmetric, possibly nonhollow matrix

```math
G_a=\begin{pmatrix}J_m-aI_m&R\\R&-J_m+aI_m\end{pmatrix}.
```

Then

```math
\boxed{\beta(G_a)=2m(m-a).}                     \tag{14}
```

Moreover every positive maximizing pair x,y has coordinatewise product
`x*y=(1_m,-1_m)`. Here beta for a nonhollow matrix retains its usual
bilinear definition; no off-diagonal Q is assigned to its diagonal.

### Proof of (14)

The corner pair `x=(1,1), y=(1,-1)` attains `B0=2m(m-a)` because
the two regular-Hadamard cross terms cancel. Write X1,Y1,X2,Y2 for the
four block spin sums of an arbitrary pair. The cross contribution has
absolute value at most `2m sqrt(m)`, while the diagonal correction has
absolute value at most `2am`. Thus a pair with value at least B0 must
satisfy

```math
\delta_1+\delta_2
:=(m^2-X_1Y_1)+(m^2+X_2Y_2)
\le2m\sqrt m+4am<m^2.
```

Both deltas are nonnegative. It follows that X1,Y1 have the same
nonzero sign, and X2,Y2 have opposite nonzero signs. Relative to the
corresponding four constant block vectors, let the numbers of flipped
coordinates be `b1,b2,b3,b4`, and let `h=sum b_i`. Each `b_i<m/2`.
The main internal loss equals

```math
\delta_1+\delta_2=2mh-4(b_1b_2+b_3b_4)\ge mh,
\quad h\le2\sqrt m+4a.
```

The diagonal correction above its corner value is at most `2ah`.
Regularity of R makes the first-order cross correction at most
`2sqrt(m)h`; the second-order correction is at most
`4(b1 b4+b2 b3)`. Since

```math
b_1b_2+b_3b_4+b_1b_4+b_2b_3
=(b_1+b_3)(b_2+b_4)\le h^2/4,
```

the total value minus B0 is at most

```math
h[-2m+2a+2\sqrt m+h]
\le h[-2m+6a+4\sqrt m]<0
```

whenever h>0, for m>=64 and a<=2. Thus all maximizers are the corner
pairs h=0, proving the value and product characterization.

Now A=G1 is a hollow full signing of order `N=2m`. For
`d=(-1_m,+1_m)`, one has `A+diag(d)=G2`, so

```math
\gamma(A)\le\beta(G_2)=2m(m-2)=\beta(A)-N.
```

The reverse inequality is the universal diagonal triangle bound. Hence
`gamma(A)=beta(A)-N` exactly. Every beta-maximizing product has m positive
and m negative coordinates, so is balanced. These are high-cap,
quadratic-scale seeds; the example does not rule out a stronger
diagonal estimate under additional near-minimality assumptions.

## 8. Further independent positive-certificate audits

The director's order-twenty symmetric-Hadamard excess certificate passes.
For a symmetric Hadamard H of order20, switched row sums r all have the
same residue0 or2 modulo4, and `sum r_i^2=400`. Residue2 implies
`sum(r_i-2)(r_i-6)>=0`, hence `sum r_i<=80`. In residue0 put r=4q.
Then `sum q_i^2=25`; Cauchy and parity give `sum q_i<=21`. Equality21
would force `q=1+e_a+e_b-e_c`. For the switched symmetric Hadamard K,
`K1=r` and `Kr=20*1`, yielding
`K(e_a+e_b-e_c)=5*1-r`. Its c-coordinate is5, impossible for a sum
of three signs. Applying the argument to -H controls the negative
energy too. The conference block construction has trace0, so hollowing
it gives cap at most40; the stored attaining witness establishes equality.

The uniform agent's universal coordinate-diamond identity also passes.
After switching b nonnegative, put `h=(b+d-1)_+`, `u=b+d-h`,
`v=b-h`, and `w=d-h`. The exact six-form identity

```math
H(a)-H(c)+b^TAd
=H(u)-H(v)+\tfrac12H(a+h)+\tfrac12H(a-h)
 -\tfrac12H(c+w)-\tfrac12H(c-w)
```

has all six arguments in the cube when `|a|+|b|<=1` and
`|c|+|d|<=1`. Its coefficient absolute sum is4. For the completed
twin clone of an order-n signing, the parent main energy is four times
this expression, and the two twin-edge payments and parent matching
cost at most `2n+2n+2n`. Therefore the claimed explicit clone-class
bound `Q(parent)<=4Q(child)+10n` is correct. It does not identify
these cloned children with near-minimizers.

The later, stronger bilinear diamond certificate also passes:
`Delta(A)<=4 beta(A)/3` for every hollow real symmetric A. Hollow
affinity reduces each coordinate diamond to its vertices; a common
switch makes `a+b=1`, leaving exactly the eight types `(i,j,r)`.
The fourteen-term certificate in
`computations/twisted_chiral_diamond_bilinear_2026_09_19.py` has positive
integer weights totaling32 and satisfies
`sum w(XY^T+YX^T)=24K`. All64 scalar coefficient identities were
independently replayed with integer loops. Arbitrary repetitions of
the eight types lift its sign vectors to any matrix order, and
`F=(1/2)<A,K>=sum(w/24)X^TAY`; hence the coefficient is exactly4/3.
No LP optimizer or floating-point equality is needed for verification.

## 9. Why this padding cannot have vanishing cap cost and fixed loss

There are two distinct notions of “small change”: the difference of two
caps, and the cap of the difference matrix. They cannot be interchanged.

### General norm continuity of the family objective

For hollow symmetric A,A' of the same order, let E=A'-A. The matching
terms cancel when the same signed permutation and diagonal are used in
both doubles. The parent difference has energy
`H_E(x)-H_E(y)+x^T g^T E g y`. Consequently

```math
|F(A')-F(A)|\le2Q(E)+\beta(E)\le6Q(E).           \tag{15}
```

This applies to every perturbation, not just padding. In particular,
changing only `o(N^(3/2))` unit edges cannot create a fixed normalized
twisted-family loss. It does not say that `Q(A')-Q(A)=o(N^(3/2))`
alone prevents such a loss: the perturbation norm may remain large.

### Uniform near-order continuity of the global minimum

Write `M_N=min_A Q(A)`. Extend an optimal order-(N-r) signing by filling
all edges incident to r new vertices with independent random signs.
There are `L=rN-r(r+1)/2` such edges. For every fixed projective spin,
their energy is a sum of L independent signs. Hoeffding's inequality
and the `2^(N-1)` projective spins give

```math
\Pr\{Q(E)>t\}\le2^N\exp[-t^2/(2L)].
```

Choosing `t=sqrt(2L(N log(2)+1))` proves the elementary extension bound

```math
0\le M_N-M_{N-r}\le\sqrt{2L(N\log2+1)}
=O(N\sqrt r).                                  \tag{16}
```

The lower inequality is principal-restriction monotonicity. No
convergence assumption is needed. This is a classical random-extension
argument; no novelty claim is made for it.

Therefore, if `Q(A_N)<=M_N+eta_N N^(3/2)` with `eta_N->0`, **every**
principal set C of size N-r satisfies

```math
Q(A_N)-\eta_NN^{3/2}-O(N\sqrt r)
\le Q(A_N[C])\le Q(A_N).                        \tag{17}
```

This is uniform over the deleted set, not merely an average or a selected
restriction.

### A fixed opposite-clique gain has a fixed near-minimizer cost

Overwrite any two disjoint m-blocks by opposite cliques, and let C be
their complement. Pick a maximizing sign of the old complement energy.
The two cliques have an internal energy of that same sign and magnitude
`floor(m^2/2)`. Independently multiply each of the three block spin
vectors by a common random sign. All internal energies stay fixed and
all cross terms average to zero. Some choice therefore has the desired
signed cross contribution nonnegative. This proves, for every block
choice and every old signing,

```math
Q(A')\ge Q(A[C])+\lfloor m^2/2\rfloor.           \tag{18}
```

For an additive-o(N^(3/2)) near-minimizer and
`m~sqrt(2 epsilon)N^(3/4)`, (17)--(18) give

```math
Q(A')\ge Q(A)+\epsilon N^{3/2}-o(N^{3/2}).       \tag{19}
```

Thus the fixed bilinear gain in the existing padding cannot be made
free by a more cleverly targeted choice of its two blocks. Conversely,
if `m=o(N^(3/4))`, its perturbation has `Q(E)=O(m^2)=o(N^(3/2))`,
and (15) forbids a fixed twisted-family change at all.

For arbitrary bounded-cap baselines, the selected padding construction
can also be made sharp: retain a Q-maximizing complement witness as a
fourth nonnegative-deficit constraint in the Markov argument. A factor
five in place of four leaves positive probability. Equations (6), (7),
and (18), together with the triangle upper bound on beta of the modified
blocks, then give

```math
Q(A')=Q(A)+\epsilon N^{3/2}+O(N^{5/4}),
\qquad
\beta(A')=\beta(A)+4\epsilon N^{3/2}+O(N^{5/4}). \tag{20}
```

The constants depend only on the bounded-cap constant and fixed epsilon.
This sharp curve explains the limitation of the opposite-clique method.
It does not rule out a one-sided perturbation using pre-existing energy
imbalance; the next section makes that distinction explicit.

There is also a general concentration corollary. The elementary block
inequality used here is already recorded in
`macroscopic_closure_block_dichotomy.md`, Section2; the present consequence
combines it with the near-order minimum bound. For any deleted set R,
write `U_R=max H_A[R]` and `V_R=-min H_A[R]`. Common flips give
`Q(A)>=Q(A[C])+min(U_R,V_R)`. Combining with (17), every r-set in an
additive-eta near-minimizer obeys

```math
\min(U_R,V_R)\le\eta_NN^{3/2}+\sqrt{2L(N\log2+1)}. \tag{21}
```

Thus an o(N)-vertex block cannot carry fixed positive **and** negative
`N^(3/2)` energies near optimality. Two disjoint microscopic blocks with
opposite hot energies are likewise excluded, by aligning their separate
endpoint energies before the common-flip argument. One-sided hot energy
alone is not excluded.

## 10. Conditional zero-cost padding from one-sided energy imbalance

Let `U(A)=max H_A` and `V(A)=-min H_A`, and suppose a bounded-cap
sequence has limits

```math
U(A_N)/N^{3/2}\to u,\quad V(A_N)/N^{3/2}\to v,
\quad\beta(A_N)/N^{3/2}\to b,
\qquad u\ge v.
```

For every fixed epsilon>0 there are full signings A' obtained by
overwriting **one** principal block of size
`m=floor(sqrt(2 epsilon)N^(3/4))` with a negative clique such that

```math
U(A')=U(A)+O(N^{5/4}),
\quad V(A')=V(A)+\epsilon N^{3/2}+O(N^{5/4}),
\quad\beta(A')=\beta(A)+2\epsilon N^{3/2}+O(N^{5/4}). \tag{22}
```

The constants depend only on the bounded-cap constant and epsilon.

**Proof.** Select a uniform m-block with small old principal beta
(equation(12), or its beta version), and retain in its complement
separate maximizing witnesses for U,V,beta. Each witness deficit is
nonnegative and has expectation `O(m sqrt(N))`. The four Markov
constraints have a simultaneous solution using factor five. The old
block has beta `O(N^(9/8))` and the three retained witnesses lose
`O(N^(5/4))`.

The negative clique has positive endpoint `floor(m/2)` and negative
endpoint `m(m-1)/2`. Subtracting the old small block gives the upper
bounds in (22). Combining the complement and clique endpoint witnesses
by an independent common flip gives the lower bounds. Bilinear
superadditivity gives the beta lower bound, while the triangle inequality
and the old block's small beta give its matching upper bound. Rounding
m costs only `O(N^(3/4))`.

If `u>v`, take epsilon=u-v. Then

```math
Q(A')/N^{3/2}\to u,\qquad
\beta(A')/N^{3/2}\to b+2(u-v).                  \tag{23}
```

Thus this one-sided balancing has **zero leading cap cost**. If the
baseline additionally satisfies

```math
b+2(u-v)>2\sqrt2\,u,                            \tag{24}
```

every twisted double of A' violates the target factor by a positive
`N^(3/2)` amount. This is a precise conditional mechanism, not an
assertion that the required baseline near-minimizers exist. If the
baselines are additive-o(N^(3/2)) near-minimizers, so are A', but the
modification does not preserve **exact** minimality.

The presently certified half-range lower bound only gives
`u+v>=2 ell`, with `ell=0.4333221116640807`, while the reported upper
range gives `u<=0.493608094`. Hence `u-v<=0.1205719646718386`.
That leaves room for a nonzero balancing gain but does not imply (24).
The trivial bound `b>=2u` is insufficient throughout this range. The
unresolved quantitative question is the **joint** behavior of beta and
one-sided imbalance, not either scalar statistic in isolation.

### Uniform beta-accounted balancing, including N-dependent imbalance

**Archive distinction.** Same-order endpoint balancing is already proved
more strongly in `transfer_director_same_order_orientation_repair_2026_09_06.md`
(gap O(N)) and `principle_synthesis_2026_09_07_global_balancing.md`
(gap at most2 with operator-norm control), both with cap cost
`O(N^(5/4))`. Thus the balancing assertion below is a weaker reproduction,
not new preparation progress. The additional feature relevant here is
the explicit **increase of beta by twice the old endpoint imbalance**,
and the resulting necessary chiral condition (24). The clique version
does not inherit the older global repair's operator-norm control.

The preceding argument does not need convergent endpoint ratios or a
fixed positive imbalance. For every fixed C and every sufficiently large
order-N full signing with `Q(A)<=C N^(3/2)`, there is a same-order full
signing A' such that

```math
Q(A')=Q(A)+O_C(N^{5/4}),
\qquad |U(A')-V(A')|=O_C(N^{5/4}),
\beta(A')=\beta(A)+2|U(A)-V(A)|+O_C(N^{5/4}).    \tag{25}
```

After negating A if needed, take `Delta=U-V>=0` and
`m=floor(sqrt(2 Delta))`. Then `m<=sqrt(2C)N^(3/4)`. If m<8,
leave A unchanged; its imbalance is bounded by32 and (25) is immediate.
Otherwise use the one-block construction and put

```math
B_m=32m^2\beta(A)/N^2+64m^{3/2},\qquad
\theta_m=1-\frac{(N-m)_2}{(N)_2}.
```

These are explicit upper bounds for the mean old-block beta and the
relative loss of each fixed complement witness. For m>=8, the Bernoulli
reservoir comparison used for (12) has success probability at least
`1-exp(-m/4)>1/2`. Four Markov constraints with factor five give a block
for which the following **finite** inequalities all hold:

```math
U+\lfloor m/2\rfloor-5\theta_m U
\le U(A')\le U+\lfloor m/2\rfloor+\tfrac52 B_m,
V+\tfrac12m(m-1)-5\theta_m V
\le V(A')\le V+\tfrac12m(m-1)+\tfrac52 B_m,
\beta(A)+m(m-1)-5\theta_m\beta(A)
\le\beta(A')\le\beta(A)+m(m-1)+5B_m.            \tag{26}
```

Now `B_m=O_C(N^(9/8))`,
`theta_m beta(A)=O_C(m sqrt(N))=O_C(N^(5/4))`, and
`|m(m-1)/2-Delta|=O(m+1)`. These estimates are uniform for the entire
permitted range of m, proving (25). If A was negated initially, negate
the final matrix back; Q, beta and the absolute endpoint imbalance are
unchanged.

In particular exact minimizing inputs admit endpoint-balanced
near-minimizing outputs with an explicit `O_C(N^(5/4))` defect and the
displayed beta bookkeeping. The existence of balanced near-minimizers
was already known from the stronger archived repairs. This version
neither preserves exact minimality nor proves a favorable twist exists
for the prepared output.

For scope, a future estimate on this **prepared near-minimizer class**,
`F(A')<=2sqrt(2)Q(A')+O(N^(3/2-delta))`, would imply

```math
M_{2N}\le2\sqrt2 M_N
 +O(N^{3/2-\min(\delta,1/4)}).
```

Such an estimate has not been proved, and an estimate restricted only
to exact minimizers would not apply to A'. The separate all-order
coverage issue would still remain.

### Finite joint-invariant checks

`computations/twisted_chiral_endpoint_balance_2026_09_19.py` recomputes
U,V,beta and `beta+2|U-V|` for the archived small-order witnesses and
the complete order-nine and order-ten minimizing censuses. It also
exhaustively checks the exact complement-deficit expectations and the
finite one-clique inequalities on both order-ten classes. All checks
pass. Both order-ten classes and all order-nine minimizing classes are
endpoint-balanced. Of these census inputs, only the previously known
order-ten beta40 class exceeds `2sqrt(2)Q` before the O(N) matching
payment. This does not supply a growing seed family satisfying (24).

## 11. A bounded-operator counterexample, still at high cap

The arbitrary-seed obstruction is not solely a consequence of an
unbounded normalized operator norm. Here is a compact explicit variant;
no attempt is made to optimize its cap constant.

Let `N=4^k`, k>=3, and H be the Sylvester Hadamard. Write
`H=H4 tensor H_(N/4)`. In its first factor choose
`y0=(1,-1,1,1)` and `x0=(1,1,-1,1)`, so `H4 y0=2x0` and
`x0 dot y0=0`. Tensor both with the usual Boolean positive eigenvector
of the remaining even-dimensional Sylvester factor. This gives Boolean
x,y with `Hy=sqrt(N)x` and balanced product `t=x*y`. Set
`A0=H-diag(H)`; then `Q(A0)<=N^(3/2)/2` and
`x^T A0 y=N^(3/2)` (the diagonal contribution vanishes by the remaining
factor's zero trace).

Partition each of the four first-factor fibers into coordinate blocks
of size `m=2sqrt(N)`, varying a fixed set of lower binary coordinates.
Every old block is a signed Sylvester Hadamard of order m, minus its
diagonal, hence has operator norm at most `sqrt(m)+1`. On each block
replace the old edges by `t_i x_i x_j`; t is constant on a block. Half
the vertices receive positive switched cliques and half negative ones,
with the same block sizes. Call the resulting full signing A.

The new block-diagonal clique matrix has cap exactly `Nm/4`, and
bilinear evaluation at the old x,y equal to `N(m-1)`. Subtracting all
old blocks costs at most `N(sqrt(m)+1)/2` in Q and
`N(sqrt(m)+1)` in that bilinear evaluation. Therefore

```math
Q(A)\le(1+o(1))N^{3/2},\qquad
\beta(A)\ge(3-o(1))N^{3/2},\qquad
\|A\|_{op}\le(3+o(1))\sqrt N.                  \tag{27}
```

Every twist consequently exceeds `2sqrt(2)Q(A)` by a positive leading
amount. This is a high-cap counterexample, not one near the optimizing
0.43--0.50 range. It says that bounded normalized operator norm by
itself cannot be substituted for near-optimality in a uniform claim.

## 12. Exact behavior under local pair quarter-turns

Order the parent coordinates as n pairs, whose cross-pair block is
`M(a,b)=[[a,b],[b,-a]]`. Let
`J=[[0,-1],[1,0]]`. Local orientation-preserving rotations that retain
the zero diagonal of a nonzero matching edge must be quarter-turns
`J^(k_i)`. Even powers are ordinary paired switches. Write
`r_i=k_i mod2`, and ignore these paired switches when discussing
triangle signs.

Since `J^(-r_i) M(a,b) J^(r_j)=M(a,b)J^(r_i+r_j)`, the new hollow
blocks satisfy, for `q_ij=a_ij b_ij`,

```math
a'_{ij}=a_{ij}(-1)^{r_i r_j}q_{ij}^{r_i\mathbin\oplus r_j},
\qquad
q'_{ij}=a'_{ij}b'_{ij}=(-1)^{r_i+r_j}q_{ij},
\qquad d'_i=(-1)^{r_i}d_i.                     \tag{28}
```

For any signing X write `tau_X(i,j,k)=X_ij X_jk X_ki`. Equation (28)
shows that the **relative triangle array**

```math
\tau_A\tau_B=\tau_q
```

is invariant under all local quarter-turns. More explicitly, put

```math
f_r(i,j,k)=(-1)^{r_ir_j+r_jr_k+r_kr_i}
 q_{ij}^{r_i\oplus r_j}q_{jk}^{r_j\oplus r_k}q_{ki}^{r_k\oplus r_i}.
```

Then both individual triangle arrays are multiplied by the same factor:
`tau_A'=tau_A f_r` and `tau_B'=tau_B f_r`.

Triangle arrays characterize switching classes exactly. Thus the new
bridge is a signed-permutation copy of the new child **if and only if**
there exists a permutation pi such that

```math
\tau_{A'}(\pi i,\pi j,\pi k)
 =[\tau_A\tau_B](i,j,k)\,\tau_{A'}(i,j,k)
\quad\text{for every triple}.                 \tag{29}
```

Once pi is known, its switching signs are recovered from the root row;
there is no need to search all `2^(n-1)` switches. If a particular
original permutation p witnessed `tau_B=tau_A composed with p`, that
**same** p still works exactly when `f_r=f_r composed with p`, an
O(n^3) check. Other permutations may work even if this check fails.

In the important switch-only subclass `B=SAS`, every relative triangle
is +1. Equation (28) then proves that **all** local quarter-turns
preserve the switch-only-copy relation. This does not hold for the
general permutation-copy family.

A smallest transparent counterexample has n=4, with upper edges
ordered `(01,02,03,12,13,23)`:

```text
A = (+,+,+,-,-,+),
B = A with vertices 1 and 2 interchanged.
```

Rotate only pair2. The resulting A' has four negative triangles and
B' has two. They cannot be switching-permutation equivalent. The full
parent remains a valid chiral signing, with unchanged cap, so this
example precisely separates the broader chiral gauge symmetry from
the restricted same-child twist relation.

The exact formulas, original-permutation survival criterion, and
switch-only closure were exhaustively checked on 24576 order-four
triples (A, signed permutation, rotation mask) by
`computations/twisted_chiral_pair_rotation_audit_2026_09_19.py`.
Its retained counterexample also reconstructs the actual signed
quarter-turn and checks all 256 parent-spin energies before and after.

## 13. Focused audit: orthogonal opposite grounds and the lower response

This angle was checked on the resumed campaign's thirty-minute budget.
It yields useful exact geometry but **no improved asymptotic lower
constant** from the presently proved marked/unmarked response theorem.
The sources read for this check are
`decisive_audit_fresh_full_lower_chain_2026_09_07.md`, especially Sections2--5,
`resumed_response_center_theorem_independent_audit_2026_09_06.md`, and
`paper_discrepancy_low_cap_response_audit_2026_09_17.md`.

### Exact additional geometry

If `J^2=-I` is a signed permutation and `J^T D J=-D`, then
`J^T=-J`. Thus a positive ground z has an exactly orthogonal negative
ground `w=Jz`: `z dot w=0` and `H_D(w)=-Q(D)`. The product `z*w`
is balanced, with one + and one - in each J-pair.

Let I be its + coordinates and J0 its - coordinates. Writing the two
energies as the sum of the two internal energies plus/minus the cut
energy proves exactly

```math
H_{D[I]}(z_I)+H_{D[J0]}(z_{J0})=0,
\qquad z_I^T D[I,J0]z_{J0}=Q(D).
```

For any principal rectangle the bilinear norm is at most Q(D), by
flipping the entire spin word on one side. Therefore this balanced
rectangle **saturates** that bound:

```math
\beta(D[I,J0])=Q(D).                            \tag{30}
```

In the all-positive ground gauge, local stabilities are nonnegative
and the opposite-ground stabilities are exactly their J-permutation.
These facts do not impose a distributional limit on the local fields.

### Why the existing marked calculation does not gain a factor

The archived mark is the independent **own-coordinate input spin** of
an injective tree field, not a selected ground state. Replacing the
product Rademacher input law by the two ground words destroys the
own-spin independence, low influences, and Gaussian tree limit needed
for the marked/unmarked identities. Orthogonality of two deterministic
words is not a substitute for those hypotheses.

Keeping iid inputs is valid, but coupling S with JS does not strengthen
the final energy comparison. For an old marked-tree field X, all its
vertices have odd degree, so the tree has an odd number of edges;
signed-permutation equivariance and anticommutation give
`X(JS)=-J X(S)`. For the existing odd F, even H response class this
implies

```math
\mu_+(JS)=J\mu_-(S),\qquad
\mu_\pm=\pm F+H\,\operatorname{sign}(BF).
```

Since JS has the same product law as S, the two expected energies are
opposite. But the proof already bounds their full difference by 2Q.
It therefore gives exactly its previous scalar reward
`E H E|K_F+tau Z|`, with no additional factor or second independent
energy payment. Treating the opposite ground as another additive
reward would double-count the same energy interval.

### The variance step is already sharp inside the chiral class

The normalized square `C=B^2` commutes with J, so paired off-diagonal
entries of C vanish. Corresponding old Gaussian fields at paired roots
are uncorrelated. Nevertheless the response theorem's odd-Schur
variance lower bound cannot acquire a uniform strict gain merely from
this structure: `C=I` saturates it and is compatible with chirality.
There are actual hollow full-sign chiral sequences approaching this
case. Let H be a Sylvester symmetric Hadamard, chiral in its first
binary tensor factor, and put `D=H-diag(H)`. Then, with
`B=D/sqrt(N-1)`,

```math
\|B^2-I\|_{op}\le\frac{2+2\sqrt N}{N-1}\longrightarrow0.
```

Both H and its removed diagonal anticommute with the same quarter-turn,
so D remains chiral. This example does not disprove a stronger chiral
cap lower theorem; it identifies why a strict gain at this particular
variance/Jensen step is unavailable.

The spectral-core deletion can preserve chirality without changing its
leading estimates: average the existing diagonal majorant with its
J-conjugate. It remains a majorant for both signs, has the same trace,
and is constant on J-pairs, so threshold deletion removes whole pairs.
This retains the hypothesis but creates no numerical improvement.

The missing step is a **new** feasible response or a new uniform
constraint on the joint ground-field law, not the mere existence of
the orthogonal opposite ground. No half-constant lower bound, no gain
over 0.433322..., and no contradiction to the finite cap40 chiral
witness follow from this audited angle. The angle is stopped here.

## 14. Closing independent identity and selector audits

The director's `twisted_chiral_director_native_gap_2026_09_19.md`, Section1,
was independently reconstructed. For a general chiral matrix
`D=[[A,C],[C,-A]]`, rotate pairs in J by
`(x_i,y_i)=(y'_i,-x'_i)`. Its new bridge, in the order I=J^c,J, is

```math
C_J=\begin{pmatrix}C[I,I]&A[I,J]\\A[J,I]&-C[J,J]\end{pmatrix}.
```

Every such bridge has beta at most Q(D). Conversely, at a parent
extremizer choose `J={i:x_i=-y_i}`. Both transformed halves equal the
old x, so the internal energies cancel and the bridge quadratic value
attains Q(D). Thus exactly `Q(D)=max_J beta(C_J)`. In the positive
all-ones ground gauge, the two single-spin stabilities are `c_i+a_i`
and `c_i-a_i`, where c and a are the bridge and internal row sums.
Both are nonnegative, proving `c_i>=|a_i|`. The stronger all-spin
inequality stated there is equivalent to the original cap assertion,
not an independently smaller sufficient state. Audit: **PASS**.

The uniform track's proposed beta-ground-product switching rule also
has a genuine quadratic-scale counterexample. Put

```math
A_m=\begin{pmatrix}J_m-I_m&J_m\\J_m&-J_m+I_m\end{pmatrix},\qquad m\ge2.
```

For a spin y with block sums u,v, its response norm is
`sum_first |u+v-y_i|+sum_second |u-v+y_i|`. If `|u|>|v|`, this is
`(2m-1)|u|+v sign(u)<=2m^2-2`; the strict magnitude difference is at
least two by parity. If `|v|>|u|`, it is
`(2m-1)|v|-u sign(v)<=2m^2-2`. At equal magnitudes it is
`(2m-1)|u|+m`, attaining `2m^2` only at macro corners. All fields at
those corners are nonzero, and their unique positive bilinear response
is x=+y or x=-y. Consequently

```math
\beta(A_m)=2m^2,
\qquad \operatorname{diag}(x_i y_i)=\pm I
\quad\hbox{for every bilinear extremizing pair.}
```

That rule therefore selects only the identity bridge. Its core parent
cap is exactly `4m^2-2m`: x=(+,+), y=(+,-) attains it. On the
two-dimensional macroconstant subspace the child eigenvalue square is
`2m^2-2m+1`, while the other eigenvalues are +/-1. Hence the spectral
parent upper bound is `2m sqrt(4m^2-4m+2)`, strictly below
`4m^2-2m+1`; the squared gap is `(2m-1)^2`. Integrality completes the
upper certificate. A matching changes the cap by at most 2m.

Thus **every beta-ground-product selector** fails the proposed
`4 beta/3+O(n)` target on this family. This is not a counterexample to
optimized switches or permutations. Audit: **PASS**, with that scope.

The uniform coarse theorem's subsequent diagonal correction also passes:
the scalar diamond maximum of `|a^2/2-c^2/2+bd|` is exactly1, so
`Delta(diag t)=sum|t_i|`. Hollowing a block average P therefore gives
`Delta(P)<=4 beta(P)/3+(7/3)sum|P_ii|`. Together with the previous
matching n cost, the final finite coarse remainder is
`(40/3)epsilon n^2+4kn+(10/3)n`, not the earlier last term n. Its
`O(n^2/sqrt(log n))` optimized scope is unchanged; it is not native
subleading-scale purification.

## 15. Bounded native counterattack: what does and does not amplify

The stronger native question

```math
F(A)\le\frac43\beta(A)+o(n^{3/2})
```

remains unresolved. The finite order-six conference child has beta=12
and minimum core cap18, exceeding16. That excludes zero-error native
purification, not an O(n) defect. The following checks delimit concrete
attempts to amplify that failure.

### Fine-fiber cloning removes, rather than preserves, the seed gap

The exact integer diamond profile of the order-six seed gives
`Delta(C_6)=12=beta(C_6)`. For even m, the weighted clone
`C_6 tensor J_m` admits a switch that is balanced inside each of its
six fibers. Its parent core cap is **exactly** `m^2 Delta(C_6)=12m^2`,
not the extrapolated `18m^2`. Indeed normalized sums and switched sums
of the two parent spin words give precisely the two coordinate
diamonds; their vertices are all realized by constant spins on the
two half-fibers. Hollow multilinearity supplies a maximizing vertex.

For example take m along powers of four and fill the six zero diagonal
macroblocks with hollow symmetric Hadamard signings. The block-diagonal
filling E has `beta(E)=O(m^(3/2))`; its contribution to the selected
parent is at most `2Q(E)+beta(E)`, and the matching costs only6m.
The resulting actual full-sign child therefore has

```math
\beta=12m^2+O(m^{3/2}),\qquad
F\le12m^2+O(m^{3/2}).
```

These are high-cap examples, but they decisively invalidate this
particular extrapolation of the microscopic factor18/12. Disjoint sums
are no substitute: even an additive fixed gadget gap gives only O(n),
and their missing cross edges cannot be silently ignored.

### Exact spin witnesses exclude two natural low-beta tensor families

Tensor beta is not multiplicative. Already the two-by-two Sylvester
Hadamard satisfies `beta(H_2)=2` but `beta(H_2 tensor H_2)=8`, not4.
The fine spin words are not restricted to product words, and the native
twist permutations need not preserve tensor fibers either.

There is a quantitative certificate for the most immediate
conference-six-derived Hadamard candidates. Define

```math
H_{12}=\begin{pmatrix}C_6-I&C_6+I\\C_6+I&-C_6+I\end{pmatrix}.
```

It is symmetric, trace zero, and `H_12^2=12I`. Exact enumeration gives
`beta(H_12)=36`, normalized `sqrt(3)/2=0.866025...`. However the retained
explicit integer spin pairs prove the following **lower bounds only**:

| matrix | order | certified bilinear value | normalized value |
| --- | ---: | ---: | ---: |
| H12 tensor H2 | 24 | 112 | 0.9525793444 |
| H12 tensor H4 | 48 | 324 | 0.9742785793 |
| H12 tensor H12 tensor H12 | 1728 | 66592 | 0.9270570295 |

The certificate is
`computations/twisted_chiral_tensor_beta_witness_2026_09_19.py`, with
spin words in the matching results JSON. The small witnesses are
reconstructed by plain scalar integer sums; the tensor-cube witness is
independently reconstructed by dense integer multiplication, not only
the tensor response used to find it. No optimum is claimed at these
three larger orders.

For any symmetric Hadamard H of order h, `vec(H)` is a sign eigenvector
of `H tensor H` with eigenvalue h. Thus every even power of H12 has
normalized beta exactly1. Tensoring the retained cube witness with
copies of this flat square witness gives normalized beta at least
0.927057 for **every odd power at least three**. Likewise H4 has the
flat sign eigenvector `(1,1,1,-1)`, so the first two witnesses extend
to every odd and positive even Sylvester exponent, respectively, in
`H12 tensor H_(2^k)`.

Finally hollowing changes beta by at most the order, because the
removed diagonal has entries +/-1. Consequently neither growing
family can have asymptotic normalized beta below the approximately
0.919226 threshold that, via the independently archived parent lower
bound, would falsify native `4beta/3+o(n^(3/2))` for every twist. This is
an exclusion of **these candidate low-beta counterattacks**, not a
proof of the native inequality and not an exclusion of all Hadamard or
conference constructions. The finite-to-asymptotic native gap remains.

## 16. Actual insertion preserving a chosen same-child twisted parent

This is an insertion operation on the actual signing and its chosen
parent, not a reformulation of their joint norm. No cap hypothesis on
the input is needed.

Let A be any order-n full signing, n>=1, and fix any allowed parent

```math
D=\begin{pmatrix}A&G^T A G+\operatorname{diag}(d)\\
G^T A G+\operatorname{diag}(d)&-A\end{pmatrix},
```

where G is a signed permutation. For any integer r>=1 set

```math
N=n+r,\qquad L=nr+\binom r2,\qquad V=8nr+4\binom r2.
```

There is a full-sign extension A' on N vertices and an allowed parent
D' of A' such that A and D are their respective induced submatrices,
and **simultaneously**

```math
Q(A')\le Q(A)+\sqrt{2L(N+2)\log2},
```

```math
\beta(A')\le\beta(A)+\sqrt{8L(2N+1)\log2},
```

```math
Q(D')\le Q(D)+r+\sqrt{2V(2N+1)\log2}.             \tag{35}
```

In particular, starting with a minimizing allowed parent gives

```math
F(A')\le F(A)+r+\sqrt{2V(2N+1)\log2}.
```

The added matching signs can be prescribed arbitrarily. The r=0 case
is the unchanged input. For r=o(n), all displayed errors are
`O((n+r)sqrt(r))=o(n^(3/2))`. This does not bridge a fixed proportional
order gap or provide the missing native doubling bound.

### Construction and exact coefficient variance

Independently choose the L child edges incident to the r new vertices
as fair signs. Write `A'=A direct-sum 0_r+E` and extend the signed
permutation by `G'=G direct-sum I_r`. Retain d on old vertices and use
the prescribed new matching signs. The old child and old parent are
then induced submatrices **exactly**, and the bridge remains the
signed-permutation copy `G'^T A' G'`. Parent entries are correlated;
the independent variables used below are the new **child** edges.

Temporarily omit the r new matching edges. At a fixed parent spin
(x,y), put u=Gx and v=Gy on the old coordinates. The coefficient of
the random old-new child edge (i,j) is

```math
c_{ij}=x_j(x_i+v_i)+y_j(u_i-y_i).
```

For fixed new spins `a=x_j,b=y_j`, its entire old-coordinate vector is

```math
c_{\cdot j}=(aI+bG)x+(aG-bI)y.
```

The concatenated matrix T=[aI+bG, aG-bI] satisfies `TT^T=4I` because
G is orthogonal and a,b are signs. Therefore

```math
\sum_{i=1}^n c_{ij}^2\le4(\|x\|_2^2+\|y\|_2^2)=8n.
```

For a new-new child edge the coefficient is
`x_i x_j-y_i y_j+x_i y_j+x_j y_i`, always +2 or -2. Hence the total
squared coefficient sum is at most V. This improves the bound16L
obtained by separately bounding every coefficient by4.

The constant8 in the old-new square sum is sharp for this uniform
coefficient method: take n even, G a direct sum of signed quarter-turns,
y=-Gx, and all new x_j=y_j=1. Then
`c_(.j)=2(I+G)x` and its squared norm is8n. The new-new coefficients
still have squared magnitude4, so V itself is attained.

### One completion controls all three norms

A Rademacher sum with coefficient square sum at most s has one-sided
tail at most `exp(-t^2/(2s))`. For the child quadratic increment the
square sum is L. Modulo overall sign there are `2^(N-1)` spin words;
using both tails gives

```math
\Pr\{Q(E)>t\}\le2^N\exp[-t^2/(2L)].
```

For the child bilinear increment each coefficient is
`x_i y_j+x_j y_i`, of magnitude at most2, so the square sum is at most
4L. Its absolute maximum equals its positive maximum, and simultaneous
overall reversal of x,y leaves the value unchanged. Thus there are at
most `2^(2N-1)` one-sided tests and

```math
\Pr\{\beta(E)>t\}\le2^{2N-1}\exp[-t^2/(8L)].
```

The matching-free parent increment is chiral. Its negative and positive
extrema therefore agree in absolute value, and overall spin reversal
also preserves its energy. Again there are at most `2^(2N-1)`
one-sided tests, now with coefficient square sum V:

```math
\Pr\{Q(\text{parent core increment})>t\}
\le2^{2N-1}\exp[-t^2/(2V)].
```

At the three respective thresholds in (35), before adding the matching
r, each failure probability is at most1/4. The union bound leaves
success probability at least1/4 for the **same** random child extension.
Triangle inequalities and the deterministic matching cost r prove all
claims. This also gives a finite randomized realization with exact norm
verification if desired; no polynomial-time verification claim is made.

Restriction monotonicity additionally gives `Q(A')>=Q(A)`,
`beta(A')>=beta(A)`, and `Q(D')>=Q(D)` for this chosen parent. It does
not give `F(A')>=F(A)`, since another optimizing extended permutation
may mix old and new vertices.

### Archive attribution and exact scope

The concentration mechanism is elementary and already archived.
`minimal_all_order_action_recovery.md`, Section7.2, equation(AR.13),
proves random completion of any signing with `O(N sqrt(r))` cap cost.
`flatify_director_variance_budget_cap_comparison_2026_09_07.md`, Section2,
and its independent audit give the broader simultaneous variance-budget
insertion principle. `cut_triangle_midpoint_balancing.md`, Section7,
already controls a random chiral matrix by a Rademacher coefficient
square sum and a Boolean union bound.

Those statements do not impose the coupled same-child bridge relation
or retain a prescribed old parent: ordinary independent parent-edge
completion would break that relation. A targeted search of the chiral,
dependent-lift, padding, and insertion artifacts did not find this exact
`G direct-sum I_r` extension with its coefficient variance V. The scoped
addition here is its explicit constraint-preserving realization and
simultaneous Q/beta/parent accounting, not a new concentration principle
or an externally established novelty claim. Together with the audited
cycle-repair deletion theorem it supplies two-sided **near-order**
operations inside the same-child twisted family, not fixed-ratio
scale transfer.

### Corollary: a relative-order modulus for the optimized family

Define the finite-family optimum

```math
T_n=\min_{A\text{ an order-}n\text{ full signing}}F(A).
```

Then there is an absolute constant C such that, for every integer
`1<=r<=n/4`,

```math
|T_{n+r}-T_n|\le Cn\sqrt r.                         \tag{36}
```

First this family has a uniform crude cap bound. Choose the child
edges independently as fair signs and use the identity twist. Each
child edge has parent-core coefficient exactly +/-2 at every parent
spin, so the variance proxy is `4 binom(n,2)`. The same one-sided
projective union bound gives

```math
T_n\le n+\sqrt{4n(n-1)(2n+1)\log2}
\le K_0 n^{3/2},\qquad K_0=1+\sqrt{8\log2}.         \tag{37}
```

This bound is only used to put optimizing seeds in a bounded-cap
class; it is not a competitive original minimax upper bound.

For every seed let U=max H_A and V=-min H_A. Reversing a whole parent
half proves `Q(parent)>=U+V>=Q(A)`, hence `Q(A)<=F(A)` and
`beta(A)<=4Q(A)<=4F(A)`. An optimizing seed for T_n therefore satisfies
the bounded-cap hypothesis uniformly.

Applying (35) to such a seed gives
`T_(n+r)<=T_n+O(n sqrt(r))`. For the reverse direction choose an
optimizing seed at N=n+r. The independently audited deletion theorem
in `twisted_chiral_near_order_2026_09_19.md` supplies an n-vertex
restriction C with

```math
F(A[C])\le T_N+
300\left[\frac rN\beta(A)+N\sqrt r\right]+r.
```

Here r<=n/4 implies r<=N/4, so that theorem applies; also
`beta(A)<=4K_0 N^(3/2)` by (37). Since `r sqrt(N)<=N sqrt(r)` and
`N<=5n/4`, its error is `O(n sqrt(r))`. As `T_n<=F(A[C])`, this is
the other half of (36).

Equivalently, `T_n/n^(3/2)` has an absolute relative-order modulus
`O(sqrt(r/n))` in this range. Such a modulus allows bounded slow
oscillation and does **not** imply existence of a limit, a favorable
doubling inequality, or passage across a proportional gap at vanishing
cost. The proof combines the archived elementary random-padding
mechanism with the exact orbit-preserving insertion and audited cycle
repair; it does not assume the chosen child minimizes Q(A).

## 17. Audit of the local-stability-filtered random-twist first moment

The director's closing conditional-permanent derivation passes. This
section records its exact normalization and a degree-only upper bound;
it does not assert that the resulting first moment is asymptotically
small.

Fix parent spins x,y and write `t_i=x_i y_i`, with k plus coordinates.
Choose a uniformly random signed permutation specified by
`g e_i=s_i e_(pi(i))`, independently of uniform matching signs d, and
put z=gx,w=gy. Then v=z*w is uniform among masks with k plus
coordinates, and z is an independent uniform spin word. Conditional
on z,w, pi is uniform over the `k!(n-k)!` bijections taking each
t-sector to its corresponding v-sector. All signs are forced by
`s_i=z_(pi(i)) x_i`.

Define

```math
\alpha_i=x_i(Ax)_i,\quad \gamma_i=y_i(Ay)_i,\quad
b_j=z_j(Aw)_j,\quad c_j=w_j(Az)_j,
\qquad u_i=d_i t_i.
```

The two signed local fields at the parent coordinates paired with i
are exactly

```math
\alpha_i+b_{\pi(i)}+u_i,\qquad
-\gamma_i+c_{\pi(i)}+u_i.
```

The u_i remain independent fair signs after this conditioning. The
parent energy is

```math
E_0+\sum_i u_i,\qquad
E_0=H_A(x)-H_A(y)+z^T A w,
```

and E0 is independent of the remaining sector bijection pi. Set
`m_ij=min(alpha_i+b_j,-gamma_i+c_j)`. It is even: each of its two
summands before taking the minimum is a sum or difference of two
integers congruent to n-1 modulo2. A parent has odd degree2n-1, so
single-spin local maximality is strict positivity of these fields.
Thus u=+1 is allowed exactly when m>=0, and u=-1 exactly when m>=2.

Use a scalar indeterminate q, distinct from the product mask t, and
form the two sector matrices with entries

```math
K_{ij}(q)=\frac12\left(q\,1_{m_{ij}\ge0}
                         +q^{-1}1_{m_{ij}\ge2}\right).
```

The exact conditional Laurent generating polynomial for positive
local stability, weighted by `q^(sum u_i)`, is

```math
P(q)=\frac{\operatorname{per}K_+(q)\,
              \operatorname{per}K_-(q)}{k!(n-k)!}.          \tag{38}
```

Empty-sector permanents and 0! equal1. Consequently the conditional
probability of a locally stable parent with energy greater than L is
the sum of coefficients `[q^h]P(q)` over `h>L-E0`. For any lambda>=0
it is at most `exp[-lambda(L-E0)]P(exp(lambda))`. The strict-tail
condition and any lattice improvement should be retained in exact
computations.

### A coefficientwise Brégman bound, not just a Chernoff estimate

Within each row's own sector, let d_i^0 and d_i^2 be its numbers of
columns satisfying m>=0 and m>=2. Define `phi(0)=0` and
`phi(a)=(a!)^(1/a)` for positive integers a. The Brégman permanent
inequality is `per(M)<=product_i phi(rowdegree_i)` for a zero-one
matrix. A primary proof is
[Radhakrishnan, An Entropy Proof of Bregman's Theorem](https://doi.org/10.1006/jcta.1996.2727).

Expand a sector permanent row-by-row, selecting either its
`q 1_(m>=0)/2` row or its `q^(-1)1_(m>=2)/2` row. Apply Brégman's
inequality separately to each resulting zero-one permanent. All
coefficients are nonnegative, so this gives the **coefficientwise**
Laurent inequality

```math
\operatorname{per}K_\pm(q)
\preceq \prod_{i\text{ in that sector}}
\frac{\phi(d_i^0)q+\phi(d_i^2)q^{-1}}2.                \tag{39}
```

Multiplying the sector bounds and dividing by `k!(n-k)!` preserves
coefficientwise domination. Therefore the appropriate coefficient tail
of this simple n-factor product bounds (38)'s high-energy tail
directly; using Chernoff is optional. This still needs the actual
conditional field-degree data, not only Q or beta of the seed.

The director's alternative auxiliary-edge argument also passes:
`K(e^lambda)/cosh(lambda)` has entries1 on m>=2, entries
`p=e^lambda/(2cosh(lambda))` on m=0, and zero otherwise. Independent
Bernoulli(p) weak edges preserve the expected permanent. Applying
Brégman and independence of row degrees gives

```math
\operatorname{per}K_\pm(e^\lambda)
\le\cosh(\lambda)^{|\text{sector}|}
\prod_i \mathbb E\phi\bigl(d_i^2+
                 \operatorname{Bin}(d_i^0-d_i^2,p)\bigr).
```

One may instead correlate all weak auxiliary edges within each row,
independently between rows. A permanent monomial uses one entry per
row, so its expectation is still correct. This gives (39) evaluated
at e^lambda with only two phi values per row, avoiding binomial sums.
No unproved comparison between these two bounds is required.

For completeness, the standard entropy proof of Brégman is short.
Choose a uniform allowed perfect matching and reveal its rows in an
independent uniform order. The entropy contribution of a row is at
most the expected logarithm of its number of still-available allowed
columns. Conditional on the matching, the row's position among the
rows matched to its r allowed columns is uniform, so that expectation
is `(log(r!))/r`. Summing the entropy chain rule and exponentiating
proves the inequality. A zero row makes the permanent zero.

### First-moment cutoff and scope

Summing (38)'s conditional tails over x,y and averaging z,v with weight
`1/(2^n binom(n,k))` gives the exact expected number of positive
locally stable violations for a random allowed twist and matching.
If all spins x,y are counted, a violation has at least its two global
sign copies. Thus expectation less than2 suffices to exhibit cap<=L.
Equivalently restrict to projective representatives x_1=+1 and use
cutoff1. The old eightfold bad-pair multiplicity cannot be reused:
the chiral quarter-turn sends a positive local maximum to a negative
local minimum, not another positive locally stable violation.

This is a valid filtered sufficient certificate with an exact
conditional identity and an optional explicit row-degree relaxation.
It neither proves a uniform first-moment bound nor turns failure of
the certificate into failure of an actual matrix or twist.

The completed stability wind-tunnel source was also independently
audited: `twisted_chiral_stability_histogram_2026_09_19.cpp` and its
Python driver preserve the Gray energy and field recurrences exactly,
check strict stability in **all** coordinates including the fixed
projective coordinate, and use exact integer strict-threshold tests.
The matching is present, so all fields are odd. Known minimizing
parents are excluded from random-sample aggregates. The result contains
151 exact per-matrix histograms and six independent direct NumPy
replays; the reported sample averages remain Monte Carlo diagnostics,
not population expectations or existence certificates. Source audit:
**PASS**.

## 18. Direct tensor-Gram proof simplification and completion audit

The complete fixed-seed even-power classification is already proved in
`transfer_seed_tensor_power_instability_2026_09_06.md`. The director's
closing observation gives a substantially shorter direct Boolean proof;
it is preserved as a **proof simplification**, not a newly discovered
classification or an unrestricted tensor no-go theorem.

Let H be a symmetric full sign matrix, including its diagonal, of
fixed order q>=2. Under either consistent vectorization convention,

```math
(H\otimes H)\operatorname{vec}(H)=\operatorname{vec}(H^3),
\qquad
\operatorname{vec}(H)^T(H\otimes H)\operatorname{vec}(H)
=\operatorname{tr}(H^4).                              \tag{40}
```

The vector is itself Boolean. Its k-fold tensor therefore has full
quadratic value `(tr H^4)^k` on `H^(tensor 2k)`. Hollowing subtracts
the spin-independent trace `(tr H)^(2k)`. Set
`rho=tr(H^4)/q^3` and `tau=(tr H)^2/q^3`. This gives

```math
\frac{Q(\operatorname{hollow}(H^{\otimes2k}))}{q^{3k}}
\ge\frac12(\rho^k-\tau^k).                          \tag{41}
```

Since `tr H^2=q^2`, Cauchy--Schwarz on the squared eigenvalues gives
`tr H^4>=q^3`, with equality exactly when `H^2=qI`. Also
`tau<=1/q`. A fixed non-Hadamard seed thus has rho>1 and divergent
normalized even-power caps. In the Hadamard case the spectral upper
is `(1+tau^k)/2`, matching (41) asymptotically at one half.

For heterogeneous paired tensors `tensor_j(H_j tensor H_j)` of total
order `N=product_j q_j^2`, the same explicit product spin gives

```math
\frac{Q(\operatorname{hollow}(\bigotimes_j(H_j\otimes H_j)))}{N^{3/2}}
\ge\frac12\prod_j\frac{\operatorname{tr}(H_j^4)}{q_j^3}
       -\frac1{2\sqrt N}\ge\frac12-\frac1{2\sqrt N}. \tag{42}
```

Thus these literal paired families cannot be asymptotically sub-half,
even when their seeds vary. No statement about arbitrary external
Hadamard padding, unpaired heterogeneous tensors, optimized subsequent
restrictions, or non-tensor constructions is being inferred.

Non-Hadamard **odd** powers also diverge by appending one fixed spin
v with `w=v^T H v!=0` to the paired witness. For power2k+1 the lower
bound becomes

```math
\frac12\left[
\frac{|w|}{q^{3/2}}\rho^k-
\frac{|\operatorname{tr}H|}{q^{3/2}}\tau^k\right].
```

Such v is explicit: use the all-ones spin unless its value is zero.
In that case q is even, and flipping coordinate1 changes the value by
minus four times an odd off-diagonal row sum, hence makes it nonzero.
This proves non-Hadamard divergence along all powers, but does not
classify the odd-power Hadamard branch.

### Exact optimization over all sign diagonals

For a hollow order-q full signing D, put
`c_i=(D^3)_ii` and `Omega0=sum_(i!=j)((D^2)_ij)^2`. For any sign
diagonal h and `H=D+diag(h)`, the off-diagonal squared-Gram defect is

```math
\Omega(H)=\Omega_0+2q(q-2)+2\left(\sum_i h_i\right)^2
                     +4\sum_i c_i h_i.                  \tag{43}
```

Indeed `(H^2)_ij=(D^2)_ij+D_ij(h_i+h_j)` off the diagonal.
The cross term sums to `4 sum_i c_i h_i`, while the remaining square
sum is `2q(q-2)+2(sum h)^2`. Also `tr H^4=q^3+Omega(H)` exactly.

If exactly k entries of h are plus, the linear term is minimized by
choosing the k smallest c_i. Sorting c and sweeping its prefix sums
therefore finds the exact best diagonal in `O(q log q)` additional
time, **after** the matrix quantities c and Omega0 are computed.
No enumeration of all `2^q` diagonals is needed.

### Chiral seeds remove the diagonal correction entirely

If D anticommutes with a signed quarter-turn, the diagonal of D^3 has
paired opposite entries. Choose `h_i=-sign(c_i)` on nonzero pairs and
opposite signs on each zero pair. This simultaneously attains zero
sum and the smallest possible linear term. Thus

```math
\Omega_*:=\min_h\Omega(D+\operatorname{diag}h)
=\Omega_0+2q(q-2)-4\|c\|_1.                           \tag{44}
```

For **every** h, writing s=sum h gives
`Omega(H)>=Omega_*+2s^2`. Set `a=q^3+Omega_*`, t=s^2. The exact
even-tensor witness is at least
`(a+2t)^k-t^k>=a^k`, by binomial expansion. Therefore, for every
completion and every k>=1,

```math
\frac{Q(\operatorname{hollow}(H^{\otimes2k}))}{q^{3k}}
\ge\frac12\left(1+\frac{\Omega_*}{q^3}\right)^k.      \tag{45}
```

There is no diagonal correction in this chiral bound. A positive
defect for each fixed seed does not alone give a uniform positive
normalized defect for a growing sequence of seeds; that quantifier
distinction remains essential.

The source `twisted_chiral_tensor_gram_2026_09_19.py` was independently
audited. In addition, both saved order20 parents were reconstructed by
plain Python integer sums, without importing that verifier or NumPy.
The cubic pairing, balanced minimizing diagonal, Gram defect, and
full four-index quadratic expression in (40) all agree. The original
Hadamard-completable witness (hash beginning `d327680f`) has Omega*=0
and square witness4000. The new witness (hash beginning `247d887a`)
has Omega*=960; **every** sign-diagonal completion therefore has
order400 tensor-square cap at least4480, normalized0.56. These are
actual Boolean tensor witnesses, not only spectral certificates, and
they do not assert the exact cap of the tensor-square matrices.

## 19. Final exact-average audit and sharp projective multiplicity

The complete small-order average in
`computations/twisted_chiral_stable_average_2026_09_19.py` passes source
audit and arithmetic replay. The sector permanent numerators count
sector bijections and matching words without dividing by their
probabilities; summing over image sectors and z then covers every
signed permutation and matching exactly once for each projective
parent spin. Hence its common denominator `4^n n!` is correct.
Factorial-root upper approximations are certified by integer powers,
and their n-factor products use denominator `SCALE^n` as required.

A separate scalar implementation,
`computations/twisted_chiral_stable_average_adversary_2026_09_19.py`,
does not use the conditional fields, permanents, or NumPy. It directly
forms all distinct bridges with their exact signed-permutation
multiplicities, then all matching words and all projective spins.
Every local field is recomputed by integer sums. Its complete raw and
stable energy histograms agree with the archived JSON. Separately,
the original permanent/Brégman calculation was rerun in memory and
its entire returned records agree with the saved data.

The exact outcomes are:

| child order | complete group/matching denominator | L at the known family minimum | expected projective stable tail | first L certified by either exact stable or Brégman tail |
| --- | ---: | ---: | ---: | ---: |
| 3 | 384 | 5 | 5/4 | 7 |
| 4 | 6144 | 10 | 4/3 | 12 |

At the first successful L, the exact stable tails are respectively
1/2 and13/32; the rational Brégman upper bounds are approximately
0.5000003162 and0.4251483286, both strictly below1. Thus the filtered
criterion fails to certify the already known best caps5 and10 on
these two examples. This is a limitation of this particular sufficient
first-moment certificate, not a failure of the optimizing matrices.

### Chirality cannot improve the universal cutoff above one

There is an explicit counterexample at every child order m>=1. Set
`A=J_m-I_m`, choose the identity bridge and all-positive matching, so

```math
D_m=\begin{pmatrix}J_m-I_m&J_m\\J_m&-J_m+I_m\end{pmatrix}.
```

For a parent spin (x,y), let u=sum x and v=sum y. Its signed local
fields are `x_i(u+v)-1` and `y_i(u-v)+1`. Strict positivity of the
first group forces x to be constant, say x=a1. Since u-v is even,
the second condition is `y_i(u-v)>=0`. If y had both signs, this
would force v=u=am, impossible for a mixed word. Thus y is constant;
the opposite choice y=-a1 has negative second-block fields, while
y=a1 has all fields positive. Therefore **exactly one** positive
stable state exists modulo global reversal, namely both blocks
constant with the same sign, of energy m^2.

In particular the actual order-three child used in the average has a
parent with just one projective stable violation above L=5 (and L=7).
The quarter-turn supplies a negative local minimum, not another
positive maximum. No universal projective positive-stability
multiplicity larger than1 follows from chirality.

There is also a tail-specific witness inside the precise order-four
optimal-child family used above. For

```math
A_4=\begin{pmatrix}0&1&1&1\\1&0&1&1\\1&1&0&-1\\1&1&-1&0\end{pmatrix},
```

take identity permutation, switch `(1,-1,-1,-1)`, and matching
`(-1,-1,-1,-1)`. Exact enumeration of128 projective parent spins gives
positive stable energies `[10,10,14]`. Consequently its stable tail
above L=10 contains **one** state, represented by
`(1,1,1,1; -1,1,-1,-1)`. A cutoff2 would therefore also be invalid
for that fixed child's L=10 test. The actual matrices, stable words,
and signed local fields are retained in the adversarial replay JSON.

## 20. Final audit: cheap dominance-support feasibility

The conditional weak support in Section17 is exactly the two-coordinate
dominance relation `b_j>=-alpha_i` and `c_j>=gamma_i`. Given these
thresholds and column coordinates, a perfect matching is decided and
constructed in `O(n log n)` time as follows: process rows in decreasing
order of `a_i=-alpha_i`, insert all columns with `b_j>=a_i`, and take
the available column of smallest c_j satisfying `c_j>=gamma_i`.

The exchange proof passes. In any feasible matching extending the
previous greedy choices, let j serve the current row and h be the
greedy choice. Then `c_h<=c_j`. If h instead serves a later row ell,
exchange these two columns. Column j is feasible for ell because
`b_j>=a_i>=a_ell` and `c_j>=c_h>=gamma_ell`. Earlier choices are
untouched. Failure to find h therefore certifies infeasibility.
Ties and duplicate coordinates cause no problem if column identities
and multiplicities are retained.

The strong support uses thresholds `(2-alpha_i,2+gamma_i)` and the
same proof. More generally any specified row-wise choice of weak or
strong threshold is handled after sorting its resulting first
thresholds. A balanced successor multiset or coordinate-compressed
Fenwick tree supplies the advertised running time; ordinary linear
deletion from a Python sorted list would not. This cost is for the
matching test **after** the fields are available.

This is an imported elementary dominance-matching mechanism applied
to the actual support, not a novel general matching algorithm. It
can eliminate conditional terms with a genuine Hall obstruction that
row degrees alone miss. It neither counts the remaining perfect
matchings nor estimates the averaged stable high-energy tail.
