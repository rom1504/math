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

### Uniform balancing preparation, including N-dependent imbalance

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
near-minimizing outputs with an explicit `O_C(N^(5/4))` defect. This is
a uniform positive preparation theorem. It neither preserves exact
minimality nor proves a favorable twist exists for the prepared output.
