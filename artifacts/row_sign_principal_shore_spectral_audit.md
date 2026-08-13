# Row-sign recoupling: conference shore spectra and their limitation

## Status and theorem-first protocol

This note audits the explicit asymmetric law

```math
X\sim\operatorname{Unif}\{\pm1\}^n,
\qquad Y=\operatorname{sign}(AX),                 \tag{1}
```

through the same-spin recoupling theorem.  All energies use the doubled
normalization `Q(A)=max_z |z^T A z|`.  The exact statements in Sections
1--4 were derived before any finite experiment in Section 5.

The positive result is an exact row-law theorem and a correct-scale
nuclear/projector certificate for symmetric conference matrices.  The
negative result is that its leading constant is too small: conference
identity, interlacing, trace norm, and bounded-rank projector rounding alone
cannot make the row-sign recoupling defect subleading.  No project bound is
improved.

## 1. Exact switched-shore formulation

Assume first that every coordinate of `AX` is nonzero.  Put

```math
D_X=\operatorname{Diag}(X),\qquad C_X=D_XAD_X,
\qquad \ell=C_X\mathbf1,
\qquad g=\operatorname{sign}(\ell)=X\circ Y.      \tag{2}
```

Let `I={i:ell_i>0}` and `J={i:ell_i<0}`, and write

```math
P=\mathbf1_I^TC_X[I]\mathbf1_I,qquad
R=\mathbf1_J^TC_X[J]\mathbf1_J.                 \tag{3}
```

Switching preserves every principal eigenvalue, so the spectral
certificates of `C_X[I]` and `C_X[J]` equal those of `A[I]` and `A[J]`.
The row-sign bilinear response and agreement/disagreement identity are

```math
\boxed{
X^TAY=\sum_i|(AX)_i|=\sum_i|\ell_i|=P-R.}         \tag{4}
```

For a zero-diagonal symmetric block `D`, retain the common one-sided
nuclear certificate

```math
S(D)=\left[{\lVert D\rVert_*\over\pi}
-\left(1-{2\over\pi}\right)|D|\right]_+,         \tag{5}
```

and let `Gamma_sigma(D)` denote the stronger sign-specific positive-projector
certificate from the recoupling theorem.

### Theorem 1 (universal row-law spectral recoupling)

For every realization in (1),

```math
\boxed{
Q(A)\ge {X^TAY+S(A[I])+S(A[J])\over2}.}           \tag{6}
```

Indeed, if `PR<0`, this is exactly the averaged two-shore recoupling
theorem.  If `PR>=0`, that theorem gives `Q(A)>=X^TAY`.  Independently,
choose positive-energy Boolean witnesses on both principal shores.  Of the
two full spins obtained by changing the global sign of the second witness,
one has absolute energy at least the sum of the two internal energies;
hence `Q(A)>=S(A[I])+S(A[J])`.  The maximum of these two lower bounds is at
least their average, proving (6).

When `PR<0`, the projector defect additionally obeys the exact upper
reduction

```math
\Delta_\Gamma
\le {1\over2}\left[
X^TAY-Gamma_{\operatorname{sgn}R}(A[I])
-\Gamma_{\operatorname{sgn}P}(A[J])
\right]_+.                                      \tag{7}
```

This is just `min(u_+,v_+) <= [(u+v)/2]_+`, using
`X^TAY=|P|+|R|`.  Thus a subleading defect would follow from a
near-lossless aligned projector sum, but a merely correct-scale projector
sum is not enough.

## 2. Exact conference row-sign law

Let `C` be a symmetric conference signing of order `n`,

```math
C^2=(n-1)I,
```

and set `n-2=2m`.  Since `n-1` is odd, the row fields in (1) never vanish.
Define

```math
\eta_n={\binom{2m}{m}\over2^{2m}},
\qquad
\rho_n=
\begin{cases}
\left[\binom m{m/2}/2^m\right]^2,&m\text{ even},\\
0,&m\text{ odd}.
\end{cases}                                      \tag{8}
```

### Theorem 2 (exact first and second shore moments)

For the row-sign law on `C`, with `i=|I|` and `j=n-i`,

```math
\boxed{\begin{aligned}
\mathbb E[X^TCY]&=n(n-1)\eta_n,\\
\mathbb EP&={n(n-1)\eta_n\over2},
&\mathbb ER&=-{n(n-1)\eta_n\over2},\\
\mathbb Ei&={n\over2},
&\operatorname{Var}(i)&={n+n(n-1)\rho_n\over4}.
\end{aligned}}                                   \tag{9}
```

In particular the selected shores are balanced up to root-order
fluctuations, while the two witnessed internal energies have opposite
expectations of order `n^(3/2)`.

#### Proof

Fix an edge `ij` and put `e=c_ij X_iX_j`.  Orthogonality of rows `i,j`
implies that among the remaining `2m` coordinates, the products
`c_ij c_ik c_jk` take each sign exactly `m` times.  After an invertible
Rademacher change of variables, `e` is independent of two independent sums
`U,V` of `m` Rademachers and

```math
\ell_i=e+U+V,
\qquad \ell_j=e(1+U-V).                           \tag{10}
```

Writing `g_i=sign(ell_i)`, direct averaging over `e` gives

```math
\mathbb E[e g_i]=\eta_n,qquad
\mathbb E[e g_i g_j]=0,qquad
\mathbb E[g_i g_j]=\rho_n,qquad
\mathbb E g_i=0.                                 \tag{11}
```

The last pair-correlation formula follows because the only uncancelled
case has `U=V=0`; it is absent when `m` is odd.  Finally,

```math
P-R=\sum_{a<b}e_{ab}(g_a+g_b),
\qquad
P+R=\sum_{a<b}e_{ab}(1+g_ag_b),                  \tag{12}
```

and `i=(n+sum_a g_a)/2`.  Summing (11) proves (9).

An independent correction-pass audit enumerated the finite `(e,U,V)` law for
every `m=1,...,8` and exactly enumerated all projective spins of the saved
conference matrices of orders `6,10,14`.  It reproduced every identity in
(9), including `E P=-E R`, and the variance formula.  No correction to the
conference moments or their doubled normalization was required.

## 3. Correct-scale expected nuclear and projector certificates

Put `r=sqrt(n-1)`.  Every principal conference compression `D=C[S]` of
order `s` satisfies

```math
\lVert D\rVert_{\rm op}\le r,qquad
\lVert D\rVert_F^2=s(s-1),
```

so Schatten interpolation gives

```math
\lVert D\rVert_*\ge {s(s-1)\over r}.             \tag{13}
```

Therefore, realization by realization,

```math
S(C[I])+S(C[J])
\ge {i(i-1)+j(j-1)\over\pi r}
-\left(1-{2\over\pi}\right)n.                  \tag{14}
```

The same raw lower bound holds for every pair of sign-specific projector
certificates, since taking `theta=1` gives

```math
\Gamma_\sigma(D)\ge {\lVert D\rVert_*\over\pi}
-\left(1-{2\over\pi}\right)|D|.                 \tag{15}
```

Theorem 2 yields the exact expectation

```math
\mathbb E[i(i-1)+j(j-1)]
={n^2-n+n(n-1)\rho_n\over2}.                     \tag{16}
```

Consequently

```math
\boxed{
\mathbb E[S(C[I])+S(C[J])]
\ge {n^2-n+n(n-1)\rho_n\over2\pi\sqrt{n-1}}
-\left(1-{2\over\pi}\right)n.}                 \tag{17}
```

For the nontrivial symmetric-conference orders, `m` is even and
`rho_n=Theta(1/n)`.  Hence the leading term in (17) is exactly

```math
{1\over2\pi}n^{3/2}+O(n).                        \tag{18}
```

Combining (6), (9), and (17) gives the explicit deterministic conference
bound

```math
Q(C)\ge {1\over2}\left[
n(n-1)\eta_n
+{n^2-n+n(n-1)\rho_n\over2\pi\sqrt{n-1}}
-\left(1-{2\over\pi}\right)n
\right].                                        \tag{19}
```

Its asymptotic leading constant is

```math
{1\over2}\sqrt{2\over\pi}+{1\over4\pi}
=0.4785\ldots                                    \tag{20}
```

in doubled normalization.  This is below the project's existing doubled
lower constant.  Thus the conference identity proves that the selected
shore certificates are on the correct scale, but not that the recoupling
shortfall is subleading.

## 4. Spectral-information barrier

The coefficient in (13) cannot be improved at leading order from conference
identity, zero diagonal, Frobenius mass, interlacing, and complementary
compression pairing alone.  Here is an exact abstract realization.

Let `n=2s`, `r=sqrt(n-1)`, and

```math
q={s(s-1)\over n-1},\qquad
t=\lfloor q/2\rfloor,
\qquad a=\sqrt{(q-2t)/2}.                         \tag{21}
```

Let `Lambda` have eigenvalues `+1,-1` each with multiplicity `t`, then
`+a,-a`, and zeros in the remaining positions.  For all sufficiently large
`s` this list has length at most `s`.  It has trace zero, squared norm `q`,
and

```math
\lVert\Lambda\rVert_*=q+2(a-a^2)\le q+{1\over2}. \tag{22}
```

With `H=sqrt(I-Lambda^2)`, the block matrix

```math
U_0=\begin{pmatrix}\Lambda&H\\H&-\Lambda\end{pmatrix}
```

is a symmetric involution.  By the Schur--Horn theorem, trace zero permits
orthogonal conjugations of both diagonal blocks to zero diagonal.  Applying
those two conjugations to `U_0` produces a symmetric orthogonal `U` with
zero diagonal and balanced complementary compressions `D,F` satisfying all
the spectral identities of a normalized conference matrix, including

```math
\lVert rD\rVert_F^2=\lVert rF\rVert_F^2=s(s-1),
\quad
\lVert rD\rVert_*+\lVert rF\rVert_*
={2s(s-1)\over r}+O(r).                          \tag{23}
```

Thus the lower bound (13) is sharp up to `O(sqrt(n))` within the complete
spectral data available to interlacing and the conference identity.  The
model is not asserted to be a flat sign matrix.  Precisely for that reason,
any stronger theorem must use flat-entry arithmetic or a quantitative
correlation between the local-field selection vector and the eigenvectors;
trace norm, interlacing, or restricted invertibility applied only through
singular values cannot supply it.

The moment theorem (9) gives the tractable part of local-field selection:
balance, pair correlations, and mean internal energies.  None of these
quantities controls the eigenvectors or the nuclear norm above (13).  The
finite audit below tests whether an unproved spectral bias is nevertheless
visible on saved conference matrices; it is diagnostic only.

## 5. Cross-aware collapsed shore

The exact recoupling proof actually permits retaining the cross field.  Fix
the agreement witness `p` and put

```math
D=C[J],\qquad h=C[J,I]p,\qquad
E_J=\begin{pmatrix}D&h\\h^T&0\end{pmatrix}.       \tag{24}
```

For every `r in {+-1}^J` and `t in {+-1}`,

```math
(tp,r)^TC(tp,r)=P+(r,t)^TE_J(r,t).               \tag{25}
```

Thus the strongest collapsed-shore version replaces the internal
one-sided cap of `D` by the aligned one-sided cap of `E_J`:

```math
Q(C)\ge |P|+C_{\operatorname{sgn}P}(E_J).
```

There is a symmetric bound anchored on `q` over `J`.  This is still a
strictly smaller Boolean problem than the original one, but its entries are
weighted: the last column consists of integer local fields, not signs.
The displayed aligned form is for `P != 0`; when `P=0`, replace its last
term by `max(C_+(E_J),C_-(E_J))`, and analogously for a zero `R`.

For a conference matrix the block equations give the exact identities

```math
\boxed{\begin{aligned}
\lVert h\rVert_2^2
&=(n-1)|I|-\lVert C[I]p\rVert_2^2,\\
\lVert E_J\rVert_F^2
&=|J|(|J|-1)+2\lVert h\rVert_2^2,\\
E_J^2&=\begin{pmatrix}
(n-1)I-C[J,I]C[I,J]+hh^T&Dh\\
h^TD&\lVert h\rVert_2^2
\end{pmatrix}.
\end{aligned}}                                    \tag{26}
```

Under the row-sign selection, after switching so that `p=1_I`,

```math
h_j=\ell_j-(D\mathbf1)_j
\qquad(j\in J).                                  \tag{27}
```

The row-sign condition gives `ell_j<0`, but it imposes no sign on `h_j`
relative to the internal row sum.  In particular it supplies no immediate
deterministic lower bound on `||h||_2`.

### Weighted projector theorem and its obstruction

For an arbitrary symmetric zero-diagonal weighted matrix `E`, the projector
rounding proof gives the exact audited extension

```math
C_\sigma(E)\ge\max_{0\le\theta\le1}{2\over\pi}
\left[
\theta\operatorname{tr}((\sigma E)_+)
-\left({\pi\over2}-1\right)\theta^2
\sum_{u\ne v}|e_{uv}|\Pi_{uv}^2
\right],                                         \tag{28}
```

where `Pi` projects onto the positive eigenspace of `sigma E`.  The absolute
weights in the remainder are necessary: the proof bounds each unknown sign
of `e_uv(arcsin(theta Pi_uv)-theta Pi_uv)` separately.  For a signing they
collapse to at most `rank(Pi)`.  For (24), however, the extra contribution is

```math
2\sum_{j\in J}|h_j|\Pi_{j,*}^2,                  \tag{29}
```

and conference identity does not bound this by `O(n)`.  The generic estimate
is `O(||h||_infinity rank Pi)`, which may be `Theta(n^(3/2))` and consumes
the desired leading gain.

Trace norm has the same unresolved interface.  From (26) and
`||E_J||op`, Schatten interpolation gives only

```math
\lVert E_J\rVert_*
\ge {|J|(|J|-1)+2\lVert h\rVert_2^2
       \over\lVert E_J\rVert_{op}},              \tag{30}
```

while `||E_J||op` is not controlled by `sqrt(n-1)` because collapsing `I`
uses the unnormalized vector `p`.  The block triangle inequality and (26)
give only

```math
\lVert E_J\rVert_{op}
\le\sqrt{n-1}+\lVert h\rVert_2
\le\sqrt{n-1}+\sqrt{(n-1)|I|},
```

which is too large by `sqrt(n)` on a balanced shore.  Hence neither
(28) nor (30) presently gives a uniform `n^(3/2)` aligned certificate.

This identifies an exact positive successor lemma rather than a vague
spectral request:

> **Collapsed-shore leverage lemma.** For the row-sign law on every
> symmetric conference matrix, prove either
> `E ||E_J||_* = Omega(n^(3/2))` together with an `O(n)` weighted rounding
> remainder in (28), or directly prove an aligned one-sided cap of
> `E_J` equal to `|R|-o(n^(3/2))` (and symmetrically for `E_I`).

The first clause is not implied by the conference singular values alone;
it needs a joint eigenvector/local-field statement.  The second is exactly
the remaining recoupling shortfall, but on an augmented shore of only about
`n/2+1` vertices.

## 6. Finite audit after the theorem

Define the exact restricted row-law defect by

```math
\Delta_{\rm aug}^{\rm exact}
=\left[
X^TAY-
\max\left\{
|P|+C_{\operatorname{sgn}P}(E_J),
|R|+C_{\operatorname{sgn}R}(E_I)
\right\}
\right]_+.                                      \tag{31}
```

The positive part is essential because a recoupled common spin may have
larger energy than the original bilinear response.  The exact theorem is

```math
\boxed{Q(A)\ge X^TAY-\Delta_{\rm aug}^{\rm exact}.} \tag{32}
```

Replacing each exact cap in (31) by its weighted-projector lower bound (28)
defines `Delta_aug^proj` and preserves (32) with the larger polynomial
defect.

The restricted defect is genuinely weaker than computing the full cap.
Conditioned on `X`, each branch permits only two common values on the whole
anchor shore and arbitrary spins on the free shore.  Computing the maximum
of **both** anchor branches in (31) costs
`O(2^(|I|+1)+2^(|J|+1))=O(2^max(|I|,|J|))`, rather than `2^n`.
One may deliberately keep only the branch with the smaller free shore, at
cost `O(2^min(|I|,|J|))`; that gives a weaker but still valid defect and is
not the exact maximum displayed in (31).  The retained state is `(D,P,h)`,
where `D` is the entire free principal block, not merely `(P,h)`.  Thus the
reduction is strict but not bounded-complexity: balanced shores still leave
`2^(n/2+O(1))` states and a dense block of order about `n/2`.

Only after Sections 1--5 and (31)--(32) were fixed, the reproducible auditor
`computations/audit_row_sign_recoupling_law.py` was run.  Its full output is
`computations/results/row_sign_recoupling_law_audit.json`, with canonical
payload SHA-256
`591acb823f980cd88ad1689081b1032e9071a02f757b8f23d0f1ad0a7f570cb7`.
The table reports means divided by `n^(3/2)`.  Orders `6,10,14` are exact;
the remaining rows are reproducible Monte Carlo.

| order | mode/count | seed | row response | exact augmented defect | weighted-projector defect |
|---:|---|---:|---:|---:|---:|
| 6 | exact/32 | -- | 0.765466 | 0.085052 | 0.185527 |
| 10 | exact/512 | -- | 0.778217 | 0 | 0.206590 |
| 14 | exact/8192 | -- | 0.783775 | 0.011875 | 0.194808 |
| 18 | MC/32768 | 26081318 | 0.786885 | 0.000170 | 0.197557 |
| 26 | MC/16384 | 26081326 | 0.790251 | -- | 0.206169 |
| 98 | MC/2048 | 26081398 | 0.795852 | -- | 0.246547 |

To improve the existing doubled constant `c_*=0.6729867...` from the row
response `sqrt(2/pi)=0.7978845...`, the asymptotic defect must be below

```math
\sqrt{2/\pi}-c_*=0.1248978\ldots .              \tag{33}
```

The weighted-projector defect misses (33) by a large margin.  By contrast,
the exact augmented defect is `0` and `0.011875` at the exactly enumerated
conference orders `10` and `14`; the same auditor gives `0.02643` on the
saved exact order-10 minimizer.  These finite values are far below (33).
Thus the row-sign law and exact restricted target survive; only the present
weighted-projector surrogate is falsified.  The order-18 exact-cap column is
sampled, not a finite certificate.

## 7. Research judgment

The row-sign response is exactly

```math
\mathbb E[X^TAY]
=n\,\mathbb E\left|
\sum_{k=1}^{n-1}\varepsilon_k
\right|
=\left(\sqrt{2/\pi}+o(1)\right)n^{3/2}.         \tag{34}
```

The exact sufficient lemma is now

```math
\mathbb E_X\Delta_{\rm aug}^{\rm exact}(A,X)
\le
\left(\sqrt{2/\pi}-c_*-\eta\right)n^{3/2}
+o(n^{3/2})                                      \tag{35}
```

for some `eta>0`, uniformly over signings with project-scale cap.  Equations
(32) and (34) would then improve the current lower constant by `eta`; a
subleading defect is stronger than necessary.

Conference moments show that the hard opposite-sign branch is typical.
Nuclear/interlacing information is spectrally sharp in the abstract model of
Section 4, and the weighted projector misses the required finite threshold.
A proof of (35) must use the sign-selection law `I={ell>0}` and flat entries
to control the exact restricted cap, not merely improve a trace-norm
constant.  Conversely, a scalable family with expected exact defect above
the threshold in (33) would falsify this leading target.

## 8. Greedy dynamics: exact progress identity and present barrier

There is a natural deterministic algorithm for the free-shore maximum in
(31).  Fix the sign `sigma=sgn(P)`, absorb the collapsed coordinate sign
into the absolute field, and maximize

```math
F_\sigma(r)=\sigma r^TDr+2|h^Tr|.                \tag{36}
```

Equivalently retain a spin `t` and maximize
`G_sigma(r,t)=sigma(r^TDr+2t h^Tr)` by one-coordinate flips, starting with
`t=1, r=sign(sigma h)`.  If `L_r=Dr+t h`, flipping coordinate `r_j`
has the exact gain

```math
G_\sigma(r-2r_je_j,t)-G_\sigma(r,t)
=-4\sigma r_jL_{r,j}.                            \tag{37}
```

Hence terminal one-flip stability is precisely

```math
\boxed{\sigma r_j(Dr+t h)_j\ge0\quad(j\in J),
\qquad \sigma t h^Tr\ge0}                       \tag{38}
```

if the collapsed coordinate is also allowed to flip.  Summing the first
conditions gives the rigorous but tautological lower bound

```math
\sigma(r^TDr+t h^Tr)=\sum_j| (Dr+t h)_j|,
\quad
G_\sigma(r,t)=\sum_j|(Dr+t h)_j|+\sigma t h^Tr. \tag{39}
```

The cumulative-gain identity is equally exact.  If coordinates
`j_1,...,j_T` are flipped and `g_t>0` is the chosen gain at step `t`, then

```math
G_\sigma(r^{(T)},t^{(T)})
=G_\sigma(r^{(0)},t^{(0)})+\sum_{u=0}^{T-1}g_u. \tag{40}
```

These identities explain why the dynamics is certifying, but do not yet
bound its deficit.  Conference orthogonality bounds `||Dr+h||_2` through
the block relations behind (26), whereas (39) needs a lower `l1` bound *correlated with the terminal
signs*.  The norm implication `||v||_1>=||v||_2` supplies only order `n`,
and `||v||_1>=||v||_2^2/||v||infinity` loses the required scale when the
maximum field is order `sqrt(n)`.  Also, a coordinate may flip more than
once, so (40) has no telescoping expression in the original conference
rows.  One-flip stability by itself is therefore insufficient for an
`O(n)` loss theorem.

After this derivation, a separate reproducible computation reported that
best-improvement ascent from the stated initialization has sampled
normalized defects `0.001278, 0.000124, 0.000133` at orders `30,54,62`, and
`0.000034, 0.000009, 0` at orders `74,90,98`.  This is much stronger
empirical evidence than the raw spectral certificates, but is not a theorem:
local optimality (38) admits no present conference-specific quantitative
margin.

The exact next lemma is now algorithmic and falsifiable:

> **Greedy conference lemma.** For every row-sign shore of every symmetric
> conference matrix, best-improvement dynamics initialized by
> `r=sign(sigma h)` terminates at `(r,t)` satisfying
> `G_sigma(r,t)>=-sigma R-O(n)`; an expectation version with
> `o(n^(3/2))` loss is
> already enough.

A proof needs one additional invariant beyond energy monotonicity: for
example, a bound on repeated flips or a conference-orthogonality estimate
forcing the terminal `l1` local-field sum in (39) to dominate `-sigma R`.
Without such an invariant, (37)--(40) merely restate coordinate ascent and
do not certify the observed tiny loss.
