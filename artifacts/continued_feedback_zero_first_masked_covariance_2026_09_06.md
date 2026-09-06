# Zero-first scalar response: a masked nonlinear covariance closure

Date: 2026-09-06. Status: proposed proof, under independent adversarial
audit. This note concerns arbitrary bounded-operator hollow signings, not
involutions or presumed spectrally flat minimizers. It does not identify a
second threshold of the transported masked field.

## 1. Frozen theorem and actual-energy consequence

Let `B=A/sqrt(n-1)`, where `A` is symmetric, hollow, and has off-diagonal
entries in `{−1,+1}`, and suppose `||B||op<=L` with `L` fixed. Set
`Q=B²`, let `S` have independent Rademacher coordinates, and put `G=BS`.
Let `f` be a fixed bounded, odd, Gaussian-almost-everywhere-continuous
function of ONE scalar variable, with zero first Gaussian coefficient.
Write its normalized Hermite expansion as

```math
f=\sum_{p\ge3,\ p\text{ odd}} f_p h_p,
\qquad R_f=\sum_{p\ge3,\ p\text{ odd}}f_p^2 Q^{\circ p},
\qquad T=B R_f B.
```

Let `H` be fixed bounded, even, and Gaussian-a.e.-continuous. Let `psi`
be fixed bounded, odd, and Lipschitz. Define the actual mixed response

```math
Y=Bf(G),\qquad C_i=H(G_i)\psi(Y_i).
```

Independently sample Gaussian vectors `U~N(0,Q)` and `Z~N(0,T)`, and
define the deterministic covariance matrix

```math
K_{ij}=\mathbb E[H(U_i)H(U_j)]\,
       \mathbb E[\psi(Z_i)\psi(Z_j)].
```

The proposed closure is

```math
\frac1n\left\|\mathbb E[CC^{\mathsf T}]-K\right\|_*\longrightarrow0.       (1)
```

Thus this is also a covariance theorem for the *actual* masked transport:

```math
\frac1n\left\|\mathbb E[(BC)(BC)^{\mathsf T}]-BKB\right\|_*\to0.           (2)
```

In particular it gives both its self-energy and its average squared norm:

```math
\frac{\mathbb E[C^{\mathsf T}BC]}{2n}
 =\frac{\operatorname{Tr}(BK)}{2n}+o(1),\qquad
\frac{\mathbb E\|BC\|^2}{n}=\frac{\operatorname{Tr}(QK)}n+o(1).             (3)
```

No limiting spectral measure or limiting variance profile is assumed.
Every approximation is fixed before taking matrix order to infinity.

There is a particularly short energy form. Put `mu=E H(N)`,
`sigma_i²=T_ii`, and

```math
a_i=\begin{cases}
\mathbb E[N\psi(\sigma_iN)]/\sigma_i,&\sigma_i>0,\\
0,&\sigma_i=0.
\end{cases}
\qquad D_a=\operatorname{diag}(a_i).
```

The choice at zero is immaterial because the corresponding row of `T`
vanishes. The Lipschitz constant bounds `|a_i|`. Then

```math
\frac{\mathbb E[C^{\mathsf T}BC]}{2n}
 =\frac{\mu^2}{2n}\operatorname{Tr}(B D_a T D_a)+o(1).                     (4)
```

If additionally `H>=0`, `|f|+H<=1`, `||psi||infinity<=1`, and
`y psi(y)>=0`, the feasible
endpoints `f(G)+C` and `−f(G)+C` yield

```math
\Lambda(B)\ge
 \mu\,\frac1n\sum_i\mathbb E[\sigma_iN\psi(\sigma_iN)]
 +\frac{\mu^2}{2n}\left|\operatorname{Tr}(B D_a T D_a)\right|-o(1),        (5)
```

where `Lambda(B)=max_x |x^T Bx|/(2n)`. Formula (5) retains an actual
common feedback self-energy, rather than discarding it. It does not
claim a positive uniform lower bound for that common term.

For the discontinuous threshold `psi=sign`, smooth approximation gives
the following formula. The original proof used a fixed positive
variance floor; that assumption is now removed by the exact subset
bound and cutoff argument in
`continued_feedback_threshold_without_variance_floor_2026_09_06.md`.

```math
\Lambda(B)\ge \sqrt{2/\pi}\,\mu\,\frac1n\sum_i\sigma_i
 +\frac{\mu^2}{\pi n}
 \left|\operatorname{Tr}(B D_{1/\sigma}T D_{1/\sigma})\right|-o(1).        (6)
```

The convenient guess `T_ii>=1` fails on actual bounded-op signings
(Section 9), and the companion apex construction defeats every fixed
positive floor for appropriate fixed operator caps. Nevertheless,
the exact count `#{i:T_ii<=epsilon E f(N)²}<=epsilon n` suffices
to justify (6), with the small-variance rows removed before smoothing.

## 2. Uniform kernel bounds for a transported scalar Hermite channel

Write `b_j` for row `j` of `B`. Under Gaussian seed input, the exact
degree-`p` channel `B h_p(BN)` is a normalized multiple integral of

```math
K_{p,i}=\sum_j B_{ij}\,b_j^{\otimes p}.
```

The factors `sqrt(p!)` are accounted for by the normalized multiple
integral convention. Its covariance is exactly `B Q^{circ p} B`.
All fixed degree constants below are harmless.

The global root map has bounded operator norm, because its Gram matrix
is `B Q^{circ p}B` and the correlation Schur multiplier inequality gives
`||Q^{circ p}||op<=||Q||op<=L²`. Each fixed-root proper flattening is

```math
M_r^{\mathsf T}\operatorname{diag}(b_i)M_{p-r},\qquad1\le r<p,
```

where the rows of `M_r` are `b_j^{tensor r}`. Thus its norm is at most
`L²/sqrt(n-1)`, since `||M_r||op²=||Q^{circ r}||op<=L²`.
Row Hilbert norms are bounded. The degree-one old kernel is simply `b_i`.

Repeated marked slots can be deleted at uniformly vanishing Hilbert
cost. For the transported `p`-kernel and a specified pair of repeated
slots,

```math
K_{p,i}(a,a,\mathbf c)
 =\frac1{n-1}\sum_{j\ne a}B_{ij}
                         \prod_{v\in\mathbf c}B_{jv}.
```

The first term after adding back `j=a` is `1/(n-1)` times the
transported `(p−2)`-kernel; its squared Hilbert norm after summing `a`
is `O(1/n)`. The subtracted term has the same bound by direct counting.
The finite union of slot pairs therefore costs `O(n^{-1/2})` in row
Hilbert norm. This handles all `p>=3`, including `p=3` where the lower
kernel is the linear vector `Q_i`.

For fixed products of branch kernels, deletion of cross-branch marked
collisions also has squared Hilbert cost `O(1/n)`: each branch has
maximum marked-slot influence `O(1/n)` from the proper flattening bound
(and directly from flat entries for the degree-one branch). The sum
over a forced shared slot is bounded by a maximum influence times a
bounded total influence. There are finitely many branch pairs.

## 3. Local Wick replacement and transfer to signs

First let `f,H,psi` be polynomials with their indicated parities and
with no degree-one term in `f`. Separate the finitely many transported
input degrees of `f`, treating `Y_p=B h_p(G)` as different coordinates.
At a fixed root they are asymptotically jointly Gaussian, independent
of `G`, with variances `(B Q^{circ p} B)_ii`. Distinct `p` channels
have exactly zero covariance by original input degree.

Every proper contraction of a transported branch has Hilbert norm
`O(n^{-1/2})`; a full contraction into a different larger transported
branch is proper on the latter. A contraction with the degree-one
branch is also proper on the transported branch. The Gaussian product
formula therefore replaces local variance-Hermite monomials by pure
Wick forests with row `L²` error `O(n^{-1/2})`. Whole equal-branch
contractions are precisely those removed by the local Hermites.
Use unnormalized variance-Hermites: no division by a small variance
is needed.

The same replacement holds for Rademacher inputs. One direct route is
to delete all repeated marked slots and use replacement of one seed
at a time for the squared local error. Use a finite vector of summary
coordinates consisting of every squarefree branch AND every squarefree
forest appearing in that local error. Each summary coordinate is
multilinear in the original seeds and thus affine in the particular
seed being replaced. They have fixed degrees, bounded moments under
every Gaussian/sign hybrid law, and each seed influence `O(1/n)`.
Write the squared local error as a fixed outer polynomial of this
summary vector. Signs and Gaussians match their first three moments.
In the fourth-order Taylor remainder of the OUTER polynomial, every
derivative must hit an affine summary coordinate: its own second seed
derivative is zero. Thus the four derivatives supply four
`O(n^{-1/2})` factors; summing over `n` seeds gives `O(1/n)`.
Thus the squared local replacement error remains `O(1/n)`. Fixed-degree
moment bounds follow directly by hypercontractivity, or its elementary
finite-degree moment expansion.

The actual raw `h_p(BS)` differs locally from its squarefree pure
kernel by the same ordered replacement. Multiplication by `B` preserves
the averaged `L²` error because its operator norm is fixed. This last
replacement is NOT yet passed through an unbounded polynomial `psi`.
Sections 3--5 prove their polynomial calculation for the squarefree
transported channels, whose row moments are uniformly bounded. In
Section 6, extend first to bounded `H` and Lipschitz bounded `psi`,
and only then transfer the raw channel using averaged `L²`. For the
cubic channel the identity is especially transparent:

```math
B h_3(BS)=Y_{3,\mathrm{squarefree}}
             -\frac{\sqrt{2/3}}{n-1}Q S.
```

Finally, averaged `L²` replacement changes covariance by `o(n)` in
nuclear norm, since

```math
\|\mathbb E[UV^{\mathsf T}]\|_*
\le\sqrt{\mathbb E\|U\|^2\,\mathbb E\|V\|^2}.
```

It remains to calculate the covariance of exact injective forests.
Gaussian and Rademacher covariance of those forests agree exactly.
After making that identification, restore the undeleted Gaussian
branch-tensor forests at their already proved row `L²` cost. Charge
that cost in normalized nuclear norm. The pairing calculation below
therefore uses unrestricted Gaussian branch tensors; finite
injectivity constraints do not obstruct the exact star factorization.

## 4. The finite branch-pairing classification

A local monomial of `C=H(G)psi(Y)` contains an EVEN number of degree-one
old branches, and an ODD total number of transported branches. This
is preserved when the local variance-Hermite expansion is taken.
Transported branch degrees are odd and at least three.

For a pairing between forests at roots `i,j`, make the bipartite graph
whose vertices are whole branches and whose edges record nonempty
blocks of paired marked slots. Three component types are possible:

1. One whole-branch edge. It matches branches of the same original
   input degree, giving `Q_ij` for two old branches and the appropriate
   transported-channel covariance for two transported branches.
2. A nontrivial star. Only a transported branch can be its split center.
   It has an odd number of leaves, at least three, because its total
   degree and all leaf degrees are odd.
3. A component that is neither a whole edge nor a star. It contains
   a path on four vertices; contract the two disjoint outer edges first.
   The two inner endpoints are split transported branches, so the two
   flattenings give an entrywise `O(1/n)` bound.

Two or more nontrivial star components also give `O(1/n)` entries.
Those matrices have Frobenius norm `O(1)`, hence nuclear norm
`O(sqrt(n))`. The only case not already negligible is EXACTLY ONE
nontrivial star, with every other component a whole-branch edge.

Suppose its center is on the left; let `s` be its number of transported
leaves and `t` its number of old degree-one leaves. Matching all other
branches preserves the parity of each total count. The two local old
counts are even, so `t` is even. The two transported counts are odd,
and the star removes one on the left and `s` on the right, so `s` is
odd. Consequently

```math
s\ge1,\qquad t\ge0\text{ even},\qquad s+t\ge3.
```

This is the complete parity information needed below. For a cubic-only
channel, the lone star is impossible outright. For mixed odd degrees,
it is possible and MUST be bounded rather than declared absent.

## 5. A lone star has a small matrix operator norm

For a transported leaf of input degree `q>=3`, set

```math
S_q=Q^{\circ q}B.
```

Because `Q` is a correlation matrix and `||Q||op<=L²`,

```math
\sum_b|Q_{ab}|^q\le\sum_b Q_{ab}^2\le L^2,
\qquad |(S_q)_{aj}|\le\frac{L^2}{\sqrt{n-1}},
\qquad \|S_q\|_F=O(\sqrt n).
```

After summing the internal slots of the star, its matrix is, up to a
fixed factorial constant,

```math
B W,\qquad W_{aj}=Q_{aj}^{t}\prod_{\ell=1}^{s}(S_{q_\ell})_{aj}.          (7)
```

This formula follows directly by writing its center as
`sum_a B_ia b_a^{tensor p}` and each whole leaf as its tensor sum.
Symmetrization changes only the fixed count of equivalent slot choices.

If `s>=3`, retain one `S_q` factor in Frobenius norm and bound at least
two others entrywise. Then `||W||F=O(n^{-1/2})`.
If `s=1`, parity forces `t>=2`; both absolute row sums and absolute
column sums of `W` are `O(n^{-1/2})`, since the corresponding sums of
`|Q_aj|^t` are at most `L²`. Thus `||W||op=O(n^{-1/2})` as well.
In both cases `||BW||op=O(n^{-1/2})`.

The remaining whole-branch components multiply this matrix entrywise.
Every such covariance factor has a Gram factorization with bounded row
Hilbert norms, hence bounded Schur multiplier norm. Their fixed product
preserves the `O(n^{-1/2})` operator bound. Transposing covers a star
center on the right. This completes EVERY partial pairing class.

## 6. Whole-branch main term and bounded-function extension

The remaining pairings are exactly Gaussian Wick pairings of the local
old coordinate and the finite independent transported-degree channels.
Summing their fixed factorials gives the covariance `K` in Section 1.
Each transported degree contributes `f_p² B Q^{circ p}B`, and their
Gaussian sum has covariance `T`. This proves (1) for fixed polynomials.

For a bounded Lipschitz `psi`, first approximate on a finite partition
of the variance interval `T_ii in [delta,V]`, using fixed odd polynomial
approximants for `psi(sigma N)`. Their coefficients may depend on the
deterministic variance bin; every preceding matrix estimate permits
bounded deterministic row-dependent coefficients. Approximate `H` by
fixed even Gaussian polynomials. Local joint convergence and the fixed
moment bounds transfer these approximations to the actual variables.
For `T_ii<delta`, use `psi(0)=0` and
`E psi(Y_i)²<=Lip(psi)² E Y_i²`, so the total small-variance contribution
is uniformly negligible as `delta` decreases. Then refine bins and
polynomial accuracy, all after the matrix limit. This proves the
bounded `H`, Lipschitz `psi` extension first for the squarefree channels.
Now their averaged `L²` distance from the actual raw `Bf(BS)` tends to
zero, so bounded `H` and Lipschitz `psi` transfer the result to the raw
channel. This proves the stated extension for fixed polynomial `f`.

For bounded zero-first `f`, approximate it in Gaussian `L²` by finite
odd Hermite sums of degrees at least three. At fixed `L`,

```math
\frac1n\mathbb E\|B[f(G)-P(G)]\|^2
\le L^2\frac1n\mathbb E\|f(G)-P(G)\|^2
\longrightarrow L^2\|f-P\|_{L^2(\gamma)}^2.
```

Bounded `H` and Lipschitz `psi` transfer this to `C` in averaged `L²`.
The comparison Gaussian vectors can be coupled by using independent
Gaussian channels of covariance `B Q^{circ p}B` for every odd `p>=3`.
Changing their scalar coefficients costs at most `L^4||f-P||²` in
averaged squared norm. The same covariance Cauchy--Schwarz inequality
then transfers (1). No growing polynomial degree is used at finite
matrix order.

## 7. Reduction of the self-energy to first local response chaos

For polynomial `H,psi`, expand their Gaussian covariance using
unnormalized variance-Hermites. The constant old-mask term is `mu²`.
Every other old-mask term has a factor `Q_ij^a`, `a>=2`. The linear
new-response term is `a_i a_j T_ij`; all others have powers `T_ij^b`,
odd `b>=3`.

Both `Q` and `T` have bounded operator norm and bounded diagonal, so
their Frobenius norms are `O(sqrt(n))`. With bounded row coefficients,

```math
\frac1n\sum_{i\ne j}|B_{ij}|\,|T_{ij}|^b=O(n^{-1/2})\quad(b\ge2),
```

and terms involving `Q_ij^a T_ij^b` with `a>=2,b>=1` obey the same
bound by Frobenius Cauchy--Schwarz (retain one power from each matrix).
Hence only the constant old-mask and linear new-response term survive
the energy trace, proving (4) for polynomials. The approximation from
Section 6 transfers the energy identity by the fixed operator cap.
The Gaussian first-coefficient vector also transfers in weighted `L²`;
equivalently, keep each finite variance-bin approximation before
passing to the limiting coefficients. This proves (4) without an
unjustified termwise summation of an infinite Hermite expansion.

## 8. Endpoint identity and the retained feedback term

The local joint law used above gives independence of the old scalar
`G_i` and the Gaussianized nonlinear field, uniformly in bounded
variance profiles. Consequently

```math
\frac1n\mathbb E[f(G)^{\mathsf T}BC]
=\frac1n\mathbb E[Y^{\mathsf T}H(G)\psi(Y)]
=\mu\frac1n\sum_i\mathbb E[\sigma_iN\psi(\sigma_iN)]+o(1).
```

The old self-energy is `o(1)` after division by `n`, because `f` has
no first Gaussian chaos: its covariance is `R_f+o_* (n)` and
`Tr(B Q^{circ p})/n=O(n^{-1/2})` for `p>=3`.

The two feasible endpoint energies therefore have a common term given
by (4), and opposite cross terms given by the preceding display.
The exact identity `max(|a+b|,|a−b|)=|a|+|b|`, followed by the hollow
cube-to-Boolean rounding, proves (5). A uniform variance floor gives
Gaussian anti-concentration for threshold approximation and yields (6).

## 9. Frozen scalable variance-floor falsifier

Let

```math
C_5=\begin{pmatrix}
-1&-1&1&1&1\\
-1&1&1&-1&1\\
1&1&-1&-1&1\\
1&-1&-1&1&1\\
1&1&1&1&-1
\end{pmatrix},\qquad D=C_5/\sqrt5.
```

Direct integer arithmetic gives

```math
\operatorname{diag}\!\left(C_5[(C_5^2)^{\circ3}]C_5\right)
=(733,733,733,733,533),
```

so the cubic feedback covariance associated with `D` has diagonal
`(733,733,733,733,533)/625`. Its fifth variance is `533/625=.8528<1`.

Let `H_m` be a symmetric Sylvester Hadamard, and hollow the full sign
matrix `C_5 tensor H_m` to obtain an actual signing `A_(5m)`. Its
normalization is operator-close to `D tensor (H_m/sqrt(m))`; the
operator cap tends to `||D||op`, a fixed constant. Squaring, taking a
fixed Schur cube, and conjugating by `B` preserve the vanishing
operator error. Here Schur multiplication by bounded Gram matrices
is uniformly operator bounded. The limiting cubic covariance is

```math
\bigl(D[(D^2)^{\circ3}]D\bigr)\otimes I_m.
```

Thus a positive fifth of roots have variance tending to `533/625`.
This falsifies a pointwise `T_ii>=1` assertion on actual signings.
It does not falsify the known average-standard-deviation inequality;
the limiting average standard deviation here is approximately
`1.0510614698>1`. These trace-zero Hadamard lifts are not suitable
witnesses for nonzero common self-energy, because of their additional
sign-reversal comparator symmetry.

## 10. Scope and next falsifier

The closure is genuinely nonlinear and includes the masked field
`B[H(BS)psi(Bf(BS))]` at covariance/energy level. It requires scalar
old `f` with zero first chaos. It does NOT cover the coherent return
`B²[S K_f]` inside a threshold when the first chaos is nonzero. It
does NOT identify the joint law of the new field with the old history,
nor the energy of another arbitrary threshold of that field.

Next adversarial tests are: exact fixed sign-seed lifts for every
matrix-kernel identity; asymmetric-spectrum actual signings for the
common energy in (4); and a search for a missed lone-star configuration
under any attempted multicoordinate or nonzero-first extension. A
positive lower bound on the displayed trace is not assumed.
