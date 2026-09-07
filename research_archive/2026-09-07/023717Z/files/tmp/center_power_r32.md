# Wave 32 Route 2: low-row centers and the partial power sum

Status: exact center-star reduction, exact normalized-degree and noise
identities, and two sharply scoped walls.  The positive conclusion is a new
minimizer-specific conditional power-moment target.  The walls do **not** use
an exact signing minimizer and hence do not falsify the project.

## 1. Center incidence and completion-cylinder freedom

Fix the target pair `(A,m)` and let `F_S` be the set of projective child
labels which obey the favorable deficit allowance in (10.896).  Thus
`F_S` is a subset of `{+-1}^S/{+-1}`; all bits off `S` are free in a full
completion.  Put

```math
R_2(z)=\|Az\|_2^2,
\qquad
\mathcal C_B=\{z\in\{\pm1\}^n:R_2(z)\le Bn^2\},
\qquad \nu_B=U_{\mathcal C_B}.
```

For a center `z` define the exact fiber distance

```math
a_z(S)=\min_{[y]\in F_S}d_{\rm pr}(z_S,y).
```

Completion-cylinder freedom gives the identity

```math
\boxed{
a_z(S)=\min\{d_{\rm pr}(z,x):x\text{ is a full favorable completion for }S\}.
}
\tag{C32.1}
```

Indeed, restriction cannot increase proximity.  Conversely, orient a
minimizing child label against `z_S` and choose every free outside bit equal
to the corresponding bit of `z`.

For an integer `d>=0`, let

```math
G_z(d)=\{S:a_z(S)\le d\},
\qquad q_z(w,d)=w(G_z(d)).
```

Take the Wave 31 scales

```math
L_0=n^{3/4-c_0},\quad
k_0=\Theta(L_0/\log n),\quad
r=\left\lceil\frac{n\log(2k_0)}{L_0}\right\rceil,
\quad s=\lceil r/T\rceil,
```

where `T<=n^eta` and `eta<c_0`.  Choose

```math
d=\left\lfloor\frac{D_0}{s}\right\rfloor,
\qquad D_0=C_Dk_0.
```

For iid `S_1,...,S_s~w`, sample also `z~nu_B`.  If every `S_j` lies in
`G_z(d)`, choose its closest favorable completion and join all of them to
`z` by a star.  By (C32.1), its total projective Hamming length is at most
`sd<=D_0`.  Therefore (10.910), followed by (10.894), produces one eligible
row-good coset hitting these `s` selectors whenever `B=O(1)` and `D_0=O(k_0)`.
Consequently

```math
\boxed{
\Pr_{w^{\otimes r}}(J_{r,s})
\ge \mathbb E_{z\sim\nu_B}q_z(w,d)^s.
}
\tag{C32.2}
```

Only the first `s` batch positions were used; the other `r-s` positions are
irrelevant to `J_(r,s)`.  Thus the following is an exact sufficient lemma
for (10.907), and hence for convergence:

```math
\boxed{
\inf_w\mathbb E_{z\sim\nu_B}
\left[w\{S:a_z(S)\le \lfloor C_Dk_0/s\rfloor\}\right]^s
\ge e^{-C rL_0}.
}
\tag{C32.3}
```

A weaker and more flexible version replaces the equal radius by the exact
star budget:

```math
\boxed{
\inf_w\mathbb E_{z\sim\nu_B}
\Pr_{S_1,\ldots,S_s\sim w}
\left\{\sum_{j=1}^sa_z(S_j)\le C_Dk_0\right\}
\ge e^{-C rL_0}.
}
\tag{C32.4}
```

This is the sharp low-row-centered missing lemma.  It is a conditional
high-power statement, not an unnormalized degree statement.

## 2. Exact abundance and the normalized-degree route

Every complete signing satisfies

```math
\mathbb E_{z\sim U_n}R_2(z)=n(n-1).
```

Hence Markov gives

```math
|\mathcal C_2|\ge2^{n-1}.
\tag{C32.5}
```

Write

```math
g_S(d)=\nu_2\{z:a_z(S)\le d\}.
```

Fubini and Jensen give the exact lower bound

```math
\boxed{
\mathbb E_{z\sim\nu_2}q_z(w,d)^s
\ge\left(\mathbb E_{S\sim w}g_S(d)\right)^s.
}
\tag{C32.6}
```

Since `sT<=r+T<=2r` for all large `n`, a sufficient normalized-degree
condition would be

```math
\inf_Sg_S(d)\ge e^{-KTL_0}.
\tag{C32.7}
```

It would imply (C32.3) with constant `2K`.  This is the correct
normalization: the raw number of centers incident to `S` is meaningless.

Unfortunately, abundance alone is much too weak for (C32.7).  Let
`N_S=|F_S|` count projective favorable labels and put

```math
V(m,d)=\sum_{j=0}^d\binom mj.
```

For `d<m/2`, the union of the two oriented radius-`d` balls around each
projective label has at most `2N_SV(m,d)` restrictions.  Using (C32.5),

```math
\boxed{
g_S(d)\le
\min\left\{1,\frac{4N_SV(m,d)}{2^m}\right\}.
}
\tag{C32.8}
```

At the project scales,

```math
r=\Theta(n^{1/4+c_0}\log n),\qquad
s=\Theta(n^{1/4+c_0-\eta}\log n),
```

and

```math
d=O\left(
\frac{T L_0^2}{n(\log n)^2}
\right)
=O\left(
\frac{n^{1/2-2c_0+\eta}}{(\log n)^2}
\right).
\tag{C32.9}
```

In particular `d=o(n)` and

```math
\log V(m,d)\le d\log(em/d)
=O\left(\frac{n^{1/2-2c_0+\eta}}{\log n}\right)
=o(L_0).
\tag{C32.10}
```

Here the last exponent gap is
`(3/4-c_0)-(1/2-2c_0+eta)=1/4+c_0-eta>1/4`.
If `N_S=e^{o(n)}` (in particular if `N_S=e^{O(L_0)}`), (C32.8) is
`e^{-Theta(n)}`, whereas (C32.7) asks only for
`e^{-O(TL_0)}` and `TL_0=o(n)`.  Thus the Jensen/degree proof would require

```math
\log N_S\ge m\log2-o(n),
\tag{C32.11}
```

essentially a full exponential family of favorable child labels.  Such a
ground-state multiplicity theorem is neither known nor plausible as a
generic consequence of exactness.  A successful proof of (C32.3) must
instead exploit concentration of the incidence on common centers across
many different selectors.

## 3. Exact conditional quadratic identities

Put `Q=A^2`, so `Q_(ii)=n-1`.  Fix a child label `y` on `S` and fill the
outside coordinates uniformly.  The conditional row-square mean is

```math
\boxed{
\mu_S(y):=\mathbb E[R_2(Z)\mid Z_S=y]
=n(n-1)+\sum_{i\ne j\in S}Q_{ij}y_iy_j
=n(n-1)+y^TQ_{SS}y-m(n-1).
}
\tag{C32.12}
```

If, before the uniform outside fill, exactly `h` uniformly chosen
coordinates of `y` are flipped, then

```math
\boxed{
\mathbb E R_2(Z)
=n(n-1)+\theta_{m,h}[\mu_S(y)-n(n-1)],
\quad
\theta_{m,h}=1-\frac{4h(m-h)}{m(m-1)}.
}
\tag{C32.13}
```

For iid flips of the child coordinates with sign retention `rho`, the same
identity holds with `theta_(m,h)` replaced by `rho^2`.  These follow by
killing every monomial involving an outside coordinate and multiplying
each surviving off-diagonal child monomial by its two-coordinate noise
correlation.

If `mu_S(y)=O(n^2)`, Markov does show that a constant fraction of the exact
completion cylinder has `R_2=O(n^2)`.  But (C32.13) also exposes the wall:
for `h<=d=o(n)`, `theta_(m,h)=1-o(1)`.  A small Hamming ball does not wash a
large positive conditional excess down to the unconditional `n(n-1)`
scale.  Moreover, even a constant conditional fraction still occupies only
about `2^{-m}` of all centers and cannot prove (C32.7).

## 4. A complete-signing conditional-fiber wall at the correct coarse scales

The preceding obstruction is not merely an arbitrary set-system artifact.
There is an infinite probabilistic-method family of complete signings with

```math
Q(A)=O(n^{3/2}),\qquad \|A\|_{op}=O(n^{3/4}),
```

and a fixed-density child whose project-favorable completion cylinder is
disjoint from every `R_2=O(n^2)` center even after an `o(n^(3/4))` Hamming
thickening.  These signings are **not exact global minimizers**; that is the
essential scope of the wall.

Let `b=ceil(K n^(3/4))`, choose `B subset S` of size `b`, and write

```math
A=P_B+R,
```

where `P_B` is the adjacency matrix of the positive clique on `B`, while
`R` is zero on `B times B` and has symmetric Rademacher entries elsewhere.
A standard random-matrix norm bound supplies a deterministic realization
with `||R||op<=C sqrt(n)`.  The resulting `A` is a complete signing, and

```math
Q(A)\le n\|A\|_{op}
\le b n+C n^{3/2}=O(n^{7/4}),
```

is too crude; the correct quadratic bound uses the two pieces separately:

```math
\boxed{
Q(A)\le \max_x|x^TP_Bx|+\max_x|x^TRx|
\le b^2+C n^{3/2}=O(n^{3/2}).
}
\tag{C32.14}
```

The operator norm is at most `b+C sqrt(n)=O(n^(3/4))`, matching the only
coarse spectral scale used for exact minimizers.

For `S` of any fixed density containing `B`, fixing all spins on `B` equal
and averaging the other spins makes the `R` contribution zero.  Hence the
positive child optimum is at least `b(b-1)`.  Every negative-oriented
energy is at most `b+C n^(3/2)`.  Taking fixed `K` with `K^2>C` makes the
absolute child optimum positive for all large `n`.

More quantitatively, if a positive child label `y` has minority size `t`
on `B`, align that minority while retaining every outside bit.  The clique
gain is `4t(b-t)`, whereas

```math
|y'^TRy'-y^TRy|
\le\|R\|_{op}\|y'-y\|_2\|y'+y\|_2
\le4Cn\sqrt t.
```

Thus an exact child ground has

```math
t\le4C^2n^2/b^2=O(n^{1/2}).
\tag{C32.15}
```

The same conclusion follows directly from the energy deficit.  If
`M_B=sum_(i in B)y_i` and a positive-oriented label has deficit at most
`theta b^2`, then

```math
M_B^2\ge(1-\theta-C/K^2-o(1))b^2.
\tag{C32.15a}
```

A negative-oriented label has deficit at least
`(1-C/K^2-o(1))b^2`.  This covers the actual project allowance, not merely
exact grounds.  Indeed

```math
\max_{0<p<1}(p^{3/2}-p^2)=\frac{27}{256},
```

`p^(3/2)-p_2<=27/256+o(1)`, `q_n<=Q(A)`, and (C32.14) gives

```math
B_{n,m}+O(n^{3/2-c_0})
\le\left[\frac{27}{256}(1+C/K^2)+o(1)\right]b^2.
```

Choose the fixed `K` sufficiently large.  The bracket is, for example,
less than `1/4`, and (C32.15a) then gives `|M_B|>=c b` for every favorable
orientation.

If `z_S` is within `d=o(b)` of such a label, then its clique magnetization
still has magnitude at least `c b`.  Since

```math
\|P_Bz\|_2^2=(b-2)\left(\sum_{i\in B}z_i\right)^2+b,
\qquad \|Rz\|_2\le Cn,
```

the reverse triangle inequality gives

```math
\boxed{
R_2(z)=\|Az\|_2^2\ge \{c b^{3/2}-Cn\}^2
=\Omega(n^{9/4}).
}
\tag{C32.16}
```

Thus no `O(n^2)` center lies in the thickened completion cylinder.  Yet
(C32.5) says at least half of all spins are `R_2<=2n(n-1)` centers.  This
shows that global center abundance, free outside completion bits, the
correct `Q` and operator exponents, and exact child optimality still do not
give even one useful conditional fiber.  A proof must use the fact that
`A` globally minimizes `Q` among complete signings, beyond those coarse
consequences.

At the actual project radius, (C32.9) is `o(n^(3/4))`, and even the entire
forest budget `k_0=n^(3/4-c_0)/log n` is `o(n^(3/4))`, so the separation is
strict for every `c_0>0` and every `eta<c_0`.

## 5. A power-moment wall independent of conditional row cost

Even if every center is declared low-row, Hamming volume and local planting
do not imply (C32.3).  On the selector slice, independently assign each
`S` one uniform projective label `Y_S`, and use its radius-`d` completion
cylinder.  For every fixed center `z`, under the random labels,

```math
q_0:=\Pr\{d_{pr}(z_S,Y_S)\le d\}
=\frac{2V(m,d)}{2^m}=e^{-p(\log2)n+o(n)}.
```

Let `N=binom(n,m)=e^{H(p)n+o(n)}` and take uniform selector law.  For fixed
`z`, its hit degree is `Bin(N,q_0)`.  For any

```math
0<alpha<min\{H(p),p\log2\},
```

the binomial upper tail at `N e^(-alpha n)` followed by a union bound over
all `2^n` centers shows that some deterministic label assignment obeys

```math
\max_zq_z(w,d)\le e^{-alpha n}.
\tag{C32.17}
```

Therefore every center law satisfies

```math
\mathbb E_\nu q_z(w,d)^s\le e^{-alpha ns}
=e^{-\omega(rL_0)},
\tag{C32.18}
```

because

```math
\frac{ns}{rL_0}=\Theta\left(\frac n{TL_0}\right)
=\Theta(n^{1/4+c_0-\eta})\longrightarrow\infty.
```

This second wall is abstract: the independently prescribed labels need not
be child grounds of one signing.  It shows exactly why the remaining input
must correlate favorable fibers across different selectors, rather than
merely thicken every pointwise cylinder or count abundant centers.

## 6. Resulting frontier

The low-row-center route is now a clean alternative formulation, but not a
proof.  Its strongest exact sufficient target is the total-cost moment
(C32.4); the equal-radius power moment (C32.3) is a convenient stronger
version.  The following weaker inputs are formally insufficient:

1. global density `nu(C_2)>=1/2`;
2. one or `e^{O(L_0)}` favorable child labels per selector;
3. a constant low-row fraction inside each individual completion cylinder;
4. the bounds `Q(A)=O(n^(3/2))` and `||A||op=O(n^(3/4))`; and
5. pointwise completion-cylinder freedom.

The exact missing minimizer-specific statement is that, against every
selector law `w`, favorable child fibers put `e^{-O(rL_0)}` joint mass on
`s=r/T` labels having total distance `O(k_0)` from a common `O(n^2)`-row
center.  Equivalently, exact global signing minimality must rule out both
the planted quadratic mode in Section 4 and the independent-label
dispersion in Section 5 at the conditional high-power level.  An ordinary
first moment or normalized degree estimate cannot do this unless it also
proves essentially maximal child-ground multiplicity.
