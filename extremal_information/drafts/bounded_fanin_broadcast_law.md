# A bounded-fan-in law for state-dependent broadcast

Status: task-local rigorous draft. This isolates the exact incidence resource
used by the state-dependent Gram construction and determines its sharp scale
for unrestricted quadratic signings.

## 1. Fixed-child contextual model

Let `Z={0,1}^h` and let `Omega` be a finite configuration space. A child has
a fixed bounded-atom presentation

```math
H_z(x)=P(x)+sum_(e=1)^E c_e(z)phi_e(x),
\qquad |phi_e(x)|\le 1.                              \tag{BF.1}
```

Choose `I_e subseteq [h]` so that `c_e` depends only on `z|_(I_e)`, and put

```math
t=max_e|I_e|,
\qquad omega_e=max_z c_e(z)-min_z c_e(z).             \tag{BF.2}
```

The public futures `(C_theta)_theta` are fixed independently of `z`. Let
`Lambda_theta` be one-Lipschitz in uniform norm and define

```math
R_theta(z)=Lambda_theta(H_z+C_theta),
\qquad d(z,z')=sup_theta|R_theta(z)-R_theta(z')|.     \tag{BF.3}
```

This includes `max`, `min`, `max |.|`, and arbitrarily strong pinning. It
also includes `C_T=-H_T`: `T` chooses a member of a language fixed before
the unknown child is supplied. It does not include a hidden-dependent atom
presentation that itself changes with the query.

### Theorem BF.1 (bounded-fan-in incidence law)

For every `z in Z`,

```math
sum_(i=1)^h d(z,z+e_i)
\le sum_(e=1)^E |I_e|omega_e.                         \tag{BF.4}
```

If every neighbour has contextual distance at least `epsilon S`, then

```math
h epsilon S
\le sum_e |I_e|omega_e
\le t sum_e omega_e.                                 \tag{BF.5}
```

In particular, if `|c_e|<=B`, then

```math
h\le {2BtE\over epsilon S}.                          \tag{BF.6}
```

For distinct quadratic pair atoms on `N` Boolean variables,
`E<=binom(N,2)`. At scale `S=N^(3/2)`,

```math
h\le {Bt(N-1)\over epsilon sqrt(N)}
 < {Bt\over epsilon}sqrt(N).                         \tag{BF.7}
```

More generally, coordinate-dependent gaps `epsilon_i S` obey

```math
S sum_i epsilon_i\le sum_e |I_e|omega_e.             \tag{BF.8}
```

#### Proof

The future cancels before optimization and `Lambda_theta` is
one-Lipschitz, so

```math
d(z,z+e_i)
\le ||H_z-H_(z+e_i)||_infty
\le sum_(e:i in I_e)|c_e(z)-c_e(z+e_i)|.             \tag{BF.9}
```

Sum over `i`. One atom occurs at most `|I_e|` times and each difference is
at most `omega_e`. This proves (BF.4), and the other claims follow. `square`

The future language may be enormous. State-dependent futures, however,
invalidate the pointwise cancellation and must have their incidence charged
as part of the child--future interface.

## 2. Neighbour visibility is not automatically information

If all `2^h` response vectors are pairwise `epsilon S`-separated, any
summary with uniform response error below `epsilon S/2` needs `2^h` states,
hence `h` bits, and (BF.6) bounds that information. Neighbour separation
alone is weaker: a response cell of diameter below `epsilon S` need only be
an independent set in the `h`-cube, so the universal conclusion is one bit.

### Proposition BF.2 (exact neighbour sharpness)

Partition `gd` distinct quadratic edges into cells `E_1,...,E_g` of size
`d`. Divide `h=gt` hidden bits into groups `G_1,...,G_g` of size `t`, and put

```math
sigma_j(z)=(-1)^(sum_(i in G_j)z_i),
\qquad
H_z(x)=sum_(j=1)^g sigma_j(z)sum_({u,v}in E_j)x_ux_v. \tag{BF.10}
```

Give unused edges fixed signs and declare every continuation `-H_w`. Then
every atom has fan-in `t`, and for every `z,i`,

```math
d(z,z+e_i)=2d.                                      \tag{BF.11}
```

Equality holds in (BF.4): both sides are `2tgd`.

#### Proof

If `i in G_j`, only cell `E_j` changes, with all coefficient differences
having the same sign. In context `-H_z`, the two responses are zero and

```math
max_x 2|sum_(\{u,v\}in E_j)x_ux_v|=2d,              \tag{BF.12}
```

attained by the all-one spin. Uniform norm gives the reverse inequality.
Every atom has oscillation two and dependency size `t`. `square`

Taking `gd=Theta(N^2)` and `d=Theta(N^(3/2))` gives
`h=Theta(t sqrt(N))`. The all-plus child has cap `Theta(N^2)`, so this is
sharpness for bounded atoms, not for a spectrally flat subclass. For `t>1`
the parity construction has only `g` semantic bits, illustrating why the
incidence and information statements must be kept distinct.

## 3. A pairwise packing saturating the scale

The parity example's semantic loss is not forced by bounded fan-in. A local
switching code plus an outer code recovers the full order `t sqrt(N)`.

### Lemma BF.3 (support-to-response bound)

If a symmetric hollow matrix `D` has `m` unordered nonzero entries, each of
magnitude two, then

```math
Q(D):=max_(x in {+-1}^N)|sum_(u<v)D_(uv)x_ux_v|
\ge {m\over sqrt(2N)}.                               \tag{BF.13}
```

#### Proof

A vertex bipartition retains at least `m/2` nonzero entries across it. For
uniform signs on the right and optimally chosen signs on the left, the sharp
`p=1` Khintchine inequality gives

```math
sqrt(2)sum_u sqrt(d_u)
\ge {sqrt(2)\over sqrt N}sum_u d_u
\ge {m\over sqrt(2N)}.                               \tag{BF.14}
```

Flipping all left signs reverses the cross term and preserves both within-
part terms, so one of the two full quadratic values is at least this large
in absolute value. `square`

### Theorem BF.4 (quadratic fan-in packing)

There is an absolute `c>0` such that, for all sufficiently large `N` and
every integer `1<=t<=cN`, there exist

* `h=t floor(sqrt(N))` hidden Boolean coordinates;
* exact hollow signings `H_z` on `N` vertices for every `z in {0,1}^h`;
* the fixed continuation language `{-H_w:w in {0,1}^h}`;

such that every edge coefficient depends on at most `t` hidden coordinates,
every hypercube neighbour has response distance at least `cN^(3/2)`, and
there is a subset `C` with

```math
log_2|C|\ge c t sqrt(N),
\qquad d(z,z')\ge cN^(3/2)
\quad(z!=z',\ z,z' in C).                            \tag{BF.15}
```

Thus the `t sqrt(N)` incidence scale and its possible response-information
interpretation are both optimal up to constants for unrestricted quadratic
signings.

Here and below the response is explicitly

```math
R_w(z)=max_x|H_z(x)-H_w(x)|.                         \tag{BF.15a}
```

#### Proof

Put `g=floor(sqrt(N))`, `q=2^t`, and partition all but fewer than `g` edges
of `K_N` into `g` cells of common size

```math
d=floor({binom(N,2)\over g})=Theta(N^(3/2)).          \tag{BF.16}
```

Choose `q` spins `u_a in {+-1}^N` so that every pair has Hamming distance
in `[N/4,3N/4]`. Independent uniform spins and a union bound work for
`t<=cN`, since

```math
Pr\{|d_H(u_a,u_b)-N/2|>N/4\}\le2e^(-N/8).           \tag{BF.17}
```

Choose the edge partition so that, in every cell and for every `a!=b`, at
least `d/4` edges cross the cut `u_au_b`. Conditional on (BF.17), a random
cell has mean at least `3d/8`; its hypergeometric lower-tail probability is
at most `e^(-d/48)`. A union bound over `gq^2` choices succeeds for the same
sufficiently small `c`, because `d=Theta(N^(3/2))`.

View `z` as `g` consecutive `t`-bit symbols `a_1,...,a_g`. On edge `{u,v}`
in cell `j`, define

```math
A_z(u,v)=u_(a_j)(u)u_(a_j)(v).                       \tag{BF.18}
```

Leftover edges get sign `+1`. This lookup has fan-in at most `t`. If one
hidden coordinate changes, only one symbol changes. At least `d/4`
coefficients flip, and the old switching vector makes all their differences
agree. Thus context `-H_z` gives response at least `d/2=Theta(N^(3/2))`.

Take a `q`-ary code of length `g`, relative Hamming distance at least `1/4`,
and size at least `q^(c g)`. Greedy packing suffices: a radius-`g/4` ball
has at most `2^gq^(g/4)` words; for `q=2` use the standard positive-rate
binary Gilbert estimate. Two codewords differ in at least `g/4` cells, so
their signings differ on at least

```math
{g\over4}{d\over4}={gd\over16}=Theta(N^2)            \tag{BF.19}
```

edges. Lemma BF.3 and context `-H_z` give their response separation, and
the code size gives the information bound. `square`

These children need not have cap `O(N^(3/2))`. Simultaneous spectral
flatness is an additional resource constraint, not a consequence of fan-in.

## 4. Gram broadcast is necessarily high-fan-in

Use any basis of `Alt(F_2^r)` in `state_dependent_gram_broadcast.md`. Every
edge evaluation is a linear functional `w_e in F_2^h`, and its coefficient
is

```math
A(e)(-1)^(w_e dot z).                                \tag{BF.20}
```

The fan-in of edge `e` is `|supp(w_e)|`. The sampler has minimum distance
at least `E/4`, where `E=binom(k,2)`. Hence every basis coordinate occurs in
at least `E/4` edge functionals. Double counting yields

```math
sum_e|supp(w_e)|=sum_i #\{e:w_(e,i)=1\}\ge {hE\over4},
\qquad t=max_e|supp(w_e)|\ge {h\over4}.              \tag{BF.21}
```

Since `h=Theta(k)`, every coordinate system for the alternating-form state
has worst-case fan-in `Theta(k)` and average edge fan-in at least `h/4`.
Each coordinate affects at least `E/4=Theta(k^2)` coefficients. The generic
threshold from (BF.6) is only `Theta(k^(3/2))`; Gram broadcast exceeds it by
`Theta(sqrt(k))`. It therefore occupies the genuinely nonlocal regime in
both fan-in and influence. This conclusion is basis-independent because
every basis vector is a nonzero alternating form.

## 5. Research judgment

At the unrestricted bounded-atom level, the law is complete:

* total-scale visibility of `h` coordinates needs dependency incidence
  `Omega(hN^(3/2))`;
* `O(N^2)` unit pair atoms imply `h=O(t sqrt(N))`;
* exact quadratic signings attain both this incidence scale and an
  `exp(Omega(t sqrt(N)))` contextual packing;
* the flat Gram code escapes only through linear fan-in and quadratic
  coordinate influence.

The next nontrivial question is whether simultaneous cap bounds
`Q(H_z)=O(N^(3/2))`, together with a natural small continuation algebra,
force a stronger tradeoff. Any such theorem must use cancellation or
spectral geometry; coefficients, fan-in, and contextual nonexpansiveness
alone are exhausted.

Finite checks appear in
[`../experiments/verify_bounded_fanin_broadcast_law.py`](../experiments/verify_bounded_fanin_broadcast_law.py).
