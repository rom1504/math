# Fixed-rank bridges: exact response algebra and the exact/approximate split

Status: main-agent proof draft for independent audit.  The finite algebra and
error bounds below are elementary.  The metric-entropy paragraph imports
Bronshtein's convex-function entropy theorem and is not used in the algebraic
proof.

## 1. Featured landscapes and declared contexts

A finite `r`-featured landscape is a triple

```math
L=(X,H,\phi),\qquad H:X\to\mathbb R,\quad \phi:X\to\mathbb R^r.
```

Its linear-field response is

```math
V_L(t)=\max_{x\in X}\{H(x)+\langle t,\phi(x)\rangle\}.
                                                            \tag{FR.1}
```

Put `K_L=conv(phi(X))` and define the upper concave roof

```math
\bar H_L(u)=\max\left\{
 \sum_x p_xH(x):p\in\Delta_X,\ \sum_xp_x\phi(x)=u
 \right\}.                                                   \tag{FR.2}
```

Equivalently, the hypograph of `bar H_L` is the downward completion of
`conv{(phi(x),H(x)):x in X}`.  Formula (FR.1) is the restricted support
function

```math
V_L(t)=\max_{u\in K_L}\{\bar H_L(u)+\langle t,u\rangle\}.     \tag{FR.3}
```

The declared future contexts below include every singleton landscape whose
feature is `t`.  Therefore two landscapes are contextually equivalent for
these futures exactly when their functions `V_L` agree.  Fenchel--Moreau
duality for a closed concave function says this is equivalent to equality of
their upper roofs.  Thus the roof is the exact semantic quotient; it is not
merely a sufficient statistic chosen in advance.

## 2. Complete feature coupling

For two landscapes with the same feature dimension, define

```math
L\star K=(X\times Y,H\star G,\phi\star\psi),                  \tag{FR.4}
```

where

```math
(H\star G)(x,y)=H(x)+G(y)+\langle\phi(x),\psi(y)\rangle,
\qquad
(\phi\star\psi)(x,y)=\phi(x)+\psi(y).                        \tag{FR.5}
```

This is the correct operation for blocks with all cross interactions carried
by one common `r`-dimensional feature.

### Theorem FR.1 (exact roof algebra)

The upper roof of `L star K`, and hence every future linear-field response,
is determined by the two input roofs.  Explicitly,

```math
V_{L\star K}(t)
=\max_{u\in K_L,\,v\in K_K}
 \{\bar H_L(u)+\bar G_K(v)+\langle u,v\rangle
                         +\langle t,u+v\rangle\}.             \tag{FR.6}
```

Equivalently, on roof points use

```math
(u,h)\circ(v,k)=(u+v,h+k+\langle u,v\rangle),                 \tag{FR.7}
```

then take the upper concave hull.  The resulting operation on upper roofs is
associative.  For `m` inputs its microscopic formula is

```math
\sum_{i=1}^m H_i(x_i)
+\sum_{1\le i<j\le m}\langle\phi_i(x_i),\phi_j(x_j)\rangle,
\qquad \phi_{\rm out}=\sum_i\phi_i.                           \tag{FR.8}
```

#### Proof

The right side of (FR.6) is separately affine in the probability vector used
to realize `u` and in the one used to realize `v`.  Its maximum over the two
convexified landscapes is consequently attained at original states and is
exactly the definition of `V_(L star K)`.  Formula (FR.7) is the same
calculation before taking the upper hull.  For three roof points, either
parenthesization gives

```math
(u+v+w,
 h+k+l+<u,v>+<u,w>+<v,w>),
```

so the point operation, and therefore upper-hull composition, is
associative.  Induction gives (FR.8). `square`

If `C_L=conv{(phi(x),H(x)):x in X}` is the compact lifted convex body, the
same law can be written without roof notation as

```math
C_{L\star K}=conv\{p\circ q:p\in C_L,\ q\in C_K\}.             \tag{FR.8a}
```

Indeed `circ` is bi-affine, so a product of convex combinations is the
corresponding product convex combination of pure starred points.  The final
convex hull in (FR.8a) is essential: fixing the total feature before this
concavification can give the wrong roof.

### Corollary FR.2 (rank-factorized Boolean bridge)

Let

```math
H(x,y)=H_1(x)+H_2(y)+x^TRy,
\qquad R=UV^T,\quad rank(R)\le r.                              \tag{FR.9}
```

Taking `phi(x)=U^Tx` and `psi(y)=V^Ty` makes (FR.9) an instance of
Theorem FR.1.  Thus a rank-`r` bridge has an exact `r`-dimensional semantic
interface even when the two internal landscapes are arbitrary.  Quantitative
claims must depend on the feature radii of the chosen factorization, not on
rank alone.

This is stronger than the fixed-rank Curie--Weiss benchmark: no finite atom
or exchangeability assumption is imposed inside either block.  It is weaker
than arbitrary dense composition because all future cross-block queries must
factor through the declared common feature.

## 3. Stability and fixed-error compression

Let a roof body be represented by its compact upper graph and use Euclidean
Hausdorff distance.  Suppose

```math
||u||\le R_L,\quad ||v||\le R_K,\quad ||t||\le T.              \tag{FR.10}
```

### Proposition FR.3 (one-composition stability)

If both roof bodies are replaced by Hausdorff-`delta` approximants, their
composed responses on `||t||<=T` differ by at most

```math
delta(R_L+R_K+2T+2)+delta^2.                                  \tag{FR.11}
```

The same conclusion holds for any finite approximants, so the approximating
states need not themselves be realizable Boolean landscapes.

#### Proof

For matched points `(u,h),(u',h')` and `(v,k),(v',k')`, each feature and
height discrepancy is at most `delta`.  The two height terms cost `2delta`,
the field terms cost `2Tdelta`, and

```math
|<u,v>-<u',v'>|
\le delta R_K+(R_L+delta)delta.
```

Match an optimizer in one product to the other product and repeat in the
opposite direction. `square`

There is also a completely explicit quotient, without invoking convex-body
entropy.  Let `Q` map the feature ball to an internal `eta`-net `C`, so that
`||phi(x)-Qphi(x)||<=eta`, and store only

```math
w(c)=\max\{H(x):Q\phi(x)=c\}.                                  \tag{FR.11a}
```

Only nonempty buckets are retained (equivalently, put `w(c)=-infinity` for an
empty bucket).  Then

```math
\widehat V_L(t)=\max_{c\in C}\{w(c)+\langle t,c\rangle\}
```

obeys `|V_L(t)-Vhat_L(t)|<=T eta` for `||t||_*<=T`.  If the
feature radius is `P`, one may take

```math
|C|\le(1+2P/\eta)^r.                                          \tag{FR.11b}
```

Quantizing `w` to height mesh `zeta`, when `|H|<=M`, adds at most `zeta`
response error and needs at most

```math
(1+2P/\eta)^r\log_2(2+2M/\zeta)                               \tag{FR.11c}
```

bits, up to encoding the fixed net itself.  Quantizing both sides of one
bridge changes the cross term by at most

```math
R_K\eta_L+(R_L+\eta_L)\eta_K.                                 \tag{FR.11d}
```

These are explicit feature buckets; their centers need not themselves be
realized features.  This is not a claim that every approximating convex body
is realized by the original constrained landscape class.

For fixed `r`, bounded feature radius, bounded height, and a fixed compact
future-field domain, the response functions `V_L` are uniformly bounded,
convex, and Lipschitz.  Bronshtein's entropy theorem for that function class
gives

```math
\log N(\epsilon)=O_r(\epsilon^{-r/2}).                         \tag{FR.12}
```

For a fixed full-dimensional convex field domain with nonempty interior and
nondegenerate uniform height/slope bounds, the worst-case order is
`Theta_r(epsilon^(-r/2))`.  Lower-dimensional or singleton future domains can
have smaller entropy.  Finite max-affine responses are dense in the
full-dimensional class, so restricting to finite landscapes does not change
that exponent.  This does not assert uniform sup-norm control of roofs whose
attainable feature domains vary.  An elementary grid argument, with worse
constants/exponent, already gives a finite approximate state whose
description is independent of `|X|`.

Consequently fixed bridge rank gives uniform **approximate** contextual
compression at each declared scale.  When rank grows, even the elementary
grid description has size exponential in `r log(C/epsilon)`.  A sufficient
subextensivity regime for the number of stored cells of an `n`-spin block is
therefore

```math
r\log(C/\epsilon_n)=o(\log n),                                \tag{FR.13}
```

At constant error, `r=O(log n)` still gives a polynomial-size state, and
`r=o(log n)` gives a subpolynomial one.  The crossover with brute-force
`2^n` storage is only at `r=Theta(n)`.  Formula (FR.13) is an upper-bound
regime, not a claimed sharp threshold; if `epsilon_n=n^{-a}`, even polynomial
size from this grid estimate requires fixed `r`.

## 4. Rank one can still have exponentially many exact response atoms

### Proposition FR.4 (exact rank-one exposure)

For every `n` there is a quadratic Boolean landscape on `{-1,1}^n` with one
linear scalar feature for which all `2^n` configurations are distinct exposed
states of the exact upper roof.

#### Proof

Put

```math
a=(1,2,4,...,2^{n-1}),\qquad
u(x)={a^Tx\over 2^n-1},\qquad H(x)=-u(x)^2.                    \tag{FR.14}
```

The `2^n` feature values are the equally spaced grid from `-1` to `1`.
For a field `t`,

```math
H(x)+tu(x)=-(u(x)-t/2)^2+t^2/4.                               \tag{FR.15}
```

Choosing `t=2u(x)` exposes that configuration uniquely.  Also `H` is a
constant plus a homogeneous quadratic Boolean form.  Coupling two such
blocks by `u(x)v(y)` is a rank-one bilinear bridge. `square`

Thus rank one does not imply a small exact contextual quotient or a bounded
number of roof facets.  The exposure margins in (FR.15) are exponentially
small, which is exactly why this does not contradict fixed-error convex
compression.  Any theory of growing bridge rank must specify its response
scale and composition depth; exact state count and macroscopic
rate--distortion are different resources.

### Proposition FR.5 (growing-rank fixed-error lower bound)

For every fixed `0<rho<1/2`, there are `r`-dimensional feature landscapes
with feature radius one whose response class has

```math
2^{2^{Omega_rho(r)}}
```

members separated by a positive constant in uniform response distance.
Hence every fixed-error quotient for the unrestricted class needs
`2^{Omega(r)}` bits.  Together with (FR.11b)--(FR.11c), the constant-error
polynomial-space transition occurs at `r=Theta(log n)`, up to constants.

#### Proof

Take a code `C subset {-1,1}^r` of relative distance `rho` and size
`2^{Omega_rho(r)}`.  Use feature `p_x=x/sqrt(r)`, and for each
`b in {0,1}^C` set

```math
H_b(x)=a b_x\quad(x\in C),\qquad H_b(x)=-1\quad(x\notin C),    \tag{FR.16}
```

where `0<a<2rho`.  At the actual Boolean query `q_c=c/sqrt(r)`, the state
`c` has value `1+ab_c`; another codeword has value at most
`1-2rho+a<1`, and a noncode state has value at most zero.  Therefore

```math
V_{H_b}(q_c)=1+ab_c.                                           \tag{FR.17}
```

The `2^|C|` response functions are pairwise `a`-separated.  The usual greedy
Hamming packing supplies `|C|=2^{Omega_rho(r)}`. `square`

This lower bound uses arbitrary bounded internal landscapes.  It calibrates
the general response theory; it is not a lower bound for quadratic internal
Hamiltonians.

## 5. Director assessment and next falsifiers

This algebra is a genuine near-original benchmark: it derives the state from
declared cross-block queries, composes arbitrary internal landscapes, and
separates exact from approximate complexity.  It does **not** reconnect the
theory to the arbitrary dense signing problem because a general dense bridge
has effective rank of order the block size and no bounded feature radius at
the natural normalization.

The next discriminating questions are:

1. Can normalized balanced composition contract roof approximation error so
   that one fixed-scale convex state remains valid through unbounded depth?
2. Can exposed-face dimension or Gaussian width replace ambient `r` in
   (FR.12) for a declared family of bridge queries?
3. Does any algebraically structured full-rank bridge admit a strict quotient
   of this roof algebra, or can exposed rank-one slices already prove an
   extensive lower bound?
