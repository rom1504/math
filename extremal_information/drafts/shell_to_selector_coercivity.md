# From an augmented-cut shell law to selector coercivity

**Status:** rigorous scoped theorem/falsifier draft.  No claim about arbitrary
near-minimizers is made.  The example below respects the augmented-cut
geometry, has an exact active shell, vanishing normalized edge
bias, and a logarithmic port algebra.  It shows that the conclusion of FB.1
alone cannot imply same-sign selector coercivity: an orientation-coherence
statistic is necessary.

Throughout this note `A` is a symmetric hollow matrix and

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j={1\over2}x^TAx.
```

## 1. The exact inequality that a shell-to-selector theorem must use

Let `r>=||A||_(2->2)`, let `w_1,...,w_p in {+-1}^N`, and let `tau` be an
odd Boolean selector.  Write

```math
tau(t)=\sum_{S\in\mathcal F}\widehat\tau(S)\prod_{j\in S}t_j,
\qquad |S|\text{ odd for }S\in\mathcal F,
```

and put

```math
w_S=\bigodot_{j\in S}w_j,
\qquad c_S^epsilon=\widehat\tau(S)\prod_{j\in S}\epsilon_j,
\qquad y^epsilon=\tau(\epsilon_1w_1,...,\epsilon_pw_p)
                 =\sum_Sc_S^epsilon w_S.                 \tag{SC.1}
```

With `Z=(w_S)_(S in mathcal F)`, set

```math
G={Z^TZ\over N},\qquad R={Z^TAZ\over rN},
\qquad D_rho=G-rho R\succeq0\quad(rho\in\{+-1\}).          \tag{SC.2}
```

The selector defect has the exact identity

```math
\boxed{
1-{rho(y^epsilon)^TAy^epsilon\over rN}
=(c^epsilon)^TD_rho c^epsilon.}                           \tag{SC.3}
```

This is the minimal algebraic bridge from poles to a selector.  If a
Fourier product `w_S` occurs in an augmented near-top shell with orientation
`sigma_S`, then

```math
sigma_S w_S^TAw_S\ge 2Q(A)-2uN^{3/2}.                     \tag{SC.4}
```

Thus (SC.4) makes the **diagonal entry of `D_(sigma_S)`** small.  It does
not necessarily make the diagonal entry of one common `D_rho` small.  This
distinction cannot be repaired by multiplying `w_S` by a global sign,
because its quadratic energy is unchanged.

For completeness, suppose a set `C subseteq mathcal F` has the common
orientation `rho` and `(D_rho)_(SS)<=d` on `C`.  Since `D_rho` is positive
semidefinite, Cauchy--Schwarz entrywise gives the sharp general estimate

```math
\sqrt{(c^epsilon)^TD_rho c^epsilon}
\le \sqrt d\sum_{S\in C}|\widehat\tau(S)|
   +\sqrt2\sum_{S\notin C}|\widehat\tau(S)|.               \tag{SC.5}
```

Here `(D_rho)_(SS)<=2` was used on the uncontrolled products.  If every
Fourier product has the same orientation and defect at most `d`, this
reduces to

```math
(c^epsilon)^TD_rho c^epsilon
\le d||\widehat\tau||_1^2.                                \tag{SC.6}
```

Equations (SC.3)--(SC.6) isolate two obligations which FB.1 does not state:

1. the products in the Fourier support, not merely the original sampled
   ports, must remain near the shell;
2. those products must have one coherent orientation (or the Fourier mass
   on the other orientation must be negligible).

The first-marginal quantity in FB.1,

```math
{1\over {N\choose2}}\sum_{i<j}
  \left|\mathbb E[A_{ij}\sigma w_iw_j]\right|,             \tag{SC.7}
```

contains neither piece of information.  The next theorem shows that this is
not just a gap in (SC.5).

## 2. A Walsh shell with two incompatible low-bias laws

Let `m>=2`, `N=2^{2m}`, and index coordinates by
`x=(s,t) in F_2^m times F_2^m`.  Let

```math
K_(x,y)=(-1)^{x dot y},\qquad
A=K-diag(K),\qquad
f_u(s,t)=(-1)^{s dot t+u dot s+u dot t}           \tag{SC.8}
```

for `u in F_2^m`, and put `sigma_u=(-1)^{|u|}`.  The elementary Walsh
calculation

```math
Kf_u=sigma_u\sqrt N f_u                                \tag{SC.9}
```

follows by first summing over `t`.  Moreover `tr(K)=0`, so

```math
f_u^TAf_u=sigma_uN^{3/2}.                            \tag{SC.10}
```

For every Boolean `x`, diagonal deletion contributes
`x^Tdiag(K)x=tr(K)=0`.  Since `||K||=sqrt N`, (SC.10) therefore gives

```math
\boxed{Q(A)={N^{3/2}\over2}.}                       \tag{SC.11}
```

Consequently every augmented cut

```math
z^(u)_(xy)=sigma_u f_u(x)f_u(y)                    \tag{SC.12}
```

is exactly active, hence belongs to `S_0(A)`.  In particular this is stronger
than an eta-thick shell at the correct `N^{3/2}` scale, not an arbitrary
collection of Boolean vectors.

### Theorem SC.1 (augmented shell balance does not choose a coherent side)

Let `mu_all` be the uniform law on the augmented cuts (SC.12).  Then

```math
{1\over {N\choose2}}\sum_{x<y}
 \left|\mathbb E_(u\sim mu_all)
 [A_(xy)sigma_uf_u(x)f_u(y)]\right|
={\sqrt N\over N-1}=O(N^{-1/2}).                  \tag{SC.13}
```

This is the smallest possible normalized `l_1` bias for any law supported
on these exact active atoms.  Indeed every such law has
`sum_(x<y)m_(xy)=Q(A)=N^(3/2)/2`, hence
`E^(-1)sum|m_(xy)|>=sqrt N/(N-1)`.  The law `mu_all` attains equality
(in fact each of its nonzero edge means is `+1`).  Thus strengthening FB.1
to optimal first-moment balance would still not remove the example.

Nevertheless there are three ports from its support whose majority
selector has unit positive and negative quadratic defect.  Namely, for two
distinct coordinate vectors `e_1,e_2 in F_2^m`, use

```math
w_1=f_0,\qquad w_2=f_(e_1),\qquad w_3=f_(e_2),
\qquad y=Maj(w_1,w_2,w_3).                        \tag{SC.14}
```

Then

```math
y^TAy=0,
\qquad
1-{rho y^TAy\over(\sqrt N+1)N}=1
\quad\text{for both }rho=+-1.                    \tag{SC.15}
```

Thus a shell law can have asymptotically vanishing FB.6 bias while a
constant-size (hence `O(log N)`) selector built from shell atoms has constant
quadratic defect.

#### Proof

For `x=(s,t)` and `y=(s',t')`, the expectation in (SC.13) vanishes unless

```math
s+t+s'+t'=1^m.                                    \tag{SC.16}
```

For every `x` there are exactly `sqrt N` choices of `y`, none equal to
`x`, satisfying (SC.16).  This gives `Nsqrt N/2` exceptional unordered
edges out of `N(N-1)/2`, proving (SC.13).

The four Walsh characters

```math
f_0, f_(e_1), f_(e_2), f_(e_1+e_2)
```

are mutually orthogonal, and the majority identity gives

```math
y={f_0+f_(e_1)+f_(e_2)-f_(e_1+e_2)\over2}.       \tag{SC.17}
```

Their eigenvalue signs in (SC.9) are respectively `+,-,-,+`.  Hence
`y^TKy=0`.  Since `y` is Boolean,
`y^Tdiag(K)y=tr(K)=0`, proving (SC.15). `square`

There is an exact formula behind this example.  Suppose the subset sums
`u_S=sum_(j in S)u_j` are distinct on the Fourier support of an odd
selector (as they are in (SC.14)).  Orthogonality of the modulated bent
vectors gives

```math
{y^TAy\over N^{3/2}}
=\sum_{S\in\mathcal F}\widehat\tau(S)^2
  \prod_{j\in S}\sigma_(u_j)
=2^{-p}\sum_{xi\in\{+-1\}^p}
  \tau(xi)\tau((sigma_(u_j)xi_j)_(j=1)^p).        \tag{SC.17a}
```

Thus the missing datum in this exact shell is the selector's Fourier-noise
response to the **orientation word** `(sigma_(u_1),...,sigma_(u_p))`.
Edgewise shell balance does not control that word.  This is also a strict
compressed statistic in this model--`p` orientation bits, not the full
Boolean landscape--so the obstruction identifies a plausible augmentation
rather than merely demanding full parent maximization.

The constant three-port witness amplifies to a genuinely growing port
language.  Take an odd `m=2k+1` and the basis

```math
u_1=e_1,...,u_k=e_k,
\qquad u_(k+1)=e_1+e_(k+1),...,u_m=e_1+e_m.       \tag{SC.17b}
```

Exactly `k` basis elements have negative orientation and `k+1` have
positive orientation.  Use the `p=m=(log_2N)/2` ports `f_(u_j)` and the
`p`-bit majority selector.  All subset sums are distinct, so (SC.17a)
applies.  If `X` and `Y` are the sums of the `k+1` unflipped and `k`
flipped independent Rademacher coordinates, respectively, its right side
is

```math
E[sign(X+Y)sign(X-Y)].                            \tag{SC.17c}
```

In fact this expectation is elementary and exact.  Write the `k+1`-step
sum as `U+epsilon`, where `U` is a `k`-step Rademacher sum independent of
`epsilon`, and let `V` be an independent copy of `U`.  Symmetrize in
`U,V`.  If `|U|` and `|V|` differ, their difference is at least two, and
the two swapped signs cancel.  If `|U|=|V|>0`, averaging `epsilon=+-1`
again cancels.  Only `U=V=0` contributes.  Consequently

```math
{y^TAy\over N^{3/2}}
=1_({k\text{ even}})
 \left(2^{-k}{k\choose k/2}\right)^2
=O(m^{-1}).                                      \tag{SC.17d}
```

Consequently this exact active shell contains a nonvacuous
`Theta(log N)` port frame whose majority selector has absolute quadratic
defect `1-O((log N)^(-1))`; on the infinite subsequence `k` odd the
selector energy is exactly zero.

The obstruction is not the archived generic Gram collision.  Every vector
is an explicit quadratic Walsh/bent vector, every signed edge word is an
augmented cut, and every atom is exactly active.  The
failure is specifically cancellation between the positive and negative
Walsh eigenspaces, which the absolute near-top shell regards as equally
good.

### Theorem SC.2 (orientation-pure multiplicative closure is sufficient)

Let

```math
U_0=\{u in F_2^m:|u|=0\pmod2\}
```

and let `mu_+` be uniform on `(+,f_u)`, `u in U_0`.  It is supported on the
same exact active shell and obeys

```math
{1\over {N\choose2}}\sum_{x<y}
 \left|\mathbb E_(u\sim mu_+)[A_(xy)f_u(x)f_u(y)]\right|
={2\sqrt N-1\over N-1}=O(N^{-1/2}).              \tag{SC.18}
```

Choose any ports from this law and any odd Boolean selector `tau`.  Every
Fourier product `w_S`, `|S|` odd, is another
`f_u` with `u in U_0`.  Therefore it lies in the positive `sqrt N`
eigenspace of `K`, and so does their linear combination
`tau(w_1,...,w_p)`.  As the latter is Boolean,

```math
tau(w_1,...,w_p)^TA tau(w_1,...,w_p)=N^{3/2}.      \tag{SC.19}
```

In particular, with the valid roof `r=sqrt N+1`, every such joint selector
has defect exactly

```math
1-{\sqrt N\over\sqrt N+1}=O(N^{-1/2}).           \tag{SC.20}
```

For an explicitly nonvacuous logarithmic language, take `f_0` together
with `f_(u_1),...,f_(u_(m-1))`, where the `u_i` form a basis of the
even-parity subspace.  The projective row patterns run through the full
`2^(m-1)` character group, while `p=m=(log_2N)/2`.  Thus this gives a
nonvacuous `p=Theta(log N)` selector language whose joint
coercivity follows from an explicit algebraic statistic strictly stronger
than edge balance: an orientation-pure odd-product algebra.  Equation
(SC.18) follows as in (SC.13), now using that the orthogonal complement of
the even-parity subgroup is `{0,1^m}`.  There are `2sqrt N-1` exceptional
off-diagonal choices of `y` for each `x`. `square`

## 3. Consequence for the smallest missing lemma

Theorems SC.1 and SC.2 place the boundary sharply:

```text
eta-thick augmented-shell law + vanishing edge means
    does not imply joint selector coercivity;

eta-thick shell + common orientation + odd-product closure
    does imply logarithmic-port selector coercivity
    in the Walsh model.
```

Therefore `L_sync` cannot be proved from FB.1 by a generic convexity or
first-moment argument.  At least one statistic detecting **oriented
Fourier-product closure** (or a substitute controlling the same coherent
quadratic form) must be forced by near-minimality.  A precise next candidate
is the following strictly stronger, still finite statement:

> For one declared `p=O(log n)` port frame sampled/generated from the FB.1
> shell, almost all Fourier mass of its selector is carried by products
> lying in one oriented thick shell, with the remaining weighted mass small
> in the PSD metric (SC.3), not merely in counting measure.

This note does **not** prove or disprove that candidate for genuine
`epsilon_n ->0` near-minimizers.  The Walsh signing has cap at the conjectural
`1/2` scale, but current rigorous knowledge does not certify it to be an
`o(1)`-near-minimizer relative to `M_N`.  Thus the result is Level 4 rather
than Level 5.  What it rigorously falsifies is the arrow from the conclusion
of FB.1 alone to `L_sync`; what remains available is an essentially
near-minimizer-specific theorem forcing orientation-pure product coherence.
