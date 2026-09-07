# Wave 55: nonnegative core mixtures and the square-root threshold wall

## Status

This memo proves a sharp scoped no-go theorem and a sharper necessary partner
capacity theorem for the self-loop-free identity (10.1292).  It does **not**
prove spectral excess for the favorable families of an exact minimizing
signing, the box witness, the restriction estimate, or convergence.

The main conclusions are:

1. Every normalized nonnegative generating polynomial in the intersection
   size is exactly a convex mixture of the common-core kernels.  Its spectrum,
   diagonal, weighted excess, and spectral-extraction denominator all mix
   linearly.  Consequently a mixture never gives a stronger extraction bound
   than its best component core scale.
2. For a genuinely linear core, the exact weighted histogram can carry
   constant load only through exponentially many partners at one sharp
   entropy rate.  This strengthens the crude maximum-weight necessary bound
   (10.1293).
3. The one-threshold lemma (10.1294) is in fact asymptotically impossible even
   on the entropy-tuned surface (10.1307): after the exponential rates match,
   replacing all pair weights by the threshold weight loses a factor
   `Theta(sqrt(n))`.  Uniformly over compact linear-core ranges, the largest
   load certifiable by any one threshold is `O(n^(-1/2))`, whereas the required
   second-eigenvalue load is bounded below by a positive constant.
4. The same square-root wall holds for every nonnegative mixture supported on
   such linear cores.  The exact multilevel histogram remains viable because
   it sums `Theta(sqrt(n))` typical intersection levels.  Thus this result
   retires threshold compression, not the weighted identity itself.

All finite identities and the tuned asymptotic table are reproduced by
`tmp/weighted_histogram_r55_check.py`; its saved output is
`tmp/weighted_histogram_r55_check.out`.

## 1. Exact normalization of every nonnegative intersection polynomial

Let `Omega` be the `m`-slice and put

```math
d_\ell=\binom m\ell\binom{n-\ell}{m-\ell}.
```

The common-core transition kernel is

```math
K_\ell(S,T)=\frac{\binom{|S\cap T|}{\ell}}{d_\ell}.
\tag{R55H.1}
```

Indeed, for every fixed `S`, double counting a core contained in both sets
gives

```math
\sum_{T\in\Omega}\binom{|S\cap T|}{\ell}=d_\ell.
\tag{R55H.2}
```

Now take an arbitrary nonzero polynomial in the binomial basis with
nonnegative coefficients,

```math
Q(j)=\sum_{\ell=0}^{m-1}c_\ell\binom j\ell,
\qquad c_\ell\ge0,
```

and row-normalize it.  If

```math
Z_Q=\sum_\ell c_\ell d_\ell,
\qquad w_\ell=\frac{c_\ell d_\ell}{Z_Q},
```

then the exact identity is

```math
\boxed{
K_Q(S,T):=\frac{Q(|S\cap T|)}{Z_Q}
=\sum_\ell w_\ell K_\ell(S,T),
\qquad w_\ell\ge0,\quad\sum_\ell w_\ell=1.}
\tag{R55H.3}
```

Thus a nonnegative generating polynomial supplies no kernels beyond convex
mixtures of common-core kernels.  Conversely every such mixture has the form
(R55H.3).

On Johnson harmonic degree `j`, its spectrum is therefore

```math
\boxed{
\Lambda_j(Q)=\sum_\ell w_\ell\lambda_j(\ell),
\qquad
\lambda_j(\ell)=
\frac{(\ell)_j(n-m)_j}{(m)_j(n-\ell)_j},}
\tag{R55H.4}
```

with `lambda_j(ell)=0` for `j>ell`.  Its self-loop is

```math
\boxed{
h_Q=K_Q(S,S)=\sum_\ell w_\ell h_\ell,
\qquad
h_\ell=\frac{\binom m\ell}{d_\ell}
=\binom{n-\ell}{m-\ell}^{-1}.}
\tag{R55H.5}
```

No self-loop is hidden in the off-diagonal formulas below.

## 2. Exact excess and extraction convexity: mixtures cannot beat a scale

Retain the center law and the `r_z^2`-biased weights of (10.1292).  By
linearity, write

```math
P_Q=\sum_\ell w_\ell P_\ell,
\qquad
E_\ell=P_\ell-\lambda_2(\ell).
```

Then

```math
\boxed{P_Q-\Lambda_2(Q)=\sum_\ell w_\ell E_\ell.}
\tag{R55H.6}
```

In particular, if a mixture has saved excess at least `epsilon`, some
component in its support has excess at least `epsilon`.  This observation is
exact and does not depend on the number of component scales or on their
weights.

The sharper extraction comparison also closes exactly.  Put

```math
D_\ell=(1-\lambda_1(\ell))
+n(\lambda_1(\ell)-\lambda_2(\ell)),
```

and define `D_Q` from the mixture eigenvalues in the same way.  Then

```math
\boxed{D_Q=\sum_\ell w_\ell D_\ell.}
\tag{R55H.7}
```

The proof of the spectral-excess extraction theorem (10.1261) applies to the
positive semidefinite mixture and gives

```math
M\ge\frac{[P_Q-\Lambda_2(Q)]_+}{D_Q}.
```

Whenever the numerator is positive,

```math
\boxed{
\frac{P_Q-\Lambda_2(Q)}{D_Q}
=\sum_\ell\frac{w_\ell D_\ell}{D_Q}
\frac{E_\ell}{D_\ell}
\le\max_{\ell:w_\ell>0}\frac{E_\ell}{D_\ell}.}
\tag{R55H.8}
```

Hence nonnegative scale mixing cannot improve the degree lower bound over
the best component.  It may be a convenient way to *prove that some scale
works*, but it is not a weaker final proof obligation than existential
single-scale weighted excess.  This is the exact analogue, at the current
second-eigenvalue target, of the earlier mixture walls (10.1196) and
(10.1263).

## 3. Exact maximum load for a prescribed number of distinct partners

For a fixed base `S`, the number of partners at intersection `j` is

```math
A_j=\binom mj\binom{n-m}{m-j}.
\tag{R55H.9}
```

Write

```math
q_Q(j)=K_Q(S,T)\quad\text{when }|S\cap T|=j.
```

Every binomial coefficient `binom(j,ell)` is nondecreasing in `j`, so
`q_Q(j)` is nondecreasing.  It follows by a direct exchange argument that,
among all sets `A` of `r` distinct partners `T != S`, the exact maximum of

```math
\sum_{T\in A}q_Q(|S\cap T|)
```

is obtained by filling the largest-intersection shells first.  More
explicitly, if `s` is chosen so that

```math
\sum_{j>s}A_j\le r\le\sum_{j\ge s}A_j,
```

then the exact capacity is

```math
\boxed{
\mathcal C_Q(r)=
\sum_{j>s}A_jq_Q(j)
+\left(r-\sum_{j>s}A_j\right)q_Q(s).}
\tag{R55H.10}
```

The term `j=m` is omitted throughout, so (R55H.10) contains no self-loop.
This is the sharp law-free partner-load optimization for every nonnegative
generating polynomial.

## 4. Sharp exponential capacity of one linear common-core kernel

Assume

```math
\frac mn\to p,\qquad \frac\ell n\to\alpha,
\qquad 0<\alpha<p<1,
```

with the parameters in a fixed compact subset of this open region.  Under
`K_ell(S,.)`, first choose a uniform `ell`-core of `S` and then complete it.
Consequently

```math
J-\ell\sim
\operatorname{Hypergeom}(n-\ell,m-\ell,m-\ell),
\qquad J=|S\cap T|.
\tag{R55H.11}
```

Its normalized mean is

```math
\boxed{
\beta_*(p,\alpha)
=\alpha+\frac{(p-\alpha)^2}{1-\alpha}
=\frac{p^2+\alpha(1-2p)}{1-\alpha}.}
\tag{R55H.12}
```

Equivalently,

```math
\beta_*-p^2=\frac{\alpha(1-p)^2}{1-\alpha}>0,
\qquad
p-\beta_*=\frac{(p-\alpha)(1-p)}{1-\alpha}.
```

Let

```math
v_*(p,\alpha)=
pH\!\left(\frac{p-\beta_*}{p}\right)
+(1-p)H\!\left(\frac{p-\beta_*}{1-p}\right).
\tag{R55H.13}
```

This is the entropy rate of the Johnson shell at the transition-typical
intersection.  Hypergeometric concentration and uniform Stirling estimates
give the following sharp **partner-capacity theorem**: for every fixed
`eta>0`, uniformly on compact parameter sets,

```math
\boxed{
\begin{aligned}
r\le e^{n(v_*-\eta)}&\quad\Longrightarrow\quad
\mathcal C_\ell(r)=e^{-\Omega_\eta(n)},\\
\text{some }r\le e^{n(v_*+\eta)}&\quad\Longrightarrow\quad
\mathcal C_\ell(r)=1-o(1).
\end{aligned}}
\tag{R55H.14}
```

For the first line, discard the exponentially small event
`|J/n-beta_*|>delta`, and note that every point in the remaining typical
window has probability at most `exp{-n(v_*-eta/2)}` after choosing `delta`
small enough.  The second line takes all partners in that typical window;
their number is at most `exp{n(v_*+eta)}` and their transition mass tends to
one.  Since (R55H.10) is the maximal load, these two bounds are sharp at the
exponential scale.

Also

```math
\lambda_1(\ell)\to\theta,
\qquad \lambda_2(\ell)\to\theta^2,
\qquad
\theta=\frac{\alpha(1-p)}{p(1-\alpha)}\in(0,1).
\tag{R55H.15}
```

Thus any fixed positive fraction of the second-eigenvalue load requires
`exp{n v_*+o(n)}` distinct partners.  This is sharper than bounding every
pair by its near-diagonal maximum, as in (10.1293).

There is an aggregate consequence directly for (10.1292).  Draw a center
with the `omega_C` law and then a uniform base `S` in its favorable family,
and define its self-loop-free load

```math
L(z,S)=\sum_{\substack{T\in F_z\\T\ne S}}K_\ell(S,T).
```

If `P_ell-lambda_2(ell)>=epsilon_n>=0`, then

```math
\mathbb E L=P_\ell-h_\ell
\ge\lambda_2(\ell)-h_\ell+\epsilon_n=\theta^2+o(1).
```

Since `0<=L<=1`, a fixed positive fraction of these base incidences have,
say, `L>=theta^2/4`.  By (R55H.14), every such base has

```math
\boxed{|F_z\setminus\{S\}|\ge\exp\{nv_*(p,\alpha)-o(n)\}.}
\tag{R55H.16}
```

This remains a necessary condition, not a signing theorem.  Unlike the old
one-threshold condition, it does not require all of those partners to lie on
one side of a fixed overlap threshold.

## 5. The tuned one-threshold surface still loses `sqrt(n)`

For `s<=m-1`, let the entire distinct Johnson ball be

```math
\mathcal B_{n,m}(s)=\sum_{j=s}^{m-1}A_j.
```

No favorable family can have more than this many partners of intersection at
least `s`.  Multiplying the maximal possible partner count by the per-pair
weight used in (10.1294) gives

```math
\mathcal B_{n,m}(s)K_\ell(s)
=\mathcal B_{n,m}(s)\frac{\binom s\ell}{d_\ell}.
\tag{R55H.17}
```

The prior entropy comparison found equality only when

```math
\alpha=\frac{\beta-p^2}{1+\beta-2p}.
```

Solving for `beta` shows that this equality surface is exactly

```math
\boxed{\beta=\beta_*(p,\alpha).}
\tag{R55H.18}
```

So the unique entropy-tuned threshold is simply the mean intersection of the
`K_ell` transition.  Exponential equality is not enough.  Uniformly for
linear cores in a compact subset of `0<alpha<p<1`, one has the stronger
bound

```math
\boxed{
\sup_{\ell\le s\le m-1}
\mathcal B_{n,m}(s)\frac{\binom s\ell}{d_\ell}
=O(n^{-1/2}).}
\tag{R55H.19}
```

Here is a short proof.  The gap `beta_*-p^2` is uniformly positive.  Below a
fixed point strictly between these two values, the entropy comparison is
uniformly strict, so (R55H.17) is exponentially small.  Above that point,
the Johnson shell sizes decrease at a uniform geometric rate, because

```math
\frac{A_{j+1}}{A_j}
=\frac{(m-j)^2}{(j+1)(n-2m+j+1)}<1-c.
```

Hence `B(s)<=C A_s`, and

```math
\mathcal B_{n,m}(s)K_\ell(s)
\le C A_sK_\ell(s)
=C\Pr_{K_\ell}\{J=s\}.
```

The hypergeometric law (R55H.11) has all four proportions bounded away from
zero and one and variance `Theta(n)`.  Uniform Stirling bounds give
`sup_s Pr{J=s}=O(n^(-1/2))`, proving (R55H.19).  At the tuned mean, the bound
has the right order.  Conceptually, the exact histogram sums
`Theta(sqrt(n))` typical levels, each carrying `Theta(n^(-1/2))` transition
mass; replacing all their weights by the boundary weight retains only one
level's scale.

Since

```math
\frac{b_\ell-\binom m\ell}{d_\ell}
=\lambda_2(\ell)-h_\ell\longrightarrow\theta^2>0,
```

(R55H.19) proves the sharpened **one-threshold no-go theorem**:

```math
\boxed{
\mathcal B_{n,m}(s)<
\frac{b_\ell-\binom m\ell+d_\ell\epsilon_n}
{\binom s\ell}
\quad\text{for every }s\ge\ell}
\tag{R55H.20}
```

for all sufficiently large `n` and every `epsilon_n>=0`.  Thus (10.1294)
is asymptotically impossible for *every* fixed linear `(p,alpha,beta)`,
including the tuned surface.  Away from (R55H.18) the deficit is exponential;
on it the remaining deficit is of order `sqrt(n)`.

For the exact tuned example

```math
p=\frac35,\qquad \alpha=\frac9{25},\qquad \beta_*=\frac9{20},
```

the checker finds the maximizing threshold exactly at `s=beta_* n` for all
tested multiples.  The quantity in (R55H.19) times `sqrt(n)` tends
numerically to about `2.66`, while the required baseline tends to
`theta^2=9/64=0.140625`.  The full-slice threshold already fails from
`n=400` onward in the tested sequence; asymptotic failure follows from
(R55H.19), not from the numerics.

## 6. No nonnegative linear-core mixture repairs threshold compression

Let a possibly `n`-dependent mixture be supported on linear core sizes in a
fixed compact subset of `0<alpha<p<1`, and write

```math
q_Q(s)=\sum_\ell w_\ell K_\ell(s).
```

The constants in (R55H.19) are uniform, so

```math
\boxed{
\sup_s\mathcal B_{n,m}(s)q_Q(s)
\le\sum_\ell w_\ell
\sup_s\mathcal B_{n,m}(s)K_\ell(s)
=O(n^{-1/2}).}
\tag{R55H.21}
```

Meanwhile `Lambda_2(Q)` is bounded below by a positive constant and `h_Q`
is exponentially small.  Therefore every mixture analogue that selects one
overlap threshold and replaces all contributing weights by `q_Q(s)` is
impossible, even for the whole slice.  Allowing a separate threshold for
each component still gives a weighted sum of `O(n^(-1/2))` terms and has the
same wall.

The compact linear-core hypothesis matters.  Near-identity cores with
`m-ell=o(n)` can have a transition supported on a much thinner range, and
the local anti-concentration scale need not be `n^(-1/2)`.  This memo makes
no no-go claim for that regime.  The exact convexity statements
(R55H.3)--(R55H.8) hold at every scale.

## 7. Research judgment

The exact weighted identity (10.1292) is genuinely weaker than its
one-threshold sufficient condition.  For example, on the full slice its
self-loop-free load is exactly `1-h_Q`, so it exceeds
`Lambda_2(Q)-h_Q`; nevertheless every one-threshold lower bound in the
linear-core regime tends to zero by (R55H.21).  This separation is not a
renaming: it is the missing sum over `Theta(sqrt(n))` transition-typical
intersection levels.

At the same time, nonnegative mixtures do not create a new spectral route.
By (R55H.6)--(R55H.8), any successful mixture already contains a component
core with at least as strong an exact excess/extraction certificate.  The
clean remaining overlap target is therefore:

- work at one or several genuinely non-direct core scales;
- retain the complete self-loop-free weighted histogram, with no minimum-
  weight threshold replacement; and
- use exact-minimizer/controlled-row structure to prove that its weighted
  load exceeds the second-eigenvalue baseline by the saved amount.

For fixed linear cores, any such proof must also create, on a positive
fraction of the biased base incidences, the sharp partner capacity
`exp{n v_*+o(n)}` from (R55H.16).  The full typical-shell rate is feasible as
an abstract Johnson histogram, so this is not a counting contradiction; the
unproved content is why exact minimizing signings should produce it at one
controlled project-row center class.
