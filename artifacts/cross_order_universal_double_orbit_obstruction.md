# The universal double is an induced-clique orbit pressure

Status: **exact all-order recurrence, exact actual-minimizer falsifier for
the zero-error tensor inequality, and a scalable linear floor for every
proof that uses only scalar pressure minimality and the Frobenius
normalization**.  The scalable floor is a method-class obstruction in a
weighted block model, not a lower bound on the true signing defect.

For a hollow symmetric signing `A` of order `r`, write

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
\phi_A(u)=\log\mathbb E_x\cosh(uH_A(x)),
```

and

```math
P_r(\beta)=\min_A\phi_A(\beta/\sqrt r),
\qquad
E_{r,r}(\beta)=P_{2r}(\beta)-2P_r(\beta).
```

The zero-matching core and the integral universal double are

```math
K_0(A)=\begin{pmatrix}A&A\\A&-A\end{pmatrix},
\qquad
K_\eta(A)=
\begin{pmatrix}A&A+\eta I\\A+\eta I&-A\end{pmatrix},
\quad \eta\in\{+1,-1\}.                              \tag{1}
```

The `r` entries joining the two copies of the same vertex are zero in
`K_0`; `K_eta` fills them and is a valid order-`2r` signing.

## 1. Exact orbit formula

For `J subseteq [r]`, let `A^J` be obtained by negating exactly the edges
with both endpoints in `J`:

```math
a^J_{ij}=a_{ij}(-1)^{\mathbf1_{\{i,j\subseteq J\}}}. \tag{2}
```

**Lemma 1 (configuration identity).**  For Boolean `x,y`, put
`p_i=x_i y_i` and `J={i:p_i=-1}`.  Then

```math
\boxed{
H_{K_\eta(A)}(x,y)
=2H_{A^J}(x)+\eta(r-2|J|).}                          \tag{3}
```

**Proof.**  Since `y_i=p_i x_i`, the coefficient of
`a_ij x_i x_j` in the three nonmatching blocks is

```math
1-p_ip_j+p_i+p_j
=2(-1)^{\mathbf1_{\{i,j\subseteq J\}}}.
```

The matching edges contribute `eta sum_i p_i=eta(r-2|J|)`. `square`

The energy law of either matrix in (1) is symmetric: the bijection
`(x,y) -> (y,-x)` negates the energy.  Its normalized cosh partition is
therefore its normalized exponential partition.

Let `D_J` have diagonal `-1` on `J` and `+1` off `J`.  Directly from (2),

```math
A^{J^c}=-D_JA^JD_J.                                  \tag{4}
```

Consequently the positive exponential partition of `A^{J^c}` is the
negative exponential partition of `A^J`.  Averaging complementary subsets
in (3) proves the exact zero-matching identity

```math
\boxed{
\phi_{K_0(A)}(t)
=\log\left\{2^{-r}\sum_{J\subseteq[r]}
 \exp\bigl(\phi_{A^J}(2t)\bigr)\right}.}            \tag{5}
```

Thus the universal core is not a function of the one child pressure or of
the child energy histogram.  It is a log-sum-exp over the entire
induced-clique-negation orbit at a temperature larger by `sqrt(2)`.

For completeness, put

```math
Z_B^\pm(u)=\mathbb E_xe^{\pm uH_B(x)},
\qquad
v_B(u)={Z_B^+(u)-Z_B^-(u)\over Z_B^+(u)+Z_B^-(u)}.
```

If `mathcal R` contains one member of each complementary pair
`{J,J^c}`, then (3)--(4) also give the exact integral formula

```math
\boxed{
e^{\phi_{K_\eta(A)}(t)}
=2^{1-r}\sum_{J\in\mathcal R}
e^{\phi_{A^J}(2t)}
\left[
 \cosh(t(r-2|J|))
 +\eta v_{A^J}(2t)\sinh(t(r-2|J|))
\right].}                                           \tag{6}
```

No orientation term was discarded in (6).

## 2. Exact pressure-to-defect implication

Fix `beta>0` and put

```math
s={\beta\over\sqrt r},
\qquad
t={\beta\over\sqrt{2r}}={s\over\sqrt2}.
```

Let `A` be an actual minimizer of `P_r(beta)`.  More generally, let
`K_b(A)` fill the `r` missing matching edges of `K_0(A)` by arbitrary
signs `b_i`.  Averaging independent fair choices of the `b_i` before
taking the logarithm gives

```math
\mathbb E_b e^{\phi_{K_b(A)}(t)}
=e^{\phi_{K_0(A)}(t)}(\cosh t)^r.
```

Thus some valid sign parent has pressure at most
`phi_(K_0(A))(t)+r log cosh(t)`.  Combining this with (5) gives the
sharpened direct recurrence

```math
\boxed{
\begin{aligned}
E_{r,r}(\beta)
&\le
 \log\left\{2^{-r}\sum_J
 e^{\phi_{A^J}(\sqrt2s)}\right\}
 -2\phi_A(s)+r\log\cosh t\\
&\le
 \log\left\{2^{-r}\sum_J
 e^{\phi_{A^J}(\sqrt2s)}\right\}
 -2\phi_A(s)+{\beta^2\over4}.
\end{aligned}}                                      \tag{7}
```

For comparison, choosing the same matching sign `eta` at every coordinate
only gives the cruder pointwise estimate

```math
\boxed{
|\phi_{K_\eta(A)}(t)-\phi_{K_0(A)}(t)|
\le tr=\beta\sqrt{r/2}.}                            \tag{8}
```

The annealed filling argument is deterministic by the probabilistic
method and improves this `O_beta(sqrt(r))` payment to `O_beta(1)`.

The needed new theorem is therefore completely explicit.  It is enough
that, at every large `r`, **one** actual minimizer (which may be selected
from the minimizing fibre) satisfy, for some `delta>0`,

```math
\log\left\{2^{-r}\sum_J
 e^{\phi_{A^J}(\sqrt2s)}\right\}
\le2\phi_A(s)+O_\beta(r^{1-\delta}),                \tag{9}
```

then

```math
\boxed{
E_{r,r}(\beta)
=O_\beta\!\left(r^{1-\delta}+1\right)
=o_\beta(r).}                                       \tag{10}
```

This is the exact `P => E` implication for the universal double.  The
matching `I` is already harmless at the requested scale; all substantive
content is (9).

Exact minimality alone points in the opposite direction.  It says
`phi_{A^J}(s)>=phi_A(s)` for every `J`, while (9) needs an upper bound on an
exponential orbit average at the larger temperature `sqrt(2)s`.

## 3. An analytic actual-minimizer falsifier

Take the order-four signing

```math
A=\begin{pmatrix}
0& 1&-1&-1\\
1& 0&-1&-1\\
-1&-1&0&-1\\
-1&-1&-1&0
\end{pmatrix}.                                      \tag{11}
```

Up to switching, permutation, and global negation, order four has two
cosh-energy classes.  Their normalized partitions are

```math
a(s)={4+8\cosh(2s)+4\cosh(4s)\over16}
={1\over2}\cosh(2s)(1+\cosh(2s))                   \tag{12}
```

for (11), and

```math
b(s)={8+6\cosh(2s)+2\cosh(6s)\over16}.              \tag{13}
```

Writing `c=cosh(2s)` gives

```math
b(s)-a(s)={1\over2}(c-1)^2(c+1)>0                  \tag{14}
```

for every `s>0`.  Hence (11) is an actual exact order-four pressure
minimizer at every positive temperature.

The exact core and integral-double partitions are

```math
k_0(t)=
{72+124\cosh(4t)+56\cosh(8t)+4\cosh(12t)\over256}, \tag{15}
```

and

```math
k_+(t)=\frac1{256}\left[
40+64\cosh(2t)+56\cosh(4t)+48\cosh(6t)
+24\cosh(8t)+16\cosh(10t)+8\cosh(12t)
\right].                                            \tag{16}
```

At scaled temperature `beta=4`, one has `s=2` and `t=s/sqrt(2)=sqrt(2)`.
Both hoped-for inequalities fail:

```math
\boxed{
\begin{aligned}
\log k_0(\sqrt2)-2\log a(2)
 &=0.252257576116984\ldots>0,\\
\log k_+(\sqrt2)-2\log a(2)
 &=1.019529787938264\ldots>0.
\end{aligned}}                                      \tag{17}
```

The first sign in (17) has a short exact proof, independent of the decimal.
Put `c=cosh 4`.  The elementary bounds `e>8/3` and `e^8>256` give

```math
c>25,
\qquad
c<{257\over512}e^4,
\qquad
a(2)<{13\over25}\left({257\over512}\right)^2e^8. \tag{18}
```

The two positive-energy `12` states in the full core cube give
`k_0(sqrt(2))>=e^{12sqrt(2)}/128`.  Moreover,

```math
e^{12\sqrt2-16}>e^{4/5}
>{1389\over625}
>{737257497769\over335544320000}
=128\left({13\over25}\right)^2
 \left({257\over512}\right)^4.                     \tag{19}
```

Equations (18)--(19) prove `k_0(sqrt(2))>a(2)^2`.  The integral double has
twice as many positive-energy `12` states, so the same argument proves its
strict inequality as well.

More generally,

```math
\begin{aligned}
\log k_0(s/\sqrt2)-2\log a(s)
  &=(6\sqrt2-8)s-\log2+o(1),\\
\log k_+(s/\sqrt2)-2\log a(s)
  &=(6\sqrt2-8)s+o(1)
\end{aligned}                                       \tag{20}
```

as `s` tends to infinity.  Thus the failure is analytic and persists on a
whole temperature ray.  Equation (17) is a floor on this construction
certificate, not a lower bound on `E_{4,4}`: the true order-eight minimizer
may use another parent.

### The archived order-six double has a larger finite floor

The order-six conference signing with root-gauge code `220` has

```math
c_6(s)={20\cosh(3s)+12\cosh(5s)\over32}.          \tag{20a}
```

It is an exact pressure minimizer for every `s>0`, not only at the audited
grid points.  Root gauge leaves `2^10=1024` order-six signings and produces
nine distinct absolute-energy histograms.  Besides `{3:20,5:12}`, the other
eight are

```text
{1:12,3:8, 5:8,7:4}
{1:12,3:10,5:6,7:3,9:1}
{1:12,3:12,5:4,7:2,9:2}
{1:14,3:8, 5:7,7:2,11:1}
{1:14,3:10,5:5,7:1,9:1,11:1}
{1:14,3:11,5:4,7:2,13:1}
{1:15,3:10,5:6,15:1}
{1:24,7:6,9:2}.
```

For each row `h`, substitution of `u=e^s` gives

```math
64u^{15}\{Z_h(s)-c_6(s)\}
=(u-1)^4(u+1)^4(u^2+1)^3R_h(u),                  \tag{20b}
```

where `R_h` is a nonzero polynomial with nonnegative coefficients.  Thus
every difference is strictly positive for `u>1`.  The exact histogram
fibres are included in the JSON certificate cited below, so (20b) is a
finite polynomial check rather than a floating-point ordering claim.

The core and integral parent both have cap `18`, while the child cap is
`5`.  Their exact histograms at scaled temperature `beta=4` give

```math
\boxed{
\begin{aligned}
\phi_{K_0(A)}(4/\sqrt{12})-2\log c_6(4/\sqrt6)
 &=2.444721543\ldots,\\
\phi_{K_+(A)}(4/\sqrt{12})-2\log c_6(4/\sqrt6)
 &=2.717572032\ldots.
\end{aligned}}                                    \tag{20c}
```

At large `beta`, both certificate defects have positive slope

```math
{18\over\sqrt{12}}-{10\over\sqrt6}
=3\sqrt3-{5\sqrt6\over3}>0.                       \tag{20d}
```

This is still a finite-order construction floor, not a scalable lower bound
on the true defect.

### An exact order-eight thermal collision

Order eight already shows that even the **full child energy law** does not
determine the orbit pressure.  The two switching/permutation classes with
canonical root-gauge codes `44280` and `111980` have the same projective
signed-energy histogram

```text
{-10:4,-8:10,-6:12,-4:16,-2:16,0:12,
   2:16,4:16,6:12,8:10,10:4}.
```

Both classes are actual pressure minimizers at every positive temperature,
not merely at the audited grid points.  Here is a finite exact certificate.
The common absolute histogram is

```math
h_* =\{0:12,2:32,4:32,6:24,8:20,10:8\}.            \tag{20e}
```

For any of the `96` absolute histograms `h` in the order-eight exhaustive
certificate, put

```math
L_h(u)=\sum_e h_e\bigl(u^{e/2}+u^{-e/2}\bigr),
\qquad u=e^{2s}.
```

Exact integer long division gives, for every `h != h_*`,

```math
u^{14}\bigl(L_h(u)-L_{h_*}(u)\bigr)
=(u^2-1)^4R_h(u).                                   \tag{20f}
```

For `93` of the `95` differences, `R_h` is a nonzero polynomial with
nonnegative integer coefficients.  The remaining two quotients are,
respectively,

```math
\begin{aligned}
R_1(u)&=2u^8(1-u+2u^2-u^3+u^4),\\
R_2(u)&=2u^8(2-u+4u^2-u^3+2u^4).
\end{aligned}                                       \tag{20g}
```

They are positive for `u>1`: in the first bracket group
`u^3(u-1)+(2u^2-u+1)`, and in the second group
`u^3(2u-1)+(4u^2-u+2)`.  Since
`L_h(u)=256 Z_h(s)`, (20f)--(20g) prove strict thermal minimality of
`h_*` for every `s>0`.

The doubles nevertheless have different exact histograms and caps:

```text
                         code 111980       code 44280
core cap                      32                40
integral +I cap               32                40
beta=4 core defect       1.565740894       7.959683657
beta=4 integral defect   1.986203979       7.967889534
beta=8 core defect       5.338579158      19.948258639
beta=8 integral defect   5.409820802      19.948282357
```

Thus scalar minimality, and even the complete signed scalar energy
histogram, cannot imply (9) for an arbitrarily chosen minimizer.  This is
not a collision of labeled geometry: the two induced-clique orbits differ.
It also does **not** obstruct the universal-double route, because a
recurrence may select the favorable minimizer.  Concretely, because the
integral doubles are valid parents, the favorable class gives the finite
consequences

```math
E_{8,8}(4)\le1.986203979\ldots,
\qquad
E_{8,8}(8)\le5.409820802\ldots .                    \tag{20h}
```

There is no asymptotic conclusion for `E_{r,r}` from this collision alone.

## 4. Scalable scalar-minimality obstruction

The finite example can be amplified while preserving the exact pressure
minimizer property and the signing Frobenius normalization.  This identifies
precisely what an attempted proof of (9) must use beyond scalar minimality.

Let `q>=1`, `r=4q`, and

```math
w_q=\sqrt{{4q-1\over3}}.
```

Consider the weighted block class

```math
\mathcal C_q=\left\{
w_q\operatorname{diag}(C_1,\ldots,C_q):
C_i\text{ is an order-four signing}
\right\}.                                           \tag{21}
```

Every member is hollow and has the exact dense-signing second moment

```math
\operatorname{Tr}B^2
=12q w_q^2
=r(r-1).                                             \tag{22}
```

The entries are weighted and cross-block entries vanish, so (21) is not the
complete-signing class.

Let `B_q=w_q diag(A,...,A)` with `A` from (11).  Then `B_q` is an exact
pressure minimizer in `mathcal C_q` at every positive temperature.  To see
this without a numerical claim, put `u=e^{2s}`.  For the other order-four
class, the two oriented exponential partitions are

```math
p(s)={6u^{-1}+8+2u^3\over16},
\qquad
q(s)={6u+8+2u^{-3}\over16}.
```

The partition `a(s)` of (11) is orientation-symmetric, and direct
factorization gives

```math
p(s)q(s)-a(s)^2
={(u-1)^4(u+1)^4\over32u^4}\ge0.                    \tag{23}
```

For arbitrary blocks, the full cosh partition is

```math
{\prod_i Z_{C_i}^+(s)+\prod_i Z_{C_i}^-(s)\over2}
\ge\sqrt{\prod_iZ_{C_i}^+(s)Z_{C_i}^-(s)}
\ge a(s)^q.                                         \tag{24}
```

This proves exact minimality of `B_q` in (21).

Now fix

```math
\beta_0=2\sqrt3,
\qquad
s_q={\beta_0w_q\over\sqrt r}
=2\sqrt{1-{1\over4q}}\longrightarrow2.             \tag{25}
```

After reordering coordinates, `K_0(B_q)` is the direct sum of `q` copies of
`w_qK_0(A)`.  Its pressure therefore factorizes exactly, and its universal-
double certificate has defect

```math
D_q^{(0)}
=q\left[\log k_0(s_q/\sqrt2)-2\log a(s_q)\right].  \tag{26}
```

By (17) and continuity,

```math
\boxed{
\lim_{q\to\infty}{D_q^{(0)}\over r}
={0.252257576116984\ldots\over4}
=0.063064394029246\ldots.}                          \tag{27}
```

Filling the matching edges changes log pressure by at most (8), so

```math
D_q^{(+)}
\ge D_q^{(0)}-\beta_0\sqrt{r/2}.
```

Consequently the integral-double certificate has the same linear floor:

```math
\boxed{
\liminf_{q\to\infty}{D_q^{(+)}\over r}
\ge0.063064394029246\ldots.}                        \tag{28}
```

Equations (22), (24), and (28) are a scalable method-class obstruction:
exact scalar pressure minimality, the correct `1/sqrt(r)` temperature, and
even the exact Frobenius identity of a complete signing do not imply (9).
Any positive theorem for actual signings must use their coordinate-flat,
complete-support geometry to exclude the positive-density localized blocks
in (21).  The obstruction does **not** prove a linear floor for actual
complete-signing children.

## 5. Reproducible finite audit and verdict

[`audit_universal_double_actual_minimizers.py`](../computations/audit_universal_double_actual_minimizers.py)
exhausts all switching-gauged signings through child order eight, keeps
global negations separate, and records every pressure-minimizing
switching/permutation class on the scaled-temperature grid.  The exact
integer histograms and evaluated defects are frozen in
[`universal_double_actual_minimizers_n8.json`](../computations/results/universal_double_actual_minimizers_n8.json).
The computation independently verifies (5) to absolute error below
`9e-16` in every recorded case.

The permanent target is therefore **not proved**.  What survives is the
single exact obligation (9).  The matching edges cost only `O_beta(1)`
after annealed sign completion; spectral or scalar-pressure information
cannot control the induced-clique orbit, and the weighted block floor
proves that exact minimality plus variance cannot do so either.
