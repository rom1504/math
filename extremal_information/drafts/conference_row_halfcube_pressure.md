# Exact rowwise antipodal halfcubes retain the conference pressure

**Status.**  Task-local theorem note.  This resolves the rowwise
antipodal-halfcube candidate identified after the speed-`r` basin audit.
It makes no canonical edits.

The family has exactly the desired speed-`r` entropy and its naive nearest
repair can require `Theta(r^(3/2))` bit flips.  Nevertheless its conditioned
pressure has the uniform-bridge rate `h_beta`, not a smaller rate.  The
mechanism is an exact gauge identity: rowwise canonicalization is an
adaptive switching of the first conference child, and the
operator-regular `r^2`-speed concentration theorem survives a union over
all `2^r` possible switchings.

## 1. Setup

Let `A_r` be a symmetric conference signing, fix
`epsilon in {+-1}`, and put

```math
f_{\epsilon,r}(B)
=\log\left[2^{-2r}\sum_{x,y}
 \cosh\left\{{\beta\over\sqrt{2r}}
 \big(H_{A_r}(x)+\epsilon H_{A_r}(y)+x^TBy\big)\right\}
 \right].
\tag{RH.1}
```

Assume

```math
0<\beta<{\sqrt2\over6}.
\tag{RH.2}
```

For a uniform sign bridge `U_r`, the audited conference theorem gives,
separately for both orientations,

```math
{f_{\epsilon,r}(U_r)\over r}\longrightarrow
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}
\tag{RH.3}
```

in probability and in mean.  The smaller same-temperature child target is

```math
\tau_\beta=2\psi(\beta)=h_\beta-\gamma(\beta),
\qquad \gamma(\beta)>0.
\tag{RH.4}
```

## 2. The exact halfcube family

Fix `u,v in {+-1}^r`.  Choose a deterministic odd tie rule

```math
\chi:\{R:\langle R,v\rangle=0\}\longrightarrow\{+-1\},
\qquad \chi(-R)=-\chi(R).
\tag{RH.5}
```

For row `i`, let `H_i` contain the unique member of every antipodal pair
`{R,-R}` satisfying

```math
u_i\langle R,v\rangle>0,
\tag{RH.6}
```

or, at a tie, `u_i chi(R)=+1`.  Thus

```math
|H_i|=2^{r-1}.
\tag{RH.7}
```

Define the rowwise halfcube

```math
\mathcal H_r(u,v)=\{B:B_{i,*}\in H_i\text{ for every }i\}.
\tag{RH.8}
```

Its cardinality is exactly

```math
\boxed{|\mathcal H_r(u,v)|=2^{r^2-r}.}
\tag{RH.9}
```

Equivalently it has probability `2^{-r}` under the full uniform bridge
law.  This exact count includes the even-`r` tie mass; no central-binomial
correction is hidden.

If `W` is a uniform bridge, define `s_i(W) in {+-1}` so that

```math
B_{i,*}=s_i(W)W_{i,*}\in H_i.
\tag{RH.10}
```

Then

```math
B=D_{s(W)}W
\tag{RH.11}
```

is uniform on `H_r(u,v)`: every output has exactly `2^r` preimages, obtained
by independently changing the signs of the input rows.

The elementary repair which flips only enough entries to change a typical
wrong row's majority can cost `Theta(sqrt(r))` flips per affected row and
therefore `Theta(r^(3/2))` in total.  Equation (RH.11), however, exposes a
different exact symmetry and avoids estimating that repair.

## 3. Uniform concentration over the switching orbit

For a deterministic `s in {+-1}^r`, changing variables `x -> D_sx`
gives the exact identity

```math
f_{A_r,\epsilon A_r}(D_sW)
=f_{D_sA_rD_s,\epsilon A_r}(W).
\tag{RH.12}
```

The switched child `D_sA_rD_s` is again a conference signing.  More
importantly, conjugating the full parent by `diag(D_s,I)` identifies its
operator norm with that of the parent having fixed children and bridge
`D_sW`.

Choose `delta>0` and `kappa<1/2` with

```math
{\beta(3+\delta)\over\sqrt2}<\kappa.
\tag{RH.13}
```

The standard rectangular Rademacher estimate gives

```math
\Pr\{\|W\|_{op}>(2+\delta)\sqrt r\}\le e^{-c_0r}.
\tag{RH.14}
```

On the complementary event, simultaneously for every `s`,

```math
\left\|{\beta\over\sqrt{2r}}
 \begin{pmatrix}D_sA_rD_s&W\\W^T&\epsilon A_r\end{pmatrix}
\right\|_{op}
\le\kappa.
\tag{RH.15}
```

Indeed, both conference diagonal blocks have norm `sqrt(r-1)` and row
switching does not change `||W||_op`.

For each fixed `s`, the convex extension from
`conference_regular_conditioned_all_tilts.md` has a dimension-free
Frobenius Lipschitz constant.  Talagrand convex-Lipschitz concentration is
two-sided, so for every fixed `eta>0`, all large `r` satisfy

```math
\Pr\left\{
 \left|f_{D_sA_rD_s,\epsilon A_r}(W)-h_\beta r\right|>\eta r,
 \ \|W\|_{op}\le(2+\delta)\sqrt r
 \right\}
\le2e^{-c_{\beta,\eta}r^2}.
\tag{RH.16}
```

The center is uniform in `s`: after the change of bridge variables
`W -> D_sW`, its law and the regular event are exactly those of the
unswitched conference calculation.  The `o(r)` error in the mean is
therefore the same for every `s`, rather than merely pointwise in `s`.

Union over the `2^r` switch vectors yields

```math
\Pr\left\{
 \max_s\left|f_{D_sA_rD_s,\epsilon A_r}(W)-h_\beta r\right|>\eta r
 \right\}
\le e^{-c_0r}+2^{r+1}e^{-c_{\beta,\eta}r^2}.
\tag{RH.17}
```

Thus even a selector that sees the entire bridge cannot find a different
linear pressure rate inside the regular switching orbit.

## 4. Main theorem

### Theorem RH.1 (the halfcube-conditioned pressure rate is unchanged)

Let `B_r` be uniform on `H_r(u,v)`, for arbitrary deterministic sign
vectors `u=u_r,v=v_r`.  Then, separately for the two orientations,

```math
\boxed{
{f_{\epsilon,r}(B_r)\over r}\longrightarrow h_\beta
\quad\hbox{in probability and in }L^1.}
\tag{RH.18}
```

More quantitatively, for every fixed `eta>0`, there are positive constants
`c_0,c_1` such that

```math
\boxed{
\Pr_{B_r\in\mathcal H_r(u,v)}
 \{|f_{\epsilon,r}(B_r)-h_\beta r|>\eta r\}
\le e^{-c_0r}+e^{-c_1r^2}.}
\tag{RH.19}
```

Consequently

```math
{1\over|\mathcal H_r(u,v)|}
\#\{B:f_{\epsilon,r}(B)
 \le\tau_\beta r+o(r)\}\longrightarrow0.
\tag{RH.20}
```

In particular, the rowwise antipodal halfcube is not a speed-`r`
same-temperature favorable basin.

**Proof.**  Generate the uniform halfcube bridge by (RH.11).  Identity
(RH.12) and the uniform event (RH.17) apply to the random adaptive choice
`s=s(W)`, proving convergence in probability and (RH.19).  Every sign
parent has `0<=f=O_beta(r^(3/2))`; multiplying this bound by the exceptional
probability `e^{-c_0r}+e^{-c_1r^2}` proves uniform integrability and hence
the `L^1` assertion.  Finally use the fixed gap
`h_beta-tau_beta=gamma(beta)>0`. `square`

## 5. Scope and archive comparison

1. This is a genuine exact speed-`r` candidate: its cardinality is
   `2^(r^2-r)`, not an asymptotic entropy estimate.
2. The result does not follow from the small-Hamming retraction theorem.
   Majority repair can have the critical `r^(3/2)` scale.  It is the
   switching covariance plus a simultaneous orbit concentration theorem
   that closes the family.
3. The switching-gauge quotient theorem for optimized bridges identifies
   equivalent child labels under a *minimum over bridges*.  RH.1 instead
   proves a quenched pressure theorem for a nonlinear canonical bridge
   fibre and an adaptive selector; it is a different consequence.
4. The regular-conditioned source states its concentration for fixed
   conference children.  The new step is uniformization over all `2^r`
   switched children.  The `r^2` tail survives because the switching orbit
   has only exponential size.
5. The theorem does not exclude an `e^{-Theta(r)}` exceptional subfamily
   carried by the operator-irregular event.  It proves that the declared
   halfcube constraint itself has the typical pressure rate and supplies no
   target-reaching certificate.
6. A genuinely different remaining rowwise candidate must be invariant
   under row sign--for example a constraint on
   `|<B_{i,*},v>|`--and must correlate operator irregularity with lower
   pressure.  A one-sided majority or any other antipodal selector is now
   closed.
