# Wave 39: selector--cut information transports parent row to favorable Gram

## Outcome

There is a verified averaged inequality which removes the internal-Gram term
as an independent selection requirement, provided the favorable labels are
restrictions of a low-information family of low-parent-row complement
certificates.

Let `P(S,D)` be any joint law, where `S` is an `m`-selector and the full
projective word represented by `D` is `z^D`.  Put

```math
\mathsf H
=D(P_S\Vert U_m)+I_P(S;D),
\qquad
\overline R=\mathbb E_P R_2(z^D).
```

Then, for an exact minimizer `A`,

```math
\boxed{
\begin{aligned}
\mathbb E_P\lVert A[:,S]z^D_S\rVert_2^2
\le{}&p^2\overline R+p(1-p)n(n-1)\\
&+C\sqrt{q_n(n^2+\overline R)
                 (\mathsf H+\log(n+1))}\\
&+Cq_n(\mathsf H+\log(n+1)).
\end{aligned}}
\tag{R39.G1}
```

The constant is absolute and the estimate is uniform in `m`.  At

```math
\overline R=O(n^{9/4-c}),
\qquad
\mathsf H=O(n^{3/4-c}),
```

every term in (R39.G1) is `O(n^(9/4-c))`.  Since internal Gram is bounded
pointwise by full-column Gram, this proves

```math
\mathbb E_P\lVert A[S]z^D_S\rVert_2^2
=O(n^{9/4-c}). \tag{R39.G2}
```

Thus the exact remaining selection problem is no longer (10.1080) plus an
unexplained internal-row clause.  In the high-ratio window `p_2>=1/2`, it is
to select an affordable coherent kernel of complement certificates with:

```math
\boxed{
\begin{aligned}
D(P_S\Vert U_m)+I(S;D)&=O(n^{3/4-c}),\\
\mathbb E R_2(D)&=O(n^{9/4-c}),\\
\text{anchored expected conflict}&\le\text{the project conflict budget}.
\end{aligned}}
\tag{R39.G3}
```

Support on complement incidences makes every restricted label favorable by
(10.1045), regardless of certificate deficit.  Equation (R39.G1) gives the
whole Gram quantity directly, including the formerly open internal half.  A
finite probabilistic argument derandomizes the kernel while preserving all
averaged clauses up to absolute constants.

This is genuine progress but not a convergence proof: no construction of a
kernel satisfying all three lines of (R39.G3) is known.  The remaining issue
is now a joint low-information, low-row, coherent congestion lemma.  It
cannot be obtained from an unrestricted deletion witness or independent
noise.

## 1. Fixed-word quadratic selector transport

Fix a full word `z` and set

```math
G_z=D_zA^2D_z,
\qquad
K_S(z)=\lVert A[:,S]z_S\rVert_2^2.
```

If `xi` is the selector indicator, then exactly

```math
K_S(z)=\xi^TG_z\xi. \tag{R39.G4}
```

Let independent `xi_i~Ber(p)`, with `p=m/n`, and put `eta=xi-p1`.  Expansion
gives

```math
\begin{aligned}
K_\xi(z)-\mathbb E K_\xi(z)
={}&2p\langle G_z\mathbf1,\eta\rangle\\
&+\{\eta^TG_z\eta-\mathbb E(\eta^TG_z\eta)\},
\end{aligned}
\tag{R39.G5}
```

with baseline

```math
\mathbb E K_\xi(z)
=p^2R_2(z)+p(1-p)n(n-1). \tag{R39.G6}
```

The standard bounded-variable Hanson--Wright mgf estimate, the subgaussian
linear mgf, and Cauchy--Schwarz between the two exponentials imply

```math
\log\mathbb E\exp\{\lambda(K_\xi-\mathbb EK_\xi)\}
\le C\lambda^2
\left(\lVert G_z\rVert_F^2+\lVert G_z\mathbf1\rVert_2^2\right)
\tag{R39.G7}
```

whenever `0<=lambda<=c/||G_z||_op`.  No independence between the linear and
quadratic terms is asserted; Cauchy--Schwarz is why their separate mgfs may
be combined.

The exact-minimizer spectral estimate `||A||_op^2<=2q_n` gives

```math
\begin{aligned}
\lVert G_z\rVert_{op}&=\lVert A\rVert_{op}^2\le2q_n,\\
\lVert G_z\rVert_F^2
&=\operatorname{tr}A^4
\le2q_n n(n-1),\\
\lVert G_z\mathbf1\rVert_2^2
&=z^TA^4z
\le2q_nR_2(z).
\end{aligned}
\tag{R39.G8}
```

Now condition the Bernoulli sample on `sum_i xi_i=m`.  This is exactly
`U_m`.  Since `m` is a mode of `Bin(n,m/n)`, its probability is at least
`1/(n+1)`.  Hence conditioning adds only `log(n+1)` to the logarithmic mgf.
Entropy duality and optimization in `lambda` prove, for every selector law
`w`, with `H_w=D(w||U_m)+log(n+1)`,

```math
\boxed{
\begin{aligned}
\mathbb E_wK_S(z)
\le{}&p^2R_2(z)+p(1-p)n(n-1)\\
&+C\sqrt{q_n(n^2+R_2(z))H_w}+Cq_nH_w.
\end{aligned}}
\tag{R39.G9}
```

For comparison, the exact uniform-slice mean is (10.1050).  Using the
Bernoulli baseline costs only a harmless lower-order difference and avoids
any unproved transfer of pair marginals under `w`.

Equation (R39.G9) is the upper-transport counterpart to the lower projection
estimate (10.1054).  It is possible here because `K_S(z)` is a quadratic
selector polynomial with matrix `D_zA^2D_z`; the earlier internal statistic
`L_S(z)` is cubic in the selector and required degree-three information.

## 2. Joint selector--certificate version

Apply (R39.G9) conditionally on `D=d`, using the selector law `P(S|D=d)`.
The KL chain rule gives

```math
\mathbb E_DD(P_{S|D}\Vert U_m)
=D(P_S\Vert U_m)+I(S;D)=\mathsf H. \tag{R39.G10}
```

Averaging the linear terms and using Cauchy--Schwarz on the square-root term
gives

```math
\begin{aligned}
&\mathbb E_D
\sqrt{q_n(n^2+R_2(D))
 [D(P_{S|D}\Vert U_m)+\log(n+1)]}\\
&\qquad\le
\sqrt{q_n(n^2+\overline R)
 [\mathsf H+\log(n+1)]}.
\end{aligned}
```

This proves (R39.G1).  Notice that the certificate can depend on the
selector; the price is exactly mutual information rather than the number of
all possible parent cuts.

At the project scales,

```math
\sqrt{q_n(n^2+\overline R)(\mathsf H+\log n)}
=O(n^{9/4-c}),
\qquad
q_n(\mathsf H+\log n)=O(n^{9/4-c}),
```

so there is no exponent loss.

## 3. Complement certificates and conflict

Suppose `P_S` is uniform on an anchored family `mathcal G`, and a kernel
`pi_S(d)` is supported on complement-flip certificates for every `S`.  Write
the projective word of `d` with the common anchor equal to `+1` and set
`y^{S,d}=z^d_S`.

For independent `(S,D)` and `(T,D')` drawn from this anchored joint law,
define the kernel conflict

```math
\mathcal C_v(\pi)
=\mathbb E\sum_{i\in(S\cap T)\setminus\{v\}}
\frac{\mathbf1\{y_i^{S,D}\ne y_i^{T,D'}\}}{p_i}.
\tag{R39.G11}
```

In the high-ratio window, the support condition and (10.1045) make every
label favorable.  For reference, (10.1079) also gives pointwise

```math
q_n-Q(A[S])\le\frac{\Delta_D}{2}, \tag{R39.G12}
```

but neither cap persistence nor a deficit estimate is needed here because
(R39.G1) controls the full-column Gram directly.

Sample one certificate independently from `pi_S` for every selector in
`mathcal G`.  The expected average parent row and Gram are the corresponding
joint averages.  The expected deterministic corrected conflict
is at most (R39.G11): when the two sampled selectors happen to be identical,
the deterministic assignment contributes zero, whereas (R39.G11) may still
sample two different certificates.  Applying the probabilistic method to the
sum of the three normalized nonnegative costs gives one deterministic
assignment with each cost within a fixed factor of its budget.

Therefore (R39.G3) plus affordability bypasses (10.1080), proves the full
favorable Gram estimate (10.1057) directly, and then (10.1056) supplies the
desired global row cut.  No generic noise and no deletion of the term
`q_n-q_m` is used.

For a single parent cut `d`, take the kernel to be deterministic on `d` and
`mathcal G subseteq I_d`.  Then `I(S;D)=0`, the anchored conflict is zero, and
(R39.G1) shows that the internal Gram condition follows automatically from

```math
R_2(d)=O(n^{9/4-c}),\quad
U_m(\mathcal G)\ge e^{-O(n^{3/4-c})}. \tag{R39.G13}
```

This single-state version is especially transparent, but it is stronger than
the kernel formulation.

## 4. Fractional covers automatically pay the information term

The mutual-information clause has an exact fractional-cover interpretation.
Let nonnegative weights `w_d` cover every `S` in a family `mathcal G`:

```math
Z_S:=\sum_{d:S\in I_d}w_d\ge1,
\qquad W:=\sum_dw_d.
```

Define

```math
\pi_S(d)=\frac{w_d\mathbf1_{\{S\in I_d\}}}{Z_S},
\qquad r(d)=\frac{w_d}{W}.
```

Then exactly

```math
\mathbb E_{S\sim U(\mathcal G)}
D(\pi_S\Vert r)
=\mathbb E_S\log\frac W{Z_S}
\le\log W. \tag{R39.G14}
```

Since divergence to any reference law decomposes as

```math
\mathbb E_SD(\pi_S\Vert r)
=I(S;D)+D(P_D\Vert r),
```

we obtain

```math
\boxed{I(S;D)\le\log W.} \tag{R39.G15}
```

Consequently an affordable family and a project-weight fractional cover
already supply the information scale in (R39.G3).  If the cover is supported
on cuts with the desired parent row cap, (R39.G1) supplies Gram regularity
automatically.  What is **not** automatic is the anchored conflict of the
induced kernel.  The exact surviving fractional statement is:

> Construct a complement-incidence fractional cover of weight
> `exp{O(n^(3/4-c))}`, supported (or averaged) at parent-row
> `O(n^(9/4-c))`, whose induced anchored kernel has project-scale conflict.

This is a joint Lagrangian/congestion target, not a first-moment average of
the flipped matrices.  It retains the three quantities that earlier
arguments lost separately.

If this cover is required over the entire selector slice with a pointwise
parent-row cap, its weight bound already implies the row-good fractional
cover theorem (10.1047), and hence is not an easier proof of convergence.
The added value of (R39.G1) is in the agreement/exceptional-center setting:
on one merely affordable anchored family, a low-information auxiliary
certificate kernel no longer needs a separate internal-Gram theorem.  This
scope distinction prevents the reduction from being advertised circularly
as a proof of the full-cover target.

## 5. Exact remaining lemma and falsification scope

A clean hard sufficient lemma is the following.  For every target pair in a
fixed high-ratio window, find an anchored family `mathcal G` of relative
density `beta` and a kernel supported on complement incidences such that

```math
\begin{aligned}
\log\beta^{-1}+I(S;D)&\le Cn^{3/4-c},\\
\mathbb E R_2(D)&\le Cn^{9/4-c},\\
\mathcal C_v(\pi)&\le C_{\rm conflict}(n),
\end{aligned}
\tag{R39.G16}
```

where `C_conflict(n)` is the already prescribed project decoder budget.
Because `P_S=U(mathcal G)` and `mathcal G` is anchored,
`D(P_S||U_m)=log beta^{-1}+O(1)`.  Equations (R39.G1)--(R39.G12) and
derandomization then prove the full favorable Gram package.

This implementation is falsified, for fixed parameters, by an unbounded
exact-minimizer family on which every affordable anchored complement kernel
with project information has either

```math
\mathbb E R_2(D)=\omega(n^{9/4-c})
\quad\text{or}\quad
\text{project-scale anchored conflict fails.}
```

That would not falsify (10.1057), because favorable labels need not arise
from complement witnesses, and it would not falsify the bare arbitrary-cut
tail.

## 6. Verification

`tmp/internal_gram_transport_r39_check.py` checks, on `A_6,A_8,A_9`, the
exact uniform Gram mean, all three matrix identities in (R39.G8), and their
finite minimizer inequalities.  It also builds the unit-weight complement
cover on `A_6,m=5` and verifies numerically to `2e-12` the KL chain rule and
fractional-cover information identity (R39.G14)--(R39.G15).

Run:

```bash
.venv/bin/python tmp/internal_gram_transport_r39_check.py
```

The current run ends with `PASS internal_gram_transport_r39_check`.
