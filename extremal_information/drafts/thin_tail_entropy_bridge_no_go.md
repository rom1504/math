# The conference entropy-tilt wall survives thin-tail conditioning

**Status.** Task-local proved no-go, awaiting independent audit.  The
centerpiece is an actual-signing theorem: the archived fixed-small-tilt
Paley-conference obstruction survives after conditioning on a set of bridge
outputs every one of which satisfies the new bounded-cap thin-tail theorem.
All constants, normalizations, and the conditioning estimate are explicit.

The result is deliberately narrow.  It does **not** rule out a large fixed
tilt, a growing tilt, a different bridge law, or a microcanonical state
retaining energy-resolved overlap/root information.  It says that the new
fixed-rate tail in spin space does not repair the already isolated
fixed-small-disorder-temperature bridge mechanism.  Conference children
are not known to minimize the contracted-temperature pressure, so this is
also not a falsifier of a minimizer-specific compensation theorem.

## 1. Exact setup and archived input

For a hollow signing `S` of order `N`, use

```math
H_S(z)=\frac12z^TSz=\sum_{i<j}s_{ij}z_iz_j,
\qquad
Q(S)=\max_z|H_S(z)|,
\tag{CT.1}
```

and the normalized `cosh` partition function

```math
\overline Z_N(S,t)=2^{-N}\sum_z\cosh(tH_S(z)).
\tag{CT.2}
```

Let `A_r` be a symmetric conference signing,

```math
A_r^2=(r-1)I.
\tag{CT.3}
```

For a bridge `B in {+-1}^{r times r}` and a relative orientation
`epsilon in {+-1}`, put

```math
S_{\epsilon,B}
=\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix},
\qquad
t={\beta\over\sqrt{2r}},
\qquad
L_{\epsilon,B}=\log\overline Z_{2r}(S_{\epsilon,B},t).
\tag{CT.4}
```

Write `U_r` for the uniform law on the `2^(r^2+1)` pairs `(epsilon,B)`.
For `lambda>0`, its bridge soft minimum is

```math
\mathcal R_{\lambda,r}
=-{1\over\lambda}\log
  \mathbb E_{U_r}e^{-\lambda L_{\epsilon,B}}.
\tag{CT.5}
```

The already audited Theorem 6.1 of
`artifacts/two_temperature_bridge_audit.md` proves the following.  Define

```math
\psi(c)={1\over4}\left[
 \sqrt{1+4c^2}-1-
 \log\left({1+\sqrt{1+4c^2}\over2}\right)\right],
\tag{CT.6}
```

```math
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4},
\qquad
\gamma(\beta)=h_\beta-2\psi(\beta)>0.
\tag{CT.7}
```

For every

```math
0<\beta<{\sqrt2\over6},
\tag{CT.8}
```

there is `lambda_0(beta)>0` such that, along the Paley conference
sequence, every fixed `0<lambda<lambda_0(beta)` satisfies

```math
\mathbb E_{U_r}e^{-\lambda L_{\epsilon,B}}
=\exp\{-\lambda h_\beta r+o(r)\},
\tag{CT.9}
```

while the same-temperature child target is

```math
{1\over r}\,2\log\overline Z_r
  (A_r,\beta/\sqrt r)\longrightarrow2\psi(\beta).
\tag{CT.10}
```

Thus the unconditioned fixed-tilt defect is
`gamma(beta)r+o(r)`.  The question here is whether the subsequently proved
bounded-cap thin-tail theorem removes the outputs responsible for that
defect.

## 2. Almost every bridge output has the uniform thin tail

Define

```math
\mathcal E_r
=\left\{B:\max_{x,y\in\{+-1\}^r}|x^TBy|
                    \le2r^{3/2}\right\}
\tag{CT.11}
```

and put

```math
c_0=2-2\log2>0.
\tag{CT.12}
```

### Lemma CT.1 (explicit bounded-cap bridge event)

One has

```math
\boxed{\Pr_B(\mathcal E_r^c)\le2e^{-c_0r}.}
\tag{CT.13}
```

For every `B in E_r` and both `epsilon=+-1`, the complete order-`2r`
signing in (CT.4) obeys

```math
\boxed{
Q(S_{\epsilon,B})
\le r\sqrt{r-1}+2r^{3/2}
\le3r^{3/2}<2(2r)^{3/2}.}
\tag{CT.14}
```

**Proof.**  For fixed Boolean `x,y`, `x^TBy` is a sum of `r^2`
independent signs.  Hoeffding gives

```math
\Pr_B\{|x^TBy|>2r^{3/2}\}
\le2\exp\left\{-{(2r^{3/2})^2\over2r^2}\right\}
=2e^{-2r}.
\tag{CT.15}
```

There are at most `2^(2r)` pairs, so the union bound is

```math
2^{2r}\,2e^{-2r}=2e^{-(2-2\log2)r},
```

which proves (CT.13).  Since `||A_r||_(2 to 2)=sqrt(r-1)`,

```math
|H_{A_r}(x)|\le{r\sqrt{r-1}\over2}.
\tag{CT.16}
```

The two internal blocks together cost at most `r sqrt(r-1)`, and the
cross term is `x^TBy`.  This proves (CT.14). `square`

Theorem 36.26 now applies with the fixed choice `C=2` to **every** output
retained by `E_r`.  Since its explicit constant is

```math
d_C={1\over200000\max\{C,1\}},
\tag{CT.17}
```

we may take

```math
d_2={1\over400000}.
\tag{CT.18}
```

After the harmless two-sided union adjustment already included in Theorem
36.26, there is a fixed `kappa_2>0` such that every retained parent satisfies

```math
\boxed{
\#\left\{z:
Q(S_{\epsilon,B})-|H_{S_{\epsilon,B}}(z)|
 <d_2(2r)^{3/2}\right\}
\le\exp\{(\log2-\kappa_2)2r\}.}
\tag{CT.19}
```

No exceptional retained output and no choice of orientation is omitted.

## 3. Main theorem: conditioning does not change the negative-moment rate

Let `U_r^E` be the law obtained by conditioning `U_r` on `B in E_r`, and
define

```math
\mathcal R^E_{\lambda,r}
=-{1\over\lambda}\log
 \mathbb E_{U_r^E}e^{-\lambda L_{\epsilon,B}}.
\tag{CT.20}
```

### Theorem CT.2 (fixed-tilt wall inside the uniformly thin-tail class)

Fix `beta` in (CT.8).  Set

```math
\boxed{
\lambda_*(\beta)
=\min\left\{\lambda_0(\beta),
             {2-2\log2\over2h_\beta}\right\}>0.}
\tag{CT.21}
```

Along the Paley conference sequence, every fixed
`0<lambda<lambda_*(beta)` satisfies

```math
\boxed{
{\mathcal R^E_{\lambda,r}\over r}\longrightarrow h_\beta,
\qquad
{\mathcal R^E_{\lambda,r}
 -2\log\overline Z_r(A_r,\beta/\sqrt r)\over r}
\longrightarrow\gamma(\beta)>0.}
\tag{CT.22}
```

Every output over which the left side is averaged satisfies the common
bounded-cap and thin-tail statements (CT.14) and (CT.19).  Nevertheless the
conference same-temperature defect is linear rather than
`O(r^(1-delta))`.

**Proof.**  First, `L_(epsilon,B)>=0`, because the normalized average of
`cosh` is at least one.  Therefore Lemma CT.1 gives the exact discarded-
moment estimate

```math
0\le
\mathbb E_{U_r}\left[
  1_{\mathcal E_r^c}e^{-\lambda L_{\epsilon,B}}\right]
\le2e^{-c_0r}.
\tag{CT.23}
```

On the other hand, (CT.9) and (CT.21) give

```math
\mathbb E_{U_r}e^{-\lambda L_{\epsilon,B}}
=e^{-\lambda h_\beta r+o(r)},
\qquad
c_0-\lambda h_\beta>{c_0\over2}.
\tag{CT.24}
```

Consequently the ratio of (CT.23) to the full moment in (CT.24) is at most

```math
2\exp\{-(c_0-\lambda h_\beta)r+o(r)\}=o(1).
\tag{CT.25}
```

Also `U_r(E_r)=1-o(1)`.  Hence

```math
\begin{aligned}
\mathbb E_{U_r^E}e^{-\lambda L}
&={\mathbb E_{U_r}[1_{\mathcal E_r}e^{-\lambda L}]
   \over U_r(\mathcal E_r)}\\
&=\exp\{-\lambda h_\beta r+o(r)\}.
\end{aligned}
\tag{CT.26}
```

Taking `-(1/lambda)log` proves the first limit in (CT.22).  Subtracting
(CT.10) and using (CT.7) proves the second. `square`

## 4. Why a spin tail does not supply a bridge basin

Theorem CT.2 is an actual-signing falsifier for the fixed-small-tilt route.
There is also a simple abstract separation explaining the quantifiers.

Let a bridge output cube have `J_N=Theta(N^2)` bits.  Fix
`0<kappa<log2`, `c>d>0`, and `beta c>log2`.  Give one output an even
landscape with two antipodal states of absolute energy `cN^(3/2)` and all
other states at zero.  Give every other output an even landscape with

```math
K_N=2\left\lfloor{e^{(\log2-\kappa)N}\over2}\right\rfloor
\tag{CT.27}
```

antipodally closed top states at the same energy and all other states at
zero.  Every output has the same cap and satisfies the same upper-tail
condition

```math
\#\{x:Q-|H(x)|<dN^{3/2}\}
\le e^{(\log2-\kappa)N}.
\tag{CT.28}
```

Their normalized `cosh` pressures satisfy

```math
L_{\rm bad}-L_{\rm good}
=(\log2-\kappa)N+o(N).
\tag{CT.29}
```

If the one good output has uniform mass `p_N=2^(-J_N)`, then its bridge
soft minimum obeys, with `Delta_N=L_bad-L_good`,

```math
0\le L_{\rm bad}-\mathcal R_{\lambda,N}
\le p_N\Delta_Ne^{\lambda\Delta_N}.
\tag{CT.30}
```

Thus for every positive `lambda_N=o(N)`,

```math
\mathcal R_{\lambda_N,N}=L_{\rm bad}+o(N).
\tag{CT.31}
```

The inequality follows directly from

```math
\mathcal R_{\lambda,N}
=L_{\rm bad}-{1\over\lambda}
 \log\{1+p_N(e^{\lambda\Delta_N}-1)\}
\tag{CT.32}
```

and `(e^u-1)/u<=e^u`.  Detecting the isolated low-pressure output needs
`lambda_N=Omega(N)`, the disorder-temperature scale at which a
`Theta(N^2)`-bit output cube is resolved to order-`N` pressure accuracy.

This construction is a rigorous abstract bridge model, not a complete-
signing realization and not a counterexample to a signing-specific basin
theorem.  Its exact lesson is that a tail rate in spin space imposes no
basin rate in bridge-disorder space.

## 5. Scalar thermal consequence and microcanonical scope

The contemporaneous task-local note
`drafts/bounded_cap_thin_tail_pressure_consequence.md` proves the complete
scalar consequence of Theorem 36.26:

```math
{\beta Q(A)\over n^{3/2}}-\log2
\le {1\over n}\log\overline Z_n(A,\beta/\sqrt n)
\le {\beta Q(A)\over n^{3/2}}
 -\min\{\kappa_C,\beta d_C\}+O(1/n).
\tag{CT.33}
```

It also gives a finite weighted quadratic countermodel satisfying the
fixed-rate thin tail, exact variance, monotonicity, adjacent-order
regularity, and exact centered subadditivity while both the normalized
maxima and every fixed-temperature diagonal pressure oscillate.  Those
results subsume the scalar pressure calculation; it is not repeated here.

There are two distinct microcanonical objects:

1. the spin density of states inside one parent;
2. the lower-tail rate of parent pressure as bridge disorder varies.

Theorem 36.26 controls one point of the first object.  The exact disorder-
counting theorem in `artifacts/microcanonical_disorder_counting_composition.md`
still contracts the child temperature and does not locate the second
object's support edge at order-`N` resolution.  The abstract model in
Section 4 shows there is no implication between the two rates without
additional structure, while Theorem CT.2 shows that the fixed-tilt wall
persists on an actual-signing output class uniformly certified by the new
tail theorem.

Therefore a microcanonical interpolation can still help only if it retains
information coupling the two spaces, such as an energy-resolved
overlap/rooted response rate or a direct bridge-disorder support-edge
theorem.  A scalar cap plus one fixed spin-tail deficit cannot produce a
summable same-temperature defect.

## 6. Archive comparison and stopping judgment

1. **Two-temperature collision.**
   `artifacts/two_temperature_bridge_audit.md` already proves the exact
   Renyi identity, linearly rare basin criterion, and unconditioned
   fixed-small-tilt conference wall.  The new increment is precisely
   Lemma CT.1 plus Theorem CT.2: intersect with an explicit
   `1-exp(-Theta(r))` event on which every output satisfies the later-proved
   Theorem 36.26, and show the negative-moment wall is unchanged.
2. **Scalar collision.**  The separate thin-tail pressure note proves a
   stronger scalar-axiom no-go.  This draft does not claim its elementary
   pressure sandwich as new.
3. **Positive tail use is compatible.**  Theorem BCL.0 uses a fixed-rate
   spin tail positively by retaining a chosen root and constructing a
   specially matched bridge.  That result amplifies rooted geometry;
   scalar entropy tilting does not.  There is no contradiction.
4. **Remaining route.**  The only entropy-weighted target not closed here
   is the one already isolated in the archive: prove a signing-specific
   `exp(-O(N))` low-pressure bridge basin or a coupled spin/disorder rate
   theorem, in particular for actual contracted-temperature pressure
   minimizers.  Merely increasing the tilt without such a theorem approaches
   full bridge minimization.  Since conference children are not known to be
   those minimizers, CT.2 is a benchmark obstruction rather than a universal
   no-go for the sufficient criterion (4.4) in the two-temperature audit.

**Classification.**  Lemma CT.1 and Theorem CT.2 are rigorous conditional
only on the already audited Theorem 36.26 and conference fixed-tilt theorem.
Section 4 is a rigorous abstract sharpness example, explicitly not a
quadratic-signing construction.  No rigorous interval for
`M_n/n^(3/2)`, cross-order recurrence, or convergence theorem is improved.
