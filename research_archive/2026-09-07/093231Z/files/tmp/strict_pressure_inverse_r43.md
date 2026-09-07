# Wave 43 memo: strict mixed-pressure inverse extraction

Status: the extraction theorem below is **Verified**.  Its scalar
minimizer premise is an **Open target**, not a consequence of exact
minimality presently proved.  The finite tables are **Numerical**, from exact
enumeration of the integer deficit data followed by floating exponentiation.

## 1. Setup

Fix an exact-minimizer signing `A` of order `n`, a slice size `m`, and

```math
L_0=n^{3/4-c},\qquad 0<c<1/4,\qquad
H=B_{n,m}+t>0,\qquad \theta=\frac{bL_0}{H},
```

where `b>0` is fixed and `t=O(n^{3/2-c})`.  For a full spin `z` and a
selector `S`, define

```math
C_S(z)=z_S^{\mathsf T}A[S]z_S,\qquad
d_S(z)=Q(A[S])-|C_S(z)|\ge0,
```

and

```math
q(z)=\mathbb E_{S\sim U_m}e^{-\theta d_S(z)},\qquad
K_z^{\rm abs}(\theta)=\theta H+\log q(z).
```

The absolute sign is harmless at the final hard event: after a center is
chosen, split the favorable selectors according to the sign of `C_S(z)` and
retain the larger of the two classes.

## 2. A global child-free-energy lemma

Let `U_z` be the uniform spin law (using projective spins gives the same
quantities), and put

```math
\mathcal Z_A(\theta)
=\mathbb E_{z\sim U_z}q(z)
=\mathbb E_{S\sim U_m}\left[
 e^{-\theta Q(A[S])}
 \mathbb E_{y\sim U_{\{\pm1\}^m}}
 e^{\theta|y^{\mathsf T}A[S]y|}
\right].
\tag{R43.1}
```

Thus `mathcal Z_A` is a scalar average of exact absolute child partition
functions; it does not mention a center.

**Verified extraction theorem.**  Suppose fixed constants `0<alpha<b`
satisfy

```math
-\log\mathcal Z_A(bL_0/H)\le \alpha L_0.
\tag{R43.2}
```

Then there are a full spin `z` and one fixed orientation
`sigma in {+1,-1}` such that

```math
R_2(z)=O(n^{9/4-c}),
\tag{R43.3}
```

and

```math
U_m\left\{S:Q(A[S])-\sigma C_S(z)\le H\right\}
\ge \frac18e^{-\alpha L_0}
\tag{R43.4}
```

for all sufficiently large `n`.  By (10.861), every selector in (R43.4)
has `widehat ell(S,(sigma,z))<=t`.  Hence (R43.3)--(R43.4) are precisely the
project-scale arbitrary-cut certificate (10.795), and a uniform proof of
(R43.2) on the active density window would prove convergence.

### Proof

Tilt the joint uniform law of `(S,z)` by `e^{-theta d_S(z)}` and call the
result `P`.  Its exact relative entropy is

```math
D(P\Vert U_m\otimes U_z)
=-\theta\mathbb E_Pd_S(z)-\log\mathcal Z_A(\theta)
\le-\log\mathcal Z_A(\theta).
\tag{R43.5}
```

Data processing therefore gives, for its center marginal `P_z`,

```math
\mathcal H:=D(P_z\Vert U_z)\le\alpha L_0.
\tag{R43.6}
```

Apply the entropy-duality/Hanson--Wright inequality underlying
(10.865) and (10.870), now directly to
`R_2(z)=z^T A^2z`:

```math
\mathbb E_{P_z}R_2(z)
\le \operatorname{tr}(A^2)
+O\left(\lVert A^2\rVert_{\rm F}\sqrt{\mathcal H}
+\lVert A^2\rVert_{\rm op}\mathcal H\right).
\tag{R43.7}
```

Here `tr(A^2)=n(n-1)`.  For an exact minimizer,
`||A||_op<=sqrt(2q_n)`, so

```math
\lVert A^2\rVert_{\rm F}
\le\lVert A\rVert_{\rm op}\lVert A\rVert_{\rm F}=O(n^{7/4}),
\qquad
\lVert A^2\rVert_{\rm op}=O(n^{3/2}).
```

At `mathcal H=O(L_0)`, the two errors in (R43.7) are

```math
O(n^{17/8-c/2})\quad\hbox{and}\quad O(n^{9/4-c}).
```

The first is `O(n^(9/4-c))` exactly when `c<=1/4`; the trace term is also
absorbed for `c<1/4`.  This proves

```math
\mathbb E_{P_z}R_2(z)=O(n^{9/4-c}).
\tag{R43.8}
```

Let `G={z:R_2(z)<=2 E_(P_z)R_2}`.  Markov gives `P_z(G)>=1/2`.  Since
`P_z(z)=U_z(z)q(z)/mathcal Z_A`,

```math
\mathbb E_{U_z}[q(z)1_G]
=\mathcal Z_A P_z(G)\ge\frac12\mathcal Z_A.
```

Consequently some `z in G` obeys

```math
q(z)\ge\frac12\mathcal Z_A\ge\frac12e^{-\alpha L_0},
\tag{R43.9}
```

and it satisfies (R43.3).  In particular

```math
K_z^{\rm abs}(\theta)
\ge(b-\alpha)L_0-\log2.
\tag{R43.10}
```

This is at least `aL_0` for, say, `a=(b-alpha)/2` and large `n`.  The cap
soft-to-hard inequality (10.1132), or the same one-line calculation, gives

```math
U_m\{d_S(z)\le H\}
\ge\frac{e^{K_z^{\rm abs}(\theta)}-1}{e^{bL_0}-1}
\ge\frac14e^{-\alpha L_0}.
\tag{R43.11}
```

Assign selectors with `C_S(z)=0` to the positive class and split the event
in (R43.11) into `C_S(z)>=0` and `C_S(z)<0`.  One class has at least half
the mass.  On it one fixed sign `sigma` turns `d_S(z)` into
`Q(A[S])-sigma C_S(z)`, proving (R43.4).

This extraction supplies only the project cap (R43.3), not the older hard
class `mathcal C_2={R_2<=2n(n-1)}`.  That distinction is real, but (10.795)
needs only the project cap, so no convergence implication is lost.

## 3. Exact open premise and its honest falsifier

The new scalar sufficient lemma is:

> Uniformly on one fixed high-ratio window and all relevant target pairs,
> choose a target-specific exact minimizer `A` and prove that for some fixed
> `b>alpha>0`, (R43.2) holds at `theta=bL_0/H`.

Equivalently, with

```math
\Gamma_A(b)=-\frac1{L_0}\log\mathcal Z_A(bL_0/H),
\tag{R43.12}
```

one needs `Gamma_A(b)<b` with a uniform constant gap.  This is an exact
nonlinear child free-energy result.  Minimality alone has not established it:

- `0<=d_S(z)<=Q(A[S])<=q_n` gives only
  `-log mathcal Z_A<=theta q_n=b(q_n/H)L_0`, whose coefficient exceeds `b`;
- retaining one child ground gives only `mathcal Z_A>=2^(1-m)`, costing
  `Theta(n)` rather than `O(L_0)`;
- a small-tilt variance correction has size at most the already recorded
  wrong-scale correction and does not cancel the leading deficit.

An honest asymptotic falsifier for this aggregate lemma is an unbounded
exact-minimizer family for which

```math
\Gamma_A(b)\ge b-o(1)
\tag{R43.13}
```

for every admissible fixed `b`.  This would falsify (R43.2), but not the
exceptional-center route: a rare center can have positive strict pressure
while contributing too little to the first center moment.

The sharper aggregate formulation preserves such rare centers.  With
`nu_2` uniform on `mathcal C_2`, set

```math
\mathcal P_r(b)
=\frac1r\log\mathbb E_{z\sim\nu_2}
 e^{rK_z^{\rm abs}(bL_0/H)},
\qquad
r=\left\lceil\kappa\frac n{L_0}\right\rceil.
\tag{R43.14}
```

Then `max_(z in mathcal C_2)K_z>=mathcal P_r`.  Conversely,

```math
\mathcal P_r(b)
\ge\max_{z\in\mathcal C_2}K_z^{\rm abs}(bL_0/H)
-\frac{\log|\mathcal C_2|}{r}
\ge\max_zK_z^{\rm abs}-\left(\frac{\log2}{\kappa}+o(1)\right)L_0.
\tag{R43.15}
```

Thus proving `mathcal P_r(b)>=aL_0` is a sharp sufficient low-row replica
lemma, and by increasing the fixed `kappa` it is exponent-equivalent to the
existence of a strict-pressure center.  It is the faithful continuation if
the first-moment premise (R43.2) proves too strong; (R43.14) is a formulation,
not evidence that minimizers satisfy it.

## 4. Finite exact-minimizer audit

The checker `strict_pressure_inverse_r43_check.py` uses `t=0` and normalizes
`theta H=1`.  It verifies (R43.1), (R43.5), KL data processing, Markov
extraction, (10.1132), and the fixed-sign split.  Its principal output is:

| signing, `m` | `H` | `mathcal Z_A` | `-log mathcal Z_A` | `max_z K_z` | best absolute tail |
|:---|---:|---:|---:|---:|---:|
| `A_6,5` | 0.940591 | 0.625076 | 0.469882 | 0.817719 | 0.833333 |
| `A_8,6` | 2.276095 | 0.159979 | 1.832711 | 0.203325 | 0.392857 |
| `A_9,7` | 2.462453 | 0.058229 | 2.843369 | 0.041747 | 0.361111 |

For `A_8` and `A_9`, the normalized first-moment condition
`-log mathcal Z_A<theta H=1` fails although an old-`mathcal C_2` center has
strictly positive pressure and a substantial hard tail.  All centers in
these three finite examples happen to meet the old cap (maximum row squares
`30,72,128` versus caps `60,112,144`).  This is a finite, not asymptotic,
wall, but it demonstrates exactly that averaging centers can erase the
exceptional-center signal.  It favors the replica formulation (R43.14) over
treating (R43.2) as the likely final minimizer theorem.

## 5. Frontier recommendation

The verified contribution is a clean global-partition-to-project-row
inverse theorem: a scalar absolute child free energy with a constant
shortfall extracts one fixed oriented cut and proves the bare tail.  It does
not prove that scalar bound for minimizers.  Finite exact examples show the
first center moment is already strictly stronger than positive exceptional
pressure, so the primary strict-pressure target should remain the low-row
high-replica inequality (R43.14), with (R43.2) retained as a simpler strong
sufficient lemma.  Do not infer either premise from minimality, unweighted
moments, or selector self-bounding.
