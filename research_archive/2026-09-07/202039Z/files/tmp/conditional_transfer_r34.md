# Wave 34 Route A: penalized block response and conditional witness transfer

## Status and verdict

The minimax identities and transfer bounds below are **Verified**.  The
decisive `A_9` certificates use only integer response tables and explicit
primal/dual witnesses; they are checked by
`tmp/conditional_transfer_r34_check.py`.  The checker also exhausts the
stated finite `A_6` claims.

The naive bridge is false: exact global minimality does **not** force the
zero-loss block-minimax law to put any mass in a prescribed favorable child
fiber.  On the `3+6` block of the exact `A_9`, the unpenalized response value
is `q_9=24`, but every optimal dual law gives zero mass to the exact-child
ground fiber.  In the switched optimal `A_6` replacement orbit, every
maximizing witness lies outside that fiber with a gap of at least eight.

There is nevertheless an exact positive interface.  A favorable-fiber
penalty gives a convex response profile whose retained value, or whose
right derivative at zero, is precisely the missing hypothesis.  A strict
gap between the full response game and the bad-state-only game quantitatively
forces favorable mass.  Global minimality controls the value of the
unpenalized game up to the rounding error, but controls neither this gap nor
the penalty derivative.

These are finite mechanism walls and conditional theorems.  They do **not**
give an asymptotic minimizer family, disprove convergence, or falsify the
weaker uniform-selector target in `tmp/uniform_partial_coset_r34.md`.

## 1. Exact penalized response game

Let `A` be an exact order-`n` minimizer and write `q=Q(A)=q_n`.  Fix a block
`U`, let `E=E(U)`, `M=|E|`, and let `Omega_U` be the `2^|U|` oriented
projective local states.  For `u in Omega_U`, retain

```math
d_e(u)=\sigma x_ix_j,
\qquad
Z(u)=\max_{x_{U^c}}\sigma\{2x_U^{\mathsf T}A_{U,U^c}x_{U^c}
+x_{U^c}^{\mathsf T}A_{U^c,U^c}x_{U^c}\}.
```

For a fractional internal replacement `z in [-1,1]^E`, put

```math
R_z(u)=Z(u)+2\sum_{e\in E}z_ed_e(u).
```

Let `a=(a_e)_(e in E)` be the actual internal signing, let

```math
q_U=Q(A[U]),
\qquad
\delta_U(u)=q_U-2\sum_ea_ed_e(u),
```

and define the parent slack after taking the best outside response by

```math
\Delta(u)=q-R_a(u)\ge0.
```

Thus the favorable child fiber at cutoff `D` is
`F_D={u:delta_U(u)<=D}`.  Let `p:Omega_U->[0,infinity)` be any penalty; the
two important choices are `p=delta_U` and
`p=1_(Omega_U setminus F_D)`.  Define

```math
\boxed{
V_p(\tau)=\min_{z\in[-1,1]^E}\max_{u\in\Omega_U}
\{R_z(u)-\tau p(u)\},\qquad \tau\ge0.
}
\tag{R34.A1}
```

Finite minimax and minimization over the cube give the exact dual

```math
\boxed{
V_p(\tau)=\max_{\mu\in\Delta(\Omega_U)}
\left\{\mathbb E_\mu Z
-2\sum_{e\in E}|m_e|-\tau\mathbb E_\mu p\right\},
\qquad m_e=\mathbb E_\mu d_e.
}
\tag{R34.A2}
```

If `mu_tau` is any optimizer, comparison with the actual block gives the
key **exact accounting identity**

```math
\boxed{
q-V_p(\tau)
=\mathbb E_{\mu_\tau}\Delta
+4\sum_{e\in E}(a_em_e)_+
+\tau\mathbb E_{\mu_\tau}p.
}
\tag{R34.A3}
```

Indeed, the dual equality says
`E Z=V_p+2sum|m_e|+tau E p`, while
`E Delta=q-E Z-2sum a_em_e`, and
`2(|m|+am)=4(am)_+` for `a in {+-1}`.  Every term on the right of (R34.A3)
is nonnegative.

Consequently, if `V_p(tau)>=q-epsilon`, one law simultaneously satisfies

```math
\mathbb E\Delta\le\varepsilon,
\qquad
4\sum_e(a_em_e)_+\le\varepsilon,
\qquad
\mathbb E p\le\frac{\varepsilon}{\tau}.
\tag{R34.A4}
```

For `p=1_(F_D^c)`, `tau>0`, and `t>0`, this yields the concrete conditional
transfer bound

```math
\boxed{
\mu_\tau\{u\in F_D:\Delta(u)\le t\}
\ge1-(q-V_p(\tau))
\left(\frac1\tau+\frac1t\right).
}
\tag{R34.A5}
```

This is the desired near-active common law if the right side is positive.
It keeps the orientation, rather than silently replacing an oriented state
by an unoriented cut.

## 2. The exact missing derivative and a testable sufficient gap

Write

```math
J(\mu)=\mathbb E_\mu Z-2\|\mathbb E_\mu d\|_1,
\qquad V_0=\max_\mu J(\mu).
```

The function `V_p` is convex and nonincreasing in `tau`.  Its one-sided
derivatives have the exact envelope interpretation

```math
V'_{p,+}(\tau)
=-\min\{\mathbb E_\mu p:\mu\text{ optimizes }V_p(\tau)\},
\tag{R34.A6}
```

with the left derivative using the maximum instead.  In particular, for
`p=1_(F_D^c)`,

```math
\boxed{
\max_{\mu\in\operatorname{argmax}J}\mu(F_D)
=1+V'_{p,+}(0).
}
\tag{R34.A7}
```

Thus the favorable mass is exactly a one-sided response susceptibility.
The rounding argument of (10.945) gives only
`V_0>=q-eta_|U|`; the automatic Lipschitz bound
`V_p(tau)>=V_0-tau` gives no positive lower bound in (R34.A5).  Controlling
`V_0` does not control (R34.A7).

There is a useful stronger, directly falsifiable sufficient hypothesis.
Put `B=Omega_U setminus F_D` and

```math
V_B=\min_{z\in[-1,1]^E}\max_{u\in B}R_z(u)
=\max_{\mu:\operatorname{supp}\mu\subset B}J(\mu),
```

and set

```math
W_F^+=\max_{\nu:\operatorname{supp}\nu\subset F_D}
\{\mathbb E_\nu Z+2\|\mathbb E_\nu d\|_1\}.
```

If `g=V_0-V_B>0`, every full-game optimal law with favorable mass `alpha`
obeys

```math
\boxed{
\alpha\ge\frac{g}{W_F^+-V_B}.
}
\tag{R34.A8}
```

To prove it, decompose the law as
`mu=(1-alpha)mu_B+alpha mu_F` and use the reverse triangle inequality:

```math
J(\mu)\le(1-\alpha)V_B+\alpha W_F^+.
```

The denominator is finite and
`W_F^+<=q+4M`, because the parent cap gives `Z(u)<=q+2M`.
Combining (R34.A8) with (R34.A3) at `tau=0`, for `t>0`, gives

```math
\mu\{u\in F_D:\Delta(u)\le t\}
\ge\frac{V_0-V_B}{W_F^+-V_B}-\frac{q-V_0}{t}.
\tag{R34.A9}
```

Hence a bad-state response gap exceeding the fractional-rounding loss by a
quantified amount is a concrete sufficient extra lemma.  The gap is not
necessary: cancellation can permit both bad-only and favorable-containing
optimal laws.  The exact derivative (R34.A7) is the minimal formulation;
(R34.A8) is the more robust testable sufficient condition.

## 3. Exact `A_9` terminal wall: zero favorable mass

Take the displayed exact `A_9`, `q_9=24`, and

```text
U = {3,4,5,6,7,8}.
```

Its induced block has `q_U=14`.  At the target-pair normalization,

```math
B_{9,6}=\frac{16\sqrt6}{3}-10=3.0639\ldots,
```

so the zero-slack favorable fiber is exactly `F_0={delta_U=0}`; it has four
oriented states.  In lexicographic local-edge order, the fractional block

```text
z_* = (-1,1,1,1,-1,1,0,0,1,-1,1,-1,-1,0,0)
```

satisfies

```math
\max_uR_{z_*}(u)=24,
\qquad
\max_{u\in F_0}R_{z_*}(u)=16.
\tag{R34.A10}
```

The matching dual certificate uses the local spin `x=++++--` in its two
opposite absolute orientations.  Their data are

```text
(sigma, delta_U, internal energy, Z, parent score)
(-1,       16,              -2, 26,           24)
(+1,       12,               2, 22,           24).
```

Their edge vectors cancel and their mean external response is `24`.
Therefore, for **every** fractional or complete internal replacement, one
of these two bad states scores at least `24`.  Together with (R34.A10),

```math
\boxed{V_0=V_B=24.}
\tag{R34.A11}
```

More strongly, any optimal dual law must be supported on states active at
`z_*`: its mean score at `z_*` is at least its minimax value `24`, while
every pointwise score is at most `24`.  Since no favorable state is active,

```math
\boxed{
\max_{\mu\in\operatorname{argmax}J}\mu(F_0)=0.
}
\tag{R34.A12}
```

For `p=1_(F_0^c)`, the same certificates give the exact initial penalty
profile

```math
\boxed{V_p(\tau)=24-\tau\qquad(0\le\tau\le8),}
\tag{R34.A13}
```

so `V'_(p,+)(0)=-1`, in exact agreement with (R34.A7).  This rigorously
falsifies a theorem asserting that global block minimality alone supplies a
zero-loss common response law with positive mass in every prescribed
favorable child fiber.  It does not say that `F_0` contains no parent
ground; it says favorable grounds cannot participate in a law retaining the
full block-minimax value.

The witness migration persists in exactly the switched optimal replacement
family of (10.702)--(10.705).  Replace the terminal block by every
switched/complemented copy of the displayed exact `A_6`.  Across the 64
translates, the full replacement norms have histogram

```text
28^6, 32^20, 36^26, 40^12.
```

For every translate, every maximizing state lies outside `F_0`; the best
`F_0` score trails the maximum by at least eight.  Thus even the strict
translate-cover margin in (10.705) is carried entirely by migrated child
labels.

## 4. Exact one-orientation shield

Take instead

```text
U = {0,1,2,6}.
```

Here `q_U=8` and `B_(9,4)=28/9`.  Two local states, both with
`sigma=+1`, have spins `+-++` and `+++-`; each has child deficit `12`,
internal energy `-4`, outside response `28`, and original parent score
`24`.  Their mean edge vector in lexicographic local order is

```text
(0,1,0,0,-1,0).
```

Hence their mean replacement score is

```math
28+2(z_{02}-z_{16})\ge24
\qquad(z\in[-1,1]^6).
\tag{R34.A14}
```

It follows that the bad-state continuous and integral values equal `24`
for every favorable cutoff through child deficit eight.  This is not merely
the cancellation of opposite absolute orientations: one orientation alone
supports the obstruction.  The signed internal-energy control in (10.945)
therefore cannot be converted into orientation-free favorable-fiber control.

## 5. Uniform-selector weakening and exact `A_6` scope

The uniform-selector reduction in `tmp/uniform_partial_coset_r34.md` is a
real logical weakening: the ultimate restriction argument need only find a
large enough uniform-selector set, not transfer a witness for every fixed
selector.  The single `A_9` block above does not falsify that target.

However, selector averaging alone does not repair the response theorem.  On
the displayed exact `A_6`, for every proper `U` of sizes three, four, and
five, an explicit opposite-orientation pair outside the exact-child ground
fiber has cancelling edge vectors and mean response `q_6=10`.  Therefore

```math
V_B=V_0=q_6=10
```

for all `20+15+6=41` such selectors, in both the continuous and integral
games.  At zero target slack, each `B_(6,m)` lies below the first positive
child deficit, so this is the actual favorable fiber.  This is an exact
finite, zero-slack, proper-selector wall.  It is **not** an asymptotic family
and does not rule out a softened cutoff, a rare good selector set at large
order, or a proof of the uniform center moment by a different mechanism.

Moreover, per-selector favorable response laws would still not by themselves
produce the one common cut/center needed by the uniform moment.  A useful
continuation must prove a selector-averaged derivative or bad-gap estimate
together with cross-selector overlap, rather than merely replace `for every
U` by an average of the same failed local statement.

## 6. Surviving open targets

1. **Penalized retention.**  For a softened asymptotic favorable fiber, prove
   `q-V_p(tau)=o(tau)` (or the quantitatively weaker version needed in
   (R34.A5)) using structure beyond block minimality.
2. **Bad-fiber gap.**  Prove `V_0-V_B>0` with a uniform quantitative bound on
   a sufficiently large uniform-selector set, then use (R34.A8)--(R34.A9).
3. **Orientation coupling.**  Supply a relation between positive and negative
   child orientations that excludes the same-orientation shield (R34.A14).
4. **Cross-selector overlap.**  Even successful local transfer must be made
   coherent across enough selectors to produce a shared low-row cut or
   center; the response theorem alone does not do this.

The correct conclusion is therefore negative for the proposed direct bridge
but positive as a reduction: the missing input is exactly a penalty-slope or
bad-response-gap theorem, and `A_9` gives sharp finite falsification tests for
any candidate proof.
