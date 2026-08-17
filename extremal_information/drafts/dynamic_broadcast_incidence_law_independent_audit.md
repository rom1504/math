# Independent audit of `dynamic_broadcast_incidence_law.md`

## Verdict

**PASS WITH A LOCAL PROOF REPAIR AND AN EXPLICIT ATOM MODEL.**

The finite-depth decomposition (DB.3), exact rate identity (DB.4), incidence
conclusions (DB.6)--(DB.9), and their constants are valid under the intended
bounded-dictionary-atom interpretation.  There is no missing terminal or
transient factor.  In fact DB.6 is slightly looser than the direct charge.

The written proof of DB.7 uses
`C_i<=sum_e|Delta m_e^i|` and then charges a dictionary atom only once.  That
step is unjustified if one atom can have nonzero components on many visible
edges, as the current prose permits.  The theorem does not need edge-local
atoms: replace that step by the stronger inequality `C_i<=M_i` and state the
atom expansion explicitly.  This repairs the proof without changing any
displayed theorem.

The finite experiment
[`../experiments/verify_dynamic_broadcast_incidence.py`](../experiments/verify_dynamic_broadcast_incidence.py)
passes, although it exercises the three-channel decomposition rather than
the general incidence bookkeeping.

## 1. Exact interpretation of the atom assumption

A sufficient precise version of Section 2 is

```math
u_q^z=\sum_\alpha c_\alpha(z)\phi_{\alpha,q},
\qquad
a_e^z=\sum_\beta d_\beta(z)\psi_{\beta,e},             \tag{A.1}
```

with

```math
\max_q\|\phi_{\alpha,q}\|_2\le1,
\qquad
\max_e\|\psi_{\beta,e}\|_2\le1.                     \tag{A.2}
```

The coefficient is invariant under flipping bit `i` outside its dependency
set, and a flip inside it changes the coefficient by at most its stated
oscillation.  Atoms may be global dictionaries over all vertices or edges;
edge locality is not required.

For a probability-space vector `f`, orthogonal mean/centred decomposition
gives

```math
|\pi f|\le\|f\|_2,
\qquad
\|f-(\pi f)1\|_2\le\|f\|_2.                         \tag{A.3}
```

Thus one coefficient oscillation can be charged once to its mean projection
and once to its centred projection, each with coefficient one.  No hidden
factor two occurs inside either channel.

## 2. Audit of DB.3

For the difference dictionary, Theorem 33.1 gives the exact orthogonal
decomposition

```math
D_p=\left(\Delta\bar u_{q_L}^i+
          \sum_{s=1}^L\Delta m_{e_s}^i\right)1+Z_p,
```

where

```math
\|Z_p\|_2\le\rho^LR_i+{1-\rho^L\over1-\rho}B_i.     \tag{A.4}
```

Delete directed cycles from the visible path.  Deleted cycles have total
absolute scalar weight at most their total length times `C_i`; the remaining
simple directed path has at most `|Q|-1` edges, each of magnitude at most
`M_i`.  Hence

```math
\left|\sum_s\Delta m_{e_s}^i\right|
\le LC_i+(|Q|-1)M_i.                                 \tag{A.5}
```

Adding the terminal scalar bound `U_i` and (A.4) proves DB.3 exactly as
stated.  There is no omitted `U_i` multiplier, no extra transient cycle
factor, and no terminal `R_i` term without `rho^L`.

For the degenerate one-vertex graph with no edge, arbitrarily long paths do
not exist.  Either exclude that convention by assuming the strongly
connected carrier has a directed cycle, or inherit Theorem 33.1's convention
that the asymptotic assertion is vacuous.  This does not affect any
nontrivial carrier.

## 3. Audit of DB.4 and the limsup

The upper bound follows from DB.3 after division by `L`.  For a directed
simple cycle attaining `C_i`, repeat it.  The scalar part grows by the signed
cycle sum, while terminal and centred terms stay bounded.  Constant and
centred channels are orthogonal, so centred cancellation cannot reduce the
scalar lower bound.  Taking the subsequence of cycle-multiple lengths proves

```math
\limsup_L d_L(z,z^i)/L=C_i.
```

Periodicity creates no issue because only a limsup is asserted.

There is also no invalid interchange of a limsup and a growing sum in DB.7:
`h` is finite and DB.4 is first applied separately to every `i`, giving

```math
\sum_i\limsup_L d_L(z,z^i)/L=\sum_i C_i.             \tag{A.6}
```

## 4. Audit and sharpening of DB.6

Let

```math
A_i^U=\sum_{\alpha:i\in I_\alpha}\omega_\alpha,
\qquad
A_i^E=\sum_{\beta:i\in I_\beta}\omega_\beta.
```

Equations (A.1)--(A.3) imply

```math
U_i\le A_i^U,quad R_i\le A_i^U,
\qquad
B_i\le A_i^E,quad M_i\le A_i^E.                    \tag{A.7}
```

Every cycle average is bounded by the maximum edge magnitude, so crucially

```math
C_i=\chi_G(\Delta m^i)\le M_i\le A_i^E.             \tag{A.8}
```

This is the correct replacement for the edge sum in the draft proof.  On
summing over bits,

```math
\sum_i A_i^U=J_U,
\qquad
\sum_i A_i^E=J_A.                                    \tag{A.9}
```

Substitution into DB.3 actually proves the sharper bound

```math
\sum_i d_L(z,z^i)
\le(1+\rho^L)J_U+
\left({1-\rho^L\over1-\rho}+L+|Q|-1\right)J_A.       \tag{A.10}
```

Since `1+rho^L<=2`, A.10 implies DB.6.  The additional `+1` multiplying
`J_A` in DB.6 is harmless slack; it is not needed for scalar/centred double
charging.  The latter is already represented by the separate geometric and
scalar coefficients.

## 5. Audit of DB.7--DB.9

By (A.6), (A.8), and (A.9),

```math
\sum_i\limsup_L d_L(z,z^i)/L
=\sum_i C_i\le J_A,
```

which proves DB.7 for global as well as edge-local atoms.  If all `h`
neighbours have rate at least `epsilon`, DB.8 follows immediately.  Finally,
for `E_0` scalar atoms, fan-in at most `t`, and coefficient oscillation at
most `2B`,

```math
J_A\le\sum_{\beta=1}^{E_0}t(2B)=2BtE_0,
```

so DB.9 has the correct constant and direction.

## Required repair before canonicalization

1. Insert an explicit dictionary expansion equivalent to (A.1)--(A.2), or
   state unambiguously that “component norm” means the maximum `L^2` norm
   across the whole vertex/edge dictionary.
2. Replace `C_i<=sum_e|Delta m_e^i|` in the proof by `C_i<=M_i`.
3. Optionally replace DB.6 by the sharper A.10; keeping the displayed DB.6
   is mathematically valid.
4. Clarify the vacuous no-cycle singleton convention.

No theorem-level weakening is required.

## Resolution check

The task draft now includes the explicit global dictionary expansion
(DB.5a), uses `C_i<=M_i` in the proof, and states the no-recurrent-cycle
convention.  These changes implement the required repair.  The displayed
DB.6 intentionally retains its harmless slack.  After the repair, the audit
verdict is **PASS**.
