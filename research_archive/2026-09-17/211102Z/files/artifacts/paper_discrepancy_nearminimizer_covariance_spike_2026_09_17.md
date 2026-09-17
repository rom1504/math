# Near-minimizing full signings can force every ground-law covariance to diverge

2026-09-17. Discrepancy-track strengthening of an archived structural
counterexample. The diverging-covariance consequence below was derived
independently; a subsequent archive check found the same planted-clique
mechanism already in the September 6 orientation-gap construction.
Section 4 gives the strongest same-order consequence of that earlier
argument, with better cap error. Sections 1--3 retain a self-contained
random-bridge alternative, not a claim of a new basic construction.
The enlarged signings are quantitatively near-minimizing, NOT
claimed to be exact global minimizers. That distinction is essential.

## 1. Finite pinned-block theorem

Let B be a full signing of order n with `max H_B=Q(B)=M_n`; an exact
minimizing core can always be globally negated to have this polarity.
Append an m-vertex ferromagnetic block F=J_m-I_m and any full sign
bridge C. Write

```
A = [[B,C],[C^T,F]],       N=n+m,
L_C=max_(x,y Boolean)|x^T C y|,       f_m=binom(m,2).
```

Then the resulting genuine full signing satisfies

```
M_n+f_m <= Q(A) <= M_n+f_m+L_C <= M_N+f_m+L_C.       (1)
```

The first inequality follows by taking a positive core maximizer,
setting every appended spin equal to one, and reversing all core spins
if necessary to make the bridge nonnegative. The upper bound is the
triangle inequality. Finally M_n<=M_N follows by principal restriction
and averaging the omitted spins.

The negative absolute sector is bounded separately:

```
max_(x,y) -H_A(x,y) <= M_n+floor(m/2)+L_C.           (2)
```

Indeed the ferromagnetic energy is
`H_F(y)=((sum_j y_j)^2-m)/2 >=-floor(m/2)`.
Consequently, for EVERY window

```
0<=T<f_m-floor(m/2)-L_C,                            (3)
```

the complete absolute nearcode E_A(T) lies in the POSITIVE energy
sector. Every one of its words (x,y) satisfies

```
(sum_j y_j)^2 >= m^2-2(L_C+T),                      (4)
H_B(x) >= M_n-L_C-T.                               (5)
```

For (4), compare `H_A>=Q(A)-T>=M_n+f_m-T` with
`H_A<=M_n+H_F(y)+L_C`. For (5), instead bound H_F by f_m.
Both assertions concern the WHOLE absolute nearcode, not a selected
ground family or one favorable energy polarity silently imposed.

Thus ANY probability law mu supported on E_A(T), including EVERY
maximizing dual law for the unrestricted physical-column response game,
has

```
||E_mu (x,y)(x,y)^T||op >= m-2(L_C+T)/m.             (6)
```

Test its covariance on the unit vector equal to 1/sqrt(m) on the
appended coordinates and zero elsewhere. No entropy assumption and no
choice of an unusually ill-conditioned optimal dual is involved.

## 2. A near-optimal family at the campaign's existing error scale

Independent fair bridge entries, followed by the elementary sign-MGF
bound and a union over all 2^N spin pairs, give existence of C with

```
L_C <= sqrt(2nm[(N+2)log 2]) <= C_0 N sqrt(m).       (7)
```

For example the union probability above this threshold is at most
1/2. All entries remain literal signs. Choose an even m asymptotic to

```
m=N^(2/3)(log N)^(1/12),          n=N-m,
```

and take an exact minimizing core of that order. Then

```
Q(A_N) <= M_N+O(N^(4/3)(log N)^(1/6)),
L_C/m^2=O((log N)^(-1/8)),
||Sigma_mu||op >= (1-o(1))m -> infinity              (8)
```

for EVERY ground-supported law mu. The same covariance conclusion
holds on every full window T_N=o(m^2). Condition (3) holds eventually
for each such window. The construction works at every sufficiently
large order N and assumes neither existence nor nonexistence of the
normalized cap limit.

Its cap error is smaller than the campaign's
`N^(4/3)(log N)^(1/3)` same-order preparation/comparison budget.
Therefore neither asymptotic near-optimality nor membership in that
cap-error class forces a bounded-covariance optimal query law. The
signing itself also has operator norm at least m-1, by testing the
appended all-one direction. Bounded spectral norm cannot be inferred
from those cap-error conditions either.

This does NOT disprove either property for EXACT global minimizers.
The ferromagnetic block has a positive subleading cost and may prevent
exact optimality. Establishing or refuting a structural theorem using
exact minimality remains a separate obligation.

## 3. Why this spike is removable, not a hard-response theorem

The example also identifies an important distinction for attempted
spectral trimming. The high covariance is concentrated in only m=o(N)
coordinates, and (5) projects its full microscopic nearcode into a
nearcode of an ACTUAL exact minimizing core. Under the positivity-window
hypothesis (3), in particular T=o(m^2), the loss in signed energy is
only L_C+T=o(N^1.5). Merely imposing T=o(N^1.5) without (3) would not
justify the positive-polarity conclusion.

If a physical law nu_old bounds the response on that projected code
by b, extend it with independent fair signs on the m deleted coordinates.
For every original word,

```
E|h_old dot x+h_new dot y| <= b+sqrt(m).             (9)
```

This follows from the triangle inequality and the exact variance m
of the independent new sum. Centering, isotropy, and a common linear
subGaussian proxy are preserved if they held for nu_old (the proxy
becomes its maximum with one). Since sqrt(m)=o(sqrt(N)), the covariance
spike causes no leading scalar-response obstruction by itself.

Accordingly the construction falsifies a literal bounded-covariance
shortcut, but not a more flexible theorem allowing removal of negligible
coordinates while preserving high energy. A generic spectral projector
need not correspond to such coordinates; the example does not resolve
that harder reduction for exact minimizing signings.

Reproducibility:
`computations/paper_discrepancy_2026_09_17_nearminimizer_covariance_spike.py`
enumerates the entire signed parent cube for four finite constructions
and checks (1)--(6) as exact integer inequalities. Imported core
global-minimality labels remain archive provenance.

## 4. Archive collision and a stronger same-order covariance consequence

The earlier
`artifacts/transfer_adversary_nearmin_clique_orientation_gap_2026_09_06.md`
already plants a clique inside a low-norm principal core of a SAME-ORDER
exact minimizer. Its stated conclusions concern orientation and exclusion
of approximately isotropic laws. Its finite inequalities also imply the
stronger divergence of the covariance norm, as follows.

Start with an exact minimizer A of order N, switch and negate it so
H_A(1)=M_N, and choose S of size r in the archived half-size spectral
core. Let u=Q(A_S). The simultaneous diagonal-majorant theorem gives

```
u<=4 K_G (r/N) M_N=O(r sqrt(N)),
K_G=pi/[2 asinh(1)].
```

Its complete elementary majorant proof was reread in Section 1 of
`artifacts/resumed_bound_audit_minimal_proof_2026_09_06.md`: polarization
bounds the bilinear sign norm by 4Q, odd tensor/Gaussian rounding bounds
the diagonal-majorant SDP by K_G times that norm, and deleting the
largest half of its diagonal gives the stated core. No lower-bound or
optimizer-regularity conjecture enters that dependency.

Overwrite only A_S by a positive clique, obtaining A'. The archived
finite inequalities are

```
Q(A')>=M_N+binom(r,2)-u,
max -H_A'<=M_N+r/2+u,
Q(A')<=M_N+binom(r,2)+u.
```

For every T below the resulting positive/negative sector gap, every
word of the COMPLETE absolute T-nearcode therefore satisfies

```
((sum_(i in S) x_i)^2-r)/2
  =H_A'(x)-H_A(x)+H_A_S(x)
 >=binom(r,2)-2u-T.
```

Consequently every supported law mu has

```
||E_mu xx^T||op >= r-4u/r-2T/r.                     (10)
```

Taking r=floor(N^(2/3)) and T=o(r^2) gives

```
Q(A'_N)=M_N+O(N^(4/3)),
inf_(mu supported on the full T-nearcode) ||E_mu xx^T||op
       >=(1-o(1))N^(2/3) -> infinity.               (11)
```

This is stronger than the random-bridge cap-error scale in (8). It is
an explicit new consequence of the old construction, not a new near-
minimizer mechanism. Neither version settles the covariance question
for exact global minimizers. The removable-core response statement of
Section 3 is specific to that independent appended-block construction;
it is NOT silently transferred to the overwritten-block version here.
