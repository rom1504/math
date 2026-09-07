# Direct E stopping: minimal dependencies for the new cap upper bound

Date: 2026-09-07. Independent audit PASS. The exact fixed-point identity
H=E is NOT needed for the original upper bound below E.

The proof is the direct-T stopping proof in
`transfer_reconstruction_standalone_2026_09_06.md`, equations (25)--(26),
with E in place of T. The continuity and compactness inputs are also
proved in Sections1--4 of
`continued_convergence_terminal_gap_reduction_2026_09_06.md`.

Fix t, write `G(nu)=g_t(m_2(nu))`, `Phi=-F_t/2`,
`f_r=B^r Phi`, `D=Phi-E`, `Gamma=Phi-B Phi`, and `K=-g_t`.
The sufficient inputs are precisely:

- `B Phi<=Phi`, `G<=f_r`, `G<=E<=Phi<=0`;
- `B E<=E`, supplied by the precision/Schur theorem;
- Phi and E are weakly continuous on every bounded-second-moment ball,
  and Gamma is weakly lower semicontinuous there;
- Gamma vanishes only at centered Gaussian laws, where `Phi=G=E`;
- `K(s)/s->0` as s grows.

The lower `G<=f_r` uses iid pairings, or `B G>=G` followed by
monotonicity. It does not use E<=f_r or any conditional-copy construction.

For epsilon>0 and C>m_2(nu0), compactness gives

```math
\kappa=\inf\{\Gamma(\nu):m_2(\nu)\le C,
                             \Phi(\nu)-E(\nu)\ge\epsilon\}>0
```

unless the set is empty. Put `E0=m_2(nu0)` and
`B0=Phi(nu0)-G(nu0)`. Then the exact sufficient upper estimate is

```math
f_r(\nu_0)-E(\nu_0)
\le\epsilon+E_0\sup_{s>C}\frac{K(s)}s
       +\frac{K(C)B_0}{r\kappa}.                    (1)
```

To verify the budget, choose a depth-r policy within zeta of f_r. Its
value is at least G(nu0)-zeta. Along a uniformly chosen tree branch,
the expected sum of actual Phi drifts is Phi(root) minus the policy
value, at most B0+zeta; each actual drift is at least Gamma. The second
moment is a nonnegative martingale by the exact parallelogram identity.
Stop on D<=epsilon or moment>C, otherwise at depth r. Surviving branches
have probability at most `(B0+zeta)/(r kappa)`. Continuation values are
at most Phi, whose excess over E is bounded by epsilon, K(moment), or
K(C) in the three cases. Bounded optional stopping controls the energy
exit, and the stopped E supersolution controls the prefix. Letting zeta
decrease to zero proves (1).

Choose C large, then epsilon small, then finite r large. This proves
`limsup_r B^r Phi(nu0)<=E(nu0)`. It is only an upper comparison: no
lower E comparison is silently added. In particular, when the certified
ternary phase gives `E_4(nu_(31/32))=g_4(1)`, the existing fixed-depth
all-order weave construction yields the new cap upper bound directly.
The choice of finite depth precedes the matrix-order limit.

Thus the original upper-bound dependency chain can omit:
the construction of H, H's Bellman fixed-point property, the
conditional-copy proof H>=E, and any two-sided boundary identification.
Those theorems remain useful analytical statements, but are not required
for this consequence. Gaussian zero-drift rigidity, weak continuity on
moment balls, the Schur supersolution, the certified scalar phase, and
the actual finite-depth/all-order weave realization remain required.
