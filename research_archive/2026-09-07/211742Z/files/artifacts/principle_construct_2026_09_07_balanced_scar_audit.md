# Independent audit: macroscopic balanced scars and local margins

2026-09-07. **PASS.** Complete read of
`principle_invent_2026_09_07_macroscopic_balanced_scars.md`.

For each negative selected target in a rank-one balanced cross block,
opposite-sign companion rows and columns outside the selected sets remain
available: at least q/2-ell-(ell-1)>0 choices exist. Selected targets and
unused companions ensure that all four entries are still original when
the switch is made. Their pattern is exactly [-1,+1;+1,-1]. The switch
preserves every row and column sum and changes no other selected-selected
entry. Applying the transposed switch preserves physical symmetry.

The Hamiltonian change is 2(x_a-x_u)(y_c-y_v), bounded by eight. Thus the
safe total edit budget is 8ell binom(m,2), giving the corrected normalized
balanced-face cost 4r/sqrt(p). No probabilistic assumption about how many
targets need a switch is used in this cap estimate.

For independently uniform ordered selections, both the original selected
energy and each target's original sign have mean zero by exact block
balance. The expected scarred selected energy is therefore precisely
ell binom(m,2). Some selection attains this. Exact balance then gives
H(1-2z)=4H(z). The centered variance is 4mell(1-ell/q), so the normalized
variance coefficient is at least 1/[2sqrt(p)(1-r)]. Fixed positive r
ensures this survives any additive o(N^(3/2)) allowance.

## Explicit actual mixed gain despite uniform seed margins

There is a useful direct consequence for the suggested local-margin
criterion. Fill each diagonal fibre by J_q-I_q, and denote this block
diagonal seed component by A0. Every ground word of A0 is constant on
each fibre, with arbitrary signs between fibres. All those grounds have
the SAME uniform positive local field q-1, and

    Q(A0)=m q(q-1)/2,
    Q(A0)-H_A0(x)=2sum_i ell_i(q-ell_i),

where ell_i is the minority count in fibre i. In particular the seed has
a global linear barrier to its ground code, not merely one-spin stability.
The scarred balanced bulk vanishes on all these constant-fibre words.

For the scar word with the same ell minority vertices in every fibre,
the full ACTUAL signing A0+B' has energy above those constant grounds
by at least

    4ell binom(m,2)-2mell(q-ell)
      =2ell m[(m-1)-q(1-ell/q)].

For q/m->p<1 and ell/q->r>0, dividing by N^(3/2) gives the strictly
positive lower limit

    2r[(1-p)+pr]/sqrt(p).

Thus exact bulk balance, a strict balanced-face bulk certificate, and
even uniform local margins plus a linear global seed barrier do not
alone prevent an actual mixed excursion. The quantitative bulk response
must be smaller than the seed's actual deficit at the relevant densities.

A0 is the explicitly separated block-diagonal constant-mode component,
not an actual full-sign optimal child. A0+B' is a full hollow signing.
This consequence does not rule out favorable unscarred ensembles or
globally rewritten optimal-child constructions.

## Exhaustive finite regression

The independent checker
`computations/principle_construct_2026_09_07_balanced_scar_check.py`
and its JSON result exhaust m=q=4, ell=1: all 256 ordered selected-set
choices and all 1296 fibre-balanced Boolean words. All exact balance,
symmetry, four-edge switch, energy, variance, cap-edit, and ferro-completion
identities PASS. The mean scarred selected energy is exactly six. The
finite base balanced cap is 48 and scar caps range from 32 to 64; this
finite base is not claimed to obey the asymptotic strict-face premise.
At these parameters the ferro gain is exactly zero, in agreement with
the finite formula; the positive asymptotic regime uses p<1 as stated.
