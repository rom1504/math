# Thick near-extreme levels already have linear factorization complexity

2026-09-07. **Proved caution.** This prevents overinterpreting the
Gram--Schmidt augmentation criterion as a new special property of
minimizers when it is applied only to a fixed-width entire superlevel.

Let B be symmetric hollow of order n, with Q(B)<=C n. For eta>0 set

    G_eta={x in {+-1}^n: |H_B(x)|>=Q(B)-eta n}.

Then G_eta contains an axis-aligned Boolean face of dimension at least

    floor(min(n/2, eta n/(80C))).                        (1)

Consequently its VC dimension, and therefore the square of any
unit-column factorization row bound gamma, is at least (1). This holds
for EVERY bounded-cap matrix, not only an asymptotic minimizer.

Proof. The simultaneous diagonal Grothendieck majorant from
`decisive_audit_fresh_full_lower_chain_2026_09_07.md`, Section 1, gives

    D>=B,-B,    Tr D<=K_G beta(B)<=4K_G Q(B)<=8C n,

using the elementary bound K_G<2. Choose an orientation sigma and
a maximizing spin g with H_(sigma B)(g)=Q(B). Its signed local fields

    h_i=g_i(sigma Bg)_i

are nonnegative and sum to 2Q(B)<=2Cn. Hence at most n/4 coordinates
have h_i>8C, and at most n/4 have D_ii>32C. Their common good set I
has size at least n/2, and ||B_I||op<=32C.

For ANY S contained in I let z_i=g_i 1_(i in S), and let g^S flip those
coordinates. The exact half-energy change is

    Q(B)-H_(sigma B)(g^S)
       =2 sum_(i in S)h_i-2sigma z^T Bz
       <=16C |S|+64C |S|=80C |S|.                       (2)

Fix any subset J of I with cardinality given by (1). Every one of the
2^|J| choices S contained in J then lies in G_eta. This proves the face
claim. A face shatters its freely varying coordinates, and the margin-one
VC argument in `principle_invent_2026_09_07_balanced_valley_augmentation.md`
gives VC dimension at most gamma^2 for any such factorization.

Without invoking the diagonal majorant, the same argument is immediate
under a fixed operator bound ||B||op<=L: use only the bounded-local-field
subset, and replace 80C by 16C+2L.

Thus a statement merely asserting gamma2(G_eta)^2=Omega_eta(n) for
liminf minimizers is mathematically true but structurally generic. The
robust-valley augmentation remains different: it factors reference
CENTERS and separately pays Hamming excursions through an energy
coercivity function. It does not factor the whole thick superlevel,
which necessarily contains the large face above.
