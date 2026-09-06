# Arbitrary vector controlled-index translations from phase modules

Date: 2026-09-06. Status: a general positive nonlocal realization
theorem. Every map f:F_2^r -> F_2^k, with no degree restriction,
gives a legitimate full Boolean-profile R test with Fourier-index
permutation (a,c)->(a,c+f(a)). This is a simultaneous construction,
not assumed composition of separately available gates.

## 1. Precise signed-module approximation supplied by reflections

For every fixed m and sign function epsilon:F_2^m->{+1,-1},
the binary phase-closure construction supplies normalized allowed
Hadamard outers U_j, Boolean carriers h_j,h*_j, and label maps
Phi_j,Psi_j into F_2^m such that their label laws approach uniform
and, for every index alpha,

    ||U_j[h_j chi_alpha(Phi_j)]
          -epsilon(alpha) h*_j chi_alpha(Psi_j)||_2 ->0.  (1)

The norm here is the uniform probability L2 norm on physical
coordinates. Since m is fixed, the convergence can be taken
uniformly over all indices. Independent signed row/column gauges
are permitted in this bilinear outer representation.

Here are the necessary approximation details. A uniform N=2^m
point reflection Q_N is approximated by exact weighted reflection
quotients Q_p, with p tending to uniform. For any character chi_a,

    ||Q_p chi_a-Q_uniform chi_a||_(L2(p))
              <=2 ||p-uniform||_1.

The same bound holds after the Boolean character gauges used to
flip one selected frequency. Tensor the finitely many required
gauged reflections, then use the sum label map. Telescoping the
tensor product on character inputs bounds the L2 error by the
sum of the individual errors: every operator is an isometry and
every Boolean input/output factor has L2 norm one. The sum-label
law approaches uniform by finite-product continuity. This proves
(1), including its common Boolean carriers and the L2 mode of
approximation; mere convergence of one optimized norm is not
being substituted for a signed-module statement.

By linearity, (1) extends to every fixed real profile on the
finite label space. For a Boolean scalar profile the total error
is at most sqrt(2^m) times the largest character error, by
Parseval and Cauchy--Schwarz. This finite-dimensional factor is
harmless because m is fixed before j tends to infinity.

## 2. Promotion of an arbitrary vector-valued control function

Fix ANY map f:F_2^r -> F_2^k. Apply (1) with m=r+k and phase

    epsilon(a,s)=(-1)^(s dot f(a)).                        (2)

Split the old labels as Phi=(Phi_1,Phi_2), Psi=(Psi_1,Psi_2),
with dimensions r and k. Add two independent k-bit physical
selector variables s,t, and tensor U_j with the normalized
Walsh operator H_(2k). The new input carrier and labels are

    h'(x,s,t)=h(x)(-1)^[s dot (t+Phi_2(x))],
    Phi'(x,s,t)=(Phi_1(x),s).                             (3)

For the input character indexed by (a,c), its full Boolean
exponent beyond the old carrier is

    s dot t+c dot s+a dot Phi_1+s dot Phi_2.

At output selector coordinates (u,b), summing over t forces
s=b. Applying the old component identity at index (a,b)
then gives, up to the L2 error in (1),

    h*(y)(-1)^[b dot u+c dot b+a dot Psi_1
                                +b dot Psi_2+b dot f(a)].

Define

    h*'(y,u,b)=h*(y)(-1)^[b dot (u+Psi_2(y))],
    Psi'(y,u,b)=(Psi_1(y),b).                             (4)

The new component identity is exactly

    (H_(2k) tensor U_j)[h' chi_(a,c)(Phi')]
              ~=h*' chi_(a,c+f(a))(Psi').                 (5)

Both new label laws approach uniform: the first label block is
a marginal of an asymptotically uniform old label, and the
second block is an independent uniform selector. No further
direct-sum equidistribution or zero-channel repair is needed.

The new error for a fixed (a,c) has squared L2 norm equal to
the average over b of the old squared component errors at
(a,b). Thus the maximum error does not increase. The calculation
remains valid for the limiting reflection modules, not only
exact finite bent identities.

## 3. The resulting full profile operator and norm inequality

The index map

    sigma_f(a,c)=(a,c+f(a))                               (6)

is an involutive permutation on F_2^(r+k). Equation (5), the
uniform label laws, and finite Fourier expansion prove that
for every symmetric seed B and every Boolean profile array F,

    R(B)>=(1/2) E_z ||[T_f F](z) B||_1,
    T_f=H_(r+k) P_(sigma_f) H_(r+k).                      (7)

All profile columns may be unbalanced. The common carrier is
what prevents a missing constant channel. The outer remains
in the prescribed regularized family because the added Walsh
factor has even dimension 2k; its standard form is equivalent
by independent Boolean row/column gauges and permutations to
an R4 tensor power.

For k=1 there are 2^(2^r) possible control functions, and hence
that many distinct actual index permutations. This is already
exponential in the label-space size. In particular, this theorem
is not confined to quadratic inverse pencils or bounded-degree
cubic parameter families.

## 4. Exact norm reduction: one promoted gate adds no phase-norm power

There is a stronger structural limitation than its restricted
matching graph. In the PRIMAL target coordinate y in F_2^k,
T_f is block diagonal. On the y-block its action on x is

    U_y=H_r diag_a[(-1)^(y dot f(a))] H_r.                (8)

This follows by Fourier-transforming the target translation
c->c+f(a): its eigenvalue at the target character y is
(-1)^(y dot f(a)). Consequently, if

    N(U;B)=(1/2) max_(F Boolean) E ||(UF) B||_1,

then independent optimization on the primal blocks gives the
EXACT identity

    N(T_f;B)=2^(-k) sum_y N(U_y;B).                      (9)

Every U_y is already an available binary Fourier-sign module.
Thus a single controlled translation cannot improve on the
supremum of those phase-only seed tests, even though its formal
index alphabet has exponentially many possibilities. This
distinction is essential for the original-problem application.

## 5. Scope: translations inside fibers are not arbitrary matchings

Within each fixed a-fiber, sigma_f translates the target block
c by one vector f(a). It does NOT arbitrarily permute that
fiber. The exponential number of scalar controlled flips only
chooses whether to swap each of many disjoint pairs.

Reversible circuits can use networks of such controlled
translations, but no closure under arbitrary overlapping
networks has been proved here. In particular, this theorem
does not yet supply the arbitrary data-dependent permutations
used in the Fourier-Gaussian matching theorem for T(B).

The positive advance is the exact controlled-index realization
(7) with arbitrary truth tables and full L2 signed-module
control. A network or other genuinely non-block-diagonal
construction is needed to obtain norm power beyond the
phase-only class. The original convergence question is not
settled by this theorem.

The exact integer selector verifier
`computations/resumed_convergence_controlled_promotion_verify_2026_09_06.py`
passes all 16 scalar control truth tables for r=2,k=1, all 256
vector tables for r=2,k=2, and 17 additional r=3,k=2 examples.
