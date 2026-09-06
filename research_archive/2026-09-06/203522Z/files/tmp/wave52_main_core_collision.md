# Wave 52 main audit: core-collision sandwich

Status: independently derived exact theorem, pending agent cross-audit.

Let `N=C(n,m)`, let `F_z` have size `r_z=a_z N`, and for every `ell`-core
`R` put

    c_z(R) = |{S in F_z : R subset S}|,
    mu_z(R) = c_z(R)/r_z.

Double counting ordered selector pairs through their common cores gives

    p_ell(z) = sum_R c_z(R)^2 / (r_z d_ell)
             = a_z C(n,ell)/C(m,ell)^2 sum_R mu_z(R)^2,

where `d_ell=C(m,ell)C(n-ell,m-ell)` and
`sum_R mu_z(R)=C(m,ell)`. Cauchy and `0<=mu_z(R)<=1` therefore prove

    a_z <= p_ell(z) <= a_z/rho_ell,
    rho_ell = C(m,ell)/C(n,ell) = (m)_ell/(n)_ell.       (A)

After weighting centers by `1_C a_z^2`, (A) becomes the exact aggregate
collision sandwich

    r_C <= P_ell(C) <= r_C/rho_ell,
    r_C = E[1_C a_z^3]/E[1_C a_z^2] <= M.               (B)

At `ell=1`, this recovers (10.1264). At `ell=2`, the common-core second
eigenvalue is

    lambda_2(2) = 2(n-m)(n-m-1) /
                  [m(m-1)(n-2)(n-3)].

Consequently, provided `D_C>0`, the mere non-strict threshold

    P_2(C) >= lambda_2(2)                                (C)

already implies, without using the small excess in (10.1261),

    M >= r_C >= rho_2 P_2(C)
      >= 2(n-m)(n-m-1)/[n(n-1)(n-2)(n-3)]
       = Theta(n^-2)                                    (D)

uniformly on a compact fixed-density window. This is much stronger than the
required `exp{-O(n^(3/4-c))}` selector degree. Equality in (C) is sufficient.

Thus the controlled-cap coarea target can be sharpened again: a box witness
plus a cap `R0<=R<=C(R0+n^2)` satisfying (C) proves (10.795) and convergence.
The result does not prove (C); singleton favorable families have exponentially
small `P_2<lambda_2`. Its value is that no quantitative spectral-excess gap is
needed at level two.

More generally, if `ell=O(n^(3/4-c))`, then `rho_ell=exp{-O(n^(3/4-c))}`
on a compact density window, while `lambda_2(ell)` is at worst polynomially
small for `ell>=2` in the sublinear range. Hence `P_ell>=lambda_2(ell)` also
gives saved degree directly through (B). The exact level-two statement (D)
is the cleanest target and avoids asymptotic qualifications.
