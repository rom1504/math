# Wave 52 main audit: parent-cap pairing removes the quadratic noise

Status: independently derived exact reduction, pending agent cross-audit and
an abundance/profile theorem.

For a fixed local state, write as in Wave 51

    Z(w)=L(w)+Q_T(w),   Z(-w)=-L(w)+Q_T(w),
    b=q-e.

The parent energy cap says `Z(w)<=b` for every outside completion. Hence on
each projective pair

    Q_T(w)+|L(w)| = max{Z(w),Z(-w)} <= b,

and therefore

    min{Z(w),Z(-w)} = Q_T(w)-|L(w)| <= b-2|L(w)|.

For every `r>=0` this proves the exact structural lower bound

    P_w{Z(w)<=-r}
      >= (1/2) P_w{|L(w)| >= (b+r)/2}.                  (A)

Because `L` is symmetric, the right side is at least the corresponding
one-sided linear tail (with harmless endpoint conventions). Unlike the Wave
51 union-bound estimate, (A) pays no separate probability for controlling
the dependent quadratic part: the global parent cap controls it pointwise on
the same rare linear-tail event.

Let `H=n^(3/4-c)` and use Montgomery-Smith's primary theorem for Rademacher
sums: for universal constants `c0,C0>0`,

    P{L >= c0 K_{1,2}(beta,sqrt(H))} >= c0 exp(-C0 H),

with

    K_{1,2}(beta,sqrt(H))
      asymp sum_{i<=H} beta_i^*
            +sqrt(H)(sum_{i>H}(beta_i^*)^2)^(1/2).      (B)

Consequently the exact profile condition

    b+r <= c1 K_{1,2}(beta,sqrt(H))                    (C)

gives a saved conditional completion tail at `-r`. If `r` is not confined to
`O(n^(3/2-c))`, this can escape the local Hanson--Wright recurrence wall.

The condition is nonvacuous at the right scale. If `s^2=sum beta_i^2` and
`B_inf=max |beta_i|`, the Holmstedt expression in (B) obeys

    K >= c min{sqrt(H) s, s^2/B_inf}.                  (D)

Indeed, if `E` is the squared mass of the largest `H` coefficients, their
`l1` mass is at least `E/B_inf`, and the remaining term is
`sqrt(H)sqrt(s^2-E)`; this concave function is minimized at an endpoint.
Here `B_inf<=2m=O(n)`. Macroscopic cross variance `s^2>=eta n^(5/2)` thus
gives `K>=c_eta n^(3/2)`, large enough in principle for genuinely macroscopic
negative margins.

What remains open is decisive: exact-minimizer averaging must force saved
local-state mass on which the constant in (C) beats `b+r`, or give a
complementary route when the cross K-functional is smaller. Neither (A) nor
(B) proves that abundance. A profile with `K << b+r` on every far-margin
state falsifies this implementation without falsifying the full annealed CDF.

Primary source audited: S. J. Montgomery-Smith, "The Distribution of
Rademacher Sums", Proc. AMS 109 (1990), 517--522, theorem and Holmstedt
formula; author preprint https://stephenmontgomerysmith.github.io/preprints/tail.pdf.
