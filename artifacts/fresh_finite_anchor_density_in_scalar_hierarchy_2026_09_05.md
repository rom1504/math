# Finite ancestor anchors are dense in the scalar marked hierarchy

Date: 2026-09-05. This is a precise interpretation of the scalar-family
supremum, not a new bound on the original Boolean optimum.

Let `Z_T` be the independent Gaussian coordinates indexed by the old
rooted marked trees, and let `U h_T=Z_T`, where `h_T` is the normalized
Hermite product of the child coordinates of `T`. The edge has `h_edge=1`.

Consider a finitely supported unit first-chaos variable

`V=sum_T v_T Z_T`, `sum_T v_T^2=1`.

Suppose its support includes a nonedge tree. Choose a support tree `T*`
of maximal height with `v_(T*)!=0`. Take as anchors all other support
trees, all proper child-descendant trees of `T*`, and the closure under
child descendants of this finite set. This is a finite ancestor-closed
anchor set and does not contain `T*`: any occurrence of `T*` as a proper
descendant would require a larger-height support tree.

Set `s=|v_(T*)|`, put the other support coefficients in the anchor vector
`rho` (zero on extra anchors), and define

`h(anchors,z)=sign(v_(T*)) h_(T*)(anchors)`.

This response is independent of the innovation argument `z`. It is
centered, has Gaussian squared norm one, and is orthogonal to every
anchor input feature `h_T`, since those normalized Hermite features are
orthogonal and `T*` is not an anchor. Its conditional innovation
Dirichlet energy is exactly zero.

Thus it is an admissible response in the general finite-anchor
fixed-point theorem. The innovation map is constant:

`z -> U h = sign(v_(T*)) Z_(T*)`.

Its unique fixed point gives precisely

`rho dot anchors+s z = V`.

If a version of the theorem requires `0<s<1`, variables supported on only
one coordinate are obtained by an arbitrarily small addition of another
coordinate. Edge-only variables are treated in the same way. Every unit
first-chaos variable is an `L2` limit of finitely supported unit variables.

For a fixed central threshold, masks converge in `L2`, their `U` images
converge by isometry, and `E|UH|(1-H)` converges by Cauchy--Schwarz.
Consequently the supremum over the general finite ancestor-anchor
construction equals `C_scalar` from
`fresh_uniform_scalar_hierarchy_escape_2026_09_05.md`.

The response class here is the general admissible class in the anchor
theorem; it is not restricted to the finite Hermite resolvent ansatz used
in a numerical optimization. Hence the uniform strict-escape theorem
separates the original universal bound from the entire scalar finite-
anchor architecture, not merely from one optimizer's parameterization.
There is no corresponding assertion about arbitrary multivariate final
masks or the enlarged unmarked hierarchy.
