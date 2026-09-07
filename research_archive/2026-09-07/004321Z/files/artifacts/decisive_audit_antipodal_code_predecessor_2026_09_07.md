# Exact-family coding predecessor: source status and asymptotic scope

Date: 2026-09-07. Bibliographic and abstract audit only. No inaccessible
full-text theorem or third-party enumeration is imported.

The publisher's [official Volume 78 record](https://combinatorialpress.com/um/vol78/)
and the [author's publication list](https://www.ece.uvic.ca/~mesmaeil/publications.html)
confirm M. Esmaeili and A. Zaghian, *On covering radius of a family of
codes C_m union (1+C_m) with maximum distance between C_m and 1+C_m*,
Utilitas Mathematica78(2009),151--158. Neither accessible page supplied
the full paper.

The paper's [accessible abstract](https://www.researchgate.net/publication/266941646_On_covering_radius_of_a_family_of_codes_C_m_1C_m_with_maximum_distance_between_C_m_and_1C_m)
defines a parity-check matrix whose columns are the distinct weight-two
binary vectors of length `m`. Its row span is exactly the cut code of
`K_m`: row `i` is the incidence vector of the star at `i`. Adjoining the
all-one word therefore gives the exact antipodal cut code corresponding
to the absolute quadratic objective, not merely a related cycle code.

The abstract conjectures, for `m>=7`, equality between the covering radius
of this antipodal cut code and that of the ORDINARY cut code at order
`m-1`; it reports verification only for orders seven through ten. The
same abstract states that the all-one word is furthest from the ordinary
cut code. Its distance is the number of complete-graph edges minus the
largest cut size. Thus the proposed formula would imply

```math
 M_{2r+1}=3r,\qquad M_{2r}=3r-2.                          (1)
```

Indeed `rho(C_cut(K_s))=binom(s,2)-floor(s^2/4)` under the stated
ordinary-code assertion, and `M_n=binom(n,2)-2rho(C_antipodal(K_n))`.
Equation (1) has LINEAR growth. It is incompatible with every positive
universal `n^(3/2)` lower bound, and already disagrees with the independently
established `M_11=17` (it predicts15). Therefore this source is a direct
historical predecessor, but its conjecture cannot supply an asymptotic
convergence theorem for the original sequence.

A newly indexed public repository, CurtisAccelerate/antipodal-cut-code-k11,
also states the finite failure. Its search result led to the primary
abstract, but none of that repository's code, numerical claims, or novelty
assertions are needed or imported here. We did not retrieve the2009 full
text and do not claim to have audited its proofs beyond the algebraic
consequences of the displayed abstract.
