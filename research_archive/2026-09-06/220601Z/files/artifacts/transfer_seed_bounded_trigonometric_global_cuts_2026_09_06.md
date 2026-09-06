# Bounded trigonometric responses retain all global Boolean cuts

2026-09-06. Verified exact finite-cube lemma from the rich-feedback
extension audit; the director and adversarial agent independently checked
its matrix remainder and cover factorization. This avoids polynomial-in-coherent-variable approximation
for GLOBAL cuts. It does not yet prove small proper cuts after a flat
transport, a mixed nuclear theorem, or the rich returned-query closure.
The separate verified high-degree transport theorem is in
`transfer_seed_trigonometric_high_degree_transport_2026_09_06.md`.

## 1. A useful entrywise exponential matrix bound

For every real rectangular matrix J and real t,

```
|| [exp(it J_ab)-1]_(a,b) ||op
 <= |t| ||J||op + (t^2/2) ||J||op^2.                    (1)
```

Indeed subtract itJ. The remainder is bounded entrywise in absolute
value by `(t^2/2) J_ab^2`. Every absolute row and column sum of the
remainder is at most `(t^2/2)||J||op^2`; the Schur test gives (1).
This is a dimension-free quadratic bound, NOT an exponential bound in
the possibly polylogarithmic norm of J. Complex outputs cause no issue;
the arguments J_ab are real.

## 2. Exact cover formula for higher Boolean derivatives

Let W_i(S) be a real polynomial on independent signs. For distinct seed
sets U use `Delta_U=product_(a in U) Delta_a`, where
`Delta_a W=(W(S)-W(S^a))/(2S_a)`. For any T,

```
W_i(S^T)=W_i(S)+sum_(empty!=V subset T)
                    (-2)^|V| chi_V(S) Delta_V W_i(S).
```

Define, for each nonempty V,

```
A_i,V=exp(it (-2)^|V| chi_V Delta_V W_i)-1.
```

Then for `s=|U|>=1` the following is EXACT:

```
Delta_U exp(it W_i)
 =(-1)^s 2^(-s) chi_U exp(it W_i)
   sum_(F subset (2^U\{empty}), union_(V in F) V=U)
                                    product_(V in F) A_i,V.    (2)
```

To prove it, insert the first formula into
`2^-s chi_U sum_(T subset U)(-1)^|T| exp(it W_i(S^T))`.
Expand the FINITE product of `1+A_i,V`. The sum over T vanishes unless
the selected family covers U, in which case it is (-1)^s.
There is no infinite Taylor expansion, no exponential moment assumption,
and no passage through polynomials in the numerical value W_i.

## 3. Operator bounds and all global original Walsh cuts

Let `J_W^(r)(S)` have rows i and columns the increasing r-tuples V,
with entries Delta_V W_i. Fix s. The matrix with entries A_i,V for
|V|=r<=s is the entrywise exponential from (1), applied to
`(-2)^r J_W^(r) diag(chi_V)`. Hence

```
||A^(r)||op <= 2^r |t| ||J_W^(r)||op
              +2^(2r-1)t^2 ||J_W^(r)||op^2.             (3)
```

For each fixed relative cover family in (2), take the tensor product
of these A^(r) matrices, compress their output roots to equality, and
embed a full U-tuple into the subtuples prescribed by the cover. The
column embedding is injective precisely because the selected subsets
COVER U. Thus it is an isometry. Repeated labels between different
factors are retained by this embedding; no diagonal restoration is
assumed. The root phase exp(itW_i) and column character chi_U have norm
one. Consequently (2) bounds `||J_exp(itW)^(s)||op` by a fixed finite
sum of products of the right sides of (3), times 2^-s.

If, for every fixed r and p, the original derivative matrices obey

```
|| ||J_W^(r)||op ||_(L^p) <= C_(r,p) log(n+1)^C_(r,p),
```

then the SAME form of bound holds for exp(itW), for every fixed s,p.
Hölder is applied only to finitely many factors. A finite real
trigonometric polynomial of a fixed-dimensional vector W has the same
property, by applying this argument to each real linear combination
t dot W. Its amplitude can remain uniformly bounded even when the
coherent W has no exponential moment.

Fourier orthogonality, exactly as in Section 4 of
`continued_feedback_all_global_walsh_cuts_2026_09_06.md`, identifies the
global degree-q Walsh cut with row `(root, q-s residual marks)` and
column s marks as the L2 matrix of the degree-(q-s) component of
J_C^(s). Therefore every fixed positive original Walsh degree of this
bounded trigonometric response has all global cuts polylogarithmically
bounded. Transports by a fixed bounded-op matrix preserve the derivative
operator estimates exactly.

This extension is for fixed trigonometric stages. Passing an arbitrary
bounded response to an operator-norm global-cut statement by averaged
L2 approximation would be invalid. It is only safe inside an argument
whose final tested quantity is continuous in the appropriate averaged
L2 norm.

## 4. Why it does not finish the rich feedback proof

A bounded continuous function of finitely many literal query variables
can be approximated on a large compact box by a finite trigonometric
polynomial with controlled global amplitude (periodize a continuous
extension and use Fejer approximation). Tightness controls the part
outside the box. This avoids the genuine moment-determinacy obstruction
to polynomial-in-W approximation.

However the coherent boundary-graph small-cut theorem used a finite
polynomial source in primitives of bounded original degree. The
trigonometric source has infinite original degree, and its fixed-degree
projections do not automatically inherit that selected-diagram structure.
All GLOBAL cuts alone do not prove small proper cuts of B C, nor full
probe contractions into B C. In particular one cannot replace the missing
full-return condition by orthogonality to individual W primitives.

The lemma above therefore removes one concrete analytic obstruction for
the bounded-trigonometric route but leaves the actual higher-degree
returned-query estimate open. No rich-frame cap gain or new original
lower constant is claimed from it.

## 5. Finite checks

`computations/transfer_seed_trigonometric_cut_checks_2026_09_06.py`
passed 192 rectangular entrywise-exponential operator tests, 3200 scalar
instances of the exact Boolean cover formula, and 96 modified derivative
matrix bounds. Its explicit finite-cube computations use a degree-three
real polynomial W, retain all partial-label characters, and compare the
cover expansion with direct seed flips. These are algebra checks only.

## 6. Dimension-free bounded-trigonometric circuits

The director pointed out a stronger useful specialization. A fixed
finite circuit consisting only of real bounded trigonometric coordinate
gates, bounded-op real linear transports, and arbitrary deterministic
offsets has DIMENSION-FREE pointwise bounds for every fixed derivative
matrix order. Start with J_seed^(1)=I and higher derivatives zero.
Transports multiply the derivative matrices by their fixed operator
norms. At a trigonometric gate, (2)--(3) express each derivative norm
as a fixed finite polynomial of previously bounded derivative norms.
No coordinate maximum or exponential moment is used. Deterministic
offsets contribute zero derivatives and only unit-modulus root phases.

This is distinct from an initial unbounded polynomial circuit, whose
random derivative operator norms may be polylogarithmic. Neither version
alone supplies small proper cuts after the next transport.
