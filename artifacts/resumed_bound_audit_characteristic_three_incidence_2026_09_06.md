# Independent audit: characteristic-three Cayley incidence bound

Date: 2026-09-06. Outcome: PASS for the final sharp version of
`resumed_convergence_characteristic_three_incidence_gain_2026_09_06.md`.

For every sequence of hollow symmetric additive-Cayley signings on
F_3^r, q=3^r tending to infinity, the proved lower bound is

    liminf Q(A)/q^(3/2) >= 2080/(9 sqrt(269441))
                         =0.44523467985944279335138857059...

The decimal 0.4452346798594428 is rounded, not a rigorous truncated
lower endpoint. There is no flat-spectrum or row-sum assumption.

## 1. Spectral factors and exact second moment

The three-point Boolean square wave has constant Fourier weight 1/9
and total weight 8/9 on its two conjugate nonzero frequencies. Hence
with K=2Q/q^(3/2), t=lambda(0)/sqrt(q), and
h(a)=(8lambda(a)+lambda(0))/(9sqrt(q)), the exact Boolean tests give
|h(a)|<=K and |t|<=K. Evenness makes h a function on projective points.

Hollowness and the signing condition imply sum lambda=0 and
sum lambda^2=q(q-1). Expanding over all nonzero characters gives

    sum_(a!=0) (8lambda(a)+lambda(0))^2
       =64q(q-1)+(q-81)lambda(0)^2.

Dividing by 81q(q-1) proves the stated exact moment identity, with
gamma_q=(q-81)/(81(q-1)). No spectral independence is used.

## 2. The sharp projective-plane inequality

The final strengthened claim 3N_0+N_2>=3 in PG(2,3) checks
independently. If no line is monochromatic, each color is a blocking
set with no full line. Three or fewer points meet at most twelve lines.
For a four-point blocking set, line/incidence/pair counts force a
negative number of two-point lines. For five points they force exactly
three three-point lines. Their intersections within the five-point set
would have to be pairwise at most one. Two triples already cover the
set and meet once; no third triple can meet each in at most one.
Thus both colors have at least six points, so their sizes are six and
seven. The line-pair identity then gives N_2=3. If N_0>=1, the desired
inequality is immediate.

Averaging over all three-dimensional subspaces is legitimate because
every projective line is contained in equally many of them. With
L=(q-1)/2 and total color sum M, the exact pair average is
r_pair=(M^2-L)/(L(L-1))>=-1/(L-1). Every point pair lies on one line,
giving p_2=3p_0-3r_pair. Combining with 3p_0+p_2>=3/13 yields

    p_0>=1/26-1/(2(L-1))=(q-29)/(26(q-3))=:rho_q.

All inequalities have the required direction. Colors at zero h may be
chosen arbitrarily; h=s|h| remains true there.

## 3. The nine-point Boolean test and deficit direction

Any five-plus/four-minus profile on F_3^2 has squared mean 1/81 and
remaining Fourier mass 80/81. Averaging invertible linear coordinate
changes assigns mass 20/81 to each of its four projective frequencies.
The quotient pullbacks are actual Boolean vectors on F_3^r. Their
average normalized energy is exactly

    (t+20 sum_line lambda(a)/sqrt(q))/81
             =(5/18) sum_line h(a)-t/9.

It has absolute value at most K. On a monochromatic line of sign s,
putting r(a)=K-|h(a)| gives

    sum_line r(a)>=(2/5)(K-s t)>=(2/5)(K-|t|).

The final factor is nonnegative by the constant Boolean test. Uniform
line incidence and nonnegativity on all other lines give
E r>=rho_q(K-|t|)/10. The square deficit obeys

    D=K^2-E h^2=E[(K-|h|)(K+|h|)]>=K E r.

Substitution therefore yields

    64/81 <= K^2(1-rho_q/10)
                  +rho_q K|t|/10-gamma_q t^2.

For q>81, gamma_q is positive. Completing the square bounds the final
two terms above by rho_q^2 K^2/(400 gamma_q). This proves the finite
inequality without requiring a separate bound on t; dropping its
constraint only weakens this upper estimate.

## 4. Exact arithmetic and replay

The limiting denominator factor is exactly

    1-1/260+81/(4*260^2)=269441/270400.

Since sqrt(270400)=520 and Q/q^(3/2)=K/2, the constant is exactly
2080/(9sqrt(269441)). The earlier weaker parity-only constant is
superseded by this final strengthened one.

The updated verifier
`computations/resumed_convergence_characteristic_three_incidence_verify_2026_09_06.py`
was independently replayed. It passed all 256 profiles up to global
sign on F_3^2 and all 8192 projective-plane point colorings, including
the sharp inequality and the two equality types. Fraction arithmetic
independently checked the final rational factor. These computations
corroborate, but do not replace, the elementary proof above.

The theorem excludes low-cap characteristic-three Cayley targets below
this constant. It neither applies to arbitrary hollow signings nor
establishes convergence of the original minimum sequence.
