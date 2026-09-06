# A uniform finite-query escape and a strict query-depth hierarchy

Date: 2026-09-06. The coefficient-peeling proof below was independently
reconstructed by the director and response agent. The response agent
supplied the shared elementary inverse modulus and the branchwise escape
corollary. Scope: the previously audited finite paired-query Haar class
and its fixed-rule exact-incoherent-involution transfer. This is not an
original-signing universal lower bound.

## 1. Statement

Fix m. A paired history consists of at most m orthonormal pairs
(q_j,p_j), where p_j=N_j is a fresh independent standard Gaussian and
q_j and earlier pairs are measurable before N_j is exposed. Arbitrary
finite initial side information is independent of these fresh normals.
The limiting involution swaps each pair. Assume the usual joint oddness
for coordinate policies and spin registers.

There is an explicit positive number eta_m, independent of all the
history functions and side-information dimension, with the following
property. For every odd |f|<=1, set

    c=E[H f],  k=H^T Jc,
    Edef=1-||c||^2,
    delta=E[|k|-f k].

Then Edef and delta cannot both be at most eta_m. Consequently every
such policy has a depth-only positive actual-energy improvement using
at most one additional query. In the population variational closure the
gain is at least

    Delta_m = min(eta_m^2, G(sqrt(eta_m))^2)/8 > 0,
    G(x)=E|sqrt(1-x^2)+xN|-sqrt(1-x^2),  0<=x<=1.             (1)

Fixed Lipschitz/gate approximants realize, for example, half this margin
with the same query budget. No constant is claimed uniform in m.

## 2. Peeling one newest pair

Write the newest two coefficients of c as (a,b), so

    H^T c=M+bN,    k=L+aN,

with M,L measurable in the filtration before this N. All coefficients
are deterministic. Bessel's identity gives

    ||f-H^T c||_2^2=E f^2-||c||^2 <= Edef.

Since |f|<=1, the squared distance of H^T c from [-1,1] is at most
Edef. The convex even function x -> (|x|-1)_+^2 has its Gaussian-shift
average minimized at zero. Conditioning on the past therefore gives

    T(|b|)<=Edef,    T(x)=E[(|xN|-1)_+^2].                    (2)

Let fbar=E[f|past]. Orthogonal projection of the preceding approximation
identity off the past yields

    ||f-fbar||_2 <= sqrt(Edef)+|b|.                            (3)

This remains valid after later pairs have been deleted, even if f still
depends on those later normals: the error term is only being projected
in L2, not assumed adapted.

Because |fbar|<=1, ||L||_2^2=||c||^2-a^2<=1-a^2, and N is independent
of L, the convex-in-L^2 Jensen bound gives

    delta >= E[Psi(L,a^2)-|L|]-|a| ||f-fbar||_2
          >= G(|a|)-sqrt(Edef)-|b|.

Thus

    G(|a|) <= delta+sqrt(Edef)+|b|.                            (4)

Both T and G are continuous, strictly increasing on [0,1], and vanish
only at zero. For G, increasing its noise parameter while decreasing
sqrt(1-x^2) strictly increases the Gaussian absolute-value excess.

Delete this pair from H and c, keeping f unchanged. If r^2=a^2+b^2,
the new quantities satisfy exactly or by Cauchy--Schwarz

    Edef_new=Edef+r^2,
    delta_new<=delta+2r.                                     (5)

Indeed ||k-k_new||_2=r and |f|<=1. Also r<=1. Repeating deletion will
end with Edef_final=1, irrespective of f.

## 3. Fully explicit elementary modulus

For 0<x<=1, restrict the Gaussian integrals in T and G to
2/x<=N<=3/x. The respective squared excess and positive excess are at
least one there. Consequently both functions obey

    T(x), G(x) >= 2 phi(3/x)
                    >= (2/3) exp(-9/(2x^2)).                 (6)

The last coefficient uses sqrt(2/pi)>2/3. Define omega(0)=0 and

    omega(r)=min(1, 3/sqrt(2 log(2/(3r))))   for 0<r<2/3,
    omega(r)=1                             for r>=2/3.

Then T(x)<=r or G(x)<=r implies x<=omega(r). This function is
continuous, nondecreasing, and tends to zero at zero. If
w>=max(Edef,delta), equations (2)--(5) imply

    |b|<=omega(w),
    |a|<=omega(w+sqrt(w)+omega(w))=:A(w),
    max(Edef_new,delta_new)<=P(w):=w+4 A(w).                   (7)

Here r<=sqrt(2)A(w), r^2<=2A(w)^2, and A(w)<=1. The loose constant
four conveniently controls both updates. P is continuous and
nondecreasing, with P(0)=0.

An explicit positive starting threshold can be given without any
inverse functions or a numerical search. Put

    h(x)=(2/3) exp(-9/(2x^2)),
    z(w)=h(w/8)/3,
    tau(w)=min(w/4, z(w)^2, h(z(w)))  for 0<w<=1/2,
    eta_m=tau composed m times, evaluated at 1/2.             (8)

Every eta_m is positive, at most 1/2, and decreases with m. If
v=tau(w), then v<=w/4, sqrt(v)<=z(w), and omega(v)<=z(w).
Since v<=sqrt(v),

    v+sqrt(v)+omega(v)<=3z(w)=h(w/8),
    A(v)<=w/8,
    P(v)<=3w/4 < w.

Therefore P composed m times at eta_m is <1/2 when m>=1; for m=0 the
initial value is 1/2. If both Edef and delta were <=eta_m, repeated
peeling of r<=m pairs would force the final Edef to be at most
P composed r times at eta_m <=1/2, contradicting its exact value one.
This proves the uniform exclusion in Section 1. The explicit constants
are extraordinarily small; no claim of practical numerical usefulness
is made.

## 4. Escape without optimizing the old moment body

The actual half-energy of a policy is e(f)=c(f)^T Jc(f)/2.
There are two exhaustive cases.

If delta>eta_m, put v=sign(k), using zero at k=0, and move within the
old history to f_theta=(1-theta)f+theta v. Its energy satisfies

    e(f_theta)-e(f)>=theta delta-2theta^2.

Taking theta=eta_m/4 yields gain at least eta_m^2/8. A fixed softsign
approximates this gain by its pointwise regret bound, so zero atoms of
k cause no problem. No new matrix query is needed in this branch.

Otherwise Edef>eta_m. Use fresh independent even/odd seed gates to
round f to a Boolean u with E[u|old history]=f. This preserves c and
its actual limiting energy. One genuine matrix query now has response

    Bu=k+sqrt(Edef)N_new,

with independent fresh N_new. Its actual stability gap is at least
G(sqrt(eta_m)). The feasible damped best-response endpoint with step
G(sqrt(eta_m))/4 gains at least G(sqrt(eta_m))^2/8. Independent final
Boolean rounding preserves normalized energy in the hollow or
vanishing-diagonal involution class.

Every map is approximated at a fixed tolerance before taking the matrix
limit. The strict margin permits fixed finite Lipschitz approximants
with at least Delta_m/2 gain, rather than only an iterated-limit formal
policy. Side-information gates are free but finite; the second branch
adds exactly one physical query. Exposing the final energy for analysis
does not add an input to the output rule.

## 5. A strict hierarchy of query budgets

Let C_m be the supremum of limiting positive actual energy over fixed
equivariant algorithms exposing at most m nondegenerate query pairs,
with arbitrary finite independent seed information. Zero-residual
queries can be removed in normalized L2 before a fixed Lipschitz
continuation, as proved in the separate finite-history extension audit.

Take policies with values tending to C_m and apply Section 4. The new
query is nondegenerate in the second branch. Its appended residual is
measurable before its fresh Gaussian reply, so it is an admissible
(m+1)-pair history. Approximation errors may be made arbitrarily small
without changing query count. Hence, in terms of the population margin
(1),

    C_(m+1) >= C_m + Delta_m > C_m.                            (9)

No global optimization or compactness over all histories is required
for this deduction. The previously proved fixed-GFOM ceiling gives
C_m<=sqrt(15)/8, and therefore C_m<=sqrt(15)/8-Delta_m at every finite
m. The increasing sequence C_m has a supremum, but neither its value
nor any useful upper convergence rate is determined here. In particular
this hierarchy does not prove convergence of M_n/n^(3/2), and it does
not identify the Boolean optimum of Hadamard matrices.
