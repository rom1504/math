# Arithmetic-phase catalyst tests: director reconstruction

Date: 2026-09-06. Status: exact finite-field reduction below, independently
reconstructed by the proof auditor from the primary source; its use for a
general optimal seed is OPEN. This is a nonlocal construction, not another
Gaussian response iteration.

## 1. A consequence of the primary Gauss-sum formula

I read Lemmas 4--7 of Kai-Uwe Schmidt, *Asymptotically optimal Boolean
functions*, [author manuscript](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Diskrete_Mathematik/Publications-schmidt/pw.pdf).
For characters of order seven over F_(2^(3s)), Lemma 6 gives normalized
Gauss sums -(-1)^s zeta^(+/-s), where
zeta=(-1+i sqrt(7))/sqrt(8). Frobenius conjugation makes the sign constant
on each of {1,2,4} and {3,5,6}; complex conjugation interchanges them.
The argument of zeta is irrational modulo pi, as in Lemma 7. Even positive
s therefore give a dense phase circle as well as odd s.

The earlier application selected phases approaching one. Here allow any
fixed phase. On the seven cosets, write P_(j,k)=1_{j+k=0}, let E denote
uniform averaging, and let C_(j,k)=chi(j+k)/sqrt(7), with chi the quadratic
character modulo seven and chi(0)=0. On the zero-mean subspace,

    U_theta=cos(theta) P+sin(theta) C.                 (1)

Changing the initial multiplicative character can reverse theta. Direct
character sums give C1=0, C^2=I-E, and PC=-CP. Thus (1) is a real
self-adjoint isometry on that subspace, for EVERY real theta.

For completeness, if f is a real function on the cosets with sum f=0,
extend it to the field by f(0)=0 and constancy on each nonzero coset.
The normalized additive Fourier transform on a nonzero coset j is

    (1/7) sum_k f(k) sum_(l=1)^6 g_l exp(-2 pi i l(j+k)/7),

where g_l is the corresponding normalized Gauss sum. When g_l tends
to exp(i chi(l) theta), the elementary quadratic Gauss identity makes
this expression tend to (1). The transform at zero is exactly zero.
At fixed seven-point f, convergence is uniform over the nonzero cosets.
All of these statements are linear and apply jointly to any fixed finite
list of functions. The total mass of the omitted additive zero is 1/2^(3s).

## 2. Exact lower tests for the existing regularization

Use the already proved bilinear identity

    2R4(B)=sup_a ||H_(4^a) tensor B||_(infinity -> 1)/(4^a)^(3/2)

where R4 is the regularization using only the four-point generator.
The larger regularization R allowing order 144 satisfies R>=R4. Here H_(4^a) is the regular
four-point tensor family; its two-sided signed equivalence to the even
Sylvester transform preserves the bilinear norm. Select even s above,
so 3s is even and the field Fourier matrix has precisely such an order.

For any vectors of seven-point functions f=(f_i), g=(g_i), with every
component in [-1,1], set P0=I-E. This proves the stronger lower test

    R(B) >= (1/2) | E_j sum_(i,k) B_(i,k)
                                  g_i(j) (U_theta P0 f_k)(j) |. (2)

For nonzero-mean f, its constant coset part has an additive-zero output
of order sqrt(field size), but its contribution against a BOUNDED
bilinear output is only O(1/sqrt(field size)) after averaging. Its
nonzero-coset contribution is also O(1/sqrt(field size)). Thus it can
be discarded in this bilinear lower test, not in a claimed L2 isometry.
Output functions g need no centering at all. The finite-order bilinear
maximum dominates every cube-valued input; equivalently independently
round the two cube-valued inputs and take the expected energy. Limit order:
B, theta, f, g fixed first; the even field extension subsequence next.

This is a continuum of exact limiting tests inside the PRESCRIBED outer
family. It is not an enlargement of that family to all orthogonal matrices.

## 3. Falsification and scope

For the weaker restriction to zero-mean inputs, the cube vertices have
one zero, three +1, and three -1 entries (140 possibilities). The stronger
test (2) instead permits all 128 Boolean inputs; the best output is
their ordinary sign response after applying the kernel and seed.
This gives a small finite optimization at each fixed theta. It can discriminate
whether (2) improves a particular seed witness, without solving the full
outer Boolean optimization.

The limiting kernel on the FULL seven-point space is the contraction
U_theta P0, not an orthogonal probability-algebra realization. Its lost
constant part is the leading L2 spike just described. Likewise, composing two
noncommuting U_theta on overlapping interfaces is NOT justified by tensor
closure alone. A successful polar-coupling theorem needs a separate
realization argument.

Higher index 7^e produces phases +/-theta,+/-7theta,...,+/-7^(e-1)theta,
with orientation fixed within each Frobenius character orbit. Lemma 6
alone does not identify the relative signs across distinct orders.
Those signs require another proof or an exact character-sum verification;
the phase values are not independent freely chosen parameters. No claim
that (2) attains T(B), improves the
original upper bound, or proves convergence has been made.
