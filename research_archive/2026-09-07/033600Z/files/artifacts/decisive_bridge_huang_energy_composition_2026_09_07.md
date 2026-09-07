# Boolean energy composition for Huang's signed hypercube

Scope: this is a sparse, weighted diagnostic, not a theorem about the dense
signing optimum. The four-ary energy construction below was proposed by the
director and independently verified here.

Put H_0=0 and H_d=[[H_(d-1),I],[I,-H_(d-1)]], and let Gamma_d
be diagonal vertex parity. Then H_d^2=dI and Gamma_d H_d=-H_d Gamma_d.
Write Q_d=max_f |f^T H_d f|/2 over Boolean f, and e_d=Q_d/2^d.
Parity switching reverses the energy, so the absolute maximum equals the
positive maximum. The spectral upper bound is e_d<=sqrt(d)/2.

## 1. Exact four-block composition

For any positive d_1,...,d_4,

    e_(d_1+...+d_4) >= (e_d1+...+e_d4)/2.                 (1)

Choose positive-maximizing Boolean functions f_i on the four blocks. Let
h=(1,1,-1,1,1,-1,-1,1,1,1,-1,-1,-1,-1,-1,1) in binary
vertex order; direct integer multiplication gives H_4 h=2h. Write
H_4=sum_j C_j, where C_j flips bit j with the preceding-parity sign.
The Clifford anticommutator {C_j,H_4}=2I gives

    <h,C_j h> / 16 = 1/2                              (2)

for every j. This is stronger than just knowing their sum.

For x=(x_1,...,x_4), let p_i be the parity bit of x_i and set
F(x)=prod_i f_i(x_i) h(p_1,...,p_4). Under the natural block order,
H_d=sum_i Gamma_(d1) tensor ... tensor Gamma_(d(i-1)) tensor H_di
tensor I. Each summand flips just p_i. The conditional expectation of
f_i(x_i)(H_di f_i)(x_i), given either parity, equals its unconditional
expectation 2e_di: every undirected edge contributes the same signed
product once to each parity side. All parity bits are independent fair
under uniform x, so the i-th summand in E[F H_d F] is exactly
2e_di times the directional expectation in (2). Thus the resulting
energy density is sum_i e_di/2, proving (1).

This argument requires no flatness of the inner f_i.

## 2. Explicit all-d lower bound

Starting with the d=1 positive eigenvector, build a four-ary tree. A full
depth-r tree has d=4^r leaves and energy density 2^(r-1). Split any j of
its leaves once, 0<=j<=4^r. The new dimension is d'=4^r+3j, and (1) gives

    e_d' >= (4^r+j)/(2*2^r) = (d'+2*4^r)/(6*2^r).     (3)

For 4^r<=d'<=4^(r+1), the ratio of the last expression to sqrt(d')
is minimized at d'=2*4^r and is sqrt(2)/3.

One may ignore extra coordinate factors: if f uses d' coordinates and is
constant on the remaining d-d' coordinates, its old energy density is
unchanged and the extra summand averages to zero because its old-block
parity factor has mean zero. Therefore e_d>=e_d' for d>=d'>=1.
Choose r with 4^r<=d<4^(r+1), and j=floor((d-4^r)/3), so 0<=d-d'<=2.
Consequently

    Q_d/(2^d sqrt(d)) >= (sqrt(2)/3) sqrt(1-2/d),
    liminf_d Q_d/(2^d sqrt(d)) >= sqrt(2)/3 = .47140452079... . (4)

At d=4^r equality in the spectral upper bound holds, so the limsup is
exactly 1/2. This does NOT establish that the limit is 1/2.

## 3. Multiplicative flat dimensions

More generally, if Boolean f,h satisfy H_a f=sqrt(a)f and H_b h=sqrt(b)h,
then F(x_1,...,x_b)=prod_i f(x_i) h(p_1,...,p_b) satisfies
H_(ab)F=sqrt(ab)F. Indeed H_a[f g(p)]=sqrt(a)f g(1-p), by
anticommutation and the two-dimensional span of 1 and parity.
This proves flatness in all dimensions 4^r directly.

The tempting second-generator strategy is actually IMPOSSIBLE. If a
Boolean f satisfies H_d f=sqrt(d)f, the same anticommutator calculation
gives <f,C_j f>/2^d=1/sqrt(d) for every j. The left side is dyadic
rational. Also sqrt(d) must be an integer, since H_d f is integer-valued.
Consequently sqrt(d) divides a power of two, hence d=4^r. Conversely,
the construction above attains every 4^r. We therefore have the exact
characterization

    Q_d=2^(d-1) sqrt(d) if and only if d is a power of four. (5)

In particular H_9 f=3f is impossible: its anticommutator identity would
require 3 to divide 512. A preliminary 60-second CP-SAT search returned
UNKNOWN, but the arithmetic proof settles nonexistence without search.
Because every signed edge sum for d=9 is even, Q_9<=766, strictly below
the spectral upper 768. This does not rule out asymptotic approach to
1/2 at non-flat dimensions; the arithmetic gap is subleading.

## 4. Weighted diagnostic and literature boundary

Set N=2^d and W_d=sqrt(N/d) H_d. Every row has squared norm N and
max_ij |W_ij|/sqrt(N)=1/sqrt(d)->0. Nevertheless the rigorously proved
cap coefficient is at least sqrt(2)/3 asymptotically. Thus this family
cannot provide the proposed coefficient below .433322; it also does
not prove a general weighted-to-dense comparison.

The exact same fully frustrated hypercube is studied in Marinari--Parisi--
Ritort, https://arxiv.org/pdf/cond-mat/9410089. Their d=9 annealing result
falls short of the spectral bound; numerical failure is not nonexistence.
The 2026 primary manuscript https://arxiv.org/html/2603.29127, Open Problem
7, explicitly leaves d=9 spectral attainment open. Equation (5) above
answers that particular question negatively by a direct independent proof;
we have not checked whether the observation appears in older literature.
Its unrestricted
C4-free graph constructions must not be imported as odd-square/signing
witnesses without checking the extra hypothesis. This proof does not
depend on those literature claims or on their computational certificates.

Reproducibility: decisive_bridge_huang_flat_eigenvectors_2026_09_07.py
checks the exact H4 witness and searched H9; the separate composition
checker verifies the energy identities for arbitrary inner Boolean functions.

## 5. Exact plateau one dimension above every flat dimension

For every r>=0, put s=2^r and d=s^2+1. Then

    Q_d = 2^(d-1) s.                                  (6)

Indeed each coordinate z=(H_d f)_v is an integer congruent to d modulo
two. Thus u=|z| lies outside the open interval (s-1,s+1), and
(u-(s-1))(u-(s+1))>=0. Equivalently

    u <= (u^2+s^2-1)/(2s).

The identity H_d^2=dI gives E u^2=d=s^2+1, whence E u<=s.
For every Boolean f, f^T H_d f<=sum_v |(H_d f)_v|, proving the
upper bound. Ignoring the extra coordinate in the flat s^2-dimensional
witness gives the matching lower bound. This also includes d=2,s=1.

## 6. Precise multiplicative limit criterion, and why (6) is insufficient

Define rho_d=2Q_d/(2^d sqrt(d)), so 0<rho_d<=1. The block construction
actually gives, for ANY positive dimensions a,b,

    rho_(ab) >= rho_a rho_b.                           (7)

Use b identical maximizing a-dimensional inner functions and a maximizing
b-dimensional outer function. The parity-conditioned inner Rayleigh
quotient is constant as in section 1; summing the outer directional
expectations gives exactly the product of the two Rayleigh quotients.
No eigenvector assumption is needed for this energy identity.

A sufficient criterion for rho_d->1 is a sequence q_j=4^(r_j) exp(delta_j)
of integer dimensions with delta_j>0 tending to zero and

    -log(rho_qj)=o(delta_j).                            (8)

Proof: fix j and use dimensions 4^a q_j^b. For every sufficiently large
target d choose 0<=b<=ceil(log(4)/delta_j) such that b delta_j, reduced
modulo log(4), lies immediately below log(d) modulo log(4), with cyclic
gap at most delta_j. Such a choice exists because the initial multiples
0,delta_j,...,floor(log(4)/delta_j)delta_j have gaps at most delta_j.
Choose the integer a>=0 so d'=4^a q_j^b<=d and log(d/d')<=delta_j.
Then ignoring coordinates and (7) yield

    rho_d >= exp(-delta_j/2) rho_qj^ceil(log(4)/delta_j).

First take liminf as d->infinity for fixed j, then j->infinity; (8)
makes the right side tend to one. This is a sufficient theorem, not
evidence that its hypothesis holds.

The exact plateau dimensions q_j=4^j+1 FAIL the needed small-loss test:
delta_j=log(1+4^(-j)) and rho_qj=exp(-delta_j/2) exactly. Their losses
accumulate at a nonzero rate per logarithmic phase shift. Repeated
composition of these particular guaranteed witnesses therefore cannot
by this criterion prove the all-d limit. This is not an upper bound
on rho for the composite dimensions, where better witnesses may exist.

For comparison, at d=s^2+2 (s=2^r>=2), the same lattice-moment argument
gives E|H_d f| <= s+1/(s+1). Its ratio to sqrt(s^2+2) differs from one
by O(s^-3), whereas the logarithmic phase is Theta(s^-2). Thus actual
witnesses within o(1/s) in Rayleigh quotient of the spectral bound at
these dimensions would satisfy (8). The moment upper bound does NOT
provide such witnesses. Exact attainment of that moment bound is itself
usually excluded by its odd denominator; rounding errors in energy are
negligible here, but existence near the bound remains unproved. We stop
at this concrete construction gap rather than infer a limit from an
upper bound.
