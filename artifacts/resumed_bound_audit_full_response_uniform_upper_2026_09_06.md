# Independent whole-functional ceiling: strictly below 0.45

Date: 2026-09-06. The response agent's structural proposal passes
independent reconstruction. An exact rational adaptive certificate
proves the entire final ONE-response variational functional is strictly
less than 0.45. This is NOT an upper bound on the signing optimum M_n,
nor on arbitrary local algorithms, iterated center updates, or nonlocal
methods. It is a ceiling on this particular lower-bound functional.

## 1. Statement, including arbitrary nonternary responses

Let F be odd, H even, 0<=H, |F|+H<=1 on the Gaussian creation space.
Write P1 for first-chaos projection, U for the creation isometry, and

    K=U*P1F, r=||P1F||_2=||K||_2,
    t0=||F||_2^2-r^2,
    J(F,H)=E[H Psi(K,t0)],
    Psi(k,t)=E|k+sqrt(t)N|.

Set p=E|F| and h=1-p. Since Psi is nonnegative and increasing in t,

    J(F,H)<=E[(1-|F|) Psi(K,p-r^2)].                    (1)

Thus an upper bound does not need an actual purification construction:
both substitutions in (1) are numerical monotonicity statements.
This also avoids any finite-frame or unused-coordinate issue.

For 0<p<1 define

    alpha=Phi^{-1}(1-p/2), m(p)=2 phi(alpha).

The Gaussian moment body gives r<=m(p). Indeed, if r>0, then
G=P1F/r is standard Gaussian and E[FG]=r. Therefore

    r<=E[|F||G|]<=E[|G| 1{|G|>alpha}]=m(p),

where the middle extremum follows by scalar upper-tail rearrangement
under 0<=|F|<=1 and E|F|=p. One has m(p)^2<p. The endpoint cases
p=0 or p=1 have J=0 and can be treated separately.

## 2. Weighted square concavity and monotonicity in the first chaos

For t>0 put s(k,t)=2Phi(k/sqrt(t))-1. On k>=0,

    d/dx Psi(sqrt(x),t)=s(sqrt(x),t)/(2sqrt(x)).

This derivative is nonincreasing because

    s(k,t)=2 integral_0^(k/sqrt(t)) phi(v)dv
           >=2(k/sqrt(t))phi(k/sqrt(t)).

Hence x->Psi(sqrt(x),t) is concave (with its continuous endpoint
extension). Apply Jensen under the probability weight (1-|F|)/h
and use E[(1-|F|)K^2]<=r^2 to obtain

    J(F,H)<=h Psi(r/sqrt(h),p-r^2).                    (2)

For fixed p, put u=r^2 and k=sqrt(u/h). The derivative of the
right side of (2) is

    s(k,p-u)/(2k)-h phi(k/sqrt(p-u))/sqrt(p-u)
      >=p phi(k/sqrt(p-u))/sqrt(p-u)>0.

The k=0 value follows by continuity. Consequently

    J(F,H)<=B(p):=(1-p) Psi(m(p)/sqrt(1-p),p-m(p)^2).  (3)

The floating maximum of this scalar envelope is about .4496994559.
That floating optimization is not used in the uniform certificate.

## 3. Rigorous monotone rectangle certificate

Define the three-variable function

    Fcal(h,p,u)=h Psi(sqrt(u/h),p-u),
    0<h<=1, 0<=u<p.

It is separately increasing in h,p,u. Writing k=sqrt(u/h) and
tau=sqrt(p-u), its derivatives have signs determined by

    d_h Fcal=Psi(k,tau^2)-(k/2)s(k,tau^2)>=0,
    d_p Fcal=h phi(k/tau)/tau>0,
    d_u Fcal=s(k,tau^2)/(2k)-h phi(k/tau)/tau
               >=(1-h)phi(k/tau)/tau>=0.

Parameterize (3) by alpha>=0, with p(alpha)=2Phi(-alpha) and
m(alpha)=2phi(alpha). If alpha lies in a rational interval [a,b],
then monotonicity gives the exact envelope

    B(p(alpha))
      <=Fcal(1-p(b),p(a),m(a)^2).                      (4)

The comparison path first increases h, then p, then u; it stays in the
domain because u<=m(a)^2<p(a). This is a range bound for every alpha,
not an interpolation between endpoint evaluations.

For alpha>=3, Cauchy gives B(p(alpha))<=sqrt(p(alpha))<=sqrt(p(3)).
This tail is less than .052. Near alpha=0, the same Cauchy envelope
avoids evaluating a very large Gaussian argument.

The certificate recursively bisects [0,3] into rational intervals until
the outward interval value in (4), or its Cauchy upper bound, is less
than 9/20. It verifies exact adjacency and coverage of all retained
intervals and includes the entire infinite tail.

## 4. Exact result and replay

Source:
`computations/resumed_bound_audit_full_response_upper_certificate_2026_09_06.py`.
Output:
`computations/results/resumed_bound_audit_full_response_upper_certificate_2026_09_06.json`.

Replay:

    .venv/bin/python computations/resumed_bound_audit_full_response_upper_certificate_2026_09_06.py
      --output computations/results/resumed_bound_audit_full_response_upper_certificate_2026_09_06.json

All arithmetic uses the previously audited 60-digit outward rational
interval class. Pi is enclosed by Machin's formula; density and Gaussian
integrals use finite Taylor series with explicit remainder bounds.
There is no floating-point arithmetic in the certificate.

The run accepted 109 intervals after 217 envelope evaluations and proves

    sup_(F,H) J(F,H)
      <= .449998927491931736970926031607336182570891990067547942316383
      < 9/20.

This removes any possibility of reaching 0.45 by further optimization
of the same final one-response functional, even with infinitely many
creation coordinates and arbitrary feasible measurable F,H. It does
not address whether original M_n/n^(3/2) converges.
