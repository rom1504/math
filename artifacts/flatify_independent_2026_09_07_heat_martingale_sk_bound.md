# An explicit heat-martingale Parisi lower bound and its dense-bridge consequence

Date: 2026-09-07. Status: proved from the cited Parisi/control theorem;
independently reconstructed by the adversary. This improves the stopped-exit
control without numerical optimization of the Parisi functional.

## 1. The variance-linearized heat martingale

Let W be standard Brownian motion on [0,1], sigma=sign(W_1), and

```math
v_t=E[sigma | F_t]
   =2 Phi(W_t/sqrt(1-t))-1,  0<=t<1.
```

This bounded martingale has terminal value sigma and diffusion coefficient

```math
k_t=sqrt(2/pi)(1-t)^(-1/2)
       exp[-W_t^2/(2(1-t))],
dv_t=k_t dW_t.
```

Direct Gaussian integration gives E k_t=sqrt(2/pi). Two conditionally
independent terminal Gaussians sharing W_t have correlation t, so the
Gaussian sign identity also gives

```math
a(t)=E v_t^2=(2/pi)arcsin(t).
```

Alternatively differentiating by Ito gives
E k_t^2=2/[pi sqrt(1-t^2)], which integrates to the same formula.

Use the deterministic inverse time change

```math
t(s)=sin(pi s/2),  t'(s)=(pi/2)cos(pi s/2),
u_s=v_(t(s)).
```

Define B_s=int_0^(t(s)) [t'(a(r))]^(-1/2)dW_r. Its quadratic variation is
s, and deterministic inversion preserves the Brownian filtration before
the endpoint. The singularity at r=1 is square integrable, since the full
quadratic variation equals 1. Thus B extends to a standard Brownian motion
on [0,1], u is adapted to it, |u|<=1, u_1=sigma, and

```math
E u_s^2=s,
du_s=sqrt(t'(s)) k_(t(s)) dB_s.
```

## 2. A common control valid for every Parisi parameter

For SK covariance Nq^2/2, use
[Auffinger--Chen, Theorem 1 and Corollary 2](https://arxiv.org/html/1606.05335v2).
The zero-temperature variational functional equals a supremum over adapted
controls |u|<=1 of

```math
E|B_1+int_0^1 gamma(s)u_s ds|
 -(1/2)int_0^1 gamma(s)E u_s^2 ds
 -(1/2)int_0^1 s gamma(s)ds.
```

For the particular control above, |z|>=sigma z and the martingale identity
E[sigma u_s]=E u_s^2=s cancel BOTH gamma penalties. This works for every
nonnegative nondecreasing integrable gamma, including its endpoint limits.
Hence the infimum over gamma, namely the SK constant P, satisfies

```math
P>=E[sigma B_1]
  =int_0^1 sqrt(t'(s)) E k_(t(s)) ds
  =int_0^1 sqrt(cos(pi s/2))ds
  =Gamma(3/4)/[sqrt(pi) Gamma(5/4)].                 (1)
```

The beta-integral evaluation in the final line is elementary. Numerically
the constant is .762759763501813188...; the proof uses its exact integral,
not an asserted numerical evaluation of the SK constant.

The key direction is that a SINGLE bounded martingale with linear variance
provides a lower bound on EVERY gamma functional. Choosing a trial gamma
alone would instead provide an upper bound and would not prove (1).

## 3. A simple rational check, independent of special-function numerics

The function f(s)=sqrt(cos(pi s/2)) is concave on [0,1], by differentiation
on [0,1) and continuity. Its 8-panel trapezoidal sum is therefore a lower
bound on its integral. At the seven interior nodes j/8, rational lower
bounds for f are respectively

```math
(990,961,911,840,745,618,441)/1000.
```

Each follows from pi<22/7 and the alternating cosine lower polynomial
sum_(l=0)^5 (-1)^l z^(2l)/(2l)! at z=22j/(7*16), by comparing to the
square of the stated rational. Thus already

```math
P>=int_0^1 f(s)ds > 1/16+(990+961+911+840+745+618+441)/8000
                 =3003/4000 > 3/4.
```

This rational check only verifies a convenient weaker threshold; the exact
integral bound (1) is the theorem.

The exact replay is
`computations/flatify_independent_2026_09_07_heat_martingale_certificate.py`;
its JSON stores every positive rational square margin. Combining the
integral threshold with sqrt(2)<99/70 gives the fully rational bridge
threshold 637/1200>53/100.

## 4. Actual iid balanced bridges

The independently audited Gaussian comparison, Bernoulli Lindeberg transfer,
and pointwise child-cancellation argument in
`flatify_independent_2026_09_07_iid_bridge_control_audit.md` apply unchanged.
For an m by m iid symmetric Bernoulli bridge J, with probability tending
to 1, every full signing having J as its balanced cross block satisfies

```math
Q(parent)/(2m)^(3/2)
 >= Gamma(3/4)/[sqrt(2pi) Gamma(5/4)]-o(1)
 = .539352601188379356...-o(1).                      (2)
```

The diagonal child blocks may be chosen arbitrarily and adaptively after J.
Even the elementary rational threshold gives a floor >3/(4sqrt(2))>.53,
well above the original all-order upper <.493608094.

This is a strong failure theorem for iid dense bridges. It does not rule out
correlated dense bridges and does not settle original convergence.
