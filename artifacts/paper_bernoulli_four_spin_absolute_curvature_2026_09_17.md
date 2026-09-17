# Full absolute four-spin child: variance concavity fails

2026-09-17, final bounded director proposal, independently reconstructed
by the Bernoulli track. Unlike the two-spin signed-branch example, this
one retains BOTH absolute polarities. It applies to every actual unit-sign
child of order4, including minimizing children.

## 1. Exact full absolute maximum

Fix any hollow full signing D of order4. Take independent Gaussian fields
g1~N(0,v) and g2,g3,g4~N(0,10000). Global reversal preserves H_D, so

    max_y |H_D(y)+g.y|=max_y[|H_D(y)|+g.y].          (1)

For z in{+-1}^3 put w_z=|H_D(1,z)|, which lies in[0,6], and

    a_+=max_z[w_z+g_other.z],
    a_-=max_z[w_z-g_other.z],
    S=(a_+-a_-)/2.

Equation(1) is exactly

    (a_++a_-)/2+|g1+S|.                            (2)

The constant and S depend only on the other three fields, independently
of g1. The intercept difference also gives |S|<=3, since
a_-=max_z[w_-z+g_other.z] and |w_z-w_-z|<=6.

## 2. The high-field event and its exact sign distribution

Let E be the event that every |g_i|>3 for i=2,3,4. On E the unique
maximizers defining a_+ and a_- are z=sign(g_other) and -z respectively:
any changed coordinate loses more than6 in its linear term, exceeding
the whole possible gain from w. Write

    r=sum_(j=2)^4 D_1j z_j,
    h=sum_(2<=i<j<=4) D_ij z_i z_j.

Both r and h are odd sums of three signs, so their absolute values
are1 or3. Hence on E,

    S=(|h+r|-|h-r|)/2=sign(hr)min(|h|,|r|),
    |S| in{1,3}.                                   (3)

The conditional signs z remain uniform and independent because E
depends only on absolute field magnitudes. For every D, exactly6 of
the8 sign patterns have |r|=1. Thus at least3/4 of the good event has
|S|=1, regardless of the internal triangle signs.

## 3. Strict positive curvature, with a rational certificate

For F(v)=E|sqrt(v)G+S|,

    F''(v)=E[phi(S/sqrt(v))(S^2-v)]/(2v^(5/2)).      (4)

Differentiation is legitimate on positive compact intervals; S is
bounded. Put v=1/4. Using phi(0)<2/5, the Gaussian union bound gives

    P(E^c)<3*(6/100)*(2/5)=9/125=.072.

On E intersect{|r|=1}, the numerator integrand in(4) is
phi(2)*3/4>3/80. On the remainder of E it is nonnegative.
Everywhere, including E^c, it is at least-phi(0)/4>-1/10.
Consequently

    E[phi(2S)(S^2-1/4)]
       >(3/4)*(116/125)*(3/80)-(9/125)*(1/10)
       =189/10000.

Since2(1/4)^(5/2)=1/16,

    F''(1/4)>189/625=.3024.                        (5)

The constants phi(0)<2/5 and phi(2)>1/20 have elementary proofs using
25/8<pi<22/7 and e<11/4. The v-independent constant in(2) does not
affect curvature. Thus the expected FULL absolute child maximum is
strictly convex near v=1/4, for every actual D of order4. A sufficiently
small independent symmetric two-point mixing of that variance increases
the expectation.

This is a finite failure of automatic variance concavity under complete
child optimization, not merely a cavity-family diagnostic. The example
uses unequal field variances10000 and1/4; no
unconditional isotropic-bridge or asymptotic minimizing-parent obstruction
is inferred. It does not contradict the positive theorem for independently
averaged uniform-scale-mixture offsets.

Exact replay enumerates ALL64 actual order4 children:48 have absolute
cap4 and16 have cap6. All512 conditional sign-pattern checks and3584
full absolute-maximum identities pass. The lower curvature constant
is replayed as the exact rational189/625. Python compilation passes.
