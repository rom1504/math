# Conditional endpoint gain from the actual high-degree feedback probe

2026-09-06. Independent endpoint audit. The implication below is proved
CONDITIONALLY on the literal-query stable Gaussian comparison proposed in
Section 4 of `transfer_seed_high_degree_actual_feedback_probe_2026_09_06.md`.
This note does not prove or assume as audited that Boolean closure step.
It keeps the ACTUAL first history and the ACTUAL regressed return.

## 1. Precise quantitative statement and required comparison

Let `B=A/sqrt(n-1)` be a symmetric hollow signing with `||B||op<=L`, where
L is fixed. Necessarily L>=1. Retain the actual fields of the source note:

```math
 F=f(G,Y)\in\{-1,0,1\}^n,\quad H=1-F^2,\quad
 T=BF=V+Z,\quad C=H\operatorname{sign}(T),\quad K=BC,
 u_+=F+C,\qquad u_-=-F+C.
```

The two u's are Boolean, and all these definitions are literal, not
Gaussian substitutes. Define

```math
 e_n(x)=\frac{x^TBx}{2n},\quad
 q_n=\max_{x\in\{-1,1\}^n}|e_n(x)|
     =\frac{Q_{\rm abs}(A)}{n\sqrt{n-1}},\quad
 j_n=\frac1n\mathbb E\,C^TBF
     =\frac1n\sum_i\mathbb E H_i|T_i|.
```

Let E denote one of the FIXED FINITE catalog of probes E_k supplied by
Section 3 of the source. E is exactly centered because its polynomial
is odd. C and K need NOT be centered at finite n: for example the hard
convention `sign(0)=+1` breaks oddness on ties. No centering of those
fields is used below. Set

```math
 v_i=\mathbb E E_i^2,\quad c_i=\mathbb E K_iE_i,\quad
 w_i=\mathbb E K_i^2,\quad t_i^{\rm old}=\mathbb E T_i^2.
```

Thus c_i is still the covariance because E_i is centered, while w_i and
t_i^old are RAW second moments. The regression identity used below is
valid for these raw moments even when K and R have nonzero means.

Assume the established moment/correlation conclusions

```math
 \liminf_n\frac1n\sum_i c_i\ge a_*>0,\qquad
 \frac1n\sum_i v_i\le L^2+o(1),\qquad
 \frac1n\sum_i(w_i+t_i^{\rm old})\le2L^2,             (1)
```

and deterministic ideal variances `0<=vhat_i<=M_0`, with fixed M_0 and
`n^-1 sum_i |v_i-vhat_i|->0`. If k varies with n in a fixed finite set,
take M_0 to be the maximum of the corresponding finite ideal row caps.
The last bound in (1) is exact: `|C_i|,|F_i|<=1` and bounded-op transport
give each of its two summands an average at most L^2.

For positive v_i define the ACTUAL regression

```math
 \beta_i=c_i/v_i,\qquad R_i=K_i-\beta_i E_i.          (2)
```

The CONDITIONAL INPUT is the following joint stable comparison, on fixed
`v_i>=delta` and fixed root second-moment caps, with bounded beta_i:

> E_i can be replaced by sqrt(v_i) N, where N is independent of the
> entire jointly retained literal list `(W_i, Z_i, R_i)`, in averaged
> bounded tests, also under deterministic root restrictions and bounded
> coefficient families with the same fixed caps. Here Z means the ordered
> bounded limit of the fixed-degree
> old-channel lists; its joint law with W and R remains the actual law.

For the endpoint application this must hold also after passing to the
old-measurable bounded functions F,H,sign(T),u_+,u_-. The source's old
small-ball bound and ordered bounded approximation are sufficient for
that passage IF its joint comparison is valid. This is not a claim that
separate comparisons with each member of the list imply the joint one.

Under this conditional input there exists an explicit gamma>0, depending
only on a_*,L,M_0, such that

```math
 \boxed{\displaystyle\liminf_n(q_n-j_n)\ge\gamma>0.}  (3)
```

Thus if the established first-query benchmark gives
`liminf j_n>=j_first`, then `liminf q_n>=j_first+gamma`. No existence
of a limit of j_n is needed. This is a fixed-L assertion; it does not
compare gamma with a spectral-core deletion loss O(1/L).

## 2. Deterministic extraction of a positive set of usable roots

Here are conservative explicit constants. Put

```math
 a=a_*/2,\quad M=2\max(1,M_0),\quad
 \delta=\left(\frac{a}{8L}\right)^2,\quad
 R=\frac{128ML^4}{a^2},\quad
 \rho=\frac{a}{4\sqrt{MR}},\quad
 s_0=\frac{a}{4\sqrt M},\quad K_0=2\sqrt R.           (4)
```

For all large n, the positive average covariance is at least a. The roots
with v_i>M cost only o(1) in average covariance: on that set
`v_i<=2|v_i-vhat_i|`, hence `average_bad(v_i)=o(1)`, and Cauchy--Schwarz
with `average(w_i)<=L^2` gives `average_bad |c_i|=o(1)`.

The roots with v_i<delta cost at most

```math
 \frac1n\sum_{v_i<\delta}|c_i|
 \le\sqrt{\delta\,\frac1n\sum_iw_i}\le L\sqrt\delta=a/8.
```

Among the remaining roots, those with `w_i+t_i^old>R` have fraction at
most `2L^2/R`. Since v_i<=M there, their covariance cost is at most

```math
 \sqrt{\left(\frac{2ML^2}{R}\right)L^2}
 =L^2\sqrt{2M/R}=a/8.
```

Consequently the remaining covariance sum is at least `a n/2` for all
large n. Remove also roots with `c_i<a/4`; their SIGNED covariance sum
is at most `a n/4` (negative terms only improve this removal estimate).
On the final deterministic root set I_n,

```math
 \delta\le v_i\le M,\quad w_i+t_i^{\rm old}\le R,\quad
 c_i\ge a/4,\quad |I_n|/n\ge\rho.                   (5)
```

The last assertion uses `c_i<=sqrt(v_iw_i)<=sqrt(MR)` and the remaining
covariance sum at least `a n/4`. In particular

```math
 s_i:=c_i/\sqrt{v_i}\in[s_0,\sqrt R],\qquad
 |\beta_i|\le\sqrt{R/\delta},\qquad
 \mathbb E R_i^2=w_i-c_i^2/v_i\le w_i.               (6)
```

All cutoffs in (4) are fixed BEFORE the signing limit. We do not need
to remove them afterwards: this one retained positive fraction suffices.

The other literal queries also have adequate root caps on this SAME set.
For completeness the coloring bound (6) of
`continued_feedback_marked_local_noise_separation_2026_09_06.md` gives
`||Y_i||_2<=(27/sqrt(2))L` and
`||(QD)_i||_2<=(27/sqrt(2))L^2`; here rows of B have Euclidean norm one
and rows of Q have norm at most L. Also `||S_i||_2=||G_i||_2=1` and
`||(QS)_i||_2<=L`. Since `b0^2+b1^2<=1`, it follows conservatively that

```math
 \mathbb E V_i^2\le C_V:=L^2+(729/2)L^4,\qquad
 \mathbb E Z_i^2\le2R+2C_V\quad(i\in I_n).
```

Thus W has a uniform L-dependent row second-moment cap, Z has a fixed
cap on I_n, and R has cap R. Original-degree projections and their
partial sums are L2-contractive, so fixed original-degree truncations
of Z and R inherit these caps. No further root deletion depending on
the increasing old-channel cutoff is necessary for these moment bounds.

## 3. Both endpoint orientations have a fixed instability density

The exact endpoint fields are

```math
 Bu_+=T+K=T+R+\beta E,\qquad
 Bu_-=-T+K=-T+R+\beta E.                            (7)
```

The comparison in Section 1 replaces these, jointly with their OLD
endpoint spins, by

```math
 T_i+R_i+s_iN,\qquad -T_i+R_i+s_iN,                 (8)
```

respectively. It does NOT independently Gaussianize T or R. Since u_+
and u_- are old-measurable, the conditional symmetry of N remains valid
after multiplication by either endpoint spin.

We need positive-energy ascent from u_+ and negative-energy ascent from
u_-. Write `(u,sigma)=(u_+,+1)` or `(u_-,-1)`. The negative part of its
oriented local field is

```math
 [-\sigma u_i(Bu)_i]_+.
```

On every retained root, both actual coherent shifts `m_i^+=T_i+R_i`
and `m_i^-=-T_i+R_i` have second moment at most 2R. Thus
`Pr(|m_i^+|<=K_0)>=1/2` and `Pr(|m_i^-|<=K_0)>=1/2`.
For all real |m|<=K_0, every sign z, and every s>=s_0,

```math
 \mathbb E_N[-z(m+sN)]_+
 \ge\mathbb E_N(s_0N-K_0)_+
 =:h_0
 =s_0\phi(K_0/s_0)-K_0\Phi(-K_0/s_0)>0.             (9)
```

This follows from Gaussian symmetry, monotonicity in the worst shift,
and monotonicity of Gaussian positive-part expectation in s. Therefore
the ideal expected negative part on each retained root is at least h_0/2,
for EACH of the two required orientations.

The actual endpoint field has row second moment at most 2R on I_n.
Its ideal counterpart has second moment at most 3R. Truncation at a
fixed large level thus upgrades averaged bounded comparison to this
linear-growth negative-part test, with tail error O(R/level). The old
endpoint discontinuities are handled by the comparison assumption's
ordered old-function passage, not by an independence assertion made
after observing their values.

More explicitly, replace an endpoint spin u by a bounded continuous
old-query approximant u_app. For either actual or ideal field X,
`|[-sigma u X]_+-[-sigma u_app X]_+|<=|u-u_app||X|`.
Averaged Cauchy--Schwarz bounds this replacement by its averaged L2
error times sqrt(3R). The audited old marginal approximation and
small-ball passage provide arbitrarily small such endpoint errors.
This licenses the discontinuous endpoint test without declaring every
bounded measurable test automatically continuous in a joint weak limit.

Consequently the ACTUAL instability densities satisfy

```math
 \liminf_n d_{+,n}\ge\rho h_0/2,\qquad
 \liminf_n d_{-,n}\ge\rho h_0/2,                     (10)
```

where

```math
 d_{+,n}=\frac1n\sum_i\mathbb E[-u_{+,i}(Bu_+)_i]_+,
 \qquad
 d_{-,n}=\frac1n\sum_i\mathbb E[+u_{-,i}(Bu_-)_i]_+.
```

In particular the second density is about POSITIVE physical local
fields at u_-; maximizing minus-energy reverses the instability sign.

## 4. Exact finite spin-ascent certificate and normalization

For any Boolean u and orientation sigma, let
`I_i=1_{sigma u_i(Bu)_i<0}`. Independently conditional on u, flip each
eligible spin with probability p, leaving the others fixed. The resulting
vector x is Boolean. Hollowness gives the exact conditional expectation

```math
 \mathbb E[\sigma e_n(x)-\sigma e_n(u)\mid u]
 =\frac{2p}{n}\sum_i[-\sigma u_i(Bu)_i]_+
   +\frac{2p^2\sigma}{n}(u\circ I)^TB(u\circ I)
 \ge\frac{2p}{n}\sum_i[-\sigma u_i(Bu)_i]_+-2Lp^2.   (11)
```

The interaction term is bounded by the operator norm; no independence
between the eligible spins or the actual field is assumed. The only
independence here is the fresh thinning coins conditional on the actual
endpoint. In particular the formula still holds for very correlated
eligible sets.

Set

```math
 d_0=\min\{L,\rho h_0/4\}>0,\qquad
 p=d_0/(2L),\qquad \gamma=d_0^2/(2L)>0.              (12)
```

Then p<=1/2. By (10), for all sufficiently large n both expected
instability densities are at least d_0. Applying (11) separately yields
actual randomized Boolean x_+,x_- with

```math
 \mathbb E e_n(x_+)\ge\mathbb E e_n(u_+)+\gamma,
 \mathbb E[-e_n(x_-)]\ge\mathbb E[-e_n(u_-)]+\gamma.  (13)
```

Finally put `b_n=E[e_n(F)+e_n(C)]`. Exact symmetry gives

```math
 \mathbb E e_n(u_+)=b_n+j_n,\qquad
 \mathbb E e_n(u_-)=b_n-j_n.
```

Since the deterministic q_n bounds the expectation of every oriented
Boolean energy, (13) implies

```math
 q_n\ge\max\{b_n+j_n+\gamma,-b_n+j_n+\gamma\}
      =j_n+\gamma+|b_n|\ge j_n+\gamma.
```

This proves (3). The common self-energy is retained and canceled only
through the two separately improved orientations. No ordinary block
cancellation or convergence assumption enters the argument.

## 5. Regression and independence audit: exact boundaries

1. Regression alone is NOT independence. Equation (2) proves only
   `E[E_i R_i]=0` and the variance bound in (6). For an elementary
   abstract warning, let E be standard Gaussian and U an independent
   uniform sign, and take K=(1+U)E. Then c=v=1 and R=UE. E is independent
   of U and uncorrelated with R, yet it is not independent of the JOINT
   pair (U,R), which determines E. This is not presented as an original
   signing counterexample; it explains why the full stable module is
   indispensable.

2. The old field BF=T must be retained in the independent query list.
   It contains old Z channels of arbitrarily high original degree.
   Their probe covariance cannot be dismissed merely from the probe's
   high degree. Section 4 of the source explicitly claims the separate
   full-contraction estimate needed here; this note does not audit it.
   Without that input, writing T as an independent coherent shift in
   (8) is invalid, and cancellation in one endpoint can occur.

3. Subtracting beta E from K does not alter the endpoint spins: those
   depend only on the OLD T,F,H. Once the JOINT comparison with (W,Z,R)
   is supplied, independence from these spins follows by old-measurability.
   It is not inferred from the algebraic act of subtraction.
   Means of K and R are retained. A deterministic choice on hard ties
   does not affect this algebra, and need not be falsely declared odd.

4. Both orientations retain the same coefficient beta of the probe in
   (7); changing endpoint or energy orientation only changes a coherent
   shift or a sign. Centered Gaussian symmetry therefore supplies both
   tails. There is no second unexplained innovation assumption.

5. All variance floors, row caps and deterministic regression coefficients
   remain fixed while taking n to infinity. The source's fixed primitive
   degrees P and finite probe-degree catalog are chosen first. Polynomial
   and old-channel cutoffs are then handled in their audited order:
   fixed approximation stage, n first, then approximation error and
   old-channel tail. No increasing probe degree or n-dependent floor
   enters this endpoint implication.

Status: the quantitative endpoint implication PASSES under the explicitly
stated joint stable comparison. Its Boolean-contraction premise remains
the separate audit task; this artifact does not certify that premise.

## 6. Subsequent supply of the separate premise

After this endpoint proof was frozen, the independent reconstruction
`transfer_adversary_high_degree_actual_feedback_audit_2026_09_06.md`,
Sections 2--4, passed the actual mixed covariance and JOINT literal-query
comparison at the required scope. It explicitly restores the whole old
BF through the ordered Z tails, and permits keeping the fixed positive
fraction of capped roots. Its Section 5 also independently checks the
endpoint orientations and normalization here. The director separately
reconstructed the constants and caught the finite hard-tie centering
qualification, which is now incorporated above.

Accordingly the premise of this conditional implication is supplied by
that separate audit for the fixed admissible f, fixed-L first marked
history. The assembled actual fixed-L positive-gain conclusion is
available for integration. The separation of proof obligations is kept
explicit: this note proves the endpoint implication, not the Boolean
contraction estimates. All unrestricted-bound and convergence exclusions
stated above remain in force.
