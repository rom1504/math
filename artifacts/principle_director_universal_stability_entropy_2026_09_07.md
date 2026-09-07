# Uniform entropy loss for almost-stable states of actual low-cap signings

2026-09-07. **Verified, including the total-violation version, by independent
reconstruction.** See `principle_construct_2026_09_07_universal_stability_audit.md`.
This is a structural theorem about EVERY actual bounded-cap full signing,
not merely a random construction or a spectrally flat near-minimizer.

## 1. The statement, including arbitrary external fields

For a symmetric hollow full signing A of order n and any b in R^n define

    I_(A,b)(x)=sum_i [-x_i((Ax)_i+b_i)]_+.

For every fixed C there are delta(C)>0 and n0(C) such that, whenever
Q(A)<=C n^(3/2), n>=n0(C), and b is ANY deterministic external field,

```math
\#\{x\in\{\pm1\}^n:I_{A,b}(x)\le n^{3/2}/1024\}
 \le 2^n\exp[-\delta(C)n].                         (1)
```

In particular this bounds all weak one-spin-stable states, including ties.
It applies separately to A and -A, hence to both possible energy polarities.
The field b can depend on A and n and have arbitrarily large coordinates.
No operator-norm assumption is imposed on A. Exact minimizers and all
additive near-minimizers are included through their bounded normalized cap.

The constant1/1024 is deliberately nonoptimal. All logarithms are natural.

## 2. A uniform-field core lemma

Let n/2<=k<=n, let B be hollow symmetric of order k with off-diagonal
entries +/-1/sqrt(n), and suppose ||B||op<=L, with L>=1. For any b in R^k,
and uniform independent signs X, there is delta_L>0 such that

    P(sum_i[-X_i((BX)_i+b_i)]_+ <= k/512)
       <=exp(-delta_L k)                             (2)

for all sufficiently large k. The constants are uniform over B and b.

Put T=128L and J={i:|b_i|>T}, h=|J|.

### Many large fields

If h>=k/4, a mismatch X_i != sign(b_i) for i in J either has
|(BX)_i|>T/2 or contributes at least T/2 to the total violation. Since
||BX||_2^2<=L^2 k, on the event in (2) the number of mismatches is at most

    4L^2 k/T^2 + 2(k/512)/T
      <= k/4096+k/(32768L) < k/1024.

Consequently its uniform probability is bounded by

    2^(-h) sum_(j<=k/1024) binom(h,j)
       <=exp[-k(log2/4-h2(1/1024))+o(k)].             (3)

Here the sum is enlarged to a Hamming ball in k coordinates. The fixed
coefficient in parentheses is positive. This argument does not assume
independence of the interacting local fields.

### Few large fields

If h<k/4, condition on X_J and put R=[k]\J, k'=|R|>3k/4. On the
remaining independent spins the field is b'=b_R+B_(R,J)X_J, and

    ||b'||_2 <= T sqrt(k)+L sqrt(h) <=129L sqrt(k).

The remaining matrix B_R is hollow, has the same entry magnitudes, and
operator norm at most L. Define the convex function

    F(x)=||B_R x+b'||_1.

Its Euclidean Lipschitz constant is at most L sqrt(k). Every coordinate
of B_R X has a symmetric distribution, so E|Z+b_i'|>=E|Z|. The elementary
Rademacher Khintchine lower bound gives

    E F >= k' sqrt((k'-1)/(2n)) >= k/8               (4)

for large k. The coarse final bound follows from k'>3k/4 and n<=2k.

Apply the convex-zero inequality of
`principle_director_stability_orthant_penalty_2026_09_07.md` to
(F-k/16)_+. Its mean is at least k/16 and its Lipschitz constant is
at most L sqrt(k), whence

    P(F<=k/16)<=exp[-k/(4096L^2)].                   (5)

If the total violation on R is at most k/512, the EXACT identity is

    X^T B_R X+b' dot X
       = F(X)-2sum_(i in R)[-X_i((B_R X)_i+b_i')]_+.

Outside (5), this is greater than k/16-k/256>k/32. Therefore at least
one of X^T B_R X and b' dot X exceeds k/64. The latter has probability
at most exp[-k/(2*64^2*129^2 L^2)] by the elementary product-cosh MGF.
The former has probability at most

    2exp[-c min(k/(64^2),k/(64L))]                 (6)

for a universal positive c, by Hanson--Wright. Indeed the independent
Rademacher coordinates have a uniformly bounded subgaussian norm,
E X^T B_R X=Tr B_R=0, and ||B_R||_F^2=k'(k'-1)/n<=k.

Adding (5)--(6) and the linear tail proves a bound K exp(-c_L k);
enlarging k0 and halving c_L gives (2). It is uniform in the values
conditioned on J, so averaging them completes this case.

Imported tail theorem: Rudelson--Vershynin, *Hanson--Wright inequality and
sub-gaussian concentration*, Theorem1.1,
[primary paper](https://arxiv.org/pdf/1306.2872). Its Frobenius and operator
normalizations, mean subtraction, and independent-coordinate hypotheses
are precisely those checked above. The convex-zero estimate is reconstructed
from Talagrand in the companion artifact, not from a concave-tail shortcut.

## 3. Remove the global operator bound without losing arbitrary fields

Set B=A/sqrt(n). Hollow polarization gives beta(B)<=4Q(A)/sqrt(n)<=4Cn.
The real bilinear Grothendieck relaxation G(B) is at most K_G beta(B), with
K_G<2. Its finite SDP dual gives diagonal D_u,D_v with

    [[D_u,-B],[-B,D_v]] positive semidefinite,
    (Tr D_u+Tr D_v)/2=G(B).

Strong duality follows from a strictly feasible scalar-diagonal dual point.
Since B is symmetric, interchange the two block coordinates and average.
For D=(D_u+D_v)/2 this gives D-B>=0 and D+B>=0, with

    Tr D<=8Cn.

An arbitrarily small dual slack can be absorbed in C. This derivation and
the exact bilinear normalization are independently recorded in
`flatify_construct_2026_09_07_spectral_fourth_moment.md`; the fourth-moment
conclusion is NOT needed for the present theorem.

Write C0=max(C,1) and choose S={i:D_ii<=32C0}. Then |S|>=3n/4 and
||B_S||op<=32C0. Condition on all signs outside S. Stability or total
violation on S is now the core problem with the ARBITRARY external field

    b_S/sqrt(n)+B_(S,S^c)X_(S^c).

If I_(A,b)(X)<=n^(3/2)/1024, its normalized violation on S is at most
n/1024<=|S|/512. The uniform core lemma applies, independently of the
conditioned signs and their resulting field. Averaging yields (1), with
delta(C) a fixed fraction of delta_(32C0). No control of that field's norm
or of the full matrix's largest eigenvalue has been assumed.

## 4. Consequences and limits

Any deterministic OR randomized operation taking a uniform Boolean input
to an output supported on the almost-stable set in (1) has output Shannon
entropy at most n(log2-delta(C)). Thus it loses at least delta(C)n entropy
from the uniform input. This extends the earlier weighted-Haar falsifier
to the actual full-sign bounded-cap class, including actual minimizers.
It does NOT imply the same entropy loss for arbitrary Gibbs inputs, whose
input entropy may already be much smaller.

Moreover choose any fixed rho>0 with h2(rho)<delta(C)/2. A Hamming-ball
union bound shows that a uniform input is farther than rho n from EVERY
almost-stable output in (1), with probability at least1-exp(-delta(C)n/2)
for large n (slightly shrink the strict entropy margin to absorb o(n)).
Thus any such cleanup, regardless of its implementation, changes a linear
number of coordinates on almost all uniform starts. This is stronger than
an entropy-loss assertion and does not assume a particular ascent dynamics.

The result is stronger in this low-cap class than the generic quadratic
strict-local-maxima bound binom(n,floor(n/2)): that bound has only a
polynomial entropy deficit and does not count ties. See Spink,
[Theorem1](https://arxiv.org/pdf/1310.1570). No claim of literature novelty
is made; the new conclusion here is derived with explicit low-cap scope.

This is not a fixed small set answering every external-field query. As b
varies, every spin x can be made strictly stable by choosing b to pin x.
It is a uniform bound on the size for EACH field, not a simultaneous union
bound over all fields. Nor is it an all-order cap construction, a recurrence,
or a proof of convergence. It removes a concrete ambiguity: an entropy-free
cleanup from uniformly distributed spins is impossible even in the actual
near-minimizing class, rather than only in a wider weighted surrogate.

## 5. Sharp universal lower-tail threshold

**Stronger theorem, independently audited by both other roles.** The constant1/1024
in (1) can be replaced by ANY fixed v<1/sqrt(2pi): for every C<infinity
and 0<=v<1/sqrt(2pi) there are delta(C,v)>0 and n0(C,v) such that

```math
\sup_{A:Q(A)\le Cn^{3/2}}\ \sup_{b\in\mathbb R^n}
2^{-n}\#\{x:I_{A,b}(x)\le vn^{3/2}\}
 \le \exp[-\delta(C,v)n] .                         (7)
```

The suprema include all actual full signings of every sufficiently large
order; b need not be bounded. The endpoint is optimal for such a uniform
statement, as proved below. This is a lower-deviation theorem, not two-sided
exponential concentration of the entire field distribution.

Choose fixed epsilon,theta>0 sufficiently small that, with
rho=(1-epsilon)(1-theta),

    sqrt(2/pi) rho^(3/2)>2v.

The diagonal-majorant construction in Section3 gives a principal set S
of size k>=(1-epsilon)n and operator norm at most L=8C0/epsilon, with
arbitrarily small harmless slack if necessary. Condition on the outside
spins and denote the resulting arbitrary normalized field by b.
We work with B=A_S/sqrt(n) and total violation at most vn.

Let J={i:|b_i|>T}. If h=|J|>=theta k, then every sign mismatch on J
either satisfies |(BX)_i|>T/2 or contributes at least T/2 violation.
Their number is at most

    4L^2 k/T^2+2vn/T.

Choose the fixed T large enough that this is at most theta k/4; for
example it suffices that 4L^2/T^2<=theta/8 and
2v/[(1-epsilon)T]<=theta/8. Thus the mismatch fraction among the h
large-field coordinates is at most1/4, and the probability is at most

    exp[-h(log2-h2(1/4))] <=exp[-c theta k].          (8)

If h<theta k, condition on these heavy spins too. The remaining size
k'=k-h satisfies k'>=rho n. Its normalized field has norm at most
(T+L)sqrt(n). With F=||B_R X+b'||_1, symmetry of each row sum gives

    E F >= k' E|S_(k'-1)|/sqrt(n)
        = [sqrt(2/pi)(k'/n)^(3/2)+o(1)]n,            (9)

where S_j is a sum of j independent fair signs. The o(1) is uniform for
k'>=rho n. One elementary exact formula is

    E|S_(2ell)|=2ell binom(2ell,ell)/4^ell,
    E|S_(2ell+1)|=(2ell+1)binom(2ell,ell)/4^ell;

Stirling proves (9), without a joint field Gaussian approximation. Let
a>0 be fixed so that E F>=(2v+3a)n for all large n. The convex-zero
lemma gives

    P(F<=(2v+2a)n)<=exp[-a^2 n/(16L^2)].             (10)

On the event of total violation at most vn, the exact identity from
Section2 then forces X^T B_R X+b' dot X>2an outside (10). Hence either
the quadratic or the linear term exceeds an. Hanson--Wright, with hollow
mean zero and Frobenius squared at most n, and the product linear MGF,
with ||b'||^2<=(T+L)^2 n, give exp(-c_(C,v)n) bounds. All estimates
are uniform in both rounds of conditioning. This proves (7) after absorbing
the fixed prefactors. Parameter order: first C,v, then epsilon,theta,L,T,a,
then n; none depend on a particular signing or external field.

### Why the threshold cannot be larger

For b=0 and uniform X, X_i is independent of (AX)_i because A is hollow.
Moreover (AX)_i has exactly the law of S_(n-1), regardless of the signing.
Consequently, for EVERY A,

    E I_(A,0)(X)=(n/2)E|S_(n-1)|
       =[1/sqrt(2pi)+o(1)]n^(3/2).                 (11)

If v>1/sqrt(2pi), Markov's inequality gives

    P(I_(A,0)<=vn^(3/2))
       >=1-[1/sqrt(2pi)+o(1)]/v>c_v>0.

Thus no exponentially small bound of the form (7) is possible above the
stated endpoint, even when restricted to exact minimizers. At the endpoint
itself no exponential claim is made.

The exact threshold is independent of the unknown optimum. Its role is to
quantify how much local unfavorable response is typical in every low-cap
signing and how atypical a near-stable state must be. An adaptation that
starts from a correlated low-entropy response law may escape the uniform-
input obstruction; (7) must not be silently applied to such a law.

## 6. Quantitative zero-field law without a global operator bound

Let c0=1/sqrt(2pi) and C0=max(C,1). There are universal positive constants
K,c such that for every actual Q(A)<=Cn^(3/2), uniformly in 0<eta<c0/2,

```math
\Pr\{I_{A,0}(X)/n^{3/2}\le c_0-\eta\}
 \le K\exp[-c\eta^3n/C_0^2].                       (12)
```

After increasing K, the inequality holds also when its right side exceeds
one, so small n or very small eta present no problem. In particular

```math
\mathbb E\left|I_{A,0}(X)/n^{3/2}-c_0\right|
 =O(C_0^{2/3}n^{-1/3}),                             (13)
```

uniformly over the actual bounded-cap class. This includes every exact
minimizer; no claim that its minimizing constant has a limit is involved.

Here is the improved quantitative argument. Use D from Section3 and choose
S={D_ii<=8C0/epsilon}, so k>=(1-epsilon)n. The contraction factorization

    B=D^(1/2) T D^(1/2),  ||T||op<=1

follows from D+-B>=0. In particular B_S D_S^(-1) B_S<=D_S. For every
s in [-1,1]^S,

    ||B_S s||^2<=dmax s^T D_S s<=dmax Tr D_S
       <=64C0^2 n/epsilon.                          (14)

Thus F(x)=||B_S x+b||_1 has squared Euclidean Lipschitz constant bounded
by (14), not by the weaker ||B_S||op^2 k. Also, after conditioning outside
spins, b=B_(S,S^c)X_(S^c) satisfies

    ||b||^2<=dmax Tr D_(S^c)<=64C0^2 n/epsilon.       (15)

To see (15), restrict D^(1/2) T D^(1/2) on its output to S and bound its
input D_(S^c)^(1/2)X_(S^c) by its exact squared norm Tr D_(S^c).
This is special to zero ORIGINAL field; an arbitrary added field is not
bounded by (15). The qualitative arbitrary-field result still uses Section5.

Take epsilon=eta/(8c0). The exact signed-sum mean gives

    EF/n>=2c0(1-epsilon)^(3/2)-O(1/n)
         >=2c0-eta/2

once n is sufficiently large compared with1/eta. On the event in (12),
the normalized violation on S is at most(c0-eta)n. Unless F is below
(2c0-3eta/2)n, the quadratic-plus-linear term in the identity of Section2
is at least eta n/2. The convex-zero inequality and (14) bound the former
exception by exp(-c eta^3n/C0^2); (15) gives the same rate for the linear
tail at eta n/4. Hanson--Wright gives the stronger quadratic rate, since
||B_S||_F^2<=n and ||B_S||op<=8C0/epsilon. A fixed union factor gives
(12). If n eta is too small for the signed-sum error estimate, then
n eta^3 is bounded by a universal constant, and a larger K makes (12)
trivial; the constants remain uniform.

Finally integrate (12) over eta to bound the expected negative part below
c0 by O(C0^(2/3)n^(-1/3)). For larger eta use the endpoint eta=c0/2.
The exact mean (11) differs from c0 by O(1/n), so positive and negative
deviations have the same expectation up to that error. This proves (13).

As a further immediate consequence,

    E| ||AX||_1/n^(3/2)-sqrt(2/pi)|=O_C(n^(-1/3)),

because ||AX||_1=2I_(A,0)(X)+X^T AX and
E|X^T AX|<=sqrt(2n(n-1))=O(n). This is an annealed random-spin response
law uniform over deterministic bounded-cap signings, not an extremal
response law. It does not justify replacing the exponentially rare
optimizing spins by independent Gaussian samples.

## 7. Stronger variance law, retained separately

The independently reconstructed convex-Poincare/sign-covariance argument in
`principle_construct_2026_09_07_instability_variance.md` now proves

    Var I_A <=K[n^2+Tr(A^4)/n].

Together with the actual-sign fourth-moment theorem this improves the
normalized L2 error to O_C(n^(-1/2)), and even gives a universal typical
instability law under Q(A)=o(n^2). Thus (13) is no longer the best
mean-square/mean-absolute rate. The exponential lower tail (12), its sharp
threshold, and the arbitrary-external-field theorem remain independent
stronger tail information and are not superseded by Chebyshev's inequality.
