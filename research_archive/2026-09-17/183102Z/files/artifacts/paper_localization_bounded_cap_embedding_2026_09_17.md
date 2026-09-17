# Linear block sharpness with bounded normalized full-sign cap

2026-09-17. New deduction, self-audited; independent audit pending. This strengthens
critical_full_cap_embedding.md: the deterministic child bias is now
only of order m^(-1/2), so the compared caps are O(N^(3/2)), not N^2.
The constants are deliberately nonoptimal and need not be close to the
minimizing constant.

## 1. Input row law and a spectral-moment lemma

The director's tilted row law in dimension m has exact mean zero,
covariance I_m, an absolute Euclidean L-subGaussian constant, and

    E|sum_j epsilon_j|-E|sum_j Z_j| >= eta sqrt(m)          (1)

for one fixed eta>0 and all sufficiently large m. The ordinary sign
law also has these properties with L=1. Use a common L>=1 for both.

If B consists of n independent rows of either law, then

    E ||B||op^2 <= 32 L^2 (n+m).                         (2)

Indeed, for fixed unit x,y, x^T B y is L-subGaussian. Two 1/4-nets
have combined size at most 9^(n+m), and their bilinear maximum is at
least ||B||op/2. Hence

    P(||B||op^2 > 8 L^2[(n+m)log9+log2+s]) <= exp(-s).

Integrating gives (2). No independence within a row is used.

## 2. Deterministic weakly planted right child

Fix constants K>=1 and lambda>0, set m=Kn (integer K), and let
rho=lambda/sqrt(m), which lies in [0,1] for all large m. There exists
a symmetric hollow full-sign D such that

    D = rho(J-I)+R,       ||R||op <= C sqrt(m), C=8.       (3)

To prove existence, take independent upper-triangle signs of mean rho.
For each unit u, the centered quadratic u^T R u has MGF <=exp(t^2).
A 1/4-net of size 9^m and the symmetric-net inequality give positive
probability of ||R||op<=8sqrt(m). This is a direct bounded-entry
argument, not an imported random-matrix premise.

Write beta=H_D(1). From (3),

    |beta-rho m(m-1)/2| <= (C/2)m^(3/2).                 (4)

Let A be ANY n-vertex full-sign child with Q(A)<=C_A n^(3/2).
It may in particular be an exact minimizer. Set

    W(B) = [ A  B ; B^T D ],     N=n+m.

Every entry of W(B) is a sign. D and A are the SAME for the two row
laws being compared.

## 3. Positive maximum: quantitative approximate pinning

For a right word y, use simultaneous reversal of x,y to arrange that
its negative set S has size k<=m/2; write s=1_S. Then

    H_D(y)-beta
       = -2rho k(m-k)-2s^T R 1+2s^T R s
       <= -lambda sqrt(m)k+2Cm sqrt(k)+2C sqrt(m)k.

For lambda>=4C this is at most

    -(lambda/2)sqrt(m)k+2Cm sqrt(k).

The bridge change from y=1 is at most
2sqrt(n)||B||op sqrt(k). Optimizing the resulting quadratic in sqrt(k)
therefore proves, for P(B)=max_(x,y) H_W(x,y) and r=B1,

    beta+sum_i|r_i|-Q(A) <= P(B)
      <= beta+sum_i|r_i|+Q(A)
         +4C^2 m^(3/2)/lambda
         +4n||B||op^2/(lambda sqrt(m)).                 (5)

The lower bound simply chooses y=1 and x_i=sign(r_i).

## 4. Negative polarity costs no leading uncontrolled term

Let N_-(B)=max_(x,y) -H_W(x,y). Since
-rho[(sum y)^2-m]/2<=rho m/2,

    N_-(B) <= Q(A)+rho m/2+(C/2)m^(3/2)
                    +sqrt(nm)||B||op.

For m>=4 and lambda>=16(C+C_A+1), (4) implies

    beta-2Q(A)-rho m/2-(C/2)m^(3/2)
                                      >= (lambda/16)m^(3/2).

Using P(B)>=beta-Q(A) and (u-t)_+<=u^2/(4t), we obtain

    (N_-(B)-P(B))_+ <= 4n||B||op^2/(lambda sqrt(m)).     (6)

Consequently the ABSOLUTE, unrestricted full-sign cap satisfies

    beta+sum_i|r_i|-Q(A) <= Q(W(B))
      <= beta+sum_i|r_i|+Q(A)
          +4C^2m^(3/2)/lambda
          +8n||B||op^2/(lambda sqrt(m)).                (7)

All bounds before taking expectation hold pointwise in the sign bridge.

## 5. Leading-order separation and its scope

Let W_0 use iid fair sign rows, and W_1 use the tilted rows. Combining
(1), (2), and (7), with m=Kn, gives

    E Q(W_0)-E Q(W_1)
      >= eta n sqrt(m)-2C_A n^(3/2)
         -[4C^2+256L^2(K+1)/K^2]m^(3/2)/lambda.         (8)

The last displayed coefficient follows exactly from
8n/(lambda sqrt(m))*32L^2(n+m).

First choose fixed integer K so that 2C_A/sqrt(K)<=eta/4.
Then choose fixed lambda sufficiently large that it meets Section 4
and the last term in (8) is at most (eta/4)n sqrt(m). Thus

    E Q(W_0)-E Q(W_1)
        >= (eta/2)n sqrt(m)
        = [eta sqrt(K)/(2(1+K)^(3/2))] N^(3/2).         (9)

Both models use n independent random blocks of size m=Theta(N),
with EXACT matching mean and covariance I_m and a common uniform
subGaussian constant. All child entries are common deterministic
full signs. Therefore (9) disproves an o(N^(3/2)) all-offset cap
universality statement at linear block size under those assumptions.
Comparison of both models to their common Gaussian reference would
contradict (9), so no direct Gaussian pinning analysis is required.

In contrast to the earlier hard-ferromagnet embedding,

    E Q(W_a)=O_(K,lambda,L,C_A)(N^(3/2)),  a=0,1.

The same order bound holds with probability 1-exp(-cN), by the net
tail for ||B||op and the deterministic inequalities above. This is a
bounded-normalized-cap example in expectation and with high probability;
the full support can contain exceptional high-cap bridges.

The fixed constants K and lambda are large. The construction is NOT
near-minimizing and does not refute universality confined to optimizing
or quantitatively near-optimizing parent matrices. The left child can
be an actual minimizer, but the right child is deliberately planted.

## 6. Explicit uniform constants and replay

The row law is the three-cell tilt in Section 5 of
[the block theorem](paper_block_universality_2026_09_17.md). No decimal
normal integral is needed. Its limiting middle-cell probability is
greater than 2/45, because it is at least `(2/5)phi(11/10)` and
`sqrt(2pi)<3`, `exp(121/200)<3`. The denominator in delta is at most
`2(21/10)^2=441/50`. Thus `lim delta>20/3969`; the cell inequality
`D<=-117/250` proves a limiting absolute-mean decrease greater than
`26/11025`. One may therefore fix `eta=13/11025` in (1) for all large m.

The limiting probabilities of the first and third cells exceed 1/45
and 2/405 respectively: bound their Gaussian densities below by
`1/(3e)` and `1/(3e^3)`, and use e<3. Eventually p0>1/90 and p2>1/405.
The tilt weights obey `w0<=1+1/(2p0)<=46` and
`w2<=1+1/(2p2)<=407/2<204`, while w1=1/2 and w=1 elsewhere.
The symmetric bounded-density moment argument therefore permits the
explicit common subGaussian choice `L^2=204` for all large m.

An exact minimizing left child has `C_A=1` available elementarily:
an iid full signing has at most 2^n signed Boolean witnesses, and
Hoeffding at threshold n^(3/2) gives union probability at most
`2^n exp(-n^2/(n-1))<1` for n>=2. Thus no archived upper constant is
needed. For example, choose an integer

    K >= (8 C_A/eta)^2

and then choose lambda at least

    max{16(C+C_A+1),
        (4K/eta)[4C^2+256L^2(K+1)/K^2]}, C=8.

All these are FIXED before n tends to infinity. Eventually m>=lambda^2
and m>=4, so rho<=1 and every previous hypothesis holds. Large constants
are harmless for the asymptotic counterexample; no numerical practicality
or small cap constant is claimed.

The replay
`computations/paper_localization_2026_09_17_bounded_cap.py` performs:

- exact symbolic verification of the threshold and completed-square
  inequalities, and the coefficient in (8);
- 240 finite full-cap enumerations through N=16, using
  both a hard-ferromagnetic child and a perturbed one, and checking
  every pointwise inequality in (5)--(7);
- a separate exact pinned-word identity for the hard-ferromagnetic case.

The finite caps are integers and are fully enumerated. Numerical spectral
norms only check the finite inequality hypotheses and are diagnostics,
not proof premises. These small examples test the algebra, not the
asymptotic weak-planting regime. The existence of D in that regime is
the analytic net argument in Section 2.

The successful output is
`tmp/paper_portfolio_2026_09_17/localization/bounded_cap_audit.json`,
seed 2026091703. Two test-setup/representation fixes are retained in
`bounded_cap_validation_notes.md` in the same directory. The independent
exact binomial row-law audit is `critical_block_audit.json`; it checked
normalization and covariance with rational arithmetic at m=100,400,1000,
and 10000. The current proof uses no numerical diagnostic as a premise.
