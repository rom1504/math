# Large operator normalization makes the opposite-spectral bridge iid-like

2026-09-07. A quantitative scope limit for one actual Gaussian-sign law, using the newly proved quenched universality theorem. This does not exclude other correlated constructions or global old-edge changes.

Let A,D be actual hollow full-sign children of equal order n with Q(A),Q(D)<=C n^(3/2). Put

    p=||A||op ||D||op,
    R=I-rho A tensor D/p,

where 0<=rho<1 is fixed. This is a unit-diagonal correlation matrix with fixed ellipticity and operator bounds. For C=sign(G), Cov(G)=R, its exact covariance is

    Cov(vec C)=I-kappa A tensor D,
    kappa=(2/pi)arcsin(rho/p)=O(1/p).

Let Y be the covariance-matched Gaussian bridge and Z the iid standard Gaussian bridge. Then their expected absolute parent caps, including the SAME deterministic children A,D, differ by at most

    O_C(n^2/sqrt(p)).                                 (1)

Together with quenched sign-to-Gaussian universality, the actual sign bridge has the same expected normalized parent cap as the iid Gaussian bridge whenever p/n tends to infinity.

## Uniform covariance comparison

Index the absolute parent maximum by omega=(sigma,x,y), where sigma is a sign and x,y are Boolean. Its deterministic offset is sigma(H_A(x)+H_D(y)); its random feature is sigma x^T Y y. For another omega'=(sigma',x',y'), the covariance difference from iid is

    -sigma sigma' kappa (x^T A x')(y^T D y').

Since beta(A)<=4Q(A), beta(D)<=4Q(D), its magnitude is uniformly at most

    M=kappa beta(A) beta(D)=O_C(n^3/p).                (2)

This controls cross covariances, not only variances at individual spin pairs.

## Softmax proof of (1)

For a finite indexed Gaussian process with arbitrary deterministic offsets, compare two covariance matrices whose entries differ by at most M. Interpolate the Gaussian laws in the pressure log sum exp(t(offset+process)). Gaussian integration by parts expresses its derivative as one half t^2 times the difference between a diagonal covariance average and a two-replica covariance average. Both are bounded in magnitude by M, so the pressure difference is at most t^2 M.

Divide by t. Each soft maximum differs from its true maximum by at most log|Omega|/t. Thus

    |E max(offset+Y)-E max(offset+Z)|
       <=t M+2log|Omega|/t.

Optimize in t to obtain O(sqrt(M log|Omega|)). Here |Omega|=2^(2n+1), so (2) gives (1). This proof is unaffected by the deterministic child offsets or dependence among the indexed Gaussian variables.

The normalized comparison cost is O_C(sqrt(n/p)). The additional sign-to-Gaussian cost from `flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md` is O_(rho)(n^(-1/6)sqrt(log n)), because epsilon=1-rho and K=1+rho are fixed.

## Consequence and scope

The iid Gaussian parent absolute cap is at least its bridge bilinear cap, for every bridge realization. Indeed whole-child spin reversal leaves both internal energies unchanged and reverses the cross term, giving |H_A(x)+H_D(y)|+|x^T Z y| after maximizing that reversal. The rigorously audited iid-bridge lower theorem in `flatify_adversary_2026_09_07_iid_bridge_floor.md` therefore applies.

Accordingly this opposite-spectral sign ensemble cannot provide a favorable EXPECTED-cap bound when p/n diverges: in that regime its expectation is iid-like at leading order. This does not rule out atypically favorable realizations of the ensemble. It does show that the newly paid quenched expectation transfer, followed by expectation-based existence, needs p=O(n) for this particular covariance to have a nonvanishing leading effect.

The fourth-moment theorem does not force p=O(n), even for general near-minimizers. Its planted-clique examples allow p much larger. Whether selectable exact optimizers have a stronger usable property, or a different covariance avoids this denominator, remains unresolved here.
