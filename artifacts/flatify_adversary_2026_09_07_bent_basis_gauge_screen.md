# A prescribed bent basis does not force joint child energy

2026-09-07. Independent actual-child screen. This is not a construction of a low-cap parent.

Let H be an order-n sign Hadamard, and suppose Boolean pairs (x_j,y_j), 1<=j<=n, satisfy x_j^T H y_j=n^(3/2). Let A,D be ANY fixed hollow symmetric sign children, including actual minimizers. There exist row and column sign gauges L,R such that, simultaneously for every j,

    |H_A(L x_j)|, |H_D(R y_j)| <= 15^(1/4) n^(5/4),

for every n for which these pairs exist. Here H_A(x)=x^T A x/2. The gauged bridge LHR still saturates on these gauged pairs.

## Fourth moment and gauge selection

For independent uniform signs epsilon, put Z=sum_{i<j} a_ij epsilon_i epsilon_j. Then

    E Z^4 <= (15/4)[tr(A^2)]^2 <= (15/4)n^4.

One elementary proof compares with independent standard Gaussians g. In the fourth-power expansion, surviving edge multisets are a simple four-cycle or repeated edges. A four-cycle has the same Gaussian and Rademacher moment, including its possibly negative coefficient. Every repeated-edge surviving multiset has positive coefficient and Gaussian moment at least its Rademacher moment. Thus E Z^4 <= E (g^T A g/2)^4. Diagonalizing A and using tr A=0 gives the latter exactly

    (3/4)[tr(A^2)]^2 + 3 tr(A^4)
      <= (15/4)[tr(A^2)]^2.

Each L x_j is uniformly distributed when L is a uniform diagonal sign matrix; independence across j is unnecessary. The union bound over both children and all n pairs, at threshold 15^(1/4)n^(5/4), bounds failure probability by 1/2. This proves existence. Alternatively threshold n^(11/8) gives failure at most (15/2)n^(-1/2), hence success with probability tending to one. Neither optimality nor a cap bound on the children is needed.

More generally, any prescribed family of M=o(n²) Boolean pairs can have all its child energies made o(n^(3/2)) by gauges: the failure bound at threshold epsilon n^(3/2) is at most (15/2)M/(epsilon^4 n²), and epsilon may tend to zero sufficiently slowly. Bridge saturation is needed only for interpreting these pairs as a proposed obstruction, not for this energy estimate.

## Stronger actual-low-cap bound by elementary decoupling

If Q(A),Q(D)<=C n^(3/2), the prescribed family may have size

    log M=o(n^(3/4)),

while retaining the conclusion that all child energies are o(n^(3/2)). In particular a prescribed basis has all child energies O(n sqrt(log n)) with high probability. The following derivation supplies the needed quadratic-chaos bound without an external concentration theorem.

First lambda=||A||op satisfies lambda²<=beta(A)<=4Q(A): each row of A is a cube vector, so every absolute row sum of A² is at most beta(A). For uniform Rademacher epsilon write Z=H_A(epsilon). Choose an independent random vertex bipartition and let B be its rectangular cross matrix. Each edge is cross with probability 1/2, so Z is the partition average of 2u^T Bv. Jensen, followed by the scalar Rademacher bound cosh z<=exp(z²/2), gives

    E exp(tZ) <= E_partition E_v exp(2t²||Bv||²).

Introduce an independent standard Gaussian g to linearize this positive quadratic exponential. Applying the same scalar Rademacher bound to v gives

    E_v exp(2t²||Bv||²)
      <= det(I-4t² BB^T)^(-1/2)
      <= exp(2t²||B||F²/(1-4t²lambda²)).

This holds for 4t²lambda²<1, uniformly in the partition. Since ||B||F²<=n²/4, for |t|<=1/(sqrt(8)lambda) the log moment generating function is at most t²n². Applying Chernoff to both signs with t=min(s/(2n²),1/(sqrt(8)lambda)) yields

    P(|Z|>=s) <= 2 exp[-min{s²/(4n²), s/(4sqrt(2)lambda)}].

Use lambda<=2sqrt(C)n^(3/4) and union over both children and M prescribed pairs. At s=epsilon n^(3/2), the exponent is at least

    min{epsilon² n/4, epsilon n^(3/4)/(8sqrt(2C))}.

When log M=o(n^(3/4)), epsilon can tend to zero slowly enough to dominate log M in both terms. For M=n, choosing s=4n sqrt(log n) gives failure at most 4n^(-3) for sufficiently large n (depending only on C). All this remains a selection of favorable gauges, not an assertion that a fixed arbitrary gauge works.

### The entire elementary permutation-based bent class is covered

Let n=q² with q=2^s and use the normalized Sylvester Walsh matrix on pairs (a,b) in F_2^s x F_2^s. For every permutation pi of F_2^s and every Boolean g on F_2^s, set

    y(a,b)=(-1)^[a dot pi(b)+g(b)].

Its Walsh transform at (u,v) is exactly

    q (-1)^[g(pi^(-1)(u))+v dot pi^(-1)(u)],

by first summing over a. Thus x=Hy/sqrt(n) is Boolean and (x,y) is a saturating bridge pair. There are at most q! 2^q such pairs, whose logarithm is O(sqrt(n) log n)=o(n^(3/4)). Therefore ONE favorable pair of Hadamard gauges neutralizes BOTH child energies on this ENTIRE class, not just on one bent basis. In fact the tail estimate gives a simultaneous bound O_C(n^(5/4) log n): choose a sufficiently large constant multiplying n^(5/4)log n to dominate log(q!2^q) in the linear branch of the tail exponent. Including all invertible affine changes of the 2s input coordinates and affine output additions multiplies the count by only exp(O((log n)²)), so the same conclusion holds for that explicitly enlarged class.

This is an elementary counting-and-selection statement about these explicit bent functions. It is not a bound on the number of all bent functions, and it gives no control of witness families chosen adaptively after the gauges.

## Persistence under low-rank orthogonal modification

Suppose additionally X=[x_1 ... x_n] and Y=[y_1 ... y_n] are orthogonal Boolean bases. Then X=HY/sqrt(n). Replace H by a real T with T T^T=nI and rank(T-H)<=R. Define nonnegative deficits d_j=n^(3/2)-x_j^T T y_j. Exactly

    sum_j d_j = (sqrt(n)/2)||T-H||_F^2 <= 2R n^(3/2).

The equality follows from YX^T=sqrt(n)H^T; the inequality uses ||T-H||op<=2sqrt(n). Therefore all but at most 2R/epsilon of the n basis pairs have bridge value at least (1-epsilon)n^(3/2). The same calculation applies after choosing the gauges above, and T may be chosen after those gauges. A full sign recovery C satisfying beta(C-T)=o(n^(3/2)) preserves these values uniformly.

Consequently, for R=o(n), almost all prescribed basis witnesses can retain a saturated bridge while BOTH child energies remain o(n^(3/2)). Their parent energies are merely n^(3/2)+o(n^(3/2)), not the larger balanced parent budget relevant to the recurrence. This prevents a no-go proof based solely on survival of this prescribed basis. It does not control all exponentially many Boolean Hadamard-saturating pairs, nor the cap of the parent outside the basis.

The permutation-based class above is closed under multiplication of y by every Walsh character: this changes pi by a constant output translation and adds a linear term to g. It therefore splits into free character orbits, each an orthogonal Boolean basis of saturating pairs. Applying the deficit identity in every orbit shows that at least a fraction 1-2R/(epsilon n) of the entire class retains bridge value at least (1-epsilon)n^(3/2). All its child energies can simultaneously remain O_C(n^(5/4)log n) by the favorable gauge selection. Thus retaining many saturation witnesses, even throughout this explicit family, is compatible with energetically neutral actual children on those witnesses.

Walsh orders n=4^s provide such bases: take a bent sign word y, multiply it by all Walsh characters, and use H y_j/sqrt(n) for x_j. Character multiplication preserves bentness, and both resulting column families are orthogonal.
