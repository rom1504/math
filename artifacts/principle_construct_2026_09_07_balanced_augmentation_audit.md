# Independent audit: actual Gram--Schmidt augmentation

2026-09-07. Audit of `principle_invent_2026_09_07_balanced_valley_augmentation.md`, the director's fixed-width simplification, and `principle_invent_2026_09_07_thick_level_complexity_caution.md`.

**PASS, with the substantive scope caution in Section 3.** These are actual full-sign extension arguments, not weighted relaxations. They do not prove convergence.

## 1. Primary input and robust-center theorem

I independently opened the primary [Gram--Schmidt Walk publication](https://theoryofcomputing.org/articles/v015a021/v015a021.pdf). Its Definition on printed page 3 uses the MGF convention exp(sigma^2 ||theta||^2/2); Theorem 1.4 on printed page 5 gives sigma=sqrt(40) for full sign colorings, with arbitrary norm-at-most-one input vectors and initial coloring zero allowed. Thus the imported variance constant is exactly 40, not 40 squared. The guarantee is simultaneous over all deterministic test directions.

For G=LV with column norms of V at most one and row norms of L at most gamma, apply that theorem to w_i=(V_i,e_i)/sqrt(2). For a reference g and Boolean x at distance d, the test vector sqrt(2)(L_g,x-g) represents v dot x and has squared norm at most 2gamma^2+8d. This gives the claimed variance proxy 80gamma^2+320d. Independent full-sign samples v give actual bridge columns; independence across their signs is neither needed nor assumed.

The VC estimate VC(G)<=gamma^2 follows by averaging the margin-one shattering identity over the shattered patterns. Sauer's bound then gives log|G|=o(n) when gamma^2=o(n). Counting a Hamming shell by |G| binom(n,d) is an upper bound even with nonunique nearest references. Union over shells, both cross-energy signs, and all new child spins gives the displayed simultaneous threshold in the audited paper. The confidence term is logarithmic in n, not linear in n.

For q=floor(epsilon n), after n tends to infinity its normalized bridge envelope is at most

    sqrt(640 epsilon r [h(r)+epsilon log 2]),  r=d/n.

The uniform convergence assertion follows from boundedness of h and continuity of square roots on a fixed compact interval. Splitting this into O(sqrt(epsilon) r sqrt(log(e/r))+epsilon sqrt(r)) and paying half the negative coercivity kappa r^alpha against each term gives the two exponents in the paper. Both incremental-cost terms are o(epsilon) for 1<alpha<2. Independent-edge rounding supplies a new child cap at most q^(3/2), adding epsilon^(3/2). Absolute triangle inequalities pay the entire cross channel, so no child-reversal cancellation is hidden.

The liminf contradiction has the correct order of limits: first choose a small fixed epsilon, then let the selected original orders grow. It applies to liminf-realizing signings, not automatically to every exact-minimizer order unless their normalized caps approach the liminf.

## 2. Fixed-width simplification: valid actual construction

Let G_eta contain ALL x with |H_A(x)|>=Q(A)-eta n^(3/2). Use the same augmented feature vectors. Two test directions give simultaneous variance proxies:

    x in G_eta: 80gamma^2, using sqrt(2)(L_x,0);
    arbitrary x: 80n, using sqrt(2)(0,x).

For q independent bridge columns, the following two bounds hold together with positive probability, uniformly in old and new spins:

    sup_(x,y) |x^T C y| <= U,
    sup_(x in G_eta,y) |x^T C y| <= T,

    U=sqrt(160 q n log(8*2^(n+q))),
    T=sqrt(160 q gamma^2 log(8 |G_eta| 2^q)).

Each union bound fails with probability at most 1/4. In particular these bounds concern the SAME actual bridge. For g=gamma^2/n<1/2, Sauer's estimate gives

    U/n^(3/2) <= sqrt(160 epsilon(1+epsilon)log 2)+o(1),
    T/n^(3/2) <= sqrt(160 epsilon g[h(g)+epsilon log 2])+o(1).

Choose epsilon so that the first quantity plus epsilon^(3/2) is below eta. Outside G_eta, the parent energy then stays below Q(A). Inside G_eta the increment is at most T+q^(3/2). If g tends to zero, its normalized limsup is at most epsilon^(3/2). Small fixed epsilon contradicts liminf dilution for a liminf-realizing sequence because c_inf>0. By keeping strict margins, this also gives a uniform positive lower bound on g for signings in a sufficiently small fixed normalized-cap window above c_inf.

No coercivity premise or exponent alpha is needed for this valid implication. However, it has the generic limitation below.

## 3. The whole thick level already has linear complexity

The invent caution is independently valid. Put B=A/sqrt(n), so Q(B)<=Cn. The simultaneous diagonal Grothendieck majorant supplies D>=B,-B with Tr D<=8Cn. Orient a maximizing ground positively. Its local fields are nonnegative with sum at most 2Cn. At least n/2 coordinates simultaneously have field at most 8C and D_ii at most 32C. On this coordinate set I, ||B_I||op<=32C.

For a subset S of I, and z equal to the ground spin restricted to S, the EXACT energy drop after flipping S is

    2 sum_(i in S) h_i - 2 sigma z^T B z <=80C |S|.

Consequently G_eta contains an axis-aligned Boolean face of dimension floor(min(n/2,eta n/(80C))). It therefore has VC dimension, and hence squared unit-column factorization norm, at least this large. This conclusion holds for every bounded-cap matrix, not merely near-minimizers.

Thus the fixed-width simplification is a correct augmentation proof but does not by itself add optimizer-specific structure. The robust-center theorem remains different: it factors selected reference centers and pays excursions through coercivity, rather than factoring the entire thick level that already contains a large Boolean face.
