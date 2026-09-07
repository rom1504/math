# Bounded regular weights can require dense sign changes to homogenize

2026-09-07. A concrete obstruction to amplitude-only homogenization, not to arbitrary global flatification. The input family below is actual, symmetric, hollow, uniformly bounded in amplitude, exactly row-square-regular, and has cap O(n^(3/2)). It is not asserted to consist of exact or additive near-minimizers.

## 1. Explicit low-cap weighted matrices with a biased sign pattern

Let q=2^d, d>=2, and n=q^2. On the additive group F_q times F_q choose alpha in F_q outside {0,1}. Define

    F(u,v)=(Tr(uv), Tr(alpha uv)) in F_2^2.

Let D be the Cayley adjacency matrix whose connection set consists of z with F(z)=(1,0). Since F(0)=0, D is hollow; characteristic two makes D symmetric. Every row has degree

    d_q=(q^2-q)/4.

To verify this and the spectral estimate without a design-existence import, use

    1_(F(z)=c)=1/4 sum_(s in F_2^2) (-1)^(s.(F(z)-c)).

For every nonzero tau in F_q and additive character indexed by (a,b), direct summation over v gives

    sum_(u,v) (-1)^Tr(tau uv+a u+b v)
       =q (-1)^Tr(ab/tau).

All three nonzero F_2-linear combinations of 1 and alpha are nonzero. At the trivial character the three coefficients for c=(1,0) sum to -1, giving d_q. At every nontrivial character the constant term vanishes and the remaining three Walsh sums each have magnitude q. Hence every nonconstant eigenvalue of D has magnitude at most 3q/4.

Set

    theta_q=d_q/(n-1)=q/[4(q+1)],
    a_q=sqrt(theta_q/(1-theta_q)),
    b_q=sqrt((1-theta_q)/theta_q)=sqrt(3+4/q),
    W=a_q(J-I)-(a_q+b_q)D.

Thus every nonedge of D has weight +a_q, every edge of D has weight -b_q, and W is hollow. Each row has EXACT mean zero and squared norm n-1:

    (1-theta_q)a_q-theta_q b_q=0,
    (1-theta_q)a_q^2+theta_q b_q^2=1.

Its maximum amplitude b_q is at most 2 and tends to sqrt(3). The constant vector is an eigenvector of eigenvalue zero. On its orthogonal complement,

    ||W||op <= a_q+(a_q+b_q)3q/4=O(sqrt(n)).

Consequently Q(W)<=n||W||op/2=O(n^(3/2)). Its sign pattern nevertheless has asymptotic positive fraction 3/4 and negative fraction 1/4. In fact

    Q(sign W) >= H_(sign W)(1)
       =n[(n-1)-2d_q]/2=(1/4+o(1))n^2.

Plain entrywise sign rounding is therefore quadratically expensive on this bounded row-square-regular family.

## 2. Every uniform sign-preserving amplitude damping fails at leading scale

Consider any transformation that retains this sign pattern and maps its two original amplitudes uniformly to +a'_q and -b'_q, then normalizes the rows back to squared norm n-1. This includes arbitrary odd entrywise amplitude functions followed by common row normalization. Suppose its maximum amplitude is at most K'<sqrt(3), where K' is fixed as q grows.

Row normalization requires

    (1-theta_q)(a'_q)^2+theta_q(b'_q)^2=1.

Since b'_q<=K'<b_q eventually, necessarily a'_q>a_q. Its normalized row mean is bounded below by

    (1-theta_q)a'_q-theta_q b'_q
       >=sqrt((1-theta_q)(1-theta_q K'^2))-theta_q K'.

The right side tends to a strictly positive constant for every K'<sqrt(3). Therefore the transformed cap is Omega_(K')(n^2). The same argument rules out any fixed decrease in maximum amplitude from b_q, not only flattening all the way to one.

This conclusion uses exact row-square accounting. Merely scaling W down preserves its zero row mean and lowers cap, but it does not preserve the required row energy and hence is not an amplitude-homogenization solution.

## 3. A stronger polarity-preserving barrier and a dense-edit lower bound

Allow the new amplitudes to vary arbitrarily from edge to edge, with maximum K'. Let V remain hollow, symmetric and row-square-regular. Then

    sum_(i<j) V_ij^2 = n(n-1)/2,
    sum_(i<j) |V_ij| >= n(n-1)/(2K').

If f original positive edges are changed to negative, the total number of negative edges of V is at most n d_q/2+f. Therefore

    H_V(1) >= n(n-1)/(2K')-K' n d_q-2K'f.             (1)

In particular, with NO sign changes and K'<sqrt(2), (1) is a positive constant times n^2. So even fully nonuniform sign-preserving amplitude changes cannot reach such a cap while retaining low cap and row energy.

More generally, if Q(V)<=C' n^(3/2), then (1) forces

    f >= n(n-1)/(4K'^2)-n d_q/2-O_(C',K')(n^(3/2)).   (2)

For K'=1 this becomes

    f >= (1/8-o(1)) n^2.

Thus every low-cap full-sign output must change the POLARITY of a positive fraction of all original edges. The lower bound permits every other edge and amplitude to change arbitrarily; it only uses the resulting low cap, row energy, and amplitude bound.

## Exact scope

This proves that a natural sign-preserving amplitude damping cannot serve as a universal row-regular flatifier, and quantifies the dense sign changes required on an explicit actual input family. It does not rule out a genuinely global operation that makes those changes with correlated cap control. The constructed inputs have merely bounded normalized cap, not a proved optimum-level cap, so no claim about selectable exact minimizers is made.

## 4. Implantation at power-saving distance from actual optimal cap

The obstruction to a COMMON entrywise amplitude map can be implanted in weighted inputs whose cap is within a power-saving additive error of M_N. Let N=2^(5d), s=2^(4d)=N^(4/5), and choose an actual optimal full signing A_N. The archived Grothendieck-Pietsch regular-core theorem gives a principal core of at least N/2 vertices and operator norm O(sqrt(N)), with an absolute constant because Q(A_N)=O(N^(3/2)). Choose S of size s in that core and replace A_N[S] by the explicit weighted matrix W_s above; leave every other edge unchanged. Denote the resulting global weighted matrix by V_N.

Every row still has squared norm N-1: the replaced block has the same local row energy s-1 as the old sign block. Amplitudes remain bounded by 2. Principal spectral monotonicity on the core gives Q(A_N[S])=O(s sqrt(N)), while Q(W_s)=O(s^(3/2)). Therefore

    |Q(V_N)-M_N| <= Q(A_N[S])+Q(W_s)=O(N^(13/10)).

Now apply a common sign-preserving scalar amplitude map, with any common normalization, and require exact global row-square regularity. Every vertex outside S sees only original unit amplitudes, so the final unit-amplitude image must be exactly one. The outside and cross edges consequently stay unchanged, and the mapped S block must again have local row energy s-1. Section 2 applies inside S. If the new maximum amplitude is below a fixed K'<sqrt(3), it creates a positive S-block energy of order s^2=N^(8/5). The remaining edges have cap O(N^(3/2)), because their matrix is A_N with A_N[S] removed. Hence the output cap is Omega(N^(8/5)).

This is a failure on bounded-amplitude weighted POWER near-minimizers relative to the original actual value, not merely on the high-cap Cayley family by itself. The operator class is important: a vertex-dependent or globally correlated transformation is not a common scalar amplitude map.

Any low-cap actual-sign output must also change Omega(s^2)=Omega(N^(8/5)) edge polarities. Indeed the old A_N[S] has negative-edge count s(s-1)/4+O(s sqrt(N)), while the replacement sign pattern has negative-edge count (1/8+o(1))s^2. Thus sign(V_N) has a negative-edge deficit (1/8+o(1))s^2 relative to the globally almost-balanced sign count forced by Q(A_N)=O(N^(3/2)). Every low-cap output has globally half its edges negative up to O(N^(3/2)), which is negligible compared with s^2.

## 5. Exact global variance-compensation bypass, and why it is not progress

The preceding example is NOT a counterexample to global homogenization. For completeness, an explicit intermediate amplitude-reduction operation exists. Suppose all nonunit amplitudes are inside S, its local row energy is s-1, and its largest amplitude is b>1. Set alpha=1/b, multiply the S block by alpha, all cross edges by gamma, and the outside block by beta, where

    gamma^2=1+(1-alpha^2)(s-1)/(N-s),
    beta^2=1-(1-alpha^2)s(s-1)/[(N-s)(N-s-1)].

For s=o(N), beta is real. Direct row counting verifies exact row energy N-1 in both vertex classes. The new maximum amplitude is gamma=1+O(s/N), while the S-block maximum is one. The cap change is bounded by

    (1-alpha)Q(W[S])+O(s/N)Q(W).

Here the cross-cut bilinear norm and outside principal cap are each bounded by Q(W); whole-shore reversal proves the former. Scaling globally by 1/gamma and independent sign rounding has variance O(Ns), yielding an additional O(N sqrt(s)) cap error.

However this operation is DOMINATED for the implanted example by an elementary block replacement: keep every edge outside S and put any good full-sign matrix B_s inside S. Then

    Q(output)<=Q(W)+Q(W[S])+Q(B_s).

When Q(W[S])=O(s^(3/2)), this costs only O(s^(3/2)), better than the compensation construction. Thus no cross-order progress is attributed to variance compensation or to the localized example. The actual macroscopic two-block amplitude problem remains open.
