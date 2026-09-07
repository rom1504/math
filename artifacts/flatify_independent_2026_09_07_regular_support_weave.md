# A positive correlated bridge: direct-E weaving on regular supports

2026-09-07. Extension of the already audited row theorem. The adversarial
researcher independently reconstructed the defect, graph contraction,
dimension uniformity, and all-order bipartite normalization.

Let

    U=[97/20+(24/25)log2-5151/6250]/[(97/10)sqrt(24/25)]
      =0.49360809358874863... .

The existing direct-E argument extends to loopless regular fibre-support
graphs, rather than only a complete graph of fibres. In particular, defining

    R_n=min_{B in {+1,-1}^{n by n}} max_{x,y in {+1,-1}^n}|x^T B y|,

the extension gives the all-order bound

    limsup_n R_n/n^(3/2) <= 2U < 0.987216188.

This is a genuine dense correlated bridge with no iid assumption. It does
not by itself transfer an arbitrary child's cap or prove convergence of
the original M_n. Triangle filling of two child blocks remains insufficient.

## 1. Actual construction on a d-regular support

Let G be any loopless simple d-regular graph on m fibre vertices. At each
vertex i take an independent d by d Hadamard basis H_i from the finite-depth
recursive ensemble in the direct-E proof, with a fresh uniform permutation
of its output columns. Label its d columns by the neighbors of i. Retain
k rows in each fibre, k/d -> p. For each support edge {i,j}, choose an
independent fair sign S_ij, and put

    W[(i,a),(j,b)]=S_ij H_i[a,j] H_j[b,i].

All these entries are actual signs. Nonedges, including every within-fibre
block, are zero. The resulting order is N=mk and every row has degree dk.
No unsupported coordinate is silently filled at this stage.

For x=(x_i), put h_i=H_i[T_i,:]^T x_i. Orthogonality gives

    sum_i ||h_i||^2=m d k.

For sigma=+/-1 define the directed support-edge defect

    D_sigma=sum_i sum_{j~i}(h_i(j)-sigma S_ij h_j(i))^2
           =2(m d k-sigma x^T W x).

Every undirected edge occurs twice. Thus averaging its sign in
exp[-t D_sigma/(2k)] produces exactly the same folded Gaussian kernel

    K_t(a,b)=[exp(-t(a-b)^2)+exp(-t(a+b)^2)]/2,

at a=h_i(j)/sqrt(k), b=h_j(i)/sqrt(k). There is no diagonal defect here and
no lost factor of two.

## 2. Why the row bound does not require the support to be complete

The kernel K_t is a positive-semidefinite Gram kernel: use the even part of
the Gaussian Fock feature map. For a fixed row vector v_i in R^d, average
the tensor product of its d feature vectors over output permutations. Call
the resulting symmetric tensor T_i. Its squared Hilbert norm is precisely

    ||T_i||^2=P_t(v_i)
      =E_{signed permutation g} exp(-t||v_i-gv_i||^2).

After averaging the independent output permutations, the support-edge
kernel product is the contraction of the tensors T_i along the edges of G.
Every tensor index occurs in exactly two factors. Repeated Cauchy--Schwarz,
equivalently the graph form of generalized Holder with exponent two at
every vertex, bounds this contraction by

    product_i ||T_i|| = product_i L_t(v_i).

One elementary version successively merges adjacent tensors, contracting
ALL their common indices at once: after reshaping, each step is matrix
multiplication, with ||AB||_F<=||A||_F||B||_F. Contracting all common indices
prevents unaccounted self-loop trace factors. This argument applies to any
finite loopless support graph. Finite feature
truncation proves it in finite dimension; convergence of the Gaussian feature
expansion passes to the full kernel. Completeness of G is not used. Also,
unlike the complete-fibre proof with self coordinates, no deleted-coordinate
factor is needed.

The previously proved recursive row theorem applies in dimension d:

    limsup_d d^(-1) log E sum_{x_i in {+1,-1}^k}
        L_t(H_i[T_i,:]^T x_i/sqrt(k))
        <=p log2+B^r Phi_t(nu_p)

at every fixed depth r. The row estimate is uniform in the terminal
Hadamards and independent of the number m of fibres. Independence across
fibres therefore yields, for every fixed strict exponent margin,

    P(exists x,sigma: D_sigma<=2 gamma m d k)
       <=2 exp[m d {t gamma+p log2+B^r Phi_t(nu_p)+o_d(1)}].

Choose a sufficiently large but fixed r using B^r Phi -> E, then let d grow.
The error is o(md), even if m also grows. Choosing gamma just below
-(p log2+E_t(nu_p))/t gives an actual signing with

    Q(W) <= [C(p,t)+o(1)] sqrt(d/m) N^(3/2),
    C(p,t)=[t+p log2+E_t(nu_p)]/(2t sqrt(p)).

Equivalently Q(W)<=[C(p,t)+o(1)] N sqrt(dk), the natural row-degree scale.
The certified rational point makes C(p,t)<=U.

## 3. Bipartite specialization and all orders

Take G=K_{d,d}, so m=2d. The actual matrix is

    W=[[0,B],[B^T,0]],

where each side has n=d k vertices and B is a FULL square sign matrix.
Its original quadratic cap is exactly the bilinear norm of B. Hence

    ||B||_{infinity to one}
       <=[U+o(1)](1/sqrt2)(2n)^(3/2)
       =[2U+o(1)]n^(3/2).

At fixed recursive depth, the allowed basis dimensions d=2^r 2^a 12^b have
successive multiplicative gaps tending to one, by the already checked
H2/H12 order-filling argument. Then n=d floor(pd) has the same property.
For an arbitrary target order n0 take the next such n>=n0, with n/n0->1,
and restrict B to any n0 rows and n0 columns. Bilinear cap cannot increase:
extend the retained spins by independent unbiased omitted spins and average.
This proves the displayed all-order limsup, without a growing-depth theorem.

Before row retention, the bipartite block B is itself a Hadamard matrix of
order d^2. Indeed, orthogonality in each right fibre first kills products
between distinct left fibres, and orthogonality in the remaining left basis
then gives BB^T=d^2 I. The retained bridge is an n by n submatrix of this
larger Hadamard, with n approximately p d^2; it need not be Hadamard at its
retained order. Thus this positive operation genuinely leaves the finite
Hadamard-bridge class excluded in the order-eight child theorem.

There is also an exact check on the scale. Match the d left fibres bijectively
with the d right fibres. In each matched pair take its spins to be the
retained partner columns, choosing relative signs to satisfy S_ij. Give each
pair an independent common random reversal. Matched edges contribute k^2;
all unmatched fibre interactions have mean zero. Therefore EVERY such
retained bipartite weave has

    ||B||_{infinity to one} >= d k^2
                              =sqrt(k/d) n^(3/2).

At the certified p=24/25 this class floor tends to sqrt(24/25)=.9797958...,
below the proved upper2U=.98721618... . At p=1 the matching lower bound and
Hadamard spectral upper bound agree exactly at n^(3/2). This explains why
nontrivial retention, not merely a different full Hadamard bridge, is the
positive mechanism.

## 4. Balanced multipartite holes

For G complete L-partite with equal fibre classes of size h, m=Lh and
d=(L-1)h. The same operation gives a full sign matrix between L equal macro
blocks, zero within those blocks, and

    Q(cross matrix)<=[U+o(1)]sqrt(1-1/L) N^(3/2).

This claim is along dimensions for which the required d-dimensional bases
exist; L=2 has the explicit all-order realization above. General fixed L
needs the corresponding harmless arithmetic compatibility checked before
asserting all-order realization.

These cross signs are highly correlated through their fibre bases. The
independent balanced-bridge SK floor does not apply to them. Conversely,
their low cross cap alone does not pay the child's intermediate energy
shells: the elementary bound after inserting children is still
Q(parent)<=Q(cross)+sum Q(child), which is too weak for convergence.

## Dependencies and remaining audit

The only analytic inputs beyond the displayed support-graph calculation are
the fixed-depth row theorem, its scalar closure B^r Phi -> E, and the
rational E certificate already reconstructed in
`decisive_audit_standalone_direct_E_upper_2026_09_07.md` and
`flatify_independent_2026_09_07_ternary_upper_proof.md`.
No new source-channel optimization or unproved phase equality is used.

## 5. Brief rectangular corollary (including an explicit 2:1 primitive)

On K_{d,d}, retain k_L rows per left fibre and k_R rows per right fibre,
with k_s/d -> p_s in (0,1]. The two side orders are n_s=d k_s. Using the
normalized edge defect

    sum_{i,j}(h_L,i(j)/sqrt(k_L)-sigma h_R,j(i)/sqrt(k_R))^2
       =2d^2-2sigma H_B/sqrt(k_L k_R)

and applying the same row theorem separately to both sides gives

    ||B|| <= [1+{(p_L+p_R)log2+E_t(nu_pL)+E_t(nu_pR)}/(2t)+o(1)]
              d^3 sqrt(p_L p_R).

This is a full rectangular signing. For p_L>=p_R, division by the usual
rectangular scale n_L sqrt(n_R) gives coefficient

    [1+{(p_L+p_R)log2+E_t(nu_pL)+E_t(nu_pR)}/(2t)]/sqrt(p_L).

For a cheap explicit 2:1 example set p_L=24/25,p_R=12/25,t=97/20. Besides
the existing E_left<=-5151/6250 certificate, use E_t(nu)<=Phi_t(nu)=-F_t(nu)/2.
In any self-coupling (X,Y) of the ternary source, e=P(X!=Y) has distortion
at least e/p and conditional entropy at most h(e)+e log2. Hence

    F_t(nu_p)>=H(nu_p)-log(1+2 exp(-t/p)),
    E_t(nu_p)<=-[H(nu_p)-log(1+2 exp(-t/p))]/2.

This follows by maximizing h(e)+e log2-(t/p)e over e in [0,1]; no source
optimizer is assumed. Directed50-digit arithmetic then certifies a bridge
coefficient less than99/100 (the displayed nonessential approximation is
.9850025401). The script is
`computations/flatify_independent_2026_09_07_rectangular_bridge_certificate.py`
with its same-stem result JSON.

Choose k_R=floor((12/25)d), k_L=2k_R. The side sizes are exactly2n and n,
where n=d k_R. Dense allowed d orders, followed by restriction to2n0 rows
and n0 columns, give this bridge bound for EVERY sufficiently large n0.
Again this controls the bridge alone, not its alignment with child shells.
