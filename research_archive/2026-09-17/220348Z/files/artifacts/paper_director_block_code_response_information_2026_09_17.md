# Sharp information cost of cheap responses on block codes

2026-09-17. **Proved; independent localization/discrepancy reconstructions PASS.** This is a quantitative
consequence of the physical feature-realization theorem, not an improvement
to the original signing bound. Classical ingredients are the scalar
Khintchine inequality, the binomial maximal-atom bound, and entropy duality.
The upper bound additionally requires actual sign realization with exact
covariance, not just a Gaussian target law.

## 1. Operational quantity and theorem

Partition n=br coordinates into r blocks, each of size b. Let C_{b,r}
contain all 2^r Boolean vectors constant on each block. For 0<epsilon<1/4,
let K_{b,r}(epsilon) be the infimum of D(nu||uniform cube) among centered,
exactly isotropic, full-support physical sign laws satisfying

    sup_(x in C_{b,r}) E_nu |h.x|/sqrt(n)<=epsilon.          (1)

If there is no such law, the infimum is infinity. There is a universal
constant C for which, whenever b>=epsilon^(-2),

    K_{b,r}(epsilon)>=r log(1/epsilon)-Cr.                 (2)

Conversely, the stronger finite construction in Appendix A proves for
EVERY r>=1 and b>=4096 epsilon^(-2),

    K_{b,r}(epsilon)<=r log(1/epsilon)+24r.                 (3)

Thus the leading information price is r log(1/epsilon) in a fully finite
regime, uniformly in the number of blocks.
The lower bound is finite-order and does not require isotropy or centering.
The upper bound enforces both exactly, and all physical signs have
positive probability. No solver or unknown optimizing signing is used.

## 2. Finite lower bound from all declared queries

For h define the r normalized block magnetizations

    S_j(h)=b^(-1/2) sum_(i in block j) h_i.

Every block-constant query corresponds to z in {+-1}^r and satisfies
h.x/sqrt(n)=z.S/sqrt(r). Averaging (1) over independent uniform z and
using the sharp lower L1 Khintchine inequality gives

    E_nu ||S||_2 <=sqrt(2r) epsilon,
    E_nu sum_j |S_j| <=sqrt(2)r epsilon.                   (4)

Under the uniform physical cube, S_1,...,S_r are independent copies of
b^(-1/2) times a sum of b uniform signs. Its maximal atom is at most
2/sqrt(b), and the lattice spacing is 2/sqrt(b). Summing the geometric
series on either side of zero therefore gives, for every t>0,

    E_uniform exp(-t|S_j|)<=C0[1/t+1/sqrt(b)],             (5)

with a universal C0, for both parities of b. The estimate is valid even
when its right side exceeds one; that harmless case is not used.

The finite Gibbs variational inequality, applied to
-t sum_j|S_j|, yields

    D(nu||uniform)
      >=-t sqrt(2)r epsilon
         -r log(C0[1/t+1/sqrt(b)]).

Taking t=1/(sqrt(2)epsilon) and b>=epsilon^(-2) proves (2), after
absorbing 1+log[C0(sqrt(2)+1)] into C. No Gaussian approximation, rare
event heuristic, differential entropy, or feature-space density assumption
enters this lower bound.

## 3. First physical upper bound; the finite Appendix A supersedes its range

Let U have columns equal to the normalized indicators of the blocks.
Then U^TU=I_r and every leverage equals1/b=r/n. Every declared query
lies EXACTLY in its feature subspace, so t_x=||U^Tx/sqrt(n)||^2=1.
The hypothesis r^2/b->0 is precisely r^3/n->0, the regime of the
growing-rank variance-realization theorem.

Put kappa=sqrt(2/pi) and fix

    v=epsilon^2/(16 kappa^2),
    cold variance v,       hot variance 1/v,
    cold weight 1/(1+v),   hot weight v/(1+v).

The weighted variance is exactly one. Its target absolute response on
the feature space is

    kappa[ sqrt(v)/(1+v) + (v/(1+v))/sqrt(v)]
       =2 kappa sqrt(v)/(1+v)<=epsilon/2.

The physical construction therefore satisfies (1) for all sufficiently
large orders, with the strict factor-two margin absorbing its uniform
response error and exact-isotropy correction.

Its entropy is bounded by

    (r/2) [(1-v)/(1+v)] log(1/v)+o_epsilon(r)
      <=r log(1/epsilon)+r log(4 kappa)+o_epsilon(r).

This proves the earlier asymptotic version of the upper bound in its
stated range. Appendix A proves the stronger finite assertion (3).
The covariance repair in this earlier version has vanishing mass O_epsilon(r^2/b)
and bounded information cost. Replacing the physical construction by an
informal Gaussian feature mixture would NOT prove (3); exact diagonal
and off-diagonal covariance and the scalar Boolean response are required.

## 4. What the combination explains

This gives a matching operational information law for a nontrivial family
with exponentially many declared queries, including growing rank. It
distinguishes the number of queries (2^r), the entropy of the physical
law relative to independent signs (r log(1/epsilon)), and a law's explicit
description or sampling complexity; these are not identified.

For comparison, a query law with covariance bounded by a fixed multiple
of I requires extensive relative entropy to obtain a fixed normalized
discount (the convex-query concentration theorem). The block code has
query covariance bP, of norm b=n/r, so its r-dimensional price is
consistent with that obstruction. Approximate isotropy alone would not
establish either upper law.

This model is a repetition-code/fixed-aggregate interface, not an actual
minimizer ground code. It supplies an independently useful benchmark and
a quantitative distinction between low-rank and diffuse macroscopic query
geometry. It does not show that original minimizers fall into the cheaper
regime, and no convergence consequence is asserted.

## Appendix A. Exact finite shared-phase construction, with no rank restriction

The localization track supplied and proved this stronger physical upper
bound after auditing Sections1--3. It uses the special block structure
directly, not a growing-rank approximation. For EVERY integer r>=1,
0<epsilon<1, and integer

    b>=4096 epsilon^(-2),

there is a centered, exactly isotropic, full-support physical sign law
with

    sup_(x in C_(b,r)) E|h.x|/sqrt(br)<=7epsilon/8,
    D(nu||uniform cube)<=r log(1/epsilon)+24r.        (A1)

Together with the finite lower bound in Section2, this proves the same
leading information price without an asymptotic aspect-ratio condition
or an o(r) error. The larger absolute constant in (A1) pays fully for
both rare activation and exact physical realization.

### A.1 An elementary binomial local lower bound

Let S=b^(-1/2)sum_(i<=b)epsilon_i. For b>=16 and any point s in its
lattice with |s|<=sqrt(b)/2,

    P(S=s)>=[e^(-1)/(4sqrt(b))] exp(-5s^2).         (A2)

Here is a direct proof avoiding a uniform local-CLT assumption. The
central binomial atom is at least1/(4sqrt(b)). For even b=2m this
follows by induction from
binom(2m,m)4^(-m)>=1/(2sqrt(m)): its successive ratio
(2m-1)/(2m) is at least sqrt((m-1)/m). For odd b the adjacent central
atom differs by the factor (2m+1)/(2m+2)>=1/2, giving the stated safe
bound as well.

By reflection it suffices to move ell places down from k_0=floor(b/2).
Writing b=2k_0+d with d in{0,1}, the probability ratio is

    prod_(j=0)^(ell-1) (k_0-j)/(k_0+d+j+1).

The constraint |s|<=sqrt(b)/2 ensures k_0-ell>=b/4. Each ratio is
1-u_j with u_j<=2/3. Using log(1-u)>=-3u and k_0>=b/3 gives log of
the product at least -9ell(ell+1)/b. Since
|s|=(2ell+d)/sqrt(b), this is at least -5s^2-1. Combining with the
central atom proves (A2), for either parity.

### A.2 Two conditional block laws and their exact moments

Put a=epsilon/4 and t=4/epsilon. Define a cold block law by conditioning
the uniform b-sign cube on |S|<=a, and a hot block law by conditioning
on t<=|S|<=t+1. Write their event probabilities as P_0,P_1 and their
second moments as v_0,v_1. Both laws are invariant under permutations
and global sign reversal within their block. In particular every
individual coordinate has mean zero and

    v_0<=a^2<1,    t^2<=v_1<=(t+1)^2.

The assumed b ensures both bands lie in (A2)'s range and contain
sufficient lattice points. The cold band contains at least a sqrt(b)/2
points, and the positive half of the hot band contains at least
sqrt(b)/4. Consequently

    P_0>=epsilon exp(-21/16)/32,
    P_1>=exp[-1-5(t+1)^2]/16,
    log(1/P_0)<=log(1/epsilon)+5,
    log(1/P_1)<=9t^2.                              (A3)

These are finite estimates, not limits. The last inequality uses t>=4.

Draw ONE common hot/cold label, with hot probability

    pi=(1-v_0)/(v_1-v_0)<=2/t^2.

Conditional on this label, sample the r blocks independently from the
corresponding conditional block law. Denote the resulting physical law
by nu_*. Its coordinates are centered. For two different coordinates
within one block, permutation symmetry gives

    E h_i h_j=(E S^2-1)/(b-1).

Since (1-pi)v_0+pi v_1=1, this correlation vanishes after the common
label is averaged. Cross-block correlations vanish already conditional
on the label, because the blocks are independent and centered. Thus

    E_(nu_*) hh^T=I EXACTLY.                       (A4)

The blocks are not unconditionally independent; their shared phase is
essential. No Gaussian rounding or approximate covariance is involved.

### A.3 Uniform response, information cost and full support

For ANY block-constant query x, conditional block independence and
centering show that (h.x)/sqrt(br) has second moment v_j in phase j.
Cauchy--Schwarz therefore gives

    E_(nu_*) |h.x|/sqrt(br)
       <=(1-pi)sqrt(v_0)+pi sqrt(v_1)
       <=a+2/t=3epsilon/4.                         (A5)

The intermediate bound follows from
pi sqrt(v_1)<=sqrt(v_1)/(v_1-a^2)<=2/sqrt(v_1)<=2/t.
It holds uniformly in the 2^r queries and for EVERY r, without a CLT.

A conditional block law has entropy log(1/P_j) relative to its uniform
block cube. Tensorization within each phase and convexity across the
shared label give

    D(nu_*||uniform)
       <=r[(1-pi)log(1/P_0)+pi log(1/P_1)]
       <=r log(1/epsilon)+23r,                     (A6)

using (A3) and pi<=2/t^2. This explicitly pays for the hot phase's
small probability under the reference cube; it is not treated as a
cost-free variance reservoir.

Finally let delta=epsilon/8 and set

    nu=(1-delta)nu_*+delta Uniform.

This law has full physical support, is centered, and remains exactly
isotropic. The uniform law's normalized response is at most1, so
(A5) gives a response at most7epsilon/8. Entropy convexity gives
D(nu||Uniform)<=(1-delta)D(nu_*||Uniform), proving (A1) with room in
its displayed constant24.

This construction is a genuinely non-Gaussian shared finite-phase law.
It proves a sharp physical information benchmark for block-constant
queries, not a low-response theorem for unknown minimizing nearcodes.
Dimension-free arbitrary-linear subGaussianity is not inferred from its
exact covariance; the following separate proof supplies it.

### A.4 A uniform subGaussian proxy for the same finite law

An elementary general fact is useful here. If a globally even physical
law has density at most K>=1 relative to independent fair signs, then
its linear subGaussian proxy is at most

    8e log(2eK).                                  (A7)

To prove it, normalize a real linear query to ||theta||_2=1. The ordinary
Rademacher MGF/tail bound gives
P_nu(|theta.h|>=s)<=min(1,2K exp(-s^2/2)). Put
A=sqrt(2log(2K)); then |theta.h| is stochastically dominated by A+R,
where R has Rayleigh tail P(R>=s)=exp(-s^2/2). With L=log(2eK)>=1,
Minkowski and E R^(2k)=2^k k! give

    ||theta.h||_(2k)<=sqrt(2L)+sqrt(2k)
                         <=2sqrt(2Lk).

Using k^k<=e^k k!<=e^k(2k-1)!!, its even moments are at most
[8eL]^k(2k-1)!!. Odd moments vanish by global evenness. Summing the
exponential series proves (A7).

Each conditional BLOCK law in A.2 is even and has density at most
K_j=1/P_j relative to its independent b-sign block. By (A3), both
therefore have proxy at most C epsilon^(-2), with a universal C,
independent of b. Within each phase, block independence tensors this
same proxy for arbitrary real queries on all br coordinates. Mixing
the two centered phases retains the maximum proxy, as does the final
uniform mixture. Thus the very same finite exactly isotropic full-
support law in (A1) is C epsilon^(-2)-subGaussian, uniformly in b,r.

There is a sharper anisotropic bound using the independently reconstructed
[fixed-slice swap lemma](paper_discrepancy_kernel_slice_repair_2026_09_17.md),
Section2. That lemma says that, conditional on ANY block magnetization,
the MGF of a block-perpendicular real query theta is at most
exp(||theta||^2). Conditional on either phase in A.2, the normalized
magnetization S is symmetric and bounded by t+1. Hoeffding therefore
bounds its parallel MGF by exp((t+1)^2 u^2/2). Tensorize within the
phase and then average the common phase to obtain, with P the
block-constant orthogonal projector,

    E_nu exp(theta.h)
       <=exp{theta^T[2(I-P)+(4/epsilon+1)^2 P]theta/2}. (A8)

The final uniform mixture preserves this matrix proxy because it
dominates I. Thus a common scalar proxy (4/epsilon+1)^2 is valid for
EVERY r, without a growing-rank or block-size limit. The discrepancy
researcher supplied this refinement; this track read and independently
reconstructed the full slice proof before applying it.

Finite replay:
`computations/paper_localization_2026_09_17_block_code_information.py`
checks40,137 local binomial inequalities,10 cold/hot cases (both block
parities, epsilon down to1/16), and40 rank-scaled covariance/entropy
checks, all PASS. Output is
`tmp/paper_portfolio_2026_09_17/localization/block_code_information.json`.

The discrepancy researcher independently read and reconstructed A.1--A.3
in full, including the parity/lattice estimates and exact isotropy: PASS.
A separate full physical-marginal composition discriminator is proved
in the
[shared-phase comparison](paper_bernoulli_shared_phase_block_composition_2026_09_17.md):
replacing the common phase by independent block phases preserves EVERY
block marginal and global covariance, but restores the leading response
kappa as r grows. Its extensive total-correlation cost is retained
exactly. This track independently audited that proof in full.
