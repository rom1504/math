# Independent rank-two tile audit

2026-09-07. Algebraic operation PASS; no asymptotic Boolean half-floor or
uniform improving cap bound is asserted.

Let k=2m, N=mk=2m^2, and let F_i be arbitrary real full Hadamard matrices
of order k. Label their columns (j,b), j in [m], b in {0,1}. Write
f_i(a,j)=(F_i[a,(j,0)],F_i[a,(j,1)]) and H=[[1,1],[1,-1]]. Define

    S[(i,a),(j,c)] = (1/2) f_i(a,j)^T H f_j(c,i).

## Full signs, symmetry, and the exact operator

For every two-bit sign vector v, Hv has one zero coordinate and one
coordinate of magnitude two. Consequently every entry of S is a sign,
including diagonal entries. Symmetry follows from H=H^T.

Let U=diag(F_i/sqrt(k)). On microcoordinates (i,j,b), let K transpose i,j
and act by H/sqrt(2) on b. Thus K is symmetric orthogonal. Directly,

    S = sqrt(N) U K U^T,       S^2=N I.

Moreover tr(K)=m tr(H/sqrt(2))=0, hence **tr(S)=0 exactly**. For the
hollow full signing A=S-diag(S), therefore

    H_A(x)=(1/2)x^T Sx,       Q(A)=(1/2)max_x |x^T Sx|.

No O(N) diagonal subtraction is needed for this particular completion.
Changing the diagonal is irrelevant; replacing whole internal fibre
blocks is a different operation with a different error budget.

The spectral upper is Q(A)<=N^(3/2)/2. Equality at finite m is impossible:
a maximizing Boolean eigenvector would satisfy Sx=±sqrt(2)m x, whereas
Sx is integral. This strict finite inequality supplies no asymptotic gap.

## Correct joint tile geometry

Two arbitrary Hadamard columns need NOT have all four two-bit patterns
balanced: an all-ones column provides an immediate counterexample. Their
agreement and disagreement classes each have size k/2, which is enough.
For U_0 consisting of those two columns, p=U_0^T x obeys

    |p_1+p_2|<=k, |p_1-p_2|<=k,
    equivalently |p_1|+|p_2|<=k.

All four vertices ±k e_1,±k e_2 are attained by column spins. The same
holds on the other side of a tile. Hence, without separating four scalar
terms, its bilinear cap is exactly

    max_(x,y) |(1/2)p^T Hq| = k^2/2.

A perfect matching of distinct fibres using pure columns gives normalized
cap 1/(2 sqrt(2)), not one half. The rank-one reciprocal-column floor
therefore cannot simply be reused.

## Exact and sharp joint exponential-moment bounds

Randomize each column of every frame independently by a sign. Conditional
on the unrandomized transformed-spin magnitudes, distinct unordered fibre
pairs then involve disjoint random signs and are independent. Let a,b be
the nonnegative magnitudes at (i,j), and c,d those at (j,i), in normalized
coordinates U^T x. The contribution to Z=x^T Sx/sqrt(N) of this pair is

    Z_ij=sqrt(2)(ac eps_1 eta_1 + ad eps_1 eta_2
                 +bc eps_2 eta_1 - bd eps_2 eta_2).

The exact moment generating function is

    E exp(t Z_ij)
      = product_(w in {ac,ad,bc,bd}) cosh(sqrt(2)t w)
        - product_(w in {ac,ad,bc,bd}) sinh(sqrt(2)t w).

Expand each exponential into cosh times (1+sign*tanh). The only nonempty
even subgraph of K_(2,2) is its four-cycle, with negative coupling product.
This proves the formula, including its sign, without an approximation.

There is also a sharp envelope using only the two pair norms. Put
r^2=a^2+b^2 and s^2=c^2+d^2. The law is symmetric, |Z_ij|<=2rs by the
orthogonality of H/sqrt(2), and E Z_ij^2=2r^2s^2. The elementary inequality

    cosh(tz) <= 1 + [z^2/M^2](cosh(tM)-1),   |z|<=M,

follows termwise from the power series. At M=2rs it gives

    E exp(t Z_ij) <= (1+cosh(2trs))/2 = cosh(trs)^2.

The inequality is sharp: if (a,b)=(r,0) and
(c,d)=(s/sqrt(2),s/sqrt(2)), then Z_ij takes values -2rs,0,2rs with
probabilities 1/4,1/2,1/4. This identifies the reciprocal one-hot/balanced
geometry which a row-profile entropy argument must pay for jointly.
When all four magnitudes equal one, the exact mgf instead simplifies to
cosh(2 sqrt(2)t), strictly smaller than this envelope for t!=0.

An angular refinement avoids losing this distinction. For nonzero r,s set

    U=((a^2-b^2)/(a^2+b^2))^2,
    V=((c^2-d^2)/(c^2+d^2))^2,
    D=1-U-V+2UV = UV+(1-U)(1-V),  0<=D<=1.

Direct fourth-moment expansion gives

    E Z_ij^4/(r^4s^4) = 8-4D.

The negative four-cycle is essential to this formula. If u=t rs, a sharp
two-moment envelope is

    E exp(t Z_ij)
       <= [cosh(u sqrt(2D)) + (1-D)cosh(2u)]/(2-D).       (*)

To prove it, put W=Z_ij^2/(r^2s^2). Then 0<=W<=4, EW=2 and
EW^2=8-4D. The function f(w)=cosh(u sqrt(w)) has f'''(w)>=0, seen from
its power series, including at zero. Its quadratic Hermite interpolant
which is tangent at a=2D and meets f at b=4 majorizes f on [0,4]: the
interpolation remainder has the sign of (w-a)^2(w-b). The distribution
with mass 1/(2-D) at 2D and mass (1-D)/(2-D) at 4 has precisely the same
first two moments, so evaluating that quadratic proves (*).

This envelope is optimal among distributions with these moments. At D=0
it is cosh(u)^2; at D=1 it is cosh(sqrt(2)u), agreeing with the exact tile
law for balanced/balanced and one-hot/one-hot pairs. No claim that every
intermediate extremal moment distribution is realized by an actual tile
is needed. If r*s=0, the mgf is simply one.

A simpler, weaker penalty following from the fourth moment and bounded
support is

    E exp(t Z_ij)
       <= cosh(u)^2 - (D/4)(cosh(2u)-1-2u^2).

These angular penalties alone still do not close an entropy bound: an
unconstrained directed edge-profile relaxation can make one side one-hot
and the other balanced at every pair, giving D=0 everywhere. Actual
Boolean/Hadamard row-profile compatibility is the missing aggregate input.

The diagonal fibre's corresponding pair has mgf

    exp(t(a^2-b^2)/sqrt(2)) cosh(sqrt(2)t ab).

Such fibre terms cannot be discarded uniformly on spiky profiles: they
are not merely the scalar matrix diagonal removed in defining A.

Finally, independent column signs erase a scalar seed sign multiplying
each reciprocal H tile: absorb that seed sign into both column signs at
one endpoint. Thus this independent-pair moment route may yield a new
universal sign construction, but does not by itself retain actual child
information. Correlating column signs can retain it, at the expense of
the independence used above. No aggregate entropy bound is claimed here.

## A finite improving upper valid for every choice of frames

At m=2, N=8, spectral and energy parity imply Q(A)<=10; the known exact
M_8=10 makes every output an actual minimizer.

At m=4, N=32, one has the stronger universal upper Q(A)<=88, giving
88/32^(3/2)=.486135912... . Indeed for any Hadamard S of order 32 and any
spin x, z=Sx has even integral entries, all in the SAME residue class
modulo four. Two rows differ in 16 positions, so their dot products with
x differ by a multiple of four. Also sum z_i^2=32^2=1024.

- In residue class zero, |z_i|<=(z_i^2+32)/12 (tight at magnitudes 4,8),
  so sum |z_i|<=170+2/3.
- In residue class two, |z_i|<=(z_i^2+12)/8 (tight at magnitudes 2,6),
  so sum |z_i|<=176.

Thus |x^T Sx|<=176. This discrete improvement vanishes in the general
large-order estimate and is not a proof of an asymptotic constant below
one half.

## Canonical common Walsh frames: exact scope and tensor lift

If every F_i=H_m tensor H, with m a power of two and the common column
ordering (j,b), then

    S_m = R_m tensor H,
    R_m[(i,a),(j,c)] = H_m[a,j] H_m[c,i].

For m=rL (powers of two), rearranging indices gives

    S_(rL) = R_L tensor S_r.

The sign vector h[(p,q)]=H_L[q,p] satisfies R_L h=Lh, directly by column
orthogonality. Therefore every Boolean witness x for S_r lifts to h tensor
x with quadratic energy multiplied by L^3. Since the order is multiplied
by L^2, its normalized cap lower bound is unchanged. This is an exact
scalable embedding of finite witnesses for the COMMON WALSH family only.
It does not cover independent arbitrary frame choices or column labels.

The canonical exponent is a·j+c·i+de over F_2. Its nondegenerate symmetric
bilinear form is a sum of hyperbolic planes and one anisotropic line. A
hyperbolic pair u,v and unit vector z are replaced by the mutually
orthogonal unit vectors u+z,v+z,u+v+z. Induction makes the form the standard
dot product. Hence S_m is permutation-congruent to the ordinary Sylvester
Hadamard matrix of order 2^(2 log_2(m)+1).

This identifies the canonical case with the odd-variable Boolean Rayleigh
quotient studied by Carlet, Danielsen, Parker, and Solé, *Self-Dual Bent
Functions*. Their Theorem 3 gives a 1/sqrt(2) spectral-fraction lower bound
and the text leaves the odd-variable maximum open; it does not establish
an asymptotic half-floor. Primary source: the paper as reproduced in
[Danielsen's thesis, printed pp.185–187](https://www.codetables.de/larsed/phdthesis.pdf).
This is a statement about what that source proves, not a claim that no
subsequent literature resolves any variant.

Replay: `computations/flatify_adversary_2026_09_07_rank_two_tile_check.py`
and its same-stem result JSON in `computations/results/`.
