# Adversarial statistical mechanics near the cap: a switching no-go for overlap order

## Scope and verdict

I tested one nonlocal, standard statistical-mechanical object: the
**measure-valued Parisi overlap trajectory**
\(\mathscr P(A)=(\operatorname{Law}_{G_{A,\beta}^{\otimes2}}R_{12})
_{\beta\in[0,\infty]}\). It is a trajectory of probability measures on the
one-dimensional interval \([-1,1]\), hence genuinely sub-landscape. The answer
is negative for arbitrary growing dense interfaces, already on **exact
minimizers**. In fact, the proof gives the stronger collision of the entire
replica-overlap array law:

> the complete overlap-array law is exactly invariant along a switching orbit,
> while a fully free, complete bipartite exact-sign interface distinguishes two
> members of that orbit by at least \(M_n=\Omega(n^{3/2})\).

Thus overlap order by itself is a neighboring-replica / disorder-averaged order
parameter, not a reusable contextual state. At least one gauge-relative,
single-replica “conventional” order parameter must be added. The result below
does **not** say that overlaps are useless for random-disorder thermodynamics,
and it does **not** rule out overlap-based control under an additional low-cap
restriction on the parent interface.

## Definitions and normalization

Let \(\mathcal S_n\) be the hollow symmetric matrices with
\(a_{ij}\in\{-1,1\}\) for \(i\ne j\). In the repository normalization,

\[
 H_A(x)=\sum_{i<j}a_{ij}x_ix_j=\frac12x^{\mathsf T}Ax,
 \qquad
 Q(A)=\max_{x\in\{-1,1\}^n}|H_A(x)|,
\]

\[
 M_n=\min_{A\in\mathcal S_n}Q(A),\qquad
 \mathcal N_n(\epsilon)=
 \{A\in\mathcal S_n:Q(A)\le M_n+\epsilon n^{3/2}\}.
\]

For \(\beta\in[0,\infty)\), use the SK-scale Gibbs measure

\[
 G_{A,\beta}(x)=Z_{A,\beta}^{-1}
 \exp\!\left(\frac{\beta}{\sqrt n}H_A(x)\right).
\]

The factor \(n^{-1/2}\) maps our unnormalized fixed-sign Hamiltonian to the
standard mean-field extensive-energy scale: \(H_A/\sqrt n\) is order \(n\) at
its extremum, and \(Q(A)/n^{3/2}\) is its ground-state energy density. None of
the invariance argument depends on this choice of factor.

For i.i.d. replicas \(X^1,X^2,\ldots\sim G_{A,\beta}\), set

\[
 R_{\ell k}=\frac1n\sum_{i=1}^nX_i^\ell X_i^k,
 \qquad
 \mathscr O_\beta(A)=
 \operatorname{Law}\bigl((R_{\ell k})_{\ell,k\ge1}\bigr).
\]

At \(\beta=\infty\), \(G_{A,\infty}\) means the uniform measure on the
maximizers of \(H_A\). Using the complete array, rather than merely
\(\operatorname{Law}(R_{12})\), makes the falsifier of the candidate
\(\mathscr P\) strictly stronger.

For an \(n\times m\) exact-sign bridge \(B\), let
\(K_m=J_m-I_m\) be the all-\(+1\) hollow child signing, and define the
all-spins-free response

\[
 \mathscr R_B(A)=
 \max_{x\in\{-1,1\}^n,\ y\in\{-1,1\}^m}
 \{H_A(x)+x^{\mathsf T}By+H_{K_m}(y)\}.
\]

Thus the full block matrix
\(\left(\begin{smallmatrix}A&B\\B^{\mathsf T}&K_m\end{smallmatrix}\right)\)
is itself a complete exact signing. Every parent and child spin is optimized.

## The one candidate

### **PROVED STATEMENT — Exact-minimizer switching no-go for Gibbs overlap order**

For every \(n\ge2\), there exist \(A,A'\in\mathcal N_n(0)\) and a complete
\(n\times n\) exact-sign bridge \(B\) such that

\[
 \mathscr O_\beta(A)=\mathscr O_\beta(A')
 \quad\text{for every }\beta\in[0,\infty],                 \tag{1}
\]

but

\[
 |\mathscr R_B(A)-\mathscr R_B(A')|
 \ge M_n
 \ge \frac{n\sqrt{n-1}}{6\sqrt2}.                         \tag{2}
\]

The same equality in (1) holds jointly for every finite collection of
temperatures. Consequently, no function of the complete Gibbs overlap-array
laws—finite-temperature, zero-temperature, or their whole temperature
trajectory—can approximate every growing dense-interface response on exact
minimizers to \(o(n^{3/2})\).

#### Proof

Choose an exact minimizer \(A_0\in\mathcal S_n\). Replacing \(A_0\) by
\(-A_0\), if necessary, gives an exact minimizer \(A\) and a spin \(u\) with

\[
 H_A(u)=Q(A)=M_n.                                        \tag{3}
\]

The uniform average of \(H_A\) on the Boolean cube is zero, so choose \(v\)
with \(H_A(v)\le0\). Let

\[
 D=\operatorname{diag}(u_1v_1,\ldots,u_nv_n),
 \qquad A'=DAD.
\]

Switching preserves hollowness, exact signs, and \(Q\), so
\(A'\in\mathcal N_n(0)\). It also sends \(u\) to \(v\).

For every \(x\), \(H_{A'}(x)=H_A(Dx)\). Hence \(x\mapsto Dx\) pushes
\(G_{A',\beta}\) to \(G_{A,\beta}\). Applying \(D\) simultaneously to all
replicas leaves every overlap fixed:

\[
 \frac1n(DX^\ell)^{\mathsf T}(DX^k)
 =\frac1n(X^\ell)^{\mathsf T}X^k.
\]

This proves (1), including \(\beta=\infty\).

Now take the rank-one, but entrywise dense, exact-sign bridge

\[
 B=u\mathbf1_n^{\mathsf T}\in\{-1,1\}^{n\times n}.
\]

Writing \(a=u^{\mathsf T}x\) and \(t=\mathbf1_n^{\mathsf T}y\), the child
contribution is \(at+(t^2-n)/2\). This is convex in \(t\in[-n,n]\), so its
maximum occurs at \(t=\pm n\). Optimizing the child therefore gives

\[
 \mathscr R_B(C)=\frac{n^2-n}{2}
   +\max_x\{H_C(x)+n|u^{\mathsf T}x|\}.                   \tag{4}
\]

We record the pinning calculation directly. Given \(x\), choose
\(p\in\{u,-u\}\) so that \(d=d_H(x,p)\le n/2\). Only the cut between the
flipped and unflipped coordinates changes in the quadratic energy, whence

\[
 |H_C(x)-H_C(p)|\le2d(n-d).                              \tag{5}
\]

The bridge loss relative to \(p\) is \(2nd\). For \(d>0\), this strictly
exceeds (5); for \(d=0\) there is equality at the pole. Since a quadratic
Hamiltonian is even, (4) is therefore pinned exactly at \(x=\pm u\):

\[
 \mathscr R_B(C)=\frac{3n^2-n}{2}+H_C(u).                 \tag{6}
\]

Apply (6) to \(A\) and \(A'\). Because \(Du=v\),

\[
 \mathscr R_B(A)-\mathscr R_B(A')
 =H_A(u)-H_{A'}(u)
 =M_n-H_A(v)\ge M_n.                                    \tag{7}
\]

It remains only to make the scale in (2) self-contained. For an arbitrary
\(C\in\mathcal S_n\), write

\[
 L(C)=\max_z|z^{\mathsf T}Cz|=2Q(C),\qquad
 K(C)=\max_{x,y}|x^{\mathsf T}Cy|.
\]

If \(X\) is uniform and \(Y_i=\operatorname{sign}((CX)_i)\), then

\[
 \mathbb E[Y^{\mathsf T}CX]
 =n\,\mathbb E|\varepsilon_1+\cdots+\varepsilon_{n-1}|
 \ge n\sqrt{(n-1)/2}.                                   \tag{8}
\]

For fixed \(x,y\), independently choose \(Z_i\) with
\(\mathbb E Z_i=(x_i+y_i)/2\). Hollowness gives

\[
 \mathbb E[Z^{\mathsf T}CZ]
 =\frac14(x+y)^{\mathsf T}C(x+y),
\]

so

\[
 x^{\mathsf T}Cy
 =2\mathbb E[Z^{\mathsf T}CZ]
   -\frac12x^{\mathsf T}Cx-\frac12y^{\mathsf T}Cy.
\]

Thus \(K(C)\le3L(C)=6Q(C)\). Combining with (8) yields
\(Q(C)\ge n\sqrt{n-1}/(6\sqrt2)\), and minimizing over \(C\) proves (2).

\(\square\)

## What is imported, and what is not

### **IMPORTED RESULTS**

1. Panchenko formulates the mean-field order parameter as the overlap array of
   i.i.d. replicas from a random Gibbs measure and proves ultrametricity under
   the Ghirlanda–Guerra identities. This verifies that the one-overlap
   distribution and the stronger array audit used here are standard replica
   objects, rather than surrogates invented for the counterexample:
   [Panchenko, *The Parisi
   ultrametricity conjecture*, Annals of Mathematics 177 (2013)](https://annals.math.princeton.edu/2013/177-1/p08).

2. Auffinger–Chen obtain the zero-temperature Parisi variational formula for
   Gaussian mixed \(p\)-spin ground-state energy. Their normalization has an
   extensive Hamiltonian and ground-state energy per spin; this corresponds to
   \(H_A/\sqrt n\) and \(Q(A)/n^{3/2}\) here. Their proof uses Gaussian
   integration by parts and a zero-temperature limit of the Parisi PDE, neither
   of which holds for a coupling matrix selected adversarially by minimizing
   \(Q\): [Auffinger–Chen, *Parisi formula for the ground state energy in the
   mixed p-spin model* (2016)](https://arxiv.org/abs/1606.05335).

3. Guerra–Toninelli's interpolation proves thermodynamic limits for quenched
   mean-field disorder by smoothly interpolating a large random system with
   two independent random subsystems. The disorder expectation and covariance
   structure supply the sign of the interpolation derivative. They do not
   supply a deterministic inequality for the selected minimizer \(A\), and a
   fixed bridge carries a relative gauge absent from that interpolation:
   [Guerra–Toninelli, *The Thermodynamic Limit in Mean Field Spin Glass Models*
   (2002)](https://arxiv.org/abs/cond-mat/0204280).

4. Huang–Sellke's optimized-spin-glass result concerns algorithms optimizing
   the spins of a random mixed \(p\)-spin Hamiltonian and assumes suitably
   Lipschitz dependence on random disorder coefficients. It does not optimize
   over the disorder coefficients themselves, as \(M_n\) does. This is the
   precise optimized/adversarial-disorder mismatch relevant here:
   [Huang–Sellke, *Tight Lipschitz Hardness for Optimizing Mean Field Spin
   Glasses* (2021)](https://arxiv.org/abs/2110.07847).

5. Chen's deterministic-interaction extension explicitly introduces a
   single-replica “conventional order parameter” in addition to replica
   overlap. This is consistent with, but not used to prove, the switching
   obstruction above: [Chen, *Free energy in spin glass models with
   conventional order* (2024)](https://arxiv.org/abs/2401.10223).

6. The only imported inequality in the proof is the sharp \(p=1\) Khintchine
   lower constant \(1/\sqrt2\) used in (8): [Haagerup, *The best constants in
   the Khintchine inequality*, Studia Mathematica 70
   (1981)](https://www.impan.pl/en/publishing-house/journals-and-series/studia-mathematica/all/70/3/103223/the-best-constants-in-the-khintchine-inequality).

### **PROVED HERE**

All switching identities, exact-minimizer membership, exact dense-interface
pinning in a complete \(2n\)-spin exact-sign parent, response separation, and
the reduction from bilinear to quadratic Boolean norm are proved above. No
random-disorder theorem is used in those steps.

### **CONJECTURES**

None. In particular, I do not conjecture that the overlap array plus one fixed
magnetization suffices: the theorem only proves that an orientation-aware
component is necessary.

## Frontier interpretation

The object is genuinely incomparable with simple contextual shell first
moments. The full overlap array records all multi-replica relative geometry
and possible ultrametric organization, whereas a shell first moment is a
one-replica statistic relative to a declared external pole. Conversely, the
switching pair above has identical full overlap geometry while its pole-relative
response differs at leading scale. A finite family of shell first moments
cannot in general reconstruct the joint law of all replica overlaps, while the
overlap array cannot reconstruct absolute orientation against an interface.

This is an exact-near-minimizer (\(\epsilon=0\)) falsifier for Route R1 if \(P\)
is taken to be a standard Gibbs/Parisi overlap object alone. Its limitation is
equally sharp: the total common pinned baseline is \((3n^2-n)/2\), and the
bridge has rank one as a linear operator, although it has all \(n^2\) nonzero
exact-sign edges, the full \(2n\)-spin parent is an exact signing, and all spins
are free. Therefore the result does not kill a theorem restricted to
parents of cap \(O(n^{3/2})\), to interfaces with operator norm \(O(\sqrt n)\),
or to an overlap object augmented by a gauge-relative conventional parameter.
Those restrictions are exactly where a further statistical-mechanics route
would have to live.

**Benchmark level:** Level 4 physical exact-sign, all-spins-free falsifier.
The child is an exact minimizer, but the \(\Theta(n^2)\) parent baseline keeps
this from certifying a target-budget near-minimizer composition obstruction.
