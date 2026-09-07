# An exponential stability penalty from diffuse tilted row variance

2026-09-07. **Verified by independent reconstruction.** This is a
quantitative addition to the ACTUAL weave counting operation. It does not
assume that all pressure-relevant configurations satisfy its hypotheses.

## 1. A direct convex-violation inequality

Let X have independent coordinates in [-1,1], with arbitrary, possibly
different and biased laws. If L is a nonnegative convex Euclidean-Lipschitz
function with Lipschitz constant K, then

```math
\Pr\{L(X)=0\}\le \exp\left[-\frac{(\mathbb E L(X))^2}{16K^2}\right].
```

If K=0, the assertion is interpreted directly. For K>0 let A={L=0} on
the product support. At any z take a subgradient g of norm at most K.
For every y in A, convexity gives
L(z)<=g.(z-y)<=2 sum_j |g_j| 1_{z_j!=y_j}.
The variational formula for Talagrand convex distance therefore gives
d_T(z,A)>=L(z)/(2K). The product-space convex-distance theorem yields
E exp(d_T(X,A)^2/4)<=1/Pr(A). Jensen now proves the displayed inequality.
Empty A is harmless. The argument does not replace a lower-tail theorem
for convex functions by an invalid application to a concave function.

Imported theorem: Talagrand, *Concentration of Measure and Isoperimetric
Inequalities in Product Spaces*, Theorem4.1.1 and Lemma4.1.2,
[primary paper](https://arxiv.org/pdf/math/9406212). Its coordinate induction
works unchanged with different marginal measures. The arbitrary-product
form is also stated in Ledoux, Section3.2,
[primary exposition](https://www.numdam.org/item/SPS_1999__33__120_0.pdf).
The only imported assertion here is the convex-distance exponential moment;
the violation functional and the application below are derived explicitly.

## 2. Orthogonal row theorem, including arbitrarily heavy coordinates

Fix constants p0,C,V,kappa>0 and D>=0. Let H be a k-by-m sign matrix
with HH^T=m I_k and k>=p0 m. Fix signs x_a, real v_j with
sum_j v_j^2<=C m, and d in R^k with ||d||_2<=D sqrt(k).
Let S_j be independent signs of means mu_j. Put J={j:|v_j|<=V} and assume

```math
\sum_{j\in J}v_j^2(1-\mu_j^2)\ge\kappa m.
```

Define the normalized signed physical fields

```math
Z_a=x_a\sum_j\frac{H_{aj}}{\sqrt m}v_jS_j+d_a.
```

Set M=sqrt(2)(sqrt(C/p0)+D), and

```math
g=\sqrt\kappa\,\varphi(M/\sqrt\kappa)
       -M\,\Phi(-M/\sqrt\kappa)>0.
```

For every sufficiently large m, uniformly over all the stated data,

```math
\Pr\{Z_a\ge0\text{ for all }a\}
 \le \exp\left[-\frac{g^2}{256V^2}k\right].                 (1)
```

Proof. Condition on all heavy signs j outside J. Orthogonality bounds the
conditional mean vector by sqrt(Cm)+D sqrt(k), since every coefficient in
its input vector has modulus at most |v_j|. At least k/2 of its coordinates
have absolute mean at most M. Every coordinate has the SAME light variance
s^2=m^-1 sum_J v_j^2(1-mu_j^2), lying in [kappa,C].

Each centered summand in a light coordinate is bounded by2V/sqrt(m), and
the sum of absolute third moments of the summands and their matched
Gaussian replacements is O(VC/sqrt(m)). For completeness one can smooth
the hinge (-z)_+ by (sqrt(z^2+delta^2)-z)/2. The uniform approximation
error is at most delta/2 and its third derivative is O(delta^-2).
The elementary independent Lindeberg telescoping expansion through second
order then has error O(VC/(sqrt(m)delta^2)), uniformly in the coordinate
mean and in the conditioned heavy signs. Fix delta sufficiently small
relative to g, then let m tend to infinity.

For a Gaussian with mean at most M and variance at least kappa, expected
negative part is at least g: it decreases with the mean and increases
with the standard deviation. Thus each of the k/2 selected coordinates
has E(-Z_a)_+>=g/2 for all sufficiently large m.

The convex function L(S_J)=sum_a(-Z_a)_+ therefore has conditional mean
at least gk/4. Its Euclidean Lipschitz constant is at most V sqrt(k):
H/sqrt(m) is a contraction, diag(v_J) has norm at most V, and the sum of
coordinate negative parts is sqrt(k)-Lipschitz. Section1 gives (1), uniformly
in the heavy signs, and averaging them removes the conditioning.

The hypotheses permit concentrated heavy coordinates. What they forbid is
the COMPLETE loss of diffuse thermal variance. No high-dimensional normal
approximation of an exponentially small orthant probability was used.

## 3. Exact mapping to the full-sign weave

Use order-m Hadamards H_i, k retained physical rows, N=mk, independent
fair off-diagonal S_ij=S_ji, and fixed diagonal S_ii in {+-1}. For fixed
full spins x put u_i(j)=sum_a H_i(a,j)x_i(a)/sqrt(k); each row has squared
norm m. The actual full matrix has entries
W_(i,a),(j,b)=S_ij H_i(a,j)H_j(b,i); let C be W with its physical diagonal removed.

Under the exponential moment exp[(t/k)sigma x^T C x], the edge signs are
still independent with means

```math
\mu_{ij}=\tanh(2t\sigma u_i(j)u_j(i)).
```

The exact condition for one-spin stability in fibre i, divided by sqrt(N),
is Z_a>=0 with light/heavy coefficients v_j=u_j(i), j!=i, row sign sigma x_i(a),
and deterministic term

```math
d_a=\sigma S_{ii}x_i(a)H_i(a,i)u_i(i)/\sqrt m
       -\sigma S_{ii}/\sqrt{mk}.
```

Since |u_i(i)|<=sqrt(k), we have ||d||_2<=2sqrt(k), uniformly even on
maximally spiked diagonal coordinates. Deleting column i from H preserves
its contraction bound; zeroing that coefficient lets Section2 apply as
written. A fibre is covered if its incoming energy is at most Cm and
its bounded incoming coordinates retain at least kappa m of tilted variance.

Write A_i for its actual stability event. Each independent edge sign
belongs to just two fibre events, so graph Holder/Finner gives

```math
\Pr_t\{\cap_i A_i\}\le\prod_i\Pr_t(A_i)^{1/2}.
```

Consequently theta m covered fibres give an additional factor
exp[-theta g^2 m k/(512V^2)] in the ACTUAL fixed-spin moment. It is a
strict leading entropy saving at order N, not a formal response relabeling.

## 4. Exact upper-bound consequence and remaining obligation

Suppose the existing row-moment argument bounds the full unconditioned
spin sum by exp[(K+o(1))m^2], where K=t+p log2+E_t(nu_p), k/m->p.
Suppose, in addition, the configurations with fewer than theta m covered
fibres have unconditioned weighted spin sum at most exp[(K-eta+o(1))m^2].
Then counting only stable extrema improves that exponent by

```math
\Delta=\min\{\eta,\theta p g^2/(512V^2)\}>0.
```

The same full-spin/both-polarity Markov argument gives all-order cap bound
(K-Delta)/(2t sqrt(p)), with the existing finite-depth and all-order limits
unchanged once the additional bad-profile estimate is uniform.

The BAD-PROFILE estimate is OPEN. Total incoming energy alone controls
overlarge columns by Markov, but it does not supply diffuse thermal variance:
large paired coordinates can freeze the edge signs. Establishing a penalty
for those profiles requires the actual recursive row-code constraints or
a different estimate. Neither Section2 nor an ordinary Gaussian CLT proves
that missing entropy statement. No convergence or optimal-seed transfer
is asserted here.
