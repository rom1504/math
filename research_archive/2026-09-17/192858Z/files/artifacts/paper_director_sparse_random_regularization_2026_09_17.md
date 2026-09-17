# Sparse random-star regularization of every full signing

2026-09-17. **Verified:** director construction, independently rebuilt
by all three researchers. This is the strongest simplified regularizer
from the paper-combination campaign. It supersedes the less efficient
replicated-field deployment, whose proof remains valid and preserved.
No external-priority claim is established.

## 1. Main theorem at the original order

There are universal constants C,c>0 such that the following holds.
Let A be ANY hollow symmetric full sign matrix of order N. Let
1<=q, m=N-q, and q^3<=m. Keep an arbitrary principal m by m
block B unchanged. Replace all edges between it and the q remaining
vertices by independent fair signs, and choose ANY full signing D
on those q vertices with Q(D)<=q^(3/2). For all sufficiently large N,
with positive probability the resulting full signing W satisfies

```math
\begin{split}
Q(W)&\le Q(A)+C N\sqrt q,\\
w\bigl(E_W(N/\sqrt q)\bigr)&\le C N q^{-1/2},\\
\log|E_W(N/\sqrt q)|&\le C N q^{-1/2}\sqrt{\log(e q)}.
\end{split}                                                    \tag{1}
```

Here E_W(T)={z:Q(W)-|H_W(z)|<=T}, and w is ordinary Gaussian
width. The last inequality is intended with a sufficiently large
universal C; the trivial N log2 bound covers small q. The construction
changes at most qN edges and has EXACT signs and zero diagonal.
It never assumes that B is optimal, spectrally flat, or has a simple
energy landscape. Labels may be selected for the given B; a single
random bridge working for all B simultaneously is not claimed.

In particular q=ceil(log N) gives, at EVERY large order,

```math
Q(W)\le Q(A)+O(N\sqrt{\log N}),\qquad
\log|E_W(N/\sqrt q)|
 =O\left(N\sqrt{\frac{\log\log N}{\log N}}\right)=o(N),
```

while only O(N log N) edges are rewritten. Any diverging q within
the stated range gives subexponential near-ground cardinality.
Taking q=floor((N/2)^(1/3)) instead gives cap cost O(N^(7/6)),
physical near-level window and width of order N^(5/6), and entropy
O(N^(5/6)sqrt(log N)). The exact window in (1) is N/sqrt(q);
implicit asymptotic constants do not enlarge that window for free.

## 2. Exact reduction to finitely many random external fields

Write C for the m by q iid sign bridge. For y in{+-1}^q set

```math
F_y(x)=\frac{|H_B(x)|}{\sqrt m}+\frac{x^TCy}{\sqrt m},
\qquad M_y=\max_xF_y(x),\qquad M_* =\max_y M_y.
```

Ignoring D, the absolute parent cap is EXACTLY sqrt(m) M_*:
the outer polarity can be absorbed into y and one may then maximize
the old polarity to obtain |H_B(x)|. Restoring D changes this cap
by at most Q(D). If an actual parent word(x,y) belongs to E_W(N/sqrt(q)),
there is a pattern y' (possibly the outer sign times y) such that

```math
F_{y'}(x)\ge M_{y'}-b,
\qquad b=\frac{N/\sqrt q+2Q(D)}{\sqrt m}
       =O(\sqrt{m/q}).                                  \tag{2}
```

Thus the old-spin projection of the full nearcode lies in the union
of ALL2^q branch nearcodes C_y={x:F_y(x)>=M_y-b}. No selected
field pattern, polarity, or new spin is omitted.

## 3. Reconstructed one-dimensional Stein replacement

For independent fair signs and real coefficients with sum a_j^2=1,

```math
d_W\left(\sum_j a_j\epsilon_j,G\right)\le\sum_j|a_j|^3. \tag{3}
```

We use the classical scalar Stein solution estimate: for a1-Lipschitz
test function h, the solution f'-zf=h-Eh(G) has ||f''||infinity<=2.
Its exact hypotheses are stated in Nathan Ross,
[Fundamentals of Stein's method, Lemma2.5 and Theorem3.1](https://arxiv.org/html/1109.1880).
This derivative estimate is an imported classical lemma, not a new theorem.
The application is reconstructed as follows. Put S=sum a_j epsilon_j
and S_j=S-a_j epsilon_j. Independence gives E[a_j epsilon_j f(S_j)]=0.
Taylor-expand f(S_j) ABOUT S; since sum(a_j epsilon_j)^2=1
pointwise, the linear terms sum exactly to E f'(S). The remainders
have total absolute expectation at most sum|a_j|^3. Stein's identity
proves(3). Smooth approximation treats nonsmooth Lipschitz tests.

For a fixed y, the field coordinates (Cy)_i/sqrt(m) are independent,
each with variance t=q/m. Equation(3) gives scalar Wasserstein error
at most1/sqrt(m) against sqrt(t)G. The map

```math
z\longmapsto\max_x[H0(x)+z\cdot x],\qquad
H0(x)=|H_B(x)|/\sqrt m,
```

is1-Lipschitz in EACH coordinate, uniformly in all offsets H0.
Replacing independent coordinates one by one therefore costs at most
sqrt(m) in expectation. The SAME statement holds after adding an
independent Gaussian field, by conditioning on that field. No nearcode
or random optimum is used as an offset in this comparison.

The sharper estimate that gives the main theorem is

```math
\sup_{a\in\mathbb R}
 \left|\mathbb E\left|q^{-1/2}\sum_{j=1}^q\epsilon_j-a\right|
            -\mathbb E|G-a|\right|\le 3/q.              \tag{3a}
```

Its full proof is [the shifted-absolute lemma](paper_bernoulli_shifted_absolute_2026_09_17.md).
Briefly the exact symmetric Stein trapezoid identity has a smooth
contribution O(1/q) and a kink contribution bounded by q^(-1/2)
times the largest atom of the (q-1)-sign sum, also O(1/q).
This is NOT an O(1/q) Wasserstein theorem. For a Boolean maximum,
conditioning on all other field coordinates gives precisely
max(A+delta z,B-delta z)=const+delta|z-a|, where A and B retain
all old energy offsets. Thus replacing the m field coordinates costs
at most 3sqrt(m/q) in expectation, including after conditioning on
an additional independent Gaussian field. The comparison is used
for each fixed pattern y, NEVER directly for the joint maximum over y.

The older generic Wasserstein comparison (3) proves a weaker,
still valid regularizer: window N, width O(Nq^(-1/4)), and range
q^2<=m. Its proof is preserved in the dated research archive.
The sharper window-width tradeoff here exploits actual Boolean
two-slope geometry rather than merely a Lipschitz bound.

## 4. A secondary field bounds the width of the realized nearcode

Put t=q/m, s=1/m, and define for each fixed y

```math
\Delta_y(C)=\mathbb E_g\max_x[F_y(x)+\sqrt s\,g\cdot x]-M_y.
```

For every realized C,

```math
\sqrt s\,w(C_y)\le b+\Delta_y(C).                         \tag{4}
```

This simply restricts the secondary maximum to C_y and pays its
deficit b. Define G(u)=E max_x[H0(x)+sqrt(u)g dot x]. Coupling
the same Gaussian vector at two amplitudes gives

```math
G(t+s)-G(t)\le\sqrt{2/\pi}\,m(\sqrt{t+s}-\sqrt t)
 \le\sqrt{2/\pi}\,ms/(2\sqrt t)=\sqrt{2/\pi}\,\sqrt{m/q}/2.
```

Applying the sharper shifted-absolute estimate (3a) below separately
to the TWO ordinary maxima defining E Delta
therefore gives, uniformly in y and B,

```math
\mathbb E\Delta_y\le(6+\sqrt{2/\pi}/2)\sqrt{m/q}.         \tag{5}
```

Changing one entry of C changes Delta_y by at most4/sqrt(m).
There are mq independent entries, so bounded differences gives
centered subGaussian proxy4q. In particular

```math
\Pr\{\Delta_y>\mathbb E\Delta_y+4\sqrt{m/q}\}
 \le e^{-2m/q^2}.
```

Union over all2^q patterns costs exp(q log2). Since q^3<=m,
with failure probability at most exp[-(2-log2)m/q^2] we have
Delta_y=O(sqrt(m/q)) simultaneously for ALL patterns. Equations(2),(4)
then give

```math
\max_y w(C_y)\le C m q^{-1/2}.                            \tag{6}
```

For any fixed family of Boolean codes, Gaussian concentration of each
supremum (Lipschitz constant sqrt m) gives

```math
w\left(\bigcup_y C_y\right)
 \le\max_y w(C_y)+\sqrt{2mq\log2}.
```

The union overhead and the q free new coordinates both fit inside
O(mq^(-1/2)) when q^2<=m, which follows from our range. This
proves the width bound in(1).
The Gaussian field is a proof tool; the final matrix is entirely signed.

## 5. Cap control and probability of simultaneous success

For each fixed pair(x,y), x^TCy has subGaussian proxy mq. A union
over2^(m+q) pairs and both signs shows, for a universal constant,

```math
\max_{x,y}|x^TCy|\le2m\sqrt q
```

except with probability exp(-c m), for all large m and q<=sqrt m.
Hence Q(W)<=Q(B)+2m sqrt(q)+Q(D). Averaging omitted spins of A
shows Q(B)<=Q(A). The cap event and Section4's event have positive
intersection (indeed probability tending to one). This proves the first
two conclusions together, for the same genuine parent.

## 6. Width-to-entropy by Gaussian information and decoding

For C subset{+-1}^d let X be uniform on C, g independent standard
Gaussian, and Y_u=sqrt(u)X+g. Its posterior mean belongs to conv(C).
Differentiating the finite posterior and integrating by parts gives

```math
\sqrt u\,\operatorname{MMSE}(u)
 =\mathbb E\,g\cdot\mathbb E[X\mid Y_u]\le w(C).
```

The Gaussian I-MMSE identity, reconstructed in the localization track,
then yields I(X;Y_u)<=w(C)sqrt(u). Equivalently the observation
X+a g has mutual information at most w(C)/a.

A nearest-codeword decoder satisfies
||Xhat-X||^2<=2a g dot(Xhat-X), whence its expected Hamming error
is at most a w(C)/2. Coordinatewise conditional binary entropy gives

```math
\frac{\log|C|}{d}\le\frac{\alpha}{a}
                  +h(a\alpha/2),\qquad \alpha=w(C)/d,
                                                                  \tag{7}
```

provided a alpha/2<=1/2. Taking a=sqrt(2/log(1/alpha)) proves
log|C|/d<=O(alpha sqrt(log(e/alpha))). Hamming balls show the
shape of this general bound is sharp up to constants. This is a
classical width/entropy phenomenon, with a directly reconstructed
information proof; it is not claimed as a new external discovery.
Together with(6), it proves the final line of(1).

## 7. Exact original-problem consequence, without overstating it

Let R_N consist of full signings satisfying(1)'s nearcode width bound
with q=ceil(log N) and its fixed universal constant. This definition
does NOT refer to the unknown optimum. The theorem gives

```math
M_N\le\min_{A\in R_N}Q(A)\le M_N+C N\sqrt{\log N}.       \tag{8}
```

So asymptotic optimization may be restricted, at every order, to this
regularized class at a vanishing normalized cost. The arbitrary large
principal block is still present; neither efficient optimization over
R_N nor a compositional law for R_N is proved. Thus this is not claimed
to resolve convergence or to make all remaining work easy.

The controlled normalized energy window is1/sqrt(Nq), NOT a fixed
positive eta. This order of limits is essential. Explicit small-width
Hadamard ground subcodes obstruct cheap isotropic response laws, and
their obstruction survives near-cap-preserving small principal extensions
at every fixed eta. Even a pure Gaussian-field landscape fails the
finite Hamming-cover slope needed by an earlier extension criterion.
Neither counterexample refutes(1); both forbid an unjustified bootstrap
from(1) to the original convergence theorem.

## 8. Proof provenance and next obligation

The discovery combined the campaign's Gaussian localization/information
view with all-offset sign/Gaussian replacement and exact full-spin
accounting. Independent audit then found the simpler, sharper scalar
Stein implementation above. We retain that simplification rather than
claiming that the common-Gibbs machinery is necessary for this gadget.

Root and all three researchers independently checked normalization,
the absolute cap identity, all2^q patterns, simultaneous concentration,
same-order restriction, and the shrinking-window scope. The earlier
replicated-field artifact remains a verified but superseded deployment.

The precise useful next question is whether a geometry-preserving operation
can use this ACTUALLY REALIZABLE regularized class without paying a fixed
leading bridge mean or requiring its full fixed-eta energy landscape.
That is open; a new cap bound or convergence claim is not manufactured.
