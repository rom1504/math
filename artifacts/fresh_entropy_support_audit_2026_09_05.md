# Entropy support audit: ordinary exponential scale is sufficient

Status: **proved correction to an archived interpretation; no large-deviation
limit or convergence theorem for the signing problem is proved here**.

Write E_n=binom(n,2), Q(A)=max_x |sum_(i<j) a_ij x_i x_j|,
g_n(A)=Q(A)/n^(3/2), and a_n=min_A g_n(A). Let

```math
Z_n(c)=\#\{A:g_n(A)\le c\},\qquad
I_n(c)=-\frac1{n^2}\log\bigl(2^{-E_n}Z_n(c)\bigr),
\qquad S_n(c)=\frac1{n^2}\log(1+Z_n(c)).
```

An empty event has I_n(c)=+infinity. It does **not** have the maximal finite
rate (log 2)/2. This elementary distinction invalidates the blanket claim
in Section 5 of `good_signing_entropy_threshold.md` that the two cases are
indistinguishable at ordinary speed n^2. A single switching orbit and many
subexponentially numerous orbits do have the same finite leading rate.

## 1. Exact lower-tail limits detect the minimum

**Theorem 1.** Suppose a_n is bounded and I_n(c) has an extended limit for
every c in a dense subset of an interval containing all cluster values of
a_n. Then a_n converges.

**Proof.** If liminf a_n<a<limsup a_n, the event defining Z_n(a) is empty
for infinitely many n and nonempty for infinitely many n. In the latter
case its probability is at least 2^(-E_n), so

```math
0\le I_n(a)\le\frac{E_n}{n^2}\log2\le\frac{\log2}{2}.
```

In the former case I_n(a)=+infinity. No extended limit exists. Choose a
from the stipulated dense set to obtain the contradiction. QED.

The theorem needs no uniformity in c, differentiability, or exclusion of
a finite maximal-rate phase. Proving those lower-tail limits is still an
unresolved mathematical task, not an automatic benefit of this criterion.

## 2. Uniform positive-entropy thickening

The following strengthens deterministic Hamming-ball thickening at a fixed
normalized relaxation. It reconstructs the planted-noise theorem in
`entropic_franz_parisi_bernoulli.md` with an explicit finite-n tail bound.

Fix 0<p<1/2. Given **any** seed signing A, independently reverse each
edge with probability p, obtaining B. Then

```math
H_B(x)=(1-2p)H_A(x)+W_x,
\qquad \operatorname{Var}(W_x)=4p(1-p)E_n.
```

Every summand in W_x is centered, has absolute value at most 2, and the
summands are independent. Bernstein and the 2^n signed/projective tests
give, for t>0,

```math
\Pr\left(\max_x|W_x|\ge t n^{3/2}\right)
\le 2^n\exp\left\{
-\frac{t^2n^3}{8p(1-p)E_n+(4/3)t n^{3/2}}
\right\}.
```

This tends to zero exponentially whenever
t>2 sqrt(p(1-p) log 2), uniformly over all seed signings.
The number K of reversed edges lies within E_n^(3/4) of p E_n with
probability tending to one. Every outcome in that shell has probability
at most

```math
\exp\left\{-E_n h(p)
+E_n^{3/4}\left|\log\frac p{1-p}\right|\right\},
```

where h(p)=-p log p-(1-p)log(1-p). Counting the good outcomes therefore
proves, uniformly over seeds with g_n(A)<=c,

```math
\log Z_n\bigl((1-2p)c+t\bigr)
\ge E_n h(p)-o(n^2).
```

**Corollary.** For every epsilon>0 there is eta(epsilon)>0 such that,
uniformly over c>=0 and all sufficiently large n,

```math
Z_n(c)>0\quad\Longrightarrow\quad
S_n(c+\varepsilon)\ge\eta(\varepsilon).
```

Choose fixed p sufficiently small that 2 sqrt(p(1-p) log2)<epsilon/2,
then choose t between that bound and epsilon. The term -2pc is favorable.
For example, sufficiently large n permits eta=h(p)/4.

This statement concerns many signings after a fixed relaxation, not a
smaller cap than the seed or entropy at the exact minimum.

## 3. Ordinary nonnegative entropy limits also suffice

**Theorem 2.** Suppose S_n(c) has a limit for every c in a dense set
containing points between any two possible cluster values of a_n. Then
a_n converges.

**Proof.** If a_-=liminf a_n<a_+=limsup a_n, choose
a_-<c_0<c_1<a_+, with c_1 in the dense set. Along a subsequence
Z_n(c_0)>0, so the corollary gives
limsup S_n(c_1)>=eta(c_1-c_0)>0. Along another subsequence
Z_n(c_1)=0 and S_n(c_1)=0. This contradicts existence of its limit. QED.

Thus speed n^(3/2) log n is **not necessary** for an entropy convergence
criterion. It was sufficient using deterministic Hamming balls; ordinary
speed n^2 suffices after the already available random-noise thickening is
included. Neither existence of the ordinary entropy limits nor a useful
formula for them has been established.

## 4. Genuine LDP versus a bulk calculation

There is also a general finite-atom statement. Let finite probability laws
nu_n have every nonzero atom at least exp(-C v_n), for a fixed C and
v_n->infinity. If their real-valued observables X_n have bounded minima
and satisfy a large-deviation upper bound on closed lower intervals and
lower bound on open lower intervals at speed v_n, with one rate function,
then their minima converge.

Indeed, a subsequence attaining a closed interval below a putative limsup
has probability at least exp(-C v_n). The LDP upper bound gives a point
of finite rate in that interval (or arbitrarily close finite-rate points).
The LDP lower bound on a slightly larger open interval then forces eventual
nonemptiness, a contradiction. Compactness/attainment of the infimum is not
needed: an approximate minimizer of the finite infimum suffices.

For uniform edge signings, v_n=n^2 and C=(log2)/2 work. One may clip
g_n above a fixed value greater than 1/2 when formulating a compact-space
LDP; the clipping does not affect the eventual support minimum. This avoids
imposing an unnecessary good-rate/exponential-tightness assertion on the
far upper tail at speed n^2.

In contrast, a replica calculation valid only where entropy is already
positive, a limsup rate, or a rate function obtained by discarding sparse
phases is not the LDP hypothesis above. Those restricted calculations may
still miss the minimum. Existing counterexamples to fixed-replica and
Gaussian replacement arguments are not contradicted by this correction.

## Research consequence

The revived precise question is existence of the **ordinary labeled
good-signing entropy limit**, or an adequate genuine lower-tail LDP,
for actual Rademacher signings at this scale. No change of scale to refined
entropy, no extra boundary orbit classifier, and no row-cavity formula is
logically required merely to imply convergence. These might still be
needed to prove the limit; no such necessity is currently established.

The director derived these criteria from first principles and the source
researcher independently checked the atom-floor argument and uniform
Bernstein/shell calculation. This removes an erroneous restriction from
the proof-obligation map. It does not itself determine the original limit.
