# Designer pressure audit and loss of entropy under local limits

Status: verified elementary statements; no new convergence theorem. This
records the initial independent proposal, its archival duplication, and a
precise limitation on a possible asymptotic variational principle.

## 1. Initial proposal and archive comparison

Before reading the synthesis I proposed, with `g_n(A)=Q(A)/n^(3/2)` and
`m_n=min_A g_n(A)`,

```math
F_n(b)=-{1\over b n^2}\log\mathbb E_U e^{-b n^2g_n(A)}.
```

The exact bounds are `m_n <= F_n(b) <= m_n+binom(n,2)log(2)/(bn^2)`.
Limits of this pressure at an unbounded collection of fixed positive b
would imply convergence. The proposal is already precisely present in
`signing_free_energy_limit.md`; it is NOT new. That artifact establishes
the changing-temperature obstruction to straightforward Shearer restriction.
The more recent `fresh_entropy_support_audit_2026_09_05.md` correctly shows
that a genuine n²-speed lower-tail LDP suffices, including sparse phases.
There is no residual need to demand a finer speed merely to detect minima.

## 2. Every bounded-cap law has the same finite-edge limit

Let A be any hollow signing with `Q(A)<=C n^(3/2)`. Randomly relabel its
vertices, uniformly. For every fixed k, the induced labeled signing on k
vertices converges in total variation to independent fair edge signs,
uniformly over all such A.

Here is a direct proof, including its scale limitation. Write
`q(z)=z^T A z`. Since A is hollow, extrema of q on the cube occur at
vertices, so `|q(z)|<=2Q(A)` for z in [-1,1]^n. For u,v in that cube,
polarization gives

```math
u^TAv=q((u+v)/2)-q((u-v)/2),\qquad |u^TAv|\le4Q(A).
```

Represent `(1+A_ij)/2` off the diagonal as a graphon W on n equal
intervals, and set the n diagonal blocks to 1/2. Its cut distance from
the constant graphon 1/2 is at most `2Q(A)/n²`, hence at most `2C/sqrt(n)`.
For any specified induced graph on k labeled vertices, telescope the
product of its edge/nonedge factors against the constant graphon. In each
term, condition on the other k-2 integration variables. The factors
involving the two endpoints split as f(x)g(y), with f,g in [0,1]; all
remaining factors are bounded by one. The layer-cake formula bounds the
remaining integral by the cut norm. Consequently its pattern probability
differs from `2^(-binom(k,2))` by at most `2C binom(k,2)/sqrt(n)`.
Sampling distinct vertices instead of independent intervals costs at most
`binom(k,2)/n`. Thus, conservatively,

```math
\mathrm{TV}\le 2^{\binom{k}{2}-1}
\left({2C\binom{k}{2}\over\sqrt n}+{\binom{k}{2}\over n}\right).
```

This proves the fixed-k assertion; the displayed estimate only gives a
growing-k conclusion on approximately the sqrt(log n) scale, not a
macroscopic or mesoscopic entropy theorem. The diagonal-block convention
has no effect on distinct-vertex samples.

The same conclusion holds for any vertex-exchangeable probability law
supported on bounded-cap signings, by conditioning on A. More generally,
it holds for laws with uniformly bounded expected normalized cap, by the
same bound averaged over A (or truncation).

In particular this applies to the ACTUAL designer Gibbs law
`P_(n,b)(A) proportional to exp(-b n²g_n(A))` at every fixed b>0.
The Gibbs identity is

```math
F_n(b)=\mathbb E_{P_{n,b}}g_n+
{D(P_{n,b}\Vert U_n)\over bn^2}.
```

Thus its expected cap is at most `m_n+log(2)/(2b)`, uniformly bounded.
The law is vertex-exchangeable. Its finite-edge limit is therefore the
same iid law for EVERY b>0, regardless of any pressure nonconvergence.
An ordinary exchangeable-array state space collapses the entire actual
designer-temperature family to one point, not only the two examples below.

## 3. Different entropy densities within that identical local limit

Choose any deterministic bounded-cap sequence A_n. Let P_n be its uniform
random-relabeling law. Then

```math
\operatorname{Ent}(P_n)\le\log(n!)=o(n^2).
```

Fix p in (0,1/2). Independently flip each edge of the relabeled matrix with
probability p, yielding R_n. Conditional on the relabeling its entropy is
exactly `binom(n,2) h(p)`. Entropy cannot decrease on averaging the
conditional laws, and cannot exceed conditional entropy plus label
entropy. Hence

```math
{\operatorname{Ent}(R_n)\over n^2}\longrightarrow {h(p)\over2}.
```

The standard independent-noise union bound gives, for every fixed
`t>2sqrt(p(1-p)log2)`,

```math
Q(B)/n^{3/2}\le (1-2p)C+t
```

with probability `1-exp(-Omega(n))`, uniformly over the initial A_n.
Condition R_n on that bounded-cap event to obtain R'_n. Its entropy differs
from that of R_n by o(n²): for an event of probability 1-epsilon, the
entropy chain rule bounds the difference by
`h(epsilon)+epsilon binom(n,2)log2`, up to the harmless normalization of
the conditional entropy. Thus R'_n has entropy density h(p)/2 and is
supported on bounded-cap signings.

Both P_n and R'_n have exactly the same iid-fair finite-edge limit, but
their entropy densities differ by h(p)/2. Taking p small makes the cap
relaxation arbitrarily small while retaining a strictly positive entropy
density difference. Therefore the limiting finite-dimensional exchangeable
array, or the collection of fixed induced-pattern densities, does not
determine the entropy term in the designer variational principle—even
inside an arbitrarily small fixed cap relaxation of a minimizing sequence.

This is a no-go for a specific choice of variational state, NOT a no-go for
all variational principles, an entropy-limit theorem, or nonconvergence.
It indicates an exact missing obligation: any successful compact state
must retain enough growing-order information to recover the n² entropy
without losing the normalized cap or all-order realizability.

## 4. Primary literature check

Chen, Guionnet, Ko, Lacroix-A-Chez-Toine, Mourrat, *One-sided large
deviations for the ground-state energy of spin glasses*, March 6, 2026:
https://perso.ens-lyon.fr/jean-christophe.mourrat/gs_ldp.pdf
(also arXiv:2603.06368).

Theorem 1.2 is for deviations ABOVE the typical maximum, at speed N,
for the Gaussian covariance model (1.1). Theorem 2.1 uses fractional
moments with exponent s in (0,1). Neither theorem gives the Bernoulli
LOWER-tail n²-speed pressure here. The paragraph following (1.8), page 3,
explicitly describes n² lower-tail speed in zero external field as
expected, and points to spherical results, not an Ising theorem. The
absolute two-sided cap is also not its one-sided observable. No import is
made from this source. The theorem statement and explanatory paragraph
were read directly from the author-hosted PDF.

## 5. Failed lower-tail-speed shortcut

The sequential greedy spin construction for an iid signing has energy
distributed as `sum_(i=1)^(n-1) |S_i|`, where S_i is an i-step simple random
walk and these variables are independent across newly exposed rows. Its
mean is asymptotic to `(2/3)sqrt(2/pi)n^(3/2)`, above 1/2. A lower deviation
of order n^(3/2), however, controls only n independent row increments and
n-scale exponential cost. It does not certify the n²-scale lower tail.
This shortcut was considered and rejected; no numerical evidence was
used to bridge that missing factor n.
