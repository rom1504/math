# A necessary square-root susceptibility of actual near-minimizers

Status: proved deduction, pending independent audit. This is a structural
statement about actual asymptotic minimizers, not a convergence theorem.
It makes no linear-well or confinement assumption. Date: 2026-09-07.

Let A_n be actual hollow signings with Q(A_n)/n^(3/2)-m_n ->0, and take a
subsequence on which m_n=M_n/n^(3/2) -> c>0. All limits below are along it.
For independent standard Gaussian upper entries G define

```math
S_n(e)=\{(s,x):sH_{A_n}(x)\ge Q(A_n)-e n^{3/2}\},\qquad
w_n(e)=n^{-3/2}\mathbb E\max_{(s,x)\in S_n(e)}sH_G(x),
\quad W(e)=\limsup_n w_n(e).
```

Here s=+/-1 and x is Boolean. The set is nonempty even at e=0. Its width
is nonnegative (it contains a centered Gaussian variable), is nondecreasing
in e, and is bounded by sqrt(log2)+o(1), by the Gaussian exponential-moment
bound over at most 2^(n+1) variables of variance n(n-1)/2.

## Theorem

```math
\boxed{\limsup_{e\downarrow0}\frac{W(e)^2}{e}\ge 2c.}
```

In particular W(e)=o(sqrt(e)) is impossible. This is a statement with the
order n -> infinity FIRST, then e ->0. It does not claim uniform finite-n
control at shrinking windows, or exponentially many exact ground states.

## Proof, including the interchange of maxima and expectations

The independently proved signing Gaussian-deformation inequality says,
uniformly in a signing A and for every fixed t>=0,

```math
n^{-3/2}\mathbb E Q(A+tG)
\ge\sqrt{1+t^2}(m_n-10n^{-1/6}).
```

Its proof replaces each sign by an independent biased sign of mean
A_ij/sqrt(1+t^2), then uses centered-entry Lindeberg replacement and softmax
at inverse temperature n^(1/6). See
`decisive_independent_gaussian_stability_2026_09_06.md`. This theorem uses
global signing optimality; it is not a Gaussian comparison for arbitrary
weighted matrices.

Partition signed configurations by their normalized deficit into intervals
[j d,(j+1)d], using fixed d>0 and a fixed number of bins covering [0,2C]
where Q(A_n)/n^(3/2)<=C eventually. For the bin beginning at j d,

```math
s(H_A+tH_G)\le Q(A)-j d\,n^{3/2}
                  +t\max_{S_n((j+1)d)}sH_G.
```

Each Gaussian maximum on the right is sqrt(n(n-1)/2)-Lipschitz in the
independent Gaussian entries. The Gaussian mgf/concentration inequality
therefore bounds the expectation of the maximum of the finitely many bin
values by the largest of their expectations plus
t sqrt(n(n-1) log J), where J is the number of bins. Dividing by n^(3/2)
makes this error vanish at fixed d. Consequently

```math
c(\sqrt{1+t^2}-1)
\le\max_j\{-j d+t W((j+1)d)\}.                         (1)
```

Suppose the displayed theorem were false. Choose K with K^2<2c and e0>0
such that W(e)<=K sqrt(e) for every 0<e<=e0. For small fixed t, all bins
with j d>=e0/2 have negative right-hand value, since W is uniformly bounded.
Choose d sufficiently small; on all the remaining bins the right side is
at most

```math
\sup_{e\ge0}\{-e+tK\sqrt e\}+d=K^2t^2/4+d.
```

Let d decrease to zero in (1), then t decrease to zero. The result
c/2<=K^2/4 contradicts K^2<2c. This proves the theorem.

## What is new, and what this does not say

The earlier confined-well theorem required a uniform linear deficit bound
relative to specified centers. The present consequence places a necessary
critical exponent on the ACTUAL deficit-conditioned Gaussian width with no
such premise. Cardinality alone would not give it: nearby spin clouds can
have many points but strongly correlated Gaussian responses. Conversely the
inequality still does not control a prescribed deterministic bridge, nor
identify the near-ground overlap law. A bridge can annihilate a rich cloud
unless an additional transversality statement is proved.

This does not remove the original all-order gap. A proposed use must prove
that the susceptibility forces a deterministic composition or a universal
variational law; simply naming W is not such a proof. The result is retained
as an exact optimizer-specific constraint, not promoted to a solved route.
