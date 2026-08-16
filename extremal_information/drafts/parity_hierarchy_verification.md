# Independent verification of the parity hierarchy

Date: 2026-08-16

Verdict: **CORRECT**. Theorem 5.2 is mathematically correct. Its final
interpretation should be strengthened by stating an explicit Cartesian-power
corollary: the theorem as written rules out exact recovery of covering radius
from any fixed unrooted replica order, while the power corollary is what rules
out recovery of the *normalized* radius to vanishing error along sequences.

## 1. Construction and proper-subconfiguration identity

Let $r=2s+1$, $M=r+1$, $N=2^{r-1}$, and

\[
P_\epsilon=\{v\in\mathbb F_2^r:{\bf1}\cdot v=\epsilon\}.
\]

The code $A_\epsilon\subseteq\mathbb F_2^N$ has labels
$0,a_1,\ldots,a_r$, with $a_i(v)=v_i$. Every nonzero function among
$a_i$ and $a_i+a_j$ is nonconstant on $P_\epsilon$, so every two
distinct codewords are at distance $N/2$.

For a proper label set
$S\subsetneq\{0,a_1,\ldots,a_r\}$, translate one selected label to zero.
The other selected words are evaluations of a linear space $W$ of forms on
$v$. The parity form $\ell(v)={\bf1}\cdot v$ is not in $W$:

- if $0\in S$, at least one coordinate form is absent, so the all-one
  coefficient vector is not spanned;
- if $0\notin S$, the translated forms have even coefficient weight,
  whereas the all-one coefficient vector has odd weight.

Therefore every fiber of the selected-pattern map meets $P_0$ and $P_1$
equally. The translated column-pattern multisets, hence the selected cube
subconfigurations, are isometric for $A_0$ and $A_1$. I independently
enumerated all such label subsets for $r=3,5,7$; every comparison agreed.

There is a small proof detail worth making explicit in the source report.
For a prescribed membership pattern in $T_t$, inclusion--exclusion over
the positions declared outside the code reduces the count to counts with a
set of positions constrained to codewords. If their tuple uses fewer than
$M$ distinct codewords, the preceding cube isometry identifies the number
of ambient completions. If it uses all $M$ codewords, then $t=M$ and
there are no unconstrained positions; the recorded metric is just the common
equilateral metric. This proves

\[
T_t(A_0)=T_t(A_1)\qquad(1\le t\le M).
\]

Thus the claimed ambient census equality is stronger than merely equality
of intrinsic codeword moments and is valid.

## 2. Radius calculation

At a coordinate $v$ of weight $w$, the $M$ codeword bits contain
$w$ ones and $M-w$ zeros. Consequently every root contributes at most
$\max\{w,M-w\}$ to the sum of its $M$ distances, proving

\[
\rho(A_\epsilon)\le
\left\lfloor S_\epsilon/M\right\rfloor,
\qquad
S_\epsilon=\sum_{w\equiv\epsilon(2)}
\binom rw\max\{w,M-w\}.
\]

Put $p=s\bmod2$. The tie layer has weight
$M/2=s+1\equiv1-p\pmod2$, so it is absent from $P_p$. Choosing the root
bit opposite the strict majority maps the coordinates bijectively to the
parity-$p$ subsets of the $M$ codeword labels having size greater than
$M/2$. This family is permutation invariant. Hence every label occurs
equally often and all root distances equal $S_p/M$. In particular $S_p/M$
is an integer and

\[
\rho(A_p)=S_p/M.
\]

Pairing weights $w$ and $2s+1-w$ gives

\[
S_0-S_1
=\sum_{w=0}^s(-1)^w\binom{2s+1}{w}
=(-1)^s\binom{2s}{s}.
\]

Therefore

\[
S_p-S_{1-p}=\binom{2s}{s}>0,
\]

and the displayed upper bound for $A_{1-p}$ is strictly below $S_p/M$.
The radii differ exactly as claimed. Independent exact checks gave

| $r$ | $N$ | $p$ | $(S_0,S_1)$ | verified radii $(\rho(A_0),\rho(A_1))$ |
|---:|---:|---:|---:|---:|
| 3 | 4 | 1 | (10, 12) | (2, 3) |
| 5 | 16 | 0 | (66, 60) | (11, 10) |
| 7 | 64 | 1 | (316, 336) | (39, 42) |

The first two rows were checked by exhaustive root enumeration; the third
was checked by an exact binary linear optimization. The proof itself does
not rely on these computations or on equality in the upper bound for the
nonselected parity.

## 3. Quantifiers and the normalized consequence

For every fixed $k$, choose odd $r$ with $M=r+1\ge k$. The preceding
argument gives $T_t(A_0)=T_t(A_1)$ for every $t\le k$ and unequal radii.
This already proves that no fixed unrooted replica hierarchy determines the
exact covering radius universally.

To obtain the stronger asymptotic statement relevant to normalized extrema,
one should add the following corollary. For fixed $r$, take Cartesian
powers $A_0^m,A_1^m\subseteq Q_{mN}$. A block $t$-tuple has a membership
vector $b\in\{0,1\}^t$ and a pair-distance vector $d$. Under products,
these combine by

\[
(b,d)\star(b',d')=(b\wedge b',d+d').
\]

Thus equality of the base $T_t$ distributions implies equality of their
$m$-fold convolutions, so

\[
T_t(A_0^m)=T_t(A_1^m)\qquad(t\le k,\ m\ge1).
\]

Meanwhile covering radius is additive under Cartesian products:

\[
\rho(A_\epsilon^m)=m\rho(A_\epsilon).
\]

Consequently the normalized radii remain separated by

\[
\frac{|\rho(A_0)-\rho(A_1)|}{N}>0
\]

for all $m$. This is the precise theorem showing that no fixed collection
of full unrooted $t$-replica censuses preserves normalized covering radius
uniformly to $o(1)$.

## 4. Exact correction requested

No formula or construction in Theorem 5.2 needs correction. Before using its
last sentence as an asymptotic information lower bound, add the
Cartesian-power corollary above. Also expand the one-line
“ambient-completion and inclusion--exclusion” step by separating the case of
all $M$ distinct constrained codewords, as in Section 1; this removes a
minor logical omission without changing the result.
