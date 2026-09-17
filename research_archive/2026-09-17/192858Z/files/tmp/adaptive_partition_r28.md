# Wave 28 structural attack: adaptive partitions and the global collision wall

## Outcome

The adaptive-partition idea has a clean **local** theorem: any one aligned
low-row cut can be planted in an entire balanced block coset whose row cap is
of the same asymptotic order. Thus adapting the partition costs no power for
one selector.

The global quantifier in (10.837) is substantially stronger. I prove an exact
signature-class criterion for many selector witnesses to share one planted
block coset. It exposes two additional requirements:

1. the witness family must use at most \(k\) coordinate sign profiles (with
   the precise balanced version below); and
2. the **whole blockwise closure** of the witnesses, not merely the planted
   witnesses, must obey the row cap.

The second issue is real even for an audited exact minimizer. In \(A_8\),
three individually favorable exact-child completions have row squares
\(64,64,40\), but every block coset containing all three also contains their
ternary product, whose row square is \(72\). This is a rigorous scoped wall
to choosing low-row child completions separately and then grouping their sign
profiles. It does not rule out a more correlated choice of completions.

The checker is `tmp/adaptive_partition_r28_check.py`. All finite indexing is
zero-based.

## 1. A local planted-partition theorem

### Verified proposition

Let \(A\) be any symmetric zero-diagonal sign matrix, let
\(x\in\{\pm1\}^n\), and let \(2\le k\le n\), with
\(n/k\gg\log k\). There is a partition \(P\) into \(k\) nonempty blocks,
each of size at most \(2n/k\), such that, for \(D=\operatorname{diag}(x)\),

\[
\max_{z\in\{\pm1\}^k}\|ADPz\|_2^2
\le
2R_2(x)
+8u\max\{k\|A\|_{\rm op}^2,n(n-1)\}
+\frac{16}{9}knu^2,
\qquad
u=\log(4(n+k)).
\tag{AP.1}
\]

In particular, for an exact minimizer, fixed \(0<c<1/4\), and

\[
k\asymp\frac{n^{3/4-c}}{\log n},
\]

the right side of (AP.1) is

\[
2R_2(x)+O(n^{9/4-c}).
\tag{AP.2}
\]

### Proof

Hash the coordinates independently and uniformly into \(k\) bins. Put
\(v_i=x_iAe_i\), \(g(i)\) for the bin of \(i\), and

\[
M=ADP=\sum_i v_i e_{g(i)}^{\mathsf T}.
\]

Then

\[
\mathbb EM=\frac1k(Ax){\bf1}^{\mathsf T},
\qquad
\|\mathbb EM\|_{\rm op}=\frac{\|Ax\|_2}{\sqrt k}.
\]

For
\(Z_i=v_i(e_{g(i)}-k^{-1}{\bf1})^{\mathsf T}\), the exact rectangular
Bernstein parameters are

\[
\|Z_i\|_{\rm op}\le\sqrt n,
\]

\[
\left\|\sum_i\mathbb EZ_iZ_i^{\mathsf T}\right\|_{\rm op}
=(1-1/k)\|A\|_{\rm op}^2,
\]

\[
\left\|\sum_i\mathbb EZ_i^{\mathsf T}Z_i\right\|_{\rm op}
=\frac{n(n-1)}k.
\]

The standard inverted Bernstein bound therefore gives, with probability at
least \(3/4\),

\[
\|M-\mathbb EM\|_{\rm op}
\le \sqrt{2\nu u}+\frac23\sqrt n\,u,
\quad
\nu=\max\{(1-1/k)\|A\|_{\rm op}^2,n(n-1)/k\}.
\]

Chernoff and a union bound give nonempty bins of size at most \(2n/k\)
with probability tending to one. The two events intersect. Finally,
\(\|z\|_2^2=k\), the triangle inequality, and
\((a+b)^2\le2a^2+2b^2\) give (AP.1). For exact minimizers use the recorded
\(\|A\|_{\rm op}^2\le2q_n=O(n^{3/2})\). This proves (AP.2).

### Consequence for one selector

If \(x\) extends a child Boolean-norm optimizer on \(S\), then

\[
Q(P^{\mathsf T}DH_SDP)
\ge |x^{\mathsf T}H_Sx|
\ge Q(A[S])-p_2q_n
=Y_A(S)+B_{n,m}.
\tag{AP.3}
\]

Thus any such completion with
\(R_2(x)=O(n^{9/4-c})\) can be embedded in a whole balanced, row-good
coset and aligns its selector with margin \(B_{n,m}\).

**Quantifier warning.** Here \(P,D\) depend on \((S,x)\). This proves
neither (10.837), where one \(P,D\) must work for all selectors, nor the
uniform hit probability (10.838).

## 2. Exact criterion for a family to share a planted block coset

### Verified signature theorem

Let

\[
\mathcal F=\{x^{(0)},x^{(1)},\ldots,x^{(r-1)}\}
\subset\{\pm1\}^n
\]

be spin representatives. Associate to coordinate \(i\) its relative
signature

\[
\tau_i=
\bigl(x_i^{(1)}x_i^{(0)},\ldots,
x_i^{(r-1)}x_i^{(0)}\bigr)\in\{\pm1\}^{r-1}.
\tag{AP.4}
\]

Changing a projective representative only relabels signature values, so the
equality classes are intrinsic. Let the nonempty signature classes have
sizes \(N_1,\ldots,N_L\).

For a partition \(P\), there exists a diagonal \(D\) whose block coset
contains \(\mathcal F\) if and only if every block of \(P\) is contained in
a signature class. Consequently:

- without a block-size constraint, a coset with at most \(k\) blocks
  contains \(\mathcal F\) if and only if \(L\le k\);
- with maximum block size \(b\), the minimum possible number of blocks is

  \[
  K_b(\mathcal F)=\sum_{\ell=1}^L\left\lceil\frac{N_\ell}{b}\right\rceil.
  \tag{AP.5}
  \]

  Thus membership in some at-most-\(k\)-block coset is equivalent to
  \(K_b(\mathcal F)\le k\).

To prove this, suppose \(x^{(a)}=DPz^{(a)}\). Absorb the block signs
\(z^{(0)}\) into \(D\), so that \(D=\operatorname{diag}(x^{(0)})\).
On each block,
\(x_i^{(a)}x_i^{(0)}=z_{g(i)}^{(a)}z_{g(i)}^{(0)}\) is constant.
This proves necessity. Conversely, use the signature classes, put
\(D=\operatorname{diag}(x^{(0)})\), and assign to each block the
corresponding component of (AP.4). Splitting classes proves (AP.5).

The signature-class partition is the coarsest partition whose coset contains
\(\mathcal F\). Every other containing block coset is a refinement and
therefore contains this minimal coset as a subset. Hence refining to enforce
balance cannot remove a bad row-square word.

### Verified row-closure conditions

Every block coset is ternary closed:

\[
x,y,w\in\mathcal C(D,P)
\quad\Longrightarrow\quad
x\odot y\odot w\in\mathcal C(D,P).
\tag{AP.6}
\]

Indeed \(D_i^3=D_i\), and the product of the three block signs is again
block-constant. Thus individual row bounds on the planted family do not
control the row bound of a containing coset; at a minimum, the entire affine
closure under (AP.6) must be checked.

There is also a deterministic field obstruction. Normalize
\(D=\operatorname{diag}(x^{(0)})\), and let a containing partition split
signature class \(C_\ell\) into \(k_\ell\) blocks. With

\[
V_\ell=A D{\bf1}_{C_\ell},
\]

uniform averaging over all block signs and then Cauchy--Schwarz give

\[
\begin{aligned}
\max_z\|ADPz\|_2^2
&\ge\sum_{B\in P}\|AD{\bf1}_B\|_2^2\\
&\ge\sum_{\ell=1}^L\frac{\|V_\ell\|_2^2}{k_\ell}
\ge\frac1k\left(\sum_{\ell=1}^L\|V_\ell\|_2\right)^2,
\end{aligned}
\tag{AP.7}
\]

whenever the partition has at most \(k\) blocks. Therefore a row cap \(C\)
requires the exact signature-field total-variation condition

\[
\sum_{\ell=1}^L
\|A\operatorname{diag}(x^{(0)}){\bf1}_{C_\ell}\|_2
\le\sqrt{kC}.
\tag{AP.8}
\]

This condition is only necessary; the exact test is the Boolean row maximum
of a balanced refinement of the signature partition.

## 3. Exact \(A_8\) child-completion collision wall

### Verified finite counterexample to independent planting

For the audited exact minimizer \(A_8\), take the spin words

\[
x^{(0)}=(+------+),\qquad
x^{(1)}=(+-----+-),\qquad
x^{(2)}=(+-+---+-).
\]

Their exact row squares are

\[
\bigl(R_2(x^{(0)}),R_2(x^{(1)}),R_2(x^{(2)})\bigr)=(64,64,40).
\tag{AP.9}
\]

They are not arbitrary eligible cuts. On the respective four-sets

\[
S_0=\{0,1,2,6\},\qquad
S_1=\{0,1,2,7\},\qquad
S_2=\{0,1,2,3\},
\tag{AP.10}
\]

their restrictions attain the exact child Boolean norm \(Q(A_8[S_a])=8\).
Thus each is an exact-child completion and is loss-favorable by (10.825).

But their ternary product is

\[
x^{(0)}\odot x^{(1)}\odot x^{(2)}
=(+-+----+),
\qquad R_2=72.
\tag{AP.11}
\]

By (AP.6), **every** block coset containing all three completions violates
the row cap \(64\). More explicitly, their three signature classes are

\[
\{0,1,3,4,5\},\qquad\{2\},\qquad\{6,7\},
\]

and the minimal containing coset has exact row-square values
\(\{40,64,72\}\).

This is a scoped falsification:

- it disproves the rule “choose individually low-row child completions and
  then group equal sign profiles; the resulting coset is automatically
  row-good”;
- it does **not** show that these selectors lack a different correlated
  choice of completions, and it does not falsify (10.837). Indeed the earlier
  finite pair-block audit shows other \(A_8\) cosets can work at cap \(64\).

## 4. Independent rounding does not create global collisions

### Verified probabilistic wall

For \(r\) independent uniform full spins, the relative coordinate signatures
in (AP.4) are iid uniform on

\[
N=2^{r-1}
\]

possible values. Projectivizing the spins merely relabels all signatures and
leaves the number of occupied values unchanged. If \(L\) is the number of
occupied signatures, then for \(k<N\)

\[
\Pr\{L\le k\}
\le {N\choose k}\left(\frac{k}{N}\right)^n.
\tag{AP.12}
\]

Choose \(r=1+\lceil\log_2(2k)\rceil\), so \(2k\le N<4k\). If
\(k=o(n)\), then

\[
\Pr\{L\le k\}
\le (eN/k)^k2^{-n}=e^{-\Omega(n)}.
\tag{AP.13}
\]

Thus only \(O(\log n)\) independently randomized witnesses already require
more than \(k\asymp n^{3/4-c}/\log n\) blocks with overwhelming probability.
Independent outside rounding from (10.826) can be useful for finding one
low-row completion, but independence itself supplies the wrong collision
geometry. Actual child-conditioned completions are not uniform full spins,
so (AP.13) is a scoped wall to the naive random-completion mechanism, not a
falsification of every correlated selector rule.

## 5. Exact global reformulation and remaining target

For each selector \(S\), choose a spin \(x_S\) satisfying

\[
|x_S^{\mathsf T}H_Sx_S|\ge Y_A(S)-t.
\tag{AP.14}
\]

Then (10.837) is exactly the problem of choosing this exponentially indexed
family so that:

1. its signature classes admit a balanced refinement with at most
   \(k=O(n^{3/4-c}/\log n)\) blocks; and
2. the **entire** resulting block coset, not just the selected \(x_S\)'s,
   has \(R_2=O(n^{9/4-c})\).

The local theorem proves this when the family contains one witness. The
finite wall and (AP.12) show why separately chosen completions do not scale
to the required global family.

**Open structural target.** A genuine child-field construction must choose
the completions jointly so that their relative sign vectors take only
\(O(k)\) coordinate profiles *and* the signature-field obstruction (AP.8)
and the full Boolean row maximum are controlled. Equivalently, it must
produce a low-complexity correlated completion rule, not merely a low-row
completion for every selector. A weaker and perhaps more realistic route is
still (10.838): prove that row-good cosets obtained without planting have
stretched-exponential alignment probability, so that natural selector
coverage supplies the collisions.

No convergence proof or asymptotic counterexample is obtained.
