# Wave 29 route 1: coherent Johnson selections and signature closure

## Status

- **Verified (analytic):** an exact spanning-tree formula identifies the
  signature classes in (10.856), and gives a global-variation sufficient
  condition for \(K_b(\mathcal F)\le k\).
- **Verified (exact finite exhaustive computation):** for the named exact
  minimizer \(A_9\), a triangle in \(J(9,6)\) is pairwise coverable by
  balanced, cap-\(80\), exact-child block cosets, but is not jointly
  coverable by any such coset. The least common cap is \(88\).
- **Falsified mechanism:** adjacency-wise Lipschitz selection, pairwise Hall
  feasibility, or the bare signature count \(K_b\le k\) does not by itself
  control the row cost of the required full coset closure.
- **Open:** an asymptotic jointly selected family satisfying both the global
  signature budget and full low-row closure. The finite wall below does not
  falsify (10.837).

The computation is checked by
`tmp/coherent_signatures_r29_search.py`; it uses only integer arithmetic
and no floating-point operations.

## 1. Exact spanning-tree description of the signatures

Let \(G=(V,E)\) be a connected graph indexing coherently oriented full
spins \(x^v\in\{\pm1\}^n\). (For projective spins, choose a root
orientation and orient successively down a spanning tree.) Fix a spanning
tree \(T\), root \(v_0\), and for \(e=uv\in T\) define

\[
 F_e=\{i:x_i^u\ne x_i^v\},\qquad
 \eta_i=\bigl({\bf1}_{i\in F_e}\bigr)_{e\in T}.
\tag{R29.1}
\]

Then the coordinate signature classes of (10.856) are *exactly* the level
sets of \(i\mapsto\eta_i\). In particular, if
\(N_\eta=|\{i:\eta_i=\eta\}|\), then

\[
 J=\bigl|\{\eta_i:i\in[n]\}\bigr|,
 \qquad
 K_b(\mathcal F)=\sum_{\eta:N_\eta>0}
     \left\lceil\frac{N_\eta}{b}\right\rceil .
\tag{R29.2}
\]

**Proof.** Along the unique root-to-\(v\) path in \(T\),

\[
 x_i^v x_i^{v_0}
 =(-1)^{\sum_{e\in[v_0,v]}{\bf1}_{i\in F_e}}.
\tag{R29.3}
\]

Thus equal edge profiles imply equal relative signatures. Conversely,
the change between the relative signs at the endpoints of each tree edge
recovers \({\bf1}_{i\in F_e}\), so equal signatures imply equal edge
profiles. This proves (R29.2).

Put

\[
 U_T=\bigcup_{e\in T}F_e,\qquad
 u_T=|U_T|,\qquad D_T=\sum_{e\in T}|F_e|.
\]

All coordinates outside \(U_T\) have the zero profile, while the nonzero
classes contain a total of \(u_T\) coordinates. Consequently

\[
 J\le1+u_T\le1+D_T,\qquad
 K_b(\mathcal F)
 \le\left\lceil\frac{n-u_T}{b}\right\rceil+u_T
 \le\left\lceil\frac n b\right\rceil+D_T.
\tag{R29.4}
\]

For projective spins one may orient the child of each tree edge to make
\(|F_e|=\min\{d_H(x^u,x^v),n-d_H(x^u,x^v)\}\); there is no cycle
constraint on a tree. Thus (R29.4) is a genuine sufficient signature
criterion. It is global: a bound \(d_H(x^u,x^v)\le h\) on every edge
only gives \(D_T\le h(|V|-1)\), which is useless when \(V\) contains the
whole selector slice unless the selection is constant on large regions.

There is a sharp mechanism counterexample on the Johnson graph. For
\(n\ge3\) and \(1\le m<n\), set

\[
 x_i^S=\begin{cases}+1,&i\in S,\\-1,&i\notin S.\end{cases}
\qquad(S\in\tbinom{[n]}m).
\tag{R29.5}
\]

Adjacent selectors give spins at Hamming distance exactly \(2\), but the
\(n\) coordinate signatures are all distinct. For any \(i\ne j\), the
parity of membership of \(i,j\) is not constant on the slice: one selector
contains exactly one, while another contains both (if \(m\ge2\)) or neither
(if \(m\le n-2\)). Thus \(x_i^Sx_j^S\) varies with \(S\), which is exactly
what rules out equality of the two relative signatures. Hence
\(J=K_b=n\), for every
\(b\). A uniform local Lipschitz constant therefore does not imply the
sublinear signature budget in (10.856). This example concerns signature
geometry only; it does not assert that the spins in (R29.5) are child
optimizers for an exact minimizer.

## 2. An exact pairwise-but-not-global wall in \(A_9\)

For a selector \(S\), define the projective exact-child completion set

\[
 {\cal W}_C(S)=
 \{x\in\{\pm1\}^9/\{\pm1\}:
   |x_S^{\mathsf T}A_9[S]x_S|=Q(A_9[S]),\
   R_2(x)\le C\}.
\tag{R29.6}
\]

Consider the three pairwise adjacent selectors, with common five-set
\(\{0,1,3,4,5\}\),

\[
\begin{aligned}
S_0&=\{0,1,2,3,4,5\},\\
S_1&=\{0,1,3,4,5,6\},\\
S_2&=\{0,1,3,4,5,7\}.
\end{aligned}
\tag{R29.7}
\]

In all three child matrices \(Q(A_9[S_a])=18\). With projective
representatives normalized to start with `+`, exhaustive enumeration
gives

| selector | all words in \({\cal W}_{80}(S)\), with \(R_2\) |
|:--|:--|
| \(S_0\) | `+-+++++++` (72), `+-+---++-` (80) |
| \(S_1\) | `+--+++-++` (80), `+--+++--+` (80) |
| \(S_2\) | `+--++++-+` (72), `+--+++--+` (80) |

Every edge of this Johnson triangle is feasible in the full finite analogue
of the block architecture: five blocks of sizes \(1,2,2,2,2\), and every
word of the coset has row square at most \(80\). Explicit certificates
are:

| edge | partition | coset representative | exact-child witnesses |
|:--|:--|:--|:--|
| \(S_0S_1\) | `(0) (1,5) (2,6) (3,8) (4,7)` | `+-+++++++` | `+-+++++++`, `+--+++-++` |
| \(S_0S_2\) | `(1) (0,6) (2,4) (3,8) (5,7)` | `++-++-+++` | `+-+---++-`, `+--++++-+` |
| \(S_1S_2\) | `(0) (1,2) (3,4) (5,8) (6,7)` | `+++++++++` | `+--+++--+` for both |

Direct enumeration of the words in each displayed coset verifies the
cap \(80\).

Nevertheless:

\[
\boxed{\text{No five-block, maximum-block-size-two coset of row cap \(80\)
contains a member of every \({\cal W}_{80}(S_a)\).}}
\tag{R29.8}
\]

The exhaustive universe has

\[
\frac{9!}{(2!)^4\,4!}=945
\]

unlabelled \(1+2+2+2+2\) partitions and \(2^{9-5}=16\) diagonal
cosets modulo block flips for each partition, hence \(15120\) cases.
Exactly \(452\) cosets have all their words at row square at most \(80\).
Among the \(504\) Johnson star triangles having a common five-set, \(122\)
are pair-covered by these cosets and six, including (R29.7), are not
triple-covered. All coverage tests compare
integer child energies and integer row squares.

The gap is exact (the complete \(A_9\) row-square spectrum jumps from
\(80\) to \(88\)):

\[
\boxed{\min_{\substack{\text{five max-size-two block cosets}\\
\text{hitting exact grounds on }S_0,S_1,S_2}}
\ \max_{x\text{ in the coset}}R_2(x)=88.}
\tag{R29.9}
\]

For example, the witnesses

\[
 x^{S_0}=\texttt{+-+++++++},\qquad
 x^{S_1}=\texttt{+--+++-++},\qquad
 x^{S_2}=\texttt{+--+++--+}
\]

fit the partition

```text
(7) (0,3) (1,5) (2,6) (4,8)
```

and its full coset has maximum row square \(88\).

The obstruction is stronger than a failed signature count. Choosing

\[
x^{S_0}=\texttt{+-+++++++},\qquad
x^{S_1}=x^{S_2}=\texttt{+--+++--+}
\tag{R29.10}
\]

produces only two signature classes, of sizes \(6\) and \(3\), so

\[
K_2=\lceil6/2\rceil+\lceil3/2\rceil=5.
\]

Its projective ternary span consists only of the two displayed low-row
words. Yet every required maximum-size-two refinement introduces another
block word of row square at least \(88\). Thus even the exact signature
budget and low-row ternary span do not guarantee low-row *full coset*
closure.

## 3. What this rules out, and what remains open

The local statement is at the single finite tuple
\((A,n,m,C,b,k)=(A_9,9,6,80,2,5)\). It proves that a coupling theorem
whose hypotheses only check adjacent pairs (even by exhibiting a fully
row-good balanced coset on every edge) cannot conclude one common
row-good coset. Higher-order/refinement control is indispensable.

It does **not** falsify (10.837):

1. (R29.8) asks the common coset to hit an *exact child ground* at every
   selector. Alignment in (10.837) permits arbitrary oriented cuts that
   meet its positive tolerance and retained \(B_{n,m}\) allowance.
2. The cap \(80\) is a finite normalization; cap \(88\) already removes
   this example.
3. One finite exact minimizer does not address the required asymptotic
   ratio window.
4. It says nothing against the multi-coset hit law (10.838).

The global deterministic target has the following quantifiers. One needs
fixed \(\rho\in[1/2,1)\), \(0<c<1/4\), and uniform constants such that, for every
sufficiently large \(n\) and every integer \(m\in[\rho n,n)\), one may
choose a target-specific exact order-\(n\) minimizer \(A\), an
\((A,m)\)-adapted balanced partition \(P\), and one diagonal \(D\), with

\[
k=O(n^{3/4-c}/\log n),\quad
b=O(n^{1/4+c}\log n),\quad
\max_z R_2(DPz)=O(n^{9/4-c}),
\tag{R29.11}
\]

such that the alignment inequality (10.837), with
\(t=O(n^{3/2-c})\), holds for **every** \(m\)-selector. The constants,
\(\rho\), and \(c\) cannot depend on \(n,m,S\); \(A,P,D\) may depend on
\((n,m)\), and the block word and orientation witnessing alignment may
depend on \(S\).

A clean successor would therefore need both:

1. a globally coupled witness selection whose spanning-tree change union
   makes (R29.4) fit the block budget (local edge bounds alone do not); and
2. a direct full-refinement row theorem, with slack, for the resulting
   signature classes (the \(A_9\) wall shows that low-row witnesses and
   their ternary span are insufficient).

No such theorem is proved here.
