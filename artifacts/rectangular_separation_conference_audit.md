# Conference audit of the proposed rectangular-separation lemma

## Setup

Write a symmetric signing in blocks
\[
A=\begin{pmatrix}B&C\\ C^\top&D\end{pmatrix},
\qquad |S|=m,\quad |S^c|=N-m.
\]
Choose a positive extremizer \(u\) and a negative extremizer \(v\) of \(B\), and set
\[
R_C(u,v)=\max\!\left\{\|C^\top(u-v)\|_1,\,
                         \|C^\top(u+v)\|_1\right\}.
\]
The proposed scale-transfer lemma asked for some \(S,u,v\) with
\[
R_C(u,v)\ge 2\bigl(1-(m/N)^{3/2}\bigr)M(A)-o(N^{3/2}).
\tag{1}
\]

For each symmetric Paley conference matrix below, I exhaustively enumerated:

1. every \(m\)-element principal subset \(S\);
2. every spin vector on \(S\) (fixing one spin by global sign symmetry);
3. all positive and negative extremizers of \(B=A[S]\);
4. every extremizer pair \((u,v)\), retaining the largest \(R_C(u,v)\).

The implementation is `rectangular_separation_enum.cpp`.

## Exact results

| conference order \(N\) | \(m\) | \(M(A)\) | required RHS of (1) | maximum \(R_C\) over all \(S,u,v\) | conclusion |
|---:|---:|---:|---:|---:|:---|
| 6 | 3 | 5 | 6.46447 | 6 | fails |
| 14 | 7 | 21 | 27.1508 | 26 | fails |
| 18 | 9 | 33 | 42.6655 | 34 | fails |

The full order-14 and order-18 parents are exactly centered:
\[
P(A)=Q(A)=21,\qquad P(A)=Q(A)=33,
\]
respectively.  At order 14 all 3432 half-size children satisfy
\(|P(B)-Q(B)|\le2\).  At order 18, 44676 of the 48620 half-size children
are exactly centered, and the best separation among those centered children is
still only \(34\).  Thus midpoint balance does not repair (1).

At \(N=18\), the weakest optimized subset has \(R_C=18\), and no subset reaches
the required \(42.6655\).  The best possible separation \(34\) is not a parity
rounding of the target; the deficit is \(8.6655\).

## Constant-level obstruction

Let \(I=\{i:u_i=v_i\}\) and \(J=\{i:u_i\ne v_i\}\).  These sets partition \(S\),
and exactly
\[
R_C(u,v)=2\max\!\left\{
 \|C_{I,S^c}^{\top}u_I\|_1,\,
 \|C_{J,S^c}^{\top}u_J\|_1
\right\}.
\tag{2}
\]
For a pseudorandom cross block and balanced \(|I|\sim|J|\sim m/2\), the
natural scale in (2) is
\[
R_C\sim \frac{2}{\sqrt\pi}(1-\alpha)\sqrt{\alpha}\,N^{3/2},
\qquad \alpha=m/N.
\tag{3}
\]
If \(M(A)\sim cN^{3/2}\), the desired scale is
\[
2c(1-\alpha^{3/2})N^{3/2}.
\tag{4}
\]
As \(\alpha\uparrow1\), (3) has first-order coefficient \(2/\sqrt\pi
\approx1.12838\), whereas (4) has coefficient \(3c\).  For the relevant
conference/ROM range \(c\approx0.4841\) to \(0.5\), this is \(1.452\) to \(1.5\).
Thus the exact finite failures agree with a genuine constant mismatch for
balanced positive/negative extrema.

## Verdict

The proposed rectangular-separation inequality is false even after:

- optimizing the principal subset;
- optimizing over every pair of child extrema;
- imposing exact midpoint balance on the parent; and
- restricting to exactly centered children.

The coupled max-extension inequality remains correct, but convergence cannot be
obtained by demanding the whole scale deficit from this single
\(\ell_1\)-separation statistic.  Any surviving variant must use more than one
positive/negative state (an energy-layer or multistate transport quantity), or
combine \(R_C\) with a separate quantitative contribution from \(D\).

