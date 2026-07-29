# Nested restriction / paving route

## Goal

For a symmetric zero-diagonal signing \(A\) of order \(N\), write

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
P(A)=\max_x H_A(x),\qquad
Q(A)=-\min_x H_A(x),
\]
\[
M(A)=\max(P(A),Q(A)),\qquad
W(A)=\frac{P(A)+Q(A)}2.
\]

The desired scale-preserving restriction statement is

\[
M(A[S])\le
\left(\frac{|S|}{N}\right)^{3/2}M(A)+o(N^{3/2})
\tag{R}
\]

for some \(S\) of each prescribed proportional size.  If (R) held uniformly for
every fixed \(\alpha=|S|/N\), applying it to a liminf subsequence of near-optimal
signings would transfer the liminf value to all intermediate orders and prove
convergence of \(M_n/n^{3/2}\).

The exact finite form of (R) is false, already for Paley conference matrices.
The asymptotic \(o(N^{3/2})\) form is not decided by these examples.

## 1. Exact partition inequality for the centered range

Let \(S_1,\ldots,S_k\) be any partition of the vertices. Then

\[
\boxed{\sum_{r=1}^k P(A[S_r])\le P(A),\qquad
\sum_{r=1}^k Q(A[S_r])\le Q(A).}
\tag{1}
\]

Indeed, choose on every block a positive ground state \(x^{(r)}\), multiply
each block independently by a sign \(\sigma_r\), and join the block vectors.
The internal energy is always \(\sum_rP(A[S_r])\), while every cross-block
edge has expectation zero over the \(\sigma_r\)'s. Some choice of the block
signs therefore has total energy at least \(\sum_rP(A[S_r])\). The negative
statement is identical, using negative ground states.

Consequently

\[
\boxed{\sum_{r=1}^k W(A[S_r])\le W(A).}
\tag{2}
\]

For equal blocks, at least one block obeys \(W(A[S_r])\le W(A)/k\). Since
\(W\le M\le2W\), this gives only

\[
M(A[S_r])\le\frac{2M(A)}k.
\tag{3}
\]

For \(|S_r|=N/k\), (3) has coefficient \(\alpha\), whereas the required one is
\(\alpha^{3/2}\). Thus the exact block-randomization argument loses a factor
\(\alpha^{-1/2}\).

The exponent in (2) cannot be improved from scalar partition information
alone: the proof uses only the ability to randomize one global sign per block,
and it is exact whenever cross-block ground-state couplings can be made zero.

More formally, iterate (2) through an arbitrary rooted tree of equitable
partitions. At the leaves one still has

\[
\sum_{\text{leaves }L}W(A[L])\le W(A);
\]

if every leaf has relative size \(\alpha\), this says only
\(\min_LW(A[L])\le\alpha W(A)\), independent of the number of levels.
The abstract set function \(w(S)=c|S|\) saturates every inequality in this
calculus, so no nonnegative combination of partition inequalities can produce
the extra factor \(\sqrt\alpha\).

Vertex switching is not a way around this loss. For every diagonal sign matrix
\(D\),

\[
M(DAD)=M(A),\qquad W(DAD)=W(A).
\tag{3a}
\]

The same is true after restriction. Random block signs are therefore a proof
device for (1), not a construction that can lower the norm of a block or of a
union of blocks.

## 2. Exact random-restriction moments

Let \(S\) be a uniformly random \(m\)-subset of \([N]\), and set

\[
p_j=\frac{(m)_j}{(N)_j}.
\]

For a fixed \(x\in\{\pm1\}^N\), put \(h=H_A(x)\) and define its switched local
fields

\[
r_i=x_i(Ax)_i=\sum_{j\ne i}a_{ij}x_ix_j.
\]

Direct classification of pairs of edges as equal, adjacent, or disjoint gives

\[
\boxed{\mathbb E_S H_{A[S]}(x_S)=p_2h,}
\tag{4}
\]

and

\[
\boxed{
\begin{aligned}
\mathbb E_S H_{A[S]}(x_S)^2
={}&p_4h^2+(p_2-p_4)\binom N2\\
&+(p_3-p_4)
\left(\sum_{i=1}^N r_i^2-N(N-1)\right).
\end{aligned}}
\tag{5}
\]

Formula (5) exposes the obstruction to a direct concentration/union-bound
proof: the fluctuation depends on the full local-field square
\(\sum_i r_i^2\), not only on \(M(A)\). Even the bootstrap
\(\|A\|_{\mathrm{op}}=O(N^{5/6})\) for a near-minimizer permits this quantity
to be as large as \(N\|A\|_{\mathrm{op}}^2=O(N^{8/3})\). For a conference
matrix it is \(O(N^2)\), but taking a maximum over exponentially many spin
vectors restores fluctuations on the leading \(N^{3/2}\) scale.

## 3. Exact balanced-cut identity

Encode \(S\) by \(s_i=1\) on \(S\) and \(s_i=-1\) off \(S\). For every full
spin vector \(x\),

\[
\boxed{
4H_{A[S]}(x_S)
=H_A(x)+H_A(xs)+\sum_{i=1}^N s_i x_i(Ax)_i.}
\tag{6}
\]

The first two terms are bounded by \(M(A)\), but the last term is a signed
local-field discrepancy. A triangle inequality therefore again gives a
linear-size restriction estimate, not the desired \(3/2\)-homogeneous one.
Any proof through (6) must choose \(s\) so as to control the local-field term
simultaneously for every \(x\) whose restriction is a ground state.

## 4. Exact conference stress tests

For \(q\equiv1\pmod4\), the symmetric Paley conference matrix of order
\(N=q+1\) is

\[
C_{\infty,j}=C_{j,\infty}=1,\qquad
C_{ij}=\chi(i-j)\quad(i,j\in\mathbb F_q,\ i\ne j),
\]

and satisfies \(C^2=qI\).

Exhaustive enumeration over all principal subsets gives:

### Order \(6\) (\(q=5\))

\[
M(C_6)=5.
\]

| \(m\) | \(\min_{|S|=m}M(C_6[S])\) | \((m/6)^{3/2}M(C_6)\) | ratio |
|---:|---:|---:|---:|
| 2 | 1 | 0.962 | 1.039 |
| 3 | 3 | 1.768 | 1.697 |
| 4 | 4 | 2.722 | 1.470 |
| 5 | 4 | 3.804 | 1.052 |

### Order \(14\) (\(q=13\))

\[
M(C_{14})=21.
\]

| \(m\) | minimum | maximum | target | min/target |
|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1.134 | 0.882 |
| 3 | 3 | 3 | 2.083 | 1.440 |
| 4 | 4 | 6 | 3.207 | 1.247 |
| 5 | 4 | 8 | 4.482 | 0.892 |
| 6 | 5 | 9 | 5.892 | 0.849 |
| 7 | 9 | 11 | 7.425 | 1.212 |
| 8 | 10 | 14 | 9.071 | 1.102 |
| 9 | 12 | 14 | 10.824 | 1.109 |
| 10 | 15 | 17 | 12.677 | 1.183 |
| 11 | 19 | 19 | 14.626 | 1.299 |
| 12 | 20 | 20 | 16.665 | 1.200 |
| 13 | 20 | 20 | 18.791 | 1.064 |

In particular,

\[
\min_{|S|=7}M(C_{14}[S])=9>
\left(\frac12\right)^{3/2}21=7.4246\ldots .
\tag{7}
\]

Thus the coefficient in (R), with no additive error, is false even for exact
conference matrices. These finite examples do not refute an error
\(o(N^{3/2})\).

### Order \(18\) (\(q=17\))

An independent exhaustive computation gives

\[
M(C_{18})=P(C_{18})=Q(C_{18})=W(C_{18})=33.
\]

| \(m\) | minimum | maximum | target | min/target |
|---:|---:|---:|---:|---:|
| 5 | 4 | 8 | 4.831 | 0.828 |
| 6 | 5 | 11 | 6.351 | 0.787 |
| 7 | 9 | 13 | 8.003 | 1.125 |
| 8 | 10 | 14 | 9.778 | 1.023 |
| 9 | 12 | 18 | 11.667 | 1.029 |
| 10 | 15 | 19 | 13.665 | 1.098 |
| 11 | 19 | 21 | 15.765 | 1.205 |
| 12 | 20 | 24 | 17.963 | 1.113 |
| 13 | 22 | 24 | 20.254 | 1.086 |
| 14 | 25 | 27 | 22.636 | 1.104 |
| 15 | 29 | 29 | 25.104 | 1.155 |
| 16 | 32 | 32 | 27.656 | 1.157 |
| 17 | 32 | 32 | 30.289 | 1.057 |

At the half-restriction, all three tested conference orders miss the literal
target by exactly the unavoidable energy lattice:

| \(N\) | target | minimum \(M\) | minimum \(W\) |
|---:|---:|---:|---:|
| 6 | 1.768 | 3 (next admissible odd value) | 2 (next integer) |
| 14 | 7.425 | 9 (next admissible odd value) | 8 (next integer) |
| 18 | 11.667 | 12 (next admissible even value) | 12 (next integer) |

The normalized excess

\[
\frac{\min_{|S|=N/2}M(C_N[S])-2^{-3/2}M(C_N)}{N^{3/2}}
\]

is \(0.0838,0.0301,0.00436\) at \(N=6,14,18\), respectively. This is not an
asymptotic argument, but the exact data give no leading-order conference
obstruction; at one-half they are consistent even with an \(O(1)\) correction.

## 5. Spectral paving cannot preserve a constant below \(1/2\)

A paving theorem could seek a block \(S\) with small
\(\|A[S]\|_{\mathrm{op}}\), followed by

\[
M(A[S])\le\frac{|S|}{2}\|A[S]\|_{\mathrm{op}}.
\tag{8}
\]

But every signing \(B\) of order \(m\) satisfies

\[
\|B\|_{\mathrm{op}}^2\ge\frac{\operatorname{tr}B^2}{m}=m-1.
\tag{9}
\]

Hence (8), even with an optimally flat paved block, has the unavoidable
asymptotic floor

\[
\frac12m^{3/2}.
\tag{10}
\]

It cannot preserve a hypothetical liminf constant \(c<1/2\). Kadison--Singer
or interlacing can recover the correct exponent for conference matrices, but
not the near-optimal constant required for a convergence proof.

## 6. An energy-layer-sensitive restriction lemma

There is a clean sufficient condition that uses precisely the information
discarded by \(M(A)\), \(W(A)\), and the spectrum.

For \(0\le t<L<M(A)=K\), define the two full-cube layers

\[
\mathcal L_t^+=\{x:H_A(x)\ge t\},\qquad
\mathcal L_t^-=\{x:H_A(x)\le-t\},
\]

and put

\[
\delta=\frac{L-t}{K-t}.
\]

Then the fraction of \(m\)-subsets \(S\) with \(P(A[S])\ge L\) is at most

\[
\frac{|\mathcal L_t^+|}{\delta\,2^{N-m}},
\tag{11}
\]

and the fraction with \(Q(A[S])\ge L\) is at most

\[
\frac{|\mathcal L_t^-|}{\delta\,2^{N-m}}.
\tag{12}
\]

Consequently,

\[
\boxed{
|\mathcal L_t^+|+|\mathcal L_t^-|
<\delta\,2^{N-m}
\quad\Longrightarrow\quad
\exists\,|S|=m:\ M(A[S])<L.}
\tag{13}
\]

**Proof.** If \(P(A[S])\ge L\), choose a witnessing \(y_S\). Extend \(y_S\)
by independent uniform spins on \(S^c\). Both the cross term and the internal
term on \(S^c\) have mean zero, so

\[
\mathbb E[H_A(y_S,X_{S^c})]=H_{A[S]}(y_S)\ge L.
\]

The random energy lies in \([-K,K]\). If \(p\) is the probability that it is
at least \(t\), then

\[
L\le pK+(1-p)t,
\]

and hence \(p\ge\delta\). Thus each positive-bad support has at least
\(\delta2^{N-m}\) extensions in \(\mathcal L_t^+\). Count pairs
\((S,x)\). A fixed \(x\) is paired with at most \(\binom Nm\) supports, which
proves (11). Apply the same argument to \(-H_A\) for (12), then use the union
bound over bad supports.

For fixed \(\alpha\in(0,1)\), (13) gives the following asymptotic conditional
form. If \(m=\alpha N+o(N)\) and for some \(\eta>0\)

\[
\limsup_{N\to\infty}\frac1N
\log_2\#\left\{x:
|H_A(x)|\ge(\alpha^{3/2}-\eta)M(A)\right\}
<1-\alpha,
\tag{14}
\]

then for all large \(N\) there is an \(m\)-subset with

\[
M(A[S])\le\alpha^{3/2}M(A).
\tag{15}
\]

This is a genuine energy-layer criterion, rather than a scalar norm
inequality. It is also demanding: the exact \(n=5,6\) examples with a positive
fraction of the entire cube at the ground energy show that no such entropy
hypothesis follows formally from small \(M(A)\). Removing (14), or proving an
alternative when (14) fails, is the natural dichotomy left by this route.

## 7. Current verdict

The literal restriction theorem is false at finite scale. Three natural
repairs currently stop at precise barriers:

1. centering by the energy range gives the exact hereditary inequality (2),
   but only exponent \(1\);
2. random restriction gives the exact moment formula (5), whose local-field
   term is not controlled strongly enough uniformly over the Boolean cube;
3. spectral paving has the intrinsic constant-\(1/2\) floor (10).

The surviving possibility is a nonlinear, ground-state-sensitive paving:
choose \(S\) by simultaneously balancing the local-field discrepancy in (6)
only on the relevant high-energy layers. This is strictly more information
than \(M(A)\), \(W(A)\), or the spectrum, and it is the same missing
energy-layer statistic encountered in proportional gluing and vertex
insertion. Lemma (13) proves the restriction theorem under the explicit
entropy condition (14); the high-layer alternative remains uncontrolled.
