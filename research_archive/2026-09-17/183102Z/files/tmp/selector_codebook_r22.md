# Wave 22 selector codebooks and exact one-deletion incidence programs

## 1. An operational sufficient lemma

Let (S\sim U_m), let (c_A(S,d)=\langle A[S],d[S]\rangle), and let
(q_S=Q(A[S])).  Suppose a finite codebook
(\mathcal C\subseteq\mathcal D_n) satisfies

\[
 \mathbb E_{S\sim U_m}
 \left[q_S-\max_{d\in\mathcal C}c_A(S,d)\right]\leq\delta.
\]

Choose, with fixed tie breaking, a maximizing codeword (D=D(S)).  Then
(D) has the required average distortion and

\[
 I(S;D)\leq H(D)\leq\log |\mathcal C|.
\]

Consequently (R_{U_m}(\delta)\leq\log|\mathcal C|).  In the notation of
the fixed-slice log-mgf theorem (10.690)--(10.704), it is therefore enough
at fixed density to construct codebooks with

\[
 \log|\mathcal C|=O(n^{1/2-2c}),\qquad
 \delta=O(n^{3/2-c}),\qquad 0<c<1/4.
\]

This is a purely combinatorial shared-prior target: find a subexponential
family of full cuts whose best restrictions approximate the optimized
principal sections on average.

## 2. Exact one-deletion incidence program

For (m=n-1), associate to an oriented full cut (d) its coverage pattern

\[
 P(d)=\{i:d[-i]\text{ is an exact ground state of }A[-i]\}.
\]

For a distribution (\nu) on full cuts put
(z_i=\nu\{d:i\in P(d)\}).  Exact target landing has rate

\[
 R_{U_{n-1}}(0)=\min_\nu -\frac1n\sum_{i=1}^n\log z_i.
\]

The finite convex program has a useful exact optimality certificate.  If
(\pi_i=1/n), a candidate (\nu) is optimal when

\[
 h(P):=\sum_{i\in P}\frac{\pi_i}{z_i}\leq1
\]

for every realizable pattern (P), with equality on every pattern in the
support of (\nu).  This is just the simplex KKT condition, and it can be
checked with rational arithmetic after enumerating all oriented cuts.

The verifier gives the following exact certificates for the three stored
minimizers.

* **A6.**  Put mass (1/6) on each of the six parent grounds whose pattern
  omits one coordinate.  Then every (z_i=5/6), all 26 realizable nonempty
  patterns obey (h(P)\leq1), and
  (R_{U_5}(0)=\log(6/5)=0.1823215568\) nats.
* **A8.**  Put mass (1/8) on each of the eight size-three parent-ground
  patterns listed in the verifier.  Then every (z_i=3/8), all 20 patterns
  obey the KKT inequality, and
  (R_{U_7}(0)=\log(8/3)=0.9808292530\) nats.
* **A9.**  A rational optimum supported on seven parent-ground patterns is

  \[
  \begin{array}{c|ccccccc}
  P&(3,6,8)&(3,5,7)&(0,4,8)&(1,3,8)&(2,3,5)&(1,2,7)&(0,6,7)\\
  \nu(P)&1/10&1/20&1/4&1/20&1/5&1/5&3/20.
  \end{array}
  \]

  It yields
  (z=(2/5,1/4,2/5,2/5,1/4,1/4,1/4,2/5,2/5)).  All 22 realizable
  patterns satisfy (h(P)\leq1), equality holds on the support, and

  \[
  R_{U_8}(0)
   =-\frac19\left(5\log\frac25+4\log\frac14\right)
   =1.1251812338\text{ nats}.
  \]

Uniform parent grounds are not optimal for A9 (they give about 1.4475
nats), so the incidence optimization is genuinely useful.

## 3. Scope

These computations show that exact optimized child landing can require
only constant information in nontrivial one-deletion examples, and that a
small mixture of parent grounds can outperform the uniform ground prior.
They do **not** give the fixed-density asymptotic codebook required by the
operational lemma.  The next mathematical question is whether minimizer
structure forces a family of at most
(\exp(O(n^{1/2-2c}))) full cuts to cover most optimized (m)-sections up
to total deficit (O(n^{3/2-c})), or whether a minimizer family has an
exponentially rich collection of mutually incompatible section grounds.

Run `./.venv/bin/python tmp/selector_codebook_r22.py` from the repository
root for the exact enumeration and rational KKT checks.
