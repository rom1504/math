# Action-limit revision and moving-stabilizer experiment

Date: 2026-08-16.  I read the restriction/action verification, coding-toolkit
§2.5, and the requested \(L_{\rm drift}\) section only.

## 1. Retain/withdraw judgment on spectral regularization

Retain fixed-\(C\) spectral regularization as a standalone class-B structural
lemma, but **withdraw it as a convergence architecture**:

\[
\boxed{L_{\rm reg}}\qquad
\exists C<\infty,\ \epsilon_n\downarrow0\quad\forall n\gg1\quad
\exists A_n\in\mathfrak A_n:
\quad Q(A_n)\le M_n+\epsilon_n n^{3/2},\quad
\|A_n\|_{\rm op}\le C\sqrt n. \tag{SR}
\]

This is strictly stronger than the archived two-limit tradeoff
\[
Q(A'_n)\le M_n+O(K^{-1/2}n^{3/2}),\qquad
\|A'_n\|_{\rm op}=O(K\sqrt n),
\]
because fixing the operator bound there leaves a fixed leading error.  It
does not recover a Boolean maximizer or a coset histogram.

It nevertheless cannot be strengthened into a strict convergence reduction
using action compactness alone.  Fix any witnesses in (SR), and let
\(\mathscr M_C\) be their action cluster set.  Compactness and the proved
continuity of \(\Phi\) give the exact identity
\[
\left\{\frac{\Phi(T)}2:T\in\mathscr M_C\right\}
=\operatorname{Clust}\left\{\frac{M_n}{n^{3/2}}\right\}. \tag{1}
\]
Indeed every scalar cluster subsequence has an action-convergent
subsubsequence, and every action cluster has the corresponding scalar limit.
Thus “all minimizing action limits have the same \(\Phi\)” is equivalent to
the desired convergence.  Turning one cluster limit into competitors at all
other orders is precisely the archived AR obligation.  SR removes neither
alternative.

**Decision:** keep (SR) as an independently worthwhile regularization target;
delete \(L_{\rm FA}=\mathrm{SR}+\mathrm{AR}\) from the candidate convergence
list.  No revised sufficient action lemma survives the collision audit.

The exact falsifier for (SR) is:
\[
\forall C<\infty\ \exists\delta_C>0\ \exists^\infty n:\quad
\min_{\substack{A\in\mathfrak A_n\\\|A\|_{\rm op}\le C\sqrt n}}Q(A)
\ge M_n+\delta_Cn^{3/2}. \tag{2}
\]

## 2. Foreign mechanism tested: moving coset harmonics

The tested import is the moving-projection mechanism of Chapter 2 of
[*Ten Advances in Mathematics and Theoretical Computer Science*](https://cdn.openai.com/pdf/ten-proofs-oai.pdf).
For packing, a nontrivial Boolean harmonic module \(E_k\) is embedded across
Fourier levels \(k,\ldots,L\) at each moving codeword.  Multiplicity-free
branching turns overlaps into scalar two-point kernels, and a block-Jacobi
Perron vector gives the certificate.

The covering analogue would root at every ambient coset
\(U\in\mathbb F_2^{E_n}/\mathcal C_n^+\), translate \(U\) to zero, attach a
module \(E_{k,U}\), and localize the nearest-code shell.  One would hope that
growing \(k,L\) retain the terminal extreme while a small transition state
can be moved from order \(n\) to every order \(m\).

### Exact finite experiment

I enumerated every augmented-cut coset, computed
\[
r(U)=d(U,\mathcal C_n^+),\qquad
b(U)=|\{e:r(U+e)=r(U)+1\}|,
\]
and, at \(n=6\), enumerated the exact \(S_6\)-orbits on cosets.

| order | exact result |
|---|---|
| \(n=6\), \(E=15\), 512 cosets | at \(r=2\), \(b=10\) or \(13\); at \(r=3\), \(b=0\) or \(4\); 25 nondeep dead ends occur at \(r=3<\rho=5\) |
| \(n=6\), orbit check | each of \((r,b)=(3,0)\) and \((3,4)\) contains two distinct \(S_6\)-orbits |
| \(n=7\), \(E=21\), 16,384 cosets | at \(r=4\), \(b\in\{6,9,12,14,17\}\); at \(r=5\), \(b\in\{0,2,4,5,10\}\); 1,260 nondeep dead ends occur at \(r=5<\rho=6\) |

This disproves radial/terminal lumpability at these orders, not the asymptotic
\(L_{\rm drift}\): the \(n=7\) nondeep dead ends have
\(z=(21-10)/7^{3/2}=0.594\ldots\), outside the mandatory subinterval
\([0.33,0.51]\).  The candidate may choose \(I\) to exclude this value.

### Why the moving module does not yield all-order exact signs

The packing theorem moves among **codewords** under the full Hamming
symmetry.  An arbitrary covering root instead sees the translated coset
\(U+\mathcal C_n^+\).  Its nearest shell is the function
\[
c\longmapsto d(U,c),
\]
equivalently the complete Boolean energy table of the signing represented by
\(U\).  The coordinate automorphism group of the cut code is only \(S_n\);
generic cosets have tiny stabilizer, and there are at least
\(2^{E_n-n}/n!\) coset orbits.  The scalar strong-Gelfand reduction used for
packing therefore disappears at arbitrary roots.

Positive-degree localizers can distinguish the finite scatter only by
retaining root-dependent orbit variables.  At growing degree
\(\Theta(E_n)\), the \(S_n\)-orbits are graph types and the state is
exponential.  Full finite convergence of a moment/Lasserre hierarchy restores
integrality only at the level that determines an actual root-supported
measure; here that recovers the coset landscape the proposal was meant to
avoid.  Lower levels provide PSD pseudo-moments and inequalities, not a
symmetric hollow sign matrix of each prescribed order.

This matches the source boundaries:

- Schrijver and the moving-projection theorem are packing-only.
- Gijswijt--Polak add genuine all-root localizers, but optimize the size of an
  unrestricted covering code at one fixed length and radius.
- Riener--Rolfes--Vallentin give finite-space hierarchy convergence, while
  their first symmetric level collapses to the volume bound.

None supplies a cross-order flat-extension theorem for the prescribed cut
code.  The moving-stabilizer experiment therefore gives **no AR theorem
without full coset data**.

The strongest exact falsifier for any proposed moving-stabilizer AR is an
order-confined limit:
\[
\exists T\quad\forall C'<\infty\quad\exists\varepsilon_{C'}>0\
\exists^\infty m:
\quad
\inf_{\substack{B\in\mathfrak A_m\\\|B\|_{\rm op}\le C'\sqrt m}}
d_M(T_B,T)\ge\varepsilon_{C'}. \tag{3}
\]

## 3. Cross-domain critic: terminal drift

The strongest failure reason for \(L_{\rm drift}\) is that \(b_n(U)/E_n\) is
a microscopic, discontinuous derivative of distance at an exponentially rare
root.  Neighboring normalized heights differ by only
\(2/n^{3/2}\), so an action limit of the normalized coset random walk together
with the mark \(z_n\) collapses inward, flat, and outward moves to the same
limiting height.  Moreover, uniformity over every terminal coset silently
lumps exponentially many order-dependent \(S_n\)-orbits into one
\(\beta(z)\).  Its unique-zero condition excludes all terminal nondeep traps;
the finite experiment shows such traps and nonlumpability already occur.

The strongest graph-limit rescue is Backhausz--Szegedy compactness for Markov
graphops and simultaneous action convergence.  Introduce the normalized coset
walk \(P_n\) and the outward sub-Markov operator
\[
(P_n^+f)(U)=\frac1{E_n}\sum_{e:\,r(U+e)=r(U)+1}f(U+e);
\qquad P_n^+\mathbf1(U)=b_n(U)/E_n. \tag{4}
\]
A simultaneous **pointed** action limit of \((P_n,P_n^+,z_n)\) would retain
the terminal statistic, unlike \((P_n,z_n)\) alone.

This does not presently rescue the proposal.  Constructing \(P_n^+\) already
uses exact distances to the code; unpointed compactness loses exponentially
rare deepest roots; compactness is subsequential; and no theorem realizes a
chosen pointed extreme at every large order.  Proving a common law
\(P^+\mathbf1=\beta(z)\) uniformly on terminal roots is the unproved
lumpability, not a consequence of action convergence.

## 4. Final verdict and confidence

- Revised exact lemma: retain only \(L_{\rm reg}\) in (SR), with no convergence
  claim.
- Convergence architecture: **negative verdict**; AR remains unchanged, and
  the moving-stabilizer import either stays nonconstructive or reaches exact
  signs by encoding full coset data.
- \(L_{\rm drift}\): no realizable all-order limit object is supplied by known
  graph-limit theorems; uniform lumpability bears the extreme-coset burden.

Confidence: (SR) truth/tractability \(0.55/0.25\); moving-stabilizer
all-order recovery without full coset data \(0.05/0.03\); negative
architecture verdict \(0.95\); the conditional implication of (SR)+AR remains
\(0.99\).
