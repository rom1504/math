# Amplified exchangeable-block regularization of actual signings

2026-09-17. Director derivation. **Verified by independent reconstruction
by the Bernoulli and localization researchers.**
This combines exact cloned sign gadgets with overlapping-row information;
there is no Gaussian replacement error and no separately paid random-child
cap. It strengthens the controlled near-energy window, not convergence.

## 1. Quantitative theorem

Fix 0<gamma<2/3. For every sufficiently large N and every actual hollow
full signing A of order N there is another such W of the SAME order with

```math
 Q(W)\le Q(A)+C N^{3/2-\gamma/2},
 \log|\{z:Q(W)-|H_W(z)|\le N^{3/2-\gamma}\}|
       \le C_\gamma N^{1-\gamma/4}\log N.             \tag{1}
```

Only edges incident to O(N^(1/2-gamma/4)) vertices are rewritten.
The construction and constants are independent of the unknown optimum.
In particular this applies to actual exact minimizers, although the output
need not be an exact minimizer. Both polarities and EVERY physical spin
are included. No spectral-norm condition is used.

## 2. Exact physical model and compressed optimizer

Keep any principal m-block B. Append p exchangeable blocks, each of r
vertices. For each old vertex i and block a, choose one fair sign h_ia
and use it on ALL r edges. For each distinct block pair a,b choose one
fair sign e_ab and use it on ALL r^2 edges. Make every within-block edge
positive. All chosen random signs are independent. Every edge of the
resulting matrix is exactly plus or minus1.

For old spins x and block magnetizations t_a=r^(-1)sum_(j in a)z_j,

```math
 H_{W_p}(z)=H_B(x)+r\sum_a t_a h_a\cdot x
       +r^2\sum_{a<b}e_{ab}t_at_b
       +\frac12\sum_a(r^2t_a^2-r).                  \tag{2}
```

Each t_a lies in the finite set{-1,-1+2/r,...,1}; every such choice is
physically realizable. Choose an absolute maximizing compressed state
U=(x,t_1,...,t_(p+1)) uniformly among compressed maximizers. This tie
rule is equivariant under block permutations, with

```math
 H(U)\le m\log2+(p+1)\log(r+1).                     \tag{3}
```

Uniformity is over compressed states, not over all physical multiplicities.
Either rule is allowable, but(3) uses the compressed selection explicitly.

## 3. Read-two deletion estimate with amplified response

Let P_p=E Q(W_p), with B and r fixed. The random incident row R_a of
block a consists of its m independent old-edge signs and its p independent
other-block signs. Every underlying random sign is read at most twice.
The elementary overlapping-coordinate information inequality therefore
gives, by exchangeability,

```math
 I(U;R_a)\le 2[m\log2+(p+1)\log(r+1)]/(p+1).        \tag{4}
```

For fixed u, the row field h_a dot x+r sum_(b!=a)e_ab t_b is a
Rademacher sum with square coefficient sum at most V=m+pr^2. Its absolute
mean is at most sqrt(V) and its centered MGF has proxy V. The entropy
selection inequality consequently bounds its mean at the selected U by
sqrt(V)+sqrt(2V I(U;R_a)). Deleting block a at that optimizer loses at
most r times this absolute field plus r^2/2 from its internal clique
(the same bound holds for r=1). The deleted matrix has exactly the W_p law. Thus

```math
 0\le P_{p+1}-P_p
 \le r\sqrt{m+pr^2}
 \left[1+2\sqrt{\frac{m\log2}{p+1}+\log(r+1)}\right]
       +r^2/2.                                     \tag{5}
```

Nonnegativity follows by averaging deleted spins and convexity of absolute
value. This estimate includes all interblock and intrablock interactions.

For a fixed W_p let Delta_r(W_p)=E[Q(W_(p+1))|W_p]-Q(W_p)>=0.
Then E Delta_r=P_(p+1)-P_p. Write c_r=r(r-1)/2. For a fixed old
physical state z, restricting the new block to its two aligned words gives

```math
 \max_{\text{new spins}}|H_{W_{p+1}}|
 \ge |H_{W_p}(z)+c_r|
       +r|h\cdot x+r\sum_a e_a t_a|
 \ge |H_{W_p}(z)|-c_r
       +r|h\cdot x+r\sum_a e_a t_a|.                \tag{6}
```

For E(T)={z:Q(W_p)-|H_W(z)|<=T}, maximize(6) over E(T).
Averaging the new interblock signs by Jensen leaves at least
max_(x in projection E(T))|h dot x|. The projection is antipodal, so
its Bernoulli width obeys SIMULTANEOUSLY for every T>=0

```math
 b(\operatorname{proj}_{old}E(T))
       \le [T+c_r+\Delta_r(W_p)]/r.                 \tag{7}
```

The gain is the DIVISION BY r of the energy tolerance, with a completely
physical sign gadget. The unamplified vertex theorem had only T+Delta.

## 4. Jointly selecting the cap and every nearcode

Assume pr^2<=N=m+pr and m>=N/2. For any compressed state the random part
of(2) has variance proxy at most

```math
 mp r^2+\binom p2r^4\le(3/2)Npr^2.
```

There are at most2^N physical states. A union bound for both polarities
puts the random cap below C N r sqrt(p) with probability at least7/8.
The deterministic within-block term is at most pr^2/2 in absolute value,
absorbed by the same bound. Markov and(5) give a simultaneous event of
probability at least3/4 on which

```math
 Q(W_p)\le Q(B)+C N r\sqrt p,
 \Delta_r(W_p)/r\le
 C\sqrt{N+pr^2}\left[1+\sqrt{N/p+\log(r+1)}\right]+Cr.
                                                               \tag{8}
```

Only the fresh-block expectation is used; no need to sample that extra
block for the final signing. The SAME realization W_p satisfies(7) at
all tolerances. The code entropy is at most its projected entropy plus
pr log2. VC dimension of any Boolean code is at most its Bernoulli width,
so Sauer's bound gives

```math
 \log|E(T)|\le B_T\log(em/B_T)+pr\log2,
 \quad B_T=\min\{m,(T+c_r+\Delta_r)/r\}.             \tag{9}
```

The convention at zero is zero. This counts even all physically distinct
within-block words with the same magnetizations; none are omitted.

## 5. Exponents, exact order and what remains

Take p=floor(N^(gamma/2)), r=floor(N^(1/2-3gamma/4)), m=N-pr,
and keep that principal block of A. Then for all sufficiently large N,

```math
 pr=\Theta(N^{1/2-\gamma/4}),\quad
 pr^2=\Theta(N^{1-\gamma}),\quad
 N r\sqrt p=\Theta(N^{3/2-\gamma/2}),
 \Delta_r/r=O(N^{1-\gamma/4}),\quad
 N^{3/2-\gamma}/r=\Theta(N^{1-\gamma/4}).            \tag{10}
```

Since Q(B)<=Q(A), (8)--(10) prove(1). Floors cause only universal
constant changes once p,r tend to infinity. The final order is exactly N.

This is a strengthened unconditional regularization theorem: one may
choose a window arbitrarily close to N^(3/2) on a power scale, with
vanishing normalized preparation cost and subexponential cardinality.
It does NOT control a fixed normalized window eta N^(3/2), nor does it
prove favorable bridge value. For a comparable new child, a generic
subGaussian old-code selection error scales as
N sqrt(log|E|)=O(N^(3/2-gamma/8)sqrt(log N)), which exceeds both the
controlled window and the preparation cost. Consequently old-word escape
and a summable cross-order recurrence do not follow from(1).

Classical ingredients are read-k information, selected-bias bounds,
bounded differences, and Sauer's lemma. The exact complete-sign
amplification and quantitative tradeoff are the construction proved here.
External priority has not been established.

## 6. ANY prescribed shrinking normalized window

Let T_N>=0 be ANY prescribed sequence with T_N=o(N^(3/2)). Set

```math
 \bar\eta_N=\max\{T_N/N^{3/2},N^{-2/3}\},\qquad
 p=\lfloor\bar\eta_N^{-1/2}\rfloor,\quad
 r=\lfloor\sqrt N\bar\eta_N^{3/4}\rfloor.
```

For all sufficiently large N, r>=1, pr^2<=N, and pr=o(N). Applying the
same construction gives a same-order full signing W with

```math
 Q(W)\le Q(A)+C\sqrt{\bar\eta_N}\,N^{3/2},\qquad
 \log|E_W(T_N)|
 \le C N\bar\eta_N^{1/4}\log(e/\bar\eta_N)=o(N).      \tag{11}
```

Indeed r is at least half its unfloored value, as that value is >=1;
p is at least half its value once eta tends to zero. Consequently
T_N/r<=2N eta^(1/4), N/sqrt(p)<=sqrt(2)N eta^(1/4),
pr<=sqrt(N)eta^(1/4), and Nr sqrt(p)<=sqrt(eta)N^(3/2).
The extra logarithm in(8) is uniformly dominated because N/p>=N^(2/3).
Equations(7)--(9) prove(11), including T_N=0.

The choice of W_N DEPENDS on the prescribed window sequence. This is
not a single family controlling every shrinking window simultaneously.
It gives an unconditional actual near-minimizer preparation at every
order and every prescribed vanishing normalized tolerance.

## 7. Paid preparation versus parent escape

The general parameters make the tradeoff explicit. Write
d=r sqrt(p)/sqrt(N) for the certified normalized preparation scale.
At T=eta N^(3/2), (7)--(8) supply normalized core width of order

```math
 \eta\sqrt p/d+1/\sqrt p
```

up to the displayed smaller terms. Optimizing this RHS gives scale
sqrt(eta/d). Thus if a proposed extension needs d=o(epsilon) and must
control its deficit band eta~epsilon, THIS estimate is trivial. This
is a limitation of the proved certificate, not a lower bound on actual
width and not an impossibility theorem for other preparations.
The stronger prescribed-window theorem removes a genuine regularization
obligation, but does not supply a free favorable-parent estimate.

## 8. The offset mechanism is not specific to a quadratic core

Sections2--4 remain valid if H_B(x) is ANY deterministic real function on
the old Boolean cube. Define the new landscape by the right side of(2),
leaving H_B unchanged. Its maximum absolute value increases by at most
C N r sqrt(p), and the same nearcode width and entropy bounds hold.
Indeed the proof uses the old function only as a deterministic offset;
exchangeable selection, deletion, and the fresh aligned-block identity
are unchanged. If the old function is not even, its projected code need
not be antipodal, but b(C)<=E max_(x in C)|h dot x| still proves(7).

The SAME-order corollary(1) uses principal restriction monotonicity of
the original quadratic signing, so that particular formulation is not
claimed for an arbitrary offset without an appropriate restriction rule.
For arbitrary offsets this is a near-order pairwise-sign-gadget extension,
not a statement that an arbitrary function becomes a pure quadratic form.
It explains the mechanism's scope: physical weak regularization of an
arbitrary optimization landscape, with exact complete sign support on
every NEW edge, rather than an unproved rigidity property of old optima.

## 9. Independent verification and reproducibility

Both the Bernoulli and localization researchers reconstructed the complete
proof, including the compressed tie rule and physical multiplicities.
The director's exact finite replay is
`computations/paper_director_2026_09_17_cloned_regularization.py`.
It enumerates every random block signing and every physical spin in six
small parameter settings, recomputes the conditional fresh-block increment,
and checks(7) using integers at every selected tolerance. All checks pass.
The output is preserved under the campaign's dated research archive.
This finite replay checks the identities, not the asymptotic concentration
or entropy proof, which is given above.
