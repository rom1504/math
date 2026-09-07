# Coding revision after archive check and spin-glass packet

## Decision

**Withdraw `L_drift` as a genuinely weaker execution route.** Its state `(r,b)` is strictly poorer than a
coset histogram, but its asymptotic assertion is not poorer than the target. An infinite Paley trap family
forces its unique zero to be `1/2`, so `L_drift` would prove the sharp universal lower constant, not merely
convergence to an unspecified constant. I retain only the bounded disproof-first test below. `L_star` and
`L_mom` remain Class C and are not revived.

## 1. Exact Paley audit

### The recorded order-10 trap is not itself conference

For the signing `A` with the 23 negative edges listed in the verification, direct multiplication gives

```math
(A^2)_{01}=-2.
```

A conference signing of order 10 has `A^2=9I`. Switching diagonally conjugates `A^2`, and global negation
leaves it unchanged. Hence this trap is not switching-equivalent to Paley conference, although its certified
data remain `Q(A)=15`, `Q(A^e)=17` for every edge, and `b(A)=0`.

### A separate infinite Paley conference trap family does exist

Let `p` be any odd prime, `q=p^2`, `n=q+1`, and let `C` be the symmetric Paley conference matrix on
`{infinity} union F_q`:

```math
C_{infinity,a}=C_{a,infinity}=1,\qquad C_{a,b}=chi_q(a-b),\qquad C_{a,a}=0.
```

Then `C^2=qI`. Write `F_q=F_p(omega)`, with `omega^2=d` nonsquare in `F_p`. Choose
`T subset F_p`, `|T|=(p-1)/2`, and define

```math
S=\bigcup_{t\in T}(t\omega+F_p),\qquad
x_{infinity}=1,\qquad x_a=-1\ (a\in S),\quad x_a=1\ (a\notin S).
```

This is an exact Boolean eigenvector:

```math
Cx=px.                                                         \tag{P1}
```

Indeed, `sum_{F_q}x_a=p`. Within one additive `F_p`-coset all nonzero differences are squares in `F_q`;
between two distinct cosets the quadratic-character sum is `-1`. For the finite Paley core `B`,

```math
(B1_S)_a=(p+1)/2\ (a\in S),\qquad (B1_S)_a=(1-p)/2\ (a\notin S),
```

so `1+B(1-2 1_S)=p(1-2 1_S)`. The spectral bound and (P1) give

```math
Q(C)=np/2.                                                     \tag{P2}
```

The Paley two-graph has `PSL(2,q)` as a 2-transitive switching group on the projective line: group elements
lift to signed permutations `T` with `T^TCT=C`. The vector in (P1) disagrees with `C` on some edge because
`np/2<N_n`; two-transitivity transports such a `+p` Boolean eigenvector to any prescribed edge `e`. Hence
`Q(C^e)>=Q(C)+2`. An edge flip changes every energy by at most `2`, proving

```math
Q(C^e)=Q(C)+2\quad(e\in E(K_n)),\qquad b(C)=0.                 \tag{P3}
```

The switching action is standard; see [de Launey--Stafford](https://doi.org/10.1016/j.disc.2007.07.118).
Finally,

```math
z(C)={Q(C)\over n^{3/2}}={1\over2}\sqrt{{q\over q+1}}\longrightarrow {1\over2}. \tag{P4}
```

For these roots `b/N=0`. Uniformity and continuity in `L_drift` imply `beta(1/2)=0`; uniqueness forces
`c=1/2`. Applied to deepest cosets, the lemma proves `M_n/n^{3/2}->1/2`. Thus the director's implication is
correct, although the listed order-10 matrix was not the conference representative.

## 2. Exact minimal theorem left by drift

Let `I` be compact and contain `[0.33,0.51]`. The minimal trap-only consequence of `L_drift` is:

> **`L_trap(1/2)`.** For every `epsilon>0` there is `n_0` such that, for every `n>=n_0` and every coset
> `U` with `z_n(U) in I` and `b_n(U)=0`,
>
> ```math
> |z_n(U)-1/2|<epsilon.                                        \tag{T}
> ```

Every deepest coset has `b=0`, and the frontier puts it in `I`; hence (T) proves the sharp limit `1/2`.
It asks more than that limit: it collapses **all** terminal dead ends, while the original problem controls
only the smallest defect. At order 10 the deepest root has `z=0.411096...` and a nondeep dead end has
`z=0.474341...`, displaying the extra finite-order content.

Therefore terminal lumpability is weaker only informationally: actual roots with equal `(r,b)` can have
different outer histograms. As an asymptotic obligation after (P4), it contains the sharp lower theorem plus
an all-traps assertion. It is not a strict reduction of the original optimization.

## 3. Disproof-first stopping test

Before any proof attempt, seek

```math
\exists\delta>0,\ n_j\to\infty,\ U_j:\quad b_{n_j}(U_j)=0,\quad z_{n_j}(U_j)\in I,
\quad z_{n_j}(U_j)\le1/2-\delta.                               \tag{STOP}
```

Paley already supplies a trap branch tending to `1/2`; (STOP) supplies a second branch and falsifies both
(T) and unique-zero `L_drift`. Use the exact certificate

```math
b(U)=N-\left|\bigcup_{g_v\le1}N_v\right|:
```

an upper certificate for `Q(U)` plus top-two switching states whose disagreement sets cover every edge.
The sole checkpoint is to amplify the order-10 nonconference trap, or another algebraic trap, to infinitely
many orders with a fixed normalized gap. If one substantive amplification attempt fails, return the route to
hold; do not start a general lumpability proof.

## 4. Foreign mechanism: constrained-overlap metastable interpolation

The foreign mechanism is Panchenko's [exact-overlap interpolation](https://arxiv.org/abs/math/0405359):
partition replica pairs into overlap cells, compare exact constraints with windows, and obtain
near-superadditivity with an `O(sqrt n)` size defect. Coding theory ordinarily tracks distance distributions,
not the constrained Gibbs complexity of near-ground switching states.

Here the tempting metastable cluster is

```math
V_1(U)=\{v\in C_n^+:g_v=(Q-a\cdot v)/2\le1\}.
```

The one-edge identity says `b(U)=0` iff the disagreement supports of `V_1(U)` cover every edge. One might
constrain overlaps in `V_1(U)` and interpolate its pressure, hoping that a zero-complexity boundary classifies
trap energies.

This is not a genuinely new route. Edge coverage is a union property not determined by pair overlap and can
be witnessed by only `O(n^2)` states, hence at zero Gibbs complexity. The layer `g<=1` is microscopic and
zero-temperature, whereas exact-to-window interpolation controls extensive pressure. Splitting `K_n` also
leaves `Theta(n_1n_2)` adversarial cross edges, with neither independent Gaussian blocks nor a signed convex
covariance remainder. The 2026 centered-Ising Hamilton--Jacobi theorem removes convexity only for quenched
Gaussian disorder; the 2026 LDP is the upper tail at speed `n`, while minimization needs the lower tail at
speed `n^2`. A deterministic coverage-to-overlap theorem would just be a new form of the missing trap lemma.

## 5. Cross-domain critique: fixed-temperature `L_Lap`

For a coset/root `U`, spin-glass pressure is exactly a coset-enumerator evaluation:

```math
P_{n,beta}(U)={1\over beta n}\log\sum_{c\in C_n^+}
e^{beta(N_n-2d(U,c))/\sqrt n}
={1\over beta n}\log\sum_i A_i(U)e^{beta(N_n-2i)/\sqrt n}.       \tag{L}
```

Since `|C_n^+|=2^n`,

```math
z_n(U)\le P_{n,beta}(U)\le z_n(U)+{\log2\over beta}.           \tag{L1}
```

Thus limits of `min_U P_{n,beta}(U)` for all fixed `beta`, followed by `beta->infinity`, imply convergence.

**Strongest failure:** temperature smooths the inner codeword maximum but leaves the outer minimum over
arbitrary cosets. Guerra interpolation controls quenched compatible random ensembles; it does not commute
with `min_U`, and a block split has uncontrolled cross-edge signs. Equation (L) is a full rooted
coset-enumerator transform, while zero-root MacWilliams data do not control its minimum over roots.

**Strongest coding rescue:** complete regularity would make every `A_i(U)` depend only on `r(U)` and turn
(L) into a distance-regular quotient recurrence. Same-layer variation of `b(U)` already refutes exact
complete regularity here. Delsarte external distance is all-coset but controls only the endpoint, not (L).

## Confidence and recommendation

Paley `b=0` theorem: **0.98**. Unique-zero pinning at `1/2`: **0.99**. Value of one scalable-trap checkpoint:
**0.55**. Viability of constrained-overlap interpolation: **0.03**. **Withdraw the proof architecture; retain
only (STOP).**
