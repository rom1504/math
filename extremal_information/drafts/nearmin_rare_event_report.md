# Rare-event renormalization versus finite-depth contextual information

**Scoped verdict.** Yes, at benchmark Level 3.  A deterministic product
landscape with all translations allowed as adversarial contexts has an exact
renormalized extremal state: its log-moment function adds under product
composition, and its Legendre transform gives the exact exponential cost of
an upper-tail witness.  At every finite depth, however, answering all
translated threshold queries requires an exponentially large witness library.
For libraries that nearly minimize this contextual problem at exponential
scale, the same mechanism forces a precise **subextensive contextual witness
entropy**.  This conclusion is genuinely different from, and does not imply,
the signed edge first-moment balance of `FB.1`.

The model has a literal dense-edge realization, but its composition is a
product of independent edge coordinates.  It is therefore not a theorem
about dense sign quadratics.  The cut constraint and the new cross-block edges
in a quadratic signing destroy the exact product law.  Consequently the
result is a positive benchmark, not a Level 5 near-minimizer lemma for
`Q(A)`.

Parts 1--2 below reprove, in self-contained form, the orbit-query mechanism
already isolated in `orbit_query_large_deviations.md`.  The scoped addition
needed for the present question is Part 3: the near-minimizer witness-entropy
lemma, together with the dense-edge specialization and transfer audit.

## One benchmark theorem

Let `G` be a finite abelian group of order `q`, let `f:G->R` be nonconstant,
and put

```math
 F_n(x)=\sum_{i=1}^n f(x_i),\qquad x\in G^n.                 \tag{1}
```

Fix

```math
 \mathbb E f<a<\max f,                                      \tag{2}
```

where the expectation is under the uniform law on `G`, and define the rare
upper set and its density by

```math
 W_n(a)=\{x:F_n(x)\ge an\},\qquad
 p_n(a)=|W_n(a)|/q^n.                                      \tag{3}
```

A library `C subset G^n` answers every adversarial translated threshold query
if

```math
 \forall s\in G^n\quad \max_{c\in C}F_n(s+c)\ge an.        \tag{4}
```

Let `L_n(a)` be the minimum cardinality of such a library.  For a covering
library define its witness multiplicity at context `s` by

```math
 N_C(s)=|\{c\in C:s+c\in W_n(a)\}|.                         \tag{5}
```

Finally put

```math
 \Lambda_f(\theta)=
 \log\left(q^{-1}\sum_{z\in G}e^{\theta f(z)}\right),
 \qquad
 I_f(a)=\sup_{\theta\ge0}\{\theta a-\Lambda_f(\theta)\}.
                                                                    \tag{6}
```

### Theorem (orbit-covering renormalization and near-minimizer witness entropy)

The following statements hold.

1. **Exact finite-depth contextual cost and its rate.**

   ```math
   {1\over p_n(a)}\le L_n(a)
   \le \left\lceil {n\log q+1\over p_n(a)}\right\rceil,
   \qquad
   \lim_{n\to\infty}{1\over n}\log L_n(a)=I_f(a)>0.       \tag{7}
   ```

   Thus the exact context carrier is exponentially large.

2. **Exact renormalized composition.**  For arbitrary finite landscapes
   `H:X->R` and `J:Y->R`, with uniform reference laws, define

   ```math
   K_H(\theta)=\log\left(|X|^{-1}\sum_{x\in X}e^{\theta H(x)}\right).
                                                                    \tag{8}
   ```

   Under product composition `(H\boxplus J)(x,y)=H(x)+J(y)`,

   ```math
   \boxed{K_{H\boxplus J}=K_H+K_J}.                         \tag{9}
   ```

   In particular, `K_{F_n}=n\Lambda_f`, and the exposed upper-tail
   information rate is its Legendre dual `I_f`.  The rate state therefore
   composes while the exact finite-depth context library in (7) does not
   compress.

3. **Near-minimizer lemma.**  If `S` is uniform on `G^n`, every covering
   library satisfies the exact bounds

   ```math
   0\le \mathbb E\log N_C(S)
   \le \log\bigl(|C|p_n(a)\bigr),                           \tag{10}
   ```

   and, for every integer `k>=2`,

   ```math
   \mathbb P\{N_C(S)\ge k\}
   \le {|C|p_n(a)-1\over k-1}.                              \tag{11}
   ```

   Consequently, if `C_n` is a genuine exponential-scale near-minimizer,

   ```math
   \log|C_n|\le \log L_n(a)+\delta_n n,
   \qquad \delta_n\longrightarrow0,                        \tag{12}
   ```

   then

   ```math
   \mathbb E\log N_{C_n}(S)=o(n),                           \tag{13}
   ```

   and for each fixed `t>0`,

   ```math
   \mathbb P\{N_{C_n}(S)\ge e^{tn}\}
   \le \exp\{-(t-o(1))n\}.                                \tag{14}
   ```

   At the sharper multiplicative sphere bound
   `|C|p_n(a)<=1+eta`, one has

   ```math
   \mathbb P\{N_C(S)>1\}\le eta,
   \qquad \mathbb E\log N_C(S)\le\log(1+eta).             \tag{15}
   ```

   Hence rate-near-minimality forces subextensive average ambiguity of the
   extremal witness, even though storing witnesses sufficient for every
   individual context requires `exp(nI_f(a)+o(n))` entries.

### Proof

For each `c`, the contexts answered by `c` form the translate `W_n(a)-c`,
which has cardinality `p_n(a)q^n`.  Since these translates cover `G^n`,
counting incidences gives

```math
 |C|p_n(a)q^n\ge q^n,                                      \tag{16}
```

and proves the lower bound in (7).  Conversely choose `k` independent
uniform points of `G^n`.  A fixed context is uncovered with probability
`(1-p_n(a))^k<=exp(-kp_n(a))`.  Thus the expected number of uncovered
contexts is at most

```math
 q^n e^{-kp_n(a)}<1                                       \tag{17}
```

when `k=ceil((n log q+1)/p_n(a))`.  Some choice covers all contexts, proving
the upper bound.

It remains to identify the exponent.  If `Z_1,...,Z_n` are uniform on `G`,
then `p_n(a)=P{sum_i f(Z_i)>=an}`.  Chernoff's inequality gives

```math
 p_n(a)\le \exp\{-n(\theta a-\Lambda_f(\theta))\}
 \quad(\theta\ge0).                                      \tag{18}
```

For the reverse exponential bound, group points of `G` with equal `f` value.
A type `nu` has probability

```math
 \exp\{-nD(\nu\Vert\mu)+O(\log n)\},                     \tag{19}
```

where `mu` is the uniform pushforward law of `f`.  There are only polynomially
many types.  Minimizing `D(nu||mu)` subject to `E_nu f>=a`, approximating a
minimizer by rational types, and using finite-dimensional entropy duality
gives

```math
 -{1\over n}\log p_n(a)\longrightarrow
 \inf_{\mathbb E_\nu f\ge a}D(\nu\Vert\mu)
 =I_f(a).                                                 \tag{20}
```

The strict inequalities in (2) make `I_f(a)>0`.  Combining (20) with the two
bounds in (7) proves the rate assertion.  No probabilistic environment is
being averaged here; randomness in (17)--(20) is only a proof device for a
fully deterministic covering problem.

Equation (9) is the factorization

```math
 {|X|}^{-1}{|Y|}^{-1}\sum_{x,y}e^{\theta(H(x)+J(y))}
 =\left({|X|}^{-1}\sum_xe^{\theta H(x)}\right)
  \left({|Y|}^{-1}\sum_ye^{\theta J(y)}\right).             \tag{21}
```

For the near-minimizer statement, (4) says `N_C(s)>=1` for every `s`, while
translation invariance and double counting give

```math
 \mathbb E N_C(S)=|C|p_n(a).                               \tag{22}
```

Jensen's inequality applied to `log` proves (10).  Also

```math
 (k-1)1_{\{N_C\ge k\}}\le N_C-1,                          \tag{23}
```

and averaging proves (11).  From (7),

```math
 1\le L_n(a)p_n(a)\le n\log q+1+p_n(a),                   \tag{24}
```

so (12) implies

```math
 \log(|C_n|p_n(a))\le\delta_n n+O(\log n)=o(n).            \tag{25}
```

Equations (10)--(11), with `k=ceil(e^{tn})`, now give (13)--(14).
If `|C|p_n(a)<=1+eta`, (22) and `N_C>=1` give
`P{N_C>1}<=E(N_C-1)<=eta`, while (10) gives the other half of (15).
This completes the proof. `square`

## Dense-edge specialization

Take a complete graph on `m` vertices, let `N=binom(m,2)`, identify every
dense exact edge signing with a point of `G^N` for `G=Z/2Z`, and take
`f(0)=1`, `f(1)=-1`.  Then

```math
 F_N(B+C)=N-2d_H(B,C),\qquad \Lambda_f(\theta)=\log\cosh\theta. \tag{26}
```

For `0<a<1`, a library answers every dense edge-sign context exactly when it
is a covering code of relative radius `(1-a)/2` in the `N`-dimensional edge
cube.  Here

```math
 I_f(a)=D\left({1+a\over2}\middle\Vert{1\over2}\right),     \tag{27}
```

so an exact finite-depth library has size `exp(Theta(N))=exp(Theta(m^2))`,
whereas the renormalized extremal state is the one convex function
`log cosh(theta)`.  This is a deterministic adversarial landscape on dense
exact-sign objects, not a random-energy assertion.

The qualification is essential: disjoint union in (9) composes edge
coordinates but creates no new cross edges.  It is not the block composition
of two complete signed quadratic forms.

## Why the near-minimizer lemma is not shell first-moment balance

The shell theorem `FB.1` concerns a probability law on near-top signed cuts
and controls the signed edge marginals
`E_mu[a_e z_e]`.  Equations (10)--(15) instead concern the number of extremal
witnesses available after an adversary chooses a context.  They use no
signed average and retain no edge marginal.  The two controls are logically
independent even in elementary product examples:

* Balanced coordinate marginals do not control contextual ambiguity.  In
  the binary model take `C=G^n`.  Its nontrivial coordinate characters have
  zero mean, but `N_C(s)=|W_n(a)|` for every `s`, exponentially large for an
  interior threshold.
* Unique contextual witnesses do not force coordinate balance.  Let
  `G=Z/6Z`, `f=1` on `{0,1}` and `f=0` elsewhere, and use the boundary
  threshold `a=1`.  Then `W_n={0,1}^n` and
  `C={0,2,4}^n` give a unique decomposition for every context, so `N_C=1`
  identically, while the natural parity character `(-1)^c` equals `+1` on
  every coordinate of every library point.

These examples are not counterexamples to `FB.1`, whose witnesses live in
the cut variety of a quadratic signing.  They certify the narrower claim
needed here: witness entropy is an incomparable near-minimizer observable,
not a disguised first-moment shell balance.

## Benchmark level and transfer judgment

**Benchmark level: Level 3 (growing adversarial interface with asymptotic
compression), not Level 4 or 5.**  The dense-edge specialization has exact
signs and `Theta(m^2)` coordinates, but it does not have the constrained
quadratic response

```math
 Q(A)=\max_{x\in\{+-1\}^m}|\langle A,(x_ix_j)_{i<j}\rangle|. \tag{28}
```

There are three precise failures of transfer.

1. The context orbit in the theorem is the full edge cube.  Quadratic
   response queries only the `2^(m-1)` unsigned cut vectors `(x_i x_j)` (or
   at most `2^m` after adjoining the overall sign needed for the absolute
   value), a very thin algebraically constrained subset of the `2^N` edge
   cube.
2. Product composition in (9) assumes additive energies and independent
   coordinates.  Gluing complete quadratic forms introduces all
   `m_1m_2` cross edges; their response depends on the relative spin geometry
   and is not determined by `K_H+K_J`.
3. The near-minimization in (12) optimizes a context library for a fixed rare
   set.  The signing problem optimizes the disorder `A` itself.  No step above
   turns `Q(A)<=M_m+epsilon m^(3/2)` into (12).

There is also an architectural limitation: (7) lower-bounds the size of an
explicit witness library.  It does not rule out a short algorithmic or
algebraic description that generates its witnesses.  “Finite-depth
incompressibility” here is therefore library incompressibility, exactly as
declared in (4), not a universal data-structure lower bound.

Accordingly, the theorem answers the scoped existence question positively
and supplies a precise incomparable near-minimizer lemma, but it makes **NO
FRONTIER CHANGE** for the current Level 5 smallest missing lemma.  A transfer
would require a new statement showing that genuine low-cap cut landscapes
admit an approximately additive log-moment/rate state under bounded-cap
cross-block gluing.  That statement is neither proved nor suggested by the
product proof above.

## Literature dependence

No imported extremal-process, metastate, or Gamma-convergence theorem is
needed.  The only large-deviation input is proved directly by the finite-type
calculation (18)--(20).  This avoids borrowing random-environment conclusions
whose hypotheses fail for adversarially selected dense signings.
