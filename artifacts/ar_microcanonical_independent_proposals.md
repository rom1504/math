# Independent microcanonical/exchangeable proposals for all-order recovery

Date: 2026-08-16.

Status: independent, ledger-blind proposal freeze.  The three architectures in
Section 1 were fixed after reading only
`ar_sampling_rounding_exchangeability_toolkit.md`,
`minimal_all_order_action_recovery.md`, and
`exchangeable_recovery_obstruction.md`, and before the literature check used
below.  No archive, steering file, ledger, or active-state file was read.

## 0. Problem and a non-negotiable boundary

Let \(\Omega_n\) be the \(2^{N_n}\) symmetric hollow signings of order
\(n\), where \(N_n=\binom n2\), and put

```math
T_A=A/\sqrt n,
\qquad
\Phi(T_A)=\frac{2Q(A)}{n^{3/2}},
\qquad
a_n=\frac{M_n}{n^{3/2}},
\qquad
L=\liminf_n a_n.
```

The useful limit objects are the purified clusters \(T_\eta\) for which

```math
\|T_\eta\|_{2\to2}\le C_\eta,
\qquad
2L\le \Phi(T_\eta)\le 2L+\eta.                 \tag{0.1}
```

There is no viable projective-exchangeability proposal.  A single jointly
exchangeable infinite signing whose restrictions have
\(\|A_n\|_{\rm op}/\sqrt n\) tight is necessarily the i.i.d. Rademacher
array, and its restrictions satisfy

```math
\liminf_n Q(A_n)/n^{3/2}
\ge \frac23\sqrt{\frac2\pi}
=0.531923\ldots>\frac12.                         \tag{0.2}
```

Thus every proposal below uses a **different, vertex-exchangeable law at
each order** (or a deterministic signing followed by a uniform random
relabeling).  None may be coupled as the restrictions of one infinite array.
This loss of Kolmogorov consistency is exactly where the all-order work has
to occur.

## 1. Frozen list (at most three)

1. **Canonical directed-profile pressure**: prove a thermodynamic limit for
   a bounded mismatch Hamiltonian under the uniform signing measure.  This is
   the best proposal.
2. **Hard microcanonical profile-shell entropy**: prove existence of the
   entropy density of every fixed-radius directed-profile shell.  This is a
   clean large-deviation formulation, but it hides non-emptiness more directly
   than proposal 1.
3. **Balanced finite-type blow-up with a joint residual absorber**: realize a
   large finite prototype at all multiples by exact block-balanced signs while
   preventing the microscopic residual from creating new outer profiles.
   This is constructive but presently the least plausible.

Ordinary graphon sampling and projective graphon sampling are not a fourth
proposal: at the \(A/\sqrt n\) normalization they respectively erase the
fluctuation object and force the obstructed i.i.d. object.

## 2. Best proposal: canonical directed-profile pressure

### 2.1 A structural Hamiltonian

Fix one \(T=T_\eta\) from (0.1), and fix \(D>C_\eta\).  For
\(A\in\Omega_n\), define

```math
R_{T,D}(A)
:=\partial_1(T_A,T)
 +\min\{1,(\|T_A\|_{2\to2}-D)_+\}.              \tag{2.1}
```

This is relabeling-invariant and takes values in a fixed bounded interval.
Let

```math
m_n(T,D):=\min_{A\in\Omega_n}R_{T,D}(A),          \tag{2.2}
```

and, for an integer inverse temperature \(\beta\ge1\), let

```math
Z_n(\beta;T,D)
:=2^{-N_n}\sum_{A\in\Omega_n}
   \exp\{-\beta N_nR_{T,D}(A)\},

F_n(\beta;T,D)
:=-\frac{1}{\beta N_n}\log Z_n(\beta;T,D).       \tag{2.3}
```

The corresponding Gibbs law is finitely exchangeable at order \(n\):

```math
\mu_{n,\beta}^{T,D}(A)
\propto \exp\{-\beta N_nR_{T,D}(A)\}.            \tag{2.4}
```

It is intentionally not projective in \(n\).

### 2.2 The exact missing lemma

Only a null sequence of purified clusters is needed.

> **Lemma \(L_{\rm can}\) (profile-pressure thermodynamic limit).**  There
> are \(\eta_\ell\downarrow0\), choices of purified extremal clusters
> \(T_\ell=T_{\eta_\ell}\), and finite constants
> \(D_\ell>C_{\eta_\ell}\), such that for every \(\ell\) and every integer
> \(\beta\ge1\), the following all-order limit exists:
>
> ```math
> \boxed{
> F_\ell(\beta)
> :=\lim_{n\to\infty}F_n(\beta;T_\ell,D_\ell)
> \quad\text{exists in }\mathbb R.}               \tag{L_can}
> ```

The lemma does not ask for the value of the limit, a rate-function
minimizer, ensemble equivalence, concentration of the Gibbs law, or a
projective realization.

### 2.3 Exact known inputs plus \(L_{\rm can}\) imply convergence

The known/proved inputs are precisely:

1. **Purification and action compactness:** (0.1), witnessed by exact
   signings \(B_j\) with \(T_{B_j}\to T\) and
   \(\|T_{B_j}\|_{2\to2}\le C_\eta<D\).
2. **One-sided action continuity:** whenever the two operator norms are at
   most \(D'\) and \(\partial_1(S,T)\le\delta\),

   ```math
   \Phi(S)\le\Phi(T)+5D'\sqrt\delta+\delta.        \tag{2.5}
   ```

3. **The exact normalization:**
   \(\Phi(T_A)=2Q(A)/n^{3/2}\).
4. **Finite-volume Laplace sandwich:** for every finite set of size
   \(2^{N_n}\), (2.2)--(2.3) give

   ```math
   m_n(T,D)
   \le F_n(\beta;T,D)
   \le m_n(T,D)+\frac{\log2}{\beta}.               \tag{2.6}
   ```

Here (2.6) is elementary: bound the average by its largest summand for the
left inequality and retain one minimizing summand for the right inequality.

Along the purified subsequence, \(R_{T,D}(B_j)\to0\), hence
\(m_{n_j}(T,D)\to0\).  If \(L_{\rm can}\) holds, (2.6) therefore yields, for
each fixed \(\beta\),

```math
0\le F_\ell(\beta)\le\frac{\log2}{\beta},
\qquad
\limsup_{n\to\infty}m_n(T,D)\le F_\ell(\beta).    \tag{2.7}
```

Sending the integer \(\beta\to\infty\) gives \(m_n(T,D)\to0\).  Choose
\(C_n\in\Omega_n\) with
\(R_{T,D}(C_n)\le m_n(T,D)+1/n\).  Then

```math
\partial_1(T_{C_n},T)\to0,
\qquad
\|T_{C_n}\|_{2\to2}\le D+o(1).                  \tag{2.8}
```

Equations (0.1), (2.5), and (2.8) imply

```math
\limsup_n\frac{M_n}{n^{3/2}}
\le\limsup_n\frac{Q(C_n)}{n^{3/2}}
\le\frac12\Phi(T)
\le L+\frac\eta2.                                \tag{2.9}
```

Apply this to \(\eta=\eta_\ell\downarrow0\).  The reverse inequality is the
definition of \(L\), so \(M_n/n^{3/2}\) converges.

Thus the exact logical chain is

```math
\boxed{
\text{purified extremal cluster}
+\text{one-sided action continuity}
+L_{\rm can}
\Longrightarrow a_n\text{ converges}.}            \tag{2.10}
```

Notice the order of limits: first \(n\to\infty\) at each fixed \(\beta\),
then \(\beta\to\infty\).  This is what makes the proposal non-tautological.

### 2.4 Why \(L_{\rm can}\) uses strictly less raw information

The Hamiltonian (2.1) never uses \(Q(A)\), \(M_n\), or the map
\(x\mapsto H_A(x)\).  It retains only

- a one-sided distance from the unlabeled one-profile of one fixed cluster;
- one scalar operator-norm excess; and
- after averaging over all signings, one scalar log-Laplace value for each
  integer \(\beta\).

It discards vertex labels, reverse profile inclusion, all joint
\(k\)-profiles for \(k\ge2\), the identity of every contributing signing,
and every individual Boolean energy.  At fixed accuracy, the input alphabet
can be gridded with
\(q=O((1+D)/\epsilon^2)\), and output laws can be truncated and tested by a
finite bounded-Lipschitz net whose size is independent of \(n\).  Thus the
**state being prescribed** is a finite quotient at fixed accuracy, not a
table of \(2^{n-1}\) Boolean energies and not a target-order minimizer.

This is an information reduction, not evidence that the lemma is easier.
The outer condition still ranges over \(q^n\) colorings.  In particular, a
proof of \(L_{\rm can}\) may be just as hard as convergence unless it exploits
the canonical averaging without reconstructing those colorings one by one.

### 2.5 Plausible literature route, and the exact scale mismatch

There are two plausible entry points.

**Direct pressure interpolation.**  Guerra--Toninelli's theorem for the SK
and \(p\)-spin models proves existence of the thermodynamic free-energy
limit by a smooth interpolation between an \(n\)-site system and two
independent subsystems, producing subadditivity of the quenched free energy
([arXiv:cond-mat/0204280](https://arxiv.org/abs/cond-mat/0204280)).  The
desired analogue would establish an \(o(n^2)\)-defect comparison for the
profile-mismatch pressure (2.3).  The obstacle is concrete: \(R_{T,D}\) is a
Hausdorff-type supremum over all vertex colorings, and cross-block signs do
not yield a covariance interpolation with a known derivative sign.

**An action-profile LDP followed by Varadhan's lemma.**  Chatterjee--Varadhan
Theorem 2.3 gives an LDP for \(G(n,p)\) in graphon cut topology at speed
\(n^2\), with good rate

```math
I_p(W)=\frac12\int
\left[W\log\frac Wp+(1-W)\log\frac{1-W}{1-p}\right].       \tag{2.11}
```

([arXiv:1008.1946](https://arxiv.org/abs/1008.1946)).  Recent probability-
graphon work extends this Sanov-type statement to compact Polish edge marks
with a relative-entropy rate
([Dionigi--Zucal, arXiv:2509.14204](https://arxiv.org/abs/2509.14204)).
The Dawson--Gärtner projective-limit theorem says that good LDPs for all
finite projections yield an LDP in the projective-limit topology with rate
\(I(x)=\sup_j I_j(p_jx)\).  These results suggest:

1. grid \([-1,1]\) and truncate outputs at fixed \((D,\epsilon)\);
2. prove an \(n^2\)-speed LDP for each finite colored one-profile observable;
3. handle the set of all colorings in Hausdorff topology;
4. pass through the finite-profile projective system; and
5. apply Varadhan's lemma to (2.3).

Steps 2--3 are the new mathematics.  Existing colored-graphon LDPs concern
bounded kernels and ordinary edge densities.  Here the kernel values are
\(\sqrt n\,a_{ij}\); all centered signings have ordinary graphon limit zero,
while their action profiles and \(n^{3/2}\) ground states differ.  Therefore
neither (2.11) nor the probability-graphon extension can simply be contracted
to give \(L_{\rm can}\).  Backhausz--Szegedy's concentration of i.i.d. sign
matrices around deterministic action representatives still takes a
subsequence for convergence and does not supply rare low-profile shells
([arXiv:1811.00626](https://arxiv.org/abs/1811.00626), Section 11).

### 2.6 Finite and scalable falsifier

Use a large computed prototype \(B_k\) from a low subsequence as the finite
surrogate for \(T\).  For a grid alphabet \(G_q\), clipped outputs, and a
finite bounded-Lipschitz test net, define the computable outer-profile
penalty \(R^{(q,h)}_{B_k,D}(A)\).

- For \(n\le7\), enumerate all \(2^{N_n}\) signings and all
  \(q^n\) colorings (quotienting by relabeling and global sign), obtaining
  certified values of \(F_n^{(q,h)}(\beta)\).
- Beyond enumeration, use a counterexample-guided loop.  A master
  SAT/MILP chooses \(A\); an adversarial coloring oracle finds a violated
  profile constraint; that coloring is added as a cut.  Branch-and-bound or
  an SDP relaxation supplies certified lower bounds.  Thermodynamic
  integration or sequential Monte Carlo under (2.4) supplies non-rigorous
  estimates at larger \(n\).
- The pressure lemma is falsified if, for some fixed
  \((q,h,\beta)\), certified upper and lower intervals on two growing
  subsequences remain separated by a fixed positive gap.  The recovery
  consequence is falsified more directly if certified lower bounds on
  \(m_n^{(q,h)}\) stay above a fixed \(\epsilon\) on a growing residue-class
  subsequence, while prototype orders have penalty tending to zero.

A useful earlier test of the interpolation route is to enumerate the defect
in any proposed two-block subadditivity inequality.  One violating pair
\((n_1,n_2)\) rules out that claimed inequality and its constants, even
though it does not by itself disprove \(L_{\rm can}\).

### 2.7 Circularity traps

1. Do not prove pressure convergence using convergence of
   \(m_n(T,D)\), \(M_n\), or the existence of target-order low-cap signings.
2. Do not replace (2.1) by \(|\Phi(T_A)-\Phi(T)|\); that is the target scalar
   in disguise and loses the structural reduction.
3. Do not take \(\beta\to\infty\) at fixed \(n\) first.  That merely returns
   the target-order structural minimum \(m_n\).
4. Do not omit the norm penalty.  Directed LP proximity alone can hide a
   high-norm spike, and (2.5) then has no useful uniform error.
5. Do not infer a hard microcanonical statement from a canonical one by
   unproved ensemble equivalence.  Dense graph ensembles exhibit genuine
   equivalence breaking
   ([arXiv:1703.08058](https://arxiv.org/abs/1703.08058)).
6. Do not realize the Gibbs laws as restrictions of one exchangeable array;
   (0.2) then applies.
7. Do not contract an ordinary graphon LDP after centering.  It maps both the
   extremal candidates and i.i.d. signs to the same zero graphon.

## 3. Backup: hard microcanonical profile-shell entropy

For fixed \((T,D,\epsilon)\), let

```math
\Gamma_n(T,D,\epsilon)
:=\{A\in\Omega_n:
       \partial_1(T_A,T)<\epsilon,
       \ \|T_A\|_{2\to2}\le D\},                  \tag{3.1}
```

and use the convention \(\log0=-\infty\) in

```math
s_n(T,D,\epsilon)
:=\frac1{N_n}\log|\Gamma_n(T,D,\epsilon)|.         \tag{3.2}
```

The missing hard-shell lemma would say that, for each selected purified
cluster and some \(\epsilon_r\downarrow0\),

```math
\lim_{n\to\infty}s_n(T,D,\epsilon_r)
\quad\text{exists in }[-\infty,\log2].             \tag{L_mc}
```

Along the realizing subsequence the shell is nonempty, so
\(s_{n_j}\ge0\).  Consequently the limit in \(L_{\rm mc}\) cannot be
\(-\infty\), and every sufficiently large order has a shell member.
Diagonalizing in \(r\) gives (2.8), hence convergence by the same argument as
above.

This is the direct microcanonical/LDP architecture.  Its attraction is that
standard graphon LDPs can give hard-shell count limits for cut-continuous
bounded statistics at rate-continuity radii.  Its defect is that the convention
\(\log0=-\infty\) makes existence of the entropy limit carry much of the
desired no-gap conclusion.  A proof based on concatenating two nonempty
shells would already need to sign all cross edges without creating a new
outer profile.  That is the original obstruction in entropic notation.
Microcanonical/canonical equivalence cannot be assumed, and fixed motif or
color-density constraints do not control the \(q^n\) outer quantifier.

The same CEGIS shell search as in Section 2.6 is a falsifier: certified
emptiness of successively larger fixed-radius shells on one residue class,
contrasted with nonempty prototype shells, attacks \(L_{\rm mc}\) directly.

## 4. Backup: balanced finite-type blow-up plus joint residual absorption

Take a large finite signing \(B\) of order \(k\) from the purified cluster.
For \(n=kr\), split the target vertices into \(k\) equal blocks.  To reproduce
\(T_B=B/\sqrt k\) on block-constant functions, a row in block \(a\) must
have signed sum approximately

```math
\sum_{j\in V_b}c_{ij}\simeq b_{ab}\sqrt r             \tag{4.1}
```

toward block \(b\).  Equivalently, each \(r\times r\) block has microscopic
bias about \(b_{ab}/\sqrt r\).  Exact row/column balances and parity defects
are design-type constraints and are plausibly repairable.

The genuinely new lemma is not row balancing but the following
**sequence-level** absorber.  It cannot reasonably be demanded at arbitrary
accuracy for one fixed finite \(B\): laws on a \(k\)-point space have atom
granularity \(1/k\), while a large lift admits much finer input laws.

> **Joint residual absorber.**  For a realizing sequence
> \(B_j\) of orders \(k_j\to\infty\) with \(T_{B_j}\to T\), there are a
> uniform \(D<\infty\), errors \(\rho_j\downarrow0\), and thresholds
> \(R_j<\infty\) such that, for every \(r\ge R_j\), a symmetric hollow
> signing \(C_{k_jr}\) satisfies the block balances (4.1) up to \(O(1)\) and
>
> ```math
> \boxed{
> \|T_{C_{k_jr}}\|_{2\to2}\le D,
> \qquad
> \partial_1(T_{C_{k_jr}},T_{B_j})\le\rho_j.}       \tag{L_abs}
> ```

Choose \(j=j(n)\to\infty\) so slowly that
\(r=\lceil n/k_j\rceil\ge R_j\) and \(k_j=o(n)\).  The constructed order
\(m=k_jr\) is \(n+o(n)\); the triangle inequality through
\(T_{B_j}\to T\), one-sided continuity, and lossless deletion from \(m\) to
\(n\) then prove convergence.  The uniform norm bound in \(L_{\rm abs}\)
can be weakened to the corresponding
\(D_j\sqrt{\rho_j+d_M(T_{B_j},T)}\to0\) condition.

Known balancing/design theorems can enforce finitely many integer
statistics, and matrix discrepancy can keep a residual on the
constant-times-\(\sqrt n\) scale.  Neither implies \(L_{\rm abs}\).  A random
balanced fill leaves a Wigner-like microscopic residual, and the i.i.d.
lower bound (0.2) warns that this residual is not harmless.  A norm error
\(c\sqrt n\) also changes the objective at leading order.  The absorber must
instead prove a **same-spin, joint cancellation before the supremum**, for
every fixed-alphabet coloring at once.

This architecture has a particularly concrete falsifier.  For small
\((k,r)\), impose (4.1) in SAT/MILP, enumerate or branch over all balanced
lifts, and use a coloring oracle to maximize profile violation (with
\(Q(C)\) as the first necessary witness).  If the minimum violation among
all balanced lifts stays bounded below as \(r\) grows along the increasingly
accurate prototypes of a proposed realizing sequence, that instance of
\(L_{\rm abs}\) is false.  The orbit of any successful lift gives finite
exchangeability, but trying to nest these orbits
projectively is prohibited by (0.2).

## 5. Verdict

No cited theorem currently supplies any of the three missing lemmas.  The
projective i.i.d.-sign obstruction is decisive, not cosmetic; ordinary dense
sampling and graphon LDPs are at the wrong normalization; and fixed
microcanonical constraints do not see the universal coloring profile.

The canonical pressure lemma \(L_{\rm can}\) is the best research target.  It
is exact, all-order, invariant under relabeling, strictly coarser in raw state
than the Boolean landscape, and its implication uses only the elementary
Laplace sandwich plus the already proved directed continuity theorem.  Its
most plausible proof would be a new action-profile thermodynamic-limit or
large-deviation theorem.  Until such a theorem controls the \(q^n\) outer
quantifier at the \(A/\sqrt n\) scale, the honest conclusion is: **there is a
sharply formulated route, but no presently validated recovery architecture
in exchangeability/graphon microcanonical theory.**
