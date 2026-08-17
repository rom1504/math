# Zero-rate quadratic upper tails are carried by zero-density principal cores

**Status.** Task-local theorem draft awaiting independent audit.  This note
uses the archived Grothendieck--Pietsch factorization to identify the precise
structural branch left by the Hanson--Wright tail/spike dichotomy.  It does
not prove the exact-minimizer tail lemma.

## 1. Setup

For a real symmetric hollow order-`n` matrix `A` with `|A_(ij)|<=1`, put

```math
H_A(x)={1\over2}x^TAx,
\qquad
P(A)=\max_xH_A(x),
\qquad
Q(A)=\max_x|H_A(x)|.                              \tag{PC.1}
```

For `T subseteq[n]`, let `A[T]` be the principal submatrix and retain the
one-sided notation `P(A[T])`.  At fixed global scale `t>0`, write

```math
p_A(t)=2^{-n}\#\{x:H_A(x)\ge t n^{3/2}\}.         \tag{PC.2}
```

The distinction between `P` and `Q` matters: the tail has a declared global
orientation.

## 2. Factorized conditional concentration

### Theorem PC.1 (principal-core gap forces a linear entropy deficit)

Fix `C,t>0`, `0<epsilon<=1`, and `0<eta<t`.  There is a constant
`c=c(K_G)>0` with the following property.  If

```math
Q(A)\le Cn^{3/2},                                  \tag{PC.3}
```

then there is a set `T=T(A,epsilon)` with `|T|<epsilon n` such that

```math
P(A[T])\le(t-eta)n^{3/2}                           \tag{PC.4}
```

implies

```math
\boxed{
p_A(t)\le3\exp\left[-c n\min\left\{
 eta^2,{eta epsilon\over C},{eta^2epsilon\over C^2}
\right\}\right].}                                \tag{PC.5}
```

The set is obtained without looking at a maximizing Boolean spin.  It is the
heavy-coordinate set of one Grothendieck--Pietsch factorization of `A`.

#### Proof

Let

```math
B=\|A\|_(infinity to1)
 =\max_{x,y\in\{+-1\}^n}|x^TAy|.
```

Box polarization in the normalization (PC.1) gives `B<=4Q(A)`.  The real
Grothendieck factorization theorem supplies probability vectors `mu,nu` for
which

```math
|u^TAv|\le K_GB
 \left(\sum_i\mu_i u_i^2\right)^{1/2}
 \left(\sum_j\nu_j v_j^2\right)^{1/2}             \tag{PC.6}
```

for all real `u,v`.  Delete

```math
T=\{i:\mu_i>2/(epsilon n)\}
  \mathbin\cup\{i:\nu_i>2/(epsilon n)\},
\qquad R=[n]\setminus T.                          \tag{PC.7}
```

Then `|T|<epsilon n`, and (PC.6) gives

```math
\|A[R]\|_(2 to2)\le {2K_GB\over epsilon n}
                  \le {8K_GC\over epsilon}\sqrt n.          \tag{PC.8}
```

Fix any spin `x_T` and put `h=A_(R,T)x_T`.  A second use of (PC.6), now
with the `T`-factor bounded by total `nu`-mass one, gives

```math
\|h\|_2\le K_GB\sqrt{2\over epsilon n}
          \le {4\sqrt2K_GC\over\sqrt epsilon},n.            \tag{PC.9}
```

For uniform `X_R`, the exact conditional decomposition is

```math
H_A(x_T,X_R)=H_(A[T])(x_T)+h^TX_R+H_(A[R])(X_R).   \tag{PC.10}
```

Under (PC.4), the event on the left of (PC.2) forces either

```math
h^TX_R\ge {eta\over2}n^{3/2}
```

or

```math
H_(A[R])(X_R)\ge {eta\over2}n^{3/2}.              \tag{PC.11}
```

The Rademacher linear tail and (PC.9) bound the first probability by

```math
\exp\{-c eta^2epsilon n/C^2\}.                    \tag{PC.12}
```

For the second, Hanson--Wright, (PC.8), and
`||A[R]||_F^2<=n^2` give

```math
2\exp\{-c n\min(eta^2,eta epsilon/C)\}.           \tag{PC.13}
```

Both estimates are uniform in `x_T`.  Averaging the conditional probability
over `x_T` proves (PC.5). `square`

The only imported inputs are the real Grothendieck--Pietsch factorization
already proved and cited in
`artifacts/orientation_even_grothendieck_localization.md` and the
Rademacher Hanson--Wright inequality already used in Theorem 21.8 and TS.1.

## 3. The converse core construction

### Proposition PC.2 (a core above threshold creates a zero-rate tail)

Let `A_n` be complete signings and let `T_n subseteq[n]` satisfy

```math
|T_n|=o(n),
\qquad
P(A_n[T_n])\ge(t+eta)n^{3/2}                      \tag{PC.14}
```

for fixed `t,eta>0`.  Then

```math
p_(A_n)(t)\ge2^{-|T_n|}(1-o(1))
            =\exp\{-o(n)\}.                       \tag{PC.15}
```

#### Proof

Fix a maximizing spin on `T_n` and make the complementary spins uniform.
The remaining cross-linear plus complementary-quadratic polynomial has mean
zero.  Orthogonality of distinct Boolean Fourier monomials gives variance at
most

```math
n|T_n|^2+n^2=o(n^3).                              \tag{PC.16}
```

(The first term is the worst possible coalescence of the `|T_n|` cross
edges incident with one free vertex.)  Chebyshev shows that a `1-o(1)`
fraction of the `2^(n-|T_n|)` extensions lose less than
`eta n^(3/2)`, proving (PC.15). `square`

## 4. Sequential zero-rate characterization

### Corollary PC.3 (tail failure forces a zero-density global-energy core)

Suppose `Q(A_n)<=Cn^(3/2)` and, for fixed `t>0`,

```math
r_n=-{1\over n}\log p_(A_n)(t)\longrightarrow0.   \tag{PC.17}
```

Then there are sets `T_n` with

```math
|T_n|=o(n),
\qquad
\boxed{P(A_n[T_n])\ge(t-o(1))n^{3/2}.}            \tag{PC.18}
```

Indeed put `delta_n=max(r_n,1/n)` and choose
`epsilon_n=eta_n=delta_n^(1/4)`.  The exponent in (PC.5) is at least a
constant times `delta_n^(3/4)n`, whereas
`r_n=o(delta_n^(3/4))` and `n delta_n^(3/4)->infinity`.  Hence (PC.4)
must fail for the factorization set `T(A_n,epsilon_n)`.  This gives
(PC.18).  Necessarily

```math
|T_n|\ge(\sqrt{2t}-o(1))n^{3/4},                  \tag{PC.19}
```

because a `k`-vertex signing has one-sided cap at most `k(k-1)/2`.

Together, PC.2--PC.3 give a slack-stable characterization:

```text
zero exponential upper-tail rate
  <=>
an o(n)-vertex principal core carries the global n^(3/2) energy scale,
```

where the reverse implication requires any fixed positive energy margin and
the forward implication loses `o(n^(3/2))`.  This is not a literal equality
at the threshold boundary.

## 5. Consequence for the exact-minimizer target

Take `t_0=c_--d_0` in TS.2.  The exact structural obstruction to the
`L_tail` route is now narrower than an unspecified spectral spike:

> **Exact-minimizer no-core lemma (`L_core`, minimal sequential form).**
> For some fixed `t_0>0`, every sequence of suitably oriented exact
> minimizers `A_n` and every sequence `T_n subseteq[n]` with `|T_n|=o(n)`
> satisfy
> ```math
> \limsup_{n\to\infty}{P(A_n[T_n])\over n^{3/2}}<t_0.
> ```

A fixed-gap formulation with constants `0<eta_0<t_0` is a convenient
stronger sufficient statement, but it is not the minimal target.  The
sequential statement retains only one-sided caps of zero-density principal
submatrices.  It neither identifies a global optimizer nor compares orders,
so it is an explicit inverse-structure target rather than a renamed full
landscape.  This does **not** by itself demonstrate that proving it is
strictly easier than the original tail assertion.  If proved, PC.3 gives
`L_tail`, and Theorem 21.8 gives the conditional boundary-response packing
in Theorem 36.22.

The quantifier may be discontinuous, and is certainly unstable below the
planted energy scale.  The archived clique overwrite and multi-clique
constructions implant **subleading** one-sided principal
cores inside broad `o(n^(3/2))` near-minimizer halos.  They do not furnish a
fixed-`t_0 n^(3/2)` core at subleading cost: a clique of the minimal
`Theta(n^(3/4))` size already changes `Theta(n^(3/2))` edges.  Thus those
examples warn that no-core statements are unstable below their planted
energy scale, but they do not falsify a fixed-level near-minimizer version.

## 6. Archive comparison and frontier classification

- Grothendieck--Pietsch common-support removal is archived; PC.1's new point
  is to retain the cross field in a conditional tail bound and charge the
  remaining obstruction to the **principal one-sided core energy**.
- TS.1 proves that tail failure is spectral.  PC.3 strengthens the structural
  conclusion to a zero-density principal core carrying the global energy
  scale; it does not merely relabel the top eigenvector.
- The weighted-clique and exact two-clique examples in the TS audit realize
  the same core mechanism outside the exact-minimizer class.
- No existing exact-minimality identity in the archive bounds (PC.18).

Classification:

```text
PROVES A NEW INVERSE/CONCENTRATION COROLLARY AND REFINES AN ARROW:
  L_tail is reduced to excluding one explicit zero-density principal-core
  obstruction rather than all possible upper-tail geometries.

BENCHMARK:
  the concentration theorem is uniform for all bounded-cap dense signings;
  the exact-minimizer application remains a conditional Level-5 target.
  This is not yet a frontier reset or a proved strict reduction.

NO CONVERGENCE IMPLICATION:
  even L_tail plus Theorem 21.8 still yields conditional profiles until the
  separate low-cap scalar selector is proved.
```
