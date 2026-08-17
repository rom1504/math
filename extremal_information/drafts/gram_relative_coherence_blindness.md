# Gram-relative synchronization needs joint coherence

**Status.**  Rigorous contraction-level falsifier and rigorous symmetry
criterion.  The falsifier is not an exact-sign or hollow child construction.

The relative product-synchronization certificate

```math
D\preceq\delta G
```

is intrinsically checkable once the **full** active-product Gram--Rayleigh
pair `(G,R)` is known, because `D=G-R`.  It therefore makes no sense to ask
for two systems with the same full pair and different relative defects.
The substantive question is whether the certificate follows from the
original ports, low-degree product data, individual product Rayleigh
deficits, or an average near-top condition.  The theorem below gives a
scalable negative answer at logarithmic port arity.

## 1. A generated Walsh pole algebra

Fix an odd integer

```math
p=2m+1
```

and put `N=2^(2p)`.  Index the `N` rows by

```math
(x,t)\in\{\pm1\}^p\times[2^p].
```

The base ports are

```math
w_i(x,t)=x_i.
```

For every odd set `S subseteq [p]`, let

```math
z_S=\prod_{i\in S}w_i,
```

and let `V` list these columns.  There are

```math
q=2^{p-1}={\sqrt N\over2}
```

of them.  They form the odd affine coset of the Boolean character algebra.
Character orthogonality gives

```math
E={V\over\sqrt N},\qquad E^TE=G=I_q.             \tag{CB.1}
```

Thus this is a generated product-pole algebra with
`p=(log_2 N)/2`, not an arbitrary list of unrelated vectors.

Let `tau=Maj_p` and use normalized Boolean Fourier coefficients

```math
a_S=\widehat\tau(S).
```

All odd sets are active.  For `|S|=2k+1`, the elementary symmetric-majority
calculation gives

```math
a_S=(-1)^k a_1{\binom m k\over\binom{2m}{2k}},
\qquad
a_1={\binom{2m}m\over2^{2m}}.                    \tag{CB.2}
```

One way to obtain (CB.2) is to group the Fourier sum by the number of minus
signs inside and outside `S` and apply the alternating Vandermonde identity.
In particular,

```math
|a_S|\le a_1,
\qquad
\sum_{|S|\ \mathrm{odd}}a_S^2=1.                 \tag{CB.3}
```

The second identity is Parseval.

Remove the degree-one coefficients:

```math
J=\{S:|S|\ge3\text{ is odd}\},
\qquad
\rho_p=\sum_{S\in J}a_S^2=1-pa_1^2.             \tag{CB.4}
```

The central-binomial asymptotic gives

```math
pa_1^2\longrightarrow {2\over\pi},
\qquad
\rho_p\longrightarrow1-{2\over\pi}>0.          \tag{CB.5}
```

Define a unit vector on the high-product coordinates by

```math
u_S=\begin{cases}a_S/\sqrt{\rho_p},&S\in J,\\0,&|S|=1.
\end{cases}                                      \tag{CB.6}
```

## 2. Coherent and diagonal defects

Consider the two positive defects

```math
D_{\rm coh}=uu^T,
\qquad
D_{\rm diag}=\operatorname{diag}(u_S^2),         \tag{CB.7}
```

and put `R_*=I-D_*`.  Both `R_*` are positive contractions.  They are
genuine Gram--Rayleigh compressions of symmetric positive contractions:
on the row space set

```math
T_*=ER_*E^T                                      \tag{CB.8}
```

and let `T_*` vanish on `range(E)^perp`.  Then

```math
\|T_*\|_{op}\le1,
\qquad
{V^TT_*V\over N}=R_*,
\qquad
G-R_*=D_*.                                       \tag{CB.9}
```

The matrices `T_*` are real symmetric weighted contractions.  In general
they are neither hollow nor entrywise signs; this limitation is essential
to the evidentiary status of the result.

### Theorem CB.1 (marginal near-top data do not control joint coherence)

The two systems in (CB.7)--(CB.9) have all of the following identical data.

1. Every matrix entry involving a degree-one pole: for `|S|=1` and every
   odd `T`,

   ```math
   (D_{\rm coh})_{ST}=(D_{\rm diag})_{ST}=0.
   ```

   Hence the complete generator Gram--Rayleigh block, including its cross
   entries against all higher products, is the same.

2. Every individual product Rayleigh deficit:

   ```math
   (D_{\rm coh})_{SS}=(D_{\rm diag})_{SS}=u_S^2.
                                                               \tag{CB.10}
   ```

3. The trace and average defects:

   ```math
   \operatorname{tr}D_{\rm coh}
   =\operatorname{tr}D_{\rm diag}=1,
   \qquad {1\over q}\operatorname{tr}D_*={1\over q}.          \tag{CB.11}
   ```

Indeed every individual pole has deficit at most

```math
\max_Su_S^2\le {a_1^2\over\rho_p}=O(p^{-1}).     \tag{CB.12}
```

Nevertheless their Gram-relative defects are macroscopically different:

```math
\boxed{
\delta_{\rm coh}=\|D_{\rm coh}\|_{op}=1,
\qquad
\delta_{\rm diag}=\|D_{\rm diag}\|_{op}
\le {a_1^2\over\rho_p}=O(p^{-1}).}               \tag{CB.13}
```

The discrepancy is exposed by the actual majority selector coefficient
vector at endpoint `epsilon=(1,...,1)`, namely `a=(a_S)_S`:

```math
a^TD_{\rm coh}a=\rho_p\longrightarrow1-{2\over\pi},
\qquad
a^TD_{\rm diag}a
={\sum_{S\in J}a_S^4\over\rho_p}
\le a_1^2=O(p^{-1}).                             \tag{CB.14}
```

#### Proof

Because `u` vanishes on singleton coordinates, the first assertion is
immediate.  Both diagonal tables are `(u_S^2)_S`, and `||u||_2=1`, proving
(CB.10)--(CB.12).  A rank-one projector has operator norm one, whereas the
operator norm of a diagonal positive matrix is its largest diagonal entry;
this proves (CB.13).  Finally,

```math
a^Tu={1\over\sqrt{\rho_p}}\sum_{S\in J}a_S^2=\sqrt{\rho_p},
```

while

```math
{1\over\rho_p}\sum_{S\in J}a_S^4
\le {\max_{S\in J}a_S^2\over\rho_p}
       \sum_{S\in J}a_S^2
\le a_1^2.
```

This gives (CB.14).  `square`

The theorem is stronger than an average-defect counterexample: even the
**complete individual deficit table** is identical and uniformly small.
What differs is one coherent off-diagonal direction.  Consequently no
argument from generator data plus separate near-top product estimates can
prove Gram-relative synchronization in this general contraction class.

Equation (CB.14) concerns the prescribed selector witness in the robust
product theorem.  It need not equal the full Boolean trust optimum, because
another Boolean spin could repair some of this witness loss.  The theorem
falsifies an inference about the relative certificate; it is not a signing
counterexample or a separated-optimum theorem.

## 3. What full `(G,R)` does and does not cost

For a possibly redundant pole table, `ker G subseteq ker D`, and the least
relative constant is

```math
\delta_{rel}
=\lambda_{max}\!\left(G^{\dagger/2}DG^{\dagger/2}\right)     \tag{CB.15}
```

on `range(G)`.  Thus full active-product `(G,R)` makes the certificate a
generalized eigenvalue computation, not a Boolean maximization.  At the
present logarithmic arity, however, the unrestricted table has
`Theta(q^2)=Theta(N)` real coherences.  This remains vastly smaller than the
entire `2^N` energy landscape but is much larger than the
`O(sqrt(N) log N)` exact row-histogram carrier.  Factorization or symmetry is
therefore what makes the relative certificate genuinely economical.

The repository's independent equal-Gram--Rayleigh Boolean collision also
shows that even full `(G,R)` does not reconstruct the Boolean response.  The
relative certificate is a deliberately stronger spectral sufficient
condition on one declared pole span, not a hidden exact optimizer.

## 4. A positive symmetry criterion

The obstruction is precisely off-diagonal coherence.  A group symmetry can
remove it without enumerating Boolean spins.

Let a finite abelian two-group `Gamma` act orthogonally on a pole span `U`.
Assume that `U` has an orthonormal basis `(e_s)_{s in I}` of pairwise
distinct real characters:

```math
\rho(g)e_s=\chi_s(g)e_s,
\qquad \chi_s(g)\in\{\pm1\}.                    \tag{CB.16}
```

Let `A=P_UT|_U` be a self-adjoint compression with `A preceq I`.  Its group
twirl is

```math
\overline A={1\over|\Gamma|}\sum_{g\in\Gamma}
 \rho(g)^*A\rho(g).                              \tag{CB.17}
```

Character orthogonality makes `overline A` diagonal in `(e_s)`, with
diagonal entries `A_ss`.

### Theorem CB.2 (twirled coherence criterion)

Suppose

```math
A_{ss}\ge1-d\quad\hbox{for every }s,
\qquad
\|A-\overline A\|_{op}\le\eta.                  \tag{CB.18}
```

Then

```math
\boxed{0\preceq I-A\preceq(d+\eta)I.}           \tag{CB.19}
```

Equivalently, in an orthogonal pole presentation,

```math
D\preceq(d+\eta)G.                               \tag{CB.20}
```

Exact equivariance of `A` under `Gamma` gives `eta=0`, so individual
Rayleigh deficits are then sufficient.

More quantitatively,

```math
\|A-\overline A\|_{op}
\le {1\over|\Gamma|}\sum_g\|[A,\rho(g)]\|_{op}. \tag{CB.21}
```

If `Gamma` is generated by `p` commuting involutions and every generator
commutator has norm at most `gamma`, then

```math
\|A-\overline A\|_{op}\le {p\gamma\over2}.       \tag{CB.22}
```

#### Proof

From (CB.18),

```math
A\succeq\overline A-\eta I\succeq(1-d-\eta)I.
```

Together with `A preceq I`, this proves (CB.19).  For (CB.21), note that

```math
\|A-\rho(g)^*A\rho(g)\|_{op}=\|[A,\rho(g)]\|_{op}
```

and average.  A group element represented by a word of Hamming length `ell`
has commutator norm at most `ell gamma` by telescoping.  The average word
length in the Boolean group is `p/2`, proving (CB.22).  `square`

Theorem CB.2 isolates a legitimate route around CB.1: prove exact or
approximate equivariance of the compressed child operator on the generated
pole algebra.  Without such a structural coherence statement, separate
near-top estimates do not approach the required relative certificate.

## 5. Research consequence

For a factored tensor construction, the factor presentation already supplies
the missing coherence and the relative defect is constant-size per factor.
For an arbitrary child with `p=Theta(log N)` generated ports, CB.1 shows that
small generator deficits, all small individual product deficits, and a
near-top trace do not suffice.  A near-original use of robust product
synchronization therefore needs a theorem of one of the following forms:

1. near-minimal sign children force approximate group equivariance on the
   declared pole span;
2. their product-defect operator has a factored or otherwise compressed
   presentation; or
3. the semantic selector orbit can be controlled directly without upgrading
   to the stronger all-vector inequality `D preceq delta G`.

The current result neither proves nor refutes any of these sign-specific
possibilities.
