# Complete bridge likelihood invertibly stores the projective latent law

**Status.** Rigorous finite theorem and scope audit.  At every fixed real
`t != 0`, the complete labelled bridge-likelihood table determines an
antipodally symmetric latent law exactly.  The only kernel on arbitrary
signed laws is the antipodally odd subspace.  For rank-one `m by n` words,
this means that the table stores the full probability law on
`2^(m+n-2)` projective atoms, not merely an overlap or Gram quotient.

This is an exact injectivity statement, not a stable recovery theorem.  The
inverse amplifies high Walsh levels by powers of `|tanh t|^(-1)`.  It
therefore rules out an **exact** uniform all-bridge compression, but it does
not rule out coarse uniform approximation or the inverse-escort-typical
approximation sought in the Ghirlanda--Guerra route.

## 1. Exact Walsh diagonalization and kernel

Let

```math
G_d=\{\pm1\}^d,
\qquad
\chi_S(q)=\prod_{i\in S}q_i
\quad(S\subseteq[d]),
```

and use normalized Walsh coefficients

```math
\widehat f(S)=2^{-d}\sum_{B\in G_d}f(B)\chi_S(B).
```

For a finite signed measure `nu` on `G_d`, put

```math
m_\nu(S)=\int\chi_S(q)\,d\nu(q)
```

and define its complete cosh bridge transform by

```math
(\mathcal K_t\nu)(B)
=\int\cosh\{t\langle B,q\rangle\}\,d\nu(q).
                                                               \tag{BL.1}
```

Write `c=cosh t`, `s=sinh t`, and `rho=tanh t`.

**Theorem BL.1 (bridge-likelihood inversion).**  For every real `t`,

```math
\boxed{
\widehat{\mathcal K_t\nu}(S)
=\mathbf1_{\{|S|\ {\rm even}\}}
 c^{d-|S|}s^{|S|}m_\nu(S).}                       \tag{BL.2}
```

If `t != 0`, then

```math
\boxed{
\ker\mathcal K_t
=\{\nu:\nu(-q)=-\nu(q)\text{ for every }q\}.}    \tag{BL.3}
```

Consequently two arbitrary probability laws have the same complete table
if and only if their antipodal symmetrizations agree.  If `mu` is
antipodally symmetric, `mu(q)=mu(-q)`, then it is recovered by

```math
\boxed{
\mu(q)=2^{-d}\sum_{|S|\ {\rm even}}
 {\widehat p_\mu(S)\over c^{d-|S|}s^{|S|}}\chi_S(q),
\qquad
p_\mu=\mathcal K_t\mu.}                            \tag{BL.4}
```

Equivalently, for the unit-mean normalized likelihood

```math
P_\mu(B)={p_\mu(B)\over(\cosh t)^d},
\qquad E_UP_\mu=1,
```

one has

```math
\boxed{
\widehat P_\mu(S)
=\mathbf1_{\{|S|\ {\rm even}\}}\rho^{|S|}m_\mu(S),
\qquad
\mu(q)=2^{-d}\sum_{|S|\ {\rm even}}
 \rho^{-|S|}\widehat P_\mu(S)\chi_S(q).}          \tag{BL.5}
```

These conclusions remain valid when `mu` is supported on any antipodal
subset of `G_d`: extend it by zero to the full cube and apply (BL.4).

*Proof.*  For `z in G_d`, product expansion gives

```math
\begin{aligned}
\cosh\left(t\sum_i z_i\right)
&={1\over2}\left\{\prod_i(c+sz_i)+\prod_i(c-sz_i)\right\}\\
&=\sum_{|S|\ {\rm even}}c^{d-|S|}s^{|S|}\chi_S(z).
                                                               \tag{BL.6}
\end{aligned}
```

Substitute `z=Bq`, use
`chi_S(Bq)=chi_S(B)chi_S(q)`, and take the normalized Walsh
coefficient in `B`; this proves (BL.2).  For real `t != 0`, every displayed
even multiplier is nonzero.

Let `iota(q)=-q` and
`nu^+=(nu+iota_#nu)/2`.  Its odd Walsh moments vanish, while its even
Walsh moments equal those of `nu`.  Thus (BL.2) says that
`mathcal K_tnu=0` exactly when every Walsh coefficient of `nu^+` is zero.
Walsh inversion gives `nu^+=0`, equivalently
`iota_#nu=-nu`.  This proves (BL.3).  An antipodally symmetric `mu` has
`m_mu(S)=0` at every odd level, so ordinary Walsh inversion of its mass
function gives (BL.4), and division by `c^d` gives (BL.5). `square`

The symmetry hypothesis is essential.  If only the **support** is
antipodal, while the law need not satisfy `mu(q)=mu(-q)`, then the theorem
with “determines `mu`” is false: the cosh table sees only

```math
\mu^{+}(q)={\mu(q)+\mu(-q)\over2}.                 \tag{BL.7}
```

At `t=0` all nonconstant multipliers vanish as well, which is why that value
is excluded.

## 2. Exact rank-one projective dimension

Let

```math
\mathcal R_{m,n}
=\{xy^{\mathsf T}:x\in\{\pm1\}^m,
                   y\in\{\pm1\}^n\}
\subseteq G_{mn}.                                  \tag{BL.8}
```

This is a subgroup under entrywise multiplication.  The parametrization
`(x,y) -> xy^T` has the two-element kernel
`{(mathbf1,mathbf1),(-mathbf1,-mathbf1)}`, and global matrix negation acts
freely.  Hence

```math
|\mathcal R_{m,n}|=2^{m+n-1},
\qquad
|\mathcal R_{m,n}/\{\pm1\}|=2^{m+n-2}.             \tag{BL.9}
```

An antipodally symmetric rank-one law is therefore an arbitrary probability
law on `2^(m+n-2)` projective atoms.  Its affine simplex dimension is

```math
\boxed{2^{m+n-2}-1.}                               \tag{BL.10}
```

The ambient inversion (BL.5) contains redundant characters on this
rank-one subgroup.  A minimal explicit family is obtained as follows.  Fix
a spanning tree `T` of `K_(m,n)`.  For an edge set `S`,

```math
\chi_S(xy^{\mathsf T})
=\prod_i x_i^{\deg_S(i)}\prod_j y_j^{\deg_S(j)}.  \tag{BL.11}
```

Two even subsets of `T` give the same character on the projective rank-one
group only if their symmetric difference is Eulerian.  A tree has no
nonempty Eulerian edge set, so the `2^(m+n-2)` even subsets of `T` give all
distinct projective characters and hence the complete dual group.  If
`bar mu([q])=mu(q)+mu(-q)` is the projective law, quotient Fourier inversion
is

```math
\boxed{
\bar\mu([q])
=2^{-(m+n-2)}
 \sum_{\substack{S\subseteq T\\|S|\ {\rm even}}}
 \rho^{-|S|}\widehat P_\mu(S)\chi_S(q).}          \tag{BL.12}
```

Thus the complete bridge table reconstructs the entire rank-one projective
law, although only `2^(m+n-2)` suitably chosen Walsh channels are needed.
For comparison, if two even edge sets in the full bipartite graph have the
same vertex-degree boundary, their characters agree on every rank-one
word; their difference is Eulerian.  This is the redundancy behind the
mandatory Eulerian coefficients in the existing cycle-code audit.

## 3. Complete replica Gram data still miss the labelled likelihood

Injectivity also gives a small explicit obstruction to any coordinate lift
from replica Gram data alone.  In `G_5`, set

```math
\begin{array}{c|cccc}
 &q_0&q_1&q_2&q_3\\ \hline
\mathcal A&+++++&--+++&-+-++&+--++
\end{array}
```

and replace the last word in the second family by

```math
q'_3=-++-+.
```

Let `mu_A` and `mu_B` be uniform on the eight-point antipodal sets generated
by their respective four displayed words.

**Proposition BL.2 (finite Gram-only lifting falsifier).**  The complete
infinite replica Gram arrays under `mu_A` and `mu_B` have the same law, but,
at the common labelled bridge `B=q_3`,

```math
\begin{aligned}
p_{\mu_A}(q_3)
&={\cosh(5t)+3\cosh t\over4},\\
p_{\mu_B}(q_3)
&={\cosh(3t)+3\cosh t\over4}.                    \tag{BL.13}
\end{aligned}
```

The two values are different for every real `t != 0`.

*Proof.*  In each four-word family every word has squared norm five and
every two distinct displayed words have inner product one.  Thus both
labelled base Gram matrices have diagonal five and off-diagonal one.  An
i.i.d. replica from either antipodal law can be written
`Q^ell=sigma_ell q_(I_ell)`, where `I_ell` is uniform on four labels and
`sigma_ell` is an independent fair sign.  Every replica inner product is
therefore

```math
Q^\ell\mathbin\cdot Q^k
=\sigma_\ell\sigma_k
 \{5\mathbf1_{I_\ell=I_k}+\mathbf1_{I_\ell\ne I_k}\},
                                                               \tag{BL.14}
```

which proves equality of the complete Gram-array laws, not just equality of
one-overlap marginals.  In family `A`, the bridge `q_3` has inner products
`1,1,1,5` with the four centers.  It has inner products `1,1,1,-3` with the
four centers of family `B`, since `q_3 dot q'_3=-3`.  Cosh removes the fair
antipodal signs and gives (BL.13).  Finally
`cosh(5t)>cosh(3t)` for real `t != 0`. `square`

This is already a rank-one example by viewing the words as `1 by 5` sign
matrices.  It isolates the exact issue in the GG proposal: replica Gram data
remember the abstract Hilbert-space configuration, while a labelled bridge
queries its particular embedding in the coordinate cube.

There is a tempting but invalid diffuse amplification of this example.
Repeating every coordinate, adding fixed-crossover BSC noise, and pairing
with a similarly smoothed left factor does make the two finite-replica
overlap arrays converge to the same limit.  It also gives the factor laws
the every-subset conditional atom bound
`(1-epsilon)^|U|`.  It does **not**, however, retain the proposed leading
bridge-likelihood separation at the physical scale `t=beta/sqrt(N)`.
Fixed BSC noise has full support.  For a rank-one bridge `B=uv^T`, the
exact aligned latent word `Q=B` has probability `exp{-O(N)}`, while its
Boltzmann gain is

```math
\exp\{t\langle B,Q\rangle\}
=\exp\{\Theta(N^{3/2})\}.                          \tag{BL.14a}
```

This entropy cost is lower order, so both smoothed families have the
universal leading support maximum `t mn+O(N)`, rather than leading slopes
set by their typical alignments `1` and `3/5`.  More explicitly, fixed BSC
densities on two centers have pointwise likelihood ratios at most
`exp{O(N)}`; averaging the common positive kernel preserves that bound, so
the two log bridge likelihoods can differ by at most `O(N)`, not by the
suggested `Theta(N^(3/2))`.  A scalable diffuse Gram-only no-go would need a
different query scale or non-full support.  Proposition BL.2 is only the
finite exact obstruction asserted here.

## 4. Exact cavity responses also reconstruct the likelihood

For symmetric `mu`, (BL.5) also has the binary-channel form

```math
P_\mu(B)
=E_\mu\prod_{f=1}^d(1+\rho B_fQ_f).               \tag{BL.15}
```

Fix a bridge coordinate `e`, average out `B_e`, and define

```math
\begin{aligned}
P_{-e}(B_{-e})
&=E_\mu\prod_{f\ne e}(1+\rho B_fQ_f),\\
r_e(B_{-e})
&={E_\mu[Q_e\prod_{f\ne e}(1+\rho B_fQ_f)]
       \over P_{-e}(B_{-e})}.
\end{aligned}                                      \tag{BL.16}
```

Since `|rho|<1` at finite real `t`, all denominators are positive, and

```math
\boxed{
P_\mu(B)=P_{-e}(B_{-e})\{1+\rho B_er_e(B_{-e})\},
\qquad
{P_\mu(B_e=+,B_{-e})\over P_\mu(B_e=-,B_{-e})}
={1+\rho r_e(B_{-e})\over1-\rho r_e(B_{-e})}.}     \tag{BL.17}
```

Therefore the complete labelled table of deleted-edge cavity responses
determines all neighboring likelihood ratios on the bridge cube.  Integrate
the ratios along cube paths and impose `E_UP_mu=1`; this determines
`P_mu`, and then BL.1 determines `mu`.

The same conclusion holds if a lifting returns the full-posterior
coordinate means

```math
m_e(B)=E[Q_e\mid B].
```

Indeed one-edge Bayes algebra gives

```math
m_e(B)={r_e(B_{-e})+\rho B_e
        \over1+\rho B_er_e(B_{-e})},
\qquad
r_e(B_{-e})={m_e(B)-\rho B_e
        \over1-\rho B_em_e(B)}.                   \tag{BL.18}
```

Hence an exact coordinate-response lift that is uniform over all labelled
bridges is already an exact encoding of the complete projective latent law.

## 5. Collision with existing results

There is no contradiction with the repository, but there is a strong
algebraic collision.

1. Equations (CY.1)--(CY.2) in
   [`actual_child_cycle_code_overlap_floor.md`](actual_child_cycle_code_overlap_floor.md)
   already prove, for the actual rank-one channel,

   ```math
   \widehat P(S)=\rho^{|S|}E\chi_S(Q).
   ```

   That audit uses only the child-independent Eulerian moments.  Theorem
   BL.1 is the full converse: retain all even coefficients and apply
   projective Walsh inversion.  Thus the present theorem is essentially the
   unrecorded inversion corollary of CY.1--CY.2.

2. Theorem 21.50 in [`theorems.md`](../theorems.md) proves the same
   projective-convolution principle for the different kernel
   `z -> |sum_i z_i|`: all of its even Walsh multipliers are nonzero, so its
   complete labelled response recovers a projective row histogram.  Here
   the kernel is `cosh(t sum_i z_i)` and the multipliers are the simpler
   powers in (BL.2).

3. Theorem 37.60 in [`theorems.md`](../theorems.md) recovers a complete
   augmented child law from all values on a basis edge-flip cube, through an
   invertible subset-triangular transform followed by Walsh inversion.  It
   has the same “complete table is not a compression” conclusion, but uses
   different queries and a different latent alphabet.

The theorem therefore sharpens, rather than conflicts with, the missing
coordinate-response lifting boundary in
[`actual_child_generic_gg_perturbation_mapping.md`](actual_child_generic_gg_perturbation_mapping.md).

## 6. Exact implication and ill-conditioning caveat

A replica overlap/Gram object can be a strict exact sufficient statistic for
**all** labelled bridge likelihoods or cavity responses only if, together
with any supplied coordinate data, it is injective on the full projective
latent law.  For unrestricted symmetric rank-one laws that means retaining
the `2^(m+n-2)-1`-dimensional simplex in (BL.10).  Thus a proposed exact
uniform all-bridge lift is not a nontrivial latent-law compression.

This statement must not be promoted to a stable lower bound.  For the
normalized operator, an even level-`k` mode has eigenvalue `rho^k`.  A
uniform additive likelihood error `epsilon` can create Walsh-coefficient
error as large as `epsilon`, and inversion then amplifies that channel by
`|rho|^(-k)`.  On the full cube the largest needed even degree is

```math
d_{\rm even}=2\lfloor d/2\rfloor.
```

On the rank-one subgroup, the spanning-tree inversion improves this to

```math
k_{\rm tree}=2\left\lfloor{m+n-1\over2}\right\rfloor,
```

but the worst displayed amplification is still
`|rho|^(-k_tree)`.  At the physical scale
`t=beta/sqrt(N)` with `m+n=N`, this is
`exp{Theta(N log N)}`.  Exact injectivity consequently gives no useful
recovery from an `o(N)` log-response error or any other coarse error metric.

Moreover every finite inverse escort has full bridge support, so exact
equality “almost surely under the escort” is still equality on the complete
cube.  A genuinely compressive GG statement must instead be approximate in
an escort-weighted norm, in probability, or at the aggregate `o(N)` target
scale.  Such inverse-escort-typical approximation remains possible.  Coarse
**uniform** approximation is not logically excluded either, precisely
because the inverse is ill-conditioned; what BL.1 excludes is exact uniform
reconstruction, or approximation accurate enough to resolve every attenuated
Walsh channel.
