# A fixed-seed, cycle-sensitive weave at the old strict upper bound

2026-09-07. Status: **proved and independently audited**. The new points
are a mixed-group orbital estimate, its fixed-macro-signing compiler, and
an explicit cycle statistic proving that this compiler does not erase its
seed in law. This is **not** a numerical seed-value transfer or a proof of
convergence. All final column signs are absent: the distinction between
input and output randomization is essential.

Write `G_m` for the signed-permutation group and `P_m` for the ordinary
permutation group. For either group `J`, put

```math
P_{t,J}(v)=\mathbb E_{g\in J}e^{-t\|v-gv\|^2},
\qquad L_{t,J}(v)=P_{t,J}(v)^{1/2}.                 \tag{1}
```

The archived row theory uses `L_{t,G_m}`. The fixed-seed contraction below
instead produces `L_{t,P_m}` on the signed, not absolute, spectrum.

## 1. Mixed input and output groups

**Lemma.** Fix `t,C>0`, and put
`D=ceil(e^2(2tC+1)m)`. For every orthogonal `U` and every
`v` with `||v||^2<=Cm`,

```math
\mathbb E_{g\in G_m}L_{t,P_m}(Ugv)
\le \sqrt2\exp\left\{\frac\pi2\sqrt{2D/3}\right\}
        L_{t,G_m}(v).                              \tag{2}
```

Thus changing only the OUTPUT orbital group costs `exp(O_{t,C}(sqrt m))`,
not an exponential-in-`m` loss. The input still has a fresh signed
permutation. A random `U` is allowed conditionally on `v` if it is
independent of that fresh input permutation.

### Proof

Use the Gaussian Fock feature

```math
\phi_t(v)=e^{-t\|v\|^2}
 \bigoplus_{d\ge0}\frac{(2t)^{d/2}}{\sqrt{d!}}v^{\otimes d},
\qquad \langle\phi_t(v),\phi_t(w)\rangle
 =e^{-t\|v-w\|^2}.                                \tag{3}
```

Let `Q_P` be the projection to ordinary-permutation invariants and
`Pi_D` the degree cutoff. The signed-input orbit covariance

```math
C_v=\mathbb E_{g\in G_m}
 |\phi_t(gv)\rangle\langle\phi_t(gv)|
```

has operator norm exactly `P_{t,G_m}(v)`: its nonzero eigenvalues are
those of the finite orbit Gram matrix, which has nonnegative entries
and constant row sum `P_{t,G_m}(v)`. This step concerns the INPUT group
only. The output invariant space in degree `d` has dimension at most
the partition number `p(d)`, since its polynomials are symmetric in
the coordinates themselves, with no even-degree restriction. Therefore

```math
\operatorname{rank}(Q_P\Pi_D)
 \le\sum_{d\le D}p(d)\le e^{\pi\sqrt{2D/3}}.        \tag{4}
```

The last bound follows by applying the partition generating function
at `e^{-s}` and minimizing `sD+pi^2/(6s)`.

Orthogonal transformations preserve degree and act unitarily, so

```math
\mathbb E_g\|Q_P\phi_t(Ugv)\|^2
\le e^{\pi\sqrt{2D/3}}P_{t,G_m}(v)
 +\Pr\{\operatorname{Poisson}(2t\|v\|^2)>D\}.       \tag{5}
```

The same Poisson tail estimate as in the archived orbit proof makes
the last term at most `exp(-(4tC+1)m)`. Meanwhile
`P_{t,G_m}(v)>=exp(-4tCm)`. Taking a square root and using Jensen
proves (2). No output sign average has been inserted.

Two exact ordinary-permutation inequalities will also be used:

```math
L_{t,P_m}(v_+,v_-)
 \le L_{t,P_{m_+}}(v_+)L_{t,P_{m_-}}(v_-),
\qquad
\max_j L_{t,P_{m-1}}(v\setminus v_j)
 \le\sqrt m\,L_{t,P_m}(v).                         \tag{6}
```

The first follows because the invariant subspace for all permutations
is contained in that for within-block permutations. For the second,
the permanent terms fixing coordinate `j` have diagonal kernel value
`exp(-t(v_j-v_j)^2)=1`, and contribute `P(v\setminus v_j)/m`.

## 2. The old signed-input recursion gives the same exponent

Take EXACTLY the recursive orthogonal matrices from the archived
`continued_convergence_recursive_orbit_bound_2026_09_06.md`:

```math
U_s=\operatorname{diag}(U_q^{(1)},U_q^{(2)})
 \frac1{\sqrt2}\begin{pmatrix}I&I\\I&-I\end{pmatrix}g,
\quad s=2q,
\qquad U_q^{\rm terminal}=U_q^0g.                  \tag{7}
```

Every `g` is a fresh independent SIGNED INPUT permutation. Both child
bases are independent. Terminal matrices `U_q^0` are fixed normalized
Hadamard matrices. There is **no** signed permutation on the output
of `U_s`. Later an independent ordinary output permutation is allowed.

Set `H_s=sqrt(s) U_s^T`. Hence the top `g` acts as a PHYSICAL ROW gauge
and permutation of `H_s`, not as a final column sign pattern. Internal
input permutations occur before higher Hadamard mixing. They must not
be silently replaced by arbitrary final column signs.

For fixed depth `r`, finite initial alphabet and symmetrized initial law
`nu`, (2), the first inequality in (6), and the very same signed-pair
type calculation as in the archive give

```math
\limsup_{m\to\infty}\frac1m
 \log\mathbb E L_{t,P_m}(U_m v_m)
 \le (\mathcal B^r\Phi_t)(\nu),                    \tag{8}
```

uniformly over terminal Hadamards, whenever the symmetrized empirical
type of `v_m` tends to `nu` and `||v_m||^2=O(m)`.

For completeness, the induction uses ordinary-permutation output norms
at EVERY node. Splitting the output gives the product of the two child
ordinary-permutation norms by (6). The fresh signed input permutation
gives exactly the archived pair-table probability
`exp[-s D(pi||nu\otimes nu)/2+o(s)]`. At each terminal, (2) bounds the
ordinary output norm by the signed INPUT norm, whose finite-type exponent
is `Phi_t`. Thus no unsigned-input Bellman problem, non-symmetric child
supersolution, or implicit output sign symmetry is being assumed. The
changed constants in (2) add only `O_{r,t,C}(sqrt m)` over all terminals.

## 3. Fixed macro signs: reflection-twisted graph contraction

Fix ANY symmetric sign matrix `S` of order `m`, including arbitrary
diagonal signs. It may be an exact original optimizer, a low-cap signing,
or the all-positive signing; no hypothesis on its cap is used. Choose
independent bases `H_i` from (7), and independent ordinary column
permutations. Fix any selectors `T_i` of common cardinality `k`.
Form the rank-one-block weave

```math
W_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i),
\qquad a\in T_i,\ b\in T_j.                       \tag{9}
```

This is an actual symmetric full sign matrix of order `N=mk`. For spins
`x_i`, put `h_i=H_i[T_i,:]^T x_i`, and

```math
D_\sigma=\sum_{i,j}(h_i(j)-\sigma S_{ij}h_j(i))^2
 =2(m^2k-\sigma x^TWx).                            \tag{10}
```

Condition on each randomly assigned diagonal coordinate and drop the
nonnegative scalar-coordinate terms in `D_sigma`. The remaining edge
kernel on the SIGNED normalized spectra `v_i=h_i/sqrt k` is

```math
K_{\sigma S_{ij}}(a,b)=e^{-t(a-\sigma S_{ij}b)^2}
 =\langle\phi_t(a),R_{\sigma S_{ij}}\phi_t(b)\rangle,
                                                               \tag{11}
```

where `R_-` is the unitary reflection acting by `(-1)^d` in Fock degree
`d`, and `R_+=I`. Thus each edge is a contraction of norm one between
Gaussian feature spaces; it need not itself be a positive-semidefinite
two-variable kernel. Absorb its reflection on either endpoint tensor.
Graph Cauchy--Schwarz/Finner bounds the contraction by the product of
vertex tensor norms. The squared norm of a vertex's ordinary-permutation
average is exactly `P_{t,P_{m-1}}` of its remaining SIGNED spectrum.
This proves, for every fixed `S`,

```math
\Pr\{\exists x,\sigma:D_\sigma\le2\gamma m^2k\}
\le 2e^{t\gamma m^2}\prod_i Z_i^{P}(t),
\quad
Z_i^{P}(t)=\sum_{x_i}
 \max_j L_{t,P_{m-1}}(h_i/\sqrt k\setminus h_i(j)/\sqrt k).
                                                               \tag{12}
```

Here probability in (12) first refers to the independent ordinary column
permutations with the bases fixed. Averaging the independent recursive
bases and using (6),(8) gives, when `k/m->p`,

```math
\Pr\{\exists x,\sigma:D_\sigma\le2\gamma m^2k\}
\le \exp\{m^2[t\gamma+p\log2+
             (\mathcal B^r\Phi_t)(\nu_p)+o(1)]\}.   \tag{13}
```

The error is uniform in `S`. No macro signs are averaged in (11)--(13).

In particular EVERY old strict finite-depth row certificate is valid for
EVERY fixed macro seed under this different ensemble. If its limiting
upper coefficient is `U`, then for every fixed `epsilon>0` one can choose
a finite depth and obtain

```math
\Pr\{Q(W-\operatorname{diag}W)
       >(U+\epsilon)N^{3/2}\}\le e^{-c_\epsilon m^2}
                                                               \tag{14}
```

for all large realizable `m`, uniformly over all `S` of order `m`.
The diagonal payment is only `N/2`. The old exact ternary certificate
`flatify_independent_2026_09_07_ternary_upper_proof.md` therefore applies
with

```math
U=\frac{97/20+(24/25)\log2-5151/6250}
         {(97/10)\sqrt{24/25}}
 =0.49360809358874865\ldots .                        \tag{15}
```

Terminal order density and principal restriction give the old all-order
numerical bound. Such a final restriction need not retain the full macro
seed interface; the fixed-seed statement itself is at the actual parent
orders `N=m floor(pm)`.

## 4. A nonzero cycle statistic rules out hidden seed erasure

The following calculation uses terminal PHYSICAL Hadamards
`H_q^0=sqrt(q)(U_q^0)^T` normalized so that one row is all positive.
One may simply dephase both the first row and column. Every other
physical row is then balanced. It does not assume
Walsh algebra, so it is compatible with the arbitrary-terminal theorem.

Let `b=2^r`, `q=m/b>=4`. Before the final ordinary column permutation, a
fixed physical row of the recursive `H_m` is a concatenation of `b`
independent signed-uniform rows from the fixed terminal Hadamards. This
follows by induction from (7): each child has an independent fresh input
signed permutation, so any queried child physical row is signed-uniform,
independently of the queried index. At a terminal it is immediate.

For a signed-uniform row of a normalized order-`q` Hadamard, an average
product over four distinct coordinates equals `1/(q-3)`. Indeed the
all-positive row contributes `1/q`; each balanced row contributes
`3/[(q-1)(q-3)]`, by the coefficient of `z^4` in
`(1-z^2)^{q/2}` divided by `binom(q,4)`. Their average is `1/(q-3)`.
Distinct-coordinate second moments are zero, and odd moments are zero
because of the independent row sign. Consequently only four coordinates
from a single terminal block contribute after the final permutation.
For any four distinct macro coordinates the exact moment is

```math
\mu_4(m,r)
 =\frac{b(q)_4}{(m)_4}\frac1{q-3}
 =\frac{bq(q-1)(q-2)}{m(m-1)(m-2)(m-3)}>0.          \tag{16}
```

Fix five distinct macro fibres `J` and one fixed retained physical vertex
in each. Their ten cross-fibre entries form a `K_5`. Independence between
the five fibre bases and (16) give the EXACT identity

```math
\mathbb E\prod_{\{i,j\}\subset J}W_{(i,a_i),(j,a_j)}
 =\mu_4(m,r)^5\prod_{\{i,j\}\subset J}S_{ij}.       \tag{17}
```

In contrast, independent final column signs would make the ten effective
macro signs independent and force this moment to zero. Thus (17) audits
the absence of precisely the gauge randomization that would erase `S`.

At fixed depth, `mu_4^5=Theta_r(m^{-5})`. Conditioning on the successful
cap event in (14) changes the expectation of this bounded statistic by
at most `2e^{-c_epsilon m^2}`. For all sufficiently large `m`, the sign
of the conditioned expectation is therefore still the macro `K_5`
product. There are actual low-cap output laws which remain different
for seeds having different `K_5` products.

For `m>=7`, all macro `K_5` products determine the OFF-DIAGONAL entries
of `S` modulo switching and global negation; its arbitrary diagonal is
not encoded by this statistic. To see this, choose disjoint vertices `a,b,c,d` and a
three-element set `R`. The symmetric difference of the four `K_5`
edge sets on `Rac,Rbc,Rad,Rbd` is the four-cycle `ac,bc,ad,bd`.
Thus equality of all `K_5` products implies equality of all four-cycle
products. The ratio signing then has all triangle products equal, and
is of the form `epsilon xi_i xi_j`. The converse is immediate because
each `K_5` has even degree and ten edges.

This is a LAW distinction. It does not assert that one sampled output
permits efficient or exact reconstruction of an unknown seed.

## 5. Exact scope and next obligation

The theorem applies to **every fixed macro signing**, hence to any
selectable exact optimizer, and constructs original full sign matrices
at the old strict bound while keeping a nontrivial cycle-space dependence.
It removes macro-edge independence as a necessary ingredient of the old
row certificate. Its proof is not just the observation that a fixed `S`
can be absorbed into independent final column signs: those signs have
been excluded, and (17) would be impossible under that erased law.

The output bound is still `U+epsilon`, regardless of `Q(S)/m^(3/2)`.
In particular this does NOT prove landing at an unusually favorable
optimal value, an asymptotically lossless multiplier, width-to-cap
equivalence, convergence, or nonconvergence. The cycle statistic here
has only polynomially small expectation; it does not itself furnish
a leading `m^2` pressure gain depending on the seed.

Nor may the director's new independent-tilted-edge stability penalty
be imported automatically: the present final macro signs are fixed.
Obtaining a seed-dependent improvement from the unaveraged, reflection-
twisted contraction is a concrete next step; replacing it by its product
of vertex norms is exactly where the displayed certificate loses `S`.

## 6. Finite regression record

`computations/principle_synthesis_2026_09_07_fixed_seed_mixed_orbit_check.py`
and its namesake results JSON preserve five exact fourth-moment checks
(4440 terminal row words in total), 32 reflection-twisted triangle
contractions with positive Gaussian kernels, and 231 exact symmetric-
difference identities generating four-cycles from four `K_5` edge sets.
All PASS. These are regressions of the finite identities, not substitutes
for the asymptotic Fock/type proof.
