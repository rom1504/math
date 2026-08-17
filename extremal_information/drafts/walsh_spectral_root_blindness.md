# Walsh spectra are blind to the characteristic-root fibre

Status: rigorous task-local theorem, independently audited after the scope
repairs below. The result separates the exact rooted response state from all
spectra conditional on the tuple Gram/relation state, including the complete
synchronized-repetition spectral family.

## 1. Weyl reduction

Let `m>=1`, `V=F_2^m`, `q=2^m`, and let `F=R/q^(1/2)` be the normalized Walsh
involution on `R^V`.  Write

```math
(M_af)(x)=(-1)^(a dot x)f(x),
\qquad (T_af)(x)=f(x+a).                                 \tag{SB.1}
```

Then

```math
FM_aF=T_a,
\qquad T_bM_a=(-1)^(a dot b)M_aT_b.                     \tag{SB.2}
```

The linear-label child factor and the common edge factor are

```math
J_a=M_aFM_a=M_aT_aF,
\qquad J_0=F.                                           \tag{SB.3}
```

### Lemma SB.1 (trace of a Walsh word)

For every word `a_1,...,a_l in V`, there is a phase
`theta(a_1,...,a_l) in F_2`, given by a quadratic polynomial in the pairwise
inner products, such that, with `s=sum_j a_j`,

```math
J_(a_1)\cdots J_(a_l)
=(-1)^theta M_sT_sF^(l mod 2).                           \tag{SB.4}
```

Consequently

```math
\operatorname{tr}(J_(a_1)\cdots J_(a_l))
=\begin{cases}
(-1)^theta q,&l\text{ even and }s=0,\\
0,&\text{otherwise}.
\end{cases}                                             \tag{SB.5}
```

#### Proof

Induct using (SB.2). Multiplication by `J_a` toggles the final `F` and adds
`a` to both Weyl coordinates; every commutation contributes one bilinear
phase. This proves (SB.4) and the stated dependence of `theta`.

For an even word, `tr(M_sT_s)` vanishes unless `s=0`, when it is `q`. For an
odd word, the diagonal sum is

```math
q^(-1/2)\sum_x(-1)^(s dot x+(x+s) dot x),               \tag{SB.6}
```

up to the already recorded phase. The two `s dot x` terms cancel and
`x dot x=omega dot x` with `omega=(1,...,1) ne 0`, so the character sum is
zero. `square`

Only two data enter (SB.5): whether a specified binary combination of labels
vanishes, and the pairwise Gram values that determine the phase.

## 2. Root-blind spectrum theorem

For a graph `G` on `k` labeled vertices, let the normalized `kq by kq`
coefficient carrier be

```math
(K_(G,a))_(ii)=J_(a_i),
\qquad
(K_(G,a))_(ij)=w_(ij)J_0,                               \tag{SB.7}
```

where arbitrary real symmetric edge weights `w_ij` are allowed.

### Theorem SB.2 (Gram and relations determine every unrooted spectrum)

If two ordered label tuples have the same binary Gram matrix and relation
kernel, then for every weighted graph `G`,

```math
\operatorname{spec}K_(G,a)
=\operatorname{spec}K_(G,a')                            \tag{SB.8}
```

with multiplicity. The same holds for the full Walsh quadratic operator
`F tensor K_(G,a)` from Theorem 21.9. No characteristic-root fibre is needed.

#### Proof

Expand `tr(K^l)` over closed length-`l` block walks. Each summand is a scalar
product of edge weights times the trace of a word in the operators `J_0` and
`J_(a_i)`. By Lemma SB.1, the term is determined by:

1. whether the mod-two multiplicity vector of its nonzero labels lies in the
   relation kernel; and
2. a quadratic phase determined by their Gram matrix.

Thus all power traces agree for the two tuples. Finite real symmetric
matrices with the same power sums have the same eigenvalue multiset, proving
(SB.8). Tensoring with the fixed `F` preserves equality. `square`

### Corollary SB.3 (a scalable spectrum/response separation)

The singleton pairs in Theorem 21.13 have identical spectra under every
unrooted weighted Walsh graph made from synchronized copies of the compared
child, because the resulting constant tuples have identical Gram and
relation states. More generally the conclusion holds for any appended
context that leaves the **complete** tuple Gram/relation state matched. It
does not hold for an arbitrary common appended label: its cross pairing with
the two children can differ. Nevertheless one canonical rooted Boolean
future separates their projective responses by at least

```math
{1\over6}n^{3/2}.                                       \tag{SB.9}
```

Hence no collection of unrooted spectral moments, even of unbounded order
and across every synchronized-repetition weighted graph experiment, recovers
the rooted extremal state.

This does not prove that the characteristic-root fibre affects an unrooted
scalar Boolean maximum. It proves the sharper query-relative statement:
rooted response contains information absent from the complete unrooted
family of synchronized-repetition Walsh-graph spectra. The next semantic
question is whether unrooted Boolean graph maxima themselves collapse from
`(G,R,R_omega)` to `(G,R)`.
