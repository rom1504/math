# Flat transport of two or more globally controlled rooted tensors

Date: 2026-09-06. Director theorem, independently submitted to both feedback
and audit researchers. This is stronger than the earlier transported-power
lemma: the input tensors need not have small fixed-root cuts.

## 1. Statement

For a=1,...,p, p>=2 fixed, let K_j^a be an order-q_a tensor, q_a>=1,
indexed by a common root j=1,...,n. Assume EVERY global proper flattening
of K^a, treating the root as one additional tensor slot, has operator norm
at most C_a. The root-map cut j versus all marked slots is included.
The constants may grow polylogarithmically with n.

Let b_ij be arbitrary scalars and define

```math
L_i=\sum_j b_{ij}\,K_j^1\otimes\cdots\otimes K_j^p.
```

For every proper cut of the marked slots of this fixed-root tensor,

```math
\boxed{\|L_i^{\rm cut}\|_{op}
       \le \max_j|b_{ij}|\prod_{a=1}^p C_a.} \tag{1}
```

No bound on the input FIXED-root proper cuts is assumed. Thus a flat
signing transport has output proper cuts O(n^(-1/2) polylog n), even
when some coherent input cuts have order-one norm.

## 2. Orient the common root and factor the flattening

Fix the desired left/right cut of all marked slots. For each branch,
assign its root j to one side so that the branch's resulting global cut
is proper. A branch whose marked slots are all left must have its root
on the right; an all-right branch must have its root on the left. A
straddling branch permits either assignment.

The assignments can always be made with at least one root-left and at
least one root-right branch. If there are whole branches on both sides,
this is forced already. If there are whole branches on just one side,
properness of the total cut gives a straddler or a whole opposite branch,
which provides the other orientation. If every branch straddles, use the
two orientations on two different branches. Here p>=2 is essential.

For the root-left group, tensor the individual bounded flattenings and
compress their several output root indices to one common j. This is a
partial isometry, so the resulting operator has norm at most the product
of the group's C_a. Denote it

```math
U:\mathcal R_+\longrightarrow\mathbb R^n\otimes\mathcal L_+.
```

For the root-right group, duplicate the input root into the several
equal-root coordinates (an isometry), then apply the tensor product of
the branch operators. The result is

```math
V:\mathbb R^n\otimes\mathcal R_-\longrightarrow\mathcal L_-.
```

Up to coordinate permutations, the desired flattening is EXACTLY

```math
(I_{\mathcal L_+}\otimes V)
\,(D_{b_i}\otimes I_{\mathcal L_+\otimes\mathcal R_-})
\,(U\otimes I_{\mathcal R_-}). \tag{2}
```

The common root is output by U and consumed as input by V. Its summation
is an ordinary operator composition, not a partial trace, and introduces
no dimension factor. Multiplication by D_b has norm max|b_ij|. This proves
(1), including mixed branch orders and arbitrary nontrivial cuts.

## 3. Exact distinctness and the important collision boundary

Symmetrizing finitely many slot patterns changes only fixed constants.
Projection onto pairwise distinct marked labels preserves any local cut
norm up to a factor at most 2^(q(q-1)/2), q=sum q_a: a cross-cut equality
is a rectangular block pinching, its exclusion has norm at most two,
and same-side exclusions are row or column projections.

Consequently (1) applies to the EXACT top squarefree source component.
It does not say that this component approximates the entire Boolean
product in L². Between-branch label collisions can carry leading lower
Walsh components when the coherent inputs have high individual influence.
Those components must be retained or separately bounded. In particular,
this theorem does not permit replacing (QS)^3 by an unrestricted Gaussian
cubic source when Q=I.

## 4. Why one factor is different

The exclusion p=1 is real. Let B be a normalized symmetric conference
matrix. With three independent sign colors, set D_i=S_i^1 S_i^2 S_i^3
and K=BD. The global tensor of K has all proper global cuts bounded by
one, and its fixed-root cuts are O(n^(-1/2)). But BK=D because B²=I.
The returned fixed-root tensor has proper-cut norm one. An additional
nonlinear factor, not merely a linear transport, is what creates (1).

## 5. Application being checked, not silently assumed

The newly audited higher-Boolean-derivative theorem provides the global
cut hypotheses for exact Walsh components of centered fixed-depth
polynomial circuits. For an even centered coherent coefficient A_i(W_i)
and an odd noise component Z_i of degree at least three, the transported
top product B[A Z] therefore has small proper cuts and degree at least
five. The first marked-history primitive coherent fields have degree at
most three, so the local first-merge contraction lemma can then test this
top product against their fixed polynomial responses.

To conclude the centered-coefficient ENERGY estimate, the discarded
noise-touching product collisions must still be charged. A full noise
contraction into an aggregated high-degree A is not automatically proper;
one must expand A into its primitive coherent factors and use the
first-merge estimate there. The present tensor theorem alone does not
assert that final energy conclusion.
