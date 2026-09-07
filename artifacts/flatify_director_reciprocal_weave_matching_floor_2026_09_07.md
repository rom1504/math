# A matching floor for every reciprocal-column weave

Date: 2026-09-07. Status: root elementary proof; independent reconstruction
requested. This is a construction-class obstruction, not a universal
lower bound for full sign matrices.

## Theorem

Partition N vertices into m nonempty fibres of sizes k_i. For each fibre
let F_i be ANY k_i by m full sign matrix, and let s_ij=s_ji be arbitrary
signs for i!=j. Suppose the cross-fibre entries are

```math
C_{(i,a),(j,b)}=s_ij F_i(a,j) F_j(b,i),  i!=j.       (1)
```

The hollow full-sign matrices inside the fibres are completely arbitrary.
Then

```math
Q(C) >= max_{matchings P on [m]} sum_{ij in P} k_i k_j.   (2)
```

Neither orthogonality nor randomness nor a common seed is assumed.

## Proof, including the absolute-value polarity

Fix a matching P. For each matched pair {i,j}, initially choose
x_i=F_i(:,j) and x_j=s_ij F_j(:,i). This makes its cross energy k_i k_j.
Choose arbitrary spins on unmatched fibres. The total internal-fibre energy
T is unchanged if any entire fibre's spin is reversed.

For a common polarity sigma in {+1,-1}, replace x_j by sigma x_j in
each matched pair. Their total cross energy is sigma K, where
K=sum_P k_i k_j, while T stays the same. Choose sigma so sigma T>=0.
Now independently multiply all spins in each matched component, and each
unmatched singleton, by one unbiased common sign. Matching energies and T
are unchanged. Every cross term between different components has mean zero.
Therefore the expected sigma times the total energy is

```math
sigma T+K >= K.
```

Some deterministic assignment attains at least this expectation. It can
also be found by conditional expectation over the component signs. This
proves (2), with no assertion of cancellation between a bridge and its
own two children's internal energies.

## Uniform and nonuniform fibre consequences

For equal k_i=k,

```math
Q(C) >= floor(m/2) k^2.
```

In particular full retention k=m forces
Q(C)/N^(3/2)>=1/2-O(1/m). This includes arbitrary seed-dependent choices
of the row bases and arbitrary internal completions. Consequently changing
only the bases of a full-retention reciprocal weave cannot propagate
original near-minimizers whose normalized cap is eventually below .493609.
This does not forbid non-rank-one cross tiles or other global changes.

If k/m->p>0, the floor is sqrt(p)/2-o(1), consistent with the strict
upper construction at p=24/25. No contradiction with that construction
or its independently varying retained row sets occurs.

For even m, pair adjacent entries of the sorted sequence k_(1)<=...<=k_(m).
The corresponding matching weight is

```math
(1/2)sum_i k_i^2-(1/2)sum_j (k_(2j)-k_(2j-1))^2
 >= (1/2)sum_i k_i^2-(1/2)(k_max-k_min)^2.            (3)
```

Thus if m->infinity, k_max=O(m), and N is comparable to m^2, the leading
floor depends on the SECOND moment of fibre sizes, not just their mean.
Irregular retention cannot improve this floor at fixed total N and m:
sum k_i^2>=N^2/m, so (3) gives sqrt(N)/(2m)-o(1) after normalization.
For odd m, omit a smallest fibre and apply the same estimate, paying
O(k_max^2)=o(N^(3/2)).

## Research implication

This removes orthogonality, spectral flatness, and a common-seed assumption
from one precise weave obstruction. It does not show that arbitrary
flatification is false. To use this construction class for a putative
liminf c, uniform retention necessarily has p<=4c^2+o(1); that condition
is not asserted sufficient. The missing favorable-flatification theorem
remains open.
