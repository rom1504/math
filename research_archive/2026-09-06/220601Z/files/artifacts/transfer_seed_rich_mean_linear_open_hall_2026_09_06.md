# Fixed-stage rich mean-linear return: open k-branch Hall bound

2026-09-06. Finite-polynomial-stage lemma, independently reconstructed by
the adversarial agent. It generalizes the exact open-mark argument in
`continued_audit_next_return_open_mark_hall_2026_09_06.md`, Section 2,
from degree-at-most-three old primitives to a fixed degree bound M.
It does not by itself assert complete rich bounded-response closure.

## 1. Exact scope

Use the same literal marked-tree/forest representations, exact coherent
Boolean equality patterns and global-cut bounds as the cited proof.
Let Z_(P),j be a fixed fully injective homogeneous noise channel of
degree P>M. Its proper and fixed-slot cuts are uniformly at most

```
epsilon_n=C n^(-1/2) log(n+1)^C.
```

Let U_j be the exact degree-p disjoint k-branch main of the normalized
Hermite polynomial in Z_(P),j, where p=kP and k>=2 is fixed. Different
left noise branches have disjoint seed labels in this main. Let r_(c,d)
be an exact original-degree-d component, d>p, of a fixed polynomial in
the literal old source primitives at c, each of degree at most M.
All exact internal coherent equality patterns on the right are retained.

For bounded deterministic root diagonals D and D_a, and bounded-op B,
consider the probe E=B D U and the mean-linear return K=B D_a B r.
Then

```
(1/n) sum_i ||K_(E_i) star_p K_(K_i,d)||F^2 -> 0.      (1)
```

The mean coefficient in D_a is deterministic. This is not a claim for
an arbitrary random multiplier inside the last return.

## 2. Retain both source roots and every right free mark

For source roots j,c define the open tensor

```
J_(j,c;V)=K_(U_j) star_p K_(r_(c,d)), |V|=d-p.
```

Fix one of the finitely many exact contraction and coherent equality
patterns. Every left branch contributes P distinct contracted labels,
and different left branches contribute disjoint sets. For any subset
of k' left branches there are Pk' distinct matched labels. A right
primitive contains at most M distinct slots, including after coherent
label equalities. Hence the number of neighboring right primitives is
at least Pk'/M>k'. Hall matching supplies k distinct matched primitives.
The extra d-p right free marks cannot reduce this capacity bound.

Each assigned merge costs epsilon_n. If a proper nonempty portion of
the left branch is summed, use its proper cut. If a label stays open
or is retained at another right primitive, use its fixed-slot influence
cut. A complete unretained merge of one left branch into a single right
primitive is impossible because P>M; thus the equal-degree exceptional
primitive comparison from the original degree-three proof is absent.

The k charged merges have disjoint selected vertices, so their gains
multiply. The remaining merges are precisely the Hilbert-contractive
ones in the cited open-mark proof. Right open marks remain tensor
indices through the final Frobenius norm. Finite exact distinctness
inclusion-exclusion is handled there as well: repeats within a primitive
vanish; retained equality marks use the influence case. The prescribed-port
graph theorem supplies the bounded norms of the surviving exact coherent
diagrams. No coherent pre-Walsh tensor is declared close to its
squarefree projection. Consequently

```
sup_(j,c) ||J_(j,c;.)||F <= C epsilon_n^k,
||J||F <= C n epsilon_n^k=o(sqrt(n)).                 (2)
```

Only two gains are necessary. The full k-gain form is available under
the stated uniform cut rate; the odd Hermite feedback probe has k>=3.

## 3. Separate root transports, then take the common diagonal

Put L=B D_a B. The desired full contraction is exactly

```
sum_(j,c) B_ij D_jj L_ic J_(j,c;V).
```

Apply B D on the j axis and L on the c axis of J. Only afterwards
restrict the two resulting output roots to equality. This final
restriction is an orthogonal coordinate projection. Its Frobenius norm
is at most ||B D||op ||L||op ||J||F=o(sqrt(n)), proving (1).
Summing row absolute values first would lose this bound.

This statement concerns the exact uniformly controlled main kernels.
Any collision/raw-kernel errors must be removed in averaged L2 before
invoking (2), or passed directly through the tested contraction using
a uniformly bounded opposite kernel norm. Merely averaged small cuts
cannot be inserted into the supremum estimate (2). All approximation
stages, degree cutoffs and deterministic coefficient bounds are fixed
before n tends to infinity.

## 4. Use and remaining scope

The lemma closes the finite-stage mean-linear branch in

```
BC=Bc0(W)+B[(A(W)-a)Z+R_(>=2)]+B D_a B r,
a_i=E A(W_i),
```

relative to the original finite old primitive frame. It complements
the verified trig high-degree proper-cut theorem for Bc0 and the local
bounded coefficient/noise approximation theorem for the middle terms.
An enlarged literal variable R=BF-Z_P from a separate mixed-covariance
argument must not be silently added to the coherent frame here: its
original polynomial degree may exceed p. Restore the full old-noise
expansion relative to the original degree-M coherent frame, or explicitly
isolate its high-degree one-noise terms.

The complete ordered assembly of the rich literal returned-query law,
including source approximation, covariance comparison, all root caps
and a positive selected source block, is not claimed by this standalone
fixed-stage lemma.
