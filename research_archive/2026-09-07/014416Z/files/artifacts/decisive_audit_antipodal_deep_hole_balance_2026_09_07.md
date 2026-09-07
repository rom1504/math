# Antipodal linearity alone does not balance deepest cosets

Date: 2026-09-07. Exact finite counterexample to a possible imported
coding principle. This is NOT a counterexample for complete-graph cut
codes or for exact signing minimizers.

Let C be the binary kernel of the six-by-twelve parity-check matrix
whose columns, encoded as integers in F2^6, are

`29,63,57,12,41,49,19,39,13,40,42,27`.

The columns are distinct and nonzero, so this is not a counterexample
created by repeated columns or an unconstrained coordinate. They span
F2^6, and their xor is j=47. In particular the all-one word is not in C.
Let `D=C union (1+C)` be the antipodal augmented linear code.

Distances to C depend only on syndrome s. They are the shortest-path
distances from zero in the Cayley graph with the twelve displayed
generators. Exact breadth-first search, independently checked by all
4096 physical words, gives

```math
\max_s\min\{d_C(s),d_C(s\mathbin{\mathrm{xor}}47)\}=2,
\qquad d_C(2)=2,\quad d_C(2\mathbin{\mathrm{xor}}47)=4.
```

Thus syndrome2 is a deepest coset of D, but its two C-halves have
distances differing by TWO, not at most one. The complete exact distance
table and a physical leader are printed by
`computations/decisive_audit_antipodal_deep_hole_balance_2026_09_07.py`.

For the original signing problem, the corresponding C is instead the
cut code of a COMPLETE GRAPH. Its distance formulas give
`P=m-2d(a,C)` and `R=m-2d(a,1+C)`. Consequently a proof that exact
signing minimizers have `|P-R|<=2` must use more than linearity,
antipodality, or simple parity-check columns. This counterexample does
not address the necessary extra complete-graph structure.
