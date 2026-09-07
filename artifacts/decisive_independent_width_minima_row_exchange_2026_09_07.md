# Actual width minima: complete order-nine classification and row exchange

## Scope and archive audit

Write P(A)=max q_A, R(A)=max(-q_A), W(A)=(P+R)/2, and
W_n=min_A W(A). The original objective is M_n=min_A max(P,R).

The target M_n-W_n=o(n^(3/2)) is NOT new. It is explicit in
`asymptotic_centered_width_recovery.md`; the stronger W_n>=M_n-1
is the distinguished-coordinate target in
`augmented_cut_code_midpoint_audit.md`. That artifact already proved
W_9=M_9=12 and W_10=M_10=13 by integer programming. The archive
`global_block_replacement_midpoint.md` already proves a global block
exchange inequality and explains why diffuse midpoint bias survives it.
The present contribution is a complete classification and an exact
single-row obstruction at an ACTUAL globally width-minimal matrix.
It proves neither convergence nor failure of asymptotic width recovery.

## Exhaustive results

The C++ source `computations/decisive_independent_exact_range_scan_2026_09_07.cpp`
enumerates all 2^((n-1)(n-2)/2) switching classes, represented by first
row all positive, and evaluates all 2^(n-1) projective spin states using
integer Gray-code updates. Counts do not quotient vertex permutations
or global negation. Every energy and objective comparison is exact.

| n | M_n | W_n | number of width minima | cap histogram of width minima |
|---|---:|---:|---:|---|
|3|3|2|2|3: 2|
|4|4|4|8|4: 6; 6: 2|
|5|4|4|12|4: 12|
|6|5|5|12|5: 12|
|7|9|8|3240|9: 3240|
|8|10|10|4200|10: 4200|
|9|12|12|1970640|12: 1607760; 14: 362880|

The prior JSON file of the basic scan is retained. The augmented source
also classifies one-vertex-deleted children. It fixes the known M_n
targets when accumulating extension flags; its independent full cap scan
recovers those exact targets, so the flags do not assume an unverified
cap value. Inputs are restricted to n<=9.

Among order-nine children obtained from width-minimal parents, 399840
switching classes occur. Of these, 75600 have NO extension attaining
both width 12 and cap 12. For each of the 362880 width-minimal cap-14
parents, count vertices whose row can be replaced to reach cap 12
while keeping width 12. The exact histogram is:

| repairable vertices | parents |
|---:|---:|
|0|120960|
|4|181440|
|9|60480|

Thus choosing the row optimally does not always repair an actual width
minimizer. This is stronger than failure for a prescribed row, but does
not exclude multi-row repair, a sequence of width-preserving row moves,
or existence of a different balanced width optimizer.

## Independent finite witness

The gauge code 134624811 uses lexicographic edges (i,j), 1<=i<j<9,
with bit one meaning negative; vertex zero has all positive neighbors.
Its integer spectrum has minimum -10 and maximum 14, so W(A)=12=W_9
and Q(A)=14>M_9.

The independent Python replay
`computations/decisive_independent_width_row_replay_2026_09_07.py`
deletes each of the nine vertices in turn and enumerates all 128 row
choices modulo negating the added row. For EVERY deleted vertex:

- the minimum extension cap is 14, even with no width constraint;
- the minimum extension width is 12;
- the minimum cap among minimum-width extensions is 14.

The row replay uses the exact identity for a fixed child B and row a:

    P(B,a)=max_x [q_B(x)+|a.x|],
    R(B,a)=max_x [-q_B(x)+|a.x|].

Fixing one coordinate of a loses no spectrum because a and -a differ
by switching the added vertex. The source prints the complete matrix
and endpoint histogram for every deleted vertex, making the witness
independently reproducible without the exhaustive C++ scan.

## What remains

The selectable optimized-value equality is still consistent with every
exact order through ten. The order-nine result specifically rules out
a proof that starts from an arbitrary width optimizer and repairs it
in one chosen row. It does not weaken the original required quantifier:
there need only exist an appropriate sequence of signings whose cap is
within o(n^(3/2)) of the optimized half-width.
