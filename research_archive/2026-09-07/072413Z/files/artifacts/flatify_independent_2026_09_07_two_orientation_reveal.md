# A two-orientation reveal with an exact original-cap child endpoint

Date: 2026-09-07. Status: exact alternative interpolation; no endpoint
monotonicity or sublinear cumulative drift bound is proved.

The optimized shared-orientation absolute pressure at a two-child endpoint
can cancel the children's energy asymmetries. In particular equal-child
endpoints minimize to twice optimized half-width pressure. The following
alternative removes that stronger-than-needed endpoint requirement.

## Definition and endpoint identities

Use the posterior variance profiles v(k,r) of the companion reveal note.
Let P be the revealed-positive vertices, and write H_P for the Hamiltonian
of edges internal to P, and H_R=H-H_P for all remaining edges. Put

```math
G_{k,r}(beta)=min_A log sum_x
 cosh[beta H_P(x)/sqrt(N-1)] cosh[beta H_R(x)/sqrt(N-1)].
```

Equivalently average exp[beta(sigma_1 H_P+sigma_2 H_R)/sqrt(N-1)] over
two independent orientation signs sigma_1,sigma_2. Initially H_P=0, so

```math
G_{0,0}=Psi_N^abs(beta).
```

At the terminal profile P and its complement are disconnected, the two
orientation signs are independent, and minimization separates EXACTLY:

```math
G_{N,m}=Psi_m^abs(beta)+Psi_n^abs(beta).
```

This does not rely on global child polarity selection, and is not a
half-width endpoint. For an intermediate fixed signing its associated
zero-temperature quantity is

```math
max_{x,sigma_1,sigma_2}(sigma_1 H_P+sigma_2 H_R)
 = max_x (|H_P|+|H_R|)
 = max_{x,tau=+-1}|tau H_P+H_R|.
```

Thus it is the robust absolute cap against one allowed internal-block
polarity reversal. Its optimized value C_{k,r} obeys

```math
beta C_{k,r}/sqrt(N-1)-log4 <= G_{k,r}
 <= beta C_{k,r}/sqrt(N-1)+Nlog2.
```

At the start C=M_N, while at the end
C/sqrt(N-1)=M_m/sqrt(m-1)+M_n/sqrt(n-1).

For every fixed profile G_*>=F_*^abs: the two-orientation partition
function is the average of the ordinary absolute partition functions
before and after reversing all internal-P edge signs; both are at least
the optimized ordinary partition function.

## Count cleanup survives without a new hypothesis

At an edge optimum the flip observable is now sigma_c A_e x_i x_j, where
c is the edge's orientation group. The flip-ratio proof still gives
d_e<=tanh(lambda_e) and the same one-sided variance derivative bound.
Scaling all edges incident to a vertex is still monotone: condition on
the two orientations and all other spins, sum the selected spin, and
obtain a nonnegative weight times cosh(t rowfield). Vertex insertion by
independent row signs gives the same [log2,log2+beta^2/2] increment.

For consecutive revealed counts, delete the single known vertex whose
label changes. On the remainder BOTH the orientation-group assignments
and the known labels agree. The variance comparisons and unknown-only
diagonal contraction from Section 7 of the count-concentration note
therefore apply verbatim. Terminal G is independent of the particular
partition, so the common-terminal-value tail comparison also survives.
Consequently

```math
sup_k E|G_{k,K_k}-G_{k,r_0(k)}|<=C_delta beta^2 sqrt(N).
```

The corresponding robust caps concentrate with normalized expected error
O_delta(N^(-1/4)) at beta=N^(1/4).

## What is still missing

An actual useful one-sided endpoint estimate would be

```math
G_{0,0}(beta)<=G_{N,m}(beta)+O_delta(beta^2 sqrt(N))
```

at beta=N^(1/4), or another error permitting a sublinear normalized
energy recurrence. It would imply directly

```math
M_N/sqrt(N-1) <= M_m/sqrt(m-1)+M_n/sqrt(n-1)
                 +O_delta(N^(3/4)).
```

No such estimate is proved here. Changing orientation assignment as the
revealed-positive set grows is a real part of the new operation and must
be included in its discrete generator. One cannot simply copy a smooth
fixed-orientation variance derivative and omit that assignment change.

## A distribution-free cavity inequality that does NOT suffice

Summing an individual bath spin produces a row functional involving
E_mu cosh(U+V+W), where U and V are its two revealed-block row fields.
Even allowing a relative global sign, comparison to the split fields
sqrt(2)U or sqrt(2)V has no universal favorable sign. For independent
U,V uniform on {+-h}, and W=0, both relative-sign choices give

```math
log E cosh(U+-V)=2logcosh(h),
log E cosh(sqrt(2)U)=logcosh(sqrt(2)h).
```

The first exceeds the second for h>0 (its small-h difference is h^4/6
+O(h^6), and its large-h difference grows as (2-sqrt(2))h).
Thus ordinary cosh convexity or an arbitrary cavity-source PSD condition
does not prove the desired drift. This is NOT a counterexample with an
actual globally optimizing Gibbs source or with unrestricted optimized
row signs; those restrictions are precisely where any surviving theorem
must obtain additional information.
