# Paired max-plus switching and ties: finite wind-tunnel report

**Reproducible script:**
`/home/math/quadra/tmp/paired_maxplus_windtunnel.py`

## Executive verdict

The conjecture

> uniform bounded error iff every reachable cycle in a finite
> active-face/product automaton has zero accumulated relative drift

is **false as a literal graph-cycle statement**, but has a correct and useful
version.

Three qualifications are necessary.

1. A cycle must be **pumpable in the underlying polyhedral dynamics**, not
   merely a cycle in an existential active-cell graph.
2. “Drift” is a **twisted/transported cocycle** determined by the selector
   product, not the untransported sum of kernel differences. A reset direction
   can absorb a nonzero raw offset.
3. At a tie, selectors cannot be branched independently. The true tangent
   update takes a max over the complete active set; a valid finite automaton
   must refine by the perturbation cone or work with exact paired values.

Under the additional hypothesis that the realizable paired selector language
is represented by a finite exact regular automaton, a corrected theorem does
hold: lift it by ordered coordinate pairs (four coordinates for a pair of
systems), attach the transported relative-offset cocycle, and test every
pumpable closed walk. Zero cycle class is then equivalent to a bounded vertex
potential by ordinary finite graph cohomology. There is no additional
nonperiodic obstruction **inside that exact finite lift**.

## 1. Exact paired propagation

For a max-plus matrix `S`, with rows indexing input coordinates and columns
indexing outputs,

```math
(F_Sx)_j=max_i(x_i+S_(ij)).
```

On a cell with selected maximizers `sigma(j)`, write

```math
F_Sx=P_sigma x+a_sigma,
\qquad (P_sigma x)_j=x_(sigma(j)),
\qquad (a_sigma)_j=S_(sigma(j),j).
```

For a paired generator `(S,T)`, if `sigma` is active at `x` and `tau` is
active at `y`, the **exact** update is

```math
(x',y')=(P_sigma x+a_sigma,
          P_tau y+b_tau).                           \tag{1.1}
```

If `sigma=tau`, the error closes:

```math
e'=P_sigma e+(b_tau-a_sigma).
```

If the selectors differ, `e=y-x` does not close; the base state is genuinely
needed. This is why a single-selector residual automaton is insufficient.

There is nevertheless a finite exact observable lift once a selector-pair
word is fixed. For ordered pairs define

```math
L_(u,v;p,q)(x,y)
=(y_u-y_v)-(x_p-x_q).                              \tag{1.2}
```

Then one paired affine step gives

```math
L_(u,v;p,q)(x',y')
=L_(tau(u),tau(v);sigma(p),sigma(q))(x,y)
 +(b_u-b_v)-(a_p-a_q).                             \tag{1.3}
```

The Hilbert error is half the maximum of (1.2) over initial states
`(u,v;p,q)=(j,k;j,k)`. Thus the exact relative cocycle lives on at most
`r^4` ordered-coordinate states, crossed with the active-face and word
automata. Formula (1.3) is the minimal finite bookkeeping missing from a
separately paid left/right analysis.

For a fixed paired selector cycle `C`, its affine composite has the form

```math
A_C(z)=P_Cz+c_C,
```

and

```math
A_C^mz=P_C^mz+sum_(s=0)^(m-1)P_C^s c_C.            \tag{1.4}
```

Because `P_C` is a coordinate selector, its powers are eventually periodic.
The repeated cycle is projectively bounded exactly when the offset has zero
class on every surviving functional-cycle direction--equivalently it is a
twisted coboundary there. A raw condition `c_C=0` is too strong; after a
rank-one selector, all old projective directions have disappeared.

## 2. Tiny exact examples

All examples below use rational arithmetic. For a two-state vector normalized
as `(0,z)`, Hilbert distance is `|z-z'|/2`.

### 2.1 Unique-selector drift

The sparse diagonal kernels

```math
S=((0,-infinity),(-infinity,0)),
T=((0,-infinity),(-infinity,1))
```

act as `z->z` and `z->z+1`. Starting together at zero, their errors at depths
one through eight are

```text
1/2  1  3/2  2  5/2  3  7/2  4
```

This is the basic nonzero pumpable cocycle.

### 2.2 Each letter is bounded, but switching drifts

Let the raw action of both letters be reflection `z->-z`. Let the perturbed
letters be

```math
a:z->1-z,
\qquad b:z->-z.
```

Each letter repeated alone has period two and paired error at most `1/2`.
But the two-letter product is a translation:

```math
b circ a:z->z-1.
```

The alternating word has errors

```text
1/2, 1/2, 1, 1, 3/2, 3/2, 2, 2, 5/2, 5/2, 3, 3.
```

Exhaustive enumeration over all binary words through depth ten found the
alternating word as a maximizer, with depth-ten error `5/2`.

This kills any test which checks generators separately or only sums raw
offsets on simple base-graph loops. The selector action transports the first
offset before the second is added. One must inspect closed walks in the
selector-semigroup skew product; the word `ab` is the first obstruction.

### 2.3 A reachable symbolic cycle need not be pumpable

The all-finite matrix

```math
S_c=((0,0),(c-1,c))
```

acts as

```math
z->clip(z+c,0,1).
```

Compare `S_0` and `S_(1/4)`. Both images lie in `[0,1]`, so their paired
Hilbert error is uniformly at most `1/2` for every input and every depth.
On the middle active cell, however, the second map translates by `1/4`, and
the active-cell transition graph contains a self-loop with nonzero apparent
drift. No point can traverse that middle self-loop more than four times: it
hits the saturation face and stays there. For seed zero the exact errors are

```text
1/8, 1/4, 3/8, 1/2, 1/2, 1/2, ...
```

Thus “reachable graph cycle” is not enough for necessity. The cycle must
admit arbitrarily long repetitions satisfying all cell inequalities. If the
intended bound is `O(c)` uniformly as `c->0`, this family still correctly
acts as an obstruction--its saturation error is order one--but it is not an
obstruction to boundedness for each fixed pair. The scale of the conjecture
must therefore be stated.

### 2.4 Tied selectors can create spurious surviving directions

The matrix

```math
S_tie=((0,1),(-1,0))
```

maps **every** projective input to `z=1`. At `z=1`, both inputs tie in both
output columns. A naive selector automaton permits, among others, the
identity selector and the swap selector, neither of which is a projective
reset as a linear map.

But the exact directional derivative on a tie face is

```math
(DF_S(x)e)_j=max_(i in I_j(x))e_i.                 \tag{2.1}
```

Here `I_1=I_2={1,2}`, so both outputs equal `max(e_1,e_2)` and the true
projective derivative is identically zero. The script checks this for all
integer directions `e_2-e_1` from `-5` through `5`.

Hence independently choosing one tied maximizer per output produces false
nonreset cycles. A tie-aware automaton must include the common inequalities
which make each selector optimal for the **same** perturbation direction, or
retain (2.1) as a max-plus tangent map.

## 3. Corrected finite certificate

The following is the strongest statement supported by the wind tunnel.

### Proposed finite paired-cycle theorem

Assume:

1. the permitted generator words form a regular language;
2. there is a finite polyhedral automaton whose paths are exactly the
   realizable paired active-face itineraries (including tie-direction cones);
3. every directed graph cycle is marked pumpable only when all its repeated
   affine cell constraints are feasible; and
4. initial endpoint functionals in (1.3) are uniformly bounded, or a prior
   selector reset removes them.

Cross this automaton with the ordered-coordinate lift (1.3). Delete states
whose relevant ordered pairs have coalesced, since they carry no projective
memory. On every remaining pumpable strongly connected component, regard the
last term of (1.3) as an additive real edge label.

Then the paired projective error is uniformly bounded over all legal words
if and only if every pumpable closed walk has zero label (after endpoint
gauge/coboundary normalization). Equivalently, on each such component there
is a vertex potential whose difference equals every edge label. A nonzero
closed-walk label can be pumped to give linear drift; zero labels telescope,
and the acyclic condensation contributes only a finite transient.

This is ordinary finite graph cohomology **after** the noncommutative selector
transport has been put into the finite coordinate lift. It does not follow by
placing the raw kernel differences on the base active-cell graph.

### Proof outline

Iterate (1.3) backwards through a word. The generalized relative coordinate
equals its bounded initial endpoint term plus the sum of the edge labels on
the lifted path. There are finitely many such coordinates. If a pumpable
closed walk has nonzero sum, repeating it makes one generalized relative
coordinate grow linearly. Conversely, if all closed-walk sums vanish, choose
a root in each recurrent component and define a potential by any root-to-node
path sum; zero cycles make it path-independent. Labels telescope to endpoint
potential differences. A finite acyclic prefix/suffix adds a bounded amount.

No genuinely nonperiodic product can evade this argument once hypotheses
1--3 provide an **exact finite** path language: every long lifted path
decomposes into a bounded simple part plus closed walks. Without exact
regularity, however, a finite existential cell graph is only an overapproximation;
its nonzero cycles can be nonpumpable, as Example 2.3 shows, and an aperiodic
feasible language need not be decided by raw graph cycles.

## 4. Computational protocol

For rational kernels and small `r`, an exact wind tunnel should implement:

```text
frontier = {(normalized x_0, normalized y_0)}
for depth = 1,...,T:
    next = empty map
    for paired state and every legal generator ell:
        x' = normalize(maxplus(S_ell,x))
        y' = normalize(maxplus(T_ell,y))
        record exact active sets ArgMax(S_ell,x), ArgMax(T_ell,y)
        record d_H(x',y') and predecessor word
        insert (x',y') into next
```

For symbolic certification:

```text
1. partition paired projective space by all active comparisons;
2. intersect each proposed transition with its exact preimage inequalities;
3. refine tie faces by tangent max regions, not arbitrary selectors;
4. for a candidate cycle, solve feasibility of all iterated cell constraints;
5. attach the r^4 ordered-coordinate transition and relative offset (1.3);
6. search pumpable SCCs for a nonzero cocycle cycle;
7. return either a repeated word with linear drift or vertex potentials.
```

For a cycle composite with coordinate-selector linear part, iterated
feasibility is decidable by decomposing the finite functional graph of the
selector: its powers are eventually periodic and its offsets are eventually
linear. This reduces the invariant-cell test to finitely many rational linear
inequalities plus recession-direction checks.

## 5. Research judgment

The paired-cycle idea survives, but not in its naive form. The right finite
object is a **tie-consistent, pumpability-certified, ordered-coordinate
skew product**. Its cocycle criterion is a real extension of the current
fixed-itinerary theorem because it transports cancellation jointly through
different left/right selectors. The two-state reflection example proves
that this transport is indispensable; the tie-reset and clip examples show
exactly where an unrefined active-face graph lies.

The most useful next theorem is the proposed finite paired-cycle theorem for
tie-free rational polyhedral cells, with an explicit LP criterion for cycle
pumpability. Ties should be added only afterward via the max-plus tangent
map (2.1).
