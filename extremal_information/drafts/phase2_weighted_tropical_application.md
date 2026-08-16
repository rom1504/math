# Weighted tropical exposure on code transversals

**Status:** focused application/no-go for Theorem WE.1.  The argument below
is self-contained once the transversal block in Sheshadri's exact trellis-rank
theorem is granted.  It is not a surface-file promotion.

## 1. Verdict

The weighted tropical-exposure theorem does not turn the canonical
Sheshadri transversal into a macroscopic mean-square lower bound under either
of the two natural diffuse query laws:

1. uniform roots in the whole conditional distance table; or
2. independent uniform left and right states in the transversal block.

More generally, let the transversal have `q=2^s` states and normalize the
conditional distance table to take values in `[0,1]`.  If `p_*` is the
largest probability of a cell in the `q by q` transversal block, then every
crossing-witness graph on those anchors satisfies

```math
\boxed{
\mathfrak m_k\le \left\lfloor{q\over2}\right\rfloor p_* .}
                                                               \tag{WT.1}
```

Thus a fixed positive exposure lower bound requires an atom of order
`1/q`, whereas an independent uniform state-pair query has atoms `1/q^2`.
Equivalently, its density must be tilted by a factor `Omega(q)=Omega(2^s)`.
This is a scalable obstruction to this witness architecture, not a
counterexample to average-error tropical-rank lower bounds by other methods.

An explicit graph-code family shows the decay is real rather than an artifact
of a loose estimate: for one channel and the complete transversal witness
graph, the parameter is exactly `1/(8q)`.

## 2. Conditional code table and transversal

Let `C<=F_2^m`, split the coordinates into nonempty sets `L` and `R`, and
write

```math
W(x_L,x_R)=d((x_L,x_R),C).                          \tag{WT.2}
```

Let `C_L` and `C_R` denote the subcodes supported on the corresponding
coordinate sets and put

```math
s=\dim C-\dim C_L-\dim C_R,
\qquad q=2^s.                                      \tag{WT.3}
```

The transversal in Sheshadri's proof gives distinct row labels
`u_1,...,u_q` and distinct column labels `v_1,...,v_q` such that

```math
W(u_i,v_i)=0,
\qquad W(u_i,v_j)\ge1\quad(i\ne j).                \tag{WT.4}
```

The row and column labels really are separately injective.  If two chosen
codewords had the same left part, their difference would lie in `C_R`; if
they had the same right part, their right projections would represent the
same quotient class.  Either event contradicts the choice of distinct
trellis-state classes.

Normalize

```math
M={W\over m}.                                      \tag{WT.5}
```

Then `0<=M<=1`.  Use the `q` diagonal cells `(u_i,v_i)` as anchors and let
`E` be any graph on their indices for which the positive crossing gaps are
declared as witnesses.  Let `mu` be any probability law on the full
conditional table and define

```math
p_*=\max_{1\le i,j\le q}\mu(u_i,v_j).              \tag{WT.6}
```

Zero-mass cells are allowed, with the convention in (WE.2).

## 3. Exposure-concentration theorem

### Theorem WT.1 (a transversal needs one-dimensional query concentration)

For every channel budget `k>=1`, every witness graph `E` on the Sheshadri
anchors, and every query law `mu`,

```math
\mathfrak m_k(M,E,\mu)
\le \left\lfloor{q\over2}\right\rfloor p_*.
                                                               \tag{WT.7}
```

Consequently:

1. Under the uniform law on the whole `2^|L| by 2^|R|` table,

   ```math
   \mathfrak m_k\le 2^{s-m-1}
   \le 2^{-m/2-1}.                                 \tag{WT.8}
   ```

2. Under the uniform law on the `q by q` transversal block,

   ```math
   \mathfrak m_k\le {1\over2q}=2^{-s-1}.           \tag{WT.9}
   ```

3. If `mu` is supported on that block and has density at most `D` relative
   to its uniform law, then

   ```math
   \mathfrak m_k\le {D\over2q}.                    \tag{WT.10}
   ```

   Hence `mathfrak m_k>=delta>0` forces

   ```math
   D\ge2\delta q,
   \qquad
   H_\infty(\mu)\le s-\log_2(2\delta).             \tag{WT.11}
   ```

   Independent uniform state-pair sampling has min-entropy `2s`; a
   nonvanishing certificate of this form requires min-entropy at most
   `s+O_delta(1)`.

#### Proof

For an edge `ij`, the normalized crossing gap is

```math
G_{ij}=M(u_i,v_j)+M(u_j,v_i)\le 2,                \tag{WT.12}
```

because the two diagonal anchor values vanish and each crossed value lies in
`[0,1]`.  If one of its four witness cells has zero mass, its exposure weight
is zero.  Otherwise every one of those masses is at most `p_*`, so

```math
\sum_{c\in C_{ij}}{1\over\mu(c)}\ge {4\over p_*}.
```

Therefore

```math
w_{ij}
={G_{ij}^2\over\sum_{c\in C_{ij}}1/\mu(c)}
\le p_*.                                          \tag{WT.13}
```

Every matching on `q` indices has at most `floor(q/2)` edges.  This bounds
the weight of every monochromatic matching for every coloring.  Taking first
the maximum over matchings and then the minimum over colorings proves
(WT.7).

For the full uniform law, `p_*=2^{-m}`.  Moreover,

```math
s\le\min\{|L|,|R|\}\le m/2.                       \tag{WT.14}
```

Indeed, `C/(C_L+C_R)` injects into each appropriate quotient of a coordinate
projection, whose dimension is at most the corresponding block size.  This
proves (WT.8).  For uniform sampling in the transversal block,
`p_*=q^{-2}`, giving (WT.9).  Under a density cap `D`,
`p_*<=D/q^2`, giving (WT.10).  If the left side is at least `delta`, (WT.7)
implies `p_*>=2delta/q`; the density and min-entropy conclusions follow.
`square`

### Scope of the entropy conclusion

Equation (WT.11) is a min-entropy/maximum-atom necessity.  It does **not**
by itself imply a comparable Shannon relative-entropy cost: one unusually
heavy atom can raise the maximum density without carrying linear KL
divergence.  A KL or mutual-information obstruction would need a massive
family of simultaneously necessary exposed cells, not merely the matching
argument in WE.1.

## 4. Exact structured example

For `t>=1`, let

```math
C_t=\{(z,z):z\in\mathbb F_2^t\}
\subseteq\mathbb F_2^t\times\mathbb F_2^t.         \tag{WT.15}
```

Here `m=2t`, `C_L=C_R={0}`, `s=t`, and `q=2^t`.  Coordinatewise
minimization gives the complete conditional table

```math
W_t(x,y)=d_H(x,y).                                 \tag{WT.16}
```

Take the anchors `(z,z)`, the complete witness graph, the normalized table
`M_t=W_t/(2t)`, and the uniform law on all `q^2` table cells.

### Proposition WT.2 (exact one-channel exposure)

For this family,

```math
\boxed{
\mathfrak m_1(M_t,K_q,\mathrm{Unif})={1\over8q}
={2^{-t}\over8}.}                                  \tag{WT.17}
```

#### Proof

For anchors `x,y`,

```math
G_{xy}
={W_t(x,y)+W_t(y,x)\over2t}
={d_H(x,y)\over t}.                               \tag{WT.18}
```

Every witness cell has mass `q^{-2}`, hence

```math
w_{xy}={d_H(x,y)^2\over4t^2q^2}
\le {1\over4q^2}.                                 \tag{WT.19}
```

For `k=1` all edges are monochromatic.  A matching has at most `q/2`
edges, so its weight is at most `1/(8q)`.  Pair every `x` with its bitwise
complement.  This is a perfect matching and every edge has Hamming distance
`t`, attaining the bound. `square`

The code has exact min-plus factor rank `q=2^t`, and the complete canonical
witness graph has the largest possible set of positive transversal
crossings.  Nevertheless its WE.1 mean-square certificate for a rank-one
approximation decays as `2^{-t}` at fixed normalized distortion.  Replacing
the complete witness graph by a subgraph cannot improve `mathfrak m_1`, since
it only removes candidate matchings.

For this explicit family the vanishing certificate reflects an actual
low-rank average approximation, not merely a weakness of WE.1.  Under the
uniform cell law, the rank-one constant table `\widetilde M_t=1/4` has

```math
\mathbb E(\widetilde M_t-M_t)^2
=\operatorname{Var}\left({d_H(X,Y)\over2t}\right)
={1\over16t},                                      \tag{WT.20}
```

because `d_H(X,Y)` is `Binomial(t,1/2)`.  Thus exact exponential tropical
rank and macroscopic uniform approximation hardness can coexist with
vanishing normalized mean-square rank-one error.

## 5. What this rules out, and what remains open

The result rules out the following proposed inference:

```text
exponential exact trellis rank
  + the canonical transversal crossing gaps
  + diffuse independent root/state queries
  => a nonvanishing average-response rank obstruction.
```

The failure is geometric and scalable.  WE.1 forbids double counting by
passing to a matching, so at most `q/2` witnesses can be charged.  Diffuse
state-pair sampling assigns only order `q^{-2}` mass per four-cell witness.
The product is order `q^{-1}`.

For a general code table, the theorem does not by itself construct a good
low-rank mean-square approximation; it only says that the canonical
transversal plus pairwise crossing/matching parameter cannot prove the
opposite.  For the explicit graph-code family, (WT.20) supplies such an
approximation directly.  A successful average-error theorem for a wider
class would need at least one of:

1. a higher-dimensional witness that charges `Theta(q^2)` diffuse cells
   without invalid double counting;
2. a declared correlated query law with genuine mass on a sparse exposed
   skeleton; or
3. a global spectral, communication, or information argument not reducible
   to disjoint four-cell crossings.

The first option is the mathematically substantive open door.  Merely
choosing another graph on the same anchors cannot help, by Theorem WT.1.

## 6. Hypothesis audit

The live statement of WE.1 assumes that the anchor rows `i -> x_i` and
anchor columns `i -> y_i` are separately injective.  This is exactly what
makes the four-cell witness sets of a matching disjoint.  Sheshadri's
transversal satisfies both injectivity conditions, as checked in Section 2.
Thus Theorem WT.1 and Proposition WT.2 apply directly, without any extension
of WE.1.

## 7. Director assessment

This application is a genuine no-go, not a new positive complexity law.
It quantitatively explains why the exact `2^s` trellis-state lower bound and
its lattice-scale uniform robustness do not survive the most natural
average-query relaxation.  The surviving parameter is useful as a detector
of query concentration, but in this model it does not bridge exact tropical
rank to macroscopic mean-square extremal information.

The minimal next question is therefore not another pairwise witness graph.
It is whether a code conditional table admits a **massive joint witness**
whose loss inequality sums over a positive fraction of the state-pair table
before absolute values or squares are paid.  Without such an object, the
weighted tropical route stops at the concentration diagnosis above.
