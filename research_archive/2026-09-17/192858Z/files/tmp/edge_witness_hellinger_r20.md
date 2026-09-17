# Wave 20: edge-flip witnesses do not control coordinate Hellinger overlap

This note tests only the information supplied by exact one-edge flip
stability.  It derives the exact fractional-cover and Gibbs-overlap formulas,
then gives two complementary scoped models showing where that information
stops.  Neither model is a counterexample to (10.617).

## 1. Exact witness incidence and its fractional cover

Let `B` be an order-`r` complete signing with `q=Q(B)`.  On the lifted
oriented state space write

```math
s_e(\omega)=\sigma b_{ij}x_ix_j\in\{\pm1\},
\qquad
e_B(\omega)=2\sum_es_e(\omega),
\qquad
\Delta_\omega=q-e_B(\omega).
```

Put

```math
L_4=\{\omega:\Delta_\omega\le4\},
\qquad
W_e=\{\omega\in L_4:s_e(\omega)=-1\}.
```

The one-edge case of (3.21) says exactly that an exact `Q`-minimizer obeys

```math
W_e\ne\varnothing\qquad\text{for every edge }e.
```

The natural fractional state cover is

```math
\tau_4(B)=
\min\left\{
\sum_{\omega\in L_4}\lambda_\omega:
\lambda_\omega\ge0,
\ \sum_{\omega\in W_e}\lambda_\omega\ge1\ (e\in E)
\right\}.
```

Its exact dual assigns edge weights `y_e>=0` subject to

```math
\sum_{e:s_e(\omega)=-1}y_e\le1
\qquad(\omega\in L_4).
```

This cover is intrinsically constant-scale.  If `N=binom(r,2)` and
`N_omega={e:s_e(omega)=-1}`, the energy identity gives

```math
\boxed{
|N_\omega|
=\frac N2-\frac q4+\frac{\Delta_\omega}{4}.
}
```

Thus every low witness covers at most `N/2-q/4+1` edges and

```math
\boxed{
\tau_4(B)
\ge
\frac{N}{N/2-q/4+1}
=2+O(r^{-1/2})
}
```

for a competitive signing.  The certificate itself therefore forces only
a constant fractional multiplicity.  Cut consistency may raise the cover
number in a particular matrix, but no growing lower bound follows from the
capacity identity.

## 2. Exact link to coordinate-flip pairs

For a state `omega`, define its negative degree and signed row field at
coordinate `i` by

```math
d_i^-(\omega)=|\{j:s_{ij}(\omega)=-1\}|,
\qquad
h_i(\omega)=(r-1)-2d_i^-(\omega).
```

Flipping coordinate `i` changes the doubled energy by `-4h_i`, hence

```math
\boxed{
\Delta_{\tau_i\omega}=\Delta_\omega+4h_i(\omega).
}
```

For Gibbs weights `w_omega=e^{-beta Delta_omega}` and
`D_beta(B)=sum_omega w_omega`, the directed geometric mass and affinity are
therefore exactly

```math
\boxed{
g_i(\omega)
=\sqrt{w_\omega w_{\tau_i\omega}}
=e^{-\beta[\Delta_\omega+2h_i(\omega)]},
\qquad
\operatorname{BC}_i
=\frac1{D_\beta(B)}\sum_\omega g_i(\omega).
}
```

This is (10.623) written state by state.

To expose the direction of a fractional incidence argument, choose for every
edge a probability distribution `eta_{e,omega}` supported on `W_e`, and put

```math
ell_i(\omega)=\sum_{j\ne i}\eta_{\{i,j\},\omega}.
```

Then

```math
\sum_{i,\omega}\ell_i(\omega)=2N,
\qquad
0\le\ell_i(\omega)\le d_i^-(\omega).
```

But `g_i(omega)` is increasing in `d_i^-(omega)`:

```math
g_i(\omega)
=e^{-\beta[\Delta_\omega+2(r-1)-4d_i^-(\omega)]}.
```

Consequently witness load can only provide **lower** bounds on geometric
mass and affinity.  It cannot give the upper bound
`BC_i<=exp(-Theta(sqrt(r)))` required by (10.624).  More basically, the
incidence cover contains no information about the denominator `D_beta(B)`
or about the coordinate partners of nonwitness states.

For each individual low witness,

```math
\sum_i h_i(\omega)=q-\Delta_\omega=\Theta(r^{3/2}),
```

so some coordinate contribution of that witness is exponentially small.
Different witnesses can place this coordinate differently, and one small
summand does not upper-bound the full normalized affinity.  This is the
precise overlap gap.

## 3. Incidence-only hard-ball model

The following model is deliberately **not** claimed to be a quadratic Gibbs
law.  It identifies exactly what raw witness incidence plus cube geometry
cannot prove.

Take lifted states `(sigma,x)`, genuine augmented-cut signs

```math
s_{ij}(\sigma,x)=\sigma x_ix_j,
```

and the true coordinate-flip involutions.  Let `B_{r,k}` be the Hamming ball
of radius `k` around the all-positive spin, take both orientation copies,
assign deficit zero on this support and infinite deficit outside, and use the
resulting normalized hard Gibbs measure.  The state `(sigma=-1,x=1)` has
`s_e=-1` for every complete edge, so it is a deficit-zero witness for every
one-edge certificate.  All sign patterns obey the augmented-cut triangle
relations.

The coordinate affinities are nevertheless

```math
\boxed{
\operatorname{BC}_i
=
\frac{2\sum_{j=0}^{k-1}\binom{r-1}{j}}
     {\sum_{j=0}^{k}\binom rj}.
}
```

For `k=floor(sqrt(r))`, this is asymptotic to `2/sqrt(r)`, so the reward is
only `O_beta(log r)=o(sqrt(r))`.

Scope is essential: these artificial deficits violate

```math
\Delta_\omega=q-2\sum_es_e(\omega)
```

for any single complete signing.  In particular the two orientation copies
of the same spin cannot both be top states when `q>0`.  The model satisfies
actual cut signs, coordinate adjacency, normalized weights, and every raw
edge-witness incidence, but not single-signing energy realizability.  It
disproves only an incidence/support-count proof.

## 4. Energy-realizable sparse quadratic model

There is a complementary asymptotic model which retains exact quadratic
energies and Gibbs normalization but loses completeness.  Let

```math
B_t=\bigoplus_{a=1}^t A_9
```

with zero interblock edges.  Since `A_9` is balanced and has norm 24,

```math
Q(B_t)=24t.
```

Flipping any present edge changes one `A_9` block.  That block still has
norm at least 24, and the remaining balanced blocks can be aligned in its
attaining orientation, so every present-edge flip has norm at least `24t`.
Equivalently, the `A_9` low witness for that edge extends by exact grounds in
the other blocks and has total deficit at most four.

The one-orientation partition polynomial of `A_9` is symmetric under energy
negation.  Hence the oriented pressure of `B_t` factorizes.  Deleting a
coordinate in one block cancels all other block factors, giving exactly

```math
\boxed{
\kappa_{\beta,i}(B_t)
=\kappa_{\beta,i}(A_9)<3
\qquad\text{for every }i,t,\beta>0.
}
```

Thus edge witnesses, the exact energy-deficit relation, and genuine
quadratic Gibbs weights still do not force square-root reward on a sparse
interaction graph.  This is not a signing counterexample: the interblock
entries are zero, and there are no certificates for the missing complete
edges.  A positive theorem must use complete-graph coupling essentially.

## 5. Actual finite complete-signing audit

The order-nine minimizer supplies the finite wall with every required axiom.
Its deficit-at-most-four layer is exactly its 25 oriented-projective grounds.
The edge-witness fractional cover has exact value

```math
\boxed{\tau_4(A_9)=4.}
```

An exact primal of total four, with weights in `{1/3,2/3}`, is recorded and
integer-checked in the verifier.  For the dual, give weight one to the four
edges

```text
(0,2), (0,4), (2,5), (4,5).
```

Every one of the 25 grounds is unfavorable on exactly one of these four
edges, so the dual value is four.  This proves equality without relying on
floating-point LP output.

The same verifier evaluates (10.623) on the full lifted state space; the
global-spin duplication cancels from numerator and denominator.  At
`beta=0.1,0.5,1` every coordinate has `kappa<3`, consistent with the
all-temperature proof (10.616).  Thus even nontrivial exact witness overlap
does not force a large reward at finite order.

## 6. Scoped conclusion and quantifiers

- The edge-flip certificates are existential covers of complete edges by a
  constant-capacity low-deficit layer.  Their natural fractional loads point
  toward **large**, not small, Hellinger affinity.
- The hard-ball model closes arguments using only witness incidence, cut-word
  support, and coordinate adjacency.  It is not energy-realizable.
- The block-`A_9` model closes the same argument for genuine quadratic Gibbs
  weights on sparse graphs.  It is not a complete signing.
- Actual `A_9` is a finite complete-signing wall, not an asymptotic
  counterexample.

Therefore none of these results falsifies literal (10.617), nor its weaker
target-specific convergence quantifier.  A surviving positive proof must use
the joint fact that **all** complete-graph edges arise from one signing,
relate the whole near-ground Gibbs denominator to coordinate partners, and
produce one common coordinate with affinity `exp(-Theta(sqrt(r)))`.  Raw
multiplicity, overlap, or fractional edge coverage does not supply that fact.

All finite and combinatorial claims are checked in
`edge_witness_hellinger_r20.py`.
