# Independent audit: matroid contraction and projective congruences

**Scope.** This audit reconstructs MQ.1, TC.1, PMQ.1, PMQ.2, the
state-count asymptotics, and the response-packing comparison in
`phase3_projective_matroid_projection_quotient.md`.  It does not assume the
theorem-builder's narrative.  I also exhaustively enumerated every join
congruence of `L(F_2^w)` through `w=3` and tested the proposed decomposition.

## Verdict

All substantive claims are correct under the declared deterministic scalar
closed-summary model.  No hidden use of the full energy landscape or target
order occurs.  The fixed-subspace quotient is a strict, exactly composable
quotient, and the zero-separating remainder in PMQ.2 is real rather than a
vacuous placeholder: TC.1 supplies nontrivial examples.

PMQ.2 is a classification of the **form** of projective join congruences,
not a solution of their optimal state-count problem.  It leaves precisely
the advertised question of how much a zero-separating remainder can merge.

## 1. Reconstruction of MQ.1

For flats `X,Z,W` of a finite matroid,

```math
h_W(X\vee Z)=(X\vee Z)\vee W
             =(X\vee W)\vee(Z\vee W),
```

so `h_W` is a join homomorphism, and each `Y>=W` is hit by `Y`.  If
`h_W(X)=Y`, then `X<=Y` and

```math
r(Y)=r(X\vee W)\le r(X)+r(W)=r(X)+d.
```

Consequently

```math
R-r(Y)\le F_M(X)\le R-r(Y)+d.
```

The printed decoder is the midpoint of this interval.  The bottom flat and
`W` lie in the fibre over `W` and have ranks `0,d`, so the quotient's worst
fibre width really is `d`; `d/2` is not merely a loose upper bound.  Since
all future joins take place in the exact quotient, this error is decoded
once and cannot accumulate.

This reasoning remains valid with loops: the bottom flat may be nonempty as
a set, but it still has rank zero and lies below every flat.

## 2. Reconstruction of TC.1

The only nontrivial congruence check begins with distinct `U,V>=X` satisfying
`U\vee Y=V\vee Y`.  For every `Z`, both joins remain above `X` and

```math
(U\vee Z)\vee Y=(U\vee Y)\vee Z
                 =(V\vee Y)\vee Z=(V\vee Z)\vee Y.
```

Thus the conditional relation is a genuine join congruence, including when
joining crosses the trigger boundary.  Its classes above `X` are exactly the
nonempty fibres indexed by `P>=Y`; all other flats are singletons.  Hence the
state count is

```math
|L(M)|-|[X,E]|+|[Y,E]|.
```

For `U\vee Y=P`, submodularity and `U\wedge Y>=X` give

```math
r(P)-d\le r(U)\le r(P),
\qquad d=r(Y)-r(X).
```

The class indexed by `Y` contains both `X` and `Y`, so the global width is
exactly `d` even for a general matroid.  When `r(X)>0`, zero is outside the
trigger interval and remains a singleton.  When `Y<E`, the class containing
`X,Y` moves after joining a flat not contained in `Y`; it is therefore not
the absorbing class of a Rees quotient.

In projective geometry every fibre, not only the worst fibre, has width
exactly `d`: in `P/X`, choose a complement to `Y/X` and lift it to a subspace
`U` with `U\cap Y=X` and `U+Y=P`.

## 3. Reconstruction of PMQ.1 and the state count

If `Y` is a subspace of `V/W` and `h_W(X)=Y`, rank--nullity gives

```math
\dim X=\dim Y+\dim(X\cap W).
```

The intersection dimension ranges over both endpoints `0,d`: take a linear
section of `Y`, then add `W`.  The residual-rank fibre is therefore exactly

```math
[w-\dim Y-d,\;w-\dim Y],
```

which verifies the decoder `w-dim(Y)-d/2` and exact error `d/2`.

For `n=w-d`, the standard bounds

```math
q^{j(n-j)}\le {n\brack j}_q\le C_q q^{j(n-j)}
```

show that the largest Gaussian coefficient has exponent
`floor(n^2/4)`, while summing over `n+1` dimensions costs `O(log n)` bits.
Thus PMQ.6--PMQ.8 are correct for fixed `q`.  For fixed trigger rank
`k>=1`,

```math
N_q(w-k)/N_q(w)=q^{-\Theta(w)},
```

so `N_q(w)-N_q(w-k)+N_q(w-k-d)` retains the full
`(w^2/4)log_2(q)+o(w^2)` logarithmic state count, as claimed.

## 4. Reconstruction of PMQ.2

Let `theta` be a join congruence.  Its zero class is join-closed.  If `W` is
the join of all its members, finiteness gives `W theta 0`.  For every
`U<=W`, congruence compatibility applied to `W theta 0` gives

```math
W=W+U\mathrel\theta U=0+U.
```

Hence the zero class is exactly `L(W)`.  The same compatibility gives
`X+W theta X`, so equality after projection modulo `W` forces equivalence.
The projection-fibre congruence is therefore contained in `theta`, and the
standard quotient relation on `L(V/W)` is well-defined.  Its zero class is
trivial: a quotient subspace equivalent to zero would lift to a subspace in
the original zero class, hence one contained in `W`.

Conversely, pulling any zero-separating quotient congruence back along
`X -> (X+W)/W` has zero class exactly `L(W)`.  These constructions are
inverse, which proves uniqueness and completeness of the decomposition.

For a quotient subspace `Y`, its full projection fibre has original
residual ranks

```math
F_{V/W}(Y),F_{V/W}(Y)+1,\ldots,F_{V/W}(Y)+d.
```

Both endpoints occur.  Taking the union over one induced congruence class
therefore gives the exact identity

```math
\operatorname{osc}_{F_V}=d+
\operatorname{osc}_{F_{V/W}},
```

not just an inequality.  CSC.1 then yields the final `2 eta` fibre-width
criterion.

## 5. Exhaustive adversarial check

The reproducible script
`experiments/verify_phase3_matroid_quotients.py` performs the following.

- It verifies every fixed-subspace homomorphism and every exact fibre width
  through `w=4` (67 flats at `w=4`).
- It verifies every triggered contraction through `w=3`, including all
  congruence compatibility triples, state counts, exact nontrivial fibre
  widths, zero separation, and nonabsorption.
- It enumerates **all** join congruences, rather than a sample, through
  `w=3`.  The counts are `1,2,12,3616` for widths `0,1,2,3`.  Enumeration is
  complete because the search begins at equality and closes every possible
  additional pair under all joins; every finite congruence is generated by
  some sequence of its pairs.
- For all 3,616 width-three congruences it checks the zero-class formula,
  containment of every projection fibre, recovery from the induced
  zero-separating congruence, and the exact oscillation identity class by
  class.

The saved result is
`experiments/phase3_matroid_quotients_results.json`; all checks pass.

## 6. Packing comparison and scope

For middle-dimensional binary subspaces the complete response distance is
the injection distance.  Taking separation

```math
s=\lfloor2\epsilon w\rfloor+1
```

in the Grassmann packing makes sharing an `epsilon w`-accurate message
impossible.  Substitution in the existing greedy packing bound yields

```math
\log_2 K_{\rm closed}(\epsilon w)
\ge ((1/2-2\epsilon)^2-o(1))w^2
```

for fixed `epsilon<1/4`.  The projection construction with
`d=floor(2 epsilon w)` gives the stated upper coefficient.  Thus
PMQ.15--PMQ.16 are normalized correctly, but they do not match at positive
distortion.

This is a genuine validation outside the syndrome word-length model: the
observable is matroid residual rank and the exact response algebra is the
flat lattice.  It should nevertheless be described as a refinement of the
existing HRC matroid model, not as yet another independent model beyond it.

## 7. Minor corrections before promotion

No theorem correction is required.  Three editorial clarifications would
improve the final version.

1. Equation PMQ.9 is duplicated verbatim.
2. The sentence using `d=floor(2 eta)` should declare `eta>=0` and
   `floor(2 eta)<=w`; PMQ.8 already supplies the relevant macroscopic range
   `0<epsilon<1/2`.
3. In the final PMQ.2 characterization, “the residual-rank oscillation” of
   the zero-separating quotient should explicitly mean the oscillation of
   **every** congruence class.

These are endpoint/wording issues only.  They do not alter any displayed
rate or proof.
