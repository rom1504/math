# A diffuse third witness from a mesoscopic bipartite agreement core

Date: 2026-08-17.

Status: **proved draft; no canonical edit**.

This note attacks the surviving PP.4 branch.  Exact multi-edge optimality
does force a third near-ground projective direction.  A localized entropy
argument makes that direction diffuse on a constant fraction of the common
correct part of the bipartite core, while paying only `O(n log^2 n)` shell
deficit.  The resulting projective separation is
`(1/2-o(1))M_n`, however, not a fixed fraction of `binom(n,2)`.  The proof
also identifies why this architecture stops there: its complementary-distance
estimate is exactly AO.20.

## 1. The common-correct reservoir

Let `A` be an exact minimizer of order `n`, put

```math
E={n\choose2},\qquad M=Q(A)=M_n,
```

and use the PP.4 gauge in which

```math
z_0=\mathbf1,
\qquad z_1=-c(v_S)
```

are positive signed-cut words of deficits at most `2s`.  Thus

```math
\langle a,z_j\rangle\ge M-2s\qquad(j=0,1).         \tag{MB.1}
```

Let `P=delta(S)`, `D=|P|=|S|(n-|S|)`, and

```math
C=\sum_{e\in P}a_e\ge M-2s.                       \tag{MB.2}
```

Both poles equal `+1` on `P`.  Split that core into

```math
Z=\{e\in P:a_e=+1\},
\qquad B=\{e\in P:a_e=-1\},
\qquad p=|Z|={D+C\over2}.                          \tag{MB.3}
```

In the code dictionary, `Z` is exactly the set of core coordinates correct
for both error supports, while `B` is wrong for both.  Since `D>=C>=0`,

```math
\boxed{p\ge C\ge M-2s.}                            \tag{MB.4}
```

PP.4 is used below only through this large common-correct reservoir.  The
collapsed augmented-cut branch supplies the reservoir, but the argument
below in fact uses only its size and not the complete-bipartite identity.
No arbitrary rectangular bridge is introduced or paid separately.

## 2. Localized multi-edge optimality

For an augmented cut `z`, write

```math
d_A(z)=M-\langle a,z\rangle
```

for its positive-shell deficit.  The following exact observation is the
multi-edge engine.

### Lemma MB.1 (flipping common-correct coordinates)

For every `F subseteq Z` with `|F|=r`, there is an augmented cut `z` such
that, with `d=d_A(z)`,

```math
0\le d\le2r,
\qquad
|\{e\in F:z_e=-1\}|\ge {r\over2}+{d\over4}.        \tag{MB.5}
```

#### Proof

Flip exactly the edges in `F`, obtaining `A^F`.  Exact minimality of `A`
among all order-`n` signings gives `Q(A^F)>=M`.  Choose the orientation of
an augmented cut `z` so that

```math
\langle a^F,z\rangle\ge M.
```

Every edge of `F` has `a_e=+1`.  If `q` of their `z_e` values are `-1`,
then

```math
M\le\langle a^F,z\rangle
=M-d-2\sum_{e\in F}z_e
=M-d-2r+4q.                                        \tag{MB.6}
```

Hence `q>=r/2+d/4`; since `q<=r`, also `d<=2r`. `square`

For one prescribed `F`, Lemma MB.1 gives only `r/2` disagreements.  The
point of the next theorem is that `F` can be chosen before the response so
that no augmented cut concentrated on a small part of `Z` can supply those
disagreements.

## 3. The third-witness theorem

### Theorem MB.2 (localized entropy forces a diffuse third witness)

Fix `0<theta<1/2`.  Suppose `1<=r<=p`, `2r<M`, and

```math
2^n\exp\{-2(1/2-\theta)^2r\}<1.                   \tag{MB.7}
```

Then there is a positive signed-cut word `z` of deficit `d<=2r` such that

```math
|\{e\in Z:z_e=-1\}|>\theta p.                     \tag{MB.8}
```

Moreover, simultaneously for `j=0,1`,

```math
\boxed{
d_{\rm P}(z,z_j)
\ge
\min\left\{\theta p, M-s-{d\over2}\right\}
\ge
\min\{\theta p, M-s-r\}.}                       \tag{MB.9}
```

Consequently the three words `z_0,z_1,z` lie in the common deficit
`2max{s,r}` shell and satisfy

```math
d_{\rm P}(z_0,z_1)\ge M-2s,
\qquad
d_{\rm P}(z,z_j)\ge
\min\{\theta(M-2s),M-s-r\}.                       \tag{MB.10}
```

The hypothesis `2r<M` ensures that the response orientation constructed in
Lemma MB.1 has strictly positive original energy, so no reorientation caveat
is needed.

#### Proof

Choose `F` uniformly among the `r`-subsets of `Z`.  Fix an augmented cut
`w` having at most `theta p` negative coordinates on `Z`.  Hoeffding's
inequality for hypergeometric sampling gives

```math
\Pr_F\{|F\cap\{e\in Z:w_e=-1\}|\ge r/2\}
\le\exp\{-2(1/2-\theta)^2r\}.                      \tag{MB.11}
```

There are at most `2^n` augmented cuts.  By (MB.7) and a union bound, some
`F` has the property that **every** augmented cut with at most `theta p`
negative core coordinates has fewer than `r/2` of them in `F`.

Apply Lemma MB.1 to this `F`.  Its response word has at least `r/2+d/4`
negative coordinates in `F`, so it cannot be one of the concentrated words.
This proves (MB.8) and `d<=2r`.

Let `h_j=d_E(z,z_j)` be the actual signed-edge distance.  Both poles are
`+1` on `Z`, so (MB.8) gives

```math
h_j>\theta p.                                      \tag{MB.12}
```

On the other hand, two positive energy words can collect their summed
energy only on coordinates where they agree:

```math
\langle a,z\rangle+\langle a,z_j\rangle
\le2(E-h_j).
```

Using (MB.1) and `\langle a,z\rangle=M-d` yields

```math
E-h_j\ge M-s-d/2.                                  \tag{MB.13}
```

Since `d_P(z,z_j)=min{h_j,E-h_j}`, equations
(MB.12)--(MB.13) prove the first bound in (MB.9).  The second uses `d<=2r`,
and (MB.4) gives (MB.10).  Finally PP.4 already gives
`d_P(z_0,z_1)=D>=M-2s`. `square`

### Corollary MB.3 (an `O(n)`-deficit constant and an asymptotically sharp form)

Assume the intended PP.4 regime

```math
M\ge c_0n^{3/2},\qquad s=o(M).                     \tag{MB.14}
```

1. Taking `theta=1/4` and

   ```math
   r=\lfloor8n\log2\rfloor+1                     \tag{MB.15}
   ```

   gives, for all sufficiently large `n`, a third word of deficit `O(n)`
   and projective distance at least

   ```math
   (1/4-o(1))M                                    \tag{MB.16}
   ```

   from both original poles.

2. More sharply, take

   ```math
   \theta_n={1\over2}-{1\over\log n},
   \qquad
   r_n=\left\lceil{(n\log2+1)\log^2n\over2}\right\rceil. \tag{MB.17}
   ```

   Then `r_n=O(n log^2n)=o(M)`, condition (MB.7) holds, and the third word
   has deficit `O(n log^2n)=o(M)` with

   ```math
   \boxed{
   d_{\rm P}(z,z_j)\ge(1/2-o(1))M
   \quad(j=0,1).}                                  \tag{MB.18}
   ```

Thus a mesoscopic two-cap core cannot consist of only two isolated poles:
an asymptotically vanishing shell enlargement forces a third direction at
the full `M_n` edge scale.

#### Proof

For fixed `theta=1/4`, the exponent in (MB.7) is `r/8`; (MB.15) makes it
strictly larger than `n log2`.  Equations (MB.4), (MB.9), and (MB.14) give
(MB.16).  For (MB.17), the exponent in (MB.7) is
`2r_n/log^2n>n log2`, while `r_n=o(M)` and `p>=M-2s`; (MB.9) gives
(MB.18). `square`

## 4. Signed interface diffusion versus projective diffusion

The theorem actually proves the stronger signed-reservoir statement

```math
d_E(z,z_j)\ge\theta p
\ge{\theta\over2}(D+M-2s).                         \tag{MB.19}
```

Hence, if the PP.4 shore has `|S|ggsqrt n`, the new response disagrees with
both poles on a macroscopic fraction of the entire bipartite interface,
potentially much more than `M` edges.  Projectivization loses that gain:
the response is allowed to lie close to the opposite signed lift, and
(MB.13) protects only `M-s-d/2` agreement edges.

This is not a technical loss in the displayed proof.  Equation (MB.13) is
exactly the AO.20 positivity inequality applied to `z` and `z_j`.  In the
intended `s=o(M)` regime (and up to integer rounding), optimizing the
elementary one-set version of Lemma MB.1 gives only

```math
\max_r\min\{r/2,M-s-r\}={M-s\over3},               \tag{MB.20}
```

whereas the localized entropy argument upgrades its first term to roughly
`p/2` and reaches the ceiling `(1/2-o(1))M`.  Neither estimate can turn the
second term into `Theta(E)`.

Since `M=Theta(n^(3/2))` and `E=Theta(n^2)`, (MB.18) corresponds to AO
overlap gap

```math
{2d_{\rm P}\over E}=(1-o(1)){M\over E}
=Theta(n^{-1/2}),                                  \tag{MB.21}
```

not fixed `gamma`.  Even when `D` is much larger than `M`, the present
multi-edge/positivity architecture therefore stalls at the energy scale.
A fixed-projective-gap theorem needs a new fact ruling out a fresh
near-antipodal lift of each new response, or a recursive argument showing
that the resulting succession of agreement cores cannot persist.  The
displayed one-reservoir entropy argument plus AO.20 has reached its
certifiable ceiling; this is not a no-go theorem for every nonrecursive
multi-edge argument.

## 5. Is this FB.3 or AO.20 in disguise?

The answer is precise.

* **MB.1 is the local form of the exact certificate underlying FB.3.**  In
  code notation, for every `r`-set `F` a closest word to the perturbed deep
  hole has error support `N` satisfying

  ```math
  |F\cap N|\ge {r+d_A(z)/2\over2}.
  ```

  FB.3 applies a large-deviation union bound over a hypothesized small
  shell to force shell cardinality.  MB.2 instead applies the union bound
  over **all** `2^n` augmented cuts, but only after PP.4 supplies a reservoir
  `Z` of size at least `M-2s`.  It forces one geometrically diffuse witness.
  Thus MB.2 is a new localized geometric consequence of the same primitive,
  not an independent optimality principle and not merely FB.3's cardinality
  statement.

* **The complementary half of MB.9 is exactly AO.20.**  It is the sole
  reason signed diffusion becomes only `Theta(M)` projective diffusion.
  Therefore this result refines PP.4 and FB.3 but does not evade the existing
  absolute-overlap ceiling.

* **No arbitrary bridge is paid.**  The perturbation flips actual
  common-correct edges of the exact signing, and the witness is selected by
  exact order-`n` minimality.  After PP.4 provides the two poles, the proof
  uses no cut-specific fact: the same implication holds for any response
  family of size at most `2^n` with a common-correct reservoir of this size.

The result is consequently a rigorous conditional third-witness theorem and
a sharp audit of what the generic flip/one-reservoir/AO.20 ingredients alone
certify.  The generic PP.2 countermodel prevents upgrading those premises to
fixed scale.  A positive upgrade must add a cut-specific fact or recursive
control of the new near-antipodal agreement cores.
