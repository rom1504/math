# Independent audit: normalized syndrome-response rate--distortion

**Scope.** This report audits
`phase2_normalized_code_rate_distortion.md`. It reconstructs the block
fragments, NRD.1--NRD.3, and the Shannon consequence. It makes no surface
edits.

## Verdict

**Accept after one interpretation correction.** The common-length/full-rank
construction is valid, the response formula and induced metric are exact,
the Gilbert--Varshamov constants give the stated `epsilon<1/8` regime, and
the deterministic and Fano bounds have the right factors.

The important correction is about the matching upper bound and “strict
quotient” language:

- the `q=w/L`-bit exact upper bound holds for the **restricted block source
  family** `(H_a)_a`, while answering **all unrestricted appended-fragment
  queries** on that family;
- it is not an `O(w)` upper bound for arbitrary syndrome supports in
  `F_2^w`; and
- if one canonical fragment is fixed for each `a`, then `a` identifies that
  fragment and can generate its whole root-distance table. The normalized
  lower bound alone does not prove strict nonreconstruction on this canonical
  source. Strict forgetting follows only after the source is explicitly
  enlarged to multiple multiplicity patterns with the same support, or from
  the earlier CR.1 examples.

With that wording repaired, the result merits concise promotion as one
normalized code rate--distortion theorem. It is a genuine advance over the
raw `<1/2` lattice-scale packing, although its proof is an elementary direct-
sum/GV construction and its scope remains the structured block family.

## 1. Common length and full rank

For `w=Lq`, write

```math
G=V_1\oplus\cdots\oplus V_q,
\qquad \dim V_j=L.
```

Every sparse support `B_j` and every dense support
`D_j=V_j\setminus\{0\}` contains the fixed basis `B_j`. Their union over all
blocks therefore spans `G`, so every `H_a` has row rank `w`.

The common per-block length can be `2^L-1`: use every type once in a dense
block and repeat arbitrary basis types in a sparse block until that length is
reached, while keeping every basis type. Duplicate binary columns do not
change the coset-leader profile. Hence every state fragment has the common
length

```math
n=q(2^L-1).
```

The environment support `T_P` also contains a basis in every direct summand,
whether its block is `B_j` or `D_j`. Thus every `E_P` is full rank. The same
padding makes all environments length `n`. For fixed `L` this is
`Theta_L(w)`, and after choosing `L=L(epsilon)` it is
`Theta_epsilon(w)` exactly as claimed.

## 2. Exact response formula and metric

For a support that splits across the direct sum, every representation of
`s=(s_1,...,s_q)` splits uniquely into block representations. Therefore

```math
\lambda(s)=\sum_j\lambda_j(s_j),
\qquad
\rho=\sum_j\rho_j.
```

A basis of `F_2^L` has covering radius `L`, attained by the all-basis-vector
sum, while the complete nonzero support has radius one. In the composite
support `S_a union T_P`, block `j` remains only a basis precisely when
`j in P` and `a_j=0`; every other block contains `D_j`. Consequently

```math
\mathcal R_{H_a}(E_P)
=q+(L-1)|\{j\in P:a_j=0\}|,
```

so NRD.7 is exact.

Let `Z_a={j:a_j=0}`. After subtracting two responses and dividing by `L-1`,
the restricted response distance is

```math
\max_P\bigl||P\cap Z_a|-|P\cap Z_b|\bigr|.
```

If

```math
A=Z_a\setminus Z_b,
\qquad B=Z_b\setminus Z_a,
```

then the maximum is `max{|A|,|B|}`: choosing `P=A` or `P=B` attains the two
directed counts, and no subcounts can have a larger absolute difference.
This proves NRD.8. Since `|A|+|B|=d_H(a,b)`, NRD.9 follows with the necessary
factor `1/2`.

The lower bound uses only the valid block environments `E_P`; it therefore
applies a fortiori when the decoder must answer every appended fragment.

## 3. Deterministic packing and constants

A maximal binary code of minimum distance at least `d` has Hamming balls of
radius `d-1` covering the cube. Hence the greedy/GV bound used in NRD.18 is

```math
|\mathcal A|\ge
{2^q\over\sum_{i=0}^{d-1}\binom qi}.
```

Two packing words differ in at least `d` positions, so one directed
difference has size at least `ceil(d/2)`. NRD.8 therefore gives response
separation

```math
\Delta_*=(L-1)\left\lceil{d\over2}\right\rceil.
```

If two sources shared one deterministic summary, their common decoded
response would put their true responses at distance at most `2 eta`. Thus
the strict condition `2 eta<Delta_*` forces injectivity on the packing and
gives NRD.16. Both the ceiling and factor two are correct.

Take `d=ceil(delta q)` with fixed `0<delta<1/2`. Then

```math
\log_2\sum_{i=0}^{d-1}\binom qi
\le q(h_2(\delta)+o(1)),
```

and the summary lower bound is

```math
\left({1-h_2(\delta)\over L}-o(1)\right)w.
```

For `eta=epsilon w=epsilon Lq`, the asymptotic separation condition is

```math
2\epsilon Lq<(L-1){\delta q\over2},
```

equivalently

```math
\epsilon<{\delta(L-1)\over4L}.
```

The strict margin handles the ceilings. There exists a
`delta<1/2` satisfying this inequality precisely when

```math
{4\epsilon L\over L-1}<{1\over2},
```

or `L(1-8 epsilon)>1`. Thus for every `epsilon<1/8`, choosing an integer

```math
L>{1\over1-8\epsilon}
```

and then a `delta` in the displayed interval proves NRD.23--NRD.25. The
result holds along widths divisible by this fixed `L`, hence at infinitely
many widths, as stated. The threshold is a limitation of this binary block
packing, not a converse above `1/8`.

## 4. Exact upper bound and its scope

Given `a` and the declared block apparatus, the support `S_a` is known.
The earlier syndrome-response algebra computes

```math
\rho(\ker[H_a\ E])
```

for every appended fragment `E` from `S_a` and the support of `E`; column
multiplicities are irrelevant. Thus the `q=w/L` latent bits give an exact
summary for **all queries on this block family**, not merely the test family
`(E_P)_P`. Combined with the lower bound, the deterministic response
complexity of this restricted source is `Theta(w)` in the stated distortion
regime.

Nothing here gives an `O(w)` summary for the full class of the
`2^(2^w-1)` possible syndrome supports. The draft's final open question
correctly leaves a subexponential-in-`2^w` approximate quotient for that
unrestricted class unresolved.

There is also a distinction between raw table size and information. If one
canonical padded matrix `H_a` is fixed per `a`, then the `q` bits determine
that finite landscape, however expensive its explicit table is to enumerate.
To make “strict quotient” literal, define the source family to include
different common-length multiplicity patterns with the same `S_a`; all such
fragments have identical future-radius responses while their codes/root
tables can differ. Otherwise replace “strict quotient” by “a succinct
generative response state for the restricted family.”

## 5. Shannon/Fano claim

The packing responses are separated by `Delta_*` in sup norm. Nearest-
neighbor decoding from `Rhat_Z` is correct whenever

```math
\|\widehat R_Z-\mathcal R_{H_A}\|_\infty
<\Delta_*/2.
```

Markov's inequality therefore gives

```math
p_e\le {2D\over\Delta_*}
```

under the expected sup-error hypothesis NRD.27. For a nontrivial packing,
Fano gives

```math
I(A;Z)
\ge\log_2|\mathcal A|-h_2(p_e)
-p_e\log_2(|\mathcal A|-1).
```

When `p=2D/Delta_*<=1/2`, the Fano penalty is increasing on the relevant
range, so replacing `p_e` by `p` yields NRD.29. A positive-rate packing and
any fixed `p<1` give a linear lower bound; the stated `p<=1/2` condition is a
clean sufficient regime.

This is an expected **global sup-response** statement. It does not imply an
information lower bound for mean squared error under a diffuse distribution
on individual environments, and the draft states that exclusion correctly.

## 6. Promotion recommendation

Promote one concise theorem after the exact syndrome algebra:

> For each fixed `epsilon<1/8`, there are fixed `L` and
> `c_epsilon>0` and, for every sufficiently large `w` divisible by `L`, a
> common-length full-rank block family such that uniform error `epsilon w`
> on all appended-fragment radius queries costs at least `c_epsilon w` bits,
> while `w/L` bits answer every query exactly on that family.

Keep the subset-count derivation, detailed GV optimization, and Fano variant
in the draft. Surface wording must say “restricted source family, unrestricted
query family” and must not extrapolate the `O(w)` upper bound to arbitrary
syndrome supports. With those qualifications, this is the first promoted code
result in the program that survives fixed normalized additive distortion.
