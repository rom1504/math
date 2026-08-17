# Independent audit: deep-hole and farthest-shell consequences

**Scope.**  This audit checks
[`nearmin_deep_hole_shell_report.md`](nearmin_deep_hole_shell_report.md)
against the proved FB.1/FB.3 statements and against the hypotheses of the
cited coding results.  It makes no canonical edit.

## Overall verdict

The two new propositions pass.  Their constants and asymptotic scales are
valid as stated:

| Claim | Verdict | Exact scope |
|---|---|---|
| signed-barycenter overlap identity and near-half-edge diameter | **PASS** | theorem-level consequence of the FB.6 clause of FB.1; it is more than cardinality but is not an absolute-overlap or response theorem |
| `exp(Omega(r/n))` packing at `Omega(r/log n)` Hamming separation | **PASS** | theorem-level consequence of FB.3 when `r/(n log n) -> infinity`; it is a carrier-metric packing, not a contextual packing |
| augmented-cut/deep-hole dictionary | **PASS** | exact, including the factor `1/2` between energy deficit and excess code radius |
| UPWS shortcut | **CORRECTLY REJECTED** | the exact equality hypothesis fails by `Omega(n^(3/2))`; no UPWS conclusion is imported |
| `n=7` finite example | **PASS, independently reproduced** | `rho=6`, `s(C_7^perp)=8`, and deep-coset leader multiplicities are exactly `3,4,7` |

No result in the report supplies contextual or compositional information in
the campaign's response sense.  DH.1 supplies genuine collective geometry,
and DH.2 supplies genuine Hamming metric entropy, so neither is *mere
cardinality*.  But neither declares continuations/queries whose response
values remain separated, controls selector products, or proves that a small
generative state is impossible.

## 1. Exact augmented-code dictionary

Write `N=binom(n,2)`, `y_e=(1-a_e)/2`, and
`c_e=(1-z_e)/2`.  Then exactly

```math
 \langle a,z\rangle=N-2d_H(y,c).
```

The set of such `c` is the cut space of `K_n` enlarged by the all-one word.
For `n>=3` it is a binary linear `[N,n]` code of size `2^n`.  Consequently

```math
 Q(a)=N-2d_H(y,C_n),
 \qquad M_n=N-2\rho(C_n).
```

An `epsilon`-near-minimizer has

```math
 d_H(y,C_n)\ge\rho(C_n)-{\epsilon\over2}n^{3/2},
```

and an energy shell `S_u(a)` is the code shell

```math
 d_H(y,c)\le d_H(y,C_n)+{u\over2}n^{3/2}.
```

Thus FB.3's energy deficit `2r` is exactly excess code radius `r`.  There is
no missing factor of two in the report.

The signed-graph citation is also mapped correctly.  Solé--Zaslavsky's
Lemmas 1 and 2 concern the ordinary cut/cocycle code: switching classes are
its cosets and ordinary frustration is coset minimum weight.  Enlarging by
the all-one word changes the distance to

```math
 \min\{\ell(\Sigma),\ell(-\Sigma)\},
```

so their ordinary covering-radius results cannot silently be substituted for
the augmented radius.

For `n>=5`, the reported parameters `[N,n,n-1]_2` are correct.  The dual is
exactly the set of Eulerian subgraphs with an even number of edges: cut-space
orthogonality imposes even degrees and all-one orthogonality imposes even
edge cardinality.

## 2. DH.1: overlap moment and diameter

Let

```math
 R(z,z')={1\over N}\sum_e z_ez'_e
         =1-{2d_H(c,c')\over N}
```

and let `Z,Z'` be independent draws from the FB.1 measure.  Since
`m_e=E[a_eZ_e]` and `a_e^2=1`,

```math
 E[Z_eZ'_e]=(E Z_e)^2=m_e^2.
```

Therefore the exact identity and bound are

```math
 0\le E R(Z,Z')={1\over N}\sum_e m_e^2
 \le {1\over N}\sum_e|m_e|\le\delta.             \tag{A.1}
```

Some support pair has `R<=E R`, hence

```math
 \operatorname{diam}_H(\operatorname{supp}\mu)
 \ge {1-\delta\over2}N.                           \tag{A.2}
```

The shell itself has at least this diameter.  For `delta<1` the selected pair
is automatically distinct.  The tail estimate in DH.3 also follows exactly
from `R>=-1`:

```math
 \Pr\{R>t\}\le {1+\delta\over1+t},\qquad -1<t<1.
```

The exact FB.1 value used here is

```math
 \delta={2(\epsilon+\eta_n(\kappa))
             \over\kappa(1-1/n)}.                 \tag{A.3}
```

For fixed sufficiently small `epsilon` and `kappa=sqrt(epsilon)`, the code
radius excess is `sqrt(epsilon)n^(3/2)` and

```math
 \delta={2\sqrt\epsilon\over1-1/n}+o_n(1).
```

Thus the displayed `(1-O(sqrt(epsilon)))N/2` diameter scale is correct (with
the usual order: fix `epsilon`, then let `n` grow).  For an exact minimizer
and `kappa=n^(-1/6)`,

```math
 \eta_n(\kappa)=C(n^{-1/3}+n^{-1/2}),
 \qquad
 \delta={2C(n^{-1/6}+n^{-1/3})\over1-1/n},
```

while the energy deficit is `2n^(4/3)` and the excess code radius is exactly
`n^(4/3)`.  Hence the claimed
`(1-O(n^(-1/6)))N/2` diameter is correct.

**Classification.**  This is a theorem-level corollary/proposition, not an
empirical observation.  It gives a replica-overlap first moment and an
existential far pair.  Calling it a "signed overlap moment" is safe only if
the text retains the definition: it does **not** control `E|R|`, a vertex
spin overlap, selector correlations, or the distribution of `R`.  A slightly
less ambiguous phrase would be "replica edge-overlap moment induced by the
signed shell barycenter."

## 3. DH.2: separated FB.3 packing

Put `s=c_1r/n` and choose a `sigma` fibre containing at least `exp(s)/2`
shell words.  In that fibre, projective vertex distance `k<=n/2` gives the
exact edge distance

```math
 d_H(c,c')=k(n-k).                                 \tag{A.4}
```

Let `L=log(en)` and `h=floor(s/(8L))`.  Since
`s/L -> infinity`, eventually `h>=1`, `h=o(n)`, and a projective radius-`h`
ball has size

```math
 \sum_{j=0}^h\binom nj
 \le \exp(hL)\le\exp(s/8).                         \tag{A.5}
```

Greedy deletion therefore retains at least

```math
 {e^s\over2e^{s/8}}={e^{7s/8}\over2}\ge e^{s/2}  \tag{A.6}
```

points for all sufficiently large `n`.  Selected projective distances exceed
`h`.  Once `s/(8L)>=2`, `h>=s/(16L)`, so (A.4) gives the slightly stronger
intermediate bound

```math
 d_H(c,c')\ge {hn\over2}
 \ge {c_1r\over32\log(en)}.
```

The report's denominator `40` is therefore safe.  The FB.3 assumption
`r=o(n^(3/2))` gives `s=o(sqrt n)`, which justifies all uses of `h=o(n)`.
The additional DH.2 hypothesis `s/log n -> infinity` is equivalent, up to
the fixed positive `c_1`, to

```math
 {r\over n\log n}\longrightarrow\infinity.
```

For `r=floor(n^(3/2)/log n)`, this yields exactly

```math
 |P|\ge\exp\!\left(\Omega\left({\sqrt n\over\log n}\right)\right),
 \qquad
 d_{\min}(P)=\Omega\!\left({n^{3/2}\over\log^2 n}\right).
```

The claim that the logarithmic loss is intrinsic should be read as a
worst-case statement **from cardinality alone**.  Indeed a projective Hamming
ball of radius `Theta(s/log(n/s))` already has `exp(Theta(s))` points, and in
the FB.3 regime `s=o(sqrt n)` one has `log(n/s)=Theta(log n)`.  Such a set has
edge diameter only `O(ns/log n)=O(r/log n)`, so FB.3's count cannot by itself
force `Omega(r)` edge separation.

**Classification.**  DH.2 is theorem-level and genuinely stronger than
FB.3 cardinality: it proves raw Hamming metric entropy
`log|P|=Omega(r/n)`.  It is nevertheless not a contextual response packing.
The metric is distance between codewords themselves, not a supremum over
allowed continuations, and every word still has the ordinary `O(n)` vertex
label description `(sigma,[x])`.  No independent writable coordinates or
sublinear composable response state follows.

## 4. UPWS and external-distance audit

The published implication is used with the correct hypotheses and only as a
rejected route:

```text
C is UPWS  <=>  rho(C)=s(C^perp)
              => all weight-rho cosets have one full weight distribution
              => equal nearest-leader multiplicity at deep holes.
```

Here `s(C^perp)` is the number of distinct nonzero weights in `C^perp`.
The equality, not approximate equality and not merely linearity, is essential.

The report's maximum dual weights are correct:

```math
W_n=
\begin{cases}
N-n/2,&n\text{ even},\\
N,&n\equiv1\pmod4,\\
N-3,&n\equiv3\pmod4.
\end{cases}                                             \tag{A.7}
```

For even `n`, the complement of an Eulerian graph has all degrees odd, so it
has at least `n/2` edges, attained by a perfect matching; the parity condition
is then also satisfied.  For odd `n`, `K_n` is Eulerian.  It has even size in
the `1 mod 4` case, while in the `3 mod 4` case the smallest odd-size
Eulerian deletion is a triangle.

The Bryant--Horsley--Pettersson application is valid, with one congruence
detail that should be made explicit in a polished proof.  Their cycle
decomposition theorem applies to `K_n` for odd `n` and to `K_n` minus a
perfect matching for even `n`, for every list of cycle lengths in `[3,n]`
having the required total.  For every even `w` with `4<=w<=W_n-4`, write
both `w` and `W_n-w` as sums of integers in `[3,n]`.  If `n` is even or
`n=1 mod 4`, decompose the relevant base graph according to those two
concatenated sublists.  If `n=3 mod 4`, append one additional `3` to the
list before decomposing `K_n`; that triangle accounts for
`N-W_n=3`.  Taking the cycles in the `w` sublist gives an even-edge Eulerian
graph of weight `w`, while taking all cycles except the extra triangle also
realizes `W_n`.  Hence

```math
 s(C_n^perp)\ge W_n/2-2=N/2-O(n).                 \tag{A.8}
```

On the other hand `M_n>=c_0n^(3/2)` gives

```math
 \rho(C_n)={N-M_n\over2}
 \le {N\over2}-{c_0\over2}n^{3/2}.                \tag{A.9}
```

Combining (A.8)--(A.9),

```math
 s(C_n^perp)-\rho(C_n)=\Omega(n^{3/2})>0
```

eventually.  Thus the code is not UPWS and the equal-deep-coset theorem is
unavailable.  No conclusion in the report relies on a nonexistent
"near-UPWS" stability theorem.

There is one bibliographic labeling nuance, not a mathematical defect:
Theorem 1.4 of the cited 2025 Davydov--Marcugini--Pambianco manuscript is a
modern restatement that itself attributes the UPWS equivalence to earlier
papers.  It is primary for its own Theorem 5.2 formulation, but should not be
described as the original primary source of the equivalence.  The report's
separate attribution of the external-distance bound to Delsarte is correct.

## 5. Other coding-literature mappings

The remaining nonapplications are correctly scoped.

1. **MCF/APMCF.**  These notions assume or certify multiplicity at the actual
   covering radius.  FB.3 counts the interval of coset weights
   `rho,...,rho+r`; it gives no lower bound at weight exactly `rho`.
2. **MDS deep cosets.**  The quoted theorem requires an MDS code with
   `R=d-1`.  Here `d=n-1`, the Singleton value would be `N-n+1`, and
   `rho=N/2-Theta(n^(3/2))`; neither hypothesis is remotely present.
3. **Average-coset results.**  The augmented cut dual contains a 4-cycle, so
   the large bilateral-dual-distance route is absent.  More importantly, an
   average over cosets cannot be conditioned on the rare deepest cosets
   without an additional theorem.
4. **Coset-leader algorithms.**  An output-sensitive enumeration method is
   not a multiplicity, overlap, or description-length bound.

The independently repeated `n=7` enumeration found `|C_7|=128`, `16384`
cosets, covering radius `6`, dual weight support

```text
0,4,6,8,10,12,14,16,18,
```

and hence external distance `8`.  The deep cosets have precisely the three
nearest-leader multiplicities `3,4,7`; the three representatives printed in
the report realize those values.  This validates the finite obstruction but,
as the report says, does not replace the asymptotic non-UPWS proof.

## 6. Theorem-status and information-status decision

The report can retain DH.1 and DH.2 as named propositions.

| Result | Theorem status | Information supplied | What remains absent |
|---|---|---|---|
| DH.1 | rigorous theorem-level corollary of FB.1 | collective replica moment and a macroscopic far pair | absolute overlap, selector/Fourier control, continuation stability |
| DH.2 | rigorous theorem-level corollary of FB.3 | `Omega(r/n)` raw carrier bits at `Omega(r/log n)` Hamming resolution | response-query separation, independent bits, composable state lower bound |
| non-UPWS gap | rigorous asymptotic obstruction | rules out one published multiplicity shortcut | says nothing positive about actual deep-hole multiplicities |
| `n=7` enumeration | exact finite computation | disproves automatic constant multiplicity from symmetry | no asymptotic classification |

Accordingly, the bottom-line sentence that neither conclusion is yet a
composable response state is correct.  DH.1 is valuable Level-5 necessary
geometry, and DH.2 is a valid packing baseline, but neither crosses the
contextual/compositional frontier.
