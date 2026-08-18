# Independent audit: growing-Renyi Hadamard cluster edge counterexample

**Frozen source:**
`extremal_information/drafts/growing_renyi_hadamard_cluster_edge_counterexample.md`

**SHA-256:**
`3e6b5b8e758d5f6a46dac8c5640cb32f41a3030c3c715a492fcf8b30e40ecb5a`

**Verdict:** **PASS.**  The construction is a valid exact-sign, central,
identity-second-moment law; the Renyi-two formula and its asymptotic constant
are correct; the occupancy lemma remains true after an arbitrary
deterministic deletion of `o(r)` latent directions; and every deterministic
rank-`o(r)` column projection leaves operator norm at least
`(L-o_Pr(1))sqrt(r)`.

No source repair is required.  One sentence attributes the rowwise union
bound merely to `c_r^2 r -> infinity`; that condition alone would not be
enough.  The displayed value in HC.1 satisfies the stronger estimate needed,
and the source's displayed failure bound is indeed `o(1)`, so the proof
itself is complete.

## 1. Exact density and second moment

For component `(j,s)`, coordinate `a` has probability

```math
\Pr(R_a=x\mid j,s)={1+s c_rh_{j,a}x\over2}.
```

Dividing by the uniform coordinate mass `1/2` and multiplying independent
coordinates gives HC.3 exactly.  Replacing `x` by `-x` replaces `s` by
`-s`, so the equal signed mixture is central and has zero mean.

For `a != b`, conditional independence gives

```math
\mathbb E(R_aR_b\mid j,s)=c_r^2h_{j,a}h_{j,b}.
```

The diagonal is one.  Since Hadamard column orthogonality says

```math
{1\over r}\sum_jh_{j,a}h_{j,b}=\mathbf1_{a=b},
```

averaging over `j` gives `E RR^T=I` exactly.  There is no population spike
and the canonical mean/covariance peel is empty.

## 2. Renyi-two formula and leading constant

Uniform-cube coordinate independence yields

```math
\langle g_{j,s},g_{k,t}\rangle_U
=\prod_a(1+stc_r^2h_{j,a}h_{k,a}).
```

For `j=k`, the two equal-sign ordered pairs contribute
`(1+c_r^2)^r` and the two opposite-sign ordered pairs contribute
`(1-c_r^2)^r`.  Across all `j`, their contribution after division by
`(2r)^2` is

```math
{(1+c_r^2)^r+(1-c_r^2)^r\over2r}.
```

For `j != k`, the Hadamard product has `r/2` plus and `r/2` minus entries.
All four sign pairs therefore have inner product `(1-c_r^4)^(r/2)`.
There are `4r(r-1)` such ordered component pairs, giving

```math
{r-1\over r}(1-c_r^4)^{r/2}.
```

This proves HC.6 with no missing factor two.

With

```math
c_r^2=L^2{\log\log r\over\log r},
```

one has

```math
r\log(1+c_r^2)
=L^2{r\log\log r\over\log r}
+O\left({r(\log\log r)^2\over(\log r)^2}\right).
```

The error, `log(2r)`, the decaying `(1-c_r^2)^r` term, and the at-most-one
cross term are all negligible on the leading scale.  Hence the exact
coefficient in HC.7 is `L^2`, and `log K_(2,r)=o(r)`.

## 3. Occupancy after an arbitrary deterministic deletion

The source uses the unsigned labels `J_i`, so this is `r` balls in `r` bins,
not `r` balls in `2r` bins.  This is correct: the later singular-value lower
bound squares row responses, so the random component signs do not cancel.

For

```math
t=(1-\epsilon){\log r\over\log\log r}+O(1),
```

the one-bin exact probability is

```math
{r\choose t}r^{-t}(1-r^{-1})^{r-t}
=r^{-(1-\epsilon)+o(1)}.
```

Thus any deterministic `G_r` with `|G_r|=r-o(r)` has

```math
\mathbb EZ_r=r^{\epsilon+o(1)}.
```

For two distinct retained bins, the exact joint probability is HC.12.
Dividing it by the square of the one-bin probability gives `1+o(1)` because
`t=o(sqrt r)`.  Therefore

```math
{\operatorname{Var}Z_r\over(\mathbb EZ_r)^2}
\le {1\over\mathbb EZ_r}+o(1)=o(1).
```

The second-moment argument is uniform in the deterministic set: it uses
only its cardinality.  For every fixed `epsilon>0`, a retained bin of load
at least `(1-epsilon)log(r)/loglog(r)` exists with probability tending to
one.  This is exactly the `1-o_Pr(1)` lower bound in HC.10.

## 4. An arbitrary deterministic rank-`o(r)` projection leaves enough basis vectors

Let `P` have rank `k=o(r)` and `V=I-P`.  For the orthonormal Hadamard basis,

```math
\sum_j\|Pu_j\|_2^2=\operatorname{tr}P=k.
```

Choose `eta=sqrt(k/r)` when `k>0` (and any vanishing positive value when
`k=0`).  Markov counting gives

```math
\#\{j:\|Pu_j\|_2^2>\eta\}
\le{k\over\eta}=\sqrt{kr}=o(r).
```

For every remaining `j`,

```math
a_j=\|Vu_j\|_2\ge\sqrt{1-\eta}=1-o(1).
```

This verifies the arbitrary-`P` residual-basis step.  It requires no
commutation of `P` with the Hadamard basis.  Since `P` is deterministic, the
retained set is eligible for HC.3 before the latent occupancies are sampled.

For `w_j=Vu_j/a_j`, orthogonal projection gives `Vw_j=w_j`.  Moreover

```math
\langle h_j,w_j\rangle
=\sqrt r{\langle u_j,Vu_j\rangle\over a_j}
=\sqrt r{\|Vu_j\|_2^2\over a_j}
=\sqrt r\,a_j.
```

Thus HC.18 has the correct residual factor `a_j`, rather than `a_j^2` or
its reciprocal.

## 5. Conditional mean, noise, and the union bound

Conditional on `(J_i,S_i)=(j,s)`,

```math
\mathbb E\langle R_i,w_j\rangle
=s c_r\sqrt r\,a_j.
```

The centered coordinates are independent, and the squared coefficients of
`w_j` sum to one.  Hoeffding therefore gives, uniformly in `j,s,P`,

```math
\Pr\left(
 |\langle R_i,w_j\rangle-s c_r\sqrt r\,a_j|
 >\zeta c_r\sqrt r\,a_j
 \mid j,s\right)
\le2e^{-c\zeta^2c_r^2r a_j^2}.
```

On the good set `a_j=1-o(1)`.  With

```math
\zeta_r=(c_r^2r)^{-1/4},
```

the exponent is of order

```math
\sqrt{c_r^2r}
=L\sqrt{{r\log\log r\over\log r}}
\gg\log r.
```

Hence the source's failure bound

```math
2r\exp\{-c\sqrt{c_r^2r}\}=o(1)
```

is correct.  This verifies simultaneous control for every sampled row whose
own latent label lies in the good set.  The random maximally occupied label
can then be selected without a post-selection gap.

The conditional signs `S_i` cause no cancellation: HC.19 is an absolute
value statement, and the operator lower bound uses the Euclidean square sum
of the row responses.

## 6. Final edge constant

For the retained label of occupancy

```math
N_j\ge(1-o_\Pr(1)){\log r\over\log\log r},
```

restricting `B_rw_j` to those rows gives

```math
\begin{aligned}
{\|B_rV_r\|_{op}^2\over r}
&\ge N_j(1-o(1))c_r^2a_j^2\\
&\ge(1-o_\Pr(1))
 {\log r\over\log\log r}
 L^2{\log\log r\over\log r}\\
&=L^2-o_\Pr(1).
\end{aligned}
```

Taking square roots proves HC.14 with the claimed exact lower constant `L`.
The estimates above are uniform for each deterministic projection sequence
of rank `o(r)`.  They do not claim simultaneous survival against a
sample-dependent projection, and the source explicitly excludes that
tautological setting.

## 7. Scope of the counterexample

The theorem rigorously falsifies a sharp-edge extension under the broad
condition `log K_(2,r)=o(r)`, even after deterministic sublinear-rank
peeling.  It does not itself produce a low-pressure conference bridge; a
large singular direction could instead increase pressure.  The source
states that limitation correctly.
