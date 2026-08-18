# A subexponential-Renyi isotropic row law with a super-Bai--Yin edge

**Status.**  Task-local counterexample theorem.  This shows that the
bounded-`L2` row theorem cannot be extended to arbitrary
`K_r=exp(o(r))`.  The law is central and has second moment exactly `I`, yet
its iid `r`-row matrix has a fixed singular-value excess above `2sqrt(r)`.
The excess survives every deterministic `o(r)`-rank column peel.

This is a random-matrix obstruction to the projected-coupling route.  It
does not assert that the resulting bridge law lowers conference pressure.

## 1. Construction

Let `r=2^m` tend to infinity through Walsh--Hadamard orders.  Let
`h_1,...,h_r in {+-1}^r` be the rows of a Hadamard matrix, so

```math
\langle h_j,h_k\rangle=r\mathbf1_{j=k},
\qquad u_j={h_j\over\sqrt r}

```

is an orthonormal basis.  Fix a constant `L>2` and put

```math
c_r=L\sqrt{{\log\log r\over\log r}}.
\tag{HC.1}
```

For all sufficiently large `r`, `0<c_r<1`.  Generate one exact-sign row
`R` as follows:

1. choose `J` uniformly from `[r]`;
2. choose `S` uniformly from `{+-1}`;
3. conditional on `(J,S)`, choose the coordinates independently with

   ```math
   \mathbb E[R_a\mid J,S]=Sc_rh_{J,a}.
   \tag{HC.2}
   ```

Equivalently, the component density relative to the uniform row cube is

```math
g_{j,s}(x)=\prod_{a=1}^r(1+sc_rh_{j,a}x_a),
\tag{HC.3}
```

and the row density is

```math
g_r(x)={1\over2r}\sum_{j=1}^r\sum_{s=+-1}g_{j,s}(x).
\tag{HC.4}
```

Let `B_r` have `r` independent rows with this marginal law.  One may realize
them using independent latent labels `(J_i,S_i)`; after the labels are
forgotten, the rows still have exactly the iid marginal (HC.4).

## 2. Exact mean and second moment

### Proposition HC.1 (the law is central and exactly isotropic in second moment)

```math
g_r(-x)=g_r(x),
\qquad \mathbb E R=0,
\qquad \boxed{\mathbb E RR^T=I_r.}
\tag{HC.5}
```

**Proof.**  Negating `x` interchanges `g_(j,+)` and `g_(j,-)`, proving
centrality and zero mean.  Conditional on `(j,s)`, coordinate independence
gives, for `a!=b`,

```math
\mathbb E[R_aR_b\mid j,s]
=c_r^2h_{j,a}h_{j,b},
```

while every diagonal entry is one.  Hadamard column orthogonality says

```math
{1\over r}\sum_{j=1}^rh_{j,a}h_{j,b}=\mathbf1_{a=b}.
```

Averaging over `j` proves (HC.5). `square`

Thus the mean/covariance spectral peel is empty at every threshold: the
counterexample is not a hidden population spike.

## 3. Exact Renyi-two cost

### Proposition HC.2 (subexponential but superpolynomial row density cost)

The squared density has the exact value

```math
\boxed{
K_{2,r}:=\mathbb E_Ug_r^2
={ (1+c_r^2)^r+(1-c_r^2)^r\over2r}
+{r-1\over r}(1-c_r^4)^{r/2}.}
\tag{HC.6}
```

Consequently

```math
\boxed{
\log K_{2,r}
=(L^2+o(1)){r\log\log r\over\log r}=o(r).}
\tag{HC.7}
```

In particular `K_(2,r)=exp(o(r))`, although it is not polynomial.

**Proof.**  Independence of the uniform cube coordinates gives

```math
\langle g_{j,s},g_{k,t}\rangle_U
=\prod_{a=1}^r(1+stc_r^2h_{j,a}h_{k,a}).
\tag{HC.8}
```

For `j=k`, this is `(1+c_r^2)^r` when `s=t` and
`(1-c_r^2)^r` when `s=-t`.  For `j!=k`, the sign vector
`h_jh_k` has equally many plus and minus coordinates, so every sign pair
gives `(1-c_r^4)^(r/2)`.  Summing the `4r^2` ordered component pairs proves
(HC.6).

Now `c_r^2->0`, `rc_r^2=L^2r loglog(r)/log(r)`, and
`rc_r^4=o(rc_r^2)`.  The first term in (HC.6) dominates on the logarithmic
scale, while `log r=o(rc_r^2)`, proving (HC.7). `square`

## 4. Occupancy creates the outlier

Let

```math
N_j=\#\{i:J_i=j\}
\tag{HC.9}
```

be the latent occupancy numbers.

### Lemma HC.3 (large occupancy persists after deleting `o(r)` bins)

If `G_r subseteq[r]` is any deterministic set with
`|G_r|=r-o(r)`, then

```math
\boxed{
\max_{j\in G_r}N_j
\ge(1-o_\Pr(1)){\log r\over\log\log r}.}
\tag{HC.10}
```

**Proof.**  Fix `epsilon in (0,1)` and put

```math
t_r=\left\lfloor{(1-\epsilon)\log r\over\log\log r}\right\rfloor,
\qquad
Z_r=\sum_{j\in G_r}\mathbf1_{N_j=t_r}.
```

The exact first moment is

```math
\mathbb EZ_r
=|G_r|{r\choose t_r}r^{-t_r}(1-r^{-1})^{r-t_r}
=r^{\epsilon+o(1)}.
\tag{HC.11}
```

For distinct bins the joint probability is

```math
{r!\over t_r!^2(r-2t_r)!}
r^{-2t_r}(1-2/r)^{r-2t_r}.
\tag{HC.12}
```

Since `t_r=o(sqrt r)`, the ratio of (HC.12) to the square of the one-bin
probability tends to one.  Hence

```math
{\operatorname {Var}Z_r\over(\mathbb EZ_r)^2}
\le {1\over\mathbb EZ_r}+o(1)\longrightarrow0.
```

Thus `Z_r>0` with probability tending to one.  Since `epsilon` was
arbitrary, (HC.10) follows. `square`

## 5. The excess survives every deterministic sublinear-rank peel

### Theorem HC.4 (super-edge after arbitrary deterministic `o(r)` rank)

Let `P_r` be any deterministic orthogonal projection with

```math
\operatorname {rank}P_r=o(r),
\qquad V_r=I-P_r.
\tag{HC.13}
```

Then

```math
\boxed{
\liminf_{r\to\infty}{\|B_rV_r\|_{op}\over\sqrt r}
\ge L
\quad\hbox{in probability}.}
\tag{HC.14}
```

In particular, choosing any fixed `L>2` gives a fixed Bai--Yin-edge excess
despite exact identity second moment and subexponential row Renyi cost.

**Proof.**  Write `k_r=rank(P_r)`.  Since `(u_j)` is an orthonormal basis,

```math
\sum_{j=1}^r\|P_ru_j\|_2^2=k_r.
\tag{HC.15}
```

Choose numbers `eta_r downarrow0` with `k_r/(eta_rr)->0`, for example
`eta_r=sqrt(k_r/r)` when `k_r>0`, with an arbitrary vanishing replacement
when `k_r=0`.  The good-direction set

```math
G_r=\{j:\|P_ru_j\|_2^2\le\eta_r\}
\tag{HC.16}
```

has size `r-o(r)`.  For `j in G_r`, put

```math
a_j=\|V_ru_j\|_2\ge\sqrt{1-\eta_r}=1-o(1),
\qquad
w_j={V_ru_j\over a_j}.
\tag{HC.17}
```

Conditional on a row label `(J_i,S_i)=(j,s)`, coordinate independence and
(HC.2) give

```math
\mathbb E[\langle R_i,w_j\rangle\mid j,s]
=sc_r\langle h_j,w_j\rangle
=sc_r\sqrt r\,a_j.
\tag{HC.18}
```

The centered linear form is uniformly subgaussian because its independent
coordinate coefficients have squared sum one.  Since
`c_r^2r->infinity`, a union bound over all `r` rows gives

```math
|\langle R_i,w_{J_i}\rangle|
\ge(1-o(1))c_r\sqrt r\,a_{J_i}
\tag{HC.19}
```

simultaneously for every row whose label lies in `G_r`, with probability
tending to one.  For example, use relative error
`zeta_r=(c_r^2r)^(-1/4)`; the failure probability is at most
`2r exp[-c zeta_r^2c_r^2r]=o(1)`.

By Lemma HC.3, with probability tending to one there is `j in G_r` with

```math
N_j\ge(1-o(1)){\log r\over\log\log r}.
```

For that direction, `w_j in ran(V_r)`, so `B_rV_rw_j=B_rw_j`, and (HC.19)
gives

```math
\begin{aligned}
\|B_rV_r\|_{op}^2
&\ge\|B_rw_j\|_2^2\\
&\ge N_j(1-o(1))c_r^2r a_j^2\\
&\ge(L^2-o(1))r.
\end{aligned}
\tag{HC.20}
```

This proves (HC.14). `square`

The projection may be chosen optimally from the row law, including every
mean/covariance response direction.  It may not depend on the realized
sample; a sample-dependent rank-one projection could of course delete the
top right singular vector tautologically.

## 6. What the counterexample proves

1. **Bounded Renyi-two is a genuine threshold in the current proof, not a
   cosmetic convenience.**  Replacing a fixed bound by the broad condition
   `log K_(2,r)=o(r)` makes the sharp projected edge false.
2. **Mean/covariance information is insufficient.**  Here the mean is zero
   and the second moment is exactly identity; the outlier is created by a
   sample-level collision among many latent directions.
3. **Every deterministic sublinear-rank projection misses the collision.**
   The Hadamard directions distribute their squared mass evenly across any
   fixed `o(r)`-rank projection, while occupancy finds a large cluster among
   the remaining `r-o(r)` directions.
4. **The missing resource is a Renyi--occupancy tradeoff.**  A row law can
   spend only `exp(o(r))` likelihood to plant a weak bias in each of `r`
   possible directions; iid sampling amplifies the most occupied direction
   to leading operator scale.

The construction does not show a favorable conference-pressure basin.
Its super-edge places it outside the regular-bulk comparison, but the
resulting finite-rank sample outlier could raise rather than lower pressure.
What it decisively falsifies is a subexponential-density extension of the
sharp-edge/no-gain proof through deterministic mean/covariance peeling.
