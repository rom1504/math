# Sparse stability in the actual balanced bulk

2026-09-07. **Exact flip identities and a positive phase-exclusion theorem;
director and independent construction-agent audits PASS.** A sharp undiluted Walsh example preserves
the limitation of zero-row-sum reasoning. No endpoint interpolation or
whole mixed-profile cap theorem is asserted.

## 1. Actual balanced weave and the ratio completion

There are `m` physical fibres of even size `q`, with `q>=p0 m` and
`q=O(m)`, and `N=mq`. At fibre `i` let `T_i` be a `q` by `m` real
transform whose active columns are balanced, and suppose

```
T_i^T 1=0,       ||T_i||op<=L sqrt(m).                   (1)
```

The application uses actual sign columns produced by balanced-column
repair; zero self/masked columns are harmless. Conditional on these
transforms and a fixed symmetric macro-edge mask, use independent fair
signs `S_ij=S_ji` on the unmasked off-diagonal macro edges. Put

```
W_ij=S_ij T_i(:,j) T_j(:,i)^T,       W_ii=0.             (2)
```

This is the actual masked bulk, not a Gaussian replacement. It is hollow
and annihilates EVERY fibre-constant vector. The transforms and mask must
be chosen independently of the remaining fair edge signs, or have that
conditional law; an arbitrary seed-dependent edge conditioning is not
silently permitted.

For a spin word `x`, write `s_i=sum_a x_i(a)`, `a_i=s_i/q`, and
`P_i=I-J/q`. A variance-weighted bulk inequality at coefficient `b>0`
is exactly equivalent to BOTH inequalities

```
F_(sigma,b)(x):=sigma H_W(x)+(b sqrt(N)/q) sum_i s_i^2
               <=b N^(3/2),       sigma=+-1.             (3)
```

Indeed `sum_i||P_i x_i||^2=N-sum_i s_i^2/q`. This is an actual
ferromagnetic completion: its added off-diagonal within-fibre coupling
is `2b sqrt(N)/q`. The difference between the expression in (3) and
the hollow quadratic form of that completion is the fixed constant
`b m sqrt(N)`.

## 2. Exact one-spin correction, including the `1/q` term

Set `f=Wx`. Flipping coordinate `x_i(a)` changes (3) by exactly

```
-2 sigma x_i(a) f_i(a)
       -(4b sqrt(N)/q) x_i(a) s_i+4b sqrt(N)/q.
```

Thus every local maximum of `F_(sigma,b)` satisfies

```
x_i(a)[sigma f_i(a)/sqrt(N)+2b a_i]>=2b/q.               (4)
```

Dropping the right side without comment would lose the physical
single-spin correction. Equivalently, (4) is ordinary stability for
the hollow ferromagnetic completion.

If a nonconstant fibre has majority sign `epsilon_i`, minority set
`R_i`, fraction `r_i=|R_i|/q<=1/2`, then `epsilon_i a_i=1-2r_i`.
At EVERY minority coordinate, (4) requires

```
sigma epsilon_i f_i(a)/sqrt(N)
       <=-2b(1-2r_i)-2b/q.                              (5)
```

In particular the negative threshold in (5) has magnitude at least `b`
when `r_i<=1/4`. Unlike ordinary mean-zero stability, it does not tend
to zero with the minority fraction.

For comparison, if one maximizes only `sigma H_W`, a stable fibre field
has the exact cone constraint

```
||f_i||_1<=2sqrt(|R_i|)||f_i||_2.                        (6)
```

The signs of `sigma epsilon_i f_i` are nonnegative off `R_i` and
nonpositive on `R_i`, while its total sum is zero by (1). Hence its
positive and negative masses agree, and Cauchy--Schwarz on `R_i`
proves (6). A fully constant stable fibre must have `f_i=0` exactly.
These ordinary-stability statements do NOT replace (4) for the ratio.

## 3. A uniform nearconstant stability penalty under an incoming-tail bound

Fix a spin word and all the transforms. Define the actual incoming
coefficients at fibre `i` by

```
v_j=q^(-1/2) T_j(:,i)^T x_j,
```

and set masked/self coefficients to zero. Then
`f_i/sqrt(N)=(T_i/sqrt(m))(S_ij v_j)_j` exactly. Fix `A>0` and

```
0<epsilon<=b^2 p0/(16L^2).
```

Suppose `0<r_i<=r0<=1/4` and

```
sum_(j:v_j^2>A r_i) v_j^2<=epsilon r_i m.               (7)
```

Then the probability, over the actual fair edge signs incident to this
fibre, that its minority coordinates satisfy (5) is at most

```
exp(-c q),       c=b^2/(64L^2 A).                       (8)
```

This bound is UNIFORM as `r_i` tends to zero through positive allowed
fractions. It needs no positive variance lower bound, no normal
approximation, and no high-dimensional Gaussian orthant comparison.

To prove it, condition on all heavy signs with `v_j^2>A r_i`. Their
total conditional mean on the minority set, in normalized field units,
has absolute value at most

```
(1/sqrt(m)) ||T_i^T 1_(R_i)|| ||v_heavy||
 <=L sqrt(|R_i|) sqrt(epsilon r_i m)
 <=L sqrt(epsilon/p0)|R_i| <=(b/4)|R_i|.                 (9)
```

Let `alpha_i=2b(1-2r_i)+2b/q>=b`. The convex violation function

```
V_i(S_light)=sum_(a in R_i)
       (sigma epsilon_i f_i(a)/sqrt(N)+alpha_i)_+
```

therefore has conditional expectation at least `b|R_i|/2`, by
Jensen and (9). Its Euclidean Lipschitz constant is at most

```
L sqrt(A r_i |R_i|).
```

The elementary convex-violation version of Talagrand's product theorem,
proved in `principle_director_stability_orthant_penalty_2026_09_07.md`,
says `Pr(V_i=0)<=exp(-(E V_i)^2/(16 Lip(V_i)^2))`.
It gives (8), uniformly in the heavy signs, so their conditioning may
be removed. If the Lipschitz constant vanishes, the positive expectation
already makes the zero event empty.

## 4. A genuine global phase-exclusion statement

Call a nonconstant fibre GOOD for a given word if it has `r_i<=r0`
and satisfies (7); otherwise call it BAD. A constant fibre is neither.
Choose `r0>0` sufficiently small and `delta in (0,1)` sufficiently small
that

```
h(r0)+delta log(2)<c(1-delta)/4,                        (10)
```

where `h` is binary entropy. These are fixed constants before the order
grows. Then, with probability `1-exp(-c' q)` over the actual fair macro
signs, EVERY nonconstant GLOBAL maximum of either `F_(+,b)` or
`F_(-,b)` has at least a fraction `delta` of its nonconstant fibres BAD.
Here `c'>0` depends only on the displayed fixed parameters.

For the proof, let `ell` be the number of nonconstant fibres. In a
global maximum, all constant fibres may be reset to the positive
constant without changing the value of `F`: `W` annihilates each such
fibre constant and its squared magnetization is unchanged. The resulting
word is STILL a global maximum and hence satisfies (4). This is why
there is no artificial `2^m` charge for constant-fibre polarities.

There are at most `binom(m,ell)` choices of active fibres. A good one
has at most `2 sum_(h=1)^(floor(r0q)) binom(q,h)` possible words, whose
logarithm is `q h(r0)+O(log q)`. A bad one has at most `2^q` words.
Specifying their designations costs at most `2^ell`.

For any fixed candidate with `g` good fibres, each independent edge sign
occurs in at most two minority-stability events. Graph Holder/Finner and
(8) give probability at most `exp(-c q g/2)`. The classification (7)
depends on the fixed word and transforms, not the outer signs, so it
can be imposed before this estimate.

If fewer than `delta ell` fibres are bad, the log of the resulting
union bound, divided by `ell q`, is at most

```
h(r0)+delta log2-c(1-delta)/2+O(log(mq)/q),
```

strictly negative by (10). Sum the resulting geometric bound over
`ell>=1` and both `sigma`. This proves the stated phase exclusion.

It is a theorem about the actual fair-edge construction conditional on
EVERY transform array obeying (1). It is not a cap bound for every biased
slice: the maximum inside a constrained slice need not satisfy (4).
Nor does it control a global maximum whose active fibres have a positive
fraction of intermediate biases or coherent incoming tails. Those are
the precisely remaining branches.

## 4a. Stronger minority-mass dichotomy, with an `exp(-cN)` failure bound

The light cutoff need not use the same scale as the local minority
fraction. For ANY `rho>0`, replace the threshold `A r_i` in (7) by
`A rho`, but KEEP its right side `epsilon r_i m`. The heavy-mean
bound (9) is unchanged. Only the Lipschitz constant changes, to
`L sqrt(A rho |R_i|)`. The same proof therefore gives

```
Pr(fibre i satisfies minority stability)
       <=exp[-c q r_i/rho],                             (14)
```

with the SAME `c=b^2/(64L^2 A)`. There is no lower bound on `rho`
or on a positive `r_i`. If there are no light variables, the zero event
is empty by the heavy-mean bound, as before.

For a nonconstant word put

```
rho=rbar=(1/m)sum_i r_i>0.
```

Call fibre `i` MASS-GOOD if

```
0<r_i<=1/4,
sum_(j:v_j^2>A rbar) v_j^2<=epsilon r_i m.               (15)
```

Fix any `delta in (0,1)`. There is a sufficiently small fixed `r0>0`
such that, with probability `1-exp(-c' N)`, neither ferromagnetic
completion in (3) has a nonconstant LOCAL maximum satisfying BOTH

```
0<rbar<=r0,
sum_(i mass-good) r_i >=(1-delta) sum_i r_i.             (16)
```

Thus every very near-fibre-constant local maximum must put at least
a fixed fraction of its TOTAL minority mass in intermediate-bias
fibres or incoming-heavy-tail fibres. This strengthens the active-fibre
version: it does not merely require one exceptional small fibre.

For each fixed word, apply Finner to (14) over its mass-good fibres.
Their total stability probability is at most

```
exp[-(c q/(2rbar)) sum_(i mass-good) r_i]
       <=exp[-c(1-delta)N/2].                          (17)
```

The classification (15) is independent of the fair outer signs. All
words with `rbar<=r0` can be described by one majority sign per fibre
and a minority subset of the `N` sites of size at most `r0 N`. Therefore
their count is at most

```
2^m sum_(s=0)^(floor(r0N)) binom(N,s)
       <=exp[N h(r0)+m log2+O(log N)].                 (18)
```

Choose, for example, `h(r0)<c(1-delta)/8`. Since `m/N=1/q ->0`,
(17),(18), and a union over the two polarities give the asserted
`exp(-c'N)` bound. This time there is no need to canonically reset
constant-fibre polarities or to require a GLOBAL maximum.

The theorem is still not a bound on a constrained biased-slice maximum.
It is a positive uniform local-stability theorem for the actual
ferromagnetic completion, with the remaining heavy/intermediate
minority mass quantified. In particular no assertion that (15) holds
for every response of the repaired recursive ensemble is being inserted.

## 4b. Entropy-adaptive sharpening: almost all sparse minority mass is exceptional

There is a useful finite-order strengthening of Section 4a. For every
nonconstant word let

```
f_good=(sum_(i mass-good) r_i)/(sum_i r_i),
eta_N(r)=(4/c)[h(r)+(log2)/q+2(log N)/N].               (19)
```

With probability at least `1-2^(1-m) N^(-3)`, simultaneously every
nonconstant local maximum of either completion satisfies

```
f_good<eta_N(rbar) whenever eta_N(rbar)<=1.             (20)
```

To prove this, fix the total minority count `s=1,...,floor(N/2)`.
There are at most `2^m binom(N,s)` possible words. Their classifying
majority signs and minority subsets need not be unique; overcounting
is harmless. For each word with `f_good>=eta_N(s/N)`, (14) and Finner
bound stability by

```
exp[-cN eta_N(s/N)/2]
 =exp[-2N h(s/N)-2m log2-4log N].
```

Multiplying by the word count, using
`binom(N,s)<=exp[N h(s/N)]`, summing over at most `N` possible counts
and both polarities, gives failure probability at most
`2^(1-m)N^(-3)`. If the threshold exceeds one the statement is vacuous.

Consequently, along ANY sequence of local maxima with `rbar->0`, the
fraction of minority mass on mass-good fibres tends to zero in this
same simultaneous high-probability event. Almost all minority mass
must be on fibres with `r_i>1/4` or

```
sum_(j:v_j^2>A rbar) v_j^2>epsilon r_i m.
```

This strengthens the fixed-fraction alternative without assuming a
bounded ratio between the different positive minority densities. It
does not exclude concentration of almost all mass on a few
intermediate-bias fibres, nor does it show that coherent-tail rows fall
in one bounded-ratio band.

## 4c. Mean-weighted version: noncoherent mass must be nearly balanced

The cutoff `r_i<=1/4` is convenient but not essential. Put
`d_i=1-2r_i=|a_i|`, and call a fibre MEAN-GOOD if `r_i>0`, `d_i>0`,
and

```
sum_(j:v_j^2>A rbar) v_j^2<=epsilon d_i^2 r_i m.        (21)
```

The same proof applies: the minority threshold in (5) is at least
`2b d_i`; the absolute heavy aggregate mean is at most
`(b/4)d_i|R_i|`; and the light hinge Lipschitz constant is still at most
`L sqrt(A rbar |R_i|)`. In particular its expectation is at least
`(b/2)d_i|R_i|`. Thus, with the original conservative constant `c`,

```
Pr(minority stability at i)<=exp[-c q d_i^2 r_i/rbar].  (22)
```

No lower bound on a positive `d_i` is needed. A balanced fibre has
`d_i=0` and supplies no penalty by this estimate; its exact `2b/q`
threshold has not been mistaken for a leading field.

Define the normalized mean-weighted good mass

```
f_mean=(sum_(i mean-good) d_i^2 r_i)/(sum_i r_i).
```

Repeating Section 4b verbatim, with the SAME `eta_N`, proves that with
probability at least `1-2^(1-m)N^(-3)`, every nonconstant local maximum
of either completion has

```
f_mean<eta_N(rbar) whenever eta_N(rbar)<=1.             (23)
```

For any fixed `d0>0`, (23) bounds the fraction of total minority mass on
mean-good fibres with `|a_i|>=d0` by `eta_N(rbar)/d0^2`. Consequently,
as `rbar->0`, almost all minority mass is either on fibres with
`|a_i|<d0` (near the balanced face) or on fibres violating the
mean-weighted tail budget (21). There is no separate arbitrary
intermediate-bias exception in this version.

This is still not a whole mixed-face energy bound: a certificate on
balanced words cannot simply be applied to the balanced part of a
mixed word without paying its cross interaction. The statement
identifies the surviving local-maximizer phases; it does not add their
separately optimized energies.

## 5. A sharp zero-row-sum boundary: sparse stable Walsh eigenstates

Zero row sums alone cannot remove the coherent branch in Section 4.
Here is an exact example in the UNDILUTED balanced weave, `q=m=2^d`.
Index physical and macro coordinates by `F_2^d`, and take

```
T_i(a,j)=(-1)^(a dot(i+j)),       j!=i,
T_i(a,i)=0.
```

These are full Walsh frames with the marked constant self-port removed;
`T_i T_i^T=mI-J=mP`. Choose a physical subspace `B` of fixed codimension
`d0>=1`, put `r=2^(-d0)`, and let `K=B^perp`, of size `1/r`. In every
fibre choose the SAME nearconstant word

```
x_i(a)=1-2 1_B(a),        mean(x_i)=1-2r.
```

Its nonconstant port transform is exactly

```
T_i(:,j)^T x_i=-2rm  if i+j in K\{0},
T_i(:,j)^T x_i=0     otherwise.                         (11)
```

Set the outer macro signs to `+1` on differences in `K\{0}`. These
edges form disjoint cliques on the cosets of `K`, with only
`m(1/r-1)/2=O_r(m)` edges. On every other edge the macro signs are
arbitrary. In particular one may start with ANY bounded-cap macro
signing and force these edges, at cap cost only `O_r(m)`.

Character orthogonality and (11) give the EXACT identity

```
Wx=m P_all x.                                          (12)
```

Indeed the field at `a` is
`-2rm sum_(k in K\{0})(-1)^(a dot k)
 =m[-2 1_B(a)+2r]`. Thus `x` is a STRICT positive one-spin-stable
word at every arbitrarily small fixed minority fraction `r=2^(-d0)`.

The signed-swap factorization gives `||W||op<=m=sqrt(N)` and
`W=P_all W P_all`. Consequently

```
|H_W(y)|<=sqrt(N)||P_all y||^2/2   for every y,
H_W(x)=sqrt(N)||P_all x||^2/2.                          (13)
```

The variance ratio is therefore EXACTLY `1/2`, attained by this sparse
strictly stable word. For every `b<1/2`, (3) is violated. Under independent
fair macro signs the event used in (11) has probability
`2^(-m(1/r-1)/2)=exp(-O_r(m))`, not `exp(-c N)`.

This is an actual masked balanced weave and can have a bounded-cap
macro seed. It demonstrates why balance alone cannot give a universal
leading stability penalty or a strict variance-ratio theorem. It does
NOT falsify the independently selected, strictly diluted `p<1`
recursive repaired ensemble. In the example every active incoming
profile places its energy on only `1/r-1` coordinates, and explicitly
violates (7). Thus it is consistent with, and explains the need for,
the positive tail-sensitive phase theorem above.

## 6. Exact finite verification and scope ledger

`computations/principle_synthesis_2026_09_07_balanced_bulk_stability_check.py`
checks eight Walsh examples of orders `m=4,8,16`, at minority fractions
`1/2,1/4,1/8` where available. It verifies symmetry and the hollow
diagonal, the kernel of EVERY fibre-constant vector, the exact identity
`Wx=m P_all x`, strict ordinary one-spin stability, and the exact
variance ratio `1/2`. All examples violate (3) at the test coefficient
`b=7/16`. It also checks 800 rational-valued instances of the complete
single-spin difference in Section 2. All checks PASS; their exact
outputs are preserved in the corresponding `_check_results.json`.

The initial checker invocation stopped before any mathematical test
because this workspace uses Python 3.9, without `int.bit_count()`.
Replacing that call by binary-string population counting was the sole
compatibility correction before the successful run. These finite tests
check the identities and boundary example, not the asymptotic
Talagrand or Finner theorem.

Sections 3--4a hold conditional on every transform array satisfying
(1), provided the remaining outer signs are independent and fair.
They do not require an optimizer, a chosen macro seed, or a Gaussian
surrogate. Their output is an actual local-maximizer exclusion with a
quantified remaining phase, not a uniform bound on all spin words or
on constrained-slice maximizers. The coherent and intermediate-bias
phase must still be analyzed under the favorable recursive/selector
law before a whole mixed-profile inequality can follow.

The complementary positive packet theorem is
`principle_invent_2026_09_07_stratified_marked_selector.md`, especially
its clipping argument and bounded-ratio band conclusion. It controls
all nonconstant rows with minority densities in `[r,Kr]` under an
actual one-node stratified selector law, with a strict variance
coefficient for the one-hole case. Its band assumption and our
incoming-tail hypothesis are different: the present theorem permits
arbitrarily many density scales but leaves coherent/intermediate mass;
the packet theorem treats the coherent packet geometry within a
bounded-ratio band. Neither statement by itself, nor their formal
juxtaposition, handles every multiscale mixed profile. No such closure
is asserted here.
