# Wave 38 S2: favorable full-column Gram estimate

## Verdict

There is an exact signing-specific bridge from **near-parent complement-flip
certificates** to the external half of (10.1057), but it does not control the
internal half.  The bridge is

```math
2t(t-1)+4\lVert A[S^c,S]y\rVert_2^2
\le w(y)^2-a(y)^2
\le q_n^2-e_A(y)^2, \qquad t=n-m. \tag{R38.1}
```

For a complement-flip certificate with parent deficit `Delta`, after orienting
the signing so its child energy is `e`, (10.1044) gives

```math
e\ge q_n-\frac{\Delta}{2}.
```

Consequently

```math
\boxed{
4\lVert A[S^c,S]y\rVert_2^2+2t(t-1)
\le q_n\Delta-\frac{\Delta^2}{4}
\le q_n\Delta .}
\tag{R38.2}
```

Thus an affordable anchored family of such certificates with
`E Delta=O(n^(3/4-c'))` automatically has external Gram average
`O(n^(9/4-c'))`, and its restricted labels are genuinely favorable.  This is
a useful positive reduction.

It is also sharply limited.  The same certificate necessarily satisfies

```math
\Delta\ge 2\{q_n-Q(A[S])\}. \tag{R38.3}
```

Hence a project-scale-deficit family can exist only on selectors whose child
cap is within `O(n^(3/4-c'))` of the full parent cap.  If the natural
order-`n^(3/2)` cap drop persists on all affordable favorable families, this
bridge cannot reach (10.1057).  Moreover, completion width contains no control
of `||A[S]y||^2`: exact `A_9` child grounds with identical `(w,a,b)` have
different internal Gram energies, including a favorable exact ground with
zero external Gram and internal Gram `86`.

The surviving sufficient theorem is therefore a joint selection statement:
affordable anchored thinning, project-scale conflict, near-parent
complement-flip deficit, **and a separate internal-Gram average**.  Parseval
removes the external clause but does not supply that last clause.

## 1. Exact Parseval identity for the completion polynomial

Fix `S`, put `T=S^c` and `t=|T|`, and fix a projective child word `y`.  If a
response orientation `sigma` is present, replace `A` by `sigma A`; this
preserves every Gram norm, `q_n`, and completion width.  Write

```math
e=y^TA[S]y,
\qquad
f(x)=x^TA[T]x+2y^TA[S,T]x,
\qquad x\in\{\pm1\}^T.
```

Let

```math
Z_+=\max_x f(x),\quad Z_-=-\min_x f(x),\quad
w=\frac{Z_++Z_-}{2},\quad a=\frac{Z_+-Z_-}{2},\quad b=q_n-w.
```

Under uniform outside spins, the quadratic and linear Walsh levels are
orthogonal and `E f=0`.  The quadratic coefficients are `2a_ij`, while the
linear coefficient vector is `2A[T,S]y`.  Therefore

```math
\boxed{
\mathbb E_x f(x)^2
=2t(t-1)+4\lVert A[T,S]y\rVert_2^2 .}
\tag{R38.4}
```

Since `f` is mean zero and lies in `[-Z_-,Z_+]`, the elementary
Bhatia--Davis bound (equivalently, `(Z_+-f)(f+Z_-)>=0` averaged over `x`)
gives

```math
\mathbb E f^2\le Z_+Z_-=w^2-a^2. \tag{R38.5}
```

Every full completion has energy `e+f(x)` in `[-q_n,q_n]`.  Hence

```math
Z_+\le q_n-e,\qquad Z_-\le q_n+e,
```

and multiplication proves the second half of (R38.1).  This independently
recovers and strengthens the useful part of (10.1058): since `w=q_n-b`, the
endpoint inequalities imply

```math
|a+e|\le b,
```

so more precisely

```math
\boxed{
2t(t-1)+4\lVert A[T,S]y\rVert_2^2
\le (q_n-b)^2-(|e|-b)_+^2 .}
\tag{R38.6}
```

When `e>=b`, its right side is

```math
(q_n-e)(q_n+e-2b)\le 2q_n(q_n-e). \tag{R38.7}
```

All statements above are deterministic and apply to any signing, selector,
and word.

## 2. What width controls, and in which direction

Equation (R38.4)--(R38.5) gives the exact useful direction

```math
\boxed{
\lVert A[T,S]y\rVert_2^2
\le\frac{w^2-a^2-2t(t-1)}4 .}
\tag{R38.8}
```

Thus the relevant completion statistic is the **endpoint product**
`Z_+Z_-=w^2-a^2`, not width alone.  A large width can coexist with a small
product when the center is close to one endpoint.

This direction does not by itself help the paired-response target (10.1038).
That target asks for small width deficit `b`, hence `w` close to `q_n`.
Without simultaneous information forcing `|a|` close to `w`, (R38.8) is only
`O(q_n^2)=O(n^3)`, worse than the generic `O(n^(5/2))` child bound.  The
identity `|a+e|<=b` shows exactly what extra information is needed: `e` must
be close to `q_n`, not merely close to the child cap.

There is no useful reverse inequality from width to external Gram.  The
outside quadratic Fourier level already contributes `2t(t-1)` to (R38.4),
and an endpoint may occur on exponentially little outside mass.  Even the
full data `(w,a,b)` can have zero external level-one energy.  A reverse bound
would require a separate anti-concentration theorem for endpoint mass; it is
not a Parseval consequence.

## 3. Complement-flip certificates give the missing endpoint asymmetry

Let an oriented full state `d=(sigma,x)` have

```math
E_d=q_n-\Delta,
```

and suppose it certifies `S` through (10.1044).  Orient `A` by `sigma` and
put `y=x_S`.  Then

```math
e=y^T(\sigma A[S])y\ge q_n-\frac\Delta2.
```

Applying (R38.1), using `0<=Delta<=2q_n`, and squaring
`e>=q_n-Delta/2>=0` proves (R38.2).  Two further consequences should be
recorded with the Gram bound:

1. Since `e<=Q(A[S])<=q_n`,

   ```math
   Q(A[S])-e\le q_n-e\le\frac\Delta2. \tag{R38.9}
   ```

   Thus the restricted word is a genuine near-ground child label.  In the
   high-ratio window `p_2>=1/2`, (10.1045) already implies the stronger
   project-specific favorable inequality.

2. Since `e<=Q(A[S])`, the same inequality rearranges to (R38.3).  Near-parent
   deficit is therefore a strong cap-persistence demand, not a generic
   feature of complement certificates.

For a law `P` on an anchored affordable family `G`, choose one complement
certificate per selector.  Averaging (R38.2) yields

```math
\mathbb E_P\lVert A[S^c,S]y^S\rVert_2^2
\le\frac{q_n}{4}\mathbb E_P\Delta_S. \tag{R38.10}
```

As `q_n=O(n^(3/2))`, the hypothesis

```math
\mathbb E_P\Delta_S=O(n^{3/4-c'}) \tag{R38.11}
```

has exactly the required external exponent `O(n^(9/4-c'))`.  Markov thinning
can turn (R38.11) into a pointwise `Delta` cap at only constant or polynomial
loss once the average is known, which is negligible relative to the allowed
`exp{-O(n^(3/4-c'))}` density.  Thinning cannot create (R38.11) from the
generic `Delta=O(n^(3/2))` scale.

## 4. Exact sufficient selection package

The following is a clean sufficient successor to (10.1057).  For some fixed
`c'>0`, find an anchor `v`, a family

```math
G\subseteq\{S:|S|=m,\ v\in S\},
\qquad
U_v(G)\ge\exp\{-O(n^{3/4-c'})\},
```

and, for every `S in G`, an oriented complement-flip certificate `d_S` with
restricted label `y^S`, such that

```math
\begin{aligned}
\mathcal C_v(\mathbf y)&\le(d+1)/2,\\
\mathcal C_v(\mathbf y)+\log\beta^{-1}+\log n
&=O(n^{3/4-c'}),\\
\mathbb E_{U(G)}\Delta_{d_S}&=O(n^{3/4-c'}),\\
\mathbb E_{U(G)}\lVert A[S]y^S\rVert_2^2
&=O(n^{9/4-c'}).
\end{aligned}
\tag{R38.12}
```

Here `d` is the already chosen project decoder radius, and the second line is
the exact row-transfer scale because exact minimizers have
`||A||_op^2=O(n^(3/2))`.  The labels are favorable by (10.1045), (R38.10)
supplies the external Gram bound, and the last line supplies the internal
Gram bound.  Therefore

```math
\mathbb E_{U(G)}\lVert A[:,S]y^S\rVert_2^2
=O(n^{9/4-c'}),
```

which is precisely (10.1057); (10.1056) then gives the desired global row
target for the plurality decoder.

This package keeps the support, conflict, and thinning requirements explicit.
It does **not** use the unrestricted deletion minimizer from (10.1059), which
the `A_9` favorable-support wall invalidates.  A fractional-cover version
could put a Lagrange price on `Delta` and internal Gram simultaneously, but
the algebra above supplies only the `Delta -> external Gram` implication; it
does not prove an affordable cover with either priced moment.

## 5. The remaining internal-Gram wall

Decompose

```math
\lVert A[:,S]y\rVert_2^2
=\underbrace{\lVert A[S]y\rVert_2^2}_{L_S(y)}
+\underbrace{\lVert A[S^c,S]y\rVert_2^2}_{E_S(y)}.
```

The completion polynomial depends on `A[S]` only through the scalar internal
energy `e`; its Fourier coefficients contain no internal row vector
`A[S]y`.  Therefore no estimate based only on `(w,a,b,e)` can control
`L_S(y)`.  Exact child-ground stability still gives only

```math
L_S(y)\le(m-1)Q(A[S])=O(n^{5/2}),
```

missing (R38.12) by `n^(1/4+c')`.

A generic independent-noise smoothing does not repair this.  If each sign of
`y` is retained with correlation `theta`, then exactly

```math
\mathbb E\lVert A[:,S]Y\rVert_2^2
=\theta^2\lVert A[:,S]y\rVert_2^2
+(1-\theta^2)m(n-1),
\qquad
\mathbb E[Y^TA[S]Y]=\theta^2 e. \tag{R38.13}
```

Uniformly suppressing a possible `Theta(n^(5/2))` Gram value to
`O(n^(9/4-c'))` requires `theta^2=O(n^(-1/4-c'))`; this destroys a
`Theta(n^(3/2))` child energy down to `O(n^(5/4-c'))`.  Thus ordinary noise
cannot simultaneously guarantee Gram saving and a nontrivial favorable
energy threshold.  This is a generic-smoothing obstruction, not an
exact-minimizer falsifier.

## 6. Finite exact audit

The checker `tmp/favorable_gram_r38_check.py` verifies (R38.1), (R38.4),
(R38.6), and (R38.2) by integer enumeration on `A_6,A_8,A_9`.  It also records
the following actual exact-child-ground data after orienting each ground to
positive energy:

- On `A_8,m=5`, `(w,a,b)=(8,0,12)` has fixed external Gram `3`, while its
  internal Gram ranges over `24,32,40`.
- On `A_9,m=6`, `(w,a,b)=(8,-2,16)` has external Gram `4` or `8` and internal
  Gram `30,62,70,78`.
- Most sharply, `(w,a,b)=(4,-2,20)` on `A_9,m=6` has external Gram exactly
  `0`, while internal Gram is `46` or `86`.  One `86` witness is

  ```text
  S=(0,1,2,3,5,7), y=(1,-1,1,-1,-1,-1),
  q_9=24, Q(A_9[S])=e=22.
  ```

Every listed word is an exact positive child ground and hence a favorable
label.  These finite examples do not asymptotically falsify (10.1057); they
do prove that completion width and external Parseval data do not determine
the missing internal Gram term.

Run:

```bash
.venv/bin/python tmp/favorable_gram_r38_check.py
```

The current run ends with `PASS favorable_gram_r38_check`.

## 7. Falsification criteria and next useful attack

The near-parent-certificate implementation is falsified for a fixed saving
`c'>0` if, along an unbounded exact-minimizer sequence, every affordable
anchored favorable family supported on complement certificates has either

```math
\mathbb E\Delta=\omega(n^{3/4-c'})
```

or project-scale conflict failure.  By (R38.3), it already fails if every
such affordable family is supported on selectors with
`q_n-Q(A[S])=Omega(n^(3/2))`.  This would not falsify (10.1057), because low
Gram favorable labels need not arise from near-parent complement witnesses.

After (R38.2), the sharp remaining positive question is internal and
signing-specific:

> Can exact edge-sign minimality force an affordable coherent family of
> favorable child labels with
> `E ||A[S]y^S||^2=O(n^(9/4-c'))`, possibly after pricing this internal
> statistic directly in the complement-flip fractional cover?

Repeating completion-width, unrestricted deletion, or generic noise arguments
cannot answer this question.  A viable successor must use simultaneous
minimality of the parent signing (or prove a cap-persistence/row-concentration
falsifier), because the completion polynomial has now exhausted its direct
Gram content.
