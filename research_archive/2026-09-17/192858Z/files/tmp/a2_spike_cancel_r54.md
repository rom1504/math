# Wave 54: `A^2` completion cancellation and a hidden-spike obstruction

## Status and outcome

No exact-minimizer box theorem is proved.  The attack does produce three
exact results.

1. A large internal Gram term is exactly positive curvature in the *frozen*
   principal block of `A^2-(n-1)I`; its only direct outside trace is the
   leakage vector `d`.  There is no trace identity transferring that
   curvature to a negative direction of `H_T`.
2. There are explicit product-bias certificates that quantify exactly when
   `d` or a delocalized negative eigenvector of `H_T` cancels the internal
   term.  A spike-column identity shows the precise cross-correlation terms
   that can erase the apparent contribution `k r_i`.
3. Those erasures are real at the full obstruction scale.  A scalable
   fixed-density family of complete signings has competitive
   `Q(A)=O(n^(3/2))`, an exact child ground with `X=0` and
   `I=Theta(n^(5/2))`, but identically `d=0` and `H_T=0`.  Every outside
   completion then has row square `Theta(n^(5/2))`.

The scalable family is **not asserted to consist of exact minimizers**.  It
therefore falsifies a generic block-algebra or operator-norm implication, not
a theorem using exact discrete signing minimality.  Stored exact A9 and A10
instances give finite pointwise versions of the same hiding phenomenon, but
they are not asymptotic counterexamples.

All finite identities are reproduced by
`tmp/a2_spike_cancel_r54_check.py`; its saved output is
`tmp/a2_spike_cancel_r54_check.out`.

## 1. The spike first appears as frozen curvature, not outside curvature

Write the signing in blocks

```math
A=\begin{pmatrix}P&B^{\mathsf T}\\ B&C\end{pmatrix},
\qquad P=A[S],\quad B=A[T,S],\quad C=A[T],
```

where `|S|=m`, `|T|=k`, and let `y` be an oriented child ground.  Put

```math
u=By,\qquad I=\lVert Py\rVert_2^2,\qquad X=\lVert u\rVert_2^2,
\qquad \mathsf H=A^2-(n-1)I_n.
```

The blocks relevant to completion are exactly

```math
d=\mathsf H[T,S]y=BPy+Cu,
\qquad
\mathsf H_T=BB^{\mathsf T}+C^2-(n-1)I_k.
\tag{W54.1}
```

If `mathsf H_S=mathsf H[S]`, then direct multiplication gives the
**Verified frozen-curvature identity**

```math
\boxed{
I+X=m(n-1)+y^{\mathsf T}\mathsf H_Sy,
\qquad d=\mathsf H[T,S]y.}
\tag{W54.2}
```

Equivalently, for `v=(y,0)`,

```math
v^{\mathsf T}\mathsf Hv=I+X-m(n-1),
\qquad
\mathsf Hv=(\mathsf H_Sy,d).
```

Thus

```math
\boxed{
\lVert\mathsf H_Sy\rVert_2^2+\lVert d\rVert_2^2
\ge \frac{\{I+X-m(n-1)\}^2}{m}.}
\tag{W54.3}
```

This is the strongest immediate curvature conclusion.  It does **not**
lower-bound `d`: all of (W54.3) can remain in the frozen block.  Moreover,
both principal blocks of `mathsf H` have zero trace separately because
`diag(A^2)=(n-1)1`.  A positive Rayleigh quotient of `mathsf H_S` therefore
creates negative curvature inside `S`, but trace balance does not move that
negative curvature to `mathsf H_T`, the block accessible to completion.

## 2. Exact linear and negative-curvature completion certificates

Let

```math
\mu=I+X+k(n-1).
```

For any `z in [-1,1]^T`, independently choose outside signs with
`E w=z`.  Since `diag(mathsf H_T)=0`, the exact centered polynomial gives

```math
\boxed{
\mathcal V(S,y)
\le \mu+2d^{\mathsf T}z+z^{\mathsf T}\mathsf H_Tz.}
\tag{W54.4}
```

Conditional expectation derandomizes this bound.  Two useful specializations
are stronger than the uniform pairing bound when the relevant direction is
structured.

### 2.1 Signed-`d` bias

Choose any sign vector `s` with `s_j=sign(d_j)` on the nonzero coordinates
and put `h_d=s^T mathsf H_T s`.  Taking `z=-theta s` proves

```math
\boxed{
\mathcal V(S,y)
\le \mu-
\max_{0\le\theta\le1}
\{2\theta\lVert d\rVert_1-\theta^2h_d\}.}
\tag{W54.5}
```

Unlike (10.1288), this can extract an `l1` rather than `l2` linear gain,
provided the same sign direction does not pay too much positive quadratic
curvature.

### 2.2 Delocalized negative eigenvector

Let `mathsf H_T e=lambda e`, `||e||_2=1`, with `lambda<0`.  Use one of the
two mean vectors `z=+e/||e||_infty` and `z=-e/||e||_infty`, choosing the
sign against `d`.  Then

```math
\boxed{
\mathcal V(S,y)
\le \mu-
\frac{|\lambda|}{\lVert e\rVert_\infty^2}
-\frac{2|d^{\mathsf T}e|}{\lVert e\rVert_\infty}.}
\tag{W54.6}
```

In particular, a negative mode cancels an internal obstruction `I` if it is
both sufficiently negative and sufficiently delocalized that
`|lambda|/||e||_infty^2` is of order `I`.  Merely knowing that the trace-zero
matrix `mathsf H_T` has a negative eigenvalue supplies neither requirement.

Equations (W54.5)--(W54.6) give a concrete sufficient continuation.  On a
low-cross child-ground incidence, it would suffice to prove from exact
minimality that one of their displayed gains is at least
`I-O(n^(9/4-c))`.  No such implication is currently known.

## 3. Exact spike-column identity and its cancellation defect

Orient the internal fields as in (10.1289):

```math
r_i=\sigma y_i(Py)_i\ge0,
\qquad Py=\sigma\operatorname{diag}(y)r.
```

Let `b_i=Be_i` be cross column `i` and define the sign vector

```math
c_i=\sigma y_i b_i\in\{\pm1\}^T,
\qquad a_i=c_i^{\mathsf T}d,
\qquad h_i=c_i^{\mathsf T}\mathsf H_Tc_i.
```

Directly expanding (W54.1) gives the **Verified spike-column identities**

```math
\boxed{
\begin{aligned}
a_i
&=kr_i+
\sum_{j\ne i}y_iy_jr_j\langle b_i,b_j\rangle
+\sigma y_i\langle b_i,Cu\rangle,\\
h_i
&=k^2+
\sum_{j\ne i}\langle b_i,b_j\rangle^2
+\lVert Cb_i\rVert_2^2-k(n-1).
\end{aligned}}
\tag{W54.7}
```

Applying (W54.4) with `z=-theta c_i` yields

```math
\boxed{
\mathcal V(S,y)
\le\mu+\min_{0\le\theta\le1}
\{-2\theta a_i+\theta^2h_i\}.}
\tag{W54.8}
```

When `a_i>=0`, the certified gain is

```math
G_i=
\begin{cases}
2a_i-h_i,&h_i\le a_i,\\
a_i^2/h_i,&h_i>a_i.
\end{cases}
\tag{W54.9}
```

The first term in (W54.7) is the hoped-for contribution of a spike.  The
other cross columns can cancel it with either sign, however, and low `X`
only bounds the final term by

```math
|\langle b_i,Cu\rangle|
\le\lVert Cb_i\rVert_2\sqrt X.
```

There is also a useful aggregate form.  Put `g=BPy`.  Since `d=g+Cu`,

```math
\boxed{
\sum_{i\in S}r_i a_i
=g^{\mathsf T}d
=\frac12\{
\lVert d\rVert_2^2+
\lVert g\rVert_2^2-
\lVert Cu\rVert_2^2\}.}
\tag{W54.10}
```

Thus a large observable cross gradient `g`, relative to
`||C||op sqrt(X)`, forces a positive weighted spike-column coefficient.
But (W54.10) does not lower-bound `g` from `r`: the cross block can annihilate
the entire paired field vector.  The next construction does exactly that.

## 4. Scalable competitive signing with a completely hidden spike

This section is a **Verified scalable generic obstruction**, not an
exact-minimizer construction.

Take an unbounded sequence of primes `q=1 mod 4` and set `k=q+1`; Dirichlet's
theorem supplies infinitely many such primes.  Let `C` be the symmetric
Paley conference signing of order `k`, so

```math
C^2=(k-1)I_k.
```

Let `s` be the least power of two at least `k`.  Then `k<=s<2k`; a Sylvester
Hadamard matrix `W` of order `s` exists.  Put `h=floor(sqrt(s))`.  Choose a
competitive signing `G` of order `s-h`, with
`Q(G)=O(s^(3/2))` and `||G||op=O(sqrt(s))`, and switch and orient it so that
`1` is a positive absolute ground.  Such `G` exist simultaneously by the
standard random-signing union and spectral-norm bounds.

Adjoin `h` positive hubs to `G`, all joined positively to each other and to
the core, producing an order-`s` signing `D`.  Triangle inequality is exact
at the all-one ground:

```math
Q(D)=Q(G)+2h(s-h)+h(h-1)=O(s^{3/2}).
\tag{W54.11}
```

Each hub field is `s-1`.

Now pair-duplicate `D`.  The child signing `P` has vertices `i^+,i^-`, twin
edge `+1`, and all four edges between distinct pairs `i,j` equal to `D_ij`.
For a spin word `x`, set
`z_i=(x_(i^+)+x_(i^-))/2 in {-1,0,1}`.  Then

```math
x^{\mathsf T}Px
=4z^{\mathsf T}Dz+2\sum_i x_{i^+}x_{i^-}.
```

Randomly completing the zero coordinates of `z` proves
`|z^TDz|<=Q(D)`, and equality is attained at the all-one word.  Therefore

```math
\boxed{Q(P)=4Q(D)+2s=O(s^{3/2}).}
\tag{W54.12}
```

The two vertices over a base vertex have equal fields.  In particular the
`2h` duplicated hubs have field `2s-1`, so for `y=1_(2s)`,

```math
I=\lVert Py\rVert_2^2
\ge2h(2s-1)^2=\Theta(s^{5/2}).
\tag{W54.13}
```

Take any `k` rows `W_0` of `W` and define the cross block

```math
B=(W_0\; -W_0).
```

Then, with `m=2s`,

```math
BB^{\mathsf T}=mI_k,
\qquad By=0,
\qquad BPy=0,
\tag{W54.14}
```

where the last identity uses equality of each pair of fields.  Form the
complete signing

```math
A=\begin{pmatrix}P&B^{\mathsf T}\\B&C\end{pmatrix},
\qquad n=m+k.
```

Its density obeys

```math
\frac23\le\frac mn<\frac45.
```

For every full Boolean word, block Cauchy--Schwarz and
`||B||op=sqrt(m)` give

```math
\boxed{
Q(A)\le Q(P)+Q(C)+2m\sqrt k=O(n^{3/2}).}
\tag{W54.15}
```

The operator estimate is direct as well.  The positive hub join has norm
`O(s^(3/4))`, so `||D||op=O(s^(3/4))`.  In pair ordering,
`P=D tensor J_2+I_s tensor K_2`, whence
`||P||op<=2||D||op+1`.  Together with
`||B||op=sqrt(m)` and `||C||op=sqrt(k-1)`, the block triangle inequality
gives `||A||op=O(n^(3/4))`.  Nevertheless, (W54.14) and the conference
identity give, **exactly**,

```math
X=0,
\qquad d=BPy+CBy=0,
\qquad
\mathsf H_T=BB^{\mathsf T}+C^2-(n-1)I_k=0.
\tag{W54.16}
```

Consequently every one of the `2^k` completions has the same row square:

```math
\boxed{
\mathcal V(S,y)=I+k(n-1)=\Theta(n^{5/2}).}
\tag{W54.17}
```

For every `c>0`, this is asymptotically larger than the desired
`n^(9/4-c)` box scale.  The construction does not claim that this selector
is the best selector of `A`, and it does not claim that `A` minimizes `Q`
among order-`n` signings.  Its exact conclusion is narrower and decisive:

> Child-ground optimality, low cross energy, competitive cut cap, the known
> operator bound, and the raw block identities cannot force a spike into
> either `d` or negative `H_T` curvature.  A positive theorem must use
> discrete exact-minimizer structure or a different selector-existence
> mechanism.

The checker realizes this template at `(n,m,k)=(22,16,6)` using A6 as both
the conference exterior and the core below two hubs.  It obtains

```math
r_{\max}=15,
\quad I=1840,
\quad X=\lVert d\rVert_2^2=\lVert\mathsf H_T\rVert_F^2=0,
\quad \mathcal V=1966.
```

## 5. Stored exact-minimizer finite evidence

Exhaustive integer enumeration gives three sharply scoped checks.

- On A9 at `m=8`, there is a child ground with fields
  `(7,3,7,1,1,1,3,1)`, `I=120`, `X=0`, `d=0`, and `H_T=0`.  Its only two
  completions (an antipodal pair) both have row square `128`.
- On stored A10 at `m=8`, there is a child ground with fields
  `(1,1,1,7,1,3,7,3)`, `I=120`, `X=0`, `d=0`, and `H_T=0`; every completion
  has row square `138`.
- The best stored A10 box at `m=6` exhibits the opposite mechanism:
  `mu=74`, `H_T=0`, and `d=(-8,8,8,-8)`.  Choosing
  `w=-sign(d)` makes (W54.5) exact and gives `mathcal V=74-64=10`.

The first two are **Verified exact finite pointwise obstructions** to a
universal claim that exact minimality makes every field spike visible in
`d` or `H_T`.  Their codimensions are one and two, and their field sizes are
finite, so they do not falsify a fixed-density asymptotic existence theorem.
The third verifies that the signed-`d` route can supply all of the decisive
finite cancellation; in that example no negative quadratic curvature is
present.

## 6. Exact remaining step and research judgment

At present a low-cross field spike cannot be proved to force `d` or a
negative direction of `H_T` at any useful scale from the generic hypotheses.
The obstruction above makes the failure exact: both quantities can vanish
while the spike produces `Theta(n^(5/2))` internal mass.

The cleanest positive continuation is therefore minimizer-specific:

> **Open localized-curvature exclusion.**  Show that every required exact
> order-`n` minimizer has some fixed-density selector and child ground in the
> low-cross branch for which either the signed-`d` gain in (W54.5), or the
> delocalized negative-curvature gain in (W54.6), is at least
> `I-O(n^(9/4-c))`.

Together with (10.1287), this proves the child-ground box witness.  A weaker
route may first show that exact edge minimality rules out the hidden-curvature
configuration `BPy approx 0`, `H_T approx 0` when `I` is above project scale.
The paired construction identifies what such a theorem must detect; trace,
one-spin optimality, cut cap, and operator norm do not detect it.

An asymptotic family of exact minimizers reproducing (W54.16)--(W54.17) for
every usable low-cross child-ground incidence would falsify this
implementation.  No such exact-minimizer family is known.  Accordingly the
box route remains open, but direct `A^2` cancellation should not be credited
without a genuinely new discrete-minimality input.
