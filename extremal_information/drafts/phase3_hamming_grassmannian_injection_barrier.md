# Low-weight sum codes and the injection-metric barrier

**Status.** The finite statements and asymptotic comparison below are proved.
They give a general-`k` Hamming--Hausdorff packing construction not requiring
one common separated host.  The same calculation then proves that this
construction can never improve the elementary common-host Gilbert exponent.
An explicit `(D,k,t)=(5,2,2)` example shows that this is a limitation of the
argument, not an exact optimality theorem for common hosts.

Throughout, `B_N(t)=sum_(i=0)^t binom(N,i)`, with `B_N(t)=2^N` when
`t>=N`.  For two `k`-subspaces `C,C'<=F_2^D`, write

```math
d_I(C,C')=k-\dim(C\cap C')
```

for their injection distance, and let `d_Hs` be Hausdorff distance induced by
ambient Hamming distance.

## 1. Hausdorff closeness forces many low words in the sum code

### Lemma HI.1 (low-weight sum-code certificate)

Let `C,C'` be `k`-subspaces, put `r=d_I(C,C')`, and put `L=C+C'`.  If

```math
d_Hs(C,C')\le t,
```

then

```math
\boxed{
2^r\le |L\cap B_D(0,t)|\le B_{k+r}(t).
}                                                    \tag{HI.1}
```

Consequently,

```math
B_{k+r}(t)<2^r
\quad\Longrightarrow\quad
d_Hs(C,C')>t.                                      \tag{HI.2}
```

#### Proof

The image of `C` in `L/C'` has dimension `r`, so it has `2^r` cosets.  For
one representative `c` of each coset, Hausdorff closeness supplies
`c' in C'` with `wt(c+c')<=t`.  The resulting words `c+c'` are distinct
modulo `C'`, hence distinct in `L`.  This proves the first inequality.

The code `L` has dimension `k+r`.  Choose an information set `I` of that
size.  Coordinate projection `L->F_2^I` is injective and cannot increase
weight.  It therefore injects `L cap B_D(0,t)` into the radius-`t` ball in
`F_2^(k+r)`, proving the second inequality. `square`

The upper count in (HI.1) is best possible from `dim L` alone: the coordinate
subspace `F_2^(k+r) times {0}` contains exactly `B_(k+r)(t)` words of weight
at most `t`.  Thus improving (HI.2) requires the rooted placement of `C` and
`C'` inside `L`, not a sharper universal weight-enumerator estimate for an
unrooted `(k+r)`-dimensional sum code.

This certificate uses the whole sum code, not a scalar distance between two
chosen bases.  It is nevertheless only a weight-enumerator certificate; the
comparison in Section 3 identifies its asymptotic limit.

## 2. A general injection-code lower bound

The ratio

```math
{B_{k+s}(t)\over 2^s}                              \tag{HI.3}
```

is nonincreasing in `s`.  Indeed,

```math
B_{n+1}(t)\le2B_n(t).
```

Thus one check at the minimum injection distance controls all more distant
pairs.

### Theorem HI.2 (injection-to-Hausdorff transfer)

Suppose

```math
1\le r\le\min\{k,D-k\},
\qquad
B_{k+r}(t)<2^r.                                    \tag{HI.4}
```

Then there is a family `A subseteq Gr_k(F_2^D)` with pairwise Hamming
Hausdorff distance greater than `t` and

```math
|A|
\ge
{ {D\brack k}_2
 \over
 \displaystyle\sum_{s=0}^{r-1}
 2^{s^2}{k\brack s}_2{D-k\brack s}_2 }           \tag{HI.5}
```

and hence

```math
\boxed{
|A|\ge
{2^{k(D-k)-(r-1)(D-r+1)}\over16r}.
}                                                  \tag{HI.6}
```

#### Proof

Greedily pack the Grassmannian in injection distance at least `r`.  The
number of `k`-subspaces at injection distance exactly `s` from a fixed one is

```math
2^{s^2}{k\brack s}_2{D-k\brack s}_2.             \tag{HI.7}
```

This follows by choosing the codimension-`s` intersection, the `s`
new quotient directions, and the graph map between the two `s`-dimensional
quotients.  Dividing the total Grassmannian by the closed radius-`r-1` ball
gives (HI.5).

For binary Gaussian coefficients,

```math
{n\brack j}_2\le4\,2^{j(n-j)},
\qquad
{D\brack k}_2\ge2^{k(D-k)}.
```

Thus the summand in (HI.5) is at most `16*2^(s(D-s))`.  Since
`s<=r-1<=(D/2)`, this exponent is increasing through the required range,
which proves (HI.6).  Every selected pair has injection distance at least
`r`; monotonicity of (HI.3) and Lemma HI.1 give Hamming Hausdorff distance
greater than `t`. `square`

### Corollary HI.3 (asymptotic exponent)

Let

```math
{k_D\over D}\longrightarrow\kappa,
\qquad
{r_D\over D}\longrightarrow\rho,
\qquad
{t_D\over D}\longrightarrow\delta,
```

where `0<rho<=kappa<=1/2`.  If

```math
{\delta\over\kappa+\rho}<\frac12
```

and

```math
(\kappa+\rho)
H_2\left({\delta\over\kappa+\rho}\right)<\rho,  \tag{HI.8}
```

then

```math
\liminf_{D\to\infty}{1\over D^2}
\log_2\operatorname{Pack}
(\operatorname{Gr}_{k_D}(F_2^D),d_Hs,t_D)
\ge
\kappa(1-\kappa)-\rho(1-\rho).                 \tag{HI.9}
```

#### Proof

The standard Hamming-ball asymptotic turns (HI.8) into (HI.4) for all
sufficiently large `D`.  Divide (HI.6) by `D^2` and pass to the limit.
`square`

This family generally uses many different sum codes; it is not obtained by
taking every carrier inside one fixed good linear code.

## 3. The construction cannot beat the common-host Gilbert exponent

The preceding transfer looks like a possible way around the separated-host
restriction.  At the leading exponent, it is not.

### Theorem HI.4 (entropy domination)

Under the hypotheses of Corollary HI.3,

```math
\boxed{
\kappa(1-\kappa)-\rho(1-\rho)
\le
\kappa\bigl(1-H_2(\delta)-\kappa\bigr).
}                                                  \tag{HI.10}
```

The right side is exactly the elementary common-host exponent: a binary
Gilbert code of rate `1-H_2(delta)-o(1)` and relative distance greater than
`delta` contains `2^(k(s-k)+o(D^2))` different `k`-subspaces.

#### Proof

Put

```math
a=\kappa+\rho,
\qquad
p={\delta\over a},
\qquad
h=H_2(p).
```

Condition (HI.8) says `ah<rho`.  Since `rho<=kappa`, one has
`rho<=a/2<=1/2`.  For fixed `a,p`, define

```math
F(x)=x(1-x)-(a-x)H_2(ap).
```

On `[0,1/2]`,

```math
F'(x)=1-2x+H_2(ap)>0.                            \tag{HI.11}
```

It is therefore enough to evaluate at `x=ah`.  Write

```math
g(a)=h(1-ah)-(1-h)H_2(ap).
```

Because `p<1/2`,

```math
g'(a)=-h^2-(1-h)pH_2'(ap)<0,
```

while `g(1)=0`.  Hence `g(a)>=0`, and

```math
F(ah)=a g(a)\ge0.
```

Monotonicity and `rho>ah` now give

```math
\rho(1-\rho)\ge\kappa H_2(\delta),
```

which is equivalent to (HI.10). `square`

Thus the following entire proof architecture is asymptotically sterile for
improving the known lower exponent:

1. infer separation from the number of low-weight words in `C+C'`;
2. impose a minimum injection distance;
3. use a Grassmannian Gilbert packing.

Any improvement beyond a common host has to use more than the scalar
low-weight count of the sum code--for example, the arrangement of its low
words relative to the two distinguished subcodes.

## 4. Common hosts are nevertheless not exactly optimal

The entropy barrier above should not be misread as an optimality theorem.
In `F_2^5`, consider the four two-dimensional subspaces

```math
C_1=\langle01100,10010\rangle,
\quad
C_2=\langle00100,10001\rangle,
```

```math
C_3=\langle00011,01000\rangle,
\quad
C_4=\langle00111,11001\rangle.                 \tag{HI.12}
```

Their directed Hausdorff distances, with the row giving the direction from
the lower-indexed subspace to the higher-indexed one, are

| pair | forward | reverse | Hausdorff |
|---|---:|---:|---:|
| `(1,2)` | 3 | 3 | 3 |
| `(1,3)` | 3 | 3 | 3 |
| `(1,4)` | 2 | 3 | 3 |
| `(2,3)` | 3 | 3 | 3 |
| `(2,4)` | 2 | 3 | 3 |
| `(3,4)` | 2 | 3 | 3 |

Hence they form a size-four packing at threshold `t=2`.  On the other hand,
a binary `[5,s,3]` linear code obeys the radius-one Hamming bound

```math
2^s(1+5)\le2^5,
```

so `s<=2`.  A common distance-three host therefore contains only one
two-dimensional carrier (and `C_4` itself shows that dimension two is
attainable).  Arbitrary carriers can beat the common-host
construction even though Theorem HI.4 shows that the most direct
injection-metric amplification of this fact cannot beat its asymptotic
Gilbert exponent.

## 5. Relation to the coding barrier

The exact `k=1` reduction to unrestricted binary codes and the systematic
chart upper bound are proved separately in
[`phase3_hamming_grassmannian_coding_barrier.md`](phase3_hamming_grassmannian_coding_barrier.md).
The present note adds two facts not contained there:

* a rigorous all-`k` lower construction using the low-weight geometry of
  pairwise sum codes rather than one common host;
* a matching entropy calculation proving why that construction cannot close
  the Hamming code--anticode gap.

Together with the four-carrier example, this isolates a more precise middle
invariant: the low-weight words of `C+C'` must be tracked **relative to both
distinguished subcodes**, not merely counted and not merely replaced by the
injection distance.
