# Walsh divisibility gives an exponentially rare selector class

Date: 2026-09-06. Independently audited exact finite statements and scopes.
The finite statements are followed by their scoped use in
the restricted weave. No near-plateau stability theorem is assumed.

Let m=2^d, and index the Walsh matrix by F_2^d, with
H(y,x)=(-1)^{x dot y}. For a ternary vector xi, put h=H xi and
T=support(xi). Logarithms in asymptotic rates are natural.

## 1. Divisibility forces low degree of the support, not of its signs

If 2^r divides every h(y), with 1<=r<=d, the Boolean indicator 1_T,
viewed as an F_2 polynomial, has algebraic degree at most d-r.

Proof. Modulo 2, xi equals 1_T. The coefficient of the squarefree monomial
indexed by J is the sum modulo 2 of xi over the coordinate subcube where
coordinates outside J vanish. Fourier inversion gives this integer sum as

```math
  \sum_{x:\ x_{J^c}=0}\xi(x)
    =2^{|J|-d}\sum_{y:\ y_J=0}h(y).
```

If |J|>d-r the right side is divisible by 2. Thus all such algebraic
coefficients vanish. This proves the statement for EVERY choice of the
signs on T. It is only a necessary condition on T, not a sufficient
condition for a signing with the prescribed divisibility.

## 2. A finite low-weight count from a translated information set

Let RM(s,d) be the binary degree-at-most-s evaluation code, of dimension

```math
 K=\sum_{j=0}^s\binom dj.
```

The Hamming ball I={x: |x|<=s} is an information set of size K. Indeed
values at these points determine all algebraic coefficients of degree
at most s by Boolean Mobius inversion. Translation of the arguments
preserves degree, so every I+v is also an information set.

The number of codewords of weight at most delta*m, for 0<=delta<1/2,
is bounded by

```math
   m\sum_{j\le\lfloor\delta K\rfloor}\binom Kj.          (1)
```

For any one such word, the mean number of its ones on I+v, averaged over
all m translations, is at most delta*K. At least one translation has
at most floor(delta*K) ones. The values there uniquely specify the whole
word; union over v and possible information words proves (1).

Consequently, for p>1/2 and a uniform k=pm selector (integer when used),

```math
 Pr_T\{\exists\xi:\operatorname{supp}\xi=T,
                   \ 2^r\mid (H\xi)(y)\ \forall y\}
 \le \frac{m\sum_{j\le(1-p)K}\binom Kj}{\binom mk},
 \qquad s=d-r.                                         (2)
```

Apply (1) to the complement of T, which has degree at most s whenever
1_T does. If s=d/2+O(1), then K/m->1/2, and (2) has exponent at most

```math
  -\tfrac12 h(p)m+o(m).                                (3)
```

This proves an exponentially rare necessary selector class for every
spin row whose entire unnormalized spectrum has a power-of-two divisor
of order sqrt(m). No count of the sign choices on an accepted selector
has been inferred.

## 3. Hamming neighborhoods can also be removed

Let B_m be the set of k-selectors in (2), or any subset of the degree-s
support class. If G_m excludes all k-selectors obtainable from B_m by
at most alpha*m deletions and the same number of insertions, then

```math
 Pr(G_m^c)\le
 \frac{m\sum_{j\le(1-p)K}\binom Kj}{\binom mk}
 \sum_{a\le\lfloor\alpha m\rfloor}\binom ka\binom{m-k}a.
                                                               (4)
```

For fixed alpha<p(1-p), the extra exponential rate is

```math
 p h(alpha/p)+(1-p)h(alpha/(1-p)).                       (5)
```

Thus, if s=d/2+O(1) and (5)<h(p)/2, this entire fixed-radius neighborhood
still has exponentially small probability. Independent fibre selectors
may be conditioned to avoid it, without changing the exponential rate
in the one-row permanent criterion: log Pr(G_m)/m -> 0.

This is a real, testable condition on selectors alone. It does not inspect
the unknown parent cap and does not store every child spin response.

## 4. Why exact divisibility is not a spectral stability theorem

If a nonzero supported xi has every Walsh coefficient divisible by 4,
flip one nonzero coordinate. Then

```math
 h'(y)=h(y)-2\xi(x_0)(-1)^{x_0\cdot y}.
```

Every coefficient of h' is congruent to 2 modulo 4. Nevertheless, with
k=|T|, its normalized empirical squared spectral displacement is

```math
 \frac1m\sum_y\left|\frac{h'(y)-h(y)}{\sqrt k}\right|^2
       =4/k\longrightarrow0.
```

The support has not changed, so this observation does NOT defeat the
selector filter: it still rejects that T because the original xi exists.
It does defeat any inference that a nearby spin spectrum itself retains
the exact divisibility. More generally, the present proof does not show
that every approximately plateaued spectrum has a support close to the
low-degree support class. Such a repair theorem would require a new
argument; neither Fourier inversion modulo 2 nor the count (1) gives it.

## 5. Operational conclusion and remaining obligation

Intersecting the good event here with the anti-sparsity event in
`continued_feedback_selector_antisparsity_2026_09_06.md` is legitimate.
The bad probabilities add and remain exponentially small. This removes
several rigorously described exceptional selector families at negligible
conditioning cost.

It does NOT establish the needed upper bound on the truncated one-row
partition sum. Dangerous nondivisible, nonsparse spectral profiles and
their multiplicities remain uncontrolled. The original convergence
question and its rigorous numerical interval are unchanged.
