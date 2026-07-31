# Switching-template obstruction

Status: proved. This note is a self-contained proof artifact for the second
sustained computational--composition campaign.

## Definition

Call an order-`n` signing `A` a **switching `K`-template** if there are

- a partition `V=V_1 union ... union V_r` with `r<=K`;
- switching signs `d_i in {+1,-1}`; and
- symmetric template signs `sigma_ab in {+1,-1}` for `1<=a<=b<=r`

such that, for distinct `i in V_a` and `j in V_b`,

~~~math
a_{ij}=d_i d_j\sigma_{ab}.
~~~

Thus the definition permits the best partition, switching, and signed
complete quotient; it does not refer to `M_n`.

## Theorem

For `n>=4`, every switching `K`-template satisfies

~~~math
\operatorname{cap}(A)
:=\max_{x\in\{+1,-1\}^n}
\left|\sum_{i<j}a_{ij}x_ix_j\right|
\ge \frac{n^2}{8K}.                                \tag{T1}
~~~

Consequently, if `F_n(K_n)` is the family of all switching `K_n`-templates
and `K_n=o(n^(1/2))`, then

~~~math
\frac{1}{n}\left(
 \min_{A\in F_n(K_n)}\operatorname{cap}(A)^{2/3}
 -M_n^{2/3}\right)\longrightarrow+\infty.          \tag{T2}
~~~

In particular this entire family has a superlinear, rather than `o(n)`,
landing gap in the `b_n=M_n^(2/3)` scale.

## Proof of (T1)

Write `s_a=|V_a|` and let `t=max_a s_a`. First, cap is monotone under
principal restriction. Indeed, fix spins on a subset and choose all outside
spins independently and uniformly. The conditional mean of the full energy
is the energy on the subset, so one extension has absolute full energy at
least the absolute restricted energy. The largest cell becomes a constant
signing after switching, hence

~~~math
\operatorname{cap}(A)\ge {t\choose2}.              \tag{T3}
~~~

There is also a quotient-variance bound. Restrict the Boolean spins to
`x_i=d_i y_a` on `V_a`, where the `y_a` are independent uniform signs. The
resulting energy is

~~~math
H(y)=W+\sum_{a<b}\sigma_{ab}s_as_b y_ay_b,
\qquad
W=\sum_a\sigma_{aa}{s_a\choose2}.
~~~

The distinct degree-two characters `y_a y_b` are orthogonal. Therefore

~~~math
\mathbb E H(y)^2
=W^2+\sum_{a<b}s_a^2s_b^2,
~~~

and hence

~~~math
\operatorname{cap}(A)
\ge\left(\sum_{a<b}s_a^2s_b^2\right)^{1/2}.        \tag{T4}
~~~

Empty cells may be discarded, so assume `K<=n`. Put
`S_2=sum_a s_a^2` and `S_4=sum_a s_a^4`. If
`t^2>=n^2/(2K)`, then `t>=2` and (T3) gives

~~~math
\operatorname{cap}(A)\ge\frac{t^2}{4}
\ge\frac{n^2}{8K}.
~~~

If instead `t^2<n^2/(2K)`, Cauchy--Schwarz gives `S_2>=n^2/K`, while
`S_4<=t^2S_2`. Consequently

~~~math
\sum_{a<b}s_a^2s_b^2
=\frac{S_2^2-S_4}{2}
\ge\frac12S_2\left(S_2-\frac{n^2}{2K}\right)
\ge\frac{n^4}{4K^2}.
~~~

Equation (T4) then gives `cap(A)>=n^2/(2K)`, which is stronger than (T1).

## Proof of (T2)

The elementary random-signing argument gives `M_n=O(n^(3/2))`: for each
fixed Boolean spin vector, the energy of a uniformly random signing is a sum
of `binom(n,2)` independent signs, and a Hoeffding bound followed by a union
bound over the `2^n` spin vectors supplies a universal constant `C` with
`M_n<=C n^(3/2)`.

On the other hand, (T1) gives

~~~math
\min_{A\in F_n(K_n)}\operatorname{cap}(A)^{2/3}
\ge 8^{-2/3}\frac{n^{4/3}}{K_n^{2/3}}.
~~~

After division by `n`, the right side tends to infinity when
`K_n=o(n^(1/2))`, whereas `M_n^{2/3}/n=O(1)`. This proves (T2).

## Research meaning and limit

The result falsifies switching-step states of bounded or sub-square-root
complexity as a landing family. It also explains why merely giving a fixed
signed quotient, bounded number of orbit types, or bounded equitable
partition cannot solve structured landing.

It does **not** rule out states with at least order `n^(1/2)` types or states
whose cross-cell blocks carry substantial pseudorandom structure. Allowing
arbitrary internal signings alone is addressed by the balanced extension
below.

## Balanced arbitrary-child extension

The cross-block part of the proof does not require constant signs *inside*
the cells. Suppose there are `K>=2` cells, every cell has size at least
`alpha n/K`, the signing inside each cell is arbitrary, and every cross-cell
block has the switching-rank-one form

~~~math
a_{ij}=d_i d_j\sigma_{ab}qquad
(i\in V_a, j\in V_b, a\ne b).
~~~

On the restricted spins `x_i=d_i y_a`, each internal-cell energy is a
constant independent of `y_a`; the cross-cell degree-two characters remain
orthogonal. The same second-moment calculation therefore gives

~~~math
\operatorname{cap}(A)
\ge\left(\sum_{a<b}s_a^2s_b^2\right)^{1/2}
\ge {\alpha^2 n^2\over 2K}.                          \tag{T7}
~~~

Thus any such balanced construction with cap `O(n^(3/2))` needs
`K=Omega(sqrt(n))`, even if every diagonal block is itself an arbitrary exact
minimizer. In particular, a bounded number of comparable good children
cannot be composed with switching-rank-one bridge blocks. A viable recursive
composition must put genuine pseudorandom/Boolean complexity in the bridges,
not merely in its child signings.

## A separate universal edit-net obstruction

There is another elementary limit on a possible structured family. If two
signings differ on `h` edges, then for every Boolean spin vector their
energies differ by at most `2h`. Therefore

~~~math
|\operatorname{cap}(A)-\operatorname{cap}(B)|\le2h. \tag{T5}
~~~

Thus an `o(n^(3/2))`-radius Hamming approximation to an optimizer would be a
sufficient landing mechanism. It cannot, however, be obtained from a small
**universal** edit net. Put `L=binom(n,2)`. A family covering every signing at
Hamming radius `r=o(n^(3/2))` must have size at least

~~~math
\frac{2^L}{\sum_{j=0}^r{L\choose j}}
=2^{L-o(n^2)}.                                      \tag{T6}
~~~

Indeed, the usual Hamming-ball estimate gives
`log(sum_(j<=r) binom(L,j))=O(r log(eL/r))=o(n^2)`.
Equation (T6) rules out generic edit rounding to a low-entropy catalogue. It
does not rule out a much smaller family specially arranged to contain a
nearby point for at least one optimizer at every order; proving that special
intersection is exactly the non-universal information a landing theorem must
supply.
