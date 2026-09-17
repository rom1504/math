# Low-entropy Boolean codes admit a favorable critical sign frame

2026-09-17. **Verified:** director proof and independent discrepancy-track
reconstruction of every finite step and normalization. This theorem
combines finite-code information with the common-Gibbs fourth-order
replacement theorem. No external novelty claim is established. The
result concerns restricted maxima, not the unrestricted parent cap.

## 1. Finite frame-selection theorem

Let n=kp+ell, 0<=ell<k, where a real k by k Hadamard matrix exists.
Let C be a nonempty subset of {+-1}^n, and let H=log|C|. Choose
q mode labels h_j from its rows by balanced repetition, so no row is
used more than b=ceil(q/k) times. For every u>0, with probability at
least 1-exp(-u), a uniformly random diagonal sign matrix D satisfies
SIMULTANEOUSLY for all x in C

```math
 W_D(x):=\sum_{j=1}^q\sum_{t=1}^p
          \langle h_j,(Dx)_{:,t}\rangle^4+q\ell
 \le 4b k^3\{H+p\log(\sqrt2 k)+u\}+q\ell.       \tag{1}
```

The same D works for the entire code. Its entries may depend on C;
there is no claim of one frame working for all codes simultaneously.
The operation changes the bridge frame, NOT the internal child edges.

### Proof

For each column t and fixed x, z_t=(Dx)_{:,t} is an independent uniform
sign vector of length k. Orthogonality gives
sum_h <h,z_t>^2=k^2. Set L_t=max_h <h,z_t>^2/k^2. Then

```math
 \sum_h\langle h,z_t\rangle^4\le k^4 L_t.          \tag{2}
```

For any row h, S=<h,z_t> has MGF E exp(aS)<=exp(ka^2/2).
An auxiliary scalar Gaussian consequently gives

```math
 \mathbb E e^{S^2/(4k)}
   =\mathbb E_G\mathbb E e^{G S/\sqrt{2k}}
   \le\mathbb E_G e^{G^2/4}=\sqrt2.
```

Taking the maximum over the k rows costs at most their sum, hence
E exp(k L_t/4)<=sqrt(2)k. Independence of the p columns and
Markov give

```math
 \Pr\left\{\sum_tL_t>\frac4k
       [H+p\log(\sqrt2k)+u]\right\}\le e^{-H-u}.
```

Union over the code, use (2), and pay at most b repeated copies of
each mode row. Leftover coordinates have coefficient magnitude one.
This proves (1), with exact sign entries and no Gaussian physical edges.

## 2. All-offset restricted-parent universality

Take arbitrary full sign children A_n and B_q. Choose D as in (1),
and construct a physical bridge whose jth core column is
D(h_j tensor g_j), where g_j=(g_(j,1),...,g_(j,p)) has independent
scalar fair-sign coordinates, independently for all j;
leftover entries are also independent fair signs. The covariance bound
is unchanged by D, and every bridge entry is a sign.

Restrict the old spin x to C, retain ALL new spins y and both absolute
polarities, and retain the exact offsets H_A(x)+H_B(y). Replacing only
the scalar bridge drivers by independent standard Gaussians changes the
expected restricted maximum by at most

```math
 8\{4b k^3[H+p\log(\sqrt2k)+u]+q\ell\}^{1/4}
       [(q+1)\log2+H+\log2]^{3/4}.                \tag{3}
```

This follows directly from the proved truncated common-Gibbs theorem
in [the universality artifact](paper_symmetric_frame_universality_2026_09_17.md).
There are K=2^(q+1)|C| witnesses; the logarithm is log(2K).
For each witness, multiplication by its new spin and polarity preserves
the fourth coefficient sum W_D(x). There is no separate payment of
left/right channels and no optimization-dependent omitted offset.

When q=Theta(n), k=Theta(sqrt n), and H=o(n), choose u=o(n) tending
to infinity. The right side of (3) is o(n^(3/2)). Powers of two supply
Hadamard orders with k between fixed multiples of sqrt n at EVERY n;
the ell<k leftovers are already included. Thus subexponential old code
size is SUFFICIENT to choose a critical-dimensional physical sign frame
with vanishing universality error on that entire code.

The formerly exhibited coherent ground word only falsifies uniformity
over ALL frames. It does not falsify favorable frame selection (1).

## 3. Stability under Hamming coverings

Suppose C lies within Hamming distance r n of a center code F. Choose
D by (1) for F. Minkowski in the finite l4 coefficient space gives

```math
 \sup_{x\in C}W_D(x)^{1/4}
 \le \sup_{f\in F}W_D(f)^{1/4}
       +[16(q+k)k^2rn+16q\ell]^{1/4}.              \tag{4}
```

For completeness, if x and f differ at r_t core coordinates in column
t, their difference v has |<h,v>|<=2r_t and sum_h<h,v>^2=4kr_t.
Thus sum_h<h,v>^4<=16k r_t^3. Repetition gives at most
16b k sum_t r_t^3<=16(q+k)k^2 sum_t r_t, since r_t<=k.
Leftovers contribute at most16qell. The difference coefficients are
linear before taking l4, so (4) includes them by ordinary Minkowski.

Consequently a code covered at radius r_n n, r_n->0, by exp(o(n))
centers also admits vanishing critical-scale universality error.
This statement neither claims cheap mean absolute response nor uses it.

## 4. Application and exact remaining obligations

The independently reconstructed same-order sign-star regularizer makes
an arbitrary signing asymptotically cap-preserving while producing a
subexponential microscopic nearcode. Therefore (1)--(3) apply to actual
regularized exact minimizers, at every order, rather than only a chosen
algebraic family. Using the q0~n^(1/3) regularizer gives

```math
 H=O(n^{5/6}\sqrt{\log n}),\quad
 W=O(n^{17/6}\sqrt{\log n}),\quad
 \text{comparison error}=O(n^{35/24}(\log n)^{1/8}).
```

The last exponent is 17/24+3/4=35/24<3/2. Here q0 is the small
regularizing star, while q=Theta(n) is the prospective new child.
They are different parameters and cannot be interchanged.

This removes a PARTICULAR coefficient-concentration obligation on the
actual microscopic nearcode: a good frame can be selected, not assumed.
It does not control states outside that window and does not bound the
Gaussian restricted-parent maximum by the desired child optimum.
The full hybrid certificate still needs those two value controls.
No convergence, recurrence, or improved original constant follows yet.
