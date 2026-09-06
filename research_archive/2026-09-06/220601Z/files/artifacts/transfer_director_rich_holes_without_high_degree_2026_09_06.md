# Actual bounded-operator signings: holes with exactly zero high degree

2026-09-06. Proved scalable counterexample to one possible overextension of
the first-marked theorem. This concerns actual hollow sign matrices, but
NOT near-minimizing matrices. It does not contradict the theorem for a
fixed response f(G,Y) of the two Gaussian-marginal first-marked variables.

Let k tend to infinity through powers of two, and take the symmetric
Sylvester Hadamard matrix H_k, so H_k^2=kI. Put n=2k and

```math
K=J_2\otimes H_k,\quad D=\operatorname{diag}(K),\quad
A=K-D,\quad B={A\over\sqrt{2k-1}},\quad Q=B^2,
\quad T=J_2\otimes I_k.
```

A is an exact symmetric hollow signing, and

```math
\|B\|_{op}\le{2\sqrt k+1\over\sqrt{2k-1}}=\sqrt2+o(1).
```

Since K^2=2kT and D^2=I,

```math
Q-T={T-KD-DK+I\over2k-1},\qquad
\|Q-T\|_{op}\le{4\sqrt k+3\over2k-1}=O(k^{-1/2}).    (1)
```

For uniform independent signs S, let i' be the matching coordinate in
the other k-block. Then `(TS)_i=S_i+S_i'`. Consider the FIXED odd ternary
response of the literal coherent query QS:

```math
F_i=\operatorname{sign}((QS)_i)\,\mathbf1_{|(QS)_i|>1},
\qquad U_i={S_i+S_{i'}\over2}.
```

U is an exact degree-one Boolean polynomial taking values in {-1,0,1}.
Whenever `|((Q-T)S)_i|<1`, F_i=U_i. Consequently

```math
{1\over n}\mathbb E\|F-U\|_2^2
\le {4\over n}\|Q-T\|_F^2
\le4\|Q-T\|_{op}^2=O(k^{-1}).                       (2)
```

This also handles equality at the threshold by charging the exceptional
event with `>=1`. Orthogonality of Boolean Fourier degrees now gives

```math
{1\over n}\sum_i\sum_{|E|\ge2}\widehat F_i(E)^2
\le {1\over n}\mathbb E\|F-U\|^2=O(k^{-1}),          (3)
\qquad {1\over n}\mathbb E\|F\|^2\longrightarrow{1\over2},
\quad {1\over n}\sum_i\mathbb E(1-F_i^2)
                                      \longrightarrow{1\over2}.
```

The last assertions follow from (2) and EU_i^2=1/2. Indeed both F and U
are ternary, so `|F_i^2-U_i^2|<=|F_i-U_i|`, and Cauchy--Schwarz suffices.

Thus a fixed rich-frame response can retain a positive density of both
marks and holes while its ENTIRE high-original-degree Fourier mass tends
to zero. The coherent query QS asymptotically exposes a pair of discrete
seed signs; it is not a Gaussian coordinate. A bounded-nonconstant-
polynomial impossibility argument on Gaussian space cannot be applied to it.

## Stronger exact Sylvester identity

For this particular family no error is needed: F=U IDENTICALLY at every
k=2^r>=2. Label rows of H by bit strings a and put d_a=H_aa=(-1)^|a|.
Exactly k/2 values of d_a are each sign. Set
`t_a=S_(1,a)+S_(2,a)`. Direct multiplication gives, for i=(ell,a),

```math
(A^2S)_i=(2k-2)t_a+S_i-R_a,
\quad R_a=\sum_{b\ne a}H_{ab}(d_a+d_b)t_b,
\quad |R_a|\le4(k/2-1)=2k-4.                         (4)
```

Only b with d_b=d_a contribute, and the diagonal term already removed is
`H_aa(2d_a)t_a=2t_a`. If t_a=0, (4) implies
`|(A²S)_i|<=2k-3`. If t_a=2s, then S_i=s and
`s(A²S)_i>=4k-3-(2k-4)=2k+1`. These lie strictly below and above the
threshold 2k-1 in F, respectively. Therefore

```math
F_i={S_i+S_{i'}\over2}\quad\hbox{for EVERY spin input},
\qquad \mathbb E F_i^2=\mathbb E(1-F_i^2)=1/2,
\qquad \widehat F_i(E)=0\quad(|E|\ge2).                (5)
```

The earlier perturbation proof is retained because it is robust to small
operator perturbations; the exact theorem (5) is the sharp counterexample.
`computations/transfer_director_rich_holes_exact_2026_09_06.py` checks all
spins at n=4,8,16 with integer arithmetic and rational Fourier masses.
Those finite checks are regressions, not the proof of the infinite family.

The non-near-minimizing scope can also be certified. Sylvester H_k has
trace zero for k>=2. Writing H_k^0 for its hollowing,
`H_A((z,z))=4H_(H_k^0)(z)`. Hence
`Q(A)/(2k)^(3/2)>=sqrt(2) M_k/k^(3/2)`. The already proved lower endpoint
above 1/(2sqrt(2)), together with the all-order upper below 1/2, separates
this family from asymptotic minimizers. No new lower result is needed.

The conditional mixed-nuclear theorem can still be correct when a
positive-variance high-degree channel P is explicitly present. This example
shows that selecting such a channel is another obligation for arbitrary
rich responses, separate from mixed covariance and returned-query closure.
It does not assert that the current optimized rich policy has this defect,
or that exact minimizers admit this particular macroscopic twin structure.
