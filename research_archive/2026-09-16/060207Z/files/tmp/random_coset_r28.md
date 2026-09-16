# Wave 28: random block-coset hits for one fixed selector

## Scope and conclusion

Fix an exact order-`n` minimizer `A`, an `m`-selector `S`, and

```math
p_2=\frac{m(m-1)}{n(n-1)},\qquad
H=H_S=A\circ(\xi_S\xi_S^{\mathsf T}-p_2\mathbf1\mathbf1^{\mathsf T}).
```

Put `T=(Y_A(S)-t)_+`.  This memo studies only the probability

```math
\Pr\{Q(P^{\mathsf T}DHD P)\ge T\}
```

under random vertex signs `D`, and optionally a random block partition `P`.
There are three rigorous conclusions.

1. For a fixed partition the hit event is exactly the saturation of the
   scalar tail set by the block-sign subgroup.  A `2^k`-word coset can amplify
   the sign-symmetric scalar tail by at most `2^(k-1)`, and randomizing the
   partition does not change this bound.
2. Hanson--Wright, with the exact norms of `H_S`, shows that uniform random
   diagonals cannot give the desired `exp(-O(n^(3/4-c)))` hit probability when
   `T=omega(n^(3/2-c))`.  In particular a macroscopic
   `T >= eta n^(3/2)` gives only `exp(-Omega(n^(3/4)))`, even after the whole
   coset is searched.  The boundary `T=O(n^(3/2-c))` is not ruled out: at that
   boundary the Hanson--Wright exponent is exactly the desired scale.
3. For one fixed scalar-hit cut `x` which already has
   `R_2(x)=O(n^(9/4-c))`, a random hash partition makes the *entire* coset
   row-good with high probability.  Thus, for one selector, a sufficiently
   dense joint scalar set `{hit and low row}` can be converted into a row-good
   random-coset hit.  This does not prove that this set is dense, and does not
   yield one partition/law uniformly for every selector.

The result is a sharp obstruction to the **uniform random-diagonal law**, not
a falsification of (10.838): (10.838) permits an `A`-adapted law supported on
row-good cosets, common to all selectors, and such a law need not resemble the
uniform quotient law.

## 1. Exact quotient identity

For a partition `P` into `k` blocks, let

```math
K_P=\{Pz:z\in\{\pm1\}^k\}\subset\{\pm1\}^n.
```

Products are coordinatewise.  The subgroup `K_P` has size `2^k`, and a
uniform diagonal `D=diag(epsilon)` chooses a uniform coset `epsilon K_P`.
Define the scalar tail set

```math
E_T=\{x\in\{\pm1\}^n:|x^{\mathsf T}Hx|\ge T\},
\qquad \pi_T=2^{-n}|E_T|.
```

Directly from (10.832),

```math
\boxed{
\Pr_D\{Q(P^{\mathsf T}DHD P)\ge T\}
=\frac{|E_TK_P|}{2^n}.
}
\tag{R28.1}
```

Indeed, the event says exactly that `epsilon Pz` lies in `E_T` for some `z`.
Consequently

```math
\boxed{
\pi_T\le \Pr_D\{Q(P^{\mathsf T}DHD P)\ge T\}
\le \min\{1,2^{k-1}\pi_T\}.
}
\tag{R28.2}
```

The factor is `2^(k-1)`, rather than `2^k`, because `E_T=-E_T` and the
two global-sign translates coincide.  Both constants are exact: the lower
bound is equality if `E_T` is already a union of cosets, while the upper
bound is equality when all effective translates before saturation are
disjoint.  Since (R28.2) holds for every `P`, it also holds after averaging
over any law of balanced or unbalanced partitions.

Let `G_P` be the row-good event used in (10.833).  It too is constant on
`K_P`-cosets, because right multiplication of `ADP` by a diagonal block sign
does not change its operator norm or the set of codewords.  The proof of
(10.833) gives `Pr_D(G_P)>=1/2`.  Therefore

```math
\Pr_D(hit\mid G_P)\le 2\Pr_D(hit),
\tag{R28.3}
```

but the only general lower bound is
`Pr(hit | G_P) >= max(0,2 Pr(hit)-1)`.  For a stretched-exponentially rare
hit this is zero.  Thus conditioning the uniform diagonal law on row goodness
can annihilate all hits; the half-probability row theorem alone gives no
positive rare-event intersection theorem.

## 2. Exact norms and the random-diagonal upper bound

Let `N=n(n-1)` and `a=m(m-1)=Np_2`.  Since `A` is zero-diagonal with
off-diagonal entries of modulus one,

```math
\boxed{
\|H_S\|_F^2
=a(1-p_2)^2+(N-a)p_2^2
=Np_2(1-p_2).
}
\tag{R28.4}
```

Writing `J_S` for coordinate projection onto `S`,

```math
H_S=J_SAJ_S-p_2A,
\qquad
\|H_S\|_{op}\le(1+p_2)\|A\|_{op}.
\tag{R28.5}
```

For an exact minimizer, the previously verified bound
`||A||_op^2 <= 2q_n` gives

```math
\|H_S\|_{op}\le(1+p_2)\sqrt{2q_n}=O(n^{3/4}),
\qquad \|H_S\|_F\le n.
\tag{R28.6}
```

Use the following normalized Hanson--Wright inequality: there is one
universal numerical constant `c_HW>0` such that, for a Rademacher vector
`epsilon` and a symmetric zero-diagonal matrix `H`,

```math
\Pr\{|\epsilon^{\mathsf T}H\epsilon|\ge u\}
\le 2\exp\left[-c_{HW}\min\left{
\frac{u^2}{\|H\|_F^2},\frac{u}{\|H\|_{op}}
\right}\right].
\tag{R28.7}
```

Here `||.||_F` counts both orientations, exactly as in (R28.4), and the
quadratic form has no factor `1/2`.  Applying (R28.7) to every one of the at
most `2^k` block words gives the fully normalized bound

```math
\boxed{
\Pr_D\{Q(P^{\mathsf T}DH_SDP)\ge T\}
\le 2^{k+1}\exp\left[-c_{HW}\min\left{
\frac{T^2}{n(n-1)p_2(1-p_2)},
\frac{T}{(1+p_2)\sqrt{2q_n}}
\right}\right].
}
\tag{R28.8}
```

This holds for every fixed partition, so also for jointly random `(P,D)`.
Conditioning on the row-good event from (10.833) changes its right side by at
most a factor two, by (R28.3).

Now fix `0<c<1/4`, put `L_n=n^(3/4-c)`, and take (10.835), so
`k=O(L_n/log n)=o(L_n)`.  Uniformly in the active ratio window, if

```math
T=\omega(n^{3/2-c}),
```

then both exponents in (R28.8) are `omega(L_n)`; for the quadratic exponent
use `1-2c > 3/4-c`, which is exactly where `c<1/4` enters.  Hence

```math
\Pr_D(hit)\le\exp\{-\omega(L_n)\}.
\tag{R28.9}
```

This is incompatible with a lower bound `exp(-O(L_n))`.  In particular,
for fixed `eta>0`, `T>=eta n^(3/2)` gives
`Pr_D(hit)<=exp(-Omega(n^(3/4)))`.

At the sharp boundary `T=K n^(3/2-c)`, however, the coarse uniform spectral
bound guarantees only a `Theta(KL_n)` lower bound on the linear exponent,
while `k log 2=o(L_n)`.  Thus (R28.8), at the level justified without
further structure on `H_S`, is compatible with `exp(-O(L_n))`; it neither
proves nor disproves the desired lower bound there.  Random-diagonal success
needs a minimizer-specific **lower** tail at the exact aligned level
`T=(Y_A(S)-t)_+`, not another upper concentration theorem.

Equation (R28.2) also makes the needed entropy explicit.  A uniform random
coset hit probability at least `exp(-C L_n)` necessarily requires

```math
\boxed{
\pi_{(Y_A(S)-t)_+}
\ge 2^{-k}e^{-CL_n}=e^{-O(L_n)}.
}
\tag{R28.10}
```

Conversely this scalar density gives the same unconditional lower bound for
every partition, because the word `z=1` is in the coset.  It does **not** by
itself survive row-good conditioning.  Notice that (R28.10) is explicitly
aligned with `Y_A(S)`; a sign-tail at an unrelated centered scale would say
nothing.

## 3. Random hashing repairs row goodness once a low-row hit exists

There is a useful fixed-selector joint statement.  Fix a spin `x` and assign
each vertex independently to one of `k` labeled blocks.  Let `P` be the
incidence matrix and `v_i=x_iAe_i`.  Then

```math
AD_xP=M+W,
\qquad
M=\frac1k(Ax)\mathbf1_k^{\mathsf T},
\qquad
W=\sum_i v_i(e_{g(i)}-k^{-1}\mathbf1_k)^{\mathsf T}.
```

The summands of `W` are independent and centered.  Their two rectangular
variance matrices and uniform norm bound are exactly

```math
\begin{aligned}
\sum_i\mathbb E Z_iZ_i^{\mathsf T}&=(1-k^{-1})A^2,\\
\sum_i\mathbb E Z_i^{\mathsf T}Z_i
&=n(n-1)(k^{-1}I-k^{-2}J),\\
\|Z_i\|_{op}&\le\sqrt{n-1}.
\end{aligned}
\tag{R28.11}
```

Set

```math
v=\max\{\|A\|_{op}^2,n(n-1)/k\},
\qquad u=\log((n+k)/\eta).
```

Rectangular matrix Bernstein gives, with probability at least `1-eta`,

```math
\|W\|_{op}\le \sqrt{2vu}+\frac23\sqrt n\,u.
\tag{R28.12}
```

Since `||M||_op^2=R_2(x)/k`, every word in this random coset then obeys

```math
\boxed{
\max_zR_2(D_xPz)
\le 2R_2(x)+8kvu+\frac{16}{9}knu^2.
}
\tag{R28.13}
```

For (10.835), exact-minimizer spectral control implies `v=O(n^(3/2))`.
Taking `eta=n^-3` gives

```math
R_2(x)=O(n^{9/4-c})
\quad\Longrightarrow\quad
\max_zR_2(D_xPz)=O(n^{9/4-c})
\tag{R28.14}
```

with probability `1-o(1)`.  Standard binomial Chernoff bounds simultaneously
give nonempty blocks and maximum block size at most `2n/k` with probability
`1-o(1)`, because `n/k >> log n`.  Thus this random hash is balanced to the
order required in (10.835).

For the fixed selector, define the genuinely aligned joint set

```math
E_{T,R}=\{x:|x^{\mathsf T}H_Sx|\ge T,
R_2(x)\le C_Rn^{9/4-c}\}.
```

If `2^-n |E_{T,R}| >= exp(-O(L_n))`, averaging (R28.14) first over `P`
and then over `x` shows that some balanced-order partition has row-good hit
mass `exp(-O(L_n))`.  Conditioning its uniform quotient law on row-good
cosets retains at least that mass.  This is a valid fixed-selector conversion.

What remains open is precisely the hard part:

- no proved minimizer property gives the required lower bound on
  `|E_{T,R}|` at `T=(Y_A(S)-t)_+`;
- the partition obtained by averaging may depend on `S`; polynomial failure
  probability in (R28.14) cannot be union-bounded over all
  `binom(n,m)=exp(Theta(n))` selectors without enlarging the row budget; and
- an arbitrary `A`-adapted law in (10.838) could outperform the uniform
  quotient law, so (R28.8)--(R28.9) do not falsify the strategic lemma.

## 4. Proved, numerical, and open status

**Proved:** (R28.1)--(R28.14), subject only to the standard universal-constant
Hanson--Wright and rectangular matrix Bernstein inequalities in the explicitly
stated normalizations.  The quotient identities and exact norms were also
checked by exhaustive small random instances in `tmp/random_coset_r28_check.py`.

**Numerical:** only the small-instance checker; no asymptotic conclusion rests
on it.

**Open:** a minimizer-specific lower bound for the joint set `E_{T,R}`, or a
single `A`-adapted partition and row-good law giving that bound for every
selector.  Uniform random diagonals have now been sharply delimited: they can
only be relevant once the aligned threshold is already at most the
`n^(3/2-c)` scale, and even there an aligned reverse tail is indispensable.
