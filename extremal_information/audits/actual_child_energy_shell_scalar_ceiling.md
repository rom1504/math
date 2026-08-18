# Scalar energy data do not control negative-path shell retuning

Status: **rigorous generic obstruction, exact actual finite response
falsifier, and narrowed optimizer-specific observable**.  This note is
separate from the exact shell decomposition theorem.  It asks whether the
shell KL can be bounded sublinearly using only scalar child pressure, cap,
energy histogram, and Hamming-sphere spread data.

The answer is negative at that level of information.  Those data give only
a linear ceiling.  A product binary-channel construction has polynomially
many energy shells, exponential atom spread, and nevertheless linear shell
KL.  More sharply, two certified actual finite pressure minimizers with the
same complete signed energy histogram have different negative-path shell
responses.  The exact missing observable is the shell-resolved inverse
output moment, not another scalar child-energy statistic.

**Scope warning.** Shell retuning is only one coarse diagnostic of the
negative-disorder interaction. Even a theorem with `bar p=p` (and hence
zero shell KL) would not by itself control the conditional law inside each
shell or prove that `J_N-I_N^{\leftarrow}` is sublinear. In particular, a coherent
row-factor retuning mechanism can in principle preserve the shell marginal
while changing the conditional factor geometry extensively. That stronger
no-go is being audited separately; nothing below promotes sublinear shell KL
to the full balanced-product conclusion.

## 1. The scalar inequalities stop at `O(N)`

Use the notation of Theorem 37.62.  Its pointwise surprise estimate is

```math
p_e\ge e^{-C_\beta N},
\qquad C_\beta=\log2+{\beta^2\over4}.                 \tag{SC.1}
```

For an arbitrary retuned shell law `bar p`, this gives only

```math
\boxed{
D(\bar p\Vert p)
=E_{\bar p}[-\log p_e]-H(\bar p)
\le C_\beta N.}                                      \tag{SC.2}
```

The quadratic spread theorem does not improve this direction.  It upper-
bounds individual prior atoms and narrow-cap masses, whereas (SC.2) needs a
stronger lower bound on every shell mass or a theorem preventing the
negative output tilt from selecting rare shells.  The scalar energy
histogram contains no information about that selection.

## 2. A scalable generic obstruction

The obstruction already occurs in a central two-bit binary channel.  Let

```math
q_0^\pm=\mathord\pm(1,1),
\qquad q_1^\pm=\mathord\pm(1,-1),
```

put `E(q_0^\pm)=0`, `E(q_1^\pm)=1`, and take the centrally symmetric prior

```math
\mu(q_0^+)=\mu(q_0^-)={2\over5},
\qquad
\mu(q_1^+)=\mu(q_1^-)={1\over10}.                    \tag{SC.3}
```

Relative to fair `B\in\{\pm1\}^2`, use the product binary channel

```math
k(B\mid Q)=\prod_{j=1}^2\left(1+{1\over2}B_jQ_j\right). \tag{SC.4}
```

Writing `z=B_1B_2`, the output density and the posterior probability of
shell one are

```math
p(B)=1+{3\over20}z,
\qquad
P(E=1\mid B)=
\begin{cases}
3/23,&z=1,\\
5/17,&z=-1.
\end{cases}                                           \tag{SC.5}
```

For the negative escort `dq_(-1)/dU proportional p^(-1)`, one has

```math
q_{-1}\{z=1\}={17\over40},
\qquad
q_{-1}\{z=-1\}={23\over40}.
```

Therefore the averaged posterior shell mass is exactly

```math
\boxed{
\bar w={17\over40}{3\over23}
       +{23\over40}{5\over17}
       ={439\over1955}\ne {1\over5}=w.}              \tag{SC.6}
```

Take `k` independent copies, declare the total energy to be the number of
one-shell blocks, and retain the product channel.  Then:

1. there are only `k+1` combined-energy shells;
2. the prior and averaged posterior are uniform conditional on each shell;
3. the largest latent atom is `(2/5)^k`, so min-entropy is linear;
4. the negative escort and averaged posterior factor blockwise; and
5. the two shell laws are respectively `Bin(k,w)` and `Bin(k,bar w)`.

Consequently

```math
\boxed{
D(\operatorname {Bin}(k,\bar w)
  \Vert\operatorname {Bin}(k,w))
=k\,d_0,
\qquad
d_0=D(\operatorname {Ber}(439/1955)
      \Vert\operatorname {Ber}(1/5))
=0.001829590523564\ldots>0.}                         \tag{SC.7}
```

This is a generic inverse binary-channel model, not an actual dense
quadratic child.  Its role is exact: polynomial shell state, central
symmetry, shell uniformity, and exponential atom spread do not imply
sublinear shell KL.  Some optimizer-specific relation between energy shells
and channel geometry is indispensable.

## 3. An actual finite scalar-data falsifier

Let `A_0,A_1` be the two certified order-eight pressure-minimizer classes
from FC.22--FC.24.  They have the same exact signed energy histogram

```text
-10:4, -8:10, -6:12, -4:16, -2:16, 0:12,
  2:16, 4:16, 6:12, 8:10, 10:4,
```

the same cap `10`, and the same scalar pressure at every temperature.  The
exhaustive cap classification proves that both are actual pressure
minimizers for every raw `t>=3`.  Attach the unique order-two child, take
orientation `epsilon=1`, raw temperature `t=3` (scaled parent temperature
`beta=3sqrt(10)`), and enumerate all `2^16` bridges.

The two combined-energy priors are identical, but under
`q_(-1) proportional p^(-1)` their averaged shell laws differ.  Numerically,

```math
\begin{array}{c|cc}
 &A_0&A_1\\ \hline
\bar p(-11)&0.21477821453936\ldots&0.24439189667794\ldots\\
D(\bar p\Vert p)&3.86163734362958\ldots&2.85355036639597\ldots
\end{array}                                           \tag{SC.8}
```

The response separation has an exact certificate, not only a floating-
point one.  Put `z=e^t`.  For disorder exponent `-1`, every shell response
is a rational function in `z` with rational coefficients, because
`cosh(kt)=(z^k+z^(-k))/2`.  Exact bridge grouping and rational arithmetic at
`z=2` give different values for the `e=-11` response.  Hence the two
rational functions are distinct.  Since `e^3` is transcendental, their
values cannot coincide at the actual `t=3`.  The KL values in (SC.8) are
stable long-double numerical evaluations; the shell-response inequality is
exact.

This finite witness does not disprove possible additional rigidity in the
asymptotic contracted regime `t->0`.  It does prove, inside the actual
finite minimizer class, that the complete signed energy histogram, cap,
pressure, and the Hamming-sphere scalar bounds do not determine the shell
response.

The computation is reproducible in
[`experiments/actual_child_energy_shell_response_falsifier.py`](../experiments/actual_child_energy_shell_response_falsifier.py),
with the recorded output in
[`computations/results/actual_child_energy_shell_response_falsifier.json`](../../computations/results/actual_child_energy_shell_response_falsifier.json).

## 4. The exact extra observable

Let `k_t(B|Q)` be the forward binary-channel density and define the
shell-output components

```math
K_e(B)=\sum_{Q:E(Q)=e}\mu(Q)k_t(B\mid Q),
\qquad p(B)=\sum_eK_e(B).                             \tag{SC.9}
```

The scalar shell histogram records only

```math
E_UK_e(B)=p_e.                                        \tag{SC.10}
```

For `q_a proportional p^aU`, however, the exact averaged posterior shell
law is

```math
\boxed{
\bar p_a(e)
={E_U[K_e(B)p(B)^{a-1}]\over E_Up(B)^a}.}             \tag{SC.11}
```

Thus the missing object is the vector of **shell-resolved inverse output
moments**

```math
R_a(e)={\bar p_a(e)\over p_e}
={E_{B\sim P_e}p(B)^{a-1}\over E_Up(B)^a},
\qquad {dP_e\over dU}={K_e\over p_e}.                \tag{SC.12}
```

For fixed `a`, this is only `O(N^2)` scalar values, but it contains channel
geometry absent from the child energy histogram.  Equivalently, it is a
shell-resolved overlap/likelihood profile.  Formula (SC.11) shows exactly why
the actual finite pair separates.

A sufficient optimizer-specific theorem would be

```math
\max_e\log R_a(e)=o(N)                               \tag{SC.13}
```

uniformly on the required negative path; then
`D(bar p_a||p)<=max_e log R_a(e)=o(N)`.  The scalable generic construction
shows that (SC.13) cannot follow from shell count, shell uniformity, or atom
spread alone.

## 5. Narrowed conclusion

No optimizer-specific sublinear shell-KL bound follows from the presently
available scalar pressure/cap/Hamming-sphere data.  The exact extra
obligation is not another energy moment: it is negative association between
rare combined-energy shells and inverse output likelihood, quantitatively
captured by (SC.12)--(SC.13).  Proving that response balance for actual
contracted-temperature minimizers would close the shell branch of Theorem
37.62.  Without it, the generic product example permits linear shell KL and
the actual finite pair proves scalar non-determinacy.
