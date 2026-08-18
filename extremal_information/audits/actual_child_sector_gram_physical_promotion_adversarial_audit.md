# Adversarial audit of the sector--Gram physical-promotion draft

## Verdict

**PASS WITH REQUIRED CORRECTIONS.**  Lemmas SP.1--SP.2 and the inequality in
SP.3 are mathematically sound.  In particular, bounded density relative to a
row cube, together with central symmetry, really does give the dimension-free
joint subgaussian estimate needed by the off-block Hanson--Wright argument;
coordinate independence inside a row is not required.  The physical
normalization is also correct:

```math
\theta=-\lambda t^2,
\qquad
t={\beta\over\sqrt N},
\qquad
\theta^2K_\epsilon
={\lambda^2\beta^4\over N^2}K_\epsilon.
```

There are nevertheless four required corrections or qualifications:

1. the cumulant-series hypothesis must quantify over **every bridge word**;
2. the proof of SP.2 should use the exact random-cut Jensen argument below,
   rather than the phrase "a constant multiple";
3. SB.2--SB.3 remove the bounded-sector-bias hypothesis, after making an
   existential balanced orientation and transpose choice;
4. the cluster tail is a child-only scalar query, but it has not been shown
   to be a lower-information or easier obligation than the physical
   interaction cumulant.  It therefore does not by itself justify a strict
   SML reset or a Level-6 claim.

Two unambiguous source typos were repaired together with this audit: the
corrupted parity subscript in (SP.2), and the fragment "Apply The" in the
proof of SQ.4.

## 1. Exact audit of SP.1

For a fixed bridge word `B`, the binary-channel identity is

```math
p_u(B)
=E_{\nu_\epsilon}\exp\left(u\sum_eB_eQ_e\right)
 (\cosh u)^{-mn}.
```

The row marginal has the identical formula restricted to the edges of that
row, with `(cosh u)^(-n)`.  Hence all channel-normalization constants cancel
in

```math
h_u=\log p_u-\sum_i\log p_{i,u}.
```

Subtracting the row series removes exactly the cumulant tuples whose left
endpoints all agree.  The transformation `X -> -X` preserves the zero-bridge
prior and sends every `Q_e` to `-Q_e`, so every odd joint cumulant vanishes.
For order two, the ordered tuples with distinct rows occur in the two orders;
the factor `1/2!` therefore leaves exactly

```math
u^2\sum_{i<k}\sum_{j,\ell}
\Gamma_{ik;j\ell}^{(\epsilon)}B_{ij}B_{k\ell}
=u^2H_2(B).
```

Thus there is no missing factor of two in (SP.3).  Repeated-edge tuples may
produce constants because `B_e^2=1`; placing those terms in `c_u` only
improves the oscillation bound.  The triangle inequality gives
`osc(R_u)<=2 C_(>=4)(u)` with the stated ordered-tuple normalization.

The hypothesis should say explicitly that, for **each** `B` and every real
`v` between `0` and `u`, the Taylor cumulant series at zero converges to

```math
\log E_{\nu_\epsilon}\exp\left(v\sum_eB_eQ_e\right).
```

Absolute finiteness of the displayed aggregate at `u` controls the
coefficient sums, but the equality with the log MGF (including at the
endpoint) should remain an explicit hypothesis.  At finite order the MGF is
positive on the real axis, yet complex zeros can make its zero-field Taylor
radius much smaller than the physical real amplitude, so this condition is
not automatic.

## 2. SP.2 is valid, with an explicit proof

Put `K_0=e^C` and let `f=dP/dU_n`.  Central symmetry first gives, for every
`v`,

```math
\begin{aligned}
E_Pe^{\langle v,R\rangle}
&=E_P\cosh\langle v,R\rangle\\
&=1+E_U f\{\cosh\langle v,R\rangle-1\}\\
&\le1+K_0\{E_U\cosh\langle v,R\rangle-1\}\\
&\le1+K_0\{e^{\|v\|^2/2}-1\}
 \le e^{K_0\|v\|^2/2}.
\end{aligned}                                      \tag{A.SP.1}
```

The last step is `(1+x)^(K_0)>=1+K_0x`.  Thus (SP.8) holds with
`sigma_C^2=K_0`, uniformly in the row dimension.  Independence of the row
blocks gives the same linear-functional estimate for any concatenation.

For completeness, let `delta_i` be iid fair bits and let `H_delta` retain
only the pairs cut by the random bipartition.  Since
`E_delta H_delta=H/2`, convexity gives the exact decoupling

```math
e^{\theta H}
\le E_\delta e^{2\theta H_\delta}.                 \tag{A.SP.2}
```

For a fixed cut, write `H_delta=R_S^TBR_T`.  Conditioning successively on
the two independent shores and using (A.SP.1) plus a standard Gaussian
linearization yields

```math
E e^{2\theta R_S^TBR_T}
\le
\det(I-4\sigma_C^4\theta^2B^TB)^{-1/2}.            \tag{A.SP.3}
```

Here `||B||op<=||M||op` and `||B||F^2<=V`.  If

```math
|\theta|\,\|M\|_{op}
\le {1\over2\sqrt2\,\sigma_C^2},                 \tag{A.SP.4}
```

then `-log(1-x)<=2x` in (A.SP.3) gives

```math
\log E e^{\theta H}
\le4\sigma_C^4\theta^2V.                          \tag{A.SP.5}
```

Thus one may take, for example,
`a_C=1/(2 sqrt(2)e^C)` and `b_C=4e^(2C)`.  This also verifies that there is
no dimension-dependent prefactor hidden in SP.2.

For the canonical row inverse law, central symmetry needed here follows
from `p_row(-b)=p_row(b)`, obtained by the global flip `Y -> -Y`.  The rows
are independent by definition of `r_u`.  Therefore SP.2 applies to `H_2`,
with `V=K_epsilon` and no missing diagonal term.

## 3. Audit and sharpened statement of SP.3

The centered identity

```math
\mathcal J_t
=\log E_{r_t}e^{-\lambda(h_t-E_{r_t}h_t)}
```

combined with SP.1 implies pointwise

```math
-\lambda(R_t-E_{r_t}R_t)
\le\lambda\operatorname {osc}R_t.
```

Since `E_(r_t)H_2=0`, SP.2 at `theta=-lambda t^2` gives exactly

```math
\mathcal J_t
\le b_C\lambda^2t^4K_\epsilon
   +2\lambda\mathfrak C_{\ge4}(t).                 \tag{A.SP.6}
```

Consequently, if `K_epsilon<=kappa N^2` and
`lambda beta^2 sqrt(2 kappa)<=a_C`, its contribution is at most the constant
`b_C lambda^2 beta^4 kappa`.  This explicit formulation is preferable to
`K=O(N^2) (with the constant satisfying ...)`, because the smallness of the
hidden constant is essential.  Rearranging (A.SP.6) verifies (SP.14), with
the advertised factor `1/(2 lambda)` and an `O(1)` remainder.

The bounded-bias premise in the draft is unnecessary.  Given exact children
`C,D`, relabel them so that `|gamma_C|<=|gamma_D|`, use `D` as the erased-row
base, and choose the relative orientation which cancels their bias signs.
SB.3 gives

```math
D_\infty(r_{C\to D}^\epsilon\Vert U_n)
\le\lambda\{\delta_n(t)+\log2\}
\le\lambda\{\beta^2/2+\log2\}.                    \tag{A.SP.7}
```

Thus SP.3 remains valid for this balanced presentation with

```math
C=\lambda(\beta^2/2+\log2),
```

uniformly in both child biases.  The resulting `K_epsilon` and cluster tail
must, of course, be evaluated in that same chosen orientation.  This is an
existential orientation-and-transpose theorem.  Although that orientation
minimizes the zero-bridge augmented pressure, no current theorem proves that
it remains target-reaching after the bridge is optimized; any convergence
application must still price this orientation choice.

## 4. Information and frontier judgment

The theorem genuinely replaces a parent bridge/external-field table by a
condition on the actual zero-bridge child law.  Its input is a constant-size
sector--Gram carrier plus one nonnegative scalar functional.  That makes it
a useful **named child-only closure criterion** and not a conference or
generic-row surrogate.

It is not yet demonstrably a lower-information reduction in the stronger
operational sense used by the frontier.  The scalar
`mathfrak C_(>=4)(t)` sums absolute joint cumulants of every order and every
cross-row edge tuple.  Certifying it can require the complete high-order
child Gibbs law; an exact real scalar has no finite-bit complexity bound.
Moreover, absolute values discard cancellations, so the sufficient condition
`mathfrak C_(>=4)=o(N)` can be strictly stronger than the needed conclusion
`mathcal J=o(N)`.  A large absolute tail does not by itself prove physical
gain, row total correlation, or reverse-product information.

Accordingly, the rigorous conclusion is:

```text
balanced row support + small quadratic tangent + sublinear absolute
connected cluster tail  =>  sublinear canonical interaction cumulant;

linear canonical interaction cumulant + small quadratic tangent
=> linear absolute connected cluster tail.
```

This is a valid physical-scale conditional theorem, but selecting which
branch actual optimizing children occupy remains open.  It narrows the
mathematical form of one sufficient condition; it does not yet establish
that the remaining condition is simpler than the balanced product-phase
problem or create a credible Level-5-to-6 route on its own.
