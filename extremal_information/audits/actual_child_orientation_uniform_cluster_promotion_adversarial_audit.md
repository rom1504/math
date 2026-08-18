# Adversarial audit: orientation-uniform actual-child cluster promotion

## Verdict

**OU.1--OU.2 PASS, with one wording correction.**  A centrally symmetric
row law with bounded Renyi-two divergence from the fair cube is uniformly
vector-subgaussian.  The canonical erased-row inverse escort has exactly
this symmetry and a dimension-free Renyi-two bound in every relative child
orientation and in either transpose direction.  Substituting that bound in
the random-row-cut proof of SP.2 makes the physical sector--Gram/cluster
promotion theorem orientation-uniform.

The result does **not** prove that the bias-canceling orientation of
Theorem 37.32 is target-optimal or within `o(N)` of the target-optimal
orientation.  It makes that comparison unnecessary: the same theorem can
be applied directly in whichever actual orientation reaches the target.
Thus the accurate claim is that the *target-orientation obligation is
removed*, not that the balanced orientation has been proved target-relevant.

## 1. Audit of the Renyi-two subgaussian lemma

Let `f=dP/dU_n`, put `X=<v,R>` and `a=||v||_2^2`, and suppose
`D_2(P||U_n)<=C`.  Central symmetry gives

```math
 E_Pe^X=E_P\cosh X.
```

Since `E_Uf^2<=e^C`, Cauchy--Schwarz gives

```math
E_P\cosh X
\le 1+e^{C/2}\{E_U(\cosh X-1)^2\}^{1/2}.             \tag{A.1}
```

For independent fair signs,

```math
 E_U\cosh(2X)=\prod_j\cosh(2v_j)\le e^{2a},
 \qquad E_U\cosh X\ge1+{a\over2}.
```

Using `cosh^2X=(1+cosh(2X))/2` therefore yields

```math
 E_U(\cosh X-1)^2
 \le {e^{2a}-1-2a\over2}
 \le a^2e^{2a}.                                      \tag{A.2}
```

The last step is the standard exponential Taylor remainder at `2a`.
Also

```math
 ae^a\le e^{2a}-1
```

because `e^{2a}-1-ae^a=e^a(2sinh a-a)>=0`.  With
`K=e^(C/2)>=1`, convexity of `y mapsto y^K` at `y=1` gives

```math
1+K(e^{2a}-1)\le e^{2Ka}.
```

Combining these inequalities proves

```math
\boxed{E_Pe^{\langle v,R\rangle}
       \le e^{2e^{C/2}\|v\|_2^2}.}                 \tag{A.3}
```

In the convention `Ee^(<v,R>)<=exp(sigma^2||v||^2/2)`, the proxy is

```math
\boxed{\sigma_C^2=4e^{C/2}.}                        \tag{A.4}
```

All constants in OU.1 are valid.  Central symmetry is essential: a
bounded-`D_2` law may otherwise have nonzero mean, precluding a centered
quadratic MGF bound of this form.

## 2. The canonical row satisfies the hypotheses in every orientation

For fixed orientation `epsilon`, the erased-row output likelihood is

```math
p_{\rm row}(b)
=E_Y\prod_j(1+\tanh(u)b_jY_j).
```

Every sector-biased child spin law is invariant under `Y mapsto -Y`, since
the child Hamiltonian is quadratic.  Hence `p_row(-b)=p_row(b)`, and its
inverse escort `r_epsilon proportional p_row^(-lambda)U_n` is centrally
symmetric.  This remains true for arbitrary scalar sector bias.

Flipping one row bit changes `log p_row` by at most `2u`.  Applying the
conditional escort cube lemma gives exactly CR.8,

```math
D_2(r_\epsilon\Vert U_n)
\le n\log\{1+\tanh^2(\lambda u)\}
\le\lambda^2u^2n.                                    \tag{A.5}
```

At `u=beta/sqrt(N)`, one may use

```math
C_2=\lambda^2u^2n\le\lambda^2\beta^2,
\qquad \sigma_*^2=4e^{\lambda^2\beta^2/2}.          \tag{A.6}
```

No optimizer envelope, sector-bias bound, choice of orientation, or choice
of transpose direction enters (A.5).  The row product in CR.3 is independent
by definition, so the remaining hypothesis of the row-cut argument is also
satisfied.

## 3. Check of the quadratic-chaos constants

The proof of SP.2 needs only the linear-functional MGF bound and independence
between row blocks.  If that bound has proxy `sigma^2`, its two-shore
Gaussian determinant calculation gives

```math
\log Ee^{\theta H}\le4\sigma^4\theta^2V
\quad\hbox{when}\quad
|\theta|\|M\|_{op}\le {1\over2\sqrt2\,\sigma^2}.    \tag{A.7}
```

Substitution of (A.6) gives the conservative valid constants

```math
a_*={1\over8\sqrt2e^{\lambda^2\beta^2/2}},
\qquad b_*=64e^{\lambda^2\beta^2}.                   \tag{A.8}
```

These are not the constants from the balanced `D_infinity` theorem.  Its
cutoff is `1/(2sqrt(2)e^(lambda(beta^2/2+log2)))`; neither cutoff dominates
the other for every `beta,lambda`, and the Renyi-two cutoff is much smaller
at many large target parameters.  Orientation-uniformity therefore trades
away the target comparison at the price of a different, sometimes more
restrictive, fixed normalized Gram threshold.  It does not transport every
instance covered by SP.3 verbatim.

For the sector--Gram chaos, `V=K_epsilon` and
`||M||op<=sqrt(2K_epsilon)`.  Taking
`theta=-lambda beta^2/N` turns the determinant condition into

```math
\lambda\beta^2\sqrt{2\kappa}\le a_*
\quad\text{when}\quad K_\epsilon\le\kappa N^2.      \tag{A.9}
```

Independent symmetric rows have zero mean, so `E H_2=0`; no centering term
was omitted.  Combining (A.7) with the already audited connected-cluster
remainder `osc R<=2mathfrak C_(>=4)^epsilon` proves OU.10 and OU.11 with the
displayed constants.

## 4. Exact orientation/target identities

The distinction between balanced orientation and target orientation can be
made completely explicit.  Put

```math
T=\log\overline Z_A(t)+\log\overline Z_D(t),
\quad d=mn,
\quad a=\tanh\gamma_A,
\quad b=\tanh\gamma_D.
```

Let `p_epsilon=dPi_epsilon/dU_B` be the forward binary-channel output
density from the exact zero-bridge child prior in orientation `epsilon`.
The sector partition and channel factorizations give, for every bridge,

```math
\boxed{
L_\epsilon(B)
=T+\log(1+\epsilon ab)+d\log\cosh t+\log p_\epsilon(B).} \tag{A.10}
```

Consequently the exact negative-disorder value is

```math
\boxed{
V_\epsilon-T
=d\log\cosh t+\log(1+\epsilon ab)
-D_{1+\lambda}(U_B\Vert\Pi_\epsilon).}             \tag{A.11}
```

Thus the zero-bridge advantage of the bias-canceling orientation can be
reversed after bridge optimization precisely when the opposite orientation
has a larger reverse-Renyi output response by more than

```math
\log{1+|ab|\over1-|ab|}=2\operatorname {artanh}|ab|. \tag{A.12}
```

The response in (A.11) uses the full bridge output law, so (A.11) is an
exact diagnostic identity, not an admissible low-information solution of
the balanced-product SML.

There is nevertheless no target loss in applying OU.2.  If

```math
V_{\rm mix}
=-{1\over\lambda}\log\left{
 {e^{-\lambda V_+}+e^{-\lambda V_-}\over2}\right},
```

then

```math
\boxed{
\min_\epsilon V_\epsilon
\le V_{\rm mix}
\le\min_\epsilon V_\epsilon+{\log2\over\lambda}.}  \tag{A.13}
```

Hence joint-orientation target reach implies that an actual deterministic
orientation reaches at least as well, while a preselected optimal
orientation trivially remains target-reaching.  OU.2 applies in that same
orientation.  Transposition is a measure-preserving reindexing of bridges,
and OU.2 is valid in either row direction, so it introduces no hidden
transpose loss either.

## 5. Sharp finite scope check

A complete actual-child enumeration shows why OU.2 should not be described
as a proof that the balanced orientation itself is optimal.  At

```text
N=10, split 3+7, beta=4, t=4/sqrt(10),
```

the order-three child is the unique signed-permutation/global-sign class,
and exhaustive enumeration of all `32,768` rooted-gauge order-seven
signings selects one minimizing class with a `100`-digit pressure gap
`0.1726269103...` to the next absolute-energy histogram.  For the
representatives with hashes

```text
f71efe16fa6b412d486f6d4f29d78174c49117aed7862ace6cebb9c67e762995
93191291e27943d7bda3c5fdb4a4fee10a8b95fa2d4d34dec01e90ec5d7c559e
```

their sector biases are approximately `.7239803` and `-.5840622`, so
`epsilon=+1` is bias-canceling.  Complete enumeration of all `2^21`
bridges gives

```math
\begin{array}{c|cc}
 &\epsilon=+1\text{ (balanced)}&\epsilon=-1\\ \hline
\min_B L_\epsilon(B)&13.8749219632&13.3492885313\\
V_{\epsilon,\lambda=2}&16.9491125147&16.9200683898\\
V_{\epsilon,\lambda=4}&16.1329693311&15.8112524421.
\end{array}                                          \tag{A.14}
```

The convolution was checked against direct spin sums and the largest
reported log-pressure discrepancy was below `3e-15`.  This is a robust
finite actual-minimizer counterexample to pointwise or exact
orientation-optimality.  It is **not** an asymptotic counterexample to an
`o(N)` orientation-loss theorem, and neither orientation reaches the child
target in this example.  Its correct role is to rule out deriving target
relevance from the zero-bridge identity alone.

## 6. Frontier consequence

Within the explicit OU.9 Gram regime, the target-orientation clause can be
deleted from the physical-promotion SML.  In the orientation which actually
reaches the target, one must still decide from an operationally smaller
actual-child observable whether

```math
J=o(N),\qquad I^\leftarrow=\Omega(N),\qquad
\text{or}\qquad J-I^\leftarrow=\Omega(N).
```

OU.2 does not decide this trichotomy, prove cluster tightness, or show that
the all-order absolute cluster scalar is operationally simpler than the
complete child Gibbs law.  The correct classification is therefore a
**target-relevance RESET at Level 5**, not a Level-6 recurrence and not a
solution of `L_balanced-product-phase` itself.
