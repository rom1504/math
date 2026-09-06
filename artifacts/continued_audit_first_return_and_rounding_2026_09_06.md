# Independent first-return, covariance, and rounding audit

Date: 2026-09-06. This audit reconstructs the arguments, rather than importing
the positive verdicts in previous audits. It concerns the continued campaign;
it does not prove convergence or change the banked universal coefficient.

## 1. Normalization and independently checked mechanism

Let `A` be symmetric, hollow, and sign-valued off the diagonal, and put
`B=A/sqrt(n-1)`, `Q=B^2`, and `Lambda(B)=max_x |x^T Bx|/(2n)`.
Exactly, `Q_ii=1`, `|Q_ij|<=1`, and, under `||B||op<=L`,
`Tr(Q^2)<=L^2 n`. The original normalization is
`Lambda(B)=Qabs(A)/(n sqrt(n-1))`, not `Qabs(A)/n^(3/2)` at finite n.

The full nonlinear covariance proof has a sound and essential norm choice.
After local Hermite/Wick replacement, a partial matching of two old forests
has either a non-star component containing a four-vertex path, two nontrivial
stars, or one star accompanied by a whole-branch edge. Two disjoint proper
flattening gains give `O(1/n)` entries in the first two cases. In the last
case, an `O(n^-1/2)` entry factor times one old covariance matrix of
Frobenius norm `O(sqrt(n))` gives Frobenius norm `O(1)` as well. Hence every
partial-pairing error has nuclear norm `O(sqrt(n))=o(n)`.

This does NOT prove vanishing operator norm. The omitted distinction is
material: a matrix all of whose entries are `1/n` has operator norm one.
Whole-branch matches yield the Schur main
`sum_(k odd,k>=3) w_k Q^(circ k)` with the local normalized Hermite factorials
canceling exactly. Local replacement errors pass to normalized nuclear norm
through `||E[UV^T]||_* <= sqrt(E||U||^2 E||V||^2)`.

For old odd `F`, the first/nonlinear energy cross term is controlled by the
endpoint Frobenius estimate, not by a nonexistent operator estimate for the
mixed covariance. The endpoint parity graph has forest-root degree at least
three, excluding the dangerous single surviving edge. This reconstructs

```math
\frac{\mathbb E F^{\mathsf T}BF}{2n}
=\frac{\|P_1F\|_2^2}{2}\frac{\operatorname{Tr}B^3}{n}+o(1).
```

The first-chaos part of the NEXT field is nevertheless `Q[S K_F]`.
Cubic self-energy cancellation is not a statement that `Q` is identity.

## 2. An actual-signing, bounded smooth first-return falsifier

Let `H_m` be the symmetric Sylvester matrix for `m=2^r`, `r>=1`, and set

```math
C=J_2\otimes H_m,\quad D=\operatorname{diag}C,\quad A=C-D,
\quad n=2m,\quad B=A/\sqrt{2m-1}.
```

This is a dense hollow signing, with no weighted or missing off-diagonal
entries. Write `P=J_2 tensor I_m` and
`B_0=(J_2/sqrt(2)) tensor (H_m/sqrt(m))`. Then

```math
\|B-B_0\|_{\rm op}=O(m^{-1/2}),\qquad B_0^2=P,
\qquad \|B\|_{\rm op}\longrightarrow\sqrt2.
```

The cubic trace vanishes EXACTLY, not merely asymptotically. Indeed
`Tr H_m=0`, `C^2=2m P`, `D^2=I`, and `Tr D=0`; cyclicity gives

```math
\operatorname{Tr}(C-D)^3
=\operatorname{Tr}C^3-3\operatorname{Tr}C^2D
 +3\operatorname{Tr}CD^2-\operatorname{Tr}D^3=0.
```

For independent Boolean `S`, the actual return obeys

```math
\frac1n\mathbb E\|B^2S-PS\|^2
=\frac1n\|B^2-P\|_F^2\longrightarrow0.
```

At twin coordinates, `(PS)_i=S_i+S_i'`. Thus its variance tends to two,
and its fourth moment tends to eight, versus one and one for `S_i`.
There is also a bounded, smooth test with a gap of one:

```math
\frac1n\sum_i\mathbb E\left[S_i
 \sin\left(\frac\pi2(B^2S)_i\right)\right]\longrightarrow0,
\qquad
\mathbb E\left[S_i\sin\left(\frac\pi2 S_i\right)\right]=1.
```

The first limit follows immediately by the Lipschitz bound for sine and
the normalized-L2 estimate, since `sin(pi(S_i+S_i')/2)=0` pointwise.
For finite checks its exact product formula is
`sin(pi Q_ii/2) product_(j!=i) cos(pi Q_ij/2)`.

The example disproves a nonlinear exceptional-return transfer under a
bounded operator cap even with exactly vanishing cubic trace. It does not
claim that these high-cap twin matrices are minimizers.

## 3. Why a one-root law does not specify a randomly masked covariance

This section is an information-theoretic Gaussian counterexample, not a
signing construction or a normalized-energy counterexample.

Take jointly Gaussian `(R1,R2,X1,X2)` with covariance

```math
\begin{pmatrix}
1&0&0&a\\0&1&a&0\\0&a&1&b\\a&0&b&1
\end{pmatrix},\qquad a=b=\tfrac12.
```

It is positive definite because its Schur complement is
`[[1-a^2,b],[b,1-a^2]]`. Each same-root pair `(Ri,Xi)` is independent
standard normal; `Cov(R)=I`, and `Cov(X)` is fixed. Compare with a second
model having the same two marginal covariance matrices but independent
vectors `R` and `X`. All these one-root data are identical.

For the bounded even mask `H(x)=(1+cos(tx))/2`, Gaussian integration by parts
gives, in the crossed model,

```math
\mathbb E[R_1H(X_1)R_2H(X_2)]
=\frac{a^2t^2}{4}e^{-t^2}\sinh(bt^2)>0.
```

It is zero in the independent model. Repeating the two-coordinate blocks
produces an order-n nuclear covariance discrepancy. Thus one-root joint
laws plus the separate old/new covariance matrices do not determine a
randomly masked covariance. A new cross-root contraction theorem, or a
different inequality genuinely avoiding it, is necessary.

An actual-signing parity warning is also exact: for `G=BS`,
`E[S_i G_i S_j G_j]=1/(n-1)` at `i!=j`, for EVERY hollow signing. Hence
the own-spin covariance-identity theorem cannot drop its even-response
hypothesis. No probabilistic matrix ensemble is needed for this warning.

## 4. Fourfold improvement of the finite slack theorem

The previously recorded finite correlated-rounding proof is correct.
Its profile step, however, loses a factor of two before squaring.
The sharp linear-in-slack profile inequality is

```math
\phi\!\left(\Phi^{-1}\!\left(\frac{1+u}{2}\right)\right)
\ge\frac{1-u^2}{\sqrt{2\pi}},\qquad -1\le u\le1. \tag{1}
```

Here is a calculus proof. Let the left side be `p(u)` and put
`g(u)=p(u)-phi(0)(1-u^2)`. By symmetry it suffices to use `[0,1]`.
The endpoint values are `g(0)=g(1)=0`, and `g'(0)=0`. On `(0,1)`,
with `t=Phi^-1((1+u)/2)`,

```math
g''(u)=2\phi(0)-\frac1{4\phi(t)}.
```

This is strictly decreasing, starts positive because `pi<4`, and tends to
minus infinity. Therefore `g'` first increases and then decreases, with
one positive-to-negative crossing; `g` increases and then decreases to its
zero endpoint. This proves (1). Its coefficient `phi(0)` is optimal since
equality holds at `u=0`.

Choose an orientation `B'=+B` or `-B`, and sample a Gaussian with correlation
matrix `I+B'/L`. This is positive semidefinite even at the spectral boundary.
Thresholding coordinate i at `t_i=Phi^-1((1+u_i)/2)` gives Boolean means
exactly `u_i`; at endpoints use infinite thresholds. The first threshold
Hermite coefficient is `-2 phi(t_i)`. All degrees at least two satisfy

```math
|R_{ij}|\le\frac{\sqrt{(1-u_i^2)(1-u_j^2)}}{L^2(n-1)}.
```

Writing `d=1-||u||^2/n`, the higher-chaos energy remainder is at most
`d/[2L^2 sqrt(n-1)]`. The exact first-chaos gain is

```math
\frac{2}{nL(n-1)}\left[
  \left(\sum_i\phi(t_i)\right)^2-\sum_i\phi(t_i)^2\right].
```

By (1), the first square is at least `n^2 d^2/(2pi)`. Consequently

```math
\Lambda(B)\ge |e_B(u)|+\frac{d^2}{\pi L}
-\frac1{\pi L(n-1)}-\frac1{2L^2\sqrt{n-1}}. \tag{2}
```

The gain coefficient is four times the prior recorded one. A useful
slack-sensitive finite variant follows from the elementary Hermite bound
`4 phi(t_i)^2<=1-u_i^2`:

```math
\Lambda(B)\ge |e_B(u)|+\frac{n d^2}{\pi L(n-1)}
-\frac{d}{2L(n-1)}-\frac{d}{2L^2\sqrt{n-1}}. \tag{3}
```

Thus for `|e_B(u)|>=Lambda(B)-epsilon`, letting

```math
c=\frac{n}{\pi L(n-1)},\qquad
a=\frac1{2L(n-1)}+\frac1{2L^2\sqrt{n-1}},
```

one has `d<=(a+sqrt(a^2+4c epsilon))/(2c)`. In particular the certified
error is now `O_L(sqrt(epsilon)+n^-1/2)`, rather than the cruder
`O_L(sqrt(epsilon)+n^-1/4)` from a slack-independent error.

There is no universal-bound improvement here: the banked endpoint is
Boolean. Damping a Boolean energy `e=.433322...` loses `ed`; even the
improved guaranteed gain `d^2/(pi L)<=d/pi` cannot pay this loss.

## 5. Operator deletion and countable-to-finite orders

The spectral-cap removal reconstructs independently. Cube rounding and
polarization give `beta(A)<=4 Qabs(A)`. The tensor/hyperplane proof of real
Grothendieck yields a nonnegative diagonal majorant `D>=+/-A` with
`Tr D<=K beta(A)`. Deleting entries above `K beta(A)/(eta n)` retains at
least `(1-eta)n` coordinates with normalized operator cap `L_eta<infinity`
for any low-cap signing sequence. Exact principal monotonicity holds by
averaging independent omitted spins.

Fix eta FIRST. For this fixed operator cap, fix the finite coordinates,
polynomials, smoothing widths, and Gaussian approximation accuracies needed
by the desired response. Send the retained matrix order to infinity. Only
then remove response approximation and smoothing errors, and finally send
eta to zero. The resulting loss is exactly `(1-eta)^(3/2)`.

No rate uniform in eta, growing response depth, or original operator cap is
being inferred. Finite conditional expectations preserve joint feasibility
`|F|+H<=1`; applying the same Gaussian Markov smoothing to both functions
also preserves it. The contraction of `U^-1 P1` gives continuity of the
response functional in Gaussian L2. This licenses countable-space closures
as limits of FIXED finite constructions, not as growing-depth matrix
algorithms.

## 6. Reproducible checks and scope

Run

```bash
.venv/bin/python computations/continued_audit_first_return_and_rounding_2026_09_06.py
```

The script verifies exact cubic traces for the actual twin signings,
integer/Fraction return errors, the bounded sine witness, the covariance
model, and a numerical profile grid. For n=512 its exact normalized squared
return errors are `1017/261121` to the twin return and `261116/261121` to
identity. The sine witness is approximately `.00305924`, versus one for
the involution law. These computations supplement the proofs above; the
grid is not the proof of the Gaussian-profile inequality.

No existing untracked artifact, global ledger, steering file, or stored
certificate was modified by this audit.
