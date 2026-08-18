# Independent audit of the actual-child energy-shell decomposition

Status: **exact identities and constants passed, with an operational-scope
caveat**.  The shell marginal is a genuine polynomial-cardinality quotient
of the latent law, and the KL split is exact for the negative-path posterior
average.  The quotient is not dynamically closed: obtaining its retuned
marginal or its within-shell entropy deficits may still require the full
bridge geometry.  Two sharp abstract examples below show that either term
can remain linear while all displayed static spread constraints hold.

## 1. Signed-word fibres, uniformity, and shell count

A signed rank-one word `Q` has exactly two factorizations,
`Q=xy^T=(-x)(-y)^T`.  Hence it determines `[x]` and `[y]`, and

```math
E_\epsilon(Q)=H_A(x)+\epsilon H_D(y)                 \tag{VES.1}
```

is well-defined because both Hamiltonians are invariant under global spin
flip.  The four augmented preimages used before conditioning on `Q` have
already been summed in PA.11; no additional factor belongs in ES.2.
Consequently every signed word in a fixed shell has the same mass
`2^(2-m-n) cosh(te)/mathcal Z_epsilon`, proving exact shell uniformity.

If `d_m=binom(m,2)` and `d_n=binom(n,2)`, then
`H_A(x)+epsilon H_D(y)` lies between `-(d_m+d_n)` and `d_m+d_n` and has one
fixed parity.  It therefore assumes at most `d_m+d_n+1` values.  More
explicitly, if `g_A(a)` and `g_D(b)` are the full-cube signed energy
histograms, then

```math
|\mathcal S_e|
={1\over2}\sum_{a+\epsilon b=e}g_A(a)g_D(b).          \tag{VES.2}
```

The factor `1/2` is exactly the simultaneous-sign fibre.  Thus the prior
shell law really is determined by two `O(N^2)` histograms and a convolution.

## 2. Surprise ceiling

The sector normalizer satisfies

```math
\mathcal Z_\epsilon
\le (Z_A^++Z_A^-)(Z_D^++Z_D^-)
=4Z_A(t)Z_D(t)
\le4(\cosh t)^{d_m+d_n}.                             \tag{VES.3}
```

Using ES.2, `cosh(tE)>=1`, `t=beta/sqrt(N)`, and
`d_m+d_n<=N^2/2` gives

```math
\mu_\epsilon(Q)
\ge2^{-N}(\cosh t)^{-(d_m+d_n)}
\ge\exp\{-(\log2+\beta^2/4)N\}.                     \tag{VES.4}
```

Every nonempty shell contains a word, so the same lower bound holds for
`p_e`.  The constant `C_beta=log2+beta^2/4` in ES.5c is correct.

If `D(bar p||p)>=cN`, put `X=-log p_E`.  Then
`D=E_(bar p)X-H(bar p)<=E_(bar p)X`, while (VES.4) gives
`0<=X<=C_beta N`.  Splitting at `cN/2` yields

```math
\bar p\{X\ge cN/2\}\ge {c\over2C_\beta-c}.          \tag{VES.5}
```

The union bound over at most `O(N^2)` shells gives the second half of ES.9c.
The constants in ES.9c--ES.9d and Theorem 37.62 are therefore correct.

## 3. KL chain rules on the actual negative path

For any `bar mu<<mu_epsilon`, the ordinary chain rule under the deterministic
map `Q->E_epsilon(Q)` is

```math
D(\bar\mu\Vert\mu_\epsilon)
=D(\bar p\Vert p)
 +\sum_e\bar p_eD(\bar\mu_e\Vert U_{\mathcal S_e}),  \tag{VES.6}
```

and the conditional term equals `log|S_e|-H(bar mu_e)`.

For the bridge posterior, its density relative to the augmented prior is a
function of `Q` alone.  This remains true after averaging against any bridge
law, including every actual negative escort `q_a`.  Thus both the pointwise
and averaged posterior conditional laws on the fibre `z->Q` equal the prior
conditional, and the conditional KL vanishes:

```math
D(\bar\nu\Vert\nu)=D(\bar\mu_Q\Vert\mu_\epsilon).    \tag{VES.7}
```

Applying the same chain rule to the joint law of `(B,E)` gives ES.10a;
`I(B;E)<=H(E)<=log O(N^2)=O(log N)`.  Hence ES.6, ES.10, and ES.10a all
pass without a hidden preimage or fibre charge.

## 4. Is this a strict low-information branch?

Semantically, yes, for the scalar retuning KL: the two shell marginals and
one conditional entropy per shell are only `O(N^2)` real numbers and do not
determine the conditional tables.  In particular, many exponentially
different within-shell laws have the same entries in ES.6.

Operationally, no closure theorem has yet been proved.  The prior marginal
`p` is determined by the child histograms, but the retuned marginal `bar p`
and the conditional entropies are produced by the bridge likelihood.  They
are not currently computable from the two child histograms alone.  Thus
Theorem 37.62 is a strict *output quotient*, not yet a reusable
low-information evolution law.  This distinction is consistent with the
source's explicit statement that neither row lifetime nor target reach is
controlled.

The finite bridge collision in
[`../../artifacts/phase2j_augmented_cut_gram_response_audit.md`](../../artifacts/phase2j_augmented_cut_gram_response_audit.md)
is a useful warning: with the positive order-three children, bridge codes
`78` and `85` have identical conditional second moments in every internal
energy-shell pair but exact parent caps `11` and `9`.  This is finite rather
than scalable and does not falsify ES.6, but it proves that even substantial
shell-conditioned data need not close the geometric response.

## 5. Sharpness: either branch can be linear

The static conclusions of Theorems 37.61--37.62 alone cannot make either
term sublinear.  Here are two exact information-theoretic witnesses.

For the radial branch, take two shells of `exp(aN+o(N))` atoms each, put

```math
p_0=e^{-cN},\qquad p_1=1-e^{-cN},
```

and make the prior uniform inside each shell.  Let `bar mu` be uniform on
the first shell.  Then

```math
D(\bar p\Vert p)=cN,\qquad
\sum_e\bar p_eD(\bar\mu_e\Vert U_{\mathcal S_e})=0.  \tag{VES.8}
```

Choose `0<a<log 2`, `0<c`, and `a+c<C_beta`.  There are then enough
rank-one labels for the abstract atoms, and every prior atom is at least
`exp(-C_beta N)` and at most `exp(-aN+o(N))`, so this example simultaneously
has a linear min-entropy rate and satisfies the surprise ceiling.

For the geometric branch, take one shell of `exp(aN+o(N))` atoms with its
uniform prior, and let `bar mu` be uniform on a subset of
`exp((a-d)N+o(N))` atoms, where `0<d<a`.  Then

```math
D(\bar p\Vert p)=0,\qquad
D(\bar\mu\Vert U_{\mathcal S})=dN+o(N),              \tag{VES.9}
```

while the posterior itself still has an exponential effective support.

These witnesses are abstract shell-uniform laws, not asserted to arise from
actual quadratic children or from their negative escort.  Their precise
role is to show that shell count, atom diffuseness, and the KL chain rule by
themselves cannot close either branch.  An actual-child theorem must use a
new optimizer or bridge identity.

## 6. Verdict

The combined-energy-shell theorem is mathematically correct and is a strict
compression of the *value* of latent retuning KL.  It narrows the SML to a
real dichotomy, but does not solve either side and is not a closed child
state.  The sharp remaining alternatives are exactly linear coherent
retuning onto exponentially rare energy shells or linear entropy loss spread
diffusely inside a large shell.
