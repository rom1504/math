# Positive residual variance is essential to the threshold passage

Date: 2026-09-06. Exact hypothesis check and actual-signing finite stress
test. The nonvanishing asymptotic hard-feedback energy suggested by the
test is NOT proved here.

## 1. A fixed regular response invisible to Gaussian L2

Let A_N be the exact binary-Steiner signing and delete its last vertex.
The resulting order n=N-1 is odd. Set B=A/sqrt(n-1). Principal compression
preserves a uniform operator cap, with the harmless normalization factor
sqrt((N-1)/(N-2)). Thus this is an actual symmetric hollow dense signing
sequence in the theorem's bounded-op regime.

For Boolean S define the genuine old fields G=BS, D=S h_2(G), Y=BD.
Take the FIXED functions

```math
f(g,y)={\bf1}_{g=0}\operatorname{sign}(y),\qquad
H(g,y)=1-|f(g,y)|.
```

Use sign(0)=0. Then f is bounded and odd under simultaneous reversal,
H is bounded, nonnegative and even, and |f|+H=1. Both are continuous
off the Gaussian-null line g=0. Therefore both meet the old Gaussian-
a.e.-continuity hypothesis. Nevertheless f is exactly zero Gaussian-a.e.,
so b_0=b_1=0 and its entire Gaussian residual variance tau^2 is ZERO.

On the actual discrete cube, each row of A has an even number n-1 of
signs. Consequently

```math
\Pr\{G_i=0\}=2^{-(n-1)}\binom{n-1}{(n-1)/2}
=O(n^{-1/2}).
```

It follows rigorously that E||f(G,Y)||^2/n=O(n^-1/2) and, by the
operator cap, E||Bf(G,Y)||^2/n=O(n^-1/2). For any fixed Lipschitz odd
bounded psi, the resulting feedback has normalized second moment and
absolute expected energy tending to zero. In particular the smooth
theorem is fully consistent with this example.

For the discontinuous response sign, that conclusion does not follow.
An arbitrarily small nonzero returned coordinate is still sent to a
unit sign. The ideal Gaussian comparator has variance zero and is
identically zero, so neither Gaussian small-ball control nor the inverse-
variance lemma applies. This explains the required tau>0 condition in
the ordered hard-threshold corollary.

## 2. Actual feasible hard feedback: reproducible evidence

The script
`computations/continued_audit_zero_variance_threshold_probe_2026_09_06.py`
uses exact integer row sums to recognize G_i=0. It forms

```math
C_i=H(G_i,Y_i)\operatorname{sign}((Bf)_i)
```

and checks both feasible endpoints f+C and -f+C samplewise. There is
no approximate zero test on G and no Gaussian substitution for the
original seeds. At n=27,119,495, with 6000 independent samples each:

```text
n                          27          119         495
E ||f||²/n                 .152123     .073395     .035837
hard-feedback energy       .220000     .242237     .235759
hard-energy standard error .001420     .000606     .000302
tanh-feedback energy       .044434     .024747     .012450
E ||Bf||²/n                .180260     .076554     .036407
```

The hard-feedback energy remains substantial in these finite tests while
the raw variance and smooth feedback decay. This is evidence against
silently extending the threshold identity to tau=0, but it is not an
asymptotic counterexample theorem. Establishing the limiting hard energy
would require a sparse local-limit/conditional-covariance analysis beyond
the fixed-Gaussian-L2 response theorem.

Reproduction:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_audit_zero_variance_threshold_probe_2026_09_06.py --samples 6000
```
