# Leading-order slack in the global coherent-state upper bound

Status: exact adversarial check. The finite-D truncation theorem remains
valid, but its coherent supremum is NOT asymptotically exact in general.

Let R be a real symmetric orthogonal map and v a real feature vector. Put
`c=||v||²` and `b=v^T R v>=0`. Give every vertex the same pure tensor
`T_i=v^(tensor(n-1))`. The exact complete-graph contraction is `b^(n(n-1)/2)`.

Decompose v=v_++v_- into the two eigenspaces of R. A complex unit vector of
the form u=a+i d, with a in the positive eigenspace and d in the negative
eigenspace, satisfies `u^T R u=1`. Choosing u=v_+/||v_+|| (when b>=0)
already gives `|<u,v>|²=||v_+||²=(c+b)/2`. Thus the coherent supremum in the
global theorem is at least `[(c+b)/2]^(n(n-1)/2)`, despite the exact value
being b to that power. The frame's polynomial prefactor cannot close this
exponential gap when c>b. No claim is made that this simple choice is optimal.

For an explicit positive Gaussian kernel example take q=2,
R=diag(1,-1), x=(r,r)/sqrt(2), and v=phi_L(x) with L>=1. Then

`b=exp(-2t r²)` and
`c=exp(-2t r²) sum_(j=0..L) (2t r²)^(2j)/(2j)! > b`.

The equality for b follows because x^T R x=0; all positive tensor degrees
vanish in the bilinear evaluation. Every source array is deterministically
x in all slots, exchangeable, and satisfies the energy condition with
C=r². Taking B>r deletes no edges. Therefore this is a counterexample
inside the exact hypotheses, not a tail pathology. For L=1 and tr²=1,
c=3 exp(-2), so the simple coherent choice already exceeds the actual
contraction by `2^(n(n-1)/2)`.

Consequently the new global reduction is a legitimate seed-retaining UPPER
bound framework, but it is not an organizing variational equality for the
ensemble. Any useful cap result needs an independently controlled coherent
supremum and must tolerate this leading-order cancellation loss. Positivity
of the original kernels does not remove the complex cancellation introduced
by coherent-state resolution.
