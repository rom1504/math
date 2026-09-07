# Actual finite-rank sign surgery and optimizer-derived covariance balance

Date: 2026-09-07. Status: **Initial proof independently audited; subsequent
rank-one and interpolation improvements independently reconstructed**. This is a positive
same-order construction followed by a consequence of global optimality. It
does not prove a cross-order recurrence or assume an optimizer ramp.

The stronger implementation in
`flatify_adversary_2026_09_07_target_contraction_surgery.md` now supersedes
the masking/bias bound: arbitrary nuclear perturbation budget L gives cap
error O(n sqrt(L)+n), with no input cap assumption. The construction and
failed routes below remain preserved as the derivation record. The fixed
feature-space minimax consequences remain valid, now with the improved
O(n^(-1/4)) error for fixed rank and penalty.

## 1. Uniform implementation of a signed finite-rank perturbation

Let A be any hollow full signing of order n, Q(A)<=C n^(3/2). Let U be an
n-by-r matrix with orthonormal columns, let T be any real symmetric r-by-r
matrix of operator norm at most1, and put G=UTU^T. Fix theta>0. We construct
an actual full signing A' with

```math
 Q(A')\le Q(A-\theta\sqrt n\,\operatorname{offdiag}G)
                 +n^{3/2} e_n.                           (1)
```

The error is uniform over U and T, with no incoherence assumption. Choose
mu>0 such that theta sqrt(n) mu<=1. With a=(n+2)log2, one valid error is

```math
 e_n=\theta\sqrt{r/(\mu n)}+\frac{\theta r}{2n}
       +\theta C\mu\sqrt n
       +2\sqrt{\theta r a}\,n^{-3/4}
       +\frac{4a}{3n^{3/2}}.                             (2)
```

One may replace r/(mu n) by min(1,r/(mu n)). The displayed diagonal payment
accounts for masking. When later interpreting the full quadratic form
x^TGx, add at most another theta r/(2n), as explained below.

### 1.1 Masking large leverage scores

Let Pi=UU^T. Delete coordinates with Pi_ii>mu, let Z be the remaining
diagonal mask, and let b be the number deleted. Then b<=r/mu. Since
||G||op<=1,

```math
 |x^TGx-(Zx)^TG(Zx)|\le2\sqrt{nb}\quad(x\text{ Boolean}).  (3)
```

For off-diagonal quadratic energies, the diagonal difference contributes
at most Tr(Pi)<=r, because
sum_i |G_ii-(ZGZ)_ii|<=r. Therefore replacing offdiag(G) by offdiag(ZGZ)
costs at most theta n sqrt(b)+theta sqrt(n)r/2.

This is the theta r/(2n) term already included in (2), so no diagonal
cancellation is silently assumed.

### 1.2 Valid matched-sign flip probabilities

Put v_i=Z_ii sqrt(Pi_ii) and

```math
 G_0=ZGZ,\qquad K=vv^T.
```

Then K is PSD, K_ii<=Pi_ii<=mu on retained coordinates and zero elsewhere,
and |(G_0)_ij|<=K_ij<=mu. Write t=theta sqrt(n), and independently flip edge
ij with probability

```math
                 p_{ij}=\frac t2[K_{ij}+A_{ij}(G_0)_{ij}]. (4)
```

These probabilities lie in [0,t mu] subset[0,1]. The exact off-diagonal mean
is

```math
                 \overline A=A-t\operatorname{offdiag}G_0
                                      -t(A\circ K).      (5)
```

The extra Schur term must be paid. Since A circ K=diag(v) A diag(v),
v_i<=sqrt(mu), and a hollow quadratic form attains its absolute maximum
over [-1,1]^n at a Boolean vertex, directly

```math
 Q(A\circ K)\le\mu Q(A).                                (6)
```

The rank-one majorant and elementary cube-contraction improvement were
suggested by the independent adversarial audit. The initial proof used
K=sum_l |lambda_l||Zv_l||Zv_l|^T and the Grothendieck diagonal-majorant
lemma, giving the weaker constant 2K_G in (6). Both are valid; the
rank-one majorant removes that external ingredient entirely. Unlike the
previous midpoint-ramp application, no ramp inequality is assumed here.

### 1.3 Full-cube rounding and actual edit count

Since (sum_i v_i)^2<=n sum_i v_i²<=nr,

```math
 \sum_{i<j}p_{ij}\le tnr/2.
```

The centered edge coefficients have absolute bound2 and total variance at
most2tnr. Bernstein and the union bound over all spins give simultaneous
energy error at most sqrt(4tnr a)+(4/3)a with probability at least1/2.
Markov gives probability at most1/4 of more than2tnr changed edges. Thus one
realization has both bounds. Combining (3),(5),(6), and the explicit diagonal
correction gives (1) with

```math
 e_n\le\theta\sqrt{r/(\mu n)}+\frac{\theta r}{2n}
       +\theta C\mu\sqrt n
       +2\sqrt{\theta r a}\,n^{-3/4}
       +\frac{4a}{3n^{3/2}}.                             (7)
```

At most2theta r n^(3/2) edges change. This is an actual sign construction,
not rounding of coefficients outside the cube without a bias payment.

Put alpha=r/sqrt(n), and choose mu=alpha^(1/3)/sqrt(n). For theta alpha^(1/3)<=1,

```math
 e_n\le\theta(1+C)\alpha^{1/3}
       +\frac{\theta r}{2n}
       +2\sqrt{\theta(a/n)\alpha}
       +\frac{4a}{3n^{3/2}}.                             (8)
```

In particular r=o(sqrt(n)) yields e_n=o(1), uniformly over the entire
operator-norm unit ball of signed perturbations on ANY rank-r subspace.
For fixed r and theta this gives O(n^(-1/6)) normalized error.

### 1.4 Improved fixed-r rate by elementary interpolation

The director supplied a stronger fixed-budget estimate. Complexifying the
real bilinear norm costs at most4 by separating real and imaginary parts;
thus ||A||_(infinity->1,complex)<=4 beta(A)<=16Q(A). Since
||A||_(1->infinity,complex)=1, midpoint Riesz--Thorin gives

    ||A||op<=4 sqrt(Q(A)).

For the same rank-one K=vv^T,

    Q(A circ K)<=.5 ||A||op sum_i v_i²<=2 sqrt(Q(A)) r.

Choose mu=1/(theta sqrt(n)), so the flip probabilities remain valid and
b<=theta r sqrt(n). The construction then has normalized error at most

    [theta^(3/2)sqrt(r)+2theta sqrt(C)r
       +2sqrt(theta r a/n)]n^(-1/4)
       +theta r/(2n)+4a/(3n^(3/2)).

For fixed r,theta this is O(n^(-1/4)), improving the original n^(-1/6)
rate. It is not stronger uniformly throughout r=o(sqrt(n)); use the
minimum of the two valid error bounds. The same estimate holds under a
nuclear budget Tr|G|<=r, without a literal rank bound; that convex
extension is developed by the independent track.

## 2. Global minimizers are stationary under every such finite-rank edit

Suppose A is an exact global minimizer, Q(A)=M_n. Then (1) immediately gives

```math
 Q(A-\theta\sqrt n\operatorname{offdiag}(UTU^T))
                    \ge M_n-e_n n^{3/2}                 (9)
```

for every symmetric T with ||T||op<=1. For a near-minimizer of normalized
excess delta, the right side is Q(A)-(delta+e_n)n^(3/2). This consequence
uses actual optimality only after the valid full-sign construction.

## 3. Minimax yields a signed near-ground covariance law

For Boolean x write h(x)=H_A(x)/n^(3/2), q=Q(A)/n^(3/2), and
p(x)=U^T x/sqrt(n). Then ||p(x)||<=1. Ignoring the off-diagonal deletion
costs at most theta r/(2n), since |Tr T|<=r. Let

```math
 \varepsilon=\delta+e_n+\theta r/(2n).
```

Minimize over the compact convex operator-norm unit ball of T and maximize
over the finite signed spin set. Finite-dimensional minimax, and duality
between operator norm and nuclear norm, give a probability law mu on pairs
(sigma,x), sigma in {+1,-1}, with

```math
 \mathbb E_\mu[q-\sigma h(x)]
   +\frac\theta2
       \left\|\mathbb E_\mu[\sigma p(x)p(x)^T]\right\|_*
                       \le\varepsilon.                 (10)
```

This is optimizer-derived and uniform over the chosen rank-r feature space;
it is not an assumed isotropic Gibbs law. Both terms on the left are
nonnegative. In particular mean ground slack is at most epsilon and the
signed projected covariance has nuclear norm at most2epsilon/theta.

For any eta>epsilon, conditioning on q-sigma h(x)<=eta removes probability
at most epsilon/eta. Since every p(x)p(x)^T has nuclear norm at most1, the
conditioned signed covariance has norm at most

```math
 \frac{2\varepsilon/\theta+\varepsilon/\eta}
             {1-\varepsilon/\eta}.                       (11)
```

Thus fixed r,theta and exact optimality yield a law supported within
O(n^(-1/12)) normalized ground slack whose signed rank-r covariance is
O(n^(-1/12)). The unconditioned bounds in (10) are the sharper
O(n^(-1/6)) statements.

## 4. Exact scope: no automatic insertion or global isotropy

Equation (10) balances the positive and negative covariance contributions
on a selected subspace. It does NOT say that the two polarity probabilities
are equal: a law supported on states almost orthogonal to that subspace can
already have tiny projected covariance. Nor does it supply an isotropic
unconditional spin law, or control absolute bridge response.

An added row v requires max_x(|H_A(x)|+|v dot x|); signed covariance balance
does not upper-bound that expression. A further argument is therefore needed
before this stationarity becomes an insertion or composition theorem. No
summable recurrence is claimed here.

## 5. One-sided near-ground escape (independently audited consequence)

Suppose the negative cap is at most (q-gamma)n^(3/2), gamma>0. Then every
negative-polarity pair has slack at least gamma, so (10) gives
mu{sigma=-1}<=epsilon/gamma. Taking the trace of the signed covariance,

```math
 E[1_{sigma=+1}||p(x)||^2]
       <=2 epsilon/theta+epsilon/gamma.
```

After conditioning on positive polarity and slack at most eta, the
corresponding projected second moment is bounded above by

```math
 (2 epsilon/theta+epsilon/gamma)
       /(1-epsilon/gamma-epsilon/eta),
```

provided the denominator is positive. Thus a fixed one-sided cap gap forces
positive near-ground states to escape every chosen o(sqrt(n))-dimensional
feature space. It does not exclude a high-dimensional one-sided landscape.
The full surgery and minimax proof was independently checked in
`flatify_adversary_2026_09_07_finite_rank_stationarity_audit.md`.

## 6. Construction attempts that remain unsuccessful

### Insertion by a quadratic penalty

For a new row v, the exact cost is max_x(|H_A(x)|+|v dot x|). The elementary
completion of squares |z|-a z^2<=1/(4a) suggests a rank-one penalty along v,
but this suppresses the positive cap and enlarges the negative cap. The
signed covariance balance above does not supply a quadratic penalty with
the opposite signs on the two extremal landscapes. No insertion estimate
follows from that completion of squares.

### Cloning a sparse fraction of vertices

Inherited interactions between clone groups scale the old extremal energy
as (1+epsilon)^2, whose linear coefficient is 2 rather than the desired
3/2 for n^(3/2) scaling. Independent edge flips with density epsilon have
full-cube fluctuation of order sqrt(epsilon)n^(3/2), which is larger than
the desired linear-order correction as epsilon tends to zero. This is a
failed estimate, not a theorem excluding correlated cloning/surgery.

### Trying to enlarge the surgery rank to o(n)

The current unsigned PSD envelope K has diagonal leverage r/n in the
incoherent case; its paid Schur bias becomes order r/sqrt(n). Merely
replacing the spectral basis by another factorization does not reduce that
diagonal for a projector. A uniform entrywise envelope K=tau J, after
discarding entries |G_ij|>tau, also gives no improvement from the available
bounds: sum G_ij^2<=r yields discarded l1 at most r/tau, so its normalized
penalty is order r/(n tau), while the Schur bias is order sqrt(n)tau.
Balancing these again requires r=o(sqrt(n)).

For a rank-r cluster projector with equal clusters of size n/r, the desired
correction theta sqrt(n)G has within-cluster entry theta r/sqrt(n), exceeding
the coefficient range when r is much larger than sqrt(n). Its total
quadratic penalty can be order n^(3/2), even though the total ACTUAL signing
energy inside all those clusters is at most n^2/(2r)=o(n^(3/2)). A favorable
descent direction need not itself be regular: adding such an unnecessary
penalty to an already strong descent may preserve some net decrease.
Accordingly, a useful extension would have to produce an alternative
implementable descent, not assume every favorable direction is incoherent.
No such replacement theorem is currently proved.
