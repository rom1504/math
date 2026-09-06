# Exact equality criterion for the conditional-variance supersolution

Date: 2026-09-06. Status: verified analytic theorem and exact finite-source
separation certificate. The equality proof was independently scrutinized
by the campaign's transfer-seeds researcher. This concerns the analytic
Gaussian-boundary Bellman certificate, not actual ensemble-pressure equality
or convergence of the original minimax.

Fix t>0. Use the original binary Bellman operator B, its Gaussian
potential g_t, the classwise envelope

```math
T(\nu)=\sup_L\{\mathbb E g_t(\operatorname{Var}(X\mid L))-I(X;L)\},
```

and the older average-variance envelope

```math
E(\nu)=\sup_L\{g_t(\mathbb E\operatorname{Var}(X\mid L))-I(X;L)\}.
```

Let H denote the common increasing/decreasing Gaussian-boundary Bellman
limit established in the audited dependencies. We know `E<=H<=T`,
`B H=H`, and `B T<=T`.

## 1. The equality theorem

**Theorem.** For every finite symmetric real source nu,

```math
(BT)(\nu)=T(\nu)
\quad\Longleftrightarrow\quad E(\nu)=T(\nu)
\quad\Longleftrightarrow\quad H(\nu)=T(\nu).
```

In the equality case there exists a T-optimal channel with constant
conditional variance. Consequently this constant-variance attainment is
also equivalent to the three conditions above.

### Attainment and safe symmetrization

Fix the finite source alphabet. Its child alphabets of normalized sums
and differences are finite and fixed. The posterior representation of T
is a concave-envelope optimization of the continuous function
`H(posterior)+g_t(Var(posterior))`, minus source entropy. A finite number
of posterior labels suffices by Caratheodory, and the maximum is attained.
T is continuous on the probability simplex: compact posterior-mixture
representations give upper semicontinuity, while each fixed finite forward
channel gives a continuous objective and hence lower semicontinuity.

The admissible pair polytope is compact, and its entropy cost relative to
the fixed source product is continuous. Thus BT attains its maximum.
Average that pair under simultaneous sign reversal and input interchange.
These operations preserve the two output absolute laws individually and
cannot increase relative entropy. The resulting actual signed marginals
are BOTH nu, and its cost is I(A;B). This step is essential; an averaged
absolute-marginal condition alone would not justify equal parent values.

### All equality slacks

Suppose BT(nu)=T(nu). Choose such an attaining pair (A,B), define
`U=(A+B)/sqrt(2)`, `V=(A-B)/sqrt(2)`, and choose optimal finite child
channels M|U and N|V independently given (U,V). Put L=(M,N).
Use parent labels `(L,B)` for A and L for B.

Let R_parent and R_child be their respective total averaged g rewards,
and I_parent and I_child their respective total channel informations.
The exact information identity is

```math
I_{\rm parent}=I_{\rm child}+I(A;B)-I(M;N).
```

The supersolution proof gives `R_parent>=R_child`. Thus

```math
2T(\nu)\ge R_{\rm parent}-I_{\rm parent}
\ge R_{\rm child}-I_{\rm child}-I(A;B)=2(BT)(\nu).
```

Equality at the endpoints forces BOTH parent channels to be optimal,
`I(M;N)=0`, and equality in every nonnegative covariance, Jensen, and
label-refinement slack in the reward comparison.

### Strict covariance comparison

Fix L=l and write its conditional covariance as Sigma=(a,c;c,b).
For positive definite Sigma, strict concavity of g_t(exp(s)) forces
equality between its Schur-pivot reward and eigenvalue reward only if
the two multisets of pivots and eigenvalues agree. This requires b to
be an eigenvalue and hence c=0. The Hadamard-rotated diagonal is then
((a+b)/2,(a+b)/2). Strict ordinary convexity of g_t forces a=b.
Therefore equality implies Sigma=v(l)I.

There is no nonzero singular exception. If b>0 and det Sigma=0, the
Schur-pivot pair is (0,b) and the eigenvalue pair is (0,a+b). Strict
decrease of g makes that comparison strict unless a=0. If a=0, c=0 and
the next comparison is `g(0)+g(b)>2g(b/2)` for b>0. If b=0,a>0 the same
last strict inequality applies with a. The zero matrix is allowed.

Strict Jensen and optimal linear-prediction equality additionally imply
`Var(A|B,L)=v(L)` almost surely. The second parent has
`Var(B|L)=v(L)` by definition.

### Independence eliminates heterogeneous residual variances

Equality in label refinement uses the two inequalities

```math
\mathbb E[g_t(\operatorname{Var}(U\mid M,N))\mid M]
\ge g_t(\mathbb E[\operatorname{Var}(U\mid M,N)\mid M])
\ge g_t(\operatorname{Var}(U\mid M)).
```

Strict convexity and strict decrease force the conditional variances to
agree, so `v(M,N)=Var(U|M)`. The analogous V argument gives
`v(M,N)=Var(V|N)`. But M and N are independent, since their mutual
information vanished. A function of M equal almost surely to a function
of independent N must be constant. Hence v(L) is a constant v_0.

The optimal parent channels found above therefore have constant posterior
variance v_0. Their objective values are unchanged when evaluated in E
instead of T. Thus E(nu)>=T(nu), proving E=T.

Conversely, E=T implies H=T by E<=H<=T. Then monotonicity and BH=H give
`BT>=BH=H=T`, while the supersolution gives the opposite inequality.
If H=T directly, the same argument gives BT=T, closing the equivalences.
Finally any constant-variance T-optimal channel gives E=T directly, while
the forward implication above constructs one. This proves the theorem.

## 2. An explicit finite source where the supersolution is strictly loose

Take t=1, R=2048, m=2R+1=4097, and let Z be uniform on the integers
`{-R,...,R}`. Its variance is `v=R(R+1)/3=1398784`. Define

```math
\nu={3\over4}\delta_0+{1\over4}\mathcal L(Z).
```

The following bounds prove `T(nu)-E(nu)>3/20` without optimizing any
channel or solving any Bellman problem numerically.

### Lower bound on T

Reveal the mixture component and use no additional label within either
component. This costs at most h(1/4), and yields
`T(nu)>=g_1(v)/4-h(1/4)`.

For v>=1, test the Gaussian self-transport problem with correlation
`rho=1-1/(2v)`. Its cost is
`1+(1/2)log v-(1/2)log(1-1/(4v))`.
Since g is minus one half the optimal cost,

```math
g_1(v)\ge-{1\over2}-{\log v\over4}-{\log(4/3)\over4}.
```

Consequently

```math
T(\nu)\ge-h(1/4)-{1\over8}-{\log v\over16}-{\log(4/3)\over16}.
```

### Upper bound on E

Write `J_lambda(X)=inf_L[I(X;L)+lambda E Var(X|L)]`. Source concavity
gives `J_lambda(nu)>=J_lambda(Z)/4`. It follows directly from the
reproduction-distribution representation: J is an infimum of linear
source expectations.

For any Z-channel let D=E Var(Z|L). Append an independent uniform
[-1/2,1/2] variable W. Disjoint unit intervals give
`h(Z+W|L)=H(Z|L)`. Gaussian maximum entropy, followed by concavity of log,
gives

```math
H(Z\mid L)\le{1\over2}\log(2\pi e(D+1/12)).
```

Thus for 0<lambda<=1, minimizing the resulting scalar bound over D>=0
(the optimum is `D=1/(2lambda)-1/12>=0`) yields

```math
J_\lambda(Z)\ge\log m+{1\over2}\log(\lambda/\pi)-{\lambda\over12}.
```

Use `E(nu)=sup_(0<lambda<=1)[c_1(lambda)-J_lambda(nu)]` and
`c_1(lambda)<=log(2lambda)/4`. The coefficient of log lambda is 1/8,
which is nonnegative, while log lambda<=0 and lambda<=1. Therefore

```math
E(\nu)\le-{\log m\over4}+{\log2\over4}
                  +{\log\pi\over8}+{1\over48}.
```

### Exact outward certificate

The same-prefix certificate script uses the already preserved rational
log/entropy interval primitives. It replaces pi by the classical upper
bound 22/7. Its exact output is

```math
T(\nu)\ge-{25436158232087\over16000000000000},\qquad
E(\nu)\le-{10453445417251\over6000000000000},
```

and therefore

```math
T(\nu)-E(\nu)\ge{7319088641747\over48000000000000}
>{3\over20}.
```

The equality theorem now proves `BT(nu)<T(nu)` and `H(nu)<T(nu)`.
It does NOT give a numerical lower bound on T-BT or T-H from the displayed
E/T gap. In particular the classwise envelope is provably not the deep
Bellman fixed point on all finite sources.

## 3. Relation to other possible conclusions

This result does not refute the older temperature-alignment conjecture
`BE<=E`: it leaves open whether H=E generally. It explains precisely when
the larger, successfully certified envelope T is already exact. Symmetric
binary sources satisfy the equality condition; the explicit lattice
mixture above does not.

For the original ternary root, the accompanying finite-grid BT/T tests
suggest a substantial strict gap, but the explicit finite separation here
is a different source. No numerical statement about the ternary root is
promoted by this theorem, and no improved original-signing cap is claimed.

## 4. Exact strictness at the actual ternary root

For p=31/32,t=4, use two reflected posterior orbits with
`(z,s)=(1363/1500,0)` and `(374/375,1)`, with first-orbit weight 49/152.
Their barycentre is exactly the ternary root. Their conditional variances
are respectively 10904/11625 and 11968/4359375. The same-prefix ternary
strict-certificate script evaluates this feasible channel outward in
rational arithmetic. In root-offset notation `offset(F)=p log2+F+4(1-sqrt p)`,
it gives

```math
\operatorname{offset}(T)\ge
-{63833371773687289\over4240800000000000000}.
```

The old E certificate was independently replayed to the identical endpoint
`offset(E)<=-19678127864847/800000000000000`. Therefore

```math
T(\nu_p)-E(\nu_p)\ge
{20240192018933329\over2120400000000000000}>{9\over1000}.
```

Thus `BT(nu_p)<T(nu_p)` and `H(nu_p)<T(nu_p)` are proved at the ACTUAL
construction root, not only at the separate wide-lattice example. This
does not quantify T-BT; the displayed lower bound is for T-E.

## 5. Subsequent director breakthrough

After this equality investigation, the director supplied a precision-side
Schur argument proving the previously open `BE<=E`. It passed independent
reconstruction here and in the transfer-seeds track; see
`decisive_bridge_precision_alignment_independent_audit_2026_09_06.md`.
Consequently H=E is now proved, and the exact T-E gaps above ARE exact
lower bounds on T-H. Statements above describing temperature alignment as
open record the chronology before that breakthrough; they are superseded
by this update. No original-minimizer convergence follows automatically.
