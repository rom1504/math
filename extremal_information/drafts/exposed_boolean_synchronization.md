# Exposed Boolean synchronization for Hadamard trust carriers

**Status.** Rigorous quantitative recovery theorem, a no-go theorem for the
stronger whole-subspace condition, and a nontrivial scalable regular-Walsh
family satisfying the exposed condition.  The result concerns the SA.3
spherical trust relaxation and its exact Boolean counterpart.

The right condition is not that an entire eigenspace be close to the cube.
That is impossible in dimension at least two.  It is enough that one
**exposed spherical optimizer** be asymptotically flat in coordinate
magnitude.

## 1. Trust response and exposed flatness

Let `H` be symmetric with `H^2=r^2I` and trace zero.  Let
`w_1,...,w_p in {+-1}^n` be ports, each repeated `m` times.  For a channel
`(sigma,epsilon) in {+-1} x {+-1}^p`, put

```math
z_\epsilon=\sum_{a=1}^p\epsilon_aw_a,
```

and define

```math
F_{\sigma,\epsilon}(u)
={\sigma\over2}u^THu+mz_\epsilon^Tu.                 \tag{BS.1}
```

The exact Boolean and spherical trust responses are

```math
\mathcal B=\max_{x\in\{+-1\}^n,\sigma,\epsilon}
                 F_{\sigma,\epsilon}(x),
\qquad
\mathcal S=\max_{\|u\|_2^2=n,\sigma,\epsilon}
                 F_{\sigma,\epsilon}(u).             \tag{BS.2}
```

For a spherical vector define its coordinate-flatness deficit

```math
\phi(u)=1-{\|u\|_1\over n}.                          \tag{BS.3}
```

Because `||u||_2^2=n`, Cauchy--Schwarz gives `phi(u)>=0`.  Coordinatewise
sign rounding `x=sgn(u)` satisfies the exact identity

```math
{\|x-u\|_2^2\over n}=2\phi(u),                       \tag{BS.4}
```

with arbitrary signs at zero coordinates.

## 2. Quantitative Boolean recovery

### Theorem BS.1 (flat exposed optimizer recovery)

Let

```math
c={mp\over r}.                                       \tag{BS.5}
```

Suppose some channel and spherical vector `u`, `||u||_2^2=n`, obey

```math
F_{\sigma,\epsilon}(u)\ge\mathcal S-\xi rn,
\qquad \phi(u)\le\varphi.                            \tag{BS.6}
```

Then

```math
\boxed{
0\le\mathcal S-\mathcal B
\le rn\left[\xi+(1+c)\sqrt{2\varphi}\right].}       \tag{BS.7}
```

In particular, at bounded total port mass `c=O(1)`, an exposed optimizer
with `xi=o(1)` and `phi(u)=o(1)` gives

```math
\mathcal S-\mathcal B=o(rn).                         \tag{BS.8}
```

#### Proof

Let `x=sgn(u)` and `delta=||x-u||_2/sqrt(n)`.  Since `||H||=r` and both
vectors have norm `sqrt(n)`,

```math
\left|{1\over2}x^THx-{1\over2}u^THu\right|
\le rn\delta.                                        \tag{BS.9}
```

Also `||z_epsilon||_2<=p sqrt(n)`, so

```math
m|z_\epsilon^T(x-u)|
\le mpn\delta=crn\delta.                             \tag{BS.10}
```

Evaluate the same channel at `x`, use (BS.6), and substitute
`delta=sqrt(2phi(u))` from (BS.4).  This proves (BS.7). `square`

The condition is checkable without any Boolean optimization table.  Solve
the spherical trust problem, choose one exposed optimizer, and evaluate one
`l_1` norm.  A robust version permits a near optimizer through `xi`.  The
condition stores one rounded witness, not the responses of all `2^n` cube
points.

It is also close to the weakest condition available to coordinatewise sign
rounding: (BS.4) says that `phi` is exactly its squared Euclidean distortion.

## 3. Why a whole-subspace Boolean net is the wrong condition

One might require every normalized vector in the port/eigenspace feature
span to lie near the cube.  This collapses to dimension one.

### Proposition BS.2 (no uniformly Boolean linear sphere)

Let `U subset R^n` be a `d`-dimensional subspace.  Then some
`u in U`, `||u||_2^2=n`, satisfies

```math
{\|u\|_1\over n}
\le\gamma_d
:=\sqrt d\,{\Gamma(d/2)\over
                    \sqrt\pi\,\Gamma((d+1)/2)}.       \tag{BS.11}
```

For every `d>=2`, `gamma_d<1`.  Consequently the sphere of `U` has Boolean
covering radius at least

```math
\sup_{\substack{u\in U\\\|u\|^2=n}}
 \min_{x\in\{+-1\}^n}{\|u-x\|_2\over\sqrt n}
\ge\sqrt{2(1-\gamma_d)}>0.                           \tag{BS.12}
```

#### Proof

Let `Q` be an `n by d` matrix with orthonormal columns and let `theta` be
uniform on the unit sphere in `R^d`.  Set `u=sqrt(n)Qtheta`.  If `q_i` is
row `i`, rotational invariance gives

```math
\mathbb E|q_i^T\theta|
={\Gamma(d/2)\over\sqrt\pi\Gamma((d+1)/2)}\|q_i\|_2.
```

Now `sum_i||q_i||^2=d`, so Cauchy--Schwarz gives

```math
\mathbb E\|u\|_1
\le n\gamma_d.
```

Some `u` attains at most the average.  Strict Cauchy--Schwarz in the random
coordinate `theta_1` gives
`E|theta_1|<sqrt(E theta_1^2)=1/sqrt(d)` for `d>=2`, hence `gamma_d<1`.
Finally use (BS.4). `square`

Thus a uniform Boolean net for an entire fixed two-dimensional trust span
already demands a fixed error.  Recovery must use **exposure** or
synchronization, not uniform subspace flatness.

## 4. Two correlated top ports: exact formula

The exposed condition can hold nontrivially when distinct poles synchronize.

### Lemma BS.3 (correlated top-pole response)

Suppose `a,b` are Boolean `+r` eigenvectors with

```math
{a^Tb\over n}=\rho\in[0,1].                          \tag{BS.13}
```

For two `m`-wide ports,

```math
\mathcal B={rn\over2}+mn(1+\rho),                    \tag{BS.14}
```

and

```math
\mathcal S={rn\over2}+mn\sqrt{2(1+\rho)}.            \tag{BS.15}
```

The spherical optimum is exposed by

```math
u_*={a+b\over\sqrt{2(1+\rho)}}.                      \tag{BS.16}
```

#### Proof

For Boolean endpoints,

```math
\|a+b\|_1=n(1+\rho),
\qquad \|a-b\|_1=n(1-\rho).                          \tag{BS.17}
```

The first is the larger support.  Together with the separate child bound
`rn/2`, it proves the upper bound in (BS.14), attained at `x=a`.

On the sphere, the larger field norm is
`||a+b||_2=sqrt(2n(1+rho))`.  Cauchy--Schwarz and the child spectral bound
give (BS.15), and (BS.16) attains both because it remains in the positive
eigenspace. `square`

The exact normalized gap is

```math
{\mathcal S-\mathcal B\over rn}
={m\over r}\left[\sqrt{2(1+\rho)}-(1+\rho)\right].   \tag{BS.18}
```

It vanishes when `rho->1`, even though `a,b` remain distinct.

## 5. A scalable regular-Walsh synchronized pair

We now construct distinct top poles whose disagreement support has the exact
scale `sqrt(n)`.  No optimal uncertainty claim is needed here.

Let `d>=2` be even, put `q=2^d`, `n=q^2`, and index the Sylvester Walsh
matrix `W` by

```math
V=\mathbb F_2^d\times\mathbb F_2^d.
```

Define the self-dual bent sign sequence

```math
y_0(x,z)=(-1)^{x\cdot z}.                             \tag{BS.19}
```

Direct Walsh summation gives

```math
Wy_0=qy_0.                                            \tag{BS.20}
```

Let `M` be the fixed-point-free coordinate-pairing permutation on
`F_2^d`; it is symmetric, has zero diagonal, and obeys `M^2=I`.  Set

```math
L=\{(x,Mx):x\in\mathbb F_2^d\}.                      \tag{BS.21}
```

Then `L=L^perp`, `|L|=q`, and `x dot Mx=0`.  Consequently

```math
W\mathbf1_L=q\mathbf1_L,
\qquad y_0|_L=1.                                     \tag{BS.22}
```

Therefore

```math
y_1=y_0-2\mathbf1_L                                  \tag{BS.23}
```

is another Boolean `+q` Walsh eigenvector, differing from `y_0` on exactly
`q` coordinates.

Regularize by

```math
H=D_{y_0}WD_{y_0}.                                   \tag{BS.24}
```

This is symmetric Hadamard, `H^2=nI`, `tr H=0`, and
`H1=q1`.  Its two Boolean `+q` eigenvectors

```math
a=\mathbf1,
\qquad b=D_{y_0}y_1=\mathbf1-2\mathbf1_L             \tag{BS.25}
```

obey

```math
\rho=1-{2\over q}.                                   \tag{BS.26}
```

This is a genuine two-dimensional port span at every order, but its small
Gram eigenvalue is `2/q`: the ports synchronize asymptotically rather than
becoming literally equal.

### Theorem BS.4 (nontrivial synchronized recovery family)

Use the ports (BS.25) with width `m=q/2`.  Then the total port mass is

```math
c={2m\over r}=1,                                     \tag{BS.27}
```

and the exposed optimizer (BS.16) has

```math
{\|u_*\|_1\over n}=\sqrt{1-{1\over q}},
\qquad
\phi(u_*)=1-\sqrt{1-{1\over q}}=o(1).                \tag{BS.28}
```

Hence Theorem BS.1 already gives `S-B=o(rn)`.  More sharply, the exact gap
from BS.3 is

```math
\boxed{
{\mathcal S-\mathcal B\over rn}
=\sqrt{1-{1\over q}}-\left(1-{1\over q}\right)
=O(q^{-1})=O(n^{-1/2}).}                             \tag{BS.29}
```

This is strictly beyond a single Boolean pole.  It also identifies the
mechanism precisely: a sparse eigendirection changes
one pole on `sqrt(n)` coordinates, and the exposed normalized sum remains
asymptotically flat.

Any exact-sign completion on the total `2m=q` auxiliary vertices perturbs
both Boolean and spherical responses by `O(q^2)=O(n)=o(qn)`, so the recovery
conclusion survives completion, although the exact coefficient in (BS.29)
need not.

## 6. Information content and boundary

Once a globally exposed channel is known, the exposed-flatness condition
does **not** encode an exponential response object:

- the spherical optimizer is obtained from the Gram/trust problem;
- its flatness is one `l_1` statistic;
- its rounded Boolean witness is obtained coordinatewise;
- no Boolean maximization or response histogram is assumed.

It is weaker than uniform approximation of the trust span, which BS.2 rules
out.  Identifying a globally exposed channel may itself require searching
the `2^(p+1)` trust table unless additional structure supplies it, so this is
query-local recovery rather than a reusable carrier theorem.  Nor is
flatness automatically information beyond `(G,R)`: for the correlated
two-pole family it equals `sqrt((1+rho)/2)` and is already Gram-determined.
The orthogonal case shows instead that a small spherical state can decode
the wrong discrete value by a fixed amount.

The remaining boundary is now crisp.  A reusable Boolean carrier needs a
mechanism forcing at least one globally exposed trust direction to have
vanishing coordinate-flatness deficit.  Asymptotic pole synchronization is
one such mechanism.  Fixed-angle independent poles are not.

## 7. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_exposed_boolean_synchronization.py
```

The verifier checks BS.1 numerically on finite channels, exhausts the exact
order-16 Boolean response, constructs the self-dual subspace family at
orders `16` and `256`, and verifies every Walsh, eigenvector, flatness, and
gap identity exactly where algebraic.
