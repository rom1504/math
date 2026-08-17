# Orthogonal scout: robust Gibbs--Sion convexification

Status: literature-grounded theorem synthesis and benchmark.  Portfolio
judgment: **keep warm under a strict convexity gate**.

## Exact robust free-energy duality

Let `X` be finite and let `K` be nonempty compact convex.  Suppose
`A -> H_A(x)` is continuous and convex for every state.  Define

```math
V_infinity=\min_{A\in K}\max_{x\in X}H_A(x),

V_beta=\min_{A\in K}{1\over\beta}
 \log\sum_xe^{\beta H_A(x)}.                                  \tag{AS.1}
```

The Gibbs variational identity followed by Sion's minimax theorem gives

```math
V_infinity=\max_{p\in\Delta(X)}\min_{A\in K}\mathbb E_pH_A,

V_beta=\max_{p\in\Delta(X)}
 \left\{{S(p)\over\beta}+\min_{A\in K}\mathbb E_pH_A\right\}. \tag{AS.2}
```

Moreover,

```math
0\le V_beta-V_infinity\le{\log|X|\over\beta}.                 \tag{AS.3}
```

For an affine `d`-dimensional coupling family, an optimal zero-temperature
mixed certificate can be supported on at most `d+1` active states: at a
saddle, the active lifted affine data lie in one affine hyperplane, and
Caratheodory preserves their barycentre.

This is not a contextual state quotient.  It compresses an adversarial
certificate by convex minimax exchange.

Primary source for the exchange theorem:
[Sion, *On general minimax theorems*](https://msp.org/pjm/1958/8-1/pjm-v8-n1-p14-s.pdf).

## Zero-temperature transfer

For a sequence with `log|X_N|<=sN`, let `e_N` be the optimized ground value
divided by `N` and let `f_N(beta)` be the optimized free energy divided by
`N`.  If `f_N(beta)` converges for every fixed `beta>0`, then

```math
f(\beta)-{s\over\beta}
\le\liminf e_N\le\limsup e_N\le f(\beta).                    \tag{AS.4}
```

Sending `beta` to infinity proves

```math
\lim_Ne_N=\lim_{\beta\to\infty}f(\beta).                     \tag{AS.5}
```

No positive ground-state entropy is assumed; the lower soft-max bound uses a
single maximizing point mass.  The same normalization is standard in the
spin-glass ground-state passage, for example
[Auffinger--Chen](https://arxiv.org/abs/1606.05335).  The real missing theorem
in a robust model is fixed-temperature convergence under the outer minimum.

## Exact benchmark

Let the designer distribute total antiferromagnetic edge weight `N` on the
complete graph:

```math
K_N=\{J_{ij}\ge0:\sum_{i<j}J_{ij}=N\},
\qquad
H_J(\sigma)=\sum_{i<j}J_{ij}{1-\sigma_i\sigma_j\over2}.       \tag{AS.6}
```

Permutation averaging makes the uniform coupling optimal, while the uniform
distribution on balanced spin configurations cuts every edge with the same
probability.  Hence

```math
{1\over N}\min_J\max_\sigma H_J(\sigma)
={\lfloor N^2/4\rfloor\over\binom N2}\longrightarrow{1\over2}. \tag{AS.7}
```

At finite temperature the same symmetry reduces the pressure to a one-
dimensional binomial sum, and Laplace's principle gives

```math
f(\beta)=\max_{0\le\rho\le1}
\left[2\rho(1-\rho)+{h_2(\rho)\over\beta}\right]
={1\over2}+{\log2\over\beta}.                                \tag{AS.8}
```

This benchmark is genuinely adversarial, but convex design and permutation
symmetry make it much easier than discrete signing.

## Decisive falsifier

Convexity of the allowed coupling set cannot be relaxed.  With one spin
product `x in {-1,1}`, fixed-magnitude coupling `A in {-1,1}`, and `H_A(x)=Ax`,

```math
\min_{A=\pm1}\max_{x=\pm1}Ax=1,
\qquad
\max_p\min_{A=\pm1}A\mathbb E_px=0.                           \tag{AS.9}
```

The finite-temperature gap also tends to one.  Convexifying a discrete sign
alphabet therefore changes the problem rather than proving a relaxation with
vanishing loss.

## Promotion judgment

Keep the theory warm because it imports a real zero-entropy-safe mechanism
and gives a small mixed certificate in convex affine models.  Do not promote
it to a primary route for fixed-magnitude signings unless a new replacement
lemma proves that convexifying the coupling class costs subleading energy.
Equation (AS.9) is the kill switch against assuming that lemma.
