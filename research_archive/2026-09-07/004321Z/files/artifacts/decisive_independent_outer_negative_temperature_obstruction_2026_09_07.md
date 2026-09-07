# Why the outer negative-temperature disorder partition does not import SK universality

Status: exact normalization and a proved Gaussian replacement obstruction.
This records a failed direct thermodynamic-limit route, not a new convergence
criterion being promoted as progress.

For the two-sided normalized Ising pressure define

```math
P_{n,\beta}(A)={1\over n}\log\sum_{s=\pm1,x\in\{\pm1\}^n}
 \exp\{s\beta q_A(x)/\sqrt n\},
\qquad p_{n,\beta}=\min_A P_{n,\beta}(A).
```

The discrete outer partition

```math
\Psi_{n,\beta,\gamma}
=-{1\over\gamma n^2}\log\mathbb E_{A\ {m iid\ signs}}
 \exp\{-\gamma n^2P_{n,\beta}(A)\}
```

satisfies the exact finite estimate

```math
p_{n,\beta}\le\Psi_{n,\beta,\gamma}
\le p_{n,\beta}+{\binom n2\log2\over\gamma n^2}.
```

Thus a fixed-(beta,gamma) thermodynamic limit, followed by gamma and beta
limits, would suffice. But the outer exponent equals
`-gamma n log Z`, not a fixed negative replica number. The disorder entropy
has order n² while the inner log partition has order n. A fixed replica
interpolation therefore has the wrong scale for approximation of the pure
minimum.

## Gaussian replacement changes the variational problem

Replace the outer sign disorder by independent standard Gaussian entries,
and denote the resulting quantity by Psi^G. Every hollow real matrix has
two-sided pressure at least `(1+1/n) log2`, by Jensen under uniform spins
and the objective sign. Therefore Psi^G has the same lower bound.

Use the Gaussian law with entry variance sigma² as an outer trial measure.
The finite entropy variational principle gives

```math
\Psi^G_{n,\beta,\gamma}
\le\mathbb E P_{n,\beta}(\sigma G)
 +{\binom n2\over2\gamma n^2}
    (\sigma^2-1-\log\sigma^2).
```

The elementary Gaussian maximum bound supplies

```math
\mathbb E P_{n,\beta}(\sigma G)
\le(1+1/n)\log2
 +\beta\sigma\sqrt{(1+1/n)\log2}.
```

Choose sigma=1/gamma for gamma>=1. It follows uniformly in n that

```math
0\le\Psi^G_{n,\beta,\gamma}-(1+1/n)\log2
\le {\beta\sqrt{2\log2}\over\gamma}
 +{\gamma^{-2}-1+2\log\gamma\over4\gamma}.
```

Hence the Gaussian outer variational problem collapses to zero interaction
energy as gamma tends to infinity. This is not the sign problem: for signs,
`p_(n,beta)>=beta M_n/n^(3/2)`, and the banked positive lower cap constant
separates this from log2 at sufficiently large fixed beta. The Gaussian
trial law can shrink every edge variance; the sign constraint cannot.

This rules out uniform Gaussian replacement at the outer optimization scale.
It does not rule out a genuinely discrete interpolation or a non-Gaussian
second-order graph-limit theorem. Neither was obtained here.
