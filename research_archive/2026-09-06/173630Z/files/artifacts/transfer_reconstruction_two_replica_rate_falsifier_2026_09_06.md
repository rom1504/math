# Near-identical replicas falsify the proposed annealed Gaussian rate inequality

Date: 2026-09-06. This is a scoped Gaussian rate-function falsifier. It is
not an evaluation of an actual signing second moment; rare joint row
profiles remain a separate obligation.

## 1. Covariance optimization and its endpoints

Let `C_q=[[1,q],[q,1]]`, `0<=q<1`, and require both corresponding
endpoint-replica correlations to be `rho in (0,1)`. Gaussian mutual
information is minimized over the cross-covariance matrix. Convexity of
minus log determinant allows averaging under transposition and replica
interchange, giving `D=[[rho,r],[r,rho]]`. The resulting information is

```math
I_2(q,\rho)=\min_r\left\{\log(1-q^2)
-\frac12\log\left([(1+q)^2-(\rho+r)^2]
                  [(1-q)^2-(\rho-r)^2]\right)\right\}.
```

The feasible interval is
`rho+q-1<r<1-|q-rho|`. After diagonalizing both matrices, this is
`j((rho+r)/(1+q))+j((rho-r)/(1-q))`, where
`j(z)=-log(1-z^2)/2`. Strict convexity of `j` makes the minimizer unique;
its stationary equation is

```math
r^3-(1+q^2+\rho^2)r+2q\rho=0.
```

At `q=0` it is `r=0`, giving `I_2(0,rho)=2j(rho)`. At the degenerate
endpoint `q=1`, the two coordinates of each row coincide and
`I_2(1,rho)=j(rho)`. The same value is the limit from `q<1`, by data
processing from one coordinate and the feasible trial below.

## 2. A feasible trial gives an exact positive lower witness

At the annealed relation `alpha=1-rho^2=2^(-4p)`, define

```math
F_p(q)=-p[\log2-h((1+q)/2)]
+\frac12[I_2(0,\rho)-I_2(q,\rho)].
```

Then `F_p(1)=0`. Put `x=(1-q)/2`, and choose the feasible trial
`r=rho`, requiring `rho<1-x`. Its antisymmetric-mode correlation is
zero, while its symmetric-mode correlation is `rho/(1-x)`. Since this
trial upper-bounds the minimizing information,

```math
F_p(1-2x)\ge p h(x)
+\frac14\log\frac{\alpha-2x+x^2}{\alpha(1-x)^2}.
\tag{1}
```

No stationary point or numerical optimizer is used in (1).

Choose `x=1/2000`, so `q=999/1000`. For every
`p in [31/32,1]`, `alpha>=1/16`; the logarithm in (1) increases with
`alpha`. Also `(1-x)^2>15/16`, proving trial feasibility uniformly.
Thus

```math
F_p(999/1000)
\ge\frac{31}{32}h(1/2000)
+\frac14\log\frac{3936016}{3996001}
\ge\frac{24619155871531}{64000000000000000}>0.
\tag{2}
```

The final lower value is about `0.000384674310492672`. For `p=1` alone
the same trial proves the stronger lower value
`1038119009001/2000000000000000`, about `0.0005190595045005`.
The reproducible exact verifier is

```sh
.venv/bin/python computations/transfer_reconstruction_two_replica_rate_falsifier_2026_09_06.py
```

It uses the established outward rational logarithm and entropy intervals.
All comparisons in (2), including covariance feasibility, are exact.

## 3. The obstruction occurs for every fixed positive retention

For fixed `p>0`, expansion of the explicit trial (1) gives

```math
F_p(1-2x)\ge
p x\log(1/x)+
\left[p-\frac{\rho^2}{2(1-\rho^2)}\right]x+O_p(x^2).
\tag{3}
```

The leading `x log(1/x)` term is positive. Thus the proposed global
inequality `F_p(q)<=0` fails near `q=1` for every fixed `p>0`, not only
at the two retention fractions in (2). The optimizing information has
the same linear endpoint expansion: its small antisymmetric-mode
correlation is `O(x)`, and
`I_2(1-2x,rho)=j(rho)+rho^2 x/(1-rho^2)+O(x^2)`.

The mechanism is the entropy of near-identical replicas. At the annealed
endpoint the rate is exactly balanced, but the information penalty varies
linearly in the mismatch density while the spin-pair entropy gains
`x log(1/x)`. This invalidates the proposed Gaussian rate inequality at
criticality. It neither proves an actual second-moment lower bound nor
excludes a cluster-normalized argument or a subcritical energy analysis.
