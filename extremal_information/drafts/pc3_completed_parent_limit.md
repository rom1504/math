# A completed-parent limit for the growing product-algebra family

**Status.** Rigorous near-original restricted theorem, independently
audited.  The product-algebra response state predicts an exact asymptotic cap
for a dense exact-sign family with a growing `sqrt(n)`-vertex interface.  It
does not optimize over signings and its value `3/4` is not competitive with
the motivating upper bound.

## 1. Seed row law

Use the order-16 regular-Hadamard triple `a,b,c` from PC.3.  Its normalized
Gram data are

```math
\langle a,b\rangle/n=1/2,\qquad
\langle a,c\rangle/n=0,\qquad
\langle b,c\rangle/n=-1/2.                       \tag{PL.1}
```

Take `a` as the affine-coset anchor and put

```math
X=a\odot b,\qquad Y=a\odot c.                    \tag{PL.2}
```

At a uniformly sampled seed coordinate,

```math
E X=1/2,\qquad E Y=0,\qquad E(XY)=-1/2.           \tag{PL.3}
```

These three moments determine the two-bit law:

```math
P((X,Y)=(1,1),(1,-1),(-1,1),(-1,-1))
=(1/4,1/2,1/4,0).                                \tag{PL.4}
```

At tensor level `j`, projectivize every port row by the anchor port.  The
resulting row pattern is exactly

```math
(1,X_1,Y_1,\ldots,X_j,Y_j),                      \tag{PL.5}
```

where the pairs are independent with law (PL.4).  Thus `p_j=2j+1`.

## 2. Exact maximum support asymptotics

### Lemma PL.1

Let

```math
L_j=\max_{\epsilon_0,\alpha_t,\beta_t\in\{+-1\}}
E\left|\epsilon_0+\sum_{t=1}^j
             (\alpha_tX_t+\beta_tY_t)\right|.    \tag{PL.6}
```

Then

```math
{j/2+1\over2j+1}\le {L_j\over2j+1}
\le {j/2+1+\sqrt{11j/4}\over2j+1},               \tag{PL.7}
```

and in particular

```math
\boxed{L_j/p_j\longrightarrow1/4.}               \tag{PL.8}
```

#### Proof

For `Z_t=alpha_tX_t+beta_tY_t`, (PL.3) gives

```math
E Z_t=\alpha_t/2.                                 \tag{PL.9}
```

Moreover `Z_t` is bounded by two, and direct use of (PL.4) gives
`Var(Z_t)<=11/4`.  If `S` is the expression inside the absolute value,
independence and Cauchy--Schwarz give

```math
E|S|\le |ES|+E|S-ES|
\le 1+j/2+\sqrt{11j/4}.                           \tag{PL.10}
```

Taking every `alpha_t=1` and `epsilon_0=1`, Jensen gives
`E|S|>=|ES|=1+j/2`.  Divide by `2j+1`. `square`

If `W_j` is the PC.3 port matrix, its rows have the uniform tensor law
(PL.5).  Therefore

```math
{1\over n_j}\max_\epsilon||W_j\epsilon||_1=L_j. \tag{PL.11}
```

## 3. A dense exact-sign completed-parent limit

Let

```math
n_j=16^j,\qquad r_j=4^j=\sqrt{n_j},\qquad
m_j=\lfloor r_j/p_j\rfloor.                      \tag{PL.12}
```

Delete the diagonal from the trace-zero symmetric Hadamard tensor `H_j`.
This is a hollow exact signing with the same Boolean quadratic energy.
Attach `m_j` identical vertices for each of the `p_j` ports, and fill all
edges among the `q_j=p_jm_j<=r_j` new vertices by **any** hollow exact
signing `C_j`.  Call the completed signing `P_j`; its order is
`N_j=n_j+q_j`.

### Theorem PL.2 (restricted completed-parent thermodynamic limit)

Every such public completion obeys

```math
\boxed{{Q(P_j)\over N_j^{3/2}}\longrightarrow {3\over4}.}    \tag{PL.13}
```

#### Proof

Odd product closure and the same-selector theorem give the exact pre-
completion cap

```math
\max_\epsilon B_\epsilon
={r_jn_j\over2}+m_jn_jL_j.                       \tag{PL.14}
```

Since `m_j/r_j=1/p_j+O(1/r_j)`, Lemma PL.1 implies

```math
{1\over r_jn_j}\max_\epsilon B_\epsilon
={1\over2}+{m_j\over r_j}L_j\longrightarrow {3\over4}.      \tag{PL.15}
```

The completion changes the cap by at most

```math
Q(C_j)\le {q_j\choose2}=O(r_j^2)=O(n_j)=o(r_jn_j).           \tag{PL.16}
```

Finally `q_j/n_j<=r_j/n_j=1/r_j` tends to zero and
`r_jn_j=n_j^(3/2)`, so replacing `n_j` by total order `N_j` does not change
the limit. `square`

## 4. Robust extension and scope

The same conclusion holds for a sequence satisfying the hypotheses of the
robust selector theorem with intrinsic defects `Delta_j=o(1)`: (PL.14) then
has additive error at most `Delta_j r_jn_j/2`.  In particular, the Cartesian
relative theorem supplies this when factorwise pole defects have total
`sum_i delta_i=o(1)` and the seed row law is unchanged.

This is a theorem about a declared structured family, not about its optimality.
Its `3/4` limit is worse than the known `1/2`-scale constructions, the orders
are `16^j+O(4^j)`, and no all-order recovery or statement about `M_n` follows.
Its role is a proof-of-value benchmark: the response algebra selected a
growing dense interface, computed its macroscopic response from a two-bit
factor law, and remained valid after arbitrary exact-sign microscopic
completion.
