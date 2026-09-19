### Wave 23, Route 2: exact variation identities and a rare-tail obstruction

**Status.** All identities and the $A_9$ calculations below are **Verified**.
The checker is `tmp/conditional_free_energy_r23.py`. The resulting bounds do
not prove (10.723) or (10.727), but they isolate the missing uniform tail
estimate and rule out a base-variance shortcut.

Write $T=[n]\setminus S$ and represent an oriented cut as
$d_{ij}=\sigma x_ix_j$. Put

```math
e_S(d)=\langle A,d\rangle-c_S(d),\qquad
K_{\beta,S}(y)=\sum_{d:d[S]=y}e^{\beta e_S(d)},\qquad
F_S(y)=\log K_{\beta,S}(y).
```

#### Child-state variation and the conditional KL chain

Flip one child spin $i\in S$, writing the resulting state as $y^i$, and set

```math
h_{i,T}(d)=\sum_{v\in T}a_{iv}d_{iv}.
```

Since $e_S(d^i)-e_S(d)=-4h_{i,T}(d)$, exact exponential tilting gives

```math
\boxed{
F_S(y^i)-F_S(y)
=\log\mathbb E_{P_y}e^{-4\beta h_{i,T}},
}
\tag{R23.1}
```

where $P_y=\nu_\beta(\cdot\mid D[S]=y)$. Moreover,

```math
\boxed{
D_{\rm KL}(P_y\Vert P_{y^i})
=4\beta\mathbb E_{P_y}h_{i,T}+F_S(y^i)-F_S(y).
}
\tag{R23.2}
```

Let $T_i$ flip coordinate $i$ on the full cut, let $\nu_{\beta,S}$ be the
parent marginal on $S$, and put
$H_i(d)=\sum_{j\ne i}a_{ij}d_{ij}$. Chain rule under the deterministic map
$d\mapsto d[S]$ gives the stronger exact decomposition

```math
\boxed{
4\beta\mathbb E_{\nu_\beta}H_i
=D_{\rm KL}(\nu_{\beta,S}\Vert T_i\nu_{\beta,S})
+\mathbb E_{y\sim\nu_{\beta,S}}
 D_{\rm KL}(P_y\Vert P_{y^i}).
}
\tag{R23.3}
```

Both terms on the right are nonnegative. Averaging over a uniform $m$-set
and a uniform $i\in S$ therefore yields

```math
\mathbb E_{S,i,y}D_{\rm KL}(P_y\Vert P_{y^i})
\le {4\beta\over n}\mathbb E_{\nu_\beta}\langle A,D\rangle
\le {4\beta Q(A)\over n}.
\tag{R23.4}
```

For an exact minimizer and
$\beta=\Theta(n^{-1/2+c})$, this is $O(n^c)$. This is a real
minimizer-scale gain, but it is under the parent marginal
$\nu_{\beta,S}$, whereas (10.723) compares the child Gibbs law to that
marginal. Summing a coordinatewise estimate over $\Theta(n)$ flips would
also cost $O(n^{1+c})$.

#### Selector-neighbor identity

For adjacent selectors $S=C\cup\{u\}$ and $S'=C\cup\{v\}$ define

```math
G_S(y)=\beta c_S(y)+F_S(y)
=\log Z_A(\beta)+\log\nu_{\beta,S}(y).
```

If $D\sim\nu_\beta$, Bayes' rule gives

```math
\boxed{
G_S(D[S])-G_{S'}(D[S'])
=\log\Pr(D_u\mid D[C])-\log\Pr(D_v\mid D[C]).
}
\tag{R23.5}
```

Each new oriented spin is binary conditional on $D[C]$. Since
$p\log^2(1/p)\le4/e^2$ and
$(a-b)^2\le2a^2+2b^2$,

```math
\boxed{
\mathbb E_{D\sim\nu_\beta,\,S\sim S'}
[G_S(D[S])-G_{S'}(D[S'])]^2
\le {32\over e^2}.
}
\tag{R23.6}
```

Equivalently, if $W=C\cup\{u,v\}$ and
$w_{ab}=e^{\beta c_W(y_{ab})}K_{\beta,W}(y_{ab})$, then
$e^{G_S(y_a)}=\sum_bw_{ab}$ and
$e^{G_{S'}(y_b)}=\sum_aw_{ab}$. Thus adjacent selector values are row and
column log-sums of one positive $2\times2$ table.

The desired potential is $F$, not $G$. For
$s_{ij}=a_{ij}D_{ij}$ and $k=m-1$, $N=n-2$, a uniform
$k$-set $C\subset[n]\setminus\{u,v\}$ satisfies the exact identity

```math
\boxed{
\begin{aligned}
\mathbb E_C(c_{C+u}-c_{C+v})^2
=4\bigg[&{k(N-k)\over N(N-1)}
 \sum_{i\ne u,v}(s_{ui}-s_{vi})^2\\
&+{k(k-1)\over N(N-1)}(H_u-H_v)^2\bigg].
\end{aligned}}
\tag{R23.7}
```

Pair the full Gibbs states under the flip $T_i$. Conditional on a pair with
$|H_i|=a$, its mean field is $a\tanh(2\beta a)$. Since
$a/\tanh(2\beta a)$ is increasing,

```math
\boxed{
\mathbb E_{\nu_\beta}\sum_iH_i^2
\le C_{\beta,n}\mathbb E_{\nu_\beta}\langle A,D\rangle
\le C_{\beta,n}Q(A),qquad
C_{\beta,n}={n-1\over\tanh(2\beta(n-1))}.
}
\tag{R23.8}
```

For $m\ge4$, averaging (R23.7) over $u,v$, using
$\mathbb E_{u\ne v}(H_u-H_v)^2
=2\sum_i(H_i-\bar H)^2/(n-1)$, gives the explicit estimate

```math
\mathbb E(c_S-c_{S'})^2
\le {16k(N-k)\over N-1}
+{8k(k-1)C_{\beta,n}Q(A)\over N(N-1)(n-1)}.
```

At fixed density this is $O(n+Q(A))$ whenever
$\beta n\to\infty$. Combining this with
$F_S-F_{S'}=(G_S-G_{S'})-\beta(c_S-c_{S'})$ yields

```math
\boxed{
\begin{aligned}
\mathbb E(F_S-F_{S'})^2
\le{}&{64\over e^2}+{32\beta^2k(N-k)\over N-1}\\
&+{16\beta^2k(k-1)C_{\beta,n}Q(A)
\over N(N-1)(n-1)}\\
=O\!\left(1+\beta^2[n+Q(A)]\right).
\end{aligned}
}
\tag{R23.9}
```

For an exact minimizer and
$\beta=\Theta(n^{-1/2+c})$, (R23.9) is
$O(n^{1/2+2c})$. The raw bounded-difference estimate is
$O(\beta^2n^2)=O(n^{1+2c})$, so (R23.9) saves a square root. It is still
larger by a factor $n^{1+4c}$ than the
$O(n^{-1/2-2c})$ adjacent scale needed by the Johnson inequality (10.729).
Thus parent selector smoothness does not close the rate bound after the
child-energy term is restored.

#### What a Jensen-gap proof must control

For the two-temperature law in (10.721), put

```math
\ell_S(y)=F_S(y)+(\beta-\gamma)c_S(y),qquad
\mu_t(y)\propto\mu_{\gamma,S}(y)e^{t\ell_S(y)}.
```

Twice differentiating the log-mgf gives the exact identity

```math
\boxed{
\log\mathbb E_{\mu_\gamma}e^{\ell_S}
-\mathbb E_{\mu_\gamma}\ell_S
=\int_0^1(1-t)\operatorname{Var}_{\mu_t}(\ell_S)\,dt,
}
\tag{R23.10}
```

where

```math
\mu_t(y)\propto
e^{[(1-t)\gamma+t\beta]c_S(y)}K_{\beta,S}(y)^t.
```

Therefore variance only at $t=0$ is insufficient; one needs uniform control
through the whole interpolation, or an equivalent exponential-tail bound.

The exact-minimizer energy cap supplies only

```math
(n-m)\log2\le F_S(y)
\le(n-m)\log2+\beta[Q(A)-c_S(y)].
\tag{R23.11}
```

The lower bound is Jensen over uniform outside extensions, and the upper
bound is the pointwise inequality
$\langle A,d\rangle\le Q(A)$. Hence

```math
\operatorname{osc}(\ell_S)
\le2Q(A)[\beta+|\beta-\gamma|].
\tag{R23.12}
```

At $\beta,\gamma=\Theta(n^{-1/2+c})$, this is only
$O(n^{1+c})$, versus the required $O(n^{1/2-2c})$. More importantly, the
upper bound in (R23.11) permits a child state of deficit $\Delta$ to have an
outside weight $e^{\beta\Delta}$ which exactly cancels its child-Gibbs
penalty. This is the rare-tail mechanism which a positive proof must exclude
using information beyond the scalar energy cap.

Indeed, the parent marginal is exactly

```math
\nu_{\beta,S}(y)
=\frac{\mu_{\gamma,S}(y)e^{\ell_S(y)}}
       {\mathbb E_{\mu_{\gamma,S}}e^{\ell_S}}.
```

Thus for the near-ground set $N_{S,t}$,

```math
\boxed{
\Omega_{\beta,S}(t)
=\frac{\mathbb E_{\mu_{\gamma,S}}
 [e^{\ell_S}\mathbf1_{N_{S,t}}]}
 {\mathbb E_{\mu_{\gamma,S}}e^{\ell_S}}.
}
\tag{R23.13}
```

High probability of $N_{S,t}$ under the child Gibbs law is therefore not
enough: it must contain almost all exponentially tilted mass.

#### Exact $A_9$ obstruction

For the hidden-optimal deletion law, the matched-temperature Jensen gap and
the variance at the base child law are:

| $\beta$ | Jensen gap | $\operatorname{Var}_{\mu_\beta}F$ |
|---:|---:|---:|
| 0.05 | 0.00124020 | 0.00239306 |
| 0.10 | 0.01477733 | 0.02639884 |
| 0.25 | 0.18993504 | 0.26812528 |
| 0.50 | 0.66275766 | 0.68662565 |
| 1.00 | 1.33013934 | 0.62872423 |
| 2.00 | 1.57391279 | 0.05868751 |
| 4.00 | 1.58304246 | 0.00008648 |

This has an exact limiting explanation. Every exact child ground has two
extensions and $F_S=\log2$, so
$\operatorname{Var}_{\mu_\beta}F_S\to0$. But only
$k_i=(4,8,4,8,4,4,8,8,8)$ of the 25 positive parent grounds restrict to a
child ground after deleting $i$. Consequently

```math
\lim_{\beta\to\infty}
\left(\log\mathbb E_{\mu_\beta}e^{F_S}
-\mathbb E_{\mu_\beta}F_S\right)
=\log{25\over k_i},
```

and the $\pi_*$-average tends
$1.5830484787467298\ldots$. Thus no universal estimate of the Jensen gap by
$C\operatorname{Var}_{\mu_\gamma}\ell$, or by any function of that variance
which vanishes at zero, can hold. Rare non-ground child states carry the
missing tilted mass.

#### Relation to the regular-ground codebook condition

The overlap identities do **not** imply the row cap required in
`tmp/regular_ground_codebook_r23.md`. For a parent ground $d$, its row fields
$r_i=H_i(d)$ satisfy only $r_i\ge0$ and $\sum_i r_i=q_n$. If
$S=[n]\setminus\{i\}$, then exactly

```math
Q(A[S])-c_S(d)
=2r_i-[q_n-Q(A[S])].
\tag{R23.14}
```

Hard near-ground conditioning therefore bounds the deleted row only. For a
general deleted set $T=[n]\setminus S$,

```math
q_n-c_S(d)=2\sum_{i\in T}r_i-c_T(d),
\tag{R23.15}
```

which is an aggregate identity with an internal signed cancellation and gives
no maximum-row bound. In $A_9$, $Q(A_9-i)=q_9=24$, so exact child-ground
conditioning forces $r_i=0$; nevertheless the selected parent grounds have
maximum other-row field 6 or 8 for deletions 0, 1, 2, 4, 5, 6, and 7, and
can also have maximum 8 for deletions 3 and 8. This is a finite separation,
not an asymptotic counterexample.

The moment estimate (R23.8) also gives only
$\sum_i r_i^2\le(n-1)q_n=O(n^{5/2})$ on a ground state. It cannot force
$\max_i r_i=O(n^{3/4-c})$. A positive overlap-to-codebook implication would
need a new uniform control of the child-state gradients (R23.1), or a
nonmigration theorem for heavy rows; neither follows from hard overlap or
the scalar values of $F_S$.

#### Disposition and falsification criteria

The useful proved gains are the dimension-free selector identity (R23.6),
the low-$Q$ second-moment improvement (R23.9), and the conditional-channel KL
budget (R23.4). None controls the exponential interpolation (R23.10).

A proposed variance proof is falsified by any sequence for which the base
variance is small but the tilted tail remains large; $A_9$ already falsifies
all universal base-variance-only inequalities. More quantitatively, for any
event $B$ and $X=\ell_S-\mathbb E\ell_S$,

```math
\log\mathbb E e^X
\ge\log\mu_{\gamma,S}(B)+\inf_{y\in B}X(y).
```

Thus an $\omega(n^{1/2-2c})$ lower bound on the right, averaged over every
admissible near-uniform selector law and every allowed temperature mixture,
falsifies (10.723) for the proposed $c$. Equivalently for the hard form, one
must prove that every target-specific exact minimizer and admissible selector
has average $-\log\Omega=\omega(n^{1/2-2c})$ at every permitted tolerance.
The positive missing lemma is precisely uniform tilted-tail control in
(R23.10) or (R23.13), not another scalar variance estimate.
