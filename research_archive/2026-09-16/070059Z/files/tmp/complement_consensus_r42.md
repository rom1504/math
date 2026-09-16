# Wave 42: incidence-preserving complement consensus

## Verdict

There is an exact finite convex formulation, but it must be lifted from local
selector kernels to distributions over **global incidence-preserving maps**.
The reason is structural: for a local kernel the independent anchored
conflict is concave, not convex.  The lifted primal is a finite LP whose exact
dual is a three-price incidence-constrained consensus problem.

Shared random priorities provide a concrete family of points in this lifted
polytope.  Their pair law has an exact weighted MinHash formula.  The relevant
quantity is not ordinary intersection/Jaccard codegree, but a **colored
codegree** which records coordinate agreement between exclusive witnesses and
the other fiber.  A second exact cost appears: derandomizing one shared race
pays code length, equal to fractional-cover information plus the conditional
entropy inside each fiber.

Finite exact audits show a substantial MinHash gain on `A5,A6`, but only a
`6--9%` gain on `A8,A9`; row pruning slightly worsens consensus.  These are
mechanism data, not an asymptotic wall.  The precise surviving theorem is a
project-scale weighted colored-codegree bound for actual complement incidence
fibers, simultaneously with row and code-length budgets.

## 1. Actual complement incidence fibers

Let `G` be a finite anchored selector family with law `u`, and let `D` be the
finite set of oriented projective parent cuts, all put in the common anchor
gauge.  For an exact minimizer `A`, define

```math
B_S=-A+2P_SAP_S,
\qquad
\mathcal I_S=\{d\in\mathcal D:\langle B_S,d\rangle\ge q_n\}.
\tag{R42.1}
```

Equation (10.1044) says every `I_S` is nonempty.  No approximate or synthetic
witness is added below.  If desired, replace `I_S` by its intersection with a
row cutoff, provided it remains nonempty.

Let

```math
\mathcal F=\prod_{S\in\mathcal G}\mathcal I_S
```

be the set of deterministic incidence-preserving selections.  For
`f in F`, put

```math
\begin{aligned}
\mathsf R(f)&=\mathbb E_{S\sim u}R_2(f(S)),\\
\mu_f(d)&=\Pr_{S\sim u}\{f(S)=d\},\\
\mathsf H(f)&=H(\mu_f)=I(S;f(S)),\\
\mathsf C(f)&=\mathbb E_{S,T\sim u}
\sum_{i\in(S\cap T)\setminus\{v\}}
\frac{\mathbf1\{y_i^{f(S)}\ne y_i^{f(T)}\}}{p_i},
\qquad p_i=\Pr_u\{i\in S\}.
\end{aligned}
\tag{R42.2}
```

These are exactly the row, deterministic information, and anchored-conflict
costs needed in (10.1089).

## 2. Why the local-channel formulation is not convex

For a local incidence kernel `pi_S(d)`, define

```math
m_i=\mathbb E[y_i^D\mid i\in S].
```

Two independent selector--certificate draws have exactly

```math
\boxed{
\mathcal C_v(\pi)
=\frac12\sum_{i\ne v}p_i(1-m_i^2).}
\tag{R42.3}
```

Thus `C_v` is concave in the channel marginals.  Row is linear and mutual
information is convex, but adding a positive conflict price does **not** give
a convex local-kernel primal.  Pairwise couplings chosen independently for
each selector pair would be an inconsistent relaxation unless they belong to
the global marginal polytope.

The exact convex repair is to put a probability law on whole maps
`f in F`.  Conflict then becomes a linear atom cost.

## 3. Finite convex primal and exact minimax dual

Fix positive budgets `R_*`, `H_*`, and `C_*`, and normalize

```math
r_f=\mathsf R(f)/R_*,\qquad
h_f=\mathsf H(f)/H_*,\qquad
c_f=\mathsf C(f)/C_*.
```

The convexified consensus value is the finite LP

```math
\boxed{
\begin{aligned}
\Theta_*:=\min_{\eta\in\Delta(\mathcal F),\,t}\quad&t\\
\text{subject to}\quad
&\sum_f\eta_fr_f\le t,\\
&\sum_f\eta_fh_f\le t,\\
&\sum_f\eta_fc_f\le t.
\end{aligned}}
\tag{R42.4}
```

Finite-dimensional minimax, equivalently LP duality, gives the exact dual

```math
\boxed{
\Theta_*
=\max_{\lambda\in\Delta_3}
\min_{f\in\prod_S\mathcal I_S}
\left\{\lambda_Rr_f+\lambda_Hh_f+\lambda_Cc_f\right\}.}
\tag{R42.5}
```

This is the separation oracle: for every nonnegative three-price vector,
find an actual-incidence assignment minimizing row, output information, and
pairwise anchored disagreement.  There is no relaxation of the fibers.

The entropy term also has an exact code-length form:

```math
\boxed{
H(\mu_f)=min_{r\in\Delta(\mathcal D)}
\mathbb E_{S\sim u}\log\frac1{r(f(S))}.}
\tag{R42.6}
```

Consequently the inner problem in (R42.5) is exactly

```math
\min_{\substack{f(S)\in\mathcal I_S\\r\in\Delta(\mathcal D)}}
\left\{
\frac{\lambda_R}{R_*}\mathbb E_SR_2(f(S))
+\frac{\lambda_H}{H_*}\mathbb E_S\log\frac1{r(f(S))}
+\frac{\lambda_C}{C_*}\mathsf C(f)
\right\}.
\tag{R42.7}
```

This is a genuine incidence-constrained consensus/codebook problem, not the
exceptional-center event under another name.  The latter chooses one parent
word close to many favorable child fibers; (R42.7) jointly selects actual
complement witnesses and prices their output code and pairwise labels.

If `Theta_*<=1/3`, some atom in an optimal mixture has all three original
costs within budget: the expected sum of normalized costs is at most one.
More generally `Theta_*<=t` yields one atom with every cost at most `3t`
times its budget.  Thus the convexification loses only a fixed constant in
the project estimates.

## 4. Shared weighted priorities: exact pair law

Choose global rates `w_d>=0`, not all zero, with

```math
Z_S=\sum_{d\in\mathcal I_S}w_d>0.
```

Give every response one shared independent exponential clock
`T_d~Exp(w_d)`, and set

```math
f_T(S)=\mathop{\rm argmin}_{d\in\mathcal I_S}T_d.
\tag{R42.8}
```

This is weighted MinHash and is incidence preserving for every realization.
For two fibers `I=I_S`, `J=I_T`, let

```math
A=I\setminus J,\qquad B=J\setminus I,
\qquad U=W(I\cup J),\qquad W(X)=\sum_{d\in X}w_d.
```

Memorylessness gives the complete joint law:

```math
\begin{aligned}
\Pr\{f(S)=f(T)=d\}&=\frac{w_d}{U},
&&d\in I\cap J,\\
\Pr\{f(S)=d,f(T)=e,\ d\text{ is first}\}
&=\frac{w_dw_e}{U Z_T},
&&d\in A,\ e\in J,\\
\Pr\{f(S)=d,f(T)=e,\ e\text{ is first}\}
&=\frac{w_dw_e}{U Z_S},
&&d\in I,\ e\in B.
\end{aligned}
\tag{R42.9}
```

For `delta_i(d,e)=1{y_i^d ne y_i^e}`, the exact expected disagreement at a
common coordinate is therefore

```math
\boxed{
\begin{aligned}
\chi_i(S,T;w)=\frac1U\bigg[&
\frac1{Z_T}\sum_{d\in A,\ e\in J}w_dw_e\delta_i(d,e)\\
&+\frac1{Z_S}\sum_{d\in I,\ e\in B}w_dw_e\delta_i(d,e)
\bigg].
\end{aligned}}
\tag{R42.10}
```

The shared-race conflict is exactly

```math
\boxed{
\mathsf C_{m race}(w)=
\mathbb E_{S,T}\sum_{i\in(S\cap T)\setminus\{v\}}
\frac{\chi_i(S,T;w)}{p_i}.}
\tag{R42.11}
```

Ordinary weighted MinHash gives only

```math
\Pr\{f(S)=f(T)\}=\frac{W(I_S\cap I_T)}{W(I_S\cup I_T)}.
```

Replacing every `delta_i` by one yields the Jaccard bound.  Formula
(R42.10) is strictly sharper: distinct witnesses may still agree on most
overlap coordinates.  It is the required **colored codegree** statistic.

## 5. Exact row, information, and derandomization costs

The one-selector race marginal is

```math
\pi_S(d)=\frac{w_d\mathbf1_{d\in\mathcal I_S}}{Z_S}.
```

Hence

```math
\boxed{
\mathsf R_{m race}(w)=
\mathbb E_S\sum_{d\in\mathcal I_S}\frac{w_d}{Z_S}R_2(d).}
\tag{R42.12}
```

With `W=sum_d w_d` and reference law `r(d)=w_d/W`, the fractional-cover
identity is

```math
\mathbb E_SD(\pi_S\Vert r)=
\mathbb E_S\log\frac W{Z_S}
=I(S;D)+D(P_D\Vert r).
\tag{R42.13}
```

There are two different information quantities here, and they cannot be
interchanged.  Let `J` denote the entire shared clock realization and put
`D=f_J(S)`.  Since `J` is independent of `S`, the chain rule gives exactly

```math
\boxed{
\mathbb E_J I(S;D\mid J)
=I(S;D)+I(S;J\mid D).}
\tag{R42.14}
```

The marginal race channel has the inexpensive information bound (R42.13),
but two independent draws from that channel use independent clocks and have
the independent-kernel conflict (R42.3).  The consensus gain (R42.11) uses
one common `J`; retaining or fixing it makes (10.1087) pay the left side of
(R42.14), including the nonnegative extra term `I(S;J|D)`.  Cheap marginal
information and shared-race conflict therefore cannot be combined for free.

For a **single shared-race realization**, use the code cost

```math
\mathsf L_r(f)=\mathbb E_S\log\frac W{w_{f(S)}}.
```

It always dominates the deterministic information by (R42.6), and its race
expectation is exactly

```math
\boxed{
\begin{aligned}
\mathsf L_{m race}(w)
&=\mathbb E_S\sum_{d\in\mathcal I_S}
\frac{w_d}{Z_S}\log\frac W{w_d}\\
&=\mathbb E_S\log\frac W{Z_S}
+\mathbb E_SH(\pi_S).
\end{aligned}}
\tag{R42.15}
```

Thus fractional-cover weight/information controls only the first term.  The
conditional fiber entropy is the exact additional term in this particular
code-length bound for extracting one global consensus map.  The full code
length is a sufficient upper bound, not an identity with the extra
information in (R42.14):

```math
\boxed{
\mathbb E_J I(S;D\mid J)
\le \mathsf L_{\rm race}(w).}
```

Indeed, for every fixed `J`, (R42.6) bounds the deterministic output entropy
by cross entropy against `r`, and averaging gives the claim.

Equations (R42.11), (R42.12), and (R42.15) give a concrete sufficient lemma:
if some global rates satisfy

```math
\boxed{
\frac{\mathsf R_{m race}(w)}{R_*}
+\frac{\mathsf L_{m race}(w)}{H_*}
+\frac{\mathsf C_{m race}(w)}{C_*}
\le1,
}
\tag{R42.16}
```

then one clock realization is an actual-incidence deterministic selection
satisfying all three clauses of (10.1089).  Constant-factor variants follow
by Markov.  This is the exact project-scale weighted colored-codegree theorem
to prove.  It simultaneously asks for:

1. low-row mass in every actual incidence fiber;
2. a short global response code, including conditional fiber entropy; and
3. small colored exclusive-fiber disagreement (R42.10).

Nonempty fibers, small fractional-cover weight, cheap marginal race
information, or large uncolored pair intersections alone do not establish
(R42.16).

## 6. Finite audit on exact minimizers

The checker uses the full anchored slice and the actual incidence fibers
`<B_S,d>>=q_n`.  It also audits the exact `B_S`-ground fibers.  Priorities are
uniform, so each selector marginal is uniform on its fiber.  Conflict and
MinHash/Jaccard quantities are computed exactly; information values below are
floating evaluations of exact finite sums of logarithms.

| matrix/fiber | fiber sizes | mean row | local `I(S;D)` nats | independent conflict | shared MinHash conflict | Jaccard bound |
|---|---:|---:|---:|---:|---:|---:|
| `A5`, incidence = ground | 10 | 20.80 | 0.9704 | 1.4400 | 0.8593 | 1.7778 |
| `A6`, incidence = ground | 12 | 30.00 | 0.8291 | 2.0000 | 1.2000 | 2.4000 |
| `A8`, incidence | 8--10 | 60.74 | 2.5041 | 2.4331 | 2.2080 | 4.5264 |
| `A8`, exact ground | 2--4 | 64.00 | 2.7145 | 2.3048 | 2.1683 | 4.5820 |
| `A9`, incidence = ground | 5--15 | 89.72 | 2.8035 | 2.6549 | 2.4939 | 5.5981 |

The exact rational triples `(independent conflict, shared conflict, Jaccard
bound)` are

```text
A5:          (36/25, 116/135, 16/9)
A6:          (2, 6/5, 12/5)
A8 incidence:(10949/4500, 83253679/37705500, 1706708/377055)
A8 ground:   (242/105, 683/315, 866/189)
A9:          (141621383/53343360,
              56595209687329/22693598928000,
              93358253/16676660).
```

Decimals are retained in the table only for comparison; the checker carries
all conflict and codegree sums as exact `Fraction` values.  Only the logarithmic
information column is numerical.

Shared priorities reduce conflict by about `40%` on `A5,A6`, but only
`9.3%` on `A8` incidence fibers and `6.1%` on `A9`.  The Jaccard-only bound
is about twice the exact colored value, so uncolored intersection data lose
important label agreement while still failing to predict a strong consensus
gain.

Hard incidence-preserving row pruning shows the tradeoff directly:

```text
A8, row cap 64: mean row 60.74 -> 58.06,
                 MinHash conflict 2.2080 -> 2.2608,
                 information 2.5041 -> 2.6594.

A9, row cap 96: mean row 89.72 -> 81.46,
                 MinHash conflict 2.4939 -> 2.5342,
                 information 2.8035 -> 2.9501.
```

The finite data neither prove nor falsify the asymptotic colored-codegree
theorem.  They do show that choosing low-row witnesses and choosing colliding
labels are distinct objectives, exactly as the three-price dual predicts.

## Verification

`tmp/complement_consensus_r42_check.py` verifies the unweighted specialization
of (R42.9)--(R42.11) against exhaustive priority permutations on a synthetic
fiber pair, then computes every finite table entry from the actual
`A5,A6,A8,A9` incidence relations.  It ends with

```text
PASS complement_consensus_r42_check
```
