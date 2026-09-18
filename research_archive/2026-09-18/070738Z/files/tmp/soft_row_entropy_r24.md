# Wave 24 root memo: soft row-square/entropy interface

## Status

The two inequalities below are proved.  They strictly relax the hard support
cap in (10.741) to an output-weighted row-square condition.  They do not
construct the required channel, so they do not prove convergence.

Keep the notation of (10.690)--(10.742).  For a full oriented cut `d`, let

```math
r_i(d)=\sum_{j\ne i}\sigma a_{ij}x_ix_j,
\qquad R_2(d)=\sum_i r_i(d)^2,
```

and let `Lambda_d(lambda)` be the uniform-slice centered log-mgf from
(10.691).

## 1. A domain-free linear bound

For independent Bernoulli-`p` selectors, the exact decomposition is

```math
T_\xi-p^2W=L+Q,
\qquad L=p\sum_i r_i(\xi_i-p).
```

Hoeffding's lemma uses only that every centered Bernoulli variable has range
length one and gives, for every real `u`,

```math
\log\mathbb E e^{uL}\le \frac{u^2p^2R_2}{8}.
```

Use the same random-bipartition/Gaussian argument as in Wave 21 for `Q`, and
then Cauchy--Schwarz with parameter `2t`.  This yields

```math
\boxed{
\log\mathbb E e^{t(T_\xi-p^2W)}
\le \frac{t^2p^2R_2}{4}
+\frac{t^2F_\times(\mathsf W)}
{4(1-t^2\|A\|_{\rm op}^2)}.
}
```

The only domain condition is `|t|<1/||A||_op`; there is no
`R_infinity` denominator.  Conditioning on the slice and putting
`t=2 lambda` proves

```math
\boxed{
\Lambda_d(\lambda)
\le \chi_{n,m}
+\lambda\epsilon_{n,m}|\langle A,d\rangle|
+\lambda^2p^2R_2(d)
+\frac{\lambda^2F_\times(\mathsf W_d)}
{1-4\lambda^2\|A\|_{\rm op}^2}.
}
```

This is weaker than Bernstein for one regular cut, but it can be averaged
over cuts with no common row-maximum assumption.

## 2. Exact output-adapted reference optimization

Let `S` be uniform on the `m`-slice and let `D=f(S)` be any deterministic
full-cut channel with output law `p_d`.  Put

```math
Z(S,d)=c_A(S,d)-p_2\langle A,d\rangle.
```

For any full-support reference law `nu` on the output alphabet, entropy
duality gives

```math
\lambda\mathbb E Z(S,D)
\le
\mathbb E\log\frac1{\nu(D)}
+\log\sum_d\nu(d)e^{\Lambda_d(\lambda)}.
```

This reference can be optimized exactly.  If
`w_d=nu(d)e^{Lambda_d}/sum_e nu(e)e^{Lambda_e}`, the right side is

```math
\mathbb E[\Lambda_D-\log w_D].
```

It is minimized by `w=p`, so

```math
\boxed{
\inf_\nu\left\{
\mathbb E\log\frac1{\nu(D)}
+\log\mathbb E_\nu e^{\Lambda_d}
\right\}
=H(D)+\mathbb E\Lambda_D.
}
```

An optimizer is proportional to
`p_d exp(-Lambda_d)`.  Thus using `log|C|+sup_C Lambda` is not intrinsic;
the exact cost of a deterministic codebook assignment is output entropy plus
the output-weighted log-mgf.

If the channel has average selector loss at most `delta`, then
`V_ad<=E Q(A[S])<=E c_A(S,D)+delta`, while
`<A,d><=q_n`.  Hence

```math
\boxed{
V_{\rm ad}(A,m)-p_2q_n
\le \delta+
\frac{H(D)+\mathbb E\Lambda_D(\lambda)}{\lambda}.
}
```

The same calculation works for a stochastic channel and an arbitrary selector
marginal `pi`.  If `p_D` is its output law, then

```math
D(P_{SD}\Vert U_m\otimes\nu)
=I(S;D)+D(\pi\Vert U_m)+D(p_D\Vert\nu),
```

and optimizing `nu` gives

```math
\boxed{
V_{\rm ad}(A,m)-p_2q_n
\le \delta+
\frac{I(S;D)+D(\pi\Vert U_m)
+\mathbb E_{p_D}\Lambda_D(\lambda)}{\lambda}.
}
```

Indeed, after the same substitution by `w`, the reference-dependent part is
`D(p_D||w)+E Lambda_D`, minimized at `w=p_D`.  Thus entropy in the
deterministic statement is exactly mutual information, not a proof artifact.

On the halved spectral domain
`lambda<=1/(2 sqrt(2)||A||_op)`, average the domain-free slice bound and put

```math
J=I(S;D)+D(\pi\Vert U_m)+\chi_{n,m},
\qquad
V_{\rm soft}=p^2\mathbb E R_2(D)+\frac{n^2}{2}.
```

Then the exact soft analogue of (10.693) is

```math
\boxed{
V_{\rm ad}(A,m)-p_2q_n
\le \delta+\epsilon_{n,m}q_n
+\frac{J}{\lambda}+\lambda V_{\rm soft}.
}
```

If `sqrt(J/V_soft)` lies in the domain, optimization gives
`delta+epsilon q+2 sqrt(J V_soft)`.  Otherwise a fixed sufficiently small
multiple of `n^{-3/4}` gives the same power scales under the hypotheses below.

## 3. Soft row-square sufficient lemma

At fixed selector ratios, suppose a target-specific exact minimizer admits a
possibly stochastic channel taking values in its parent grounds and, for fixed
`0<c<1/4`,

```math
\boxed{
\begin{aligned}
\mathbb E\ell(S,D)&=O(n^{3/2-c}),\\
I(S;D)+D(\pi\Vert U_m)&=O(n^{3/4-c}),\\
\mathbb E R_2(D)&=O(n^{9/4-c}).
\end{aligned}
}
```

Take `lambda=eta n^{-3/4}` with a sufficiently small fixed `eta`.  The exact
minimizer bounds `||A||_op^2<=2q_n=O(n^{3/2})`,
`F_x<=n^2/4`, and `<A,D>=q_n` give

```math
\mathbb E\Lambda_D
=O(\log n)+O(n^{-1/4})
+O(n^{3/4-c})+O(n^{1/2}).
```

Consequently

```math
V_{\rm ad}(A,m)-p_2q_n=O(n^{3/2-c}),
```

which is the same power-saving restriction edge as (10.741)--(10.742).

The old hypothesis implies this one: if every output ground has
`R_infinity<=O(n^{3/4-c})`, then parent-ground positivity gives
`R_2<=R_infinity q_n=O(n^{9/4-c})`, and a deterministic codebook has
`I=H(D)<=log|C|`.  The converse is false for abstract channels: a small output
mass may use spiky grounds while the mean row-square remains subcritical.

## 4. Exact joint-overlap target

This has a particularly clean global common-prior form.  Let `nu` be a law on
parent grounds, put

```math
G=\{(S,d):\ell(S,d)\le t\},
\qquad Z=(U_m\otimes\nu)(G),
```

and condition the single product law `U_m x nu` on `G`.  The resulting joint
law has distortion at most `t` and

```math
D(P_{SD}\Vert U_m\otimes\nu)=-\log Z.
```

By the chain rule,
`I(S;D)+D(pi||U_m)<=-log Z`; the omitted term is
`D(p_D||nu)>=0`.  Its row-square cost is exactly

```math
\mathbb E_P R_2(D)
=\frac{\mathbb E_{U_m\otimes\nu}
[R_2(D)\mathbf1_G]}{Z}.
```

It is therefore enough to prove, for `t=O(n^{3/2-c})`, one parent-ground prior
satisfying the exact joint average conditions

```math
\boxed{
\begin{aligned}
-\log Z&=O(n^{3/4-c}),\\
\mathbb E[R_2(D)\mid G]&=O(n^{9/4-c}).
\end{aligned}
}
```

Equivalently, the two nontrivial costs enter through the product
`(-log Z+chi)(p^2 E[R_2|G]+n^2)` in the optimized inequality; the displayed
separate exponents are a robust sufficient allocation of the required
`O(n^{3-2c})` product budget.

This is a genuinely joint condition: regularity is measured after conditioning
a ground and selector to be compatible.  It asks only for average overlap,
not the capped lower selector tail in (10.739), because `V_ad` is a minimum and
the selector KL already charges concentration.  The `A_9` fact that the most
regular grounds miss some selectors therefore does not obstruct the averaged
pair.  A delta prior recovers the still weaker-looking one-cut target: one
ground with sufficiently large uniform coverage and subcritical `R_2`.

For the uniform law on the 25 oriented-projective `A_9` parent grounds and
uniform one-deletion selectors, the exact global overlaps at ledger tolerances
`t=0,4,8` are respectively `56/225`, `132/225`, and `194/225`.  The
corresponding conditional row-square means are `752/7`, `3472/33`, and
`10032/97`.  Thus the soft joint statistic remains benign in the finite
example where the hardest row cap and coverage separate.  This is an audit,
not asymptotic evidence.

## 5. Exact next combinatorial target

The hard-cap learner can therefore be replaced by a softer joint problem:
construct an average-distortion channel while controlling the entropy and
the output-weighted second row moment.  Row regularity and selector coverage
still must be proved jointly, but they need not hold word by word.  A useful
weak-learning theorem may charge a spiky learner by the selector mass it
actually captures, rather than exclude it from the action set.

This does not follow from the universal ceiling
`R_2<=2(n-1)q_n=O(n^{5/2})`; that ceiling loses `n^{1/4+c}`.  Nor does it
remove the Wave 23 converse: a leading restriction gap can still force the
needed channel complexity to the critical scale.  The gain is precisely that
finite examples such as `A_9`, where the most regular grounds miss selectors,
no longer force a hard separation.
