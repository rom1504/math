# Two-level uniform deletion and adaptive descent (Wave 12)

## 1. Verdict

There is an exact two-level jackknife calculus, but every **linear** identity
obtained from it is a scalar tautology.  It partitions the fresh-ground
replenishment slack and does not control the terminal excess.

There are two useful positive conclusions.

1. If `A` is an order-`n` minimizer and `T` is uniform of order `m=n-h`,
   then

   ```math
   0\le \mathbb E\bigl[Q(A[T])-q_m\bigr]
   \le h(2n-h-1).
   \tag{D1}
   ```

   Thus, for example, `h=floor(n^{1/3})` gives a macroscopic child with
   mean terminal excess `O(n^{4/3})=o(n^{3/2})`.  This is a genuine uniform
   estimate, but the steps are too short for the errors to be summable along
   a long descent.

2. The all-pairs tail condition (10.500) can be replaced by an exact
   **shortest-path/adaptive tail** condition.  It is enough to certify one
   cheap path of normalized comparisons from a future liminf order to each
   target order.  A vanishing potential-drop bound on the selected edges
   makes the errors telescope, rather than sum.  This is stated precisely in
   Section 6.

Neither variance nor generic set-function curvature supplies the missing
nonlinear remainder.  The certified order-nine minimizer contains an order-8
child with excess `4` for which all eight order-7 children also have excess
`4`; the norms have zero variance at both levels.  In addition,
`S -> Q(A[S])` is neither submodular nor supermodular already by order four.

## 2. Exact two-level jackknife

For a fixed order-`n` signing `A`, define

```math
F_A(k)=\mathbb E_{|S|=k}Q(A[S]),
\qquad
\beta_{r,k}=\frac{(k)_2}{(r)_2}
\quad(2\le k\le r\le n).
\tag{D2}
```

For `r` and `k` define the fresh-ground restriction slack

```math
R^A_{r\to k}=F_A(k)-\beta_{r,k}F_A(r).
\tag{D3}
```

This is nonnegative.  Indeed, condition on a uniform `r`-set `U`, orient a
ground state of `A[U]` positively, and then choose a uniform `k`-set
`S subset U`.  The expected energy of the restricted ground is
`beta_{r,k}Q(A[U])`, while `Q(A[S])` is at least that oriented energy.
Average over `U`.

The falling-factorial coefficients compose:

```math
\beta_{n,\ell}=\beta_{m,\ell}\beta_{n,m}.
```

Consequently, for every `2 <= ell <= m <= n`, there is the exact nonnegative
cocycle

```math
\boxed{
R^A_{n\to\ell}
=R^A_{m\to\ell}
+\beta_{m,\ell}R^A_{n\to m}.
}
\tag{D4}
```

Here `R^A_{m->ell}` means the average of the local `m->ell` slack over the
uniform `m`-set.  Equation (D4) is just

```math
F_A(\ell)-\beta_{n,\ell}F_A(n)
=F_A(\ell)-\beta_{m,\ell}F_A(m)
+\beta_{m,\ell}
  [F_A(m)-\beta_{n,m}F_A(n)].
```

Now put

```math
\varepsilon_k^A=F_A(k)-q_k.
\tag{D5}
```

Substitution into (D3) gives

```math
\boxed{
R^A_{m\to\ell}
-\varepsilon_\ell^A
+\beta_{m,\ell}\varepsilon_m^A
=q_\ell-\beta_{m,\ell}q_m.
}
\tag{D6}
```

Thus all linear two-level combinations eliminate to scalar minima and
terminal excesses.  There is no second independent equation.

There is an equivalent replenishment statement.  Compare direct restriction
from the original order-`n` ground with restriction after choosing a fresh
ground on the intermediate `m`-set.  If their mean replenishments at order
`ell` are respectively `G_dir` and `G_fresh`, then

```math
\boxed{
G_{\rm dir}-G_{\rm fresh}
=\beta_{m,\ell}R^A_{n\to m}.
}
\tag{D7}
```

So intermediate reoptimization merely transfers the nonnegative first-stage
slack between the two levels.  It does not remove it.

## 3. A uniform near-top terminal-excess estimate

Let `A` be an exact order-`n` minimizer.  Principal monotonicity gives
`Q(A[T]) <= Q(A)=q_n` for every `T`: extend a ground state on `T` by
independent mean-zero signs and choose a completion with at least the same
absolute oriented energy.

Conversely, start with an order-`m` minimizer and fill all edges incident to
the `n-m` new vertices arbitrarily.  The doubled contribution of those edges
has absolute value at most

```math
2\left[\binom n2-\binom m2\right]
=n(n-1)-m(m-1).
```

Hence

```math
q_n-q_m\le n(n-1)-m(m-1).
\tag{D8}
```

Combining the two observations proves (D1).  For `m=n-h`, its normalized
form is

```math
\frac{\mathbb E[Q(A[T])-q_m]}{m^{3/2}}
\le
\frac{h(2n-h-1)}{(n-h)^{3/2}}.
\tag{D9}
```

It is `o(1)` whenever `h=o(sqrt(n))`.  It is not a descent theorem: if steps
of size `h(k)` are repeated, the edge-count bound has cumulative scale

```math
\sum \frac{h(k)}{\sqrt{k}}
\asymp \int \frac{dk}{\sqrt{k}},
\tag{D10}
```

because the number of steps per unit decrease is `1/h(k)`.  The block size
cancels, and the worst-case tail diverges.  Thus (D1) cannot replace
(10.500).

## 4. Exact finite obstruction to variance corrections

Use the certified order-nine minimizer in (10.298), with doubled norm
`q_9=24`.  Exhaustive enumeration of all principal restrictions gives

| order `m` | `F_A(m)` | `q_m` | mean excess |
|---:|---:|---:|---:|
| 4 | `184/21` | `8` | `16/21` |
| 5 | `248/21` | `8` | `80/21` |
| 6 | `118/7` | `10` | `48/7` |
| 7 | `21` | `18` | `3` |
| 8 | `24` | `20` | `4` |
| 9 | `24` | `24` | `0` |

At order eight every one of the nine restrictions has norm `24`, so its
variance and range are both zero despite mean excess `4`.

More strongly, delete zero-based vertex `7` from the displayed order-nine
matrix and call the resulting order-eight signing `B`.  Then

```math
\boxed{
Q(B)=24=q_8+4,
\qquad
Q(B[-i])=22=q_7+4
\quad(i\in V(B)).
}
\tag{D11}
```

Thus the excess is exactly constant across the whole `8 -> 7` jackknife and
both the parent norm and all child norms have zero variance.  Numerically in
the exact cocycle,

```math
R^B_{8\to7}=22-\frac34\,24=4,
```

and (D6) reads

```math
4-4+\frac34\,4
=18-\frac34\,20
=3.
```

This defeats:

- contraction of terminal excess with no creation term (the exact parent has
  excess zero and all order-eight children have excess four);
- strict contraction on near-minimizer descendants (the excess remains four
  from `B` to every child);
- any correction whose only nonlinear input vanishes with the variance or
  range of first- and second-level principal norms.

For completeness, the exact law-of-total-variance decomposition contains no
hidden sign.  If `U` is uniform of order `m`, `T` is uniform of order `ell`
inside `U`, and

```math
\mu(U)=\mathbb E[Q(A[T])\mid U]
=\beta_{m,\ell}Q(A[U])+R_U,
```

then

```math
\begin{aligned}
\operatorname{Var}Q(A[T])
={}&\mathbb E\operatorname{Var}(Q(A[T])\mid U)
+\beta_{m,\ell}^2\operatorname{Var}Q(A[U])\\
&+\operatorname{Var}R_U
+2\beta_{m,\ell}\operatorname{Cov}(Q(A[U]),R_U).
\end{aligned}
\tag{D12}
```

The covariance has no universal sign, and (D11) already makes all the most
obvious dispersion corrections vanish while the excess persists.

The other standard nonlinear route also fails generically.  Exhaustive
enumeration gives:

- submodularity fails already at order three (two distinct singleton sets
  have defect `-2`);
- supermodularity fails at order four for edge signs
  `(-1,-1,-1,-1,-1,+1)` and subset masks `7,11`, with defect `+2`.

Thus `S -> Q(A[S])` has neither curvature sign needed for a generic nonlinear
jackknife theorem.

## 5. Exact checker

The script `tmp/audit_adaptive_deletion_r12.py` independently evaluates all
spins on every principal restriction of the order-nine matrix.  It verifies
the table, (D4) for every `2 <= ell < m <= 9`, and the flat child (D11).
The separate script `tmp/check_q_setmod_r12.py` verifies the submodularity and
supermodularity failures.  Both use exact integer or rational arithmetic.

## 6. Adaptive shortest-path bridge (strictly weaker than (10.500))

This is the clean surviving theorem.  Put

```math
a_n=\frac{q_n}{n^{3/2}}.
```

Suppose a directed edge `n -> m`, where
`ceil(n/2) <= m < n`, records a proved comparison

```math
a_m\le a_n+\eta_{n,m},
\qquad \eta_{n,m}\ge0.
\tag{D13}
```

For `M >= N`, let `d(M,N)` be the infimum of `sum eta` over all directed
paths

```math
M=n_0>n_1>\cdots>n_L=N;
```

put `d(M,N)=infinity` if no such path exists.

### Adaptive bridge theorem

If

```math
\boxed{
\Omega_{\rm ad}(N)
:=\sup_{M\ge N}d(M,N)
\longrightarrow0,
}
\tag{D14}
```

then `a_n` converges.

Indeed, let `M_j -> infinity` be a subsequence on which
`a_{M_j} -> L=liminf a_n`.  Telescope (D13) along a cheapest path from
`M_j` to `N`:

```math
a_N\le a_{M_j}+d(M_j,N)
\le a_{M_j}+\Omega_{\rm ad}(N).
```

Let `j -> infinity` and then `N -> infinity`.  This gives
`limsup a_N <= L`, hence convergence.

Condition (D14) is strictly weaker than (10.500): (10.500) controls one
prescribed canonical path using an error bound on every allowed pair, whereas
(D14) optimizes over whatever sparse set of comparisons has actually been
proved.  The logically weakest version only asks, for one fixed liminf
subsequence `M_j`, that

```math
\liminf_{j\to\infty}d(M_j,N)\le\omega(N),
\qquad \omega(N)\to0.
\tag{D15}
```

The same proof applies.

### Potential-drop corollary

A particularly checkable sufficient condition avoids summing independent
errors.  Let `V(n)>=0` with `V(n)->0`.  If for every target `N` and every
`M>=N` one can select a path to `N` such that every edge satisfies

```math
\boxed{
\eta_{n,m}\le V(m)-V(n),
}
\tag{D16}
```

then

```math
\sum_{(n,m)\text{ on the path}}\eta_{n,m}
\le V(N)-V(M)\le V(N),
\tag{D17}
```

so (D14) holds.  If, away from the final landing window, the selected child
also satisfies `m <= theta n` for a fixed `theta<1`, the path has
`O(log(M/N))` steps and cannot stall.  Only the selected edges need (D16).

There is an analogous terminal-excess lemma for nonoptimal intermediate
restrictions.  Write

```math
e(B)=\frac{Q(B)-q_{|B|}}{|B|^{3/2}}.
```

If every encountered order-`r` restriction `B` admits an induced child `C`
of deterministic order `m` conditional on `B`, with

```math
\boxed{
\mathbb E[e(C)\mid B]
\le e(B)+V(m)-V(r),
}
\tag{D18}
```

then `e(B_t)-V(|B_t|)` is a supermartingale.  Starting from an exact
minimizer and stopping at a bounded landing time `L` gives
`E e(B_L)\le E V(|B_L|)-V(n)`, not a sum of local excesses.  When the
terminal order is the deterministic value `N`, this becomes
`E e(B_L)\le V(N)-V(n)`.  For a random terminal order the expectation on
`V(|B_L|)` must be retained, or bounded over a declared landing window.
This controls terminal excess only; it does not by itself imply (D13),
(D14), or convergence.  It is the weakest useful
replacement target suggested by the jackknife audit: prove a one-sided
potential drop for one adaptively selected child, with a landing rule.  A
variance-only term that vanishes on flat profiles cannot serve as that
potential because of (D11).

## 7. Equation-ready conclusions

1. **Verified two-level no-go:** (D4), (D6), and (D7).  Linear uniform
   deletion/jackknife identities are exact cocycles and scalar tautologies.
2. **Verified partial positive result:** (D1)/(D9) gives `o(n^{3/2})` mean
   terminal excess for a near-top child `m=n-o(sqrt(n))`, but (D10) shows why
   this does not sum along descent.
3. **Verified nonlinear obstruction:** the exact plateau (D11), plus failure
   of both submodularity signs, rules out generic variance/range/curvature
   repairs.
4. **Verified weaker descent theorem:** replace the all-pairs canonical tail
   by the adaptive shortest-path tail (D14), or prove the stronger local
   potential certificate (D16).  For terminal excess on inherited
   near-minimizers, the corresponding precise open lemma is (D18).
