# Complemented half-stopping gives an unconditional `K <= 4`

This memo uses only the verified endpoint identities and the verified
path-cover recursion in the ledger.  It does **not** use the false inequality
`B(U) >= I(U)` from (10.435).

## 1. Notation at one endpoint node

Let `U` be a nontrivial node of a fixed endpoint partition tree, with
endpoint shores `S,T`.  Write

- `R(U)=P(U)+N(U)`;
- `I(U)=|P(U)-N(U)|`;
- `I(X)=|P(X)-N(X)|` for `X=S,T`;
- `D(U)=I(S)+I(T)`;
- `h_S,h_T` for the inherited internal endpoint energies;
- `b_X=max_sigma (mu_X^sigma)_+` and `B(U)=b_S+b_T`.

The exact exposure identity (10.432), and then (10.433), give

```math
b_X=\frac{R(U)}2-Q(X)+|h_X|.
```

Consequently

```math
B(U)=R(U)-Q(S)-Q(T)+|h_S|+|h_T|.                 \tag{A}
```

For every block `X`, purely from the definitions,

```math
Q(X)=\frac{R(X)+I(X)}2.                          \tag{B}
```

Separate `P/N` block-superadditivity gives

```math
R(S)+R(T)\le R(U).                               \tag{C}
```

Finally (10.432) gives

```math
I(U)=2|h_S+h_T|,
```

so `|h_S|+|h_T| >= I(U)/2`.  Substituting (B)--(C) in (A) gives the
complement inequality

```math
\begin{aligned}
B(U)
&=R(U)-\frac{R(S)+R(T)+D(U)}2+|h_S|+|h_T|\\
&\ge \frac{R(U)+I(U)-D(U)}2.
\end{aligned}                                    \tag{D}
```

Since `R(U) >= I(U)`, this implies

```math
\boxed{B(U)+\frac{D(U)}2\ge I(U).}               \tag{E}
```

This is stronger than merely combining (10.407) with (10.434): it uses the
exact formula (10.433), the identity (B), and range superadditivity.  It also
does not imply the falsified `B(U)>=I(U)`; the missing local capacity in that
false inequality is supplied by half-sized recursive imbalance reservoirs.

## 2. Continuous reservoir lemma

**Lemma.**  Fix any endpoint tree below `U`.  Let `0<=t<=1` and let `U`
receive an obligation

```math
0\le w_U\le t I(U).
```

There is a conserved allocation below `U` whose path-cover value satisfies

```math
\theta_U\le2t.                                   \tag{F}
```

**Proof.**  At `U` make four nonnegative reservoirs:

```math
t b_S,\qquad t b_T,\qquad \frac t2 I(S),\qquad
\frac t2 I(T).                                   \tag{G}
```

The first two are allocation reservoirs: on the edge toward `X`, use only
one orientation attaining `b_X`.  The last two are child-obligation
reservoirs.  By (E), their total is at least `t I(U)`, hence at least `w_U`.
Because these are continuous fractional reservoirs, choose nonnegative
amounts bounded by (G) whose sum is exactly `w_U`.  Denote the two child
amounts by `w_S,w_T`.  Then

```math
w_X\le\frac t2 I(X).                              \tag{H}
```

Recurse in `X` with parameter `t/2`.  A singleton has `I=0`, so its incoming
obligation is zero and is the induction base.  The tree is finite because
every nonzero-range nontrivial endpoint cut strictly reduces the vertex set.

On the edge from `U` toward `X`, the local directed-edge load is at most `t`:
only the selected best bucket is used, and it receives at most `t b_X`.
If `b_X=0`, its reservoir and allocation are both zero, so the load is zero
under the convention in (10.411).  By induction and (H),
`theta_X<=2(t/2)=t`.  Therefore the exact recursion (10.412) gives

```math
\theta_U
=\sum_{X=S,T}\max\{\ell_{UX},\theta_X\}
\le t+t=2t.
```

This proves (F), including zero capacities and zero-imbalance children.

## 3. Root construction

At the root `A`, let

```math
C_0=\sum_{X=S,T}\sum_{\sigma=\pm}(\mu_X^\sigma)_+,
\qquad D_0=I(S)+I(T).
```

The verified coverage inequality (10.407) says

```math
C_0+D_0\ge R(A).
```

For `R(A)>0`, set

```math
t_0=\frac{R(A)}{C_0+D_0}\le1.                   \tag{I}
```

(If `R(A)=0`, take the zero allocation.)  At the root allocate exactly
`t_0` times every nonzero directed capacity and pass the two obligations

```math
w_X=t_0 I(X).
```

Their sum is exactly `t_0(C_0+D_0)=R(A)`, so root conservation holds.
On either root edge there are at most two nonzero orientation buckets, hence

```math
\ell_{AX}\le2t_0.                                \tag{J}
```

Apply the reservoir lemma to each child with parameter `t_0`; it gives
`theta_X<=2t_0`.  One final use of (10.412) yields

```math
\boxed{
K=\theta_A
\le\sum_{X=S,T}\max\{2t_0,2t_0\}
=4t_0\le4.}                                      \tag{K}
```

Thus every fixed endpoint partition tree admits a conserved allocation of
the full root range with path-cover mass at most four.  No endpoint-tie
optimization is required.

## 4. Harvesting consequence and scope

By the exact path-cover theorem (10.412), the allocation in (K) is dominated
coordinatewise by a finite measure on oriented root-to-leaf chains of total
mass at most four.  Applying the verified all-successor bound (10.341) to
each chain and integrating gives the same implication stated after (10.412):

```math
3KQ(A)\le12Q(A).
```

Therefore this supplies the previously missing universal `12Q` path-cover
harvesting bound within the ledger's reduction.  It does not rehabilitate
(10.435), which remains false, and it makes no claim that the optimized
constant four is sharp.  The exact lower witness (10.446) still shows that
this framework cannot have a universal constant below `10/3`.

## 5. Exact finite checks

`tmp/half_stopping_k4_r10.py` implements the construction with rational
arithmetic, checking all capacities, conservation equations, complement
slacks, and the recursion (10.412).  On the named fixed trees it returns:

| tree | `t_0` | constructed `theta_A` |
|---|---:|---:|
| positive triangle | `2/3` | `2` |
| matrix (10.420) | `8/11` | `32/11` |
| `A_7` (10.443) | `5/7` | `20/7` |
| `A_8` (10.445) | `5/6` | `10/3` |
| reset matrix (10.403), lexicographic tree | `7/8` | `7/2` |
| rank-two family (10.416), `m=4` | `2/3` | `8/3` |

The `A_8` construction coincides with its exact optimum `10/3`.  On the
strict weighted counterexample (10.448), `D=I`, so (E) has ample positive
slack even though `B<I`; this identifies exactly why that counterexample
does not obstruct the recursive argument.
