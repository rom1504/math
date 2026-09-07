# Wave 17: mixed `A_9` products and a growing-Hadamard no-go

Status labels in this memo are literal.  The algebraic identities and the
growing-Hadamard theorem are **Verified**.  The only finite input to the latter
is a complete integer enumeration of `512*10*15=76,800` cross-Gram entries,
independently strengthened by exhaustion of all order-eighteen Boolean states.
The Paley-conference optimization data at the end are **Numerical** only.

## 1. Classification of separable Boolean products

Let `A=A_9`, let `R` be any symmetric zero-diagonal order-`m` signing, and
choose arbitrary diagonal sign matrices

```math
D=\operatorname{diag}(d_1,\ldots,d_9),\qquad
E=\operatorname{diag}(e_1,\ldots,e_m).
```

Put `P=A+D` and `H=R+E`.  The natural separable Kronecker construction is

```math
\boxed{
\mathcal T(R,E;D)
=H\otimes P-\operatorname{diag}(H\otimes P)
=R\otimes P+E\otimes A.
}
\tag{1}
```

It is a symmetric zero-diagonal `\{\pm1\}` signing of order `9m`.  Conversely,
among genuinely separable rules `t_{(a,i),(b,j)}=h_{ab}p_{ij}`, every rule
whose entry at distinct seed and factor coordinates is `r_{ab}a_{ij}` has
this form: the only missing signs are those at `a=b` and `i=j`, and they are
exactly the two diagonal completions `E,D`.
This is different from the clone and lexicographic products already excluded
in (10.583)--(10.584).

For a Boolean state written as fibres `x_a\in\{\pm1\}^9`, (1) gives the exact
finite-state formula

```math
x^{\mathsf T}\mathcal T x
=\sum_{a\ne b}r_{ab}x_a^{\mathsf T}Px_b
{}+\sum_ae_a x_a^{\mathsf T}Ax_a.
\tag{2}
```

Choose 256 projective representatives `u_t` of `\{\pm1\}^9/\{\pm1\}` and
write `x_a=\sigma_a u_{t_a}`.  With

```math
K_D(t,u)=u_t^{\mathsf T}Pu_u,\qquad
L(t)=u_t^{\mathsf T}Au_t,
```

equation (2) is an exact signed 256-state model:

```math
\boxed{
Q(\mathcal T)=
\max_{\sigma,t}\left|
\sum_{a\ne b}r_{ab}\sigma_a\sigma_bK_D(t_a,t_b)
{}+\sum_ae_aL(t_a)
\right|.
}
\tag{3}
```

The second sum is at most `24m`.  Thus for a growing exact conference factor,
the leading problem is precisely the old finite-channel conference problem
with channel matrix `P/3`; indeed
`\operatorname{diag}((P/3)^2)=\mathbf1`.  The existing finite-type lower
theorem/conjecture at constant `c_2=0.7833875\ldots` does not by itself prove
nonminimality, which here requires a constant strictly above the Paley upper
constant one.

There is also an exact product-state bound.  Write
`\delta=\operatorname{tr}D`, `\varepsilon=\operatorname{tr}E`.  For seed and
factor states `u,v`, respectively,

```math
(v\otimes u)^{\mathsf T}\mathcal T(v\otimes u)
=(u^{\mathsf T}Au+\delta)(v^{\mathsf T}Rv+\varepsilon)
-\delta\varepsilon.
\tag{4}
```

Because `A_9` has both `+24` and `-24` grounds, (4) implies

```math
\boxed{
Q(\mathcal T)\ge(24+|\delta|)Q(R)-24|\varepsilon|.
}
\tag{5}
```

On the arithmetic-resonant Paley subsequence from Section 1.12, where
`Q(R)/(m\sqrt m)\to1`, this alone excludes every completion with
`|\delta|\ge5`, since `24+|\delta|>27`.  It does not exclude the balanced
cases `|\delta|=1,3`.

## 2. What happens to the `3+6` wall

Write the seed partition as

```math
A=C_0+(A_3\oplus A_6),
```

where `C_0` is the seed cross-only matrix.  Under the enlarged partition
`[m]\times3\;\sqcup\;[m]\times6`, equation (1) splits exactly as

```math
\begin{aligned}
\mathcal C_m&=(R+E)\otimes C_0,\\
\mathcal D_{3,m}&=R\otimes(A_3+D_3)+E\otimes A_3,\\
\mathcal D_{6,m}&=R\otimes(A_6+D_6)+E\otimes A_6.
\end{aligned}
\tag{6}
```

If one only replaces `A_3,A_6` separately inside each of the `m` intact seed
fibres, the perturbation is block diagonal and

```math
\begin{aligned}
Q(\text{perturbation})
&\le m\{Q(G_3-A_3)+Q(G_6-A_6)\}\\
&\le m\{q_3+Q(A_3)+q_6+Q(A_6)\}=36m.
\end{aligned}
\tag{7}
```

Here `(Q(A_3),Q(A_6))=(6,14)` and `(q_3,q_6)=(6,10)`.  In particular the
four-unit seed internal excess becomes only the within-fibre replicated excess
`4m=o((9m)^{3/2})`.  This is the
mixed-product specialization of (10.582): arbitrary inter-fibre conference
coupling does not make a fibrewise seed selector leading.

To obtain a leading replacement effect one must also change the
`R\otimes(A_r+D_r)` terms, hence `\Theta(m^2)` edges inside the enlarged
shores.  Replacing these by the analogous products built from seed minimizers
is admissible in the response problem only if those order-`3m` and order-`6m`
products are themselves global minimizers.  Membership
`G_r\in\mathcal M_r` gives no such conclusion.  This is the exact missing
tensor-compatibility hypothesis; the six-state order-nine certificate cannot
simply be tensored.

There is nevertheless a genuine leading **internal-excess diagnostic** for
conference products.  The displayed `A_3` is an unbalanced triangle with
energies `-6` and `2`.  If `\delta_3=\operatorname{tr}D_3<0`, (4) gives

```math
Q(\mathcal D_{3,m})
\ge(6+|\delta_3|)Q(R)-6|\varepsilon|.
\tag{8}
```

For any exact symmetric-conference sequence, the verified depth-two theorem
gives `Q(R)\ge(c_2-o(1))m\sqrt m`.  Since

```math
7c_2=5.4837127\ldots>3\sqrt3=5.1961524\ldots,
```

even `\delta_3=-1` makes (8) exceed the asymptotic Paley upper bound for
`q_{3m}` by a positive multiple of `m^{3/2}`.  Thus these candidates really
do carry leading captured excess in their `3m` shore.  This does **not** say
the order-`9m` product is a minimizer, and it does not provide compatible
order-`3m` and order-`6m` minimizing replacements.

## 3. A rigorous growing-Hadamard no-go

The Hadamard-like branch can be closed for a natural infinite family and for
**every** seed diagonal completion.

Let

```math
H_2=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad P=A_9+D,
```

and form the order-eighteen signing

```math
B_D=H_2\otimes P-\operatorname{diag}(H_2\otimes P).
\tag{9}
```

The `A_9` positive and negative projective ground sets have respectively 10
and 15 elements.  Exact enumeration of their 150 cross pairs for all 512
choices of `D` proves

```math
\boxed{
\min_D\ \max_{p^{\mathsf T}A_9p=24,\ n^{\mathsf T}A_9n=-24}
|p^{\mathsf T}(A_9+D)n|=15.
}
\tag{10}
```

The histogram of the inner maximum over the 512 completions is

```text
15: 7, 17: 82, 19: 204, 21: 165, 23: 48, 25: 6.
```

For any pair in (10), the two states `(p,n)` and `(p,-n)` in (9) have
energies

```math
48\pm2p^{\mathsf T}Pn.
```

Consequently

```math
\boxed{Q(B_D)\ge78\quad\text{for every }D.}
\tag{11}
```

Full exhaustion of all `2^{17}` projective order-eighteen states independently
checks (11); the exact minimum over `D` is 78.  Already
`78>18\sqrt{17}`, so every order-two mixed product is nonminimal.

Now put `s=4^r`.  Let

```math
H_s=H_4^{\otimes r},\qquad H_4=H_2\otimes H_2,
\qquad w=(1,1,1,-1).
```

Since `H_4w=2w`, the Boolean vector `y=w^{\otimes r}` satisfies

```math
H_sy=\sqrt s\,y,qquad y^{\mathsf T}H_sy=s\sqrt s.
\tag{12}
```

Take the growing full Hadamard factor

```math
H_m=H_s\otimes H_2,qquad m=2s,
```

its zero-diagonal signing `R_m=H_m-\operatorname{diag}H_m`, and precisely that
Hadamard diagonal as `E_m`.  The mixed signing is

```math
\mathcal T_m=H_m\otimes P-\operatorname{diag}(H_m\otimes P),
\qquad |\mathcal T_m|=9m=18s.
\tag{13}
```

Write `K_D=H_2\otimes P`.  Associativity turns the full matrix in (13) into
`H_s\otimes K_D`.  Moreover

```math
\operatorname{tr}K_D=\operatorname{tr}H_2\operatorname{tr}P=0,
```

so deleting the diagonal changes neither `z^{\mathsf T}K_Dz` nor the product
state below.  Choose `z` from (11).  Equations (11)--(12) give the exact
Boolean witness

```math
\boxed{
Q(\mathcal T_m)\ge78s\sqrt s.
}
\tag{14}
```

Since

```math
\frac{78s\sqrt s}{(18s)^{3/2}}
=\frac{78}{18\sqrt{18}}
=1.0213764617\ldots>1,
\tag{15}
```

while principal submatrices of Paley conference matrices give
`q_N\le(1+o(1))N^{3/2}`, (13) is not globally minimizing for all sufficiently
large `r`, for every one of the 512 choices of `D`.  This is a genuine
Boolean-witness theorem, not an operator-norm or heuristic claim.

The theorem closes the most natural Sylvester/Hadamard mixed amplification.
It does not cover every symmetric Hadamard family, every conference factor,
or nonseparable/nonuniform gadgets.

## 4. Conference data (numerical only)

For orientation, coordinate ascent was run on symmetric Paley conference
factors of orders

```text
6, 14, 18, 30, 38, 42, 54, 62, 74, 102, 150, 174, 198
```

with `E=I`.  At outer order 198 the heuristic lower ratios
`Q(\mathcal T)/(9m\sqrt{9m-1})` were

```text
D=+I:       1.2656
D=-I:       1.3158
trace(D)=1: 1.1832
trace(D)=-1:1.1512
alternating:1.2116.
```

All tested values were above one, including balanced completions, but these
are only local-search lower bounds and are not used in any proof.  The precise
remaining conference target is to prove, for the 256-state kernel (3),

```math
\liminf_m\frac{Q(\mathcal T_m)}{27m^{3/2}}>1,
\tag{16}
```

or to find a counterexample.  The existing finite-channel lower constant
`c_2` is insufficient for (16).

## Reproducibility

- `tmp/check_mixed_order2_r17.py`: exact ground-pair lemma, complete
  order-eighteen audit, and the Sylvester constant.
- `tmp/mixed_product_r17.py`: Paley construction and heuristic coordinate
  ascent.  Its large-order values are explicitly numerical only.
