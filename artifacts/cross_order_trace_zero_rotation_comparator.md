# Trace-zero Hadamard cores: the exact independent-row comparator

Status: **exact reduction, conditional rather than an exponent improvement**.
For an explicit symmetric Hadamard order `k` with `r<=k<4r`, the rotated
core defect is exactly an orthogonal-cube replacement defect minus a
nonnegative one-sided-asymmetry credit.  The pressure and cap formulas below
retain that credit and give direct arrows to `P_(kr)-kP_r` and
`b_(kr)-kb_r`.  What is not proved is the required subleading replacement
bound for growing-order actual minimizers.

Throughout,

```math
H_A(x)={1\over2}x^{\mathsf T}Ax,
\qquad
\phi_A(s)=\log\mathbb E_x\cosh(sH_A(x)),
\qquad s={\beta\over\sqrt r}.                       \tag{1}
```

## 1. A balanced symmetric Hadamard at a comparable order

Put

```math
H_2=\begin{pmatrix}1&1\\1&-1\end{pmatrix},
\qquad R_4=J_4-2I_4,
\qquad R_q=R_4^{\otimes d},\quad q=4^d,             \tag{2}
```

and let

```math
T_k=H_2\otimes R_q,
\qquad k=2q.                                        \tag{3}
```

Then `T_k` is symmetric, has only sign entries, and

```math
T_k^2=kI_k,
\qquad \operatorname{Tr}T_k=0.                     \tag{4}
```

For every `r>=2`, choosing the smallest member of the sequence
`2,8,32,...` which is at least `r` gives

```math
r\le k<4r.                                          \tag{5}
```

In an orthogonal diagonalization

```math
T_k=\sqrt k\,U\Sigma U^{\mathsf T},
\qquad \Sigma=\operatorname{diag}(\sigma_1,\ldots,\sigma_k),          \tag{6}
```

equation (4) says that exactly `k/2` of the `sigma_p` are `+1` and
exactly `k/2` are `-1`.

The direct implication of (5) is that the completion errors below are only
`O_beta(k)=O_beta((kr)^(1/2))` for pressure and `O(k sqrt(r))=o(kr)` in the
`b` normalization.  Thus selecting a trace-zero outer does not create a new
leading completion cost.

## 2. Exact pressure decomposition

For a Boolean `k` by `r` matrix `X`, put `Y=U^T X`.  Its columns are
independent with law

```math
\nu_U=(U^{\mathsf T})_\#\operatorname{Unif}\{\mathord\pm1\}^k.
```

The diagonalized-core identity is

```math
{\beta\over\sqrt{kr}}H_{T_k\otimes A}(X)
=s\sum_{p=1}^k\sigma_pH_A(y_p).                    \tag{7}
```

Define the rotated and independent-row partition functions

```math
\begin{aligned}
Z_{\rm rot}(T_k,A)
&=\mathbb E_{Y\sim\nu_U^{\otimes r}}
  \cosh\!\left(s\sum_p\sigma_pH_A(y_p)\right),\\
Z_{\rm ind}(\Sigma,A)
&=\mathbb E_{\epsilon\in\{\mathord\pm1\}^{k\times r}}
  \cosh\!\left(s\sum_p\sigma_pH_A(\epsilon_p)\right),              \tag{8}
\end{aligned}
```

and

```math
Z_+(A)=\mathbb E_xe^{sH_A(x)},
\qquad Z_-(A)=\mathbb E_xe^{-sH_A(x)},
\qquad
\delta_A={1\over2}\log{Z_+(A)\over Z_-(A)}.         \tag{9}
```

### Proposition 2.1 (balanced comparator identity)

For every hollow symmetric signing `A`,

```math
\boxed{
Z_{\rm ind}=(Z_+Z_-)^{k/2}}
                                                               \tag{10}
```

and the exact core defect is

```math
\boxed{
\phi_{T_k\otimes A}\!\left({\beta\over\sqrt{kr}}\right)
-k\phi_A\!\left({\beta\over\sqrt r}\right)
=\bigl(\log Z_{\rm rot}-\log Z_{\rm ind}\bigr)
-k\log\cosh\delta_A.}                              \tag{11}
```

Proof.  Independence of the Boolean rows in the second line of (8) gives

```math
\begin{aligned}
Z_{\rm ind}
={1\over2}\left[
 \prod_p Z_{\sigma_p}+\prod_p Z_{-\sigma_p}\right]
=(Z_+Z_-)^{k/2},
\end{aligned}
```

because the signature is balanced.  On the other hand,

```math
e^{\phi_A(s)}={Z_++Z_-\over2}
=\sqrt{Z_+Z_-}\cosh\delta_A.                       \tag{12}
```

Equation (7) says that the first pressure in (11) is `log Z_rot`.
Taking logarithms in (10)--(12) proves (11).  \(\square\)

The term `k log cosh(delta_A)` is nonnegative.  It is the exact credit for
an asymmetric positive/negative child partition function; discarding it
strengthens the needed rotation theorem and can lose a leading term.

### Direct arrow to the permanent pressure

Let `A` be an actual order-`r` pressure minimizer, so that
`phi_A(s)=P_r(beta)`.  Annealed filling of the
`binom(k,2)r` zero matching-coordinate edges gives

```math
\boxed{
\begin{aligned}
P_{kr}(\beta)-kP_r(\beta)
\le{}&\log Z_{\rm rot}-\log Z_{\rm ind}
-k\log\cosh\delta_A\\
&+{\beta^2(k-1)\over4}.
\end{aligned}}                                      \tag{13}
```

Indeed, each filled edge costs `log cosh(beta/sqrt(kr))`, and

```math
{k(k-1)r\over2}\log\cosh{\beta\over\sqrt{kr}}
\le {\beta^2(k-1)\over4}.                           \tag{14}
```

Consequently the precise missing one-sided estimate is

```math
\boxed{
\log Z_{\rm rot}-\log Z_{\rm ind}
\le k\log\cosh\delta_A+o_\beta(kr)}.                \tag{15}
```

For the comparable choice (5), (15) implies
`P_(kr)(beta)-kP_r(beta)=o_beta(kr)`.  Quantitatively, an
`O_beta((kr)^(1-delta))` remainder in (15) gives

```math
P_{kr}(\beta)-kP_r(\beta)
=O_\beta\!\left((kr)^{1-\min\{\delta,1/2\}}\right). \tag{16}
```

No estimate of the form (15) is proved here.

## 3. Exact cap decomposition

Put

```math
M_+(A)=\max_xH_A(x),
\qquad M_-(A)=-\min_xH_A(x),
\qquad M(A)=\max\{M_+(A),M_-(A)\}.                  \tag{17}
```

Define

```math
\begin{aligned}
C_{\rm rot}(T_k,A)
&=\max_{Y\in\operatorname{supp}(\nu_U^{\otimes r})}
 \left|\sum_p\sigma_pH_A(y_p)\right|,\\
C_{\rm ind}(\Sigma,A)
&=\max_{\epsilon\in\{\mathord\pm1\}^{k\times r}}
 \left|\sum_p\sigma_pH_A(\epsilon_p)\right|,\\
\Gamma_{k,r}(T_k,A)&=C_{\rm rot}-C_{\rm ind}.       \tag{18}
\end{aligned}
```

### Proposition 3.1 (balanced cap comparator)

For every `A`,

```math
\boxed{
C_{\rm ind}={k\over2}(M_++M_-)}                     \tag{19}
```

and

```math
\boxed{
M(T_k\otimes A)-k^{3/2}M(A)
=\sqrt k\left[
 \Gamma_{k,r}(T_k,A)-{k\over2}|M_+(A)-M_-(A)|
 \right].}                                         \tag{20}
```

Proof.  To maximize `sum_p sigma_p H_A(epsilon_p)`, independently put a
positive-sign row at a positive maximizer and a negative-sign row at a
negative minimizer.  This gives `k(M_++M_-)/2`; minimizing gives its
negative.  This proves (19).  Equation (7) without the temperature gives

```math
M(T_k\otimes A)=\sqrt k\,C_{\rm rot}.               \tag{21}
```

Finally,

```math
kM-C_{\rm ind}
={k\over2}\bigl(2M-M_+-M_-\bigr)
={k\over2}|M_+-M_-|,
```

which together with (18)--(21) proves (20).  \(\square\)

### Direct arrow to `b_(kr)-kb_r`

Let `A` now be an actual cap minimizer, `M(A)=M_r`, and put

```math
\mathcal E_{k,r}(T_k,A)
=\Gamma_{k,r}(T_k,A)-{k\over2}|M_+(A)-M_-(A)|,
\qquad
C_{k,r}=\sqrt{k(k-1)r(kr+2)\log2}.                  \tag{22}
```

Random sign completion and concavity of `u^(2/3)` give the explicit
implication

```math
\boxed{
b_{kr}-kb_r
\le {2\over3}(k^{3/2}M_r)^{-1/3}
 \bigl(\sqrt k\,\mathcal E_{k,r}+C_{k,r}\bigr)_+.} \tag{23}
```

Thus the precise cap-side missing estimate is

```math
\boxed{
\Gamma_{k,r}(T_k,A)
\le {k\over2}|M_+(A)-M_-(A)|+o(k r^{3/2})}          \tag{24}
```

along some sequence of actual minimizing children.  Using the project's
uniform lower bound `M_r>=c r^(3/2)`, (23) shows more quantitatively that

```math
\mathcal E_{k,r}=O(k r^{3/2-\delta})
\quad\Longrightarrow\quad
b_{kr}-kb_r=O(k r^{1-\delta}+k\sqrt r).             \tag{25}
```

For (5), the right side is `o(kr)` for every `delta>0`.  Again, (24) is an
optimizer-specific rotation theorem, not a conclusion proved here.

## 4. The credit cannot be replaced by zero

The tempting universal inequalities

```math
C_{\rm rot}\le C_{\rm ind}
\quad\hbox{or}\quad
\log Z_{\rm rot}\le\log Z_{\rm ind}                \tag{26}
```

are false even for the first balanced outer and an actual cap-minimizing
child.  Take

```math
A=\begin{pmatrix}
0&-1&1&-1&1\\
-1&0&-1&1&1\\
1&-1&0&1&-1\\
-1&1&1&0&-1\\
1&1&-1&-1&0
\end{pmatrix},
\qquad A^2=5I-J.                                    \tag{27}
```

For Boolean `x`, the odd number `1^T x` and (27) give

```math
|H_A(x)|\le {1\over2}\sqrt5\sqrt{25-(1^Tx)^2}
\le\sqrt{30}.
```

All energies are even, hence `M(A)<=4`.  The spins

```text
(-1,-1, 1, 1,-1),  (-1,-1,-1, 1, 1)
```

have energies `4` and `-4`, respectively.
Conversely every order-five signing has
`E_x H_A(x)^2=binom(5,2)=10`; since its energies are even, its cap is at
least `4`.  Thus this `A` is an actual minimizer and

```math
M_+(A)=M_-(A)=M_5=4.                                \tag{28}
```

For `T_2=H_2`, put `K=H_2 tensor A` and write a Boolean state as two
blocks `z_1,z_2`.  Since each block sum `s_a=1^Tz_a` is odd,

```math
\|Kz\|_2^2
=2\sum_{a=1}^2(25-s_a^2)\le96.
```

Thus `|H_K(z)|<=sqrt(10)sqrt(96)/2=4sqrt(15)<16`.  The nonzero-support
graph of `K` has even degree eight, and `H_K(1)=0` because `A1=0`.
Successive spin flips therefore show that every energy is a multiple of
four.  Hence `M(K)<=12`.  The two blocks

```text
(-1,-1,-1,-1, 1),  (-1,-1, 1, 1,-1)
```

attain energy `-12`, so

```math
M(H_2\otimes A)=12.                                 \tag{29}
```

Equations (18), (19), and (29) now give

```math
C_{\rm rot}=6\sqrt2,
\qquad C_{\rm ind}=8,
\qquad
\Gamma_{2,5}=6\sqrt2-8>0,                           \tag{30}
```

while the asymmetry credit in (20) is zero.  This disproves the cap
inequality in (26).  The zero-temperature slopes of the two finite
partition functions are `C_rot` and `C_ind`, respectively, so (30) also
disproves the pressure inequality in (26) for all sufficiently large
`beta`.  This finite witness only delimits the replacement theorem; it is
not a growing-order obstruction.

## 5. Why the immediate norm relaxation does not close (24)

For completeness, define the bipartite vector relaxation

```math
G(A)={1\over2}\max_{\|u_i\|=\|v_j\|=1}
 \left|\sum_{i,j}a_{ij}\langle u_i,v_j\rangle\right|.               \tag{31}
```

Writing a parent state by its Boolean columns `c_i`, setting
`u_i=c_i/sqrt(k)` and `B=T_k/sqrt(k)`, and then taking `v_j=Bu_j` proves

```math
M(T_k\otimes A)\le k^{3/2}G(A).                     \tag{32}
```

Hence the norm route would imply directly

```math
b_{kr}-kb_r
\le {2\over3}(k^{3/2}M_r)^{-1/3}
 \bigl(k^{3/2}[G(A)-M_r]+C_{k,r}\bigr)_+.           \tag{33}
```

In particular `G(A)-M_r=O(r^(3/2-delta))` would give the same conclusion
as (25).  But this premise already contains essentially the sharp scalar
answer: choosing `u_i=e_i` and

```math
v_j={1\over\sqrt{r-1}}(a_{1j},\ldots,a_{rj})
```

gives, for every signing,

```math
G(A)\ge {r\sqrt{r-1}\over2}.                        \tag{34}
```

Thus `G(A)-M_r=o(r^(3/2))`, combined with the standard random-signing upper
bound, already forces `M_r/r^(3/2)->1/2`.  The unrestricted norm relaxation
therefore replaces the rotated-core theorem by a statement at least as
strong as determining the extremal constant; it supplies no independent
closure.

## 6. Boundary

Equations (11) and (20) are exact and work for every child.  Equations (13)
and (23) are the promised direct composition implications for actual
minimizers.  They show that the correct target is a **net** rotation defect:
the raw rotated-versus-independent discrepancy may be leading, provided it
is paid for by the child's positive/negative asymmetry.

No bound of the form (15) or (24) is currently proved for growing-order
actual minimizers.  A successful continuation must use optimizer-specific
control of the orthogonal-cube replacement defect while retaining the
asymmetry credit; covariance matching, `T^2=kI`, or the bipartite vector norm
alone does not provide it.
