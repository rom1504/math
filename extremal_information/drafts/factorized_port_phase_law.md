# A factorized port phase law and its exact-sign obstruction

**Status.** Rigorous theorem, scalable counterexample, and independent
audit. On a Cartesian product of finite pole algebras, the large-arity
labelled trust response is controlled by a signed first-moment phase,
uniformly over all endpoint labels. Combined with relative product
synchronization, this gives a thermodynamic limit for a restricted dense
exact-sign parent class. Exact product closure alone does **not** give a
limit: an explicit order-16 regular-Hadamard construction with zero relative
defect can retain an oscillating empirical factor phase.

This concerns explicit signing sequences, not the minimizing values `M_n`.

## 1. Cartesian pole interfaces

For factor `t`, let `H_t` be a real symmetric `n_t by n_t` matrix with

```math
\|H_t\|_{op}\le r_t.                                  \tag{FP.1}
```

Choose a Boolean base pole `a_t` and Boolean relative generators
`g_(t,1),...,g_(t,q_t)`. Let

```math
U_t=\operatorname {span}\left\{
a_t\bigodot_{j\in S}g_{t,j}:S\subseteq[q_t]\right\}. \tag{FP.2}
```

At depth `L`, put

```math
H^{(L)}=\bigotimes_{t=1}^LH_t,\quad
n_L=\prod_{t=1}^Ln_t,\quad r_L=\prod_{t=1}^Lr_t,
\quad p_L=1+\sum_{t=1}^Lq_t.                         \tag{FP.3}
```

The base port is `w_0=otimes_t a_t`. For every `(t,j)`, the corresponding
port is `w_(t,j)=w_0` multiplied coordinatewise by `g_(t,j)` in factor `t`.
All odd products of these `p_L` ports lie in

```math
U^{(L)}=\bigotimes_{t=1}^LU_t.                       \tag{FP.4}
```

For an endpoint label `epsilon` and integer shore width `m_L`, write

```math
z_\epsilon=\epsilon_0w_0+
 \sum_{t,j}\epsilon_{t,j}w_{t,j}                    \tag{FP.5}
```

and define the labelled Boolean trust response

```math
B_L(\epsilon)=
\max_{x\in\{+-1\}^{n_L},\ \sigma\in\{+-1\}}
\left\{{\sigma\over2}x^TH^{(L)}x+m_Lz_\epsilon^Tx\right\}.
                                                               \tag{FP.6}
```

For every `L`, fix once and for all one antipodally odd tie-broken majority
selector `tau_L:{+-1}^(p_L)->{+-1}`. It agrees with the sign of the
coordinate sum off the zero layer; on ties use any fixed antipodally
consistent rule. The same selector is used for every endpoint label at that
level.

Let `alpha_t` be a uniform row coordinate in factor `t` and set

```math
X_{t,j}=g_{t,j}(\alpha_t),\qquad
\mu_{t,j}=\mathbb E X_{t,j}.                       \tag{FP.7}
```

The vectors `X_t=(X_(t,j))_(j<=q_t)` are independent across factors; no
independence within one block is assumed. Gauging every row by its base
sign gives the exact support identity

```math
{\|z_\epsilon\|_1\over n_L}
=\mathbb E\left|\epsilon_0+
       \sum_{t,j}\epsilon_{t,j}X_{t,j}\right|.      \tag{FP.8}
```

## 2. Uniform first-moment collapse

Put

```math
V_L=\sum_{t=1}^Lq_t^2,
\qquad
d_L(\epsilon)={\epsilon_0+
 \sum_{t,j}\epsilon_{t,j}\mu_{t,j}\over p_L},
\qquad
\theta_L={1+\sum_{t,j}|\mu_{t,j}|\over p_L}.       \tag{FP.9}
```

### Lemma FP.1 (factorized support collapse)

Uniformly over all `2^(p_L)` endpoint labels,

```math
0\le {\|z_\epsilon\|_1\over n_Lp_L}
       -|d_L(\epsilon)|
\le {\sqrt {V_L}\over p_L}.                        \tag{FP.10}
```

Consequently,

```math
\left|{1\over n_Lp_L}\max_\epsilon\|z_\epsilon\|_1
      -\theta_L\right|
\le {\sqrt {V_L}\over p_L}.                        \tag{FP.11}
```

#### Proof

For fixed `epsilon`, let

```math
Y_t=\sum_{j=1}^{q_t}\epsilon_{t,j}X_{t,j}.
```

The `Y_t` are independent and `Var(Y_t)<=q_t^2`. Jensen and Cauchy--Schwarz
therefore give

```math
0\le \mathbb E\left|\epsilon_0+\sum_tY_t\right|
 -\left|\epsilon_0+\sum_t\mathbb EY_t\right|
\le \sqrt{\sum_t\operatorname {Var}(Y_t)}
\le\sqrt {V_L}.                                    \tag{FP.12}
```

This proves (FP.10). The largest absolute mean is exactly
`1+sum_(t,j)|mu_(t,j)|`: choose every endpoint sign to align the means and
choose `epsilon_0` with their common sign. Maximizing (FP.10) proves
(FP.11). `square`

Thus the exact projective histogram, which can have `2^(p_L-1)` bins, has a
strict asymptotic quotient for this one max-over-label query. If

```math
{\sqrt {V_L}\over p_L}\longrightarrow0,             \tag{FP.13}
```

then only the signed first moments in (FP.9) survive at leading order.
Uniformly bounded nonzero `q_t` imply (FP.13).

## 3. Relative synchronization and the limit theorem

Let `T_L=H^(L)/r_L`. Define the relative defect on the represented pole
span by

```math
e_L=\sup_{0\ne u\in U^{(L)}}
 {\langle u,(I-T_L)u\rangle\over\|u\|_2^2}.         \tag{FP.14}
```

It is nonnegative under (FP.1). If the factor defects are `delta_t`, where

```math
\delta_t=\sup_{0\ne u\in U_t}
 {\langle u,(I-H_t/r_t)u\rangle\over\|u\|_2^2},    \tag{FP.14a}
```

the Cartesian relative-synchronization theorem gives

```math
e_L\le\sum_{t=1}^L\delta_t;                        \tag{FP.15}
```

For positive compressed contractions on the full tensor span (FP.4), it
gives the sharper exact identity
`e_L=1-prod_t(1-delta_t)`. The same statements apply row by row to a
triangular family whose factor presentations improve with target depth.
Exact pole closure is the important special case `e_L=0`.

Put

```math
c_L={m_Lp_L\over r_L}.                              \tag{FP.16}
```

### Theorem FP.2 (restricted Cartesian trust limit)

For every endpoint label,

```math
\left|{B_L(\epsilon)\over r_Ln_L}
 -\left({1\over2}+c_L|d_L(\epsilon)|\right)\right|
\le {e_L\over2}+c_L{\sqrt {V_L}\over p_L}.         \tag{FP.17}
```

Moreover,

```math
\left|{\max_\epsilon B_L(\epsilon)\over r_Ln_L}
 -\left({1\over2}+c_L\theta_L\right)\right|
\le {e_L\over2}+c_L{\sqrt {V_L}\over p_L}.         \tag{FP.18}
```

Hence, if

```math
c_L\to c,\qquad e_L\to0,
\qquad {\sqrt {V_L}\over p_L}\to0,
\qquad \theta_L\to\theta,                         \tag{FP.19}
```

then

```math
{\max_\epsilon B_L(\epsilon)\over r_Ln_L}
\longrightarrow {1\over2}+c\theta.                \tag{FP.20}
```

For a prescribed sequence of endpoint labels it is enough that
`|d_L(epsilon^(L))|` converge. Conversely, in the exact case `e_L=0`, if
`c_L->c>0` and (FP.13) holds, convergence of the maximum is equivalent to
convergence of `theta_L`. Thus (FP.19) isolates the minimal scalar phase
law in this product class.

#### Proof

The majority-selector witness has all its active Fourier products in
`U^(L)`: antipodal oddness puts its Walsh support on odd products, and the
tie rule fixed after (FP.6) makes one common witness definition valid for
every endpoint label. Relative robust synchronization and (FP.14) give

```math
0\le {1\over2}+{m_L\|z_\epsilon\|_1\over r_Ln_L}
 -{B_L(\epsilon)\over r_Ln_L}\le {e_L\over2}.       \tag{FP.21}
```

Since `m_L/r_L=c_L/p_L`, combine (FP.21) first with (FP.10), then with
(FP.11). This proves (FP.17)--(FP.18), and the limit statements follow.
`square`

Ordinary summability `sum_t delta_t<infinity` for one fixed infinite factor
sequence is **not** the same as `e_L->0`: the robust theorem accumulates
defect from the root and can leave a fixed leading-width interval. Thus a
claim with merely summable fixed defects requires an additional theorem
showing convergence of the actual selector deficit. Factorwise relative
synchronization alone does not supply that second phase law.

## 4. Hollow exact-sign completion

Suppose now that every `H_t` is an entrywise sign matrix, the tensor product
has trace zero, and

```math
{r_L\over\sqrt {n_L}}\to\rho,qquad
d_L^{aux}:=m_Lp_L=o(n_L).                           \tag{FP.22}
```

Delete the diagonal of `H^(L)`, append one shore of width `m_L` for each
port, and fill all auxiliary--auxiliary edges by an arbitrary public hollow
signing `C_L`. Call the completed parent `P_L` and its order
`N_L=n_L+d_L^(aux)`. Trace zero makes diagonal deletion invisible on the
Boolean cube. The completion estimate is

```math
\left|Q(P_L)-\max_\epsilon B_L(\epsilon)\right|
\le Q(C_L)\le {d_L^{aux}\choose2}.                 \tag{FP.23}
```

If also `(d_L^(aux))^2=o(r_Ln_L)` and (FP.19) holds, Theorem FP.2 yields

```math
{Q(P_L)\over N_L^{3/2}}
\longrightarrow
\rho\left({1\over2}+c\theta\right).               \tag{FP.24}
```

This is a genuine dense exact-sign restricted thermodynamic limit. It is
not a statement about minimizing over all signings.

## 5. The repeated order-16 pole algebra

The regular-Hadamard seed of PC.3 has order `16`, top eigenvalue `4`, trace
zero, and Boolean positive poles `a,b,c,a*b*c` with

```math
{a^Tb\over16}={1\over2},\qquad {a^Tc\over16}=0.    \tag{FP.25}
```

Use in each tensor factor the two relative generators `a*b` and `a*c`.
Then `q_t=2`,

```math
\mathbb E(a*b)={1\over2},\qquad \mathbb E(a*c)=0, \tag{FP.26}
```

and all affine products are exact positive poles. Thus

```math
p_L=2L+1,\quad e_L=0,\quad
\theta_L={1+L/2\over2L+1}\to{1\over4}.             \tag{FP.27}
```

With `m_L=floor(4^L/(2L+1))`, one has `c_L->1`, while
`n_L=16^L`, `r_L=4^L=sqrt(n_L)`. Every public exact-sign completion in
Section 4 therefore satisfies

```math
\boxed{{Q(P_L)\over |P_L|^{3/2}}\longrightarrow {3\over4}.} \tag{FP.28}
```

This is near the original scale but remains a deliberately structured
logarithmic-interface class.

## 6. Zero-defect phase obstruction

The phase hypothesis in (FP.19) cannot be deleted. Use the same `H_16` in
every factor, but take only one relative generator per factor:

* an `X` factor uses `a*b`, whose mean is `1/2`;
* a `Y` factor uses `a*c`, whose mean is `0`.

Both choices have exact affine pole closure, so every factor and every
Cartesian product has relative defect zero. Let the factor types occur in
alternating blocks. If `N_(k-1)` factors have already been chosen, choose
the next block length at least `k^2N_(k-1)`, using `X` for odd `k` and `Y`
for even `k`. At odd block endpoints the density of `X` factors tends to
one; at even endpoints it tends to zero.

Here `p_L=L+1`, and if `K_L` denotes the number of `X` factors, then

```math
\theta_L={1+K_L/2\over L+1}.                       \tag{FP.29}
```

Take `m_L=floor(4^L/(L+1))` and complete exactly as above. Equations
(FP.18), (FP.23), and `sqrt(V_L)/p_L=sqrt(L)/(L+1)->0` give

```math
\lim_{\substack{L\to\infty\\L\text{ odd-block end}}}
 {Q(P_L)\over|P_L|^{3/2}}=1,
\qquad
\lim_{\substack{L\to\infty\\L\text{ even-block end}}}
 {Q(P_L)\over|P_L|^{3/2}}={1\over2}.              \tag{FP.30}
```

Thus exact relative synchronization, exact Cartesian histogram evolution,
sublinear completion cost, and the correct `n^(3/2)` scale do not force a
limit. The missing datum is the empirical factor phase (FP.29), not an
uncontrolled Boolean landscape.

This differs from the Walsh-prefix obstruction of Theorem 30.1. Here the
failure already occurs **at tensor endpoints**, because the factor-type
empirical phase is deliberately nonconvergent. In Theorem 30.1 one fixed
factor has convergent geometric tensor endpoints, while interpolation to
all integer orders creates a continuous mantissa phase. Consequently a
restricted all-order theorem needs two logically separate controls:

1. convergence of the factor/histogram phase at tensor depths; and
2. constancy (or synchronization) of the scale phase used to fill the gaps
   between those depths.

## 7. Verification

The companion script reconstructs the order-16 seed, checks its pole
relations and feature law, exhausts the first-moment inequalities on small
Cartesian products, verifies the `3/4` normalization, and checks the
alternating-block phase bounds:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_factorized_port_phase_law.py
```

An independent proof and normalization check is recorded in
[`factorized_port_phase_law_audit.md`](factorized_port_phase_law_audit.md).
