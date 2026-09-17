# Wave 45 memo: a hard common-core branching criterion

Status: all identities and the branching theorem below are **Verified**.
The required minimizer-specific collision inequality remains **Open**. The
finite counts are exact. No soft pressure or high-replica functional is used.

## 1. Hard incidence and its exact two-selector ratio

Fix a project row class

```math
\mathscr P=\{[z]:R_2(z)\le C_Rn^{9/4-c}\}
```

and let `nu=U_(mathscr P)` be its uniform projective law. The trace bound
makes this class nonempty at the project scale. The proof also works for any
explicitly supplied selector-independent prior, but choosing a prior after
knowing the best center would be circular. Put

```math
f_z(S)=\mathbf1\{Q(A[S])-|z_S^{\mathsf T}A[S]z_S|\le H\},
\qquad a_z=\mathbb E_{S\sim U_m}f_z(S),
\qquad J_1=\mathbb E_{z\sim\nu}a_z.
\tag{R45.1}
```

Assume `J_1>0`; centers with `a_z=0` then receive zero weight in every
size-biased formula below. Direct double counting gives

```math
J_2=\mathbb E_za_z^2
=\Pr\{f_Z(S)=f_Z(T)=1\},
\qquad
\boxed{\max_{z\in\mathscr P}a_z\ge\frac{J_2}{J_1}.}
\tag{R45.2}
```

Thus the exact averaged-collision sufficient lemma is

```math
\sum_{S,T}\nu\{z:f_z(S)=f_z(T)=1\}
\ge e^{-O(L_0)}\binom nm
\sum_S\nu\{z:f_z(S)=1\}.
\tag{R45.3}
```

It permits arbitrary empty selector pairs. It is stronger than the desired
existential conclusion and is not established by first-center averaging;
an exceptionally rare center can still be diluted in both `J_1` and `J_2`.

## 2. PSD common-core transport

For `0<=ell<m`, define the slice down-up kernel `K_ell`: from `S`, choose a
uniform `ell`-subset `R` of `S`, then choose a uniform `m`-set `T` containing
`R`. If `D_(ell,m)f(R)=E[f(S)|S contains R]`, then

```math
K_\ell=D_{\ell,m}^{\!*}D_{\ell,m}.
\tag{R45.4}
```

Hence `K_ell` is self-adjoint, Markov, and positive semidefinite on the
uniform slice. Its Johnson harmonic eigenvalues are

```math
\lambda_j
=\frac{(\ell)_j(n-m)_j}{(m)_j(n-\ell)_j}
\quad(0\le j\le\ell),
\qquad \lambda_j=0\quad(j>\ell).
\tag{R45.5}
```

The ratios of consecutive nonconstant eigenvalues are below one, so the
nonconstant operator norm is

```math
\boxed{\lambda_\ell=\lambda_1
=\frac{\ell(n-m)}{m(n-\ell)}.}
\tag{R45.6}
```

For an `ell`-core define its hard extension density

```math
b_z(R)=\Pr_{S\sim U_m}\{f_z(S)=1\mid S\supset R\}.
```

The factorization (R45.4) gives the exact common-core square identity

```math
\langle f_z,K_\ell f_z\rangle_{U_m}
=\mathbb E_{R\sim U_\ell}b_z(R)^2.
\tag{R45.7}
```

Define the incidence-conditioned retention

```math
\boxed{
\rho_\ell
=\frac{\mathbb E_{z\sim\nu}\langle f_z,K_\ell f_z\rangle}
       {\mathbb E_{z\sim\nu}a_z}
=\frac{\mathbb E_{z,R}b_z(R)^2}{J_1}.}
\tag{R45.8}
```

Equivalently: sample `(Z,S)` from `nu times U_m` conditioned on hard
incidence, choose a uniform `ell`-core of `S`, branch to a uniform
`m`-superset `T`, and ask whether the same center remains hard-favorable.
This is an averaged regeneration probability, not a uniform Helly premise.

## 3. Exact hard branching theorem

For a center with `a_z>0`, put

```math
r_t(z)=\frac{\langle f_z,K_\ell^tf_z\rangle}{a_z}.
```

Because the spectral measure of `K_ell` lies in `[0,1]`, Jensen gives
`r_t(z)>=r_1(z)^t`. A second Jensen inequality under the center law
size-biased by `a_z` gives

```math
\mathbb E_{a\nu/J_1}r_t(z)\ge\rho_\ell^t.
\tag{R45.9}
```

On the other hand, writing `f_z=a_z+g_z` and using (R45.6),

```math
r_t(z)
\le a_z+\lambda_\ell^t(1-a_z)
\le a_z+\lambda_\ell^t.
\tag{R45.10}
```

Combining (R45.9)--(R45.10) proves, for every integer `t>=1`,

```math
\boxed{
\max_{z\in\mathscr P}a_z
\ge\rho_\ell^t-\lambda_\ell^t.}
\tag{R45.11}
```

In particular the new hard common-core sufficient lemma is simply

```math
\boxed{
\rho_\ell-\lambda_\ell\ge e^{-O(L_0)}.}
\tag{R45.12}
```

This proves the hard exceptional-center degree directly at `t=1`. For a
fixed macroscopic core fraction, a fixed positive gap in (R45.12) would
already force constant degree and is therefore much stronger than needed.
Powers do not create positivity when `rho_ell<=lambda_ell`.

They can only amplify an existing small positive excess by a polynomial
mixing factor. To see the relevant near-identity scale, write
`ell=m-s`. Exactly,

```math
1-\lambda_{m-s}=\frac{sn}{m(n-m+s)}.
\tag{R45.13}
```

If `s=Theta(n/L_0)`, then `1-lambda=Theta(1/L_0)`. If

```math
\rho_{m-s}=1-\frac\beta{L_0}+o(L_0^{-1}),
\qquad
\lambda_{m-s}=1-\frac\alpha{L_0}+o(L_0^{-1}),
\qquad0\le\beta<\alpha,
\tag{R45.14}
```

then `t=Theta(L_0)` makes (R45.11) a positive constant, while
`t=Theta(L_0^2)` gives a difference of order `e^{-O(L_0)}`. More generally,
within the comparable-rate regime (R45.14), optimizing powers enlarges an
infinitesimal excess by only a polynomial `O(L_0)` factor. Thus (R45.11) is a precise branching mechanism,
not an exponential amplification miracle; the substantive missing theorem
is still the hard spectral excess (R45.12), or the rate gap (R45.14).

## 4. Finite exact-minimizer audit

The checker `hard_center_branching_r45_check.py` independently verifies the
kernel spectrum, (R45.7), (R45.11) for powers one through six, and all hard
incidence counts. It uses the old class `R_2<=2n(n-1)`, which is a stronger
row qualification and happens to contain every projective center in these
three examples.

| signing, `m` | `max a_z` | `J_1` | `J_2/J_1` | empty pairs | empty triples |
|:---|---:|---:|---:|---:|---:|
| `A_6,5` | `5/6` | `5/8` | `2/3` | `0/15` | `0/20` |
| `A_8,6` | `11/28` | `25/224` | `11/50` | `22/378` | `1700/3276` |
| `A_9,7` | `13/36` | `7/192` | `157/864` | `362/630` | `6070/7140` |

The exact common-core spectral gaps `rho_ell-lambda_ell` are:

| signing | successive `ell=0,...,m-1` gaps |
|:---|:---|
| `A_6,m=5` | `2/3, 16/25, 3/5, 8/15, 2/5` |
| `A_8,m=6` | `11/50, 33/175, 836/5625, 121/1250, 11/375, -11/225` |
| `A_9,m=7` | `157/864, 5275/32928, 4945/37044, 293/2940, 659/11760, -1/21168, -55/882` |

Thus the nearest-core criterion fails on exact `A_8`, and both one- and
two-replacement versions fail on exact `A_9`, even though a large hard-degree
center exists. This is a finite scoped falsifier to overly local regeneration,
not to a macroscopic common-core theorem and not to the asymptotic hard tail.
The many empty arbitrary pairs and triples coexist with positive averaged
gaps at less-local core sizes, exactly as (R45.8) permits.

## 5. Frontier recommendation

Retain (R45.12) as an exact hard, clipped, project-row conditional theorem.
It replaces a false uniform-overlap request by one averaged common-core
collision ratio and has a precise spectral threshold. It is not yet a
minimizer theorem and is strictly stronger than the existential degree;
finite exact minimizers show that core size is essential and rule out nearest-
neighbor regeneration. Any continuation should seek a global-minimality
lower bound on `rho_ell-lambda_ell` for a nonlocal or near-identity core scale,
while keeping the hard indicator and row-supported base law explicit.
