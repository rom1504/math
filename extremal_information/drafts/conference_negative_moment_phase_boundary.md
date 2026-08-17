# The conference negative-moment phase boundary

**Status.** Task-local theorem draft.  This note does not modify the frozen
thin-tail-conditioned source.  It proves an explicit fixed-tilt wall beyond
the exact-small-tilt asymptotic, and it identifies an exact necessary-and-
sufficient lower-tail theorem for extending the wall to every fixed tilt.
It also audits the March 2026 spin-glass large-deviation result against that
missing theorem.

The result does not prove that a finite-tilt phase transition exists.  It
shows precisely what would create one: a speed-`r` lower deviation of bridge
pressure.  Conversely, superexponential rarity at speed `r` is equivalent to
the absence of any finite-tilt transition.

The evidentiary status of the three statements is different:

| statement | exact range | status / increment |
|---|---|---|
| `R^E_lambda/r -> h_beta` | `0<lambda<lambda_*(beta)` | frozen CT theorem; not reproved here |
| `liminf (R^E_lambda-T_r)/r >= gamma-lambda beta^2/4` | every fixed `lambda>0`; positive for `lambda<4gamma/beta^2` | new transfer of the archived global MGF bound through uniform thin-tail conditioning |
| `R^E_lambda/r -> h_beta` for every fixed `lambda` | all fixed `lambda>0` | open; NP.2 reduces it exactly to the superexponential lower-tail lemma |

Theorem NP.2 is an elementary exponential-moment/tail equivalence, not a
new large-deviation estimate.  Its value is quantitative research steering:
it proves that improving constants in speed-`r` concentration can never
establish the all-fixed-tilt statement.  A genuinely faster lower-tail speed
is necessary and sufficient.

## 1. Conference output notation

Retain the normalization of
`thin_tail_entropy_bridge_no_go.md`.  For a symmetric conference signing
`A_r`, a uniform sign bridge `B`, and `epsilon=+-1`, put

```math
S_{\epsilon,B}
=\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix},
\qquad
t={\beta\over\sqrt{2r}},
\tag{NP.1}
```

```math
L_{\epsilon,B}
=\log\left[
 2^{-2r}\sum_z\cosh\{tH_{S_{\epsilon,B}}(z)\}
 \right].
\tag{NP.2}
```

Here `H_S(z)=z^TSz/2`.  Write `U_r` for the uniform joint output law and

```math
\mathcal R_{\lambda,r}
=-{1\over\lambda}\log\mathbb E_{U_r}e^{-\lambda L}.
\tag{NP.3}
```

For

```math
0<\beta<{\sqrt2\over6},
\tag{NP.4}
```

the audited conference theorem gives, separately for both orientations,

```math
{1\over r}\mathbb E_B L_{\epsilon,B}\longrightarrow h_\beta,
\qquad
{T_r\over r}\longrightarrow\tau_\beta,
\tag{NP.5}
```

where

```math
T_r=2\log\overline Z_r(A_r,\beta/\sqrt r),
\tag{NP.6}
```

```math
\tau_\beta=2\psi(\beta),
\qquad
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4},
\qquad
\gamma(\beta)=h_\beta-\tau_\beta>0.
\tag{NP.7}
```

It also gives `L/r -> h_beta` in probability.  The event `E_r` from the
frozen thin-tail note has

```math
U_r(\mathcal E_r^c)\le2e^{-(2-2\log2)r},
\tag{NP.8}
```

and every output retained by `E_r` satisfies the common Theorem 36.26 tail.
Denote the conditioned law and soft minimum by `U_r^E` and
`R^E_(lambda,r)`.

## 2. An explicit wall for a larger tilt interval

### Theorem NP.1 (global transport wall, preserved by thin-tail conditioning)

For every fixed `lambda>0`,

```math
\boxed{
\liminf_{r\to\infty}
 {\mathcal R^E_{\lambda,r}-T_r\over r}
\ge \gamma(\beta)-{\lambda\beta^2\over4}.}
\tag{NP.9}
```

Consequently every fixed

```math
\boxed{0<\lambda<{4\gamma(\beta)\over\beta^2}}
\tag{NP.10}
```

has a strictly positive linear same-temperature conference defect even
after every retained output is required to have the common bounded-cap thin
tail.

**Proof.**  Fix one orientation and write `f(B)=L_(epsilon,B)`.  Flipping
one bridge bit changes every Hamiltonian value by `2` and therefore changes
`f` by at most

```math
c_e=2t={2\beta\over\sqrt{2r}}.
\tag{NP.11}
```

There are `r^2` bridge bits, so the bounded-difference exponential lemma
gives

```math
\log\mathbb E_B
 \exp\{-\lambda(f-\mathbb E_Bf)\}
\le{\lambda^2\over8}\sum_ec_e^2
={\lambda^2\beta^2r\over4}.
\tag{NP.12}
```

Both orientation means are `h_beta r+o(r)` by (NP.5).  Averaging their two
negative moments therefore yields

```math
\mathcal R_{\lambda,r}
\ge h_\beta r-{\lambda\beta^2r\over4}+o(r).
\tag{NP.13}
```

Conditioning can only increase this lower bound up to its negligible
normalization correction:

```math
\begin{aligned}
\mathbb E_{U_r^E}e^{-\lambda L}
&\le {\mathbb E_{U_r}e^{-\lambda L}\over U_r(\mathcal E_r)},\\
\mathcal R^E_{\lambda,r}
&\ge\mathcal R_{\lambda,r}
  +{1\over\lambda}\log U_r(\mathcal E_r)
 =\mathcal R_{\lambda,r}+o(1).
\end{aligned}
\tag{NP.14}
```

Subtract (NP.5)--(NP.7) to obtain (NP.9). `square`

This is an explicit **positive-wall** interval.  The frozen theorem is
stronger, namely `R^E/r -> h_beta`, on a possibly different nonzero interval
`lambda<lambda_*(beta)`.  No ordering between the unpublished constant
`lambda_*(beta)` and the explicit endpoint in (NP.10) is claimed.

### Entropy interpretation and sharpness of the method

For a law `q` on fixed-orientation bridges, the same bounded-difference MGF
and the entropy variational formula give

```math
\mathbb E_U L-\mathbb E_qL
\le\beta\sqrt{rD(q\Vert U)}.
\tag{NP.15}
```

The orientation contributes only `o(r)` in the joint law.  If
`d=D(q||U)/r`, then every entropy-regularized competitor consequently has

```math
{1\over r}\left{
\mathbb E_qL+{1\over\lambda}D(q\Vert U)
\right}
\ge h_\beta-\beta\sqrt d+{d\over\lambda}+o(1).
\tag{NP.16}
```

Minimizing the displayed quadratic in `sqrt(d)` gives

```math
\inf_{d\ge0}\left\{-\beta\sqrt d+{d\over\lambda}\right\}
=-{\lambda\beta^2\over4}.
\tag{NP.17}
```

Thus (NP.10) is the optimal conclusion obtainable from the global
bounded-difference/transport inequality alone.  Extending it requires a
strictly stronger lower-tail theorem, not a different optimization of the
same MGF estimate.

## 3. Exact equivalence for all fixed tilts

The next theorem is elementary but decisive.  It identifies the missing
probabilistic statement without assuming an LDP.

### Theorem NP.2 (all-fixed-tilt annealing iff lower deviations are superexponential)

Let `mu_r` be probability laws and let `X_r>=0` satisfy

```math
{X_r\over r}\longrightarrow h>0
\quad\hbox{in }\mu_r\hbox{-probability}.
\tag{NP.18}
```

The following are equivalent.

1. For every fixed `lambda>0`,

   ```math
   -{1\over\lambda r}\log
     \mathbb E_{\mu_r}e^{-\lambda X_r}\longrightarrow h.
   \tag{NP.19}
   ```

2. For every fixed `a<h`,

   ```math
   \boxed{
   {1\over r}\log\mu_r\{X_r\le ar\}\longrightarrow-\infty.}
   \tag{NP.20}
   ```

In words: every fixed disorder tilt remains on the typical branch exactly
when every fixed lower deviation is superexponentially rare at speed `r`.

**Proof.**  Assume (NP.20).  For fixed `epsilon>0`, convergence in
probability gives

```math
\mathbb E e^{-\lambda X_r}
\ge e^{-\lambda(h+\epsilon)r}
   \mu_r\{X_r\le(h+\epsilon)r\}
=e^{-\lambda(h+\epsilon)r+o(r)}.
\tag{NP.21}
```

On the other hand,

```math
\mathbb E e^{-\lambda X_r}
\le e^{-\lambda(h-\epsilon)r}
 +\mu_r\{X_r<(h-\epsilon)r\}.
\tag{NP.22}
```

The second term is smaller than `e^(-Kr)` for every fixed `K` eventually,
so (NP.19) follows by sending `epsilon` to zero.

Conversely, for `a<h`,

```math
\mu_r\{X_r\le ar\}e^{-\lambda ar}
\le\mathbb E e^{-\lambda X_r}
=e^{-\lambda hr+o(r)}.
\tag{NP.23}
```

Thus

```math
\limsup_r{1\over r}\log\mu_r\{X_r\le ar\}
\le-\lambda(h-a)
\tag{NP.24}
```

for every fixed `lambda>0`.  Send `lambda` to infinity. `square`

Apply NP.2 with `X_r=L_(epsilon,B)` and either `mu_r=U_r` or
`mu_r=U_r^E`.  The conference theorem and (NP.8) verify (NP.18).  Therefore

```math
\boxed{
\forall\lambda>0:\ {\mathcal R^E_{\lambda,r}\over r}\to h_\beta
\quad\Longleftrightarrow\quad
\forall\delta>0:\
U_r^E\{L\le(h_\beta-\delta)r\}=e^{-\omega_\delta(r)}.}
\tag{NP.25}
```

Here `omega_delta(r)/r -> infinity`.  In particular, a speed-`r^2` lower
tail would suffice, while a genuine speed-`r` lower tail can create a
finite-disorder-temperature phase transition.

## 4. Exact Renyi/chaos form of the same obligation

For the unconditioned output law, put

```math
a_r=\mathbb E_{U_r}e^{L},
\qquad
\Pi_r(o)={e^{L_o}\over|\Omega_r|a_r}.
\tag{NP.26}
```

The exact noisy-rank-one output identity from the archive gives

```math
\boxed{
\mathcal R_{\lambda,r}
=\log a_r-D_{1+\lambda}(U_r\Vert\Pi_r).}
\tag{NP.27}
```

Moreover `log(a_r)/r -> h_beta`.  Hence Theorem NP.2 is equivalently

```math
\boxed{
\forall\lambda>0:\quad
{1\over r}D_{1+\lambda}(U_r\Vert\Pi_r)\longrightarrow0}
\tag{NP.28}
```

if and only if the superexponential lower-tail statement in (NP.25) holds
without conditioning.  This is the exact moment/chaos formulation.  It is
not enough to bound the ordinary reverse KL (`lambda downarrow0`): a finite
phase transition is precisely the onset of an extensive higher-order
reverse Renyi divergence.

## 5. If a speed-`r` LDP exists, its slope is the phase boundary

The previous equivalence does not assume an LDP.  If, additionally,
`L/r` obeys a lower-tail LDP at speed `r` with good rate `I_beta`, unique
zero `I_beta(h_beta)=0`, then the Laplace principle gives

```math
\boxed{
\lim_r{\mathcal R_{\lambda,r}\over r}
=\inf_{a\ge0}\left\{a+{I_\beta(a)\over\lambda}\right\}.}
\tag{NP.29}
```

The typical branch survives exactly when

```math
I_\beta(a)\ge\lambda(h_\beta-a)
\qquad(0\le a<h_\beta).
\tag{NP.30}
```

Thus its critical tilt is

```math
\boxed{
\lambda_{\rm ann}(\beta)
=\inf_{0\le a<h_\beta}
 {I_\beta(a)\over h_\beta-a}.}
\tag{NP.31}
```

If `I_beta(a)=infinity` below `h_beta` at speed `r`, then
`lambda_ann=infinity`; this is exactly the speed-faster-than-`r` case of
Theorem NP.2.  If some fixed lower deviation has a finite speed-`r` rate,
then `lambda_ann` is finite.

For the same-temperature target `tau_beta=h_beta-gamma(beta)`, the exact
criterion is

```math
\inf_{a\ge0}\left\{a+{I_\beta(a)\over\lambda}\right\}
\le\tau_\beta.
\tag{NP.32}
```

In particular, only `a<=tau_beta` can realize it.  The isolated algebraic
double in the archive has probability `exp(-Theta(r^2))` for fixed child
representatives, so it supplies no finite value of the speed-`r` rate and
does not locate `lambda_ann`.

Equations (NP.25), (NP.28), and (NP.31) are three exact versions of the same
missing theorem: bridge-disorder lower-tail, output reverse Renyi, and rate-
function slope.

## 6. Audit of the March 2026 spin-glass large-deviation theorem

Chen, Guionnet, Ko, Lacroix-A-Chez-Toine, and Mourrat,
[*One-sided large deviations for the ground-state energy of spin
glasses*](https://arxiv.org/abs/2603.06368) (2026), prove for Gaussian mixed
`p`-spin Ising Hamiltonians an explicit speed-`N` LDP for deviations of the
normalized maximum **above** its typical value.  Their proof starts from a
Parisi formula for positive fractional moments `E[Z_N^s]`, `0<s<1`, and
then passes to a positive Laplace transform of the maximum.

The paper explicitly separates the opposite direction: with zero external
field, deviations of the maximum **below** its typical value are expected at
speed `N^2`; it cites a proof for spherical models, not for the Ising model.
That phenomenology is exactly favorable to (NP.25), but it is not an
applicable theorem here.

There are four non-removable hypothesis/direction mismatches.

1. Their disorder is a full Gaussian mixed-`p`-spin process.  Here only the
   bipartite cross block is random, it is Rademacher, and both deterministic
   conference blocks remain present.
2. Their observable is the zero-temperature maximum.  Here it is the
   finite-temperature normalized `cosh` pressure.
3. Their theorem controls the upper disorder tail of the maximum.  We need
   the lower disorder tail (NP.25).
4. Their fractional moment is a positive power `Z^s`.  Equation (NP.27)
   requires negative powers `Z^(-lambda)`, or equivalently reverse Renyi
   divergence.

Gaussian interpolation and the Parisi formula therefore cannot be imported
by changing signs or temperatures.  Nor does a spectral LDP alone suffice:
rare bridges can violate the open-walk/delocalization hypotheses under which
conference pressure is a spectral functional.  No verified primary theorem
found in this audit supplies a speed-`r^2` lower LDP for the deformed
bipartite Rademacher Ising pressure.

The precise imported-theorem target suggested by the paper is now:

> **Conference bridge lower-LDP lemma.**  For every fixed `beta` in (NP.4)
> and every `delta>0`, prove
>
> ```math
> \Pr_{\epsilon,B}\{L_{\epsilon,B}
>       \le(h_\beta-\delta)r\}
> \le\exp\{-r\,\omega_{\beta,\delta}(r)\},
> \qquad \omega_{\beta,\delta}(r)\to\infty.
> \tag{NP.33}
> ```
>
> The conjectural spin-glass scale would strengthen the right side to
> `exp(-c_(beta,delta)r^2)`.

By Theorem NP.2, (NP.33) is exactly sufficient for every fixed bridge tilt
to remain annealed.  It is a lower-deviation theorem for one explicit random
conference-bipartite model, not the original signing optimization.

## 7. Archive comparison and judgment

1. Equation (NP.12) and its unconditioned consequence already appear as
   (5.5)--(5.6) in `artifacts/two_temperature_bridge_audit.md`.  The new
   theorem-level increment in NP.1 is to carry the **full explicit interval**
   through the later thin-tail conditioning, rather than only the smaller
   exact-rate interval in the frozen CT theorem.
2. The exact Renyi identity (NP.27) is archived.  The new synthesis is the
   equivalence NP.2/NP.25: all finite Renyi orders are subextensive if and
   only if every lower pressure deviation is superexponential at speed `r`.
3. The LDP formula (NP.29) is a standard Laplace principle, used here to
   locate the exact phase slope and to prevent a vague request for “better
   concentration.”
4. The March 2026 paper supports the **direction and conjectural speed** of
   NP.33, but proves none of its hypotheses in this model.  It should be
   cited as motivation, not as progress on the conference lemma.

**Stopping judgment.**  The fixed-small-tilt obstruction is now extended to
the explicit interval (NP.10), and the all-fixed-tilt question is reduced
exactly to (NP.33).  Further bounded-difference optimization cannot cross
the endpoint (NP.10).  Continue only through a genuine lower-LDP method for
deformed bipartite Rademacher pressure (fractional negative moments,
enumeration of low-pressure bridge types, or a theorem controlling rare
traffic-irregular outputs).  Generic bridge sampling or another scalar
pressure rearrangement is closed.
