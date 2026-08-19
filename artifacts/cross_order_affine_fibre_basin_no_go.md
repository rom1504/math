# A successful switch quotient requires a thick affine good basin

Status: **proved exact all-order characterization and conditional scalable
no-go for actual optimizing children**.  A labelled Fourier quotient is not
merely a truncated Fourier estimate: its certificate is exactly the
logarithmic mean of the true fixed-child parent pressures over one affine
switch fibre.  Consequently a quotient of effective dimension
`O(sqrt(N))` can give a sublinear cross-order defect only if essentially an
entire fibre, of size `exp((log 2)N-O(sqrt(N)))`, consists of individually
sublinear-good switches.

This theorem does not prove that such fibres exist or do not exist for
large optimizing children.  It sharpens the surviving labelled-quotient
route into a basin-entropy dichotomy and gives an immediate quantitative
arrow to the permanent cross-order defect.

## 1. Exact effective quotient

Retain the notation of
[`cross_order_labelled_fourier_quotient_no_go.md`](cross_order_labelled_fourier_quotient_no_go.md).
Thus `N=m+n`, `t=beta/sqrt(N)`, `A,D` are exact own-scale pressure
minimizers, and

```math
 P=P_m(\beta)+P_n(\beta),\qquad
 X_g=L_{\epsilon,g}(B)-P.                         \tag{1.1}
```

Normalize the actual child and bridge kernels by

```math
 a={w_\epsilon\over\mathbb Ew_\epsilon},\qquad
 b={k_B\over\mathbb Ek_B},\qquad f=a*b,           \tag{1.2}
```

and put

```math
 \mathcal C_{\epsilon,B}
 =\log\mathbb Ew_\epsilon+\log\mathbb Ek_B-P.     \tag{1.3}
```

Let `W` be the `(N-2)`-dimensional even--even character space.  The exact
Fourier formula for the child kernel gives

```math
 \widehat f(\chi)=0\quad(\chi\notin W).            \tag{1.4}
```

For an arbitrary character subspace `V`, set

```math
 U=V\cap W,\qquad r=\dim U.                        \tag{1.5}
```

The number `r`, rather than the nominal dimension of `V`, is the effective
quotient dimension.

**Theorem 1.1 (affine-fibre formula).**  Conditional expectation onto `V`
retains exactly the frequencies in `U`, so, after lifting quotient functions
back to the switch group,

```math
 f_V=f_U.                                          \tag{1.6}
```

If `F_u=pi_U^(-1)(u)`, then

```math
 |F_u|=2^{N-1-r}                                   \tag{1.7}
```

and

```math
 \boxed{
 \mathcal C_{\epsilon,B}+\log f_U(u)
 =\log\left{{1\over|F_u|}\sum_{g\in F_u}e^{X_g}\right}.}
                                                               \tag{1.8}
```

In particular the exact quotient certificate is

```math
 \boxed{
 Q_{V,\epsilon,B}
 =\min_u\log\left{{1\over|F_u|}
            \sum_{g\in F_u}e^{L_{\epsilon,g}(B)-P}\right}.}  \tag{1.9}
```

*Proof.*  Fourier conditional expectation keeps precisely the characters
in `V`.  Equation (1.4) deletes those outside `W`, proving (1.6).  Every
fibre of the rank-`r` homomorphism `pi_U` has the cardinality (1.7).  The
exact switch-convolution identity gives

```math
 e^{L_{\epsilon,g}(B)}
 =\mathbb Ek_B\,\mathbb Ew_\epsilon\,f(g).         \tag{1.10}
```

Average (1.10) over `F_u`, take logarithms, and subtract `P`; this is
(1.8), hence (1.9).  `square`

Parent minimization and the fact that a minimum is at most an exponential
mean give the direct defect arrow

```math
 \boxed{
 Q_{V,\epsilon,B}\le\omega_N
 \quad\Longrightarrow\quad
 E_{m,n}(\beta)\le\omega_N.}                       \tag{1.11}
```

Thus the concrete sufficient lemma

```math
 {1\over|F|}\sum_{g\in F}e^{L_{\epsilon,g}(B)-P}
 \le e^{C_\beta N^{1-\delta}}                     \tag{1.12}
```

for one actual-child affine fibre immediately implies

```math
 \boxed{E_{m,n}(\beta)\le C_\beta N^{1-\delta}.}   \tag{1.13}
```

There is no unlinked intermediate state in (1.12)--(1.13).

## 2. Quotient success forces a thick good fibre

**Theorem 2.1 (affine good-basin necessity).**  If
`Q_(V,epsilon,B)<=omega_N`, then, for every `s>0`, a fibre `F` of the
effective quotient satisfies

```math
 \boxed{
 \#\{g\in F:X_g\le\omega_N+s\}
 \ge(1-e^{-s})2^{N-1-r}.}                          \tag{2.1}
```

*Proof.*  Choose the minimizing fibre in (1.9).  Its exponential mean is at
most `e^(omega_N)`.  If a fraction `p` of the fibre had
`X_g>omega_N+s`, that mean would exceed `p e^(omega_N+s)`.  Hence
`p<=e^(-s)`.  `square`

For a fixed bridge and orientation define the directly relevant count

```math
 M_{B,\epsilon}(u)
 =\#\{g:L_{\epsilon,g}(B)-P\le u\}.                \tag{2.2}
```

Combining (1.11) and (2.1) gives

```math
 \boxed{
 Q_{V,\epsilon,B}\le\omega_N
 \Longrightarrow
 \begin{cases}
 E_{m,n}(\beta)\le\omega_N,\\[2mm]
 r\ge N-1-\log_2M_{B,\epsilon}(\omega_N+s)
       +\log_2(1-e^{-s}).
 \end{cases}}                                      \tag{2.3}
```

Two scalable consequences are immediate.

1. If `r<=K sqrt(N)`, `omega_N=o(N)`, and
   `s_N->infinity`, `s_N=o(N)`, then quotient success requires

   ```math
   M_{B,\epsilon}(o(N))
   \ge(1-o(1))2^{N-1-K\sqrt N}.                    \tag{2.4}
   ```

   Thus a critical quotient needs an almost full-entropy affine flat of
   individually good switches.  One exceptional switch or even
   `exp(o(N))` good switches cannot suffice.

2. If, uniformly over the bridge/orientation family under consideration,

   ```math
   M_{B,\epsilon}(\omega_N+s)
   \le e^{\gamma N+o(N)},\qquad\gamma<\log2,       \tag{2.5}
   ```

   then every successful quotient has

   ```math
   \boxed{r\ge(1-\gamma/\log2)N-o(N).}             \tag{2.6}
   ```

   In particular an `exp(o(N))` good landscape forces
   `r=N-o(N)`, essentially the full `(N-2)`-dimensional even--even
   landscape.

If the pre-cancellation term also has a positive linear floor
`C_(epsilon,B)>=cN-o(N)`, then every switch counted in (2.4) satisfies

```math
 \log f(g)\le-cN+o(N).                              \tag{2.7}
```

The needed cancellation is therefore a linear depression across almost an
entire affine fibre, not a localized Fourier valley.

## 3. Exhaustive actual-child wind tunnel

[`audit_switch_global_quotient_basin.py`](../computations/audit_switch_global_quotient_basin.py)
exhausts all `512` row/column-switching classes of order-four bridges, both
orientations, and every effective character subspace of the requested
dimension.  It then validates the winning quotient by directly averaging
the actual parent pressures on its decoded affine fibre.

For the exact order-four pressure minimizer, the global finite results are:

| `beta` | effective `r` | best certificate | bridge mask | fibre anatomy |
|---:|---:|---:|---:|:---|
| 4 | 3 | `-0.056514042676` | 245 | 8 switches at `-1.090324951270`, 8 at `0.440831717307` |
| 8 | 3 | `0.581268468549` | 245 | 8 at `-3.015191152099`, 8 at `1.260610479637` |
| 8 | 4 | `-3.015191152262` | 245 | all 8 switches at `-3.015191152099` |

Here each listed fibre value is `X_g=L_(epsilon,g)-2P_4`; the representative
orientation is `epsilon=+1`.  At `beta=4`, the critical dimension
`r=3` is already successful because good switches fill half of a
16-element affine fibre.  At `beta=8`, the same fibre's log-mean-exp is
positive; one additional effective bit isolates an eight-element fibre on
which every switch is optimal.  These are exhaustive floating-point finite
diagnostics, not asymptotic claims or interval-certified transcendental
inequalities.

## 4. Scope and cross-order verdict

Theorems 1.1 and 2.1 apply to every finite order, every bridge and
orientation, and the actual selected pressure-minimizing children.  They
prove a genuine basin-versus-information dichotomy for the only surviving
mesoscopic Fourier window.

They do **not** prove either side of that dichotomy uniformly at large
orders: no theorem currently supplies (1.12), and no theorem bounds
`M_(B,epsilon)` by `exp(o(N))` for every bridge selected from actual
children.  Consequently the strongest unconditional comparable-order
defect remains `O_beta(N)`.  The result is a conditional method-class no-go,
not an exponent improvement and not a lower bound on the true optimized
parent defect.
