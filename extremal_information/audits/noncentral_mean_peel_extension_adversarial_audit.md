# Adversarial audit: noncentral mean-peel extension

**Proposal audited.**  Let `mu_r` be an arbitrary law on `{+-1}^r` with
density `g_r=dmu_r/dU_r` satisfying `E_U g_r^2<=K`.  Put

```math
m_r=\mathbb E_{\mu_r}R,
\qquad \Sigma_r=\mathbb E_{\mu_r}RR^T.
```

Peel the bad eigenspace of `Sigma_r` at window
`[1-r^(-1/4),1+r^(-1/4)]`, join it with `span(m_r)`, and project rows onto
the orthogonal complement `V_r`.

**Verdict:** **PASS.**  Central symmetry is not needed after the mean
direction is explicitly peeled.  The projected iid-row operator edge is
still two, joining the mean cannot destroy the spectral window, and convex
restoration pays `o(r)` pressure.  I found no biased-halfcube, mixture, or
high-order-fibre counterexample.

The same conclusion holds if the bad eigenspace is defined from the centered
covariance `C_r=Sigma_r-m_rm_r^T`, provided its regular compression is the
one used for whitening.  On `m_r^perp`, the compressions of `C_r` and
`Sigma_r` coincide.

## 1. Fourier information and peel rank

Parseval gives separately

```math
\|m_r\|_2^2
=\sum_i\langle g_r,\chi_i\rangle_U^2\le K-1,
```

and

```math
\|\Sigma_r-I\|_F^2
=2\sum_{i<j}\langle g_r,\chi_{ij}\rangle_U^2
\le2(K-1).
```

Thus the second-moment bad eigenspace has rank `O_K(sqrt r)` at threshold
`delta_r=r^(-1/4)`.  Joining `span(m_r)` increases rank by at most one.

If one instead peels the centered covariance, then

```math
\|C_r-I\|_F
\le\|\Sigma_r-I\|_F+\|m_rm_r^T\|_F
\le\sqrt{2(K-1)}+(K-1),
```

so its bad eigenspace has the same `O_K(sqrt r)` rank.

## 2. Joining the mean preserves the spectral window

Let `P_0` be the bad spectral projection of `Sigma_r`, let `P_r` project
onto `ran(P_0)+span(m_r)`, and set `V_r=I-P_r`.  Although `V_r` need not
commute with `Sigma_r`, its range is a subspace of the original good
spectral space.  Hence every `x in ran(V_r)` satisfies

```math
(1-\delta_r)\|x\|_2^2
\le x^T\Sigma_rx
\le(1+\delta_r)\|x\|_2^2.
```

Therefore the compression

```math
S_r=V_r\Sigma_rV_r|_{\operatorname {ran}V_r}
```

still has all eigenvalues in the desired window.  This Rayleigh-quotient
argument does not require `V_r` to be a spectral projection itself.

Moreover `V_rm_r=0`, so

```math
X_r=S_r^{-1/2}V_rR
```

is centered and exactly isotropic.  The lack of central symmetry has been
removed at the only place where Chafaï--Tikhomirov requires centering.

## 3. The projected edge remains two

For every whitened rank-`k` projection `Q`, define as before

```math
M=V_rS_r^{-1/2}QS_r^{-1/2}V_r.
```

Then

```math
\|M\|_{op}\le(1-\delta_r)^{-1},
\quad \|M\|_F^2\le k(1-\delta_r)^{-2},
\quad |\operatorname {tr}M-k|
\le{\delta_r\over1-\delta_r}k.
```

These facts use only the compression `S_r`; commutation with `Sigma_r` is
unnecessary.  Indeed

```math
\mathbb E_{\mu_r}R^TMR
=\operatorname {tr}(\Sigma_rM)=\operatorname {tr}Q=k.
```

Cauchy--Schwarz transfers every unconditioned Hanson--Wright event:

```math
\mu_r(A)\le\sqrt K\,U_r(A)^{1/2}.
```

Consequently the proof of uniform STP and Yaskov (A1) in Corollary CE.2 is
unchanged.  With `d_r=rank(V_r)=r-O_K(sqrt r)`, the iid row matrix obeys

```math
{\|B_rV_r\|_{op}\over\sqrt r}\longrightarrow2
```

in probability.  Thus no projected-edge counterexample exists within the
bounded-`L^2` class.

## 4. Transport and restoration costs

Jensen under `mu_r` gives

```math
D(\mu_r\|U_r)\le\log\mathbb E_Ug_r^2\le\log K.
```

The Hamming transport argument therefore couples a `mu_r` row to a uniform
row with expected edit distance `O_K(sqrt r)`, without any symmetry
assumption.  For independent row copies,

```math
\mathbb E\|(B_r-W_r)V_r\|_F=O_K(r^{3/4})=o(r).
```

Let `k_r=rank(P_r)=O_K(sqrt r)`.  Since
`tr Sigma_r=E||R||_2^2=r`,

```math
\mathbb E\|B_rP_r\|_*
\le\sqrt{k_r}\,(r\operatorname {tr}(P_r\Sigma_r))^{1/2}
\le r\sqrt{k_r}=O_K(r^{5/4})=o(r^{3/2}).
```

The iid removed component costs only `k_r sqrt(r)=O_K(r)`.  The repaired
projected-coupling theorem now applies.  Its convex supporting-line argument
restores the full `B_rP_r`, including the mean direction, with expected
downward pressure loss

```math
O(r^{-1/2})\,\mathbb E\|B_rP_r\|_*=O_K(r^{3/4})=o(r).
```

The mean itself is even cheaper: `||m_r||=O_K(1)`, so the deterministic
rank-one mean matrix `1m_r^T` has nuclear norm `O_K(sqrt r)`.

## 5. Explicit hostile laws

1. **Pinned-coordinate halfcube.**  Conditioning on `R_1=1` has `K=2`,
   `m=e_1`, and `Sigma=I`.  Peeling `e_1` leaves an ordinary iid Rademacher
   block, whose edge is exactly two.  The restored pinned column has nuclear
   norm `sqrt r`.
2. **Majority halfcube.**  For odd `r`, conditioning on positive total sum
   selects one point from every antipodal pair.  Every even observable,
   including `RR^T`, retains its uniform expectation, so `Sigma=I`, while
   `m` is a constant-norm multiple of the all-ones direction.  The mean peel
   again removes the only noncentral obstruction.
3. **Diffuse product bias.**  Taking independent coordinate means
   `m_i=Theta(r^{-1/2})` has bounded `L^2` density when
   `sum_i m_i^2=O(1)`.  Its finite population spike is caught by the
   second-moment/mean peel; the bulk remains regular.
4. **Mixtures and high-order fibres.**  These can create arbitrary
   high-order dependence, but bounded `L^2` likelihood still transfers every
   quadratic Hanson--Wright tail with only a square-root probability loss.

None produces either an excess projected edge or a leading downward
restoration cost.

## 6. Minimal statement discipline

The extension should specify that:

- `P_r` is the orthogonal projection onto the *sum* of the mean span and the
  chosen bad eigenspace;
- `S_r` is the compression `V_rSigma_rV_r` rather than an assumed invariant
  restriction;
- whitening uses this compression; and
- rows of `B_r` are independent copies of `mu_r`.

Under those explicit conventions, the noncentral extension is rigorous.

## 7. Frozen implementation verification

I verified the implementation
`extremal_information/drafts/bounded_l2_noncentral_row_extension.md` at
SHA-256
`94d9d676f2d348584093b79f333bd503a2f10830111a6871a30cdb2117178a46`.

**Verdict:** **PASS with one non-theorem prose repair.**  BL.1--BL.6 and
Corollary BL.3 are rigorous and match Sections 1--6 above.  The source even
improves the generic removed-component estimate: BL.11 gives
`||Sigma_r||op<=L_K`, hence

```math
\operatorname {tr}(P_r\Sigma_r)\le L_K k_r
```

and therefore

```math
\mathbb E\|B_rP_r\|_*
\le\sqrt{k_r}(rL_Kk_r)^{1/2}
=O_K(r),
```

rather than merely `O_K(r^(5/4))`.  This is correct even though `P_r` does
not commute with `Sigma_r`.

The uniform quantifier needed by BL.3 is present.  If

```math
\epsilon_r(K)=
\sup_{\mu:\,\mathbb E_U(d\mu/dU)^2\le K}
\mathbb E_{\mu^{\otimes r}}
 \left(h_\beta-f_r(B)/r\right)_+
```

did not tend to zero, selecting an almost-maximizing law at each violating
order would form a triangular sequence contradicting BL.2.  Thus every
conditional component in BL.42 has shortfall at most the same
`epsilon_r(K)`, and Tonelli/integration proves BL.44.  The component
projection may depend on the latent state because it is only a proof
witness, not part of the bridge law or the integrated observable.

The one prose repair is the sentence after BL.3 saying the corollary can
cover `Theta(r) or larger` row total correlation.  Under the standard
multi-information definition and fixed `K`, superlinear total correlation
is impossible.  If `bar mu_r=int mu_(r,z) pi_r(dz)`, then

```math
\begin{aligned}
D(q_r\|U_r^{\otimes r})
&=D(q_r\|\bar\mu_r^{\otimes r})
  +rD(\bar\mu_r\|U_r),\\
D(q_r\|U_r^{\otimes r})
&\le\int D(\mu_{r,z}^{\otimes r}\|U_r^{\otimes r})\,\pi_r(dz)\\
&=r\int D(\mu_{r,z}\|U_r)\,\pi_r(dz)
\le r\log K.
\end{aligned}
```

Hence

```math
D(q_r\|\bar\mu_r^{\otimes r})\le r\log K.
```

The bound can be `Theta(r)` (for instance, one latent bit choosing between
two disjoint uniform halfcubes), so the intended and correct wording is
“linear row total correlation, and arbitrarily large latent support.”  This
does not alter Corollary BL.3 or any pressure conclusion.

### Final source confirmation

The final source with SHA-256
`e4e9a9e83e369bafabe3896e98efb8a95e8e9d49f4f70778d5ddf5b57568e282`
replaces the overstatement by “linear row total correlation and arbitrarily
large latent support.”  I verified that repair and rechecked the surrounding
BL.3 statement.  **Final verdict: PASS.**
