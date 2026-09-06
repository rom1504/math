# Independent primary-theorem audit of the marked/unmarked AMP coupling

Date: 2026-09-05. Status: the stated finite construction passes, subject to
the matrix and regularity hypotheses below. This audit does not extend it
to every bounded-operator signing.

Final archival scope update: the coupling audited here is restricted to
the imported all-power delocalization conditions, including conference
sequences. It is not a prerequisite for the subsequent universal
actual-sign weighted projection/gain theorem, and that later theorem
does not remove the hypotheses from this full AMP energy formula.
This audit imports no additional universal lower bound or original
minimum convergence claim. It does not certify the separate numerical
optimization decimal in the coupling note.

## 1. Primary statements checked directly

I read [Wang--Zhong--Fan, Proposition 2.7(b2), Theorem 2.8, and Assumption
2.1](https://arxiv.org/pdf/2206.13037) directly, including the displayed
condition (2.9). Uniform signed-permutation conjugation of a deterministic
matrix is covered when every fixed power has diagonal entries within
`n^(-1/2+epsilon)` of its normalized trace and off-diagonal entries of that
order. The limiting spectral law and bounded operator norm are additional
hypotheses. Theorem 2.8 requires continuous polynomial-growth responses,
history-Lipschitzness, admissible independent initialization/side information,
and nonsingular queried Gaussian covariance matrices.

I also read [Fan, Section 4.1, equations (4.4)--(4.7)](https://arxiv.org/pdf/2008.11892)
directly. With the response Jacobian lower triangular, the Onsager
coefficient matrix is the transpose of the free-cumulant polynomial in
that Jacobian. The covariance is the sum over both left and right Jacobian
powers. These are precisely the conventions used below.

Conference matrices satisfy the required fixed-power condition: normalized
even powers are the identity and odd powers equal the original flat hollow
matrix. Its spectral mean is zero, variance one, and fourth free cumulant
is minus one. Signed permutation conjugation preserves the Boolean
quadratic absolute norm exactly.

## 2. Earlier fields and the complete final block calculation

Let the earlier input family consist of `S a_j(X)` and `h_l(X)`, where the
`a_j` are triangular and orthonormal, the `h_l` are orthonormal and
orthogonal to all Gaussian linear coordinates, and `S` is an independent
sign. The earlier input Gram matrix is exactly the identity. Every earlier
mean derivative is zero: the marked inputs have a factor of `S`, and the
unmarked inputs have zero derivative means by Gaussian integration by
parts. The prescribed earlier Gaussian covariance is therefore the identity,
independent of the initial sign and every auxiliary sign. This checks the
induction used to produce the finite history.

Append `f=F(X,Z)+S H(X,Z)`, and put `alpha=E grad F`,
`c=E[u_prior f]`, and `d0=E(F²+H²)`. The full limiting matrices are

\[
 \Delta=\begin{pmatrix}I&c\\c^T&d_0\end{pmatrix},
 \qquad
 \Phi=\begin{pmatrix}0&0\\\alpha^T&0\end{pmatrix}.
\]

Thus `Phi²=0`. With spectral mean zero and variance one, the exact Fan
prescriptions reduce to

\[
 \operatorname{Onsager}=\Phi,
 \qquad
 \Sigma=\Delta+\kappa_3(\Phi\Delta+\Delta\Phi^T)
                    +\kappa_4\Phi\Delta\Phi^T.
\]

The last output `y=Wf-sum alpha_j u_j` has cross covariance
`Cov(y,(X,Z))=c+kappa3 alpha`. Its final variance is
`d0+2kappa3 alpha dot c+kappa4||alpha||²`. Its Schur complement against
the previous identity block is exactly

\[
 d_0-\|c\|^2+(\kappa_4-\kappa_3^2)\|\alpha\|^2.
\]

There is no missing higher-cumulant term: every such term contains at least
two successive Jacobian factors on one side and is zero. In particular,
nilpotence does not remove the `kappa4 Phi Delta Phi^T` variance term,
which is retained above.

## 3. Energy and nonsingularity

Joint Gaussian integration by parts gives
`E f y=alpha dot(c+kappa3 alpha)`; independence of `S` makes its contribution
to this expression zero. The explicit Onsager part contributes `alpha dot c`.
Therefore

\[
 {1\over2n}E f^T Wf\longrightarrow
 \sum_j E\partial_{X_j}F\,E[a_jH]
 +\sum_l E\partial_{Z_l}F\,E[h_lF]
 +{\kappa_3\over2}\|E\nabla F\|^2.
\]

This confirms the coupling artifact's formula (9), including its factor
one-half. Replacing `B` by `-B` changes the sign of `kappa3` but not the
preceding standard Gaussian model, so the favorable sign is available for
an absolute quadratic norm.

The earlier inputs and all Gaussian coordinate linears are mutually
orthonormal in scalar `L²`. Bessel's inequality gives
`d0≥||c||²+||alpha||²`. For a centered spectral variable of variance one,
`E(lambda²-m3 lambda-1)²≥0` gives `m4≥1+m3²`, hence
`kappa4-kappa3²≥-1`. The Schur complement is nonnegative.

If it vanishes, add `epsilon T` only to the final query, where `T` is a new
independent auxiliary sign. This leaves `alpha,c` unchanged and increases
`d0` by `epsilon²`. All queried covariances are then nonsingular. The
bounded operator norm controls the final energy difference by
`O(epsilon+epsilon²)`, so the original cube-valued `f` is recovered after
the dimension limit. Query regularization need not itself lie in the cube.

## 4. Regularity and scope clarifications

When defining responses on all real auxiliary arguments, extend the signs
by a bounded continuous clipping function that agrees with the identity at
`+1,-1`. Literal multiplication by an arbitrary real auxiliary argument
would not be uniformly Lipschitz in the history variables. Clipping changes
no actual iterate and supplies the required uniform Lipschitz extension.

Finite smooth approximations of the old triangular Hermite inputs may be
orthonormalized in their chronological order; subtracting earlier inputs
preserves triangular dependence. Every approximation and history length is
fixed before the matrix dimension tends to infinity. The final bounded
response and the operator bound provide uniform integrability for the
normalized energy, so the stated expectation limit follows from empirical
state evolution.

The independent audit found no flaw in this finite coupling. The unresolved
restriction is substantive: bounded operator norm, or even bounded operator
norm plus ordinary entry coherence, does not alone establish the imported
all-power delocalization condition. No universal signing lower bound or
convergence statement is obtained by dropping that condition.
