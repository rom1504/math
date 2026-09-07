# Independent audit: cap-controlled quadratic Gibbs fluctuations

Status: PASS by direct reconstruction. Canonical theorem:
`decisive_director_delocalized_quadratic_fluctuations_2026_09_07.md`.
This note adds a matrix/determinant consequence and records the precise
scope of the primary literature comparison. It does not prove convergence.

## 1. Conditioning and constants

Write H_J(x)=sum_{i<j}J_ij x_i x_j and Q(J)=max|H_J|. For D=diag(d_i)
with D>=+/-J, the matrices K_s=D+sJ are positive semidefinite. Starting
with the original quenched Gibbs law of (s,x), then sampling
Y|(s,x,g)~N(K_s x,K_s), gives conditionally independent spins given
(s,Y,g), with means tanh(Y_i). Diagonal addition changes both branch
partition functions by the SAME factor exp(tr(D)/2), so no unintended
branch reweighting occurs. This statement remains exact for singular K_s:
Gaussian exponential tilting translates by K_s x inside its range.

For each fixed g, variance decomposition first conditions on (s,Y).
Only after this inequality is established do we average in g. Thus the
left side is E_g Var_{mu_g}, as required for the second derivative of
quenched log partition functions. It is not necessary to add the variance
of the conditional mean over g.

If y_i=E|Y_i|, the posterior Gaussian law gives

    sum_i y_i <= tr(D)+E||Jx||_1+sqrt(2N tr(D)/pi).

The pointwise bound ||Jx||_1<=4Q(J) uses cube polarization and hollow
multilinear extension; the latter bounds |H_J(z)| by Q(J) on [-1,1]^N.
Independently, conditional expectation over x_i gives

    E[s x_i (Jx)_i]=E[(Jx)_i tanh((Jx)_i)].

Summation and |u|<=u tanh u+1 give E||Jx||_1<=2Q(J)+N.
This uses only the Gibbs conditional spin law, not optimality. Arbitrary
externally chosen branch weights leave it unchanged.

For completeness, the simultaneous-majorant SDP dual is

    max Tr[J(X-Y)] : X,Y>=0, diag(X+Y)=1.

Gram vectors (a_i,b_i) and (a_i,-b_i) are unit vectors. Real
Grothendieck and polarization bound this objective by 4 K_G Q(J), with
K_G=pi/(2 asinh(1)); finite strict SDP feasibility justifies duality.
There is no extra factor two. Thus, writing q=Q(J)/N, one may take

    K=min(4q,2q+1)+4K_G q+sqrt(8K_G q/pi),
    sum_i y_i<=NK.

## 2. Matrix-level statement and a determinant consequence

Index the edge-statistic vector by unordered pairs:

    T(x)=(x_i x_j)_{i<j},    T_s(s,x)=s T(x).

For either vector, let C be E_g Cov_{mu_g} of that vector. Conditional
product chaos and then Jensen prove the simultaneous Loewner inequality

    C >= diag_{i<j}( exp[-2(y_i+y_j)] ).                 (A)

Indeed, for any coefficient vector b, the degree-two centered product
part of sum b_ij x_i x_j has squared norm sum b_ij^2 v_i v_j, where
v_i=sech^2Y_i. It is orthogonal to the constant and linear parts. This
works for arbitrary real b and for s times that observable because s is
held fixed in this conditional variance. Jensen is applied to |Y_i|+|Y_j|;
independence of Y_i,Y_j is neither used nor asserted.

If V=sum b_ij^2 and eta=N max_i(sum_j b_ij^2)/V, (A) implies

    b^T C b >= V exp(-2 eta K).

For a complete flat m-by-n bridge, eta=(m+n)/min(m,n), exactly as in
the canonical theorem.

The additional determinant consequence is

    log det C >= -2(N-1) sum_i y_i >= -2N(N-1)K.        (B)

This follows by conjugating (A) by the inverse square root of its positive
diagonal right side and taking determinants. In dimension L=binom(N,2),
the geometric mean of the eigenvalues of C is therefore at least exp(-4K).
The assertion does NOT imply a dimension-free smallest eigenvalue: the
individual y_i can be large. It is a finite-temperature Fisher-information
volume bound for either choice of sufficient statistics.

## 3. Pressure integration and minimizing a bridge afterwards

For a fixed parent J and a prescribed cut, let J_0 be its block-diagonal
part, B its cross part, and U the diagonal sign matrix switching one block.
For 0<=t<=1,

    J_0+tB=((1+t)/2)J+((1-t)/2)UJU,
    Q(J_0+tB)<=Q(J).

Consequently the same cap-controlled curvature constant is valid along
the ENTIRE bridge-deletion path. Switching makes its quenched pressure
even in t (the global orientation field is unaffected), so p'(0)=0.
Integration yields

    p(1)>=p(0)+(1/2)c V_B.

It is legitimate to establish this separately for every eligible parent
and then minimize the endpoint over bridge signs. The endpoint p(0) is
independent of those bridge signs. This does not require curvature of the
minimum envelope; such curvature is generally not available.

For J=beta A/sqrt(N), the flat bridge has V_B=beta^2 mn/N, giving an
order-N pressure gain at comparable splits and fixed beta. This is a
genuine interacting gain. The inherited child temperatures and the common
orientation at the factorization endpoint still obstruct the desired
original-value recurrence.

## 4. Primary-source comparison

The directly relevant prior anti-concentration theorem is Theorem 4.2,
printed pp. 13--15, of Diakonikolas--Kane--Stewart--Sun,
[Outlier-Robust Learning of Ising Models Under Dobrushin's Condition](https://cseweb.ucsd.edu/~dakane/RobustIsing.pdf).
It bounds the variance of any shifted hollow quadratic form below by a
constant times its squared Frobenius norm for an (M,alpha)-bounded model.
Definition 2.1 means maximum absolute interaction-row sum at most M and
maximum external field at most alpha. Importantly, this particular lower
bound allows any FIXED M, not just the Dobrushin regime M<1. Its proof
uses two independent samples, agreement-set decoupling, and linear-form
anti-concentration. For J=beta A/sqrt(N), the row sum is
beta(N-1)/sqrt(N), so that hypothesis does not provide a uniform constant.

Theorem 1, printed pp. 3--4, of Eldan--Koehler--Zeitouni,
[A Spectral Condition for Spectral Gap](https://arxiv.org/pdf/2007.08200),
assumes a positive-semidefinite representative 0<=J<I and proves a
Poincare upper variance bound through the Glauber Dirichlet form. This is
a spectral high-temperature condition and a different inequality. The
new cap-controlled argument does not assume a bounded spectral width or
mixing. These comparisons identify why the retrieved theorems do not
directly yield the present dense conclusion; they are not a claim of an
exhaustive novelty search.

## 5. Reproducible finite checks

Run `.venv/bin/python computations/decisive_audit_quadratic_hubbard_2026_09_07.py`.
The checker exhausts switching classes at orders 3--5 at three inverse
temperatures, with and without a quenched global orientation field, then
tests random order-6--8 matrices. It computes the complete edge-statistic
covariance, subtracting its mean at each Gaussian quadrature node BEFORE
averaging. The Gaussian posterior absolute means y_i are computed by the
closed normal absolute-moment formula, not by sampling.

Result: PASS, 492 Gibbs systems, 984 full quadratic covariance/determinant
tests, and 3984 complete-cut directions. These floating diagnostics check
normalization and conditioning; the proof above supplies the theorem.
