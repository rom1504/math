# All-order realization of bounded-operator signed action limits

The purification theorem in
`concentration_compactness_boolean_profiles.md` reduces minimization, up to
an arbitrarily small leading error, to signings \(A_n\) with
\[
\|A_n\|_{\rm op}=O(\sqrt n).
\]
For these sequences the normalized operators \(T_n=A_n/\sqrt n\) are
action-compact and the Boolean objective is continuous.  To prove
convergence one would still need to realize every subsequential extremal
limit at every sufficiently large order.

## 1. Principal sampling only handles nearby orders

If \(S\) has size \(m=\alpha n\), then
\[
\frac{A[S]}{\sqrt m}
=\alpha^{-1/2}
P_S\frac A{\sqrt n}P_S.
\]
Thus a proportional principal sample is a rescaled compression, not the
same action object.  Sampling preserves the limit only when
\(\alpha\to1\).  It bridges \(o(n)\) gaps but not lacunary subsequences or a
fixed-ratio change of order.

## 2. Exact variance obstruction to sign block lifts

Suppose \(m=kn\) and replace every base entry \(a_{ij}\) by a
\(k\times k\) sign block \(L_{ij}\).  To preserve the action of \(A/\sqrt n\)
on functions constant on every fiber, each row of \(L_{ij}\) must have sum
\[
r_{ij}=(1+o(1))a_{ij}\sqrt k.
\tag{1}
\]
Assume the row sum is exact, and decompose
\[
L_{ij}=\frac{r_{ij}}kJ_k+R_{ij},
\qquad R_{ij}\mathbf1=0.
\]
Orthogonality of the two terms in Frobenius inner product gives
\[
\boxed{
\|R_{ij}\|_F^2
=\|L_{ij}\|_F^2-r_{ij}^2
=k^2-r_{ij}^2
=(1-o(1))k^2.}
\tag{2}
\]
Almost the entire sign-block energy is therefore forced into
high-frequency fiber modes.

This residual cannot be corrected away:

- independent residuals create a new Wigner/Gaussian action component;
- Hadamard or orthogonal residuals create a tensor-fiber component and can
  enlarge the Boolean objective;
- constant blocks have no residual but amplify the base action by
  \(\sqrt k\), because their row sum is \(k\), not \(\sqrt k\).

Equation (2) is an exact obstruction to every naive rational blow-up,
sampling-plus-discrepancy correction, or low-energy block replacement.
An all-order realization theorem would have to show that an extremal limit
is *absorbing* under this forced microscopic sign noise.

## 3. Arithmetic obstruction was not found

Exact conference matrices have order restrictions, but those restrictions
do not yield a non-realizable action limit: principal blocks of Paley
conference matrices at orders \(N=n+o(n)\) remove parity/congruence defects
at vanishing action cost.  Similar finite divisibility constraints wash out
under asymptotic approximation.

No explicit uniformly \(2\to2\)-bounded signed action limit was found which
is provably realizable along one infinite order set and bounded away from
all sign matrices on another infinite order set.

## 4. Verdict

The realization step remains open, not disproved.  Its precise missing
statement is:

> Every extremal bounded-operator action limit arising from symmetric
> off-diagonal sign matrices is approximable, with objective convergence,
> by such sign matrices at every sufficiently large order.

The variance identity (2) shows why existing graphon sampling, rational
blow-ups, independent random lifts, and Hadamard lifts do not prove it.

