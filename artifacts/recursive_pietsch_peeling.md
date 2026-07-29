# Recursive Grothendieck--Pietsch peeling: exact output and flat-core wall

Let \(S\) be a current principal block of size \(m\), write
\[
B_S=\|A[S]\|_{\infty\to1},
\]
and let \(K_G\) denote the real Grothendieck constant.

## 1. Symmetric Pietsch weight and one peeling step

Grothendieck factorization gives probability weights \(\mu,\nu\) with
\[
|u^\top A[S]v|
\le K_GB_S\|u\|_{L_2(\mu)}\|v\|_{L_2(\nu)}.
\]
Putting \(\rho=(\mu+\nu)/2\) and using
\(\mu,\nu\le2\rho\) coordinatewise gives the symmetric form
\[
\boxed{
|u^\top A[S]v|
\le2K_GB_S\|u\|_{L_2(\rho)}\|v\|_{L_2(\rho)}.}
\tag{1}
\]

Fix \(\tau>1\), and split
\[
H=\{i\in S:\rho_i>\tau/m\},\qquad R=S\setminus H.
\]
Then
\[
|H|<m/\tau
\tag{2}
\]
and, for vectors supported in \(R\),
\[
\|u\|_{L_2(\rho)}\le\sqrt{\tau/m}\|u\|_2.
\]
Consequently
\[
\boxed{
\|A[R]\|_{\rm op}\le\frac{2K_GB_S\tau}{m}.}
\tag{3}
\]

This is a clean hub/regular decomposition.  Recursing on \(H\) produces a
nested sequence
\[
S=H_{-1}\supset H_0\supset H_1\supset\cdots,
\qquad |H_j|<m\tau^{-(j+1)},
\]
and regular layers \(R_j=H_{j-1}\setminus H_j\).

## 2. Why the resulting energy bounds do not sum at the right scale

For the large regular remainder, (3) gives only
\[
Q(A[R])
\le |R|\|A[R]\|_{\rm op}
\le 2K_GB_S\tau.
\tag{4}
\]
Since \(B_S\le2Q(A[S])\), the right side of (4) is a factor
\(4K_G\tau\) *larger* than the current Boolean budget.  Monotonicity gives
the better bound \(Q(A[R])\le Q(A[S])\), but no scale reduction at all.

The same issue affects cross blocks.  From (1), sign vectors supported on
sets \(I,J\subset S\) satisfy
\[
|x_I^\top A[I,J]y_J|
\le2K_GB_S\sqrt{\rho(I)\rho(J)}.
\tag{5}
\]
Thus the natural scalar quantity attached to a union \(U\) is
\(\rho(U)\), and the bound scales linearly:
\[
Q(A[U])\le2K_GB_S\rho(U).
\tag{6}
\]
For a cardinal fraction \(|U|/m=\alpha\) in a flat block,
\(\rho(U)=\alpha\).  Equation (6) gives \(\alpha B_S\), whereas
proportional restriction needs \(\alpha^{3/2}B_S\).  The missing factor
\(\sqrt\alpha\) is precisely cancellation information not contained in
Pietsch factorization.

## 3. Explicit obstruction: a flat orthogonal core

Suppose \(C\) is a symmetric conference matrix:
\[
C^2=(m-1)I.
\]
For all \(u,v\),
\[
|u^\top Cv|\le\sqrt{m-1}\|u\|_2\|v\|_2.
\]
Taking the uniform probability weight \(\rho_i=1/m\), this becomes
\[
|u^\top Cv|
\le m\sqrt{m-1}\,
\|u\|_{L_2(\rho)}\|v\|_{L_2(\rho)}.
\tag{7}
\]
Its bilinear Boolean norm is of order \(m^{3/2}\), so (7) is a valid
Pietsch factorization up to an absolute constant.

For every threshold \(\tau>1\), the hub set
\[
\{i:\rho_i>\tau/m\}
\]
is empty.  Hence recursive high-weight peeling terminates immediately on
the whole conference core.  If one forcibly partitions the coordinates,
the factorization still gives only
\[
Q(C[U])\lesssim B_C\,\frac{|U|}{m},
\tag{8}
\]
not the needed \(B_C(|U|/m)^{3/2}\).

This is an explicit obstruction to **pure** recursive
Grothendieck--Pietsch peeling.  It does not disprove a hierarchy that also
uses genuinely new information about the flat core (entrywise rigidity,
high even-cycle data, or an Ising-specific entropy profile), but it proves
that factor weights and scalar block masses alone cannot supply the
missing square-root cancellation.

## 4. Verdict

Grothendieck--Pietsch completely resolves the common-support deletion
problem for spectral thresholds \(L\gg\sqrt n\), as recorded in
`orientation_even_grothendieck_localization.md`.  Its recursive use has a
canonical stopping state: a diffuse \(O(\sqrt m)\)-operator-norm core.
That core already contains the original hard conference/ROM regime.
Therefore the next useful input must act on the flat core itself; further
Pietsch peeling or scalar \(W\)-partitioning cannot prove convergence.

