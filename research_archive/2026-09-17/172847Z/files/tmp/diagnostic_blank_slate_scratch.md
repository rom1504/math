# Blank-slate convergence candidates

Timestamp: 2026-07-30T15:26:08Z (UTC)

Constraint note: this report was written before opening `ledger.md` or consulting any
route vocabulary from it.  I use only the definition

\[
D(A):=\max_{x\in\{\pm1\}^n}\left|q_A(x)\right|,
\qquad q_A(x):=\sum_{i<j}a_{ij}x_ix_j,
\qquad M_n:=\min_A D(A),
\]

and elementary union bounds/conditioning.  There are exactly three candidates.

## Candidate 1: block gluing and a (2/3)-power Fekete inequality

**Exact map.**  Given signings (A) on (m) vertices and (B) on (r)
vertices, form

\[
G_C=\begin{pmatrix}A&C\\C^T&B\end{pmatrix}
\]

(zero diagonal understood), where (C\in\{\pm1\}^{m\times r}).  For random
independent entries in (C), a union bound over the (2^{m+r}) spin pairs gives

\[
\Pr\left\{\max_{u,v}|u^TCv|>t\right\}
 \le 2^{m+r+1}e^{-t^2/(2mr)}.
\]

Consequently there is a deterministic bridge with

\[
D(G_C)\le D(A)+D(B)+
\sqrt{2mr(m+r+1)\log2},                                      \tag{1}
\]

and hence the same inequality for (M_{m+r}).  Together with the elementary
random-sign bound (M_s\le \sqrt{s(s-1)(s+1)\log2}=O(s^{3/2})), (1) implies

\[
\frac{M_{m+r}}{(m+r)^{3/2}}
\le \left(\frac m{m+r}\right)^{3/2}\frac{M_m}{m^{3/2}}
 +O\left((r/m)^{3/2}+\sqrt{r/m}\right).                       \tag{2}
\]

Thus a low value propagates forward through any (r=o(m)), but (1) has a
leading-order bridge cost when (r\asymp m).

**Candidate theorem.**  The bridge can be chosen *using (A,B)* so that, uniformly
in (m,r\to\infty),

\[
M_{m+r}^{2/3}\le M_m^{2/3}+M_r^{2/3}+o(m+r).                  \tag{C1}
\]

Equivalently, at the original scale,

\[
D(G_C)\le (D(A)^{2/3}+D(B)^{2/3})^{3/2}+o((m+r)^{3/2})
\]

for optimal (A,B) and some bridge (C).  The proposed map is exactly the
two-block matrix above; the needed improvement over (1) is a correlated bridge,
not a relabelling.  (C1), with a uniform (o(N)) defect, is an approximate
Fekete lemma for (b_n=M_n^{2/3}), and forces (b_n/n), hence
(M_n/n^{3/2}), to converge.  A precise falsifier is a fixed
\(\alpha\in(0,1)\) and infinitely many (N) such that every bridge between
optimal (or (o(N^{3/2}))-optimal) blocks of sizes \(\lfloor\alpha N\rfloor\)
and (N-\lfloor\alpha N\rfloor) exceeds the right side by \(\Omega(N^{3/2})\).

## Candidate 2: sharp downward sampling / exchange

**Exact map.**  From a signing (A) of (K_N), choose an injection
\(\phi:[n]\hookrightarrow[N]\) and use the principal signing
(A[\phi([n])]).  The elementary conditioning identity is

\[
q_{A[S]}(x)=\mathbb E\big[q_A(x,Z)\mid x\text{ fixed on }S\big],             \tag{3}
\]

where the spins (Z) off (S) are independent uniform signs.  Thus
(D(A[S])\le D(A)) for every (S), but this loses the desired
((n/N)^{3/2}) scaling.

**Candidate theorem.**  There is an absolute (C) such that every signing (A)
on (N\) vertices and every (n\le N) admit a principal (n)-set (S) with

\[
\frac{D(A[S])}{n^{3/2}}
\le \frac{D(A)}{N^{3/2}}
   +C\left(n^{-1/4}+\sqrt{n/N}\right).                       \tag{C2}
\]

The exponent (1/4) is not sacred, but a definite (o_n(1)), uniform in
(N\), is essential; a constant additive sampling error is useless.  Apply
(C2) to an optimal (A_N), send (N/n\to\infty) along a liminf subsequence,
and obtain

\[
M_n/n^{3/2}\le \liminf_N M_N/N^{3/2}+Cn^{-1/4};
\]

this forces convergence.  The proposed proof mechanism is random injection plus
an exchange/alteration step controlled by the full-family bound (D(A)).
A precise falsifier is a sequence (A_N) with (D(A_N)=O(N^{3/2})) and some
(n(N)\to\infty, n(N)/N\to0), for which every principal (n(N))-submatrix
has normalized discrepancy at least
(D(A_N)/N^{3/2}+c) for a fixed (c>0).

## Candidate 3: finite-temperature interpolation with uniform zero-temperature error

For a signing (A), put

\[
Z_A(\beta)=2^{-n}\sum_x
 \exp\left(\frac{\beta|q_A(x)|}{\sqrt n}\right),\qquad
P_n(\beta)=\frac1n\min_A\log Z_A(\beta).
\]

There is an exact, uniform zero-temperature squeeze

\[
\frac{M_n}{n^{3/2}}-\frac{\log2}{\beta}
\le \frac{P_n(\beta)}\beta
\le \frac{M_n}{n^{3/2}}.                                    \tag{4}
\]

**Candidate theorem.**  For every fixed \(\beta<\infty\), block interpolation
and a deterministic choice of the cross-edge signs give

\[
\left|(m+r)P_{m+r}(\beta)-mP_m(\beta)-rP_r(\beta)\right|
\le C_\beta\sqrt{m+r}                                       \tag{C3}
\]

(an (o(m+r)) error would suffice).  The exact composition map is again
(G_C), now choosing (C) sequentially to control the log-partition potential;
the candidate exchange step rounds each cross sign according to the first
derivative of \(\log Z\), while the sum of second-order losses must be
(o(m+r)).  Approximate additivity yields convergence of (P_n(\beta)) for
each fixed \(\beta\); then (4), first (n\to\infty) and then
\(\beta\to\infty), yields convergence of the original normalized minima.

The obvious obstruction is quantitative and therefore testable: naive random
bridging has accumulated curvature

\[
\sum_{i\le m<j}O(\beta^2/(m+r))=O_\beta(mr/(m+r))=O_\beta(m+r),
\]

not (o(m+r)).  Thus C3 has content only if signs can be exchanged/balanced so
that this extensive curvature cancels or is charged to the two block pressures.
A falsifier is a fixed \(\beta\) and balanced splits for which the optimal
cross-edge pressure gap remains (c_\beta(m+r)).

## Initial ranking before ledger consultation

Candidate 2 would be the cleanest convergence theorem but is vulnerable to a
sampling-noise floor.  Candidate 3 has the strongest uniform passage to the
ground state, but its bridge curvature is visibly extensive.  Candidate 1 has
the only unconditional new-scale estimate here, namely (2); its decisive issue
is whether the (2/3)-power bridge target is geometrically compatible with the
set of all spin maximizers.

---

## Post-ledger correction (2026-07-30T15:36Z)

The claims above that a merely uniform `o(N)` defect in C1 or C3 would force
convergence are too optimistic.  A scalar slowly varying countermodel is

\[
b_n=n\{1+\varepsilon\sin(\log\log(n+n_0))\}.
\]

It has no normalized limit, while

\[
\sup_{m+r=N}\frac{|b_N-b_m-b_r|}{N}\longrightarrow0.
\]

Indeed, when `min(m,r)<=delta N`, its contribution is `O(delta N)+o(N)`;
when both parts are at least `delta N`, slow variation is uniform on
`[delta N,N]`.  First send `N` to infinity and then `delta` to zero.  Thus the
meaningful versions of C1/C3 need, for example, an `O(N^(1-delta))` defect (or
a stated Hammersley/dyadic-summable modulus), not an unnamed little-oh.  The
original timestamped candidates are intentionally left visible rather than
silently repaired.
