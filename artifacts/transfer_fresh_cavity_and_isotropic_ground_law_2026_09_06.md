# Fresh cavity attempt and a Paley extension obstruction

Date: 2026-09-06. All energies in this note count each unordered edge once:

`H_A(x)=x^T A x/2`, `Q(A)=max_x |H_A(x)|`, `M_n=min_A Q(A)`.

## 1. Independently frozen convergence target; archive reconciliation

Before reading route assessments, I isolated

`E(A)=min_{b in {+-1}^n} max_x (|H_A(x)|+|b dot x|)`.

The target was: there are absolute `C<infinity`, `eta>0` such that, for
at least one exact minimizing signing A at every sufficiently large order n,

`E(A) <= M_n + 3 M_n/(2n) + C n^(1/2-eta)`.                 (1)

It would imply convergence. Indeed, writing `c_n=M_n/n^(3/2)` and using
the boundedness of c_n, expansion of `(n+1)^(-3/2)` gives

`c_(n+1) <= c_n + O(n^(-1-eta)+n^(-2))`.

The positive error is summable, so adding its remaining tail produces an
eventually nonincreasing bounded sequence differing from c_n by o(1).
The exact proof obligation in (1) is a signed-row discrepancy theorem for
the *whole labelled slack landscape* `g_A(x)=Q(A)-|H_A(x)|`.

Archive inspection then found essentially this same target in
`continued_convergence_proportional_thinning_2026_09_06.md`, Section 1,
and `ar_action_independent_proposals.md`, Section 2. It is not new progress.
The sharp all-star optimality identity is also already in
`near_cap_insertion.md`, Proposition 2.1. In particular, every row of an
exact minimizer optimally extends its deleted core, but the latter need
not minimize at its own order; the exact compensation identity there
prevents treating that observation as a recurrence for M_n.

No finite zero-defect counterexample, nor a nonoptimal conference family,
refutes the exact-minimizer/power-saving formulation (1).

## 2. An exact covariance criterion for extension obstruction

Let mu be any probability law on Boolean spins, let

`K=E_mu xx^T`, `g(x)=Q(A)-|H_A(x)|`, `m=E_mu g`,
`v=Var_mu(g)`, and `kappa=lambda_min(K)`.

For any signed row b, put

`Delta_A(b)=max_x (|b dot x|-g(x))`.

This quantity is nonnegative, since it can be tested on an exact ground
state. Pointwise `|b dot x| <= Delta_A(b)+g(x)`. Squaring and averaging,

`n kappa <= b^T K b <= (Delta_A(b)+m)^2+v`.

Consequently

`Delta_A(b) >= max(0, sqrt(max(0,n kappa-v))-m)`.           (2)

In particular, an **isotropic exact ground-state law**, meaning

`mu{ |H_A(x)|=Q(A) }=1`, `E_mu xx^T=I`,

forces `Delta_A(b)>=sqrt(n)` for EVERY signed row b. If n is even the
dot products are even integers, giving the parity improvement

`Delta_A(b)>=2 ceil(sqrt(n)/2)`.                            (3)

This is an actual extension obstruction, not a statement about an iid
row, a separately chosen row for each state, or a scalar cap profile.

## 3. Every Paley conference has an isotropic ground-state law

The square-field construction below was the first proof found, but the
covariance conclusion in fact holds for ALL Paley conferences. No
Boolean eigenvector hypothesis is required.

More generally, suppose a hollow signing A has a group G of signed
permutation symmetries `R^T A R=A`, with the underlying permutations
transitive on unordered vertex pairs. Suppose also that a signed
permutation T satisfies `T^T A T=-A`. Then A has an isotropic law on its
exact absolute ground states.

**Proof.** The anti-symmetry makes the two oriented extrema equal. Start
with a positive ground x and average `R x` over the finite group G.
Preservation of A and pair transitivity make
`a_ij E (Rx)_i(Rx)_j` constant on all unordered edges. Its sum is Q(A),
so with `D=binom(n,2)` the positive ground covariance is

`K_+=I+(Q(A)/D) A`.

Applying T to this law gives a negative ground law with covariance
`K_-=I-(Q(A)/D) A`. Their equal mixture is I. This uses exact symmetry
of one matrix, not permutation averaging across different minimizers.

Now let q be ANY odd prime power with `q=1 mod 4`, let F=GF(q), and let
C be the Paley conference on `F union {infinity}` with finite entries
`chi(z-w)` and border entries one. Write `n=q+1`. The required symmetries
are explicit:

* translations `z -> z+a`, fixing infinity, preserve C;
* inversion interchanges zero and infinity and sends nonzero z to 1/z.
  Give the inversion signed permutation column signs `d_0=d_infinity=1`
  and `d_z=chi(z)` elsewhere. The identity
  `chi(1/z-1/w)=chi(z-w)chi(z)chi(w)` proves preservation of C;
* multiplication by a nonsquare, together with a sign switch at
  infinity, sends C to -C.

Translations and signed inversion are pair-transitive: translate one
of two finite points to zero, invert it to infinity, then translate
the other to zero. Pairs already containing infinity need only a
translation. Thus the preceding symmetry argument applies.

Consequently EVERY Paley conference satisfies the exact theorem

`E(C) >= Q(C)+2 ceil(sqrt(n)/2)`.                         (P)

Since `C^2=(n-1)I` gives `Q(C)<=n sqrt(n-1)/2`,

`E(C)-(1+1/n)^(3/2) Q(C) >= (1/4-o(1)) sqrt(n)`.          (P')

The error in (P') is uniform over these prime powers; it follows by a
one-term Taylor expansion. The cap Q(C) need not be known, need not
saturate its spectral bound, and can be below the 1/2 scale at a given
order. Thus this theorem is stronger than an eigenbasis-only argument.
It is still NOT an asymptotic statement about exact minimizing parents:
only some small Paley members have certified global minimality.

## 4. An independent explicit law for square-field Paley conferences

Let r be an odd prime power, F=GF(r^2), K=GF(r), and let chi be the
quadratic character of F, with chi(0)=0. The symmetric Paley conference
matrix C, indexed by F union {infinity}, is

`C_(z,w)=chi(z-w)` for distinct finite z,w,
`C_(infinity,z)=C_(z,infinity)=1`, with zero diagonal.

Write `n=r^2+1`. The standard character identity gives `C^2=r^2 I`.
For completeness, finite off-diagonal entries of C^2 vanish because
`sum_z chi(a-z)chi(z-b)=-1` for a!=b, cancelling the border contribution.
The diagonal entries are r^2, and border off-diagonal entries vanish by
`sum_z chi(z)=0`.

Fix `t in F minus K`. For each square `u in F*`, write `uz=a+bt`, with
a,b in K. Choose f:K->{+-1} uniformly subject to `sum_b f(b)=1`, and set

`x_infinity=1`, `x_z=f(b)` whenever `uz=a+bt`.

These are Boolean positive eigengrounds: `Cx=r x`. Indeed, within one
affine K-fibre the off-diagonal character sum is r-1. Between two
distinct fibres the sum is -1. The latter follows by writing the
character as the K-character of the norm and using
`sum_a chi_K((a+c)^2-d)=-1` for nonzero d. Therefore, at a finite vertex
in fibre b,

`(Cx)_z=1+(r-1)f(b)-sum_(d!=b)f(d)=r f(b)`.

At infinity the sum is `r sum_b f(b)=r`.

Let mu_+ be the law from uniform square u and the independent uniform
balanced f above. Its covariance is EXACTLY

`E_(mu_+) xx^T=I+C/r`.                                    (4)

The diagonal is one, and border covariances equal `E f(b)=1/r`. For
different fibres, `E f(b)f(d)=-1/r`. A nonsquare finite difference is
never mapped into K* by a square u, so its covariance is -1/r. A square
finite difference is mapped into K* with probability `2/(r+1)`, giving

`2/(r+1) + (r-1)/(r+1) (-1/r)=1/r`.

Choose a nonsquare v in F*. Multiplication of finite labels by v changes
every finite character sign. Combining that permutation with a sign
switch at infinity gives a signed permutation R satisfying
`R C R^T=-C`. Pushing mu_+ forward by R therefore gives a Boolean
negative eigenstate law mu_- with covariance `I-C/r`. The equal mixture

`mu=(mu_++mu_-)/2`

is isotropic and supported on exact absolute ground states. Spectral
boundedness and the displayed eigenstates give

`Q(C)=q=rn/2`.

The positive covariance (4) is the same explicit law used for the
edge-balanced trap in `cross_order_exposed_shell_sparse_repair_no_go.md`,
Section 3. The two-orientation isotropic mixture and its consequence for
EVERY appended row are the points retained here. A targeted archive
search for insertion/eigenbasis/isotropic-ground combinations did not
locate that consequence; this is a novelty check, not a claim about all
mathematical literature.

## 5. Explicit square-field consequence

Combining (3) with n=r^2+1 and odd r gives

`E(C) >= rn/2+r+1`.                                      (5)

Thus every one-vertex extension of these parents has excess at least

`(r+1)-3q/(2n)=r/4+1`

above the leading derivative `3q/(2n)=3r/4`. More exactly,

`E(C)-(1+1/n)^(3/2) Q(C) >= r/4+1-O(1/r)`.

Along the unbounded prime-power family this rules out a universal
conference-parent version of (1), even allowing `o(sqrt n)` error and
even optimizing the new row using the full parent. It does NOT refute
(1) for selected exact minimizers. These parent caps have asymptotic
constant 1/2, and the project's strict all-order upper bound below 1/2
in fact makes the large members nonoptimal.

## 6. Why the criterion cannot silently become a block theorem

Isotropy alone gives, for an arbitrary n-by-h SIGNED bridge B,

`max_(ground x) ||B^T x||_1 >= sqrt(nh)`,

not `h sqrt(n)`. The same spin must serve all bridge columns. The
Paley positive ground law above is not uniformly subgaussian: if b is
one of the fibre-constant positive eigengrounds, its correlations can
concentrate in the rare matching fibre direction. Therefore replacing
the second-moment argument by an unjustified average absolute-value
bound would create a false mesoscopic conclusion.

## 7. Status of the original target

The only convergence-target claim remains unproved (1). Existing exact
star replacement identities do not close it. The concrete results (P)
and (5) exclude broader Paley-parent statements; they do not establish
an asymptotic exact-minimizer counterexample. No new bound
on liminf or limsup c_n, and no convergence/nonconvergence conclusion,
is asserted.

## 8. Exact computational replay and primary-source scope

`python3 computations/transfer_fresh_paley_isotropic_extension_verify.py`
checks the conference identity, both signed eigenstate identities, and
both covariance identities by integer arithmetic for r=3,5,7. The
positive laws have respectively 12,120,840 equally weighted atoms,
counted with their construction multiplicities.

For r=3, n=10, it also exhausts all 512 old spins and all 512 appended
rows modulo their independent global signs. It finds `Q(C)=15`,
`E(C)=19`, and 20 optimal appended rows modulo sign. Thus (5) is sharp
at this finite member. At r=5,7 the computation verifies the laws and
the resulting lower bounds 71 and 183, not exact extension optima.

Primary-source queries located the original construction in
[Paley, *On Orthogonal Matrices* (1933)](https://onlinelibrary.wiley.com/doi/10.1002/sapm1933121311),
the conference spectral setting in
[Haemers and Parsaei Majd, *Spectral symmetry in conference matrices*](https://arxiv.org/abs/2004.05829),
and the related square-field eigenfunction problem in
[*On eigenfunctions and maximal cliques of generalised Paley graphs of square order*](https://arxiv.org/abs/2203.16081).
The search results do not establish the best signed-row extension
inequality (5), nor a theorem for exact minimizing signings. No such
import is used here: the conference, fibre, covariance and extension
calculations required for (5) are all proved explicitly above.
