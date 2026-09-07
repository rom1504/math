# Exact finite-temperature weak duality, and the annealed saddle obstruction

2026-09-07. The finite-temperature inequality below is a proved consequence
of the audited marked lower policy and actual weave moment bound. It is not
a min-max equality. The strict Haar obstruction is independently audited in
`principle_director_haar_stability_audit_2026_09_07.md`.

## 1. Definitions and a pointwise Gibbs inequality

For a symmetric hollow matrix B of order N put

    H_B(x)=x^T B x/2,
    Z_beta^sigma(B)=sum_(x in {+-1}^N) exp(beta sigma H_B(x)),
    P_beta(B)=[log Z_beta^+(B)+log Z_beta^-(B)]/(2N).

Let F_i,H_i be any real, possibly randomly generated, fields satisfying
H_i>=0 and |F_i|+H_i<=1. Define

    d_i=H_i sign((BF)_i),       mu_i^+=F_i+d_i, mu_i^-=F_i-d_i.

Both means belong to [-1,1]. The Gibbs variational inequality applied to
the two product measures with these means gives, pointwise in all fields,

    P_beta(B) >= beta/N sum_i H_i |(BF)_i|
        + 1/(2N) sum_i [h2((1+F_i+H_i)/2)
                       +h2((1+F_i-H_i)/2)].                 (1)

Here h2(q)=-q log q-(1-q)log(1-q), with 0 log 0=0. To check the energy
factor, H_B(F+d)-H_B(F-d)=2 F^T B d. Division by 2N therefore yields
exactly the first term in (1), not twice that term. The entropy pair is
unchanged when sign((BF)_i) exchanges the two means. Hollowness ensures
the product-measure expected Hamiltonian is H_B(mu), with no diagonal
correction. No independence between F and B is needed.

## 2. Passing through the successful marked policy

Use B=A/sqrt(N), where A is a hollow full signing in the fixed-operator
regime. In the old marked Gaussian frame let F be odd, H be even,
H>=0 and |F|+H<=1. Let P1 be first-chaos projection and U the audited
isometry from even L2 to first chaos. Define

    K_F=U^(-1)P1F,
    tau^2=||F-P1F||_2^2,
    J(F,H)=E[H E_Z |K_F+tau Z|],
    s(F,H)=(1/2)E[h2((1+F+H)/2)+h2((1+F-H)/2)],             (2)

where Z is an independent standard Gaussian. The marked return theorem
and (1) imply

    liminf_N P_beta(A/sqrt(N)) >= beta J(F,H)+s(F,H).         (3)

The finite polynomial/finite frame approximations are taken first, then
N tends to infinity, then the bounded gate approximations are removed.
For the entropy, boundedness and uniform continuity of h2 on [0,1]
suffice; no derivative bound at 0 or 1 is required. This statement only
uses the fixed-operator form of the marked theorem. That is enough for
the actual weave below. No additional pressure-preserving spectral
regularization is being claimed here.

The source marked inequality and its ordered approximation are in
`decisive_audit_fresh_full_lower_chain_2026_09_07.md`, Sections 3--5,
and `decisive_audit_certified_minimum_width_lower_2026_09_07.md`.

## 3. Exact conversion of the weave parameters

The actual weave has m fibres with k retained physical rows each,
N=mk, p=k/m. Its full signed matrix W has norm m, because the conjugated
signed coordinate swap is an involution. Hollowing changes its operator
norm by at most one. Thus ||W_hol/sqrt(N)||op<=1/sqrt(p)+o(1).

The audited defect identity is

    D_sigma=2(m^2 k-sigma x^T W x).

Hence the defect tilt exp[-tD_sigma/(2k)] equals

    exp[-tm^2+(t/k)sigma x^T W x].

The temperature corresponding to the normalized half-Hamiltonian is

    beta=2t sqrt(N)/k=2t/sqrt(p).                           (4)

For every fixed depth r, the actual all-spin row bound gives

    limsup_N (1/N) log E Z_beta^sigma(W_hol/sqrt(N))
       <= [t+p log 2+B^r Phi_t(nu_p)]/p,

where nu_p=(1-p)delta_0+(p/2)(delta_(1/sqrt(p))+
delta_(-1/sqrt(p))). Sending finite depth r to infinity gives

    A(p,t)=[t+p log 2+E_t(nu_p)]/p.                        (5)

Hollowing costs o(N) in log partition, since its unnormalized Boolean
energy shift is O(N), while beta is fixed and the normalization is
sqrt(N). These are moment UPPER bounds for an actual ensemble; E_t is
not asserted to equal its limiting quenched or annealed pressure.

Jensen's inequality, (3), and (5) now give the correctly normalized
weak duality

    (2t/sqrt(p)) J(F,H)+s(F,H)
        <= [t+p log 2+E_t(nu_p)]/p                         (6)

for every feasible marked policy and every p in (0,1], t>0 to which the
weave closure applies. In particular the zero-entropy Boolean gate
policies recover the familiar cap inequality after division by beta.
Equation (6) is a genuine theorem connecting the two successful objects,
but its proof includes Jensen and the upper orbit/graph relaxations.
Neither operation carries an equality conclusion.

## 4. Why the shared Gaussian number does not prove an exact saddle

At the Gaussian boundary the upper scalar function satisfies

    a(t)=t+g_t(1)
        =sup_(-1<r<1) [tr+(1/4)log(1-r^2)].                (7)

This is exactly the limiting log moment-generating function of the
Rayleigh variable r=x^T U x/N for a Haar balanced involution U. Counting
all Boolean vectors gives

    inf_(t>0) [a(t)+log 2]/(2t)=sqrt(15)/8.

The permanent/Bellman lower floor at this same number is already proved
in `transfer_director_exact_permanent_floor_2026_09_06.md`. It is a floor
for that upper certificate, not for the actual signing cap.

The new local-stability argument proves for the actual Haar matrix

    limsup_N E max_x |x^T U x|/(2N)
       <= (1/2)sqrt(1-exp(1/25000)/16) < sqrt(15)/8.         (8)

Its proof and exact constants are in
`principle_invent_2026_09_07_haar_stability_strict_gap.md`.
Wherever the archived fixed-GFOM comparison to Haar applies, (8)
strictly lowers its ceiling. Thus the present Bellman/permanent upper
certificate and the present fixed-GFOM lower control class cannot be
matching exact primal/dual classes. The equality of their OLD boundary
numbers came from the same all-vector annealed tail, not from an exact
ground-state duality. This does not falsify the possibility of a different
min-max theory for original signs, and (8) is not a sign-matrix upper.

## 5. A separate concrete local-state incompleteness test

The universal one-root marked Gaussian frame is not by itself an exact
per-matrix ground-state descriptor. A symmetric Hadamard H, hollowed to
A=H-diag(H), has Q(A)<=(1/2+o(1))N^(3/2). In contrast, iid symmetric
off-diagonal signs have, with probability tending to one,

    Q(A)>=[(2/3)sqrt(2/pi)-o(1)]N^(3/2),

whose constant is about .531923. For the latter statement choose spins
sequentially, setting x_i=sign(sum_(j<i) A_ij x_j). Conditional on the
revealed past the increment has the law |epsilon_1+...+epsilon_(i-1)|,
independent of that past. The sum therefore has mean
[(2/3)sqrt(2/pi)+o(1)]N^(3/2) and variance O(N^2), giving the assertion
by Chebyshev. An elementary sphere-net bound gives ||A||op=O(sqrt(N))
with high probability, so this lies within the same fixed-operator
marked-frame regime.

This comparison disproves a per-model exact-value reading of only the
local Gaussian state. It does NOT disprove an exact saddle after taking
the infimum over A. The missing information in this example is concrete:
ordered exposure includes the conditional law of the unexposed matrix.
The iid variance clock survives that exposure; Hadamard orthogonality
does not permit the same conditional independence argument.

The next actual upper operation is correspondingly to retain local-field
stability under the weave's edge tilt, rather than identify the entropy
envelope with the pressure. Its exact finite statement is recorded in
`principle_invent_2026_09_07_stability_aware_weave.md`.
