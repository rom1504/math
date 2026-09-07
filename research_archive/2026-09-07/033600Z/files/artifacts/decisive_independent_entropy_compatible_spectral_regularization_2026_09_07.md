# Entropy-compatible spectral regularization of the continuous outer model

Status: independently audited PASS by transfer_seeds. This removes an
operator-norm hypothesis from the continuous variational program; it
does NOT prove the required bounded-operator thermodynamic limit.

The deterministic core/refill module is already archived in
`fresh_range_and_spectral_regularization_2026_09_05.md`. The contribution
here is a disorder-law relative-entropy comparison for an ADAPTIVELY
selected core. This is unrelated to the spin-shell entropy no-go in
`joint_entropy_spectral_regularization.md`.

## 1. Definitions and main inequality

Let d=n(n-1)/2, n>=2, let U_d be uniform probability on the cube of
hollow symmetric matrices with entries in [-1,1], and define

    q_n(B)=Q(B)/n^(3/2),
    delta(B)=d^(-1) sum_(i<j)(1-B_ij^2),
    F_tau(B)=q_n(B)+tau delta(B),
    L_n(alpha,tau)=-(alpha d)^(-1) log integral exp(-alpha d F_tau) dU_d.

Let L_n^(L)(alpha,tau) be the same expression with the integral
restricted to ||B||op<=L sqrt(n). The restricted reference measure is
NOT renormalized. Fix alpha,tau>0 and 0<epsilon<=1/2. Set

    C=sqrt(log2/3)+(2/3)tau+1,
    K=pi/[2 asinh(1)],
    L=4KC/epsilon+8.

Then

    0 <= L_n^(L)-L_n
      <= 2sqrt(epsilon)+4tau epsilon
          +(n+1)log2/(alpha d)
          -log(1-exp(-alpha d))/(alpha d).                (1)

Thus for every fixed alpha,tau the unrestricted outer free energy is
uniformly approximated, as n tends to infinity and then epsilon tends
to zero, by a fixed-normalized-operator model. In particular, proving
an all-order limit for each fixed-L restricted model would give the
unrestricted continuous-model limit, and then the soft-flatness bridge
would give the original limit. The fixed-L limit is still unproved.

## 2. Deterministic core, also valid for real cube coefficients

For any real hollow symmetric B, polarization gives beta(B)<=4Q(B),
where beta is the Boolean infinity-to-one norm. The elementary real
Grothendieck inequality and the diagonal SDP dual give D diagonal with

    D >= B, D >= -B, D_ii>=0, Tr D<=K beta(B)<=4KQ(B).

These statements do not use flat entries. Delete the coordinates with
D_ii>4KC sqrt(n)/epsilon. If q_n(B)<=C, fewer than epsilon n vertices
are deleted, and the retained principal block satisfies

    ||B_R||op<=4KC sqrt(n)/epsilon.

For measurability one need not choose an SDP solution: enumerate the
finitely many deletion sets of size less than epsilon n and select the
first whose retained block meets this operator bound. Existence follows
from the preceding argument, and the selection is measurable because
operator norm is continuous. Write S(B) for the selected deleted set.

Principal monotonicity Q(B_R)<=Q(B) follows by completing the missing
spins independently and averaging. It too is valid for real hollow B.

## 3. Uniform good-refill probability

For a fixed deletion set S, independently refill all h edges incident
to S by uniform[-1,1] variables, and put zero on all other edges of
the noise matrix E. Here h<=epsilon n^2. Hoeffding and the projectively
reduced two-sided spin family give

    Pr{Q(E)>2sqrt(epsilon)n^(3/2)}
        <=exp[-(2-log2)n].                               (2)

For a unit vector x, the quadratic form x^T E x has independent
summands 2E_ij x_i x_j and sum of squared coefficient magnitudes at
most 2. Thus

    Pr{|x^T E x|>t}<=2exp(-t^2/4).

A 1/4-net of the unit sphere with size at most 9^n and the standard
||E||op<=2 max_net |x^T E x| estimate give

    Pr{||E||op>8sqrt(n)}<=2exp[-(4-log9)n].                (3)

For n>=2 the sum of (2) and (3) is below 1/2. (If S is empty, E=0.)
Condition the refill law on both good events. Its relative entropy
against unconditioned uniform refill is at most log2, uniformly in S.

The output B'=B_R padded by zero, plus E, satisfies deterministically

    ||B'||op<=L sqrt(n),
    q_n(B')<=q_n(B)+2sqrt(epsilon),
    delta(B')<=delta(B)+h/d<=delta(B)+4epsilon.

Therefore F_tau(B')<=F_tau(B)+2sqrt(epsilon)+4tau epsilon.

## 4. Adaptive-mask relative-entropy comparison

Let P be ANY law supported on {q_n(B)<=C}, absolutely continuous with
respect to U_d. Form the joint law of B, its deterministic mask S(B),
and the good-conditioned refill noise. Compare it with the reference
law in which:

- B has law U_d;
- S is independently uniform over all 2^n vertex subsets;
- conditional on S, the refill coordinates are independent uniform.

The relative entropy of the actual joint law against this reference
is at most

    D(P||U_d)+n log2+log2.                               (4)

The n log2 term is exactly the price of recording an arbitrary adaptive
mask; the refill conditioning contributes at most log2. Applying the
map that keeps B on the retained block and writes the refill entries
elsewhere, the reference output is again U_d. Data processing therefore
proves for the output law P'

    D(P'||U_d)<=D(P||U_d)+(n+1)log2.                     (5)

The selected core may depend on the entire input matrix. No false
independence assumption is made about the selected core itself.

## 5. Applying the transformation to the outer Gibbs law

Jensen and the strong-variance smoothmax bound give

    L_n<=E_(U_d) F_tau
        <=sqrt(log2/3)+(2/3)tau=C-1.                     (6)

Let P_G be the full outer Gibbs law. Since F_tau>=q_n, its mass outside
{q_n<=C} is at most exp(-alpha d). Restricting and renormalizing that law
therefore raises the variational free energy by at most

    -log(1-exp(-alpha d))/(alpha d).                     (7)

Apply sections 2--4 to this conditional law. The Gibbs variational
identity, with reference U_d, gives

    L_n^(L)=inf_(P:op support<=Lsqrt(n))
        [E_P F_tau+D(P||U_d)/(alpha d)].

Equations (5), (7), and the energy bound in section 3 prove (1).

## 6. Scope for the original-limit program

The n^2 disorder-entropy scale is respected: the adaptive mask uses
only O(n) bits. The refill is same-order and its cap cost is O(sqrt(epsilon))
after normalization. This justifies bounded-operator restriction for
the ENTIRE continuous outer variational problem, not only for one
deterministically selected minimizing matrix.

It still does not supply a large-deviation principle for the action
profiles, nor an all-order recovery theorem for them. Ordinary dense
graph limits do not resolve the n^(3/2) observable; ordinary spectral
limits do not determine the same-spin Boolean cap. Those missing limit
theorems must not be inferred from (1).

In an action-profile compactification, the variance occupancy must also
be tracked explicitly, or its continuity proved separately. We make no
claim that action convergence alone determines delta(B). For the outer
integral, even all-order realization of a SINGLE matrix profile would
not suffice: one needs the corresponding exponential-volume recovery
at speed d, or another direct proof of the restricted pressures' size
limit. Compactness only supplies subsequences and does not identify
their microstate-volume rates.

## 7. Discrete microcanonical counting corollary

The same proof applies with the fair-sign product measure instead of
the uniform continuous cube, using fair-sign refill. The concentration
estimates and the adaptive-mask entropy accounting are unchanged.
Let

    N_n(C)=#{flat A: q_n(A)<=C},
    N_n(C,L)=#{flat A: q_n(A)<=C, ||A||op<=Lsqrt(n)}.

Count labeled signings, without quotienting switching or permutations.
For every C>0, 0<epsilon<=1/2, and L=4KC/epsilon+8,

    N_n(C+2sqrt(epsilon),L) >= 2^(-(n+1)) N_n(C).        (8)

If N_n(C)=0 there is nothing to prove. Otherwise take P uniform on its
low-cap set and apply the discrete regularization channel. Since KL
against the fair-sign product measure is d log2 minus Shannon entropy,
(5) gives

    H(P') >= log N_n(C)-(n+1)log2.

The output support is contained in the left-hand set in (8), and its
entropy is at most the logarithm of that set's cardinality. This proves
the count inequality.

Consequently bounded-operator regularization preserves the leading
n^2-speed low-cap disorder entropy, allowing only an arbitrarily small
fixed cap-threshold increase. This is stronger than applying a
deterministic existence theorem to a single minimizing signing. It
still does not assert that either counting entropy has an all-order
limit, or that one fixed L captures an o(1)-near-minimizing sequence.
