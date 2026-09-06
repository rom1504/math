# Independently programmable projective Fourier phases

Date: 2026-09-06. Status: exact auxiliary quotient theorem, transferred
to the prescribed R norm by the preceding resonant-carrier theorem.
This gives a large explicit family of nonlocal tests and an exact
remaining Boolean phase-feasibility condition. It does not prove that
condition can be met or that R=T.

## 1. An orthogonal projective character matrix

Fix r>=1, let G=F_7^r, q=7^r, and t=(q-1)/6. Choose one nonzero
representative from every one-dimensional F_7 subspace, independently
for rows a and columns l. Put

    C_(a,l)=chi_7(a dot l),       chi_7(0)=0.

Then the t-by-t matrix C obeys the exact identity

    C C^T=7^(r-1) I_t.                                     (1)

Indeed, for distinct row lines the two linear forms a dot x and
b dot x are independent, so the sum of their quadratic-character
product over x in G is zero. For the same representative the sum is
6*7^(r-1). Each nonzero projective line contributes six identical
products because chi_7(c)^2=1 for c!=0. Dividing by six proves (1).

Consequently every real function psi on G satisfying

    psi(0)=0,       psi(ca)=chi_7(c) psi(a) for c!=0        (2)

has an exact expansion

    psi(a)=sum_l theta_l chi_7(a dot l),
    theta=7^(1-r) C^T psi_representatives.                 (3)

This is all of the space in (2), not just a particular symmetric
one-parameter subfamily.

## 2. One tensor quotient realizes the entire phase family

The resonant-carrier theorem makes every full seven-point circulant
rotation V_theta available as a limiting spin-profile test for R.
Its Fourier multipliers are

    lambda(0)=1,       lambda(j)=exp(i theta chi_7(j)).

Tensor the t independently chosen operators V_(theta_l) from (3).
Let L:F_7^t -> F_7^r have these projective representatives as columns.
It is onto, so its fibers have equal sizes, and pullback along L is
a unital Boolean-preserving isometry of probability spaces.

The character indexed by a pulls back to the character indexed by
(a dot l)_l. The tensor multiplier on it is therefore

    product_l exp(i theta_l chi_7(a dot l))=exp(i psi(a)).

Thus the full algebra of functions of L is exactly invariant under
this auxiliary tensor operator. Its quotient is the real orthogonal
circulant operator U_psi on G with multipliers

    lambda(0)=1,       lambda(a)=exp(i psi(a)).             (4)

For any fixed finite seed B this proves the lower theorem

    R(B)>=(1/(2q)) max_(F,G Boolean q-by-k)
                                   |Tr[G^T U_psi F B]|.   (5)

The physical realizations of the seven-point factors are signed spin
modules, but tensoring those modules and then choosing the Boolean
profiles pulled back along L is legitimate. No multiplication of
overlapping quotient gates is used.

The exact phase restriction must be retained. Multipliers are equal
on {a,2a,4a}, and conjugate on {-a,-2a,-4a}. The theorem does not give
independent phases to all nonzero frequencies of G.

## 3. Exact elimination of all phase variables

Use the probability-normalized Fourier transform

    fhat(a)=(1/q) sum_x f(x) exp(-2 pi i a dot x/7).

For real profiles F,G put

    s(a)=sum_(i,j) conj(Ghat_i(a)) B_ij Fhat_j(a),
    S_l=sum_(c in {1,2,4}) s(ca_l).

The complementary three frequencies contribute the complex conjugate.
The constant term is real. Since the projective phases are freely
and independently programmable, their exact optimized objective is

    max_psi |(1/q)Tr[G^T U_psi F B]|
       =|mean(G)^T B mean(F)|+2 sum_l |S_l|.              (6)

Therefore (5) strengthens to the completely explicit finite test

    2R(B)>=max_(F,G Boolean)
                 [|mean(G)^T B mean(F)|+2 sum_l |S_l|].   (7)

In particular there is no continuous phase optimization left once
the Boolean profiles are selected. If G=F and B is real symmetric,
each s(a) is real, and phases zero or pi suffice in (6).

## 4. A precise remaining polar condition for the weighted seed

For B=[-1,2;2,4], let

    P=sqrt(2) diag(1,4),
    M=P^(-1/2) B P^(-1/2)=[-1,1;1,1]/sqrt(2).

Then M is a real symmetric orthogonal involution, Tr(P)=5 sqrt(2),
and T(B)=Tr(P)/2. For Boolean F,G and any orthogonal U,

    Tr(P)-(1/q)Tr[G^T U F B]
      =(1/(2q)) ||U F P^(1/2)-G P^(1/2) M||_F^2.        (8)

Thus reaching T(B) in the programmable phase family is equivalent
to making the displayed normalized squared defect tend to zero.
The phase variables in that defect can also be eliminated exactly.

Write A(a)=Fhat(a)P^(1/2) and D(a)=Ghat(a)P^(1/2)M, as complex
row vectors. On one QR triple l, collect the three rows into
matrices A_l,D_l. The minimum defect in (8) is

    (1/2)||A(0)-D(0)||^2
       +sum_l [||A_l||_F^2+||D_l||_F^2
                                -2|<D_l,A_l>_F|].       (9)

Every summand is nonnegative. Consequently this family reaches the
majorant value if and only if there are larger and larger Boolean
profiles for which the means match as in the first term and the
total phase-collinearity defects of the QR triples tend to zero.
Each triple's two 3-by-2 matrices must become proportional by a
single unit complex scalar in the aggregate L2 sense.

This is a concrete Boolean Fourier feasibility problem, not a claim
that arbitrary orthogonal polar actions have already been obtained.
It also exposes why independently programmable phases alone are
not an automatic proof of R=T: a scalar phase cannot change the
relative feature directions within a Fourier triple.

## 5. Full Paley rotations of arbitrarily large odd field degree

For odd r, identify G additively with F_q. Its quadratic character
satisfies chi_q(ca)=chi_7(c)chi_q(a) for c in F_7^*. Thus

    psi(a)=theta chi_q(a)

satisfies (2), and is programmable by (3). The resulting U_psi is
the full real Paley rotation

    E_q+cos(theta)(I-E_q)+sin(theta) A_q,

where A_q, up to an overall sign convention, has entries
chi_q(x-y)/sqrt(q), is skew-symmetric, and squares to -(I-E_q).
This follows directly from the elementary quadratic Gauss transform:
the transform of chi_q has modulus sqrt(q), and is purely imaginary
because q=3 mod4. Alternatively (3) proves the phase realization
without evaluating its spatial kernel.

There is also a uniform-angle signed-simplex description. Choose
one l per F_7 projective line with chi_q(l)=+1. Such representatives
exist since r is odd. The trace-fiber quadratic-character identity
gives

    sum_l chi_7(Tr(a l))
       =sigma 7^((r-1)/2) chi_q(a),       sigma in {+1,-1}.

For completeness, if H(t)=sum_(Tr x=t) chi_q(x), scalar changes show
H(0)=0 and H(t)=chi_7(t)H(1). The Gauss transforms give
|H(1)|=sqrt(q/7), and H(1) is a real integer. Summing over the six
nonzero t and dividing the six representatives per projective line
proves the formula. Tensor angle sigma theta/7^((r-1)/2) on every
line therefore gives the same full Paley action.

## 6. Computation is only exploratory

`tmp/resumed_convergence_programmable_phases.py` implements (6) via
finite-dimensional FFTs and alternating Boolean/phase best responses.
Its first bounded q=49 weighted-seed run found no improvement over
the elementary 7/2 baseline. This is neither an upper bound nor a
positive phase-defect certificate. The new proved results are (1),
(3), (5), (7), and the exact remaining condition (9).
