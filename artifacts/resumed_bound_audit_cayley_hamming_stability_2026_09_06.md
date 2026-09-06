# Ramsey witnesses give Hamming stability around flat Cayley signings

Date: 2026-09-06. Status: the director and convergence agent independently
reconstructed the mechanism and constants. This is a theorem about
actual signings, stronger than saturation of the Cayley family itself.
It does not prove convergence of the unrestricted minima.

## 1. Theorem and precise hypotheses

Let q tend to infinity through odd prime powers. Let A_q be real
symmetric additive-convolution matrices on F_q, and write lambda_q(a)
for their real, even additive Fourier eigenvalues. Assume that for
every fixed eta in (0,1),

    q^(-1) #{a!=0: |lambda_q(a)|<(1-eta)sqrt(q)} ->0.   (1)

No uniform operator bound on A_q is needed for the theorem below.
Let B_q be arbitrary real symmetric HOLLOW matrices satisfying

    ||B_q-A_q||_F=o(q).                                (2)

For Q(B)=max_(x in {+1,-1}^q) |x^T B x|/2, one has

    liminf Q(B_q)/q^(3/2)>=1/2.                         (3)

In particular, if A_q is an additive-Cayley hollow signing with
||A_q||op<=(1+o(1))sqrt(q), and B_q is ANY signing differing from it
on o(q^2) unordered edges, then (3) holds. There is no operator-norm
assumption on B_q.

The spectral assumption for this signing corollary follows from
q^(-1)Tr[(A_q/sqrt(q))^2]=1-1/q and the stated operator bound. The
Hamming assumption gives ||B_q-A_q||_F^2=8 d_q=o(q^2).

## 2. A bounded host for every fixed accuracy

Fix epsilon in (0,1) and eta in (0,1), BEFORE sending q to infinity.
The independently audited Ramsey/squarefree theorem supplies a finite
host H_q in F_q^*, with |H_q|<=L_epsilon uniformly in q, as follows.
For sufficiently large characteristic it is the integer Deuber host
and its negative; in each of the finitely many remaining characteristics
it is the nonzero part of a fixed-dimensional vector subspace.

For every coloring of H_q by two signs, this host contains frequencies
producing a real polynomial P and a Boolean f such that

    ||f-P||_2<epsilon, E P=0,
    Fourier support(P) is contained in H_q,

and every frequency supporting P has one common color s. Norms here
and below use uniform probability on F_q. The frequency configuration
and the chosen profile may depend on the coloring. No fixed profile
under every subsequent dilation is asserted.

Let Bad_q be the set in (1). A nonzero field scalar b is GOOD if
b H_q avoids Bad_q. The good proportion g_q satisfies

    g_q>=1-|H_q| |Bad_q|/(q-1) ->1.                    (4)

This is a union bound over the fixed host, since b h is uniform on
F_q^* for every h!=0. For each good b, color b H_q by the signs of
lambda_q and choose one resulting P_b,f_b,s_b. These obey

    ||f_b-P_b||_2<epsilon, ||P_b||_2<=1+epsilon,
    support(P_b) subset b H_q,
    s_b P_b^T A_q P_b
      >=(1-eta)(1-epsilon)^2 q^(3/2).                  (5)

The last statement uses only the lower spectral condition on the
selected support. Large exceptional eigenvalues outside it do not
enter the argument.

## 3. Averaging dilations and translations flattens the covariance

Choose b uniformly among good scalars and t uniformly in F_q. Set
P_(b,t)(x)=P_b(x+t), f_(b,t)(x)=f_b(x+t), and define the real symmetric
oriented covariance matrix

    M_xy=E_(b,t)[s_b P_b(x+t)P_b(y+t)].                (6)

Use probability-normalized Fourier coefficients
Phat_b(a)=q^(-1)sum_x P_b(x) conjugate(chi_a(x)).
Translation averaging makes M a convolution matrix. Its eigenvalues
are q mu_a, where

    mu_a=E_b[s_b |Phat_b(a)|^2].

Let nu_a=E_b|Phat_b(a)|^2. Then |mu_a|<=nu_a and

    sum_a nu_a=E_b||P_b||_2^2<=(1+epsilon)^2,
    nu_a<=(1+epsilon)^2 |H_q|/[(q-1)g_q].             (7)

For the second inequality, a fixed nonzero a belongs to bH_q for at
most |H_q| scalars b. Each coefficient square is at most the entire
squared norm. Also nu_0=0. This counting argument remains valid even
though the selected Ramsey configuration depends arbitrarily on b.

By Parseval for the convolution matrix,

    ||M||_F^2=q^2 sum_a mu_a^2
      <=q^2(1+epsilon)^4 |H_q|/[(q-1)g_q].             (8)

Consequently ||M||_F=O_epsilon(sqrt(q)). Averaging (5) and using (8)
shows, with D=B_q-A_q,

    E[s_b P_(b,t)^T B_q P_(b,t)]
      >=(1-eta)(1-epsilon)^2 q^(3/2)-||D||_F ||M||_F. (9)

Under (2), the final term is o_epsilon(q^(3/2)).

## 4. Boolean transfer for an arbitrary hollow target matrix

This is where an unproved operator bound on B_q must NOT be inserted.
Translation averaging gives the uniform coordinatewise estimates

    E[(f_(b,t)(x)-P_(b,t)(x))^2]<=epsilon^2,
    E f_(b,t)(x)^2=1,
    E P_(b,t)(x)^2<=(1+epsilon)^2.                    (10)

Let beta(B)=max_(u,v Boolean)|u^T B v|. Real Grothendieck applied to
the Hilbert-space vectors indexed by (b,t), absorbing s_b into the
left vector, gives

    |E[s_b(f^T B f-P^T B P)]|
       <=K_G beta(B) epsilon(2+epsilon).              (11)

Indeed expand the difference as (f-P)^T B f+P^T B(f-P), and use (10)
for the two vector families in each term. This is the ordinary real
bilinear Grothendieck inequality, not a matrix-operator estimate.
The elementary tensor-power/Gaussian-rounding proof with
K_G<=pi/[2 asinh(1)] is already given in Section 4 of
`fresh_range_and_spectral_regularization_2026_09_05.md`.

For completeness, hollow symmetry gives beta(B)<=4Q(B). Put
u=a+b, v=a-b with a=(u+v)/2 and b=(u-v)/2 in [-1,1]^q. Then

    u^T B v=a^T B a-b^T B b.

Since the hollow quadratic form is multilinear, averaging independent
Boolean coordinates of means a or b shows |a^T B a|,|b^T B b|<=2Q(B).
This proves the claimed bilinear bound. Hollow is an explicit and
necessary hypothesis for this particular comparison.

Every oriented Boolean quadratic energy is at most 2Q(B), so (11)
implies

    E[s_b P^T B P]
       <=[2+4K_G epsilon(2+epsilon)]Q(B).             (12)

Combining (9) and (12), then q->infinity, yields

    liminf Q(B_q)/q^(3/2)
      >=(1-eta)(1-epsilon)^2/[2+4K_G epsilon(2+epsilon)].

Now send eta and epsilon to zero. This proves (3). No rate in q
uniform as epsilon tends to zero is used or needed.

## 5. Finite Hamming form and minimizer consequence

If A,B are signings differing on d unordered edges, the proof gives

    [2+4K_G epsilon(2+epsilon)] Q(B)/q^(3/2)
      >=(1-eta)(1-epsilon)^2
        -(1+epsilon)^2 sqrt(8 d |H_q|/[q(q-1)g_q]).    (13)

All constants in the host depend only on the fixed accuracy, uniformly
in characteristic after the finite threshold split. Thus a fixed
deficit Q(B)/q^(3/2)<=1/2-delta forces a positive Hamming fraction
d/q^2, after choosing epsilon,eta sufficiently small in terms of delta.
The fraction may be extremely small because the Ramsey host is large,
but it is a genuine positive constant independent of q.

Therefore any actual minimizing family with a persistent sub-1/2
coefficient must remain a constant Hamming fraction away from every
spectrally near-Frobenius-optimal finite-field additive-Cayley family.
The same assertion holds after arbitrary vertex switchings, simultaneous
row/column permutations, and a global matrix sign, because all these
operations preserve Q, Frobenius distance, and Hamming distance.
This is a structural restriction on possible minimizers and nearby
upper constructions. It does not show that all minimizers are close
to Cayley signings, provide a fixed-block insertion estimate, or settle
convergence of M_n/n^(3/2).
