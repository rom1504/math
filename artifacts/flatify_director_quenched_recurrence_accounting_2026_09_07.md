# Exact recurrence payment supplied by quenched sign universality

2026-09-07. Proved conditional implication; its Gaussian hypothesis is OPEN.
This records what the new tool pays, not another claim of flatification.

Let N=m+l be a comparable split, and fix selectable actual optimizing
children A_m,D_l. Write their actual unscaled sum as h(x,y). Let

    T = sqrt((N-1)/(m-1)) M_m + sqrt((N-1)/(l-1)) M_l.

Suppose one constructs an ml-dimensional unit-diagonal Gaussian covariance
R with epsilon I<=R<=KI, where epsilon>0,K are uniform in all orders.
Let S=(2/pi)arcsin[R] entrywise, and Y be the rectangular Gaussian bridge
with covariance S. The exact hypothesis needed for this particular tool is

    E max_(x,y) |h(x,y)+x^T Y y| <= T + C N^(3/2-delta).       (G)

Then a hollow full signing parent exists, preserving both children, with

    Q(parent) <= T + C N^(3/2-delta)
                         + C_(epsilon,K) N^(4/3)sqrt(log N).  (1)

Proof: apply the quenched theorem to d=ml, features sigma x_i y_j,
offset beta sigma h(x,y)/sqrt(N), and lambda=beta/sqrt(N).
The pressure error is O(beta^3 sqrt(N) log(N)^(3/2)). Dividing by
lambda converts to cap error O(beta^2 N log(N)^(3/2)); replacing the
soft maximum costs O(N^(3/2)/beta). Set beta=N^(1/6)/sqrt(log N).
These calculations are uniform over comparable splits. A sign realization
with cap no larger than its expectation exists, proving (1).

After division by sqrt(N-1), the u_n=M_n/sqrt(n-1) recurrence defect is

    O(N^(1-delta) + N^(5/6)sqrt(log N)).                       (2)

This has summable geometric normalized error. For a pure power statement,
every fixed delta'<min(delta,1/6) gives O(N^(1-delta')). The verified
balanced-tree recurrence therefore proves convergence if (G) holds at
all comparable sufficiently large orders. An equal-size-only comparison
is not silently enough.

For a moving covariance gap epsilon_N=N^(-a), the same argument with
beta=N^((1/2-2a)/3)/sqrt(log N) yields normalized error

    O(N^(-(1/2-2a)/3)sqrt(log N)),

provided 0<=a<1/4 and K remains bounded. Any separate change in Gaussian
covariance must also be paid. The endpoint construction in the canonical
universality proof takes a=1/7 and pays that change by O(sqrt(epsilon_N)).

Nothing here proves (G), bounds its state complexity, or ensures existence
of a useful R. The actual-near-minimizer counterexample shows that one
natural global operator normalization cannot meet it uniformly by typical
sampling. The implication nevertheless removes a full sign-rounding and
exponential-configuration payment: the hypothesis is a Gaussian PROCESS
maximum with exact sign covariance, not a one-source variance bound.

Dependencies: `flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`
and the campaign's audited balanced-tree almost-subadditivity theorem.
