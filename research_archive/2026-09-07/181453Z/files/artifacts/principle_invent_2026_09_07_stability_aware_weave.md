# A stability-aware moment operation for the actual Hadamard weave

2026-09-07. **Exact finite identity and Finner upper bound proved.** No
uniform orthant-probability estimate or improved signing cap is asserted
in this checkpoint. The physical diagonal is kept explicitly throughout.

## 1. Actual matrix, tilt, and physical one-spin stability

For each i in [m], let H_i be k selected rows of an order-m Hadamard
matrix, so H_i H_i^T=m I_k. Let S be a symmetric m by m sign matrix,
whose off-diagonal unordered entries are independent fair signs. The
diagonal S_ii is fixed in advance. On physical vertices (i,a), a in [k],
define the full sign matrix

    W_((i,a),(j,b))=S_ij H_i(a,j) H_j(b,i),

and hollow it as C=W-diag(W). Its physical diagonal is exactly S_ii
on fibre i. For fixed Boolean row words x_i set

    h_i=H_i^T x_i,       u_i=h_i/sqrt(k).

Then

    x^T C x/(2k)=sum_(i<j) S_ij u_i(j)u_j(i)
                   +(1/2)sum_i S_ii[u_i(i)^2-1].           (1)

For an orientation sigma in {+-1}, x is one-spin stable for sigma C iff
sigma x_(i,a)(Cx)_(i,a)>=0 for every physical vertex. Equivalently, the
following event A_i in the signs incident to i holds for each fibre:

    sigma x_(i,a) [sum_(j!=i) S_ij H_i(a,j)u_j(i)
          +S_ii H_i(a,i)u_i(i)-S_ii x_(i,a)/sqrt(k)]>=0
                           for every a in [k].             (2)

This includes both the within-fibre rank-one field and the deletion of
the physical diagonal. Ties are allowed. Every Boolean global maximizer
of sigma x^T C x is stable in this sense.

## 2. Independent tilted edges and the exact Finner correction

Fix H_i and x_i throughout this section. Set

    a_ij=sigma u_i(j)u_j(i),
    D_sigma(u)=exp[t sigma sum_i S_ii(u_i(i)^2-1)].

Define a product probability measure P_(t,u,sigma) on unordered edges by

    P_(t,u,sigma)(S_ij=s)
       =exp(2t a_ij s)/(2 cosh(2t a_ij)),  s in {+-1}.      (3)

In particular its mean is tanh(2t a_ij); the tilts are not generally
small or unbiased. The elementary change of measure gives the identity

    E_S [exp((t/k)sigma x^T C x) 1_(all A_i)]
      =D_sigma(u) product_(i<j)cosh(2t a_ij)
            P_(t,u,sigma)(intersection_i A_i).             (4)

Under (3), each independent edge variable appears in precisely two of
the events A_i. The graph form of generalized Holder/Finner, with vertex
weights 1/2, therefore implies

    P_(t,u,sigma)(intersection_i A_i)
       <=product_i P_(t,u,sigma)(A_i)^(1/2).                (5)

Consequently the old all-vector moment can be replaced by the exact
stability-aware upper integrand

    D_sigma(u) product_(i<j)cosh(2t a_ij)
                   product_i p_i(u,H,x,sigma)^(1/2),       (6)

where p_i is the probability of ALL k inequalities (2) under independent
biased signs with law (3). One proof of (5) is the product-space Finner
inequality for the family of incidence sets, because every edge has sum
of assigned exponents equal to one. Equivalently it is the usual graph
Cauchy--Schwarz inequality for nonnegative vertex functions of incident
independent edge variables. This includes arbitrary biases and repeated
values among u coordinates.

For any cap threshold q, the stable-vector first moment now gives

    P_S(Q(C)>=q) <=sum_(sigma=+-1)sum_x exp(-2tq/k)
        D_sigma(u) product_(i<j)cosh(2t a_ij)
                       product_i p_i^(1/2).               (7)

Indeed a maximizer witnessing Q(C)>=q is stable for its orientation.
The exact same inequality can subsequently be averaged over the bases.
No assertion that the p_i factors decouple under that second average
is implicit in (7).

For N=mk, B=C/sqrt(N), the equivalent inverse temperature is
beta=2t sqrt(N)/k=2t/sqrt(k/m). This checks (7) against the audited
Bellman normalization. The diagonal prefactor D_sigma cannot be dropped
pointwise: a coordinate u_i(i) can carry order-m energy.

## 3. What the new operation does and does not accomplish

The row state is now an actual biased Hadamard orthant probability. It
depends on incoming magnitudes u_j(i), outgoing coordinates u_i(j), the
physical word x_i, and the selected row matrix H_i, not only on a
one-coordinate spectrum type. The operation is therefore a specific
piece of missing feedback, not an arbitrary unnamed control law.

Uniform p_i<=exp(-c k) over all rows is false: incoming fields can be
zero, or strongly aligned with a single Hadamard column, so all the
inequalities can hold simultaneously with probability bounded away from
zero. A useful theorem must either restrict a verified pressure-relevant
region or charge these exceptional configurations through a separate
entropy/energy estimate. An absolute-error high-dimensional central
limit approximation cannot by itself estimate an exp(-c k) event.

A promising sufficient condition is positive total tilted variance and
bounded incoming magnitudes. In that region the sum of negative parts
of the normalized physical fields is a convex Lipschitz function of
independent bounded signs. A positive linear expectation, combined with
convex-product concentration, can yield an exponential orthant cost
without a rectangle CLT. The director is independently deriving this
row lemma and auditing its degeneracy conditions.

## 4. Stronger exact fibre-best-response condition

With all other fibres fixed, write z=H_i(:,i) and

    f_i=sum_(j!=i) S_ij H_i(:,j) h_j(i).

The part of the half-Hamiltonian depending on x_i is exactly

    x_i^T f_i+(S_ii/2)(z^T x_i)^2-S_ii k/2.               (8)

A global maximizer must maximize this rank-one binary quadratic problem
in every fibre. Thus (2) may be replaced by a stronger event enforcing
the full fibre best response, still a function only of the incident edge
signs, and (4)--(7) remain valid. The scalar overlap z^T x_i takes k+1
values, giving an exact finite threshold/overlap description of the
conditional fibre optimum. This alone is NOT a bound of polynomial size
on the number of global fixed points: the external fields f_i depend
on every other fibre. No global entropy saving is claimed from this
parameterization.
