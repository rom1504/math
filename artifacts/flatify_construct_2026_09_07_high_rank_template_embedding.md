# High-rank quadratic families contain every fixed Boolean-compatible template

Date: 2026-09-07. Status: **Pending independent audit**. This elementary finite
linear-algebra theorem was suggested by the director and reconstructed by the
construction agent. It is an obstruction mechanism for generic quadratic-MUB
families, not a favorable original-value recurrence.

## 1. Prescribing all polar forms on one common subspace

Let W=F2^n, and let B_1,...,B_k be alternating bilinear forms on W. Suppose
every nonzero binary combination has rank at least R. Fix an integer d>=1
and suppose

```math
                  R>2k(d-1)+k.                           (1)
```

Then for arbitrary alternating forms S_1,...,S_k on F2^d, there is an
injective linear map L:F2^d->W with L^*B_i=S_i. Moreover L can be selected
so that the kd covectors B_i(L e_j,.) are linearly independent.

**Proof.** Choose u_j=L e_j inductively. At stage j, let F be the span of
the k(j-1) old covectors; they are independent by induction, so r=dim F=k(j-1).
The desired values B_i(u_l,u_j)=S_i(e_l,e_j) for l<j prescribe r independent
linear equations, leaving exactly2^(n-r) choices for u_j.

The new k covectors fail to be independent modulo F precisely if, for some
nonzero c in F2^k, the covector B_c(u_j,.) lies in F. The map u->B_c(u,.) has
rank at least R, so the inverse image of F has size at most2^(n-R+r).
There are fewer than2^k bad nonzero combinations. By (1),

```math
       2^k 2^{n-R+r}<2^{n-r}.
```

Thus some solution to the prescribed equations avoids every bad event.
The enlarged family of kj covectors is independent. This also makes u_j
independent of the old vectors, since otherwise every B_i(u_j,.) would
already lie in F. Alternation supplies all diagonal zeros and the reverse
pair values. Induction proves the claim.

## 2. One affine translation prescribes every linear part

Let Q_i be quadratic functions with polars B_i, and let q_i be arbitrary
quadratic functions on F2^d with polars S_i. After section1, the difference
between Q_i(Lz) and q_i(z) is affine linear. Since all kd covectors
B_i(u_j,.) are independent, choose one p in W solving the kd equations
which correct all k linear parts simultaneously. Therefore

```math
             Q_i(p+Lz)=q_i(z)+c_i                         (2)
```

for some constants c_i. Constants are immaterial for Boolean flatness.
There is no demand that p lie outside the image of L.

## 3. Exact Boolean-template embedding

Put D=2^d and t=2^n. Let H_D,H_t be Walsh matrices with no normalization.
For arbitrary Boolean f_i in {+1,-1}^D define unit vectors

```math
 w_i=D^{-1}\operatorname{diag}((-1)^{q_i})H_D f_i.
```

Their unit norms follow from H_D H_D^T=D Id. Embed the w_i at the common
affine coordinates p+Lz of W, setting all other coordinates to zero, and
write the embedded vectors as v_i. Equation (2) gives

```math
 x_i=H_t\operatorname{diag}((-1)^{Q_i})v_i
                   \in\{+1,-1\}^{t}.                    (3)
```

Indeed restriction of an ambient Walsh character to p+L(F2^d) gives a
constant sign times a character of F2^d; all D local characters occur, and
the local Walsh sum in (3) is f_i at that character. Every pair Gram product
of the v_i equals the original pair Gram product of the w_i exactly.

If U_i=H_t diag((-1)^Q_i)/sqrt(t) are MUBs and A is any hollow seed signing,
the actual cross blocks sqrt(t)a_ij U_iU_j^T therefore admit Boolean energy

```math
          t^{3/2}\sum_{i<j}a_{ij}\langle w_i,w_j\rangle.   (4)
```

This statement requires the actual finite Boolean-compatible templates
above. It does NOT substitute an arbitrary spherical or SDP vector family.

## 4. All within-fibre completion terms disappear under an exact finite law

Apply the same translation and character modulation
(T_(a,b)v)(z)=(-1)^(b dot z)v(z+a) to every v_i. Translation changes the
quadratic restrictions only by affine linear functions; modulation does
the same. In the local Walsh sum, these changes only permute the Boolean
f_i and change an overall sign. Thus every transformed physical vector
remains Boolean. Pair Gram products remain unchanged.

Uniform averaging over (a,b) makes each spectral covariance Id/t, and
therefore each physical covariance Id. Every hollow within-fibre quadratic
completion has zero expected energy. Some actual Boolean assignment in
this orbit consequently attains at least (4), no matter how those full
diagonal-fibre signs were chosen.

## 5. Scope and relation to the explicit Kerdock family

For the explicit field forms in
`flatify_construct_2026_09_07_kerdock_plateaued_witness.md`, the sum over any
r selected parameters differs from B_(sum parameters) by a matrix of rank
at most r. Indeed all remaining terms are combinations of the rank-one
covectors x->Tr(a_i x). If the parameters are binary-linearly independent,
every nonempty sum is nonzero, so the combined polar rank is at least n-k.
For fixed k,d, condition (1) follows by choosing n sufficiently large.

Thus growing independent-parameter Kerdock families contain all fixed
Boolean-compatible quadratic-Walsh templates. Merely arranging small bias
for every fixed nontrivial phase product does not eliminate adapted
finite-dimensional spectral witnesses.

The separate pentagon theorem is stronger for that particular seed: its
two-dimensional witness requires only nonsingular pair differences, with
no high-rank condition on larger sums. Neither theorem claims an upper
bound or original nonconvergence.
