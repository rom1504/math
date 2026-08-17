# Scout: finite-speed lexicographic count carriers

Status: orthogonal-theory scout; one scoped theorem proved below.  The result
adds genuine subexponential resolution to the one-speed microcanonical
hypograph, but the growing-descriptor falsifier shows that it is **not** a
universal compact state.  Portfolio recommendation: **keep warm; do not
promote** without a tangent-mass composition theorem.

## 1. Literature boundary

Three established theories meet at this question, but none by itself supplies
the desired finite multiscale landscape state.

- Anzellotti--Baldo and
  [Braides--Truskinovsky](https://doi.org/10.1007/s00161-008-0072-2)
  develop variational problems successively at scales
  `1 >> f_1(epsilon) >> ...`.  Braides--Truskinovsky, Sections 1 and 8, also
  identify *locking*, *choking*, nonuniformity in parameters, and the need for
  tables of blown-up Gamma limits near singular parameters.  Thus an iterated
  Gamma limit is not automatically one uniform higher-order state.
- [Akian](https://doi.org/10.1090/S0002-9947-99-02153-4) gives the
  idempotent-measure/large-deviation correspondence at one exponential speed.
  [Mariani](https://arxiv.org/abs/1204.0640), Theorems 3.3--3.5, makes the
  LDP/Gamma-convergence equivalence precise through scaled relative entropy
  and compact-set functionals.  His Theorem 6.1 is a genuine second-order
  Sanov theorem, but its second speed is produced by a particular triangular
  product model; it is not a general finite tower of rate functions.
- The algebra that a finite tower *would* have is standard.  In recent primary
  literature, Friedenberg--Mincheva explicitly use the semifield
  `T^(L)=R^L union {-infinity}` with lexicographic maximum and coordinatewise
  addition; see their 2024
  [Definition 2.3](https://doi.org/10.1007/s40687-024-00467-6).

The theorem below says exactly when this lexicographic algebra is a valid
log-count carrier.  Its sharp boundary is the number of competing
decompositions at the *smallest retained speed*.

## 2. Finite-speed state

Fix an integer `L>=1` and speeds

```math
a_{1,n}\gg a_{2,n}\gg\cdots\gg a_{L,n}\longrightarrow\infty,
\qquad {a_{j+1,n}\over a_{j,n}}\longrightarrow0.       \tag{MS.1}
```

Write `T_L=R^L union {-infinity}`.  Its operations are

```math
u\oplus v=\max_{lex}\{u,v\},\qquad
u\odot v=u+v,                                         \tag{MS.2}
```

with `-infinity` absorbing for `odot`.  A nonnegative sequence `w=(w_n)`
has valuation `nu(w)=u in R^L` if

```math
\log w_n=\sum_{j=1}^L a_{j,n}u_j+o(a_{L,n}),          \tag{MS.3}
```

and `nu(w)=-infinity` if it is eventually zero.  For a finite descriptor
alphabet `K`, the state of multiplicities `A_n:K->N_0` is the profile
`sigma_A(q)=nu(A_n(q))`; equivalently one may store its finite
lexicographic hypograph.

The word *finite* matters topologically.  Lexicographic order is not closed
in the ordinary product topology:
`(1/n,-1)>_lex(0,0)` for every `n`, while the limit `(0,-1)` is smaller than
`(0,0)`.  Thus the scalar compact-hypograph theorem does not automatically
lift by replacing one ordinate with `R^L_lex`.

This state answers not only the leading count question.  If a query has a
matching expansion

```math
V_n(q)=\sum_{j=1}^L a_{j,n}v_j(q)+o(a_{L,n}),          \tag{MS.4}
```

then it retains the response successively at every declared speed.

## 3. The theorem

### Theorem MS.1 (finite lexicographic count algebra, recovery, and sharp uniform branching boundary)

Let the speeds satisfy (MS.1).

1. **Valuation algebra.**  Whenever the displayed valuations exist,

   ```math
   \nu(uv)=\nu(u)\odot\nu(v),\qquad
   \nu(u+v)=\nu(u)\oplus\nu(v).                       \tag{MS.5}
   ```

2. **Exact finite composition.**  Let `K_1,K_2,K` be fixed finite sets and
   `m:K_1 times K_2->K`.  If `A_n,B_n` have valuation profiles `sigma_A,
   sigma_B`, then

   ```math
   C_n(z)=\sum_{m(x,y)=z}A_n(x)B_n(y)                 \tag{MS.6}
   ```

   has the exact asymptotic state

   ```math
   \boxed{
   \sigma_C(z)=
   \mathop{\max_{lex}}_{m(x,y)=z}
       \{\sigma_A(x)+\sigma_B(y)\}.}                 \tag{MS.7}
   ```

   Likewise, for a finite query alphabet,

   ```math
   \nu\!\left(\sum_q A_n(q)e^{V_n(q)}\right)
   =\max_{lex,q}\{\sigma_A(q)+v(q)\}.                \tag{MS.8}
   ```

3. **All-order abstract recovery.**  Conversely, every profile
   `sigma:K->T_L` on a fixed finite `K` is realized, modulo one common
   additive normalization, by integer multiplicities at every sufficiently
   large `n`.  More precisely, choose `C` so that
   `C+sigma_1(q)>0` for every finite entry and put

   ```math
   A_n(q)=
   \left\lfloor
      \exp\!\left(a_{1,n}(C+\sigma_1(q))
             +\sum_{j=2}^L a_{j,n}\sigma_j(q)\right)
   \right\rfloor,                                    \tag{MS.9}
   ```

   with `A_n(q)=0` at `-infinity`.  Then
   `nu(A_n(q))=sigma(q)+(C,0,...,0)`.

4. **Sharp uniform branching boundary.**  For arbitrary positive terms
   `(w_{n,i})_(i in I_n)`,

   ```math
   0\le
   \log\sum_{i\in I_n}w_{n,i}-\max_i\log w_{n,i}
   \le\log|I_n|.                                     \tag{MS.10}
   ```

   Consequently the lexicographic maximum law remains valid, uniformly, for
   growing fibres whenever `log|I_n|=o(a_{L,n})` and the term expansions are
   uniform to `o(a_{L,n})`.  This condition is sharp without further
   structure: `exp(c a_{L,n})` equal maximal terms shift the last coordinate
   by `c`.  This is a worst-case boundary, not a necessary condition for each
   structured family; a much larger fibre can be harmless if almost all of
   its terms are sufficiently suppressed.

#### Proof

For products, add (MS.3).  For sums, let `w,z` have valuations `u,v`,
suppose `u>_lex v`, and let `j` be the first coordinate at which they
differ.  Then

```math
\log w_n-\log z_n
=a_{j,n}(u_j-v_j)+o(a_{j,n})\longrightarrow+\infty,  \tag{MS.11}
```

because every later speed is `o(a_{j,n})`.  Hence
`log(w_n+z_n)=log w_n+o(1)`.  If `u=v`, the maximum of the two
`o(a_{L,n})` remainders is still `o(a_{L,n})`, and the extra `log 2` is
negligible.  This proves (MS.5).

The fibres in (MS.6) and the query alphabet in (MS.8) are fixed and finite,
so iterating (MS.5) proves (MS.7)--(MS.8).  In (MS.9), the common positive
leading coefficient makes the exponent tend uniformly to `+infinity` on the
finite effective domain.  Thus replacing the exponential by its floor changes
its logarithm by `o(1)`, proving recovery.  Finally (MS.10) is the
largest-summand inequality.  Equal summands attain its upper scale, proving
sharpness. `square`

### What is genuinely new relative to the committed one-speed theorem?

At `L=2`, with speeds `(n,sqrt(n))`, the carrier distinguishes one extremal
state from `exp(theta sqrt(n))` extremal states and answers tilts of size
`sqrt(n)`.  The speed-`n` hypograph identifies those two cases.  The exact algebra and
integer recovery persist at all declared levels as long as the number of
decompositions is subexponential at the *last* speed, not merely at the first.

This is nevertheless an algebraic refinement of a response roof, not a new
kind of contextual state.

## 4. Positive benchmark: a finite-type two-speed GREM skeleton

Take a fixed finite acyclic rooted type tree (fixed depth and finitely many
path types).  Give each edge `e` an energy increment
`h_e` and a branching multiplicity

```math
b_{e,n}=\left\lfloor
   \exp\{n s_e+\sqrt n\,t_e\}
\right\rfloor,
\qquad s_e>0.                                         \tag{MS.12}
```

The number of leaves of total energy mark `z` is

```math
N_n(z)=\sum_{p:\,\sum_{e\in p}h_e=z}\prod_{e\in p}b_{e,n}. \tag{MS.13}
```

The tree and its type set are fixed, so Theorem MS.1 gives

```math
\nu(N_n(z))=
\max_{lex,p:\,h(p)=z}
\left(\sum_{e\in p}s_e,\sum_{e\in p}t_e\right).      \tag{MS.14}
```

Product or gluing of two such finite-type skeletons is again governed by
(MS.7).  Thus the state independently recovers the familiar hierarchical
selection rule: maximize the extensive entropy first, then the subextensive
entropy among leading ties.  Its size is linear in the finite type graph, not
in the number of leaves.  Two skeletons with the same `s_e` and different
`t_e` are identical at speed `n` but have different `sqrt(n)` responses.

This is the deterministic count skeleton of a finite-level GREM, not a claim
about the Gaussian extremal process.  Derrida's original
[REM](https://doi.org/10.1103/PhysRevLett.45.79) and
[GREM](https://doi.org/10.1051/jphyslet:01985004609040100) add Gaussian energy
randomness; at an entropy-zero edge their order-one decorations are described
by extremal-process results such as
[Bovier--Klimovsky](https://arxiv.org/abs/0805.1478), beyond any diverging
finite list of speeds.

## 5. Decisive falsifier: an ordinary coefficient hypograph forgets saddle mass

The fixed-alphabet hypothesis cannot be replaced merely by ordinary
product/Hausdorff convergence of coefficient graphs on a compact descriptor
space.

Take `K=[-1,1]`, speeds `(n,log n)`, and two grids

```math
G_n=\{k/n:-n\le k\le n\},\qquad
\widetilde G_n=\{k/\sqrt n:|k|\le\lfloor\sqrt n\rfloor\}. \tag{MS.15}
```

At every supported `q`, assign the integer multiplicity

```math
A_n(q)=\lfloor e^{n(1-q^2)}\rfloor.                   \tag{MS.16}
```

Both grids become dense in `K`; in ordinary product/Hausdorff topology both
coefficient graphs converge to exactly the same smooth two-speed profile

```math
\sigma(q)=(1-q^2,0).                                  \tag{MS.17}
```

Now collapse every descriptor to one point.  On `G_n`,

```math
\sum_{q\in G_n}A_n(q)
=e^n\left(\sqrt{\pi n}+o(\sqrt n)\right),            \tag{MS.18}
```

whereas on `widetilde G_n`,

```math
\sum_{q\in\widetilde G_n}A_n(q)
=e^n\left(\sum_{k\in\mathbb Z}e^{-k^2}+o(1)\right). \tag{MS.19}
```

The floors change either sum only below the displayed scale.  The collapsed
states therefore have valuations `(1,1/2)` and `(1,0)`, respectively, despite
identical limiting coefficient hypographs.  What is missing is the tangent
counting mass of descriptor points in the `n^(-1/2)` saddle window.

This does not rule out a genuinely scale-aware topology: one that demands
descriptor recovery fine enough to make the leading exponent accurate to
`o(log n)` distinguishes the two meshes.  It does rule out an unscaled
limiting coefficient graph as a sufficient second-speed state.

The classical Vandermonde identity gives the same failure in a canonical
counting model without that topological ambiguity.  Fix `p in (0,1)` and take
an integer subsequence on which `pn` is integral.  Stirling's formula gives

```math
\nu\binom n{pn}=(h(p),-1/2),                           \tag{MS.20}
```

but maximizing the pointwise child vectors in

```math
\binom{2n}{2pn}=\sum_k\binom nk\binom n{2pn-k}         \tag{MS.21}
```

predicts second coordinate `-1`, while the exact left side has `-1/2` at
speeds `(n,log n)`.  The missing `+(1/2)log n` is the mass of the
`Theta(sqrt n)` near-saddle decompositions.

Thus even a smooth strictly concave leading profile does not make a bare
pointwise finite-speed hypograph composable by lexicographic maximization.  A
candidate Morse-class repair would retain a local mesh/tangent measure and
Hessian, but no sufficiency or closure theorem for that decoration is proved
here.  Arbitrary landscapes can produce
quadratic, quartic, flat, or fractal exposed fibres, so no fixed finite list of
scalar speeds is closed without a structural regularity class.  This is the
counting analogue of the nonuniform higher-order Gamma tables emphasized by
Braides--Truskinovsky.

## 6. Director judgment

**Do not promote as an orthogonal theory.**  The finite-speed carrier is a
clean and useful theorem, but on its valid class it is precisely the current
response algebra with coefficients upgraded from `T` to `T^(L)`.  Its abstract
recovery theorem, like the one-speed recovery theorem, says nothing about
realization inside codes, graphs, or quadratic signings.

The falsifier isolates one nontrivial possible next theorem:

> Define a tangent-mass decoration on each exposed face and prove a
> composition/recovery theorem for a natural class closed under saddle
> convolution.

Promotion would be warranted only if that decoration is (i) smaller than the
full local counting measure, (ii) closed under repeated composition, and (iii)
validated on a constrained model.  Until then, retain MS.1 as a scoped
subexponential response tool and keep rare-event compactness warm rather than
allocating a second full theory program.
