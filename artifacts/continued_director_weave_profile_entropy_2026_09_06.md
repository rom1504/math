# Spectrum-type entropy in a random exact Hadamard weave

Date: 2026-09-06. Director derivation. The exact counting inequality has
been independently reconstructed by the convergence agent. The positive-gap
extension and primary-theorem mapping have also passed the audit agent. No upper bound for
the unrestricted Boolean cap follows from the scoped applications.

## 1. The exact signing and its equality condition

Let H be any order-m real Hadamard matrix. Independently choose uniform
column permutations P_i, and put H_i=H P_i, for 1<=i<=m. Independently
choose signs S_ij=S_ji for i<j, with any fixed diagonal signs. The matrix

```math
K_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i)
```

is symmetric, has sign entries, and satisfies K²=m²I. For block spins
x_i in {−1,1}^m let h_i=H_i^T x_i. Direct expansion gives

```math
x^T Kx=\sum_{i,j}S_{ij}h_i(j)h_j(i),\qquad
\sum_i\|h_i\|_2^2=m^3,
```

and, for sigma in {−1,1},

```math
2(m^3-\sigma x^TKx)
=\sum_{i,j}\big(h_i(j)-\sigma S_{ij}h_j(i)\big)^2.       (1)
```

Thus exact saturation is an equality of actual transformed coordinates,
not a separately paid estimate on channels. Hollowing K changes every
Boolean quadratic value by at most m², lower order relative to m³.
This construction nevertheless always has matching Boolean eigenvectors;
the present theorem eliminates OTHER declared classes, not those witnesses.

## 2. An exact entropy bound for any fixed block spins

Fix x independently of the random permutations and S. For g_i=H^T x_i,
let c_(i,a) be the multiplicities of its distinct absolute values. Define

```math
t_i=|\{j:g_i(j)\ne0\}|,\qquad
N_i=\frac{m!}{\prod_a c_{i,a}!},\qquad
k_i=|\{a:c_{i,a}>0\}|.
```

Then

```math
\Pr\{|x^TKx|=m^3\}
\le 2\prod_{i=1}^m
 \frac{\sqrt{2k_i}}{\sqrt{N_i\,2^{t_i}}}.                (2)
```

The right side can exceed one; it is then simply uninformative.

Proof. The absolute row words |h_i| are independent and uniform among
N_i permutations of their multisets. First fix their diagonal symbols d_i.
The possible off-diagonal row words number N_i c_(i,d_i)/m. The number of
symmetric off-diagonal magnitude arrays with these row restrictions is at
most the square root of the product of these numbers. To see this, take
the uniform distribution on such arrays and apply the entropy chain rule
to unordered edges. The sum of row entropies is at least twice the joint
edge entropy: each edge occurs in exactly two rows, and conditioning on
earlier edges outside a row can only reduce entropy. Each row entropy is
at most the logarithm of its allowed number of words.

Sum over diagonal symbols. The sum factors and
sum_d sqrt(c_(i,d)/m)<=sqrt(k_i), so the number of compatible magnitude
arrays is at most prod_i sqrt(N_i k_i). Divide by prod_i N_i.

On any such array, the number of nonzero unordered off-diagonal edges is
at least (sum_i t_i-m)/2. Each imposes one independent constraint on S_ij
for a fixed sigma, regardless of the signs already attached to h_i.
The conditional probability is at most 2^{-(sum t_i-m)/2}. Finally union
over the two sigma values. This proves (2).

For any declared row family X_m, a direct union bound consequently gives

```math
\Pr\{\exists x_i\in X_m:\ |x^TKx|=m^3\}
\le 2\left[
 \sum_{u\in X_m}\sqrt{\frac{2k(u)}{N(u)2^{t(u)}}}
 \right]^m.                                            (3)
```

This is a one-row counting obligation on m spins, not full optimization
on m² spins. It does not imply that the bracket is small for all rows:
the 2m signed Hadamard rows alone prevent that conclusion.

## 3. Positive gaps for exact plateaued row classes

Suppose H is the Walsh matrix of order m=2^d, and X_m consists of
s-plateaued functions: every transformed coordinate is 0 or ±sqrt(m/q),
where q=2^{-s}. Thus t=qm, N=binomial(m,qm), and k<=2. Assume s is fixed
and d+s is even. Write |X_m|<=2^{(c_q+o(1))m} and let

```math
g_q=\frac{H_2(q)+q}{2}-c_q.
```

If g_q>0, then there is a fixed eta>0 such that, with probability tending
to one, EVERY x with all rows in X_m satisfies

```math
|x^TKx|<(1-\eta)m^3.                                   (4)
```

Here is a quantitative sufficient choice. For q<1 put k=2; for q=1 put
k=1. Any eta with 2 eta q<1/2 and

```math
\frac12H_2(2\eta q)+\eta q\log_2(2k^2)<g_q             (5)
```

works, with exponentially small failure probability on the m² scale.

Proof. A failed unordered equality in (1) contributes at least 2m/q to
its right side. An eta-near-saturator therefore has at most eta q m² bad
unordered edges. Fix a set E of e such edges. On all other edges magnitudes
must agree. Expose the two endpoint magnitudes on E and all diagonal
magnitudes, at cost at most k^{2e+m}. On the remaining undirected graph
the row-entropy argument above gives at most prod_i sqrt(N_i) arrays.
This bound is deliberately loose but uniform. Dividing by prod_i N_i,
and imposing the independent signs on the nonzero good edges, bounds the
probability by

```math
k^m2^{m/2}\prod_i(N_i2^{t_i})^{-1/2}(2k^2)^e.
```

The number of constrained nonzero good edges is at least
(sum t_i-m-2e)/2. Sum over the at most binomial(m,2) possible bad-edge
locations and union over all |X_m|^m row choices and both sigma. The
logarithm of the resulting probability, divided by m², is at most

```math
-g_q+\frac12H_2(2\eta q)
       +\eta q\log_2(2k^2)+o(1).
```

This proves (4)-(5). The row class is EXACTLY plateaued. A spin row close
in an unspecified spectral sense to this class is not covered.

## 4. Exact mapping of a primary counting theorem

Potapov, *Upper bounds on the numbers of binary plateaued and bent
functions*, arXiv:2303.16547v3 (18 November 2024), Theorems 1(a) and 2,
uses d Boolean variables, unnormalized Walsh coefficients, and binary
logarithms. Thus his landscape length is our m=2^d. His results give

```math
c_1=11/32,
\qquad c_q=\alpha/8+(H_2(q)+q)/4\quad(s>0\text{ fixed}),
\qquad\alpha=1+(3/8)\log_2 6.
```

For fixed s, the Hamming-ball term in Theorem 1(a), divided by 2^d,
tends to 1/8; the relevant d+s parity must be retained. These imported
bounds make g_q positive at q=1,1/2,1/4. They do not do so at q=1/8.
In particular (4) applies separately to bent rows, to 1-plateaued rows
at odd d, and to 2-plateaued rows at even d.

One common explicit choice is eta=1/1000. This does not depend on floating
logarithm evaluation: alpha<2 and log_2 3<8/5 give g_(1/4)>1/80,
g_(1/2)>1/8, while g_1=5/32. The elementary bound
H_2(z)<=z log_2(e/z), with e<3, bounds the respective penalties in (5)
by 1/250, 3/400, and 3/250, all strictly smaller. Thus the scoped gaps
can be taken as 0.001 of the full spectral energy, at all sufficiently
large allowed Walsh orders, without a numerical asymptotic claim.

[Primary theorem and full proof](https://arxiv.org/html/2303.16547v3#S5).
The counting proof relies on exact Walsh divisibility and algebraic
degree bounds. No quantitative near-bent stability theorem is imported.

## 5. Relation to the original question

This theorem explains quantitatively why dense exact plateaued witnesses
need not survive random weaving, even though matching witnesses always
do. It is not an upper cap certificate and does not improve M_n: sparse
spectra, nonplateaued spectra, and near-saturators outside these classes
are unbounded here. For restricted weaves the coordinates arise from
partially specified spin rows, not the exact plateaued classes above.

The next relevant obligation is a stable spectrum-profile count covering
ALL high-energy rows, with its entropy cost strong enough to control
the squared defect (1). Eliminating a list of named spectral classes is
not a substitute for that obligation.
