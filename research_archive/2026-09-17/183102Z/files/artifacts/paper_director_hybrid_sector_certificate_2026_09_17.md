# A full-parent hybrid fourth-mass certificate

2026-09-17. Director derivation, **pending independent audit**. This
combines the common-Gibbs/truncation comparison with the capped
fourth-mass Laplace envelope. Neither alone gives the stated hybrid
certificate. It is a finite theorem about genuine full sign parents,
not a theorem that its right side is favorable for exact minimizers.

## 1. Exact physical model and statement

Let A,D be arbitrary hollow full sign matrices of orders n,q. Write
n=kp+ell, 0<=ell<k. Fix ANY list h_1,...,h_q in {+-1}^k.
The j-th bridge column, restricted to the k by p array, is h_j tensor
g_j, where all entries of g_j are independent fair signs. Its ell
remaining entries are also independent fair signs. All drivers are
independent across columns. Thus every physical bridge entry is an
exact sign; it is not independently rounded across all physical edges.

For an old spin array x define

```math
V(x)=\sum_{j=1}^q\sum_{t=1}^p(h_j^Tx_{:,t})^2+q\ell,
\qquad
W(x)=\sum_{j=1}^q\sum_{t=1}^p(h_j^Tx_{:,t})^4+q\ell.
```

These are the exact second/fourth coefficient sums for EVERY new
spin y and either polarity. Put psi(u)=u^2/2-log cosh u. Partition
the old words into L={x:W(x)<=w0} and U=its complement. Empty
sectors may be omitted. Let Q_L^G be the maximum over both polarities,
x in L and all y of the same parent energy after replacing the scalar
drivers by independent standard Gaussians. The children remain fixed.
Define, for t>0,

```math
\begin{split}
G_L&=\mathbb E Q_L^G+
 8w_0^{1/4}[(n+q+2)\log2]^{3/4},\\
B_U(t)&=\frac1t\log\sum_{s=\pm1}
 \left(\sum_{y\in\{\pm1\}^q}e^{tsH_D(y)}\right)
 \left(\sum_{x\in U}
 e^{tsH_A(x)+t^2V(x)/2-[W(x)/k^4]\psi(tk)}\right),\\
L_*&=qp k^2+q\ell.
\end{split}
```

There exists a realization of the exact full-sign bridge with

```math
\boxed{
Q\begin{pmatrix}A&C\\C^T&D\end{pmatrix}
 \le \max\{G_L,\inf_{t>0}B_U(t)\}
       +\sqrt{2L_*\log2}.}                         \tag{1}
```

The error joining the two sectors is O(sqrt(qnk)). In particular,
q=Theta(n), k=O(sqrt(n)) makes it o(n^(3/2)). Both children and all
their energies are retained, and the infimum is over one common
temperature for the ENTIRE coherent sector. There is no choice of
a separate favorable temperature or signing for each query.

## 2. Proof

For each witness (s,x,y), its scalar coefficient at driver (j,t) is
s y_j h_j^T x_{:,t}. The leftover coefficients have magnitude one.
Its offset is exactly s[H_A(x)+H_D(y)]. The diffuse-sector comparison
therefore gives E Q_L<=G_L, with at most 2^(n+q+1) witnesses. The
Gaussian model is the SAME physical frame, not iid edge disorder.

For the coherent sector, the capped fourth-mass inequality gives

```math
\mathbb E e^{t s x^TCy}
 \le e^{t^2V(x)/2-[W(x)/k^4]\psi(tk)}.
```

The maximum is bounded by log-sum-exp; Jensen moves expectation
outside its logarithm. Summing y and the two polarities, without
altering the child offsets, gives E Q_U<=B_U(t). Taking the infimum
after this expectation bound is valid.

Changing one array driver changes either restricted maximum by at
most 2k, and changing a leftover driver changes it by at most2.
The bounded-differences MGF inequality therefore gives

```math
\log\mathbb E e^{a(Q_S-\mathbb E Q_S)}\le a^2L_*/2
\quad(S=L,U).
```

No independence between sectors is needed. Log-sum-exp over these
two centered quantities gives
E max(Q_L,Q_U)<=max(E Q_L,E Q_U)+sqrt(2L_*log2).
Since the full parent cap equals this maximum pointwise, some actual
bridge realizes at most its expectation. This proves (1).

## 3. Critical normalization and exact remaining gap

For balanced Hadamard lists at comparable orders and k~lambda sqrt(n),
V(x)<=(q+k)n and L_*=O(n^(5/2)). If w0=omega n^3,
the diffuse comparison costs O(omega^(1/4)n^(3/2)); the joining cost
is O(n^(5/4)). Coherent words pay the negative pressure correction
W(x)psi(beta k/sqrt(n))/k^4 at t=beta/sqrt(n).
Unlike the uncapped one-spike envelope, this correction is Theta(n)
when W(x)=Theta(n^3). No divisibility hypothesis is used in (1).

To get a favorable seed transfer, one must bound BOTH the actual
Gaussian-sector maximum and the coherent weighted child partition sum
by the desired seed value plus a sufficiently small increment. Neither
quantity is bounded here by simply assuming a favorable W value.
The correction may be substantial, but there is no claimed uniform
improvement for exact optimizing children. Global child reversal still
gives |H_A+H_D|+|x^TCy|; the certificate respects it by retaining every
y and both polarities.

## 4. Dependencies and interpretation

- Common-Gibbs/truncation theorem:
  paper_symmetric_frame_universality_2026_09_17.md, Sections1--3.
- Capped fourth-mass theorem:
  paper_director_capped_fourth_mass_2026_09_17.md, Sections1--2.
- Elementary bounded differences; the constants are derived above.

The original full-sign constraint is exact. The intermediate Gaussian
model is only a comparison device. There is no claim that covariance
matching preserves the Boolean maximum on the coherent sector.
