# Elementary dephased-Hadamard bridge: all fixed-center energy windows

2026-09-07. Root's stronger observation independently reconstructed PASS.
This supersedes the need for Gaussian-profile or exact-support assumptions
in the previous FIXED-CENTER CHILD-ENERGY-WINDOW conclusions. Those earlier
kernel theorems remain correct, but are unnecessarily restrictive for
this particular operation. No original convergence theorem follows.

## Exact elementary identity

Let H be any normalized Hadamard matrix of order n. Its first row and
column are all ones, so H1=H^T1=n e1. For arbitrary physical centers
x0,y0, define the full sign bridge C=diag(x0) H diag(y0). Its operator
norm is sqrt(n), and x0^T C y0=n=o(n^(3/2)).

For arbitrary target spins put u=diag(x0)x, v=diag(y0)y and
alpha=<x,x0>/n, beta=<y,y0>/n. Write r=u-alpha*1, s=v-beta*1. Exactly,

    u^T H v = r^T H s+n*(alpha*v1+beta*u1-alpha*beta).

Since ||r||^2=n(1-alpha^2), ||s||^2=n(1-beta^2),

    |x^T C y| <= n^(3/2)*sqrt((1-alpha^2)(1-beta^2))+3n.       (1)

This holds simultaneously for ALL physical x,y, with no probabilistic
ensemble, profile approximation, entropy count, or optimizer hypothesis.
It also handles both bridge polarities and arbitrary signs of the overlaps.

## Exact original-energy consequence, including unequal noise levels

Suppose the ORIGINAL child energies of a candidate satisfy

    |H_A(x)| <= c_A*alpha^2*n^(3/2)+o(n^(3/2)),
    |H_D(y)| <= c_D*beta^2*n^(3/2)+o(n^(3/2)),

with c_A,c_D<=1/2. The errors must be uniform over the stated windows.
For a=alpha^2,b=beta^2, arithmetic-geometric mean gives

    c_A*a+c_D*b+sqrt((1-a)(1-b))
      <= c_A*a+c_D*b+1-(a+b)/2 <= 1.

Therefore the actual parent expression is at most

    |H_A(x)+H_D(y)|+|x^T C y| <= n^(3/2)+o(n^(3/2)).          (2)

For equal seed constant c this is strictly below 2sqrt(2)c*n^(3/2)
whenever c>1/(2sqrt(2)), in particular throughout the currently possible
actual-minimizer interval. It covers unequal noise levels as well as
equal ones. Actual child optimality is not needed for (1)--(2); it is used
only to identify the relevant cap scale and typical noisy-energy windows.

## The pinned rank-two ensemble already has the same property

In the pinned ensemble used in the previous hybrid proofs, every right
center pin points to one of the two reserved LEFT fibres. Consequently
C y0 is supported on those 2k physical coordinates. Similarly C^T x0 is
supported on the two reserved right fibres. Since ||C||op=sqrt(n),

    ||C y0||2=||C^T x0||2=n,
    ||C y0||1,||C^T x0||1 <= n*sqrt(2k).

The center bilinear form is exactly zero. Splitting x,y into center
components and orthogonal residuals now gives for EVERY output of that
ensemble, not merely a favorable random output,

    |x^T C y| <= n^(3/2)*sqrt((1-alpha^2)(1-beta^2))
                  +2n*sqrt(2k).

As k=sqrt(2n), the extra term is O(n^(5/4)). Thus the entire Gaussian and
G/P coupling apparatus was unnecessary to obtain its fixed-center
energy-window saving. It is the localized center fields that pay those
windows. This observation does not invalidate the exact kernel identities
or the separately proved exceptional-profile counts.

## All original orders

Normalized Hadamard orders h=2^a12^b are multiplicatively dense. For every
integer n choose h>=n with h/n->1. Pad arbitrary centers and candidates
by fixed signs, use the dephased normalized Hadamard bridge of order h,
and restrict to the original n-by-n rectangle. The uniform bilinear
restriction error is at most 2h*sqrt(h-n)=o(n^(3/2)). Padded overlaps
differ from original overlaps by o(1); continuity of the square-root
factor is uniform on [-1,1]^2. Hence (1) holds at ALL orders with o(n^(3/2))
in place of 3n, and (2) follows for the original, unchanged children.
No change to any child edge or child energy is made.

## Exact remaining scope and why this is not recurrence

For one chosen center, a general physical spin need not have child energy
bounded by the center-overlap-squared window. Other extrema can have large
child energy at small overlap. Typical noisy-energy concentration around
the chosen center does not bound every such exceptional spin, and a
separate bridge for every center does not give one bridge for all centers.

Thus the remaining full-parent obligation is an energy/center covering
problem, not a need to remove Fourier-profile assumptions within these
windows. The rare high-bridge examples from unconcentrated rank-two
ensembles do not contradict (1), since those bridges were not dephased
around the relevant centers and lacked localized center fields. The
broader weighted flatification task, which allows internal-edge changes,
remains open as well.
