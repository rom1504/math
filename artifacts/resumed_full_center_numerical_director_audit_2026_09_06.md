# Independent reconstruction of the 0.4310374587290013 certificate

Date: 2026-09-06. The director independently read the new script, its
interval-arithmetic dependencies, the canonical finite-anchor definition,
and the full nonlinear theorem, and reran the exact script successfully.

The canonical covariance c is an exactly specified expression for a
particular degree-200 anchored fixed point, not merely a lower bound on
an unknown optimizing covariance. The finite grouped Hermite sums,
deleted anchor terms, and contraction ratio match that construction.

Write W=cV+dN, d=sqrt(p-c^2), with independent standard Gaussians V,N,
and F=sign(W)1{|V|>alpha}. Conditioning and integration by parts give

    a=E VF=2phi(alpha)(2Phi(c alpha/d)-1)
               +4c phi(0)Phi(-alpha sqrt(p)/d)/sqrt(p),
    b=E NF=4d phi(0)Phi(-alpha sqrt(p)/d)/sqrt(p).

No other Gaussian direction contributes to P1F. Hence tau^2=1-p-a^2-b^2
and J=ca+db. Creation gives E H K=J. Signed weighted Jensen for the convex
function Gamma_tau(k)=E|k+tau Z| under density H/p gives p Gamma_tau(J/p),
without assuming K>=0 or knowing its distribution.

I checked outward rational-grid operations, Machin pi with its alternating
remainder, the finite Gaussian integral and explicit tail, monotone endpoint
evaluation for Phi, and the density's range reduction and five squarings.
Every argument falls in the asserted ranges. Tail subtraction retains
enclosing intervals rather than relying on floating-point cancellation.

The independently reproduced lower endpoints are

    old value:
    .430658179405528602724053804634711026327173238336325190455580
    gain:
    .000379279323472719152390980507389095679106795070238950997741
    sum:
    .431037458729001321876444785142100122006280033406564141453321.

Thus liminf M_n/n^(3/2)>=.4310374587290013 is verified. The JSON saves
both endpoints. Fixed finite construction, matrix limit, approximation,
and principal deletion occur in the full theorem's prescribed order.
Convergence and a matching upper construction remain open; the upper is 1/2.
