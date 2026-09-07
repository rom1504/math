# Independent audit: iterated-selector packet tests

2026-09-07. **PASS**, Sections1--5 of
`principle_invent_2026_09_07_iterated_selector_packet_falsifier.md`.

For the full tensor transform, a pure rare packet of size s occurs with
probability theta=kr/s. Its coarse output amplitudes are -2C_j/sqrt(k),
so the clipping sum W gives normalized improvement
p(W-1)/(4 tau s). A singleton has W=4tau/p before saturation, giving
1-p/(4tau). A product of positive half-supports has d=2^J nonzero
coefficients of magnitude s=L/d. Its saturated response is exactly
p d(d-1)/(4tau L). For two32 stages the singleton lower threshold lies
well beyond product saturation; the necessary compatibility fails since
32>93/8. The finite selector likelihood of a chosen admissible rare
packet contributes only O(r), not another r log(1/r) term.

For the DC-spine orthogonal transform F, the output amplitudes instead
are -2 F1_S/sqrt(p). Thus lambda=4tau/p and the response is precisely
(W_lambda-1)/(lambda s). At stage ell the product-half packet's
coefficient magnitude is sqrt(prod_(h<=ell)L_h)/2^ell. The multiplicity
of its one non-DC character is prod_(h>ell)L_h/2^(J-ell). The last stage
adds its DC and one ordinary coefficient. This reconstructs d_star.

In the two32 case the first stage contributes16 atoms and the last
stage2, all clipped for lambda>=1; s=256. Singleton energy is1 with
maximum squared coefficient1/32. Therefore any strict-subhalf common
temperature must have lambda>32/31 and lambda<17/8. Their ratio is
527/256 exactly.

These are lower witnesses for the finite block ENVELOPE supremum, not
Boolean cap lower bounds. They test the specified tensor and DC-spine
selector models, not unspecified feedback or fresh-randomized iterations.
