# Independent audit of the sharp logarithmic port theorem

2026-09-07. **PASS.** I reconstructed both the aggregate matching estimate
and the explicit counterexample in
`principle_director_sharp_port_universality_2026_09_07.md` independently.
This audit does not extend the theorem to actual parent maxima.

For the matching estimate, independently assigning vertex `i` to a bin
with probability `e_ib/(Cm)` makes the selected residual sets disjoint.
For a selected set of size `r>=4D_b+2`, maximality gives internal signed
energy at most `rD_b-r(r-1)/2<=-r^2/4`. Independent global flips on the
disjoint words, and independent spins outside, eliminate every cross
term in expectation. Thus the sum of their squared sizes is at most
`4Q(S)`. Cauchy--Schwarz, followed by expectation of the assigned vertex
count, gives precisely the director's bound

```
E_bad/m^2 <= 4C^2/[V^2(1-(1+delta)^(-2))]
             +2CB/m+2C sqrt(BQ(S))/m.
```

This justifies the improvement to a fixed seed envelope
`q_m log(m)/m^2 -> 0`. The bounded-profile dependency needs only
`Q=o(m^2)`, so no narrower hypothesis is reintroduced there.

For sharpness, take `B=4^b`, `s=4^(B+2)`, and `m=Bs`. The exponent `B+2`
is even. Consequently the positive graph of the hollow block `-H_s`
has even degree `(s-sqrt(s))/2`. The specified sequence
`h_g=2*4^(B-g)` consists of even integers, with `h_1=s/32` below this
degree. Euler orientation and regular bipartite matchings give an
`h_g`-factor without any divisibility assumption beyond evenness.

The hollow perturbation from `T` to `S` is exactly `H_s-J_s` on each
diagonal block: its diagonal is zero. Its all-positive energy is
`B(s^(3/2)-s^2)/2`, and the original `T` energy is `O(m^(3/2))`.
Since `sqrt(m)/B -> infinity`, this proves
`Q(S)=(1/2+o(1))m^2/B`, not merely an upper bound. Also
`B=Theta(log m)` and `Q(T)<=m^(3/2)/2+m/2`.

The row profiles have exact squared energy `h_g a_g^2=m`. Equal nonzero
amplitudes occur only inside their group, where `S=-1`. Different
amplitudes have ratio at least two, so
`(a-b)^2 >= (a^2+b^2)/5`. Summing edge defects counts all directed-port
squares once, proving `D_S>=m^2/5` for every placement. The positive
factors give an actual simultaneous zero-defect placement for `T`.

The probability of the specified factor port sets is exactly the product
of `binom(m-1,h_g)^(-s)`. Its negative logarithm is at most

```
s sum_g h_g log(e(m-1)/h_g) = O(s^2 log B) = o(m^2),
```

because `h_g=(s/32)4^(-(g-1))` and both geometric sums with weights
`1` and `g` converge. Thus the same example separates fixed-temperature
pressures, not only zero-temperature minima. All scales and normalizations
in the stated counterexample pass.
