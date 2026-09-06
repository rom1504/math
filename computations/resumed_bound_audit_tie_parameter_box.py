"""Exact rational enclosures for the finite causal zero-gradient tie example."""

from fractions import Fraction as F
from math import factorial

# t=(1/4) integral_0^2 exp(-u^2/2)du. The integrated terms decrease.
ts = [sum((F((-1) ** k * 2 ** (k + 1), (2 * k + 1) * factorial(k))
           for k in range(n + 1)), F()) / 4 for n in (16, 17)]
# The tails of the alternating exponential series are decreasing here.
es = [sum((F((-2) ** k, factorial(k)) for k in range(n + 1)), F())
      for n in (18, 19)]
tl, th = min(ts), max(ts)
el, eh = min(es), max(es)
assert F(29907, 100000) < tl < th < F(29908, 100000)
assert F(13533, 100000) < el < eh < F(13534, 100000)

rl = (1 - 3 * th) / (1 + (1 - el) / (4 * tl))
rh = (1 - 3 * tl) / (1 + (1 - eh) / (4 * th))
pl, ph = tl + rl, th + rh
al, ah = rl / (4 * th), rh / (4 * tl)
bl, bh = 1 - ah, 1 - al
assert F(358, 1000) < pl < ph < F(360, 1000)
assert F(377, 1000) < pl / bh < ph / bl < F(379, 1000)
assert F(256, 1000) ** 2 < rl / bh ** 2 < rh / bl ** 2 < F(259, 1000) ** 2

# phi(0)>.398 and phi(0)<.4 follow from classical rational pi bounds.
assert F(398, 1000) ** 2 < F(7, 44)  # pi < 22/7
assert F(23, 24) * F(398, 1000) > F(379, 1000)
assert F(1, 2) * F(379, 1000) < F(256, 1000)
assert F(3, 4) * F(377, 1000) > F(259, 1000)

# J>0.31, whereas the pure-noise cap is <0.308.
assert (4 * F(398, 1000) * F(358, 1000)) ** 2 * F(299, 1000) > F(31, 100) ** 2
assert F(8, 27) / F(157, 50) < F(308, 1000) ** 2  # pi > 3.14
print("PASS: exact parameter, annulus-existence, and strict pure-noise separation boxes.")
