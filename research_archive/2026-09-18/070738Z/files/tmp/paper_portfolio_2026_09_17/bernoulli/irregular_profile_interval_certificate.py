"""Exploratory rigorous interval lower bound for an abstract profile functional.

No actual-signing realizability claim. mpmath interval arithmetic encloses all
rational operations; the final square root only covers the phase grid.
"""

import mpmath as mp

mp.iv.dps = 40
r = mp.iv.mpf(101) / 100
half = 70
period = 2 * half
lo = -15 * period
hi = 15 * period
eta = []
e = []
for i in range(lo, hi + 1):
    cycles, rem = divmod(i, period)
    exponent = -half * cycles - max(rem - half, 0)
    eta.append(r ** (-i))
    e.append(r ** exponent)
weights = []
for i in range(len(e) - 1):
    lower = (e[i] - e[i + 1]).a
    weights.append(max(mp.iv.mpf(0), lower))

minimum = None
min_phase = None
for phase in range(period):
    eps = r ** (-phase)
    best = mp.iv.mpf(0)
    lower_sum = mp.iv.mpf(0)
    for idx, weight in enumerate(weights):
        s = eta[idx] / e[idx]
        b = 2 * eps * s / (eta[idx] + eps) ** 2
        best = max(best, b.a)
        lower_sum += best * weight
    lower = lower_sum.a
    if minimum is None or lower < minimum:
        minimum = lower
        min_phase = phase

all_phase_lower = (minimum / mp.iv.sqrt(r)).a
assert all_phase_lower > mp.iv.mpf(9) / 8
print("phase_grid_lower", minimum)
print("worst_grid_phase", min_phase)
print("all_phase_lower", all_phase_lower)
print("target", mp.iv.mpf(9) / 8)
print("PASS: abstract periodic profile resource is uniformly above 9/8")
