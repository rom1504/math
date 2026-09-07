"""Finite exact-energy regression for the continuous block-rotation identity."""
import json
import math
import numpy as np

rng = np.random.default_rng(20260907)
records = []
for n in range(2, 9):
    spins = 1-2*((np.arange(1 << n)[:, None] >> np.arange(n)) & 1)
    for trial in range(3):
        a = rng.choice([-1, 1], size=(n, n))
        a = np.triu(a, 1)
        a += a.T
        energy = np.einsum('bi,ij,bj->b', spins, a, spins)//2
        cross = spins @ a @ spins.T
        k = math.sqrt((2*n-1)/(n-1))
        initial = k*(energy.max()-energy.min())
        for theta in (0, math.pi/16, math.pi/8):
            alpha, beta = k*math.cos(2*theta), k*math.sin(2*theta)
            values = alpha*(energy[:, None]-energy[None, :])+beta*cross
            cap = np.abs(values).max()
            assert cap+1e-9 >= math.cos(2*theta)*initial
            assert cap <= (1+math.sin(2*theta))*initial+1e-9
            ix, iy = np.unravel_index(np.abs(values).argmax(), values.shape)
            x, y = spins[ix], spins[iy]
            u, v = (x+y)//2, (x-y)//2
            hu, hv, z = u@a@u/2, v@a@v/2, u@a@v
            assert abs(values[ix, iy]-2*(beta*(hu-hv)+alpha*z)) < 1e-9
            if n == 2:
                assert abs(cap-math.cos(2*theta)*initial) < 1e-9
            records.append({"n": n, "trial": trial, "theta_over_pi": theta/math.pi,
                            "cap_ratio_DISPLAY": cap/initial})
print(json.dumps({"status": "ROTATION_IDENTITY_AND_BOUNDS_PASS",
                  "tests": len(records), "records": records}, indent=2))
