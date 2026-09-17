"""Exploratory, not a theorem: running precision on oscillatory power profiles."""
import numpy as np


def evaluate(period, phase, dt=0.001):
    t = np.arange(-30.0, 30.0 + dt / 2, dt)
    # Extend the integral of the 0/1 square-wave slope to negative t.
    z = (t + phase) / period
    cycles = np.floor(z)
    remainder = (t + phase) - cycles * period
    primitive = cycles * period / 2 + np.maximum(remainder - period / 2, 0)
    log_e = -primitive
    log_ref = -((np.floor(phase / period) * period / 2)
                + max(phase % period - period / 2, 0))
    ratio = np.exp(log_e - log_ref)
    u = np.exp(-t)
    scaled_b = 2 * u / (u + 1) ** 2 / ratio
    running = np.maximum.accumulate(scaled_b)
    increments = np.diff(np.concatenate(([0.0], running)))
    return float(ratio @ increments)


for p in (0.05, 0.1, 0.3, 0.6, 1, 2, 4, 6, 10):
    vals = [evaluate(p, phase) for phase in np.linspace(0, p, 80, endpoint=False)]
    print(p, min(vals), sum(vals) / len(vals), max(vals))
