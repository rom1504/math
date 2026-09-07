"""Floating exploration of the exact arbitrary-center Gaussian shell rate.

No floating extremum is a directed certificate or a global shell theorem.
"""
import json
import math
from pathlib import Path
from scipy.optimize import brentq, minimize_scalar


def rate(alpha, bridge):
    if alpha == 0 or bridge == 0:
        z = 0.0
    else:
        z = brentq(lambda v: v**3-(1+alpha**2)*v+alpha*bridge,
                   0, math.sqrt((1+alpha**2)/3))
    d = 1-alpha**2
    if d-z*z <= abs(bridge-2*alpha*z):
        return float('inf')
    determinant = (d-z*z)**2-(bridge-2*alpha*z)**2
    return (-.5*math.log(determinant/d**2)
            if determinant > 0 else float('inf'))


def threshold(delta):
    alpha = 1-2*delta
    entropy = -delta*math.log(delta)-(1-delta)*math.log1p(-delta)
    bridge = brentq(lambda b: rate(alpha,b)-2*entropy, 0, .999999999)
    return {'delta': delta, 'alpha': alpha, 'bridge_threshold': bridge,
            'seed_threshold': bridge/(2*math.sqrt(2)-2*alpha**2)}


grid = [threshold(j/2000) for j in range(1,1001)]
maximum = minimize_scalar(lambda delta: -threshold(delta)['seed_threshold'],
                          bounds=(.075,.175), method='bounded',
                          options={'xatol': 1e-13})
result = {'status': 'FLOATING EXPLORATION ONLY',
          'local_maximum': threshold(float(maximum.x)),
          'delta_one_tenth': threshold(.1), 'grid': grid}
out = Path('computations/results/flatify_adversary_2026_09_07_conditional_shell_rate.json')
out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({key: val for key,val in result.items() if key != 'grid'}, indent=2))
