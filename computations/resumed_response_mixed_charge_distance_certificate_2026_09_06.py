"""Exact rectangle lower bound on two-step distance from the old frame."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

from fresh_limit_rooted_lower_certificate import phi, Phi
from fresh_limit_mask_ascent_certificate import phi_large


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    mu = 2*Phi(F(1, 2))-1
    b = 2*phi(F(1, 2))
    tau = (1-mu-b*b).sqrt()
    sd = (mu-mu*mu).sqrt()
    assert F(38, 100) < mu.lo < mu.hi < F(39, 100)
    assert F(70, 100) < b.lo < b.hi < F(71, 100)
    assert F(32, 100) < tau.lo < tau.hi < F(37, 100)
    assert F(48, 100) < sd.lo < sd.hi < F(49, 100)
    density_z = phi(F(1, 4))
    density_w = phi_large(F(73, 50))
    density_gate = phi_large(F(5, 2))
    assert density_z.lo > F(3, 8)
    assert density_w.lo > F(1, 8)
    assert density_gate.lo > F(1, 60)
    # Rectangle |Z|<=1/4, 1/2<=W<=3/5 implies H=1 and A=W.
    standardized_w_max = (F(3, 5)+F(39, 100)/4)/F(48, 100)
    assert standardized_w_max < F(73, 50)
    assert 1/F(49, 100) > 2
    # c=b+1/2 is in (1.20,1.21). The Gaussian gate interval
    # ((-A-c)/tau,(A-c)/tau) contains [-5/2,-9/4].
    assert (-F(1, 2)-F(6, 5))/F(37, 100) < -F(5, 2)
    assert (F(1, 2)-F(121, 100))/F(32, 100) > -F(9, 4)
    assert F(3, 5)-F(6, 5) < 0  # gate probability <1/2
    rectangle_mass_lower = F(1, 2)*F(3, 8)*F(1, 10)*F(1, 4)
    gate_probability_lower = F(1, 4)*F(1, 60)
    distance_squared_lower = rectangle_mass_lower*gate_probability_lower/2
    assert distance_squared_lower == F(1, 102400)
    result = dict(
        verified=True, method='exact_fraction_rectangle_density_bounds',
        alpha0='1/2', alpha1='1/2',
        mu=mu.json(), b=b.json(), tau=tau.json(), conditional_W_sd=sd.json(),
        density_z=density_z.json(), density_w=density_w.json(),
        density_gate=density_gate.json(),
        rectangle_mass_lower=str(rectangle_mass_lower),
        gate_probability_lower=str(gate_probability_lower),
        squared_L2_distance_strict_lower=str(distance_squared_lower),
        scope='two-step actual limiting odd-charge half versus complete old canonical sigma-field',
    )
    encoded = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
