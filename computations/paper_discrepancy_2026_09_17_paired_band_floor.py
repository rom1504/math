"""Exact matching-pair sharpness replay and diagnostic spectral thresholds."""
from fractions import Fraction as F
import json
import math


def convolution_step(coefficients, kernel):
    output = [0]*(len(coefficients)+2)
    for i, value in enumerate(coefficients):
        for j, weight in enumerate(kernel):
            output[i+j] += value*weight
    return output


plus = [1]
minus = [1]
records = []
targets = {1,2,4,8,16,32,64,128,256,512}
for m in range(1,max(targets)+1):
    plus = convolution_step(plus,[1,1,1])
    minus = convolution_step(minus,[1,4,1])
    assert sum(plus)==3**m and sum(minus)==6**m
    assert plus==plus[::-1] and minus==minus[::-1]
    assert sum((j-m)**2*v for j,v in enumerate(plus)) == 2*m*3**(m-1)
    assert 3*sum((j-m)**2*v for j,v in enumerate(minus)) == m*6**m
    if m in targets:
        eplus = F(sum(abs(j-m)*v for j,v in enumerate(plus)),3**m)
        eminus = F(sum(abs(j-m)*v for j,v in enumerate(minus)),6**m)
        # Physical sum equals twice this lazy walk; average the two sectors.
        exact_unnormalized = eplus+eminus
        records.append({"n":2*m,"mean_abs_exact":str(exact_unnormalized),
                        "normalized_diagnostic":float(exact_unnormalized)/math.sqrt(2*m)})
kappa=math.sqrt(2/math.pi)
floor=(math.sqrt(2)+1)/math.sqrt(3*math.pi)
# Exact, coarse proof that this floor exceeds 3/4:
# sqrt(2)>7/5 and pi<22/7 imply floor^2>203/330>9/16.
assert F(203,330)>F(9,16)
thresholds=[]
for cap in (.4333221116640807,.45,.493608094,.5):
    ratio=1.5*cap/kappa
    contrast=2*ratio*math.sqrt(1-ratio*ratio)
    amplitude=math.sin(math.pi*contrast/2)
    thresholds.append({"cap_diagnostic":cap,"target_slope_diagnostic":1.5*cap,
                       "necessary_latent_band_radius_diagnostic":amplitude})
print(json.dumps({"status":"PASS exact lazy-walk distributions and moments",
                  "sharp_half_band_floor_diagnostic":floor,
                  "exact_coarse_floor_squared_lower":"203/330 > 9/16",
                  "matching_pair_responses":records,
                  "general_band_thresholds_diagnostic_only":thresholds},indent=2))
