"""Pointwise checks of the exact cavity algebra and proved error bound."""
import math
import random

def lc(x):
    return abs(x)+math.log1p(math.exp(-2*abs(x)))-math.log(2)

rng=random.Random(703)
worst=0.0
for _ in range(100000):
    h=rng.uniform(-15,15)
    v=rng.choice((-1,1))*10**rng.uniform(-2,1)
    m=math.tanh(h+v)
    bracket=2*v*m-2*v*v*(1-m*m)-(lc(h+v)-lc(h-v))
    bound=(8/3)*abs(v)**3
    assert abs(bracket)<=bound+2e-13
    worst=max(worst,abs(bracket)/bound)
print('PASS 100000 pointwise cavity remainder checks; maximum bound ratio',worst)
print('The checks supplement, and do not replace, the uniform Taylor proof.')
