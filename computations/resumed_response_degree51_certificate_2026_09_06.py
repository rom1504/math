"""Exact rational checks for the degree-51 arbitrary-mask direction."""
from fractions import Fraction as F
import json

j,a,t=F(43,100),F(141,1000),F(39,100)
d=j-a
linear=1+j*j-d*d
ft=t**3-2*t*t+linear*t-j*j
fprime=3*t*t-4*t+linear
curvature=1-2*t
maximum_upper=ft+fprime*fprime/(4*curvature)
assert maximum_upper==F(-18965469,80000000000)<0
D=51
root_recip_upper=F(14003,100000)
assert D*root_recip_upper**2>1
margin=a-root_recip_upper
assert margin==F(97,100000)>0
projection_lower=margin/(2**25)
assert projection_lower>F(1,40000000000)
assert projection_lower/5>F(1,200000000000)
# Independent direct verification of the double-factorial bound used here.
odd=lambda k: __import__('functools').reduce(lambda x,y:x*y,range(1,k+1,2),1)
ratio=F(odd(2*D-1),odd(D)**2)
assert ratio<=2**(D-1)
assert F(6,25)*(1-F(6,25))<j*j
print(json.dumps({'J_threshold':str(j),'edge_coefficient_lower':str(a),'uniform_cubic_polynomial_upper':str(maximum_upper),'degree':D,'nonlinear_projection_norm_lower':str(projection_lower),'individual_coefficient_lower':str(projection_lower/5),'verified':True},indent=2))
