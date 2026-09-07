"""Outward interval check of the explicit selector-count entropy gap."""
import mpmath as mp

mp.iv.dps = 50
iv = mp.iv
one = iv.mpf(1)
h = lambda z: -z*iv.ln(z)-(one-z)*iv.ln(one-z)
p = iv.mpf(24)/25
eps = one/10**6
v = iv.mpf(100)
xi = h(one/v**2)+iv.ln(one+2/iv.sqrt(eps))/v**2+h(16*eps)+16*eps*iv.ln(2)
a = p*iv.ln(2)-iv.mpf(5151)/6250
d = h(p)+a
assert xi < iv.mpf(2)/1000
assert d > iv.mpf(9)/1000
assert d-xi > iv.mpf(7)/1000
print('SELECTOR_ENTROPY_GAP_INTERVAL_PASS')
print('xi:', xi)
print('h(p)+a:', d)
print('h(p)+a-xi:', d-xi)
