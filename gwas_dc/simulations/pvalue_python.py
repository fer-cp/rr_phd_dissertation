from mpmath import mpmathify, appellf1, exp, sqrt, log, mp

mp_exp = exp
mp_sqrt = sqrt
mp_log = log

del log, exp, sqrt

def pvalue_python(xax,l21,v, digits):
  mp.dps = (digits+8) * 4
  xax = mpmathify(xax)
  l21 = mpmathify(l21)
  v = mpmathify(v)
  return float(1-appellf1((v+2)/2, 1/2, 1, 2,xax*(1-l21),xax, maxterms = 10**9)*mp_exp(v/2*mp_log(1-xax))*xax*v/2*mp_sqrt(l21))
