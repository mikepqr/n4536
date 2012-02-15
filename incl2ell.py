import math

def incl2ell(incl, edgeon = 0.2):
    inclrad = incl * math.pi / 180.
    return 1 - math.sqrt( (1 - edgeon**2) * math.cos(inclrad)**2 + edgeon**2)

def ell2incl(ell, edgeon = 0.2):
    q = 1 - ell
    inclrad = math.acos( math.sqrt((q**2 - edgeon**2)/(1 - edgeon**2)) )
    return inclrad * 180. / math.pi
