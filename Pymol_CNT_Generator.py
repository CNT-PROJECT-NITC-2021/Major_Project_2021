
##NIT CALICUT MAJOR PROJECT 2021
import cmd
from math import cos,sin,pi,ceil,floor, acos
from chempy import models, cpv
 
class cell:
  def _init_(self,N,M):
    l = 1.42
    a = 2.*l*cos(pi*30./180.)
    self.a1 = [a, 0., 0.]
    self.a2 = [a*cos(pi*60./180.), a*sin(pi*60./180.), 0.]
    self.crds = [[0., 0., 0.],
                 [l*cos(pi*30./180.), l*sin(pi*30./180.), 0]]
 
    self.Ch = cpv.add(cpv.scale(self.a1, N), cpv.scale(self.a2, M))
    self.lenCh = cpv.length(self.Ch)
    self.Ch = cpv.normalize(self.Ch)
    self.T = cpv.normalize(cpv.cross_product(self.Ch, cpv.cross_product(self.a1, self.a2)))
    self.radius = self.lenCh/math.pi/2.
 
  def get_crds(self, n=0, m=0):
    d = cpv.add(cpv.scale(self.a1,n),cpv.scale(self.a2,m))
    crds = []
    for crd in self.crds:
      v = cpv.add(crd, d)
      ang = 2. * math.pi * cpv.dot_product(v,self.Ch) / self.lenCh
      r = cpv.dot_product(v, self.T)
      crds.append([self.radius*math.cos(ang), self.radius*math.sin(ang),r])
    return crds
 
def new_at(crd, nm, sym, typ, charge=0.0):
  at = chempy.Atom()
  at.charge = charge
  at.name = nm
  at.symbol = sym
  at.type = typ
  at.coord = crd
  at.hetatm = False
  at.resn = 'CNT'
  at.resi = '1'
  at.resi_number = 1
  at.bonds = []
  return at
 
def ntgen(obj, sN,sM,sL, save = None):
  """
  Usage: ntgen obj_name, n, m, l
  """
  bdist = 2.0
  rCH = 1.1
  N,M,L = int(sN), int(sM), int(sL)
  if (M > N):
    _N = N
    N = M
    M = _N
  rN,rM,rL = float(N), float(M), float(L)
  C = cell(N,M)
## Generate atoms
  nt = models.Indexed()
#  for j in range(L):
#    for i in range(N):
#      n,m = i + (j*M)/N, (i*M)/N -j
  iat = 0
  for n in range(0,N+M*L/N):
    if M > 0:
      m_start =  max(int(ceil(-rL*(1 + rM*2/rN*2) + rM*n/rN)), int(floor(-rN*n/rM)) )
      m_stop = min(int(floor(rM*n/rN)), int(floor(rM+rN**2/rM-rN*n/rM)))
    else:
      m_start, m_stop = -L, 0
    for m in range(m_start, m_stop):
      for crd in C.get_crds(n,m):
        nt.add_atom(new_at(crd, "C%d" % (iat+1), 'C', 'opls_145'))
        iat += 1
## Add bonds
  for i in range(1,len(nt.atom)):
    for j in range(i):
      d = cpv.distance(nt.atom[i].coord, nt.atom[j].coord)
      if d <= bdist:
        b = chempy.Bond()
        b.index = [j,i]
        nt.add_bond(b)
        nt.atom[i].bonds.append(j)
        nt.atom[j].bonds.append(i)
        if d <= 0.5:
          print "WARNING: Atoms #%d and #%d are too close" % (i,j)
## Adding protons
  nat = iat
  for iat in range(len(nt.atom)):
    if len(nt.atom[iat].bonds) == 2:
      r0 = nt.atom[iat].coord
      r1 = cpv.sub(nt.atom[nt.atom[iat].bonds[0]].coord, r0)
      r2 = cpv.sub(nt.atom[nt.atom[iat].bonds[1]].coord, r0)
      rh = cpv.add(r0, cpv.scale(cpv.normalize(cpv.add(r1,r2)), -rCH))
      nt.add_atom(new_at(rh, "H%d" % (nat+1), 'H', 'opls_146', 0.115))
      nt.atom[iat].charge = -0.115
      b = chempy.Bond()
      b.index = [iat,nat]
      nt.add_bond(b)
      nt.atom[iat].bonds.append(nat)
      nt.atom[nat].bonds.append(iat)
      nat += 1
    if len(nt.atom[iat].bonds) < 2:
      print "WARNING: Lone carbon detected, id #%d" % iat
# Center system
 
#  for iat in range(len(nt.atom)): print iat, nt.atom[iat].bonds
  C = [0.,0.,0.]
  for iat in range(len(nt.atom)): C = cpv.add(C, nt.atom[iat].coord)
  C = cpv.scale(C, -1./len(nt.atom))
  for iat in range(len(nt.atom)): nt.atom[iat].coord = cpv.add(nt.atom[iat].coord, C)
 
  if save:
## Save GRO
    f = open(save+'.gro', 'w')
    print >>f, "SWNT %d-%d-%d" % (N,M,L)
    print >>f, "%5d" % len(nt.atom)
    for iat in range(len(nt.atom)):
      print >>f, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f" % (1, 'CNT', nt.atom[iat].name, iat+1, nt.atom[iat].coord[0]/10., nt.atom[iat].coord[1]/10., nt.atom[iat].coord[2]/10., 0.0, 0.0, 0.0)
    print >>f, "0.0 0.0 0.0"
    f.close()
 
## Generate topology
    f = open(save+'.itp', 'w')
## Atoms
    print >>f, """
[ moleculetype ]
; name  nrexcl
CNT         3
 
[ atoms ]
;   nr    type   resnr  residu    atom    cgnr  charge"""
    for iat in range(len(nt.atom)):
      print >>f, "%-5d %-10s %-5d %-10s %-10s %-5d %8.3f" % (iat+1, nt.atom[iat].type, 1, 'CNT', nt.atom[iat].name, iat+1, nt.atom[iat].charge)
## Bonds
    print >>f, """
[ bonds ]
; i j"""
    for iat in range(len(nt.atom)):
      for ib in nt.atom[iat].bonds:
        if ib < iat:
          print >>f, "%-5d %-5d 1" % (iat+1, ib+1)
## Angles
    print >>f, """
[ angles ]
; i j k 1"""
    for j in range(len(nt.atom)):
      for i in nt.atom[j].bonds:
        for k in nt.atom[j].bonds:
          if k < i:
            print >>f, "%-5d %-5d %-5d 1" % (i+1, j+1, k+1)
## Dihedrals
    print >>f, """
[ dihedrals ]
; i j k l 2 Theta k"""
    for b in nt.bond:
      j,k = b.index[0], b.index[1]
      if len(nt.atom[j].bonds) >= 2 and len(nt.atom[k].bonds) >= 2:
        for i in nt.atom[j].bonds:
          for l in nt.atom[k].bonds:
            if i != k and l != j:
              rjk = cpv.sub(nt.atom[k].coord, nt.atom[j].coord)
              rji = cpv.sub(nt.atom[i].coord, nt.atom[j].coord)
              rkl = cpv.sub(nt.atom[l].coord, nt.atom[k].coord)
              nijk = cpv.normalize(cpv.cross_product(rji, rjk))
              njkl = cpv.normalize(cpv.cross_product(rkl, rjk))
              cosa = cpv.dot_product(nijk, njkl)
              if cosa > 1.0: cosa = 1.0
              if cosa < -1.0: cosa = -1.0
#              print i,j,k,l,cosa
              a = 180*acos(cosa)/pi
#              print a
              if a > 90.0:
                print >>f, "%-5d %-5d %-5d %-5d 2 %-8.3f 20.0" % (i+1, j+1, k+1, l+1, a)
    f.close()
 
## Load object
#  for iat in range(len(nt.atom)): nt.atom[iat].id = iat+1
#  nt.update_index()
  cmd.delete(obj)
  cmd.load_model(nt, obj)
 
cmd.extend("ntgen", ntgen)
