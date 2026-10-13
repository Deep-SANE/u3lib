"""Freeze deterministic integration fixtures before any participant runs."""
import itertools as it
import json,sys,subprocess,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'python'))
import u3lib as u

def rowstate(g):
 return [getattr(g,n)() for n in ('n13','n23','n33','n12','n22','n11')]
def irr(x): return u.irrep(*x)
def rho(a,b,c): return u.Rho(irr(a),irr(b),irr(c))
def channels(a,b,c,t,left):
 out=[]
 # Complete finite bound: lambda+mu of any product output <= sum of input lambda+mu.
 bound=sum(a)+sum(b)
 for x in it.product(range(bound+1),repeat=2):
  if sum(x)>bound: continue
  r=rho(a,b,x); s=rho(x,c,t) if left else rho(c,x,t)
  if r*s:
   assert max(x)<=5,(a,b,c,t,x)
   out += [list(x)+[i,j] for i in range(r) for j in range(s)]
 return out

def generate():
 d={str(i):[] for i in range(1,7)}
 for la,mu in it.product(range(6),repeat=2):
  a=u.irrep(la,mu)
  branches=[[L,u.Kappa(a,L)] for L in range(la+mu+1) if u.Kappa(a,L)]
  assert sum((2*L+1)*k for L,k in branches)==a.dim()
  d['1'].append(dict(labels=[la,mu],dimension=a.dim(),branches=branches,states=[rowstate(g) for g in u.gelfandBasis(a)]))
 u.initialize()
 try:
  for aa,bb,cc in [((0,0),(3,1),(3,1)),((1,0),(0,1),(0,0)),((2,2),(2,2),(2,2)),((3,1),(1,1),(3,1)),((3,2),(2,3),(2,2))]:
   states=[]
   for labels in (aa,bb,cc):
    states.append(sorted(set(tuple(rowstate(g)[:5]) for g in u.gelfandBasis(irr(labels)))))
   count=0
   for a,b,c in it.product(*states):
    eps=lambda x:2*sum(x[:3])-3*sum(x[3:])
    js=[x[3]-x[4] for x in (a,b,c)]
    if eps(a)+eps(b)!=eps(c) or not abs(js[0]-js[1])<=js[2]<=js[0]+js[1] or sum(js)%2: continue
    v=np.asarray(u.wigU3U2db().get(*a,*b,*c))
    assert len(v)==rho(aa,bb,cc) and np.isfinite(v).all()
    d['2'].append(dict(request=[a,b,c],value=v.tolist()))
    count+=1
    if count==12: break
  for aa,bb,cc in [((0,0),(3,1),(3,1)),((2,0),(2,0),(4,0)),((2,2),(1,1),(2,2)),((3,1),(1,1),(3,1))]:
   count=0
   for la,lb,lc in it.product(range(sum(aa)+1),range(sum(bb)+1),range(sum(cc)+1)):
    ks=[u.Kappa(*x,L) for x,L in zip((aa,bb,cc),(la,lb,lc))]
    if not all(ks) or not abs(la-lb)<=lc<=la+lb: continue
    obj=u.wigU3SO3db().get_U3SO3(*aa,la,*bb,lb,*cc,lc)
    shape=ks+[rho(aa,bb,cc)]; v=np.zeros(shape)
    for ix in np.ndindex(*shape): v[ix]=obj(*ix)
    assert np.isfinite(v).all()
    d['3'].append(dict(request=[[*aa,la],[*bb,lb],[*cc,lc]],value=v.tolist()))
    count+=1
    if count==12: break
  for aa,bb,cc,t in [((1,0),(1,0),(1,0),(1,1)),((2,0),(0,2),(2,1),(2,1)),((1,1),(1,1),(1,1),(1,1)),((3,1),(1,1),(1,1),(3,1))]:
   for kind in ('U','Z'):
    left=channels(aa,bb,cc,t,True)
    right=channels(bb,cc,aa,t,False) if kind=='U' else channels(aa,cc,bb,t,True)
    v=np.zeros((len(left),len(right))); A,B,C,T=map(irr,(aa,bb,cc,t))
    for i,(xl,xm,r,s) in enumerate(left):
     x=(xl,xm)
     for j,(yl,ym,p,q) in enumerate(right):
      y=(yl,ym)
      if kind=='U': o=u.u3u6db().get_U3U6(A,B,T,C,irr(x),rho(aa,bb,x),rho(x,cc,t),irr(y),rho(bb,cc,y),rho(aa,y,t))
      else: o=u.u3z6db().get_U3Z6(B,A,T,C,irr(x),rho(aa,bb,x),rho(x,cc,t),irr(y),rho(aa,cc,y),rho(y,bb,t))
      v[i,j]=o(r,s,p,q)
    assert len(left)==len(right)
    assert np.linalg.norm(v@v.T-np.eye(len(left)))<1e-9
    d['4'].append(dict(request=[kind,aa,bb,cc,t],left=left,right=right,value=v.tolist()))
  requests=[[(2,0),(0,0),(0,0),(2,0),(2,0),(2,0),(2,0),(2,0),(4,0)]]
  pair=[(0,0),(1,1),(3,0),(0,3),(2,2)]
  for ab,cd,ac,bd in it.product(pair,repeat=4):
   A=B=C=D=(1,1); T=(1,1)
   rs=[rho(A,B,ab),rho(C,D,cd),rho(ab,cd,T),rho(A,C,ac),rho(B,D,bd),rho(ac,bd,T)]
   if all(rs) and (rs[2]!=rs[5] or ab==cd==ac==bd==(1,1)):
    requests.append([A,B,C,D,ab,cd,ac,bd,T])
    if len(requests)==9: break
  requests += [requests[0],requests[2]]
  for req in requests:
   A,B,C,D,AB,CD,AC,BD,T=map(irr,req)
   rs=[u.Rho(A,B,AB),u.Rho(C,D,CD),u.Rho(AB,CD,T),u.Rho(A,C,AC),u.Rho(B,D,BD),u.Rho(AC,BD,T)]
   o=u.u3u9db().get_U3U9(A,B,C,D,AB,CD,AC,BD,T,rs[3],rs[4],rs[2],rs[0],rs[1],rs[5])
   entries=[]
   for a,b,l,c,e,r in it.product(*(range(x) for x in rs)): entries.append([a,b,l,c,e,r,float(o(c,e,l,a,b,r))])
   d['5'].append(dict(request=req,shape=rs,entries=entries))
  raw_t6 = [
    ('C', [1, 0, 0, 1, 0,  1, 0, 0, 1, 0,  2, 0, 0, 2, 0]),
    ('P', [1, 0, 1,  0, 1, 1,  1, 1, 1]),
    ('C', [2, 1, 0, 2, 1,  2, 1, 0, 2, 0,  3, 2, 1, 3, 2]),
    ('P', [2, 0, 2,  2, 0, 2,  4, 0, 2]),
    ('C', [2, 0, 0, 2, 0,  1, 1, 0, 1, 0,  2, 1, 1, 2, 1]),
    ('P', [2, 2, 2,  1, 1, 2,  2, 2, 2]),
  ]
  for tag, labels in raw_t6:
   req_str = f"{tag} " + " ".join(map(str, labels))
   if tag == 'C':
    row = u.wigU3U2db().get(*labels)
    rhomax = len(row)
    entries = [[rho, float(row[rho])] for rho in range(rhomax)]
    d['6'].append(dict(tag='canon', request=req_str, shape=[rhomax], entries=entries))
   else:
    obj = u.wigU3SO3db().get_U3SO3(*labels)
    shape = [obj.kappaAmax(), obj.kappaBmax(), obj.kappaCmax(), obj.rhomax()]
    entries = []
    for ka, kb, kc, r in it.product(*(range(x) for x in shape)):
     entries.append([ka, kb, kc, r, float(obj(ka, kb, kc, r))])
    d['6'].append(dict(tag='phys', request=req_str, shape=shape, entries=entries))
 finally: u.finalize()
 return d
if __name__=='__main__':
 target=Path(sys.argv[2] if len(sys.argv)>2 else sys.argv[1] if len(sys.argv)>1 else Path(__file__).parent/'fixtures-pilot.json')
 data=generate()
 target.write_text(json.dumps(data,indent=2)+'\n')
 print({k:len(v) for k,v in data.items()}); print('sha256',hashlib.sha256(target.read_bytes()).hexdigest())
