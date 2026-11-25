"""Finite-order original J-midpoint kernel in a compact radius square, 0<=kappa<=1/16.
Integer interval certifier adapted from verify_continuum_shallow.py. No floating point.
Arguments as for the continuum file. See finite-shallow-square.md for derivation.
"""
from verify_continuum_shallow import (I,P,ival,frac,narrow,series,rt,ceildiv)
from fractions import Fraction
import sys
B=I(0,frac(1,256).h)
def corr(v):return 1/(1-B*v.sq()/20)
def phi_b(q):
    z=q*q;n=48
    f=2/(1-B);fz=ival(0);fzz=ival(0)
    for j in range(2,n+1):
        d=j*j-z; g=ival(j*j)/(j*j-B)
        f=f+g*(2*(1-z))/(j*d)
        fz=fz-g*ival(2*(j*j-1))/(j*d*d)
        fzz=fzz-g*ival(4*(j*j-1))/(j*d*d*d)
    g=ival((n+1)**2)/((n+1)**2-B)
    f=f+I(0,(g*frac(1,n*(n+1))).h)
    fz=fz+I(-(g*frac(1,n*(n+1))).h,0)
    fzz=fzz+I(-(g*frac(2,n*(n+1)*((n+1)**2-1))).h,0)
    return (f,2*q*fz,2*fz+4*z*fzz)
def fvals_b(q,v):
    t=q*v
    C=series('C',v);S=series('S',v);fac=corr(v)
    H=(series('cosh',t)-(S/C)*series('sinhc',t))/C
    return (v*series('E',t)/C*fac,v.sq()*series('E1',t)/C*fac,
            (v.sq()*v)*series('E2',t)/C*fac,H)
def endpoint(num,den,bpu):
    assert num*bpu%den==0 and P%bpu==0
    n=num*bpu//den;dv=frac(1,bpu);q=I(0,frac(1,16).h);Kmax=phi_b(q)[0].h
    rows=[]
    for i in range(n):
        v=I(i*dv.l,(i+1)*dv.h);t=q*v
        C=series('C',v);S=series('S',v);ec=series('E',t)/C*corr(v)
        # v*partial_v of every summand of Psi is nonnegative. These are
        # lower bounds for f, d f/dv (do not discard a negative coefficient).
        rows.append(((v*ec).l,((1-S/C)*ec).l))
    # Range-minimum tables over the closed radial cells.
    tabs=[rows]
    while 2**len(tabs)<=n:
        last=tabs[-1];jump=2**(len(tabs)-1)
        tabs.append([tuple(min(a,b) for a,b in zip(last[i],last[i+jump])) for i in range(len(last)-jump)])
    def lo(i,m):
        j=(m-i+1).bit_length()-1
        return tuple(min(a,b) for a,b in zip(tabs[j][i],tabs[j][m+1-2**j]))
    minval=10**100;bad=0;M=16
    for i in range(n):
        for m in range(i,n):
            w=t=ival(0)
            for j in range(M):
                a=i+(m-i)*j//M;c=i+ceildiv((m-i)*(j+1),M)
                wl,tl=lo(a,c)
                w=w+I(wl)/M
                t=t+I(tl)*frac(2*min(j,M-1-j)+1,4*M*M)
            # Use only the lower endpoint; w has nonnegative lower endpoint.
            assert w.l>=0 and Kmax>0
            v=t+(1-B)*w.sq()/(8*I(Kmax))
            minval=min(minval,v.l)
            bad+=int(v.l<0)
    print('cutoff slope coefficient boxes='+str(n*(n+1)//2),
          'failed='+str(bad),'lower='+str(minval)+'/'+str(P),flush=True)
    return bad==0
def run(num,bpu,nbq,den=1, targets=None,cap_num=None,cap_den=1):
    assert 0<num<=4*den and num*bpu%den==0 and P%bpu==0 and P%nbq==0
    n=num*bpu//den;dv=frac(1,bpu);dq=frac(1,nbq)
    best=10**100;which=None;failed=[]; passed_min=10**100; visited=0
    active=None if targets is None else {v[2] for v in targets}
    maxi=n if cap_num is None else min(n,cap_num*bpu//cap_den)
    assert cap_num is None or cap_num*bpu%cap_den==0
    byk={}
    if targets is not None:
        for i,m,k in targets:byk.setdefault(k,[]).append((i,m))
    for k in range(nbq):
        if active is not None and k not in active: continue
        q=I(k*dq.l,(k+1)*dq.h); z=q*q
        ph0,ph1,ph2=phi_b(q)
        # b/q extends through q=b=0 on the physical set b<=q^2.
        bh=min(q.h,frac(1,16).h); bq=I(0,bh)
        if q.l>0:bq=narrow(bq,B/q)
        f=[fvals_b(q,I(i*dv.l,(i+1)*dv.h)) for i in range(n)]
        su=[(ival(0),ival(0),ival(0))]
        for row in f:su.append(tuple(su[-1][j]+dv*row[j] for j in range(3)))
        roots=[]
        for i in range(n):
            U,V,X=(I(su[i][j].l,su[i+1][j].h) for j in range(3))
            roots.append(rt((ph0-(1-z)*U,ph1+2*q*U-(1-z)*V,ph2+2*U+4*q*V-(1-z)*X)))
        tabs=[f]
        while 2**len(tabs)<=n:
            last=tabs[-1];jump=2**(len(tabs)-1)
            tabs.append([tuple(I(min(last[i][j].l,last[i+jump][j].l),max(last[i][j].h,last[i+jump][j].h)) for j in range(4)) for i in range(len(last)-jump)])
        def rmq(lo,hi):
            j=(hi-lo+1).bit_length()-1;left=tabs[j][lo];right=tabs[j][hi+1-(1<<j)]
            return [I(min(x.l,y.l),max(x.h,y.h)) for x,y in zip(left,right)]
        pre=[[0,0,0,0,0,0]]
        for row in f:pre.append([pre[-1][j]+(row[j//2].l if j%2==0 else row[j//2].h) for j in range(6)])
        pairs=((i,m) for i in range(maxi) for m in range(i,n)) if targets is None else sorted(byk[k])
        for i,m in pairs:
                fs=rmq(i,m)
                visited+=1
                P0=fs[3]/8
                W=fs[:3]
                if m>=i+2:
                    for j in range(3):
                        # All intervening full radial cells, plus subsets of the two end cells
                        l=pre[m][2*j]-pre[i+1][2*j];h=pre[m][2*j+1]-pre[i+1][2*j+1]
                        W[j]=narrow(W[j],I(l//(m-i+1),ceildiv(h+f[i][j].h+f[m][j].h,m-i-1)))
                # F(q,v)<=v/3 implies W<=(A+B)/6
                W[0]=narrow(W[0],I(0,(frac(i+m+2,6*bpu)*corr(frac(num,den))).h))
                S,S1,S2=(roots[i][j]+roots[m][j] for j in range(3))
                L=1/(2*S*S); l1=-2*S1/S; l2=6*(S1/S).sq()-2*S2/S
                def lower(p0,weights):
                    w,v,x=weights
                    Y=L*w.sq(); Y1=L*(2*w*v+w.sq()*l1);Y2=L*(2*v.sq()+2*w*x+4*w*v*l1+w.sq()*l2)
                    return p0+(6*(1+B)-20*z)*Y+(6*(1+B)*q-10*q*z-2*bq)*Y1+((1+B)*z-z*z-B)*Y2
                val=lower(P0,W)
                if val.l<=0 and m>i:
                    # W is an ordinary uniform average under v=A+(B-A)y.
                    # P0 is the average of H against min(y,1-y)/2.
                    M=16; wt=[ival(0),ival(0),ival(0)];pp=ival(0)
                    for slot in range(M):
                        lo=i+(m-i)*slot//M
                        hi=i+ceildiv((m-i)*(slot+1),M)
                        box=rmq(lo,hi)
                        for j in range(3):wt[j]=wt[j]+box[j]/M
                        hh=box[3]
                        pp=pp+hh*frac(2*min(slot,M-1-slot)+1,4*M*M)
                    val=lower(narrow(P0,pp),[narrow(x,y) for x,y in zip(W,wt)])
                if val.l<best: best,which=val.l,(i,m,k)
                if val.l<=0:failed.append((i,m,k))
                else:passed_min=min(val.l,passed_min)
    print('edge='+str(Fraction(num,den)),'radial_bins='+str(n),'phase_bins='+str(nbq),
          'visited='+str(visited),'unresolved='+str(len(failed)),
          'certified_lower='+str(best)+'/'+str(P),'worst_cell='+str(which),flush=True)
    if failed:
        print('certified_strips_in_original_domain first_depth<='+str(Fraction(min(w[0] for w in failed),bpu))+' OR second_depth<='+str(Fraction(min(w[1] for w in failed),bpu)),flush=True)
    return failed,passed_min
def check(num=2,bpu=16,nbq=32,den=1,levels=4,cap_num=None,cap_den=1):
    print('finite kappa in [0,1/16], b in [0,1/256]; original-kernel midpoint curvature',flush=True)
    if not endpoint(num,den,max(32,bpu)):
        print('endpoint not proved');return False
    if cap_num is not None:print('first radial depth capped at',str(Fraction(cap_num,cap_den)),flush=True)
    targets=None; minbound=10**100
    for _ in range(levels):
        f,p=run(num,bpu,nbq,den,targets,cap_num,cap_den);minbound=min(minbound,p)
        if not f:
            print('PASS: all original boxes certified; minimum across accepted boxes',str(minbound)+'/'+str(P),flush=True);return True
        targets={(2*i+x,2*m+y,2*k+z) for i,m,k in f for x in [0,1] for y in [0,1] for z in [0,1] if 2*i+x<=2*m+y}
        bpu*=2;nbq*=2
    print('INCOMPLETE; no conclusion on unresolved boxes');return False
if __name__=='__main__':
    args=list(map(int,sys.argv[1:]));sys.exit(0 if check(*args) else 1)
