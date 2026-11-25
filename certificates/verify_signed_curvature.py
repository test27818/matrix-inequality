"""All-phase proof on finite original J-kernel by signed-curvature cosine envelope.
No floating-point acceptance, and never assumes that D'' is positive.
See signed-curvature.md; compares to the same finite Hermite minorant.
"""
import sys
from fractions import Fraction
from functools import lru_cache
from math import factorial
from verify_continuum_shallow import I,P,ival,frac,narrow,rt,ceildiv,cp
import verify_finite_shallow as fs
from verify_finite_shallow import B,phi_b,corr
MAXEDGE=8
# On [0, MAXEDGE] each named Taylor coefficient is <=1/n!, so
# successive scalar tail terms have ratio <=8/82<1/2 after degree 80.
assert 2*MAXEDGE**81*P<factorial(81)
NAMES=('cosh','sinhc','C','S','E','E1','E2')
RATIOS=[tuple((u*factorial(n),factorial(k)) for u,k in (cp(f,n) for f in NAMES)) for n in range(81)]
@lru_cache(maxsize=120000)
def ser_tuple(x):
    assert 0<=x<=MAXEDGE*P
    tl=th=P;rlo=[0]*7;rhi=[0]*7
    for n,coefs in enumerate(RATIOS):
        for j,(a,d) in enumerate(coefs):
            if a:
                rlo[j]+=tl*a//d;rhi[j]+=ceildiv(th*a,d)
        tl=tl*x//((n+1)*P);th=ceildiv(th*x,(n+1)*P)
    return tuple(zip(rlo,[z+1 for z in rhi]))
def series(f,x):
    j=NAMES.index(f)
    return I(ser_tuple(x.l)[j][0],ser_tuple(x.h)[j][1])
fs.series=series
fvals_b=fs.fvals_b
def endpoint(num,den,bpu):
    assert num*bpu%den==0 and P%bpu==0
    n=num*bpu//den;dv=frac(1,bpu);q=I(0,frac(1,16).h);Kmax=phi_b(q)[0].h
    rows=[];kup=[];u=ival(0)
    for i in range(n):
        v=I(i*dv.l,(i+1)*dv.h);t=q*v
        C=series('C',v);S=series('S',v);ec=series('E',t)/C*corr(v)
        # v*partial_v of every summand of Psi is nonnegative. These are
        # lower bounds for f, d f/dv (do not discard a negative coefficient).
        rows.append(((v*ec).l,((1-S/C)*ec).l))
        kup.append((I(Kmax)-(1-B)*u).h)
        u=u+dv*I((v*ec).l)
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
            S=I(kup[i]).sqrt()+I(kup[m]).sqrt()
            v=t+(1-B)*w.sq()/(2*S.sq())
            minval=min(minval,v.l)
            bad+=int(v.l<0)
    print('cutoff slope coefficient boxes='+str(n*(n+1)//2),
          'failed='+str(bad),'lower='+str(minval)+'/'+str(P),flush=True)
    return min(0,(2*B*I(minval)).l)
def profiles(num,bpu,nbq,den=1,targets=None):
    assert 0<num<=MAXEDGE*den and num*bpu%den==0 and P%bpu==0 and P%nbq==0
    n=num*bpu//den;dv=frac(1,bpu);dq=frac(1,nbq)
    pairs=list((i,m) for i in range(n) for m in range(i,n)) if targets is None else sorted(targets)
    prof={pair:[] for pair in pairs}
    for k in range(nbq):
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
            K=ph0-(1-z)*U
            roots.append(rt((K,ph1+2*q*U-(1-z)*V,ph2+2*U+4*q*V-(1-z)*X)) if K.l>0 else None)
        tabs=[f]
        while 2**len(tabs)<=n:
            last=tabs[-1];jump=2**(len(tabs)-1)
            tabs.append([tuple(I(min(last[i][j].l,last[i+jump][j].l),max(last[i][j].h,last[i+jump][j].h)) for j in range(4)) for i in range(len(last)-jump)])
        def rmq(lo,hi):
            j=(hi-lo+1).bit_length()-1;left=tabs[j][lo];right=tabs[j][hi+1-(1<<j)]
            return [I(min(x.l,y.l),max(x.h,y.h)) for x,y in zip(left,right)]
        pre=[[0,0,0,0,0,0]]
        for row in f:pre.append([pre[-1][j]+(row[j//2].l if j%2==0 else row[j//2].h) for j in range(6)])
        for i,m in pairs:
                fs=rmq(i,m)
                P0=fs[3]/8
                W=fs[:3]
                if m>=i+2:
                    for j in range(3):
                        # All intervening full radial cells, plus subsets of the two end cells
                        l=pre[m][2*j]-pre[i+1][2*j];h=pre[m][2*j+1]-pre[i+1][2*j+1]
                        W[j]=narrow(W[j],I(l//(m-i+1),ceildiv(h+f[i][j].h+f[m][j].h,m-i-1)))
                # F(q,v)<=v/3 implies W<=(A+B)/6
                W[0]=narrow(W[0],I(0,(frac(i+m+2,6*bpu)*corr(frac(num,den))).h))
                if roots[i] is None or roots[m] is None:
                    prof[(i,m)].append(None);continue
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
                coeff=(q*val).l
                # For q<kappa, D is extended by zero. We may still use a
                # negative lower bound; a positive one is truncated to zero.
                if q.l<frac(1,16).h:coeff=min(0,coeff)
                prof[(i,m)].append(coeff)
    return prof
@lru_cache(maxsize=150000)
def trig_s(q, kind):
    # Trig denoting cos(q) or sinc(q). Halve until |q|<=1; then signed
    # alternating Taylor intervals. sinc(2t)=sinc(t)*cos(t).
    assert isinstance(q,Fraction) and q>=0
    if q>1:
        c=trig_s(q/2,'c')
        return 2*c.sq()-1 if kind=='c' else trig_s(q/2,'s')*c
    j=0 if kind=='c' else 1
    sq=q*q
    term=Fraction(1); val=term
    for i in range(1,13):
        term=-term*sq/((2*i-1+j)*(2*i+j));val+=term
    tail=-term*sq/((25+j)*(26+j))
    def fi(t):return frac(t.numerator,t.denominator)
    return I(fi(val+tail).l,fi(val).h)
@lru_cache(maxsize=30000)
def weights(n,w0,w1,small=False):
    # On cell p in [i/n,(i+1)/n] enclose the cell average of
    # 1-cos(w p), for every w in [w0,w1]. In the small band divide by w^2.
    out=[]
    for i in range(n):
        if small:
            # For 0<=w<=1, 1/2-(wp)^2/24 <= (1-cos(wp))/(wp)^2 <=1/2.
            a=Fraction(3*i*i+3*i+1,3*n*n)
            hi=a/2;lo=a*(Fraction(1,2)-Fraction(1,24)*Fraction(i+1,n)**2)
            out.append((frac(lo.numerator,lo.denominator).l,frac(hi.numerator,hi.denominator).h));continue
        w=(w0+w1)/2;t=Fraction(2*i+1,2*n)
        c=trig_s(w*t,'c');sc=trig_s(w/(2*n),'s')
        v=1-c*sc; eps=(w1-w0)*t/2
        e=frac(eps.numerator,eps.denominator).h
        out.append((max(0,v.l-e),min(2*P,v.h+e)))
    return out
def frequency_margin(reversed_coeffs,atom,step=Fraction(1,2),cap=128):
    if any(c is None for c in reversed_coeffs):return None
    # Coeffs returned in q order; cosine variable p=1-q.
    a=reversed_coeffs[::-1];n=len(a)
    # All weights positive; this fast path is legitimate even if min(a)=0 at cutoff.
    if min(a)>=0 and atom==0:return (0,0)
    mass=sum(a)+2*n*atom
    if mass<=0:return None
    # Piecewise-constant curvature lower profile has TV = sum of jumps.
    # Boundary sine at p=0 vanishes, so no atom |a[0]| is needed.
    tv=abs(a[-1])+sum(abs(a[i+1]-a[i]) for i in range(n-1))
    omega=max(1,ceildiv(n*tv,mass))
    if omega>cap:return None
    def lower(weights):
        return sum(c*(l if c>=0 else h) for c,(l,h) in zip(a,weights))
    low=lower(weights(n,Fraction(0),Fraction(1),True))+n*atom*(P//2)
    if low<=0:return None
    w=Fraction(1)
    while w<omega:
        h=min(Fraction(omega),w+step)
        v=lower(weights(n,w,h))+2*n*atom*P
        if v<=0:return None
        w=h
    return (omega,low)
def check(edge=4,bpu=8,nq=16,levels=5):
    assert edge in range(1,MAXEDGE+1) and P%bpu==0 and P%nq==0
    atom=endpoint(edge,1,max(32,bpu));print('negative cutoff allowance='+str(atom)+'/'+str(P),flush=True)
    targets=None
    for level in range(levels):
        prof=profiles(edge,bpu,nq,targets=targets)
        bad=[];seenmax=0
        for ij,coeff in prof.items():
            proof=frequency_margin(coeff,atom)
            if proof is None: bad.append(ij)
            else:seenmax=max(seenmax,proof[0])
        print(f'edge={edge} bpu={bpu} phase={nq} pairs={len(prof)} unresolved={len(bad)} max_high_frequency_cut={seenmax}',flush=True)
        if not bad:
            print('PASS: finite kappa signed-curvature all-phase transform certified',flush=True);return True
        print('first_unresolved',bad[:7],flush=True)
        targets={(2*i+x,2*m+y) for i,m in bad for x in (0,1) for y in (0,1) if 2*i+x<=2*m+y}
        bpu*=2;nq*=2
    print('INCOMPLETE');return False
if __name__=='__main__':sys.exit(0 if check(*map(int,sys.argv[1:])) else 1)
