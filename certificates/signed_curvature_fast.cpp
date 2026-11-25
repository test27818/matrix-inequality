// Integer alternative implementation for the finite signed-curvature kernel.
// Start with run_signed_curvature.py, which supplies the rigorous cutoff atom.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <set>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using I64=int64_t; using W=__int128_t;
constexpr int bits=40;constexpr I64 P=I64(1)<<bits;
I64 ck(W x){if(x<INT64_MIN||x>INT64_MAX)throw std::overflow_error("dyadic integer");return I64(x);}
I64 down(W a,I64 b){if(b<=0)throw std::runtime_error("denominator");if(b==P) return ck(a<0?-((-a+P-1)>>bits):a>>bits);return ck(a/b-(a%b<0));}
I64 up(W a,I64 b){return ck(-W(down(-a,b)));}
struct I{I64 l,h;I():l(0),h(0){}I(I64 n):l(ck(W(n)*P)),h(l){}
 static I raw(I64 a,I64 b){if(a>b)throw std::runtime_error("empty");I x;x.l=a;x.h=b;return x;}
};I operator+(I a,I b){return I::raw(ck(W(a.l)+b.l),ck(W(a.h)+b.h));}
I operator-(I a){return I::raw(ck(-W(a.h)),ck(-W(a.l)));}
I operator-(I a,I b){return a+(-b);}
I operator*(I a,I b){W ll=W(a.l)*b.l,lh=W(a.l)*b.h,hl=W(a.h)*b.l,hh=W(a.h)*b.h;
return I::raw(down(std::min({ll,lh,hl,hh}),P),up(std::max({ll,lh,hl,hh}),P));}
I operator/(I a,I b){if(b.l<=0)throw std::runtime_error("interval division");return I::raw(down(W(a.l)*P,a.l<0?b.l:b.h),up(W(a.h)*P,a.h<0?b.h:b.l));}
I operator/(I a,I64 n){return I::raw(down(a.l,n),up(a.h,n));}
I hull(I x,I y){return I::raw(std::min(x.l,y.l),std::max(x.h,y.h));}
I narrow(I x,I y){return I::raw(std::max(x.l,y.l),std::min(x.h,y.h));}
I sq(I a){if(a.h<0)a=-a;if(a.l<0)return I::raw(0,up(W(std::max(-a.l,a.h))*std::max(-a.l,a.h),P));return a*a;}
I64 iroot(I64 z){if(z<0)throw std::runtime_error("negative root");W n=W(z)*P;I64 r=I64(std::sqrt((long double)n));
while(W(r)*r>n)--r;while(W(r+1)*(r+1)<=n)++r;return r;}
I root(I a){return I::raw(iroot(a.l),iroot(a.h)+1);}
using F4=std::array<I,4>;using F3=std::array<I,3>;using F7=std::array<I,7>;
std::unordered_map<I64,F7> cache;
F7 const& series(I64 x){
 auto it=cache.find(x);if(it!=cache.end())return it->second;
 if(x<0||x>8*P)throw std::runtime_error("Taylor range");
 I z=I::raw(x,x),t(1);F7 out{};
 // Each coefficient at degree k bounded by 1/k!. For x<=8 the tail beyond
 // 80 is bounded by 2*8^81/81! < 2^-40; recurrence terms rounded out.
 for(int k=0;k<=80;++k){if(k%2==0){
  out[0]=out[0]+t;out[1]=out[1]+t/(k+1);
  auto c=t/((k+1)*(k+2));out[2]=out[2]+c;out[3]=out[3]+I(k)*c;
  out[4]=out[4]+t/((k+1)*(k+2)*(k+3));
  out[6]=out[6]+t/((k+3)*(k+4)*(k+5));
 }else out[5]=out[5]+t/((k+2)*(k+3)*(k+4));
 t=(t*z)/(k+1);}
 for(auto& a:out)a.h++;
 return cache.emplace(x,out).first->second;
}
F7 both(I x){F7 a=series(x.l),b=series(x.h);for(int j=0;j<7;++j)a[j].h=b[j].h;return a;}
const I B=I::raw(0,P/256);
I corr(I v){return I(1)/(I(1)-B*sq(v)/20);}
F3 phi(I q){auto z=sq(q);I f=I(2)/(I(1)-B),fz(0),fzz(0);const int N=48;
for(int j=2;j<=N;++j){auto inv=I(1)/(I(j*j)-z),g=I(j*j)/(I(j*j)-B);
 f=f+g*I(2)*(I(1)-z)*inv/j;fz=fz-g*I(2*(j*j-1))*sq(inv)/j;fzz=fzz-g*I(4*(j*j-1))*sq(inv)*inv/j;}
 auto g=I((N+1)*(N+1))/(I((N+1)*(N+1))-B);auto t=g/(N*(N+1));f.h+=t.h;fz.l-=t.h;
 fzz.l-=(g*I(2)/(N*(N+1)*((N+1)*(N+1)-1))).h;
 return {f,I(2)*q*fz,I(2)*fz+I(4)*z*fzz};}
F4 inputs(I q,I v){F7 a=both(v),t=both(q*v);auto iv=I(1)/a[2],u=v*iv*corr(v);
return {u*t[4],u*v*t[5],u*sq(v)*t[6],(t[0]-a[3]*iv*t[1])*iv};}
F3 rootJet(F3 t){auto r=root(t[0]);return {r,t[1]/(I(2)*r),t[2]/(I(2)*r)-sq(t[1])/(I(4)*r*r*r)};}
using Pair=std::pair<int,int>;using Job=std::array<int,3>;
using Table=std::vector<std::pair<I64,I64>>;
std::map<std::tuple<int,int>,Table> wtcache;
I trig(I u, bool sinc){
 if(u.l<0||u.h>256*P)throw std::runtime_error("trig range");
 if(u.h>P){auto c=trig(u/2,false);return sinc?trig(u/2,true)*c:I(2)*sq(c)-I(1);}
 int j=sinc?1:0;I term(1),v(1),U=sq(u);
 for(int k=1;k<=10;++k){term=(-term*U)/((2*k-1+j)*(2*k+j));v=v+term;}
 auto tail=(-term*U)/((21+j)*(22+j));return hull(v,v+tail);
}
Table const& weights(int n,int k){
 auto key=std::make_tuple(n,k);auto ext=wtcache.find(key);if(ext!=wtcache.end())return ext->second;
 Table tab;for(int j=0;j<n;++j){
  if(k==-1){auto a=I(3*j*j+3*j+1)/(3*n*n),right=I(j+1)/n;
      tab.emplace_back((a*(I(1)/2-sq(right)/24)).l,(a/2).h);
  }else{
    // frequency [k/2,(k+1)/2], centre (2k+1)/4 and halfwidth 1/4
    auto t=I(2*j+1)/(2*n),mu=I(2*k+1)/4;
    auto avg=I(1)-trig(t*mu,false)*trig(mu/(2*n),true);auto er=t/4;
    tab.emplace_back(std::max<I64>(0,avg.l-er.h),std::min<I64>(2*P,avg.h+er.h));
  }
 }
 return wtcache.emplace(key,std::move(tab)).first->second;
}
W wadd(W a,W b){W c;if(__builtin_add_overflow(a,b,&c))throw std::overflow_error("signed budget");return c;}
int freq(std::vector<I64> const&a,I64 atom) {
 int n=a.size();if(n>4096||atom>0)throw std::runtime_error("freq input");
 I64 mv=*std::min_element(a.begin(),a.end());if(mv==INT64_MIN)return -1;
 if(mv>=0&&atom==0)return 0;
 W mass=W(2)*n*atom,tv=std::abs(a.front()); // sin(0)=0, q-array reversed
 for(int i=0;i<n;++i){mass=wadd(mass,a[i]);if(i)tv=wadd(tv,std::abs(ck(W(a[i])-a[i-1])));}
 if(mass<=0||tv*n>128*mass)return -1;
 int omg=std::max(1,int((tv*n+mass-1)/mass));
 auto chk=[&](int k){W sum=W(n)*atom*(k==-1?P/2:2*P);auto const&w=weights(n,k);
   for(int i=0;i<n;++i){auto [l,h]=w[n-1-i];sum=wadd(sum,W(a[i])*(a[i]>=0?l:h));}return sum>0;};
 if(!chk(-1))return -1;
 for(int k=2;k<2*omg;++k)if(!chk(k))return -1;
 return omg;}
int main(int argc,char**argv){
 if(argc!=7){std::cerr<<"edge cap radial_bins_per_unit phase_bins levels atom_lower_numerator\n";return 2;}
 int edge=atoi(argv[1]),cap=atoi(argv[2]),bpu=atoi(argv[3]),nq=atoi(argv[4]),levels=atoi(argv[5]);
 if(!(edge>=1&&edge<=8&&cap>=1&&cap<=edge&&bpu>0&&nq>0&&levels>0))return 2;
 std::vector<Pair> jobs;I64 atom=std::stoll(argv[6]);
 if(atom>0)throw std::runtime_error("cutoff allowance sign");
 for(int lev=0;lev<levels;++lev,bpu*=2,nq*=2){
  int n=edge*bpu; if(P%bpu||P%nq||n>8192||nq>4096)throw std::runtime_error("grid");I64 dv=P/bpu,dq=P/nq;
  if(!lev)for(int i=0;i<cap*bpu;++i)for(int m=i;m<n;++m)jobs.push_back({i,m});
  std::vector<std::vector<I64>> profiles(jobs.size());std::vector<Pair> failed;I64 visited=0,localMin=INT64_MAX;int mi=n,mm=n;
  for(int k=0;k<nq;++k){
   // q slices are independent. Clear cached Taylor evaluations periodically
   // so large refinement runs maintain a bounded memory footprint.
   if(cache.size()>250000)cache.clear();
   I q=I::raw(k*dq,(k+1)*dq),z=sq(q),one(1);F3 ph=phi(q);I bq=I::raw(0,std::min(q.h,P/16));if(q.l>0)bq=narrow(bq,B/q);
   std::vector<F4> ff(n);std::vector<F3> pre(n+1),roots(n);std::vector<I> pr(n+1);
   for(int i=0;i<n;++i){I vi=I::raw(i*dv,(i+1)*dv);ff[i]=inputs(q,vi);F3 u;
    pr[i+1]=pr[i]+vi*both(q*vi)[1]/both(vi)[2];
    for(int j=0;j<3;++j){pre[i+1][j]=pre[i][j]+ff[i][j];u[j]=I::raw(pre[i][j].l,pre[i+1][j].h)/bpu;}
    auto T=ph[0]-(one-z)*u[0];
    // The positive lower bound T_s >= j=1 integral contribution:
    // T_s=(1-q^2) J_q(s)/q^3 >= (1-q^2) exp(-s) (s^3+3s^2+6s+6)/3,
    // not used; skip ill-conditioned coarse cells rather than assume positivity.
    if(T.l>0)roots[i]=rootJet({T,ph[1]+I(2)*q*u[0]-(one-z)*u[1],
      ph[2]+I(2)*u[0]+I(4)*q*u[1]-(one-z)*u[2]});
   }
   std::vector<std::vector<F4>> rm{ff};for(int j=1;(1<<j)<=n;++j){
    int len=1<<(j-1);std::vector<F4> nxt(n-2*len+1);for(int i=0;i<(int)nxt.size();++i)
     for(int v=0;v<4;++v)nxt[i][v]=hull(rm.back()[i][v],rm.back()[i+len][v]);rm.push_back(std::move(nxt));}
   auto range=[&](int a,int b){int j=31-__builtin_clz((unsigned)(b-a+1)),h=b+1-(1<<j);F4 x;
      for(int v=0;v<4;++v)x[v]=hull(rm[j][a][v],rm[j][h][v]);return x;};
   auto eval=[&](int i,int m,int idx){++visited;auto S=roots[i][0]+roots[m][0];I val;
    if(roots[i][0].l>0&&roots[m][0].l>0){
     auto r=range(i,m);F3 W{r[0],r[1],r[2]};auto PP=r[3]/8;
     if(m>i+1)for(int j=0;j<3;++j){
       auto lo=pre[m][j].l-pre[i+1][j].l,hi=pre[m][j].h-pre[i+1][j].h;
       W[j]=narrow(W[j],I::raw(down(lo,m-i+1),up(__int128_t(hi)+ff[i][j].h+ff[m][j].h,m-i-1)));
     }
     W[0]=narrow(W[0],I::raw(0,((I(i+m+2)/(6*bpu))*corr(I(edge))).h));
     auto S1=roots[i][1]+roots[m][1],S2=roots[i][2]+roots[m][2];
     auto L=one/(I(2)*sq(S)),l1=I(-2)*S1/S,l2=I(6)*sq(S1/S)-I(2)*S2/S;
     auto lower=[&](I p0,F3 wv){auto [w,v,x]=wv;auto w2=sq(w),Y=L*w2;
      auto Y1=L*(I(2)*w*v+w2*l1),Y2=L*(I(2)*sq(v)+I(2)*w*x+I(4)*w*v*l1+w2*l2);
      return p0+(I(6)*(one+B)-I(20)*z)*Y+(I(6)*(one+B)*q-I(10)*q*z-I(2)*bq)*Y1+((one+B)*z-sq(z)-B)*Y2;};
     if(m-i>=2){
       // Integrate the explicit primitive H(v) = d/dv of v*sinhc(qv)/C(v).
       // P0 is the second divided difference of its first radial primitive.
       int a=(i+m)/2,b=(i+m+3)/2; // lower floor, upper ceil; b<=n
       I U=I::raw(pr[i].l,pr[i+1].h),V=I::raw(pr[m].l,pr[m+1].h),T=I::raw(pr[a].l,pr[b].h);
       I d=I::raw((m-i-1)*dv,(m-i+1)*dv);
       PP=narrow(PP,((U+V)/2-T)/bpu/sq(d));
     }
     val=lower(PP,W);
     if(val.l<=0&&m>i){int M=(bpu>=128&&m-i>bpu)?64:32;F3 av{};I pp(0);
      for(int j=0;j<M;++j){auto box=range(i+(m-i)*j/M,i+((m-i)*(j+1)+M-1)/M);
       for(int v=0;v<3;++v)av[v]=av[v]+box[v];pp=pp+box[3]*I(2*std::min(j,M-1-j)+1);}
      for(int v=0;v<3;++v)av[v]=narrow(W[v],av[v]/M);
      val=lower(narrow(PP,pp/(4*M*M)),av);}
    }
    bool ok=roots[i][0].l>0&&roots[m][0].l>0;
    I64 c=ok?(q*val).l:INT64_MIN;
    if(q.l<P/16)c=std::min<I64>(c,0);profiles[idx].push_back(c);
   };
   for(size_t j=0;j<jobs.size();++j)eval(jobs[j].first,jobs[j].second,int(j));
  }
  int maxfreq=0;
  for(size_t j=0;j<jobs.size();++j){int f=freq(profiles[j],atom);auto [i,m]=jobs[j];
   if(f<0){failed.push_back(jobs[j]);mi=std::min(mi,i);mm=std::min(mm,m);}else maxfreq=std::max(f,maxfreq);}
  std::cout<<"edge="<<edge<<" cap="<<cap<<" bpu="<<bpu<<" nq="<<nq<<" pairs="<<jobs.size()<<" unresolved="<<failed.size()
    <<" max_cut="<<maxfreq<<" certified_first="<<mi<<"/"<<bpu<<" certified_second="<<mm<<"/"<<bpu<<std::endl;
  if(failed.empty()){std::cout<<"PASS all original radial pairs, physical parameters and frequencies"<<std::endl;return 0;}
  jobs.clear();for(auto [i,m]:failed){for(int x=0;x<=1;++x)for(int y=0;y<=1;++y)if(2*i+x<=2*m+y)jobs.push_back({2*i+x,2*m+y});}
 }
 std::cout<<"INCOMPLETE; unresolved radial pairs never counted as passed"<<std::endl;return 1;
}
