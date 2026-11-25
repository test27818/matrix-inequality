// Rigorously outward-rounded binary64 interval proof, original finite kappa.
// Extract exactly the first exponential alias from J. All other aliases
// are positive even-power series with radius >=2 around q=0.
// Compile: c++ -O3 -std=c++20 -ffp-contract=off -fno-fast-math short_zero_residual.cpp -o /tmp/szr
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
using std::max;using std::min;using std::array;using std::vector;
using Key=std::tuple<int,int,int,int,int>;
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53);
double up(double x){
 if(!std::isfinite(x))throw std::runtime_error("range");
 if(x==0)return std::numeric_limits<double>::denorm_min();
 auto b=std::bit_cast<uint64_t>(x);b+=(x>0?1ULL:-1ULL);return std::bit_cast<double>(b);
}
double dn(double x){return -up(-x);}
struct I{double l,h; constexpr I(double x=0):l(x),h(x){}constexpr I(double a,double b):l(a),h(b){}};
I operator+(I a,I b){return {dn(a.l+b.l),up(a.h+b.h)};}
I operator-(I a){return {-a.h,-a.l};}I operator-(I a,I b){return a+(-b);}
I operator*(I a,I b){
 if(a.l>=0){if(b.l>=0)return {dn(a.l*b.l),up(a.h*b.h)};if(b.h<=0)return {dn(a.h*b.l),up(a.l*b.h)};return {dn(a.h*b.l),up(a.h*b.h)};}
 if(a.h<=0)return -((-a)*b);
 if(b.l>=0)return {dn(a.l*b.h),up(a.h*b.h)};if(b.h<=0)return {dn(a.h*b.l),up(a.l*b.l)};
 return {dn(min(a.h*b.l,a.l*b.h)),up(max(a.l*b.l,a.h*b.h))};
}
I operator/(I a,I b){if(b.l<=0)throw std::domain_error("nonpositive denominator");return a*I(dn(1./b.h),up(1./b.l));}
I sq(I a){if(a.l>=0)return {max(0.,dn(a.l*a.l)),up(a.h*a.h)};if(a.h<=0)return sq(-a);
return {0,up(max(a.h*a.h,a.l*a.l))};}
I sqrtI(I a){if(a.l<=0)throw std::domain_error("sqrt lower <=0");
 // proposal never trusted: round outward then test the dyadic square enclosures
 double l=std::sqrt(a.l),h=std::sqrt(a.h);
 for(int i=0;i<8 && up(l*l)>a.l;i++)l=dn(l);
 for(int i=0;i<8 && dn(h*h)<a.h;i++)h=up(h);
 if(up(l*l)>a.l||dn(h*h)<a.h)throw std::runtime_error("root validation");
 return {l,h};}
I rawexp(double v){
 if(v==0)return I(1);
 if(v<-100){ // in coefficient sums exp(-jt) times long partial sums; application avoids this branch for large jt
  int n=1;double t=v;while(t<-64){t*=.5;n*=2;}auto w=rawexp(t);while(n>1){w=sq(w);n/=2;}return w;
 }
 if(v<0)return I(1)/rawexp(-v);
 int n=0;while(v>.5){v*=.5;n++;}I x(v),z(1),res(1);
 for(int k=1;k<=18;++k){z=z*x/I(k);res=res+z;}
 auto rem=I(2)*z*x/I(19);res.h=up(res.h+rem.h);while(n-->0)res=sq(res);return res;
}
I Exp(I x){return {rawexp(x.l).l,rawexp(x.h).h};}
struct J{I f,d,s;J(I v=0,I z=0,I e=0):f(v),d(z),s(e){}};
J operator+(J a,J b){return {a.f+b.f,a.d+b.d,a.s+b.s};}
J operator-(J a,J b){return {a.f-b.f,a.d-b.d,a.s-b.s};}
J operator*(J a,J b){return {a.f*b.f,a.f*b.d+a.d*b.f,a.f*b.s+a.d*b.d+a.s*b.f};}
J operator/(J a,J b){auto r=a.f/b.f,t=(a.d-r*b.d)/b.f;return {r,t,(a.s-r*b.s-t*b.d)/b.f};}
J cj(I a,J b){return {a*b.f,a*b.d,a*b.s};}
J EJ(J a){auto u=Exp(a.f);return {u,u*a.d,u*(a.s+sq(a.d)*.5)};}
J SJ(J a){auto u=sqrtI(a.f);auto w=a.d/(I(2)*u);return {u,w,(a.s-sq(w))/(I(2)*u)};}
constexpr int N=36;constexpr int MAXDEPTH=24;
constexpr double SCALE=16384.,BMAX=1./256;
constexpr double tail=0x1p-41; // 96*(N+2)^3/4^(N+1) < 2^-41; see documentation.
using Co=array<I,N+1>;using Pol=array<J,N+1>;
struct Pt{Co R,A,H;};
std::map<int,Pt> coeffcache;
Pt const& point(int T){
 auto it=coeffcache.find(T);if(it!=coeffcache.end())return it->second;
 double t=T/SCALE;Pt p{};I xi(t);
 // individual j aliases, then upper bound unretained j by integral test OR j=1 column.
 int M=t==0?512:std::min(512,std::max(3,int(40/t))); // no proof depends on this choice
 for(int j=1;j<=M;++j){
  I y=I(j)*xi;I w=Exp(-y),a=w;I jj=I(1)/I(j);I jq=jj*jj*jj;
  for(int i=1;i<=2*N+3;i++){
   w=w*y/I(i);a=a+w;
   if(i%2==0)continue;
   if(j==1&&i<=2*N+1)p.H[(i-1)/2]={max(0.,a.l),min(1.,a.h)};
   if(i<3)continue;
   auto term=I(2)*jq*I(max(0.,a.l),min(1.,a.h));int n=(i-3)/2;
   if(j==1)p.A[n]=term;else p.R[n]=p.R[n]+term;
   jq=jq*jj*jj;
  }
 }
 I m=I(1)/I(M);I zt=m*m,ft;
 if(t>0){auto d=I(1)-Exp(-xi);ft=I(M+1)*Exp(-I(M)*xi)/sq(d);}
 for(int n=0;n<=N;n++){
  double h=(zt/I(n+1)).h;
  if(t>0)h=min(h,(ft*p.A[n]).h);
  p.R[n].h=up(p.R[n].h+h);
  p.A[n]=p.A[n]+p.R[n];
  zt=zt*m*m;
 }
 return coeffcache.emplace(T,p).first->second;
}
struct Rad{I t,h;Co r,a;};
Rad radial(int l,int h,I b){
 auto &L=point(l),&H=point(h);Rad t;t.t=I(l/SCALE,h/SCALE);
 I bp(1),hv(0);for(int n=0;n<=N;n++){
  t.r[n]={H.R[n].l,L.R[n].h};t.a[n]={H.A[n].l,L.A[n].h};
  hv=hv+bp*I(H.H[n].l,L.H[n].h);bp=bp*b;
 }
 t.h=I(2)*I(hv.l,up(hv.h+(bp/(I(1)-b)).h));return t;
}
Pol polys(I q,I b){Pol a;J z(sq(q),I(-2)*q,I(1));a[0]=J(I(1));I bp(1);
 for(int k=1;k<=N;++k){bp=bp*b;a[k]=z*a[k-1]+J(bp);}
 return a;}
J series(Co const& v,Pol const& f){
 J a;for(int k=0;k<=N;++k)a=a+cj(v[k],f[k]);a.f.h=up(a.f.h+tail);
 a.d.l=dn(a.d.l-tail);a.s.h=up(a.s.h+tail);return a;}
J P(I q,I b){return {q*(sq(q)-b),b-I(3)*sq(q),I(3)*q};}
J Tjet(Rad const&t,I p,I b,J Pj,Pol const&f){
 I q=I(1)-p,z=I(1)/(I(1)+q),s=t.t,u=Exp(-(I(1)+q)*s)*z;
 J e(u,u*(s+z),u*(sq(s)*.5+s*z+sq(z)));
 return e+J(q*t.h,-t.h)-Pj*series(t.r,f);
}
J W(Rad const&t,I pp,J T){J p(pp,I(1));return J(I(1))-p*EJ(cj(t.t,p))*T;}
struct Out{I f,curv;};
struct Box{int a,b,c,d,e,g;};
struct Phase{I pp,q;Pol f;J Pj;};
std::map<std::tuple<int,int,int>,vector<Phase>> pc;
std::vector<Phase> const& makephase(int n,int b0,int b1){
 auto key=std::make_tuple(n,b0,b1);auto it=pc.find(key);if(it!=pc.end())return it->second;
 I b(BMAX*b0/64,BMAX*b1/64);std::vector<Phase> vec;vec.reserve(n);
 for(int i=0;i<n;i++){I pp(double(i)/n,double(i+1)/n),q(1-pp.h,1-pp.l);vec.push_back({pp,q,polys(q,b),P(q,b)});}
 return pc.emplace(key,std::move(vec)).first->second;
}
double atom(int x0,int x1,int y0,int y1,const Rad&x,const Rad&y,const Rad&m,I b){
 if(b.h==0)return 0;
 // at q^2=b, coefficients are (n+1)b^n, last remainder bounded by our usual tail
 I bp(1),xv(0),yv(0),mv(0);for(int k=0;k<=N;k++){
  auto w=I(k+1)*bp;xv=xv+w*x.a[k];yv=yv+w*y.a[k];mv=mv+w*m.a[k];bp=bp*b;}
 xv.h=up(xv.h+tail);yv.h=up(yv.h+tail);
 auto d=mv-sqrtI(xv*yv);
 return min(0.,(I(2)*b*d).l); // signed mass D'' atom lower min(0,-D'(L-))
}
using Table=vector<I>;
std::map<std::pair<int,int>,Table> weightcache;
I trig_raw(I x,bool sinc){
 if(x.l<0)throw std::runtime_error("trig negative");
 if(x.h>.5){auto c=trig_raw(x*.5,false);return sinc?trig_raw(x*.5,true)*c:I(2)*sq(c)-I(1);}
 int j=sinc?1:0;I term(1),v(1),U=sq(x);
 for(int k=1;k<=10;++k){term=-term*U/((2*k-1+j)*(2*k+j));v=v+term;}
 auto sm=v-term*U/((21+j)*(22+j));return {min(v.l,sm.l),max(v.h,sm.h)};
}
Table const& weights(int n,int k){auto key=std::make_pair(n,k);auto it=weightcache.find(key);if(it!=weightcache.end())return it->second;
 Table tt;for(int i=0;i<n;++i){if(k<0){I qq=I(3*i*i+3*i+1)/(3*n*n),e=I(i+1)/n;
   tt.emplace_back((qq*(I(.5)-sq(e)/24)).l,(qq/2).h);}
   else{I t=I(2*i+1)/(2*n),mu=I(2*k+1)/4;
   I a=I(1)-trig_raw(t*mu,false)*trig_raw(mu/(2*n),true);tt.emplace_back(max(0.,(a-t*.25).l),min(2.,(a+t*.25).h));}}
 return weightcache.emplace(key,std::move(tt)).first->second;}
int freq(vector<Out> const&pr,double atm){
 int n=pr.size();I mass=I(2*n)*atm,tv=std::abs(pr.back().curv.l);
 bool positive=atm==0;
 for(int i=0;i<n;i++){mass=mass+pr[i].curv.l;if(pr[i].curv.l<0)positive=false;if(i>0){auto s=I(pr[i].curv.l)-I(pr[i-1].curv.l);tv=tv+max(-s.l,s.h);}}
 if(positive)return 1;
 if(mass.l<=0)return 0;
 auto om=(tv*I(n)/mass).h;if(om>128)return 0;
 auto test=[&](int k){
   auto &ww=weights(n,k);I s=I(n)*atm*(k==-1?.5:2);
   for(int i=0;i<n;i++)s=s+I(pr[i].curv.l)*(pr[i].curv.l>=0?ww[i].l:ww[i].h);return s.l>0;};
 if(!test(-1))return 0;for(int k=2;k<2*std::ceil(om);k++)if(!test(k))return 0;
 return 2;
}
static int acceptedExp; static double margin;
int eval(Box B,int n){
 int x0=B.a,x1=B.b,y0=B.c,y1=B.d,b0=B.e,b1=B.g;
 if((x0+y0)%2||(x1+y1)%2)throw std::runtime_error("midpoint");
 I b(BMAX*b0/64,BMAX*b1/64);auto x=radial(x0,x1,b),y=radial(y0,y1,b),m=radial((x0+y0)/2,(x1+y1)/2,b);
 I c=(y.t-x.t)*.5;vector<Out> pr;pr.reserve(n);
 auto const& ph=makephase(n,b0,b1);
 for(auto const&t:ph){
  J D; if(t.pp.h<=.5){
   auto tx=Tjet(x,t.pp,b,t.Pj,t.f),ty=Tjet(y,t.pp,b,t.Pj,t.f),tm=Tjet(m,t.pp,b,t.Pj,t.f);
   auto u=SJ(W(x,t.pp,tx)),v=SJ(W(y,t.pp,ty)),G=tx/(J(I(1))+u);
   J ep(t.pp,I(1));auto co=EJ(cj(-c,ep)),eo=EJ(cj(c,ep));
   D=co*G-tm+eo*u*ty/(J(I(1))+v);
  } else{
   auto X=series(x.a,t.f),Y=series(y.a,t.f),M=series(m.a,t.f);
   D=t.Pj*(M-SJ(X*Y));}
  I f=D.f,w=I(2)*D.s;
  // zero extension where q<kappa; b0/b1 are rational, no irrational comparison
  if(t.q.l*t.q.l < b.h){f=I(min(0.,f.l),max(0.,f.h));w=I(min(0.,w.l),max(0.,w.h));}
  pr.push_back({f,w});
 }
 double at=atom(x0,x1,y0,y1,x,y,m,b);
 // Union: parameter box can be accepted by any of these Helmholtz multipliers
 for(int k=-3;k<=4;++k){I ff=std::ldexp(1.,k);I ll=sq(c*ff);
  I sum=I(2*n)*at;
  for(auto const &p:pr){auto l=ll*p.f;sum=sum+min(l.l,(I(2)*p.curv-l).l);}
  if(sum.l>0){acceptedExp=k;margin=(sum/I(n)).l;return 3;}}
 return freq(pr,at);
}
int S(double v){int i=int(v*SCALE);if(double(i)!=v*SCALE)throw std::runtime_error("non-dyadic");return i;}
int main(int ac,char**av){
 if(ac!=8&&ac!=10){puts("x0 x1 y0 y1 n_initial depth maxvisits [b0 b1] ; b indices in 0..64");return 2;}
 if(std::fegetround()!=FE_TONEAREST)throw std::runtime_error("rounding mode");
 volatile double tiny=std::numeric_limits<double>::denorm_min(); if(tiny+tiny<=tiny)throw std::runtime_error("subnormal mode");
 const char* traceName=std::getenv("PSD_ZERO_TRACE");FILE*trace=traceName?std::fopen(traceName,"w"):nullptr;if(traceName&&!trace)throw std::runtime_error("write trace");if(trace)std::fprintf(trace,"# x0 x1 y0 y1 b0 b1 (denominator=16384 throughout) n method multiplier_exponent positive_lower_dyadic_hex\n");
 Box root{S(atof(av[1])),S(atof(av[2])),S(atof(av[3])),S(atof(av[4])),ac==10?atoi(av[8]):0,ac==10?atoi(av[9]):64};
 int n0=atoi(av[5]),dep=atoi(av[6]),mx=atoi(av[7]);
 if(!(0<=root.a&&root.a<=root.b&&root.b<=3*SCALE&&4*SCALE<=root.c&&root.c<=root.d&&root.d<=24*SCALE&&0<=root.e&&root.e<=root.g&&root.g<=64&&n0>=4&&n0<=1024&&(n0&(n0-1))==0&&dep>=0&&dep<=24))throw std::runtime_error("range");
 std::vector<std::pair<Box,int>> todo{{root,0}};int vis=0,passes=0,method[4]{};
 while(!todo.empty()&&vis++<mx){
  auto [v,d]=todo.back();todo.pop_back();int ok=0;
  // Refinement in radial parameters eventually refines spectral mesh.
  int n=min(1024,n0*(1<<min(5,max(0,(d-4)/3))));
  try{ok=eval(v,n);}catch(std::domain_error const&){ok=0;}
  if(ok){passes++;method[ok]++;if(trace)std::fprintf(trace,"%d %d %d %d %d %d %d %d %d %a\n",v.a,v.b,v.c,v.d,v.e,v.g,n,ok,ok==3?acceptedExp:0,ok==3?margin:0.);continue;}
  if(d>=dep){printf("INCOMPLETE x=[%.8f,%.8f] y=[%.8f,%.8f] b=[%d,%d]/16384 n=%d visits=%d passed=%d\n",v.a/SCALE,v.b/SCALE,v.c/SCALE,v.d/SCALE,v.e,v.g,n,vis,passes);return 1;}
  auto L=v,R=v;int xw=v.b-v.a,yw=v.d-v.c,bw=v.g-v.e;
  if(xw>=yw/3 && xw>=8*bw && xw>=4){R.a=L.b=(v.a+v.b)/2;}
  else if(yw>=24*bw&&yw>=4){R.c=L.d=(v.c+v.d)/2;}
  else if(bw>=2){R.e=L.g=(v.g+v.e)/2;}
  else if(xw>=4){R.a=L.b=(v.a+v.b)/2;}
  else if(yw>=4){R.c=L.d=(v.c+v.d)/2;}else{puts("INCOMPLETE resolution");return 1;}
  todo.push_back({R,d+1});todo.push_back({L,d+1});
  if(vis%200==0)fprintf(stderr,"visits %d accepted %d x=%.6g y=%.6g d=%d\n",vis,passes,v.a/SCALE,v.c/SCALE,d);
 }
 if(trace)std::fclose(trace);
 if(!todo.empty()){puts("INCOMPLETE visit cap");return 1;}printf("PASS x=[%.5g,%.5g] y=[%.5g,%.5g] b=[%d,%d]/16384 visits=%d leaves=%d methods=%d/%d/%d\n",root.a/SCALE,root.b/SCALE,root.c/SCALE,root.d/SCALE,root.e,root.g,vis,passes,method[1],method[2],method[3]);return 0;
}
