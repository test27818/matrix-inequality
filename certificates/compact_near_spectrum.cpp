// Compact near-diagonal integer/rational-input ball verifier for kappa=0.
// Outward IEEE-754 interval arithmetic, transcendental endpoints from MPFR.
// See compact-near-spectrum.md. No fastmath, fused operations, or sampled sign.
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>
#ifndef COMPACT_PORTABLE
#include <gmp.h>
extern "C" {
struct mpfr_obj {long prec; int sign; long exp; mp_limb_t *d;};
void mpfr_init2(mpfr_obj *,long); void mpfr_clear(mpfr_obj *);
int mpfr_set_d(mpfr_obj *,double,int);double mpfr_get_d(const mpfr_obj *,int);
int mpfr_exp(mpfr_obj *,const mpfr_obj *, int);
int mpfr_sin(mpfr_obj *,const mpfr_obj *, int);
int mpfr_sqrt(mpfr_obj *,const mpfr_obj *, int);
}
#endif
using std::max;using std::min;
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53);
struct I{double l,h;I(double l,double h):l(l),h(h){if(!(l<=h&&std::isfinite(l)&&std::isfinite(h)))throw std::range_error("interval");} I(double z=0):I(z,z){}};
double adj(double x,bool pos){
 if(!std::isfinite(x))throw std::range_error("not finite");
 if(!x)return pos?std::numeric_limits<double>::denorm_min():-std::numeric_limits<double>::denorm_min();
 auto b=std::bit_cast<uint64_t>(x);
 return std::bit_cast<double>((pos==(x>0))?b+1:b-1);
}
double up(double x){return adj(x,true);} double dn(double x){return adj(x,false);}
I operator+(I a,I b){return {dn(a.l+b.l),up(a.h+b.h)};}
I operator-(I a){return {-a.h,-a.l};}
I operator-(I a,I b){return a+-b;}
I operator*(I a,I b){if(a.l>=0 && b.l>=0)return {dn(a.l*b.l),up(a.h*b.h)};double ll=a.l*b.l,lh=a.l*b.h,hl=a.h*b.l,hh=a.h*b.h;return {dn(std::min({ll,lh,hl,hh})),up(std::max({ll,lh,hl,hh}))};}
I operator/(I a,I b){if(b.l<=0)throw std::range_error("division");return a*I(dn(1/b.h),up(1/b.l));}
I sq(I t){if(t.l>=0)return {max(0.,dn(t.l*t.l)),up(t.h*t.h)};if(t.h<=0)return sq(-t);return {0,up(max(t.l*t.l,t.h*t.h))};}
I hull(I a,I b){return {min(a.l,b.l),max(a.h,b.h)};}
I isect(I a,I b){return {max(a.l,b.l),min(a.h,b.h)};}
double mag(I a){return max(-a.l,a.h);}
#ifndef COMPACT_PORTABLE
struct M {mpfr_obj x,v; M(){mpfr_init2(&x,128);mpfr_init2(&v,128);}~M(){mpfr_clear(&v);mpfr_clear(&x);}};
I trans(double t,int j){static M m;static std::array<std::unordered_map<double,I>,3> tab;
 auto &tt=tab[j];auto it=tt.find(t);if(it!=tt.end())return it->second;if(tt.size()>250000)tt.clear();
 mpfr_set_d(&m.x,t,0);auto F=(j==0?mpfr_exp:(j==1?mpfr_sqrt:mpfr_sin));
 F(&m.v,&m.x,3);double l=mpfr_get_d(&m.v,3);F(&m.v,&m.x,2);double h=mpfr_get_d(&m.v,2);
 return tt.emplace(t,I(l,h)).first->second;}
#else
// Auditable alternative to MPFR: range-reduced Taylor, square root certified
// by outward squaring. All functions here are evaluated at a binary64 point.
I elementary(double x,int j){
 if(j==1){
  if(!(x>0))throw std::range_error("positive square root");
  double l=std::sqrt(x),h=l;
  for(int i=0;i<16&&up(l*l)>x;i++)l=dn(l);
  for(int i=0;i<16&&dn(h*h)<x;i++)h=up(h);
  if(l<0||up(l*l)>x||dn(h*h)<x)throw std::range_error("sqrt guard");
  return {l,h};
 }
 if(j==0 && x<0)return I(1)/elementary(-x,0);
 if(j==2 && x<0)return -elementary(-x,2);
 if(x<0||x>700)throw std::range_error("positive Taylor input");
 int k=0;while(x>.25){x*=.5;++k;}I u(x);
 if(j==0){I v=1,t=1;for(int i=1;i<=16;i++){t=t*u/i;v=v+t;}
  v=I(v.l,(v+2*(t*u/17)).h);while(k-->0)v=sq(v);return v;
 }
 I z=sq(u),co=1,si=u,tc=1,ts=u;
 for(int i=1;i<=8;i++){tc=-tc*z/((2*i-1)*(2*i));ts=-ts*z/((2*i)*(2*i+1));co=co+tc;si=si+ts;}
 co=hull(co,co-tc*z/(17*18));si=hull(si,si-ts*z/(18*19));
 while(k-->0){si=2*si*co;co=2*sq(co)-1;}return si;
}
I trans(double z,int j){
 static std::array<std::unordered_map<double,I>,3> tab;auto &t=tab[j];auto it=t.find(z);if(it!=t.end())return it->second;if(t.size()>250000)t.clear();
 return t.emplace(z,elementary(z,j)).first->second;
}
#endif
I ex(I t){return {trans(t.l,0).l,trans(t.h,0).h};}I rt(I t){if(t.l<=0)throw std::range_error("sqrt");return {trans(t.l,1).l,trans(t.h,1).h};}
// Taylor jet f,p derivative, p^2 derivative /2.
struct J{I a,b,c;J(I a=0,I b=0,I c=0):a(a),b(b),c(c){}J(double z):a(z),b(0),c(0){}};
J operator+(J u,J v){return {u.a+v.a,u.b+v.b,u.c+v.c};}J operator-(J u){return {-u.a,-u.b,-u.c};}J operator-(J u,J v){return u+-v;}
J operator*(J u,J v){return {u.a*v.a,u.a*v.b+u.b*v.a,u.a*v.c+u.b*v.b+u.c*v.a};}
J jsq(J v){return {sq(v.a),2*v.a*v.b,2*v.a*v.c+sq(v.b)};}
J operator/(J u,J v){auto a=u.a/v.a,b=(u.b-a*v.b)/v.a;return {a,b,(u.c-a*v.c-b*v.b)/v.a};}
J je(J v){auto a=ex(v.a);return {a,a*v.b,a*(v.c+sq(v.b)/2)};}
J jr(J v){auto a=rt(v.a),b=v.b/(2*a);return {a,b,(v.c-sq(b))/(2*a)};}
// second independent variable c, Taylor convention again
struct K{J a,b,c;K(J a=0,J b=0,J c=0):a(a),b(b),c(c){} K(I v):a(v){} K(double v):a(v){}};
K operator+(K u,K v){return {u.a+v.a,u.b+v.b,u.c+v.c};}K operator-(K u){return {-u.a,-u.b,-u.c};}K operator-(K u,K v){return u+-v;}
K operator*(K u,K v){return {u.a*v.a,u.a*v.b+u.b*v.a,u.a*v.c+u.b*v.b+u.c*v.a};}
K operator/(K u,K v){J a=u.a/v.a,b=(u.b-a*v.b)/v.a;return {a,b,(u.c-a*v.c-b*v.b)/v.a};}
K ke(K v){J a=je(v.a);return {a,a*v.b,a*(v.c+jsq(v.b)/2)};}
K kr(K v){J a=jr(v.a),b=v.b/(2*a);return {a,b,(v.c-jsq(b))/(2*a)};}
I err(double a){return {-a,a};}
const J eJ(err(0x1p-33),err(0x1p-33),err(0x1p-33));
#ifdef EXTENDED_COMPACT_SPECTRUM
int LAST=4; constexpr double min_x_permitted=3,max_half_gap=12;
#else
constexpr int LAST=3; constexpr double min_x_permitted=16,max_half_gap=4;
#endif
struct TJ{J f,dt,dtt;}; // radial derivatives; all 3 p slots are Taylor coeffs.
TJ lowT(I t,J p){
 J q=1-p;
 TJ r{};
 for(int j=1;j<=LAST;++j){
  J a=j+q,b=j+p,sa=1-a,sb=1-b;
  auto eh=ex((1-j)*t); I L=t+I(1)/j;
  J h=q*J(2*eh); J ha=h*J(L),hd=h*(1+J((1-j)*L)),hdd=h*(2*(1-j)+J((1-j)*(1-j)*L));
  J A=(j*je(sa*J(t)))/a,B=((j+1)*je(sb*J(t)))/b;
  r.f=r.f+ha+A-B;r.dt=r.dt+hd+A*sa-B*sb;r.dtt=r.dtt+hdd+A*jsq(sa)-B*jsq(sb);
 }
 r.f=r.f+eJ;r.dt=r.dt+eJ;r.dtt=r.dtt+eJ;
 return r;
}
// e^t J_t(q)/q^3 and first two t-derivatives for min_x_permitted<=t<=32, 0<=q<=3/32.
TJ highS(I t,J q){
 TJ r{};J zz=jsq(q);
 for(int j=1;j<=LAST;++j){
  const int K=60; // q-series n=0..30, last radial degree 63; additionally tail by Cauchy estimate.
  I e=ex((1-j)*t),v=j*t,ti=1,mom=1,prev=1,pprev=1;
  J fr=0,dr=0,rr=0; // Horner on q^2/j^2 over precomputed moments
  std::array<I,31> arr,ad,add;
  for(int k=1;k<=K+3;++k){ti=ti*v/k;pprev=prev;prev=mom;mom=mom+ti;
   if((k&1)&&k>=3){
    int n=(k-3)/2;
    arr[n]=e*mom;ad[n]=e*((1-j)*mom+j*prev);
    add[n]=e*((1-j)*(1-j)*mom+2*(1-j)*j*prev+j*j*pprev);
   }
  }
  J z=zz/J(j*j);
  for(int n=30;n>=0;--n){fr=fr*z+J(arr[n]);dr=dr*z+J(ad[n]);rr=rr*z+J(add[n]);}
  double j3=j*j*j;
  r.f=r.f+2*fr/J(j3);r.dt=r.dt+2*dr/J(j3);r.dtt=r.dtt+2*rr/J(j3);
 }
 r.f=r.f+eJ;r.dt=r.dt+eJ;r.dtt=r.dtt+eJ;return r;
}
K lift(TJ t,int s){return {t.f,s*t.dt,(s*s)*t.dtt/2};}
struct Table{std::map<std::pair<double,double>,TJ> vals;};
TJ value(I t,J p,bool low,Table& cache){
 auto key=std::make_pair(t.l,t.h);auto hit=cache.vals.find(key);if(hit!=cache.vals.end())return hit->second;
 auto v=low?lowT(t,p):highS(t,1-p);return cache.vals.emplace(key,v).first->second;
}
// true e^m D and second c derivative/2, for fixed x (m=x+c, y=x+2c).
K full(I x,I ci,J p,bool low,Table& cache){
 I m=x+ci,y=x+2*ci; // exact outward expansion harmless; clamped to theorem domain for alias estimates
 K c=K(J(ci),1),t=K(x),u=K(J(m),1),v=K(J(y),2);
 K Tx=lift(value(x,p,low,cache),0),Ty=lift(value(y,p,low,cache),2),Tm=lift(value(m,p,low,cache),1);
 if(!low){K q=1-K(p);return q*q*q*(Tm-kr(Tx*Ty));}
 K q=1-K(p),wx=kr(1-K(p)*ke(-q*t)*Tx),wy=kr(1-K(p)*ke(-q*v)*Ty);
 return ke(q*c)*Tx/(1+wx)-Tm+ke(-q*c)*wx*Ty/(1+wy);
}
// intersect direct and Peano averages with exact triangular weights, no finite-difference at c=0.
J radial(I x,I c,J p,bool low,Table& cache,int nu){
 J d;bool have=false;
 if(c.l>0){d=full(x,c,p,low,cache).a/J(sq(c));have=true;}
 if(c.l==0 || c.l<1.25 || !have){
  J h=0;for(int i=0;i<nu;++i){I arm=I(double(i)/nu,double(i+1)/nu)*c;
    J b=full(x,arm,p,low,cache).c;
    h=h+J(I(2*nu-2*i-1)/(nu*nu))*b;
  }
  if(have){return {isect(d.a,h.a),isect(d.b,h.b),isect(d.c,h.c)};}return h;
 }
 return d;
}
I cosav(int i,int np,int g,int nfreq){ // average over p cell at freq center (2g+1)/(2nfreq)
 double w=double(2*g+1)/(2*nfreq),a=double(i)/np,b=double(i+1)/np;
 return (trans(b*w,2)-trans(a*w,2))/(w/np);
}
const int NF=4;int maxomega=0;
double global_margin=std::numeric_limits<double>::max();
struct SC{bool yes=false; double E=0,delta=0,Mf=0,goal=0; int reason=0;};
SC sign(std::vector<J> const& prof,I delta){
 int n=prof.size();double E=0,V=0,mag_der=0,Ev=0,Mv=0,prev=0;
 std::vector<double> f,d;for(auto z:prof){I w=2*z.c;double wi=(w.l+w.h)/2,fi=(z.a.l+z.a.h)/2;
  f.push_back(fi);d.push_back(wi);
  E=up(E+max(up(wi-w.l),up(w.h-wi)));Ev=up(Ev+max(up(fi-z.a.l),up(z.a.h-fi)));
  if(f.size()>1)V=up(V+up(std::abs(wi-prev)));prev=wi;
 }
 E=up(E/n);Ev=up(Ev/n);V=up(V+std::abs(prev));
 double rem=dn(delta.l-E);if(rem<=0)return {false,E,delta.l,0,0,1};
 if(V>rem*240)return {false,E,delta.l,V,rem,2};
 int om=max(1,1+int(std::ceil(up(V/rem))));
 double cert=dn(rem-up(V/om)); if(!(cert>0))return {false,E,delta.l,0,0,5};
 // zero to om via either direct or twice-integrated triangle bound. Lipschitz in omega times |piecewise coefficient| times cell centre.
 for(int i=0;i<n;++i){double mid=(i+.5)/n;mag_der=up(mag_der+up(std::abs(d[i])*mid));Mv=up(Mv+up(std::abs(f[i])*mid));}
 auto errd=(I(E)+I(mag_der)/(2*NF*n)).h,errf=(I(Ev)+I(Mv)/(2*NF*n)).h;
 static std::map<std::pair<int,int>,std::vector<I>> wt;
 for(int g=0;g<NF*om;++g){
  auto key=std::make_pair(n,g);auto it=wt.find(key);if(it==wt.end()){
   std::vector<I>w;for(int i=0;i<n;++i)w.push_back(cosav(i,n,g,NF));it=wt.emplace(key,std::move(w)).first;
  }
  I sd=0,sf=0;for(int i=0;i<n;++i){sd=sd+it->second[i]*d[i];sf=sf+it->second[i]*f[i];}
  double good=max(dn((delta-sd/I(n)).l-errd),dn((sf/I(n)).l-errf));
  if(good>0){cert=min(cert,good);continue;}
  return {false,E,delta.l,errd,g/(double)NF,3};
 }
 global_margin=min(global_margin,cert);maxomega=max(maxomega,om);return {true};
}
I delta(I x,I c,int nu){ // e^m(-D'(0))/c^2 = (2-A) e^m A/2, A=mean_[x,y] (1-C'(t)).
 I A=0,B=0;
 for(int i=0;i<nu;++i){I z=I(2.*i/nu,2.*(i+1)/nu),t=x+z*c;
 I ev=ex(-t),v=2*(t-1+ev)/sq(1-ev);
 A=A+ev*v;B=B+ex((1-z)*c)*v;}
 return (2-A/nu)*(B/nu)/2;
}
struct PBox{I x,c;int depth;};
int main(int argc,char**argv){
 if(argc!=8){std::cerr<<"xlo xhi clo chi nq maxdepth nu, dyadic inputs only\n";return 2;}
 // custom inputs use exact dyadics as parsed by strtod; all built-in roots integral.
 auto arg=[&](int k){char*e=nullptr;double v=strtod(argv[k],&e);if(!e||*e)throw std::runtime_error("input");return v;};
 I x(arg(1),arg(2)),c(arg(3),arg(4));int n=atoi(argv[5]),dep=atoi(argv[6]),nu=atoi(argv[7]);
 if(!(x.l>=min_x_permitted&&x.h<=24&&c.l>=0&&c.h<=max_half_gap&&x.h+2*c.h<=32&&n>=32&&n<=4096&&(n&(n-1))==0&&dep>=0&&dep<=24&&nu>=1&&nu<=64&&(nu&(nu-1))==0))return 2;
 #ifdef EXTENDED_COMPACT_SPECTRUM
 LAST=x.l>=16?3:x.l>=12?4:x.l>=8?7:x.l>=6?9:x.l>=4?15:23;
#endif
 std::vector<PBox> st{{x,c,0}};int visited=0,passed=0,maxD=0;
 // cache safe: phase cells never share their t data.
 std::vector<Table> cache(n);
 while(!st.empty()){
  PBox t=st.back();st.pop_back();++visited;maxD=max(t.depth,maxD);
  std::vector<J> pvec;SC s;char const*guard=nullptr;
  try{
   for(int i=0;i<n;++i){if(cache[i].vals.size()>15000)cache[i].vals.clear();
    J p(I(double(i)/n,double(i+1)/n),1);bool low=i<29*n/32;
    pvec.push_back(radial(t.x,t.c,p,low,cache[i],nu));
   }
   s=sign(pvec,delta(t.x,t.c,16));
  }catch(const std::range_error&){s.reason=4;guard="interval";}
  if(s.yes){++passed;continue;}
  if(t.depth>=dep){std::cout << "UNCERTIFIED visits="<<visited<<" pass="<<passed<<" bad x="<<t.x.l<<","<<t.x.h<<" c="<<t.c.l<<","<<t.c.h<<" fail="<<s.reason<<" E="<<s.E<<" delta="<<s.delta<<" remaining="<<s.Mf<<" phase="<<s.goal<<"\n";return 1;}
  if(t.x.h-t.x.l>=2*(t.c.h-t.c.l)){double m=(t.x.l+t.x.h)/2;st.push_back({I(m,t.x.h),t.c,t.depth+1});st.push_back({I(t.x.l,m),t.c,t.depth+1});}
  else{double m=(t.c.l+t.c.h)/2;st.push_back({t.x,I(m,t.c.h),t.depth+1});st.push_back({t.x,I(t.c.l,m),t.depth+1});}
 }
 std::cout<<"CERTIFIED x=["<<x.l<<","<<x.h<<"] c=["<<c.l<<","<<c.h<<"] cells="<<n<<" peano="<<nu<<" visits="<<visited<<" leaves="<<passed<<" depth="<<maxD<<" maxfreq="<<maxomega<<" min_scaled_margin_lower="<<std::hexfloat<<global_margin<<"\n";
}
