// Certified continuum (kappa=0) signed-curvature replay on compact x>=3.
// Uses the integer frequency verifier from signed_curvature_fast.cpp unchanged.
#define main legacy_signed_main
#include "signed_curvature_fast.cpp"
#undef main
#include <bit>
#include <cfenv>
#include <iomanip>
#include <string>
namespace nearcert {
using D=double;
static_assert(std::numeric_limits<D>::is_iec559 && sizeof(D)==8 && std::numeric_limits<D>::digits==53);
D hike(D x,bool up){if(!std::isfinite(x))throw std::range_error("nonfinite");
 if(x==0) return up?std::numeric_limits<D>::denorm_min():-std::numeric_limits<D>::denorm_min();
 auto b=std::bit_cast<uint64_t>(x); return std::bit_cast<D>(up==(x>0)?b+1:b-1);}
struct R{D a,b;R():a(0),b(0){} R(D a):a(a),b(a){} R(D a,D b):a(a),b(b){if(!(a<=b)) throw std::range_error("interval");}
 static R out(D a,D b){return {hike(a,0),hike(b,1)};}};
R operator+(R x,R y){return R::out(x.a+y.a,x.b+y.b);}
R operator-(R x){return {-x.b,-x.a};}R operator-(R x,R y){return x+(-y);}
R operator*(R x,R y){if(x.a>=0){if(y.a>=0)return R::out(x.a*y.a,x.b*y.b);if(y.b<=0)return R::out(x.b*y.a,x.a*y.b);return R::out(x.b*y.a,x.b*y.b);}
 if(x.b<=0)return -((-x)*y); if(y.a>=0)return R::out(x.a*y.b,x.b*y.b);if(y.b<=0)return - (x*(-y));
 return R::out(std::min(x.a*y.b,x.b*y.a),std::max(x.a*y.a,x.b*y.b));}
R inv(R x){if(x.a<=0)throw std::range_error("division");return R::out(1/x.b,1/x.a);}
R operator/(R x,R y){return x*inv(y);}
R sqr(R x){if(x.a>=0)return R::out(x.a*x.a,x.b*x.b); if(x.b<=0)return sqr(-x);D u=std::max(-x.a,x.b);return R(0,hike(u*u,1));}
R sr(R x){if(x.a<=0)throw std::range_error("root");D l=std::sqrt(x.a),h=std::sqrt(x.b);
 for(int k=0;k<8&&hike(l*l,true)>x.a;k++)l=hike(l,false);
 for(int k=0;k<8&&hike(h*h,false)<x.b;k++)h=hike(h,true);
 if(hike(l*l,true)>x.a||hike(h*h,false)<x.b)throw std::runtime_error("root proposal failed");
 return {l,h};}
R hull(R x,R y){return {std::min(x.a,y.a),std::max(x.b,y.b)};}
R inter(R x,R y){D a=std::max(x.a,y.a),b=std::min(x.b,y.b);if(a>b)throw std::logic_error("inconsistent enclosures");return {a,b};}
// Exact real interval for positive exponential by range reduction and geometric tail.
R ep0(D t){if(t<0)return inv(ep0(-t));if(t>200)throw std::range_error("exp range");
 int n=0;while(t>.0625){t*=.5;n++;} R z(t),sum(1),term(1);
 for(int i=1;i<=12;i++){term=term*z/R(i);sum=sum+term;}
 sum.b=hike(sum.b+(term*z/R(13)*R(2)).b,1);while(n--)sum=sqr(sum);return sum;}
R ep(D t){static std::unordered_map<double,R> E; auto it=E.find(t);if(it!=E.end())return it->second;if(E.size()>150000)E.clear();R v=ep0(t);return E.emplace(t,v).first->second;}
R er(R z){return {ep(z.a).a,ep(z.b).b};}
struct J{R t,d,e;J()=default;J(double z):t(z){}J(R z):t(z){}J(R a,R b,R c):t(a),d(b),e(c){}};
J operator+(J a,J b){return {a.t+b.t,a.d+b.d,a.e+b.e};}J operator-(J a){return {-a.t,-a.d,-a.e};}J operator-(J a,J b){return a+(-b);}
J operator*(J a,J b){return {a.t*b.t,a.d*b.t+a.t*b.d,a.e*b.t+a.d*b.d+a.t*b.e};}
J operator*(J a,R s){return {a.t*s,a.d*s,a.e*s};}J operator*(R s,J a){return a*s;}
J ij(J a){R a0=inv(a.t),a1=-(a.d*a0*a0);return {a0,a1,-(a.d*a1+a.e*a0)*a0};}J operator/(J a,J b){return a*ij(b);}
J ej(J a){R t=er(a.t);return {t,t*a.d,t*(a.e+sqr(a.d)*R(.5))};}
J rj(J a){R t=sr(a.t),rr=inv(t*R(2)),d=a.d*rr;return {t,d,(a.e-sqr(d))*rr};}
J square(J a){return {sqr(a.t),a.t*a.d*R(2),sqr(a.d)+a.t*a.e*R(2)};}
J intersect(J a,J b){return {inter(a.t,b.t),inter(a.d,b.d),inter(a.e,b.e)};}
J unionj(J a,J b){return {hull(a.t,b.t),hull(a.d,b.d),hull(a.e,b.e)};}
using Vec=std::array<J,3>;
struct Rad{Vec f;J root;};
J unit_err(){R u(-0x1p-45,0x1p-45);return {u,u,u};}
// q<=.5: C_k(t)=sum_l t^(2l)/(2l+k)!, k=0..3; derivatives via positive series.
std::array<J,4> ck4(R v,J q){std::array<J,4> c{};J z=square(q*v);J u(1);
 int n=34; for(int i=0;i<=n;i++){R f(1);for(int k=0;k<4;k++){if(k)f=f/R(2*i+k);c[k]=c[k]+u*f;}
  u=(u*z)*inv(R((2*i+1)*(2*i+2)));}
 R eps(0,0x1p-45);for(auto& k:c)k=k+J(eps,eps,eps);return c;}
Rad radial(R v,J q,bool highp){
 if(v.a<3||v.b>24)throw std::range_error("v");
 J p=J(1)-q; R ev=er(-v),mult(1); if(v.a>0)ev.a=std::max(D(0),ev.a);
 int n=1;while(n*v.a<105)n++; // tail j>=n+1. n*v.a exact at model grid
 Vec a{};R one(1);
 if(!highp){J eq=ej((-q)*v),epp=ej((-p)*v);
  for(int j=1;j<=n;j++){
    R u(j),l(j-1),s(j+1);J h=q+J(l),c=q*R(2)+J(l);R z=v+inv(u);
    J f=q*R(2),g=(J(u)/(J(u)+q))*eq,k=(J(s)/(J(u)+p))*epp;
    a[0]=a[0]+(f*z+g-k)*mult;
    a[1]=a[1]+(f*(J(1)-h*z)-g*c+k*u)*mult;
    a[2]=a[2]+(f*(square(h)*z-h*R(2))+g*square(c)-k*(u*u))*mult;
    mult=mult*ev;}
  for(J &s:a)s=s+unit_err();
  J wt=J(1)-p*eq*a[0];return {a,rj(wt)};
 }else{
  if(q.t.b>.5)throw std::range_error("highp");
  auto cf=ck4(v,q);
  for(int j=1;j<=n;j++){R u(j),l(j-1),t=u*v,t2=sqr(t),t3=t2*t;
    J fac=ij(J(1)-square(q*inv(u)))*(R(2)*mult*inv(u*u*u));
    J p1=cf[0]+cf[1]*t, p2=J(1)+cf[1]*t+cf[2]*t2, p3=J(1)+J(t)+cf[2]*t2+cf[3]*t3;
    a[0]=a[0]+fac*p3;
    a[1]=a[1]+fac*(p2*u-p3*l);
    a[2]=a[2]+fac*(p1*(u*u)-p2*(R(2)*l*u)+p3*(l*l));
    mult=mult*ev;}
  for(J &s:a)s=s+unit_err();
  return {a,rj(a[0])};
 }
}
Rad join(Rad x,Rad y){for(int a=0;a<3;a++)x.f[a]=unionj(x.f[a],y.f[a]);x.root=unionj(x.root,y.root);return x;}
struct Cache{std::unordered_map<int64_t,Rad> data;int bp;J q;bool hp;
 Rad get(int i,int j){if(i>j)throw std::logic_error("cell");int64_t id=(int64_t(i)<<32)|unsigned(j);
  auto it=data.find(id);if(it!=data.end())return it->second;
  Rad r;int len=j-i;if(len<=1)r=radial(R(D(i)/bp,D(j)/bp),q,hp);
  else if((len&(len-1))==0)r=join(get(i,(i+j)/2),get((i+j)/2,j));
  else {int span=1<<(31-__builtin_clz((unsigned)len));r=join(get(i,i+span),get(j-span,j));}
  data.emplace(id,r);return r;}
};
struct Box {int i,j;};
// i,i+1 and j,j+1 are radial boxes (point option for testing).
J block(Cache& C,int i,int j,bool point,int mx){
 int d=point?0:1,b=C.bp;J q=C.q;bool hi=C.hp;R X(D(i)/b,D(i+d)/b),Y(D(j)/b,D(j+d)/b);
 R m((X.a+Y.a)*.5,(X.b+Y.b)*.5),h(std::max(0.,(Y.a-X.b)*.5),(Y.b-X.a)*.5);
 Rad x=C.get(i,i+d),y=C.get(j,j+d);
 auto psden=square(x.root+y.root); J u,hess;
 // Explicit endpoint differences are efficient for separated radii. The Peano version
 // handles h=0 without dividing by the radius difference.
 bool endpoint=(h.a>.25);
 J linear,first;
 if(endpoint){
  // dyadic outward enclosure of the midpoint in the radial grid
  int im=(i+j)/2,ip=(i+j+2*d+1)/2;
  Rad mid=C.get(im,ip);J xw=x.f[0],yw=y.f[0];
  if(!hi){xw=xw*ej(q*h);yw=yw*ej((-q)*h);}
  linear=((xw+yw)*R(.5)-mid.f[0])/J(sqr(h));
  first=(yw-xw)/J(R(2)*h);
 }
 if(!endpoint||h.b<4){
    int M=(b>=32?mx:(h.b<.6?4:(h.b<1.5?8:mx)));J hs(0),fs(0);
    for(int k=0;k<M;k++){
      // Interpolate v=(1-t)x+ty; range in original grid coordinates is
      // [i+(j-i)k/M, i+d+(j-i)(k+1)/M].
      int lo=i+(j-i)*k/M, up=i+d+((j-i)*(k+1)+M-1)/M;
      auto v=C.get(lo,up);R t((2.*k-M)/M,(2.*k+2-M)/M);
      J w=hi?J(1):ej((-q)*(h*t));
      R lam((2.*std::min(k,M-1-k)+1)/(M*M)); // normalized triangular weight integrates to 1/2 in t=-1..1 convention.
      hs=hs+(w*v.f[2])*lam;fs=fs+(w*v.f[1])*R(1./M);
    }
    if(endpoint){linear=intersect(linear,hs);first=intersect(first,fs);}else{linear=hs;first=fs;}
 }
 J ret;if(hi)ret=(square(first)*R(2)/psden-linear)*square(q)*q;
 else ret=linear+(J(1)-q)*ej((-q)*m)*square(first)*R(2)/psden;
 return ret;
}
struct Cur{std::vector<I64> a,au,v,vu;}; // arbitrary fixed denominator: 2^28; integer sign test scale-invariant.
I64 integerlo(D x){if(x<-1e9)return INT64_MIN;if(x>1e9)x=1e9;return I64(std::floor(std::ldexp(x,28)));}
R d0(int i,int j,int b,bool pt){
 const int K=16;int a=pt?0:1;
 R x(D(i)/b,D(i+a)/b),h(std::max(0.,(j-i-a)*.5/b),(j-i+a)*.5/b);
 R A=0,B=0;for(int k=0;k<K;k++){R r(2.*k/K,2.*(k+1)/K),v=x+h*r,e=er(-v);R z=R(2)*(v-R(1)+e)/sqr(R(1)-e);
  A=A+e*z;B=B+er((R(1)-r)*h)*z;}
 return (R(2)-A/R(K))*B/R(2*K);
}
int both_freq(Cur const&g,I64 delta){
 auto const&lo=g.a;int n=lo.size();std::vector<I64> c(n),d(n);W E=0,F=0,T=0;
 for(int j=0;j<n;j++){int i=n-1-j; if(lo[i]==INT64_MIN||g.au[i]==INT64_MAX||g.v[i]==INT64_MIN||g.vu[i]==INT64_MAX)return -1;
  c[j]=ck((W(lo[i])+g.au[i])/2);d[j]=ck((W(g.v[i])+g.vu[i])/2);
  E+=std::max(W(c[j])-lo[i],W(g.au[i])-c[j]);F+=std::max(W(d[j])-g.v[i],W(g.vu[i])-d[j]);
  if(j)T+=std::abs(ck(W(c[j])-c[j-1]));}
 T+=std::abs(c.back());W M=W(n)*delta-E;if(M<=0||T*n>128*M)return -1;
 int om=std::max(1,int(1+T*n/M)); if(om>128)return -1;
 // w in [0,1]: cos in [1/2,1], using interval lower sign.
 W zero=0;for(auto a:g.v)zero+=(a>=0?W(a):2*W(a));if(zero<=0)return -1;
 for(int k=2;k<2*om;k++){auto const&w=weights(n,k); W cv=-F*P,cc=M*P;
  for(int j=0;j<n;j++){auto [l,r]=w[j];cv+=W(d[j])*(d[j]>=0?P-r:P-l);
    cc-=W(c[j])*(c[j]>=0?P-l:P-r);}
  if(cv<=0&&cc<=0)return -1;
 }
 return om;
}
int run(int argc,char**argv){
 if(argc<9){std::cerr<<"usage xlo xhi ylo yhi initial_bpu nq levels mx [point]\n";return 2;}
 int xl=atoi(argv[1]),xh=atoi(argv[2]),yl=atoi(argv[3]),yh=atoi(argv[4]),b=atoi(argv[5]),n=atoi(argv[6]),levels=atoi(argv[7]),mx=atoi(argv[8]);bool point=argc>9&&std::string(argv[9])=="point";
 if(xl<3||xh>24||yl<3||yh>24||xl>xh||yl>yh||b<1||b>1024||(b&(b-1))||n<8||(n&(n-1))||mx<2||mx>512||(mx&(mx-1)))return 2;
 std::vector<Box> jobs; for(int i=xl*b;i<(point?xl*b+1:xh*b);++i)for(int j=std::max(yl*b,i);j<(point?yl*b+1:yh*b);++j)jobs.push_back({i,j});
 if(argc>9&&std::string(argv[9])=="cell")jobs={{argc>11?atoi(argv[10]):xl*b,argc>11?atoi(argv[11]):yl*b}};
 for(int l=0;l<levels;l++,b*=2,n=std::min(1024,n*2)){
  if(n>4096||b>2048)throw std::range_error("grid maximum");
  std::vector<Cur> curves(jobs.size());
  for(int k=0;k<n;k++){
   R ql(D(k)/n,D(k+1)/n);J q(ql,R(1),R(0));
   Cache Clo{{},b,q,false},Chi{{},b,q,true};
   for(size_t a=0;a<jobs.size();a++){
    I64 v=INT64_MIN,vs=INT64_MIN,ve=INT64_MIN,vse=INT64_MIN;
    for(int hp=0;hp<=1;hp++){
     bool hi=hp==1; auto job=jobs[a];double tx=ql.a*job.i/b,ty=ql.b*(job.j+1)/b;
     if((hi&&(ql.b>.5||(ty>4.25&&tx>=1.5)))||(!hi&&ql.a<=.5&&tx<.65&&ty<2.))continue;
     try{J res=block(hi?Chi:Clo,job.i,job.j,point,mx);v=std::max(v,integerlo(res.e.a*2));ve=std::max(ve,integerlo(-res.e.b*2));vs=std::max(vs,integerlo(res.t.a));vse=std::max(vse,integerlo(-res.t.b));if(jobs.size()==1&&std::getenv("TRACE")&&k%(std::max(n/16,1))==0)std::cerr<<k<<" chart="<<hi<<" curv="<<2*res.e.a<<","<<2*res.e.b<<"\n";}
     catch(std::range_error const&){/* insufficiently sharp enclosure is not counted */}
    }
    curves[a].a.push_back(v);curves[a].v.push_back(std::max<I64>(0,vs));curves[a].au.push_back(ve==INT64_MIN?INT64_MAX:-ve);curves[a].vu.push_back(vse==INT64_MIN?INT64_MAX:-vse);
   }
  }
  std::vector<Box> todo;int cut=0,fail=0,maxI=-1,maxJ=-1;Box worst{};
  for(size_t t=0;t<jobs.size();t++){auto &v=curves[t].a;auto [i,j]=jobs[t];int val=freq(v,0);if(val<0)val=both_freq(curves[t],integerlo(d0(i,j,b,point).a));
    if(point||jobs.size()==1) {D mass=0,tv=0;for(size_t k=0;k<v.size();k++){mass+=std::ldexp(D(v[k]),-28)/n;if(k)tv+=std::abs(std::ldexp(D(v[k]),-28)-std::ldexp(D(v[k-1]),-28));}std::cout<<" mass="<<mass<<" tv="<<tv<<" minimum="<<std::ldexp(D(*std::min_element(v.begin(),v.end())),-28)<<"\n";}
    if(val>=0){cut=std::max(cut,val);continue;}++fail;worst={i,j};maxI=std::max(maxI,i);maxJ=std::max(maxJ,j);
    if(point){todo.push_back({2*i,2*j});continue;}
    for(int a=0;a<2;a++)for(int z=0;z<2;z++)if(2*i+a<=2*j+z)todo.push_back({2*i+a,2*j+z});
  }
  std::cout<<"bpu="<<b<<" bins="<<n<<" boxes="<<jobs.size()<<" failed="<<fail<<" last_fail="<<worst.i<<","<<worst.j<<" unresolved_largest_i_j="<<maxI<<","<<maxJ<<" max_frequency_cut="<<cut<<std::endl;
  if(!fail){std::cout<<"PASS continuum signed curvature on entire input"<<std::endl;return 0;}
  jobs=std::move(todo);
 }
 std::cout<<"INCOMPLETE"<<std::endl;return 1;
}
}
int main(int argc,char**argv){if(std::fegetround()!=FE_TONEAREST)throw std::runtime_error("round-to-nearest mandatory");return nearcert::run(argc,argv);}
