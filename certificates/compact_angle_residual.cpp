// Exact outward-interval checker, separated windows of the original finite-κ defect.
// See compact-balanced-proof.md and cutoff-atom-logconcavity.md. Inherits
// binary64 primitives/positive-moment T,J evaluators from short_zero_residual.cpp.
#define main residual_legacy_entry
#include "short_zero_residual.cpp"
#undef main
#include <fstream>
using JJ=array<J,3>;
JJ lowform(I c0,I m0,Rad const&x,Rad const&m,Rad const&y,Phase const&t,I b){
 auto tx=Tjet(x,t.pp,b,t.Pj,t.f),ty=Tjet(y,t.pp,b,t.Pj,t.f),tm=Tjet(m,t.pp,b,t.Pj,t.f);
 auto u=SJ(W(x,t.pp,tx)),v=SJ(W(y,t.pp,ty)),G=tx/(J(I(1))+u);
 J ep(t.pp,I(1));auto co=EJ(cj(-c0,ep));
 J E=EJ(cj(c0,ep))*u*ty/(J(I(1))+v)-tm;
 return {G,E,co};
}
I meet(I a,I b){if(a.h<b.l||b.h<a.l)throw std::runtime_error("disjoint independently valid intervals");return I(max(a.l,b.l),min(a.h,b.h));}
struct Band{I dim,cor;};
int ac_eval(Box B,int n,int cap=256){
 int x0=B.a,x1=B.b,y0=B.c,y1=B.d,b0=B.e,b1=B.g;
 if((x0+y0)%2||(x1+y1)%2)throw std::runtime_error("midpoint");
 I b(BMAX*b0/64,BMAX*b1/64);auto x=radial(x0,x1,b),y=radial(y0,y1,b),m=radial((x0+y0)/2,(x1+y1)/2,b);
 I c=(y.t-x.t)*.5,c2=sq(c);
 constexpr int K=8;array<I,K> l{sq(c),sq(c*.5),sq(c*1.5),sq(c*2),I(1),I(4),I(9),I(16)};
 array<vector<Band>,K> col; for(auto& s:col) s.reserve(n);
 for(auto const&t:makephase(n,b0,b1)){
  JJ low;bool hl=false;J high;bool hh=false;
  if(t.pp.h<=.9375){
   try{low=lowform(c,m.t,x,m,y,t,b);hl=true;}catch(std::domain_error const&){}
  }
  if(t.pp.l>=.5){auto X=series(x.a,t.f),Y=series(y.a,t.f),M=series(m.a,t.f);high=t.Pj*(M-SJ(X*Y));hh=true;}
  if(!hl&&!hh)return 0;
  for(int k=0;k<K;k++){
    I dim,cor;
    if(hl){
     auto &G=low[0],&E=low[1],&co=low[2];
     auto base=co.f*G.f+E.f;dim=l[k]*base;
     // D''-λ²D = exp(-cp)*((c²-λ²)G-2cG'+G'')+E''-λ²E.
     I extra=k==0?I(0):(k>=4?c2-l[k]:(I(1)-(k==1?I(.25):k==2?I(2.25):I(4)))*c2);
     cor=co.f*(extra*G.f-I(2)*c*G.d+I(2)*G.s)+I(2)*E.s-l[k]*E.f;
    }
    if(hh){
      auto dimH=l[k]*high.f,corH=I(2)*high.s-dimH;
      if(hl){dim=meet(dim,dimH);cor=meet(cor,corH);}else{dim=dimH;cor=corH;}
    }
    // Beyond q=kappa actual function is zero; both endpoints in grid are exact dyadics.
    if(t.q.l*t.q.l < b.h){dim=I(min(0.,dim.l),max(0.,dim.h));cor=I(min(0.,cor.l),max(0.,cor.h));}
    col[k].push_back({dim,cor});
  }
 }
 for(auto const&pr:col){
   I dim,cor,easy,tv(std::abs(pr.back().cor.l));
   for(int i=0;i<n;++i){
     dim=dim+pr[i].dim.l;cor=cor+pr[i].cor.l;easy=easy+min(0.,(I(2)*pr[i].cor.l).l);
     if(i>0){I d=I(pr[i].cor.l)-I(pr[i-1].cor.l);tv=tv+max(-d.l,d.h);}
   }
   if((dim+easy).l>0)return 1;
   I mass=dim+cor;if(dim.l<=0||mass.l<=0)continue;
   double om=(I(n)*tv/mass).h;if(om>=cap)continue;
   // existing weights cover [(k/2),(k+1)/2] for ordinary 1-cos; k=0,1 now also allowed.
   bool ok=true;int last=std::max(1,(int)ceil(om));
   for(int k=0;k<2*last;k++){
     auto const& ww=weights(n,k);I s=dim;
     for(int i=0;i<n;i++){double z=pr[i].cor.l;s=s+I(z)*(z>=0?ww[i].l:ww[i].h);}
     if(s.l<=0){ok=false;break;}
   }
   if(ok)return 2;
 }
 return 0;
}
int compact_main(int argc,char**argv){
 if(argc!=8&&argc!=10&&argc!=11){puts("x0 x1 y0 y1 n_initial depth maxvisits [b0 b1] [trace]");return 2;}
 if(std::fegetround()!=FE_TONEAREST)throw std::runtime_error("rounding mode");
 Box root{S(atof(argv[1])),S(atof(argv[2])),S(atof(argv[3])),S(atof(argv[4])),argc>=10?atoi(argv[8]):0,argc>=10?atoi(argv[9]):64};
 int n0=atoi(argv[5]),dep=atoi(argv[6]),mx=atoi(argv[7]);
 if(!(3*SCALE<=root.a&&root.a<=root.b&&root.b<=16*SCALE&&3*SCALE<=root.c&&root.b<=root.c&&root.c<=root.d&&root.d<=24*SCALE&&0<=root.e&&root.e<=root.g&&root.g<=64&&n0>=4&&n0<=1024&&(n0&(n0-1))==0&&dep>=0&&dep<=28&&mx>0))throw std::runtime_error("range");
 std::ofstream o;if(argc==11){o.open(argv[10]);if(!o)throw std::runtime_error("write trace");o<<"radialScale 16384 parameterBScale 16384 roots "<<root.a<<' '<<root.b<<' '<<root.c<<' '<<root.d<<' '<<root.e<<' '<<root.g<<'\n';}
 vector<std::pair<Box,int>> todo{{root,0}};int vis=0,pa=0,meth[3]{};
 while(!todo.empty()&&vis++<mx){
  auto [v,d]=todo.back();todo.pop_back();int ok=0;
  int n=min(1024,n0*(1<<min(5,max(0,(d-5)/4))));
  try{ok=ac_eval(v,n);}catch(std::domain_error const&){ok=0;}
  if(ok){pa++;meth[ok]++;if(o)o<<"L "<<n<<' '<<ok<<'\n';continue;}
  if(d>=dep){printf("INCOMPLETE x=[%.8f,%.8f] y=[%.8f,%.8f] b=[%d,%d]/16384 n=%d visits=%d passed=%d\n",v.a/SCALE,v.b/SCALE,v.c/SCALE,v.d/SCALE,v.e,v.g,n,vis,pa);return 1;}
  auto left=v,right=v;int xw=v.b-v.a,yw=v.d-v.c,bw=v.g-v.e;char a='?';int cut=0;
  if(xw>=yw && xw>=bw*32&&xw>=4){a='x';cut=right.a=left.b=(v.a+v.b)/2;}
  else if(yw>=bw*32&&yw>=4){a='y';cut=right.c=left.d=(v.c+v.d)/2;}
  else if(bw>=2){a='b';cut=right.e=left.g=(v.e+v.g)/2;}
  else if(xw>=4){a='x';cut=right.a=left.b=(v.a+v.b)/2;}
  else if(yw>=4){a='y';cut=right.c=left.d=(v.c+v.d)/2;}
  else{puts("INCOMPLETE resolution");return 1;}
  if(o)o<<"S "<<a<<' '<<cut<<'\n';
  todo.push_back({right,d+1});todo.push_back({left,d+1});
 }
 if(!todo.empty()){puts("INCOMPLETE cap");return 1;}
 printf("CERTIFIED x=[%.5g,%.5g] y=[%.5g,%.5g] b=[%d,%d]/16384 visits=%d leaves=%d sign/phase=%d/%d\n",root.a/SCALE,root.b/SCALE,root.c/SCALE,root.d/SCALE,root.e,root.g,vis,pa,meth[1],meth[2]);return 0;
}

#ifndef COMPACT_READBACK_ONLY
int main(int argc,char**argv){return compact_main(argc,argv);}
#endif
