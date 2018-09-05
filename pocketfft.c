#include <math.h>
#include <string.h>
#include "pocketfft.h"
#include "c_utils.h"
#include "trig_utils.h"

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

static double cost_guess (size_t n)
  {
  const double lfp=1.1; // penalty for non-hardcoced larger factors
  size_t ni=n;
  double result=0.;
  size_t tmp;
  while (((tmp=(n>>1))<<1)==n)
    { result+=2; n=tmp; }

  size_t limit=(size_t)sqrt(n+0.01);
  for (size_t x=3; x<=limit; x+=2)
  while ((tmp=(n/x))*x==n)
    {
    result+= (x<=5) ? x : lfp*x; // penalize larger prime factors
    n=tmp;
    limit=(size_t)sqrt(n+0.01);
    }
  if (n>1) result+=(n<=5) ? n : lfp*n;

  return result*ni;
  }

/* returns the smallest composite of 2, 3 and 5 which is >= n */
static size_t good_size(size_t n)
  {
  if (n<=6) return n;

  size_t bestfac=2*n;
  for (size_t f2=1; f2<bestfac; f2*=2)
    for (size_t f23=f2; f23<bestfac; f23*=3)
      for (size_t f235=f23; f235<bestfac; f235*=5)
        if (f235>=n) bestfac=f235;
  return bestfac;
  }

typedef struct cmplx {
  double r,i;
} cmplx;

#define NFCT 20
typedef struct cfftp_fctdata
  {
  size_t fct;
  cmplx *tw;
  } cfftp_fctdata;

typedef struct cfftp_plan_i
  {
  size_t length, nfct;
  cmplx *mem;
  cmplx *wrk;
  cfftp_fctdata fct[NFCT];
  } cfftp_plan_i;
typedef struct cfftp_plan_i * cfftp_plan;

#define CONCAT(a,b) a ## b

#define X(arg) CONCAT(passb,arg)
#define BACKWARD
#include "pocketfft.inc"
#undef BACKWARD
#undef X
#define X(arg) CONCAT(passf,arg)
#include "pocketfft.inc"
#undef X

static void cfftp_forward(cfftp_plan plan, double c[])
  { passf_all(plan,(cmplx *)c); }

static void cfftp_backward(cfftp_plan plan, double c[])
  { passb_all(plan,(cmplx *)c); }

static void cfftp_factorize (cfftp_plan plan)
  {
  size_t length=plan->length;
  size_t nfct=0;
  while ((length%4)==0)
    { plan->fct[nfct++].fct=4; length>>=2; }
  if ((length%2)==0)
    {
    length>>=1;
    // factor 2 should be at the front of the factor list
    plan->fct[nfct++].fct=2;
    SWAP(plan->fct[0].fct, plan->fct[nfct-1].fct,size_t);
    }
  size_t maxl=(size_t)(sqrt((double)length))+1;
  for (size_t divisor=3; (length>1)&&(divisor<maxl); divisor+=2)
    if ((length%divisor)==0)
      {
      while ((length%divisor)==0)
        { plan->fct[nfct++].fct=divisor; length/=divisor; }
      maxl=(size_t)(sqrt((double)length))+1;
      }
  if (length>1) plan->fct[nfct++].fct=length;
  plan->nfct=nfct;
  }

static size_t cfftp_twsize (cfftp_plan plan)
  {
  size_t twsize=0, l1=1;
  for (size_t k=0; k<plan->nfct; ++k)
    {
    size_t ip=plan->fct[k].fct, ido= plan->length/(l1*ip);
    twsize+=(ip-1)*ido;
    l1*=ip;
    }
  return twsize;
  }

static void cfftp_comp_twiddle (cfftp_plan plan)
  {
  size_t length=plan->length;
  triggen tg;
  triggen_init(&tg,length);
  size_t l1=1;
  size_t memofs=0;
  for (size_t k=0; k<plan->nfct; ++k)
    {
    size_t ip=plan->fct[k].fct, ido= length/(l1*ip);
    plan->fct[k].tw=plan->mem+memofs;
    memofs+=(ip-1)*ido;
    for (size_t j=1; j<ip; ++j)
      {
      for (size_t i=1; i<ido; ++i)
        triggen_get(&tg,j*l1*i,&(plan->fct[k].tw[(j-1)*ido+i].i),
                               &(plan->fct[k].tw[(j-1)*ido+i].r));
      if (ip>6)
        triggen_get(&tg,j*l1*ido,&(plan->fct[k].tw[(j-1)*ido].i),
                                 &(plan->fct[k].tw[(j-1)*ido].r));
      }
    l1*=ip;
    }
  triggen_destroy(&tg);
  }

static cfftp_plan make_cfftp_plan (size_t length)
  {
  UTIL_ASSERT(length!=0,"bad FFT length");
  cfftp_plan plan = RALLOC(cfftp_plan_i,1);
  plan->length=length;
  plan->nfct=0;
  for (size_t i=0; i<NFCT; ++i)
    plan->fct[i]=(cfftp_fctdata){0,0};
  plan->mem=0;
  if (length==1) return plan;
  cfftp_factorize (plan);
  size_t tws=cfftp_twsize(plan);
  plan->mem=RALLOC(cmplx,tws+length);
  cfftp_comp_twiddle(plan);
  plan->wrk=plan->mem+tws;
  return plan;
  }

static void destroy_cfftp_plan (cfftp_plan plan)
  {
  DEALLOC(plan->mem);
  DEALLOC(plan);
  }

typedef struct rfftp_fctdata
  {
  size_t fct;
  double *tw, *tws;
  } rfftp_fctdata;

typedef struct rfftp_plan_i
  {
  size_t length, nfct;
  double *mem, *wrk;
  rfftp_fctdata fct[NFCT];
  } rfftp_plan_i;
typedef struct rfftp_plan_i * rfftp_plan;

#define WA(x,i) wa[(i)+(x)*ido]
#define PM(a,b,c,d) { a=c+d; b=c-d; }
/* (a+ib) = conj(c+id) * (e+if) */
#define MULPM(a,b,c,d,e,f) { a=c*e+d*f; b=c*f-d*e; }

#define CC(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define CH(a,b,c) ch[(a)+ido*((b)+cdim*(c))]

NOINLINE static void radf2 (size_t ido, size_t l1, const double * restrict cc,
  double * restrict ch, const double * restrict wa)
  {
  const size_t cdim=2;

  for (size_t k=0; k<l1; k++)
    PM (CH(0,0,k),CH(ido-1,1,k),CC(0,k,0),CC(0,k,1))
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      CH(    0,1,k) = -CC(ido-1,k,1);
      CH(ido-1,0,k) =  CC(ido-1,k,0);
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      double tr2, ti2;
      MULPM (tr2,ti2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1))
      PM (CH(i-1,0,k),CH(ic-1,1,k),CC(i-1,k,0),tr2)
      PM (CH(i  ,0,k),CH(ic  ,1,k),ti2,CC(i  ,k,0))
      }
  }

NOINLINE static void radf3(size_t ido, size_t l1, const double * restrict cc,
  double * restrict ch, const double * restrict wa)
  {
  const size_t cdim=3;
  static const double taur=-0.5, taui=0.86602540378443864676;

  for (size_t k=0; k<l1; k++)
    {
    double cr2=CC(0,k,1)+CC(0,k,2);
    CH(0,0,k) = CC(0,k,0)+cr2;
    CH(0,2,k) = taui*(CC(0,k,2)-CC(0,k,1));
    CH(ido-1,1,k) = CC(0,k,0)+taur*cr2;
    }
  if (ido==1) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      double di2, di3, dr2, dr3;
      MULPM (dr2,di2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1)) // d2=conj(WA0)*CC1
      MULPM (dr3,di3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2)) // d3=conj(WA1)*CC2
      double cr2=dr2+dr3; // c add
      double ci2=di2+di3;
      CH(i-1,0,k) = CC(i-1,k,0)+cr2; // c add
      CH(i  ,0,k) = CC(i  ,k,0)+ci2;
      double tr2 = CC(i-1,k,0)+taur*cr2; // c add
      double ti2 = CC(i  ,k,0)+taur*ci2;
      double tr3 = taui*(di2-di3);  // t3 = taui*i*(d3-d2)?
      double ti3 = taui*(dr3-dr2);
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr2,tr3) // PM(i) = t2+t3
      PM(CH(i  ,2,k),CH(ic  ,1,k),ti3,ti2) // PM(ic) = conj(t2-t3)
      }
  }

NOINLINE static void radf4(size_t ido, size_t l1, const double * restrict cc,
  double * restrict ch, const double * restrict wa)
  {
  const size_t cdim=4;
  static const double hsqt2=0.70710678118654752440;

  for (size_t k=0; k<l1; k++)
    {
    double tr1,tr2;
    PM (tr1,CH(0,2,k),CC(0,k,3),CC(0,k,1))
    PM (tr2,CH(ido-1,1,k),CC(0,k,0),CC(0,k,2))
    PM (CH(0,0,k),CH(ido-1,3,k),tr2,tr1)
    }
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      double ti1=-hsqt2*(CC(ido-1,k,1)+CC(ido-1,k,3));
      double tr1= hsqt2*(CC(ido-1,k,1)-CC(ido-1,k,3));
      PM (CH(ido-1,0,k),CH(ido-1,2,k),CC(ido-1,k,0),tr1)
      PM (CH(    0,3,k),CH(    0,1,k),ti1,CC(ido-1,k,2))
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
      MULPM(cr2,ci2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1))
      MULPM(cr3,ci3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2))
      MULPM(cr4,ci4,WA(2,i-2),WA(2,i-1),CC(i-1,k,3),CC(i,k,3))
      PM(tr1,tr4,cr4,cr2)
      PM(ti1,ti4,ci2,ci4)
      PM(tr2,tr3,CC(i-1,k,0),cr3)
      PM(ti2,ti3,CC(i  ,k,0),ci3)
      PM(CH(i-1,0,k),CH(ic-1,3,k),tr2,tr1)
      PM(CH(i  ,0,k),CH(ic  ,3,k),ti1,ti2)
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr3,ti4)
      PM(CH(i  ,2,k),CH(ic  ,1,k),tr4,ti3)
      }
  }

NOINLINE static void radf5(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=5;
  double tr11= 0.3090169943749474241, ti11=0.95105651629515357212,
     tr12=-0.8090169943749474241, ti12=0.58778525229247312917;

  for (size_t k=0; k<l1; k++)
    {
    double cr2, cr3, ci4, ci5;
    PM (cr2,ci5,CC(0,k,4),CC(0,k,1))
    PM (cr3,ci4,CC(0,k,3),CC(0,k,2))
    CH(0,0,k)=CC(0,k,0)+cr2+cr3;
    CH(ido-1,1,k)=CC(0,k,0)+tr11*cr2+tr12*cr3;
    CH(0,2,k)=ti11*ci5+ti12*ci4;
    CH(ido-1,3,k)=CC(0,k,0)+tr12*cr2+tr11*cr3;
    CH(0,4,k)=ti12*ci5-ti11*ci4;
    }
  if (ido==1) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2; i<ido; i+=2)
      {
      double ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3,
         dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
      size_t ic=ido-i;
      MULPM (dr2,di2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1))
      MULPM (dr3,di3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2))
      MULPM (dr4,di4,WA(2,i-2),WA(2,i-1),CC(i-1,k,3),CC(i,k,3))
      MULPM (dr5,di5,WA(3,i-2),WA(3,i-1),CC(i-1,k,4),CC(i,k,4))
      PM(cr2,ci5,dr5,dr2)
      PM(ci2,cr5,di2,di5)
      PM(cr3,ci4,dr4,dr3)
      PM(ci3,cr4,di3,di4)
      CH(i-1,0,k)=CC(i-1,k,0)+cr2+cr3;
      CH(i  ,0,k)=CC(i  ,k,0)+ci2+ci3;
      tr2=CC(i-1,k,0)+tr11*cr2+tr12*cr3;
      ti2=CC(i  ,k,0)+tr11*ci2+tr12*ci3;
      tr3=CC(i-1,k,0)+tr12*cr2+tr11*cr3;
      ti3=CC(i  ,k,0)+tr12*ci2+tr11*ci3;
      MULPM(tr5,tr4,cr5,cr4,ti11,ti12)
      MULPM(ti5,ti4,ci5,ci4,ti11,ti12)
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr2,tr5)
      PM(CH(i  ,2,k),CH(ic  ,1,k),ti5,ti2)
      PM(CH(i-1,4,k),CH(ic-1,3,k),tr3,tr4)
      PM(CH(i  ,4,k),CH(ic  ,3,k),ti4,ti3)
      }
  }

NOINLINE static void radfg(size_t ido, size_t ip, size_t l1,
  double *cc, double *ch, const double *wa, const double *csarr)
  {
  const size_t cdim=ip;
  size_t ipph=(ip+1)/ 2;

  double *tarr=RALLOC(double,4*ipph);
  double *tarr2=tarr+2*ipph;
  for(size_t k=0; k<l1; k++)
    {
    {
    double v0=tarr[0]=CC(0,k,0);
    for (size_t j=1; j<ipph; ++j)
      {
      double t10=CC(0,k,ip-j), t20=CC(0,k,j);
      v0+=tarr[2*j]=t10+t20;
      tarr[2*j+1]=t10-t20;
      }
    CH(0,0,k)=v0;
    for(size_t l=1; l<ipph; l++)
      {
      size_t aidx=2*l;
      double tx1=tarr[0]+csarr[aidx]*tarr[2];
      double tx2=csarr[aidx+1]*tarr[3];
      for(size_t j=2; j<ipph; j++)
        {
        aidx+=2*l;
        if (aidx>=2*ip) aidx-=2*ip;
        tx1+=csarr[aidx  ]*tarr[2*j];
        tx2+=csarr[aidx+1]*tarr[2*j+1];
        }
      CH(ido-1,2*l-1,k) = tx1;
      CH(0    ,2*l  ,k) = tx2;
      }
    }
    for(size_t i=2; i<ido; i+=2)
      {
      double v0a=tarr[0]=CC(i-1,k,0);
      double v0b=tarr2[0]=CC(i  ,k,0);
      for (size_t j=1; j<ipph; ++j)
        {
        size_t idij=(j-1)*ido+1+i-2, idij2=(ip-j-1)*ido+1+i-2;
        double t1=CC(i-1,k,j   ), t2=CC(i,k,j   );
        double t3=CC(i-1,k,ip-j), t4=CC(i,k,ip-j);
        double t5, t6, t7, t8;
        MULPM(t5,t6,wa[idij -1],wa[idij ],t1,t2)
        MULPM(t7,t8,wa[idij2-1],wa[idij2],t3,t4)
        PM(t1,tarr2[2*j+1],t7,t5)
        PM(t3,tarr[2*j+1],t6,t8)
        v0a+=tarr[2*j]=t1;
        v0b+=tarr2[2*j]=t3;
        }
      CH(i-1,0,k)=v0a;
      CH(i,0,k)=v0b;
      for(size_t l=1; l<ipph; l++)
        {
        size_t aidx=2*l;
        double tx1=tarr[0]+csarr[aidx]*tarr[2];
        double tx2=csarr[aidx+1]*tarr[3];
        double tx3=tarr2[0]+csarr[aidx]*tarr2[2];
        double tx4=csarr[aidx+1]*tarr2[3];
        for(size_t j=2; j<ipph; j++)
          {
          aidx+=2*l;
          if (aidx>=2*ip) aidx-=2*ip;
          tx1+=csarr[aidx  ]*tarr[2*j];
          tx2+=csarr[aidx+1]*tarr[2*j+1];
          tx3+=csarr[aidx  ]*tarr2[2*j];
          tx4+=csarr[aidx+1]*tarr2[2*j+1];
          }
        PM (CH(i-1,2*l,k),CH(ido-i-1,2*l-1,k),tx1,tx2)
        PM (CH(i  ,2*l,k),CH(ido-i  ,2*l-1,k),tx4,tx3)
        }
      }
    }
  DEALLOC(tarr);
  }

#undef CH
#undef CC
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]

NOINLINE static void radb2(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=2;

  for (size_t k=0; k<l1; k++)
    PM (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(ido-1,1,k))
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      CH(ido-1,k,0) = 2.*CC(ido-1,0,k);
      CH(ido-1,k,1) =-2.*CC(0    ,1,k);
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      double ti2, tr2;
      PM (CH(i-1,k,0),tr2,CC(i-1,0,k),CC(ic-1,1,k))
      PM (ti2,CH(i  ,k,0),CC(i  ,0,k),CC(ic  ,1,k))
      MULPM (CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),ti2,tr2)
      }
  }

NOINLINE static void radb3(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=3;
  double taur=-0.5, taui=0.86602540378443864676;

  for (size_t k=0; k<l1; k++)
    {
    double tr2=2.*CC(ido-1,1,k);
    double cr2=CC(0,0,k)+taur*tr2;
    CH(0,k,0)=CC(0,0,k)+tr2;
    double ci3=2.*taui*CC(0,2,k);
    PM (CH(0,k,2),CH(0,k,1),cr2,ci3);
    }
  if (ido==1) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      double tr2=CC(i-1,2,k)+CC(ic-1,1,k); // t2=CC(I) + conj(CC(ic))
      double ti2=CC(i  ,2,k)-CC(ic  ,1,k);
      double cr2=CC(i-1,0,k)+taur*tr2;     // c2=CC +taur*t2
      double ci2=CC(i  ,0,k)+taur*ti2;
      CH(i-1,k,0)=CC(i-1,0,k)+tr2;         // CH=CC+t2
      CH(i  ,k,0)=CC(i  ,0,k)+ti2;
      double cr3=taui*(CC(i-1,2,k)-CC(ic-1,1,k));// c3=taui*(CC(i)-conj(CC(ic)))
      double ci3=taui*(CC(i  ,2,k)+CC(ic  ,1,k));
      double di2, di3, dr2, dr3;
      PM(dr3,dr2,cr2,ci3) // d2= (cr2-ci3, ci2+cr3) = c2+i*c3
      PM(di2,di3,ci2,cr3) // d3= (cr2+ci3, ci2-cr3) = c2-i*c3
      MULPM(CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),di2,dr2) // ch = WA*d2
      MULPM(CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),di3,dr3)
      }
  }

NOINLINE static void radb4(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=4;
  static const double sqrt2=1.41421356237309504880;

  for (size_t k=0; k<l1; k++)
    {
    double tr1, tr2;
    PM (tr2,tr1,CC(0,0,k),CC(ido-1,3,k))
    double tr3=2.*CC(ido-1,1,k);
    double tr4=2.*CC(0,2,k);
    PM (CH(0,k,0),CH(0,k,2),tr2,tr3)
    PM (CH(0,k,3),CH(0,k,1),tr1,tr4)
    }
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      double tr1,tr2,ti1,ti2;
      PM (ti1,ti2,CC(0    ,3,k),CC(0    ,1,k))
      PM (tr2,tr1,CC(ido-1,0,k),CC(ido-1,2,k))
      CH(ido-1,k,0)=tr2+tr2;
      CH(ido-1,k,1)=sqrt2*(tr1-ti1);
      CH(ido-1,k,2)=ti2+ti2;
      CH(ido-1,k,3)=-sqrt2*(tr1+ti1);
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2; i<ido; i+=2)
      {
      double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
      size_t ic=ido-i;
      PM (tr2,tr1,CC(i-1,0,k),CC(ic-1,3,k))
      PM (ti1,ti2,CC(i  ,0,k),CC(ic  ,3,k))
      PM (tr4,ti3,CC(i  ,2,k),CC(ic  ,1,k))
      PM (tr3,ti4,CC(i-1,2,k),CC(ic-1,1,k))
      PM (CH(i-1,k,0),cr3,tr2,tr3)
      PM (CH(i  ,k,0),ci3,ti2,ti3)
      PM (cr4,cr2,tr1,tr4)
      PM (ci2,ci4,ti1,ti4)
      MULPM (CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),ci2,cr2)
      MULPM (CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),ci3,cr3)
      MULPM (CH(i,k,3),CH(i-1,k,3),WA(2,i-2),WA(2,i-1),ci4,cr4)
      }
  }

NOINLINE static void radb5(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=5;
  double tr11= 0.3090169943749474241, ti11=0.95105651629515357212,
     tr12=-0.8090169943749474241, ti12=0.58778525229247312917;

  for (size_t k=0; k<l1; k++)
    {
    double ti5=CC(0,2,k)+CC(0,2,k);
    double ti4=CC(0,4,k)+CC(0,4,k);
    double tr2=CC(ido-1,1,k)+CC(ido-1,1,k);
    double tr3=CC(ido-1,3,k)+CC(ido-1,3,k);
    CH(0,k,0)=CC(0,0,k)+tr2+tr3;
    double cr2=CC(0,0,k)+tr11*tr2+tr12*tr3;
    double cr3=CC(0,0,k)+tr12*tr2+tr11*tr3;
    double ci4, ci5;
    MULPM(ci5,ci4,ti5,ti4,ti11,ti12)
    PM(CH(0,k,4),CH(0,k,1),cr2,ci5)
    PM(CH(0,k,3),CH(0,k,2),cr3,ci4)
    }
  if (ido==1) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      double tr2, tr3, tr4, tr5, ti2, ti3, ti4, ti5;
      PM(tr2,tr5,CC(i-1,2,k),CC(ic-1,1,k))
      PM(ti5,ti2,CC(i  ,2,k),CC(ic  ,1,k))
      PM(tr3,tr4,CC(i-1,4,k),CC(ic-1,3,k))
      PM(ti4,ti3,CC(i  ,4,k),CC(ic  ,3,k))
      CH(i-1,k,0)=CC(i-1,0,k)+tr2+tr3;
      CH(i  ,k,0)=CC(i  ,0,k)+ti2+ti3;
      double cr2=CC(i-1,0,k)+tr11*tr2+tr12*tr3;
      double ci2=CC(i  ,0,k)+tr11*ti2+tr12*ti3;
      double cr3=CC(i-1,0,k)+tr12*tr2+tr11*tr3;
      double ci3=CC(i  ,0,k)+tr12*ti2+tr11*ti3;
      double ci4, ci5, cr5, cr4;
      MULPM(cr5,cr4,tr5,tr4,ti11,ti12)
      MULPM(ci5,ci4,ti5,ti4,ti11,ti12)
      double dr2, dr3, dr4, dr5, di2, di3, di4, di5;
      PM(dr4,dr3,cr3,ci4)
      PM(di3,di4,ci3,cr4)
      PM(dr5,dr2,cr2,ci5)
      PM(di2,di5,ci2,cr5)
      MULPM(CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),di2,dr2)
      MULPM(CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),di3,dr3)
      MULPM(CH(i,k,3),CH(i-1,k,3),WA(2,i-2),WA(2,i-1),di4,dr4)
      MULPM(CH(i,k,4),CH(i-1,k,4),WA(3,i-2),WA(3,i-1),di5,dr5)
      }
  }

NOINLINE static void radbg(size_t ido, size_t ip, size_t l1,
  double *cc, double *ch, const double *wa, const double *csarr)
  {
  const size_t cdim=ip;
  size_t ipph=(ip+1)/ 2;

  double *tarr=RALLOC(double,4*ipph);
  double *tarr2=tarr+2*ipph;
  for(size_t k=0; k<l1; k++)
    {
    {
    double v0=tarr[0]=CC(0,0,k);
    for (size_t j=1; j<ipph; ++j)
      {
      v0+=tarr[2*j]=2.*CC(ido-1,2*j-1,k);
      tarr[2*j+1]=2.*CC(0,2*j,k);
      }
    CH(0,k,0)=v0;
    for(size_t l=1; l<ipph; l++)
      {
      size_t aidx=2*l;
      double tx1=tarr[0]+csarr[aidx]*tarr[2];
      double tx2=csarr[aidx+1]*tarr[3];
      for(size_t j=2; j<ipph; j++)
        {
        aidx+=2*l;
        if (aidx>=2*ip) aidx-=2*ip;
        tx1+=csarr[aidx  ]*tarr[2*j];
        tx2+=csarr[aidx+1]*tarr[2*j+1];
        }
      PM (CH(0,k,ip-l),CH(0,k,l),tx1,tx2)
      }
    }

    for(size_t i=2; i<ido; i+=2)
      {
      double v0a=tarr[0]=CC(i-1,0,k);
      double v0b=tarr2[0]=CC(i,0,k);
      for (size_t j=1; j<ipph; ++j)
        {
        double tx1,tx2,tx3,tx4;
        PM (tx1,tx2,CC(i-1,2*j,k),CC(ido-i-1,2*j-1,k))
        PM (tx4,tx3,CC(i  ,2*j,k),CC(ido-i  ,2*j-1,k))
        v0a+=tarr[2*j]=tx1;
        v0b+=tarr2[2*j]=tx3;
        tarr[2*j+1]=tx2;
        tarr2[2*j+1]=tx4;
        }
      CH(i-1,k,0)=v0a;
      CH(i,k,0)=v0b;
      for(size_t l=1; l<ipph; l++)
        {
        size_t aidx=2*l;
        double txn1=tarr[0]+csarr[aidx]*tarr[2];
        double txn2=csarr[aidx+1]*tarr[3];
        double txn3=tarr2[0]+csarr[aidx]*tarr2[2];
        double txn4=csarr[aidx+1]*tarr2[3];
        for(size_t j=2; j<ipph; j++)
          {
          aidx+=2*l;
          if (aidx>=2*ip) aidx-=2*ip;
          txn1+=csarr[aidx  ]*tarr[2*j];
          txn2+=csarr[aidx+1]*tarr[2*j+1];
          txn3+=csarr[aidx  ]*tarr2[2*j];
          txn4+=csarr[aidx+1]*tarr2[2*j+1];
          }
        double t5, t6, t7, t8;
        PM (t5,t6,txn1,txn4)
        PM (t7,t8,txn3,txn2)
        size_t idij=(l-1)*ido+1+i-2;
        size_t idij2=((ip-l)-1)*ido+1+i-2;
        MULPM (CH(i,k,l),CH(i-1,k,l),wa[idij-1],wa[idij],t7,t6)
        MULPM (CH(i,k,ip-l),CH(i-1,k,ip-l),wa[idij2-1],wa[idij2],t8,t5)
        }
      }
    }
  DEALLOC(tarr);
  }

#undef CC
#undef CH
#undef PM
#undef MULPM
#undef WA

static void rfftp_forward(rfftp_plan plan, double c[])
  {
  if (plan->length==1) return;
  size_t n=plan->length;
  size_t l1=n, nf=plan->nfct;
  double *ch=plan->wrk;
  double *p1=c, *p2=ch;

  for(size_t k1=0; k1<nf;++k1)
    {
    size_t k=nf-k1-1;
    size_t ip=plan->fct[k].fct;
    size_t ido=n / l1;
    l1 /= ip;
    if(ip==4)
      radf4(ido, l1, p1, p2, plan->fct[k].tw);
    else if(ip==2)
      radf2(ido, l1, p1, p2, plan->fct[k].tw);
    else if(ip==3)
      radf3(ido, l1, p1, p2, plan->fct[k].tw);
    else if(ip==5)
      radf5(ido, l1, p1, p2, plan->fct[k].tw);
    else
      radfg(ido, ip, l1, p1, p2, plan->fct[k].tw, plan->fct[k].tws);
    SWAP (p1,p2,double *);
    }
  if (p1!=c)
    memcpy (c,ch,n*sizeof(double));
  }

static void rfftp_backward(rfftp_plan plan, double c[])
  {
  if (plan->length==1) return;
  size_t n=plan->length;
  size_t l1=1, nf=plan->nfct;
  double *ch=plan->wrk;
  double *p1=c, *p2=ch;

  for(size_t k=0; k<nf; k++)
    {
    size_t ip = plan->fct[k].fct,
           ido= n/(ip*l1);
    if(ip==4)
      radb4(ido, l1, p1, p2, plan->fct[k].tw);
    else if(ip==2)
      radb2(ido, l1, p1, p2, plan->fct[k].tw);
    else if(ip==3)
      radb3(ido, l1, p1, p2, plan->fct[k].tw);
    else if(ip==5)
      radb5(ido, l1, p1, p2, plan->fct[k].tw);
    else
      radbg(ido, ip, l1, p1, p2, plan->fct[k].tw, plan->fct[k].tws);
    SWAP (p1,p2,double *);
    l1*=ip;
    }
  if (p1!=c)
    memcpy (c,ch,n*sizeof(double));
  }

static void rfftp_factorize (rfftp_plan plan)
  {
  size_t length=plan->length;
  size_t nfct=0;
  while ((length%4)==0)
    { plan->fct[nfct++].fct=4; length>>=2; }
  if ((length%2)==0)
    {
    length>>=1;
    // factor 2 should be at the front of the factor list
    plan->fct[nfct++].fct=2;
    SWAP(plan->fct[0].fct, plan->fct[nfct-1].fct,size_t);
    }
  size_t maxl=(size_t)(sqrt((double)length))+1;
  for (size_t divisor=3; (length>1)&&(divisor<maxl); divisor+=2)
    if ((length%divisor)==0)
      {
      while ((length%divisor)==0)
        { plan->fct[nfct++].fct=divisor; length/=divisor; }
      maxl=(size_t)(sqrt((double)length))+1;
      }
  if (length>1) plan->fct[nfct++].fct=length;
  plan->nfct=nfct;
  }

static size_t rfftp_comp_twsize(rfftp_plan plan)
  {
  size_t twsize=0, l1=1;
  for (size_t k=0; k<plan->nfct; ++k)
    {
    size_t ip=plan->fct[k].fct, ido= plan->length/(l1*ip);
    twsize+=(ip-1)*ido;
    if (ip>5) twsize+=2*ip;
    l1*=ip;
    }
  return twsize;
  }

static void rfftp_comp_twiddle (rfftp_plan plan)
  {
  size_t length=plan->length;
  triggen tg;
  triggen_init(&tg,length);
  size_t l1=1;
  double *ptr=plan->mem;
  for (size_t k=0; k<plan->nfct; ++k)
    {
    size_t ip=plan->fct[k].fct, ido= length/(l1*ip);
    if (k<plan->nfct-1) // last factor doesn't need twiddles
      {
      plan->fct[k].tw=ptr; ptr+=(ip-1)*ido;
      for (size_t j=1; j<ip; ++j)
        for (size_t i=1; i<=(ido-1)/2; ++i)
          triggen_get(&tg,j*l1*i,&(plan->fct[k].tw[(j-1)*ido+2*i-1]),
                                 &(plan->fct[k].tw[(j-1)*ido+2*i-2]));
      }
    if (ip>5) // special factors required by *g functions
      {
      plan->fct[k].tws=ptr; ptr+=2*ip;
      for (size_t i=0; i<ip; ++i)
        triggen_get(&tg,i*(length/ip),
          &(plan->fct[k].tws[2*i+1]),&(plan->fct[k].tws[2*i]));
      }
    l1*=ip;
    }
  triggen_destroy(&tg);
  }

static rfftp_plan make_rfftp_plan (size_t length)
  {
  UTIL_ASSERT(length!=0,"bad FFT length");
  rfftp_plan plan = RALLOC(rfftp_plan_i,1);
  plan->length=length;
  plan->nfct=0;
  plan->mem=NULL;
  for (size_t i=0; i<NFCT; ++i)
    plan->fct[i]=(rfftp_fctdata){0,0,0};
  if (length==1) return plan;
  rfftp_factorize (plan);
  size_t tws=rfftp_comp_twsize(plan);
  plan->mem=RALLOC(double,tws+plan->length);
  plan->wrk=plan->mem+tws;
  rfftp_comp_twiddle(plan);
  return plan;
  }

static void destroy_rfftp_plan (rfftp_plan plan)
  {
  DEALLOC(plan->mem);
  DEALLOC(plan);
  }

typedef struct fftblue_plan_i
  {
  size_t n, n2;
  cfftp_plan plan;
  double *mem;
  double *bk, *bkf, *akf;
  } fftblue_plan_i;
typedef struct fftblue_plan_i * fftblue_plan;

static fftblue_plan make_fftblue_plan (size_t length)
  {
  fftblue_plan plan = RALLOC(fftblue_plan_i,1);
  plan->n = length;
  plan->n2=good_size(plan->n*2-1);
  plan->mem = RALLOC(double, 2*plan->n+4*plan->n2);
  plan->bk  = plan->mem;
  plan->bkf = plan->bk+2*plan->n;
  plan->akf = plan->bkf+2*plan->n2;

/* initialize b_k */
  double *tmp = RALLOC(double,IMAX(4*plan->n,2*plan->n2));
  sincos_2pibyn(2*plan->n,2*plan->n,&tmp[1],&tmp[0],2);
  plan->bk[0] = 1;
  plan->bk[1] = 0;

  size_t coeff=0;
  for (size_t m=1; m<plan->n; ++m)
    {
    coeff+=2*m-1;
    if (coeff>=2*plan->n) coeff-=2*plan->n;
    plan->bk[2*m  ] = tmp[2*coeff  ];
    plan->bk[2*m+1] = tmp[2*coeff+1];
    }

  /* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
  double xn2 = 1./plan->n2;
  plan->bkf[0] = plan->bk[0]*xn2;
  plan->bkf[1] = plan->bk[1]*xn2;
  for (size_t m=2; m<2*plan->n; m+=2)
    {
    plan->bkf[m]   = plan->bkf[2*plan->n2-m]   = plan->bk[m]   *xn2;
    plan->bkf[m+1] = plan->bkf[2*plan->n2-m+1] = plan->bk[m+1] *xn2;
    }
  for (size_t m=2*plan->n;m<=(2*plan->n2-2*plan->n+1);++m)
    plan->bkf[m]=0.;
  plan->plan=make_cfftp_plan(plan->n2);
  cfftp_forward(plan->plan,plan->bkf);
  DEALLOC(tmp);

  return plan;
  }

static void destroy_fftblue_plan (fftblue_plan plan)
  {
  DEALLOC(plan->mem);
  destroy_cfftp_plan(plan->plan);
  DEALLOC(plan);
  }

static void fftblue_fft(fftblue_plan plan, double c[], int isign)
  {
  size_t n=plan->n;
  size_t n2=plan->n2;
  double *bk  = plan->bk;
  double *bkf = plan->bkf;
  double *akf = plan->akf;

/* initialize a_k and FFT it */
  if (isign>0)
    for (size_t m=0; m<2*n; m+=2)
      {
      akf[m]   = c[m]*bk[m]   - c[m+1]*bk[m+1];
      akf[m+1] = c[m]*bk[m+1] + c[m+1]*bk[m];
      }
  else
    for (size_t m=0; m<2*n; m+=2)
      {
      akf[m]   = c[m]*bk[m]   + c[m+1]*bk[m+1];
      akf[m+1] =-c[m]*bk[m+1] + c[m+1]*bk[m];
      }
  for (size_t m=2*n; m<2*n2; ++m)
    akf[m]=0;

  cfftp_forward (plan->plan,akf);

/* do the convolution */
  if (isign>0)
    for (size_t m=0; m<2*n2; m+=2)
      {
      double im = -akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  =  akf[m]*bkf[m]   + akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }
  else
    for (size_t m=0; m<2*n2; m+=2)
      {
      double im = akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  = akf[m]*bkf[m]   - akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }

/* inverse FFT */
  cfftp_backward (plan->plan,akf);

/* multiply by b_k* */
  if (isign>0)
    for (size_t m=0; m<2*n; m+=2)
      {
      c[m]   = bk[m]  *akf[m] - bk[m+1]*akf[m+1];
      c[m+1] = bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  else
    for (size_t m=0; m<2*n; m+=2)
      {
      c[m]   = bk[m]  *akf[m] + bk[m+1]*akf[m+1];
      c[m+1] =-bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  }

static void cfftblue_backward(fftblue_plan plan, double c[])
  { fftblue_fft(plan,c,1); }

static void cfftblue_forward(fftblue_plan plan, double c[])
  { fftblue_fft(plan,c,-1); }

static void rfftblue_backward(fftblue_plan plan, double c[])
  {
  size_t n=plan->n;
  double *tmp = RALLOC(double,2*n);
  tmp[0]=c[0];
  tmp[1]=0.;
  memcpy (tmp+2,c+1, (n-1)*sizeof(double));
  if ((n&1)==0) tmp[n+1]=0.;
  for (size_t m=2; m<n; m+=2)
    {
    tmp[2*n-m]=tmp[m];
    tmp[2*n-m+1]=-tmp[m+1];
    }
  fftblue_fft(plan,tmp,1);
  for (size_t m=0; m<n; ++m)
    c[m] = tmp[2*m];
  DEALLOC(tmp);
  }

static void rfftblue_forward(fftblue_plan plan, double c[])
  {
  size_t n=plan->n;
  double *tmp = RALLOC(double,2*n);
  for (size_t m=0; m<n; ++m)
    {
    tmp[2*m] = c[m];
    tmp[2*m+1] = 0.;
    }
  fftblue_fft(plan,tmp,-1);
  c[0] = tmp[0];
  memcpy (c+1, tmp+2, (n-1)*sizeof(double));
  DEALLOC(tmp);
  }

typedef struct cfft_plan_i
  {
  cfftp_plan packplan;
  fftblue_plan blueplan;
  } cfft_plan_i;

cfft_plan make_cfft_plan (size_t length)
  {
  UTIL_ASSERT(length!=0,"bad FFT length");
  cfft_plan plan = RALLOC(cfft_plan_i,1);
  double comp1 = cost_guess(length);
  double comp2 = 2*cost_guess(good_size(2*length-1));
  comp2*=1.5; /* fudge factor that appears to give good overall performance */
  plan->blueplan=0;
  plan->packplan=0;
  if (comp2<comp1) // use Bluestein
    plan->blueplan=make_fftblue_plan(length);
  else
    plan->packplan=make_cfftp_plan(length);
  return plan;
  }

void destroy_cfft_plan (cfft_plan plan)
  {
  if (plan->blueplan)
    destroy_fftblue_plan(plan->blueplan);
  if (plan->packplan)
    destroy_cfftp_plan(plan->packplan);
  DEALLOC(plan);
  }

void cfft_backward(cfft_plan plan, double c[])
  {
  if (plan->packplan)
    cfftp_backward(plan->packplan,c);
  else if (plan->blueplan)
    cfftblue_backward(plan->blueplan,c);
  }

void cfft_forward(cfft_plan plan, double c[])
  {
  if (plan->packplan)
    cfftp_forward(plan->packplan,c);
  else if (plan->blueplan)
    cfftblue_forward(plan->blueplan,c);
  }

typedef struct rfft_plan_i
  {
  rfftp_plan packplan;
  fftblue_plan blueplan;
  } rfft_plan_i;

rfft_plan make_rfft_plan (size_t length)
  {
  UTIL_ASSERT(length!=0,"bad FFT length");
  rfft_plan plan = RALLOC(rfft_plan_i,1);
  double comp1 = 0.5*cost_guess(length);
  double comp2 = 2*cost_guess(good_size(2*length-1));
  comp2*=1.5; /* fudge factor that appears to give good overall performance */
  plan->blueplan=0;
  plan->packplan=0;
  if (comp2<comp1) // use Bluestein
    plan->blueplan=make_fftblue_plan(length);
  else
    plan->packplan=make_rfftp_plan(length);
  return plan;
  }

void destroy_rfft_plan (rfft_plan plan)
  {
  if (plan->blueplan)
    destroy_fftblue_plan(plan->blueplan);
  if (plan->packplan)
    destroy_rfftp_plan(plan->packplan);
  DEALLOC(plan);
  }

void rfft_backward(rfft_plan plan, double c[])
  {
  if (plan->packplan)
    rfftp_backward(plan->packplan,c);
  else if (plan->blueplan)
    rfftblue_backward(plan->blueplan,c);
  }

void rfft_forward(rfft_plan plan, double c[])
  {
  if (plan->packplan)
    rfftp_forward(plan->packplan,c);
  else if (plan->blueplan)
    rfftblue_forward(plan->blueplan,c);
  }
