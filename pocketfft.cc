#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <memory>

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#define WARN_UNUSED_RESULT __attribute__ ((warn_unused_result))
#define restrict __restrict__
#else
#define NOINLINE
#define WARN_UNUSED_RESULT
#define restrict
#endif
#define RALLOC(type,num) \
  ((type *)aligned_alloc(64,(num)*sizeof(type)))
#define DEALLOC(ptr) \
  do { free(ptr); (ptr)=NULL; } while(0)

using namespace std;

namespace {

class OOM {};

template<typename T> struct arr
  {
  private:
    T *p;
    size_t sz;

  public:
    arr() : p(0), sz(0) {}
    arr(size_t n) : p(RALLOC(T, n)), sz(n) {}
    ~arr() { DEALLOC(p); }

    void resize(size_t n)
      {
      if (n==sz) return;
      DEALLOC(p);
      p = RALLOC(T, n);
      sz = n;
      }

    T &operator[](size_t idx) { return p[idx]; }
    const T &operator[](size_t idx) const { return p[idx]; }

    T *data() { return p; }
    const T *data() const { return p; }

    size_t size() const { return sz; }
  };

class sincos_2pibyn
  {
  private:
    arr<double> data;

    // adapted from https://stackoverflow.com/questions/42792939/
    // CAUTION: this function only works for arguments in the range [-0.25; 0.25]!
    void my_sincosm1pi (double a, double *restrict res)
      {
      double s = a * a;
      /* Approximate cos(pi*x)-1 for x in [-0.25,0.25] */
      double r =     -1.0369917389758117e-4;
      r = fma (r, s,  1.9294935641298806e-3);
      r = fma (r, s, -2.5806887942825395e-2);
      r = fma (r, s,  2.3533063028328211e-1);
      r = fma (r, s, -1.3352627688538006e+0);
      r = fma (r, s,  4.0587121264167623e+0);
      r = fma (r, s, -4.9348022005446790e+0);
      double c = r*s;
      /* Approximate sin(pi*x) for x in [-0.25,0.25] */
      r =             4.6151442520157035e-4;
      r = fma (r, s, -7.3700183130883555e-3);
      r = fma (r, s,  8.2145868949323936e-2);
      r = fma (r, s, -5.9926452893214921e-1);
      r = fma (r, s,  2.5501640398732688e+0);
      r = fma (r, s, -5.1677127800499516e+0);
      s = s * a;
      r = r * s;
      s = fma (a, 3.1415926535897931e+0, r);
      res[0] = c;
      res[1] = s;
      }

    NOINLINE void calc_first_octant(size_t den, double * restrict res)
      {
      size_t n = (den+4)>>3;
      if (n==0) return;
      res[0]=1.; res[1]=0.;
      if (n==1) return;
      size_t l1=(size_t)sqrt(n);
      for (size_t i=1; i<l1; ++i)
        my_sincosm1pi((2.*i)/den,&res[2*i]);
      size_t start=l1;
      while(start<n)
        {
        double cs[2];
        my_sincosm1pi((2.*start)/den,cs);
        res[2*start] = cs[0]+1.;
        res[2*start+1] = cs[1];
        size_t end = l1;
        if (start+end>n) end = n-start;
        for (size_t i=1; i<end; ++i)
          {
          double csx[2]={res[2*i], res[2*i+1]};
          res[2*(start+i)] = ((cs[0]*csx[0] - cs[1]*csx[1] + cs[0]) + csx[0]) + 1.;
          res[2*(start+i)+1] = (cs[0]*csx[1] + cs[1]*csx[0]) + cs[1] + csx[1];
          }
        start += l1;
        }
      for (size_t i=1; i<l1; ++i)
        res[2*i] += 1.;
      }

    NOINLINE void calc_first_quadrant(size_t n, double * restrict res)
      {
      double * restrict p = res+n;
      calc_first_octant(n<<1, p);
      size_t ndone=(n+2)>>2;
      size_t i=0, idx1=0, idx2=2*ndone-2;
      for (; i+1<ndone; i+=2, idx1+=2, idx2-=2)
        {
        res[idx1]   = p[2*i];
        res[idx1+1] = p[2*i+1];
        res[idx2]   = p[2*i+3];
        res[idx2+1] = p[2*i+2];
        }
      if (i!=ndone)
        {
        res[idx1  ] = p[2*i];
        res[idx1+1] = p[2*i+1];
        }
      }

    NOINLINE void calc_first_half(size_t n, double * restrict res)
      {
      int ndone=(n+1)>>1;
      double * p = res+n-1;
      calc_first_octant(n<<2, p);
      int i4=0, in=n, i=0;
      for (; i4<=in-i4; ++i, i4+=4) // octant 0
        {
        res[2*i] = p[2*i4]; res[2*i+1] = p[2*i4+1];
        }
      for (; i4-in <= 0; ++i, i4+=4) // octant 1
        {
        int xm = in-i4;
        res[2*i] = p[2*xm+1]; res[2*i+1] = p[2*xm];
        }
      for (; i4<=3*in-i4; ++i, i4+=4) // octant 2
        {
        int xm = i4-in;
        res[2*i] = -p[2*xm+1]; res[2*i+1] = p[2*xm];
        }
      for (; i<ndone; ++i, i4+=4) // octant 3
        {
        int xm = 2*in-i4;
        res[2*i] = -p[2*xm]; res[2*i+1] = p[2*xm+1];
        }
      }

    NOINLINE void fill_first_quadrant(size_t n, double * restrict res)
      {
      const double hsqt2 = 0.707106781186547524400844362104849;
      size_t quart = n>>2;
      if ((n&7)==0)
        res[quart] = res[quart+1] = hsqt2;
      for (size_t i=2, j=2*quart-2; i<quart; i+=2, j-=2)
        {
        res[j  ] = res[i+1];
        res[j+1] = res[i  ];
        }
      }

    NOINLINE void fill_first_half(size_t n, double * restrict res)
      {
      size_t half = n>>1;
      if ((n&3)==0)
        for (size_t i=0; i<half; i+=2)
          {
          res[i+half]   = -res[i+1];
          res[i+half+1] =  res[i  ];
          }
      else
        for (size_t i=2, j=2*half-2; i<half; i+=2, j-=2)
          {
          res[j  ] = -res[i  ];
          res[j+1] =  res[i+1];
          }
      }

    NOINLINE void fill_second_half(size_t n, double * restrict res)
      {
      if ((n&1)==0)
        for (size_t i=0; i<n; ++i)
          res[i+n] = -res[i];
      else
        for (size_t i=2, j=2*n-2; i<n; i+=2, j-=2)
          {
          res[j  ] =  res[i  ];
          res[j+1] = -res[i+1];
          }
      }

    NOINLINE void sincos_2pibyn_half(size_t n, double * restrict res)
      {
      if ((n&3)==0)
        {
        calc_first_octant(n, res);
        fill_first_quadrant(n, res);
        fill_first_half(n, res);
        }
      else if ((n&1)==0)
        {
        calc_first_quadrant(n, res);
        fill_first_half(n, res);
        }
      else
        calc_first_half(n, res);
      }

  public:
    sincos_2pibyn(size_t n)
      : data(2*n)
      {
      sincos_2pibyn_half(n, data.data());
      fill_second_half(n, data.data());
      }

    double operator[](size_t idx) const { return data[idx]; }
  };

NOINLINE size_t largest_prime_factor (size_t n)
  {
  size_t res=1;
  size_t tmp;
  while (((tmp=(n>>1))<<1)==n)
    { res=2; n=tmp; }

  size_t limit=size_t(sqrt(n+0.01));
  for (size_t x=3; x<=limit; x+=2)
  while (((tmp=(n/x))*x)==n)
    {
    res=x;
    n=tmp;
    limit=size_t(sqrt(n+0.01));
    }
  if (n>1) res=n;

  return res;
  }

NOINLINE double cost_guess (size_t n)
  {
  const double lfp=1.1; // penalty for non-hardcoded larger factors
  size_t ni=n;
  double result=0.;
  size_t tmp;
  while (((tmp=(n>>1))<<1)==n)
    { result+=2; n=tmp; }

  size_t limit=size_t(sqrt(n+0.01));
  for (size_t x=3; x<=limit; x+=2)
  while ((tmp=(n/x))*x==n)
    {
    result+= (x<=5) ? x : lfp*x; // penalize larger prime factors
    n=tmp;
    limit=size_t(sqrt(n+0.01));
    }
  if (n>1) result+=(n<=5) ? n : lfp*n;

  return result*ni;
  }

/* returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n */
NOINLINE size_t good_size(size_t n)
  {
  if (n<=6) return n;

  size_t bestfac=2*n;
  for (size_t f2=1; f2<bestfac; f2*=2)
    for (size_t f23=f2; f23<bestfac; f23*=3)
      for (size_t f235=f23; f235<bestfac; f235*=5)
        for (size_t f2357=f235; f2357<bestfac; f2357*=7)
          for (size_t f235711=f2357; f235711<bestfac; f235711*=11)
            if (f235711>=n) bestfac=f235711;
  return bestfac;
  }

template<typename T> struct cmplx {
  T r, i;
};

#define PMC(a,b,c,d) { a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
#define ADDC(a,b,c) { a.r=b.r+c.r; a.i=b.i+c.i; }
#define SCALEC(a,b) { a.r*=b; a.i*=b; }
#define ROT90(a) { auto tmp_=a.r; a.r=-a.i; a.i=tmp_; }
#define ROTM90(a) { auto tmp_=-a.r; a.r=a.i; a.i=tmp_; }
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
#define WA(x,i) wa[(i)-1+(x)*(ido-1)]
/* a = b*c */
#define A_EQ_B_MUL_C(a,b,c) { a.r=b.r*c.r-b.i*c.i; a.i=b.r*c.i+b.i*c.r; }
/* a = conj(b)*c*/
#define A_EQ_CB_MUL_C(a,b,c) { a.r=b.r*c.r+b.i*c.i; a.i=b.r*c.i-b.i*c.r; }

#define PMSIGNC(a,b,c,d) { a.r=c.r+sign*d.r; a.i=c.i+sign*d.i; b.r=c.r-sign*d.r; b.i=c.i-sign*d.i; }
/* a = b*c */
#define MULPMSIGNC(a,b,c) { a.r=b.r*c.r-sign*b.i*c.i; a.i=b.r*c.i+sign*b.i*c.r; }
/* a *= b */
#define MULPMSIGNCEQ(a,b) { auto xtmp=a.r; a.r=b.r*a.r-sign*b.i*a.i; a.i=b.r*a.i+sign*b.i*xtmp; }

using dcmplx = cmplx<double>;

constexpr size_t NFCT=25;

class cfftp
  {
  private:

    struct fctdata
      {
      size_t fct;
      dcmplx *tw, *tws;
      };

    size_t length, nfct;
    arr<dcmplx> mem;
    fctdata fct[NFCT];

template<typename T, bool bwd> NOINLINE void pass2 (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa)
  {
  constexpr size_t cdim=2;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
  else
    for (size_t k=0; k<l1; ++k)
      {
      PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
      for (size_t i=1; i<ido; ++i)
        {
        T t;
        PMC (CH(i,k,0),t,CC(i,0,k),CC(i,1,k))
        if (bwd)
          A_EQ_B_MUL_C (CH(i,k,1),WA(0,i),t)
        else
          A_EQ_CB_MUL_C (CH(i,k,1),WA(0,i),t)
        }
      }
  }

#define PREP3(idx) \
        T t0 = CC(idx,0,k), t1, t2; \
        PMC (t1,t2,CC(idx,1,k),CC(idx,2,k)) \
        CH(idx,k,0).r=t0.r+t1.r; \
        CH(idx,k,0).i=t0.i+t1.i;
#define PARTSTEP3a(u1,u2,twr,twi) \
        { \
        T ca,cb; \
        ca.r=t0.r+twr*t1.r; \
        ca.i=t0.i+twr*t1.i; \
        cb.i=twi*t2.r; \
        cb.r=-(twi*t2.i); \
        PMC(CH(0,k,u1),CH(0,k,u2),ca,cb) \
        }

#define PARTSTEP3b(u1,u2,twr,twi) \
        { \
        T ca,cb,da,db; \
        ca.r=t0.r+twr*t1.r; \
        ca.i=t0.i+twr*t1.i; \
        cb.i=twi*t2.r; \
        cb.r=-(twi*t2.i); \
        PMC(da,db,ca,cb) \
        A_EQ_B_MUL_C (CH(i,k,u1),WA(u1-1,i),da) \
        A_EQ_B_MUL_C (CH(i,k,u2),WA(u2-1,i),db) \
        }
template<typename T> void pass3b (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa)
  {
  constexpr size_t cdim=3;
  constexpr double tw1r=-0.5, tw1i= 0.86602540378443864676;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      PREP3(0)
      PARTSTEP3a(1,2,tw1r,tw1i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      PREP3(0)
      PARTSTEP3a(1,2,tw1r,tw1i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        PREP3(i)
        PARTSTEP3b(1,2,tw1r,tw1i)
        }
      }
  }
#define PARTSTEP3f(u1,u2,twr,twi) \
        { \
        T ca,cb,da,db; \
        ca.r=t0.r+twr*t1.r; \
        ca.i=t0.i+twr*t1.i; \
        cb.i=twi*t2.r; \
        cb.r=-(twi*t2.i); \
        PMC(da,db,ca,cb) \
        A_EQ_CB_MUL_C (CH(i,k,u1),WA(u1-1,i),da) \
        A_EQ_CB_MUL_C (CH(i,k,u2),WA(u2-1,i),db) \
        }
template<typename T> NOINLINE void pass3f (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa)
  {
  const size_t cdim=3;
  const double tw1r=-0.5, tw1i= -0.86602540378443864676;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      PREP3(0)
      PARTSTEP3a(1,2,tw1r,tw1i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      PREP3(0)
      PARTSTEP3a(1,2,tw1r,tw1i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        PREP3(i)
        PARTSTEP3f(1,2,tw1r,tw1i)
        }
      }
  }

template<typename T> NOINLINE void pass4b (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa)
  {
  const size_t cdim=4;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      T t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k))
      PMC(t3,t4,CC(0,1,k),CC(0,3,k))
      ROT90(t4)
      PMC(CH(0,k,0),CH(0,k,2),t2,t3)
      PMC(CH(0,k,1),CH(0,k,3),t1,t4)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      T t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k))
      PMC(t3,t4,CC(0,1,k),CC(0,3,k))
      ROT90(t4)
      PMC(CH(0,k,0),CH(0,k,2),t2,t3)
      PMC(CH(0,k,1),CH(0,k,3),t1,t4)
      }
      for (size_t i=1; i<ido; ++i)
        {
        T c2, c3, c4, t1, t2, t3, t4;
        T cc0=CC(i,0,k), cc1=CC(i,1,k),cc2=CC(i,2,k),cc3=CC(i,3,k);
        PMC(t2,t1,cc0,cc2)
        PMC(t3,t4,cc1,cc3)
        ROT90(t4)
        dcmplx wa0=WA(0,i), wa1=WA(1,i),wa2=WA(2,i);
        PMC(CH(i,k,0),c3,t2,t3)
        PMC(c2,c4,t1,t4)
        A_EQ_B_MUL_C (CH(i,k,1),wa0,c2)
        A_EQ_B_MUL_C (CH(i,k,2),wa1,c3)
        A_EQ_B_MUL_C (CH(i,k,3),wa2,c4)
        }
      }
  }
template<typename T> NOINLINE void pass4f (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa)
  {
  const size_t cdim=4;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      T t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k))
      PMC(t3,t4,CC(0,1,k),CC(0,3,k))
      ROTM90(t4)
      PMC(CH(0,k,0),CH(0,k,2),t2,t3)
      PMC(CH(0,k,1),CH(0,k,3),t1,t4)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      T t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k))
      PMC(t3,t4,CC(0,1,k),CC(0,3,k))
      ROTM90(t4)
      PMC(CH(0,k,0),CH(0,k,2),t2,t3)
      PMC (CH(0,k,1),CH(0,k,3),t1,t4)
      }
      for (size_t i=1; i<ido; ++i)
        {
        T c2, c3, c4, t1, t2, t3, t4;
        T cc0=CC(i,0,k), cc1=CC(i,1,k),cc2=CC(i,2,k),cc3=CC(i,3,k);
        PMC(t2,t1,cc0,cc2)
        PMC(t3,t4,cc1,cc3)
        ROTM90(t4)
        dcmplx wa0=WA(0,i), wa1=WA(1,i),wa2=WA(2,i);
        PMC(CH(i,k,0),c3,t2,t3)
        PMC(c2,c4,t1,t4)
        A_EQ_CB_MUL_C (CH(i,k,1),wa0,c2)
        A_EQ_CB_MUL_C (CH(i,k,2),wa1,c3)
        A_EQ_CB_MUL_C (CH(i,k,3),wa2,c4)
        }
      }
  }

#define PREP5(idx) \
        T t0 = CC(idx,0,k), t1, t2, t3, t4; \
        PMC (t1,t4,CC(idx,1,k),CC(idx,4,k)) \
        PMC (t2,t3,CC(idx,2,k),CC(idx,3,k)) \
        CH(idx,k,0).r=t0.r+t1.r+t2.r; \
        CH(idx,k,0).i=t0.i+t1.i+t2.i;

#define PARTSTEP5a(u1,u2,twar,twbr,twai,twbi) \
        { \
        T ca,cb; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(CH(0,k,u1),CH(0,k,u2),ca,cb) \
        }

#define PARTSTEP5b(u1,u2,twar,twbr,twai,twbi) \
        { \
        T ca,cb,da,db; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(da,db,ca,cb) \
        A_EQ_B_MUL_C (CH(i,k,u1),WA(u1-1,i),da) \
        A_EQ_B_MUL_C (CH(i,k,u2),WA(u2-1,i),db) \
        }
template<typename T> NOINLINE void pass5b (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa)
  {
  const size_t cdim=5;
  const double tw1r= 0.3090169943749474241,
               tw1i= 0.95105651629515357212,
               tw2r= -0.8090169943749474241,
               tw2i= 0.58778525229247312917;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      PREP5(0)
      PARTSTEP5a(1,4,tw1r,tw2r,+tw1i,+tw2i)
      PARTSTEP5a(2,3,tw2r,tw1r,+tw2i,-tw1i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      PREP5(0)
      PARTSTEP5a(1,4,tw1r,tw2r,+tw1i,+tw2i)
      PARTSTEP5a(2,3,tw2r,tw1r,+tw2i,-tw1i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        PREP5(i)
        PARTSTEP5b(1,4,tw1r,tw2r,+tw1i,+tw2i)
        PARTSTEP5b(2,3,tw2r,tw1r,+tw2i,-tw1i)
        }
      }
  }
#define PARTSTEP5f(u1,u2,twar,twbr,twai,twbi) \
        { \
        T ca,cb,da,db; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(da,db,ca,cb) \
        A_EQ_CB_MUL_C (CH(i,k,u1),WA(u1-1,i),da) \
        A_EQ_CB_MUL_C (CH(i,k,u2),WA(u2-1,i),db) \
        }
template<typename T> NOINLINE void pass5f (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa)
  {
  const size_t cdim=5;
  const double tw1r= 0.3090169943749474241,
               tw1i= -0.95105651629515357212,
               tw2r= -0.8090169943749474241,
               tw2i= -0.58778525229247312917;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      PREP5(0)
      PARTSTEP5a(1,4,tw1r,tw2r,+tw1i,+tw2i)
      PARTSTEP5a(2,3,tw2r,tw1r,+tw2i,-tw1i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      PREP5(0)
      PARTSTEP5a(1,4,tw1r,tw2r,+tw1i,+tw2i)
      PARTSTEP5a(2,3,tw2r,tw1r,+tw2i,-tw1i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        PREP5(i)
        PARTSTEP5f(1,4,tw1r,tw2r,+tw1i,+tw2i)
        PARTSTEP5f(2,3,tw2r,tw1r,+tw2i,-tw1i)
        }
      }
  }

#define PREP7(idx) \
        T t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7; \
        PMC (t2,t7,CC(idx,1,k),CC(idx,6,k)) \
        PMC (t3,t6,CC(idx,2,k),CC(idx,5,k)) \
        PMC (t4,t5,CC(idx,3,k),CC(idx,4,k)) \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i;

#define PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,out1,out2) \
        { \
        T ca,cb; \
        ca.r=t1.r+x1*t2.r+x2*t3.r+x3*t4.r; \
        ca.i=t1.i+x1*t2.i+x2*t3.i+x3*t4.i; \
        cb.i=y1*t7.r y2*t6.r y3*t5.r; \
        cb.r=-(y1*t7.i y2*t6.i y3*t5.i); \
        PMC(out1,out2,ca,cb) \
        }
#define PARTSTEP7a(u1,u2,x1,x2,x3,y1,y2,y3) \
        PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,CH(0,k,u1),CH(0,k,u2))
#define PARTSTEP7(u1,u2,x1,x2,x3,y1,y2,y3) \
        { \
        T da,db; \
        PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,da,db) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

template<typename T>NOINLINE void pass7(size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa, const int sign)
  {
  const size_t cdim=7;
  const double tw1r= 0.623489801858733530525,
               tw1i= sign * 0.7818314824680298087084,
               tw2r= -0.222520933956314404289,
               tw2i= sign * 0.9749279121818236070181,
               tw3r= -0.9009688679024191262361,
               tw3i= sign * 0.4338837391175581204758;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      PREP7(0)
      PARTSTEP7a(1,6,tw1r,tw2r,tw3r,+tw1i,+tw2i,+tw3i)
      PARTSTEP7a(2,5,tw2r,tw3r,tw1r,+tw2i,-tw3i,-tw1i)
      PARTSTEP7a(3,4,tw3r,tw1r,tw2r,+tw3i,-tw1i,+tw2i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      PREP7(0)
      PARTSTEP7a(1,6,tw1r,tw2r,tw3r,+tw1i,+tw2i,+tw3i)
      PARTSTEP7a(2,5,tw2r,tw3r,tw1r,+tw2i,-tw3i,-tw1i)
      PARTSTEP7a(3,4,tw3r,tw1r,tw2r,+tw3i,-tw1i,+tw2i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        PREP7(i)
        PARTSTEP7(1,6,tw1r,tw2r,tw3r,+tw1i,+tw2i,+tw3i)
        PARTSTEP7(2,5,tw2r,tw3r,tw1r,+tw2i,-tw3i,-tw1i)
        PARTSTEP7(3,4,tw3r,tw1r,tw2r,+tw3i,-tw1i,+tw2i)
        }
      }
  }

#define PREP11(idx) \
        T t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; \
        PMC (t2,t11,CC(idx,1,k),CC(idx,10,k)) \
        PMC (t3,t10,CC(idx,2,k),CC(idx, 9,k)) \
        PMC (t4,t9 ,CC(idx,3,k),CC(idx, 8,k)) \
        PMC (t5,t8 ,CC(idx,4,k),CC(idx, 7,k)) \
        PMC (t6,t7 ,CC(idx,5,k),CC(idx, 6,k)) \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r+t5.r+t6.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i+t5.i+t6.i;

#define PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,out1,out2) \
        { \
        T ca,cb; \
        ca.r=t1.r+x1*t2.r+x2*t3.r+x3*t4.r+x4*t5.r+x5*t6.r; \
        ca.i=t1.i+x1*t2.i+x2*t3.i+x3*t4.i+x4*t5.i+x5*t6.i; \
        cb.i=y1*t11.r y2*t10.r y3*t9.r y4*t8.r y5*t7.r; \
        cb.r=-(y1*t11.i y2*t10.i y3*t9.i y4*t8.i y5*t7.i ); \
        PMC(out1,out2,ca,cb) \
        }
#define PARTSTEP11a(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) \
        PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,CH(0,k,u1),CH(0,k,u2))
#define PARTSTEP11(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) \
        { \
        T da,db; \
        PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,da,db) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

template<typename T> NOINLINE void pass11 (size_t ido, size_t l1, const T * restrict cc,
  T * restrict ch, const dcmplx * restrict wa, const int sign)
  {
  const size_t cdim=11;
  const double tw1r =        0.8412535328311811688618,
               tw1i = sign * 0.5406408174555975821076,
               tw2r =        0.4154150130018864255293,
               tw2i = sign * 0.9096319953545183714117,
               tw3r =       -0.1423148382732851404438,
               tw3i = sign * 0.9898214418809327323761,
               tw4r =       -0.6548607339452850640569,
               tw4i = sign * 0.755749574354258283774,
               tw5r =       -0.9594929736144973898904,
               tw5i = sign * 0.2817325568414296977114;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      PREP11(0)
      PARTSTEP11a(1,10,tw1r,tw2r,tw3r,tw4r,tw5r,+tw1i,+tw2i,+tw3i,+tw4i,+tw5i)
      PARTSTEP11a(2, 9,tw2r,tw4r,tw5r,tw3r,tw1r,+tw2i,+tw4i,-tw5i,-tw3i,-tw1i)
      PARTSTEP11a(3, 8,tw3r,tw5r,tw2r,tw1r,tw4r,+tw3i,-tw5i,-tw2i,+tw1i,+tw4i)
      PARTSTEP11a(4, 7,tw4r,tw3r,tw1r,tw5r,tw2r,+tw4i,-tw3i,+tw1i,+tw5i,-tw2i)
      PARTSTEP11a(5, 6,tw5r,tw1r,tw4r,tw2r,tw3r,+tw5i,-tw1i,+tw4i,-tw2i,+tw3i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      PREP11(0)
      PARTSTEP11a(1,10,tw1r,tw2r,tw3r,tw4r,tw5r,+tw1i,+tw2i,+tw3i,+tw4i,+tw5i)
      PARTSTEP11a(2, 9,tw2r,tw4r,tw5r,tw3r,tw1r,+tw2i,+tw4i,-tw5i,-tw3i,-tw1i)
      PARTSTEP11a(3, 8,tw3r,tw5r,tw2r,tw1r,tw4r,+tw3i,-tw5i,-tw2i,+tw1i,+tw4i)
      PARTSTEP11a(4, 7,tw4r,tw3r,tw1r,tw5r,tw2r,+tw4i,-tw3i,+tw1i,+tw5i,-tw2i)
      PARTSTEP11a(5, 6,tw5r,tw1r,tw4r,tw2r,tw3r,+tw5i,-tw1i,+tw4i,-tw2i,+tw3i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        PREP11(i)
        PARTSTEP11(1,10,tw1r,tw2r,tw3r,tw4r,tw5r,+tw1i,+tw2i,+tw3i,+tw4i,+tw5i)
        PARTSTEP11(2, 9,tw2r,tw4r,tw5r,tw3r,tw1r,+tw2i,+tw4i,-tw5i,-tw3i,-tw1i)
        PARTSTEP11(3, 8,tw3r,tw5r,tw2r,tw1r,tw4r,+tw3i,-tw5i,-tw2i,+tw1i,+tw4i)
        PARTSTEP11(4, 7,tw4r,tw3r,tw1r,tw5r,tw2r,+tw4i,-tw3i,+tw1i,+tw5i,-tw2i)
        PARTSTEP11(5, 6,tw5r,tw1r,tw4r,tw2r,tw3r,+tw5i,-tw1i,+tw4i,-tw2i,+tw3i)
        }
      }
  }

#define CX(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define CX2(a,b) cc[(a)+idl1*(b)]
#define CH2(a,b) ch[(a)+idl1*(b)]

template<typename T> NOINLINE void passg (size_t ido, size_t ip, size_t l1,
  T * restrict cc, T * restrict ch, const dcmplx * restrict wa,
  const dcmplx * restrict csarr, const int sign)
  {
  const size_t cdim=ip;
  size_t ipph = (ip+1)/2;
  size_t idl1 = ido*l1;

  arr<dcmplx> wal(ip);
  wal[0]=(dcmplx){1.,0.};
  for (size_t i=1; i<ip; ++i)
    wal[i]=(dcmplx){csarr[i].r,sign*csarr[i].i};

  for (size_t k=0; k<l1; ++k)
    for (size_t i=0; i<ido; ++i)
      CH(i,k,0) = CC(i,0,k);
  for (size_t j=1, jc=ip-1; j<ipph; ++j, --jc)
    for (size_t k=0; k<l1; ++k)
      for (size_t i=0; i<ido; ++i)
        PMC(CH(i,k,j),CH(i,k,jc),CC(i,j,k),CC(i,jc,k))
  for (size_t k=0; k<l1; ++k)
    for (size_t i=0; i<ido; ++i)
      {
      T tmp = CH(i,k,0);
      for (size_t j=1; j<ipph; ++j)
        ADDC(tmp,tmp,CH(i,k,j))
      CX(i,k,0) = tmp;
      }
  for (size_t l=1, lc=ip-1; l<ipph; ++l, --lc)
    {
    // j=0
    for (size_t ik=0; ik<idl1; ++ik)
      {
      CX2(ik,l).r = CH2(ik,0).r+wal[l].r*CH2(ik,1).r+wal[2*l].r*CH2(ik,2).r;
      CX2(ik,l).i = CH2(ik,0).i+wal[l].r*CH2(ik,1).i+wal[2*l].r*CH2(ik,2).i;
      CX2(ik,lc).r=-wal[l].i*CH2(ik,ip-1).i-wal[2*l].i*CH2(ik,ip-2).i;
      CX2(ik,lc).i=wal[l].i*CH2(ik,ip-1).r+wal[2*l].i*CH2(ik,ip-2).r;
      }

    size_t iwal=2*l;
    size_t j=3, jc=ip-3;
    for (; j<ipph-1; j+=2, jc-=2)
      {
      iwal+=l; if (iwal>ip) iwal-=ip;
      dcmplx xwal=wal[iwal];
      iwal+=l; if (iwal>ip) iwal-=ip;
      dcmplx xwal2=wal[iwal];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        CX2(ik,l).r += CH2(ik,j).r*xwal.r+CH2(ik,j+1).r*xwal2.r;
        CX2(ik,l).i += CH2(ik,j).i*xwal.r+CH2(ik,j+1).i*xwal2.r;
        CX2(ik,lc).r -= CH2(ik,jc).i*xwal.i+CH2(ik,jc-1).i*xwal2.i;
        CX2(ik,lc).i += CH2(ik,jc).r*xwal.i+CH2(ik,jc-1).r*xwal2.i;
        }
      }
    for (; j<ipph; ++j, --jc)
      {
      iwal+=l; if (iwal>ip) iwal-=ip;
      dcmplx xwal=wal[iwal];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        CX2(ik,l).r += CH2(ik,j).r*xwal.r;
        CX2(ik,l).i += CH2(ik,j).i*xwal.r;
        CX2(ik,lc).r -= CH2(ik,jc).i*xwal.i;
        CX2(ik,lc).i += CH2(ik,jc).r*xwal.i;
        }
      }
    }

  // shuffling and twiddling
  if (ido==1)
    for (size_t j=1, jc=ip-1; j<ipph; ++j, --jc)
      for (size_t ik=0; ik<idl1; ++ik)
        {
        T t1=CX2(ik,j), t2=CX2(ik,jc);
        PMC(CX2(ik,j),CX2(ik,jc),t1,t2)
        }
  else
    {
    for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)
      for (size_t k=0; k<l1; ++k)
        {
        T t1=CX(0,k,j), t2=CX(0,k,jc);
        PMC(CX(0,k,j),CX(0,k,jc),t1,t2)
        for (size_t i=1; i<ido; ++i)
          {
          T x1, x2;
          PMC(x1,x2,CX(i,k,j),CX(i,k,jc))
          size_t idij=(j-1)*(ido-1)+i-1;
          MULPMSIGNC (CX(i,k,j),wa[idij],x1)
          idij=(jc-1)*(ido-1)+i-1;
          MULPMSIGNC (CX(i,k,jc),wa[idij],x2)
          }
        }
    }
  }

#undef CH2
#undef CX2
#undef CX

template<typename T> void pass_all(T c[], double fact,
  const int sign)
  {
  if (length==1) return;
  size_t l1=1, nf=nfct;
  arr<T> ch(length);
  T *p1=c, *p2=ch.data();

  for(size_t k1=0; k1<nf; k1++)
    {
    size_t ip=fct[k1].fct;
    size_t l2=ip*l1;
    size_t ido = length/l2;
    if     (ip==4)
      sign>0 ? pass4b (ido, l1, p1, p2, fct[k1].tw)
             : pass4f (ido, l1, p1, p2, fct[k1].tw);
    else if(ip==2)
      sign>0 ? pass2<T, true>(ido, l1, p1, p2, fct[k1].tw)
             : pass2<T, false>(ido, l1, p1, p2, fct[k1].tw);
    else if(ip==3)
      sign>0 ? pass3b (ido, l1, p1, p2, fct[k1].tw)
             : pass3f (ido, l1, p1, p2, fct[k1].tw);
    else if(ip==5)
      sign>0 ? pass5b (ido, l1, p1, p2, fct[k1].tw)
             : pass5f (ido, l1, p1, p2, fct[k1].tw);
    else if(ip==7)  pass7 (ido, l1, p1, p2, fct[k1].tw, sign);
    else if(ip==11) pass11(ido, l1, p1, p2, fct[k1].tw, sign);
    else
      {
      passg(ido, ip, l1, p1, p2, fct[k1].tw, fct[k1].tws, sign);
      swap(p1,p2);
      }
    swap(p1,p2);
    l1=l2;
    }
  if (p1!=c)
    {
    if (fact!=1.)
      for (size_t i=0; i<length; ++i)
        {
        c[i].r = ch[i].r*fact;
        c[i].i = ch[i].i*fact;
        }
    else
      memcpy (c,p1,length*sizeof(T));
    }
  else
    if (fact!=1.)
      for (size_t i=0; i<length; ++i)
        {
        c[i].r *= fact;
        c[i].i *= fact;
        }
  }

#undef PMSIGNC
#undef A_EQ_B_MUL_C
#undef A_EQ_CB_MUL_C
#undef MULPMSIGNC
#undef MULPMSIGNCEQ

#undef WA
#undef CC
#undef CH
#undef ROT90
#undef SCALEC
#undef ADDC
#undef PMC

public:

template<typename T> NOINLINE
void forward(T c[], double fct)
  { pass_all(c, fct, -1); }

template<typename T> NOINLINE
void backward(T c[], double fct)
  { pass_all(c, fct, 1); }

private:

NOINLINE void factorize ()
  {
  nfct=0;
  size_t len=length;
  while ((len%4)==0)
    { if (nfct>=NFCT) throw OOM(); fct[nfct++].fct=4; len>>=2; }
  if ((len%2)==0)
    {
    len>>=1;
    // factor 2 should be at the front of the factor list
    if (nfct>=NFCT) throw OOM();
    fct[nfct++].fct=2;
    swap(fct[0].fct, fct[nfct-1].fct);
    }
  size_t maxl=(size_t)(sqrt((double)len))+1;
  for (size_t divisor=3; (len>1)&&(divisor<maxl); divisor+=2)
    if ((len%divisor)==0)
      {
      while ((len%divisor)==0)
        {
        if (nfct>=NFCT) throw OOM();
        fct[nfct++].fct=divisor;
        len/=divisor;
        }
      maxl=(size_t)(sqrt((double)len))+1;
      }
  if (len>1) fct[nfct++].fct=len;
  }

NOINLINE size_t twsize ()
  {
  size_t twsize=0, l1=1;
  for (size_t k=0; k<nfct; ++k)
    {
    size_t ip=fct[k].fct, ido= length/(l1*ip);
    twsize+=(ip-1)*(ido-1);
    if (ip>11)
      twsize+=ip;
    l1*=ip;
    }
  return twsize;
  }

NOINLINE void comp_twiddle()
  {
  sincos_2pibyn twid(length);
  size_t l1=1;
  size_t memofs=0;
  for (size_t k=0; k<nfct; ++k)
    {
    size_t ip=fct[k].fct, ido= length/(l1*ip);
    fct[k].tw=mem.data()+memofs;
    memofs+=(ip-1)*(ido-1);
    for (size_t j=1; j<ip; ++j)
      for (size_t i=1; i<ido; ++i)
        {
        fct[k].tw[(j-1)*(ido-1)+i-1].r = twid[2*j*l1*i];
        fct[k].tw[(j-1)*(ido-1)+i-1].i = twid[2*j*l1*i+1];
        }
    if (ip>11)
      {
      fct[k].tws=mem.data()+memofs;
      memofs+=ip;
      for (size_t j=0; j<ip; ++j)
        {
        fct[k].tws[j].r = twid[2*j*l1*ido];
        fct[k].tws[j].i = twid[2*j*l1*ido+1];
        }
      }
    l1*=ip;
    }
  }

public:

cfftp(size_t length_)
  {
  length=length_;
  if (length==0) throw 42;
  nfct=0;
  if (length==1) return;
  factorize();
  size_t tws=twsize();
  mem.resize(tws);
  comp_twiddle();
  }

};

struct fftblue
  {
  size_t n, n2;
  cfftp plan;
  arr<double> mem;
  double *bk, *bkf;

  fftblue(size_t length)
    : n(length), n2(good_size(n*2-1)), plan(n2), mem(2*(n+n2)), bk(mem.data()), bkf(mem.data()+2*n)
    {
/* initialize b_k */
  sincos_2pibyn tmp(2*n);
  bk[0] = 1;
  bk[1] = 0;

  size_t coeff=0;
  for (size_t m=1; m<n; ++m)
    {
    coeff+=2*m-1;
    if (coeff>=2*n) coeff-=2*n;
    bk[2*m  ] = tmp[2*coeff  ];
    bk[2*m+1] = tmp[2*coeff+1];
    }

  /* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
  double xn2 = 1./n2;
  bkf[0] = bk[0]*xn2;
  bkf[1] = bk[1]*xn2;
  for (size_t m=2; m<2*n; m+=2)
    {
    bkf[m]   = bkf[2*n2-m]   = bk[m]   *xn2;
    bkf[m+1] = bkf[2*n2-m+1] = bk[m+1] *xn2;
    }
  for (size_t m=2*n;m<=(2*n2-2*n+1);++m)
    bkf[m]=0.;
  plan.forward((dcmplx *)bkf,1.);
    }

template<typename T> NOINLINE
void fftblue_fft(T c[], int isign, double fct)
  {
  arr<T> akf(2*n2);

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
    akf[m]=0.*c[0];

  plan.forward ((cmplx<T> *)akf.data(),fct);

/* do the convolution */
  if (isign>0)
    for (size_t m=0; m<2*n2; m+=2)
      {
      T im = -akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  =  akf[m]*bkf[m]   + akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }
  else
    for (size_t m=0; m<2*n2; m+=2)
      {
      T im = akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  = akf[m]*bkf[m]   - akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }

/* inverse FFT */
  plan.backward ((cmplx<T> *)akf.data(),1.);

/* multiply by b_k */
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

template<typename T>
void backward(T c[], double fct)
  { fftblue_fft(c,1,fct); }

template<typename T>
void forward(T c[], double fct)
  { fftblue_fft(c,-1,fct); }
};

} // unnamed namespace

struct pocketfft_c
  {
  private:
    unique_ptr<cfftp> packplan;
    unique_ptr<fftblue> blueplan;

  public:
    pocketfft_c(size_t length) : packplan(nullptr), blueplan(nullptr)
      {
      if (length==0) throw 42;
      if ((length<50) || (largest_prime_factor(length)<=sqrt(length)))
        {
        packplan=make_unique<cfftp>(length);
        return;
        }
      double comp1 = cost_guess(length);
      double comp2 = 2*cost_guess(good_size(2*length-1));
      comp2*=1.5; /* fudge factor that appears to give good overall performance */
      if (comp2<comp1) // use Bluestein
        blueplan=make_unique<fftblue>(length);
      else
        packplan=make_unique<cfftp>(length);
      }

template<typename T> void backward(T c[], double fct)
  {
  if (packplan)
    return packplan->backward((dcmplx *)c,fct);
  return blueplan->backward(c,fct);
  }

template<typename T> void forward(T c[], double fct)
  {
  if (packplan)
    return packplan->forward((dcmplx *)c,fct);
  return blueplan->forward(c,fct);
  }
  };
#define maxlen 8192

void fill_random (double *data, size_t length)
  {
  for (size_t m=0; m<length; ++m)
    data[m] = rand()/(RAND_MAX+1.0)-0.5;
  }

double errcalc (double *data, double *odata, size_t length)
  {
  double sum = 0, errsum = 0;
  for (size_t m=0; m<length; ++m)
    {
    errsum += (data[m]-odata[m])*(data[m]-odata[m]);
    sum += odata[m]*odata[m];
    }
  return sqrt(errsum/sum);
  }
#include <stdio.h>
template<typename T, int n> int test_complex(void)
  {
  double data[2*n*maxlen], odata[2*n*maxlen];
  fill_random (odata, 2*n*maxlen);
  const double epsilon=2e-15;
  int ret = 0;
  double errsum=0;
  for (int length=1; length<=maxlen; ++length)
    {
    memcpy (data,odata,2*n*length*sizeof(double));
    pocketfft_c plan(length);
for (int x=0; x<1; ++x)
  {
    plan.forward((T *)data, 1.);
    plan.backward((T *)data, 1./length);
  }
    double err = errcalc (data, odata, 2*n*length);
    if (err>epsilon)
      {
      printf("problem at complex length %i: %e\n",length,err);
      ret = 1;
      }
    errsum+=err;
    }
  printf("errsum: %e\n",errsum);
  return ret;
  }

#include <iostream>
#include <vector>
#include <complex>
#include <x86intrin.h>

int main()
  {
  double *a=RALLOC(double,1234560);
  test_complex<double,1>();
  }
