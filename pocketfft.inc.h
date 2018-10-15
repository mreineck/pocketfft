/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*
 *  Auxiliary implementation file, included from pocketfft.c
 *
 *  Copyright (C) 2004-2018 Max-Planck-Society
 *  \author Martin Reinecke
 */

#define PMC(a,b,c,d) { a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
#define ADDC(a,b,c) { a.r=b.r+c.r; a.i=b.i+c.i; }
#define SCALEC(a,b) { a.r*=b; a.i*=b; }
#define CONJFLIPC(a) { double tmp_=a.r; a.r=-a.i; a.i=tmp_; }
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
#define WA(x,i) wa[(i)-1+(x)*(ido-1)]

#ifdef BACKWARD
#define PSIGN +
#define PMSIGNC(a,b,c,d) { a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
/* a = b*c */
#define MULPMSIGNC(a,b,c) { a.r=b.r*c.r-b.i*c.i; a.i=b.r*c.i+b.i*c.r; }
/* a *= b */
#define MULPMSIGNCEQ(a,b) { double xtmp=a.r; a.r=b.r*a.r-b.i*a.i; a.i=b.r*a.i+b.i*xtmp; }
#else
#define PSIGN -
#define PMSIGNC(a,b,c,d) { a.r=c.r-d.r; a.i=c.i-d.i; b.r=c.r+d.r; b.i=c.i+d.i; }
/* a = conj(b)*c */
#define MULPMSIGNC(a,b,c) { a.r=b.r*c.r+b.i*c.i; a.i=b.r*c.i-b.i*c.r; }
/* a *= conj(b) */
#define MULPMSIGNCEQ(a,b) { double xtmp=a.r; a.r=b.r*a.r+b.i*a.i; a.i=b.r*a.i-b.i*xtmp; }
#endif

NOINLINE static void X(2) (size_t ido, size_t l1, const cmplx * restrict cc,
  cmplx * restrict ch, const cmplx * restrict wa)
  {
  const size_t cdim=2;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
  else
    for (size_t k=0; k<l1; ++k)
      {
      PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
      for (size_t i=1; i<ido; ++i)
        {
        cmplx t;
        PMC (CH(i,k,0),t,CC(i,0,k),CC(i,1,k))
        MULPMSIGNC (CH(i,k,1),WA(0,i),t)
        }
      }
  }

#define PREP3(idx) \
        cmplx t0 = CC(idx,0,k), t1, t2; \
        PMC (t1,t2,CC(idx,1,k),CC(idx,2,k)) \
        CH(idx,k,0).r=t0.r+t1.r; \
        CH(idx,k,0).i=t0.i+t1.i;
#define PARTSTEP3a(u1,u2,twr,twi) \
        { \
        cmplx ca,cb; \
        ca.r=t0.r+twr*t1.r; \
        ca.i=t0.i+twr*t1.i; \
        cb.i=twi*t2.r; \
        cb.r=-(twi*t2.i); \
        PMC(CH(0,k,u1),CH(0,k,u2),ca,cb) \
        }
#define PARTSTEP3(u1,u2,twr,twi) \
        { \
        cmplx ca,cb,da,db; \
        ca.r=t0.r+twr*t1.r; \
        ca.i=t0.i+twr*t1.i; \
        cb.i=twi*t2.r; \
        cb.r=-(twi*t2.i); \
        PMC(da,db,ca,cb) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

NOINLINE static void X(3)(size_t ido, size_t l1, const cmplx * restrict cc,
  cmplx * restrict ch, const cmplx * restrict wa)
  {
  const size_t cdim=3;
  static const double tw1r=-0.5, tw1i= PSIGN 0.86602540378443864676;

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
        PARTSTEP3(1,2,tw1r,tw1i)
        }
      }
  }

NOINLINE static void X(4)(size_t ido, size_t l1, const cmplx * restrict cc,
  cmplx * restrict ch, const cmplx * restrict wa)
  {
  const size_t cdim=4;

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      cmplx t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k))
      PMC(t3,t4,CC(0,1,k),CC(0,3,k))
      CONJFLIPC(t4)
      PMC(CH(0,k,0),CH(0,k,2),t2,t3)
      PMSIGNC (CH(0,k,1),CH(0,k,3),t1,t4)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      cmplx t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k))
      PMC(t3,t4,CC(0,1,k),CC(0,3,k))
      CONJFLIPC(t4)
      PMC(CH(0,k,0),CH(0,k,2),t2,t3)
      PMSIGNC (CH(0,k,1),CH(0,k,3),t1,t4)
      }
      for (size_t i=1; i<ido; ++i)
        {
        cmplx c2, c3, c4, t1, t2, t3, t4;
        cmplx cc0=CC(i,0,k), cc1=CC(i,1,k),cc2=CC(i,2,k),cc3=CC(i,3,k);
        PMC(t2,t1,cc0,cc2)
        PMC(t3,t4,cc1,cc3)
        CONJFLIPC(t4)
        cmplx wa0=WA(0,i), wa1=WA(1,i),wa2=WA(2,i);
        PMC(CH(i,k,0),c3,t2,t3)
        PMSIGNC (c2,c4,t1,t4)
        MULPMSIGNC (CH(i,k,1),wa0,c2)
        MULPMSIGNC (CH(i,k,2),wa1,c3)
        MULPMSIGNC (CH(i,k,3),wa2,c4)
        }
      }
  }

#define PREP5(idx) \
        cmplx t0 = CC(idx,0,k), t1, t2, t3, t4; \
        PMC (t1,t4,CC(idx,1,k),CC(idx,4,k)) \
        PMC (t2,t3,CC(idx,2,k),CC(idx,3,k)) \
        CH(idx,k,0).r=t0.r+t1.r+t2.r; \
        CH(idx,k,0).i=t0.i+t1.i+t2.i;

#define PARTSTEP5a(u1,u2,twar,twbr,twai,twbi) \
        { \
        cmplx ca,cb; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(CH(0,k,u1),CH(0,k,u2),ca,cb) \
        }
#define PARTSTEP5(u1,u2,twar,twbr,twai,twbi) \
        { \
        cmplx ca,cb,da,db; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(da,db,ca,cb) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

NOINLINE static void X(5)(size_t ido, size_t l1, const cmplx * restrict cc,
  cmplx * restrict ch, const cmplx * restrict wa)
  {
  const size_t cdim=5;
  const double tw1r= 0.3090169943749474241,
               tw1i= PSIGN 0.95105651629515357212,
               tw2r= -0.8090169943749474241,
               tw2i= PSIGN 0.58778525229247312917;

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
        PARTSTEP5(1,4,tw1r,tw2r,+tw1i,+tw2i)
        PARTSTEP5(2,3,tw2r,tw1r,+tw2i,-tw1i)
        }
      }
  }

#define PREP7(idx) \
        cmplx t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7; \
        PMC (t2,t7,CC(idx,1,k),CC(idx,6,k)) \
        PMC (t3,t6,CC(idx,2,k),CC(idx,5,k)) \
        PMC (t4,t5,CC(idx,3,k),CC(idx,4,k)) \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i;

#define PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,out1,out2) \
        { \
        cmplx ca,cb; \
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
        cmplx da,db; \
        PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,da,db) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

NOINLINE static void X(7)(size_t ido, size_t l1, const cmplx * restrict cc,
  cmplx * restrict ch, const cmplx * restrict wa)
  {
  const size_t cdim=7;
  const double tw1r= 0.623489801858733530525,
               tw1i= PSIGN 0.7818314824680298087084,
               tw2r= -0.222520933956314404289,
               tw2i= PSIGN 0.9749279121818236070181,
               tw3r= -0.9009688679024191262361,
               tw3i= PSIGN 0.4338837391175581204758;

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
        cmplx t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; \
        PMC (t2,t11,CC(idx,1,k),CC(idx,10,k)) \
        PMC (t3,t10,CC(idx,2,k),CC(idx, 9,k)) \
        PMC (t4,t9 ,CC(idx,3,k),CC(idx, 8,k)) \
        PMC (t5,t8 ,CC(idx,4,k),CC(idx, 7,k)) \
        PMC (t6,t7 ,CC(idx,5,k),CC(idx, 6,k)) \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r+t5.r+t6.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i+t5.i+t6.i;

#define PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,out1,out2) \
        { \
        cmplx ca,cb; \
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
        cmplx da,db; \
        PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,da,db) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

NOINLINE static void X(11)(size_t ido, size_t l1, const cmplx * restrict cc,
  cmplx * restrict ch, const cmplx * restrict wa)
  {
  const size_t cdim=11;
  const double tw1r =        0.8412535328311811688618,
               tw1i = PSIGN  0.5406408174555975821076,
               tw2r =        0.4154150130018864255293,
               tw2i = PSIGN  0.9096319953545183714117,
               tw3r =       -0.1423148382732851404438,
               tw3i = PSIGN  0.9898214418809327323761,
               tw4r =       -0.6548607339452850640569,
               tw4i = PSIGN  0.755749574354258283774,
               tw5r =       -0.9594929736144973898904,
               tw5i = PSIGN  0.2817325568414296977114;

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

NOINLINE static int X(g)(size_t ido, size_t ip, size_t l1,
  cmplx * restrict cc, cmplx * restrict ch, const cmplx * restrict wa,
  const cmplx * restrict csarr)
  {
  const size_t cdim=ip;
  size_t ipph = (ip+1)/2;
  size_t idl1 = ido*l1;

  cmplx * restrict wal=RALLOC(cmplx,ip);
  if (!wal) return -1;
  wal[0]=(cmplx){1.,0.};
  for (size_t i=1; i<ip; ++i)
    wal[i]=(cmplx){csarr[i].r,PSIGN csarr[i].i};

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
      cmplx tmp = CH(i,k,0);
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
      cmplx xwal=wal[iwal];
      iwal+=l; if (iwal>ip) iwal-=ip;
      cmplx xwal2=wal[iwal];
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
      cmplx xwal=wal[iwal];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        CX2(ik,l).r += CH2(ik,j).r*xwal.r;
        CX2(ik,l).i += CH2(ik,j).i*xwal.r;
        CX2(ik,lc).r -= CH2(ik,jc).i*xwal.i;
        CX2(ik,lc).i += CH2(ik,jc).r*xwal.i;
        }
      }
    }
  DEALLOC(wal);

  // shuffling and twiddling
  if (ido==1)
    for (size_t j=1, jc=ip-1; j<ipph; ++j, --jc)
      for (size_t ik=0; ik<idl1; ++ik)
        {
        cmplx t1=CX2(ik,j), t2=CX2(ik,jc);
        PMC(CX2(ik,j),CX2(ik,jc),t1,t2)
        }
  else
    {
    for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)
      for (size_t k=0; k<l1; ++k)
        {
        cmplx t1=CX(0,k,j), t2=CX(0,k,jc);
        PMC(CX(0,k,j),CX(0,k,jc),t1,t2)
        for (size_t i=1; i<ido; ++i)
          {
          cmplx x1, x2;
          PMC(x1,x2,CX(i,k,j),CX(i,k,jc))
          size_t idij=(j-1)*(ido-1)+i-1;
          MULPMSIGNC (CX(i,k,j),wa[idij],x1)
          idij=(jc-1)*(ido-1)+i-1;
          MULPMSIGNC (CX(i,k,jc),wa[idij],x2)
          }
        }
    }
  return 0;
  }

#undef CH2
#undef CX2
#undef CX

WARN_UNUSED_RESULT static int X(_all)(cfftp_plan plan, cmplx c[], double fct)
  {
  if (plan->length==1) return 0;
  size_t len=plan->length;
  size_t l1=1, nf=plan->nfct;
  cmplx *ch = RALLOC(cmplx, len);
  if (!ch) return -1;
  cmplx *p1=c, *p2=ch;

  for(size_t k1=0; k1<nf; k1++)
    {
    size_t ip=plan->fct[k1].fct;
    size_t l2=ip*l1;
    size_t ido = len/l2;
    if     (ip==4)  X( 4)(ido, l1, p1, p2, plan->fct[k1].tw);
    else if(ip==2)  X( 2)(ido, l1, p1, p2, plan->fct[k1].tw);
    else if(ip==3)  X( 3)(ido, l1, p1, p2, plan->fct[k1].tw);
    else if(ip==5)  X( 5)(ido, l1, p1, p2, plan->fct[k1].tw);
    else if(ip==7)  X( 7)(ido, l1, p1, p2, plan->fct[k1].tw);
    else if(ip==11) X(11)(ido, l1, p1, p2, plan->fct[k1].tw);
    else
      {
      if (X( g)(ido, ip, l1, p1, p2, plan->fct[k1].tw, plan->fct[k1].tws))
        { DEALLOC(ch); return -1; }
      SWAP(p1,p2,cmplx *);
      }
    SWAP(p1,p2,cmplx *);
    l1=l2;
    }
  if (p1!=c)
    {
    if (fct!=1.)
      for (size_t i=0; i<len; ++i)
        {
        c[i].r = ch[i].r*fct;
        c[i].i = ch[i].i*fct;
        }
    else
      memcpy (c,p1,len*sizeof(cmplx));
    }
  else
    if (fct!=1.)
      for (size_t i=0; i<len; ++i)
        {
        c[i].r *= fct;
        c[i].i *= fct;
        }
  DEALLOC(ch);
  return 0;
  }

#undef PSIGN
#undef PMSIGNC
#undef MULPMSIGNC
#undef MULPMSIGNCEQ

#undef WA
#undef CC
#undef CH
#undef CONJFLIPC
#undef SCALEC
#undef ADDC
#undef PMC
