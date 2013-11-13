//*************************************************************************
//* ===================
//*  EXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXHit.
//* (Requires)
//* (Provides)
//*     class EXMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//

#include "EXMeasLayer.h"
#include "EXKalDetector.h"
#include "EXHit.h"
#include "TRandom.h"

ClassImp(EXMeasLayer)

Int_t EXMeasLayer::CalcCellNo(const TVector3 &xv) const
{
   // Eq for the virtual wire passing through xv
   //
   //   	x = xv + t*et
   // wehre
   //      et: wire dir vector at xv
   //      t:  distance from xv along the virtual wire 
   //
   // et is given by
   //		et = sin(alpha) * (ez x xv)/|(ez x xv)|
   //		   + cos(alpha) * ez
   // where ez is the unit vector in the z-direction
   // and alpha is the stereo angle.
   //      cos(alpha) = wd.Z();
   //      sin(alpha) = wd * (ez x we)/|(ez x we)|
   //
   // We can then solve for t = t(z = wire end) from
   // 		we.Z() = xv.Z() + t * et.Z()
   //

   static const TVector3 ez(0., 0., 1.);
   TVector3 epe = ez.Cross(fWireEnd).Unit();
   Double_t csa = fWireDir.Z();
   Double_t sna = fWireDir * epe;
   TVector3 epv = ez.Cross(xv).Unit();
   TVector3 et  = sna * epv + csa * ez;
   Double_t t   = (fWireEnd.Z() - xv.Z())/csa;
   TVector3 xve = xv + t * et;
   Double_t phi = TMath::ATan2(xve.Y(), xve.X());
   Double_t dfi = phi - (fWirePhi - fDphi/2);

   Int_t cellno =  (Int_t)(dfi > 0. ? dfi/fDphi : dfi/fDphi-1);

   return cellno;
}
                                                                                
TKalMatrix EXMeasLayer::XvToMv(const TVTrackHit &vht,
                               const TVector3   &xv) const
{
   const EXHit &ht = *dynamic_cast<const EXHit *> (&vht);
   Int_t lr     = ht.GetLR();
   Int_t cellno = ht.GetCellNo();

   return XvToMv(xv, lr, cellno);
}

TKalMatrix EXMeasLayer::XvToMv(const TVector3 &xv,
                                     Int_t    &lr,
                                     Int_t    &cellno) const
{
   // Calculate cell No and rotate wire dir and end vectors

   if (cellno < -99999) cellno = CalcCellNo(xv);
   TVector3 wd = GetWireDir(cellno);
   TVector3 we = GetWireEnd(cellno);

   // Calculate hit coordinate information:
   //	mv(0,0) = drift length
   //     (1,0) = z from charge division

   TKalMatrix mv(kMdim,1);
   Double_t t  = wd * (xv - we);
   TVector3 xw = we + t * wd;
   mv(0,0)  = (xw - xv).Mag();   // drift length
   mv(1,0)  = xw.Z();            // z from charge division
   if (!lr) {
      TVector2 xyv  = (*(TVector3 *)&xv).XYvector();
      TVector2 xyw  = xw.XYvector();
      Double_t test = xyw.X()*xyv.Y() - xyw.Y()*xyv.X();
      lr   = (test < 0. ? -1 : 1);           // left/right
   }
   mv(0,0) *= lr;
   return mv;
}


TVector3 EXMeasLayer::MvToXv(const TKalMatrix &mv,
                                   Int_t       lr,
                                   Int_t       cellno) const
{
   // Caution !!!!!
   // The resultant Xv is not exactly on the measurement surface!
   // Use this with care.
   //
   // Rotate fWireDir and fWireEnd to the hit cell
   // Notice that cellno can be -ve!

   TVector3 wd = GetWireDir(cellno);
   TVector3 we = GetWireEnd(cellno);

   // Calculate the 3-d hit position

   Double_t hlen = 0.5*GetLength();
   Double_t dlen = mv(0,0) * lr;
   Double_t zh   = mv(1,0);
   Double_t t    = (hlen + zh - GetXc().Z())/wd.Z();
   TVector3 xhv  = we + t * wd;
   Double_t tana = GetTanA();
   TVector3 ns(xhv.X(), xhv.Y(), -xhv.Z() + tana*tana);
   TVector3 ed = wd.Cross(ns).Unit();
#if 0
   Double_t rw = xhv.Perp();
#if 1
   Double_t fi = 0.5 * lr * (dlen/rw);
   ed.Rotate(fi,wd);
   TVector3 xv = xhv + dlen * lr * ed; 
#else
   TVector3 xv = xhv + dlen * lr * ed; 
   xv -= 0.5 * (dlen*dlen/rw) * ns.Unit();
#endif
#else
   TVector3 xv = xhv + dlen * lr * ed; 
#endif

   return xv;
}

TVector3 EXMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   // Caution !!!!!
   // The resultant Xv is not exactly on the measurement surface!
   // Use this with care.

   const EXHit &ht = *dynamic_cast<const EXHit *> (&vht);

   Int_t m = ht.GetDimension();
   TKalMatrix mv(m,1);
   for (Int_t i=0; i<m; i++) mv(i,0)  = ht.GetX(i);
   return MvToXv(mv, ht.GetLR(), ht.GetCellNo());
}

void EXMeasLayer::CalcDhDa(const TVTrackHit &vht,
                           const TVector3   &xxv,
                           const TKalMatrix &dxphiada,
                                 TKalMatrix &H)  const
{
   const EXHit &ht = *dynamic_cast<const EXHit *> (&vht);

   // Calculate
   //    H = (@h/@a) = (@d/@a, @z/@a)^t
   // where
   //        h(a) = (d, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //
   // -----------
   // (1) (@d/@a)
   // -----------

   TKalMatrix xx(xxv);
   TKalMatrix xw(ht.GetWireEnd());
   TKalMatrix xxw = xx - xw;
   TKalMatrix xxwt(TKalMatrix::kTransposed,xxw);

   TVector3   xwdv = ht.GetWireDir();
   TKalMatrix wdir(xwdv);
   TKalMatrix wdirt(TKalMatrix::kTransposed,wdir);
   TKalMatrix unit(3,3);
   unit.UnitMatrix();

   TKalMatrix ddda  = xxwt * (unit - wdir * wdirt) * dxphiada;
   TVector3   xxwv  = TKalMatrix::ToThreeVec(xxw);
   Double_t   dlen  = (xxwv - xxwv.Dot(xwdv) * xwdv).Mag();
   ddda *= 1/dlen;
   Int_t lr = ht.GetLR();

   // -----------
   // (2) (@z/@a)
   // -----------

   TKalMatrix dzda = wdirt * dxphiada;
   dzda *= xwdv.Z();
   
   // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
   Int_t sdim = H.GetNcols();
   Int_t hdim = TMath::Max(5,sdim-1);

   for (Int_t i=0; i<hdim; i++) {
      H(0,i) = ddda(0,i) * lr;
      H(1,i) = dzda(0,i);
   }
   if (sdim == 6) {
#ifndef __TPC__
      H(0,sdim-1) = ht.GetVdrift() * lr;
      H(1,sdim-1) = 0.;
#else
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = ht.GetVdrift();
#endif
   }
}


void EXMeasLayer::ProcessHit(const TVector3 &xx, TObjArray &hits)
{
   Int_t lr     = 0;
   Int_t cellno = -99999999;
   TKalMatrix h = XvToMv(xx, lr, cellno);
   Double_t   x = h(0, 0);
   Double_t   z = h(1, 0);

   x += fVdrift * EXEventGen::GetT0() * lr;

   x += gRandom->Gaus(0.,fDx);
   z += gRandom->Gaus(0.,fDz);

   Double_t meas [2];
   Double_t dmeas[2];
   meas [0] = x;
   meas [1] = z;
   dmeas[0] = fDx;
   dmeas[1] = fDz;

   Double_t b = EXKalDetector::GetBfield(xx);
   hits.Add(new EXHit(*this, meas, dmeas, lr, cellno, fVdrift, xx, b));
}
