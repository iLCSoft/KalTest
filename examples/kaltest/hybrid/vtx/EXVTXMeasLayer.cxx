//*************************************************************************
//* ===================
//*  EXVTXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXVTXHit.
//* (Requires)
//* (Provides)
//*     class EXVTXMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//

#include "EXVTXMeasLayer.h"
#include "EXVTXHit.h"

ClassImp(EXVTXMeasLayer)
                                                                                
EXVTXMeasLayer::EXVTXMeasLayer(TMaterial &min,
                               TMaterial &mout,
                               Double_t   r0,
                               Double_t   lhalf,
                               Double_t   sigmax,
                               Double_t   sigmaz,
                               Bool_t     type)
              : EXVMeasLayer(min, mout, r0, lhalf, type),
                fSigmaX(sigmax),
                fSigmaZ(sigmaz)
{
}

EXVTXMeasLayer::~EXVTXMeasLayer()
{
}

TKalMatrix EXVTXMeasLayer::XvToMv(const TVector3 &xv) const
{
   // Calculate hit coordinate information:
   //	mv(0,0) = r * phi 
   //     (1,0) = z

   TKalMatrix mv(kMdim,1);
   mv(0,0)  = GetR() * TMath::ATan2(xv.Y(), xv.X());
   mv(1,0)  = xv.Z();
   return mv;
}

TKalMatrix EXVTXMeasLayer::XvToMv(const TVTrackHit &,
                               const TVector3   &xv) const
{
   return XvToMv(xv);
}

TVector3 EXVTXMeasLayer::HitToXv(const TVTrackHit &vht) const
{
   const EXVTXHit &ht = dynamic_cast<const EXVTXHit &>(vht);

   Double_t phi = ht(0,0) / GetR();
   Double_t z   = ht(1,0);
   Double_t x   = GetR() * TMath::Cos(phi);
   Double_t y   = GetR() * TMath::Sin(phi);

   return TVector3(x,y,z);
}

void EXVTXMeasLayer::CalcDhDa(const TVTrackHit &vht,
                           const TVector3   &xxv,
                           const TKalMatrix &dxphiada,
                                 TKalMatrix &H)  const
{
   // Calculate
   //    H = (@h/@a) = (@phi/@a, @z/@a)^t
   // where
   //        h(a) = (phi, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   Int_t sdim = H.GetNcols();
   Int_t hdim = TMath::Max(5,sdim-1);

   Double_t xv = xxv.X();
   Double_t yv = xxv.Y();
   Double_t xxyy = xv * xv + yv * yv;

   // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
   for (Int_t i=0; i<hdim; i++) {
      H(0,i) = - (yv / xxyy) * dxphiada(0,i) 
               + (xv / xxyy) * dxphiada(1,i);
      H(0,i) *= GetR();
      H(1,i) =  dxphiada(2,i);
   }
   if (sdim == 6) {
      H(0,sdim-1) = 0.;
      H(1,sdim-1) = 0.;
   }
}
